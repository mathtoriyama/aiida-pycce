
import numpy as np
import os, sys
from glob import iglob
from scipy.optimize import curve_fit

from aiida.engine import WorkChain, ToContext, calcfunction, append_, if_, while_
from aiida.orm import Code, Str, Float, Bool, StructureData, load_group
from aiida.plugins import CalculationFactory, DataFactory

"""
Notes on updates:
    (1) Data management with temporary folder, creating/deleting folder
    (2) Monitor coherence function tail
    (3) Fix: label --> label.value
    (4) Stop calculation after 10 iterations
    (5) Heteronuclear decoupling
"""

Dict = DataFactory("core.dict")
PycceCalculation = CalculationFactory("aiida.pycce")

@calcfunction
def T2_Fit(data):
    
    # Settings for fitting
    p0 = [2, 1]
    r = 20
    maxfev = 5000

    # Data to fit coherence time
    # Note: the input "data" is a list
    data = np.array(data)
    x = data[:,0]
    y = data[:,1]

    # Hahn-echo coherence function to fit to
    def _Lfunc(t, t2, n):
        return np.exp(-(t / t2) ** n) * y[0]
    
    # Get envelope of coherence function by getting the max. value in every 20 bins
    npoints = x.size // r
    y_max = np.empty(npoints, dtype=np.float64)
    x_max = np.empty(npoints, dtype=np.float64)
    for i in range(npoints):
        ind = np.argmax(y[r * i:r * i + r])
        x_max[i] = x[r * i + ind]
        y_max[i] = y[r * i + ind]
    
    # Ensure that all coherence function values are less than 1
    x_max = np.delete(x_max, y_max > 1)
    y_max = np.delete(y_max, y_max > 1)

    # Fit coherence function
    param, cov = curve_fit(_Lfunc, x_max, y_max, p0=p0, maxfev=maxfev)

    return Dict(dict={"t2": param[0], "n": param[1]})


class Chain_PyCCE(WorkChain):

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Specification inputs
        spec.input("code_pycce", valid_type=Code, required=True)
        spec.input("calc_params_pycce", valid_type=Dict, required=True)
        spec.input("custom_metadata", valid_type=Dict, required=True)
        spec.input("structure", valid_type=StructureData, required=True)
        spec.input("label", valid_type=Str, required=True)
        spec.input("heteronuclear_decoupling", valid_type=Bool, required=False, default=Bool(False))
        
        spec.input("cf_tol", valid_type=Float, required=False, default=lambda: Float(0.01))
        spec.input("cf_tail_max", valid_type=Float, required=False, default=lambda: Float(0.01))
        spec.input("cf_min", valid_type=Float, required=False, default=lambda: Float(0.000001))  # Threshold for what's considered `0` for the coherence function

        # Workflow procedure / schedule
        spec.outline(
            cls.Initialize,
            while_(cls.PyCCE_Status)(
                if_(cls.Heteronuclear_Decoupling_Protocol)(
                    cls.Run_HT_PyCCE_byElement,
                ).else_(
                    cls.Run_HT_PyCCE,
                ),
                cls.Check_PyCCE_Calc,
                cls.Check_CoherenceFunction_Convergence,
            ),
            cls.Get_T2_Results,
            cls.Clean_Temp_Folder,
        )
        
        spec.output("results", valid_type=Dict)
        spec.exit_code(401, "ERROR_WALLTIME", message="The walltime ran out.")
        spec.exit_code(402, "ERROR_CANNOT_FIND_COHERENCE", message="Coherence_Func_Avg.txt cannot be found.")


    def Get_Inputs(self):
        return self.inputs


    def Initialize(self):
        
        # Create temporary folder to store data files (if necessary)
        self.ctx.tmp_folder = "/home/toriyama/.tmp"
        try:
            os.mkdir(self.ctx.tmp_folder)
        except:
            pass
        
        # Initialize parameters
        self.ctx.is_converged = False
        self.ctx.iteration = 1
        try:
            self.ctx.maxtime = self.inputs.calc_params_pycce.get_dict()["maxtime"]
        except:
            self.ctx.maxtime = 50.


    def PyCCE_Status(self):
        return not self.ctx.is_converged


    def Heteronuclear_Decoupling_Protocol(self):
        return self.inputs.heteronuclear_decoupling


    ##################################################################
    ########################## Utility ###############################
    ##################################################################
    
    def Get_Builder_PyCCE(
        self, 
        structure: StructureData, 
        calc_name: Str,
    ):

        # Initialize builder object
        builder = PycceCalculation.get_builder()
        builder.code = self.inputs.code_pycce
        
        # Structure (if decoupling, only structure with single element type)
        builder.structure = structure

        # Input parameters for PyCCE
        calc_params = self.inputs.calc_params_pycce.get_dict()
        calc_params["name"] = calc_name  # Name of calculation, should be descriptive (e.g. chemical formula)
        calc_params["maxtime"] = self.ctx.maxtime  # Adjust maximum time
        builder.parameters = Dict(dict=calc_params)

        # Job submission details
        builder.metadata = self.inputs.custom_metadata.get_dict()

        return builder


    def Check_Calcjob_Node(
        self,
        calcjob_node,
    ):

        # Check to see if calculation ran out of walltime
        if "DUE TO TIME LIMIT" in calcjob_node.get_scheduler_stderr():
            return self.exit_codes.ERROR_WALLTIME
        
        # Check if coherence function files can be found
        if "Coherence_Func.txt" not in calcjob_node.outputs.remote_folder.listdir():
            return self.exit_codes.ERROR_CANNOT_FIND_COHERENCE
        if "Coherence_Func_Avg.txt" not in calcjob_node.outputs.remote_folder.listdir():
            return self.exit_codes.ERROR_CANNOT_FIND_COHERENCE


    # Function for obtaining coherence function data (averaged over all processes)
    def Get_Coherence_Data(
        self,
        calcjob_node_label: Str,
    ):
        
        if self.inputs.heteronuclear_decoupling:

            # Initialize coherence function data
            coherence_times = None
            coherence_func = None

            # Loop through each element
            elements = self.inputs.structure.get_symbols_set()
            for element in elements:

                # Get coherence function, averaged over all processes, for element
                calcjob_node_label_element = calcjob_node_label+"_"+element
                calcjob_node_element = self.ctx[calcjob_node_label_element]
                remote_element = calcjob_node_element.outputs.remote_folder
                tmp_filename_element = self.ctx.tmp_folder + "/Coherence_Func_Avg_"+calcjob_node_label_element+"_"+calcjob_node_element.uuid+".txt"
                if not os.path.exists(tmp_filename_element):
                    remote_element.getfile("Coherence_Func_Avg.txt", tmp_filename_element)  # Copy to temp folder
                coherence_data_avg_element = np.loadtxt(tmp_filename_element)

                # Get times for coherence function
                if coherence_times is None:
                    coherence_times = coherence_data_avg_element[:,0]

                # Either initialize or multiply coherence function (following the principle of heteronuclear decoupling)
                if coherence_func is None:
                    coherence_func = coherence_data_avg_element[:,1]
                else:
                    coherence_func *= coherence_data_avg_element[:,1]

            # Gather data (time, coherence function)
            coherence_data_avg = np.vstack((coherence_times, coherence_func)).T

        else:

            # Get coherence function, averaged over all processes
            calcjob_node = self.ctx[calcjob_node_label]
            remote = calcjob_node.outputs.remote_folder
            tmp_filename = self.ctx.tmp_folder + "/Coherence_Func_Avg_"+calcjob_node_label+"_"+calcjob_node.uuid+".txt"
            if not os.path.exists(tmp_filename):
                remote.getfile("Coherence_Func_Avg.txt", tmp_filename)  # Copy to temp folder
            coherence_data_avg = np.loadtxt(tmp_filename)
        
        return coherence_data_avg


    ##################################################################
    ########################## Workflow ##############################
    ##################################################################
    
    def Run_HT_PyCCE(self):

        # Builder for PyCCE calculation
        calc_name = self.inputs.structure.get_formula(mode="hill_compact")
        builder = self.Get_Builder_PyCCE(structure=self.inputs.structure, calc_name=calc_name)

        # Submit calculation
        calcjob_node = self.submit(PycceCalculation, **builder)
        pycce_calculations = {}
        pycce_calculations[str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)] = calcjob_node
        
        # Ask the workflow to continue when the results are ready and store them in the context
        return ToContext(**pycce_calculations)


    def Run_HT_PyCCE_byElement(self):

        # Initialize dictionary to track PyCCE calculation of each element
        pycce_calculations = {}

        # Loop through each element
        elements = self.inputs.structure.get_symbols_set()
        for element in elements:

            # Get structure with only the chosen element
            structure_ase = self.inputs.structure.get_ase()
            atoms_to_delete = []
            for atom_index, atom in enumerate(structure_ase):
                if atom.symbol != element:
                    atoms_to_delete.append(atom_index)
            del structure_ase[atoms_to_delete]
            new_structure = StructureData(ase=structure_ase)

            # Builder for PyCCE calculation (one for each element)
            calc_name = self.inputs.structure.get_formula(mode="hill_compact")+"_"+element
            builder_by_element = self.Get_Builder_PyCCE(structure=new_structure, calc_name=calc_name)

            # Submit calculation
            calcjob_node_by_element = self.submit(PycceCalculation, **builder_by_element)
            pycce_calculations[str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)+"_"+element] = calcjob_node_by_element

        # Ask the workflow to continue when the results are ready and store them in the context
        return ToContext(**pycce_calculations)


    def Check_PyCCE_Calc(self):

        if self.inputs.heteronuclear_decoupling:
            elements = self.inputs.structure.get_symbols_set()
            for element in elements:
                calcjob_node_by_element = self.ctx[str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)+"_"+element]
                self.Check_Calcjob_Node(calcjob_node_by_element)
        else:
            calcjob_node = self.ctx[str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)]            
            self.Check_Calcjob_Node(calcjob_node)


    def Check_CoherenceFunction_Convergence(self):

        # Get coherence function
        calcjob_node_label = str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)
        coherence_data_avg = self.Get_Coherence_Data(calcjob_node_label=calcjob_node_label)
        coherence_func = coherence_data_avg[:,1]

        # Average L by bins, to get smoother function
        average_every = 10   # Cannot be zero
        coherence_func_averaged = np.array([np.mean(coherence_func[i:i+average_every]) for i in range(0, len(coherence_func), average_every)])
        coherence_func_averaged_last10 = coherence_func_averaged[int(len(coherence_func_averaged)*9/10.):]

        # Check if the coherence function tail has leveled off
        if (np.std(coherence_func_averaged_last10) > self.inputs.cf_tol) or np.all(coherence_func_averaged_last10 > self.inputs.cf_tail_max):
            # If not, extend the maxtime
            self.ctx.maxtime += 10
            self.ctx.iteration += 1
        # Check if more than a tenth of the coherence function is > 0
        elif (coherence_func_averaged > self.inputs.cf_min).sum() < int(len(coherence_func_averaged)/10):
            # If not, set the new maxtime to the time where the coherence function hits cf_min,
            #   unless L(t=0)=0 in which case set it to first index
            new_maxtime_index = np.argmin(coherence_func_averaged[coherence_func_averaged > self.inputs.cf_min]) * average_every
            if new_maxtime_index == 0:
                new_maxtime_index = average_every
            new_maxtime = coherence_data_avg[new_maxtime_index][0]
            self.ctx.maxtime = new_maxtime
            self.ctx.iteration += 1
        else:
            self.ctx.is_converged = True

        # Check if 10 iterations have been performed
        if self.ctx.iteration > 10:
            self.ctx.iteration -= 1 # Reset iteration
            self.ctx.is_converged = True


    # Retrieve, calculate, and record coherence time T2 and stretching exponent n
    def Get_T2_Results(self):

        calcjob_node_label = str(self.inputs.label.value)+"_pycce_iter"+str(self.ctx.iteration)
        coherence_data = self.Get_Coherence_Data(calcjob_node_label=calcjob_node_label)
        data = T2_Fit(coherence_data.tolist())

        self.out("results", data)


    # Delete data files from temporary folder
    def Clean_Temp_Folder(self):
        
        if self.inputs.heteronuclear_decoupling:

            # Loop through each element
            elements = self.inputs.structure.get_symbols_set()
            for element in elements:

                # Loop through each iteration
                for i in range(1, self.ctx.iteration+1):
                    calcjob_node_label = str(self.inputs.label.value)+"_pycce_iter"+str(i)+"_"+element
                    calcjob_node = self.ctx[calcjob_node_label]
                    tmp_filename = self.ctx.tmp_folder + "/Coherence_Func_Avg_"+calcjob_node_label+"_"+calcjob_node.uuid+".txt"
                    os.system("rm "+tmp_filename)

        else:

            # Loop through each iteration
            for i in range(1, self.ctx.iteration+1):
                calcjob_node_label = str(self.inputs.label.value)+"_pycce_iter"+str(i)
                calcjob_node = self.ctx[calcjob_node_label]
                tmp_filename = self.ctx.tmp_folder + "/Coherence_Func_Avg_"+calcjob_node_label+"_"+calcjob_node.uuid+".txt"
                os.system("rm "+tmp_filename)


