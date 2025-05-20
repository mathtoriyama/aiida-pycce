
import numpy as np
from scipy.optimize import curve_fit

from aiida.engine import ExitCode, calcfunction
from aiida.orm import SinglefileData, Dict
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

PyCCECalculation = CalculationFactory('aiida.pycce')

class PyCCEParser(Parser):

    def __init__(self, node):
        super(PyCCEParser, self).__init__(node)

    def parse(self, **kwargs):
        """Parse outputs, store results in database.

        :returns: non-zero exit code, if parsing fails
        """

        output_filename = self.node.get_option('output_filename')

        # Check to make sure pycce.out was created
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [output_filename]
        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error(f"Found files '{files_retrieved}', expected to find '{files_expected}'")
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # Add output file to outgoing nodes
        self.logger.info(f"Parsing '{output_filename}'")
        with self.retrieved.open(output_filename, 'rb') as handle:
            output_node = SinglefileData(file=handle)
        self.out("output_file", output_node)

        # Add calculated coherence function to outgoing nodes
        with output_node.as_path() as filepath:
            coherence_data = []
            
            # Get number of processes
            for line in open(filepath).readlines():
                if "Number of processes:" in line:
                    n_procs = int(line.split()[-1])
                    break
            try:
                n_procs
            except:
                return ExitCode(3008)
            
            # Get data from relevant lines in pycce.out
            for line in open(filepath).readlines():
                if self.Check_Line(line, n_procs+1):
                    coherence_data.append([float(x) for x in line.split()])
            coherence_data = np.asarray(coherence_data)
        
        coherence_data_times = coherence_data[:,0]
        coherence_data_avg = np.mean( coherence_data[:,1:], axis=1 )
        avg_coherence_data = np.vstack(( coherence_data_times, coherence_data_avg )).T
        
        all_data = Dict()
        #all_data["Coherence_Data_Raw"] = coherence_data
        all_data["Coherence_Data_Averaged"] = avg_coherence_data
        self.out("output_data", all_data)

        return ExitCode(0)


    def Check_Line(self, line, n_elems):
        
        line_elems = line.split()
        
        # Check if line is empty
        if len(line_elems) == 0:
            return False
        
        # Check if line has correct number of elements
        if len(line_elems) != n_elems:
            return False
        
        # Check if any element in line is not a float
        for x in line_elems:
            try:
                float(x)
            except:
                return False
        
        return True

