
import os
import numbers
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, RemoteData, Dict, StructureData, Str


def conv_to_fortran(val, quote_strings=True):
    """Convert a python value to a format suited for fortran input.

    :param val: the value to be read and converted to a Fortran-friendly string.
    """
    import numpy as np

    # Note that bool should come before integer, because a boolean matches also isinstance(..., int)
    if isinstance(val, (bool, np.bool_)):
        if val:
            val_str = '.true.'
        else:
            val_str = '.false.'
    elif isinstance(val, numbers.Integral):
        val_str = f'{val:d}'
    elif isinstance(val, numbers.Real):
        #val_str = f'{val:18.10e}'.replace('e', 'd')
        val_str = f'{val:18.10}'
    elif isinstance(val, str):
        if quote_strings:
            val_str = f"'{val!s}'"
        else:
            val_str = f'{val!s}'
    else:
        raise ValueError(
            f"Invalid value '{val}' of type '{type(val)}' passed, accepts only bools, ints, floats and strings"
        )

    return val_str



class PyCCECalculation(CalcJob):
    """AiiDA calculation plugin wrapping the diff executable."""

    _DEFAULT_INPUT_FILE = "pycce.in"
    _DEFAULT_OUTPUT_FILE = "pycce.out"
    _CRASH_FILE = "CRASH"

    @classmethod
    def define(cls, spec):
        
        """Define inputs and outputs of the calculation."""
        super().define(spec)

        spec.input('metadata.options.input_filename', valid_type = str, default = cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type = str, default = cls._DEFAULT_OUTPUT_FILE)

        spec.input('structure', valid_type = StructureData, required = False)
        spec.input('parameters', valid_type = Dict,
            help = 'Input parameters for PyCCE.',
            required = True,
        )

        # Set parser
        spec.inputs["metadata"]["options"]["parser_name"].default = "aiida.pycce"

        # Define outputs of node
        spec.output('output_file', valid_type = SinglefileData)
        spec.output('output_data', valid_type = Dict)



    def max_retrieve_list(self):
        retrieve_list = [
            self.metadata.options.output_filename,
            "Coherence_Func.txt",
            "Coherence_Func_Avg.txt",
        ]
        return retrieve_list



    def prepare_for_submission(self, folder):
        """Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files needed by
            the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        self.verify_inputs()

        input_filecontent = self._generate_PyCCEinputdata(self.inputs.parameters)
        with folder.open(self.metadata.options.input_filename, "w") as handle:
            handle.write(input_filecontent)

        if self.inputs.structure:
            structure_str = self.inputs.structure._prepare_cif()[0].decode("utf-8")
            with folder.open("structure.cif", "w") as handle:
                handle.write(structure_str)
        
        # Gather information about code
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (['-in', self.metadata.options.input_filename])
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.code_pk = self.inputs.code.pk
        
        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.codes_info = [codeinfo]    # Information about code
        #calcinfo.retrieve_list = [self.metadata.options.output_filename, self._CRASH_FILE]  # Outputs
        calcinfo.retrieve_list = self.max_retrieve_list()  # Outputs
        calcinfo.local_copy_list = []

        return calcinfo



    @classmethod
    def _generate_PyCCEinputdata(cls, parameters):
        inputfile = ""
        for key, value in parameters.items():
            inputfile += f"  {key} = {conv_to_fortran(value)}\n"
        return inputfile



    def verify_inputs(self):

        parameters = self.inputs.parameters.get_dict()
        params_needed = ["r_bath", "r_dipole", "order", "mag_field", "pulses", "mintime", "maxtime"]
        for param in params_needed:
            if param not in parameters:
                raise ValueError(f"Need to have '{param}' tag in parameters.")

        params_needed_without_structure = ["pw_infile", "gipaw_efg_outfile"]
        for param in params_needed_without_structure:
            if (not self.inputs.structure) and (param not in parameters):
                raise ValueError(f"Need to have '{param}' tag in parameters, since structure is not provided.")



