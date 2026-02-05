
# Import AiiDA-related packages
import aiida
from aiida.engine import submit
from aiida.orm import Str, StructureData, load_code

# Import helper packages
import numpy as np
import os, sys
from ase.io import read
from glob import iglob
from copy import deepcopy

# Import PyCCE work chain, to be submitted as jobs
from Chain_PyCCE import Chain_PyCCE

### -----------------------------------------------------
### This example submission script assumes that the AiiDA
### framework has been set up properly on your machine.
### 
### AiiDA can be installed by following the instructions here: 
###   https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/index.html
### -----------------------------------------------------

# Load profile
aiida.load_profile()

# Load code (aiida-pycce)
code_pycce = load_code(label="pycce2d@midway")

# Define job submission details
# See CalcJob documentation for more details (https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/calculations/usage.html#options)
custom_metadata = {
    "options": {
        "account": "pi-gagalli",
        "queue_name": "gagalli-csl2",
        #"qos": "gagalli-debug",
        "max_wallclock_seconds": 24 * 60 * 60,
        "import_sys_environment": False,
        "max_memory_kb": 192000000,
        "resources": {
            "num_machines": 1
        }
    }
}

# Define basic input parameters
parameters_pycce = {        
    "r_bath": 120,
    "r_dipole": 30,
    "order": 2,
    "mag_field": 50000,
    "pulses": 1,
    "mintime": 0,
    "maxtime": 50,
    "time_npoints": 1001,
}

# Gather crystal structures
structures = {}
for ciffile in iglob("Structures/*.cif"):
    structure = StructureData(ase=read(ciffile))  ## Important: Ensure that the structure data type is the AiiDA StructureData object
    name = ciffile.split("/")[-1].split(".cif")[0]
    structures[name] = structure


# Submit calculations (1 job per crystal structure)
for cif_code in structures.keys():

    # Gather inputs
    inputs = {
        "code_pycce": code_pycce,
        "calc_params_pycce": parameters_pycce,
        "custom_metadata": custom_metadata,
        "structure": structures[cif_code],
        "label": cif_code,
    }

    # Submit WorkChain
    chain = submit(Chain_PyCCE, **inputs)

    print(f"Submitted: {chain}")

