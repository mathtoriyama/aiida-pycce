
import aiida
from aiida.engine import submit
from aiida.orm import Str, StructureData, load_code

import numpy as np
import os, sys
from ase.io import read
from glob import iglob
from copy import deepcopy

from Chain_PyCCE import Chain_PyCCE

aiida.load_profile()
code_pycce = load_code(label="pycce2d@midway")

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


structures = {}
for ciffile in iglob("Structures/*.cif"):
    structure = StructureData(ase=read(ciffile))
    name = ciffile.split("/")[-1].split(".cif")[0]
    print(name)
    structures[name] = structure
print(structures)


# Submit calculations
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

    print(f"Submitted {chain}")

