#!/path/to/pycce_env/bin/python

import os, sys
from mpi4py import MPI
import numpy as np
import pandas as pd
from ase.io import read
import pycce as pc
from pycce.bath import common_concentrations


def Set_Params(param_filename):
    
    params_available = ["pw_infile", "gipaw_efg_outfile", 
                        "r_bath", "r_dipole", "mag_field", "mintime", "maxtime", "central_spin", "sc_size",
                        "order", "pulses", "atom_index", "time_npoints",
                        "zdir",
                        "name",
                        ]

    params = {}
    for line in open(param_filename).readlines():

        if "=" not in line:
            sys.exit(f"Line `{line}` is not in the correct format.")

        x = [i.strip() for i in line.split("=")]
        
        if not (x[0] in params_available):
            sys.exit(f"Parameter `{x[0]}` is not legitimate.")
        
        if x[0] in ["pw_infile", "gipaw_efg_outfile"]:
            params[x[0]] = x[1].strip("'")
        elif x[0] in ["r_bath", "r_dipole", "mag_field", "mintime", "maxtime", "central_spin", "sc_size"]:
            params[x[0]] = float(x[1])
        elif x[0] in ["order", "pulses", "atom_index", "time_npoints"]:
            params[x[0]] = int(x[1])
        elif x[0] in ["zdir"]:
            if len(x) != 5:
                sys.exit(f"Please have parameter `{x[0]}` in the format `x y z`.")
            params[x[0]] = [float(x[2]), float(x[3]), float(x[4])]
    return params


def Calc_PyCCE( r_bath = 60,
                r_dipole = 10,
                order = 2,
                mag_field = 500,
                pulses = 1,
                atom_index = 0,
                mintime = 0,
                maxtime = 50,
                zdir = [0, 0, 1],
                central_spin = 1,
                sc_size = 200,
                time_npoints = 1001,
                pw_infile = None,
                gipaw_efg_outfile = None,
                ):

    # Set MPI settings
    seed = 1
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    conf = rank #+ start
    np.random.seed(seed + conf)
    np.set_printoptions(suppress=True, precision=5)

    # Print number of processes to pycce.out
    print("Number of processes: ", size)

    # PyCCE coherence function calculation setup
    calc_setup = {'r_bath': r_bath, 'r_dipole': r_dipole, 'order': order}
    calc_param = {'magnetic_field': np.array([0., 0., mag_field]), 'pulses': pulses}

    # Get structure, either from QE or structure.cif
    if pw_infile:
        struc = read(pw_infile, format="espresso-in")
    elif "structure.cif" in os.listdir(os.getcwd()):
        struc = read("structure.cif")
    else:
        return "Structure not found."

    # Initialize bath
    cell_pycce = pc.bath.BathCell.from_ase(struc)

    # Define vacancy defect
    center_atom = struc[atom_index].symbol
    center_pos = struc[atom_index].scaled_position
    center = cell_pycce.to_cartesian(center_pos)
    cell_pycce.zdir = zdir
    vc_uc = np.linalg.inv(cell_pycce.cell) @ center

    # Get isotopes and concentrations, add to bath object
    # Supposed to look like:
    #       cell_pycce.add_isotopes(('183W', 0.144), ('33S', 0.0076))
    elements = list(set(struc.get_chemical_symbols()))
    conc = []
    for element in elements:
        if element == "Ce": continue
        for isotope in common_concentrations[element].keys():
            conc.append( (isotope, common_concentrations[element][isotope]) )
    cell_pycce.add_isotopes(*conc)
    
    # Create supercell, create vacancy defect
    atoms_pycce = cell_pycce.gen_supercell(sc_size, remove=(center_atom, vc_uc))

    # Only one layer
    z_len = struc.cell.lengths()[2]
    atoms_pycce_onelayer = atoms_pycce[ np.abs(atoms_pycce["xyz"] - center[np.newaxis, :])[:, 2] < z_len/2 ]

    # Define states
    alpha = np.array([0, 0, 1])   # |0> state in Sz basis
    beta = np.array([0, 1, 0])    # |1> state in Sz basis

    # Start PyCCE simulator
    calc = pc.Simulator(spin=central_spin,
                        position=center,
                        alpha=alpha,
                        beta=beta,
                        bath=atoms_pycce_onelayer,
                        **calc_setup)

    # Read quadrupole interactions from GIPAW output
    # Supposed to look like:
    #       ba = pc.read_qe("DFT/pw.in", efg="DFT/efg.out", isotopes={"W": "183W", "S": "33S"})
    if pw_infile and gipaw_efg_outfile:
        for element in elements:
            if element == "Ce": continue
            for isotope in common_concentrations[element]:
                bath_qe = pc.read_qe(pwfile=pw_infile, efg=gipaw_efg_outfile, isotopes={element: isotope})
                calc.bath["Q"][calc.bath["N"] == isotope] = bath_qe[isotope]["Q"][0][np.newaxis, :, :]


    # Set time for coherence function
    time_space = np.linspace(mintime, maxtime, time_npoints)

    # Generate clusters and compute coherence function
    clusters = calc.generate_clusters(order=calc_setup['order'], r_dipole=calc_setup['r_dipole'])
    coherence_func = calc.compute(time_space, magnetic_field=calc_param['magnetic_field'], pulses=calc_param['pulses'], as_delay=False)
    
    # Problem: If no clusters are found (clusters == {}), then the coherence_func will return 1 (integer) instead of an array
    # Solution: Just return an array of 1's instead
    if (clusters == {}) and (not isinstance(coherence_func, np.ndarray)):
        coherence_func = np.ones(len(time_space))
    
    # Wait until all MPI processes finish running
    # Gather data from each MPI process once done
    coherence_func_all = comm.gather(np.abs(coherence_func), root=0)
    coherence_func_all = comm.bcast(coherence_func_all, root=0)
    coherence_func_all.insert(0, time_space)
    if rank == 0:
        print(coherence_func_all)
    coherence_func_all = np.asarray(coherence_func_all).T
    coherence_func_all = np.nan_to_num(coherence_func_all)

    # Save to files (but only by root process)
    if rank == 0:
        np.savetxt("Coherence_Func.txt", coherence_func_all)

        # Output average coherence function (averaged over all processes)
        coherence_func_avg = np.mean( coherence_func_all[:,1:], axis=1 )
        coherence_data_avg = np.vstack(( time_space, coherence_func_avg )).T
        np.savetxt("Coherence_Func_Avg.txt", coherence_data_avg)

        # Output to terminal (to .out file)
        for line in coherence_func_all:
            print(" ".join([str(i) for i in line]))



if __name__ == "__main__":
    param_filename = sys.argv[2]
    params = Set_Params(param_filename)
    Calc_PyCCE(**params)
