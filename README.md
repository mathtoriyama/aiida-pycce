Welcome to AiiDA-PyCCE!
=======================

AiiDA-PyCCE is a plugin for running the [PyCCE](https://pycce.readthedocs.io/en/latest/) code in a high-throughput manner using the [AiiDA](https://www.aiida.net/) framework.

If you use AiiDA-PyCCE, please remember to cite:

- *Strategies to search for two-dimensional materials with long spin qubit coherence time*
[M.Y. Toriyama, et al., npj 2D Mater. Appl., 9, 108 (2025)](https://www.nature.com/articles/s41699-025-00623-8)

Please also remember to cite:

- *AiiDA 1.0, a scalable computational infrastructure for automated reproducible workflows and data provenance*
    [S.P. Huber, et al., Sci. Data, 7, 300 (2020)](https://www.nature.com/articles/s41597-020-00638-4)

- *Workflows in AiiDA: Engineering a high-throughput, event-based engine for robust and modular computational workflows*
    [M. Uhrin, et al., Comput. Mater. Sci., 187, 110086 (2021)](https://www.sciencedirect.com/science/article/pii/S0927025620305772?via%3Dihub)


Installation
============
- Download all source files
    - `git clone https://github.com/mathtoriyama/aiida-pycce`
- Navigate to the appropriate directory
    - `cd aiida-pycce`
- Install
    - `pip install .`


Code Setup
==========
Please create the code using the `PyCCE_Wrapper_2D.py` script in the `code` directory. An example configuration file `code_setup_pycce2d.yml` is also provided. Please ensure that it's tailored to your specific system. The code can be created using the command:

    verdi code create core.code.installed --config code_setup_pycce2d.yml

Please also remember to update the first line of `PyCCE_Wrapper_2D.py`.

For more information, check out the [online guide](https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/run_codes.html#how-to-create-a-code) on how to create external codes on AiiDA.


Usage
=====
Check out the `HT_2D_Hosts.py` script in the `examples` folder.


Acknowledgements
================
This work was supported by the Midwest Integrated Center for Computational Materials (MICCoM). MICCoM is part of the Computational Materials Sciences Program funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Materials Sciences, and Engineering Division through the Argonne National Laboratory, under contract No. DEAC02-06CH11357.


