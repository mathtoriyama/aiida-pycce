Welcome to AiiDA-PyCCE!
=======================

AiiDA-PyCCE is a plugin for running the [PyCCE code](https://pycce.readthedocs.io/en/latest/) in a high-throughput manner using the [AiiDA](https://www.aiida.net/) framework.

If you use AiiDA-PyCCE, please remember to cite:

- *Strategies to search for two-dimensional materials with long spin qubit coherence time*
[Add arxiv link]()

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
Please create the code using the `PyCCE\_Wrapper\_2D.py` script. An example configuration file is provided in the `code` directory. The code can be created using the command:
`verdi code create core.code.installed --config code_setup_pycce2d.yml`


Getting Started
===============



Acknowledgements
================
Add later.


