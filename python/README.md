# Python notes for Renfro PhD Dissertation

- `my_anaconda.yml`: Export of Anaconda environment used to run simulations. Can be imported with `conda env create -f my_anaconda.yml`
- `plate_creator.py`: Main python script used for creating WARP3D input files (requires local Python module `preprocess.py` and usually requires a `*config.py` to define the model parameters)
- `*config.py`: definition of model parameters for an execution of `plate_runner.py` or `plate_creator.py`
- `plate_runner.sh`: Slurm job script used to run simulations (loads conda environment, runs `plate_runner.py`)
- `plate_runner.py`: Main python script used for simulations (requires local Python modules `postprocess.py`, `preprocess.py`, `section.py`, `solve.py`, `util.py`, and usually requires a `*config.py` to define the model parameters)
- `preprocess.py`, `postprocess.py`, `solve.py`, `util.py`: functions useful for creating WARP3D models, evaluating WARP3D output, calculating model objective functions, and manipulating semi-structured text files.
- `section.py`: Dr. Drang's module for section properties (http://leancrew.com/all-this/2018/01/python-module-for-section-properties/), used to calculate initial bending stress value
- `*wrp.inp`: WARP3D input files generated from `plate_creator.py` and FEACrack, specific to a particular loading configuration, geometry, and material
- `*msh.out`: Mesh data generated from FEACrack, specific to a particular geometry
- `bend.elt`, `tens.elt`: FEACrack input file, specific to a particular loading configuration
