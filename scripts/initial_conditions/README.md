# Initial Conditions

A set of jupyter notebooks have been created for generating initial conditions for NAARC.

The jupyter notebooks can be run using the pyic environment which can be installed following the instructions here (https://github.com/NOC-MSM/pyIC). Before running the notebooks run the exports in the file export_for_pyic.txt. 

The 3 notebooks are IC.ipynb, IC_2d_vars.ipynb and reassemble_IC_chunks.ipynb. IC.ipynb generates the 3D variables temperature and salinity. The size of the dataset means a large amount of memory is required for interpolation. To serialize the interpolation and reduce memory load, the domain is divided up into smaller chunks for processing. reassemble_IC_chunks.ipynb should be run after IC.ipynb to assemble the chunks of 3D initial conditions back into a single file. IC_2d_vars.ipynb is used to create initial conditions for ice variables, siconc and sithic. The output of IC_2d_vars.ipynb does not need to be reassembled.
