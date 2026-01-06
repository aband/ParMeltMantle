A serial simulator for Partially Molten Mantle

Packages needed:
1. PETSC
2. LAPACK (with LAPACKE)
link by providing route to .pc file

Run from ./driver
1. cd to the case desired Icase ./driver/full, /driver/half, etc).
2. build using build.sh
3. run as fast.sh (or test.sh)
4. output is in ./build (*.dat)

Case setup
1. flow.cpp stores information regarding the Darcy-Stokes solver
2. porosity.cpp stores information regarding how to compute porosity (with phase/ fixed function)

Plotting:
1. Matlab plotting files in ./driver (*.m) to process *.dat files
   e.g., plot_edge_velocity(40,40,1,'half','edgevel_stokes')
   Arguments: sizeX,sizeY,timeStepNumber,caseName,variableName

2. Python plotting files in ./driver (*.py) to process *.dat files
   Plotting settings are stored in myplot.py file.
   Plot by running commend: python3 run.py in terminal.
