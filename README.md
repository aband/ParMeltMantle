A serial simulator for Partially Molten Mantle

Packages needed:
1. PETSC
2. LAPACK (with LAPACKE)

Run from ./driver
1. cd to the case desired Icase ./driver/full, /driver/half, etc).
2. build using build.sh
3. run as fast.sh (or test.sh)
4. output is in ./build (*.dat)

Matlab plotting files in ./driver (*.m) to process *.dat files
   e.g., plot_edge_velocity(40,40,1,'half','edgevel_stokes')
   Arguments: sizeX,sizeY,timeStepNumber,caseName,variableName

