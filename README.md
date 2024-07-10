# VIC-Borg-MPI
Calibration of the VIC model using the MPI version of the Borg multi-objective optimization algorithm

# Enter the 'example' directory
0. go to https://github.com/JinfengM/VIC-Borg/tree/main/example, cd example directory and download example.tar.gzaa and example.tar.gzab
   To ensure that MPICH 3.2 has been installed

# Unzip the example file 'example.tar.gz*' to obtain the 'run_lh' directory.
1. cat example.tar.gz* | tar -xzv

2. tar -zxvf example_files.tar.gz

# Create a directory at '/home/VIC'
3. mkdir /home/VIC

# To ensure that the '/home/VIC' directory exists
4.

# To copy the run_lh directory to the /home/VIC/ directory
5.cp -r run_lh /home/VIC/

# To ensure the '/home/VIC/run_lh' exists
6. 

# Download VIC-Borg-MPI project to the local directory '/home/VIC' 
7. cd /home/VIC/VIC-Borg-MPI/routMPI
# compile streamflow routing module, copy routMPI.so to /home/VIC/VIC-Borg-MPI/
8. make
   cp routMPI.so ../
# Add the current directory to the LD_LIBRARY_PATH environment variable
9. export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
# compile BorgMS at /home/VIC/VIC-Borg-MPI/
10. cd ..
make 
# ExecuteVIC_BORG_MPI.X
mpiexec -n 9 ./VIC_BORG_MPI.X -g /home/VIC/run_lh/chanliu_input.txt

# Notice
VIC model's source code is from https://github.com/UW-Hydro/VIC, users can access the source code here.
Borg algorithm's source code is from http://borgmoea.org/, users are required to complete the Google form to request access to the source code.
We thank both the developers of the VIC model and Borg algorithm for their great contribution.

