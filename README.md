1. Building (with auto-download the source codes for LAMMPS):

* git clone https://github.com/ProtKsen/new_SSAGES.git
* cd new_SSAGES/
* mkdir build/
* cd build/
* cmake .. -DLAMMPS=YES
* to '~/build/CMakeFiles/ssages.dir/link.txt' add key '-lintlc'
* make
* ls -l ~/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/liblammps_mpi.so
* rm ~/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/liblammps_mpi.so
* make

2. Description of changes

