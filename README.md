# NGA_LeWI
Course Project for COS597E: Dynamic Load Balancing in a Massively Parallel Reacting Flow Solver

## Installing `NGA`
[Library Installation](library-installation)

[Code compilation](code-installation)

## Running `NGA`
[Generating a data file](generating-a-data-file)

[Code execution](code-execution)

[Running on `TIGER`](running-on-tiger)

[Examples](examples)

## Miscellaneous

[Converting Mechanisms for Finite Chem](finite-chem-mechanism)

## Installing `LeWI`

[Library Installation](library-installation)

## Basic usage of `LeWI`

[Basic usage](basic-usage)

## [Reference Mapping](reference-mapping)

## [Citations](citations)

----
# `NGA`
Next Generation Advanced Reacting Turbulence Solver (ARTS), a multi-dimensional, structured flame and turbulence code.

# Library Installation

[NGA](home) requires five libraries: `LAPACK`, `SUNDIALS`, `FFTW`, `HYPRE`, and some version of `MPI`. These libraries must be installed prior to [compiling NGA](code-installation). The `MPI` library only needs to be installed if it has not been on the machine in use; note that all Princeton and DOE clusters already have an `MPI` library installed.

*It is very important that the same compiler version (and `MPI` version) are used when installing these libraries and `NGA`.* There are many versions installed on the clusters, which can be viewed using the `module avail` command and loaded using the `module load` command. The recommended versions are as follows:

On `tigercpu` (Old Tiger) : Intel 19.1 (default)

      module load intel
      module load intel-mpi

On `tiger3` (New Tiger) : Intel 2024.2
```
module load intel/2024.2
module load intel-oneapi/2024.2
module load intel-mkl/2024.2
module load intel-mpi/oneapi/2021.13
```
Add the appropriate module load commands to `~/.bashrc` file so the appropriate modules are loaded every time one logs in to the cluster.

The libraries can either be installed through a single automated process, as described in the next section, or individually and manually, as described in subsequent sections.

## Automated Installation

For convenience, the Bash script `NGA-Library-Installer.sh` has been provided at the bottom of this section. It will automatically install the required libraries for `NGA` (including the serial and parallel versions of `HYPRE`, excluding `MPI`).  Copy the script to the home directory, and then enable executing the installation script with: 
```
    chmod 755 NGA-Library-Installer.sh
```
Then execute the installation:
```
    ./NGA-Library-Installer.sh tiger3
```
where `tiger3` (New Tiger) can be replaced by `tigercpu` (Old Tiger) or `<workstation>`, depending on where the libraries are installed. The script uses Intel version 2024.2 on `tiger3` and Intel version 19.1 on `tigercpu`. Note that if installing on a workstation, one may need to manually install a version of `MPI` (see below) prior to running the automated installation script for the other libraries.

This script updates the `.bashrc` file, so run the command
```
     source ~/.bashrc 
```
to make sure these changes are reflected in the working terminal. If the install script completed with no errors, the libraries are now installed, and move on to [installing NGA](code-installation). Verify installation by checking in `/home/$USER/opt/` directory, which should have subdirectories where `lapack`, `hypre`, `hypre` (with `openmp`), `fftw`, and `sundials` were installed. The specific files that should be within these directories are listed at the end of each section of the manual installation instructions.

----
The script can also be used to uninstall the libraries. To do this, add `uninstall` to the end of the command above. For example,
```
    ./NGA-Library-Installer.sh tigercpu uninstall
```
An old version of the installation bash script `NGA-Library-Installer-13.sh` has also been provided, but its use is deprecated.  It works with old versions of the compilers and libraries that may no longer be available. (`intel-13.0`, `lapack-3.5.0`, `fftw-3.3.4`, `sundials-2.6.1`, `hypre-2.10.0b`).

----
[NGA-Library-Installer.sh](/uploads/8dcbebdae64d40040e640dbe7e89866d/NGA-Library-Installer.sh)
[NGA-library-installer.sh](/uploads/bf8336fd9dbcd03727663c82ce01d3d2/NGA-library-installer.sh)
[NGA-Library-Installer-13.sh](/uploads/087b2a495cdf2c897e4f21afdd733529/NGA-Library-Installer-13.sh)

----
## Manual Installation

Detailed instructions for **manually installing** these libraries are given below.  The `MPI` library only needs to be installed if it has not been on the machine in useg; note that all Princeton and DOE clusters already have an `MPI` library installed.

The compiler flags in the instructions below assume the use of Intel compilers on `TIGER` or local [CTRFL](https://ctrfl.princeton.edu/) workstations.  These instructions will need to be modified for other systems and software.

Unless otherwise indicated, all commands should be performed as the normal user.  

On `TIGER`, load the latest versions of the compiler and `MPI` libraries using the `module load` command as described above. 

For the `SUNDIALS` installation, the latest version of Python 2.x (2.7.8 as of 08/17 on `TIGER`) is also needed.  On Old `TIGER`, execute `module load python/2.7`. On New `TIGER`, execute `module load anaconda`.

View the latest versions of the above through `module avail`.

Also, installation of `g++` through `yum install gcc-c++`  may be needed if the personal workstation does not have this compiler.

----
### `LAPACK`

`LAPACK` (Linear Algebra PACKage) is a library of various numerical linear algebra routines.  The latest version of `LAPACK` (3.7.1 as of 08/17) can be downloaded from http://www.netlib.org/lapack.

Download and unpack the archive by executing the following:

    cd /home/$USER/Downloads
    wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz
    tar -zxf v3.11.tar.gz
    cd lapack-3.11

In the new directory, copy the file `make.inc.example` to `make.inc`.  Open this file and change the Fortran compiler options to the following:

    FORTRAN  = ifort
    OPTS     = -O3 -xhost -ip
    NOOPT    = 
    LOADER   = ifort
    LOADOPTS = -O3 -xhost -ip
    ARCH      = xiar

In addition, the `TIMER` option will need to be changed.  Comment out the current setting, `TIMER = INT_ETIME`, and re-enable `TIMER = EXT_ETIME`, which is compatible with the `ifort` compiler.  There are additional options for a C compiler, but these do not need to be changed since `NGA` does not use the C interfaces to `LAPACK`.

`NGA` will need the BLAS (Basic Linear Algebra Subprograms) library for simple matrix/vector math such as multiplication, etc., and the `LAPACK` library.  These libraries must be named, `libblas.a` and `liblapack.a`, respectively.  Therefore, edit the following options in `make.inc`:

    BLASLIB      = ../../libblas.a
    LAPACKLIB    = liblapack.a

To speed up compiling, compile in parallel.  Open `Makefile` and edit the instructions under the target `blaslib` to be `$(MAKE) -C BLAS -j 4` to compile with 4 cores.  Similarly, edit the instructions under the target `lapacklib` to be `$(MAKE) -C SRC -j 4`.

The compilation of the two libraries is done in two steps in the following order:

1. To compile BLAS, execute
  
        make blaslib
  
2. To compile `LAPACK`, execute
  
        make lapacklib
  
To complete the installation of these libraries, copy `libblas.a` and `liblapack.a` to a convenient location through

    sudo cp libblas.a liblapack.a /opt/lapack/lib

or

    cp libblas.a liblapack.a /home/$USER/opt/lapack/lib

for installation on a local workstation/laptop (superuser access available) or on the `TIGER` cluster (no superuser access available), respectively.

----
### `FFTW`

`FFTW` (Fastest Fourier Transform in the West) is a library for computing discrete Fourier transforms.  The latest version of `FFTW` (3.3.10 as of 04/23) can be downloaded from http://www.fftw.org or via

    wget http://www.fftw.org/fftw-3.3.10.tar.gz

Download and unpack the archive.  The compilers and compiler flags are set with environment variables.  From the command line, execute the following commands:

    export CC=icc
    export CXX=icpc
    export F77=ifort
    export AR=xiar
    export LD=xild
    export CFLAGS="-O3 -xhost -ip"
    export CXXFLAGS="-O3 -xhost -ip"
    export FFLAGS="-O3 -xhost -ip"

To automatically generate a `Makefile` for `FFTW`, run the `configure` script.  The option `prefix` will dictate where the code is installed.  For example, depending on superuser access, execute either

    ./configure --prefix=/opt/fftw --enable-openmp

or

    ./configure --prefix=/home/$USER/opt/fftw --enable-openmp

To compile the `FFTW` library in parallel (on `TIGER` use `-j 8`), execute

    make -j 4

and, to install the `FFTW` library, execute

    sudo make install

or

    make install

depending on the availability of superuser access.

Upon successful installation, the files `fftw3.h` and `fftw3.f`  should appear in the directory `/home/$USER/opt/fftw/include`, and the files `libfftw3.a` and `libfftw3_omp.a` should appear in the directory `/home/$USER/opt/fftw/lib`.

----
### `SUNDIALS` 4.0.0

First, as usual, download and unpack the tarball from LLNL's github via the following commands:

    cd /home/$USER/Downloads
    wget https://github.com/LLNL/sundials/releases/download/v4.0.0/sundials-4.0.0.tar.gz
    tar -zxf sundials-4.0.0.tar.gz

Create a `sundials` directory along with a `builddir` subdirectory:

    mkdir -p /home/$USER/opt/sundials-4/builddir
    cd /home/$USER/opt/sundials-4/builddir

Then create bash variables to store the paths to the C and Fortran compilers. In the command-line, type

    ccpath="$(which icc)"
    fcpath="$(which ifort)"

Now use the CMake to generate the `Makefile`.  Execute the following:

    cmake -DCMAKE_INSTALL_PREFIX=/home/$USER/opt/sundials \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=$ccpath \
    -DCMAKE_C_FLAGS_RELEASE="-O3 -xhost -ip" \
    -DCMAKE_Fortran_COMPILER=$fcpath \
    -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -xhost -ip" \
    -DEXAMPLES_ENABLE=OFF \
    -DFCMIX_ENABLE=ON \
    /home/$USER/Downloads/sundials-4.0.0

Execute `make install -j` to compile in parallel and install the library.  If the installation directory requires superuser permissions for writing, the latter command must be executed as superuser.

Upon successful installation, the libraries `libsundials_cvode.a`, `libsundials_fcvode.a`, `libsundials_fnvecserial.a`, and `libsundials_nvecserial.a` should be in the `lib` sub-directory of the directory in which `SUNDIALS` was installed.

Lastly, add

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/$USER/opt/sundials/lib64

to `~/.bashrc`.  Don't forget to call `source ~/.bashrc` in the working terminal.

**Warning: `SUNDIALS` 4.0.0 drops the FCVDENSE interface, making master branch `NGA` installations not work. For now, the workaround is to comment out line 1135 in `src/combustion/finitechem.f90` (master branch) since it is no longer necessary.**

----
### `MPI` (not needed if running on a cluster)

Message Passing Interface (`MPI`) is a standard for libraries that allow code parallelization on distributed memory systems.  **Only install an `MPI` implementation if one is not already installed on the machine in use.  All clusters, including `TIGER`, already have an `MPI` implementation installed.**  The `MPI` compilers are required to compile `NGA` (and the `HYPRE` library in the next step), so include the line `module load openmpi/intel-16.0` in the `~/.bashrc` file.

Install an `MPI` library on [CTRFL](https://ctrfl.princeton.edu/) workstations if it is not already installed.  `MPI` is the recommended distribution, and the latest version (2.1.1 as of 08/17) can be downloaded from http://open-mpi.org.

Download and unpack the archive.  To automatically generate a `Makefile` for `MPI`, run the `configure` script.  The option `prefix` will dictate where the code is installed.  For example, depending on superuser access, execute either

    ./configure CC=icc CXX=icpc F77=ifort FC=ifort \
    CFLAGS="-O3 -xhost -ip" CXXFLAGS="-O3 -xhost -ip" \
    FFLAGS="-O3 -xhost -ip" FCFLAGS="-O3 -xhost -ip" \
    --prefix=/opt/openmpi

or

    ./configure CC=icc CXX=icpc F77=ifort FC=ifort \
    CFLAGS="-O3 -xhost -ip" CXXFLAGS="-O3 -xhost -ip" \
    FFLAGS="-O3 -xhost -ip" FCFLAGS="-O3 -xhost -ip" \
    --prefix=/home/$USER/opt/openmpi

To compile the `MPI` library, execute (on `TIGER` use `-j 8`):

    make -j 4 all

and, to install the `MPI` library, execute

    make install

If the installation directory requires superuser permissions for writing, the latter command must be executed as superuser.

Upon successful installation, the binaries `mpicc`, `mpicxx`, `mpif77`, `mpif90`, and `mpifort` should be in the `bin` sub-directory of the directory in which `MPI` was installed.

If installing `NGA` on a workstation machine or if a copy of `MPI` is already installed, add `MPI` to the user bash profile.  Open the hidden `~/.bashrc` file and add the following:

    export PATH=$PATH:/home/$USER/opt/openmpi/bin

The directory for the `bin` path may not appear exactly as above; it depends on where `MPI` was installed.

Additionally, add the libraries in a similar fashion as above, using:

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/$USER/opt/openmpi/lib

----
### `HYPRE`

`HYPRE` (Scalable Linear Solvers and Multigrid Methods) is a library for solving linear systems in parallel with various iterative methods (Krylov and multigrid). `HYPRE` can be acquired from their github page via `wget`:

    wget https://github.com/hypre-space/hypre/archive/refs/heads/master.zip

Download and unpack the archive.  Note that an `MPI` library must be installed prior to the installation of `HYPRE`.  Once the `MPI` library is installed, execute 

    source ~/.bashrc
 
from the home directory.  The compilers and compiler flags are set with environment variables.  From the command line, execute the following commands:

    unset AR
    unset LD
    export CC=mpicc
    export CXX=mpicc
    export FC=mpif90
    export CFLAGS="-O3 -xhost -ip"
    export CXXFLAGS="-O3 -xhost -ip"
    export FCFLAGS="-O3 -xhost -ip"

Note that the first two lines ensure that AR is not set to `xiar` and LD is not set to `xild`, otherwise an error `ar: illegal option -- H` will prevent compilation.  To automatically generate a *Makefile* for `HYPRE`, run the `configure` script while in the `src` directory.  The option `prefix` will dictate where the code is installed.  For example, depending on superuser access, execute either

    cd /home/$USER/Downloads/hypre-2.11.2/src
    ./configure --prefix=/opt/hypre

or

    cd /home/$USER/Downloads/hypre-2.11.2/src
    ./configure --prefix=/home/$USER/opt/hypre

To compile and install the `HYPRE` library, execute (on `TIGER` use `-j 8`):

    make install -j 4

If the installation directory is something that requires superuser permissions for writing, the latter command must be executed as superuser.

Upon successful installation, the file `libHYPRE.a` should be in the `lib` sub-directory of the directory in which `HYPRE` was installed.

For some problems on certain systems, when `NGA` is run in a hybrid parallel mode (`MPI` + `OpenMP`), higher performance can *'sometimes*' be achieved when `HYPRE` is allowed to use the `OpenMP` threads.  To allow the latter, the same instructions as above should be followed.  However, the compiler flags should be changed to

    export CFLAGS="-O3 -qopenmp -xhost -ip"
    export CXXFLAGS="-O3 -qopenmp -xhost -ip"
    export FCFLAGS="-O3 -qopenmp -xhost -ip"

and an additional option is required for the `configure` script:

    ./configure --with-openmp --prefix=/opt/hypre-openmp

or

    ./configure --with-openmp --prefix=/home/$USER/opt/hypre-openmp

As a practical note, it is good practice to keep two working installations of `HYPRE`, one without `OpenMP` and one with `OpenMP`.  This allows the fastest installation to be used for various problems without having to constantly recompile/reinstall `HYPRE`. `HYPRE` with `OpenMP` should also work for runs without `OpenMP` without errors/failure. Keeping both is just for 5-10% performance improvement, which may not be the case any more for new version of `HYPRE`. When installing the second `HYPRE`, remember to make clean first.



# Downloading `NGA`

To download a clear version of `NGA` from the server, execute the command

    git clone git@ctrfl-internal.Princeton.EDU:ctrfl_codes/NGA.git

To download a clear version of `NGA` from the server into a directory other than "NGA", execute the command

    git clone git@ctrfl-internal.Princeton.EDU:ctrfl_codes/NGA.git my_NGA

If `NGA` is already installed and a simple version update with the most recent changes from the server is needed, execute the command

    git pull

From `STAMPEDE` and `PERLMUTTER` (or other DOE clusters) the clone command needs to include the username:

    git clone http://<username>@ctrfl-internal.Princeton.EDU/ctrfl_codes/NGA.git

----
# Compiling `NGA`

## Update `Makefile.in`  

In the `NGA` `src` directory, open the file `Makefile.in`  The first three sections of this file should not need to be modified unless running `NGA` outside of Princeton.  Otherwise, a change to the paths to the libraries is needed.  Change the four entries `BLAS_DIR`, `LAPACK_DIR`, `HYPRE_DIR`, `FFTW_DIR`, and `SUNDIALS_DIR` to the directories in which these libraries were installed. Write out the full path to the directories, e.g., `/home/username/opt/fftw` and not `~/opt/fftw`, otherwise the compiler may fail.

Additional steps may be necessary depending on where `NGA` is being installed and what compiler version is used:

### New `TIGER` (tigercpu)
The changes in this section have already been applied to `Makefile.in` in the latest version of `NGA`, but are left for posterity until after fully transitioning to `tigercpu`.

Remark: Using `mpif90` or `mpifort` depends on which one is available, but they probably perform similarly when both are available.

The Skylake processors on the new `TIGER` support 512-bit vector extensions, which could make `NGA` run faster. To compile using these extensions, set the following two lines in `Makefile.in` as:
```
    LDFLAGS  = -qopenmp -xCORE-AVX512
    OPTFLAGS = -qopenmp -O3 -xCORE-AVX512
```

Remarks:

(1) When use Intel-16, `-O3` and `-qoverride-limits` could result in segmentation fault for some cases (e.g., DNS Box), but `-O0` to `-O2` still works.

(2) `-ipo` cannot find its required files in New `TIGER` (tigercpu) and could result in warnings, but `-ip` may work. 

Also make change these fields:
```    
    AR  = xiar rcv
    DBGFLAGS= -qopenmp -O0 -g -CA -CB -CS -traceback -debug all -ftrapuv -check noarg_temp_created -WB -warn none
```

Remark:

(1) `-CB` checks for out-of-bounds array.

(2) `-ftrapuv` stops the code if see division by zero, which is very useful. However, for DNS Box case, `-ftrapuv` indicates such a bug when calling `HYPRE_StructGMRESSetup` and stops the code. Since this is an internal subroutine of the `hypre` library, we can only see the input, and thus there is no easy way to debug here and figure out why. This bug is still there waiting for future solution.

One mustalso comment out `use mpi` in `library/parallel.f90` and uncomment `include 'mpif.h'`.

### Old `TIGER` (tiger) 

AS OF 08/2016: If the above libraries were installed using the latest Intel compilers (`intel/16.0/64/16.0.2.181` on `TIGER` as of 08/16), modify the following fields in `Makefile.in`:
```
    F90 = mpifort
    LD  = mpifort
    AR  = xiar rcv
    LDFLAGS  = -qopenmp
    DBGFLAGS = -qopenmp -O0 -g -CA -CB -CS -CV -traceback -debug all -ftrapuv -check noarg_temp_created -WB -warn none
    OPTFLAGS = -qopenmp -O3 -xhost -ipo -qoverride-limits -nowarn
```
If using the `openmpi/intel-16.0` module on `TIGER`, comment out `use mpi` in `library/parallel.f90` and uncomment `include 'mpif.h'`. Otherwise an error message like the following will be received:
```
parallel.f90(116): error #6285: There is no matching specific subroutine for this generic subroutine call.   [MPI_ALLREDUCE]
```

Also add the lines `module load intel/16.0` and `module load openmpi/intel-16.0` to the `SLURM` script.



<!-- Additionally, to speed up compilation of `NGA`, open `Makefile` and edit the instructions under the target `opt` to be `@make -j 4 "FLAGS = $(OPTFLAGS)"`. Similarly, edit the instructions under the target `debug` to be `@make -j 4 "FLAGS = $(OPTFLAGS)"`. Note that on `TIGER`, specify 8 cores instead of 4 cores for faster compilation. -->

## Compiling

If on a cluster, use the `module list` command to show what modules are currently loaded. If these are different than the modules loaded when installing the libraries for `NGA`, use the `module purge` command and then `module load` the correct modules, as described in the page on [libraries](library-installation).

Compilation of `NGA` then occurs in one single step.  To compile the code for production runs (faster code, longer compilation), execute the command

    make opt

Alternatively, to compile the code for debugging (slower code, shorter compilation), for which additional diagnostics will be printed to the screen for debugging purposes, execute the command

    make debug

Unless there is a specific need, the code should always be compile in `opt` mode.

If compilation of the code was successful, the `NGA` `bin` directory should contain several executables, notably, `ARTS` and `init_flow`.

Remark: `init_flow` and `ARTS` must be generated by the same compiler flags, otherwise the reading of initial flow by `ARTS` will be misaligned and incorrect (e.g., all zeros). May need to `make distclean` and recompile.

----
# Generating a data file
When designing a grid for simulations in [NGA](home), specify the grid geometry in an `.f90` file.
For the purpose of this tutorial, the file is named `example_simulation.f90` and placed in `~/NGA/src/tools/init_flow`.

Also consider the other `.f90` files within the above directory for guidance on grid generation. If interpolating grid points using the non-linear fitting x<sub>j</sub> = ab<sup>(j-1)</sup>+c, then use the `Matlab` program `GridCoeffs.m` to calculate a, b, and c from three given points (j, x(j)), j > 0. The Matlab program `GridEval.m` can be used to evaluate x<sub>j</sub> for any value of j. Both Matlab programs are attached to this page.

In order for [NGA](home) to incorporate our grid specifications into `init_flow`, we will append `example_simulation.f90` to the end of the `F90FILES` field in the `Makefile` of `~/NGA/src/tools/init_flow`. The `Makefile` should look like:

    include ../../Makefile.in

    F90FILES = param.f90 taylor.f90 ...
               ...
               jetcart2.f90 simplejet.f90 \
               example_simulation.f90


Additionally, we need to add the case corresponding to `example_simulation.f90` before `case default` in `init_flow.f90`. The format should be similar to

    case ('desired_case_name')
       call subroutine1_grid
       call subroutine2_data
       call subroutine3_optdata


Lastly, generate the `init_flow` executable by calling `make opt` in `~/NGA/src`. Create the `config`, `data.0`, and `optdata.0` files by calling

    ~/NGA/bin/init_flow input

in the same directory where the `input` file is located. The `config`, `data.0`, and `optdata.0` files can be moved onto the `TIGER` cluster using `scp`.

[GridCoeffs.m](/uploads/fd929afa442a81e3f58f4e4b8d526d71/GridCoeffs.m)
[GridEval.m](/uploads/954ec14da47bad5987a2d3cd2949f6b9/GridEval.m)

To interpolate data from one grid to another, first make sure that the new resolution grids for the specific problem/case are already specified in the `/tools/init_flow/`

Second, generate the new `config` using the new input file with the new grid, then call

    ~/NGA/bin/interpolatedata

Enter the source `config`, the target `config`, the source `data.*`, and the target `data.*`. After execution, the target `data.*` will be generated for use. 

The `data.*` could also be `optdata.*`. `inflow` normally can be interpolated automatically, so no need to do above.

----
# Code execution

This page covers the basics of running `NGA` that are uniform across all machines where the code may be run. For machine-specific details on running `NGA`, see 
[Running on `TIGER`](running-on-tiger), and [Examples](examples). [NGA](home) is run in two distinct steps.

1. `init_flow` is executed to create a `config` file and an initial `data` (and possibly `optdata`) file.  The config file contains all of the information about the flow geometry; the data file contains the velocities and scalars.  `init_flow` is run in serial and takes a single command line argument: the name of the input file.  For example, `init_flow` can be run by executing the following command

        /home/<username>/NGA/bin/init_flow <input file>


2. The actual flow solver `ARTS` is executed.  `ARTS` is run in parallel and takes a single command line argument (in addition to the `MPI` arguments): the name of the input file.  For example, `ARTS` can be run with `mpirun` on 16 processors by executing the command

        mpirun -np 16 /home/<username>/NGA/bin/arts <input file>

At minimum, the solver will dump limited diagnostic information (time, maximum velocities, etc.) to the screen and more detailed diagnostic information (maximum velocities, convergence information, CFL numbers, timings, etc.) to a `monitor` directory.

Remark: when changing the number of nodes/cores/omp-threads, it is OK to read the same initial flow or restart data file, and there is no need to rerun the `init_flow`.

To avoid having to write out the path to the `NGA` executables when running them, do the following. Open the `~/.bashrc` file (if on a workstation) or `~/.bash_profile` file (if on a cluster). Add the line

    export PATH=$PATH:/home/<username>/<name of NGA directory>/bin

and save the file. Then, execute `source .bashrc` (replace .bashrc with .bash_profile if on a cluster). It should now be possible to run `init_flow`, `ARTS`, and other `NGA` executables with commands like `init_flow <input file>`.


## Input File Keywords

Note: This page will continue to be updated with a menu of all keywords for `ARTS`.  Keywords for `init_flow` for each example case are documented and can be used for the [NGA example](examples) sample scripts.

The following is a list of input file keywords for `ARTS`.  Note that the keywords are case sensitive.  Many of the keywords do not need to be in the input file; default values will simply be assumed (this is being phase out).

### Files

* **Configuration file**: The name of the config (geometry) file.
* **Data file to read**: The name of the data file that the simulation will be initialized from.
* **Data file to write**: The name of the data file that will be written during the simulation.  These two files should not be the same.
* **Data frequency**: The time interval (not frequency) at which a data file is written.  This interval is in terms of the global simulation time, not the time since a simulation was locally restarted.

### Partitioning

* **Processors along X**: The number of processors used to partition the x-direction.  The total number of processors is the product of the x-, y-, and z-directions.
* **Processors along Y**: The number of processors used to partition the y-direction.  The total number of processors is the product of the x-, y-, and z-directions.
* **Processors along Z**: The number of processors used to partition the z-direction.  The total number of processors is the product of the x-, y-, and z-directions.
* **OpenMP threads**: The number of `OpenMP` threads to be used.

### Chemistry

* **Chemistry model**: The option 'none' for incompressible flow, 'finitechem' for detailed chemistry, and various options for flamelet models are supported.  For the 'none' chemistry model:
  * **Density**: The density of the flow.
  * **Viscosity**: The **''dynamic**'' viscosity of the flow.
  * **Diffusivity**: The **''dynamic**'' diffusivity of the flow.  Currently, all scalars must have the same diffusivity.
  * **Temperature**: the temperature of the flow.  This is not really needed and should be removed.

### Subfilter Model

* **Subgrid Scale model**: Logical `.true.` or `.false.` for a subfilter turbulence model for LES.  Currently, only dynamic Smagorinksy-like models are included.
* **SGS averaging**: Procedure for dynamically determining the coefficient in the subfilter turbulence model.  Three types of averaging are supported: `Ad-Hoc` (local volume averaging), `Germano` (homogeneous direction averaging), and `Lagrangian`.  Generally speaking, Lagrangian averaging is the best.  See instructions on how to use Lagrangian averaging at the bottom of this page.
* **SGS Override Lag.**: Set to `.true.` for fresh run, and to `.false.` for rerun. Such setting is for SGS accuracy.

### Time Advancement

* **Timestep size**: The value of the timestep.  The timestep is set by the most restrictive of this value and the **CFL number**.
* **CFL number**: The maximum allowable explicit-direction CFL number.  The timestep is set by the most restrictive of the this value and the **timestep size**.
* **Subiterations**: For incompressible flows, two subiterations are required.  We could do a lot of analysis of this...

### End of Simulation

* **Maximum iterations**: The number of iterations before the simulation is stopped.  The simulation is stopped when the first of these three conditions is reached.
* **Maximum time**: The maximum simulation time before the simulation is stopped.  This time is the global simulation time, not the time since a simulation was locally restarted.  The simulation is stopped when the first of these three conditions is reached.
* **Maximum wall time**: The maximum (wall) clock time before the simulation is stopped.  This time is a relative time to the restart of the simulation.  Set this time to be slightly smaller than the requested simulation time to guarantee a restart file is written before the job is time-out.  The simulation is stopped when the first of these three conditions is reached.

### Velocity
* **Velocity conv scheme**: The desired order of accuracy for the convective terms in the momentum equation.  Currently, only second-, fourth-, and sixth-order schemes are explicitly coded.
* **Velocity visc scheme**: The desired order of accuracy for the viscous terms in the momentum equation.  Currently, only second-, fourth-, and sixth-order schemes are explicitly coded.
* **Implicit directions**: Directions for which the momentum equation is evolved in time with the semi-implicit Crank-Nicolson method.  In parallel, only second-order schemes can be evolved implicitly.

### Pressure
* **Pressure solver**: Iterative linear solver used to solve the elliptic pressure equation.  Options include built-in conjugate gradients solvers and Krylov and multigrid solvers from hypre.  A whole Wiki page could be devoted to the pressure solvers...
* **Pressure precond**: Preconditioner used for the solution of the pressure equation.  Options include several simple built-in preconditions (diagonal, tridiagonal, incomplete Cholesky factorization) for use with the built-in solvers and a simple diagonal or more involved multigrid preconditions for hypre.
* **Pressure cvg**: The convergence criterion for solving the pressure equation.  The basic recommendation is 1.0e-08 unless problems arise.  The iterative solver stops once this convergence criterion is met or the maximum number of iterations is reached.
* **Pressure iterations**: The maximum number of iterations allowed for the pressure solver.  The basic recommendation is a sufficiently large number that the convergence criterion is met first.
* **Pressure fft**: Logical for using discrete Fourier transforms for solving the pressure equation in homogeneous (periodic) directions.  There are very few reasons for this to ever be set to false.

## Lagrangian SGS Averaging
If doing LES and wish to use Lagrangian averaging, do the following:

1. **Creating the `optdata` file.** Run `init_flow` for the simulation case as described above. This should generate a `data` file and potentially an `optdata` file. 
  a. If an `optdata` file is generated, run `editData` (another executable that is installed with `NGA` in the `/home/<username>/NGA/bin/` directory) and choose option 1, "Print Min Max of Variable". Check that the variable list includes LM_VEL and MM_VEL. If the simulation includes any scalars, the `optdata` file should also contain LM_<SC> and MM_<SC> for each scalar (e.g., ZMIX, ZMIX2, and PROG for an FPVA combustion simulation. If the `optdata` file is missing any of these variables, add them as described in step 2. Otherwise, continue to step 3.
  b. If no `optdata` file was generated, create a blank one by copying the data file generated by `init_flow` and emptying the variable list: 

        cp <data file> optdata

     
Then, run `editData` (another executable that is installed with `NGA` in the `/home/<username>/NGA/bin/` directory) and choose option 5, "Empty variable list".

2. **Adding variables to the `optdata` file.** The `optdata` file should always contain LM_VEL and MM_VEL. If the simulation includes any scalars, the `optdata` file should also contain LM_<SC> and MM_<SC> for each scalar (e.g., ZMIX, ZMIX2, and PROG for an FPVA combustion simulation. To add these variables, run `editData` and choose option 2, "Add variable". Add each variable individually. Typically, it is best to set the LMs to 0 and the MMs to 1.

3. **Running in `NGA`.** Add the following to the `NGA` input file:

        ! Files 
        Optional data to read :         <input optdata file>
        Optional data to write:         <output optdata file>

        ! Subfilter turbulence model
        Use SGS model :                 .true.
        SGS averaging :                 Lagrangian
        SGS Override Lag. :             .true.

After running the simulation for the first time, change the value of `SGS Override Lag.` to `.false`.

----
# Running on `TIGER`

Production runs on `TIGER` are submitted through the `SLURM` scheduling system. A `SLURM` script is needed to tell the cluster how long and on how many processors the job will run. An annotated example script for new `TIGER` (aka tigercpu) is shown below. The header of the `SLURM` file tells the cluster information about the job being submitted, and the rest of the file is the commands that are needed to run the job. The script should load the modules for whatever Intel and `MPI` version used when installing `NGA`. See the [Examples](examples) page for a sample problem that uses this procedure.

## Submitting `SLURM` Scripts

Submit `SLURM` scripts (e.g., `samplescript.slurm`) by executing

```
sbatch samplescript.slurm
``` 

in the directory containing the script. The [SLURM website](https://slurm.schedmd.com/man_index.html) has a useful list of commands, notably `scancel <job id>` to cancel a submitted job and `squeue -u $USER` to check on the status of a job.

### Queuing rules

It is important to understand how jobs are queued on the cluster. Here are some highlights:

* Each job is assigned a priority based on job length, job size, group priority, and time spent in the queue. The job in the queue with the highest priority is scheduled first, and starts whenever sufficient nodes become available. Then the second highest priority job is scheduled, and so on. Smaller lower priority jobs can start before larger higher priority jobs if they can run while the the system waits for space for the higher priority job. For example, if there are 8 nodes free, a low priority 8 node job may start before a high priority 16 node job. 
* **Job length**: shorter jobs get higher priority, with test jobs (less than or equal to 1 hour) getting the highest priority. The cutoffs for priority are 1, 6, 24, 72, and 168 hours. Jobs with lengths of 1:01, 5:59, and 6:00 all have identical priority, but the 1:01 job may still start first because it is easier to fit into the schedule. Each user can have at most two test jobs running at once, and there is a limit on the total number of nodes for the other job lengths, which can be found [here](https://www.princeton.edu/researchcomputing/computational-hardware/tiger/user-guidelines/).
* **job size**: larger jobs get higher priority
* **group priority**: proportional to funding provided for the cluster. If we use more than our share, our priority goes down. If we use less, our priority goes up.
* **time in queue**: priority *slowly* increases as the job spends time in the queue

The [Princeton research computing website](https://www.princeton.edu/researchcomputing/computational-hardware/tiger/) is a good resource for the guidelines for running on `TIGER`, including limits on job size, total number of available cores, etc. The research computing technicians can answer questions at `cses@princeton.edu`.

## Note on Number of Processors

Running `NGA` requires that the number of processors in both the `SLURM` and `NGA` input files match.  This page is created to ensure that the processors are allocated appropriately and the run fully capitalizes on requested computational resources from the `TIGER` cluster.

### `NGA` Side
As stated in the [execution page](code-execution), processor partitioning can be done in three directions and into a number of `OpenMP` threads; however, in order to correctly capitalize on threading, the processors must be allocated correctly.  Processors can be allocated in each of three directions: **X**, **Y**, and **Z**, which appear as options in the `NGA` input file.  The product of the numbers of processors in all directions (x\*y\*z) yields the total number of `MPI` processes per task, not strictly the number of processors.  In addition to the allocation of `MPI` processes, `OpenMP` threads can also be be assigned to each partition.  This is given in the input file as `OpenMP threads`, unsurprisingly.  The total number of processors should equal to the product of the number of processors in all directions and `OpenMP` threads (x\*y\*z \* OpenMP).

### `SLURM` Side
Allocating the processors in the `SLURM` file can be slightly more complicated.  On (old) `TIGER`, there are 16 processors on each node to allocate, and any number of nodes can be allocated to each job. On (new) `TIGER`, there are 40 processors per node. `SLURM` provides the option to allocate a number of nodes (nnodes) and a number of tasks per node (ntasks-per-node), as well as a number of CPUs (processors) per task (cpus-per-task).  So, how should each of these nodes be divided up?  There are three rules which should be generally obeyed to ensure the job submits as intended.
1. The product of the number of nodes and number of tasks per node should be equal to the product of x-y-z partitions. This is the total number of `MPI` tasks and that should follow the `-n` in the `srun` command. The sample script provided automatically calcuates and sets this number.
2. The number of CPUs per task should be equal to the number of `OpenMP` threads. The `SLURM` file should contain the line `export OMP_NUM_THREADS = ` this number, prior to the `srun` command
3. The product of the number of tasks per node and number of CPUs per task should be equal to 16 (old `TIGER`) or 40 (new `TIGER`).
If these rules are followed, the job should submit correctly. 

## Annotated `SLURM` Script Example

```
#!/bin/bash
# HEADER for Parallel job using 80 processors:
#SBATCH --nodes=2	       	# number of nodes
#SBATCH --ntasks-per-node=10	# number of processors per node
#SBATCH --cpus-per-task=4       # number of cpus per task
#SBATCH -t 3:00:00		# run for 3 hr max
#SBATCH --mail-type=begin	# send email when process begins...
#SBATCH --mail-type=end		# ...and when it ends...
#SBATCH --mail-type=fail	# ...or when it fails.
#SBATCH --mail-user=<your-email>@princeton.edu # send notifications to this email
#SBATCH -e job.err              # Name of output file for error messages
#SBATCH -o job.out              # Name of output file for standard output

# BODY - commands to be run
# Load required modules
module load intel               
module load intel-mpi

# Set number of openmp threads and number of MPI tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASKS=$(echo "$SLURM_NNODES*$(echo $SLURM_TASKS_PER_NODE | cut -d '(' -f 1)" | bc -l)

# Print some information
echo Directory is `pwd`
echo Time is `date`
echo
echo This job runs on the following processors:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_NNODES nodes with $SLURM_TASKS_PER_NODE cores per node.
echo

# run arts
srun -n $NTASKS ~/NGA/bin/arts input_pipe
```

----
# Examples

## LES of Pipe / Annulus Flow on New `TIGER`

This example is a simulation of fully developed flow in a pipe and surrounding annulus. The fluid properties correspond to methane in the pipe and air in the annulus. An initial velocity profile is specified, but over time the flow transitions to turbulent and becomes fully developed because a periodic boundary condition is used in the *x* (axial) direction. The simulation is run using the attached input file as described in [the section on running the code](code execution). After running `init_flow`, create an `optdata` file. The process for doing this is given on the same page. Also attached is a sample script for submitting the job to the `SLURM` scheduler on `TIGER`. Copy the script into the same directory as the input file, and then submit the job using the command

    sbatch pipe.slurm


The attached script is set up to run the job on 80 processors (2 nodes, 40 cores each). In the `NGA` input file, it is specified that the domain will be split into 8 segments in the *y* (radial) direction, and each segment will be computed on four processors using `OpenMP`. After submitting the job, check its status using 

    squeue -u <username>

[input_pipe](uploads/f0696bb8e74c0d5e0f3297316a88753b/input_pipe)
[pipe.slurm](uploads/856c4ce27f35e08f50512923212850b4/pipe.slurm)

----
# Converting `Chemkin` format to `FlameMaster` and `NGA` format

Using finite rate chemistry in `NGA` requires a chemical mechanism, which typically are provided in `Chemkin` format and need to be converted to be useable by `NGA`. This process assumes there is a mechanism file, a thermodynamic data file, and a transport data file, all in chemkin format. We will refer to these files, in order, as `mech.dat`, `therm.dat`, and `tran.dat`. There are tools installed with `FlameMaster` to facilitate the conversion process.

Navigate to the tool location in `FlameMaster` (modify the path if `FlameMaster` is installed somewhere else)

```
cd ~/FlameMaster/FlameManTools/CK2FMReinh/
```
There is a `README` file here that covers the first half of this process. It is messy and poorly written, but may have more info regarding errors.

The first time using these tools, it is necessary to build a chemkin interpreter, as well as prepare two perl scripts.

```
chmod 755 maki
./maki
chmod 755 *.perl
```
Now copy the three `.dat` files into the current directory. The interpreter does not like comments in the mechanism files (specifically, those that come before the actual mechanism). Run the following command on all three of the `.dat` files to remove them.

```
sed '/^!/ d' mech.dat > mech_stripped.dat
sed '/^!/ d' therm.dat > therm_stripped.dat
sed '/^!/ d' tran.dat > tran_stripped.dat
```
Next, run `ckintrp3seiser`

```
./ckintrp3seiser
```
Four inputs are required. The first two are the `therm` and `mech` files, and the last two are names for output files.

```
thermo data file name = ?: therm_stripped.dat
input model file name = ?: mech_stripped.dat
interpreter ascii output file name = ?: a.i
CHEMKIN link file name = ?: linkfile
```
Open `a.i` and check at the very end if there are any errors. If everything is correct, near the very end should be the message:

```
NO ERRORS FOUND ON INPUT...CHEMKIN LINKING FILE WRITTEN.
```
If there are PLOGs in the mechanism, see http://www.nuigalway.ie/combustionchemistrycentre/softwaredownloads/

Execute the first `perl` script:

```
./mechi2tex.perl a.i 1 2
```
The useful output of this script is `a.mech`. 

Execute the second `perl` script:

```
./modmech.perl -t therm_stripped.dat -r tran_stripped.dat -o MECHNAME.mech a.mech
```

There will likely be errors and warnings here. Read through them.

Two important files are created. `MECHNAME.mech` is now a `FlameMaster` format mechanism. `newthermofile` is a combined transport and thermodynamic data file for `FlameMaster`. Run the following:

```
CreateBinFile -i newthermofile -o MECHNAME.thermo.bin
```

Finally, run `ScanMan`:

```
ScanMan -w -i MECHNAME.mech -t MECHNAME.thermo.bin -S > MECHNAME.out
```

Look at the `.out` file for errors and warnings. Check the `README` for error information. 

`ScanMan` produces a `.pre` file, which is the mechanism for `FlameMaster`. The `-w` flag also produces a Fortran90 mechanism file for use with `NGA`. The two important files are `MECHNAMEF90.h` and `MECHNAMEF.f90`. 

## Using the mechanism in `NGA`
First, the `.h` and `.f90` files from before should be copied over to `NGA`, specifically into:

```
NGA/src/combustion
```
Only one chemical mechanism can be compiled at any one time in `NGA`. Therefore to remove the previous mechanisms. First, open:

```
NGA/src/combustion/Makefile
```
One of the `F90FILES` should be a mechanism file (ideally labeled as a fuel species for easy identification). Replace that filename with the name of the `.f90` file. Next, within the `combustion` directory, run:

```
make clean
```
This is necessary to remove the previous mechanism. The next step involves sterilizing the chemical species names. `NGA` currently has an eight character limit on species names. Open up the mechanism `.f90` file, and search for `GETSPECIESNAMES`. This will open a subroutine where the names from the `.h` file are converted to `NGA` friendly names. Everything to the right of the equal sign needs to be converted to at most an 8 character unique name. However, make sure there are 20 characters between each quotation mark, use spaces to pad the string. Now compile `NGA` from the `src` directory.

```
cd ..
cd make opt
```
## Initializing for a finite chem simulation
When creating an initialization for the simulation, add a new case so as not to change others work. Call this case by the chemical mechanism, since specifying and initializing species will be the major difference from others. Fill a `names` array, which involves using the 8 character maximum names from before, and then initialize data for each species. A example case would look something like:
```
     nvar = 17
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     names = ''
     data = 0.0_WP

     ! Names of general variables
     names(1) = 'U'
     names(2) = 'V'
     names(3) = 'W'
     names(4) = 'P'
     names(5) = 'RHO'
     names(6) = 'dRHO'
     names(16) = 'T'
     names(17) = 'ZMIX'

     ! Names of species
     names(7) = 'N2'
     names(8) = 'H'
     names(9) = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'HO2'
     names(15) = 'H2O2'
```
Note that velocities and pressure come first in a given order, then other scalars like density, viscosity, diffusivity, etc. Next should come all species, in the same order and with the same names as from the previous section. Finally, temperature should be second from last, and mixture fraction **MUST BE LAST**. Finally, the data should be initialized as:
```
data(:,:,:,k) = whatever
```
Where the first three indices represent location in (*x*,*y*,*z*), and the fourth index matches the index from the names selected previously.


----
# `LeWI`

`LeWI` (Lend When Idle) is a dynamic library designed to speed up HPC hybrid applications (i.e., two levels of parallelism) by improving the load balance of the outer level of parallelism (e.g., `MPI`) by dynamically redistributing the computational
resources at the inner level of parallelism (e.g., `OpenMP`). at run time. This dynamism allows `LeWI` to react to different sources of imbalance: Algorithm, data, hardware architecture and resource availability among others.The algorithm redistributes the computational resources that are not being used from one process to another process inside the same shared memory node in order to speed up its execution.

# Library Installation

1. Build requirements
* A supported platform running `GNU`/`Linux` (`i386`, `x86-64`, `ARM`, `PowerPC`, or `IA64`)
* C compiler
* `Python 2.4` or higher (`Python 3` recommended)
* `GNU` Autotools, only needed if you want to build from the repository.

2. Download the `LeWI` source code:
* Either from our website: [DLB Downloads][].
* Or from a git repository. Clone `LeWI` repository from GitHub:
      ```bash
        git clone https://github.com/bsc-pm/dlb.git
      ```
* Or download from [GitHub releases][]
* Bootstrap autotools:
      ```bash
      cd dlb
      ./bootstrap
      ```
            
3. Run `configure`. Optionally, check the configure flags by running `./configure -h` to see detailed information about some features. `MPI` support must be enabled with ``--with-mpi`` and, optionally, an argument telling where `MPI` can be located.

    ```bash
    ./configure --prefix=<DLB_PREFIX> [<configure-flags>]
    ```
    
4. Build and install

    ```bash
    make
    make install
    ```
5. Optionally, add the installed bin directory to your `PATH`

    ```bash
    export PATH=<DLB_PREFIX>/bin:$PATH
    ```

## Basic usage

Choose between linking or preloading the binary with the `LeWI` shared library `libdlb.so` and configure using the environment variable `DLB_ARGS`.

1. **Example 1:** Share CPUs between `MPI` processes

    ```bash
    # Link application with DLB
    mpicc -o myapp myapp.c -L<DLB_PREFIX>/lib -ldlb -Wl,-rpath,<DLB_PREFIX>/lib

    # Launch MPI as usual, each process will dynamically adjust the number of threads
    export DLB_ARGS="--lewi"
    mpirun -n <np> ./myapp
    ```

2. **Example 2:** Share CPUs between `MPI` processes with advanced affinity
control through `OMPT`.

    ```bash
    # Link application with an OMPT capable OpenMP runtime
    OMPI_CC=clang mpicc -o myapp myapp.c -fopenmp

    # Launch application:
    #   * Set environment variables
    #   * DLB library is preloaded
    #   * Run application with binary dlb_run
    export DLB_ARGS="--lewi --ompt"
    export OMP_WAIT_POLICY="passive"
    preload="<DLB_PREFIX>/lib/libdlb.so"
    mpirun -n <np> <DLB_PREFIX>/bin/dlb_run env LD_PRELOAD="$preload" ./myapp
    ```

3. **Example 3:** Manually reduce assigned CPUs to an `OpenMP` process.

    ```bash
    # Launch an application preloading DLB
    export OMP_NUM_THREADS=4
    export DLB_ARGS="--drom"
    export LD_PRELOAD=<DLB_PREFIX>/lib/libdlb.so
    taskset -c 0-3 ./myapp &

    # Reduce CPU binding to [1,3] and threads to 2
    myapp_pid=$!
    dlb_taskset -p $myapp_pid -c 1,3

    ```

4. **Example 4:** Get a `TALP` summary report at the end of an execution

    ```bash
    export DLB_ARGS="--talp --talp-summary=pop-metrics"
    PRELOAD=<DLB_PREFIX>/lib/libdlb_mpi.so
    mpirun <opts> env LD_PRELOAD="$PRELOAD" ./app
    ```

#### User Guide
Refer to the [DLB User Guide][] for a more complete documentation. For questions, suggestions and bug reports, the support team can be contacted via e-mail at `dlb@bsc.es`.

----
# Reference Mapping

Once the compilation is successful, any case running with a standard Cartesian grid may further utilize the zonal reference mapping scheme. Set the `refmapping` as `active` in `chemistryProperties` file to use the reference mapping method (it is necessary to add an empty `refmapping{}` `dict` even if it is unused):

```
refmapping
{
    active  true;
    
    mixtureFractionProperties
    {
        oxidizerMassFractions
        {
            N2       0.77;
            O2       0.23;
        }

        fuelMassFractions
        {
            NC12H26       1.0;
        }

        #include "$CASE/constant/thermo.f90"
    }
    tolerance	1e-4;  // mixture fraction tolerance
    deltaT	2; // temperature tolerance
}
```
Reference mapping uses mixture fraction ($Z$) and maps a reference solution to reference cells satisfying a condition. The entry above sets the $Z=0$ and $Z=1$ conditions from given mass fractions. For each iteration it finds a reference solution where $Z <$`tolerance` and solves the chemistry. Subsequent cells following the same condition are mapped from this reference solution. When `deltaT` is explicitly set, the scheme also checks the temperature between reference solution and other reference cells and ensures: $|T_\textrm{cell}-T_\textrm{ref}|<$`deltaT`.

----
# Citations

* O. Desjardins, G. Blanquart, G. Balarac, and H. Pitsch. High order conservative finite difference scheme for variable density low mach number turbulent flows. *Journal of Computational Physics*, 227(15):7125–7159, 2008. DOI: https://doi.org/10.1016/j.jcp.2008.03.027.
* J. F. MacArt and M. E. Mueller. Semi-implicit iterative methods for low mach number turbulent reacting flows: Operator splitting versus approximate factorization. *Journal of Computational Physics*, 326:569–595, 2016. DOI: https://doi.org/10.1016/j.jcp.2016.09.016.
* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban, D. E. Shumaker, and C. S. Woodward. SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers. *ACM Trans. Math. Softw.*, 31(3):363–396, 2005. DOI: https://doi.org/10.1145/1089014.1089020.
* B. Tekgül, P. Peltonen, H. Kahila, O. Kaario, and V. Vuorinen. DLBFoam: An open-source dynamic load balancing model for fast reacting flow simulations in OpenFOAM. *Computer Physics Communications*, 267:108073, 2021. DOI: https://doi.org/10.1016/j.cpc.2021.108073.
* I. Morev, B. Tekgül, M. Gadalla, A. Shahanaghi, J. Kannan, S. Karimkashi, O. Kaario, and V. Vuorinen. Fast reactive flow simulations using analytical Jacobian and dynamic load balancing in OpenFOAM. *Physics of Fluids*, 34(2):021801, Feb 2022. DOI: https://doi.org/10.1063/5.0077437.
   D. Rettenmaier, D. Deising, Y. Ouedraogo, E. Gjonaj, H. De Gersem, D. Bothe, C. Tropea, and H. Marschall. Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM. *SoftwareX*, 10:100317, 2019. DOI: https://doi.org/10.1016/j.softx.2019.100317.
* W. Zhang, A. Myers, K. Gott, A. Almgren, and J. Bell. AMReX: Block-structured adaptive mesh refinement for multiphysics applications. *The International Journal of High Performance Computing Applications*, 35(6):508–526, 2021. DOI: https://doi.org/10.1177/10943420211022811.
* M. T. Henry de Frahan, J. S. Rood, M. S. Day, S. Hariswaran, S. Yellapantula, B. A. Perry, R. W. Grout, A. Almgren, W. Zhang, J. B. Bell, and J. H. Chen. PeleC: An adaptive mesh refinement solver for compressible reacting flows. *The International Journal of High Performance Computing Applications*, 37(2):115–131, 2023. DOI: https://doi.org/10.1177/10943420221121151.
* L. D. Owen, W. Ge, M. Rieth, M. Arienti, L. Esclapez, B. S. Soriano, M. E. Mueller, M. Day, R. Sankaran, and J. H. Chen. PeleMP: The Multiphysics Solver for the Combustion Pele Adaptive Mesh Refinement Code Suite. *Journal of Fluids Engineering*, 146(4):041103, Feb 2024. DOI: https://doi.org/10.1115/1.4064494.
* M. Garcia, J. Corbalan, and J. Labarta. LeWI: A runtime balancing algorithm for nested parallelism. In *International Conference on Parallel Processing, 2009. ICPP ’09.*, pages 526–533, Sep 2009. DOI: https://doi.org/10.1109/ICPP.2009.56.
* M. Garcia, J. Labarta, and J. Corbalan. Hints to improve automatic load balancing with lewi for hybrid applications. *Journal of Parallel and Distributed Computing*, 74(9):2781–2794, 2014. DOI: https://doi.org/10.1016/j.jpdc.2014.05.004.
* M. Garcia, J. Corbalan, R.M. Badia, and J. Labarta. A dynamic load balancing approach with SMPSuperscalar and MPI. In R. Keller, D. Kramer, and J.-P. Weiss, editors, *Facing the Multicore - Challenge II*, volume 7174 of *Lecture Notes in Computer Science*, pages 10–23. Springer Berlin Heidelberg, 2012. DOI: http://dx.doi.org/10.1007/978-3-642-30397-5_2.
* R. D. Blumofe and Charles E. Leiserson. Scheduling multithreaded computations by work stealing. *J. ACM*, 46(5):720–748, Sep 1999. DOI: https://doi.org/10.1145/324133.324234.
* M. Raju, M. Wang, M. Dai, W. Piggott, and D. Flowers. Acceleration of detailed chemical kinetics using multi-zone modeling for CFD in internal combustion engine simulations. In *SAE 2012 World Congress & Exhibition*. SAE International, Apr 2012. DOI: https://doi.org/10.4271/2012-01-0135.
* R. S. Barlow. Sandia H2/He flame data. https://tnfworkshop.org/data-archives/simplejet/, 2003. Release 2.0.
* R. S. Barlow and C. D. Carter. Raman/Rayleigh/LIF measurements of nitric oxide formation in turbulent hydrogen jet flames. *Combustion and Flame*, 97(3):261–280, 1994. DOI: https://doi.org/10.1016/0010-2180(94)90020-5.
* R. S. Barlow and C. D. Carter. Relationships among nitric oxide, temperature, and mixture fraction in hydrogen jet flames. *Combustion and Flame*, 104(3):288–299, 1996. DOI: https://doi.org/10.1016/0010-2180(95)00123-9.
* R. S. Barlow and J. H. Frank. Effects of turbulence on species mass fractions in methane/air jet flames. *Symposium (International) on Combustion*, 27(1):1087–1095, 1998. Twenty-Seventh Sysposium (International) on Combustion Volume One. DOI: https://doi.org/10.1016/S0082-0784(98)80510-9.
* C. E. Lacey, A. G. Novoselov, and M. E. Mueller. In-Situ Adaptive Manifolds: Enabling computationally efficient simulations of complex turbulent reacting flows. *Proceedings of the Combustion Institute*, 38(2):2673–2680, 2021. DOI: https://doi.org/10.1016/j.proci.2020.06.207.
* Z. Chen. *Studies on the Initiation, Propagation, and Extinction of Premixed Flames*. Phd thesis, Princeton University, Princeton, NJ, Jan 2009. Available at http://www2.coe.pku.edu.cn/tpic/2011812212957550.pdf.
* E. F. Toro. *Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction*. Springer Berlin Heidelberg, Berlin, Heidelberg, 2009. DOI: https://doi.org/10.1007/b79761_1.
* P. G. Tucker. *Advanced Computational Fluid and Aerodynamics*. Cambridge Aerospace Series. Cambridge University Press, 2016.

<!---
@article{SUNDIALS,
    author      = {Hindmarsh, A. C. and Brown, P. N. and Grant, K. E. and Lee, S. L. and Serban, R. and Shumaker, D. E. and Woodward, C. S.},
    title       = {{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
    year        = {2005},
    issue_date  = {September 2005},
    publisher   = {Association for Computing Machinery},
    address     = {New York, NY, USA},
    volume      = {31},
    number      = {3},
    issn        = {0098-3500},
    url         = {https://doi.org/10.1145/1089014.1089020},
    doi         = {10.1145/1089014.1089020},
    journal     = {ACM Trans. Math. Softw.},
    pages       = {363–396},
    numpages    = {34},
}

@incollection{DLB_OpenMP_SMPS,
    year        = {2012},
    isbn        = {978-3-642-30396-8},
    booktitle   = {Facing the Multicore - Challenge II},
    volume      = {7174},
    series      = {Lecture Notes in Computer Science},
    editor      = {Keller, R. and Kramer, D. and Weiss, J.-P.},
    doi         = {10.1007/978-3-642-30397-5_2},
    title       = {A Dynamic Load Balancing Approach with {SMPSuperscalar} and {MPI}},
    url         = {http://dx.doi.org/10.1007/978-3-642-30397-5_2},
    publisher   = {Springer Berlin Heidelberg},
    author      = {Garcia, M. and Corbalan, J. and Badia, R.M. and Labarta, J.},
    pages       = {10-23},
    language    = {English}
}

@INPROCEEDINGS{LeWI_ICPP09,
    author      = {Garcia, M. and Corbalan, J. and Labarta, J.},
    booktitle   = {International Conference on Parallel Processing, 2009. ICPP '09.},
    title       = {{LeWI}: A Runtime Balancing Algorithm for Nested Parallelism},
    year        = {2009},
    month       = {Sep},
    pages       = {526-533},
    doi         = {10.1109/ICPP.2009.56},
    issn        = {0190-3918}
}

@article{GARCIA2014,
    title       = {Hints to improve automatic load balancing with LeWI for hybrid applications},
    journal     = {Journal of Parallel and Distributed Computing},
    volume      = {74},
    number      = {9},
    pages       = {2781-2794},
    year        = {2014},
    issn        = {0743-7315},
    doi         = {10.1016/j.jpdc.2014.05.004},
    url         = {https://www.sciencedirect.com/science/article/pii/S0743731514000926},
    author      = {Garcia, M. and Labarta, J. and Corbalan, J.}
}

@article{PeleC,
    author      = {Henry de Frahan, M. T. and Rood, J. S. and Day, M. S. and Hariswaran, S. and Yellapantula, S. and Perry, B. A. and Grout, R. W. and Almgren, A. and Zhang, W. and Bell, J. B. and Chen, J. H.},
    title       = {{PeleC}: An adaptive mesh refinement solver for compressible reacting flows},
    journal     = {The International Journal of High Performance Computing Applications},
    volume      = {37},
    number      = {2},
    pages       = {115-131},
    year        = {2023},
    doi         = {10.1177/10943420221121151},
    url         = {https://doi.org/10.1177/10943420221121151},
    eprint      = {https://doi.org/10.1177/10943420221121151}
}

@article{AMReX,
    author      = {Zhang, W. and Myers, A. and Gott, K. and Almgren, A. and Bell, J.},
    title       = {{AMReX}: {B}lock-structured adaptive mesh refinement for multiphysics applications},
    journal     = {The International Journal of High Performance Computing Applications},
    volume      = {35},
    number      = {6},
    pages       = {508-526},
    year        = {2021},
    doi         = {10.1177/10943420211022811},
    url         = {https://doi.org/10.1177/10943420211022811},
    eprint      = {https://doi.org/10.1177/10943420211022811}
}

@article{RETTENMAIER2019,
    title       = {Load balanced {2D} and {3D} adaptive mesh refinement in {OpenFOAM}},
    journal     = {SoftwareX},
    volume      = {10},
    pages       = {100317},
    year        = {2019},
    issn        = {2352-7110},
    doi         = {https://doi.org/10.1016/j.softx.2019.100317},
    url         = {https://www.sciencedirect.com/science/article/pii/S2352711018301699},
    author      = {Rettenmaier, D. and Deising, D. and Ouedraogo, Y. and Gjonaj, E. and {De Gersem}, H. and Bothe, D. and Tropea, C. and Marschall, H.}
}

@article{PeleMP,
    author      = {Owen, L. D. and Ge, W. and Rieth, M. and Arienti, M. and Esclapez, L. and Soriano, B. S. and Mueller, M. E. and Day, M. and Sankaran, R. and Chen, J. H.},
    title       = {{PeleMP}: {The Multiphysics Solver for the Combustion Pele Adaptive Mesh Refinement Code Suite}},
    journal     = {Journal of Fluids Engineering},
    volume      = {146},
    number      = {4},
    pages       = {041103},
    year        = {2024},
    month       = {Feb},
    issn        = {0098-2202},
    doi         = {10.1115/1.4064494},
    url         = {https://doi.org/10.1115/1.4064494},
    eprint      = {https://asmedigitalcollection.asme.org/fluidsengineering/article-pdf/146/4/041103/7243623/fe\_146\_04\_041103.pdf}
}

@article{DLBFoam_1,
    title       = {{DLBFoam}: An open-source dynamic load balancing model for fast reacting flow simulations in {OpenFOAM}},
    journal     = {Computer Physics Communications},
    volume      = {267},
    pages       = {108073},
    year        = {2021},
    issn        = {0010-4655},
    doi         = {https://doi.org/10.1016/j.cpc.2021.108073},
    url         = {https://www.sciencedirect.com/science/article/pii/S0010465521001855},
    author      = {Tekg\"{u}l, B. and Peltonen, P. and Kahila, H. and Kaario, O. and Vuorinen, V.}
}

@article{DLBFoam_2,
    author      = {Morev, I. and Tekg\"{u}l, B. and Gadalla, M. and Shahanaghi, A. and Kannan, J. and Karimkashi, S. and Kaario, O. and Vuorinen, V.},
    title       = {Fast reactive flow simulations using analytical {J}acobian and dynamic load balancing in {OpenFOAM}},
    journal     = {Physics of Fluids},
    volume      = {34},
    number      = {2},
    pages       = {021801},
    year        = {2022},
    month       = {Feb},
    doi         = {10.1063/5.0077437},
    url         = {https://doi.org/10.1063/5.0077437},
    eprint      = {https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/5.0077437/19868913/021801\_1\_5.0077437.pdf}
}

@article{Work_Stealing,
    author      = {Blumofe, R. D. and Leiserson, Charles E.},
    title       = {Scheduling multithreaded computations by work stealing},
    year        = {1999},
    issue_date  = {Sept. 1999},
    publisher   = {Association for Computing Machinery},
    address     = {New York, NY, USA},
    volume      = {46},
    number      = {5},
    issn        = {0004-5411},
    url         = {https://doi.org/10.1145/324133.324234},
    doi         = {10.1145/324133.324234},
    journal     = {J. ACM},
    month       = {Sep},
    pages       = {720–748},
    numpages    = {29}
}

@article{BARLOW1998,
    title       = {Effects of turbulence on species mass fractions in methane/air jet flames},
    journal     = {Symposium (International) on Combustion},
    volume      = {27},
    number      = {1},
    pages       = {1087-1095},
    year        = {1998},
    note        = {Twenty-Seventh Sysposium (International) on Combustion Volume One},
    issn        = {0082-0784},
    doi         = {10.1016/S0082-0784(98)80510-9},
    url         = {https://www.sciencedirect.com/science/article/pii/S0082078498805109},
    author      = {Barlow, R. S. and Frank, J. H.}
}

@article{BARLOW1994,
    title       = {{R}aman/{R}ayleigh/{LIF} measurements of nitric oxide formation in turbulent hydrogen jet flames},
    journal     = {Combustion and Flame},
    volume      = {97},
    number      = {3},
    pages       = {261-280},
    year        = {1994},
    issn        = {0010-2180},
    doi         = {10.1016/0010-2180(94)90020-5},
    url         = {https://www.sciencedirect.com/science/article/pii/0010218094900205},
    author      = {Barlow, R. S. and Carter, C. D.}
}

@article{BARLOW1996,
    title       = {Relationships among nitric oxide, temperature, and mixture fraction in hydrogen jet flames},
    journal     = {Combustion and Flame},
    volume      = {104},
    number      = {3},
    pages       = {288-299},
    year        = {1996},
    issn        = {0010-2180},
    doi         = {10.1016/0010-2180(95)00123-9},
    url         = {https://www.sciencedirect.com/science/article/pii/0010218095001239},
    author      = {Barlow, R. S. and Carter, C. D.}
}

@inproceedings{Zonal,
    author      = {Raju, M. and Wang, M. and Dai, M. and Piggott, W. and Flowers, D.},
    title       = {Acceleration of Detailed Chemical Kinetics Using Multi-zone Modeling for {CFD} in Internal Combustion Engine Simulations},
    booktitle   = {SAE 2012 World Congress \& Exhibition},
    publisher   = {SAE International},
    month       = {Apr},
    year        = {2012},
    doi         = {10.4271/2012-01-0135},
    url         = {https://doi.org/10.4271/2012-01-0135}
}

@phdthesis{Chen_thesis,
    title       = {Studies on the Initiation, Propagation, and Extinction of Premixed Flames},
    author      = {Chen, Z.},
    year        = {2009},
    month       = {Jan},
    address     = {Princeton, NJ},
    note        = {Available at \url{http://www2.coe.pku.edu.cn/tpic/2011812212957550.pdf}},
    school      = {Princeton University},
    type        = {PhD thesis}
}

@article{LACEY2021,
    title   = {In-Situ Adaptive Manifolds: Enabling computationally efficient simulations of complex turbulent reacting flows},
    journal = {Proceedings of the Combustion Institute},
    volume  = {38},
    number  = {2},
    pages   = {2673--2680},
    year    = {2021},
    issn    = {1540-7489},
    doi     = {10.1016/j.proci.2020.06.207},
    url     = {https://www.sciencedirect.com/science/article/pii/S1540748920302984},
    author  = {Lacey, C. E. and Novoselov, A. G. and Mueller, M. E.}
}

@article{DESJARDINS2008,
	title	= {High order conservative finite difference scheme for variable density low Mach number turbulent flows},
	journal	= {Journal of Computational Physics},
	volume	= {227},
	number	= {15},
	pages	= {7125-7159},
	year	= {2008},
	issn	= {0021--9991},
	doi	= {https://doi.org/10.1016/j.jcp.2008.03.027},
	url	= {https://www.sciencedirect.com/science/article/pii/S0021999108001666},
	author	= {Desjardins, O. and Blanquart, G. and Balarac, G. and Pitsch, H.}
}

@article{MACART2016,
    title   = {Semi-implicit iterative methods for low Mach number turbulent reacting flows: Operator splitting versus approximate factorization},
    journal = {Journal of Computational Physics},
    volume  = {326},
    pages   = {569-595},
    year    = {2016},
    issn    = {0021-9991},
    doi     = {10.1016/j.jcp.2016.09.016},
    url     = {https://www.sciencedirect.com/science/article/pii/S0021999116304284},
    author  = {MacArt, J. F. and Mueller, M. E.}
}

@book{Tucker2016,
    place       =   {Cambridge},
    series      =   {Cambridge Aerospace Series},
    title       =   {Advanced Computational Fluid and Aerodynamics},
    doi         =   {10.1017/CBO9781139872010},
    publisher   =   {Cambridge University Press},
    author      =   {Tucker, P. G.},
    year        =   {2016},
    collection  =   {Cambridge Aerospace Series}
}

@book{Toro2009,
    author      =   {Toro, E. F.},
    title       =   {Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction},
    year        =   {2009},
    publisher   =   {Springer Berlin Heidelberg},
    address     =   {Berlin, Heidelberg},
    doi         =   {10.1007/b79761_1},
    url         =   {https://doi.org/10.1007/b79761_1}
}
-->
