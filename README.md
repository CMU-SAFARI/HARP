## HARP 

This software provides the artifacts for evaluating Hybrid Active-Reactive
Profiling (HARP) as described in our MICRO 2021 academic paper to appear [1].
The HARP tool that we provide uses a combination of Monte-Carlo simulation and a
SAT solver to study the effectiveness of different strategies for profiling for
memory errors.

*Please send questions to Minesh Patel at minesh.patelh@gmail.com*

## HARP Overview

At the high level, HARP comprises three parts:

1. C++-based Monte-Carlo simulation of on-die ECC words across different ECC
   functions and error models. These files are largely an extension of the
   open-source BEER project [2, 3], and the individual file headers are used to
   indicate those files that are adapted. All source files are set up as a
   Makefile project contained within the ```src/``` directory, and library
   dependencies are provided within ```lib/```.

2. Python-based analysis scripts that parse the output of the Monte-Carlo
   simulations. These scripts are found under ```script/```.

3. Meta-scripts capable of reproducing the experiments in our paper. These
   scripts are found under ```evaluation/```.

We use Doxygen to document the source code and provide a Doxyfile for building
HTML and LaTeX documentation. To build the documentation, simply issue:

```
$ doxygen
```
when in the project directory, or point ```doxygen``` to the provided Doxyfile.
The HTML documentation will be built under ```doxygen/html/index.html```. 

## Dependencies

### C++ Simulations

Building and running HARP requires a working C++11 toolchain (note:
C++17 is necessary to build the latest Z3 from source). The C++ simulations have
two external dependencies that must be built separately:

1. The Monte-Carlo simulation relies on the [EINSim
   simulator](https://github.com/CMU-SAFARI/EINSim) [4, 5] to simulate injecting
   data-retention errors. EINSim must be built separately and its executable
   path provided to the HARP executable via the command line.

2. HARP uses the [Z3 Solver](https://github.com/Z3Prover/z3) for various SAT
   problems throughout the simulator. Building the latest version of Z3 requires
   a working C++17 toolchain. Alternatively, a system-wide installation (e.g.,
   as provided by Linux repositories) of the Z3 C++ development library may be
   used. If neither of these are possibilities, an older version of Z3 that does
   not rely upon C++17 features may be used.

Both dependencies are included as Git submodules in this repository.

This software has been built and tested on:

   - Debian 10 using GCC 8.3.0

Known potential incompatibilities:

   - Use of std::beta() (C++17 cmath, TR1 extension) that requires using
     "std=g++11" in certain Apple-LLVM versions.

### Python Dependencies

The analysis scripts require Python 3 and have dependencies on ```matplotlib >= 3.3``` and ```scipy >= 0.14.0```.

The scripts have been tested on:

   - Debian 10 using Python 3.7.3 using ```matplotlib 3.4.2``` and ```scipy 1.7.0```.

## Step-by-Step Instructions for Building and Running in Linux

### Building All Binaries

Preparing HARP for use requires three steps:

1. Building Z3 as a library

We supply Z3 as a submodule under ```lib/z3```. To build Z3, either follow Z3's
build instructions to install Z3 within ```lib/z3```, or use the convenience
script that we provide for building Z3 ```lib/build_z3.sh```. The final
installation must have (1) the development headers installed to
```lib/z3/include/*``` and the library binary installed to
```lib/z3/lib/libz3.so```.

2. Building EINSim as an executable

We supply the EINSim source code as a submodule under ```lib/einsim```. To build
EINSim, either enter the directory and issue ```make``` (default target is
sufficient) or use the convenience script that we provide for building EINSim
```lib/build_einsim.sh```.

3. Building HARP

HARP is organized as a Makefile project. Simply use:

```
$ make [-j <# threads>] [other make options] <target>
```

The makefile has various targets, described as follows:

- ```release``` builds ```harp``` with full optimizations
- ```debug``` builds ```harp.d``` with no optimization and debug symbols
- ```all``` builds both ```release``` and ```debug```
- ```doc``` builds ```doxygen``` documentation using the provided Doxyfile
- ```clean``` cleans build and binary files for both ```release``` and ```debug```

Omitting the ```target``` argument defaults to the ```release``` configuration,
which is sufficient for running HARP.

### Running HARP

HARP run as a command-line tool with several CLI options that are shown when
running the executable without options:

```
$ ./path/to/harp
```

HARP currently provides two types of analyses: (1) profiler evaluation and (2)
post-correction probability analysis. Both analyses dump a sizable amount of
output that is later parsed by the scripts within the ```script``` directory.
These analyses and scripts can be used to reproduce all figures in our paper, as
discussed within the Artifacts appendix.

We provide a convenience script at ```evaluation/run_sanity.sh``` to run all
experiments. Please see the detailed description below.

We describe each experiment workflow individually as follows. The general flow
in each case is to first run the C++ simulations and then the analysis scripts.

#### Manually running the ```Profiler Evaluation```

##### Step 1: Running the C++ simulations

The user can directly invoke the ```harp``` binary with the following command line:

```
$ ./path/to/harp <path_to_einsim> evaluations -j <JSON_directory> -k <K> -r <R> -c <N_CODES> -w <N_WORDS>
```

The arguments are defined as follows:
- ```path_to_einsim```: the path to the compiled EINSim binary
- ```JSON_directory```: the path to a directory in which to dump JSON files that describe the simulated ECC codes (will be needed for the analysis scripts) 
- ```K```: integer describing the ECC dataword length to simulate (e.g., 32, 64)
- ```R```: integer describing the random seed to use for generating random ECC codes (e.g., 0, 1, ...)
- ```S```: integer describing the random seed to use for generating random ECC words (e.g., 0, 1, ...)
- ```N_CODES```: integer describing the total number of ECC codes to simulate
- ```N_WORDS```: integer describing the total number of ECC words to simulate per ECC code

The output will be printed to stdout and will need to be redirected to a file.
For representative values of these arguments and their expected runtimes, please
see Appendix A in our paper [1].

##### Step 2: Running the Python analysis scripts

Once the data file(s) are obtained, the analysis script for this experiment can be invoked as:

```
python3 script/figures_6to10-parse_evaluation_data.py <JSON_directory> <simulation_files> -o <output_dir>
```

The arguments are defined as follows:
- ```JSON_directory```: directory containing all JSON files dumped by the C++ simulations
- ```simulation_files```: path to one or more files (or directories containing files) containing the stdout of the C++ simulations
- ```output_dir```: directory in which to dump matplotlib figures (PDF format)

#### Manually running the ```Post-Correction Probability Analysis```

##### Step 1: Running the C++ simulations

The user can directly invoke the ```harp``` binary with the following command line:

```
$ ./path/to/harp <path_to_einsim> probabilities -j <JSON_directory> -k <K> -r <R> -c <N_CODES> -w <N_WORDS>
```

The arguments are defined as follows:
- ```path_to_einsim```: the path to the compiled EINSim binary
- ```JSON_directory```: the path to a directory in which to dump JSON files that describe the simulated ECC codes (will be needed for the analysis scripts) 
- ```K```: integer describing the ECC dataword length to simulate (e.g., 32, 64)
- ```R```: integer describing the random seed to use for generating random ECC codes (e.g., 0, 1, ...)
- ```S```: integer describing the random seed to use for generating random ECC words (e.g., 0, 1, ...)
- ```N_CODES```: integer describing the total number of ECC codes to simulate
- ```N_WORDS```: integer describing the total number of ECC words to simulate per ECC code

The output will be printed to stdout and will need to be redirected to a file.
For representative values of these arguments and their expected runtimes, please
see Appendix A in our paper [1].

##### Step 2: Running the Python analysis scripts

Once the data file(s) are obtained, the analysis script for this experiment can be invoked as:

```
python3 script/figure_4-parse_postcorrection_probabilities_data.py <simulation_files> -o <output_dir>
```

The arguments are defined as follows:
- ```simulation_files```: path to one or more files (or directories containing files) containing the stdout of the C++ simulations
- ```output_dir```: directory in which to dump matplotlib figures (PDF format)

#### Using the Convenience Wrapper Script

Instead of directly running the ```harp``` binary and the analysis
scripts, we provide a convenience script at ```evaluations/run_sanity.sh```. The
script allows the user to run both experiments without micromanaging the
command-line arguments. The script is used as follows:

```
./script <harp_executable> <einsim_executable> <output_directory> <int: num_cores>
```

- ```harp_executable```: path to the ```harp``` binary
- ```einsim_executable```: path to the EINSim binary
- ```output_directory```: directory to create and hold all output files
- ```num_cores```: number of cores to run tasks with (defaults to ```nproc --all```)

Issuing an appropriate set of arguments will run the C++ simulations and
analysis scripts for both sets of experiments and output all files to the output
directory specified. The simulation configuration that will be run has its
parameters hardcoded in the script itself. We encourage the user to modify these
hardcoded parameters:
   - ```K```, ```R```, ```S```, ```N_CODES```, ```N_WORDS```: (same as specified above)
   -```N_TASKS```: specifies how many processes to spawn with the provided configuration, which will be parallelized per ```num_cores```

For representative configuration values and their expected runtimes, please see Appendix A in our paper [1].

## Licensing

The current version of the tools are provided as-is under the MIT license.

The following header-only libraries are used and are located under ```lib``` with their own license:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- [rapidjson](https://rapidjson.org/)

This software requires the Z3 solver built as a library. ```lib``` includes the
Z3-v4.8.12 source and a build script to help build Z3 as a library in-directory.
However, you may modify the Makefile to link against a different version (e.g.,
system-wide installation).
- [Z3 Solver](https://github.com/Z3Prover/z3)

## Attribution

Please cite the following paper when using the HARP Artifacts:

\[1\] Minesh Patel, Geraldo F. Oliveira, and Onur Mutlu, "HARP: Practically and Effectively Identifying Uncorrectable Errors in Memory Chips That Use On-Die Error-Correcting Codes", in the Proceedings of the 54rd Annual ACM/IEEE International Symposium on Microarchitecture (MICRO 2021), Virtual, October 2021.

Other references:

[\[2\] Minesh Patel, Jeremie S. Kim, Taha Shahroodi, Hasan Hassan, and Onur Mutlu, "Bit-Exact ECC Recovery (BEER): Determining DRAM On-Die ECC Functions by Exploiting DRAM Data Retention Characteristics", in the Proceedings of the 53rd Annual ACM/IEEE International Symposium on Microarchitecture (MICRO 2020), Virtual, October 2020.](https://people.inf.ethz.ch/omutlu/pub/BEER-bit-exact-ECC-recovery_micro20.pdf)

[\[3\] BEER on GitHub](https://github.com/CMU-SAFARI/BEER)

[\[4\] Minesh Patel, Jeremie S. Kim, Hasan Hassan, and Onur Mutlu, "Understanding and Modeling On-Die Error Correction in Modern DRAM: An Experimental Study Using Real Devices", in the Proceedings of the 49th Annual IEEE/IFIP International Conference on Dependable Systems and Networks (DSN 2019), Portland, OR, USA, June 2019.](https://people.inf.ethz.ch/omutlu/pub/understanding-and-modeling-in-DRAM-ECC_dsn19.pdf)

[\[5\] EINSim Simulator on GitHub](https://github.com/CMU-SAFARI/EINSim)
