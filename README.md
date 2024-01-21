# Project nbody3D optimization

## Introduction
This README provides instructions on how to use the provided scripts, compile the project, and understand the project's directory structure.

## Directory Structure
The project directory is organized as follows:
- project/
- │
- ├── Makefile
- ├── nbody/
- │ ├── nbody1.c
- │ ├── nbody2.c
- │ ├── nbody3.c
- │ ├── nbody4.c
- │ ├── nbody5.c
- │ ├── nbody6.c
- │ ├── nbody7.c
- │ ├── nbody8.c
- │ └── nbody.c
- │
- ├── particle_references/
- │ └── particle_positions.dat
- │
- ├── plot_gen/
- │ ├── plot_nbody0.gp
- │ ├── plot_nbody1.gp
- │ ├── plot_nbody2.gp
- │ ├── plot_nbody3.gp
- │ ├── plot_nbody4.gp
- │ ├── plot_nbody5.gp
- │ ├── plot_nbody6.gp
- │ ├── plot_nbody7.gp
- │ └── plot_nbody8.gp
- │
- ├── results/
- │ ├── perf_plot/
- │ │ └── nbody*.png
- │ │
- │ ├── result_clang/
- │ │ ├── O1/
- │ │ ├── O2/
- │ │ ├── O3/
- │ │ └── OFast/
- │ │
- │ ├── result_gcc/
- │ │ ├── O1/
- │ │ ├── O2/
- │ │ ├── O3/
- │ │ └── OFast/
- │ │
- │ └── simplified_data/
- │ ├── clang/
- │ └── gcc/
- │
- └── scripts/
- ├── data.sh
- └── plot.sh


- **nbody/**: Contains source code files for N-body simulation.
- **particle_references/**: Contains particle reference data.
- **plot_gen/**: Contains Gnuplot scripts for generating plots.
- **results/**: Stores simulation results and performance data.
- **scripts/**: Contains scripts for benchmarking and generating plots.

## Compilation
To compile the project, you can use the provided `Makefile`. The Makefile is configured to use GCC or LLVM Clang with various optimization flags.

### Compilation Commands

- **To compile all N-body simulations**:
make all CC= GCC or CLANG OFLAGS:-OX

- **To clean up compiled files**:
make clean

- **To remove result data**:
make result_del

- **To remove generated plots**:
make plot_del

## Benchmarking and Plot Generation
Two scripts are provided in the `scripts/` directory to automate benchmarking and plot generation.

### Benchmarking Script (data.sh)
This script compiles and runs N-body simulations with different optimization levels and saves the results in the `results/` directory.

- **Usage**:
- cd scripts/
- chmod +x data.sh
- ./data.sh

### Plot Generation Script (plot.sh)
This script generates plots for N-body simulations using Gnuplot and saves them in the `results/perf_plot/` directory.

- **Usage**:
- cd scripts/
- chmod +x plot.sh
- ./plot.sh

Make sure to customize the scripts and Makefile as needed for your specific project and compiler configurations.

