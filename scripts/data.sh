#!/bin/bash

#GCC
cd ..
make clean
make CC=gcc OFLAGS=-O1

taskset -c 2 ./nbody0 > results/result_gcc/O1/nbody0_O1.dat
taskset -c 2 ./nbody1 > results/result_gcc/O1/nbody1_O1.dat
taskset -c 2 ./nbody2 > results/result_gcc/O1/nbody2_O1.dat
./nbody3 > results/result_gcc/O1/nbody3_O1.dat
./nbody4 > results/result_gcc/O1/nbody4_O1.dat
./nbody5 > results/result_gcc/O1/nbody5_O1.dat
./nbody6 > results/result_gcc/O1/nbody6_O1.dat
./nbody7 > results/result_gcc/O1/nbody7_O1.dat
./nbody8 > results/result_gcc/O1/nbody8_O1.dat


make clean
make CC=gcc OFLAGS=-O2

taskset -c 2 ./nbody0 > results/result_gcc/O2/nbody0_O2.dat
taskset -c 2 ./nbody1 > results/result_gcc/O2/nbody1_O2.dat
taskset -c 2 ./nbody2 > results/result_gcc/O2/nbody2_O2.dat
./nbody3 > results/result_gcc/O2/nbody3_O2.dat
./nbody4 > results/result_gcc/O2/nbody4_O2.dat
./nbody5 > results/result_gcc/O2/nbody5_O2.dat
./nbody6 > results/result_gcc/O2/nbody6_O2.dat
./nbody7 > results/result_gcc/O2/nbody7_O2.dat
./nbody8 > results/result_gcc/O2/nbody8_O2.dat


make clean
make CC=gcc OFLAGS=-O3

taskset -c 2 ./nbody0 > results/result_gcc/O3/nbody0_O3.dat
taskset -c 2 ./nbody1 > results/result_gcc/O3/nbody1_O3.dat
taskset -c 2 ./nbody2 > results/result_gcc/O3/nbody2_O3.dat
./nbody3 > results/result_gcc/O3/nbody3_O3.dat
./nbody4 > results/result_gcc/O3/nbody4_O3.dat
./nbody5 > results/result_gcc/O3/nbody5_O3.dat
./nbody6 > results/result_gcc/O3/nbody6_O3.dat
./nbody7 > results/result_gcc/O3/nbody7_O3.dat
./nbody8 > results/result_gcc/O3/nbody8_O3.dat


make clean
make CC=gcc OFLAGS=-Ofast

taskset -c 2 ./nbody0 > results/result_gcc/OFast/nbody0_OFast.dat
taskset -c 2 ./nbody1 > results/result_gcc/OFast/nbody1_OFast.dat
taskset -c 2 ./nbody2 > results/result_gcc/OFast/nbody2_OFast.dat
./nbody3 > results/result_gcc/OFast/nbody3_OFast.dat
./nbody4 > results/result_gcc/OFast/nbody4_OFast.dat
./nbody5 > results/result_gcc/OFast/nbody5_OFast.dat
./nbody6 > results/result_gcc/OFast/nbody6_OFast.dat
./nbody7 > results/result_gcc/OFast/nbody7_OFast.dat
./nbody8 > results/result_gcc/OFast/nbody8_OFast.dat



#----------------------------------------------------------------#
#----------------------------------------------------------------#
#----------------------------------------------------------------#
#----------------------------------------------------------------#


#LLVM CLANG
make clean
make CC=clang OFLAGS=-O1

taskset -c 2 ./nbody0 > results/result_clang/O1/nbody0_O1.dat
taskset -c 2 ./nbody1 > results/result_clang/O1/nbody1_O1.dat
taskset -c 2 ./nbody2 > results/result_clang/O1/nbody2_O1.dat
./nbody3 > results/result_clang/O1/nbody3_O1.dat
./nbody4 > results/result_clang/O1/nbody4_O1.dat
./nbody5 > results/result_clang/O1/nbody5_O1.dat
./nbody6 > results/result_clang/O1/nbody6_O1.dat
./nbody7 > results/result_clang/O1/nbody7_O1.dat
./nbody8 > results/result_clang/O1/nbody8_O1.dat


make clean
make CC=clang OFLAGS=-O2

taskset -c 2 ./nbody0 > results/result_clang/O2/nbody0_O2.dat
taskset -c 2 ./nbody1 > results/result_clang/O2/nbody1_O2.dat
taskset -c 2 ./nbody2 > results/result_clang/O2/nbody2_O2.dat
./nbody3 > results/result_clang/O2/nbody3_O2.dat
./nbody4 > results/result_clang/O2/nbody4_O2.dat
./nbody5 > results/result_clang/O2/nbody5_O2.dat
./nbody6 > results/result_clang/O2/nbody6_O2.dat
./nbody7 > results/result_clang/O2/nbody7_O2.dat
./nbody8 > results/result_clang/O2/nbody8_O2.dat


make clean
make CC=clang OFLAGS=-O3

taskset -c 2 ./nbody0 > results/result_clang/O3/nbody0_O3.dat
taskset -c 2 ./nbody1 > results/result_clang/O3/nbody1_O3.dat
taskset -c 2 ./nbody2 > results/result_clang/O3/nbody2_O3.dat
./nbody3 > results/result_clang/O3/nbody3_O3.dat
./nbody4 > results/result_clang/O3/nbody4_O3.dat
./nbody5 > results/result_clang/O3/nbody5_O3.dat
./nbody6 > results/result_clang/O3/nbody6_O3.dat
./nbody7 > results/result_clang/O3/nbody7_O3.dat
./nbody8 > results/result_clang/O3/nbody8_O3.dat


make clean
make CC=clang OFLAGS=-Ofast

taskset -c 2 ./nbody0 > results/result_clang/OFast/nbody0_OFast.dat
taskset -c 2 ./nbody1 > results/result_clang/OFast/nbody1_OFast.dat
taskset -c 2 ./nbody2 > results/result_clang/OFast/nbody2_OFast.dat
./nbody3 > results/result_clang/OFast/nbody3_OFast.dat
./nbody4 > results/result_clang/OFast/nbody4_OFast.dat
./nbody5 > results/result_clang/OFast/nbody5_OFast.dat
./nbody6 > results/result_clang/OFast/nbody6_OFast.dat
./nbody7 > results/result_clang/OFast/nbody7_OFast.dat
./nbody8 > results/result_clang/OFast/nbody8_OFast.dat
