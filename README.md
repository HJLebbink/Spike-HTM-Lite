# Spike-HTM-Lite
C++ port of [Nupic HTM](https://github.com/numenta/nupic) with the aim of being lite and fast. 

I wanted to experiment with HTM but I could't easily compile it with the Intel c++ compiler; well, to be honest, I could't compile the Nupic HTM project with Visual Studio. First I implemented HTM in Fortran, but due to technical reasons (data structures & lack of intrinsics) I made an implementation in C++ without any depedencies. Just several header files that need to be included. Perfect for experimenting with vectorization, parallelization and MPI.

### Features:
* 7 steps of history in the TP.
* 3 Layer HTM stack.
* Swarm with GA for 1, 2 and 3 layer HTM's.
* vectorized for AVX512-BW

### Updates:
* 07-11-2017: Translated Fortran code to C++ code.
* 15-11-2017: Added vectorization for AVX512-BW
