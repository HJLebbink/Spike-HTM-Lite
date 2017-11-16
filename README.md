# Spike-HTM-Lite
C++ port of [Nupic HTM](https://github.com/numenta/nupic) with the aim of being lite and fast. 

I wanted to experiment with HTM but I could't easily compile it with the Intel c++ compiler; well, to be honest, I could't compile this project with Visual Studio. First I implemented the HTM in Fortran, but due to technical reasons (data structures & lack of intrinsics) I made an implementation in C++ without any depedencies. Just several header files that need to be included. Perfect for experimenting with vectorization, parallelization and MPI.

### Updates:
* 07-11-2017: Translated Fortran code to C++ code. 
