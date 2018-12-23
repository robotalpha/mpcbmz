# mpcbmz
Compile and test the codes with Intel Compiler:
1. Entering the corresponding code directory(base,openmp and mpi), and using `make` to generate the binary file;
   e.g.  cd ./cbmz_baseline
         make
2. Runing the binary file directly and entering the path of the input data in input directory following the direction of the programe:
   e.g,  ./cbmz_base76
         ../input/3D_case.input
3. The calculation time would be printed after runing.
