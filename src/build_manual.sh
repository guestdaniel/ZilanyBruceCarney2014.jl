clang -c -fPIC -O3 -march=native model_IHC.c
clang -c -fPIC -O3 -march=native model_Synapse.c
clang -c -fPIC -O3 -march=native complex.c 
clang -c -fPIC -O3 -march=native test.c 
clang -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o test.o

