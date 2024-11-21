gcc -c -fPIC -O3 model_IHC.c
gcc -c -fPIC -O3 model_Synapse.c
gcc -c -fPIC -O3 complex.c 
gcc -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o

