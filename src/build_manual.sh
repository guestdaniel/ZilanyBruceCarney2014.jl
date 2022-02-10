gcc -c -fPIC -Ofast model_IHC.c
gcc -c -fPIC -Ofast model_Synapse.c
gcc -c -fPIC -Ofast complex.c 
gcc -c -fPIC -Ofast test.c 
gcc -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o test.o

