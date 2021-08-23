gcc -c -Wall -fPIC model_IHC.c
gcc -c -Wall -fPIC model_Synapse.c
gcc -c -Wall -fPIC complex.c 
gcc -shared -o libihc.so model_IHC.o model_Synapse.o complex.o
