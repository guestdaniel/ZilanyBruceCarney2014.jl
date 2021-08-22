gcc -c -Wall -fPIC model_IHC.c
gcc -c -Wall -fPIC model_Syanpse.c
gcc -c -Wall -fPIC complex.c 
gcc -shared -o libihc.so model_IHC.o model_synapse.o complex.o
