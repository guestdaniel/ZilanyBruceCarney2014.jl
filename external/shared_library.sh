gcc -c -Wall -fPIC -O3 model_IHC.c
gcc -c -Wall -fPIC -O3 model_Synapse.c
gcc -c -Wall -fPIC -O3 complex.c 
gcc -c -Wall -fPIC -O3 test.c 
gcc -shared -o libihc.so model_IHC.o model_Synapse.o complex.o test.o
