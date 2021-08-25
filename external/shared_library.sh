gcc -c -Wall -fPIC -O3 -ffast-math model_IHC.c
gcc -c -Wall -fPIC -O3 -ffast-math model_Synapse.c
gcc -c -Wall -fPIC -O3 -ffast-math complex.c 
gcc -shared -o libihc.so model_IHC.o model_Synapse.o complex.o
