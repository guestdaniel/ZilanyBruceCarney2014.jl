using BinaryBuilder

name = "libzbc2014"
version = v"0.1.0"
sources = [
    DirectorySource("external")
]

script = raw"""
cd ${WORKSPACE}/srcdir/
gcc -c -Wall -fPIC -O3 model_IHC.c
gcc -c -Wall -fPIC -O3 model_Synapse.c
gcc -c -Wall -fPIC -O3 complex.c 
gcc -c -Wall -fPIC -O3 test.c 
gcc -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o test.o
echo ${libdir}
mkdir ${libdir}
cp libzbc2014.so ${libdir}
"""

platforms = [supported_platforms()[2]]

products = [
    LibraryProduct("libzbc2014", :libzbc2014)
]

dependencies = []

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)