using BinaryBuilder

name = "libzbc2014"
version = v"0.1.0"
sources = [
    DirectorySource("external")
]

script = raw"""
cd ${WORKSPACE}/srcdir/
clang -c -Wall -fPIC -Ofast model_IHC.c
clang -c -Wall -fPIC -Ofast model_Synapse.c
clang -c -Wall -fPIC -Ofast complex.c 
clang -c -Wall -fPIC -Ofast test.c 
clang -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o test.o
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