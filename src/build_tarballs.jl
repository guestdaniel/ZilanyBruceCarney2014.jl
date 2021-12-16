using BinaryBuilder

name = "libzbc2014"
version = v"0.2.1"
sources = [
    DirectorySource("external")
]

# Build for Linux
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

platforms = [supported_platforms()[1],     # Linux i686 (glibc)
             supported_platforms()[2],     # Linux x64 (glibc)
             supported_platforms()[6],     # Linux i686 (musl)
             supported_platforms()[7]]     # Linux i686 (musl)



products = [
    LibraryProduct("libzbc2014", :libzbc2014)
]

dependencies = []

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)

# Build for Windows 
# script = raw"""
# cd ${WORKSPACE}/srcdir/
# gcc -c model_IHC.c
# gcc -c model_Synapse.c
# gcc -c complex.c 
# gcc -c test.c 
# gcc -shared -o libzbc2014.dll model_IHC.o model_Synapse.o complex.o test.o
# echo ${libdir}
# mkdir ${libdir}
# cp libzbc2014.dll ${libdir}
# """

# platforms = [supported_platforms()[13]]    # Windows x64

# products = [
#     LibraryProduct("libzbc2014", :libzbc2014)
# ]

# dependencies = []

# build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)