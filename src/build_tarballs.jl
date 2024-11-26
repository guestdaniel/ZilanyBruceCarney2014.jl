using BinaryBuilder

name = "libzbc2014"
version = v"0.2.1"
sources = [
    DirectorySource("external")
]

# Build for Linux
script = raw"""
cd ${WORKSPACE}/srcdir/
gcc -c -fPIC model_IHC.c
gcc -c -fPIC model_Synapse.c
gcc -c -fPIC complex.c 
gcc -shared -o libzbc2014.so model_IHC.o model_Synapse.o complex.o
echo ${libdir}
mkdir ${libdir}
cp libzbc2014.so ${libdir}
"""

# Build for Windows
script_windows = raw"""
cd ${WORKSPACE}/srcdir/
gcc -c -fPIC model_IHC.c
gcc -c -fPIC model_Synapse.c
gcc -c -fPIC complex.c 
gcc -shared -o libzbc2014.lib model_IHC.o model_Synapse.o complex.o
echo ${libdir}
mkdir ${libdir}
cp libzbc2014.lib ${libdir}
"""

# platforms = [supported_platforms()[1],     # Linux i686 (glibc)
#              supported_platforms()[2],     # Linux x64 (glibc)
#              supported_platforms()[6],     # Linux i686 (musl)
#              supported_platforms()[7]]     # Linux i686 (musl)

platforms = [supported_platforms()[[15, 16]]]

products = [
    LibraryProduct("libzbc2014", :libzbc2014)
]

dependencies = []

build_tarballs(ARGS, name, version, sources, script_windows, platforms, products, dependencies)
