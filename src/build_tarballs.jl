using BinaryBuilder, Pkg

name = "libzbc2014"
version = v"0.2.1"
sources = [
#    DirectorySource("external")
     GitSource("https://github.com/guestdaniel/ZilanyBruceCarney2014.jl.git", ""),
]

# Build for Linux
script = raw"""
cd ${WORKSPACE}/srcdir/external
gcc -c -fPIC complex.c -o complex.o
gcc -c -fPIC model_IHC.c -o model_IHC.o
gcc -c -fPIC model_Synapse.c -o model_Synapse.o
gcc -shared -o "libzbc2014.${dlext}" model_IHC.o model_Synapse.o complex.o
mkdir ${libdir}
cp "libzbc2014.${dlext}" ${libdir}
"""

platforms = [
    Platform("i686", "linux"; libc="glibc"),
    Platform("x86_64", "linux"; libc="glibc"),
    Platform("i686", "windows"),
    Platform("x86_64", "windows"),
]

products = [
    LibraryProduct("libzbc2014", :libzbc2014)
]

dependencies = []

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies)
