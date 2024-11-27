using BinaryBuilder, Pkg

name = "libzbc2014"
version = v"0.3.0"
sources = [
     GitSource("https://github.com/guestdaniel/ZilanyBruceCarney2014.jl.git", "c55d49e656c42471abf6158c187da5a92ecc187f"),
]

# Build for Linux
script = raw"""
cd ${WORKSPACE}/srcdir/ZilanyBruceCarney2014.jl/external
${CC} -c -fPIC complex.c -o complex.o
${CC} -c -fPIC model_IHC.c -o model_IHC.o
${CC} -c -fPIC model_Synapse.c -o model_Synapse.o
${CC} -shared -o "libzbc2014.${dlext}" model_IHC.o model_Synapse.o complex.o
mkdir ${libdir}
cp "libzbc2014.${dlext}" ${libdir}
"""

#platforms = [
#    Platform("i686", "linux"; libc="glibc"),
#    Platform("x86_64", "linux"; libc="glibc"),
#    Platform("i686", "windows"),
#    Platform("x86_64", "windows"),
#]
platforms = supported_platforms()

products = [
    LibraryProduct("libzbc2014", :libzbc2014)
]

dependencies = []

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat="1.6")
