# ZilanyBruceCarney2014.jl

ZilanyBruceCarney2014.jl is a Julia package that provides access to the Zilany, Bruce, and Carney (2014) auditory-nerve model in Julia. 

```
Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and expanded simulation options for a model of the auditory periphery. The Journal of the Acoustical Society of America, 135(1), 283–286. http://dx.doi.org/10.1121/1.4837815
```

## Why a new package? Why Julia?
Existing bindings for some models already exist in several MATLAB or Python packages such as [the original code release](https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/Zilany-2014-Code-and-paper.zip), [cochlea](https://github.com/mrkrd/cochlea) and [Auditory Modeling Toolbox](https://amtoolbox.org/).
This package was written in Julia to supplement existing packages and leverage several unique benefits of Julia:
- Interoperability between Julia and C is excellent. Julia has a ["no boilerplate" philosophy](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/) that results in C bindings that easy to write and maintain.
- Julia can be exceptionally fast, but it is still a high-level langauge. This allows us to extend tools and models written in low-level langauges (such as C) using a high-level language (Julia) while avoiding some of the performance penalty usually associated with writing in a high-level language. 
- Julia has a rapidly developing ecosystem of state-of-the-art packages that could be easily integrated with this package. Examples relevant to those in the auditory modeling field including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), [Flux.jl](https://github.com/FluxML/Flux.jl), and [Turing.jl](https://github.com/TuringLang/Turing.jl). 

# Installation
Installation is very easy.
Simply open your package REPL by pressing `]` while your Julia REPL is open, and then type:
```
add ZilanyBruceCarney2014
```
The Julia package manager should then automatically install the Julia source code and download the appropriate packaged C binaries for your platform.
Binaries are available for most processor architectures on Linux (using either `glibc` or `musl`), for 32- and 64-bit Windows platforms, and for 64-bit Mac platforms.

In your own scripts or packages, you can put
```
using ZilanyBruceCarney2014
```
at the top and then use the functions described below to invoke the model.

If you run into any problems, please reach out to daniel_guest@urmc.rochester.edu.

# Usage 
There are three key functions for users of the model: `sim_ihc_zbc2014`, `sim_anrate_zbc2014`, and `sim_spikes_zbc2014`. 
You can find documentation for these functions in the source code or in the REPL, accessible by first typing `?` in the Julia REPL to access help and then typing the name of the function.
(Note that the package must be loaded via `using` before you look for help — see above for details.)
Each of these functions accepts two positional arguments:
1. Vector-valued input waveform. For `sim_ihc_zbc2014`, this is the acoustic waveform; for the other functions, this is the output of `sim_ihc_zbc2014`.
2. Scalar-valued CF in Hz.

Other parameters are passed as keyword arguments. For example, simulating an IHC response at 4 kHz CF using the cat model would be done as follows:
```
ihc_output = sim_ihc_zbc2014(stimulus, 4e3; species="cat")
```

Direct bindings are available in the form of the functions `IHCAN!`, `Synapse!`, and `SingleAN!`, which emulate the behaviors of the corresponding C functions in the model source code.
Note that the exclamation marks indicate that these functions operate on (some of) their arguments in-place, just as the original functions in C do. 
Most users will not need to interact with these functions.

# Testing
Many basic response properties of the auditory-nerve simulations (e.g., responses grow in response in increasing sound level) are tested in `test/runtests.jl`. 
If you want to run these tests yourself, follow these steps:
- Clone the repository (`git clone git@github.com:guestdaniel/ZilanyBruceCarney2014.jl`)
- Change directory to the repository 
- Open a Julia REPL
- Switch to the Pkg REPL (press `]` on your keyboard)
- Instantiate the package's dependencies (`instantiate` in the Pkg REPL)
- Call `test` from the Pkg REPL

# Funding
Development of this package was supported by the following funding resources at various points in time:
- NIH R01 DC005216
- NIH F31 DC019247
- NIH R01 DC010813
- NIH F32 DC022143
- UMN College of Liberal Arts Graduate Fellowship

# License and acknowledgments
This repository is licensed under the [GNU AGPL v3 license](https://www.gnu.org/licenses/agpl-3.0.en.html). 
The underlying model code is largely the work of:
- Muhammad S. A. Zilany
- Ian C. Bruce
- Laurel H. Carney