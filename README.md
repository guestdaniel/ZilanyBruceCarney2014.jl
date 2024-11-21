# AuditoryNerveFiber.jl

AuditoryNerveFiber.jl is a Julia package that provides access to auditory-nerve models in Julia. 
Presently, this package provides access to the Zilany, Bruce, and Carney (2014) auditory-nerve model. 

## Why a new package? Why Julia?
Existing bindings for some models already exist in several Python or MATLAB packages such as [cochlea](https://github.com/mrkrd/cochlea) and [Auditory Modeling Toolbox](https://amtoolbox.org/).
This package was written in Julia to supplement existing packages and leverage several unique benefits of Julia:
- Interoperability between Julia and C is excellent. Julia has a ["no boilerplate" philosophy](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/) that results in C bindings that are less complicated than those of Python or MATLAB.
- Julia can be exceptionally fast, but it is still a high-level langauge. This allows us to extend tools and models written in low-level langauges (such as C) using a high-level language (Julia) while avoiding some of the performance penalty usually associated with writing in a high-level language. 
- Julia has a rapidly developing ecosystem of state-of-the-art packages that could be easily integrated with this package. Examples relevant to those in the auditory modeling field including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), [Flux.jl](https://github.com/FluxML/Flux.jl), and [Turing.jl](https://github.com/TuringLang/Turing.jl). 

## Implemented models

Presently, AuditoryNerveFiber.jl provides access to:
- Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and
expanded simulation options for a model of the auditory periphery. The Journal
of the Acoustical Society of America, 135(1), 283–286.
http://dx.doi.org/10.1121/1.4837815

More models to come soon!

# Installation
In this early version of AuditoryNerveFiber.jl, installation unfortunately relies on some manual steps (this will be improved soon).
Follow the steps below and reach out if you have any questions:
1. To install AuditoryNerveFiber.jl, clone the repository from GitHub, either from the download link on the GitHub page or with the following CLI git command:
```
git clone git@github.com/guestdaniel/AuditoryNerveFiber.jl
```
2. Navigate to wherever you cloned the repository and execute the compilation script `external/shared_library.sh`. 
Note that you will need to have `gcc` installed, or modify the script to suit your installed C compiler.
3. Use a text editor or your IDE to modify line 5 of `src/AuditoryNerveFiber.jl` to point to your copy of `libzcb2014.so`, which should be in the same folder as the compilation script (assuming that the script executed successfully). 
This is an unfortunate manual step that we are working to eliminate from the install process!
4. From whatever Julia envronment you want to use the package, activate the package REPL (press `]` in the regular REPL) and install with:
```
add path/to/package/folder
```

# Interface
Currently, two "levels" of functions are provided.
- The first, for models written in C, is a set of low-level bindings which emulate the function signatures of the original C functions and directly pass their inputs to a `ccall`. For the Zilany, Bruce, and Carney (2014) model, these include `IHCAN!`, `Synapse!`, and `SingleAN!` (the exclamation marks indicate that they operate on their arguments in-place, just as the original functions do). 
- The second is a set of more convenient wrappers written for end-users. These functions insulate end-users from seeing the "guts" of calling C, handling pointers, etc., and instead look and feel just like any other Julia function. These include `sim_ihc_zbc2014`, `sim_synapse_zbc2014`, and `sim_an_zbc2014`, to simulate inner-hair-cell, synapse, and auditory-nerve responses, respectively, for the Zilany, Bruce, and Carney (2014) model. Users should default to using these functions to simulate responses. 

# Testing
Many basic response properties of the auditory-nerve simulations (e.g., responses grow in response in increasing sound level) are tested in `test/runtests.jl`. 
If you want to run these tests yourself, follow these steps:
- Clone the repository (`git clone git@github.com:guestdaniel/AuditoryNerveFiber.jl`)
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

The Julia code in AuditoryNerveFiber.jl is licensed under the [GNU AGPL v3 license](https://www.gnu.org/licenses/agpl-3.0.en.html). 
The C code included in this repository can be downloaded by users and compiled per the instructions provided above.
Modification, redistribution, or reuse of the C code requires written permission from the original authors, listed below:
- Muhammad S. A. Zilany
- Ian C. Bruce
- Laurel H. Carney