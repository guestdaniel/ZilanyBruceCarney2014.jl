# caplab_zbc2014
Internal version of "ZilanyBruceCarney2014.jl" for development as a caplab package while a lot of things get migrated over.

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
- Clone the repository 
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
