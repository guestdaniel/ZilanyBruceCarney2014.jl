# AuditoryNerveFiber.jl

AuditoryNerveFiber.jl is a Julia package that provides access to several different auditory-nerve models in Julia. 
Most implementations of auditory-nerve models in the literature are written in C [[1]](#1), MATLAB, or Python, or a mixture these languages.
This package provides access to many of these models in Julia using a consistent interface so that users can model easily in Julia without worrying about the details of implementation for each model or knowing the native language of each model.

## Why a new package? Why Julia?
Existing bindings to some models do already exist in several Python or MATLAB packages such as [cochlea](https://github.com/mrkrd/cochlea), [apcmodels](https://github.com/guestdaniel/apcmodels), and [Auditory Modeling Toolbox](https://amtoolbox.org/).
This package was written in Julia to supplement the existing packages and take advantage of several unique benefits of Julia:
- Interoperability between Julia and C is excellent and significantly less complicated than interoperability between Python and C or MATLAB and C. This allows for a very "clean" bindings to be written for models natively written in C [[1]](#1) with very low overhead.
- Julia can be exceptionally fast. While MATLAB and Python can achieve similar speeds on most problems, this may require the use of performance-oriented libraries like NumPy, Numba, or Cython (in the case of Python) that can only be applied to limited problem types or may require avoiding useful abstractions or convenience features (e.g., function calls in MATLAB). Such restrictions are absent or less onerous in Julia. 
- Julia has excellent tools for documentation and reproducible computing. For example, this package's in-line documentation is automatically rendered from Markdown to HTML and is available at LINK, and this package's set of dependencies is precisely version controlled in `Project.toml`.
- Julia has a rapidly developing ecosystem of state-of-the-art packages that could be easily integrated with this package. Examples relevant to those in the auditory modeling field including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), [Flux.jl](https://github.com/FluxML/Flux.jl), and [Turing.jl](https://github.com/TuringLang/Turing.jl). 

## Implemented models

Presently, AuditoryNerveFiber.jl provides access to:
- Zilany, Bruce and Carney (2014) [[1]](#1), originally in C

<a id="1">[1]</a> Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and
expanded simulation options for a model of the auditory periphery. The Journal
of the Acoustical Society of America, 135(1), 283–286.
http://dx.doi.org/10.1121/1.4837815

## Installation

## Interoperability mechanisms

## Interface

## References
- Bruce, I. C., Erfani, Y., & Zilany, M. S. A. (2018). A phenomenological model
  of the synapse between the inner hair cell and auditory nerve: implications of
  limited neurotransmitter release sites. Hearing Research, 360(), 40–54.
  http://dx.doi.org/10.1016/j.heares.2017.12.016


- Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and
  expanded simulation options for a model of the auditory periphery. The Journal
  of the Acoustical Society of America, 135(1), 283–286.
  http://dx.doi.org/10.1121/1.4837815


- Heinz, M. G., Colburn, H. S., & Carney, L. H. (2001). Evaluating auditory
  performance limits: ii. one-parameter discrimination with random-level
  variation. Neural Computation, 13(), 2273–2316.
  http://dx.doi.org/10.1162/089976601750541813

