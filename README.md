# AuditoryNerveFiber.jl

AuditoryNerveFiber.jl is a Julia package that provides access to auditory-nerve models in Julia. 
Most implementations of auditory-nerve models in the literature are written in C, MATLAB, Python, or a mixture these languages.
This package provides access these models in Julia using a consistent interface so that users can model easily in Julia without worrying about the details of implementation for each model or knowing the native language of each model.

## Why a new package? Why Julia?
Existing bindings for some models already exist in several Python or MATLAB packages such as [cochlea](https://github.com/mrkrd/cochlea), [apcmodels](https://github.com/guestdaniel/apcmodels), and [Auditory Modeling Toolbox](https://amtoolbox.org/).
This package was written in Julia in part as a learning exercise and in part to supplement existing packages and explore several unique benefits of Julia:
- Interoperability between Julia and C is excellent. Julia has a ["no boilerplate" philosophy](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/) that results in C bindings that are less complicated than those of Python or MATLAB.
- Julia can be exceptionally fast, but it is still a high-level langauge. This allows us to extend tools and models written in low-level langauges without taking the performance hits usually associated with writing in a high-level language. 
- Julia has excellent standard tools for documentation and reproducible computing. For example, this package's documentation is automatically rendered from Markdown to HTML and is available at LINK, and this package's set of dependencies is precisely version controlled in `Project.toml`.
- Julia has a rapidly developing ecosystem of state-of-the-art packages that could be easily integrated with this package. Examples relevant to those in the auditory modeling field including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), [Flux.jl](https://github.com/FluxML/Flux.jl), and [Turing.jl](https://github.com/TuringLang/Turing.jl). 

## Implemented models

Presently, AuditoryNerveFiber.jl provides access to:
- Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and
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

