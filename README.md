# AuditoryNerveFiber.jl

AuditoryNerveFiber.jl is a Julia package that provides access to several different auditory-nerve models in Julia. 
Most implementations of auditory-nerve models in the literature are written in C [[1]](#1), MATLAB, or Python, or a mixture these languages.
This package provides access to many of these models in Julia using a consistent interface so that users can happily model in Julia without worrying about the details of implementation for each model or knowing the native language of each model.

## Implemented models

Presently, AuditoryNerveFiber.jl provides access to:
- Zilany, Bruce and Carney (2014) [[1]](#1), originally in C

<a id="1">[1]</a> Zilany, M. S. A., Bruce, I. C., & Carney, L. H. (2014). Updated parameters and
expanded simulation options for a model of the auditory periphery. The Journal
of the Acoustical Society of America, 135(1), 283–286.
http://dx.doi.org/10.1121/1.4837815
