# Reference

Here, you can see the documentation for each type and method provided by AuditoryNerveFiber.jl.
Higher-level functions designed for end users are documented first, while lower-level functions designed mostly for internal are documented later.

```@meta
CurrentModule = AuditoryNerveFiber
```

## Wrappers

Wrappers are functions that "wrap" around lower-level bindings and/or implementations and provide more convenience than the bindings themselves.
For example, when calling wrappers a user should never have to worry about handling pointers, pre-allocating arrays, etc.

```@docs
sim_ihc_zbc2014
sim_synapse_zbc2014
sim_an_zbc2014
```

## Bindings

Bindings are functions that directly interface with external implementations of models in other languages. 
All bindings here work generally the same way:
- Inputs are passed directly to external functions with minimal or no input checking or preprocessing
- No defaults are provided
- Names and types of variables match those of the original-language source code
Generally, it is expected that end-users will never need to call these functions. 
Instead, they should use a wrapper or another higher-level function to access the same functionality in a a more "user-friendly" way.

```@docs
IHCAN!
Synapse!
SingleAN!
```