# Examples

The best way to learn to use the package is to build from examples.

```@meta
CurrentModule = AuditoryNerveFiber
```
## Pure-tone response from one ANF

Simulating responses from a single auditory-nerve fiber is easy.
The first element returned from [`sim_an_zbc2014`](@ref) is the analytic firing rate approximation, so we can pass a pure-tone stimulus through [`sim_ihc_zbc2014`](@ref) to [`sim_an_zbc2014`](@ref) and then extract the first element to get the instantaneous firing rate of the auditory nerve responding to the pure tone.

```@example
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0   # frequency and CF, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 10e4       # sampling rate, Hz
level = 50.0    # level, dB SPL

# Synthesize a pure tone
x = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);
# Simulate response 
y = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, freq), freq)[1];
plot(y[1:2000], ylabel="Firing rate (sp/s)", xlabel="Samples")
```

## Iso-level tuning curve

Simulating an iso-level tuning curve requires a bit more work.
Here, we define simple functions that synthesize a pure tone and simulate a response to that pure tone.
Then, we use `map` to simulate an average response at several pure-tone frequencies. 

```@example
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Plots
using Statistics

# Define variables
cf = 1000.0     # CF, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 10e4       # sampling rate, Hz
level = 50.0    # level, dB SPL

# Define a function to synthesize a pure tone
pure_tone(freq) = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);
# Define a function to simulate a single response
anf_response(x) = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, cf), cf)[1];
# Generate a log-spaced frequency axis
freqs = ASU.LogRange(200.0, 20000.0, 50)
# Synthesize tone and simulate response at each freq
results = map(freq -> mean(anf_response(pure_tone(freq))), freqs)
plot(freqs, results, ylabel="Firing rate (sp/s)", xlabel="Frequency (Hz)", xscale=:log)
```

## Generating spike trains

[`sim_an_zbc2014`](@ref) has multiple outputs. 
The first is an analytic approximation of the underlying instantaneous firing rate.
The second is an analytic approximation of the variance of the underlying instantaneous variance.
The third is a spike train. 
Here, we can see how to extract and analyze each of these outputs.

```@example
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0   # frequency and CF, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 10e4       # sampling rate, Hz
level = 50.0    # level, dB SPL

# Synthesize a pure tone
x = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);
# Simulate response 
(mean, var, spikes) = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(x, freq), freq);

# Plot
l = @layout [a ; b]
p1 = plot(mean, ylabel="Firing rate (sp/s)", xlabel="Samples")
p2 = plot(spikes, ylabel="Firing rate (sp/s)", xlabel="Samples")
plot(p1, p2, layout=l)
```

The top row shows the analytic firing rate.
The bottom row shows a single example spike train. 

## Extending functions with multiple dispatch

The available interface may not always satisfy your needs.
For example, you may need to simulate responses at several CFs but the default method for [`sim_ihc_zbc2014`](@ref) is only defined for a scalar CF parameter.
With Julia's multiple dispatch system, you can readily define your own extensions to the methods provided by AuditoryNerveFiber.jl.
Rather than define methods for common use cases in this package, we defer to the user or packages that extend this package to define methods that suit their needs (i.e., this package is merely a thin implementation of underlying models and not a modeling toolbox). 

```@example
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Plots

# Define extension to vector-valued CF parameters
function ANF.sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1}; kwargs...)
    map(_cf -> ANF.sim_ihc_zbc2014(input, _cf; kwargs...), cf)
end

# Define variables
freq = 1000.0   # freq, Hz
cfs = [500.0, 1000.0, 1500.0]  # CFs, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 10e4       # sampling rate, Hz
level = 50.0    # level, dB SPL

# Define a function to synthesize a pure tone
pure_tone = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);
# Simulate IHC response at several CFs
results = ANF.sim_ihc_zbc2014(pure_tone, cfs; species="human")
plot([result[1:1000] for result in results], layout=3, labels="CF = " .* string.(hcat(cfs...)))
```

## Plotting neurograms

Neurograms are likewise easy to generate, so long as we organize the simulations correctly.
Here, we define the same method for [`sim_ihc_zbc2014`](@ref)

```@example
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Plots

# Define extension to vector-valued CF parameters
function ANF.sim_ihc_zbc2014(input::Array{Float64, 1}, cf::Array{Float64, 1})
    map(_cf -> ANF.sim_ihc_zbc2014(input, _cf), cf)
end

# Define variables
freq = 1000.0   # freq, Hz
cfs = collect(ASU.LogRange(200.0, 20000.0, 100))  # CFs, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 10e4       # sampling rate, Hz
level = 50.0    # level, dB SPL

# Define a function to synthesize a pure tone
pure_tone = ASU.scale_dbspl(ASU.pure_tone(freq, phase, dur, fs), level);
# Simulate IHC response at several CFs
results = ANF.sim_ihc_zbc2014(pure_tone, cfs)


heatmap(transpose(hcat(results...))[:, 1:3000], xlabel="Samples", ylabel="CF (#)")
```

