# Examples

The best way to learn to use the package is to build from examples.
Below, several example visualizations and shown using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

```@meta
CurrentModule = AuditoryNerveFiber
```
## Pure-tone response from one ANF

Simulating responses from a single auditory-nerve fiber is easy.
[`sim_anrate_zbc2014`](@ref) returns the instantaneous firing rate response, so we can pass a pure-tone stimulus through [`sim_ihc_zbc2014`](@ref) to [`sim_anrate_zbc2014`](@ref) and to get the instantaneous firing rate of the auditory nerve responding to the pure tone.

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0                # frequency and CF, Hz
phase = 0.0                  # starting phase, rads
dur = 0.2                    # duration, seconds
fs = 100e3                   # sampling rate, Hz
level = 50.0                 # level, dB SPL
t = 0:(1/fs):prevfloat(dur)  # time axis

# Synthesize a pure tone
x = scale_dbspl(pure_tone(freq, phase, dur, fs), level);

# Simulate response 
y = sim_an_zbc2014(sim_ihc_zbc2014(x, freq), freq)[1];
plot(
    t[1:2000], y[1:2000]; 
    ylabel="Firing rate (sp/s)", 
    xlabel="Time (s)", 
    label=:none
)
```

## Iso-level tuning curve

Simulating an iso-level tuning curve requires a bit more work.
Here, we define simple functions that synthesize a pure tone and simulate a response to that pure tone.
Then, we use `map` to simulate an average response at several pure-tone frequencies. 

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots
using Statistics

# Define variables
cf = 1000.0      # CF, Hz
phase = 0.0      # starting phase, rads
dur = 0.1        # duration, seconds
dur_ramp = 0.01  # duration of ramp, seconds
fs = 100e3       # sampling rate, Hz
level = 50.0     # level, dB SPL

# Define a function to synthesize a pure tone
stim(freq) = scale_dbspl!(cosine_ramp!(pure_tone(freq, phase, dur, fs), dur_ramp, fs), level);

# Define a function to simulate a single response
resp(stim) = sim_anrate_zbc2014(sim_ihc_zbc2014(stim, cf), cf);

# Generate a log-spaced frequency axis
freqs = LogRange(200.0, 20000.0, 100)

# Synthesize tone and simulate response at each freq
results = map(freq -> mean(resp(stim(freq))), freqs)
plot(
    freqs, results; 
    ylabel="Firing rate (sp/s)", 
    xlabel="Frequency (Hz)", 
    xscale=:log10, 
    label=:none
)
```

## Generating spike trains

[`sim_an_zbc2014`](@ref) has multiple outputs. 
The first is an analytic approximation of the underlying instantaneous firing rate.
The second is an analytic approximation of the variance of the underlying instantaneous variance.
The third is a spike train. 
Utility functions are provided to only return one of those outputs, with [`sim_anrate_zbc2014`](@ref) returning the firing rate only and [`sim_spikes_zbc2014`](@ref) returning the spike train only.
When `n_rep` is greater than 1, [`sim_spikes_zbc2014`] will return a PSTH for `n_rep` sequential presentations of the stimulus.

Here, we can see how to utilize these functions to extract and analyze these outputs.

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0                # frequency and CF, Hz
phase = 0.0                  # starting phase, rads
dur = 0.2                    # duration, seconds
dur_ramp = 0.01              # ramp duration, seconds
fs = 10e4                    # sampling rate, Hz
level = 50.0                 # level, dB SPL
t = 0:(1/fs):prevfloat(dur)  # time axis

# Synthesize a pure tone
x = scale_dbspl!(cosine_ramp!(pure_tone(freq, phase, dur, fs), dur_ramp, fs), level);

# Simulate response 
(rate, var, spikes) = sim_an_zbc2014(sim_ihc_zbc2014(x, freq), freq);

# Plot
l = @layout [a ; b]
p1 = plot(
    t, rate; 
    ylabel="Firing rate (sp/s)", 
    label=:none,
)
p2 = plot(
    t, spikes; 
    ylabel="Spike [y/n]", 
    label=:none,
)
plot(
    p1, p2; 
    layout=l,
    xlabel="Time (s)",
)
```

The top row shows the analytic firing rate.
The bottom row shows a single example spike train. 

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0                # frequency and CF, Hz
phase = 0.0                  # starting phase, rads
dur = 0.2                    # duration, seconds
dur_ramp = 0.01              # ramp duration, seconds
fs = 10e4                    # sampling rate, Hz
level = 50.0                 # level, dB SPL
t = 0:(1/fs):prevfloat(dur)  # time axis

# Synthesize a pure tone
x = scale_dbspl!(cosine_ramp!(pure_tone(freq, phase, dur, fs), dur_ramp, fs), level);

# Simulate PSTH response at many n_reps
psths = map(n_rep -> sim_spikes_zbc2014(sim_ihc_zbc2014(x, freq; n_rep=n_rep), freq; n_rep=n_rep), [1, 10, 30, 300])

# Generate plots
l = @layout [a ; b ; c; d]
plots = map(psths) do psth
    plot(
        t, psth; 
        ylabel="Spike count", 
        label=:none,
    )
end
plot(
    plots...;
    layout=l,
    xlabel="Time (s)",
)
```

## Extending functions with multiple dispatch

Via a set of [macros](https://docs.julialang.org/en/v1/manual/metaprogramming/) defined in `src/AuditoryNerveFiber.jl`, most functions provided by AuditoryNerveFiber have methods to handle a range of input types.
The base functions, such as [`sim_ihc_zbc2014`](@ref) are defined in terms of a single vector input and a single scalar characteristic frequency.
However, [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/) means that the same function can be defined multiple times for different combinations of input types. 
Each of these different definitions is called a "method" in Julia.
You can inspect the methods of `sim_ihc_zbc2014` with the `methods` function, and you will see that several methods are defined.
Documentation of these methods is pending, but the options should be fairly intuitive
- Vector input + scalar CF => Vector output
- Vector input + vector CF => Matrix output
- Matrix input + vector CF => Matrix output (one CF per row of matrix)

This makes it easy to compose inner hair cell and auditory nerve stages flexibly. 
For example, if we want to generate the response to a single acoustic input at multiple CFs, we can pass a vector input and a vector of CFs. 
This is shown in the neurogram demo below.

### Plotting neurograms

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots

# Define variables
freq = 1000.0   # freq, Hz
phase = 0.0     # starting phase, rads
dur = 0.2       # duration, seconds
fs = 100e3       # sampling rate, Hz
level = 50.0    # level, dB SPL
t = 0:(1/fs):prevfloat(dur)  # time axis, s
cfs = LogRange(200.0, 20000.0, 100)  # CFs, Hz

# Synthesize a pure tone
x = scale_dbspl(pure_tone(freq, phase, dur, fs), level);

# Simulate IHC response at several CFs
results = sim_anrate_zbc2014(sim_ihc_zbc2014(x, cfs), cfs)

# Plot
heatmap(
    t[1:3000], cfs, results[:, 1:3000]; 
    xlabel="Time (s)", 
    ylabel="CF (Hz)",
    yscale=:log10,
)
```

## Simulating hearing loss

We can extend our neurogram simulation above to simulate what happens in the case of broad loss of outer hair cells by passing a new value to the `cohc` parameter in [`sim_ihc_zbc2014`](@ref).
We'll do something even a bit fancier --- by wrapping the entire routine in a function, we can easily repeat the simulation for various levels of `cohc` and compare them side-by-side.

```@example
using AuditoryNerveFiber
using AuditorySignalUtils
using Plots

function viz_hearing_loss(cohc)
    # Define variables
    freq = 1000.0   # freq, Hz
    phase = 0.0     # starting phase, rads
    dur = 0.2       # duration, seconds
    fs = 100e3       # sampling rate, Hz
    level = 50.0    # level, dB SPL
    t = 0:(1/fs):prevfloat(dur)  # time axis, s
    cfs = LogRange(200.0, 20000.0, 100)  # CFs, Hz

    # Synthesize a pure tone
    x = scale_dbspl(pure_tone(freq, phase, dur, fs), level);

    # Simulate IHC response at several CFs
    results = sim_anrate_zbc2014(sim_ihc_zbc2014(x, cfs; cohc=cohc), cfs)  # here we pass cohc to sim_ihc_zbc2014

    # Plot
    heatmap(
        t[1:3000], cfs, results[:, 1:3000]; 
        xlabel="Time (s)", 
        ylabel="CF (Hz)",
        yscale=:log10,
        clim=(0, 1000),
        title="OHC loss = $((1-cohc)*100)%",
    )
end

# Run function and plot results
plots = map(viz_hearing_loss, [1.0, 0.5, 0.1, 0.01])
plot(plots...; layout=(4, 1), size=(500, 300*4))
```

As can be seen, with significant OHC loss the response becomes weaker over all (in terms of average firing rate) and becomes more broadly tuned.