# The goal of this script is to replicate key figures from Zilany et al. (2009)
# and Zilany, Bruce, and Carney (2014) to verify the correct functioning of the
# Zilany-Bruce-Carney model.
using Distributed
if nworkers() < 10
    addprocs(10 - nworkers())
end
@everywhere begin
    using Pkg
    Pkg.activate(Base.current_project())
end
@everwhere begin
    using ZilanyBruceCarney2014
    using AuditorySignalUtils
    using Statistics
end
using CairoMakie

### 2014 - Figure 1
"""
    zilanyetal2014_figure1()

Replicate the right-hand-side of Figure 1 from Zilany, Bruce, and Carney (2014).

The right-hand-side of Figure 1 from Zilany, Bruce, and Carney (2014) depicts
average rates from a cat ANF simulation at a range of sound levels and
frequencies for an input pure tone with a frequency of 1000 Hz, a duration of
50 ms, and 2.5-ms ramps. Here, we generate this figure by simulating many
repeated spiking responses.
"""
function zilanyetal2014_figure1(;
    n_cf=50,
    n_level=50,
    dur=0.05,
    n_rep=50,
    freq=500.0,
)
    # Define ranges over which we'll iterate
    cfs = LogRange(0.2e3, 20e3, n_cf)
    levels = LinRange(-15.0, 120.0, n_level)

    # Get responses
    avg_rates = pmap(Base.Iterators.product(cfs, levels)) do (cf, level)
        # Synthesize pure_tone
        stim = pure_tone(freq, 0.0, dur, 100e3)
        # Apply 2.5 ms cosine ramp
        stim = cosine_ramp(stim, 0.0025, 100e3)
        # Scale pure pure
        stim = scale_dbspl(stim, level)
        # Simulate response
        resp = sim_an_zbc2014(sim_ihc_zbc2014(stim, cf; species="cat", n_rep=n_rep), cf; n_rep=n_rep, fractional=true)[3]
        sum(resp) / dur / n_rep
    end

    # Plot
    fig = Figure(resolution=(500, 800))
    ax = Axis(fig[1, 1]; xscale=log10)
    ax.xticks = [0.5, 1.0, 2.0, 5.0, 10.0]
    ax.yticks = 0.0:20.0:120.0
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "SPL (dB)"

    hm = heatmap!(ax, cfs ./ 1000, levels, avg_rates; colormap=:jet1, colorrange=(50, 310))
    Colorbar(fig[1, 2], hm)

    # Save figure
    save("/home/daniel/ZilanyBruceCarney2014.jl/test/outputs/zilanyetal2014_figure1.png", fig)
end
zilanyetal2014_figure1()

### 2009 - Figure 4
"""
    zilanyetal2009_figure3()
"""
function zilanyetal2009_figure3(
    n_fiber=738,
    dur=30.0,
    cf=4000.0,
    type_proportions=[0.61, 0.23, 0.16],
    fs=100e3,
)
    # Figure out how many of each fiber type we're going to run
    n_fiber_per_type = Int.(floor.(type_proportions .* n_fiber))

    # Create vector to iterate over
    fiber_types = vcat([repeat([t], n) for (t, n) in zip(["high", "medium", "low"], n_fiber_per_type)]...)

    # Run simulation for many stim durations
    results = pmap(fiber_types) do fiber_type
        # Run simulation
        spikes = sim_spikes_zbc2014(sim_ihc_zbc2014(zeros(Int(floor(dur*fs))), cf; species="cat"), cf; fractional=true, fiber_type=fiber_type)
        return sum(spikes)/dur
     end

    # Visualize
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.ylabel = "Number of units"
    ax.xlabel = "Spontaneous rate (sp/s)"
    xlims!(ax, (0.0, 120.0))
    ylims!(ax, (0, 150))
    hist!(ax, results; color=:black, bins=120)

    # Save figure
    save("/home/daniel/ZilanyBruceCarney2014.jl/test/outputs/zilanyetal2009_figure3.png", fig)
end

zilanyetal2009_figure4()

### 2009 - Figure 4
"""
    zilanyetal2009_figure4()
"""
function zilanyetal2009_figure4(
    n_rep=50,
    cf=10000.0,
    fs=100e3,
    level=12.0,
)
    # Run simulation for many stim durations
    durations = [0.050, 0.10, 0.20, 0.50, 1.0]
    results = pmap(durations) do dur
        # Synthesize stimulus
        stim = pure_tone(cf, 0.0, dur, fs)
        stim = cosine_ramp(stim, 0.0025, fs)
        stim = scale_dbspl(stim, level)
        stim = [stim; zeros(20000)]

        # Run simulation
        x = sim_synapse_zbc2014(sim_ihc_zbc2014(stim, cf; n_rep=n_rep, species="cat"), cf; power_law="actual")
        x = x ./ (1 .+ 0.75e-3 .* x)  # transform synapse into firing rate
        x = reshape(x, (Int(dur*fs + 20000), n_rep))

        # Average over resp
        x = mean(x; dims=2)
        x = dropdims(x; dims=2)

        # Return zero-padded to 1.5 seconds
        return [x; zeros(Int(1.5*fs - length(x)))]
    end

    # Visualize
    t = LinRange(0, 1.5, Int(fs*1.5))
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.ylabel = "Firing rate (sp/s)"
    ax.xlabel = "Time (s)"
    xlims!(ax, (0.0, 1.2))
    ylims!(ax, (0, 500))
    lns = map(1:length(results)) do i
        lines!(ax, t, results[i])
    end
    hlines!(ax, [75]; color=:black)
    axislegend(ax, lns, string.(durations); position=:rt)

    # Save figure
    save("/home/daniel/ZilanyBruceCarney2014.jl/test/outputs/zilanyetal2009_figure4.png", fig)
end

zilanyetal2009_figure4()
