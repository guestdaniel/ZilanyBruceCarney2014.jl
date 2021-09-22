# The goal of this script is to replicate key figures from Zilany et al. (2009) and 
# Zilany, Bruce, and Carney (2014) to verify the correct functioning of the 
# Zilany-Bruce-Carney model.
using Distributed
addprocs(3)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@assert nprocs() == 4
@everywhere using AuditoryNerveFiber
@everywhere const ANF = AuditoryNerveFiber
@everywhere using AuditorySignalUtils
@everywhere const ASU = AuditorySignalUtils
@everywhere using Statistics
using Plots

"""
    zilanyetal2014_figure1()

This function replicates the right-hand-side of Figure 1 from Zilany, Bruce, and Carney (2014). 
"""
function zilanyetal2014_figure1(resolution=50)
    # Define function to synthesize stimulus at one particular freq / level
    @everywhere function synthesize_pure_tone(freq, level)
        # Synthesize pure_tone
        stim = ASU.pure_tone(freq, 0.0, 0.050, 100e3)
        # Scale pure pure 
        stim = ASU.scale_dbspl(stim, level)
        # Apply 2.5 ms cosine ramp
        stim = ASU.cosine_ramp(stim, 0.0025, 100e3)
        # Return
        return stim
    end

    # Define function that implements model 
    @everywhere function simulate_an_response(stim, cf)
        x = mean(ANF.sim_synapse_zbc2014(ANF.sim_ihc_zbc2014([zeros(500); stim; zeros(500)], cf), cf))
        #x = x ./ (1 .+ 0.75e-3 * x)  # transform synapse out into firing rate
        return x
    end

    # Define ranges over which we'll iterate
    cfs = ASU.LogRange(200, 20000, resolution)
    levels = LinRange(-15, 120, resolution)
    domain = [(cf, level) for cf in cfs, level in levels]

    # Get responses
    avg_rates = reshape(pmap(x -> simulate_an_response(synthesize_pure_tone(500.0, x[2]), x[1]), hcat(domain...)), (resolution, resolution))

    # Plot
    #heatmap(cfs, levels, transpose(avg_rates), xscale=:log10, color=:jet, clim=(50, 320))
    heatmap(cfs, levels, transpose(avg_rates), xscale=:log10, color=:jet)
    xlabel!("CF (Hz)")
    ylabel!("Level (dB SPL)")

    # Save figure
    savefig("test/outputs/zilanyetal2014_figure1.png")
end


"""
    zilanyetal2009_figurex()
"""
function zilanyetal2009_figure4(n_rep=20)
    # Define function to synthesize stimulus
    @everywhere function synthesize_pure_tone(dur)
        # Synthesize pure tone of duration dur
        stim = ASU.pure_tone(10_000.0, 0.0, dur, 100e3)
        # Set to level
        stim = ASU.scale_dbspl(stim, 15.0)
        # Set ramp
        stim = ASU.cosine_ramp(stim, 0.0025, 100e3)
        # Zero-pad
        stim = [zeros(500); stim; zeros(20_000)]
        # Return
        return stim
    end
    # Define function to run simulation
    @everywhere function simulate_an_response(dur, n_rep)
        # Synthesize stimulus 
        stim = synthesize_pure_tone(dur)
        # Run simulation
        x = ANF.sim_synapse_zbc2014(ANF.sim_ihc_zbc2014(stim, 10_000.0; n_rep=n_rep), 10_000.0; frac_noise="actual")
        x = x ./ (1 .+ 0.75e-3 * x)  # transform synapse into firing rate
        x = reshape(x, (Int(500 + dur*100e3 + 20e3), n_rep))
        # Average over resp
        x = mean(x; dims=2)
        # Return zero-padded to 1.2 seconds
        return [x; zeros(Int(1.5*100e3 - length(x)))]
    end
    # Run simulation for many levels
    results = pmap(dur -> simulate_an_response(dur, n_rep), [0.050, 0.10, 0.20, 0.50, 1.0])

    # Visualize
    t = LinRange(0, 1.5, Int(100e3*1.5))
    plot(t, results, ylim=(0, 500), label=hcat([0.050, 0.10, 0.20, 0.50, 1.0]...), legendtitle="Duration (s)")
    hline!([75], color="black", label="")
    ylabel!("Firing rate (sp/s)")
    xlabel!("Time (s)")

    # Save figure
    savefig("test/outputs/zilanyetal2009_figure4.png")
end

zilanyetal2009_figure4(20)
zilanyetal2014_figure1(40)