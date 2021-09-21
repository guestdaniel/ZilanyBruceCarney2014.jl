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

# Zilany et al. (2009): Figure 3

"""
    zilanyetal2014_figure1()

This function replicates the right-hand-side of Figure 1 from Zilany, Bruce, and Carney (2014). 
"""
function zilanyetal2014_figure1(resolution=50)
    # Define function to synthesize stimulus at one particular freq / level
    @everywhere function stim(freq, level)
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
    @everywhere function resp(stim, cf)
        return mean(ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014([zeros(500); stim; zeros(500)], cf), cf)[1])
    end

    # Define ranges over which we'll iterate
    cfs = ASU.LogRange(200, 20000, resolution)
    levels = LinRange(-15, 120, resolution)
    domain = [(cf, level) for cf in cfs, level in levels]

    # Get responses
    avg_rates = reshape(pmap(x -> resp(stim(500.0, x[2]), x[1]), hcat(domain...)), (resolution, resolution))

    # Plot
    heatmap(cfs, levels, transpose(avg_rates), xscale=:log10, color=:jet, clim=(50, 320))
end