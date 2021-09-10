# Just a fun little script to create some gifs for the README page
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using Plots
using Images

# ===== Plot #1 ======================================
# Plot iso-level tuning curves for pure-tone
# ====================================================
fs = 10e4                               # sampling rate (Hz)
cfs = ASU.LogRange(200.0, 5000.0, 100)  # characteristic frequencies (Hz)
levels = LinRange(10.0, 50.0, 5)        # sound levels (dB SPL)

# Define function to calculate tuning curves
function tuning_curve(level, cfs)
    pt = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, 0.1, fs), level)
    return map(cf -> mean(ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(pt, cf), cf)[1]), cfs)
end

# Estimate tuning curves for each level
tcs = map(level -> tuning_curve(level, cfs), levels)
labels = string.(levels)

# Plot
plot(cfs, hcat(tcs...), label=hcat(labels...), legendtitle="Level (dB SPL)", ylabel="Firing rate (sp/s)", xlabel=("CF (Hz)"))