# Just a fun little script to create some gifs for the README page
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using Plots

# ===== Plot #1 ======================================
# Plot ANF response as function of level for pure tone
# ====================================================
# Define range of levels
levels = LinRange(0.0, 60.0, 30)

# Define function that plots what we want
function plot_an(level)
    x = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, 0.1, 100_000.0), level)
    x = [zeros(500); x; zeros(500)]
    y = ANF.sim_an_zbc2014(x, 1000.0)[1]
    plot(y[250:2000], ylims=(0, 1500))
    annotate!((250, 1250, string(round(level; sigdigits=2))))
end

# Create gif
anim = @animate for level in levels
    plot_an(level)
end
gif(anim, fps=10)