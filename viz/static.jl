# Just a fun little script to create some gifs for the README page
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using Plots
using Images

# Plot 1
fs = 10e4
pt = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, 0.1, fs), 30.0)
plot(map(cf -> mean(ANF.sim_ihc_zbc2014(pt, cf)), ASU.LogRange