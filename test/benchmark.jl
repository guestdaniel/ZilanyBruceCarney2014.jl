using BenchmarkTools
using ANF
using AuditorySignalUtils
const ASU = AuditorySignalUtils

# Declare various constants that hold across all tests in this file
fs = 100000.0
dur = 1.0
pt = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, dur, fs), 50.0)

# Benchmark running a sinusoid several times
@benchmark ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0)