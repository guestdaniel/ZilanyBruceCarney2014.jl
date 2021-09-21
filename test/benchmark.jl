using BenchmarkTools
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using BenchmarkTools

# Declare various constants that hold across all tests in this file
fs = 100_000.0
dur = 1.0
pt = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, dur, fs), 50.0)

# Benchmark running a sinusoid several times
@benchmark ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0)
@benchmark map(cf -> ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), cf), ASU.LogRange(500.0, 5000.0, 50))