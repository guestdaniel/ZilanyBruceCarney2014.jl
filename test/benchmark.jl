using BenchmarkTools
using Profile
using ZilanyBruceCarney2014
const ANF = ZilanyBruceCarney2014
using AuditorySignalUtils
const ASU = AuditorySignalUtils

# Declare various constants that hold across all tests in this file
fs = 100_000.0
dur = 1.0
pt = ASU.scale_dbspl(ASU.pure_tone(1000.0, 0.0, dur, fs), 50.0)

# Benchmark running a sinusoid several times
@benchmark ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0)
@benchmark ANF.sim_synapse_zbc2014(zeros((Int(dur*fs), )), 1000.0)
@benchmark ANF.sim_an_zbc2014(zeros((Int(dur*fs), )), 1000.0)
@benchmark ANF.sim_synapse_zbc2014(ANF.sim_ihc_zbc2014(zeros(Int(dur*fs), ), 1000.0), 1000.0)


# Profile running a sinsuoid several times
Profile.clear()
@profile (for i in 1:100; ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0); end)
Profile.print()

Profile.clear()
@profile (for i in 1:100; ANF.sim_synapse_zbc2014(zeros((Int(dur*fs), )), 1000.0); end)
Profile.print()