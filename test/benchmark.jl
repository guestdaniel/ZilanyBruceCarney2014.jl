using BenchmarkTools
using Profile
using ProfileView
include("../src/ZBC2014.jl")
include("../src/Hearing.jl")

# Declare various constants that hold across all tests in this file
fs = 100000.0
dur = 1.0
pt = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)

# Benchmark running a sinusoid several times
@benchmark ZBC2014.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0)

@time ZBC2014.sim_ihc_zbc2014(pt, repeat([1000.0], 100));
@profile ZBC2014.sim_ihc_zbc2014(pt, repeat([1000.0], 100));

input = pt
species_flag = ZBC2014.translate_species_label("cat")
cohc = 1.0
cihc = 1.0
output = zeros((length(input), ))
cf = 1000.0
@time ZBC2014.IHCAN!(input, cf, Int32(1), 1/fs, Int32(length(input)), cohc, cihc, Int32(species_flag), output)

