using Test
include("../src/Hearing.jl")
include("../external/util.jl")

# Declare various constants that hold across all tests in this file
fs = 100000
dur = 1.0

# First, we try testing the direct C binding to the IHC simulation with silence to make sure it runs
@test begin
    #px = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    px = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)
    cf = 1000.0
    nrep = Int32(1)
    tdres = 1.0/fs
    totalstim = Int32(dur*fs)
    cohc = 1.0
    cihc = 1.0
    species = Int32(1)
    ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    try
        Hearing.sim_ihc(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
        true
    catch e
        false
    end
end

# Next, we try testing the direct C binding to the simulation with silence to make sure it runs
@test begin
    ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    tdres = 1.0/fs
    cf = 1000.0
    totalstim = Int32(dur*fs)
    nrep = Int32(1)
    spont = 1.0
    noiseType = 1.0
    implnt = 1.0
    sampFreq = 20000.0
    synouttmp = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    myffGn = @cfunction(ffGn, Vector{Cdouble}, (Cint, ))
    try
        Hearing.sim_an(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp, myffGn)
        true
    catch e
        false
    end
end


@test begin
    x = 1.0
    slope = 1.0
    asym = 1.0
    cf = 1000.0
    try
        Hearing.NLogarithm(x, slope, asym, cf)
        true
    catch e
        false
    end
end
