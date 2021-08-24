using Test
include("../src/Hearing.jl")
include("../external/util.jl")

# Declare various constants that hold across all tests in this file
fs = 200000
dur = 1.0

# First, we try to make sure the basic util.jl functions work correctly
#@test begin
#    tone_orig = Hearing.pure_tone(1000.0, 0.0, 1.0, 48000)
#    tone_down = zeros((24000, ))
#    decimate(pointer(tone_orig), pointer(tone_down), Int32(48000), Int32(2))
#    size(tone_down) == (24000,)
#end

# First, we try testing the direct C binding to the IHC simulation with silence to make sure it runs
@test begin
    px = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)
    cf = 1000.0
    nrep = Int32(1)
    tdres = 1.0/fs
    totalstim = Int32(dur*fs)
    cohc = 1.0
    cihc = 1.0
    species = Int32(1)
    ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    Hearing.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
    true
end

# Next, we try testing the direct C binding to the simulation with silence to make sure it runs
@test begin
    px = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)
    tdres = 1.0/fs
    cf = 1000.0
    totalstim = Int32(dur*fs)
    nrep = Int32(1)
    spont = 100.0
    noiseType = 1.0
    implnt = 0.0
    sampFreq = 50000.0
    cohc = 1.0
    cihc = 1.0
    species = Int32(1)
    ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    synouttmp = Vector{Cdouble}(zeros((Int64(dur*fs), )))

    myffGn = @cfunction(ffGn, Vector{Cdouble}, (Cint, ))
    mydecimate = @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint))

    Hearing.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
    Hearing.Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp, myffGn, mydecimate)
    true
end

# Finally, we try testing the direct C binding to the SingleAN() function 
@test begin
    px = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)
    tdres = 1.0/fs
    cf = 1000.0
    totalstim = Int32(dur*fs)
    nrep = Int32(1)
    spont = 100.0
    fibertype = 3.0
    noiseType = 1.0
    implnt = 0.0
    sampFreq = 50000.0
    cohc = 1.0
    cihc = 1.0
    species = Int32(1)
    ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    synouttmp = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    meanrate = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    varrate = Vector{Cdouble}(zeros((Int64(dur*fs), )))
    psth = Vector{Cdouble}(zeros((Int64(dur*fs), )))

    myffGn = @cfunction(ffGn, Vector{Cdouble}, (Cint, ))
    mydecimate = @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint))

    Hearing.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
    Hearing.SingleAN(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType, implnt, meanrate, varrate, psth, myffGn, mydecimate)
    true
end
