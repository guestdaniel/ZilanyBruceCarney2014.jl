"""
This testing suite handles very basic tests of the implemented auditory nerve models... Namely, it makes sure that they can be run without error using
default arguments and empty/standard inputs and it makes sure that basic response properties (e.g., spontaneous rate and rate-level function) are
approximately correct (i.e., exact spontaneous rates are not checked, but relative spontaneous rates like high > low are checked).
"""

using Test
using AuditoryNerveFiber
const ANF = AuditoryNerveFiber
using AuditorySignalUtils
const ASU = AuditorySignalUtils
using Statistics
using DSP

# Declare various constants that hold across all tests in this file
fs = 100_000.0
dur = 0.1
freq = 1000.0
pt = ASU.scale_dbspl(ASU.pure_tone(freq, 0.0, dur, fs), 50.0)
tol = 1e-2  # general tolerance on comparisons of approximate equality (should upgrade to isapprox)

# Test ffGn
@testset "Fractional Gaussian noise" begin
    # Test that we can call the function
    @test begin
        sample = ANF.ffGn(Int32(10000), 1/fs, 0.75, 1.0, 1.0)
        true
    end
    # Test that the noiseType switch behaves as expect
    @test begin
        sample = ANF.ffGn(Int32(10000), 1/fs, 0.75, 0.0, 1.0)
        sample = unsafe_wrap(Array, sample, 10000)
        all(sample .== 0.0)
    end
    # Test that raising the mean to the branch points in the code (0.5, 18.0 results in 
    # corresponding changes in sigma
    @test begin
        myfunc(mu) = std(unsafe_wrap(Array, ANF.ffGn(Int32(100000), 1/fs, 0.75, 1.0, mu), 100000))
        vars = map(myfunc, [0.4, 10.0, 20.0])
        vars[1] < vars[2] < vars[3]
    end
end

# Test upsampling function
@test begin
    pt_upsampled = ANF.upsample(pt, 5)
    maximum(abs.(pt_upsampled[1:5:length(pt_upsampled)] - pt)) < tol
end

# Start by testing the direct bindings and just make sure that they run!
@testset "C bindings: check callable" begin
  # First, we try testing the direct C binding to the IHC code (IHCAN!)
  @test begin
      px = pt
      cf = freq
      nrep = Int32(1)
      tdres = 1.0/fs
      totalstim = Int32(dur*fs)
      cohc = 1.0
      cihc = 1.0
      species = Int32(1)
      ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
      ANF.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      true
  end

  # Next, we try testing the direct C binding to the Synapse code (Synapse!)
  @test begin
      px = pt
      tdres = 1.0/fs
      cf = freq
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

      ANF.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      ANF.Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt,
                       sampFreq, synouttmp)
      true
  end

  # Finally, we try testing the direct C binding to the SingleAN!() function
  @test begin
      px = ASU.scale_dbspl(ASU.pure_tone(freq, 0.0, dur, fs), 50.0)
      tdres = 1.0/fs
      cf = freq
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

      ANF.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      ANF.SingleAN!(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType,
                       implnt, meanrate, varrate, psth)
      true
  end
end

# Next we test the "wrappers" that provide a nicer interface to the bindings, beginning with simple checks that they are callable
@testset "Wrappers: check callable" begin
    # Check if the wrapper can be evaluated
    @test begin
        ANF.sim_ihc_zbc2014(zeros((Int(dur*fs), )), freq)
        true
    end
    @test begin
        ANF.sim_synapse_zbc2014(zeros((Int(dur*fs), )), freq)
        true
    end
    @test begin
        ANF.sim_an_zbc2014(zeros((Int(dur*fs), )), freq)
        true
    end
end

# Next we check that the outputs are sensible
@testset "Wrappers: outputs" begin
    # Check that inner hair cell response shows expected hallmarks (response to sinusoid at CF, partial rectification)
    @test begin
        output = ANF.sim_ihc_zbc2014(pt, freq)
        # Check that
        #1 at least one element is non-zero,
        #2 minima is smaller abs val than maxima (consistent with partial rectification in IHC)
        any(output .!= 0) && (abs(minimum(output)) < abs(maximum(output)))
    end
    # Check that each output of auditory nerve response is present and non-zero
    @test begin
        outs = ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(pt, freq), freq)
        (sum(outs[1]) != 0) && (sum(outs[2]) != 0) && (sum(outs[3]) != 0)
    end
    # Check that auditory nerve response shows expected pattern of results with changing spontaneous rate
    @test begin
        # Check analytic firing rates
        (low, medium, high) = 
            map(fiber_type -> mean(ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(pt, freq), freq; fiber_type=fiber_type)[1]), ["low", "medium", "high"])
        analytic = (low < medium < high)
        # Check spiking firing rate (note that this approach is not great because it's stochastic and could fail randomly)
        (low, medium, high) = 
            map(fiber_type -> mean(ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(pt, freq), freq; fiber_type=fiber_type)[3]), ["low", "medium", "high"])
        spiking = (low < medium < high)
        analytic && spiking
    end
    # Check that auditory nerve shows reasonable rate-level function (e.g., range from 10 to 30 dB is constantly increasing)
    @test begin
        inputs = [ASU.scale_dbspl(pt, level) for level in [10.0, 15.0, 20.0, 25.0, 30.0]]
        ihcouts = [ANF.sim_ihc_zbc2014(input, freq) for input in inputs]
        outputs = [ANF.sim_an_zbc2014(ihcout, freq)[1] for ihcout in ihcouts]
        means = map(Statistics.mean, outputs)
        means[5] > means[4] > means[3] > means[2] > means[1]
    end
    # Check that IHC shows reasonable frequency tuning
    @test begin
        cfs = ASU.LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(ANF.sim_ihc_zbc2014(pt, cf)), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
    # Check that synapse shows reasonable frequency tuning
    @test begin
        cfs = ASU.LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(ANF.sim_synapse_zbc2014(ANF.sim_ihc_zbc2014(pt, cf), cf)), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
    # Check that auditory nerve shows reasonable frequency tuning
    @test begin
        cfs = ASU.LogRange(200.0, 5000.0, 40)
        tuning_curve = map(cf -> mean(ANF.sim_an_zbc2014(ANF.sim_ihc_zbc2014(pt, cf), cf)[1]), cfs)
        abs(cfs[findmax(tuning_curve)[2]] - freq) < 500
    end
end
