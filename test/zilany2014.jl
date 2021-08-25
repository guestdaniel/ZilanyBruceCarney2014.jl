using Test
include("../src/ZBC2014.jl")
include("../src/Hearing.jl")

# Declare various constants that hold across all tests in this file
fs = 100000.0
dur = 1.0
pt = Hearing.scale_dbspl(Hearing.pure_tone(1000.0, 0.0, dur, fs), 50.0)

# Start by testing the direct bindings and just make sure that they run!
@testset "C bindings: check callable" begin
  # First, we try testing the direct C binding to the IHC simulation with silence to make sure it runs
  @test begin
      px = pt
      cf = 1000.0
      nrep = Int32(1)
      tdres = 1.0/fs
      totalstim = Int32(dur*fs)
      cohc = 1.0
      cihc = 1.0
      species = Int32(1)
      ihcout = Vector{Cdouble}(zeros((Int64(dur*fs), )))
      ZBC2014.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      true
  end

  # Next, we try testing the direct C binding to the simulation with silence to make sure it runs
  @test begin
      px = pt
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

      ZBC2014.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      ZBC2014.Synapse!(ihcout, tdres, cf, totalstim, nrep, spont, noiseType, implnt,
                       sampFreq, synouttmp)
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

      ZBC2014.IHCAN!(px, cf, nrep, tdres, totalstim, cohc, cihc, species, ihcout)
      ZBC2014.SingleAN(ihcout, cf, nrep, tdres, totalstim, fibertype, noiseType,
                       implnt, meanrate, varrate, psth)
      true
  end
end

# Next we test the "wrappers" that provide a nicer interface to the bindings, beginning with simple checks that they are callable
@testset "Wrappers: check callable" begin
    # Check if the wrapper can be evaluated
    @test begin
        ZBC2014.sim_ihc_zbc2014(zeros((Int(dur*fs), )), 1000.0)
        true
    end
end
# Now we check that they handle arguments correctly
@testset "Wrappers: arg handling" begin
    # Check that multiple dispatch correctly handles multiple cf requests
    @test begin
        output = ZBC2014.sim_ihc_zbc2014(zeros((Int(dur*fs), )), Array{Float64}([1000.0, 2000.0]))
        size(output) == ((Int(dur*fs), 2))
    end
    # Check that multiple dispatch correctly handles multiple cf requests and the outputs actually contain something
    @test begin
        input = pt
        output = ZBC2014.sim_ihc_zbc2014(input, [1000.0, 2000.0])
        any(output .!= 0)
    end
end
# Next we check that the outputs are sensible
@testset "Wrappers: outputs" begin
    # Check that sinusoid input produces sensible response
    @test begin
        input = pt
        output = ZBC2014.sim_ihc_zbc2014(input, 1000.0)
        # Check that
        #1 at least one element is non-zero,
        #2 minima is smaller abs val than maxima (consistent with partial rectification in IHC)
        any(output .!= 0) && (abs(minimum(output)) < abs(maximum(output)))
    end
end
