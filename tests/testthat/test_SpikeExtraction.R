library(relectro)

context("spike extraction")
test_that("spike extraction",{
  ## get a raw trace to work with
  sim<-simulateRawTrace(samplingRate=20000,
                        durationSec=2,
                        noiseSD=100,
                        noiseMean=0,
                        waveformAmplitude=1500,
                        nClusters=3,
                        waveformDifferentiationSD=200,
                        maxSpikes=10000)
  spD<-detectSpikesFromTrace(data=sim$trace,
                             samplingRate=20000,
                             powerWindowSizeMs=0.4,
                             powerWindowSlideMs=0.1,
                             SDThreshold=2.5)
  da<-spikeDetectionAccuracy(spD$spikeTime,sim$spikeTime,maxJitter=1)
  
  ## all spikes should be detected on only true spikes
  expect_equal(da$detectedTrue,length(sim$spikeTime))
  expect_equal(da$trueDetected,length(sim$spikeTime))
  rm(da,sim,spD)
})
