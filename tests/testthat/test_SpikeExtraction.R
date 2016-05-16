library(relectro)

context("spike extraction")
test_that("spike extraction",{
  ## get a raw trace to work with
  sim<-simulateRawTrace(samplingRate=20000,
                        durationSec=20,
                        noiseSD=100,
                        noiseMean=0,
                        waveformAmplitude=2000, ## very clear wave forms
                        nClusters=3,
                        nChannels=4,
                        waveformDifferentiationSD=100,
                        maxSpikes=10000)
  resD<-detectSpikesTetrodes(data=sim$trace,
                             samplingRate=20000,
                             powerWindowSizeMs=0.5,
                             powerWindowSlideMs=0.1,
                             SDThreshold=3, ## high threshold
                             simultaneousSpikeMaxJitterMs=0.4)
  da<-spikeDetectionAccuracy(resD,sim$spikeTime,maxJitter=1)
  ## all spikes should be detected on only true spikes
  expect_equal(da$detectedTrue,length(sim$spikeTime))
  expect_equal(da$trueDetected,length(sim$spikeTime))
  expect_equal(da$detectedFalse,0)
  expect_equal(da$trueNonDetected,0)
  rm(da,sim,resD)
})
