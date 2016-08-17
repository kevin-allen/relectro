library(relectro)
library(testthat)
context("spike extraction")
test_that("spike extraction",{
  ## get a raw trace to work with
  sim<-simulateRawTrace(samplingRate=20000,
                        durationSec=0.5,
                        noiseSD=100,
                        noiseMean=0,
                        waveformAmplitude=1500, ## very clear wave forms, whould always be 100 % correct in test
                        nClusters=3,
                        nChannels=4,
                        waveformDifferentiationSD=100,
                        maxSpikes=10000)
  resD<-detectSpikesTetrodes(data=sim$trace,
                             samplingRate=20000,
                             powerWindowSizeMs=0.5,
                             powerWindowSlideMs=0.1,
                             SDThreshold=2.5, ## high threshold
                             simultaneousSpikeMaxJitterMs=0.4,
                             spikeDetectionRefractoryMs = 0.5)
  ## should detect all and ony real spikes
  expect_equal(length(resD$spikeTimes),length(sim$spikeTime))
  ## the spike times should be similar to real spike time,
  ## difference is caused by jitter and noise in waveforms.
  expect_lt(max(abs(resD$spikeTimes-sim$spikeTime)),4)
  rm(sim,resD)
})
