module BrainWaveExperimental 

export calcDelayedPeakPoints, calcPeaksFromExpectedLatencies, calcPeaksInflFromExpectedLatencies, calcPeaksTroughsInflFromExpectedLatencies, findABRWaves, selectLargestPeakInWindow, selectLargestTroughInWindow, selectStrongestInflectionInWindow

using BrainWave, DataFrames, DocStringExtensions

include("findABRPeaks.jl")
include("utils.jl")

end #module
