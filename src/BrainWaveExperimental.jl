module BrainWaveExperimental 

export calcDelayedPeakPoints, calcPeaksFromExpectedLatencies, calcPeaksInflFromExpectedLatencies, calcPeaksInflFromExpectedLatencies2, calcPeaksInflFromExpectedLatencies3, calcPeaksTroughsInflFromExpectedLatencies,
findABRWaves, findInflectionsDown, findInflectionsUp,
selectLargestPeakInWindow, selectLargestTroughInWindow, selectStrongestInflectionInWindow

using BrainWave, DataFrames, DocStringExtensions

include("findABRPeaks.jl")
include("utils.jl")

end #module
