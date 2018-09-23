
using DocStringExtensions
"""
Select the largest peak among a set of peak candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
peak is found the function returns `missing`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the peak is sought.
* `peakPnts::AbstractVector{Real}`: vector with candidate peak points. Peak candidates can be found with the `findPeaks` function.
* `peakTimes::AbstractVector{Real}`: vector with candidate peak times. Peak candidates can be found with the `findPeaks` function.
* `minLat::Real`: minimum latency of the peak.
* `maxLat::Real`: maximum latency of the  peak.

##### Returns

* `peakPoint::Real`: index of the point at which the highest peak is found. If no peak is found `missing` is returned.
* `peakTime`: time at which the highest peak is found, in seconds. If no peak is found `missing` is returned.

##### Examples

```julia
    using BrainWave
    sampRate = 256; nTaps=64; minLat=0.5; maxLat=2.5
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    pkPnts, pkTimes = findPeaks(rec[1,:], sampRate)
    peakPnt, peakTime = selectLargestPeakInWindow(rec[1,:], pkPnts, pkTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Peak Candidates")
    ## s3 = scatter(;x=[tArr[peakPnt]], y=[rec[1, peakPnt]], mode="markers", marker_size=10, name="Largest Peak")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""
function selectLargestPeakInWindow(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, peakPnts::AbstractVector{P}, peakTimes::AbstractVector{S}, minLat::Real, maxLat::Real) where {T<:Real, P<:Real, S<:Real}

    peakCandidates = findall((peakTimes .<= maxLat ) .& (peakTimes .>= minLat))
    nCandidates = length(peakCandidates)
    if nCandidates == 0
        peakPnt = missing
        peakTime = missing
    elseif nCandidates == 1
        peakPnt = peakPnts[peakCandidates[1]]
        peakTime = peakTimes[peakCandidates[1]]
    else
        idx = findall(sig[peakPnts[peakCandidates]] .== maximum(sig[peakPnts[peakCandidates]]))[1]
        peakPnt = peakPnts[peakCandidates[idx]]
        peakTime = peakTimes[peakCandidates[idx]]
    end

    return peakPnt, peakTime

end

"""
Select the largest trough among a set of trough candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
trough is found the function returns `missing`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the trough is sought.
* `troughPnts::AbstractVector{Real}`: vector with candidate trough points. Trough candidates can be found with the `findTroughs` function.
* `troughTimes::AbstractVector{Real}`: vector with candidate trough times. Trough candidates can be found with the `findTroughs` function.
* `minLat::Real`: minimum latency of the trough.
* `maxLat::Real`: maximum latency of the  trough.

##### Returns

* `troughPoint::Real`: index of the point at which the highest trough is found. If no trough is found `missing` is returned.
* `troughTime`: time at which the highest trough is found, in seconds. If no trough is found `missing` is returned.

##### Examples

```julia
    sampRate = 256; nTaps=64; minLat=0.5; maxLat=2.5
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    trPnts, trTimes = findTroughs(rec[1,:], sampRate)
    troughPnt, troughTime = selectLargestTroughInWindow(rec[1,:], trPnts, trTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Trough Candidates")
    ## s3 = scatter(;x=[tArr[troughPnt]], y=[rec[1, troughPnt]], mode="markers", marker_size=10, name="Largest Trough")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""
function selectLargestTroughInWindow(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, troughPnts::AbstractVector{P}, troughTimes::AbstractVector{S}, minLat::Real, maxLat::Real, maxAmp::Real=Inf) where {T<:Real, P<:Real, S<:Real}

    troughCandidates = findall((troughTimes .<= maxLat ) .& (troughTimes .>= minLat))
    #make sure trough amplitude is smaller than peak amplitude
    troughCandidates = troughCandidates[findall(sig[troughPnts[troughCandidates]] .< maxAmp)]
    nCandidates = length(troughCandidates)
    if nCandidates == 0
        troughPnt = missing
        troughTime = missing
    elseif nCandidates == 1
        troughPnt = troughPnts[troughCandidates[1]]
        troughTime = troughTimes[troughCandidates[1]]
    else
        idx = findall(sig[troughPnts[troughCandidates]] .== minimum(sig[troughPnts[troughCandidates]]))[1]
        troughPnt = troughPnts[troughCandidates[idx]]
        troughTime = troughTimes[troughCandidates[idx]]
    end

    return troughPnt, troughTime

end

"""
Select the strongest inflection point among a set of inflection point candidates
(which must be passed as an argument to the function along with their time of occurrence),
that are within a minimum and maximum latency window. If no suitable
inflection point is found the function returns `missing`.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the inflection point is sought.
* `dy::AbstractVector{Real}`: first derivative of `sig`.
* `inflectionPnts::AbstractVector{Real}`: vector with candidate inflection points. Inflection point candidates can be found with the `findInflections` function.
* `inflectionTimes::AbstractVector{Real}`: vector with candidate inflection times. Inflection point candidates can be found with the `findInflections` function.
* `minLat::Real`: minimum latency of the inflection.
* `maxLat::Real`: maximum latency of the  inflection.

##### Returns

* `inflectionPoint::Real`: index of the point at which the highest inflection point is found. If no inflection point is found `missing` is returned.
* `inflectionTime`: time at which the highest inflection point is found, in seconds. If no inflection point is found `missing` is returned.

##### Examples

```julia
    sampRate = 256; nTaps=64; minLat=1; maxLat=2
    rec, evtTab = simulateRecording(nChans=1, dur=4, sampRate=sampRate)
    filterContinuous!(rec, sampRate, "lowpass", nTaps, [2], channels=[1], transitionWidth=0.2)
    inflPnts, inflTimes = findInflections(rec[1,:], sampRate)
    dy = zeros(size(rec[1,:])); dy[2:end] = diff(rec[1,:])
    inflectionPnt, inflectionTime = selectStrongestInflectionInWindow(rec[1,:], dy, inflPnts, inflTimes, minLat, maxLat)
    ## not run
    ## using PlotlyJS
    ## tArr = collect(1:length(rec[1,:]))/sampRate
    ## s1 = scatter(;x=tArr, y=rec[1,:], name="Waveform")
    ## s2 = scatter(;x=tArr[pkPnts], y=rec[1, pkPnts], mode="markers", name="Inflection Candidates")
    ## s3 = scatter(;x=[tArr[inflectionPnt]], y=[rec[1, inflectionPnt]], mode="markers", marker_size=10, name="Strongest Inflection")
    ## shapes = rect([minLat], [maxLat], [0], [1], fillcolor="gray", opacity=0.2, line_width=0, xref="x", yref="paper")
    ## plot([s1, s2, s3], Layout(shapes=shapes, xaxis_title="Time (s)", yaxis_title="Amplitude (a.u.)"))


```

"""
function selectStrongestInflectionInWindow(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, dy::Union{AbstractMatrix{Q}, AbstractVector{Q}}, inflPnts::AbstractVector{P}, inflTimes::AbstractVector{S}, minLat::Real, maxLat::Real; maxAmp::Real=Inf) where {T<:Real, Q<:Real, P<:Real, S<:Real}


    inflectionCandidates = findall((inflTimes .<= maxLat ) .& (inflTimes .>= minLat))
    #make sure inflection amplitude is smaller than peak amplitude
    inflectionCandidates = inflectionCandidates[findall(sig[inflPnts[inflectionCandidates]] .< maxAmp)]
    nCandidates = length(inflectionCandidates)
    if nCandidates == 0
        inflectionPnt = missing
        inflectionTime = missing
    elseif nCandidates == 1
        inflectionPnt = inflPnts[inflectionCandidates[1]]
        inflectionTime = inflTimes[inflectionCandidates[1]]
    else
        #strategy 1: find the minimum of signal (useful as surrogate of finding troughs)
        #idx = find(sig[inflPnts[inflectionCandidates]] .== minimum(sig[inflPnts[inflectionCandidates]]))[1]
        # strategy 2: find point at which first derivative is closer to zero (closer to a trough)
        idx = findall(abs.(dy[inflPnts[inflectionCandidates]]) .== minimum(abs.(dy[inflPnts[inflectionCandidates]])))[1]
        inflectionPnt = inflPnts[inflectionCandidates[idx]]
        inflectionTime = inflTimes[inflectionCandidates[idx]]
    end

    return inflectionPnt, inflectionTime

end

"""
Attempt to find peaks and troughs of the ABR response for a click of a given level. The algorithm is
largely based on Bradley and Wilson (2004).

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the ABR waveform for which the peaks and troughs are sought.
* `stimLevel::Real`: the level of the click used to evoke the ABR response.
* `sampRate::Real`: the sampling rate of the ABR recording.
* `epochStart::Real`: the time, in seconds, at which the epoch starts.
* `waveVLatencyOffset::Real`: additional wave V latency delay (e.g. if stimulus is different than a standard click you may want to specify an additional delay on top of the delay computed by the formula)
* `dBRef::String`: whether the stimulus level is specified in dB SPL `SPL`, or in dB normal hearing level `nHL`.

##### Returns

A dataframe with the following columns:

* `wave::String`: wave label.
* `peakPoint::Real`: index of the point at which the peak is detected in the waveform.
* `troughPoint::Real`: index of the point at which the peak is detected in the waveform.
* `peakLatency::Real`: latency of the peak, in seconds.
* `troughLatency::Real`: latency of the trough, in seconds.
* `peakAmp::Real`: amplitude of the peak, in microvolts.
* `troughAmp::Real`: amplitude of the trough, in microvolts.
* `peakTroughAmp::Real`: peak-to-trough amplitude, in microvolts.
* `minPeakLat::Real`: minimum peak latency, in seconds, used by the algorithm to find the peak.
* `maxPeakLat::Real`: maximum peak latency, in seconds, used by the algorithm to find the peak.
* `minTroughLat::Real`: minimum trough latency, in seconds, used by the algorithm to find the trough.
* `maxTroughLat::Real`: maximum trough latency, in seconds, used by the algorithm to find the trough.



##### References

* Bradley, A. P., & Wilson, W. J. (2004). Automated Analysis of the Auditory Brainstem Response. Proceedings of the 2004 Intelligent Sensors, Sensor Networks and Information Processing Conference, 2004., 541–546. http://doi.org/10.1109/ISSNIP.2004.1417519

* Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87.

##### Examples

```julia

```

"""
function findABRWaves(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, stimLevel::Real, sampRate::Real; epochStart::Real=0, waveVLatencyOffset::Real=0, dBRef::String="SPL") where {T<:Real}
    sig = vec(sig)
    if dBRef == "SPL"
        stimLevel = stimLevel - 31
    elseif dBRef == "nHL"
        stimLevel = stimLevel
    else
        error("`stimLevel` needs to be specified in `SPL` or `nHL`")
    end

    #equation derived from Prosser, S., & Arslan, E. (1987). Prediction of auditory brainstem wave V latency as a diagnostic tool of sensorineural hearing loss. Audiology : Official Organ of the International Society of Audiology, 26(3), 179–87. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/3662941
    avPeakVLat = (-5.859e-08*stimLevel^4 + 1.274e-05*stimLevel^3 - 0.0005424*stimLevel^2 - 0.05297*stimLevel + 9.3)/1000
    avPeakVLat = avPeakVLat + waveVLatencyOffset

    minPeakVLat = avPeakVLat-(0.00024*3)
    maxPeakVLat = avPeakVLat+(0.00024*3)
    #todo
    #shorten peak II trough II latency

    ##data from Picton book
    avPeakVpeakI_IPL = 0.00395
    avPeakVpeakIII_IPL = 0.00207

    minPeakVpeakI_IPL = avPeakVpeakI_IPL - (0.00023*2)
    maxPeakVpeakI_IPL = avPeakVpeakI_IPL + (0.00023*2)
    minPeakVpeakIII_IPL = avPeakVpeakIII_IPL - (0.00019*2)
    maxPeakVpeakIII_IPL = avPeakVpeakIII_IPL + (0.00019*2)

    ## guesstimates
    minPeakITroughILat = 0.25/1000
    maxPeakITroughILat = minPeakITroughILat+0.75/1000

    minPeakIIITroughIIILat = 0.25/1000
    maxPeakIIITroughIIILat = minPeakIIITroughIIILat+0.75/1000

    minPeakVTroughVLat = 0.25/1000
    maxPeakVTroughVLat = 1.75/1000

    minPeakIITroughIILat = 0.125/1000
    minPeakIVTroughIVLat = 0.125/1000

    # find all peaks, troughs and inflection points in the waveform
    peakPnts, peakTimes = findPeaks(sig, sampRate, epochStart=epochStart)
    troughPnts, troughTimes = findTroughs(sig, sampRate, epochStart=epochStart)
    inflPnts, inflTimes = findInflections(sig, sampRate, epochStart=epochStart)
    #dy = zeros(size(sig)); dy[:,2:size(sig)[2]] = diff(sig, 2)
    dy = zeros(size(sig)); dy[2:end] = diff(sig)


    peakVPnt, peakVTime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakVLat, maxPeakVLat)
    if ismissing(peakVPnt) == true
        peakIPnt = missing
        peakITime = missing
        peakIIPnt = missing
        peakIITime = missing
        peakIIIPnt = missing
        peakIIITime = missing
        peakIVPnt = missing
        peakIVTime = missing

        troughIPnt = missing
        troughITime = missing
        troughIIPnt = missing
        troughIITime = missing
        troughIIIPnt = missing
        troughIIITime = missing
        troughIVPnt = missing
        troughIVTime = missing
        troughVPnt = missing
        troughVTime = missing
        minPeakILat = missing
        maxPeakILat = missing
        minPeakIILat = missing
        maxPeakIILat = missing
        minPeakIIILat = missing
        maxPeakIIILat = missing
        minPeakIVLat = missing
        maxPeakIVLat = missing

        minTroughILat = missing
        maxTroughILat = missing
        minTroughIILat = missing
        maxTroughIILat = missing
        minTroughIIILat = missing
        maxTroughIIILat = missing
        minTroughIVLat = missing
        maxTroughIVLat = missing
        minTroughVLat = missing
        maxTroughVLat = missing
    else
        minPeakILat = peakVTime - maxPeakVpeakI_IPL
        maxPeakILat = peakVTime - minPeakVpeakI_IPL
        minPeakIIILat = peakVTime - maxPeakVpeakIII_IPL
        maxPeakIIILat = peakVTime - minPeakVpeakIII_IPL

        minTroughVLat = peakVTime + minPeakVTroughVLat
        maxTroughVLat = peakVTime + maxPeakVTroughVLat

        #Peak I
        peakIPnt, peakITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakILat, maxPeakILat)
        #Peak III
        peakIIIPnt, peakIIITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIIILat, maxPeakIIILat)

        #Trough V
        troughVPnt, troughVTime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughVLat, maxTroughVLat, sig[peakVPnt])
        if ismissing(troughVPnt)
            troughVPnt, troughVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughVLat, maxTroughVLat, maxAmp=sig[peakVPnt])
        end

        #Trough I
        if ismissing(peakIPnt) == true
            troughIPnt = missing
            troughITime = missing
            minTroughILat = missing
            maxTroughILat = missing
        else
            minTroughILat = peakITime + minPeakITroughILat
            maxTroughILat = peakITime + maxPeakITroughILat
            troughIPnt, troughITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughILat, maxTroughILat, sig[peakIPnt])
            if ismissing(troughIPnt)
                troughIPnt, troughITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughILat, maxTroughILat, maxAmp=sig[peakIPnt])
            end
        end

        #Trough III
        if ismissing(peakIIIPnt) == true
            troughIIIPnt = missing
            troughIIITime = missing
            minTroughIIILat = missing
            maxTroughIIILat = missing
        else
            minTroughIIILat = peakIIITime + minPeakIIITroughIIILat
            maxTroughIIILat = peakIIITime + maxPeakIIITroughIIILat
            troughIIIPnt, troughIIITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIIILat, maxTroughIIILat, sig[peakIIIPnt])
            if ismissing(troughIIIPnt)
                troughIIIPnt, troughIIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIIILat, maxTroughIIILat, maxAmp=sig[peakIIIPnt])
            end
        end



        #Peak II
        minPeakIILat = missing
        maxPeakIILat = missing
        if ismissing(peakIIITime) == false
            maxPeakIILat = peakIIITime - 0.25/1000
        end
        if ismissing(troughITime) == false
            minPeakIILat = troughITime + 0.25/1000
        end
        if ismissing(maxPeakIILat) == true
            if ismissing(peakVTime) == false
                if ismissing(minPeakIILat) == false
                    maxPeakIILat = (peakVTime + minPeakIILat)/2.5
                end
            end
        end

        peakIIPnt=missing; peakIITime=missing
        if ((ismissing(minPeakIILat) == false) .& (ismissing(maxPeakIILat) == false))
            peakIIPnt, peakIITime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIILat, maxPeakIILat)
            if ismissing(peakIIPnt)
                peakIIPnt, peakIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minPeakIILat, maxPeakIILat)
            end
        end


        #Trough II
        minTroughIILat = missing
        maxTroughIILat = missing
        if ismissing(peakIIITime) == false
            maxTroughIILat = peakIIITime - 0.25/1000
        end
        if ismissing(peakIITime) == false
            minTroughIILat = peakIITime + minPeakIITroughIILat
        end

        troughIIPnt=missing; troughIITime=missing
        if ((ismissing(minTroughIILat) == false) .& (ismissing(maxTroughIILat) == false))
            troughIIPnt, troughIITime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIILat, maxTroughIILat, sig[peakIIPnt])
            if ismissing(troughIIPnt)
                troughIIPnt, troughIITime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIILat, maxTroughIILat, maxAmp=sig[peakIIPnt])
            end
        end

        #peak IV
        minPeakIVLat = missing
        maxPeakIVLat = missing
        if ismissing(peakVTime) == false
            maxPeakIVLat = peakVTime - 0.25/1000
        end
        if ismissing(troughIIITime) == false
            minPeakIVLat = troughIIITime + 0.25/1000
        end

        peakIVPnt=missing; peakIVTime=missing
        if ((ismissing(minPeakIVLat) == false) .& (ismissing(maxPeakIVLat) == false))
            peakIVPnt, peakIVTime = selectLargestPeakInWindow(sig, peakPnts, peakTimes, minPeakIVLat, maxPeakIVLat)
            if ismissing(peakIVPnt)
                peakIVPnt, peakIVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minPeakIVLat, maxPeakIVLat)
            end
        end



        #trough IV
        minTroughIVLat = missing
        maxTroughIVLat = missing
        if ismissing(peakVTime) == false
            maxTroughIVLat = peakVTime - 0.25/1000
        end
        if ismissing(peakIVTime) == false
            minTroughIVLat = peakIVTime + minPeakIVTroughIVLat
        end

        troughIVPnt=missing; troughIVTime=missing
        if ((ismissing(minTroughIVLat) == false) .& (ismissing(maxTroughIVLat) == false))
            troughIVPnt, troughIVTime = selectLargestTroughInWindow(sig, troughPnts, troughTimes, minTroughIVLat, maxTroughIVLat, sig[peakIVPnt])
            if ismissing(troughIVPnt)
                troughIVPnt, troughIVTime = selectStrongestInflectionInWindow(sig, dy, inflPnts, inflTimes, minTroughIVLat, maxTroughIVLat, maxAmp=sig[peakIVPnt])
            end
        end
    end

    peakPoints = [peakIPnt; peakIIPnt; peakIIIPnt; peakIVPnt; peakVPnt]
    peakLatencies = [peakITime; peakIITime; peakIIITime; peakIVTime; peakVTime]
    troughPoints = [troughIPnt; troughIIPnt; troughIIIPnt; troughIVPnt; troughVPnt]
    troughLatencies = [troughITime; troughIITime; troughIIITime; troughIVTime; troughVTime]

    peakAmps = missings(Float64, length(peakPoints))
    troughAmps = missings(Float64, length(peakPoints))
    #peakAmpsEstMax = zeros(length(peakPoints))
    #troughAmpsEstMax = zeros(length(peakPoints))
    for i=1:length(peakAmps)
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            #peakAmpsEstMax[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            #peakAmpsEstMax[i] = maximum(sig[,)
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
        end
    end


    peakTroughAmps = peakAmps-troughAmps

    waveLabels = ["I", "II", "III", "IV", "V"]

    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps,
                   minPeakLat=[minPeakILat, minPeakIILat, minPeakIIILat, minPeakIVLat, minPeakVLat],
                   maxPeakLat=[maxPeakILat, maxPeakIILat, maxPeakIIILat, maxPeakIVLat, maxPeakVLat],
                   minTroughLat=[minTroughILat, minTroughIILat, minTroughIIILat, minTroughIVLat, minTroughVLat],
                   maxTroughLat=[maxTroughILat, maxTroughIILat, maxTroughIIILat, maxTroughIVLat, maxTroughVLat])

    return df
end







