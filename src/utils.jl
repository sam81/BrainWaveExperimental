using DocStringExtensions



"""
Compute ABR peak points on the basis of peak latencies in another condition.
The delay is recomputed on the basis of peak V delay if peak V is found.
"""
function calcDelayedPeakPoints(sig, pkTimes, delay, sampRate; epochStart=0)
    nPeaks = length(pkTimes)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart)
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart)

    ## try to find peak V first
    minPeakVLat = pkTimes[5]+delay[1]
    maxPeakVLat = pkTimes[5]+delay[2]

    if ismissing(minPeakVLat) == false
        peakVTimes = xpksTimes[(xpksTimes .>= minPeakVLat) .& (xpksTimes .<= maxPeakVLat)]
        peakVPnts = xpksPnts[(xpksTimes .>= minPeakVLat) .& (xpksTimes .<= maxPeakVLat)]
        if length(peakVTimes) > 0
            ## find largest peak
            idx = findall(sig[peakVPnts] .== maximum(sig[peakVPnts]))[1]
            peakVTime = peakVTimes[idx]
            peakVDelay = peakVTime - pkTimes[5]
            ## recompute delay on the basis of peak V delay
            delay = (peakVDelay-0.00035, peakVDelay+0.00035)
        end
    end
    peakTimes = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks)
    for pk=1:nPeaks
        if ismissing(pkTimes[pk]) == false
            fooTimes = xpksTimes[(xpksTimes .>= pkTimes[pk]+delay[1]) .& (xpksTimes .<= pkTimes[pk]+delay[2])]
            fooPnts = xpksPnts[(xpksTimes.>=pkTimes[pk]+delay[1]) .& (xpksTimes.<=pkTimes[pk]+delay[2])]
        else
            fooTimes = (Float64)[]
            fooPnts = (Int)[]
        end

        if length(fooTimes) > 0
            #find largest peak
            idx = findall(sig[fooPnts] .== maximum(sig[fooPnts]))[1]
            peakTimes[pk] = fooTimes[idx]
            peakPoints[pk] = fooPnts[idx]
        else
            peakTimes[pk] = missing
            peakPoints[pk] = missing
        end
    end

    minPeakTroughLat = [0.25/1000, 0.12/1000, 0.25/1000, 0.12/1000, 0.25/1000]
    maxPeakTroughLatDiff = [0.75/1000, missing, 0.75/1000, missing, 1.5/1000]
    maxPeakTroughLat = minPeakTroughLat + maxPeakTroughLatDiff

    troughTimes = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks)
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        if ismissing(maxPeakTroughLat[pk]) == true
            if pk+1 <= 5
                thisMax = peakTimes[pk+1] - peakTimes[pk]
            end
        else
            thisMax = maxPeakTroughLat[pk]
        end
        minTroughLat[pk] = peakTimes[pk]+thisMin
        maxTroughLat[pk] = peakTimes[pk]+thisMax

        if (ismissing(peakTimes[pk]+thisMin+thisMax) == false)
            fooTimes = xtrsTimes[(xtrsTimes .>= peakTimes[pk]+thisMin) .& (xtrsTimes .<= peakTimes[pk]+thisMax)]
            fooPnts = xtrsPnts[(xtrsTimes.>= peakTimes[pk]+thisMin) .& (xtrsTimes.<= peakTimes[pk]+thisMax)]
        else
            fooTimes = (Float64)[]
            fooPnts = (Int)[]
        end

        if length(fooTimes) > 0
            #find largest trough
            idx = findall(sig[fooPnts] .== minimum(sig[fooPnts]))[1]
            troughTimes[pk] = fooTimes[idx]
            troughPoints[pk] = fooPnts[idx]
        else
            troughTimes[pk] = missing
            troughPoints[pk] = missing
        end
    end

    peakLatencies = peakTimes
    troughLatencies = troughTimes

    peakAmps = missings(Float64, length(peakPoints))
    troughAmps = missings(Float64, length(peakPoints))
    for i=1:length(peakAmps)
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[1, peakPoints[i]]
        else
            peakAmps[i] = missing
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[1, troughPoints[i]]
        else
            troughAmps[i] = missing
        end
    end

    peakTroughAmps = peakAmps-troughAmps

    waveLabels = ["I", "II", "III", "IV", "V"]
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps,
                   minPeakLat=vcat(pkTimes[1:4].+delay[1], minPeakVLat),
                   maxPeakLat=vcat(pkTimes[1:4].+delay[2], maxPeakVLat),
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end

"""
Given a vector of expected peak times (e.g. from a grand average),
find peak points within a time window of exp_pk_time +/- delay.

The function first finds peaks

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the peaks and troughs are sought.
* `expPeakLatencies::AbstractVector{Real}`: vector of expected latencies of the peaks sought, in seconds.
* `winStart::AbstractVector{Real}`: the start of the time windows relative to the expected peak latencies in which to find peaks, in seconds.
* `winStop::AbstractVector{Real}`: the end of the time windows relative to the expected peak latencies in which to find peaks, in seconds.
* `sampRate::Real`: sampling rate of the waveform in Hz.
* `epochStart::Real`: the time, in seconds, at which the epoch starts relative to the starting time of the stimulus

##### Returns

* ``:

##### Examples

```julia

```

"""
function calcPeaksFromExpectedLatencies(sig::Union{AbstractMatrix{T}, AbstractVector{T}},
                                        expPeakLatencies::Union{AbstractVector{S}, AbstractVector{Union{S, Missing}}},
                                        winStart::AbstractVector{P}, winStop::AbstractVector{Q}, sampRate::Real;
                                        epochStart::Real=0,
                                        minPeakTroughLat::AbstractVector{G}=[0.25, 0.12, 0.25, 0.12, 0.25]./1000,
                                        maxPeakTroughLat::AbstractVector{H}=[1, 1, 1, 0.75, 2]./1000,
                                        minAbsLat::Union{AbstractVector{S}, AbstractVector{Union{R, Missing}}}=[missing for i=1:5],
                                        maxAbsLat::Union{AbstractVector{S}, AbstractVector{Union{W, Missing}}}=[missing for i=1:5],
                                        waveLabels::AbstractVector{String} = ["I", "II", "III", "IV", "V"]) where {T<:Real, S<:Real, P<:Real, Q<:Real, G<:Real, H<:Real, R<:Real, W<:Real}
    nPeaks = length(expPeakLatencies)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart) #find all peaks in the waveform
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart) #find all troughs in the waveform

    #############################
    ## find peaks
    peakLatencies = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks)
    peakLatenciesEst = missings(Float64, nPeaks)
    peakPointsEst = missings(Int, nPeaks)

    searchWinStart = missings(Float64, nPeaks)
    searchWinStop = missings(Float64, nPeaks)
    for pk=1:nPeaks
        searchWinStart[pk] = max(expPeakLatencies[pk]+winStart[pk], ifelse(ismissing(minAbsLat[pk]), expPeakLatencies[pk]+winStart[pk], minAbsLat[pk]))
        searchWinStop[pk] = min(expPeakLatencies[pk]+winStop[pk], ifelse(ismissing(maxAbsLat[pk]), expPeakLatencies[pk]+winStop[pk], maxAbsLat[pk]))
    end
    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            candidateTimes = xpksTimes[(xpksTimes .>= searchWinStart[pk]) .& (xpksTimes .<= searchWinStop[pk])]
            candidatePnts = xpksPnts[(xpksTimes.>=searchWinStart[pk]) .& (xpksTimes.<=searchWinStop[pk])]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest peak
            idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
            peakLatencies[pk] = candidateTimes[idx]
            peakPoints[pk] = candidatePnts[idx]
            peakLatenciesEst[pk] = candidateTimes[idx]
            peakPointsEst[pk] = candidatePnts[idx]
        else
            peakLatencies[pk] = missing
            peakPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (searchWinStart[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (searchWinStop[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== maximum(sig[idxStart:idxStop])) + idxStart-1
                peakPointsEst[pk] = idx
                peakLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                peakLatenciesEst[pk] = missing
                peakPointsEst[pk] = missing
            end
        end
    end

    ###############################
    ## determine trough windows from peak latencies
    troughLatencies = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    troughLatenciesEst = missings(Float64, nPeaks)
    troughPointsEst = missings(Int, nPeaks) #has to be float to possibly be missing
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    minTroughLatEst = missings(Float64, nPeaks)
    maxTroughLatEst = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        thisMax = maxPeakTroughLat[pk]
        ## if ismissing(maxPeakTroughLat[pk]) == true
        ##     if pk+1 <= nPeaks
        ##         thisMax = peakLatencies[pk+1] - peakLatencies[pk]
        ##     end
        ## else
        ##     thisMax = maxPeakTroughLat[pk]
        ## end
        minTroughLat[pk] = peakLatencies[pk]+thisMin
        maxTroughLat[pk] = peakLatencies[pk]+thisMax
        minTroughLatEst[pk]= peakLatenciesEst[pk]+thisMin
        maxTroughLatEst[pk] =  peakLatenciesEst[pk]+thisMax

        if ismissing(peakLatencies[pk]) == false
            candidateTimes = xtrsTimes[(xtrsTimes .>= peakLatencies[pk]+thisMin) .& (xtrsTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xtrsPnts[(xtrsTimes.>= peakLatencies[pk]+thisMin) .& (xtrsTimes.<= peakLatencies[pk]+thisMax)]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end

        if length(candidateTimes) > 0
            #find largest trough
            idx = findall(sig[candidatePnts] .== minimum(sig[candidatePnts]))[1]
            troughLatencies[pk] = candidateTimes[idx]
            troughPoints[pk] = candidatePnts[idx]
            troughLatenciesEst[pk] = candidateTimes[idx]
            troughPointsEst[pk] = candidatePnts[idx]
        else
            troughLatencies[pk] = missing
            troughPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (minTroughLatEst[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (maxTroughLatEst[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== minimum(sig[idxStart:idxStop])) + idxStart-1
                troughPointsEst[pk] = idx
                troughLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                troughLatenciesEst[pk] = missing
                troughPointsEst[pk] = missing
            end
        end
    end

    #####################################
    ## measure peak and trough amplitudes
    peakAmps = missings(Float64, nPeaks)#zeros(nPeaks)
    troughAmps = missings(Float64, nPeaks)#zeros(nPeaks)
    peakAmpsEst = missings(Float64, nPeaks)#zeros(nPeaks)
    troughAmpsEst = missings(Float64, nPeaks)#zeros(nPeaks)
    for i=1:nPeaks
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            peakAmpsEst[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            if ismissing(expPeakLatencies[i]) == false
                #idxStart = round(Int, (expPeakLatencies[i]+winStart[i]-epochStart)*sampRate);
                #idxStop = round(Int, (expPeakLatencies[i]+winStop[i]-epochStart)*sampRate);
                peakAmpsEst[i] = sig[round(Int, peakPointsEst[i])] #maximum(sig[idxStart:idxStop])[1]
            else
                peakAmpsEst[i] = missing
            end
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
            troughAmpsEst[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
            if ismissing(minTroughLatEst[i]) == false
                #idxStart = round(Int, (minTroughLat[i]-epochStart)*sampRate);
                #idxStop = round(Int, (maxTroughLat[i]-epochStart)*sampRate);
                troughAmpsEst[i] = sig[round(Int, troughPointsEst[i])] #minimum(sig[idxStart:idxStop])[1]
            else
                troughAmpsEst[i] = missing
            end
        end
    end

    ####################################
    ## measure peak-trough amplitudes
    peakTroughAmps = peakAmps-troughAmps
    peakTroughAmpsEst = peakAmpsEst-troughAmpsEst

    ## #############################################################
    ## Estimate peaks and troughs by finding point that maximizes
    ## peak-trough amplitude
    ## #############################################################
    peakPointsEstPTMax = missings(Int, nPeaks)
    troughPointsEstPTMax = missings(Int, nPeaks)
    peakAmpsEstPTMax = missings(Float64, nPeaks)
    troughAmpsEstPTMax = missings(Float64, nPeaks)
    peakLatenciesEstPTMax = missings(Float64, nPeaks)
    troughLatenciesEstPTMax = missings(Float64, nPeaks)

    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            minPeakIdx = round(Int, (searchWinStart[pk]-epochStart)*sampRate)
            maxPeakIdx = round(Int, (searchWinStop[pk]-epochStart)*sampRate)

            bestDiff = 0
            pkIdx = 1
            trIdx = 1
            for i=minPeakIdx:maxPeakIdx-1
                for j=i+round(Int, minPeakTroughLat[pk]*sampRate):i+round(Int, maxPeakTroughLat[pk]*sampRate)
                    currDiff = sig[i]-sig[j]
                    if currDiff > bestDiff
                        bestDiff = copy(currDiff)
                        pkIdx = i
                        trIdx = j
                    end
                end
            end

            peakPointsEstPTMax[pk] = pkIdx
            troughPointsEstPTMax[pk] = trIdx
            peakAmpsEstPTMax[pk] = sig[pkIdx]
            troughAmpsEstPTMax[pk] = sig[trIdx]
            peakLatenciesEstPTMax[pk] = (pkIdx-1)/sampRate+epochStart
            troughLatenciesEstPTMax[pk] = (trIdx-1)/sampRate+epochStart
        end
    end
    peakTroughAmpsEstPTMax = peakAmpsEstPTMax-troughAmpsEstPTMax


    ####################################
    ## store results in a dataframe
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps, peakPointEst=peakPointsEst, peakLatencyEst=peakLatenciesEst, peakAmpEst=peakAmpsEst, troughPointEst=troughPointsEst, troughLatencyEst=troughLatenciesEst, troughAmpEst=troughAmpsEst, peakTroughAmpEst=peakTroughAmpsEst,
                   peakPointEstPTMax=peakPointsEstPTMax, peakLatencyEstPTMax=peakLatenciesEstPTMax, peakAmpEstPTMax=peakAmpsEstPTMax, troughPointEstPTMax=troughPointsEstPTMax, troughLatencyEstPTMax=troughLatenciesEstPTMax, troughAmpEstPTMax=troughAmpsEstPTMax, peakTroughAmpEstPTMax=peakTroughAmpsEstPTMax,
                   minPeakLat=expPeakLatencies+winStart,
                   maxPeakLat=expPeakLatencies+winStop,
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end

"""
Given a vector of expected peak times (e.g. from a grand average),
find peak points within a time window of exp_pk_time +/- delay.

The difference from `calcPeaksFromExpectedLatencies` is that if troughs are not found inflection points are used to mark troughs.
The function first finds peaks.

$(SIGNATURES)

##### Arguments

* `sig::Union{AbstractMatrix{T}, AbstractVector{T}}`: the waveform for which the peaks and troughs are sought.
* `expPeakLatencies::AbstractVector{Real}`: vector of expected latencies of the peaks sought, in seconds.
* `winStart::AbstractVector{Real}`: the start of the time windows relative to the expected peak latencies in which to find peaks, in seconds.
* `winStop::AbstractVector{Real}`: the end of the time windows relative to the expected peak latencies in which to find peaks, in seconds.
* `sampRate::Real`: sampling rate of the waveform in Hz.
* `epochStart::Real`: the time, in seconds, at which the epoch starts relative to the starting time of the stimulus

##### Returns

* ``:

##### Examples

```julia

```

"""
function calcPeaksTroughsInflFromExpectedLatencies(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, expPeakLatencies::Union{AbstractVector{S}, AbstractVector{Union{S, Missing}}}, winStart::AbstractVector{P}, winStop::AbstractVector{Q}, sampRate::Real; epochStart::Real=0, minPeakTroughLat::AbstractVector{G}=[0.25, 0.12, 0.25, 0.12, 0.25]./1000, maxPeakTroughLat::AbstractVector{H}=[1, 1, 1, 0.75, 2]./1000, waveLabels::AbstractVector{String} = ["I", "II", "III", "IV", "V"]) where {T<:Real, S<:Real, P<:Real, Q<:Real, G<:Real, H<:Real}
    nPeaks = length(expPeakLatencies)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart) #find all peaks in the waveform
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart) #find all troughs in the waveform
    xinflPnts, xinflTimes = findInflections(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform
    dy = diff(vec(sig))
    #############################
    ## find peaks
    peakLatencies = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks)
    peakLatenciesEst = missings(Float64, nPeaks)
    peakPointsEst = missings(Int, nPeaks)
    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            candidateTimes = xpksTimes[(xpksTimes .>= expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes .<= expPeakLatencies[pk]+winStop[pk])]
            candidatePnts = xpksPnts[(xpksTimes.>=expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes.<=expPeakLatencies[pk]+winStop[pk])]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest peak
            idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
            peakLatencies[pk] = candidateTimes[idx]
            peakPoints[pk] = candidatePnts[idx]
            peakLatenciesEst[pk] = candidateTimes[idx]
            peakPointsEst[pk] = candidatePnts[idx]
        else
            peakLatencies[pk] = missing
            peakPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (expPeakLatencies[pk]+winStart[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (expPeakLatencies[pk]+winStop[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== maximum(sig[idxStart:idxStop])) + idxStart-1
                peakPointsEst[pk] = idx
                peakLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                peakLatenciesEst[pk] = missing
                peakPointsEst[pk] = missing
            end
        end
    end

    ###############################
    ## determine trough windows from peak latencies
    troughLatencies = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    troughLatenciesEst = missings(Float64, nPeaks)
    troughPointsEst = missings(Int, nPeaks) #has to be float to possibly be missing
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    minTroughLatEst = missings(Float64, nPeaks)
    maxTroughLatEst = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        thisMax = maxPeakTroughLat[pk]
        ## if ismissing(maxPeakTroughLat[pk]) == true
        ##     if pk+1 <= nPeaks
        ##         thisMax = peakLatencies[pk+1] - peakLatencies[pk]
        ##     end
        ## else
        ##     thisMax = maxPeakTroughLat[pk]
        ## end
        minTroughLat[pk] = peakLatencies[pk]+thisMin
        maxTroughLat[pk] = peakLatencies[pk]+thisMax
        minTroughLatEst[pk]= peakLatenciesEst[pk]+thisMin
        maxTroughLatEst[pk] =  peakLatenciesEst[pk]+thisMax

        candidateTimes = (Float64)[]
        candidatePnts = (Int)[]
        if ismissing(peakLatencies[pk]) == false
            candidateTimes = xtrsTimes[(xtrsTimes .>= peakLatencies[pk]+thisMin) .& (xtrsTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xtrsPnts[(xtrsTimes.>= peakLatencies[pk]+thisMin) .& (xtrsTimes.<= peakLatencies[pk]+thisMax)]
            if length(candidateTimes) > 0
                #find largest trough
                idx = findall(sig[candidatePnts] .== minimum(sig[candidatePnts]))[1]
                troughLatencies[pk] = candidateTimes[idx]
                troughPoints[pk] = candidatePnts[idx]
                troughLatenciesEst[pk] = candidateTimes[idx]
                troughPointsEst[pk] = candidatePnts[idx]
            end
        end

        if (ismissing(troughPoints[pk]) == true) & (ismissing(peakLatencies[pk]) == false) #try inflection points
            candidateTimes = xinflTimes[(xinflTimes .>= peakLatencies[pk]+thisMin) .& (xinflTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xinflPnts[(xinflTimes.>=peakLatencies[pk]+thisMin) .& (xinflTimes.<=peakLatencies[pk]+thisMax)]
            if length(candidateTimes) > 0
                idx = findall(abs.(dy[candidatePnts .- 1]) .== minimum(abs.(dy[candidatePnts .- 1])))[1] #choose point where first derivative is closest to zero
                troughLatencies[pk] = candidateTimes[idx]
                troughPoints[pk] = candidatePnts[idx]
                troughLatenciesEst[pk] = candidateTimes[idx]
                troughPointsEst[pk] = candidatePnts[idx]
            end
        end

        if ismissing(troughPointsEst[pk]) == true
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (minTroughLatEst[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (maxTroughLatEst[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== minimum(sig[idxStart:idxStop])) + idxStart-1
                troughPointsEst[pk] = idx
                troughLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                troughLatenciesEst[pk] = missing
                troughPointsEst[pk] = missing
            end
        end
    end



    #####################################
    ## measure peak and trough amplitudes
    peakAmps = missings(Float64, nPeaks)#zeros(nPeaks)
    troughAmps = missings(Float64, nPeaks)#zeros(nPeaks)
    peakAmpsEst = missings(Float64, nPeaks)#zeros(nPeaks)
    troughAmpsEst = missings(Float64, nPeaks)#zeros(nPeaks)
    for i=1:nPeaks
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            peakAmpsEst[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            if ismissing(expPeakLatencies[i]) == false
                #idxStart = round(Int, (expPeakLatencies[i]+winStart[i]-epochStart)*sampRate);
                #idxStop = round(Int, (expPeakLatencies[i]+winStop[i]-epochStart)*sampRate);
                peakAmpsEst[i] = sig[round(Int, peakPointsEst[i])] #maximum(sig[idxStart:idxStop])[1]
            else
                peakAmpsEst[i] = missing
            end
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
            troughAmpsEst[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
            if ismissing(minTroughLatEst[i]) == false
                #idxStart = round(Int, (minTroughLat[i]-epochStart)*sampRate);
                #idxStop = round(Int, (maxTroughLat[i]-epochStart)*sampRate);
                troughAmpsEst[i] = sig[round(Int, troughPointsEst[i])] #minimum(sig[idxStart:idxStop])[1]
            else
                troughAmpsEst[i] = missing
            end
        end
    end

    ####################################
    ## measure peak-trough amplitudes
    peakTroughAmps = peakAmps-troughAmps
    peakTroughAmpsEst = peakAmpsEst-troughAmpsEst

    ## #############################################################
    ## Estimate peaks and troughs by finding point that maximizes
    ## peak-trough amplitude
    ## #############################################################
    peakPointsEstPTMax = missings(Int, nPeaks)
    troughPointsEstPTMax = missings(Int, nPeaks)
    peakAmpsEstPTMax = missings(Float64, nPeaks)
    troughAmpsEstPTMax = missings(Float64, nPeaks)
    peakLatenciesEstPTMax = missings(Float64, nPeaks)
    troughLatenciesEstPTMax = missings(Float64, nPeaks)

    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            minPeakIdx = round(Int, (expPeakLatencies[pk]+winStart[pk]-epochStart)*sampRate)
            maxPeakIdx = round(Int, (expPeakLatencies[pk]+winStop[pk]-epochStart)*sampRate)

            bestDiff = 0
            pkIdx = 1
            trIdx = 1
            for i=minPeakIdx:maxPeakIdx-1
                for j=i+round(Int, minPeakTroughLat[pk]*sampRate):i+round(Int, maxPeakTroughLat[pk]*sampRate)
                    currDiff = sig[i]-sig[j]
                    if currDiff > bestDiff
                        bestDiff = copy(currDiff)
                        pkIdx = i
                        trIdx = j
                    end
                end
            end

            peakPointsEstPTMax[pk] = pkIdx
            troughPointsEstPTMax[pk] = trIdx
            peakAmpsEstPTMax[pk] = sig[pkIdx]
            troughAmpsEstPTMax[pk] = sig[trIdx]
            peakLatenciesEstPTMax[pk] = (pkIdx-1)/sampRate+epochStart
            troughLatenciesEstPTMax[pk] = (trIdx-1)/sampRate+epochStart
        end
    end
    peakTroughAmpsEstPTMax = peakAmpsEstPTMax-troughAmpsEstPTMax


    ####################################
    ## store results in a dataframe
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps, peakPointEst=peakPointsEst, peakLatencyEst=peakLatenciesEst, peakAmpEst=peakAmpsEst, troughPointEst=troughPointsEst, troughLatencyEst=troughLatenciesEst, troughAmpEst=troughAmpsEst, peakTroughAmpEst=peakTroughAmpsEst,
                   peakPointEstPTMax=peakPointsEstPTMax, peakLatencyEstPTMax=peakLatenciesEstPTMax, peakAmpEstPTMax=peakAmpsEstPTMax, troughPointEstPTMax=troughPointsEstPTMax, troughLatencyEstPTMax=troughLatenciesEstPTMax, troughAmpEstPTMax=troughAmpsEstPTMax, peakTroughAmpEstPTMax=peakTroughAmpsEstPTMax,
                   minPeakLat=expPeakLatencies+winStart,
                   maxPeakLat=expPeakLatencies+winStop,
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end



function calcPeaksInflFromExpectedLatencies(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, expPeakLatencies::Union{AbstractVector{S}, AbstractVector{Union{S, Missing}}}, winStart::AbstractVector{P}, winStop::AbstractVector{Q}, sampRate::Real; epochStart::Real=0, minPeakTroughLat::AbstractVector{G}=[0.25, 0.12, 0.25, 0.12, 0.25]./1000, maxPeakTroughLat::AbstractVector{H}=[1, 1, 1, 0.75, 2]./1000, waveLabels::AbstractVector{String} = ["I", "II", "III", "IV", "V"]) where {T<:Real, S<:Real, P<:Real, Q<:Real, G<:Real, H<:Real}
    nPeaks = length(expPeakLatencies)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart) #find all peaks in the waveform
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart) #find all troughs in the waveform
    xinflPnts, xinflTimes = findInflections(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform

    #############################
    ## find peaks
    peakLatencies = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    peakLatenciesEst = missings(Float64, nPeaks) #in case a peak can't be found fall back to finding the maximum value
    peakPointsEst = missings(Int, nPeaks)
    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            candidateTimes = xpksTimes[(xpksTimes .>= expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes .<= expPeakLatencies[pk]+winStop[pk])]
            candidatePnts = xpksPnts[(xpksTimes.>=expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes.<=expPeakLatencies[pk]+winStop[pk])]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest peak
            idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
            peakLatencies[pk] = candidateTimes[idx]
            peakPoints[pk] = candidatePnts[idx]
            peakLatenciesEst[pk] = candidateTimes[idx]
            peakPointsEst[pk] = candidatePnts[idx]
        else
            ##try inflection points
            if ismissing(expPeakLatencies[pk]) == false
                candidateTimes = xinflTimes[(xinflTimes .>= expPeakLatencies[pk]+winStart[pk]) .& (xinflTimes .<= expPeakLatencies[pk]+winStop[pk])]
                candidatePnts = xinflPnts[(xinflTimes.>=expPeakLatencies[pk]+winStart[pk]) .& (xinflTimes.<=expPeakLatencies[pk]+winStop[pk])]
            else
                candidateTimes = (Float64)[]
                candidatePnts = (Int)[]
            end
            if length(candidateTimes) > 0
                #find largest peak
                idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
                peakLatencies[pk] = candidateTimes[idx]
                peakPoints[pk] = candidatePnts[idx]
                peakLatenciesEst[pk] = candidateTimes[idx]
                peakPointsEst[pk] = candidatePnts[idx]
            else
                peakLatencies[pk] = missing
                peakPoints[pk] = missing
                if ismissing(expPeakLatencies[pk]) == false
                    idxStart = round(Int, (expPeakLatencies[pk]+winStart[pk]-epochStart)*sampRate)+1;
                    idxStop = round(Int, (expPeakLatencies[pk]+winStop[pk]-epochStart)*sampRate)+1;
                    idx = findfirst(sig[idxStart:idxStop] .== maximum(sig[idxStart:idxStop])) + idxStart-1
                    peakPointsEst[pk] = idx
                    peakLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
                else
                    peakLatenciesEst[pk] = missing
                    peakPointsEst[pk] = missing
                end
            end
        end
    end

    ###############################
    ## determine trough windows from peak latencies
    troughLatencies = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    troughLatenciesEst = missings(Float64, nPeaks)
    troughPointsEst = missings(Int, nPeaks) #has to be float to possibly be missing
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    minTroughLatEst = missings(Float64, nPeaks)
    maxTroughLatEst = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        thisMax = maxPeakTroughLat[pk]
        ## if ismissing(maxPeakTroughLat[pk]) == true
        ##     if pk+1 <= nPeaks
        ##         thisMax = peakLatencies[pk+1] - peakLatencies[pk]
        ##     end
        ## else
        ##     thisMax = maxPeakTroughLat[pk]
        ## end
        minTroughLat[pk] = peakLatencies[pk]+thisMin
        maxTroughLat[pk] =  peakLatencies[pk]+thisMax
        minTroughLatEst[pk] = peakLatenciesEst[pk]+thisMin
        maxTroughLatEst[pk] =  peakLatenciesEst[pk]+thisMax
        if ismissing(peakLatencies[pk]) == false
            candidateTimes = xtrsTimes[(xtrsTimes .>= peakLatencies[pk]+thisMin) .& (xtrsTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xtrsPnts[(xtrsTimes.>= peakLatencies[pk]+thisMin) .& (xtrsTimes.<= peakLatencies[pk]+thisMax)]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest trough
            idx = findall(sig[candidatePnts] .== minimum(sig[candidatePnts]))[1]
            troughLatencies[pk] = candidateTimes[idx]
            troughPoints[pk] = candidatePnts[idx]
            troughLatenciesEst[pk] = candidateTimes[idx]
            troughPointsEst[pk] = candidatePnts[idx]
        else
            troughLatencies[pk] = missing
            troughPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (minTroughLatEst[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (maxTroughLatEst[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== minimum(sig[idxStart:idxStop])) + idxStart-1
                troughPointsEst[pk] = idx
                troughLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                troughLatenciesEst[pk] = missing
                troughPointsEst[pk] = missing
            end
        end
    end

    #####################################
    ## measure peak and trough amplitudes
    peakAmps = missings(Float64, nPeaks)
    troughAmps = missings(Float64, nPeaks)
    peakAmpsEst = missings(Float64, nPeaks)
    troughAmpsEst = missings(Float64, nPeaks)
    for i=1:nPeaks
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            peakAmpsEst[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            if ismissing(expPeakLatencies[i]) == false
                #idxStart = round(Int, (expPeakLatencies[i]+winStart[i]-epochStart)*sampRate);
                #idxStop = round(Int, (expPeakLatencies[i]+winStop[i]-epochStart)*sampRate);
                peakAmpsEst[i] = sig[round(Int, peakPointsEst[i])] #maximum(sig[idxStart:idxStop])[1]
            else
                peakAmpsEst[i] = missing
            end
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
            troughAmpsEst[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
            if ismissing(minTroughLatEst[i]) == false
                #idxStart = round(Int, (minTroughLat[i]-epochStart)*sampRate);
                #idxStop = round(Int, (maxTroughLat[i]-epochStart)*sampRate);
                troughAmpsEst[i] = sig[round(Int, troughPointsEst[i])] #minimum(sig[idxStart:idxStop])[1]
            else
                troughAmpsEst[i] = missing
            end
        end
    end

    ####################################
    ## measure peak-trough amplitudes
    peakTroughAmps = peakAmps-troughAmps
    peakTroughAmpsEst = peakAmpsEst-troughAmpsEst

    ####################################
    ## store results in a dataframe
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps, peakPointEst=peakPointsEst, peakLatencyEst=peakLatenciesEst, peakAmpEst=peakAmpsEst, troughPointEst=troughPointsEst, troughLatencyEst=troughLatenciesEst, troughAmpEst=troughAmpsEst, peakTroughAmpEst=peakTroughAmpsEst,
                   minPeakLat=expPeakLatencies+winStart,
                   maxPeakLat=expPeakLatencies+winStop,
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end

##difference from calcPeaksInflFromExpectedLatencies is that this function only uses
## down-going inflections to find peaks. Also, the inflection point with the lowest absolute first derivative is chosen.
function calcPeaksInflFromExpectedLatencies2(sig::Union{AbstractMatrix{T}, AbstractVector{T}}, expPeakLatencies::Union{AbstractVector{S}, AbstractVector{Union{S, Missing}}}, winStart::AbstractVector{P}, winStop::AbstractVector{Q}, sampRate::Real; epochStart::Real=0, minPeakTroughLat::AbstractVector{G}=[0.25, 0.12, 0.25, 0.12, 0.25]./1000, maxPeakTroughLat::AbstractVector{H}=[1, 1, 1, 0.75, 2]./1000, waveLabels::AbstractVector{String} = ["I", "II", "III", "IV", "V"]) where {T<:Real, S<:Real, P<:Real, Q<:Real, G<:Real, H<:Real}
    nPeaks = length(expPeakLatencies)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart) #find all peaks in the waveform
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart) #find all troughs in the waveform
    xinflPnts, xinflTimes = findInflections(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform
    xinflPntsUp, xinflTimesUp = findInflectionsUp(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform
    xinflPntsDown, xinflTimesDown = findInflectionsDown(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform

    #############################
    ## find peaks
    peakLatencies = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    peakLatenciesEst = missings(Float64, nPeaks) #in case a peak can't be found fall back to finding the maximum value
    peakPointsEst = missings(Int, nPeaks)
    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            candidateTimes = xpksTimes[(xpksTimes .>= expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes .<= expPeakLatencies[pk]+winStop[pk])]
            candidatePnts = xpksPnts[(xpksTimes.>=expPeakLatencies[pk]+winStart[pk]) .& (xpksTimes.<=expPeakLatencies[pk]+winStop[pk])]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest peak
            idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
            peakLatencies[pk] = candidateTimes[idx]
            peakPoints[pk] = candidatePnts[idx]
            peakLatenciesEst[pk] = candidateTimes[idx]
            peakPointsEst[pk] = candidatePnts[idx]
        else
            ##try inflection points
            if ismissing(expPeakLatencies[pk]) == false
                candidateTimes = xinflTimesDown[(xinflTimesDown .>= expPeakLatencies[pk]+winStart[pk]) .& (xinflTimesDown .<= expPeakLatencies[pk]+winStop[pk])]
                candidatePnts = xinflPntsDown[(xinflTimesDown.>=expPeakLatencies[pk]+winStart[pk]) .& (xinflTimesDown.<=expPeakLatencies[pk]+winStop[pk])]
            else
                candidateTimes = (Float64)[]
                candidatePnts = (Int)[]
            end
            if length(candidateTimes) > 0
                ###find largest peak
                ##idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
                ## pick point where first derivative has minimum absolute value
                dy = diff(sig)
                idx = findall(abs.(dy[candidatePnts.-1]) .== minimum.(abs.(dy[candidatePnts.-1])))[1]#sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
                peakLatencies[pk] = candidateTimes[idx]
                peakPoints[pk] = candidatePnts[idx]
                peakLatenciesEst[pk] = candidateTimes[idx]
                peakPointsEst[pk] = candidatePnts[idx]
            else
                peakLatencies[pk] = missing
                peakPoints[pk] = missing
                if ismissing(expPeakLatencies[pk]) == false
                    idxStart = round(Int, (expPeakLatencies[pk]+winStart[pk]-epochStart)*sampRate)+1;
                    idxStop = round(Int, (expPeakLatencies[pk]+winStop[pk]-epochStart)*sampRate)+1;
                    idx = findfirst(sig[idxStart:idxStop] .== maximum(sig[idxStart:idxStop])) + idxStart-1
                    peakPointsEst[pk] = idx
                    peakLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
                else
                    peakLatenciesEst[pk] = missing
                    peakPointsEst[pk] = missing
                end
            end
        end
    end

    ###############################
    ## determine trough windows from peak latencies
    troughLatencies = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    troughLatenciesEst = missings(Float64, nPeaks)
    troughPointsEst = missings(Int, nPeaks) #has to be float to possibly be missing
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    minTroughLatEst = missings(Float64, nPeaks)
    maxTroughLatEst = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        thisMax = maxPeakTroughLat[pk]
        ## if ismissing(maxPeakTroughLat[pk]) == true
        ##     if pk+1 <= nPeaks
        ##         thisMax = peakLatencies[pk+1] - peakLatencies[pk]
        ##     end
        ## else
        ##     thisMax = maxPeakTroughLat[pk]
        ## end
        minTroughLat[pk] = peakLatencies[pk]+thisMin
        maxTroughLat[pk] =  peakLatencies[pk]+thisMax
        minTroughLatEst[pk] = peakLatenciesEst[pk]+thisMin
        maxTroughLatEst[pk] =  peakLatenciesEst[pk]+thisMax
        if ismissing(peakLatencies[pk]) == false
            candidateTimes = xtrsTimes[(xtrsTimes .>= peakLatencies[pk]+thisMin) .& (xtrsTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xtrsPnts[(xtrsTimes.>= peakLatencies[pk]+thisMin) .& (xtrsTimes.<= peakLatencies[pk]+thisMax)]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest trough
            idx = findall(sig[candidatePnts] .== minimum(sig[candidatePnts]))[1]
            troughLatencies[pk] = candidateTimes[idx]
            troughPoints[pk] = candidatePnts[idx]
            troughLatenciesEst[pk] = candidateTimes[idx]
            troughPointsEst[pk] = candidatePnts[idx]
        else
            troughLatencies[pk] = missing
            troughPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (minTroughLatEst[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (maxTroughLatEst[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== minimum(sig[idxStart:idxStop])) + idxStart-1
                troughPointsEst[pk] = idx
                troughLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                troughLatenciesEst[pk] = missing
                troughPointsEst[pk] = missing
            end
        end
    end

    #####################################
    ## measure peak and trough amplitudes
    peakAmps = missings(Float64, nPeaks)
    troughAmps = missings(Float64, nPeaks)
    peakAmpsEst = missings(Float64, nPeaks)
    troughAmpsEst = missings(Float64, nPeaks)
    for i=1:nPeaks
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            peakAmpsEst[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            if ismissing(expPeakLatencies[i]) == false
                #idxStart = round(Int, (expPeakLatencies[i]+winStart[i]-epochStart)*sampRate);
                #idxStop = round(Int, (expPeakLatencies[i]+winStop[i]-epochStart)*sampRate);
                peakAmpsEst[i] = sig[round(Int, peakPointsEst[i])] #maximum(sig[idxStart:idxStop])[1]
            else
                peakAmpsEst[i] = missing
            end
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
            troughAmpsEst[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
            if ismissing(minTroughLatEst[i]) == false
                #idxStart = round(Int, (minTroughLat[i]-epochStart)*sampRate);
                #idxStop = round(Int, (maxTroughLat[i]-epochStart)*sampRate);
                troughAmpsEst[i] = sig[round(Int, troughPointsEst[i])] #minimum(sig[idxStart:idxStop])[1]
            else
                troughAmpsEst[i] = missing
            end
        end
    end

    ####################################
    ## measure peak-trough amplitudes
    peakTroughAmps = peakAmps-troughAmps
    peakTroughAmpsEst = peakAmpsEst-troughAmpsEst

    ####################################
    ## store results in a dataframe
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps, peakPointEst=peakPointsEst, peakLatencyEst=peakLatenciesEst, peakAmpEst=peakAmpsEst, troughPointEst=troughPointsEst, troughLatencyEst=troughLatenciesEst, troughAmpEst=troughAmpsEst, peakTroughAmpEst=peakTroughAmpsEst,
                   minPeakLat=expPeakLatencies+winStart,
                   maxPeakLat=expPeakLatencies+winStop,
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end

#difference from calcPeaksInflFromExpectedLatencies2 is that there are further constraint on search window (useful for example
# to make sure wave at a lower level has delayed latency compared to wave previously measure at a higher level
function calcPeaksInflFromExpectedLatencies3(sig::Union{AbstractMatrix{T}, AbstractVector{T}},
                                             expPeakLatencies::Union{AbstractVector{S}, AbstractVector{Union{S, Missing}}},
                                             winStart::AbstractVector{P}, winStop::AbstractVector{Q},
                                             sampRate::Real; epochStart::Real=0,
                                             minPeakTroughLat::AbstractVector{G}=[0.25, 0.12, 0.25, 0.12, 0.25]./1000,
                                             maxPeakTroughLat::AbstractVector{H}=[1, 1, 1, 0.75, 2]./1000,
                                             minAbsLat::Union{AbstractVector{S}, AbstractVector{Union{R, Missing}}}=[missing for i=1:5],
                                             maxAbsLat::Union{AbstractVector{S}, AbstractVector{Union{W, Missing}}}=[missing for i=1:5],
                                             waveLabels::AbstractVector{String} = ["I", "II", "III", "IV", "V"]) where {T<:Real, S<:Real, P<:Real, Q<:Real, G<:Real, H<:Real, R<:Real, W<:Real}
    nPeaks = length(expPeakLatencies)
    xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart) #find all peaks in the waveform
    xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart) #find all troughs in the waveform
    xinflPnts, xinflTimes = findInflections(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform
    xinflPntsUp, xinflTimesUp = findInflectionsUp(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform
    xinflPntsDown, xinflTimesDown = findInflectionsDown(sig, sampRate, epochStart=epochStart) #find all inflections in the waveform

    #############################
    ## find peaks
    peakLatencies = missings(Float64, nPeaks)
    peakPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    peakLatenciesEst = missings(Float64, nPeaks) #in case a peak can't be found fall back to finding the maximum value
    peakPointsEst = missings(Int, nPeaks)

    searchWinStart = missings(Float64, nPeaks)
    searchWinStop = missings(Float64, nPeaks)
    for pk=1:nPeaks
        searchWinStart[pk] = max(expPeakLatencies[pk]+winStart[pk], ifelse(ismissing(minAbsLat[pk]), expPeakLatencies[pk]+winStart[pk], minAbsLat[pk]))
        searchWinStop[pk] = min(expPeakLatencies[pk]+winStop[pk], ifelse(ismissing(maxAbsLat[pk]), expPeakLatencies[pk]+winStop[pk], maxAbsLat[pk]))
    end
    for pk=1:nPeaks
        if ismissing(expPeakLatencies[pk]) == false
            candidateTimes = xpksTimes[(xpksTimes .>= searchWinStart[pk]) .& (xpksTimes .<= searchWinStop[pk])]
            candidatePnts = xpksPnts[(xpksTimes.>=searchWinStart[pk]) .& (xpksTimes.<=searchWinStop[pk])]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest peak
            idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
            peakLatencies[pk] = candidateTimes[idx]
            peakPoints[pk] = candidatePnts[idx]
            peakLatenciesEst[pk] = candidateTimes[idx]
            peakPointsEst[pk] = candidatePnts[idx]
        else
            ##try inflection points
            if ismissing(expPeakLatencies[pk]) == false
                candidateTimes = xinflTimesDown[(xinflTimesDown .>= searchWinStart[pk]) .& (xinflTimesDown .<= searchWinStop[pk])]
                candidatePnts = xinflPntsDown[(xinflTimesDown.>=searchWinStart[pk]) .& (xinflTimesDown.<=searchWinStop[pk])]
            else
                candidateTimes = (Float64)[]
                candidatePnts = (Int)[]
            end
            if length(candidateTimes) > 0
                ###find largest peak
                ##idx = findall(sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
                ## pick point where first derivative has minimum absolute value
                dy = diff(sig)
                idx = findall(abs.(dy[candidatePnts.-1]) .== minimum.(abs.(dy[candidatePnts.-1])))[1]#sig[candidatePnts] .== maximum(sig[candidatePnts]))[1]
                peakLatencies[pk] = candidateTimes[idx]
                peakPoints[pk] = candidatePnts[idx]
                peakLatenciesEst[pk] = candidateTimes[idx]
                peakPointsEst[pk] = candidatePnts[idx]
            else
                peakLatencies[pk] = missing
                peakPoints[pk] = missing
                if ismissing(expPeakLatencies[pk]) == false
                    idxStart = round(Int, (searchWinStart[pk]-epochStart)*sampRate)+1;
                    idxStop = round(Int, (searchWinStop[pk]-epochStart)*sampRate)+1;
                    idx = findfirst(sig[idxStart:idxStop] .== maximum(sig[idxStart:idxStop])) + idxStart-1
                    peakPointsEst[pk] = idx
                    peakLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
                else
                    peakLatenciesEst[pk] = missing
                    peakPointsEst[pk] = missing
                end
            end
        end
    end

    ###############################
    ## determine trough windows from peak latencies
    troughLatencies = missings(Float64, nPeaks)
    troughPoints = missings(Int, nPeaks) #has to be float to possibly be missing
    troughLatenciesEst = missings(Float64, nPeaks)
    troughPointsEst = missings(Int, nPeaks) #has to be float to possibly be missing
    minTroughLat = missings(Float64, nPeaks)
    maxTroughLat = missings(Float64, nPeaks)
    minTroughLatEst = missings(Float64, nPeaks)
    maxTroughLatEst = missings(Float64, nPeaks)
    for pk=1:nPeaks
        thisMin = minPeakTroughLat[pk]
        thisMax = maxPeakTroughLat[pk]
        ## if ismissing(maxPeakTroughLat[pk]) == true
        ##     if pk+1 <= nPeaks
        ##         thisMax = peakLatencies[pk+1] - peakLatencies[pk]
        ##     end
        ## else
        ##     thisMax = maxPeakTroughLat[pk]
        ## end
        minTroughLat[pk] = peakLatencies[pk]+thisMin
        maxTroughLat[pk] =  peakLatencies[pk]+thisMax
        minTroughLatEst[pk] = peakLatenciesEst[pk]+thisMin
        maxTroughLatEst[pk] =  peakLatenciesEst[pk]+thisMax
        if ismissing(peakLatencies[pk]) == false
            candidateTimes = xtrsTimes[(xtrsTimes .>= peakLatencies[pk]+thisMin) .& (xtrsTimes .<= peakLatencies[pk]+thisMax)]
            candidatePnts = xtrsPnts[(xtrsTimes.>= peakLatencies[pk]+thisMin) .& (xtrsTimes.<= peakLatencies[pk]+thisMax)]
        else
            candidateTimes = (Float64)[]
            candidatePnts = (Int)[]
        end
        if length(candidateTimes) > 0
            #find largest trough
            idx = findall(sig[candidatePnts] .== minimum(sig[candidatePnts]))[1]
            troughLatencies[pk] = candidateTimes[idx]
            troughPoints[pk] = candidatePnts[idx]
            troughLatenciesEst[pk] = candidateTimes[idx]
            troughPointsEst[pk] = candidatePnts[idx]
        else
            troughLatencies[pk] = missing
            troughPoints[pk] = missing
            if ismissing(expPeakLatencies[pk]) == false
                idxStart = round(Int, (minTroughLatEst[pk]-epochStart)*sampRate)+1;
                idxStop = round(Int, (maxTroughLatEst[pk]-epochStart)*sampRate)+1;
                idx = findfirst(sig[idxStart:idxStop] .== minimum(sig[idxStart:idxStop])) + idxStart-1
                troughPointsEst[pk] = idx
                troughLatenciesEst[pk] = ((idx-1)/sampRate)+epochStart
            else
                troughLatenciesEst[pk] = missing
                troughPointsEst[pk] = missing
            end
        end
    end

    #####################################
    ## measure peak and trough amplitudes
    peakAmps = missings(Float64, nPeaks)
    troughAmps = missings(Float64, nPeaks)
    peakAmpsEst = missings(Float64, nPeaks)
    troughAmpsEst = missings(Float64, nPeaks)
    for i=1:nPeaks
        if ismissing(peakPoints[i]) == false
            peakAmps[i] = sig[round(Int, peakPoints[i])]
            peakAmpsEst[i] = sig[round(Int, peakPoints[i])]
        else
            peakAmps[i] = missing
            if ismissing(expPeakLatencies[i]) == false
                #idxStart = round(Int, (expPeakLatencies[i]+winStart[i]-epochStart)*sampRate);
                #idxStop = round(Int, (expPeakLatencies[i]+winStop[i]-epochStart)*sampRate);
                peakAmpsEst[i] = sig[round(Int, peakPointsEst[i])] #maximum(sig[idxStart:idxStop])[1]
            else
                peakAmpsEst[i] = missing
            end
        end

        if ismissing(troughPoints[i]) == false
            troughAmps[i] = sig[round(Int, troughPoints[i])]
            troughAmpsEst[i] = sig[round(Int, troughPoints[i])]
        else
            troughAmps[i] = missing
            if ismissing(minTroughLatEst[i]) == false
                #idxStart = round(Int, (minTroughLat[i]-epochStart)*sampRate);
                #idxStop = round(Int, (maxTroughLat[i]-epochStart)*sampRate);
                troughAmpsEst[i] = sig[round(Int, troughPointsEst[i])] #minimum(sig[idxStart:idxStop])[1]
            else
                troughAmpsEst[i] = missing
            end
        end
    end

    ####################################
    ## measure peak-trough amplitudes
    peakTroughAmps = peakAmps-troughAmps
    peakTroughAmpsEst = peakAmpsEst-troughAmpsEst

    ####################################
    ## store results in a dataframe
    df = DataFrame(wave=waveLabels, peakPoint=peakPoints, troughPoint=troughPoints, peakLatency=peakLatencies, troughLatency=troughLatencies, peakAmp=peakAmps, troughAmp=troughAmps, peakTroughAmp=peakTroughAmps, peakPointEst=peakPointsEst, peakLatencyEst=peakLatenciesEst, peakAmpEst=peakAmpsEst, troughPointEst=troughPointsEst, troughLatencyEst=troughLatenciesEst, troughAmpEst=troughAmpsEst, peakTroughAmpEst=peakTroughAmpsEst,
                   minPeakLat=expPeakLatencies+winStart,
                   maxPeakLat=expPeakLatencies+winStop,
                   minTroughLat=minTroughLat,
                   maxTroughLat=maxTroughLat)

    return df

end


function findInflectionsUp(y::Union{AbstractMatrix{T}, AbstractVector{T}}, sampRate::Real; epochStart::Real=0) where {T<:Real}
    inflPntsUp = (Int)[]

    y = vec(y)
    dy = diff(y)
    ddy = diff(dy)
    sgnch = ddy[1:end-1].*ddy[2:end]
    for i=1:length(sgnch)
        if signbit(sgnch[i]) == true
            if dy[i+1] > 0
                push!(inflPntsUp, i+2)
            end
        end
    end

    inflTimesUp = zeros(length(inflPntsUp))

    for i=1:length(inflPntsUp)
        inflTimesUp[i] = (inflPntsUp[i]-1)/sampRate
    end

    inflTimesUp = inflTimesUp .+ epochStart

    return inflPntsUp, inflTimesUp

end


function findInflectionsDown(y::Union{AbstractMatrix{T}, AbstractVector{T}}, sampRate::Real; epochStart::Real=0) where {T<:Real}
    inflPntsDown = (Int)[]

    y = vec(y)
    dy = diff(y)
    ddy = diff(dy)
    sgnch = ddy[1:end-1].*ddy[2:end]
    for i=1:length(sgnch)
        if signbit(sgnch[i]) == true
            if dy[i+1] < 0
                push!(inflPntsDown, i+2)
            end
        end
    end

    inflTimesDown = zeros(length(inflPntsDown))

    for i=1:length(inflPntsDown)
        inflTimesDown[i] = (inflPntsDown[i]-1)/sampRate
    end

    inflTimesDown = inflTimesDown .+ epochStart

    return inflPntsDown, inflTimesDown

end

## function calcPeakPoints(sig, pkTimes, delay, sampRate; epochStart=0)

##     xpksPnts, xpksTimes = findPeaks(sig, sampRate, epochStart=epochStart)
##     xtrsPnts, xtrsTimes = findTroughs(sig, sampRate, epochStart=epochStart)

##     peakTimes = (Float64)[]
##     peakPoints = (Float64)[]
##     for pk=1:length(pkTimes)
##         fooTimes = xpksTimes[(xpksTimes .>= pkTimes[pk]+delay[1]) .& (xpksTimes .<= pkTimes[pk]+delay[2])]
##         fooPnts = xpksPnts[(xpksTimes.>=pkTimes[pk]+delay[1]) .& (xpksTimes.<=pkTimes[pk]+delay[2])]
##         if length(fooTimes) > 0
##             #find largest peak
##             idx = findall(sig[fooPnts] .== maximum(sig[fooPnts]))[1]
##             push!(peakTimes, fooTimes[idx])
##             push!(peakPoints, fooPnts[idx])
##         else
##             push!(peakTimes, missing)
##             push!(peakPoints, missing)
##         end
##     end

##     minPeakTroughLat = [0.25/1000, 0.12/1000, 0.25/1000, 0.12/1000, 0.25/1000]
##     maxPeakTroughLatDiff = [0.75/1000, missing, 0.75/1000, missing, 1.5/1000]
##     maxPeakTroughLat = minPeakTroughLat + maxPeakTroughLatDiff

##     troughTimes = (Float64)[]
##     troughPoints = (Float64)[]
##     for pk=1:length(pkTimes)
##         thisMin = minPeakTroughLat[pk]
##         if isnan(maxPeakTroughLat[pk]) == true
##             if pk+1 <= 5
##                 thisMax = peakTimes[pk+1] - peakTimes[pk]
##             end
##         else
##             thisMax = maxPeakTroughLat[pk]
##         end
##         fooTimes = xtrsTimes[(xtrsTimes .>= peakTimes[pk]+thisMin) .& (xtrsTimes .<= peakTimes[pk]+thisMax)]
##         fooPnts = xtrsPnts[(xtrsTimes.>= peakTimes[pk]+thisMin) .& (xtrsTimes.<= peakTimes[pk]+thisMax)]
##         if length(fooTimes) > 0
##             #find largest trough
##             idx = findall(sig[fooPnts] .== minimum(sig[fooPnts]))[1]
##             push!(troughTimes, fooTimes[idx])
##             push!(troughPoints, fooPnts[idx])
##         else
##             push!(troughTimes, missing)
##             push!(troughPoints, missing)
##         end
##     end

##     peakLatencies = peakTimes
##     troughLatencies = troughTimes

##     peakAmps = zeros(length(peakPoints))
##     troughAmps = zeros(length(peakPoints))
##     for i=1:length(peakAmps)
##         if isnan(peakPoints[i]) == false
##             peakAmps[i] = sig[1, round(Int, peakPoints[i])]
##         else
##             peakAmps[i] = missing
##         end

##         if isnan(troughPoints[i]) == false
##             troughAmps[i] = sig[1, round(Int, troughPoints[i])]
##         else
##             troughAmps[i] = missing
##         end
##     end

##     peakTroughAmps = peakAmps-troughAmps
##     prms = []## Dict("minPeakLat" => [minPeakILat, minPeakIILat, minPeakIIILat, minPeakIVLat, minPeakVLat],
##              ##    "maxPeakLat" => [maxPeakILat, maxPeakIILat, maxPeakIIILat, maxPeakIVLat, maxPeakVLat],
##              ##    "minTroughLat" => [minTroughILat, minTroughIILat, minTroughIIILat, minTroughIVLat, minTroughVLat],
##              ##    "maxTroughLat" => [maxTroughILat, maxTroughIILat, maxTroughIIILat, maxTroughIVLat, maxTroughVLat])

##     return peakPoints, troughPoints, peakLatencies, troughLatencies, peakAmps, troughAmps, peakTroughAmps, prms

## end
