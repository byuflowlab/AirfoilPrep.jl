"""
    findclmax(aoa, cl)
Determines clmax as greatest lift coefficient. Returns aoaclmax, clmax
"""
function findclmax(aoa, cl)
    indclmax = argmax(cl)
    clmax = cl[indclmax]
    aoaclmax = aoa[indclmax]

    return aoaclmax, clmax
end #findclmax

"""
    findclmax(aoa, cl, npnts::Integer)
Determines clmax as when all angles greater than the clmax angle for `npnts`
data points are less than clmax. Assumes points are ordered by angle of attack.
Returns aoaclmax, clmax
"""
function findclmax(aoa, cl, npnts::Integer)
    # default to basic method if conditions not satisfied
    aoaclmax,clmax = findclmax(aoa,cl)

    # start at zero
    idxzero = argmin(abs.(aoa))
    for i = idxzero:length(cl)-npnts
        # check if all points in range are less than current point
        if count(cl[(i+1):(i+npnts)] .> cl[i]) == 0
            clmax = cl[i]
            aoaclmax = aoa[i]
            break
        end
    end

    return aoaclmax,clmax
end #findclmax

"""
    findclmax(aoa, cl, range::AbstractFloat)
Determines clmax as when all angles greater than the clmax angles for `range`
radians are less than clmax. Assumes points are ordered by angle of attack.
Returns aoaclmax, clmax
"""
function findclmax(aoa, cl, range::AbstractFloat)
    # default to basic method if conditions not satisfied
    aoaclmax, clmax = findclmax(aoa, cl)

    # start at zero
    idxzero = argmin(abs.(aoa))
    for i = idxzero:length(cl)
        # find index of point approximately `range` greater
        idx = argmin(abs.(aoa[i] .+ range .- aoa))
        if idx == i
            break
        else
            # check if all points in range are less than current point
            if count(cl[(i+1):idx] .> cl[i]) == 0
                clmax = cl[i]
                aoaclmax = aoa[i]
                break
            end
        end
    end

    return aoaclmax, clmax
end #findclmax


"""
    findclmin(aoa, cl)
Determines clmin as smallest lift coefficient.  Returns aoaclmin, clmin
"""
function findclmin(aoa, cl)
    indclmin = argmin(cl)
    clmin = cl[indclmin]
    aoaclmin = aoa[indclmin]

    return aoaclmin,clmin
end #findclmax

"""
    findclmin(aoa, cl, npnts::Integer)
Determines clmin as when all angles less than the clmin angle for `npnts` data
points are greater than clmin. Assumes points are ordered by angle of attack.
Returns aoaclmin, clmin
"""
function findclmin(aoa, cl, npnts::Integer)
    # default to basic method if conditions not satisfied
    aoaclmin, clmin = findclmin(aoa, cl)

    # start at zero
    idxzero = argmin(abs.(aoa))
    for i = idxzero:-1:(npnts+1)
        # check if all points in range are less than current point
        if count(cl[(i-1):-1:(i-npnts)] .< cl[i]) == 0
            clmin = cl[i]
            aoaclmin = aoa[i]
            break
        end
    end

    return aoaclmin, clmin
end #findclmax

"""
    findclmin(aoa, cl, range::AbstractFloat)
Determines clmin as when all points less than the clmin angle for `range`
radians are greater than clmin. Returns aoaclmin, clmin
"""
function findclmin(aoa, cl, range::AbstractFloat)
    # default to basic method if conditions not satisfied
    aoaclmin, clmin = findclmin(aoa, cl)

    # start at zero
    idxzero = argmin(abs.(aoa))
    for i = idxzero:-1:1
        # find index of point approximately `range` less
        idx = argmin(abs.(aoa[i] .- range .- aoa))
        if idx == i
            break
        else
            # check if all points in range are greater than current point
            if count(cl[(i-1):-1:idx] .< cl[i]) == 0
                clmin = cl[i]
                aoaclmin = aoa[i]
                break
            end
        end
    end

    return aoaclmin, clmin
end #findclmax

"""
    fitliftslope(aoa, cl, tol::Real=0.05, allnegfit::Bool=false)
Returns the lift slope and zero lift angle of attack as well as the data used
to compute these parameters.  The tolerance `tol` (in terms of cl) for the fit
and whether to include all data points below zero angle of attack `allnegfit`
may be specified. Returns liftslope,zeroliftangle,aoafit,clfit
"""
function fitliftslope(aoa, cl, tol::Real=0.05, allnegfit::Bool=false, center=0.0)

    if length(aoa) != length(cl)
        error("aoa and cl must have the same length")
    end

    # split into positive and negative angles of attack
    idxcenter = argmin(abs.(aoa .- center))
    idxneg = findall(aoa .< aoa[idxcenter])
    sort!(idxneg,rev = true)
    idxpos = findall(aoa .> aoa[idxcenter])
    sort!(idxpos)

    # initialize arrays of data to fit
    aoafit = [aoa[idxcenter]]
    clfit = [cl[idxcenter]]
    # create indexes for increasing and decreasing aoa
    ipos = 1
    ineg = 1
    # initialize flag to stop adding data > 0
    upperflag = false
    if ipos > length(idxpos)
        upperflag = true
    end
    # initialize flag to stop adding data < 0
    lowerflag = false
    if ineg > length(idxneg)
        lowerflag = true
    end
    # compute lift slope and zero lift angle of attack
    success = false
    liftslope = 2pi
    zeroliftangle = 0.0
    for i = 1:(length(idxpos) + length(idxneg))
        # add positive angle of attack point to fit if closest to zero aoa
        if upperflag == false
            if lowerflag == true || abs(aoa[idxpos[ipos]]) <= abs(aoa[idxneg[ineg]])
                # add point
                push!(aoafit, aoa[idxpos[ipos]])
                push!(clfit, cl[idxpos[ipos]])
                # order points
                order = sortperm(aoafit)
                aoafit = aoafit[order]
                clfit = clfit[order]
                # create least squares fit
                lsqsol = hcat(aoafit, -ones(length(aoafit))) \ clfit
                # if liftslope not defined from previous iteration, assign it here
                if ipos == 1
                    liftslope = lsqsol[1]
                end
                # get slope and error of added point
                yerror = lsqsol[1] * (aoa[idxpos[ipos:end]] .- lsqsol[2] / lsqsol[1]) .- cl[idxpos[ipos:end]]
                # slope = (clfit[end] - clfit[end-1])/(aoafit[end] - aoafit[end-1])
                # check if on linear portion of lift curve
                if any(yerror .< tol) #slope > 0.5*liftslope ||
                    liftslope = lsqsol[1]
                    zeroliftangle = lsqsol[2] / lsqsol[1]
                    success = true
                else # if not on linear portion remove from fit
                    aoafit = aoafit[1:(end-1)]
                    clfit = clfit[1:(end-1)]
                    upperflag = true
                end
                #  increment counter and check for more points greater than aoa=0
                ipos = ipos+1
                if ipos > length(idxpos)
                  upperflag = true
                end
            end
        end
        # add negative angle of attack point to fit if closest to zero
        if lowerflag == false
            if upperflag == true || abs(aoa[idxneg[ineg]]) <= abs(aoa[idxpos[ipos]])
                # add point
                pushfirst!(aoafit, aoa[idxneg[ineg]])
                pushfirst!(clfit, cl[idxneg[ineg]])
                # create least squares fit
                lsqsol = hcat(aoafit, -ones(length(aoafit))) \ clfit
                # get error of added point
                yerror = lsqsol[1] * (aoa[idxneg[ineg:end]] .- lsqsol[2]/lsqsol[1]) .- cl[idxneg[ineg:end]]
                # check if on linear portion of lift curve
                if any(abs.(yerror) .< tol) || allnegfit
                    liftslope = lsqsol[1]
                    zeroliftangle = lsqsol[2] / lsqsol[1]
                    success = true
                else # if not on linear portion remove from fit
                    aoafit = aoafit[2:end]
                    clfit = clfit[2:end]
                    lowerflag = true
                end
                #  increment counter and check for more points less than aoa=0
                ineg = ineg + 1
                if ineg > length(idxneg)
                    lowerflag = true
                end
            end
        end
        # break if no more points to add
        if upperflag && lowerflag
            break
        end
    end
    if success == false
        liftslope, zeroliftangle, aoafit, clfit = fitliftslope(aoa, cl, tol, allnegfit, aoa[idxcenter]+1.0*pi/180)
    end

    return liftslope, zeroliftangle, aoafit, clfit
end

"""
    findclmaxlinear(aoa, cl, liftslope::Real, zeroliftangle::Real; tol::Real=0.1, 
    interpolate::Bool=true)
Returns the angle of attack and cl of the maximum linear lift coefficient. `tol`
is used to determine the liftslope line tolerance. Returns aoaclmaxlinear,clmaxlinear
"""
function findclmaxlinear(aoa, cl, liftslope::Real, zeroliftangle::Real; 
                          tol::Real=0.1, interpolate::Bool=true)

    # find where difference in cl is greater than tol
    indclmaxlinear = 0
    idxzero = argmin(abs.(aoa))
    for i = idxzero:length(aoa)
        predictedcl = liftslope * (aoa[i] - zeroliftangle)
        indclmaxlinear = i
        if abs(cl[i]-predictedcl) > tol
            break
        end
    end
    if interpolate
        m = (cl[indclmaxlinear] - cl[indclmaxlinear-1]) / (aoa[indclmaxlinear] - aoa[indclmaxlinear-1])
        aoaclmaxlinear = (m*aoa[indclmaxlinear-1] - liftslope * zeroliftangle - tol - cl[indclmaxlinear-1]) / (m - liftslope)
        clmaxlinear = liftslope * (aoaclmaxlinear - zeroliftangle) - tol
    else
        aoaclmaxlinear = aoa[indclmaxlinear]
        clmaxlinear = cl[indclmaxlinear]
    end

    return aoaclmaxlinear,clmaxlinear
end #findclmaxlinear
