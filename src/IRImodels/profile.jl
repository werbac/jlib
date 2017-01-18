function _profile(m::GeneralizedLinearModel, ind::Integer)
    pred = m.pp
    Xorig = pred.X
    remainder = ones(Bool, size(Xorig, 2))
    remainder[ind] = false
    Xreduced = Xorig[:, remainder]
    start = coef(m)[remainder]
    y = m.rr.y
    d = m.rr.d
    l = Link(m)
    β₀ = coef(m)[ind]
    col = view(Xorig, :, ind)
    stderror = StatsBase.stderr(m)[ind]
    dev₀ = deviance(m)
    parvals = [β₀]
    zvals = zeros(parvals)
    newm = fit(GeneralizedLinearModel, Xreduced, y, d, l, offset = β₀ * col, start = start)
    cvals = copy(coef(newm))
    function newval(pv)
        push!(parvals, pv)
        newm.fit = false
        copy!(newm.rr.offset, pv * col)
        push!(zvals, sign(step) * sqrt(deviance(fit!(newm, start = start)) - dev₀))
        append!(cvals, coef(newm))
    end
    fact = step = 0.6
    while abs(zvals[end]) < quantile(Normal(), 0.995)
        newval(β₀ + fact * stderror)
        fact += step
    end
    fact = step = -step
    newval(β₀ + fact * stderror)
    while abs(zvals[end]) < quantile(Normal(), 0.995)
        newval(β₀ + fact * stderror)
        fact += step
    end
    p = sortperm(parvals)
    parvals[p], zvals[p], transpose(reshape(cvals, (length(start), length(p))))[p,:]
end
