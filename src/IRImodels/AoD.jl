"""
    drop1(m::DataFrameRegressionModel, terms::Vector{Symbol}=Symbol[])

Perform a `drop1` analysis of deviance on `terms` in model `m`

Model `m` is modified and refit dropping the terms in `terms` one at a time.
An analysis of deviance is applied to the refit model.
"""
function drop1(m, terms::Vector{Symbol}=Symbol[])
    form = Formula(m.mf.terms)
    basedev = deviance(m)
    deviances, degsfreedom, Χ², pvalues, nms = [basedev], [0], [0.], [NaN], ["<none>"]
    ncoef = length(coef(m))
    dat = m.mf.df
    dinst = m.model.rr.d
    if length(terms) == 0
        terms = filter(x -> isa(x, Symbol) || isa(x, Expr), form.rhs.args[2:end])
    end
    for trm in terms
        newform = dropterm(form, trm)
        newm = fit(GeneralizedLinearModel, newform, dat, dinst, Link(m))
        newdev = deviance(newm)
        newdof = ncoef - length(coef(newm))
        newΧ² = newdev - basedev
        push!(deviances, newdev)
        push!(degsfreedom, newdof)
        push!(Χ², newΧ²)
        push!(pvalues, ccdf(Chisq(newdof), newΧ²))
        push!(nms, string(trm))
    end
    CoefTable([deviances, Χ², degsfreedom, pvalues], ["Deviance", "Χ²", "degfree", "P[> Χ²]"], nms)
end

function StatsBase.counts(v::Union{CategoricalVector,PooledDataVector})
    levs = isa(v, CategoricalArray) ? CategoricalArrays.levels(v) : DataArrays.levels(v)
    counts = zeros(Int, length(levs) + 1)
    for ind in v.refs
        counts[ind + 1] += 1
    end
    counts[1] == 0 ? shift!(counts) : unshift!(levs, "#NULL")
    result = NamedArray(counts)
    NamedArrays.setnames!(result, levs, 1)
    result
end
GLM.Link{T<:GeneralizedLinearModel}(m::DataFrames.DataFrameRegressionModel{T}) = Link(m.model)

"""
    dropterm(f::Formula, trm::Symbol)

Return a copy of `f` without the term `trm`.

# Examples
```jl
julia> dropterm(foo ~ 1 + bar + baz, :bar)
foo ~ 1 + baz
```
"""
function dropterm(f::Formula, trm::Symbol)
    rhs = copy(f.rhs)
    args = rhs.args
    if !(Meta.isexpr(rhs, :call) && args[1] == :+ && (tpos = findlast(args, trm)) > 0)
        throw(ArgumentError("$trm is not a summand in `$rhs`"))
    end
    Formula(f.lhs, Expr(:call, :+, deleteat!(args, [1, tpos])...))
end

"""
    checksingularity(form::Formula, data::DataFrame, tolerance = 1.e-8)

Return a vector of terms in `form` that produce singularity in the model matrix
"""
function checksingularity(form::Formula, data::DataFrame, tolerance = 1.e-8)
    mf = ModelFrame(form, data)
    mm = ModelMatrix(mf)
    qrf = qrfact!(mm.m, Val{true})
    vals = abs.(diag(qrf[:R]))
    firstbad = findfirst(x -> x < min(tolerance, 0.5) * vals[1], vals)
    if firstbad == 0
        return Symbol[]
    end
    mf.terms.terms[view(mm.assign[qrf[:p]], firstbad:length(vals))]
end

"""
    GLMformula(dd::Dict{String,Any})

Create the GLM formula from a configuation dictionary for the model
"""
function GLMformula(dd::Dict)
    haskey(dd, "y_var") || throw(ArgumentError("dd should have a \"y_var\" key"))
    haskey(dd, "finalvars") || throw(ArgumentError("dd should have a \"finalvars\" key"))
    rhs = :(1 + x)
    pop!(rhs.args)
    append!(rhs.args, Symbol.(dd["finalvars"]))
    Formula(Symbol(dd["y_var"]), rhs)
end

"""
    GLMMform(form::Formula, raneff::Symbol)

Replace the last term in the formula by :(1 | ranef)
"""
function GLMMform(form::Formula, raneff::Symbol)
    retrm = :(1 | x)
    retrm.args[end] = raneff
    form.rhs.args[end] = retrm
    form
end

"""
    profileθ(m::GeneralizedLinearMixedModel)

Profile the deviance as a function of θ for models with a single, scalar r.e. term
"""
function profileθ{T}(m::GeneralizedLinearMixedModel{T})
    θ = getθ(m)
    isa(θ, Vector{T}) && length(θ) == 1 ||
        throw(ArgumentError("m must have a single scalar r.e. term"))
    th = zeros(T, 1)
    δ = zeros(T, 1)
    fill!(θ, 0)
    decreasing = false
    basedev = pirls!(setθ!(m, θ), true)
    for t in 0.001:0.001:2.0
        θ[1] = t
        dd = pirls!(setθ!(m, θ), true) - basedev
        sgn = sign(dd - δ[end])
        push!(δ, dd)
        push!(th, t)
        dd ≤ 4.0 || break
        decreasing && sgn > 0 && break
        decreasing = sgn < 0
    end
    th, δ
end

"""
    relevel{T,I}(v::PooledDataVector{T,I}, lev::T)

Create a new `PooledDataVector{T,I}` with `lev` as the first (i.e. reference) level
"""
function relevel{T,I}(v::PooledDataVector{T,I}, lev::T)
    pool = v.pool
    i = convert(I, findfirst(pool, lev))
    i ≠ 0 || throw(ArgumentError("\"$lev\" is not in v.pool"))
    if i == 1
        return copy(v)
    end
    orig = convert(Vector{I}, 1:length(pool))
    perm = append!([i], deleteat!(orig, i))
    PooledDataArray(DataArrays.RefArray(invperm(perm)[v.refs]), v.pool[perm])
end

function cov2cor(m)
    d = Diagonal([inv(sqrt(m[i, i])) for i in 1:size(m, 2)])
    A_mul_B!(d, m * d)
end


"""
    vif(m::DataFrameRegressionModel)

The variance inflation factors for the terms in a model.
"""
function vif(m)
    v = vcov(m)
    assign = m.mm.assign
    nms = coefnames(m.mf)
    if (ipos = findfirst(nms, "(Intercept)")) > 0
        inds = deleteat!(collect(eachindex(nms)), ipos)
        v = view(v, inds, inds)
        assign = view(assign, inds)
    else
        warning("No intercept: vifs may not be sensible.")
    end
    terms = convert(Vector{Symbol}, m.mf.terms.terms)
    n_terms = length(terms)
    if n_terms < 2
        error("model contains fewer than 2 terms")
    end
    R = cov2cor(v)
    detR = det(Symmetric(R))
    result = NamedArray(Array(eltype(R), (length(terms), 3)))
    setnames!(result, string.(terms), 1)
    setnames!(result, ["GVIF", "Df", "GVIF^(1/(2*Df))"], 2)
    inds = collect(eachindex(terms))
    ainds = collect(eachindex(assign))
    for i in inds
        subs = find(x -> x == i, assign)
        notsubs = setdiff(ainds, subs)
        result[i, 1] = det(Symmetric(R[notsubs, notsubs])) / detR
        if (result[i, 2] = length(subs)) > 1
            result[i, 1] *= det(Symmetric(R[subs, subs]))
        end
    end
    if all(result[:, 2] .== 1)
        return result[:, 1:1]
    end
    result
end
