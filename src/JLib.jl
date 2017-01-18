isdefined(Base, :__precompile__) && __precompile__()

#__precompile__()


module JLib

using JSON, Compat, DataArrays, GLM, DataFrames, Distributions, NLopt, Showoff, StatsBase, DataStructures, StatsFuns, JuMP, MixedModels, NamedArrays, NLsolve
#using JSON, Requests, HttpParser

import DataArrays
#import DBAPI, JavaCall, JDBC
import StatsBase: coef, coeftable, df, deviance, fit!, fitted, loglikelihood, model_response, nobs, vcov
import Base: cond, std
import Distributions: Bernoulli, Binomial, Poisson, Gamma
import GLM: LogitLink, LogLink, InverseLink
import DataFrames: @~

export
       @~,
       lowercase,
       grep,
       grep1,
       lowercase!,
       loadCFG,
       df2dict,
       FS_singleLevel,
       genF,
       retypeBools!,
       poolit!, unpool,
       getColswithType,
       corDFD,
       cmparr,
       vif,
       vifDF, vifDF2,
       coefDF,
       checksingularity,
       xResiduals,
       gethostModelDistribution,
       gethostworkers,
       genF,
       raneffect,
       relevel,
       profileθ,
       readModels, expandM,  #saveModels,
       read_dfx, save_dfx, read_modelsDict, save_modelsDict, read_dfd, save_dfd, read_dfr, save_dfr,
       isnumeric, isfloat, isint, isString, isBool, 
       isDict, isArray, isNum, isFloat, isInt, isSymbol,
       filler,
       ZDict, CIs_O, calcPValue_Opt, np_val,
       dict2json, dict_Sym, mergeDict,  savejson, save_cfg, read_cfg
 

import Base: ==, *

include("common.jl")
include("featureSelection.jl")
include("CIs.jl")
#include("diagnostics.jl")
##include("misc.jl")
#include("modsel.jl")




""" x
# vif(m::DataFrameRegressionModel)  --  The variance inflation factors for the terms in a model.

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
"""


function vifDF(m::RegressionModel)
    sdf=coefDF(m)
    vdf=vif(m)
    vdf[:parameter] = convert(Array{Symbol}, vdf[:variable])
    join(sdf,vdf[[:parameter,:vif]], on = :parameter, kind=:outer)
end

#dfrm = glm( exposed_flag ~ 1 + Prd_2_Qty_POS, df_in  , Poisson(), LogLink() )

function vif(dfrm::RegressionModel)
    X = dfrm.mf.df[2:end]
    rhs = extractrhs(dfrm)
    result = DataFrame(variable=rhs.rhsarray, vif=0.0)
    i = 1
    for var in rhs.rhsarray 
        lhs = parse(var)
        rhsnew = replace(rhs.rhsstring, "+"*var, "")
        rhsnew = rhsnew[2:end]
        rhsnew = parse(rhsnew)
        fnew = Formula(lhs, rhsnew)
        newfit = fit(LinearModel, fnew, X)
        r2 = rsquared(newfit)
        result[:vif][i] = 1 / (1-r2)
        i = i + 1
    end
    result
end


function rsquared(dfrm::RegressionModel)
    SStot = sum((dfrm.model.rr.y - mean(dfrm.model.rr.y)).^2)
    SSres = sum((dfrm.model.rr.y - dfrm.model.rr.mu).^2)
    return (1-(SSres/SStot))
end



type extractrhs
    rhsarray::Array
    rhsstring::AbstractString
end


# Extract the right hand side of a RegressionModel formula
function extractrhs(dfrm::RegressionModel)
    regressors = dfrm.mf.terms.terms
    rhsarray = similar(regressors, AbstractString)
    for i in 1:length(regressors)
        rhsarray[i] = string(regressors[i])
        rhsarray[i] = replace(rhsarray[i], " ", "")
    end
    rhsstring = ""
    for var in rhsarray
        rhsstring = rhsstring*"+"*var
    end
    if !dfrm.mf.terms.intercept
        rhsarray = [0, rhsarray]
        rhsstring = "0+"*rhsstring
    end
    return extractrhs(rhsarray, rhsstring)
end





end # module
