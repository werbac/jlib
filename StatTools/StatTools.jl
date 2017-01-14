isdefined(Base, :__precompile__) && __precompile__()

#__precompile__()


module StatLib

using Compat, DataArrays, GLM, DataFrames, Distributions, NLopt, Showoff, StatsBase, DataStructures, StatsFuns, JuMP, DBAPI, JavaCall, JDBC
using JSON, Requests, HttpParser

import DataArrays
import DBAPI, JavaCall, JDBC
import StatsBase: coef, coeftable, df, deviance, fit!, fitted, loglikelihood, model_response, nobs, vcov
import DataFrames: @~




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
       rest

import Base: ==, *

include("rest.jl")

end # module
