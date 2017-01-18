module IRImodels
    using CategoricalArrays, DataArrays, DataFrames, DataStructures, Distributions,
        GLM, JSON, MixedModels, NamedArrays
    using StatsBase: CoefTable, counts
    export
        checksingularity,
        dropterm,
        drop1,
        GLMformula,
        GLMMform,
        profileθ,
        relevel,
        vif

    include("AoD.jl")

end # module
