
#Generate GLM Formula
#function genFmula(iv::Array{Symbol},isLastOffset::Bool=false)  # Y=1st, if isLastOffset:: offset last
#    isLastOffset ?  genFmula(iv[1],iv,last(iv)) : genFmula(iv[1],iv)
#end
#function genFmula(y::Symbol, iv::Array{Symbol},offset::Symbol=:Missing)  
#    xvars=setdiff(iv,[y,offset])
#    eval(parse( string(y)*" ~ 1"* reduce(*, [ " + "*  string(c) for c in xvars ] )  ))
#end
#function genFmula(y::Symbol, iv::Array{Symbol})  
#    xvars=setdiff(iv,[y])
#    eval(parse( string(y)*" ~ 1"* reduce(*, [ " + "*  string(c) for c in xvars ] )  ))
#end
function genF(y::Symbol, iv::Array{Symbol},ranef::Array{Symbol}=Symbol[])  
        vars=setdiff(iv,vcat([y],ranef))
        f = string(y)*" ~ 1"
        f = length(vars) > 0 ?  f * reduce(*, [ " + "*  string(c) for c in unique(vars) ] ) : f
        f = length(ranef) > 0 ?  f * reduce(*, [ " + "*  "(1 | "*string(c)*")" for c in unique(ranef) ] )  : f
        eval(parse(f))
end

# remove single level vars
function FS_singleLevel(dfd::DataFrame,vars::Array{Symbol}=Symbol[])
    if length(vars)==0  vars=names(dfd)  end
    return [c for c in filter(x -> length(unique( dfd[x] ))==1, vars)]
end


function raneffect(m::MixedModels.GeneralizedLinearMixedModel)
    DataFrame(parameter=levels(m.LMM.trms[1].f), coef=vec(ranef(m, named=true)'[1]), stderr = sqrt.(vec(condVar(m)[1])) )  # stderrX =vec(condVar(m)[1])
end

isDict(::Dict) = true
#isDict(::OrderedDict) = true
isDict(::Any) = false
isArray(::AbstractArray) = true
isArray(::Any) = false
isNum{T<:Number}(::AbstractArray{T}) = true
isNum(::Any) = false
isFloat{T<:Float64}(::AbstractArray{T}) = true
isFloat(::Any) = false
isInt{T<:Int64}(::AbstractArray{T}) = true
isInt(::Any) = false
isString{T<:String}(::AbstractArray{T}) = true
isString(::String) = true
isString(::Any) = false
isBool{T<:Bool}(::AbstractArray{T}) = true
isBool(::Any) = false
isSymbol{T<:Symbol}(::AbstractArray{T}) = true
isSymbol(::Any) = false


function blank(dfx::DataFrame,d::Dict=Dict())
    arr = Any[]      
    for (n, v) in eachcol(dfx)
        if haskey(d,n)
            arr = vcat(arr, [d[n]])  
        else
            if isFloat(v)
                arr = vcat(arr, [0.0])
            elseif isInt(v)
                arr = vcat(arr, [0])
            elseif isString(v)
                arr = vcat(arr, [""])
            elseif isBool(v)
                        arr = vcat(arr, [false])
            elseif isSymbol(v)
                        arr = vcat(arr, [:empty])
            else
                println("ERROR : datatype Not Found!!!!!!!!")
            end
        end     
    end 
    return  arr
end 


Base.lowercase(df::DataFrame) = names!(df, map(x->Symbol(lowercase(string(x))),names(df)))
Base.lowercase(da::DataArray{String}) = [lowercase(x) for x in da]


function filler(dfx::DataFrame, rowIn::Array{Any})
    row=filler(dfx)
    for i in 1:length(rowIn)
        row[i] = rowIn[i]
    end
    return row
end

function filler(dfx::DataFrame, cols::Vector{Symbol}=Symbol[])
    a = Any[]
    if length(cols)==0 cols=names(dfx) end
    for (n, v) in eachcol(dfx)
        if n ∈ cols && isFloat(v)
            push!(a, 0.0)
        elseif n ∈ cols && isInt(v)
            push!(a, 0)
        elseif n ∈ cols && isString(v)
            push!(a, "")
        elseif n ∈ cols && isBool(v)
            push!(a, false)
        end
    end
    return a
end


function cov2cor(m)
    d = Diagonal([inv(sqrt(m[i, i])) for i in 1:size(m, 2)])
    A_mul_B!(d, m * d)
end


function vif2(m)
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

function vifDF2(m::RegressionModel)
    sdf=coefDF(m,true)
    v=vif2(m)
    vdf=DataFrame(parameter=convert(Array{Symbol},names(v)[1]), vif=convert(Array{Float64}, v)[:,1]   )
    #vdf[:parameter] = convert(Array{Symbol}, vdf[:variable])
    join(sdf,vdf[[:parameter,:vif]], on = :parameter, kind=:outer)
end



""" Doug
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



function gethostworkers()
    d = Dict()
    @sync @async for (idx, pid) in enumerate(workers())
        d[pid] = remotecall_fetch(getipaddr,pid)
    end
    iout=Dict()
    for ip in [string(ip) for ip in unique(values(d))]
        a=Int64[]
        for (key, value) in d
            if string(value) == ip
                push!(a,key)
            end
        end    
        iout[ip] = sort(a)
    end
    l = get(iout,string(getipaddr()),[])
    iout[string(getipaddr())] = convert(Array{Int64},vcat([1],l))
    return iout
end


function gethostModelDistribution(raneff::Array{Symbol})
    function allocateCore(a::Array{Int64})
        proc = a[1]
        ao = length(setdiff(a,proc)) > 0 ?  setdiff(a,proc) : a
        return proc, ao
    end
    nodehosts = gethostworkers()
    hosts = [k for k in keys(nodehosts)]
    hostcnt = length(hosts)
    #raneff = [:creative, :publisher1, :placement]
    h=1
    m=:pen
    penprocs = Dict()
    for r in raneff
        proc, nodehosts[hosts[h]] = allocateCore(nodehosts[hosts[h]])
        penprocs[r] = Dict(:proc=>proc, :host=>hosts[h])
        println(r,"  :  ", hosts[h], "  ",proc)
        h = h == hostcnt ? h=1 : h+1
    end
    
    static_nodehosts = deepcopy(nodehosts)
    
    println("Pen : ",penprocs )
    println(nodehosts)
end



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


function coefDF(m::RegressionModel, paramsAsSymbols::Bool=false)
    vout = DataFrame( #vars=vcat([:intercept],m.mf.terms.terms)  #g.model.mm.assign # m.mm.assign
               #vars=vcat([:intercept],[m.mf.terms.terms[x] for x in m.mm.assign[2:end]])
               parameter=DataFrames.coefnames(m.mf)
             , coef=GLM.coef(m)
             , stderr=stderr(m)
             , zval=GLM.coef(m)./stderr(m) 
             , pval= ccdf(FDist(1, dof_residual(m)), abs2(GLM.coef(m)./stderr(m))  )
             )
    if paramsAsSymbols == true
        vout[:parameter] = map(x->Symbol(replace(replace(x,"group: 1","group"),"(Intercept)","Intercept")), vout[:parameter])
    end
    return vout
end


"""
Relevel DF for GLM
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
#dfd[ranef] = relevel(dfd[ranef], "none")
 
 
#GLM REsiduals
function xResiduals(g::DataFrames.DataFrameRegressionModel)
    resp = g.model.rr
    #sign(resp.y - resp.mu) .* sqrt(resp.devresid)    # Original
    #sign(resp.y - resp.mu) .* sqrt(complex(resp.devresid))
    sign(resp.y - resp.mu) .* sqrt(((resp.devresid)+1)-1)
    #sign(resp.y - resp.mu) .* sqrt(abs(resp.devresid))
end


function unpool(dfd::DataFrame)
    v_factors = Symbol[]
    for c in names(dfd)
        if typeof(dfd[c]) in [ PooledDataArray{String,UInt8,1} ]
            println("         Converting : ", c )
            push!(v_factors,c)
            dfd[c] = Array(dfd[c])
        end
    end
    return v_factors
end


function poolit!(dfd::DataFrame, vars::Array{Symbol}=Symbol[])
    if length(vars) == 0   vars=names(dfd) end
    for (i, c) in enumerate(dfd.columns)
        if typeof(c) == DataVector{String}
            if dfd.colindex.names[i] in vars
                #println("[",i,"]  : ",dfd.colindex.names[i])
                #unique(dfd[i])
                dfd[i] = pool(dfd[i])
            end
        end
        #typeof(c) == DataVector{String} && println("[",i,"]  : ",dfd.colindex.names[i])
    end
end


function getColswithType(dtype::String, dfd::DataFrame, vars::Array{Symbol}=Symbol[])
    clst=Symbol[]
    typLst = Dict("num" => [DataVector{Float64},DataVector{Int64}], "str" => [DataVector{String}] )
    if length(vars) == 0   vars=names(dfd) end
    if dtype in ["num","str"]
        for (i, c) in enumerate(dfd.columns)
            if dfd.colindex.names[i] in vars
                    if typeof(c) in typLst[dtype]
                #if typeof(c) in [ DataVector{Float64}, DataVector{Int64} ] # = = DataVector{String}
                    #println("[",i,"]  : ",dfd.colindex.names[i])
                    push!(clst,dfd.colindex.names[i])
                end
             end
         end 
    else
        println("ERROR: getColswithType : Invalid dtype!")
    end
    return clst
end
 
function corDFD(dfd::DataFrame, vars::Array{Symbol}=Symbol[])
    if length(vars) == 0   
        vars=names(dfd)  
    #else 
    #    vars=convert(Array{Symbol},vars) 
    end
    cols=getColswithType("num", dfd, convert(Array{Symbol},vars) )
    carr= cor(Array(dfd[cols]))  
    tcarr = triu(carr,1)
    d = DataFrame(tcarr) 
    names!(d,cols)
    d[:vars] = cols
    return stack(d,cols)
end
#function corDFD(dfd::DataFrame, vars::Array{Any}=[])
#    varsS=Symbol[]
#    if length(vars) > 0 varsS=convert(Array{Symbol},vars) end
#    corrDF(dfd,varsS)
#end





function cmparr(a1::Array{Symbol},a2::Array{Symbol},ignore::Array{Symbol}=Symbol[])
    a1=setdiff(a1,ignore)
    a2=setdiff(a2,ignore)
    a1_only=setdiff(a1,a2)
    a2_only=setdiff(a2,a1)
    l1=length(a1_only)
    l2=length(a2_only)
    println(l1,"--",l2)
    l1>l2 ? DataFrame(a1_only=a1_only,a2_only=vcat(a2_only,fill(:x, l1-l2))) : DataFrame(DataFrame(a1_only=vcat(a1_only,fill(:x, l2-l1)),a2_only=a2_only))    
end

function cmparr(a1::Array{Any},a2::Array{Symbol},ignore::Array{Symbol}=Symbol[])
    cmparr(convert(Array{Symbol},a1), a2)
end



