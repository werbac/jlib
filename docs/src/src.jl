
# ------- R 2 Julia Conversion ---------------
#function saveModels(d::Dict,fname::String)
#    j = JSON.json([ d])
#    open(fname,"w") do f
#        write(f,j)
#    end
#end

function expandM(di::Dict)
    d2=deepcopy(di)
    if d2[:modelName]==:occ  
        d2[:dist] = Poisson()
        d2[:lnk] = LogLink()
        d2[:Buyer_Pos_P1_is1] = true
    end
    if d2[:modelName]==:dolocc  
        d2[:dist] = Gamma()
        d2[:lnk] = LogLink()
        d2[:Buyer_Pos_P1_is1] = true
    end
    if d2[:modelName]==:pen
        d2[:dist] = Bernoulli()
        d2[:lnk] = LogitLink()
        d2[:Buyer_Pos_P1_is1] = false
    end
    return d2
end



function readModels(fname::String)
    if isfile(fname)
        d=Dict()
        j = JSON.parsefile(fname)[1]
        #for (key, value) in dict
        for (k, v) in j    
            d2=Dict()
            if Symbol(k)==:mocc  
                d2[:dist] = Poisson()
                d2[:lnk] = LogLink()
                d2[:Buyer_Pos_P1_is1] = true
            end
            if Symbol(k)==:mdolocc  
                d2[:dist] = Gamma()
                d2[:lnk] = LogLink()
                d2[:Buyer_Pos_P1_is1] = true
            end
            if Symbol(k)==:mpen  
                d2[:dist] = Bernoulli()
                d2[:lnk] = LogitLink()
                d2[:Buyer_Pos_P1_is1] = false
            end
            
            if typeof(v) == Array{Any,1}
                d[Symbol(k)] = [Symbol(vx) for vx in v]
            else
                for (k1, v1) in v
                    if typeof(v1) in [ Array{Any,1} ]
                        v2=[Symbol(v1x) for v1x in v1]
                    else
                        v2=Symbol(v1)
                    end
                    d2[Symbol(k1)] = v2
                end
                d[Symbol(k)] = d2
            end
        end
        for (k,v) in d
            if k in [:idolocc, :ipen, :iocc]
                d[k] = expandM(v)
            end
        end
        return d
    else
        return "file not found!"
    end
end


read_dfx(root::String) = readtable(root*"/dfx.csv")
save_dfx(root::String,dfx::DataFrame) = writetable(root*"/dfx.csv", dfx) 
read_modelsDict(root::String) = readModels(root*"/modelsDict.json")
function save_modelsDict(root::String, d::Dict)
    j = JSON.json([d])
    open(root*"/modelsDict.json","w") do f
        write(f,j)
    end
end
read_dfd(root::String) = readtable(root*"/dfd.csv",header=true);
save_dfd(root::String,dfd::DataFrame) = writetable(root*"/dfd.csv", dfd)

read_dfr(root::String) = readtable(root*"/scored.csv",header=true);
save_dfr(root::String,dfr::DataFrame) = writetable(root*"/scored.csv", dfr)

# ------- R 2 Julia Conversion - END ---------------



#include("trd.jl")
#function Base.lowercase(df::DataFrame)
function lowercase!(df::DataFrame)
    for (i,v) in enumerate(names(df))   if v != Symbol(lowercase(string(v)))  rename!(df,v,Symbol(lowercase(string(v)))) end end 
end
function Base.lowercase(arr::Array{Symbol})    return  Symbol[Symbol(lowercase(string(s))) for s in arr] end
function Base.lowercase(arr::Array{String})    return  String[lowercase(string(s)) for s in arr] end
function Base.lowercase(arr::Array{String})    return  String[lowercase(string(s)) for s in arr] end
function Base.lowercase(s::Symbol) return Symbol(lowercase(string(s))) end
function Base.lowercase(d::OrderedDict) 
    vout=OrderedDict()
    for (k, v) in d  vout[lowercase(k)] = typeof(v) in [Bool,Float64,Int64] ? v : lowercase(v) end 
    return vout
end
  
function grep(v::String,vars::Array{Symbol})   #function grep(v::AbstractString) filter(x->ismatch(Regex("$(v)"), string(x)), vars) end
    filter(x->ismatch(Regex("$(v)"), string(x)), vars)
end
function grep1(v::String,vars::Array{Symbol})   #function grep(v::AbstractString) filter(x->ismatch(Regex("$(v)"), string(x)), vars) end
    lst = filter(x->ismatch(Regex("$(v)"), string(x)), vars)
    return length(lst)>0 ? lst[1] : nothing 
end

function grep(varr::Array{String},vars::Array{Symbol})
    arr_out=Symbol[]
    for v in varr
        res = filter(x->ismatch(Regex("$(v)"), string(x)), vars)
        arr_out = vcat(arr_out, res)
    end
    return arr_out
end 
#function grep(varr::Array{String},vars::Array{Symbol}) return grep(convert(Array{String},varr),vars) end  #grep(["M","pan"],names(df_in))


function df2dict(idf::DataFrame)
    dout=OrderedDict()
    for n in names(idf)
       dout[n] = idf[1,Symbol(n)]
    end
    return dout
end

function df2dict(row::DataFrameRow)
    dout=OrderedDict()
    for n in names(row)
       dout[n] = row[Symbol(n)]
    end
    return dout
end






"""
app.cfg
P2_Competitor = true
pvalue_lvl = 0.20
excludedKeys = key1,key2,key3
exposed_flag_var = new_exposed_flag
sigLevel = 0.2
random_demos =
random_campaigns =
dropvars = goop,snop
scoring_vars =
TotalModelsOnly=false
"""

function loadCFG(d::OrderedDict,fn::String) 
    function cvt(v::String, t::DataType)
        if t in [Bool]
            return lowercase(v) =="true" ? true : lowercase(v) =="false" ? false : NA
        elseif t in [nothing]
            return NA
        elseif t in [Array{String,1}]
            #return AbstractString[string(x) split(v,",")]
            return String[string(x) for x in split(replace(v,"\"",""),",")]
        elseif t in [Array{Symbol,1}]
            return Symbol[Symbol(x) for x in split(replace(v,"\"",""),",")]
        elseif t in [Symbol]
            return Symbol(v)
        elseif t in [String]
            return string(v)
        elseif t in [Float64]
            return parse(Float64,v)
        else
            return "__ooops__ TYPE not found!"
        end
    end
    if isfile(fn)
        f = open(fn)
        for ln in eachline(f)
            ln=replace(ln,"\n","")
            sidx=search(ln,"=").start-1
            eidx=search(ln,"=").stop+1
            if sidx>0
                k=strip(ln[1:sidx])
                v=strip(ln[eidx:end])
                if (length(k)>0) & (length(v)>0)
                    println("LINE:: ~"*k*"~"*v*"~")
                    #println("DEF::       ______ ",d[Symbol(k)],"  ----> ",typeof(d[Symbol(k)]) )
                    newk=Symbol(k)
                    dv=get(d,newk,nothing)
                    if dv!=nothing
                        newv = cvt(v,typeof(dv))
                        #println("IN::        ______ ",Symbol(k),"  ----> ",ct, " ....... ", typeof(ct))
                        #println( d[Symbol(k)],"|",v," {",typeof(d[Symbol(k)]),"}  ----->  ",newk,"|",newv,"  {",typeof(newv),"}")
                        d[newk] = newv 
                    else
                        println(newk," in app.cfg, but not defaults")
                        #d[newk] = v 
                    end
                else
                    #println("CFG ERROR :: invalid length : "*"~"*k*"~"*v*"~",typeof(d[Symbol(k)]))
                end
            else
                if length(ln)>0   println("CFG ERROR :: invalid cfg parameter : ~",ln,"~") end
            end
        end
        close(f)
    else
        println("Using Default configuration")
    end
    return d
end



# Retype DataFrame Columns
# Ensure/Convert DataFrame Column Types
function retypeBools!(dfd::DataFrame,vars::Array{Symbol}=Symbol[])
    if length(vars)==0 vars=names(dfd) end
    for c in vars
        if typeof(dfd[c]) in [DataArray{Any,1}, DataArray{String,1}]  #DataArray{Bool,1}
            println(c," ~ ",typeof(dfd[c]))
            if  typeof(dfd[c]) == DataArray{Any,1}
                if unique(map(x-> typeof(x),dfd[c])) == [Bool]
                    println("Converting Col")
                    dfd[c] = convert(Array{Bool}, dfd[c])
                else
                    println("Type = Any; but not Bool")
                end
            else
                if true in unique(isna(dfd[c]))
                    println("String Col has NA's")
                else
                    ul=convert(Array{String},unique(map(x -> lowercase(x), dfd[c] )))
                    if (length(setdiff(ul,["true","false"]))==0) & (length(intersect(ul,["true","false"]))>0)
                        println("Mapping String Col")
                        dfd[c]=  map(x-> lowercase(x)=="true" ? true :  false , dfd[c])
                        dfd[c] = convert(Array{Bool}, dfd[c])
                    else
                        println("String col not Bool")
                    end
                end
            end
        end        
    end
end

