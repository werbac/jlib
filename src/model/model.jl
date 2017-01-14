
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays


# ************ CUSTOMIZE *********************
#root="/mnt/resource/analytics/models/Jenny-o"
#root="/mnt/resource/analytics/models/NatChoice"
#root="/mnt/resource/analytics/models/rev"
#root="/mnt/resource/analytics/models/ALL#10"
#root="/mnt/resource/analytics/models/ALL#11"
#root="/mnt/resource/analytics/models/CDW#6"
#root="/mnt/resource/analytics/models/HormelChili#8"
#root="/mnt/resource/analytics/models/Rev#8"                            
#root="/mnt/resource/analytics/models/Jenny-o"
#root="/mnt/resource/analytics/models/hunts"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/rev"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/NatChoice"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/ALL10"

# ************ END CUSTOMIZE *********************
root = !isdefined(:root) ? pwd() : root

        
modelsDict = read_modelsDict(root) 
                


function runGlm(root::String, modelname::Symbol, ranef::Symbol=:empty)
    println("Loding data : ")
    dfd = read_dfd(root) 
    modelsDict = read_modelsDict(root)
    poolit!(dfd,modelsDict[:factors])
    m=modelsDict[modelname]
    return runGlm(dfd, m, ranef)
end

function runGlm(dfd::DataFrame,m::Dict, ranef::Symbol=:empty)
    if ranef==:empty
        f=genF(m[:y_var],m[:finalvars])
        cols = convert(Array{Symbol},vcat(m[:finalvars], [m[:y_var], m[:logvar]]))
        println("Ok - Running ",m[:modelName]," -  Glm on proc ",myid()," with : ",f)
    else
        f=genF(m[:y_var],setdiff(vcat(m[:finalvars],ranef),[:group]))
        dfd[ranef] = relevel(dfd[ranef], "none")
        cols = convert(Array{Symbol},vcat(m[:finalvars], [m[:y_var], m[:logvar]], ranef ))
        println("Ok - Running ",m[:modelName]," - ",ranef," Glm on proc ",myid()," with : ",f)
    end
    if m[:Buyer_Pos_P1_is1]
        return glm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk])
    else
        return return glm(f, dfd[cols] , m[:dist], m[:lnk])
    end
end        
        
                

function runGlmm(root::String, modelname::Symbol, ranef::Symbol=:empty, runθ::Bool=true)
    println("Loding data : Glmm : ",modelname," ~~~ ",ranef)
    dfd = read_dfd(root)   #dfd = readtable(mod_fname,header=true);
    modelsDict = read_modelsDict(root) #modelsDict = readModels(jmod_fname)
    poolit!(dfd,modelsDict[:factors])
    m=modelsDict[modelname]
    if ranef==:empty
        v_out=Dict()
        for r in m[:raneff]
            v_out[r] = runGlmm(dfd, m, [r], runθ)
        end
        return v_out
    else
        return runGlmm(dfd, m, [ranef], runθ)
    end
end

function runGlmm(dfd::DataFrame,m::Dict, ranef::Array{Symbol}=[:empty], runθ::Bool=true)
    ranef = ranef==[:empty] ? m[:raneff]  : ranef
    for r in ranef
        dfd[r] = relevel(dfd[r], "none")
    end
    f=genF(m[:y_var],setdiff(m[:finalvars],[:group]),ranef)
    cols = convert(Array{Symbol},vcat(m[:finalvars], [m[:y_var], m[:logvar]], ranef ))
    println("Ok - Running ",m[:modelName]," - ",ranef," Glmm on proc ",myid()," with : ",f)
    if m[:Buyer_Pos_P1_is1]
        #gmm1 = fit!(glmm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk]), false, 0)   #gmm1 = fit!(glmm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk]), nAGQ=0)
        gmm1 = fit!(glmm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk]), verbose=false, fast=true)
        if runθ if sum(raneffect(gmm1)[:coef]) == 0 profileθ(gmm1) end  end
        return gmm1
    else
        #gmm1 = fit!(glmm(f, dfd[cols], m[:dist], m[:lnk]), false, 0)   #gmm1 = fit!(glmm(f, dfd[cols], m[:dist], m[:lnk]), nAGQ=0)
        gmm1 = fit!(glmm(f, dfd[cols] , m[:dist], m[:lnk]), verbose=false, fast=true)
        if runθ if sum(raneffect(gmm1)[:coef]) == 0 profileθ(gmm1) end end
        return gmm1
    end
end

  


# *******************************************
function runModels(root::String, modelsDict::Dict)
    dfd = read_dfd(root) #dfd = readtable(mod_fname,header=true);
    poolit!(dfd,modelsDict[:factors])
    res = OrderedDict()
    q = Symbol[]
    cnt=0
    for (k,v) in OrderedDict(m=>modelsDict[m][:raneff] for m in [:ipen,:idolocc,:iocc]) 
         for r in v push!(q,k); push!(q,r); cnt=cnt+1 end 
    end
    x=reshape(q, (2,cnt))
    ml = permutedims(x, [2, 1])
    m=:empty
    for i in 1:length(ml[:,1]) 
        if m != ml[i,:][1]
            gr=Symbol("glm_"*string(ml[i,:][1]))
            res[gr] = runGlm(dfd,modelsDict[ml[i,:][1]]) #res[gr] = runGlm(mod_fname, jmod_fname, ml[i,:][1])
            m=ml[i,:][1]
        end
        r=Symbol("glmm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2]))
        res[r] = runGlmm(dfd,modelsDict[m], [ml[i,:][2]], true)   #res[r] = runGlmm(mod_fname, jmod_fname, ml[i,:][1], ml[i,:][2]) 
        covg=Symbol("covglm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2]))
        res[covg] = runGlm(dfd,modelsDict[m], ml[i,:][2])   #res[covg] = runGlm(mod_fname, jmod_fname, ml[i,:][1], ml[i,:][2]) 
    end 
    return res                        
end        
#dx = runModels(root,modelsDict) 
dx = runModels(root,modelsDict) 
  
"""
function consolidateResultsOLD(modelsDict::Dict,dx::OrderedDict)  
  sdf = DataFrame(parameter=String[], coef=Int64[], stderr=Int64[], zval=Int64[], pval=Int64[], model=Symbol[], ranef=Symbol[], modelType=Symbol[])
    for m in [:iocc, :idolocc, :ipen]
        g = dx[Symbol("glm_"*string(m))]
        dfx = coefDF(g,false)
        dfx[:model] = Symbol(replace(string(m),"i",""))
        dfx[:ranef] = :none
        dfx[:modelType] = :GLM
        sdf = vcat(sdf,dfx)
        grp = sdf[(sdf[:modelType].==:GLM)&(sdf[:model].==Symbol(replace(string(m),"i","")))&(sdf[:parameter].=="group"),:coef][1]
        err = sdf[(sdf[:modelType].==:GLM)&(sdf[:model].==Symbol(replace(string(m),"i","")))&(sdf[:parameter].=="group"),:stderr][1]
        for r in  modelsDict[m][:raneff]
            k="_"*string(m)*"_"*string(r)
            g = dx[Symbol("glmm"*k)]
            dfx = raneffect(g)     
            dfx[:zval]=0.0
            dfx[:pval]=0.0
            dfx[:model] = Symbol(replace(string(m),"i",""))
            dfx[:ranef]=r
            dfx[:modelType]=:GLMM
            sdf = vcat(sdf,dfx)
            g = dx[Symbol("covglm"*k)]
            dfx = coefDF(g,false)
            dfx[:model] = Symbol(replace(string(m),"i",""))
            dfx[:ranef] = r
            dfx[:modelType] = :RanCov  
            dfx = dfx[findin(dfx[:parameter],filter(x->contains(x,string(r)) ,Array(dfx[:parameter]))),:]
            #dfx[:coef]=dfx[:coef]-grp  
            #dfx[:stderr]=sqrt(dfx[:stderr].^2+err^2)
            dfx[:parameter] = map(x-> replace(x,string(r)*": ",""),dfx[:parameter])
            sdf = vcat(sdf,dfx)                    
        end
    end
    return sdf
end     
dx = runModels(root,modelsDict) 
   
    
function consolidateResults(modelsDict::Dict,dx::OrderedDict)  
  sdf = DataFrame(parameter=String[], coef=Int64[], stderr=Int64[], zval=Int64[], pval=Int64[], model=Symbol[], ranef=String[], modelType=String[])
    for m in [:iocc, :idolocc, :ipen]
        g = dx[Symbol("glm_"*string(m))]
        dfx = coefDF(g,false)
        dfx[:model] = replace(string(m),"i","")
        dfx[:ranef] = "none"
        dfx[:modelType] = "GLM"
        sdf = vcat(sdf,dfx)
        grp = sdf[(sdf[:modelType].=="GLM")&(sdf[:model].==replace(string(m),"i",""))&(sdf[:parameter].=="group"),:coef][1]
        err = sdf[(sdf[:modelType].=="GLM")&(sdf[:model].==replace(string(m),"i",""))&(sdf[:parameter].=="group"),:stderr][1]
        for r in  modelsDict[m][:raneff]
            k="_"*string(m)*"_"*string(r)
            g = dx[Symbol("glmm"*k)]
            dfx = raneffect(g)     
            dfx[:zval]=0.0
            dfx[:pval]=0.0
            dfx[:model] = replace(string(m),"i","")
            dfx[:ranef]=string(r)
            dfx[:modelType]="GLMM"
            sdf = vcat(sdf,dfx)
            g = dx[Symbol("covglm"*k)]
            dfx = coefDF(g,false)
            dfx[:model] = replace(string(m),"i","")
                dfx[:ranef] = string(r)
            dfx[:modelType] = "RanCov"  
            dfx = dfx[findin(dfx[:parameter],filter(x->contains(x,string(r)) ,Array(dfx[:parameter]))),:]
            #dfx[:coef]=dfx[:coef]-grp  
            #dfx[:stderr]=sqrt(dfx[:stderr].^2+err^2)
            dfx[:parameter] = map(x-> replace(x,string(r)*": ",""),dfx[:parameter])
            sdf = vcat(sdf,dfx)                    
        end
    end
    return sdf
end         
"""
    
  
function consolidateResults(modelsDict::Dict,dx::OrderedDict=OrderedDict())  
  sdf = DataFrame(parameter=String[], coef=Int64[], stderr=Int64[], zval=Int64[], pval=Int64[], model=String[], ranef=String[], modelType=String[])
    for m in [:iocc, :idolocc, :ipen]
        g = length(dx) > 0 ? dx[Symbol("glm_"*string(m))] : getw( Symbol("glm_"*string(m)) )  
        #g = takew(  Symbol("glm_"*string(m))    )
        dfx = coefDF(g,false)
        dfx[:model] = replace(string(m),"i","")
        dfx[:ranef] = "none"
        dfx[:modelType] = "GLM"
        sdf = vcat(sdf,dfx)
        grp = sdf[(sdf[:modelType].=="GLM")&(sdf[:model].==replace(string(m),"i",""))&(sdf[:parameter].=="group"),:coef][1]
        err = sdf[(sdf[:modelType].=="GLM")&(sdf[:model].==replace(string(m),"i",""))&(sdf[:parameter].=="group"),:stderr][1]
        for r in  modelsDict[m][:raneff]
            k="_"*string(m)*"_"*string(r)
            g = length(dx) > 0 ?  dx[Symbol("glmm"*k)]  : getw(Symbol("glmm"*k))
            dfx = raneffect(g)     
            dfx[:zval]=0.0
            dfx[:pval]=0.0
            dfx[:model] = replace(string(m),"i","")
            dfx[:ranef]=string(r)
            dfx[:modelType]="GLMM"
            sdf = vcat(sdf,dfx)
            g = length(dx) > 0 ?  dx[Symbol("covglm"*k)]    : getw(Symbol("covglm"*k))
            dfx = coefDF(g,false)
            dfx[:model] = replace(string(m),"i","")
                dfx[:ranef] = string(r)
            dfx[:modelType] = "RanCov"
            dfx = dfx[findin(dfx[:parameter],filter(x->contains(x,string(r)) ,Array(dfx[:parameter]))),:]
            #dfx[:coef]=dfx[:coef]-grp  
            #dfx[:stderr]=sqrt(dfx[:stderr].^2+err^2)
            dfx[:parameter] = map(x-> replace(x,string(r)*": ",""),dfx[:parameter])
            sdf = vcat(sdf,dfx)                    
        end
    end
    return sdf
end     
#dfx = consolidateResults(modelsDict)  

#dx = runModels(root,modelsDict) 
#dfx = consolidateResults(modelsDict, dx)      
#save_dfx(root,dfx)   #writetable(root*"/campaign.csv", dfx)  
dfx = consolidateResults(modelsDict, dx)      
save_dfx(root,dfx)   #writetable(root*"/campaign.csv", dfx)  
    
    
    
     
    
        