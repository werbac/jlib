
function runModels( root::String, modelsDict::Dict, dfd::DataFrame=NA)
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
           isna(dfd) ? ( @swarm gr runGlm(root, ml[i,:][1]) ) : ( res[gr] = runGlm(dfd,modelsDict[ml[i,:][1]]) ) 
           m=ml[i,:][1]
       end
       r=Symbol("glmm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2]))
       isna(dfd) ? (@swarm r runGlmm(root, ml[i,:][1], ml[i,:][2]) ) : (res[r] = runGlmm(dfd,modelsDict[m], [ml[i,:][2]], true) )
       covg=Symbol("covglm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2])   )
       isna(dfd) ? ( @swarm covg runGlm(root, ml[i,:][1], ml[i,:][2]) ) : ( res[covg] = runGlm(dfd,modelsDict[m], ml[i,:][2])  )
   end 
   if isna(dfd)
       sleep(20)
       while !iscompletew() println("No!"); sleep(5); end 
       println("DONE")
   end
   return res
end 

    
function runModels(root::String, modelsDict::Dict, )
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
        

function genDHHMeans(dfx::DataFrame)
    println("TOTAL DOLHH")
    o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    adj_mean_score0=o*y*p 
    println(".... 0 : $adj_mean_score0 ($o, $y,$p) ...")
    o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    adj_mean_score1=o*y*p
    println(".... 1 : $adj_mean_score1 ($o, $y,$p) ...")
    unadj_avg_cntrl_hh_pre = mean(dfd[ (dfd[:group] .== 0), :prd_1_net_pr_pre] )
    unadj_avg_expsd_hh_pre = mean(dfd[ (dfd[:group] .== 1), :prd_1_net_pr_pre] )
    unadj_avg_cntrl_hh_pst = mean(dfd[ (dfd[:group] .== 0), :prd_1_net_pr_pos] )
    unadj_avg_expsd_hh_pst = mean(dfd[ (dfd[:group] .== 1), :prd_1_net_pr_pos] )
    pre=["group",0.0,0.0,0.0,0.0,"dolhh","ranef","GLM"]
    vars=[0.0,0.0,adj_mean_score0,adj_mean_score1,unadj_avg_expsd_hh_pre,unadj_avg_expsd_hh_pst,unadj_avg_cntrl_hh_pre,unadj_avg_cntrl_hh_pst]
    cis = zeros(10)  # cis and pvals [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    pvals=[0.0,0.0]
    cnts=convert(Array{Int64},collect(values(getCnts(dfd))))
    push!(dfx, filler(dfx,vcat(pre,vars,cis,cnts)))
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
#            println(ranef," ~~ ",level)
            o=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            y=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            p=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            adj_mean_score0=o*y*p      
            println(".... 1 : $adj_mean_score1 ($o, $y,$p) ...")
            o=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            y=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            p=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            adj_mean_score1=o*y*p
            println(".... 0 : ",adj_mean_score0," ... 1 :",adj_mean_score1)
            unadj_avg_expsd_hh_pre = mean(dfd[ (dfd[:group].==1) & (dfd[Symbol(ranef)].==level), :prd_1_net_pr_pre] )
            unadj_avg_expsd_hh_pst = mean(dfd[ (dfd[:group].==1) & (dfd[Symbol(ranef)].==level), :prd_1_net_pr_pos] )
            unadj_avg_cntrl_hh_pre = mean(dfd[ (dfd[:group].==0) & (dfd[Symbol(ranef)].==level), :prd_1_net_pr_pre] )
            unadj_avg_cntrl_hh_pst = mean(dfd[ (dfd[:group].==0) & (dfd[Symbol(ranef)].==level), :prd_1_net_pr_pos] )         
            exposed=true # temp hack - all exposed
            if exposed
                unadj_avg_cntrl_hh_pre = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolhh") ,:unadj_avg_cntrl_hh_pre][1]
                unadj_avg_cntrl_hh_pst = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolhh") ,:unadj_avg_cntrl_hh_pst][1]     
            end   
            pre=[level,0.0,0.0,0.0,0.0,"dolhh",ranef,"GLMM"]
            vars=[0.0,0.0,adj_mean_score0,adj_mean_score1,unadj_avg_expsd_hh_pre,unadj_avg_expsd_hh_pst,unadj_avg_cntrl_hh_pre,unadj_avg_cntrl_hh_pst ]
            cis=zeros(10) #[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            cnts=convert(Array{Int64},collect(values(getCnts(dfd,ranef,level))))
            push!(dfx, filler(dfx,vcat(pre,vars,cis,cnts)))
        end
    end
end
genDHHMeans(dfx)


        


ZDict = Dict("onetail_80_pct_intrvl" => 0.84 ,"onetail_90_pct_intrvl" => 1.28, "twotail_80_pct_intrvl" => 1.28, "twotail_90_pct_intrvl" => 1.65)



v_ttl=2

function CIs_O_LB(iDict::OrderedDict, zscore::Float64, iAccuracy::Float64=0.000000001)
    M = get(iDict, :M, NA)
    Mt = get(iDict, :Mt, NA)
    Mc = get(iDict, :Mc, NA)
    N = get(iDict, :N, NA)
    Nt = get(iDict, :Mt, NA)
    Nc = get(iDict, :Nc, NA)
    B1 = get(iDict, :B1, NA)
    B2 = get(iDict, :B2, NA)
    B3 = get(iDict, :B3, NA)
    SE1 = get(iDict, :SE1, NA)
    SE2 = get(iDict, :SE2, NA)
    SE3 = get(iDict, :SE3, NA)
    o_SE0 = get(iDict, :o_SE0, 0)
    y_SE0 = get(iDict, :y_SE0, 0)
    p_SE0 = get(iDict, :p_SE0, 0)
SEsq=sqrt(SE1^2+SE2^2+SE3^2+o_SE0^2+y_SE0^2+p_SE0^2)
    o_B0 = get(iDict, :o_B0, 0)
    y_B0 = get(iDict, :y_B0, 0)
    p_B0 = get(iDict, :p_B0, 0)

    
    o_mean_score0 = get(iDict, :o_mean_score0, NA)
    o_mean_score1 = get(iDict, :o_mean_score1, NA)
    y_mean_score0 = get(iDict, :y_mean_score0, NA)
    y_mean_score1 = get(iDict, :y_mean_score1, NA)
    p_mean_score0 = get(iDict, :p_mean_score0, NA)
    p_mean_score1 =get(iDict, :p_mean_score1, NA)
Bsum=(B1+B2+B3)-(o_B0+y_B0+p_B0)
    
    ztot = Bsum-(zscore*SEsq)
    ######CONFIDENCE INTERVAL - LB ########        
    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))
    @variable(m, Bocc <= (B1-o_B0))
    @variable(m, Bdolocc <= (B2-y_B0))
    @variable(m, Bpen <= (B3-p_B0))
    @NLobjective(m, Min, ((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))* 
                                  ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))*
                                  ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
                                  )
                                  -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))*
                                    ((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))*
                                    ((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
                                    )
                                  )
                               )
    @constraint(m, (0.000000000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAccuracy)) 
    status = solve(m)
    mval=getobjectivevalue(m)
    println(mval)
    mval_out=mval/(o_mean_score0*y_mean_score0*p_mean_score0)
    return status == :Optimal ? mval_out : -Inf 
end
CIs_O_LB(od, 0.84, 0.000000001)


onetail_80 FAILED : DataFrameRow (row 112)
modelType                 GLMM
model                     dolhh
ranef                     ad_type
parameter                 Desktop
adj_dod_effct             10.02936692284249
*****onetail_80_pct_intrvl_lb  10.955272420022512
onetail_80_pct_intrvl_ub  20.010743178509298


LB onetail_80 (ad_type~Desktop):=  - 0.000000001, DataStructures.
od = OrderedDict(:M=>108092,:Mt=>10810,:Mc=>106004,:N=>594930,:Nc=>584120,:B1=>-0.000998217,:B2=>0.0583124,:B3=>0.066694,:SE1=>0.0201122,:SE2=>0.0193296,:SE3=>0.0301203,:o_SE0=>0,:y_SE0=>0,:p_SE0=>0,:o_B0=>0,:y_B0=>0,:p_B0=>0,:o_mean_score0=>1.24344,:o_mean_score1=>1.2422,:y_mean_score0=>3.3212,:y_mean_score1=>3.52063,:p_mean_score0=>0.181431,:p_mean_score1=>0.193944,:metakey=>"ad_type~Desktop",:Bsum=>0.124008,:onetail_pval=>0.998738,:twotail_pval=>0.997476,:twotail_90_pct_intrvl_lb=>0.0678382,:twotail_90_pct_intrvl_ub=>0.245724,:twotail_80_pct_intrvl_lb=>0.0867394,:twotail_80_pct_intrvl_ub=>0.224721)
0.08208352266101582
            : Confidence Interval : 0.10955272420022512


mDict=collectModels(dfx,"GLMM","ad_type","Desktop")
calcPValue_Opt(mDict)
dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100


dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].=="ad_type")&(dfx[:parameter].=="Desktop"),:]


            adjctl = md[:o_mean_score0] * md[:y_mean_score0] * md[:p_mean_score0]
            adjexp = md[:o_mean_score1] * md[:y_mean_score1] * md[:p_mean_score1]
            println(k," ~~ ",adjctl,"~",adjexp) #,md)            
            sdf[sdf[:key].==k,:adj_mean_cntrl_grp] = adjctl
            sdf[sdf[:key].==k,:adj_mean_expsd_grp] = adjexp       
            sdf[sdf[:key].==k,:adj_dod_effct] = ((adjexp - adjctl) ./ adjctl ) * 100



function collectModels(dfx::DataFrame,modelType::String,ranef::String="",level::String="")
    mDict = OrderedDict()
    if modelType=="GLM"
        o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:]
        y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:]
        p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:]
    else
        o=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        y=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        p=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
    end
    mDict[:M]=p[:M][1]
    mDict[:Mt]=p[:Mt][1]
    mDict[:Mc]=p[:Mc][1]
    mDict[:N]=p[:N][1]
    mDict[:Mt]=p[:Nt][1]
    mDict[:Nc]=p[:Nc][1]
    if modelType=="GLM"
        mDict[:B1]=o[:coef][1]   
        mDict[:B2]=y[:coef][1]
        mDict[:B3]=p[:coef][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
        mDict[:SE1]=o[:stderr][1]
        mDict[:SE2]=y[:stderr][1]
        mDict[:SE3]=p[:stderr][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100S
    else
        mDict[:B1]=o[:adj_coef][1]   
        mDict[:B2]=y[:adj_coef][1]
        mDict[:B3]=p[:adj_coef][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
        mDict[:SE1]=o[:adj_stderr][1]
        mDict[:SE2]=y[:adj_stderr][1]
        mDict[:SE3]=p[:adj_stderr][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100S
    end
    mDict[:o_SE0]=0
    mDict[:y_SE0]=0
    mDict[:p_SE0]=0
    mDict[:o_B0]=0
    mDict[:y_B0]=0
    mDict[:p_B0]=0
    mDict[:o_mean_score0]=o[:adj_mean_score0][1]
    mDict[:o_mean_score1]=o[:adj_mean_score1][1]
    mDict[:y_mean_score0]=y[:adj_mean_score0][1]
    mDict[:y_mean_score1]=y[:adj_mean_score1][1]
    mDict[:p_mean_score0]=p[:adj_mean_score0][1]
    mDict[:p_mean_score1]=p[:adj_mean_score1][1] 
    return mDict
end





function ConfidenceIntervals(dfx::DataFrame)
   # mDict=collectModels(dfx,"GLM")
#    println("Running CI for : Total")
#    mDict[:metakey] = "Total Campaign"        
#    calcPValue_Opt(mDict)
#    CIs_O(mDict)
#    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_pval] = mDict[:onetail_pval]*100
#    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_pval] = mDict[:twotail_pval]*100
#    for k in keys(ZDict)
#        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
#        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
#    end
    for ranef in ["ad_type"] #unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].=="ad_type")&(dfx[:parameter].=="Desktop"),:ranef])
        for level in ["Desktop"] #unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")
            mDict[:metakey] = ranef*"~"*level        
            calcPValue_Opt(mDict)
            CIs_O(mDict)
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_pval] = mDict[:onetail_pval]*100
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_pval] = mDict[:twotail_pval]*100
            for k in keys(ZDict)
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
            end
        end
    end
end
ConfidenceIntervals(dfx)


        
        
# ---------------------------------------------------------

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


isnumeric{T<:Number}(::AbstractArray{T}) = true
isnumeric(::Any) = false
isfloat{T<:Float64}(::AbstractArray{T}) = true
isfloat(::Any) = false
isint{T<:Int64}(::AbstractArray{T}) = true
isint(::Any) = false
isString{T<:String}(::AbstractArray{T}) = true
isString(::Any) = false
isBool{T<:Bool}(::AbstractArray{T}) = true
isBool(::Any) = false
isSymbol{T<:Symbol}(::AbstractArray{T}) = true
isSymbol(::Any) = false


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
        if n ∈ cols && isfloat(v)
            push!(a, 0.0)
        elseif n ∈ cols && isint(v)
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













# ***************************************************************************************************************************
#export JULIA_NUM_THREADS=4
ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg" 
for i in ["10.106.128.212","10.106.128.213","10.106.128.214","10.106.128.215","10.106.128.216","10.106.128.217",
          "10.106.128.218","10.106.128.219","10.106.128.220","10.106.128.221"]
             if i != string(getipaddr()) println(i); addprocs([(i, 3)]) end
end #addprocs(1) 
@everywhere insert!(Base.LOAD_CACHE_PATH, 1, "/mapr/mapr04p/analytics0001/analytic_users/jpkg/lib/v0.5")
@everywhere pop!(Base.LOAD_CACHE_PATH)
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
@everywhere using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve

root = !isdefined(:root) ? pwd() : root
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


#jmod_fname = root*"/dfd_model.json"
#mod_fname = root*"/dfd_model.csv"         
modelsDict = read_modelsDict(root) #modelsDict = readModels(jmod_fname) 
                


@everywhere function runGlm(root::String, modelname::Symbol, ranef::Symbol=:empty)
    println("Loding data : ")
    dfd = read_dfd(root) #dfd = readtable(mod_fname,header=true);
    modelsDict = read_modelsDict(root)#modelsDict = readModels(jmod_fname)
    poolit!(dfd,modelsDict[:factors])
    m=modelsDict[modelname]
    return runGlm(dfd, m, ranef)
end

@everywhere function runGlm(dfd::DataFrame,m::Dict, ranef::Symbol=:empty)
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
        
                

@everywhere function runGlmm(root::String, modelname::Symbol, ranef::Symbol=:empty, runθ::Bool=true)
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

@everywhere function runGlmm(dfd::DataFrame,m::Dict, ranef::Array{Symbol}=[:empty], runθ::Bool=true)
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






"""
SWARMING - hold on
"""
const SHELF = OrderedDict()
freew() = setdiff(workers(),[v[:worker] for (k,v) in SHELF])             
hasfreew() = length(freew()) > 0 ? true : false
nextFreew() = length(freew()) > 0 ? freew()[1]   : NA
#statusW() =  [ string(k)*"   status: "*string(v[:status]) for (k,v) in values(SHELF)] 
#statusw() =  for (k,v) in SHELF println(string(k)*"   status: "*string(v[:status])*"   "*string(v[:worker]) ) end
notcompletew() = filter((k,v)-> v[:status]!=:complete ,SHELF) 
iscompletew() = length(notcompletew()) == 0 ? true : false
keysw() = collect(keys(SHELF))
getw(k::Symbol) = SHELF[k][:results]
function takew(k::Symbol) 
    vo = deepcopy(SHELF[k][:results]) 
    delete!(SHELF, k) 
    return vo
end


type statusW function statusW()  return new() end  end 
function Base.show(io::IO, statusw::statusW)  for (k,v) in SHELF println(string(k)*"   status: "*string(v[:status])*"   "*string(v[:worker]) ) end end
statusw=statusW()              
                
function swarm(k::Symbol, thunk)
    #println("SWARMING  : ",k)
    #println("SWARM DUMP : ",typeof(thunk)==Expr)
    SHELF[k] = Dict{Symbol,Any}(:worker =>0,:expr => thunk, :status =>:ready, :channel => Channel(1))
    t = SHELF[k]
    while !hasfreew() sleep(15) end   # println("waiting for free worker");
    t[:worker] = nextFreew()
    wid = t[:worker]        
    if isna(wid)
        return false
    else
        #remotecall_fetch(runGlmm, wid, mod_fname, jmod_fname, w[:model], w[:renef] )
        t[:status] = :running
        #splice!(expr.args, 1:1,[:remotecall_fetch,expr.args[1],wid])
        #println("Starting model")  
        println("SWARMING  : ",k, " on ", wid)
        @async put!(t[:channel], fetch(Base.sync_add(remotecall(thunk, wid)))  )
        while isready(t[:channel])==false sleep(15) end  #println("NOT Ready : ",t[:worker] );      
        t[:results] = take!(t[:channel])    
        t[:status] = :complete
        t[:worker] = 0
        println("SWARM Task Completed : ",k)
        close(t[:channel]) 
        return true
    end     
end
                
                
macro swarm(p, expr)
    expr = Base.localize_vars(esc(:(()->($expr))), false)
    ex = :( @async swarm($(esc(p)), $expr))
end   
        
# ******************************************  
           
                           
function runModels( root::String, modelsDict::Dict)
   @swarm :glm_ipen runGlm(root, :ipen)
   sleep(25)
   @swarm :glm_iocc runGlm(root, :iocc)
   sleep(25)
   @swarm :glm_idolocc runGlm(root, :idolocc)    
   q = Symbol[]
   cnt=0
   for (k,v) in OrderedDict(m=>modelsDict[m][:raneff] for m in [:ipen,:idolocc,:iocc]) 
        for r in v push!(q,k); push!(q,r); cnt=cnt+1 end 
   end
   x=reshape(q, (2,cnt))
   ml = permutedims(x, [2, 1])
   m=:empty
   for i in 1:length(ml[:,1]) 
       r=Symbol("glmm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2]))
       @swarm r runGlmm(root, ml[i,:][1], ml[i,:][2]) 
       covg=Symbol("covglm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2])   )
       @swarm covg runGlm(root, ml[i,:][1], ml[i,:][2]) 
   end 
   sleep(20)
   while !iscompletew() println("Not Complete yet!"); sleep(5); end 
   println("DONE!DONE!DONE!")
end        
runModels(root,modelsDict)                   
     
             
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
dfx = consolidateResults(modelsDict)  
save_dfx(root,dfx)      #writetable(root*"/campaign.csv", dfx)                      

    
    
# *************************************************************************************************************************** 
# *************************************************************************************************************************** 
# *************************************************************************************************************************** 
# ***************************************************************************************************************************            
# ***************************************************************************************************************************           
            
          
#mv(root*"/dfx.csv",root*"/dfx.csv_prev")
#writetable(root*"/dfx_ready.csv", dfx) 
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/rev"
#modelsDict = read_modelsDict(root)  
#modelsDict[:occ]=modelsDict[:iocc]; delete!(modelsDict,:iocc)
#modelsDict[:dolocc]=modelsDict[:idolocc]; delete!(modelsDict,:idolocc)
#modelsDict[:pen]=modelsDict[:ipen]; delete!(modelsDict,:ipen)
##dfx = read_dfx(root)   
            
           
#export JULIA_NUM_THREADS=4
ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg" 
for i in ["10.106.128.212","10.106.128.213","10.106.128.214","10.106.128.215","10.106.128.216","10.106.128.217",
          "10.106.128.218","10.106.128.219","10.106.128.220","10.106.128.221"]
             if i != string(getipaddr()) println(i); addprocs([(i, 3)]) end
end #addprocs(1) 
@everywhere insert!(Base.LOAD_CACHE_PATH, 1, "/mapr/mapr04p/analytics0001/analytic_users/jpkg/lib/v0.5")
@everywhere pop!(Base.LOAD_CACHE_PATH)
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
@everywhere using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve


"""
SWARMING - hold on
"""
const SHELF = OrderedDict()
freew() = setdiff(workers(),[v[:worker] for (k,v) in SHELF])             
hasfreew() = length(freew()) > 0 ? true : false
nextFreew() = length(freew()) > 0 ? freew()[1]   : NA
#statusW() =  [ string(k)*"   status: "*string(v[:status]) for (k,v) in values(SHELF)] 
#statusw() =  for (k,v) in SHELF println(string(k)*"   status: "*string(v[:status])*"   "*string(v[:worker]) ) end
notcompletew() = filter((k,v)-> v[:status]!=:complete ,SHELF) 
iscompletew() = length(notcompletew()) == 0 ? true : false
keysw() = collect(keys(SHELF))
getw(k::Symbol) = SHELF[k][:results]
function takew(k::Symbol) 
    vo = deepcopy(SHELF[k][:results]) 
    delete!(SHELF, k) 
    return vo
end
            # NOTE : create function rmW() # remove value from SHELF
            #                        clearW() #remove all values from SHELF
            #                        waitW()  # while !iscompletew() println("Not Complete yet!"); sleep(5); end 
            
type statusW function statusW()  return new() end  end 
function Base.show(io::IO, statusw::statusW)  for (k,v) in SHELF println(string(k)*"   status: "*string(v[:status])*"   "*string(v[:worker]) ) end end
statusw=statusW()              
                
function swarm(k::Symbol, thunk)
    #println("SWARMING  : ",k)
    #println("SWARM DUMP : ",typeof(thunk)==Expr)
    SHELF[k] = Dict{Symbol,Any}(:worker =>0,:expr => thunk, :status =>:ready, :channel => Channel(1))
    t = SHELF[k]
    while !hasfreew() sleep(15) end   # println("waiting for free worker");
    t[:worker] = nextFreew()
    wid = t[:worker]        
    if isna(wid)
        return false
    else
        #remotecall_fetch(runGlmm, wid, mod_fname, jmod_fname, w[:model], w[:renef] )
        t[:status] = :running
        #splice!(expr.args, 1:1,[:remotecall_fetch,expr.args[1],wid])
        #println("Starting model")  
        println("SWARMING  : ",k)
        @async put!(t[:channel], fetch(Base.sync_add(remotecall(thunk, wid)))  )
        while isready(t[:channel])==false sleep(15) end  #println("NOT Ready : ",t[:worker] );      
        t[:results] = take!(t[:channel])    
        t[:status] = :complete
        t[:worker] = 0
        println("SWARM Task Completed : ",k)
        close(t[:channel]) 
        return true
    end     
end                
macro swarm(p, expr)
    expr = Base.localize_vars(esc(:(()->($expr))), false)
    ex = :( @async swarm($(esc(p)), $expr))
end   
 
           
            
@everywhere function MP_collectModels(dfx::DataFrame, modelType::String,ranef::String="",level::String="") 
    mDict = OrderedDict()
    if modelType=="GLM"
        o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:]
        y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:]
        p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:]
    else
        o=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        y=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        p=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
    end
    mDict[:M]=p[:M][1]
    mDict[:Mt]=p[:Mt][1]
    mDict[:Mc]=p[:Mc][1]
    mDict[:N]=p[:N][1]
    mDict[:Nt]=p[:Nt][1]
    mDict[:Nc]=p[:Nc][1]
    if modelType=="GLM"
        mDict[:B1]=o[:coef][1]   
        mDict[:B2]=y[:coef][1]
        mDict[:B3]=p[:coef][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
        mDict[:SE1]=o[:stderr][1]
        mDict[:SE2]=y[:stderr][1]
        mDict[:SE3]=p[:stderr][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100S
    else
        mDict[:B1]=o[:adj_coef][1]   
        mDict[:B2]=y[:adj_coef][1]
        mDict[:B3]=p[:adj_coef][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
        mDict[:SE1]=o[:adj_stderr][1]
        mDict[:SE2]=y[:adj_stderr][1]
        mDict[:SE3]=p[:adj_stderr][1] #fx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100S
    end
    mDict[:o_SE0]=0
    mDict[:y_SE0]=0
    mDict[:p_SE0]=0
    mDict[:o_B0]=0
    mDict[:y_B0]=0
    mDict[:p_B0]=0
    mDict[:o_mean_score0]=o[:adj_mean_score0][1]
    mDict[:o_mean_score1]=o[:adj_mean_score1][1]
    mDict[:y_mean_score0]=y[:adj_mean_score0][1]
    mDict[:y_mean_score1]=y[:adj_mean_score1][1]
    mDict[:p_mean_score0]=p[:adj_mean_score0][1]
    mDict[:p_mean_score1]=p[:adj_mean_score1][1]     
    println("Running CI for : $ranef : $level")
    mDict[:metakey] = ranef*"~"*level        
    calcPValue_Opt(mDict)
    CIs_O(mDict)              
    return mDict
end


function MP_ConfidenceIntervals(dfx::DataFrame)           
                @swarm :total MP_collectModels(dfx, "GLM")
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
                k = Symbol(ranef*"~"*level)
                @swarm k  MP_collectModels(dfx, "GLMM",ranef,level)
        end
    end
    sleep(20)
    while !iscompletew() println("Not Complete yet!"); sleep(5); end 
    mDict = takew(:total)
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_pval] = mDict[:onetail_pval]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_pval] = mDict[:twotail_pval]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:pval] = mDict[:pval]
    for k in keys(ZDict)
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
    end
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
                mDict = takew( Symbol(ranef*"~"*level) ) 
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:pval] = mDict[:pval]
                for k in keys(ZDict)
                    dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
                    dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
                end 
        end
    end    
end   
            
MP_ConfidenceIntervals(dfx)          
         
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
# *************************************************************************************************************************** 
# *************************************************************************************************************************** 
# *************************************************************************************************************************** 
# ***************************************************************************************************************************            
# ***************************************************************************************************************************           
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
    
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
            dfx[:parameter] = map(x-> replace(x,string(r)*": ",""),dfx[:parameter])
            sdf = vcat(sdf,dfx)                    
        end
    end
    return sdf
end        
#dx = runModels(root,modelsDict) 
#dfx = consolidateResults(modelsDict, dx)      
#save_dfx(root,dfx)   #writetable(root*"/campaign.csv", dfx)  

#dx = runModels(root,modelsDict) 
#dfx = consolidateResults(modelsDict, dx)      
#save_dfx(root,dfx)   #writetable(root*"/campaign.csv", dfx)  
    
    
    
 # ***************************************************************************************************************************   
    
    
    
    
#function distributeDataset()
#    #jmod_fname = root*"/dfd_model.json"
#    #mod_fname
#    #ppp="/mnt/resource/analytics/models/"; run(`ls $ppp`)
#    ips = ["10.63.36.22","10.63.36.23"]
#    for ip in ips
#            cmd=`scp $root/dfd_model.* $ip:$root/`
#            println(cmd)
#            run(cmd)
#        #run(`ssh $ip mkdir $root`)
#        #run(`scp $root/dfd_model.* $ip:$root/`)
#    end
#end    
#addprocs([("iriadmin@10.63.36.22", 1), ("iriadmin@10.63.36.23", 1)])
#addprocs(2)
        
        
    










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



