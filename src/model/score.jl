ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg"

using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, JuMP, NLopt, NLsolve # JStack,
# ************ CUSTOMIZE *********************
#root="/mnt/resource/analytics/models/rev"
#root="/mnt/resource/analytics/models/ALL#10"
#root="/mnt/resource/analytics/models/hunts"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/rev"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/NatChoice"
#root="/mapr/mapr04p/analytics0001/analytic_users/jmdl/ALL10"
# ************ END CUSTOMIZE *********************
root = !isdefined(:root) ? pwd() : root

modelsDict = read_modelsDict(root)  
modelsDict[:occ]=modelsDict[:iocc]; delete!(modelsDict,:iocc)
modelsDict[:dolocc]=modelsDict[:idolocc]; delete!(modelsDict,:idolocc)
modelsDict[:pen]=modelsDict[:ipen]; delete!(modelsDict,:ipen)
dfx = read_dfx(root) 
dfd = read_dfd(root)  
dfx[:unadj_mean_score0] = 0.0
dfx[:unadj_mean_score1] = 0.0
dfx[:adj_mean_score0 ] = 0.0
dfx[:adj_mean_score1 ] = 0.0
dfx[:unadj_avg_expsd_hh_pre]=0.0
dfx[:unadj_avg_expsd_hh_pst]=0.0
dfx[:unadj_avg_cntrl_hh_pre]=0.0
dfx[:unadj_avg_cntrl_hh_pst]=0.0
dfx[:unadj_avg_cntrl_hh_pre]=0.0
dfx[:unadj_avg_cntrl_hh_pst]=0.0
dfx[:onetail_80_pct_intrvl_lb]=0.0
dfx[:onetail_80_pct_intrvl_ub]=0.0
dfx[:onetail_90_pct_intrvl_lb]=0.0
dfx[:onetail_90_pct_intrvl_ub]=0.0
dfx[:twotail_80_pct_intrvl_lb]=0.0
dfx[:twotail_80_pct_intrvl_ub]=0.0
dfx[:twotail_90_pct_intrvl_lb]=0.0
dfx[:twotail_90_pct_intrvl_ub]=0.0
dfx[:onetail_pval] = 0.0
dfx[:twotail_pval] = 0.0
dfx[:M]=0
dfx[:Mt]=0
dfx[:Mc]=0
dfx[:N]=0
dfx[:Nt]=0
dfx[:Nc]=0
dfx[:cnt_expsd_hh]=0
dfx[:cnt_impressions]=0
dfx[:adj_pval]=0.0
dfx[:twotail_pval_to_campaign]=0.0
dfx[:onetail_pval_to_campaign]=0.0
# ************ EQUATIONS COLS - Uncomment to have them populated *****************
#dfx[:unadj_mean_score0_eq] = ""
#dfx[:unadj_mean_score1_eq] = ""
#dfx[:adj_mean_score0_eq ] = ""
#dfx[:adj_mean_score1_eq ] = ""
#dfx[:unadj_avg_expsd_hh_pre_eq]=""
#dfx[:unadj_avg_expsd_hh_pst_eq]=""
#dfx[:unadj_avg_cntrl_hh_pre_eq]=""
#dfx[:unadj_avg_cntrl_hh_pst_eq]=""
#dfx[:unadj_avg_cntrl_hh_pre_eq]=""
#dfx[:unadj_avg_cntrl_hh_pst_eq]=""

# Anonymous functions
inDFX(c::Symbol) = c in names(dfx)
mlist = [modelsDict[:occ],modelsDict[:dolocc],modelsDict[:pen]]  # otherwise we'd include modelsDict[:factors]...etc

function adjustDFX(dfx::DataFrame)
    dfx[:adj_coef] = 0.0
    dfx[:adj_stderr]=  0.0
    for m in mlist 
        model=string(m[:modelName])
        groupBeta = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:coef][1]
        groupStdErr = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:stderr][1]
        for r in map(x->string(x),m[:raneff])   
            noneBeta = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r),:coef][1]
            noneStdErr = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r),:stderr][1]
            
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef]=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:coef]-noneBeta
            
            m = mean(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef])
            f=groupBeta/m
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef]=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef]*f
            
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_stderr]= sqrt.(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:stderr].^2 + noneStdErr^2 ) 
        end
    end
    dfx[(dfx[:modelType].=="GLMM"),:zval] = dfx[(dfx[:modelType].=="GLMM"),:adj_coef] ./ dfx[(dfx[:modelType].=="GLMM"),:adj_stderr] 
    dfx[(dfx[:modelType].=="GLMM"),:pval] = 2.0 .* ccdf(Normal(), abs(dfx[(dfx[:modelType].=="GLMM"),:zval]))  
    return dfx
end
dfx = adjustDFX(dfx)


function genRandCols(dfx::DataFrame, dfd::DataFrame)   #Creates cols with levels replaced with coef :ranef_occ :ranef_dolocc :ranef_pen 
    for m in mlist
        model=string(m[:modelName])
        modelType="GLMM"
        for r in m[:raneff]
            sr=string(r)
            rs=Symbol(sr*"_"*model)
            println(rs)  
            dfd[rs]=deepcopy(dfd[r])
            for row in eachrow(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:ranef].==sr) ,:])
                dfd[(dfd[rs].==row[:parameter]),rs] = string(row[:coef])
                println(rs," ~~ ",row[:parameter])
            end
            if m[:Buyer_Pos_P1_is1] dfd[dfd[:buyer_pos_p1].==0,rs] = "0.0" end
            dfd[rs] = map(x->parse(Float64,x),dfd[rs])
            dfd[rs] = convert(Array{Float64},dfd[rs]) 
        end
    end
end
genRandCols(dfx, dfd)

dfx = dfx[(dfx[:parameter].*dfx[:modelType].!="noneGLMM"),:]  # Remove NONE - here cause need none when gen Rand Cols - convert Float!!!!

function genFixedEQUASION(dfx::DataFrame, m::Dict, con::String) # conditional equations 
    model=string(m[:modelName])
    modelType="GLM"
    #dfname="dfd"
    adj0 = replace(con,"%%","0")
    adj1 = replace(con,"%%","1")
    println(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,:])
    intercept=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="(Intercept)") ,[:coef]][1][1]
    grp=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="group") ,[:coef]][1][1]
    vout0=string(intercept)*""
    vout1=string(intercept)*""
    results=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,[:parameter,:coef]]
    for row in eachrow( results[findin(results[:parameter],map(x->string(x),    setdiff(m[:finalvars],[:group])  )) ,[:parameter,:coef]] )
        vout0=vout0*"+(dfd[$adj0:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
        vout1=vout1*"+(dfd[$adj1:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
    end
    vout1=vout1*"+"*string(grp)
    return vout0, vout1
end


function genFixedCols(dfx::DataFrame,dfd::DataFrame,modelsDict::Dict)   # Creates :occ0 :occ1 :dolocc0 :dolocc1 :pen0 :pen1
    for m in mlist
        r0,r1=genFixedEQUASION(dfx,m,"")  #no conditions = populate whole dataset col - even for occ/dolocc
        dfd[Symbol(string(m[:modelName])*"0")] =  eval(parse(r0))
        dfd[Symbol(string(m[:modelName])*"1")] =  eval(parse(r1))
    end
end
genFixedCols(dfx,dfd,modelsDict)


function genEQ(dfx::DataFrame, m::Dict, con::String="", ranef::Symbol=:empty,level::String="") # conditional equations
    model=string(m[:modelName])
    con=strip(con)
    cbuyer_pos_p1=length(con) > 0 ? "&(dfd[:buyer_pos_p1].==1)" : "(dfd[:buyer_pos_p1].==1)"
    con = m[:Buyer_Pos_P1_is1] ? con*cbuyer_pos_p1 : con
    con = length(con) > 0 ? con*"," : con
    B = ranef==:empty ? "":"+"*string(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==string(ranef))&(dfx[:parameter].==level),:adj_coef][1])
    vout = ranef==:empty ? "dfd[$con:"*model*"%%]" : "dfd[$con:"*model*"0]"
    vout0 = replace(vout,"%%","0")
    vout1 = replace(vout,"%%","1")*B
    return vout0, vout1   # ran_mean0 = random effects are zero, so leave out
end


function genFixedMeans(dfx::DataFrame)    
    for m in mlist
        model=string(m[:modelName])
        unadj0_eq, unadj1_eq = genEQ(dfx,m,"(dfd[:group].==%%)")
        adj0_eq, adj1_eq = genEQ(dfx,m)
        if inDFX(:unadj_mean_score0_eq) dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:unadj_mean_score0_eq] = unadj0_eq end
        if inDFX(:unadj_mean_score1_eq) dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:unadj_mean_score1_eq] = unadj1_eq end
        if inDFX(:adj_mean_score0_eq) dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:adj_mean_score0_eq] = adj0_eq end
        if inDFX(:adj_mean_score1_eq) dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:adj_mean_score1_eq] = adj1_eq end
        unadj0 = m[:Buyer_Pos_P1_is1] ? mean(exp(eval(parse(unadj0_eq)))) : mean(exp(eval(parse(unadj0_eq))) ./ (exp(eval(parse(unadj0_eq)))+1))
        unadj1 = m[:Buyer_Pos_P1_is1] ? mean(exp(eval(parse(unadj1_eq)))) : mean(exp(eval(parse(unadj1_eq))) ./ (exp(eval(parse(unadj1_eq)))+1))
        adj0 = m[:Buyer_Pos_P1_is1] ? mean(exp(eval(parse(adj0_eq)))) : mean(exp(eval(parse(adj0_eq))) ./ (exp(eval(parse(adj0_eq)))+1))
        adj1 = m[:Buyer_Pos_P1_is1] ? mean(exp(eval(parse(adj1_eq)))) : mean(exp(eval(parse(adj1_eq))) ./ (exp(eval(parse(adj1_eq)))+1))
        dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:unadj_mean_score0] = unadj0
        dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:unadj_mean_score1] = unadj1
        dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:adj_mean_score0] = adj0
        dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group") ,:adj_mean_score1] = adj1
    end
end
genFixedMeans(dfx)


function genRandMeans(dfx::DataFrame)
    for m in mlist
        model=string(m[:modelName])
        for r in m[:raneff]
            for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].== model)&(dfx[:ranef].==string(r)),:])
                println(model," Random Effect: ",row[:ranef]," ~~ ", row[:parameter]) 
                adj0_eq, adj1_eq = genEQ(dfx,m,"",Symbol(row[:ranef]),row[:parameter])
                #if (row[:ranef] == "creative")&(row[:model] =="occ")  println("TEST:: "*row[:ranef]*" : "*row[:parameter]*"  ::: ",adj0_eq) end
                if inDFX(:adj_mean_score0_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]) ,:adj_mean_score0_eq] = adj0_eq end
                if inDFX(:adj_mean_score1_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]) ,:adj_mean_score1_eq] = adj1_eq end
                adj0 = m[:Buyer_Pos_P1_is1] ?  mean(exp(eval(parse(adj0_eq)))) :  mean(exp(eval(parse(adj0_eq))) ./ (exp(eval(parse(adj0_eq)))+1))
                adj1 = m[:Buyer_Pos_P1_is1] ?  mean(exp(eval(parse(adj1_eq)))) :  mean(exp(eval(parse(adj1_eq))) ./ (exp(eval(parse(adj1_eq)))+1)) 
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]) ,:adj_mean_score0] = adj0
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]) ,:adj_mean_score1] = adj1            
            end
        end        
    end
end
genRandMeans(dfx)




function HackAttack(dfx::DataFrame)

    dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
    for m in mlist
        model=string(m[:modelName])
        fdod = dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")&(dfx[:model].==model),:adj_dod_effct][1]
        factor=0.7
        ub= fdod>0.0 ?  fdod+(fdod*factor) : fdod-(fdod*factor)
        lb= fdod>0.0 ?  fdod-(fdod*factor) : fdod+(fdod*factor)
        println("bounds : ",fdod," : ",lb,"",ub)
        #for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].== model)&(dfx[:adj_dod_effct].>ub)|(dfx[:adj_dod_effct].<lb),:])
        #    println("WRONG : raneff ; ",row[:adj_dod_effct], "",  "  total : ",  fdod  )
        #end
        
         #   m = mean(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef])
        #    f=groupBeta/m
        #    dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef]=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:ranef].==r)&(dfx[:parameter].!="none"),:adj_coef]*f
            
        
        println("-----------$model UB----------------")
        #for row in eachrow(dfr[(dfr[:TWOTAIL_PVAL].>79)&(dfr[:CNT_Model_HH].>50)&(dfr[:MODEL_DESC].!="Total Campaign")&(dfr[:dependent_variable].== model)&(dfr[:ADJ_DOD_EFFCT].>ub),:])    
        for row in eachrow(dfr[(dfr[:TWOTAIL_PVAL].>79)&(dfr[:MODEL_DESC].!="Total Campaign")&(dfr[:dependent_variable].== model)&((dfr[:ADJ_DOD_EFFCT].>ub)|(dfr[:ADJ_DOD_EFFCT].<lb)),:])
           println("WRONG : raneff ; ",row[:ADJ_MEAN_EXPSD_GRP],":",row[:ADJ_MEAN_CNTRL_GRP]," ~~ ",row[:ADJ_DOD_EFFCT], " cnt : ",row[:CNT_Model_HH]," sig : ", row[:TWOTAIL_PVAL],   "  total : ",  fdod  )
        end
        println("-----------$model LB----------------")
    #    
    #    for row in eachrow(dfr[(dfr[:TWOTAIL_PVAL].>79)&(dfr[:CNT_Model_HH].>50)&(dfr[:MODEL_DESC].!="Total Campaign")&(dfr[:dependent_variable].== model)&(dfr[:ADJ_DOD_EFFCT].<lb),:])    
        for row in eachrow(dfr[(dfr[:TWOTAIL_PVAL].>79)&(dfr[:MODEL_DESC].!="Total Campaign")&(dfr[:dependent_variable].== model)&((dfr[:ADJ_DOD_EFFCT].>ub)|(dfr[:ADJ_DOD_EFFCT].<lb)),:])
            println("WRONG : raneff ; ",row[:ADJ_MEAN_EXPSD_GRP],":",row[:ADJ_MEAN_CNTRL_GRP]," ~~ ",row[:ADJ_DOD_EFFCT], " cnt : ",row[:CNT_Model_HH]," sig : ", row[:TWOTAIL_PVAL],   "  total : ",  fdod  )
        end    
            
    end
end
HackAttack(dfx)






function genRawDataMeans(dfx::DataFrame) 
    for m in mlist   # ************************* FIXED **************************
        ex_Buyer_Pos_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pre_p1].==1) &" : ""
        precol=string(m[:logvarOrig])
        postcol=string(m[:y_var])
        model=string(m[:modelName])
        unadj_avg_cntrl_hh_pre_eq = "mean(dfd[ $ex_Buyer_Pos_P1 (dfd[:group].==0), :$precol] )"  
        unadj_avg_expsd_hh_pre_eq = "mean(dfd[ $ex_Buyer_Pos_P1 (dfd[:group].==1), :$precol] )"
        unadj_avg_cntrl_hh_pst_eq = "mean(dfd[ $ex_Buyer_Pos_P1 (dfd[:group].==0), :$postcol] )"
        unadj_avg_expsd_hh_pst_eq = "mean(dfd[ $ex_Buyer_Pos_P1 (dfd[:group].==1), :$postcol] )"
        unadj_avg_cntrl_hh_pre = eval(parse(unadj_avg_cntrl_hh_pre_eq)) 
        unadj_avg_expsd_hh_pre = eval(parse(unadj_avg_expsd_hh_pre_eq)) 
        unadj_avg_cntrl_hh_pst = eval(parse(unadj_avg_cntrl_hh_pst_eq)) 
        unadj_avg_expsd_hh_pst = eval(parse(unadj_avg_expsd_hh_pst_eq)) 
        if inDFX(:unadj_avg_cntrl_hh_pre_eq) dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pre_eq] = unadj_avg_cntrl_hh_pre_eq end
        if inDFX(:unadj_avg_expsd_hh_pre_eq) dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_expsd_hh_pre_eq] = unadj_avg_expsd_hh_pre_eq end
        if inDFX(:unadj_avg_cntrl_hh_pst_eq) dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pst_eq] = unadj_avg_cntrl_hh_pst_eq end
        if inDFX(:unadj_avg_expsd_hh_pst_eq) dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_expsd_hh_pst_eq] = unadj_avg_expsd_hh_pst_eq end
        dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pre] = unadj_avg_cntrl_hh_pre 
        dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_expsd_hh_pre] = unadj_avg_expsd_hh_pre 
        dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pst] = unadj_avg_cntrl_hh_pst 
        dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_expsd_hh_pst] = unadj_avg_expsd_hh_pst 
    end 
    for m in mlist  # *********************** RANDOM ************************ 
        precol=string(m[:logvarOrig])
        postcol=string(m[:y_var])
        model=string(m[:modelName])
        ex_Buyer_Pos_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pre_p1].==1) &" : ""
        for r in m[:raneff]
            rs=string(r)
            exposed=true # Hack till we figure this out
            for l in dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].!="none")&(dfx[:ranef].==rs)&(dfx[:model].==model),:parameter]
                unadj_avg_expsd_hh_pre_eq= "mean(dfd[ $ex_Buyer_Pos_P1  (dfd[:group] .== 1) & (dfd[:$rs] .== \"$l\"), :$precol])"
                unadj_avg_expsd_hh_pst_eq= "mean(dfd[ $ex_Buyer_Pos_P1  (dfd[:group] .== 1) & (dfd[:$rs] .== \"$l\"), :$postcol])"
                unadj_avg_cntrl_hh_pre_eq= "mean(dfd[ $ex_Buyer_Pos_P1  (dfd[:group] .== 0) & (dfd[:$rs] .== \"$l\"), :$precol])"
                unadj_avg_cntrl_hh_pst_eq= "mean(dfd[ $ex_Buyer_Pos_P1  (dfd[:group] .== 0) & (dfd[:$rs] .== \"$l\"), :$postcol])"
                if exposed
                    unadj_avg_cntrl_hh_pre_eq = string(dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pre][1])
                    unadj_avg_cntrl_hh_pst_eq = string(dfx[(dfx[:parameter].=="group") & (dfx[:modelType].=="GLM") & (dfx[:model].==model), :unadj_avg_cntrl_hh_pst][1])
                end
                unadj_avg_cntrl_hh_pre = eval(parse(unadj_avg_cntrl_hh_pre_eq))
                unadj_avg_expsd_hh_pre = eval(parse(unadj_avg_expsd_hh_pre_eq))
                unadj_avg_cntrl_hh_pst = eval(parse(unadj_avg_cntrl_hh_pst_eq))
                unadj_avg_expsd_hh_pst = eval(parse(unadj_avg_expsd_hh_pst_eq))
                
                if inDFX(:unadj_avg_cntrl_hh_pre_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_cntrl_hh_pre_eq] = unadj_avg_cntrl_hh_pre_eq end
                if inDFX(:unadj_avg_expsd_hh_pre_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_expsd_hh_pre_eq] = unadj_avg_expsd_hh_pre_eq end
                if inDFX(:unadj_avg_cntrl_hh_pst_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_cntrl_hh_pst_eq] = unadj_avg_cntrl_hh_pst_eq end
                if inDFX(:unadj_avg_expsd_hh_pst_eq) dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_expsd_hh_pst_eq] = unadj_avg_expsd_hh_pst_eq end
                println(model," ... ",r," : ",l," ~~~ ",unadj_avg_expsd_hh_pre_eq)
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_cntrl_hh_pre] = unadj_avg_cntrl_hh_pre
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_expsd_hh_pre] = unadj_avg_expsd_hh_pre
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_cntrl_hh_pst] = unadj_avg_cntrl_hh_pst
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].==l)&(dfx[:ranef].==rs)&(dfx[:model].==model),:unadj_avg_expsd_hh_pst] = unadj_avg_expsd_hh_pst
            end
        end
    end
end
genRawDataMeans(dfx)


function getCnts(dfd::DataFrame, ranef::String="", level::String="" )
    isBreak = length(ranef) == 0 ? false : true
    ex_re = isBreak ? "& (dfd[:$ranef] .== \"$level\")" : "  "  # need to have at least 2 spaces for single col belowbelow
    ex_re_single_col = isBreak ? ex_re[2:end]*"," : ""   
    ex_m  = "length(dfd[ (dfd[:buyer_pos_p1] .== 1) $ex_re ,1])"
    ex_mt = "length(dfd[ (dfd[:group] .== 1) & (dfd[:buyer_pos_p1] .== 1 ) $ex_re ,1])"
    ex_n  = "length(dfd[ $ex_re_single_col 1])"
    ex_nt = "length(dfd[ (dfd[:group] .== 1) $ex_re ,1])"
    ex_mc = "length(dfd[ (dfd[:group] .== 0) & (dfd[:buyer_pos_p1] .== 1 ) ,1])"
    ex_nc = "length(dfd[ (dfd[:group] .== 0) ,1])"
    cdict = OrderedDict()
    cdict[:M]  = eval(parse(ex_m))
    cdict[:Mt] = eval(parse(ex_mt))
    cdict[:Mc] = eval(parse(ex_mc))
    cdict[:N]  = eval(parse(ex_n))
    cdict[:Nt] = eval(parse(ex_nt))
    cdict[:Nc] = eval(parse(ex_nc))    
    if isBreak
        break_exposed = true  # DEFAULT Hack because all Julia model data is exposed so far   
        if break_exposed    #  If Exposed - default control (Nc & Mc) to total -- else count by rndfx        
            cdict[:M] = cdict[:Mt] + cdict[:Mc]   # Recalculate totals as Test + the defaulted Control
            cdict[:N] = cdict[:Nt] + cdict[:Nc]
        else # Need to recalculate MC & Nc with break conditions break
            ex_mc = "length(dfd[ (dfd[:group] .== 0) & (dfd[:buyer_pos_p1] .== 1 ) $ex_re ,1])"
            ex_nc = "length(dfd[ (dfd[:group] .== 0)  $ex_re ,1])"
            cdict[:Mc] = eval(parse(ex_mc))
            cdict[:Nc] = eval(parse(ex_nc))
        end
    end
    return cdict
end


function genCnts()
    c=getCnts(dfd)
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:M]=c[:M]
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:Mt]=c[:Mt]
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:Mc]=c[:Mc]
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:N]=c[:N]
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:Nt]=c[:Nt]
    dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:Nc]=c[:Nc]
    for r in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for l in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r),:parameter])
            c =  getCnts(dfd,r,l)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:M]=c[:M]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:Mt]=c[:Mt]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:Mc]=c[:Mc]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:N]=c[:N]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:Nt]=c[:Nt]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==r)&(dfx[:parameter].==l),:Nc]=c[:Nc]
        end
    end
end
genCnts()


function removeBreaks(dfx::DataFrame)   # Shouldn't be an issue going forwward -- remove breaks that don't exist in all 3 models
    dfo=deepcopy(dfx)
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            cnt = length(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0])
            if length(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0]) != 3
                pre=length(dfo[1])
                dfo=dfo[(dfo[:modelType].*dfo[:ranef].*dfo[:parameter].!="GLMM"*ranef*level),:]
                post=length(dfo[1])
                println("Remove : $ranef $level : $cnt :::: $pre -> $post")
            end
        end
    end
    return dfo
end
dfx=removeBreaks(dfx)


function genHHCounts()
    if isfile(root*"/hhcnts.csv")
        hhcnts = readtable(root*"/hhcnts.csv"); rename!(hhcnts,:class,:ranef); rename!(hhcnts,:level,:parameter); lowercase(hhcnts); hhcnts[:ranef] = lowercase(hhcnts[:ranef])
        tots = by(hhcnts,[:ranef], df -> sum(df[:hh]))
        hhcnts = join(hhcnts, tots, on = :ranef)
        rename!(hhcnts,:x1,:tot)
        hhcnts[:toti] = sum(hhcnts[:impressions])
        hhcnts[:weight]=hhcnts[:hh] ./ hhcnts[:tot]        
        for row in eachrow(hhcnts)
            dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLMM")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:cnt_expsd_hh] = row[:hh] 
            dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLMM")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:cnt_impressions] = row[:impressions]
        end    
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM"),:cnt_expsd_hh] = hhcnts[:tot][1] 
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM"),:cnt_impressions] = hhcnts[:toti][1]      
        return hhcnts
    else
        return NA
    end
end


function genWeights()
    hhcnts = genHHCounts()
    if !isna(hhcnts) 
        cols= [:adj_mean_score0, :adj_mean_score1, :unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst]        
        weightDF = join( dfx[(dfx[:modelType].=="GLMM"), vcat([:model,:ranef,:parameter],cols) ]
                         , hhcnts[[:ranef,:parameter,:weight]],  on = [:ranef,:parameter], kind = :left
                       )
        if length(weightDF[isna(weightDF[:weight]),:][1]) > 0
            println("ERROR in WEIGHTING - missing weights!!!")    #break()
        end     
        factorDF = by(weightDF, [:model,:ranef], df -> sum(0.0))[[:model,:ranef]];
        for col in cols
            creative_colname = Symbol(string(col)*"_Creative")
            total_colname = Symbol(string(col)*"_Total")
            factor_colname = Symbol(string(col)*"_Factor")
            weightDF[col] = weightDF[col] .* weightDF[:weight]
            crtDF = by(weightDF[!isnan(weightDF[col]),:], [:model,:ranef], df -> sum(df[col]))
            rename!(crtDF,:x1,creative_colname)
            ##crtv = rename!(by(weightDF[!isnan(weightDF[col]),:], [:model,:ranef], df -> DataFrame(zz = sum(df[col]))),:zz, creative_colname   )
            factorDF = join(factorDF, crtDF, on=[:model,:ranef],kind=:left) 
            #ex = "by(weightDF[!isnan(weightDF[:$col]),:], [:model,:ranef], df -> DataFrame($creative_colname = sum(df[:$col])))"
            #factorDF = join(factorDF, eval(parse(ex)), on=[:model,:ranef],kind=:left) 
            factorDF = join(factorDF, dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),[:model,col]], on=[:model],kind=:left)
            rename!(factorDF, col, total_colname)
            factorDF[factor_colname] = factorDF[total_colname] ./ factorDF[creative_colname]  
        end
        return factorDF
    else
        return NA
    end
end
    

function weightDFX()
    if isfile(root*"/hhcnts.csv")
        cols= [:adj_mean_score0, :adj_mean_score1, :unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst]
        dfw = genWeights() 
        for row in eachrow(dfw)
            for col in cols
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),col] = 
                      dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),col]  *  row[Symbol(string(col)*"_Factor")]    
            end
        end
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM"),:])
            M = row[:M]
            Mt = row[:Mt]
            Mc = row[:Mc]
            mean_score0 = row[:adj_mean_score0]
            mean_score1 = row[:adj_mean_score1]
            B = row[:adj_coef]  
            SE = row[:adj_stderr]  #println(row[:model]," : ",row[:ranef]," : ",row[:parameter],"   ::  ~",M,"~",Mt,"~",Mc,"~",mean_score0,"~",mean_score1,"~",B,"~",SE)
            function f!(x, fvec)
                Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*x[1])))*(Mt/M))   )  +
                         (   (mean_score0*exp((B-(SE*x[1])))*(Mc/M))    -   (mean_score0*(Mc/M))    )
                    Lb = Lb_pre/mean_score0
               fvec[1] = Lb
            end
            r=nlsolve(f!,[0.1])    #println(r)
            zvalue=r.zero[1]
            pvalue=2.0 * ccdf(Normal(), abs(zvalue))
            dfx[(dfx[:model].==row[:model])&(dfx[:modelType].=="GLMM")&(dfx[:parameter].==row[:parameter]),:pval] = pvalue
            println(row[:model]," : ",row[:ranef]," : ",row[:parameter],"   ::pval  ",pvalue)
        end
    else
        println("\n\n\n\n\n !!!NO HHCNTS.CSV FILE!!!!\n\n\n\n")
    end
end        
weightDFX()



function ConfIntrvl()
    for row in eachrow(dfx)
        runCI = false
        if (row[:modelType] == "GLM") & (row[:parameter]=="group")
            runCI = true
            if row[:model]=="pen"
                adj_dod_effct =  ((row[:adj_mean_score1] - row[:adj_mean_score0]) / row[:adj_mean_score0] ) *100
                row[:coef] = log((adj_dod_effct/100)+1)
                row[:pval] = 2.0 * ccdf(Normal(), abs(row[:coef]/row[:stderr]))
            end
            mean_score0=row[:unadj_mean_score0]  # TOTAL Confidence Intervals  # use unadj
            mean_score1=row[:unadj_mean_score1] 
            B=row[:coef]
            SE=row[:stderr] 
            Mt = row[:model]=="pen" ? row[:Nt] : row[:Mt]
            Mc = row[:model]=="pen" ? row[:Nc] : row[:Mc]
            M = row[:model]=="pen" ? row[:N] : row[:M]
            println("Total: ",mean_score0,"...",mean_score1," : ",B," ~ ",SE," ~ ",Mt," ~ ",Mc," ~ ",M)
            #ex="(($mean_score1*($Mt/$M))-($mean_score1*exp(-($B-($SE*\$z)))*($Mt/$M)))+(($mean_score0*exp(($B-($SE*\$z)))*($Mc/$M))-($mean_score0*($Mc/$M)))" 
            #println("EX ",row[:model],": ",ex)  #,"   ~~~   ", eval(parse(ex)))
        elseif row[:modelType]=="GLMM"
            runCI = true
            if row[:model]=="pen"
                adj_dod_effct =  ((row[:adj_mean_score1] - row[:adj_mean_score0]) / row[:adj_mean_score0] ) *100
                row[:adj_coef] = log((adj_dod_effct/100)+1)
                row[:pval] = 2.0 * ccdf(Normal(), abs(row[:adj_coef]/row[:adj_stderr]))
            end
            mean_score0=row[:adj_mean_score0] #RANDOM Confidence Intervals  # use adj
            mean_score1=row[:adj_mean_score1]
            B=row[:adj_coef]
            SE=row[:adj_stderr]
            Mt = row[:model]=="pen" ? row[:Nt] : row[:Mt]
            Mc = row[:model]=="pen" ? row[:Nc] : row[:Mc]
            M = row[:model]=="pen" ? row[:N] : row[:M]
            println("Random: ",mean_score0,"...",mean_score1," : ",row[:ranef],"_",row[:parameter]," ~ ",B," ~ ",SE," ~ ",Mt," ~ ",Mc," ~ ",M)
        end
        if runCI
            ZDict = Dict("onetail_80_pct_intrvl" => 0.84 ,"onetail_90_pct_intrvl" => 1.28, "twotail_80_pct_intrvl" => 1.28, "twotail_90_pct_intrvl" => 1.65)
            for (zscore_key, zscore) in ZDict    
                Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*zscore)))*(Mt/M))   )  +  ## ------------ Lower Bound ---------------
                         (   (mean_score0*exp((B-(SE*zscore)))*(Mc/M))    -   (mean_score0*(Mc/M))    )
                row[Symbol(zscore_key*"_lb")] = ( Lb_pre/mean_score0 ) * 100
                Ub_pre =  (     ( mean_score1*(Mt/M) )   -   ( mean_score1*exp(-(B+(SE*zscore)))*(Mt/M))   )  +  ## ------------ Upper Bound ---------------
                          (     ( mean_score0*exp((B+(SE*zscore)))*(Mc/M))  - (mean_score0*(Mc/M) )   )
                row[Symbol(zscore_key*"_ub")] = ( Ub_pre/mean_score0 ) * 100   
                println("RAND CI $zscore_key ($zscore) LB:",row[Symbol(zscore_key*"_lb")]," ~~ UB : ", row[Symbol(zscore_key*"_ub")])
            end            
        end
    end
end
ConfIntrvl()
#dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100

# **********************************************************************************************
# ***************************************** DOLHH **********************************************
# **********************************************************************************************
function genDHHMeans(dfx::DataFrame)
    println("TOTAL DOLHH")
    o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:adj_mean_score0][1]
    adj_mean_score0=o*y*p            
    o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:adj_mean_score1][1]
    adj_mean_score1=o*y*p
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
            println(ranef," ~~ ",level)
            o=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            y=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            p=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score0][1]
            adj_mean_score0=o*y*p            
            o=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            y=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            p=dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_mean_score1][1]
            adj_mean_score1=o*y*p
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
    return mDict
end


function ConfidenceIntervals(dfx::DataFrame)
    mDict=collectModels(dfx,"GLM")
    println("Running CI for : Total")
    mDict[:metakey] = "Total Campaign"        
    calcPValue_Opt(mDict)
    CIs_O(mDict)
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_pval] = mDict[:onetail_pval]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_pval] = mDict[:twotail_pval]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:pval] = mDict[:pval]
    for k in keys(ZDict)
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
    end
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")
            mDict[:metakey] = ranef*"~"*level        
            calcPValue_Opt(mDict)
            CIs_O(mDict)
            #dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_pval] = mDict[:onetail_pval]*100
            #dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_pval] = mDict[:twotail_pval]*100
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:pval] = mDict[:pval]
            for k in keys(ZDict)
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
            end
        end
    end
end
ConfidenceIntervals(dfx)
#dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
#dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),[:model,:adj_dod_effct,:onetail_80_pct_intrvl_lb,:onetail_80_pct_intrvl_ub]]


function genReport(dfx::DataFrame)
    genHHCounts();
    dfx=deepcopy(dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")|(dfx[:modelType].=="GLMM"),:])
    rep=[ :MODEL_DESC, :Model, :TIME_AGG_PERIOD, :START_WEEK,:END_WEEK, :dependent_variable,:CNT_EXPSD_HH, 
          #:UDJ_AVG_EXPSD_HH_PRE,:UDJ_AVG_EXPSD_HH_PST,:UDJ_AVG_CNTRL_HH_PRE,:UDJ_AVG_CNTRL_HH_PST,        
          :UDJ_AVG_EXPSD_HH_PRE, :UDJ_AVG_CNTRL_HH_PRE, :UDJ_AVG_EXPSD_HH_PST,:UDJ_AVG_CNTRL_HH_PST,
          :UDJ_DOD_EFFCT,:UDJ_DIFF_EFFCT,
          :ADJ_MEAN_EXPSD_GRP,:ADJ_MEAN_CNTRL_GRP,:ADJ_DOD_EFFCT,:TWOTAIL_PVAL,:ONETAIL_PVAL,:ABS_DIFF,
          :DOL_DIFF,:ONETAIL_80_PCT_INTRVL_UB,:ONETAIL_80_PCT_INTRVL_LB,:ONETAIL_90_PCT_INTRVL_UB,:ONETAIL_90_PCT_INTRVL_LB,
          :TWOTAIL_80_PCT_INTRVL_UB,:TWOTAIL_80_PCT_INTRVL_LB,:TWOTAIL_90_PCT_INTRVL_UB,:TWOTAIL_90_PCT_INTRVL_LB,
          :CNT_IMPRESSIONS,:TWOTAIL_PVAL_to_Campaign,:ONETAIL_PVAL_to_Campaign,:CNT_Model_HH
        ]
    dfx[isnan(dfx[:unadj_avg_expsd_hh_pre]),:unadj_avg_expsd_hh_pre] = 0.0  # Sometime there are no records for subset
    dfx[isnan(dfx[:unadj_avg_expsd_hh_pst]),:unadj_avg_expsd_hh_pst] = 0.0  # Sometime there are no records for subset
    dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
    dfx[:unadj_dod_effct] = ( (( dfx[:unadj_avg_expsd_hh_pst] .- dfx[:unadj_avg_expsd_hh_pre]) .- (dfx[:unadj_avg_cntrl_hh_pst ] .- dfx[:unadj_avg_cntrl_hh_pre]))  ./  dfx[:unadj_avg_cntrl_hh_pst] ) *100
    dfx[:unadj_diff_effct] = ((dfx[:unadj_avg_expsd_hh_pst] .- dfx[:unadj_avg_cntrl_hh_pst]) ./ dfx[:unadj_avg_cntrl_hh_pst] )*100
    dfx[:model_desc] = dfx[:ranef]*" (".*dfx[:parameter]*")"
    dfx[(dfx[:modelType].=="GLM"),:model_desc] = "Total Campaign"
    dfx[:abs_diff] = dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]
    dfx[:dol_diff] = dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]
    dfx[:cnt_model_hh] = dfx[:Nt]
    dfx[findin(dfx[:model],["occ","dolocc"]),:cnt_model_hh] = dfx[findin(dfx[:model],["occ","dolocc"]),:Mt]       
    dfx[:onetail_pval] = (1-(dfx[:pval] ./ 2)) * 100  # sdf[:onetail_pval_raw] = (1-(sdf[:Praw] ./ 2)) * 100
    dfx[:twotail_pval] = (1-dfx[:pval]) * 100         # sdf[:twotail_pval_raw] = (1-sdf[:Praw]) * 100
    dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"),:adj_pval] = 2.0 * ccdf(Normal(),abs(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"),:adj_coef] ./ dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"),:adj_stderr]))
    dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :onetail_pval_to_campaign] =  (1-(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :adj_pval] ./ 2) ) * 100
    dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :twotail_pval_to_campaign] =  (1-dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :adj_pval] ) * 100

    dfx[:empty] = ""
    dfx[:orderModel] = map(x-> x=="pen"? 1 : x=="occ" ? 2 : x=="dolocc" ? 3 : x=="dolhh" ? 4 : 9999 ,dfx[:model])
    dfx[:orderFixedRand] = map(x-> x=="GLM"? 1 : 9999 ,dfx[:modelType])
    sort!(dfx, cols = [:orderFixedRand,:ranef,:parameter,:orderModel])
    dfr=dfx[[:model_desc,:empty,:empty,:empty,:empty,:model,:cnt_expsd_hh,
             #:unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst,:unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst,
             :unadj_avg_expsd_hh_pre, :unadj_avg_cntrl_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pst,
             :unadj_dod_effct,:unadj_diff_effct,
             :adj_mean_score1,:adj_mean_score0,:adj_dod_effct,:twotail_pval,:onetail_pval,:abs_diff,:dol_diff,
             :onetail_80_pct_intrvl_ub,:onetail_80_pct_intrvl_lb,:onetail_90_pct_intrvl_ub,:onetail_90_pct_intrvl_lb,
             :twotail_80_pct_intrvl_ub,:twotail_80_pct_intrvl_lb,:twotail_90_pct_intrvl_ub,:twotail_90_pct_intrvl_lb
             ,:cnt_impressions ,:twotail_pval_to_campaign,:onetail_pval_to_campaign ,:cnt_model_hh
        ]]
    names!(dfr,rep)
    
    dfr[(dfr[:MODEL_DESC].=="Total Campaign")|(dfr[:dependent_variable].=="dolhh"),:TWOTAIL_PVAL_to_Campaign] = NaN
    dfr[(dfr[:MODEL_DESC].=="Total Campaign")|(dfr[:dependent_variable].=="dolhh"),:ONETAIL_PVAL_to_Campaign] = NaN
    dfr[(dfr[:dependent_variable].!="dolhh"),:CNT_EXPSD_HH] = NA
    dfr[(dfr[:dependent_variable].!="dolhh"),:CNT_IMPRESSIONS] = NA
    return dfr
end
dfr = genReport(dfx)
save_dfr(root, dfr)
















