

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
#x=genRndMDF(cnts, mocc, mdolocc, mpen, false)
#CIs_O_LB(df2dict(x[x[:key].=="estimated_hh_income (H)",:]), 0.84, float("0.001"))



function CIs_O_UB(iDict::OrderedDict, zscore::Float64, iAccuracy::Float64=0.000000001)
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
#Bsum=B1+B2+B3 + o_B0+y_B0+p_B0
Bsum=(B1+B2+B3)-(o_B0+y_B0+p_B0)
   
    ######CONFIDENCE INTERVAL - UB ########
    ztot = Bsum+(zscore*SEsq)

    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))
    @variable(m, Bocc >= (B1-o_B0))
    @variable(m, Bdolocc >= (B2-y_B0))
    @variable(m, Bpen >= (B3-p_B0))
    @NLobjective(m, Max, ( (((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))
                                 * ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))
                                 * ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
                                  )
                                 -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))
                                  *((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))
                                  *((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
                                  )
                               )
                       )
    @constraint(m, (0.0000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAccuracy)) 
    status = solve(m)
    mval=getobjectivevalue(m)
    mval_out=mval/(o_mean_score0*y_mean_score0*p_mean_score0)
    return return status == :Optimal ? mval_out : -Inf  #mval_out    
end
#x=genRndMDF(cnts, mocc, mdolocc, mpen, false)
#CIs_O_UB(df2dict(x[x[:key].=="estimated_hh_income (H)",:]), 0.84, float("0.001"))


function CIs_O(iDict::OrderedDict)
    AccArr= [ "0.000000001","0.00000001","0.0000001","0.000001","0.00001","0.0001","0.001","0.01"]
    for (zscore_key,zscore) in  ZDict          
        pref="LB "*zscore_key[1:10]*" ("*iDict[:metakey]"):= "
        preflen=length(pref)
        dkey=Symbol(zscore_key*"_lb")
        for iAcc in AccArr
           print(pref," - "*iAcc*", ",iDict)       ;pref=lpad("", preflen, " ")
           iDict[dkey] = CIs_O_LB(iDict, zscore, float(iAcc))
           if (iDict[dkey] != -Inf) & (iDict[dkey] != Inf)
               println("            : Confidence Interval : ", string(iDict[dkey]) )
               break
           end
        end
        if (iDict[dkey] == -Inf) | (iDict[dkey] == Inf) 
             println( "!!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!!!!!!!   LB ",zscore_key[1:10]," := ",iDict[:metakey]," == !!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!" ) 
        end
       
        pref="UB "*zscore_key[1:10]*" ("*iDict[:metakey]"):= "   
        dkey=Symbol(zscore_key*"_ub")
        for iAcc in AccArr
           print(pref," - "*iAcc*", ")       ;pref=lpad("", preflen, " ")
           iDict[dkey] = CIs_O_UB(iDict, zscore, float(iAcc))
           if (iDict[dkey] != -Inf) & (iDict[dkey] != Inf)
               println("            : Confidence Interval : ", string(iDict[dkey]) )
               break
           end
        end   
        if (iDict[dkey] == -Inf) | (iDict[dkey] == Inf) 
            println( "!!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!   UB ",zscore_key[1:10]," := ",iDict[:metakey]," == !!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!" ) 
        end
    end
    return iDict
end



function calcPValue_Opt(iDict::OrderedDict)
    dout = iDict 
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
    SEsq=sqrt(SE1^2+SE2^2+SE3^2)
    o_mean_score0 = get(iDict, :o_mean_score0, NA)
    o_mean_score1 = get(iDict, :o_mean_score1, NA)
    y_mean_score0 = get(iDict, :y_mean_score0, NA)
    y_mean_score1 = get(iDict, :y_mean_score1, NA)
    p_mean_score0 = get(iDict, :p_mean_score0, NA)
    p_mean_score1 =get(iDict, :p_mean_score1, NA)
    Bsum=B1+B2+B3
    dout[:Bsum] = Bsum  
    ###### PVALUE - ONE & TWO ########
    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))
    @variable(m, Bocc <= B1)
    @variable(m, Bdolocc <= B2)
    @variable(m, Bpen <= B3)
    @objective(m, Max, (((Bocc+Bpen+Bdolocc)-Bsum)/SEsq ))
    @NLconstraint(m, 0.00000 <= ((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))
	                           * ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))
	                           * ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
	                             )
	                           -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))
	                           *((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))
	                           *((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
	                             )
	                          ) 
	               <= 0.00001
                    )
    #print(m)
    status = solve(m)
    zvalue=getobjectivevalue(m)
    pvalue=2.0 * ccdf(Normal(), abs(zvalue))
    two_tail = 1-pvalue     
    one_tail = 1-(pvalue/2)
    dout[:pval] =  pvalue
    dout[:onetail_pval] = one_tail
    dout[:twotail_pval] = two_tail
    println("z-value: ", string(zvalue)," --> p-value: ",string(two_tail))
    return dout           
end





function np_val(iDict::OrderedDict)
    M = get(iDict, :M, NA)
    Mt = get(iDict, :Mt, NA)
    Mc = get(iDict, :Mc, NA)
    mean_score0 = get(iDict, :adj_mean_cntrl_grp, NA)
    mean_score1 = get(iDict, :adj_mean_expsd_grp, NA)
    B = get(iDict, :B1_combo, NA)
    SE = get(iDict, :SE1_combo, NA)
    println(get(iDict,:key, NA),"~",M,"~",Mt,"~",Mc,"~",mean_score0,"~",mean_score1,"~",B,"~",SE)
    
    function f!(x, fvec)
        Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*x[1])))*(Mt/M))   )  +
                 (   (mean_score0*exp((B-(SE*x[1])))*(Mc/M))    -   (mean_score0*(Mc/M))    )
        Lb = Lb_pre/mean_score0
        fvec[1] = Lb
    end
    r=nlsolve(f!,[0.1])
    #println(r)
    zvalue=r.zero[1]
    pvalue=2.0 * ccdf(Normal(), abs(zvalue))
    two_tail = 1-pvalue     
    one_tail = 1-(pvalue/2)
    return one_tail, two_tail
end
