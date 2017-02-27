ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg" 
@everywhere insert!(Base.LOAD_CACHE_PATH, 1, "/mapr/mapr04p/analytics0001/analytic_users/jpkg/lib/v0.5")
@everywhere pop!(Base.LOAD_CACHE_PATH)
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
@everywhere using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve

root = !isdefined(:root) ? pwd() : root
MP = !isdefined(:MP) ? false : MP

const cfgDefaults=OrderedDict( :P2_Competitor => true
                        ,:pvalue_lvl => 0.20  #pvalue_lvl = 0.20 
                        ,:excludedBreaks => String[]    #["estimated_hh_income","hh_age","number_of_children_in_living_un","person_1_gender"]
                        ,:excludedLevels => ["none"]
                        ,:excludedKeys => String[]
                        ,:exposed_flag_var => :exposed_flag
                        ,:sigLevel => "0.2"
                        ,:random_demos => [:estimated_hh_income,:hh_age,:number_of_children_in_living_Un,:person_1_gender]
                        ,:random_campaigns => []
                        ,:dropvars => [:exposed_flag,:dma,:dma_code]  
                        ,:scoring_vars => [:Prd_1_Net_Pr_PRE,:Prd_1_Net_Pr_POS,:Buyer_Pos_P0,:Buyer_Pre_P0]
                        ,:occ_y_var => :Trps_POS_P1
                        ,:occ_logvar => :Trps_PRE_P1
                        ,:dolocc_y_var => :Dol_per_Trip_POS_P1
                        ,:dolocc_logvar => :Dol_per_Trip_PRE_P1
                        ,:pen_y_var => :Buyer_Pos_P1
                        ,:pen_logvar => :Buyer_Pre_P1
                        ,:TotalModelsOnly=>false
                       )
cfgDefaults = mergeDict(dict_Sym(read_cfg(root)),cfgDefaults)
dfd=read_orig(root)
#using Feather
#Feather.write(root*"/orig.feather", dfd)
#Feather.read(root*"/orig.feather")


function getCFG(dfd::DataFrame)
    cfg=StatLib.loadCFG(cfgDefaults, pwd()*"/app.cfg")
    cfg[:allrandoms] = vcat(cfg[:random_demos],cfg[:random_campaigns]) 
    if !haskey(cfg,:ProScore)
        ps = filter(x->contains(string(x), "MODEL"), names(dfd)) 
        cfg[:ProScore] = length(ps) > 0 ? ps[1] : :MISSING_MODEL_VARIABLE_IN_DATA
    end
    cfg[:num_products] = length( grep("Buyer_Pre_P",names(dfd)))-1 
    cfg[:random_campaigns] = intersect(cfg[:random_campaigns],names(dfd)) 
    cfg[:occ_logvar_colname] = Symbol("LOG_"*string(cfg[:occ_logvar]))
    cfg[:dolocc_logvar_colname] = Symbol("LOG_"*string(cfg[:dolocc_logvar]))
    cfg[:pen_logvar_colname] = Symbol("LOG_"*string(cfg[:pen_logvar]))
    cfg[:all_mandatory_vars] = vcat( [:experian_id, cfg[:pen_y_var],cfg[:pen_logvar], cfg[:occ_y_var], cfg[:occ_logvar], cfg[:dolocc_y_var], cfg[:dolocc_logvar] ], cfg[:random_demos] )
    return cfg
end
cfg=getCFG(dfd)


function isValid(dfd::DataFrame,cfg::OrderedDict)
    function checkValid(iarr::Array{Symbol})  length(setdiff(iarr, names(dfd))) > 0 ? false : true end
    checkValid(iSym::Symbol) = iSym in names(dfd)
    !checkValid(cfg[:all_mandatory_vars]) ? error("ERROR: Not all mandatory_vars in dataset ") : println("VALID : mandatory_vars") 
    !checkValid(cfg[:scoring_vars]) ? error("ERROR: Not all scoring_vars in dataset ") : println("VALID : scoring_vars") 
    !checkValid(cfg[:ProScore]) ? error("ERROR: Invalid ProScore ") : println("VALID : ProScore") 
end
isValid(dfd,cfg)

function reworkCFG!(dfd::DataFrame,cfg::OrderedDict)
    xflags=[cfg[:exposed_flag_var],:exposed_flag,:Nonbuyer_Pre_P1]  # :Nonbuyer_Pre_P1 no longer used
    features = setdiff(names(dfd),xflags)
    if :state in names(dfd)
        cfg[:xVarsDemos]=vcat([:experian_id, :banner, :Prd_0_Qty_PRE],features[findfirst(features, :state):findfirst(features, :Mosaic)])   # exclude demos between....
    else
        cfg[:xVarsDemos]=Symbol[]
    end
    cfg[:xVarsDemos]=setdiff(cfg[:xVarsDemos],[:person_1_gender,:number_of_children_in_living_Un]) # exclude person, # from non-used demos
    cfg[:xVarsPost] = grep(["POS","Pos","Buyer_Pre_","Nonbuyer_Pre_"],features)   #exclude POST variables
    cfg[:iVarsPREPOS] = grep(["PRE_POS","Pre_Pos"],features) 
    features=setdiff(features, setdiff(  vcat(cfg[:xVarsDemos],cfg[:xVarsPost]) ,cfg[:iVarsPREPOS])  )
    cfg[:xVarsP0] =  setdiff(grep("0", features ),[cfg[:ProScore]] ) #exclude category variables(exclude P0's)  
    features=setdiff(features,cfg[:xVarsP0])
    cfg[:xVarsReports] = grep(["Perc_","Pr_per_"],features)  #exclude Reporting vars
    features=setdiff(features,cfg[:xVarsReports])
    cfg[:xvars] = vcat(xflags,cfg[:xVarsDemos], cfg[:xVarsPost],cfg[:xVarsP0], cfg[:xVarsReports],cfg[:dropvars])
    features=setdiff(features,cfg[:dropvars])
    cfg[:ivars] = vcat(cfg[:iVarsPREPOS],cfg[:all_mandatory_vars],cfg[:scoring_vars])
    cfg[:ALL_vars_to_exclude] = setdiff(cfg[:xvars], cfg[:ivars])
    dfd[:group] = dfd[cfg[:exposed_flag_var]]
    dfd[:panid] = dfd[:experian_id] 
    features = setdiff(unique(vcat(features,cfg[:iVarsPREPOS],cfg[:all_mandatory_vars],cfg[:scoring_vars],[:group,:panid])  ) ,[:experian_id] )
    cfg[:negativevars] = grep(vec(hcat([["P"*string(i),string(i)*"_"] for i=3:cfg[:num_products]]...)),features)   # get variables that need to have negative sign
    cfg[:positivevars] = grep(["P1", "1_","MODEL"], features)    # get variables that need to have positive sign
    if cfg[:P2_Competitor] == true  
        cfg[:negativevars] = unique(vcat(cfg[:negativevars],grep(["P2","2_"],features))) 
    else
        cfg[:positivevars] = unique(vcat(cfg[:positivevars],grep(["P2","2_"],features))) 
    end
    return vcat(features,[:Prd_0_Net_Pr_PRE]) 
end
DSvars = reworkCFG!(dfd,cfg)



function pre_dataPrep(dfd::DataFrame)
    dfd[:iso]=false
    for col in [cfg[:occ_y_var],cfg[:occ_logvar],cfg[:pen_y_var],cfg[:pen_logvar], :Buyer_Pre_P0, :group ]
        if typeof(dfd[col]) in [DataArray{String,1}]
            dfd[col] = [x in ["//N","\\N"] ? 0 : parse(Int64,x) for x in Array(dfd[col])]
        else
            dfd[ DataArrays.isna(dfd[col]), col] = 0
        end
    end
    for col in [cfg[:dolocc_y_var],cfg[:dolocc_logvar], :Prd_0_Net_Pr_PRE]
        if typeof(dfd[col]) in [DataArray{String,1}]
            dfd[col] = [x in ["//N","\\N"] ? 0.0 : parse(Float64,x) for x in Array(dfd[col])]
        else
            dfd[ DataArrays.isnan(dfd[col]), col] = 0.0
        end
    end
    vars = DataFrame(names=names(dfd[DSvars]),eltypes=eltypes(dfd[DSvars]))
    for c in setdiff(vars[findin(vars[:eltypes],[String]),:names],cfg[:allrandoms])
        println("Convert String: String->Numeric: ",c)
        try dfd[c] = [x in ["//N","\\N"] ? 0.0 : parse(Float64,x) for x in Array(dfd[c])]
            dfd[DataArrays.isnan(dfd[c]), c] = 0.0
        catch e 
            try dfd[c] = [x in ["//N","\\N"] ? 0 :  x==0.0 ? 0 : parse(Int64,x) for x in Array(dfd[c])] 
                dfd[ DataArrays.isna(dfd[c]), c] = 0
            catch e 
               println("  (failed)")  
            end
        end
    end
    dfd[dfd[:number_of_children_in_living_Un].>=4,:number_of_children_in_living_Un] = 4 
    vars = DataFrame(names=names(dfd[DSvars]),eltypes=eltypes(dfd[DSvars]))
    for c in setdiff(vars[findin(vars[:eltypes],[Float64]),:names],cfg[:allrandoms])  
        println("Replace Float64 NaN (0.0): ",c)
        dfd[ DataArrays.isnan(dfd[c]), c] = 0.0
    end
    for c in  setdiff(vars[findin(vars[:eltypes],[Int64]),:names],cfg[:allrandoms])    
        println("Replace Int64 NaN (0) : ",c)
        dfd[ DataArrays.isnan(dfd[c]), c] = 0
    end 
    dfd[findin(dfd[:estimated_hh_income],["U","L"]),:estimated_hh_income]="L" 
    
    dfd[dfd[:person_1_gender].=="U",:iso] = true  #..... hash
    #****!!!!!!!!!!!!!!!  --- ISSUE WHEN RAND EFFECTS ARE BAD!!!! ---  !!!!!!!!!!!!!!!!!!!!!*****#
    for r in cfg[:random_campaigns]
        dfd[findin(dfd[r],["\\N","NULL","0","NONE"])  ,r] ="none"
        dfd[ (dfd[:iso].==false) & (dfd[r].!="none") & (dfd[:group].==0) ,:iso] = true  # so all ctrl data other than 
        dfd[ (dfd[:iso].==false) & (dfd[r].=="none") & (dfd[:group].!=0) ,:iso] = true
    end
end
println("pre_dataPrep")
pre_dataPrep(dfd)

println("Pre_3STD_out")
function Pre_3STD_out(dfd::DataFrame)
    df_cat_pre = dfd[(dfd[:iso].==false)&(dfd[:Buyer_Pre_P0] .==1) , [:Prd_0_Net_Pr_PRE,:panid]]
    df_cat_pos = dfd[(dfd[:iso].==false)&(dfd[:Buyer_Pre_P0] .==0) & (dfd[:Buyer_Pos_P0] .==1) , [:panid]]
    median_df = median(df_cat_pre[:Prd_0_Net_Pr_PRE])
    df_cat_pre[:Prd_0_Net_Pr_PRE_med1] = abs(df_cat_pre[:Prd_0_Net_Pr_PRE]-median_df)
    MAD=median(df_cat_pre[:Prd_0_Net_Pr_PRE_med1])
    df_cat_pre[:Prd_0_Net_Pr_PRE_med2] = (0.6745*(df_cat_pre[:Prd_0_Net_Pr_PRE]-median_df))/MAD
    df_cat_pre_zsc = df_cat_pre[abs(df_cat_pre[:Prd_0_Net_Pr_PRE_med2]) .< 3.5,:]
    df_cat_pre_zsc_1 = df_cat_pre_zsc[:,[:panid]]
    df_cat_pre_zsc_f = vcat(df_cat_pos,df_cat_pre_zsc_1)
    #dfd_pout =  join(dfd, df_cat_pre_zsc_f, on =  :panid , kind = :inner);
    rows2remove = setdiff(dfd[dfd[:iso].==false, :panid],df_cat_pre_zsc_f[:panid])
    dfd[findin(dfd[:panid],rows2remove),:iso]=true
end
Pre_3STD_out(dfd);


function Restrict50Buyers()
    for (k,v) in Dict(r=>countmap(dfd[(dfd[:iso].==false)&(dfd[:Buyer_Pos_P1].==1), r]) for r in cfg[:random_campaigns])
            lvls=String[]
            for (kl,vl) in v 
                if vl < 50
                    println("$k :: $kl to Other")
                    dfd[dfd[k].==kl,k] = "Other" 
                end  
            end
    end
end
println("Restrict50Buyers")
Restrict50Buyers()    

println("MatchMe")       
function MatchMe(dfd::DataFrame,cfg::OrderedDict)
    df=dfd[dfd[:iso].==false,[:panid,:group,:Buyer_Pre_P1,cfg[:ProScore]]]
    df_exp     = df[df[:group].==1,:]
    df_unexp   = df[df[:group].==0,:]
    df_exp_dim   = nrow(df_exp)
    df_unexp_dim = nrow(df_unexp)
    new_unexp_dim = df_unexp_dim*(df_unexp_dim>2000000 ? 0.3 : df_unexp_dim>1000000 ? 0.4 : df_unexp_dim>750000 ? 0.6 : 0.7)
    if length(string(cfg[:ProScore])) == 0
        df_unexp_1 =  df[(df[:group].==0)&(df[:Buyer_Pre_P1].==1),:]
        df_unexp_0 =  df[(df[:group].==0)&(df[:Buyer_Pre_P1].==0),:]
        
        df_exp_1_dim = nrow(df[(df[:group].==1)&(df[:Buyer_Pre_P1].==1),:])
        df_exp_0_dim = nrow(df[(df[:group].==1)&(df[:Buyer_Pre_P1].==0),:])
        df_unexp_1_dim = nrow(df_unexp_1)
        df_unexp_0_dim =  nrow(df_unexp_0)     
        dim_sample0 = round(Int64, (new_unexp_dim-df_exp_dim  ) / (1+(df_exp_1_dim / df_exp_0_dim)) )
        dim_sample1 = round(Int64, new_unexp_dim - df_exp_dim - dim_sample0)    
        
        new_df_unexp_1 = df_unexp_1[sample(1:size(df_unexp_1,1), dim_sample1 ),:]
        new_df_unexp_0 = df_unexp_0[sample(1:size(df_unexp_0,1), dim_sample0 ),:]
        dfd_sample  = vcat(df_exp,new_df_unexp_1,new_df_unexp_0)
    elseif length(string(cfg[:ProScore])) > 0    
        sample_control_data=similar(df_unexp, 0)
        for (key, value) in countmap(df_exp[cfg[:ProScore]])
            sample_dim=round(Int64,new_unexp_dim*(value/df_exp_dim))
            temp_data = df_unexp[df_unexp[cfg[:ProScore]].==key,:]
            if sample_dim > size(temp_data,1) sample_dim = size(temp_data,1) end  #??
            samp_data = temp_data[sample(1:size(temp_data,1), sample_dim, replace=false),:]
            sample_control_data = vcat(sample_control_data,    samp_data   )
        end
        sample_data = vcat(sample_control_data,df_exp)
        sample_df_unexp = sample_data[sample_data[:group].==0,:]
        sample_df_exp   = sample_data[sample_data[:group].==1,:]
        sample_df_unexp_1 = sample_data[(sample_data[:group].==0)&(sample_data[:Buyer_Pre_P1].==1),:]  
        sample_df_unexp_0 = sample_data[(sample_data[:group].==0)&(sample_data[:Buyer_Pre_P1].==0),:]    
        sample_df_exp_1 =  sample_data[(sample_data[:group].==1)&(sample_data[:Buyer_Pre_P1].==1) ,:]
        sample_df_exp_0 =  sample_data[(sample_data[:group].==1)&(sample_data[:Buyer_Pre_P1].==0) ,:]
        sample_df_unexp_1_dim = nrow(sample_df_unexp_1)
        sample_df_unexp_0_dim = nrow(sample_df_unexp_0)
        sample_df_exp_1_dim = nrow(sample_df_exp_1)
        sample_df_exp_0_dim = nrow(sample_df_exp_0)
        dim_sampleA = round(Int64,(sample_df_exp_1_dim/sample_df_exp_0_dim)*sample_df_unexp_0_dim)
        dim_sampleB = round(Int64,(sample_df_exp_0_dim/sample_df_exp_1_dim)*sample_df_unexp_1_dim) 
        if sample_df_unexp_1_dim/sample_df_unexp_0_dim > sample_df_exp_1_dim/sample_df_exp_0_dim
            new_df_unexp_1 = sample_df_unexp_1[sample(1:sample_df_unexp_1_dim,dim_sampleA , replace=false),:]
            dfd_sample = vcat(sample_df_exp,new_df_unexp_1,sample_df_unexp_0)
        else
            new_df_unexp_0 = sample_df_unexp_0[sample(1:sample_df_unexp_0_dim, dim_sampleB, replace=false ),:]
            dfd_sample = vcat(sample_df_exp,sample_df_unexp_1,new_df_unexp_0)
        end
    end 
    rows2remove = setdiff(dfd[dfd[:iso].==false, :panid],dfd_sample[:panid])
    dfd[findin(dfd[:panid],rows2remove),:iso]=true     
    #return dfd[dfd[:iso].==false, : ]  
end
MatchMe(dfd,cfg)
    
println("filter DFD")
#x=similar(dfd,length(dfd[(dfd[:iso].==false),:panid]))
lowercase!(dfd)
cfg=lowercase(cfg)

#DSvars=setdiff(DSvars,[ :core_based_statistical_areas,:latitude,:county_code,
#                        :person_1_birth_year_and_month,:person_1_combined_age,:person_1_marital_status,:recipient_reliability_code,
#                        :household_composition,:person_1_person_type,:homeowner_combined_homeowner,:homeowner_probability_model,
#                        :dwelling_type,:length_of_residence,:dwelling_un_size, :mail_responder
#              ])
dfd=dfd[dfd[:iso].==false, lowercase(DSvars) ] 
    
    
iocc = Dict(:modelName=>:occ, :raneff=>cfg[:random_campaigns], :y_var=>:trps_pos_p1, :logvar=>:LOG_trps_pre_p1, :logvarOrig=>:trps_pre_p1 )
idolocc = Dict(:modelName=>:dolocc, :raneff=>cfg[:random_campaigns], :y_var=>:dol_per_trip_pos_p1, :logvar=>:LOG_dol_per_trip_pre_p1, :logvarOrig=>:dol_per_trip_pre_p1)
ipen = Dict(:modelName=>:pen, :raneff=>cfg[:random_campaigns], :y_var=>:buyer_pos_p1, :logvar=>:LOG_buyer_pre_p1, :logvarOrig=>:buyer_pre_p1 ) 
for m in [iocc,idolocc,ipen] dfd[m[:logvar]]=log(Array(dfd[m[:logvarOrig]]+1)) end
function genExcludeVars!(iocc::Dict,idolocc::Dict,ipen::Dict)  
    iocc[:exclude_vars] = [ :iso, :buyer_pos_p1, iocc[:logvarOrig],idolocc[:y_var],idolocc[:logvar],idolocc[:logvarOrig] ,ipen[:y_var],ipen[:logvar],ipen[:logvarOrig] ]
                              
    idolocc[:exclude_vars] = [ :iso, :buyer_pos_p1, idolocc[:logvarOrig] ,iocc[:y_var],iocc[:logvar],iocc[:logvarOrig] ,ipen[:y_var],ipen[:logvar],ipen[:logvarOrig] ]
                                     
    ipen[:exclude_vars] = [ :iso, ipen[:logvarOrig] ,iocc[:y_var],iocc[:logvar],iocc[:logvarOrig],idolocc[:y_var],idolocc[:logvar],idolocc[:logvarOrig] ]         
end
genExcludeVars!(iocc,idolocc,ipen)

println("featureSelection")
function featureSelection(dfd::DataFrame, m::Dict)
    function rmVars(v::DataArray{Any}) rmVars(convert(Array{Symbol},v)) end
    function rmVars(v::Array{Any}) rmVars(convert(Array{Symbol},v)) end
    function rmVars(v::Array{Symbol})
        return setdiff(vars,v)  
    end        
    vars=setdiff(names(dfd),   vcat(m[:exclude_vars],m[:y_var],:panid,cfg[:random_demos],cfg[:random_campaigns],cfg[:scoring_vars]) )
    upperModName = uppercase( string(  m[:modelName]  )  )
    println( upperModName*" : SingleValue") #SingleValue
    removed_SingleLevelVars=FS_singleLevel(dfd,vars)
    vars = rmVars(removed_SingleLevelVars)
    println(upperModName*" : Singularity : "*string(genF(m[:y_var],vars))) # Singularity
    singularity_x = checksingularity(genF(m[:y_var],vars), dfd)
    vars = rmVars(singularity_x)
    println(upperModName*" : PVals") #PVals
    f=genF(m[:y_var],vars)
    g1 = glm( f, dfd[convert(Array{Symbol},vcat(f.lhs,f.rhs.args[3:end]))]  , m[:dist], m[:lnk] )
    sdf = coefDF(g1,true)
    g1_x= sdf[sdf[:pval].>0.7,:parameter]
    vars = rmVars(convert(Array{Any},g1_x))
    function chkSigns(m::Dict, vars::Array{Symbol}, dfd::DataFrame, cfg::OrderedDict)  # Pvalue & Signs
        vars=unique(vcat(vars,[:group]))
        f=genF(m[:y_var],vars)
        g = glm( f, dfd[convert(Array{Symbol},vcat(f.lhs,f.rhs.args[3:end]))]  , m[:dist], m[:lnk] )
        sdf = coefDF(g1, true)
        neutralvars = setdiff(vars,vcat(cfg[:negativevars],cfg[:positivevars])) 
        neg=intersect(cfg[:negativevars],sdf[sdf[:coef].<0,:parameter])
        pos=intersect(cfg[:positivevars],sdf[sdf[:coef].>0,:parameter])
        varstokeep = intersect(vcat(neutralvars, pos,neg) ,  sdf[sdf[:pval].<cfg[:pvalue_lvl] ,:parameter] )
        varstokeep =  convert(Array{Symbol},varstokeep)
        return g, varstokeep
    end
    println(upperModName*" : SIGN Check 1") 
    (g3, initialvars) = chkSigns(m, vars, dfd, cfg)
    println(upperModName*" : SIGN Check 2") 
    (g4, vars_2) = chkSigns(m, initialvars, dfd, cfg)
    function getCorrVars(m::RegressionModel, dfd::DataFrame, vars_2::Array{Symbol})
        sdf = coefDF(m, true)
        vars_2 = convert(Array{Symbol}, setdiff(vars_2, filter(x-> typeof(dfd[x]) in [PooledDataArray{String,UInt8,1}], vars_2)) )
        rm_lst=Symbol[]
        if (length(vars_2) > 1) & (   length(getColswithType("num", dfd, vars_2 ) ) > 1  )
            stackdf = corDFD(dfd,vars_2)
            stackdf[:variable_pval] = [ sdf[sdf[:parameter].==c,:pval][1]   for c in stackdf[:variable]]
            stackdf[:vars_pval] = [ sdf[sdf[:parameter].==c,:pval][1]   for c in stackdf[:vars]] 
            stackdf[:most_Sig] = map((x,y) -> x < y ? "variable" : "vars" ,stackdf[:variable_pval],stackdf[:vars_pval])
     
            for row in eachrow(stackdf[(stackdf[:value].> 0.8) | (stackdf[:value].<-0.8),:])
                if row[:vars] == "group"
                    push!(rm_lst,row[:variable])
                elseif row[:variable] == "group"
                    push!(rm_lst,row[:vars])
                else
                    row[:most_Sig] == "variable" ? push!(rm_lst,row[:vars]) : push!(rm_lst,row[:variable])
                end
            end    
        end
        return rm_lst
    end
    println(upperModName*" : Correlation") 
    corrvars_x = getCorrVars(g4, dfd,vars_2  )
    vars_2 = setdiff(vars_2,corrvars_x)    
    (g5, finalvars) =  chkSigns(m, convert(Array{Symbol},vars_2), dfd, cfg)
    if length(finalvars) > 2    
        println( upperModName *" : Z & Vif") 
        f=genF(m[:y_var],finalvars)
        g2 = glm( f, dfd[convert(Array{Symbol},vcat(f.lhs,f.rhs.args[3:end]))]  , m[:dist], m[:lnk] )
        sdf=vifDF2(g2)
        z = sdf[abs(sdf[:zval]).<1.96,:parameter]
        v = sdf[ !DataArrays.isna(sdf[:vif])&(sdf[:vif].>15),:parameter]
        g2_x =intersect(z,v)
        finalvars = setdiff(finalvars,g2_x)
        println("vif_vars: ",g2_x)
    end
    finalvars = setdiff(finalvars,[:group])
    finalvars = convert(Array{Symbol},vcat(finalvars,[:group]) )
    return finalvars
end
     

factor_cols=vcat( [ cfg[:proscore], :group, :panid], cfg[:allrandoms] )
for c in setdiff(factor_cols,[:panid, cfg[:proscore]]) 
    if !( typeof(dfd[c]) in [  DataArray{String,1}  ]) 
        println("converting to Strings : ", c," of type : ",typeof(dfd[c]))
        dfd[c] = map(x->string(x),dfd[c])
        dfd[c] = convert(Array{String},dfd[c]) 
    end
end
poolit!(dfd,factor_cols)

iocc[:finalvars] = featureSelection(dfd[(dfd[:buyer_pos_p1].==1),:], expandM(iocc))
idolocc[:finalvars]   = featureSelection(dfd[(dfd[:buyer_pos_p1].==1),:], expandM(idolocc))
ipen[:finalvars]  = featureSelection(dfd, expandM(ipen))
#:core_based_statistical_areas, :latitude, :county_code

    
#function scaleTestControl(dfd::DataFrame) 
#        t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        f=(t/c)
#        println(f)
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] .* f
#end
#scaleTestControl(dfd)    
#function scaleTestControl(dfd::DataFrame) 
#        t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        f=(t/c)
#        println(f)
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] .* f
#        # OCC
#        o_t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_f=(o_t/o_c)
#        println(o_f)
#        #dfd[modelsDict[:occ][:y_var]] = Float64.(dfd[modelsDict[:occ][:y_var]])
#        #dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = Int64.(round.(dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f))
#end    
#scaleTestControl(dfd) 
    
    
delete!(iocc, :exclude_vars)
delete!(idolocc, :exclude_vars)
delete!(ipen, :exclude_vars)
modelsDict = Dict()
modelsDict[:occ] = iocc
modelsDict[:dolocc] = idolocc
modelsDict[:pen] = ipen
modelsDict[:factors] = unpool(dfd) 
cols = vcat([:panid, :group, :buyer_pos_p1],cfg[:random_campaigns])
cols = vcat(cols,iocc[:finalvars],[iocc[:y_var], iocc[:logvar]])
cols = vcat(cols,idolocc[:finalvars],[idolocc[:y_var], idolocc[:logvar]])
cols = unique(vcat(cols,ipen[:finalvars],[ipen[:y_var], ipen[:logvar]]) )
cols = unique(vcat(cols,[:buyer_pre_p1, :buyer_pos_p1, :trps_pos_p1,:trps_pre_p1,:dol_per_trip_pre_p1,:dol_per_trip_pos_p1, :prd_1_net_pr_pre, :prd_1_net_pr_pos]))
    
   

function scaleTestControl(dfd::DataFrame) 
        tst=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
        ctrl=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
        f=(tst/ctrl)
        println(f)
        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] .* f
        dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:dolocc][:logvarOrig]] = dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:dolocc][:logvarOrig]] .* f
        # OCC
        o_t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
        o_c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
        o_f=(o_t/o_c)
        println(o_f)
        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = Int64.(round.(dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f))
        dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:occ][:logvarOrig]] = Int64.(round.(dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:occ][:logvarOrig]] .* o_f))
end    
scaleTestControl(dfd) 
        
    
    
save_dfd(root, dfd[cols] )
save_modelsDict(root, modelsDict)
    
    

#function scaleTestControl(dfd::DataFrame) 
#        t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        f=(t/c)
#        println(f)
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] .* f
#        # OCC
#        o_t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_f=(o_t/o_c)
#        println(o_f)
#        #dfd[modelsDict[:occ][:y_var]] = Float64.(dfd[modelsDict[:occ][:y_var]])
#        #dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = Int64.(round.(dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f))
#end
#dfd = read_dfd(root)
#println(mean(dfd[(dfd[:group].==0)|(dfd[:group].=="0"),modelsDict[:dolocc][:y_var]] ))  
#println(mean(dfd[(dfd[:group].==0)|(dfd[:group].=="0"),modelsDict[:occ][:y_var]] ))  
#scaleTestControl(dfd)    
#save_dfd(root, dfd) 
    
    
#    dfd[:pre] = deepcopy(dfd[modelsDict[:dolocc][:logvarOrig]])
#    dfd[:pos] = deepcopy(dfd[modelsDict[:dolocc][:y_var]])
 
#    dfd[modelsDict[:dolocc][:logvarOrig]] = deepcopy(dfd[:pre])
#    dfd[modelsDict[:dolocc][:y_var]] = deepcopy(dfd[:pos])
"""
function t()
    Dict(:unadj_avg_cntrl_hh_pre => mean(dfd[ (dfd[:buyer_pre_p1].==1) & ((dfd[:group].=="0")|(dfd[:group].==0)), :dol_per_trip_pre_p1] ), 
         :unadj_avg_expsd_hh_pre => mean(dfd[ (dfd[:buyer_pre_p1].==1) & ((dfd[:group].=="1")|(dfd[:group].==1)), :dol_per_trip_pre_p1] ),
        :unadj_avg_cntrl_hh_pst => mean(dfd[ (dfd[:buyer_pos_p1].==1) & ((dfd[:group].=="0")|(dfd[:group].==0)), :dol_per_trip_pos_p1] ),
        :unadj_avg_expsd_hh_pst => mean(dfd[ (dfd[:buyer_pos_p1].==1) & ((dfd[:group].=="1")|(dfd[:group].==1)), :dol_per_trip_pos_p1] ))
end
""" 
#function scaleTestControl(dfd::DataFrame) 
#        tst=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        ctrl=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:dolocc][:logvarOrig]])
#        f=(tst/ctrl)
#        println(f)
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] = dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:dolocc][:y_var]] .* f
#        dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:dolocc][:logvarOrig]] = dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:dolocc][:logvarOrig]] .* f
#        # OCC
#        o_t=mean(dfd[((dfd[:group].=="1")|(dfd[:group].==1))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_c=mean(dfd[((dfd[:group].=="0")|(dfd[:group].==0))&(dfd[:buyer_pre_p1].==1),modelsDict[:occ][:logvarOrig]])
#        o_f=(o_t/o_c)
#        println(o_f)
#        dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] = Int64.(round.(dfd[(dfd[:group].=="0")|(dfd[:group].==0),modelsDict[:occ][:y_var]] .* o_f))
#        dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:occ][:logvarOrig]] = Int64.(round.(dfd[((dfd[:group].=="0")|(dfd[:group].==0)),modelsDict[:occ][:logvarOrig]] .* o_f))
#end    
#scaleTestControl(dfd) 
    
    
    
    
    
# ******************************* MODELS **********************************************************
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays
root = !isdefined(:root) ? pwd() : root
MP = !isdefined(:MP) ? false : MP
    
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
    
#dx=OrderedDict()  # defaults to empty - so that consolidate will work with seq or MP
function runModels(root::String, modelsDict::Dict)
    dfd = read_dfd(root) #dfd = readtable(mod_fname,header=true);
    poolit!(dfd,modelsDict[:factors])
    res = OrderedDict()
    for gm in [:occ,:dolocc,:pen] res[Symbol("glm_"*string(gm))] = runGlm(dfd,modelsDict[gm]) end
    q = Symbol[]
    cnt=0
    for (k,v) in OrderedDict(m=>modelsDict[m][:raneff] for m in [:pen,:dolocc,:occ]) 
         for r in v push!(q,k); push!(q,r); cnt=cnt+1 end 
    end
    x=reshape(q, (2,cnt))
    ml = permutedims(x, [2, 1])
    m=:empty
    for i in 1:length(ml[:,1]) 
        if m != ml[i,:][1]
            #gr=Symbol("glm_"*string(ml[i,:][1]))
            #res[gr] = runGlm(dfd,modelsDict[ml[i,:][1]]) #res[gr] = runGlm(mod_fname, jmod_fname, ml[i,:][1])
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
    
    
    
#function runModels( root::String, modelsDict::Dict)
#   @swarm :glm_ipen runGlm(root, :ipen)
#   sleep(25)
#   @swarm :glm_iocc runGlm(root, :iocc)
#   sleep(25)
#   @swarm :glm_idolocc runGlm(root, :idolocc)    
#   q = Symbol[]
#   cnt=0
#   for (k,v) in OrderedDict(m=>modelsDict[m][:raneff] for m in [:ipen,:idolocc,:iocc]) 
#        for r in v push!(q,k); push!(q,r); cnt=cnt+1 end 
#   end
#   x=reshape(q, (2,cnt))
#   ml = permutedims(x, [2, 1])
#   m=:empty
#   for i in 1:length(ml[:,1]) 
#       r=Symbol("glmm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2]))
#       @swarm r runGlmm(root, ml[i,:][1], ml[i,:][2]) 
#       covg=Symbol("covglm_"*string(ml[i,:][1])*"_"*string(ml[i,:][2])   )
#       @swarm covg runGlm(root, ml[i,:][1], ml[i,:][2]) 
#   end 
#   sleep(20)
#   while !iscompletew() println("Not Complete yet!"); sleep(5); end 
#   println("DONE!DONE!DONE!")
#end        
#runModels(root,modelsDict)                   
 
dx=(isdefined(:MP))&(MP)?(MPrunModels(root,modelsDict);OrderedDict()):runModels(root,modelsDict) 

        
function consolidateResults(modelsDict::Dict,dx::OrderedDict=OrderedDict())  
  sdf = DataFrame(parameter=String[], coef=Int64[], stderr=Int64[], zval=Int64[], pval=Int64[], model=String[], ranef=String[], modelType=String[])
    for m in [:occ, :dolocc, :pen]
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
dfx = consolidateResults(modelsDict,dx)  
save_dfx(root,dfx)      #writetable(root*"/campaign.csv", dfx)                      

    
# ************************************* SCORING*************************************************************************************
  
#ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg" 
#@everywhere insert!(Base.LOAD_CACHE_PATH, 1, "/mapr/mapr04p/analytics0001/analytic_users/jpkg/lib/v0.5")
#@everywhere pop!(Base.LOAD_CACHE_PATH)
#using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
#@everywhere using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
#root = !isdefined(:root) ? pwd() : root
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib, NamedArrays, JuMP, NLopt, NLsolve
root = !isdefined(:root) ? pwd() : root
MP = !isdefined(:MP) ? false : MP  
            
modelsDict = read_modelsDict(root)  
#modelsDict[:occ]=modelsDict[:iocc]; delete!(modelsDict,:iocc)
#modelsDict[:dolocc]=modelsDict[:idolocc]; delete!(modelsDict,:idolocc)
#modelsDict[:pen]=modelsDict[:ipen]; delete!(modelsDict,:ipen)
dfx = read_dfx(root) 
dfd = read_dfd(root)  
dfx[:unadj_mean_score0] = 0.0
dfx[:unadj_mean_score1] = 0.0
dfx[:adj_mean_score0 ] = 0.0
dfx[:adj_mean_score1 ] = 0.0
dfx[:adj_mean_score0_preweight ] = 0.0
dfx[:adj_mean_score1_preweight ] = 0.0
dfx[:weight0 ] = 1.0
dfx[:weight1 ] = 1.0 
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
Rlist() = by(dfx[(dfx[:modelType].=="GLMM"),[:ranef,:parameter]], [:ranef,:parameter], df -> DataFrame(N = size(df,1)))[[:ranef,:parameter]]

#unique(rlist[rlist[:ranef].=="targeting",:parameter]
#unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].=="targeting"),:parameter])
#hhcnts[(hhcnts[:ranef].=="targeting"),:parameter] 
        
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


function genRawDataMeans(dfx::DataFrame) 
    for m in mlist   # ************************* FIXED **************************
        ex_Buyer_Pos_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pos_p1].==1) &" : ""
        ex_Buyer_Pre_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pre_p1].==1) &" : ""
        precol=string(m[:logvarOrig])
        postcol=string(m[:y_var])
        model=string(m[:modelName])
        unadj_avg_cntrl_hh_pre_eq = "mean(dfd[ $ex_Buyer_Pre_P1 (dfd[:group].==0), :$precol] )"  
        unadj_avg_expsd_hh_pre_eq = "mean(dfd[ $ex_Buyer_Pre_P1 (dfd[:group].==1), :$precol] )"
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
        ex_Buyer_Pos_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pos_p1].==1) &" : ""
        ex_Buyer_Pre_P1 = m[:Buyer_Pos_P1_is1] ? "(dfd[:buyer_pre_p1].==1) &" : ""
        for r in m[:raneff]
            rs=string(r)
            exposed=true # Hack till we figure this out
            for l in dfx[(dfx[:modelType].=="GLMM")&(dfx[:parameter].!="none")&(dfx[:ranef].==rs)&(dfx[:model].==model),:parameter]
                unadj_avg_expsd_hh_pre_eq= "mean(dfd[ $ex_Buyer_Pre_P1  (dfd[:group] .== 1) & (dfd[:$rs] .== \"$l\"), :$precol])"
                unadj_avg_expsd_hh_pst_eq= "mean(dfd[ $ex_Buyer_Pos_P1  (dfd[:group] .== 1) & (dfd[:$rs] .== \"$l\"), :$postcol])"
                unadj_avg_cntrl_hh_pre_eq= "mean(dfd[ $ex_Buyer_Pre_P1  (dfd[:group] .== 0) & (dfd[:$rs] .== \"$l\"), :$precol])"
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
    if isfile(root*"/hhcounts.csv")
        rlist = Rlist()
        hhcnts = readtable(root*"/hhcounts.csv"); rename!(hhcnts,:class,:ranef); rename!(hhcnts,:level,:parameter); lowercase(hhcnts); hhcnts[:ranef] = lowercase(hhcnts[:ranef])
        hhcnts=join(rlist,hhcnts,on=[:ranef,:parameter]) #!!!!!!            
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
    

dfx[:adj_mean_score0_pw ] = dfx[:adj_mean_score0]
dfx[:adj_mean_score1_pw ] = dfx[:adj_mean_score1]

function weightDFX()
    if isfile(root*"/hhcounts.csv")
        dfw = genWeights() 
        for row in eachrow(dfw)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),:weight1] = row[:adj_mean_score1_Factor]
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),:weight0] = row[:adj_mean_score0_Factor]
        end
        dfx[:adj_mean_score0] = dfx[:adj_mean_score0_pw] .* dfx[:weight0]
        dfx[:adj_mean_score1] = dfx[:adj_mean_score1_pw] .* dfx[:weight1]
#dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ"),[:model,:ranef,:parameter,:unadj_avg_cntrl_hh_pre,:unadj_avg_cntrl_hh_pst]]
        #for col in [:unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst]
        #    dfx[col] = dfx[col] .* dfx[:weight0]
        #    dfx[col] = dfx[col] .* dfx[:weight1]   
        #end
        for col in [:unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst]
            for row in eachrow(dfw)
                dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),col] = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==row[:model])&(dfx[:ranef].==row[:ranef]),col] .* row[Symbol(string(col)*"_Factor")] 
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

function sclaeTestCtrl()   # VIJOY's SCALING
    t = dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")&(dfx[:model].=="dolocc"),:unadj_avg_expsd_hh_pre][1]
    c = dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")&(dfx[:model].=="dolocc"),:unadj_avg_cntrl_hh_pre][1]
    f = t/c 
    dfx[:scaleF]=1.0
    dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:scaleF] = f
    dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc"),:scaleF] = f
    if c > t
        for c in [:adj_mean_score0,:adj_mean_score0_preweight,:unadj_avg_cntrl_hh_pst]
            dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),c] = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),c] .* f
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc"),c] = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc"),c] .* f 
        end
    end
end
#sclaeTestCtrl()


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
                Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*zscore)))*(Mt/M))   )  +  ## ------ Lower Bound --------
                         (   (mean_score0*exp((B-(SE*zscore)))*(Mc/M))    -   (mean_score0*(Mc/M))    )
                row[Symbol(zscore_key*"_lb")] = ( Lb_pre/mean_score0 ) * 100
                Ub_pre =  (     ( mean_score1*(Mt/M) )   -   ( mean_score1*exp(-(B+(SE*zscore)))*(Mt/M))   )  +  ## ------ Upper Bound -------
                          (     ( mean_score0*exp((B+(SE*zscore)))*(Mc/M))  - (mean_score0*(Mc/M) )   )
                row[Symbol(zscore_key*"_ub")] = ( Ub_pre/mean_score0 ) * 100   
                println("RAND CI $zscore_key ($zscore) LB:",row[Symbol(zscore_key*"_lb")]," ~~ UB : ", row[Symbol(zscore_key*"_ub")])
            end            
        end
    end
end
#ConfIntrvl()
#dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100

# **********************************************************************************************
# ***************************************** DOLHH **********************************************
# **********************************************************************************************

#function blank(dfx::DataFrame,d::Dict=Dict())
#    arr = Any[]      
#    for (n, v) in eachcol(dfx)
#        if haskey(d,n)
#            arr = vcat(arr, [d[n]])  
#        else
#            if isFloat(v)
#                arr = vcat(arr, [0.0])
#            elseif isInt(v)
#                arr = vcat(arr, [0])
#            elseif isString(v)
#                arr = vcat(arr, [""])
#            elseif isBool(v)
#                        arr = vcat(arr, [false])
#            elseif isSymbol(v)
#                        arr = vcat(arr, [:empty])
#            else
#                println("ERROR : datatype Not Found!!!!!!!!")
#            end
#        end     
#    end 
#    return  arr
#end            


        
function genDHHMeans(dfx::DataFrame)      
    dfh=similar(dfx,0)
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
            
    d=Dict{Any,Any}(:parameter=>"group", :model=>"dolhh", :ranef=>"ranef", :modelType=>"GLM",
         :adj_mean_score0=>adj_mean_score0 ,:adj_mean_score1=>adj_mean_score1, 
         :unadj_avg_expsd_hh_pre=>unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst=>unadj_avg_expsd_hh_pst, 
         :unadj_avg_cntrl_hh_pre=>unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst=>unadj_avg_cntrl_hh_pst,
         :weight0=>1.0, :weight1=>1.0
        )
    d=merge(d, getCnts(dfd))        
    push!(dfh,blank(dfx, d) )
            
    #pre=["group",0.0,0.0,0.0,0.0,"dolhh","ranef","GLM"]
    #vars=[0.0,0.0,adj_mean_score0,adj_mean_score1,unadj_avg_expsd_hh_pre,unadj_avg_expsd_hh_pst,unadj_avg_cntrl_hh_pre,unadj_avg_cntrl_hh_pst]
    #cis = zeros(10)  # cis and pvals [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    #pvals=[0.0,0.0]
    #cnts=convert(Array{Int64},collect(values(getCnts(dfd))))
    #push!(dfx, filler(dfx,vcat(pre,vars,cis,cnts)))
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
                unadj_avg_cntrl_hh_pre = dfh[(dfh[:modelType].=="GLM")&(dfh[:model].=="dolhh") ,:unadj_avg_cntrl_hh_pre][1]
                unadj_avg_cntrl_hh_pst = dfh[(dfh[:modelType].=="GLM")&(dfh[:model].=="dolhh") ,:unadj_avg_cntrl_hh_pst][1]     
            end   
             
                    
            d=Dict{Any,Any}(:parameter=>level, :model=>"dolhh", :ranef=>ranef, :modelType=>"GLMM",
                            :adj_mean_score0=>adj_mean_score0 ,:adj_mean_score1=>adj_mean_score1, 
                            :unadj_avg_expsd_hh_pre=>unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst=>unadj_avg_expsd_hh_pst, 
                            :unadj_avg_cntrl_hh_pre=>unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst=>unadj_avg_cntrl_hh_pst,
                            :weight0=>1.0, :weight1=>1.0
                            )
            d=merge(d, getCnts(dfd,ranef,level))        
            push!(dfh,blank(dfx, d) )
            #println(d)
                    
            #pre=[level,0.0,0.0,0.0,0.0,"dolhh",ranef,"GLMM"]
            #vars=[0.0,0.0,adj_mean_score0,adj_mean_score1,unadj_avg_expsd_hh_pre,unadj_avg_expsd_hh_pst,unadj_avg_cntrl_hh_pre,unadj_avg_cntrl_hh_pst ]
            #cis=zeros(10) #[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            #cnts=convert(Array{Int64},collect(values(getCnts(dfd,ranef,level))))
            #push!(dfx, filler(dfx,vcat(pre,vars,cis,cnts)))
        end
    end
    return vcat(dfx,dfh)
end
dfx=genDHHMeans(dfx)


        
        
# -------------------------------------------------------------------------------------------------------------
 

@everywhere function OPTPValue(iDict::OrderedDict)
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
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=2))
    #@variable(m, Bocc <= B1)
    #@variable(m, Bdolocc <= B2)
    #@variable(m, Bpen <= B3)
    if B1 > 0 @variable(m, Bocc <= B1) else @variable(m, Bocc >= B1) end
    if B2 > 0 @variable(m, Bdolocc <= B2) else @variable(m, Bdolocc >= B2) end
    if B3 > 0 @variable(m, Bpen <= B3) else @variable(m, Bpen >= B3) end
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
    dout[:pval] = pvalue
    dout[:onetail_pval] = one_tail
    dout[:twotail_pval] = two_tail
    println("z-value: ", string(zvalue)," --> p-value: ",string(two_tail))
    return dout           
end
            


@everywhere function OPT_LB(iDict::OrderedDict, zscore::Float64, iAccuracy::Float64=0.000000001,minval::Float64=0.0)
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
    #println(M,"\n",Mt,"\n",Mc,"\n",N,"\n",Nt,"\n",Nc,"\n",B1,"\n",B2,"\n",B3,"\n",SE1,"\n",SE2,"\n",SE3,"\n",o_SE0,"\n",y_SE0,"\n",p_SE0,"\n",SEsq,"\n",
    #        o_B0,"\n",y_B0,"\n",p_B0,"\n",o_mean_score0,"\n",o_mean_score1,"\n",y_mean_score0,"\n",y_mean_score1,"\n",p_mean_score0,"\n",p_mean_score1) 
    pen_t = p_mean_score1*(Nt/N)
    pen_c = p_mean_score0*(Nc/N)
    occ_t = o_mean_score1*(Mt/M)
    occ_c = o_mean_score0*(Mc/M)
    dolocc_t = y_mean_score1*(Mt/M)
    dolocc_c = y_mean_score0*(Mc/M)
    Bsum=(B1+B2+B3)-(o_B0+y_B0+p_B0)
    ztot = Bsum-(zscore*SEsq)
    ######CONFIDENCE INTERVAL - LB ########        
    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=2))
    @variable(m, Bocc <= (B1-o_B0))
    @variable(m, Bdolocc <= (B2-y_B0))
    @variable(m, Bpen <= (B3-p_B0))
    #@NLexpression(m, expr1 ,(( ((pen_t + pen_c*exp(Bpen))*(occ_t+occ_c*exp(Bocc))*(dolocc_t+dolocc_c*exp(Bdolocc)))
    #                           /
    #                          ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
    #                          )-1)
    #                 )
    
    @NLexpression(m, expr1 ,(  ( ((pen_t + pen_c*exp(Bpen))*(occ_t+occ_c*exp(Bocc))*(dolocc_t+dolocc_c*exp(Bdolocc)))
                                   -
                                  ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
                                ) / ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
                              
                            )
                     )
    
    #@NLobjective(m, Min, (( ((pen_t + pen_c*exp(Bpen))*(occ_t+occ_c*exp(Bocc))*(dolocc_t+dolocc_c*exp(Bdolocc)))
    #                                /
    #                               ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
    #                             )-1)
    #                    )
    @NLobjective(m, Min, expr1)
    @constraint(m, (0.000000000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAccuracy)) 
    @NLconstraint(m, constr1, iDict[:minval] <= expr1)
    status = solve(m)
    mval=getobjectivevalue(m)
    println(mval)
    println(status)
    #return status == :Optimal ? mval : -Inf
    return (status != :Optimal)|(mval < -0.85 ) ? -Inf : mval
end
#OPT_LB(collectModels(dfx,"GLMM","campaign_name","Pure_Life_Promise"), 0.84, float("0.000000001"))
#zscore=0.84
#iAccuracy=float("0.000000001")
#OPT_LB(collectModels(dfx,"GLMM","site","Evite.com"), 1.65, float("0.000000001"))



@everywhere function OPT_UB(iDict::OrderedDict, zscore::Float64, iAccuracy::Float64=0.000000001,maxval::Float64=0.0)
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
    pen_t = p_mean_score1*(Nt/N)
    pen_c = p_mean_score0*(Nc/N)
    occ_t = o_mean_score1*(Mt/M)
    occ_c = o_mean_score0*(Mc/M)
    dolocc_t = y_mean_score1*(Mt/M)
    dolocc_c = y_mean_score0*(Mc/M)
    #Bsum=B1+B2+B3 + o_B0+y_B0+p_B0
    Bsum=(B1+B2+B3)-(o_B0+y_B0+p_B0)
    ######CONFIDENCE INTERVAL - UB ########
    ztot = Bsum+(zscore*SEsq)
    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=2))
    @variable(m, Bocc >= (B1-o_B0))
    @variable(m, Bdolocc >= (B2-y_B0))
    @variable(m, Bpen >= (B3-p_B0))
    @NLexpression(m, expr1 ,(( ((pen_t + pen_c*exp(Bpen))*(occ_t+occ_c*exp(Bocc))*(dolocc_t+dolocc_c*exp(Bdolocc)))
                               /
                              ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
                              )-1)
                     )
    
    #@NLobjective(m, Max, (( ((pen_t + pen_c*exp(Bpen))*(occ_t+occ_c*exp(Bocc))*(dolocc_t+dolocc_c*exp(Bdolocc)))
    #                                /
    #                               ((pen_t*exp(-Bpen)+pen_c)*(occ_t*exp(-Bocc)+occ_c)*(dolocc_t*exp(-Bdolocc)+dolocc_c))
    #                             )-1)
    #                    )
    @NLobjective(m, Max, expr1)
    @constraint(m, (0.0000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAccuracy))     
    @NLconstraint(m, constr1, expr1 <= iDict[:maxval])
    status = solve(m)
    mval=getobjectivevalue(m)
    println(mval)
    #return return status == :Optimal ? mval : -Inf  #mval_out 
    return (status != :Optimal)|(mval > 0.85 ) ? -Inf : mval
end


@everywhere function OPTCIs(iDict::OrderedDict)
    OPTPValue(iDict)
    #AccArr= [ "0.000000001","0.00000001","0.0000001","0.000001","0.00001","0.0001","0.001","0.01"]
    AccArr= [ "0.000000001","0.00000001","0.0000001","0.000001","0.00001","0.0001","0.001","0.01"]
    for (zscore_key,zscore) in  ZDict          
        pref="LB "*zscore_key[1:10]*" ("*iDict[:metakey]"):= "
        preflen=length(pref)
        dkey=Symbol(zscore_key*"_lb")
        for iAcc in AccArr
            print(pref," - "*iAcc*", ",iDict)       ;pref=lpad("", preflen, " ")
            iDict[dkey] = OPT_LB(iDict, zscore, float(iAcc))
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
           iDict[dkey] = OPT_UB(iDict, zscore, float(iAcc))
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





function collectModels(dfx::DataFrame,modelType::String,ranef::String="",level::String="")
    mDict = OrderedDict()
    if modelType=="GLM"
        o=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:]
        y=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:]
        p=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:]
        dhh=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolhh")&(dfx[:parameter].=="group"),:]
        mDict[:metakey] = "Total_Campaign" 
    else
        o=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="occ")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        y=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="dolocc")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        p=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="pen")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        dhh=dfx[(dfx[:modelType].==modelType)&(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:]
        mDict[:metakey] = ranef*"~"*level 
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
    mDict[:o_weight0]= o[:weight0][1]== 0.0 ? 1.0 : o[:weight0][1]
    mDict[:o_weight1]= o[:weight1][1]== 0.0 ? 1.0 : o[:weight1][1]
    mDict[:y_weight0]= y[:weight0][1]== 0.0 ? 1.0 : y[:weight0][1]
    mDict[:y_weight1]= y[:weight1][1]== 0.0 ? 1.0 : y[:weight1][1]
    mDict[:p_weight0]= p[:weight0][1]== 0.0 ? 1.0 : p[:weight0][1]
    mDict[:p_weight1]= p[:weight1][1]== 0.0 ? 1.0 : p[:weight1][1]
    
    # SCALING NBNBNB
#    mDict[:y_weight0] = mDict[:y_weight0] * y[:scaleF][1]
    # SCALING END
            
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
    mDict[:o_mean_score0_pw]=o[:adj_mean_score0_pw][1]  # == 0.0 ? 1.0 : o[:adj_mean_score0_preweight][1]
    mDict[:o_mean_score1_pw]=o[:adj_mean_score1_pw][1]  # == 0.0 ? 1.0 : o[:adj_mean_score1_preweight][1]
    mDict[:y_mean_score0_pw]=y[:adj_mean_score0_pw][1]  # == 0.0 ? 1.0 : y[:adj_mean_score0_preweight][1]
    mDict[:y_mean_score1_pw]=y[:adj_mean_score1_pw][1]  # == 0.0 ? 1.0 : y[:adj_mean_score1_preweight][1]
    mDict[:p_mean_score0_pw]=p[:adj_mean_score0_pw][1]  # == 0.0 ? 1.0 : p[:adj_mean_score0_preweight][1]
    mDict[:p_mean_score1_pw]=p[:adj_mean_score1_pw][1]  # == 0.0 ? 1.0 : p[:adj_mean_score1_preweight][1]
            
    mDict[:SEsq]=sqrt(mDict[:SE1]^2+mDict[:SE2]^2+mDict[:SE3]^2)
    mDict[:Bsum]=mDict[:B1]+mDict[:B2]+mDict[:B3]
    mDict[:SD]=mDict[:SEsq]*(sqrt(mDict[:N]))
    mDict[:adj_dod] = ((dhh[:adj_mean_score1][1]-dhh[:adj_mean_score0][1])/dhh[:adj_mean_score0][1])*100 
    mDict[:minval] = (mDict[:adj_dod]-(mDict[:SD]*4))/100
    mDict[:maxval] = (mDict[:adj_dod]+(mDict[:SD]*4))/100
    mDict[:minval] = mDict[:minval] < -0.6 ? -0.6 : mDict[:minval]
    mDict[:maxval] = mDict[:maxval] > 0.6 ? 0.6 : mDict[:maxval]
    return mDict
end


function ConfidenceIntervals()
    mDict=collectModels(dfx,"GLM")
    println("Running CI for : Total") 
    #OPTPValue(mDict)
    OPTCIs(mDict)
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:pval] = mDict[:pval]
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:coef] = mDict[:SEsq]
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:stderr] = mDict[:Bsum]
    for k in keys(ZDict)
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
    end
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")  
            #OPTPValue(mDict)
            OPTCIs(mDict)
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:pval] = mDict[:pval]
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_coef] = mDict[:Bsum]
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:adj_stderr] = mDict[:SEsq]
            println("!!!!!!!!!!!!",ranef,"!!!!!!!!!!!!",level,"!!!!!!!!!!!!!!!!!!!!!! ::: ",mDict[:Bsum]," ~~~ ",mDict[:SEsq] )
            for k in keys(ZDict)
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_lb")] = mDict[Symbol(k*"_lb")]*100
                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_ub")] = mDict[Symbol(k*"_ub")]*100
            end
        end
    end
end
#ConfidenceIntervals()

        
function PVals()
    mDict=collectModels(dfx,"GLM")
    println("Running CI for : Total") 
    OPTPValue(mDict)
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:pval] = mDict[:pval]
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")  
            OPTPValue(mDict)
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:pval] = mDict[:pval]
         end
    end
end
PVals()


#function MPConfidenceIntervals(dfx::DataFrame)    
#    mDict=collectModels(dfx,"GLM")
#    @swarm :total OPTCIs(mDict)
#    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
#        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
#                rDict=collectModels(dfx,"GLMM",ranef,level)
#                k = Symbol(ranef*"~"*level)
#                @swarm k OPTCIs(rDict)
#        end
#    end
#    sleep(20)
#    while !iscompletew() println("Not Complete yet!"); sleep(5); end 
#    zDict = takew(:total)
#    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:pval] = zDict[:pval]
#    for k in keys(ZDict)
#        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_lb")] = zDict[Symbol(k*"_lb")]*100
#        dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),Symbol(k*"_ub")] = zDict[Symbol(k*"_ub")]*100
#    end
#    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
#        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
#                xDict = takew( Symbol(ranef*"~"*level) ) 
#                dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:pval] = xDict[:pval]
#                for k in keys(ZDict)
#                    dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_lb")] = xDict[Symbol(k*"_lb")]*100
#                    dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),Symbol(k*"_ub")] = xDict[Symbol(k*"_ub")]*100
#                end 
#                println("*~*~*~*~*~*~~* : ",dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_80_pct_intrvl_ub][1]," ~~~~~~~~ ",xDict[:twotail_80_pct_intrvl_ub]*100)
#                println(xDict)
#        end
#    end    
#            return dfx
#end  
#MPConfidenceIntervals(dfx)

        #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh"),[:ranef,:parameter,:twotail_80_pct_intrvl_lb,:twotail_80_pct_intrvl_ub]]
    
#isdefined(:MP)&&(MP) ? MPConfidenceIntervals(dfx) : ConfidenceIntervals()

#https://github.com/JuliaStats/Distributions.jl/blob/master/doc/source/starting.rst
function Glivenko_Cantelli_dolhh(d::OrderedDict)
    pen_t = (d[:p_mean_score1]*(d[:Nt]/d[:N]))
    pen_c = (d[:p_mean_score0]*(d[:Nc]/d[:N]))
    occ_t = (d[:o_mean_score1]*(d[:Mt]/d[:M]))
    occ_c = (d[:o_mean_score0]*(d[:Mc]/d[:M]))
    dolocc_t = (d[:y_mean_score1]*(d[:Mt]/d[:M]))
    dolocc_c = (d[:y_mean_score0]*(d[:Mc]/d[:M]))
    dfm = DataFrame(Bocc=rand(Normal(d[:B1], d[:SE1]),1000000),Bdolocc=rand(Normal(d[:B2],d[:SE2]),1000000),Bpen=rand(Normal(d[:B3], d[:SE3]),1000000))
    dfm[:expd] = ( (pen_t+pen_c.*exp(dfm[:Bpen])).*(occ_t+occ_c.*exp(dfm[:Bocc])).*(dolocc_t+dolocc_c.*exp(dfm[:Bdolocc]))  )
    dfm[:ctrl] = ((pen_t.*exp(-dfm[:Bpen])+pen_c).*(occ_t.*exp(-dfm[:Bocc])+occ_c).*(dolocc_t.*exp(-dfm[:Bdolocc])+dolocc_c))
    dfm[:diff] = dfm[:expd] .- dfm[:ctrl]
    dfm[:adj_dod] = dfm[:diff] ./ dfm[:ctrl]
    quantile(dfm[:adj_dod], [0.05,0.95,0.10,0.90])
end
function GC_CI_dolhh()
    dfx[:GC_twotail_90_pct_intrvl_ub] = 0.0
    dfx[:GC_twotail_90_pct_intrvl_lb] = 0.0
    dfx[:GC_twotail_80_pct_intrvl_ub] = 0.0
    dfx[:GC_twotail_80_pct_intrvl_lb] = 0.0
    mDict=collectModels(dfx,"GLM")
    println("Running CI for : Total")
    mDict[:metakey] = "Total Campaign"        
    Glivenko_Cantelli_dolhh(mDict)
    ci = Glivenko_Cantelli_dolhh(mDict)
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_90_pct_intrvl_lb] = ci[1]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_90_pct_intrvl_ub] = ci[2]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_80_pct_intrvl_lb] = ci[3]*100
    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_80_pct_intrvl_ub] = ci[4]*100
    println(ci[1],":",ci[3],mDict[:adj_dod]," ~~ ",ci[4],":",ci[2])
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")
            mDict[:metakey] = ranef*"~"*level        
            ci = Glivenko_Cantelli_dolhh(mDict)
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_90_pct_intrvl_lb] = ci[1]*100
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_90_pct_intrvl_ub] = ci[2]*100
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_80_pct_intrvl_lb] = ci[3]*100
            dfx[(dfx[:model].=="dolhh")&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_80_pct_intrvl_ub] = ci[4]*100
            println(ci[1],":",ci[3],mDict[:adj_dod]," ~~ ",ci[4],":",ci[2])
        end
    end
end
#GC_CI_dolhh()
#https://github.com/JuliaStats/Distributions.jl/blob/master/doc/source/starting.rst
function Glivenko_Cantelli_OLD(d::OrderedDict)
    pen_t = (d[:p_mean_score1_pw]*(d[:Nt]/d[:N]))
    pen_c = (d[:p_mean_score0_pw]*(d[:Nc]/d[:N]))
    occ_t = (d[:o_mean_score1_pw]*(d[:Mt]/d[:M]))
    occ_c = (d[:o_mean_score0_pw]*(d[:Mc]/d[:M]))
    dolocc_t = (d[:y_mean_score1_pw]*(d[:Mt]/d[:M]))
    dolocc_c = (d[:y_mean_score0_pw]*(d[:Mc]/d[:M]))         
    dfm = DataFrame(Bocc=rand(Normal(d[:B1], d[:SE1]),1000000),Bdolocc=rand(Normal(d[:B2],d[:SE2]),1000000),Bpen=rand(Normal(d[:B3], d[:SE3]),1000000))
    dfm[:p0] = (pen_t.*exp(-dfm[:Bpen])+pen_c).*d[:p_weight0]
    dfm[:p1] = (pen_t+pen_c.*exp(dfm[:Bpen])).*d[:p_weight1]   
    dfm[:o0] = (occ_t.*exp(-dfm[:Bocc])+occ_c).*d[:o_weight0]
    dfm[:o1] = (occ_t+occ_c.*exp(dfm[:Bocc])).*d[:o_weight1]
    dfm[:y0] = (dolocc_t.*exp(-dfm[:Bdolocc])+dolocc_c).*d[:y_weight0]
    dfm[:y1] = (dolocc_t+dolocc_c.*exp(dfm[:Bdolocc])).*d[:y_weight1]         
    dfm[:expd] = (dfm[:p1].*dfm[:o1].*dfm[:y1])
    dfm[:ctrl] = (dfm[:p0].*dfm[:o0].*dfm[:y0])
    dfm[:adj_dod] = (dfm[:expd] .- dfm[:ctrl]) ./ dfm[:ctrl]
    dfm[:p_adj_dod] = (dfm[:p1] .- dfm[:p0]) ./ dfm[:p0]
    dfm[:o_adj_dod] = (dfm[:o1] .- dfm[:o0]) ./ dfm[:o0]
    dfm[:y_adj_dod] = (dfm[:y1] .- dfm[:y0]) ./ dfm[:y0]
    bounds_arr=[0.05,0.95, 0.10,0.90, 0.10,0.90, 0.20, 0.80]
    for (k,v) in OrderedDict(:dolhh=>:adj_dod,:pen=>:p_adj_dod, :occ=>:o_adj_dod, :dolocc=>:y_adj_dod)      
        q = quantile(dfm[v], bounds_arr)
        d[k] = OrderedDict(:lb90=>q[1]*100,:ub90=>q[2]*100,:lb80=>q[3]*100,:ub80=>q[4]*100,:lb901=>q[5]*100,:ub901=>q[6]*100,:lb801=>q[7]*100,:ub801=>q[8]*100)        
    end
    # ---------- tst ------------
    dfm[:p_mean_score1_pw] = d[:p_mean_score1_pw]
    dfm[:p_mean_score0_pw] = d[:p_mean_score0_pw]
    dfm[:o_mean_score1_pw] = d[:o_mean_score1_pw]
    dfm[:o_mean_score0_pw] = d[:o_mean_score0_pw]
    dfm[:y_mean_score1_pw] = d[:y_mean_score1_pw]
    dfm[:y_mean_score0_pw] = d[:y_mean_score0_pw]
    dfm[:N] = d[:N]
    dfm[:Nt] = d[:Nt]
    dfm[:Nc] = d[:Nc]
    dfm[:M] = d[:M]
    dfm[:Mt] = d[:Mt]
    dfm[:Mc] = d[:Mc]
    dfm[:p_weight0] = d[:p_weight0]
    dfm[:p_weight1] = d[:p_weight1]
    dfm[:o_weight0] = d[:o_weight0]
    dfm[:o_weight1] = d[:o_weight1]
    dfm[:y_weight0] = d[:y_weight0]
    dfm[:y_weight1] = d[:y_weight1]
    #writetable(root*"/dfm.csv", dfm)
    # ------------ tst end -------------- 
end

function GC_ConfidenceIntervals_OLD()
    dfx[:GC_twotail_90_pct_intrvl_ub] = 0.0
    dfx[:GC_twotail_90_pct_intrvl_lb] = 0.0
    dfx[:GC_twotail_80_pct_intrvl_ub] = 0.0
    dfx[:GC_twotail_80_pct_intrvl_lb] = 0.0
    dfx[:GC_onetail_90_pct_intrvl_ub] = 0.0
    dfx[:GC_onetail_90_pct_intrvl_lb] = 0.0
    dfx[:GC_onetail_80_pct_intrvl_ub] = 0.0
    dfx[:GC_onetail_80_pct_intrvl_lb] = 0.0
    mDict=collectModels(dfx,"GLM")
            d=mDict
            ctl=(d[:o_mean_score1_pw]*(d[:Mt]/d[:M]))*exp(-d[:B1])+(d[:o_mean_score0_pw]*(d[:Mc]/d[:M]))
            tst=(d[:o_mean_score1_pw]*(d[:Mt]/d[:M]))+(d[:o_mean_score0_pw]*(d[:Mc]/d[:M]))*exp( d[:B1]  )   
            println("CTRL :::$ctl := ",d[:o_mean_score0_pw],",    TEST ::: $tst := ",d[:o_mean_score1_pw])
                    
            
            
    println("Running CI for : Total")   
    Glivenko_Cantelli_OLD(mDict)
    for m in [:dolhh,:occ,:dolocc,:pen]
        println(m)
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_90_pct_intrvl_lb] = mDict[m][:lb90]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_90_pct_intrvl_ub] = mDict[m][:ub90]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_80_pct_intrvl_lb] = mDict[m][:lb80]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_twotail_80_pct_intrvl_ub] = mDict[m][:ub80]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_onetail_90_pct_intrvl_lb] = mDict[m][:lb901]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_onetail_90_pct_intrvl_ub] = mDict[m][:ub901]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_onetail_80_pct_intrvl_lb] = mDict[m][:lb801]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:GC_onetail_80_pct_intrvl_ub] = mDict[m][:ub801]
    end    
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")
            mDict[:metakey] = ranef*"~"*level        
            Glivenko_Cantelli_OLD(mDict)
            for m in [:dolhh,:occ,:dolocc,:pen]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_90_pct_intrvl_lb] = mDict[m][:lb90]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_90_pct_intrvl_ub] = mDict[m][:ub90]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_80_pct_intrvl_lb] = mDict[m][:lb80]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_twotail_80_pct_intrvl_ub] = mDict[m][:ub80]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_onetail_90_pct_intrvl_lb] = mDict[m][:lb901]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_onetail_90_pct_intrvl_ub] = mDict[m][:ub901]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_onetail_80_pct_intrvl_lb] = mDict[m][:lb801]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:GC_onetail_80_pct_intrvl_ub] = mDict[m][:ub801]
            end         
        end
    end
end
#GC_ConfidenceIntervals_OLD()
 
        
        

function Glivenko_Cantelli(d::OrderedDict)
    nsamp = 100000
    pen_t = (d[:p_mean_score1_pw]*(d[:Nt]/d[:N]))
    pen_c = (d[:p_mean_score0_pw]*(d[:Nc]/d[:N]))
    occ_t = (d[:o_mean_score1_pw]*(d[:Mt]/d[:M]))
    occ_c = (d[:o_mean_score0_pw]*(d[:Mc]/d[:M]))
    dolocc_t = (d[:y_mean_score1_pw]*(d[:Mt]/d[:M]))
    dolocc_c = (d[:y_mean_score0_pw]*(d[:Mc]/d[:M])) 
    dfm = DataFrame(Bocc=rand(Normal(d[:B1], d[:SE1]),nsamp),Bdolocc=rand(Normal(d[:B2],d[:SE2]),nsamp),Bpen=rand(Normal(d[:B3], d[:SE3]),nsamp))
    dfm[:p0] = (pen_t.*exp(-dfm[:Bpen])+pen_c)
    dfm[:p1] = (pen_t+pen_c.*exp(dfm[:Bpen]))   
    dfm[:o0] = (occ_t.*exp(-dfm[:Bocc])+occ_c)
    dfm[:o1] = (occ_t+occ_c.*exp(dfm[:Bocc]))
    dfm[:y0] = (dolocc_t.*exp(-dfm[:Bdolocc])+dolocc_c)
    dfm[:y1] = (dolocc_t+dolocc_c.*exp(dfm[:Bdolocc]))         
    dfm[:expd] = (dfm[:p1].*dfm[:o1].*dfm[:y1])
    dfm[:ctrl] = (dfm[:p0].*dfm[:o0].*dfm[:y0])
    dfm[:adj_dod] = ((dfm[:expd] .- dfm[:ctrl]) ./ dfm[:ctrl] ) #*100
    dfm[:p_adj_dod] = ((dfm[:p1] .- dfm[:p0]) ./ dfm[:p0] ) #*100
    dfm[:o_adj_dod] = ((dfm[:o1] .- dfm[:o0]) ./ dfm[:o0]) #*100
    dfm[:y_adj_dod] = ((dfm[:y1] .- dfm[:y0]) ./ dfm[:y0]) #*100
    pdod = ((pen_t+pen_c.*exp(d[:B3]))-(pen_t.*exp(-d[:B3])+pen_c)) / (pen_t.*exp(-d[:B3])+pen_c)  
    pdod_w = ( d[:p_mean_score1_pw] * d[:p_weight1] - d[:p_mean_score0_pw] * d[:p_weight0] ) /  (d[:p_mean_score0_pw] * d[:p_weight0])
    pdiff = pdod_w-pdod
    odod = (((occ_t+occ_c.*exp(d[:B1]))-(occ_t.*exp(-d[:B1])+occ_c))/(occ_t.*exp(-d[:B1])+occ_c))
    odod_w = ( d[:o_mean_score1_pw] * d[:o_weight1] - d[:o_mean_score0_pw] * d[:o_weight0] ) /  (d[:o_mean_score0_pw] * d[:o_weight0])
    odiff = odod_w-odod
    ydod = (((dolocc_t+dolocc_c.*exp(d[:B2]) )-(dolocc_t.*exp(-d[:B2])+dolocc_c))/(dolocc_t.*exp(-d[:B2])+dolocc_c))
    ydod_w = ( d[:y_mean_score1_pw] * d[:y_weight1] - d[:y_mean_score0_pw] * d[:y_weight0] ) /  (d[:y_mean_score0_pw] * d[:y_weight0])
    ydiff = ydod_w-ydod
    
    h0 = (pen_t.*exp(-d[:B3])+pen_c)*(occ_t.*exp(-d[:B1])+occ_c)*(dolocc_t.*exp(-d[:B2])+dolocc_c)
    h1 = (pen_t+pen_c.*exp(d[:B3]))*(occ_t+occ_c.*exp(d[:B1]))*(dolocc_t+dolocc_c.*exp(d[:B2]))
    hdod = (h1-h0)/h0
    h0w = (d[:p_mean_score0_pw]*d[:p_weight0])*(d[:o_mean_score0_pw]*d[:o_weight0])*(d[:y_mean_score0_pw]*d[:y_weight0])
    h1w = (d[:p_mean_score1_pw]*d[:p_weight1])*(d[:o_mean_score1_pw]*d[:o_weight1])*(d[:y_mean_score1_pw]*d[:y_weight1])
    hdod_w = (h1w-h0w)/h0w
    hdiff = hdod_w-hdod
    println("Diff : $odod, $ydod, $pdiff, $hdiff")
    diffDict=Dict(:dolhh=>hdiff,:pen=>pdiff,:occ=>odiff,:dolocc=>ydiff)
    bounds_arr=[0.05,0.95, 0.10,0.90, 0.10,0.90, 0.20, 0.80]
    for (k,v) in OrderedDict(:dolhh=>:adj_dod,:pen=>:p_adj_dod, :occ=>:o_adj_dod, :dolocc=>:y_adj_dod)      
        q=(quantile(dfm[v], bounds_arr)+diffDict[k])*100
        d[k] =OrderedDict(:lb90=>q[1],:ub90=>q[2],:lb80=>q[3],:ub80=>q[4],:lb901=>q[5],:ub901=>q[6],:lb801=>q[7],:ub801=>q[8])
    end
end


function GC_ConfidenceIntervals()
    mDict=collectModels(dfx,"GLM")
    println("Running CI for : Total")   
    Glivenko_Cantelli(mDict)
    for m in [:dolhh,:occ,:dolocc,:pen]
        println(m)
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_90_pct_intrvl_lb] = mDict[m][:lb90]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_90_pct_intrvl_ub] = mDict[m][:ub90]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_80_pct_intrvl_lb] = mDict[m][:lb80]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_80_pct_intrvl_ub] = mDict[m][:ub80]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_90_pct_intrvl_lb] = mDict[m][:lb901]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_90_pct_intrvl_ub] = mDict[m][:ub901]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_80_pct_intrvl_lb] = mDict[m][:lb801]
        dfx[(dfx[:model].==string(m))&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_80_pct_intrvl_ub] = mDict[m][:ub801]
    end    
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
            mDict=collectModels(dfx,"GLMM",ranef,level)
            println("Running CI for : $ranef : $level")
            mDict[:metakey] = ranef*"~"*level        
            Glivenko_Cantelli(mDict)
            for m in [:dolhh,:occ,:dolocc,:pen]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_90_pct_intrvl_lb] = mDict[m][:lb90]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_90_pct_intrvl_ub] = mDict[m][:ub90]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_80_pct_intrvl_lb] = mDict[m][:lb80]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:twotail_80_pct_intrvl_ub] = mDict[m][:ub80]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_90_pct_intrvl_lb] = mDict[m][:lb901]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_90_pct_intrvl_ub] = mDict[m][:ub901]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_80_pct_intrvl_lb] = mDict[m][:lb801]
                dfx[(dfx[:model].==string(m))&(dfx[:ranef].==ranef)&(dfx[:parameter].==level),:onetail_80_pct_intrvl_ub] = mDict[m][:ub801]
            end         
        end
    end
end
GC_ConfidenceIntervals()  




        
        
            
  

function pvalues2C(dfx)      # pval = min(2 * min(cdf(TDist(Tn), t), ccdf(TDist(Tn), t)), 1.0)
        TµO = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:coef][1]  
        TsteO = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:stderr][1]
        TµY = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:coef][1]  
        TsteY = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:stderr][1]
        TµP = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:coef][1]  
        TsteP = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:stderr][1] 
        TµD = TµO + TµY + TµP
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh"),:])
            println("Pval : ", row[:ranef]," -- ",row[:parameter])
            #OCC
            BµO = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteO = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BσO = sqrt( BsteO^2 + TsteO^2 )
            tO = (BµO-TµO) / BσO
            pvalO = 2.0 * ccdf(Normal(), abs(tO))  #pvalO = min(2 * min(cdf(TDist(TnO), tO), ccdf(TDist(TnO), tO)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalO) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalO/2)) * 100
            #DOLOCC
            BµY = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteY = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BσY = sqrt( BsteY^2 + TsteY^2 )
            tY = (BµY-TµY) / BσY
            pvalY = 2.0 * ccdf(Normal(), abs(tY))  #pvalY = min(2 * min(cdf(TDist(TnY), tY), ccdf(TDist(TnY), tY)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalY) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalY/2)) * 100
            #PEN
            BµP = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteP = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BσP = sqrt( BsteP^2 + TsteP^2 )
            tP = (BµP-TµP) / BσP
            pvalP = 2.0 * ccdf(Normal(), abs(tP))    #pvalP = min(2 * min(cdf(TDist(TnP), tP), ccdf(TDist(TnP), tP)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalP) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalP/2)) * 100
            #DOLHH  
            BµD = BµO + BµY + BµP
            BσD =  sqrt( BsteO^2 + TsteO^2 + BsteY^2 + TsteY^2 + BsteP^2 + TsteP^2)
            tD = (BµD-TµD) / BσD
            pvalD = 2.0 * ccdf(Normal(), abs(tD))   #pval = min(2 * min(cdf(TDist(n), t), ccdf(TDist(n), t)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalD) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalD/2)) * 100  
        end
    end    
pvalues2C(dfx)
#dfx[(dfx[:modelType].=="GLMM"),[:ranef,:parameter,:onetail_pval_to_campaign,:twotail_pval_to_campaign]]
# "audience (Weeknight_Warrior)"  -- "creative_name (Baby_Video)"
        
                
        
        
        
        
        
        
dfx = dfx[(dfx[:modelType].=="GLMM")|((dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")),:]
dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100            
            
function genReport(dfx::DataFrame)
    genHHCounts();
    dfr=deepcopy(dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group")|(dfx[:modelType].=="GLMM"),:])
    rep=[ :MODEL_DESC, :Model, :TIME_AGG_PERIOD, :START_WEEK,:END_WEEK, :dependent_variable,:CNT_EXPSD_HH,  
          :UDJ_AVG_EXPSD_HH_PRE, :UDJ_AVG_CNTRL_HH_PRE, :UDJ_AVG_EXPSD_HH_PST,:UDJ_AVG_CNTRL_HH_PST,
          :UDJ_DOD_EFFCT,:UDJ_DIFF_EFFCT,
          :ADJ_MEAN_EXPSD_GRP,:ADJ_MEAN_CNTRL_GRP,:ADJ_DOD_EFFCT,:TWOTAIL_PVAL,:ONETAIL_PVAL,:ABS_DIFF,
          :DOL_DIFF,:ONETAIL_80_PCT_INTRVL_UB,:ONETAIL_80_PCT_INTRVL_LB,:ONETAIL_90_PCT_INTRVL_UB,:ONETAIL_90_PCT_INTRVL_LB,
          :TWOTAIL_80_PCT_INTRVL_UB,:TWOTAIL_80_PCT_INTRVL_LB,:TWOTAIL_90_PCT_INTRVL_UB,:TWOTAIL_90_PCT_INTRVL_LB,
          :CNT_IMPRESSIONS,:TWOTAIL_PVAL_to_Campaign,:ONETAIL_PVAL_to_Campaign,:CNT_Model_HH
        ]
     dfn=[ :model_desc, :empty,:empty_1, :empty_2, :empty_3, :model, :cnt_expsd_hh, :unadj_avg_expsd_hh_pre, :unadj_avg_cntrl_hh_pre,
           :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pst, :unadj_dod_effct, :unadj_diff_effct, :adj_mean_score1, :adj_mean_score0,
           :adj_dod_effct, :twotail_pval, :onetail_pval, :abs_diff, :dol_diff, :onetail_80_pct_intrvl_ub, :onetail_80_pct_intrvl_lb,
           :onetail_90_pct_intrvl_ub, :onetail_90_pct_intrvl_lb, :twotail_80_pct_intrvl_ub, :twotail_80_pct_intrvl_lb,
           :twotail_90_pct_intrvl_ub, :twotail_90_pct_intrvl_lb, :cnt_impressions, :twotail_pval_to_campaign,
           :onetail_pval_to_campaign,:cnt_model_hh
         ]
    #--
    dfr[isnan(dfr[:unadj_avg_expsd_hh_pre]),:unadj_avg_expsd_hh_pre] = 0.0  # Sometime there are no records for subset
    dfr[isnan(dfr[:unadj_avg_expsd_hh_pst]),:unadj_avg_expsd_hh_pst] = 0.0  # Sometime there are no records for subset
    #dfr[:adj_dod_effct] = ((dfr[:adj_mean_score1] .- dfr[:adj_mean_score0]) ./ dfr[:adj_mean_score0] ) *100
    dfr[:unadj_dod_effct] = ( (( dfr[:unadj_avg_expsd_hh_pst] .- dfr[:unadj_avg_expsd_hh_pre]) .- (dfr[:unadj_avg_cntrl_hh_pst ] .- dfr[:unadj_avg_cntrl_hh_pre]))  ./  dfr[:unadj_avg_cntrl_hh_pst] ) *100
    dfr[:unadj_diff_effct] = ((dfr[:unadj_avg_expsd_hh_pst] .- dfr[:unadj_avg_cntrl_hh_pst]) ./ dfr[:unadj_avg_cntrl_hh_pst] )*100
    dfr[:model_desc] = dfr[:ranef]*" (".*dfr[:parameter]*")"
    dfr[(dfr[:modelType].=="GLM"),:model_desc] = "Total Campaign"
    dfr[:abs_diff] = dfr[:adj_mean_score1] .- dfr[:adj_mean_score0]
    dfr[:dol_diff] = dfr[:adj_mean_score1] .- dfr[:adj_mean_score0]
    dfr[:cnt_model_hh] = dfr[:Nt]
    dfr[findin(dfr[:model],["occ","dolocc"]),:cnt_model_hh] = dfr[findin(dfr[:model],["occ","dolocc"]),:Mt]       
    dfr[:onetail_pval] = (1-(dfr[:pval] ./ 2)) * 100  # sdf[:onetail_pval_raw] = (1-(sdf[:Praw] ./ 2)) * 100
    dfr[:twotail_pval] = (1-dfr[:pval]) * 100         # sdf[:twotail_pval_raw] = (1-sdf[:Praw]) * 100
    dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"),:adj_pval] = 2.0 * ccdf(Normal(),abs(dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"),:adj_coef] ./ dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"),:adj_stderr]))
    #dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"), :onetail_pval_to_campaign] =  (1-(dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"), :adj_pval] ./ 2) ) * 100
    #dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"), :twotail_pval_to_campaign] =  (1-dfr[(dfr[:modelType].=="GLMM")&(dfr[:model].!="dolhh"), :adj_pval] ) * 100
#--
    dfr[:empty] = ""
    dfr[:orderModel] = map(x-> x=="pen"? 1 : x=="occ" ? 2 : x=="dolocc" ? 3 : x=="dolhh" ? 4 : 9999 ,dfr[:model])
    dfr[:orderFixedRand] = map(x-> x=="GLM"? 1 : 9999 ,dfr[:modelType])
    sort!(dfr, cols = [:orderFixedRand,:ranef,:parameter,:orderModel])
    dfr=dfr[[:model_desc,:empty,:empty,:empty,:empty,:model,:cnt_expsd_hh,
             #:unadj_avg_expsd_hh_pre, :unadj_avg_expsd_hh_pst,:unadj_avg_cntrl_hh_pre, :unadj_avg_cntrl_hh_pst,
             :unadj_avg_expsd_hh_pre, :unadj_avg_cntrl_hh_pre, :unadj_avg_expsd_hh_pst, :unadj_avg_cntrl_hh_pst,
             :unadj_dod_effct,:unadj_diff_effct,
             :adj_mean_score1,:adj_mean_score0,:adj_dod_effct,:twotail_pval,:onetail_pval,:abs_diff,:dol_diff,
             :onetail_80_pct_intrvl_ub,:onetail_80_pct_intrvl_lb,:onetail_90_pct_intrvl_ub,:onetail_90_pct_intrvl_lb,
             :twotail_80_pct_intrvl_ub,:twotail_80_pct_intrvl_lb,:twotail_90_pct_intrvl_ub,:twotail_90_pct_intrvl_lb
             ,:cnt_impressions ,:twotail_pval_to_campaign,:onetail_pval_to_campaign ,:cnt_model_hh
        ]]
    #names!(dfr,rep)
    for i in 1:length(rep) rename!(dfr,dfn[i],rep[i]) end
                
    #dfr[(dfr[:MODEL_DESC].=="Total Campaign")|(dfr[:dependent_variable].=="dolhh"),:TWOTAIL_PVAL_to_Campaign] = NaN
    #dfr[(dfr[:MODEL_DESC].=="Total Campaign")|(dfr[:dependent_variable].=="dolhh"),:ONETAIL_PVAL_to_Campaign] = NaN
    dfr[(dfr[:MODEL_DESC].=="Total Campaign"),:TWOTAIL_PVAL_to_Campaign] = NaN
    dfr[(dfr[:MODEL_DESC].=="Total Campaign"),:ONETAIL_PVAL_to_Campaign] = NaN
    dfr[(dfr[:dependent_variable].!="dolhh"),:CNT_EXPSD_HH] = NA
    dfr[(dfr[:dependent_variable].!="dolhh"),:CNT_IMPRESSIONS] = NA
    return dfr
end
dfr = genReport(dfx)
save_dfr(root, dfr)


#dfr[(dfr[:dependent_variable].=="occ"),:]

#Dong, Zhiyuan
#pvalue < 0.05 = statistical significance level = reject null hypothesis = and conclude that no matching on trps_pre_p1
function kstest()
    #using HypothesisTests
    ks=ApproximateTwoSampleKSTest(dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==0),:trps_pre_p1],dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==1),:trps_pre_p1])
    o=sqrt(ks.n_x*ks.n_y/(ks.n_x+ks.n_y))*ks.δ 
    ks=ApproximateTwoSampleKSTest(dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==0),:dol_per_trip_pre_p1],dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==1),:dol_per_trip_pre_p1])
    y=sqrt(ks.n_x*ks.n_y/(ks.n_x+ks.n_y))*ks.δ         
    ks=ApproximateTwoSampleKSTest(dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==0),:buyer_pre_p1],dfd[(dfd[:buyer_pre_p1].==1)&(dfd[:group].==1),:buyer_pre_p1])
    p=sqrt(ks.n_x*ks.n_y/(ks.n_x+ks.n_y))*ks.δ 
    return Dict(:occ=>o,:dolocc=>y, :pen=>p)
end        
        
        
#dfr[(dfr[:dependent_variable].=="dolhh"),[:MODEL_DESC,:ADJ_DOD_EFFCT,:ONETAIL_80_PCT_INTRVL_LB,:ONETAIL_80_PCT_INTRVL_UB]]
#dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
#dfx[(dfx[:adj_dod_effct].<dfx[:onetail_80_pct_intrvl_lb])|(dfx[:adj_dod_effct].>dfx[:onetail_80_pct_intrvl_ub]),:]
        
#dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].=="campaign_name")&(dfx[:parameter].=="Halloween/Share_A_Smile"),:] 
#dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].=="creative_size")&(dfx[:parameter].=="300x250/728x90"),:]            
            
            
            
        
        
        
        
        
        
        
        
            
            
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# xyz*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- TESt EVAL -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    function pvalues2CampaignOLD(dfx::DataFrame)
    for model in ["occ","dolocc","pen"]       
        Bt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:coef][1]  
        Errt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:stderr][1]
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model),:])
            Bb = row[:adj_coef][1]
            Errb = row[:adj_stderr][1]
            t= (Bb-Bt) / sqrt(   (  ((Errb*sqrt(1000 ))^2) /1000   ) 
                               + (  ((Errt*sqrt(10000))^2)/10000  ) 
                             )
            pval = 2 * cdf(TDist(1000),abs(t))    #cdf(TDist(5.3900745416769-1),-0.7687842353483888)*2
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
        end  
    end            
    for model in ["dolhh"]
        ot = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:coef][1] 
        yt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:coef][1]
        pt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:coef][1]            
        Bt = ot*yt*pt
        ote = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:stderr][1] 
        yte = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:stderr][1]
        pte = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:stderr][1]  
        Errt = sqrt(ote^2+yte^2+pte^2)
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model),:])
            ob = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:parameter].==row[:parameter]),:coef][1] 
            yb = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:parameter].==row[:parameter]),:coef][1]
            pb = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:parameter].==row[:parameter]),:coef][1]            
            Bb = ob*yb*pb
            obe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:parameter].==row[:parameter]),:stderr][1] 
            ybe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:parameter].==row[:parameter]),:stderr][1]
            pbe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:parameter].==row[:parameter]),:stderr][1]  
            Errb = sqrt(obe^2+ybe^2+pbe^2)
            t= (Bb-Bt) / sqrt(   (  ((Errb*sqrt(1000 ))^2) /1000   ) 
                               + (  ((Errt*sqrt(10000))^2)/10000  ) 
                             )
                        pval = 2 * cdf(TDist(1000),abs(t))    #cdf(TDist(5.3900745416769-1),-0.7687842353483888)*2
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
        end  
    end         
end
#pvalues2CampaignOLD(dfx)
 
function pvalues2CampaignNEW(dfx::DataFrame)
    for model in ["occ","dolocc","pen"]
        Tµ = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:coef][1]  
        Tste = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:stderr][1]
        ncol = model=="pen" ? :N : :M
        Tn=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),ncol][1]
        Tstd = Tste*sqrt(Tn)
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model),:])
            Bµ = row[:adj_coef][1]
            Bste = row[:adj_stderr][1]
            Bn = row[ncol][1]
            Bstd = Bste*sqrt(Bn)
            #t = (Bµ-Tµ) / sqrt( ((  (Bste*sqrt(Bn))    ^2)/Bn) + (((Tste*sqrt(Tn))^2)/Tn)  )
            t = (Bµ-Tµ) / sqrt( ((Bstd^2)/Bn) + ((Tstd^2)/Tn) )
            pval = min(2 * min(cdf(TDist(Tn), t), ccdf(TDist(Tn), t)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
            #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:pval2camp] = pval
        end  
    end 
    dfd[:dolhh] = dfd[:trps_pos_p1] .* dfd[:dol_per_trip_pos_p1]
    Tµ=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolhh"),:adj_mean_score1][1]
    Tstd=std(dfd[:dolhh])
    Tn = length(dfd[:dolhh])
    for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh"),:])
        Bµ = row[:adj_mean_score1]
        Bstd = std(dfd[(dfd[ Symbol(row[:ranef]) ].==row[:parameter]),:dolhh])
        Bn = length(dfd[(dfd[ Symbol(row[:ranef]) ].==row[:parameter]),:dolhh])
        t= (Bµ-Tµ) / sqrt( ((Bstd^2)/Bn) + ((Tstd^2)/Tn) )
        pval = min(2 * min(cdf(TDist(Tn), t), ccdf(TDist(Tn), t)), 1.0)
        dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
        dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
        #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:parameter].==row[:parameter]),:pval2camp] = pval
        println(row[:parameter]," ~~ ",Bµ," :: ",Bstd," .... ",t," ~~ ",pval)
    end            
end
#pvalues2CampaignNEW(dfx)
        
        
function pvalues2Campaign(dfx::DataFrame)
    for model in ["occ","dolocc","pen"]
        Tµ = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:coef][1]  
        Tste = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),:stderr][1]
        ncol = model=="pen" ? :N : :M
        Tn=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==model)&(dfx[:parameter].=="group"),ncol][1]
        Tstd = Tste*sqrt(Tn)
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model),:])
            Bµ = row[:adj_coef][1]
            Bste = row[:adj_stderr][1]
            Bn = row[ncol][1]
            Bstd = Bste*sqrt(Bn)
            #t = (Bµ-Tµ) / sqrt( ((  (Bste*sqrt(Bn))    ^2)/Bn) + (((Tste*sqrt(Tn))^2)/Tn)  )
            t = (Bµ-Tµ) / sqrt( ((Bstd^2)/Bn) + ((Tstd^2)/Tn) )
            pval = min(2 * min(cdf(TDist(Tn), t), ccdf(TDist(Tn), t)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
            #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:pval2camp] = pval
        end  
    end 
    stds(µ::Float64,ste::Float64) = std(rand(Normal(µ, ste),100000))
    for model in ["dolhh"]
        ot = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:coef][1] 
        yt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:coef][1]
        pt = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:coef][1]            
        Tµ = ot+yt+pt
        ote = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:stderr][1] 
        yte = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:stderr][1]
        pte = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:stderr][1]  
        Tste = sqrt(ote^2+yte^2+pte^2)
                println("DOLHH Total :: $ot * $yt * $pt = $Tµ ~~~~ $ote + $yte + $pte = $Tste")
        Tstd = stds(Tµ,Tste)
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model),:])
            ob = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:parameter].==row[:parameter]),:adj_coef][1] 
            yb = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            pb = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:parameter].==row[:parameter]),:adj_coef][1]            
            Bµ = ob+yb+pb
            obe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:parameter].==row[:parameter]),:adj_stderr][1] 
            ybe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            pbe = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]  
            Bste = sqrt(obe^2+ybe^2+pbe^2)
            Bstd = stds(Bµ,Bste)
            t = (Bµ-Tµ) / sqrt( ((Bstd^2)/100000) + ((Tstd^2)/100000) )
            pval = min(2 * min(cdf(TDist(100000), t), ccdf(TDist(100000), t)), 1.0)
                    println("DOLHH Break :: $ob * $yb * $pb = $Bµ ~~~~ $obe + $ybe + $pbe = $Bste ..... pval := $pval")
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==model)&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
        end  
    end            
end
#pvalues2Campaign(dfx)        
        
       
        
        
        
        
        
 
    function pvalues2C(dfx)      
        TµO = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:coef][1]  
        TsteO = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:stderr][1]
        TnO=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="occ")&(dfx[:parameter].=="group"),:M][1]
        TstdO = TsteO*sqrt(TnO)
        TxO=((TstdO^2)/TnO)
        
        TµY = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:coef][1]  
        TsteY = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:stderr][1]
        TnY=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="dolocc")&(dfx[:parameter].=="group"),:M][1]
        TstdY = TsteY*sqrt(TnY)
        TxY=((TstdY^2)/TnY)
        
        TµP = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:coef][1]  
        TsteP = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:stderr][1]
        TnP=dfx[(dfx[:modelType].=="GLM")&(dfx[:model].=="pen")&(dfx[:parameter].=="group"),:N][1]
        TstdP = TsteP*sqrt(TnP)
        TxP=((TstdP^2)/TnP)     
            
        #Tn=TnP
        #Td =  
            
        for row in eachrow(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh"),:])
            println("Pval : ", row[:ranef]," -- ",row[:parameter])
            #OCC
            BµO = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteO = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BnO = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:M][1]
            BstdO = BsteO*sqrt(BnO)
            BxO=((BstdO^2)/BnO)
            tO = (BµO-TµO) / sqrt( BxO + TxO )
            pvalO = min(2 * min(cdf(TDist(TnO), tO), ccdf(TDist(TnO), tO)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalO) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalO/2)) * 100
            #DOLOCC
            BµY = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteY = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BnY = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:M][1]
            BstdY = BsteY*sqrt(BnY)
            BxY=((BstdY^2)/BnY)
            tY = (BµY-TµY) / sqrt( BxY + TxY )
            pvalY = min(2 * min(cdf(TDist(TnY), tY), ccdf(TDist(TnY), tY)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalY) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolocc")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalY/2)) * 100
            #PEN
            BµP = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_coef][1]
            BsteP = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:adj_stderr][1]
            BnP = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:N][1]
            BstdP = BsteP*sqrt(BnP)
            BxP=((BstdP^2)/BnP)
            tP = (BµP-TµP) / sqrt( BxP + TxP )
            pvalP = min(2 * min(cdf(TDist(TnP), tP), ccdf(TDist(TnP), tP)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pvalP) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="pen")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pvalP/2)) * 100
                
                
            #DOLHH  
            t = ((BµO+BµY+BµP)-(TµO+TµY+TµP)) / sqrt( BxO + TxO + BxY + TxY + BxP + TxP )
            n = TnP
            pval = min(2 * min(cdf(TDist(n), t), ccdf(TDist(n), t)), 1.0)
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:twotail_pval_to_campaign] = (1-pval) * 100
            dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh")&(dfx[:ranef].==row[:ranef])&(dfx[:parameter].==row[:parameter]),:onetail_pval_to_campaign] = (1-(pval/2)) * 100
            
            
        end
    end    
#pvalues2C(dfx)
dfx[(dfx[:modelType].=="GLMM"),[:ranef,:parameter,:onetail_pval_to_campaign,:twotail_pval_to_campaign]]
#dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="dolhh"),[:ranef,:parameter,:onetail_pval_to_campaign,:twotail_pval_to_campaign]]       
        
        