
#export JULIA_NUM_THREADS=4
ENV["JULIA_PKGDIR"] = "/mapr/mapr04p/analytics0001/analytic_users/jpkg" 
for i in ["10.106.128.212","10.106.128.213","10.106.128.214","10.106.128.215","10.106.128.216","10.106.128.217",
          "10.106.128.218","10.106.128.219","10.106.128.220","10.106.128.221"]
             if i != string(getipaddr()) println(i); addprocs([(i, 1)]) end
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




root = !isdefined(:root) ? pwd() : root

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

root = !isdefined(:root) ? pwd() : root
cfgDefaults = mergeDict(dict_Sym(read_cfg(root)),cfgDefaults)

#function loadDF()
#    df_data = readtable(root*"/orig.csv",header=false);
#    df_h = readtable(root*"/origHead.csv",header=false);  
#    names!(df_data, convert(Array{Symbol}, df_h[:x1]) )
#end
#df_in=loadDF()
df_in=read_orig(root)

function Pre_out(df_in::DataFrame)
     df_cat_pre = df_in[df_in[:Buyer_Pre_P0] .==1 , [:Prd_0_Net_Pr_PRE,:experian_id]]
     df_cat_pos = df_in[(df_in[:Buyer_Pre_P0] .==0) & (df_in[:Buyer_Pos_P0] .==1) , [:experian_id]]
     median_df = median(df_cat_pre[:Prd_0_Net_Pr_PRE])
     df_cat_pre[:Prd_0_Net_Pr_PRE_med1] = abs(df_cat_pre[:Prd_0_Net_Pr_PRE]-median_df)
     MAD=median(df_cat_pre[:Prd_0_Net_Pr_PRE_med1])
     df_cat_pre[:Prd_0_Net_Pr_PRE_med2] = (0.6745*(df_cat_pre[:Prd_0_Net_Pr_PRE]-median_df))/MAD
     df_cat_pre_zsc = df_cat_pre[abs(df_cat_pre[:Prd_0_Net_Pr_PRE_med2]) .< 3.5,:]
     df_cat_pre_zsc_1 = df_cat_pre_zsc[:,[:experian_id]]
     df_cat_pre_zsc_f = vcat(df_cat_pos,df_cat_pre_zsc_1)
     df_in_pout =  join(df_in, df_cat_pre_zsc_f, on =  :experian_id , kind = :inner);
end
df_in = Pre_out(df_in);


function isValid(df_data::DataFrame,cfg::OrderedDict)
    function checkValid(iarr::Array{Symbol})  length(setdiff(iarr, names(df_data))) > 0 ? false : true end
    !checkValid(cfg[:all_mandatory_vars]) ? error("ERROR: Not all mandatory_vars in dataset ") : println("VALID : mandatory_vars") 
    !checkValid(cfg[:scoring_vars]) ? error("ERROR: Not all scoring_vars in dataset ") : println("VALID : scoring_vars") 
end


function getCFG(df_in::DataFrame)
    cfg=StatLib.loadCFG(cfgDefaults, pwd()*"/app.cfg")
    cfg[:allrandoms] = vcat(cfg[:random_demos],cfg[:random_campaigns]) 
    ps = filter(x->contains(string(x), "MODEL"), names(df_in)) 
    cfg[:ProScore] = length(ps) > 0 ? ps[1] : :MISSING_MODEL_VARIABLE_IN_DATA
    cfg[:num_products] = length( grep("Buyer_Pre_P",names(df_in)))-1 
    cfg[:random_campaigns] = intersect(cfg[:random_campaigns],names(df_in)) 
    cfg[:occ_logvar_colname] = Symbol("LOG_"*string(cfg[:occ_logvar]))
    cfg[:dolocc_logvar_colname] = Symbol("LOG_"*string(cfg[:dolocc_logvar]))
    cfg[:pen_logvar_colname] = Symbol("LOG_"*string(cfg[:pen_logvar]))
    cfg[:all_mandatory_vars] = vcat( [:experian_id, cfg[:pen_y_var],cfg[:pen_logvar], cfg[:occ_y_var], cfg[:occ_logvar], cfg[:dolocc_y_var], cfg[:dolocc_logvar] ]
                                      , cfg[:random_demos]
                                   )
    return cfg
end
cfg=getCFG(df_in)

isValid(df_in,cfg)

function reworkCFG!(df_in::DataFrame,cfg::OrderedDict)
    xflags=[cfg[:exposed_flag_var],:exposed_flag,:Nonbuyer_Pre_P1]  # :Nonbuyer_Pre_P1 no longer used
    features = setdiff(names(df_in),xflags)
    if :state in names(df_in)
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
    df_in[:group] = df_in[cfg[:exposed_flag_var]]
    df_in[:panid] = df_in[:experian_id] 
    features = setdiff(unique(vcat(features,cfg[:iVarsPREPOS],cfg[:all_mandatory_vars],cfg[:scoring_vars],[:group,:panid])  ) ,[:experian_id] )
    cfg[:negativevars] = grep(vec(hcat([["P"*string(i),string(i)*"_"] for i=3:cfg[:num_products]]...)),features)   # get variables that need to have negative sign
    cfg[:positivevars] = grep(["P1", "1_","MODEL"], features)    # get variables that need to have positive sign
    if cfg[:P2_Competitor] == true  
        cfg[:negativevars] = unique(vcat(cfg[:negativevars],grep(["P2","2_"],features))) 
    else
        cfg[:positivevars] = unique(vcat(cfg[:positivevars],grep(["P2","2_"],features))) 
    end
    return df_in[features]
end
dfd = reworkCFG!(df_in,cfg)


function data_Prep(dfd::DataFrame, cfg::OrderedDict)
    dfd[:isO]=false
    dfd[dfd[:number_of_children_in_living_Un].>=4,:number_of_children_in_living_Un] = 4  # aggregate #of children for 4+ L_114
    if typeof(dfd[:group]) in [DataArray{String,1}] 
        dfd[ findin(dfd[:group],["//N","\\N"]), :group] = "0" 
        dfd[DataArrays.isna(dfd[:group]), :group]="0"
        dfd[:group] = [parse(Int64,s) for s = dfd[:group]]
    else
        dfd[ DataArrays.isnan(dfd[:group]), :group] = 0
    end
    vars = DataFrame(names=names(dfd),eltypes=eltypes(dfd))
    for c in setdiff(vars[findin(vars[:eltypes],[String]),:names],cfg[:allrandoms])
        print("Convert String: String->Numeric: ",c)
        try dfd[c] = map( x -> DataArrays.isna.(x) ?  NaN : convert(Float64, x)  , dfd[c]) catch e println("  (failed)") end  #NA to NaN
        try dfd[c] = convert(Array{Float64}, dfd[c]) catch e end 
    end
    vars = DataFrame(names=names(dfd),eltypes=eltypes(dfd))
    for c in setdiff(vars[findin(vars[:eltypes],[Float64]),:names],cfg[:allrandoms])  
        println("Replace Float64 NaN (0.0): ",c)
        dfd[ DataArrays.isnan(dfd[c]), c] = 0.0
    end
    for c in  setdiff(vars[findin(vars[:eltypes],[Int64]),:names],cfg[:allrandoms])    
        println("Replace Int64 NaN (0) : ",c)
        dfd[ DataArrays.isnan(dfd[c]), c] = 0
    end
    dfd[dfd[:person_1_gender].=="U",:isO] = true  
    dfd[findin(dfd[:estimated_hh_income],["U","L"]),:estimated_hh_income]="L" 
    for r in cfg[:random_campaigns]    
        dfd[findin(dfd[r],["\\N","NULL","0","NONE"])  ,r] ="none"
        dfd[ (dfd[:isO].==false) & (dfd[r].!="none") & (dfd[:group].==0) ,:isO] = true 
        dfd[ (dfd[:isO].==false) & (dfd[r].=="none") & (dfd[:group].!=0) ,:isO] = true
    end
    dfd[:data_NB_NE_B] = false
    dfd[ (dfd[:Buyer_Pre_P1].==0 ) & (dfd[:group].==0 ) & (dfd[:Buyer_Pos_P1].==1 ) ,:data_NB_NE_B] = true
    dfd[:data_B_E_NB] = false
    dfd[ (dfd[:Buyer_Pre_P1].==1 ) & (dfd[:group].==0 ) & (dfd[:Buyer_Pos_P1].==0 ) ,:data_B_E_NB] = true
    dfd[:pen_reduction] = false
    dfd[ (dfd[:Buyer_Pre_P1].==1) & (dfd[:group].==0) & (dfd[:Buyer_Pos_P1].==0 )  ,:pen_reduction] = true
    dfd[:occ_reduction] = false
    dfd[ (dfd[:group].==0) & (dfd[:Buyer_Pos_P1].==1) & (dfd[:Trps_POS_P1].< dfd[:Trps_PRE_P1] ) ,:occ_reduction] = true
    dfd[:dolocc_reduction] = false
    dfd[  (dfd[:group].==0) & (dfd[:Buyer_Pos_P1].==1) & (dfd[:Dol_per_Trip_POS_P1].< dfd[:Dol_per_Trip_PRE_P1] )  , :dolocc_reduction] = true
    return dfd[dfd[:isO].==false, : ]   
end
dfd = data_Prep(dfd, cfg);

function Restrict50Buyers()
    for (k,v) in Dict(r=>countmap(dfd[(dfd[:Buyer_Pos_P1].==1), r]) for r in cfg[:random_campaigns])
            lvls=String[]
            for (kl,vl) in v 
                if vl < 50
                    println("$k :: $kl to Other")
                    dfd[dfd[k].==kl,k] = "Other" 
                end  
                #if vl > 50 push!(lvls,kl) end  
            end
            #toNone = setdiff(levels(dfd[k]),lvls)
            #println("\n\n\n$k       ::::          ",toNone)
    end
end
Restrict50Buyers()    

function MatchMe(dfd::DataFrame,cfg::OrderedDict)
    df=dfd[dfd[:isO].==false,:]
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
    rows2remove = setdiff(dfd[dfd[:isO].==false, :panid],dfd_sample[:panid])
    #dfd[findin(dfd[:panid],rows2remove),:whyO]="NoMatch"
    dfd[findin(dfd[:panid],rows2remove),:isO]=true 
        

#    --------------Matching - pre-post ------------------
"""
    dfp=join(by(dfd[(dfd[:isO].==false),:],[:group,:Buyer_Pos_P1], df->DataFrame(cnt=size(df,1))) ,  
             by(dfd[(dfd[:isO].==false),:],:group,df->DataFrame(tot=size(df,1)))
             , on=:group 
            )  
    dfp[:prop] = dfp[:cnt]./dfp[:tot]
    if dfp[(dfp[:group].==1)&(dfp[:Buyer_Pos_P1].==1),:prop][1] < dfp[(dfp[:group].==0)&(dfp[:Buyer_Pos_P1].==1),:prop][1]  
            
        scnt = Int( round(( dfp[(dfp[:group].==0)&(dfp[:Buyer_Pos_P1].==0),:cnt][1]*dfp[(dfp[:group].==1)&(dfp[:Buyer_Pos_P1].==1),:prop][1] ) / (1-dfp[(dfp[:group].==1)&(dfp[:Buyer_Pos_P1].==1),:prop][1])  ))
        smp = sample(dfd[(dfd[:isO].==false)&(dfd[:Buyer_Pos_P1].==1)&(dfd[:group].==0),:panid], scnt , replace=false)
       
        oliers = setdiff(dfd[(dfd[:isO].==false)&(dfd[:Buyer_Pos_P1].==1)&(dfd[:group].==0),:panid],smp)
        
        dfd[findin(dfd[:panid], setdiff(dfd[(dfd[:isO].==false)&(dfd[:Buyer_Pos_P1].==1)&(dfd[:group].==0),:panid],smp)),:isO] =true
   end
"""
        # ---- end matching ------------
        
        
        
    return dfd[dfd[:isO].==false, : ]  #[setdiff(names(dfd),[:isO,:whyO])] 
end
dfd = MatchMe(dfd,cfg)

lowercase!(dfd)
cfg=lowercase(cfg)

iocc = Dict(:modelName=>:occ, :raneff=>cfg[:random_campaigns], :y_var=>:trps_pos_p1, :logvar=>:LOG_trps_pre_p1, :logvarOrig=>:trps_pre_p1 )
idolocc = Dict(:modelName=>:dolocc, :raneff=>cfg[:random_campaigns], :y_var=>:dol_per_trip_pos_p1, :logvar=>:LOG_dol_per_trip_pre_p1, :logvarOrig=>:dol_per_trip_pre_p1)
ipen = Dict(:modelName=>:pen, :raneff=>cfg[:random_campaigns], :y_var=>:buyer_pos_p1, :logvar=>:LOG_buyer_pre_p1, :logvarOrig=>:buyer_pre_p1 )
custom_vars = [:iso,:whyo,:data_nb_ne_b, :data_b_e_nb ,:pen_reduction,:occ_reduction,:dolocc_reduction]  
for m in [iocc,idolocc,ipen] dfd[m[:logvar]]=log(Array(dfd[m[:logvarOrig]]+1)) end
function genExcludeVars!(iocc::Dict,idolocc::Dict,ipen::Dict)  
    iocc[:exclude_vars] = vcat( custom_vars,  
                                [ :buyer_pos_p1, iocc[:logvarOrig],idolocc[:y_var],idolocc[:logvar],idolocc[:logvarOrig] ,ipen[:y_var],ipen[:logvar],ipen[:logvarOrig] ]
                              )
    idolocc[:exclude_vars] = vcat( custom_vars,  
                                   [ :buyer_pos_p1, idolocc[:logvarOrig] ,iocc[:y_var],iocc[:logvar],iocc[:logvarOrig] ,ipen[:y_var],ipen[:logvar],ipen[:logvarOrig] ]
                                 )    
    ipen[:exclude_vars] = vcat( custom_vars,
                                [ ipen[:logvarOrig] ,iocc[:y_var],iocc[:logvar],iocc[:logvarOrig],idolocc[:y_var],idolocc[:logvar],idolocc[:logvarOrig] ]
                                 )        
end
genExcludeVars!(iocc,idolocc,ipen)

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
    if length(finalvars) > 2     #VIF is only valid on multiple cols
        println( upperModName *" : Z & Vif") #Z & Vif
        f=genF(m[:y_var],finalvars)
        g2 = glm( f, dfd[convert(Array{Symbol},vcat(f.lhs,f.rhs.args[3:end]))]  , m[:dist], m[:lnk] )
        sdf=vifDF2(g2)
        z = sdf[abs(sdf[:zval]).<1.96,:parameter]
        v = sdf[ !DataArrays.isna(sdf[:vif])&(sdf[:vif].>15),:parameter]
        g2_x =intersect(z,v)
        finalvars = setdiff(finalvars,g2_x)
        println("vif_vars: ",g2_x)
    end
    finalvars = setdiff(finalvars,[:group])# reorder for group
    finalvars = convert(Array{Symbol},vcat(finalvars,[:group]) )
    return finalvars
end
    
    

factor_cols=vcat( [ cfg[:proscore], :group, :panid], cfg[:allrandoms] )
for c in setdiff(factor_cols,[:panid, cfg[:proscore]]) #cfg[:random_campaigns]
    if !( typeof(dfd[c]) in [  DataArray{String,1}  ]) 
        println("converting to Strings : ", c," of type : ",typeof(dfd[c]))
        dfd[c] = map(x->string(x),dfd[c])
        dfd[c] = convert(Array{String},dfd[c]) 
    end
end
poolit!(dfd,factor_cols)

iocc[:finalvars] = featureSelection(dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1),:], expandM(iocc))
idolocc[:finalvars]   = featureSelection(dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1),:], expandM(idolocc))
ipen[:finalvars]  = featureSelection(dfd[(dfd[:iso].==false) ,:], expandM(ipen))

    
"""
dfcor(dfr::DataFrame) = dfcor(dfr, names(dfr))

function dfcor(dfr::DataFrame, cols::Vector{Symbol})
    nms = Symbol[]
    arr = Float64[]
    for (n, v) in eachcol(dfr)
        if n ∈ cols && isnumeric(v)
            push!(nms, n)
            append!(arr, v)
        end
    end
    result = NamedArray(cor(reshape(arr, (size(dfr, 1), length(nms)))))
    NamedArrays.setnames!(result, string.(nms), 1)
    NamedArrays.setnames!(result, string.(nms), 2)
    result
end
    
dfcor1(dfr::DataFrame) = dfcor1(dfr, names(dfr))

function dfcor1(dfr::DataFrame, cols::Vector{Symbol})
    m = size(dfr, 1)
    nms = Symbol[]
    arr = ones(m)
    for (n, v) in eachcol(dfr)
        if n ∈ cols && isnumeric(v)
            push!(nms, n)
            append!(arr, v)
        end
    end
    n = length(nms) + 1
    R = full(UpperTriangular(view(qrfact!(reshape(arr, (m, n)))[:R], 
                2:n, 2:n)))
    for j in 1:size(R, 2)
        colj = view(R, 1:j, j)
        colj ./= norm(colj)
    end
    Rt = UpperTriangular(R)
    result = NamedArray(Rt'Rt)
    NamedArrays.setnames!(result, string.(nms), 1)
    NamedArrays.setnames!(result, string.(nms), 2)
    result
end
"""
# --- Removing outlliers - to ensure overall campaign uplift
function getOutliers(r::DataFrame, dfd::DataFrame, m::Dict, segment::Symbol,pct::Int64) # :data_nb_ne_b, :data_b_e_nb, :pen_reduction, :occ_reduction, :dolocc_reduction
    if m[:Buyer_Pos_P1_is1]
        println("pre=1")
        odf = dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1)&(dfd[segment].==true),:panid]
    else
        println("pre!=1")
        odf = dfd[(dfd[:iso].==false)&(dfd[segment].==true),:panid]
    end
    pctnum = Integer(round(length(odf) / 100) * pct)
    return sort(r[findin(r[:panid],odf) ,:], cols=[order(:resids, rev = true)] )[1:pctnum,:panid]
end

function genGLM(dfd::DataFrame, m::Dict)
    f=genF(m[:y_var],m[:finalvars])
    println("GLM for ",m[:modelName]," ::: ",f)
    dfd_tmp = m[:modelName] in [:occ, :dolocc] ? dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1),:] : dfd[(dfd[:iso].==false),:] 
    g = glm(f, dfd_tmp[vcat(m[:y_var],m[:finalvars])], m[:dist], m[:lnk])
    println("GLM Residuals for ",m[:modelName])
    r = DataFrame(panid=dfd_tmp[:panid], resids=StatLib.xResiduals(g))
    return r, g 
end

function genlists(dfd::DataFrame, iocc::Dict,idolocc::Dict,ipen::Dict)   
    r1, g1 = genGLM( dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1),:], expandM(iocc))
    r2, g2 = genGLM( dfd[(dfd[:iso].==false)&(dfd[:buyer_pos_p1].==1),:], expandM(idolocc))
    r3, g3 = genGLM( dfd[(dfd[:iso].==false),:], expandM(ipen))
    rlist = [r1,r2,r3]
    glist = [g1,g2,g3]
    return rlist, glist
end
rlist, glist = genlists(dfd,iocc,idolocc,ipen) 

    
rmOutliers=false
if rmOutliers
    for i in 1:5
        mlist = [expandM(iocc),expandM(idolocc),expandM(ipen)]
        #mliftβ = [ g.sdf[g.sdf[:vars].==:group,:coef][1] for g in glist]
        mliftβ = [ sdf[sdf[:parameter].==:group,:coef][1] for sdf in [coefDF(m, true) for m in glist]]
        minidx = find(x->x==minimum(mliftβ),mliftβ)[1]
        if  mliftβ[minidx] < 0.0
            #mx=mlist[minidx]
            println("WRONG : ",mliftβ," for ",mlist[minidx][:modelName])
            panids = getOutliers(rlist[minidx], dfd, mlist[minidx], :data_nb_ne_b,5)
            dfd[findin(dfd[:panid],panids) ,:iso] = true
            rlist, glist = genlists(dfd,iocc,idolocc,ipen) 
        else
            println("All good : ",mliftβ)
            break
        end
    end
end
# ---- SAVE RESULTS -----------
delete!(iocc, :exclude_vars)
delete!(idolocc, :exclude_vars)
delete!(ipen, :exclude_vars)
modelsDict = Dict()
modelsDict[:iocc] = iocc
modelsDict[:idolocc] = idolocc
modelsDict[:ipen] = ipen
modelsDict[:factors] = unpool(dfd) 
cols = vcat([:panid, :iso, :group, :buyer_pos_p1],cfg[:random_campaigns])
cols = vcat(cols,iocc[:finalvars],[iocc[:y_var], iocc[:logvar]])
cols = vcat(cols,idolocc[:finalvars],[idolocc[:y_var], idolocc[:logvar]])
cols = unique(vcat(cols,ipen[:finalvars],[ipen[:y_var], ipen[:logvar]]) )
# add required cold for Scoring agregation
cols = unique(vcat(cols,[:buyer_pre_p1, :buyer_pos_p1, :trps_pos_p1,:trps_pre_p1,:dol_per_trip_pre_p1,:dol_per_trip_pos_p1, :prd_1_net_pr_pre, :prd_1_net_pr_pos]))
save_dfd(root, dfd[(dfd[:iso].==false),cols] )
save_modelsDict(root, modelsDict)
    
    
    
    
    
    
# ******************************* MODELS **********************************************************
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
    
dx=OrderedDict()  # defaults to empty - so that consolidate will work with seq or MP
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
dfx = consolidateResults(modelsDict,dx)  
save_dfx(root,dfx)      #writetable(root*"/campaign.csv", dfx)                      

    
# ************************************* SCORING*************************************************************************************
        
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


        
        
# ************************************* OPT ****************************************************************************************    
@everywhere function CIOx(iDict::OrderedDict)
    v_ttl=2
    M = get(iDict,:M,NA)
    Mt = get(iDict,:Mt,NA)
    Mc = get(iDict,:Mc,NA)
    N = get(iDict,:N,NA)
    Nt = get(iDict,:Mt,NA)
    Nc = get(iDict,:Nc,NA)
    B1 = get(iDict,:B1,NA)
    B2 = get(iDict,:B2,NA)
    B3 = get(iDict,:B3,NA)
    SE1 = get(iDict,:SE1,NA)
    SE2 = get(iDict,:SE2,NA)
    SE3 = get(iDict,:SE3,NA)
    o_mean_score0 = get(iDict,:o_mean_score0,NA)
    o_mean_score1 = get(iDict,:o_mean_score1,NA)
    y_mean_score0 = get(iDict,:y_mean_score0,NA)
    y_mean_score1 = get(iDict,:y_mean_score1,NA)
    p_mean_score0 = get(iDict,:p_mean_score0,NA)
    p_mean_score1 =get(iDict,:p_mean_score1,NA)
    # ****** PValue ******
    SEsq=sqrt(SE1^2+SE2^2+SE3^2)
    Bsum=B1+B2+B3
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
    status = solve(m)
    zvalue=getobjectivevalue(m)
    iDict[:pval] =2.0 * ccdf(Normal(), abs(zvalue))
    println("z-value: ", string(zvalue)," --> p-value: ",string(iDict[:pval]))         
    # ****** UB-LB ******
    o_SE0 = get(iDict, :o_SE0, 0)
    y_SE0 = get(iDict, :y_SE0, 0)
    p_SE0 = get(iDict, :p_SE0, 0)
    o_B0 = get(iDict, :o_B0, 0)
    y_B0 = get(iDict, :y_B0, 0)
    p_B0 = get(iDict, :p_B0, 0)        
    SEsq=sqrt(SE1^2+SE2^2+SE3^2+o_SE0^2+y_SE0^2+p_SE0^2)        
    Bsum=(B1+B2+B3)-(o_B0+y_B0+p_B0)        
    ZDict = Dict("onetail_80_pct_intrvl" => 0.84 ,"onetail_90_pct_intrvl" => 1.28, "twotail_80_pct_intrvl" => 1.28, "twotail_90_pct_intrvl" => 1.65)
    AccArr= [ 0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01]
    ftxt="!!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!!!!!!!"
    failmsg(b::String,interval::String,meta::String) = println(ftxt*"   ",b," ",interval," := ",meta," == "*ftxt) 
    for (zscore_key,zscore) in  ZDict  
        # **********************  LB  ***************************
        pref="LB "*zscore_key[1:10]*" ("*iDict[:metakey]"):= "
        preflen=length(pref)
        dkey=Symbol(zscore_key*"_lb")
        for iAcc in AccArr
           print(pref," - "*string(iAcc)*", ",iDict)       ;pref=lpad("", preflen, " ")
           m=nothing
           m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))               
           ztot = Bsum-(zscore*SEsq)
           @variable(m, Bocc <= (B1-o_B0))
           @variable(m, Bdolocc <= (B2-y_B0))
           @variable(m, Bpen <= (B3-p_B0))
           #@NLobjective(m, Min, ((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))* 
           #                             ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))*
           #                              ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
            #                             )
            #                             -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))*
            #                               ((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))*
            #                               ((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
            #                               )
            #                             )
            #                          )
           @NLobjective(m, Min, (((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))*
                                         ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))*
                                         ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M))))
                                         /
                                         (((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))*
                                           ((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))*
                                           ((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))))-1)                     
                        )
           @constraint(m, (0.000000000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAcc)) 
           status = solve(m)
           mval=getobjectivevalue(m)
           println("lb_... ",mval)
           #mval_out=mval/(o_mean_score0*y_mean_score0*p_mean_score0)
           mval_out=mval
           if (status==:Optimal)&(mval_out!=Inf)&(mval_out!=-Inf)&(mval!=Inf)*(mval!=-Inf)
               iDict[dkey] = mval_out
               println("            : Confidence Interval : ", string(iDict[dkey]) )
               break      
           else
               iDict[dkey] = "error"
           end
        end
        if (iDict[dkey] == "error") 
             failmsg("LB",zscore_key[1:10],iDict[:metakey])
        end
        # **********************  UB  ***************************        
        pref="UB "*zscore_key[1:10]*" ("*iDict[:metakey]"):= "   
        dkey=Symbol(zscore_key*"_ub")
        for iAcc in AccArr
           print(pref," - "*string(iAcc)*", ")       ;pref=lpad("", preflen, " ")
           m=nothing
           m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))  
           ztot = Bsum+(zscore*SEsq)
           @variable(m, Bocc >= (B1-o_B0))
           @variable(m, Bdolocc >= (B2-y_B0))
           @variable(m, Bpen >= (B3-p_B0))
           #@NLobjective(m, Max, ( (((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))
            #                            * ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))
            #                            * ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
            #                             )
            #                            -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))
            #                             *((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))
            #                             *((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
            #                             )
            #                          )
            #                  )
            @NLobjective(m, Max, (((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))*
                                         ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))*
                                         ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M))))
                                         /
                                         (((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))*
                                           ((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))*
                                           ((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M))))
                                  )-1
                                 )                     
                         )
                        
           @constraint(m, (0.0000<= (((Bocc+Bpen+Bdolocc)-ztot))<= iAcc)) 
           status = solve(m)
           mval=getobjectivevalue(m)
           println("ub_... ",mval)
           #mval_out=mval/(o_mean_score0*y_mean_score0*p_mean_score0)
           mval_out=mval             
           if (status==:Optimal)&(mval_out!=Inf)&(mval_out!=-Inf)&(mval!=Inf)*(mval!=-Inf)
               iDict[dkey] = mval_out
               println("            : Confidence Interval : ", string(iDict[dkey]) )
               break      
           else
               iDict[dkey] = "error"
           end
        end
        if (iDict[dkey] == "error") 
             failmsg("UB",zscore_key[1:10],iDict[:metakey])
        end
    end
    return iDict
end                
# -------------------------------------------------------------------------------------------------------------------        
        
        

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
#ConfidenceIntervals(dfx)
#dfx[:adj_dod_effct] = ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
#dfx[(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),[:model,:adj_dod_effct,:onetail_80_pct_intrvl_lb,:onetail_80_pct_intrvl_ub]]
 
            
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
    #calcPValue_Opt(mDict)
    #CIs_O(mDict)
    CIOx(mDict)
    return mDict
end


function MP_ConfidenceIntervals(dfx::DataFrame)           
    @swarm :total MP_collectModels(dfx, "GLM")
    for ranef in unique(dfx[(dfx[:modelType].=="GLMM"),:ranef])
        for level in unique(dfx[(dfx[:modelType].=="GLMM")&(dfx[:ranef].==ranef),:parameter])
                k = Symbol(ranef*"~"*level)
                @swarm k MP_collectModels(dfx, "GLMM",ranef,level)
        end
    end
    sleep(20)
    while !iscompletew() println("Not Complete yet!"); sleep(5); end 
    mDict = takew(:total)
#    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:onetail_pval] = mDict[:onetail_pval]*100
#    dfx[(dfx[:model].=="dolhh")&(dfx[:modelType].=="GLM")&(dfx[:parameter].=="group"),:twotail_pval] = mDict[:twotail_pval]*100
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
    

  
function pvalues2Campaign(dfx::DataFrame)
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
pvalues2Campaign(dfx)
        

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
    #--
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
    #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :onetail_pval_to_campaign] =  (1-(dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :adj_pval] ./ 2) ) * 100
    #dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :twotail_pval_to_campaign] =  (1-dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].!="dolhh"), :adj_pval] ) * 100
#--
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


            
            
            
            
            
            
       
            
            
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- TESt CIs -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            
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
                Bplus=+(SE*zscore)  #B+(SE*zscore) 
                Bminus=-(SE*zscore) #B-(SE*zscore)
                score0=mean_score0*(Mc/M)
                score1=mean_score1*(Mt/M)
                Lb_pre = ( score1 - (score1*exp(-( Bminus )) ) )  +  ## ------ Lower Bound --------
                         ( (score0*exp(( Bminus )) )  -  score0    )
                row[Symbol(zscore_key*"_lb")] = ( Lb_pre/mean_score0 ) * 100
                Ub_pre =  (  score1 - ( score1*exp(-(Bplus)) )   )  +  ## ------ Upper Bound -------
                          (  ( score0*exp(Bplus) - score0) )
                row[Symbol(zscore_key*"_ub")] = ( Ub_pre/mean_score0 ) * 100   
                println("RAND CI $zscore_key ($zscore) LB:",row[Symbol(zscore_key*"_lb")]," ~~ UB : ", row[Symbol(zscore_key*"_ub")])
            end            
        end
    end
end
ConfIntrvl()
dfx[:adj_dod_effct] =  ((dfx[:adj_mean_score1] .- dfx[:adj_mean_score0]) ./ dfx[:adj_mean_score0] ) *100
dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].=="occ")&(dfx[:ranef].=="creative_theme")&(dfx[:parameter].=="Halloween"),[:adj_mean_score0,:adj_mean_score1,:Mt,:Mc,:N,:Nt,:Nc,:adj_coef,:adj_stderr,:adj_dod_effct,:onetail_80_pct_intrvl_lb,:onetail_80_pct_intrvl_ub]]

dfx[(dfx[:adj_dod_effct].<dfx[:onetail_80_pct_intrvl_lb])|(dfx[:adj_dod_effct].>dfx[:onetail_80_pct_intrvl_ub]),[:adj_mean_score0,:adj_mean_score1,:Mt,:Mc,:N,:Nt,:Nc,:adj_coef,:adj_stderr,:adj_dod_effct,:onetail_80_pct_intrvl_lb,:onetail_80_pct_intrvl_ub]]


