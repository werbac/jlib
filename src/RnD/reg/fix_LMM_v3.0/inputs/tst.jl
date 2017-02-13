


using DataFrames, MixedModels,StatLib
dfx = readtable("/mapr/mapr04p/analytics0001/analytic_users/jmdl/v3.0/inputs/mctargetgrus00115.csv")
dfs=readtable("/mapr/mapr04p/analytics0001/analytic_users/jmdl/v3.0/sas_output/fixed1_grus00115.csv")
#mctargetgrus00147.csv  mctargetgrus00247.csv  mctargetgrus00256.csv 
#Sample 1:for ppgid :115 region: US outlet: GR    
f = logvolume ~ 1 + logregprice + logpriceindex + displayonly + featuredisplay + featureonly + logseason + absprice19 + lift19 + absprice220 + lift220 + absprice178 + lift178 + absprice132 + lift132 + absprice131 + lift131 + absprice191 + lift191 + absprice95 + lift95 + absprice140 + lift140 + absprice139 + lift139 + absprice31 + lift31 + absprice233 + lift233 + absprice289 + lift289 + absprice275 + lift275 + absprice236 + lift236 + absprice230 + lift230 + absprice235 + lift235 + absprice246 + lift246 + absprice265 + lift265 + absprice290 + lift290 + absprice239 + lift239 + absprice7 + lift7 + absprice305 + lift305 + absprice322 + lift322 + absprice302 + lift302 + absprice299 + lift299 + absprice229 + lift229 + absprice309 + lift309 + absprice231 + lift231 + absprice298 + lift298 + absprice311 + lift311 + 
( ( 0 + logregprice ) | store ) + 
( ( 0 + logpriceindex ) | store ) + 
( ( 0 + displayonly ) | store ) + 
( ( 0 + featuredisplay ) | store ) + 
( ( 0 + featureonly ) | store ) + 
( ( 0 + logseason ) | store ) 

lowercase!(dfx)
for col in[:weekending, :outlet, :region] dfx[col] = pool(dfx[col]) end

m=fit!(lmm(f, dfx), true)



mx = lmm(f, dfx)
reml!(mx,true)
m = fit(mx)

@time m1 = fit!(lmm(Y ~ 1 + (1|G) + (1|H), ml1m))


fit!(lmm(f, dfx), true)





gmm1 = fit!(glmm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk]), verbose=false, fast=true)
glm(f, dfd[(dfd[:buyer_pos_p1].==1),cols] , m[:dist], m[:lnk])




modelREML = lmm(f, dfx)
reml!(modelREML,true)
lmeModel = fit(modelREML)




    
    

Sample 2:for ppgid :247 region: US outlet: GR
logvolume ~ 1 + logregprice + logpriceindex + displayonly + featuredisplay + featureonly + logseason + absprice263 + lift263 + absprice258 + lift258 + absprice259 + lift259 + absprice254 + lift254 + absprice261 + lift261 + absprice265 + lift265 + absprice248 + lift248 + absprice260 + lift260 + absprice266 + lift266 + absprice255 + lift255 + absprice276 + lift276 + absprice271 + lift271 + absprice289 + lift289 + absprice228 + lift228 + absprice236 + lift236 + absprice242 + lift242 + absprice279 + lift279 + absprice282 + lift282 + absprice239 + lift239 + absprice230 + lift230 + absprice311 + lift311 + absprice10 + lift10 + lift189 + absprice115 + lift115 + absprice20 + lift20 + absprice33 + lift33 + absprice293 + lift293 + absprice12 + lift12 + absprice298 + lift298 + lift45 + ( ( 0 + logregprice ) | store ) + ( ( 0 + logpriceindex ) | store ) + ( ( 0 + displayonly ) | store ) + ( ( 0 + featuredisplay ) | store ) + ( ( 0 + featureonly ) | store ) + ( ( 0 + logseason ) | store ) 

Sample 3: for ppgid :147 region: US outlet: GR
logvolume ~ 1 + logregprice + logpriceindex + displayonly + logseason + absprice88 + lift88 + absprice48 + lift48 + absprice95 + lift95 + absprice194 + lift194 + absprice89 + lift89 + absprice227 + lift227 + absprice49 + lift49 + absprice65 + lift65 + absprice225 + lift225 + absprice33 + lift33 + absprice235 + lift235 + absprice272 + lift272 + absprice276 + lift276 + absprice288 + lift288 + absprice278 + lift278 + absprice265 + lift265 + absprice274 + lift274 + absprice228 + lift228 + absprice322 + lift322 + absprice13 + lift13 + absprice327 + lift327 + absprice236 + lift236 + absprice305 + lift305 + absprice299 + lift299 + absprice248 + lift248 + absprice324 + lift324 + absprice275 + lift275 + absprice289 + lift289 + absprice309 + lift309 + absprice266 + lift266 + ( ( 0 + logregprice ) | store ) + ( ( 0 + logpriceindex ) | store ) + ( ( 0 + displayonly ) | store ) + ( ( 0 + logseason ) | store ) 

Sample 4: for ppgid :256 region: US outlet: GR
logvolume ~ 1 + logregprice + logpriceindex + featuredisplay + featureonly + logseason + lift265 + absprice248 + lift248 + absprice278 + lift278 + lift229 + lift238 + absprice230 + lift230 + lift281 + lift239 + absprice246 + lift246 + lift236 + absprice275 + lift275 + absprice274 + lift274 + lift33 + lift178 + lift308 + absprice95 + lift95 + absprice20 + lift20 + absprice300 + lift300 + absprice12 + lift12 + absprice305 + lift305 + absprice107 + lift107 + absprice199 + lift199 + absprice106 + lift106 + absprice7 + absprice3 + lift3 + lift189 + lift298 + absprice297 + lift297 + absprice303 + lift303 + absprice5 + lift5 + ( ( 0 + logregprice ) | store ) + ( ( 0 + logpriceindex ) | store ) + ( ( 0 + featuredisplay ) | store ) + ( ( 0 + featureonly ) | store ) + ( ( 0 + logseason ) | store ) 

