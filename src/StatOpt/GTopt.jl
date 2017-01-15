using JuMP;
using Ipopt;
using DataFrames;
my_table = readtable( "/home/iriadmin/newdata.csv");
( N, K1 ) = size( my_table );
m = Model(solver=IpoptSolver());
@defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 )
@NLobjective( m, Min, sum{ ( my_table[i,1] - sum{ my_table[i,j] * b[ j - 1] , j in 2:K1} )^2 , i in 1:N} )
statusM = solve( m )

beta = getvalue(b)


# ---------------------------- R ---------------------------------------------------------
library( Matrix )
myd <- read.csv( "C:\\CBIG\\ZD\\logistic_regression_ex2.csv" ) 
x <- as.matrix( myd[, -c( 1 )] )
y <- myd[, c(1)]
## make sure its rank = # columns
ncol(x)
rankMatrix(x)
myglm <- glm( y ~ -1 + x, family = binomial( link = logit ) )
myglm_coeff <- summary( myglm )$coefficients
rownames( myglm_coeff ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )
print( myglm_coeff )
## using optim to reproduce the estimates
Loglikelihood <- function( beta )
{ xbeta <- x %*% beta
  Li <- -1 * y * xbeta + log( 1 + exp( xbeta )) 
  sum( Li )
}
beta_start_value <- as.matrix( rep( 0.0, ncol( x ) ), ncol(x), 1 )
myoptim <- optim( beta_start_value, fn = Loglikelihood, gr = NULL, method = "BFGS", hessian = TRUE )
var_cov_m <- solve( myoptim$hessian )
beta_se <- sqrt( diag( var_cov_m ) )
zvalue <- myoptim$par / beta_se
pvalue <- 2 * ( 1 - pnorm(abs( zvalue ) ) )
coeff_outputs <- cbind( myoptim$par, beta_se, zvalue, round( pvalue, 6 ) )
colnames( coeff_outputs ) <- c( "Estimates", "Std. Error", "z value", "Pr(>|z|)" )
rownames( coeff_outputs ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )
print( coeff_outputs )
print( summary( myglm )$coefficients ) 


#######################################################################

## this is the residual deviance
2 * myoptim$value

########## since colmn 4 to column 7 are indicators, assuming they are the compaign breaks :)
##### so we want to add penalty function to make the estimates >= 0, using the penalty that Zhiyuan recommended

## assuming the penalty function is : 

a <- c( 1000, 9, 60 )

Penalized_Loglikelihood <- function( beta )
{
	xbeta <- x %*% beta
	Li <- -1 * y * xbeta + log( 1 + exp( xbeta ))
	## constraining beta[3] to beta[7]
	myPenValue <- 0
	for ( i in 4 : 7 )
	{
		myPenValue <- myPenValue + a[1] * exp( -1.0 * a[2] * exp( a[3] * beta[ i ] ) )
	}
	sum( Li ) + myPenValue
}

beta_start_value <- as.matrix( rep( 0.0, ncol( x ) ), ncol(x), 1 )

myoptim_Penalized <- optim( beta_start_value, fn = Penalized_Loglikelihood, gr = NULL, method = "BFGS", hessian = TRUE)

myoptim_Penalized

var_cov_m <- solve( myoptim_Penalized$hessian )

beta_se <- sqrt( diag( var_cov_m ) )

zvalue <- myoptim_Penalized$par / beta_se

pvalue <- 2 * ( 1 - pnorm(abs( zvalue ) ) )

coeff_outputs <- cbind( myoptim_Penalized$par, beta_se, zvalue, round( pvalue, 6 ) )

colnames( coeff_outputs ) <- c( "Estimates", "Std. Error", "z value", "Pr(>|z|)" )
rownames( coeff_outputs ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )

print( coeff_outputs )


# ------------------------------ GMT ------------------------------------
using JuMP, DataFrames
# Use nonlinear optimization to compute the maximum likelihood estimate (MLE) of the parameters of a normal distribution  aka the sample mean and variance
n = 1000
data = randn(n)
m = Model()
@variable(m, μ, start = 0.0)
@variable(m, σ >= 0.0, start = 1.0)
@NLobjective(m, Max, (n/2)*log(1/(2π*σ^2))-sum((data[i]-μ)^2 for i=1:n)/(2σ^2))
solve(m)
println("μ = ", getvalue(μ))
println("mean(data) = ", mean(data))
println("σ^2 = ", getvalue(σ)^2)
println("var(data) = ", var(data))
println("MLE objective: ", getobjectivevalue(m))
# constrained MLE?
@NLconstraint(m, μ == σ^2)
solve(m)
println("\nWith constraint μ == σ^2:")
println("μ = ", getvalue(μ))
println("σ^2 = ", getvalue(σ)^2)
println("Constrained MLE objective: ", getobjectivevalue(m))
            
            
            
# ------------------------------ GMT ------------------------------------
# ------------------------------ GMT ------------------------------------
# ------------------------------ GMT ------------------------------------
using JuMP, DataFrames
# Use nonlinear optimization to compute the maximum likelihood estimate (MLE) of the parameters of a normal distribution  aka the sample mean and variance
n = 1000
data = randn(n)
m = Model()
@variable(m, μ, start = 0.0)
@variable(m, σ >= 0.0, start = 1.0)
@NLobjective(m, Max, (n/2)*log(1/(2π*σ^2))-sum((data[i]-μ)^2 for i=1:n)/(2σ^2))
solve(m)
println("μ = ", getvalue(μ),"\nmean(data) = ",mean(data), "\nσ^2 = ",getvalue(σ)^2,"\nvar(data) = ",var(data),"\nMLE objective: ",getobjectivevalue(m) )
@NLconstraint(m, μ == σ^2) # constrained MLE?
solve(m)
println("\nWith constraint μ == σ^2:","μ = ", getvalue(μ), "σ^2 = ", getvalue(σ)^2, "Constrained MLE objective: ", getobjectivevalue(m) )

# ------------ GMT HERE HERE HERE HERE HERE ---------------
using JuMP, DataFrames, GLM, NLsolve
dfd = readtable( "/home/iriadmin/g/JStack/src/tst/logreg.csv")
dfd=dfd[setdiff(names(dfd),[:x0])]
mm = ModelMatrix(ModelFrame(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd))

g = glm(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd,  Bernoulli())   #m = glm(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd, Poisson())

a = convert(Array,dfd) #array(dfd[:, 1:end])
B = zeros(length(mm.m[1,2:end]))
   
x = a[:,2:end]
y = a[:,1]

# sum( -1 * y .* (x*beta) + log(1+exp((x*beta))) )
                
function LogLHd(beta::Array{Float64})
    xbeta = x*beta
    Li = -1 * y .* xbeta + log(1+exp(xbeta)) 
                    println(sum(Li) )                    
    sum(Li)              
end
                
#ll(B) = B*mm.m[:,2:end]
#m = Model()
#@variable(m, μ, start = 0.0)
                
using Optim
zz=optimize(LogLHd, B, LBFGS())
                
                
beta_start_value <- as.matrix( rep( 0.0, ncol( x ) ), ncol(x), 1 )
myoptim <- optim( beta_start_value, fn = Loglikelihood, gr = NULL, method = "BFGS", hessian = TRUE )
var_cov_m <- solve( myoptim$hessian )
                    
## make sure its rank = # columns
#ncol(x)
#rankMatrix(x)
#myglm <- glm( y ~ -1 + x, family = binomial( link = logit ) )                               
m = Model()
@variable(m, μ, start = 0.0)
@variable(m, σ >= 0.0, start = 1.0)
                @NLobjective(m, Max, (n/2)*log(1/(2π*σ^2))-sum((dfd[i]-μ)^2 for i=1:n)/(2σ^2))
solve(m)

# ----
using NLsolve
function f!(x, fvec)
       Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*x[1])))*(Mt/M))   )  +
                    (   (mean_score0*exp((B-(SE*x[1])))*(Mc/M))    -   (mean_score0*(Mc/M))    )
                    Lb = Lb_pre/mean_score0
                   fvec[1] = Lb
end
r=nlsolve(f!,[0.1])    #println(r)
# ---------------------------- R ---------------------------------------------------------
library( Matrix )
myd <- read.csv( "C:\\CBIG\\ZD\\logistic_regression_ex2.csv" ) 
x <- as.matrix( myd[, -c( 1 )] )
y <- myd[, c(1)]
## make sure its rank = # columns
ncol(x)
rankMatrix(x)
myglm <- glm( y ~ -1 + x, family = binomial( link = logit ) )
myglm_coeff <- summary( myglm )$coefficients
rownames( myglm_coeff ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )
print( myglm_coeff )
## using optim to reproduce the estimates
Loglikelihood <- function( beta )
{ xbeta <- x %*% beta
  Li <- -1 * y * xbeta + log( 1 + exp( xbeta )) 
  sum( Li )
}
beta_start_value <- as.matrix( rep( 0.0, ncol( x ) ), ncol(x), 1 )
myoptim <- optim( beta_start_value, fn = Loglikelihood, gr = NULL, method = "BFGS", hessian = TRUE )
var_cov_m <- solve( myoptim$hessian )
beta_se <- sqrt( diag( var_cov_m ) )
zvalue <- myoptim$par / beta_se
pvalue <- 2 * ( 1 - pnorm(abs( zvalue ) ) )
coeff_outputs <- cbind( myoptim$par, beta_se, zvalue, round( pvalue, 6 ) )
colnames( coeff_outputs ) <- c( "Estimates", "Std. Error", "z value", "Pr(>|z|)" )
rownames( coeff_outputs ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )
print( coeff_outputs )
print( summary( myglm )$coefficients ) 

# ------------------------------ GMT END ------------------------------------
# ------------------------------ GMT END ------------------------------------
# ------------------------------ GMT END ------------------------------------
