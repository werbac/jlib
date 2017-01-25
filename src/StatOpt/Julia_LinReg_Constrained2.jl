using JuMP;
using Ipopt;
using DataFrames;
my_table = readtable("/home/iriadmin/g/jlib/src/StatOpt/newdata.csv")
( N, K1 ) = size( my_table )
m = Model(solver=IpoptSolver())
@defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 ) 
@NLobjective( m, Min, sum{ ( my_table[i,1] - sum{ my_table[i,j] * b[ j - 1] , j in 2:K1} )^2 , i in 1:N} )
statusM = solve( m )
df = DataFrame( beta = getvalue( b ) )



# ----------------------- GT LOGISTIC --------------------------------
using JuMP,Ipopt,DataFrames, GLM
m=nothing
m = Model(solver=IpoptSolver())  # sum( -1 * y .* (x*beta) + log(1+exp((x*beta))) )
dfd = readtable("/home/iriadmin/g/jlib/src/StatOpt/logreg.csv")
x = Array(dfd[:,2:end])
y = Array(dfd[:,1])
(rows, cols) = size(x)
g = glm(y ~ -1 + x0 + x1 + x2 + x3 + x4 + x6 + x7, dfd,  Bernoulli(), LogitLink()) # X0 s the replacement for intercept - otherwise glm will generate a [1...] arr
@variable(m, B[z=1:cols], start=0)  # getvalue(B)  @defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 )
@NLobjective( m, Min, sum{ ( -1*y[r]*sum{x[r,c]*B[c],c in 1:cols} ) + log(1+exp(sum{x[r,cc]*B[cc],cc in 1:cols})) ,r in 1:rows} )  
#..... sum( -1 * y .* (x*B) + log(1+exp((x*B))) )

d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:Grad,:HessVec,:Hess])
statusM = solve( m )  #statusM = solve( m, IpoptOptions=[("tol",1e-6)] )
df = DataFrame( beta = getvalue( B ) )

d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:Grad,:HessVec,:Hess])
hess = MathProgBase.hesslag_structure(d)
diag(hess)
beta_se <- sqrt( diag( var_cov_m ) )

# -----------------------------solution ------------------------
using JuMP,Ipopt,DataFrames, GLM, ForwardDiff
m=nothing
m = Model(solver=IpoptSolver())  # sum( -1 * y .* (x*beta) + log(1+exp((x*beta))) )
dfd = readtable("/home/iriadmin/g/JLib/src/StatOpt/logreg.csv")
x = Array(dfd[:,2:end])
y = Array(dfd[:,1])
(rows, cols) = size(x)
g = glm(y ~ -1 + x0 + x1 + x2 + x3 + x4 + x6 + x7, dfd,  Bernoulli(), LogitLink())
@variable(m, B[z=1:cols], start=0)  # getvalue(B)  @defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 )
@NLobjective( m, Min, sum{ ( -1*y[r]*sum{x[r,c]*B[c],c in 1:cols} ) + log(1+exp(sum{x[r,cc]*B[cc],cc in 1:cols})) ,r in 1:rows} )

statusM = solve( m )
df = DataFrame( beta = getvalue( B ) )

function mynllf(beta::Vector)
	nll = 0.0
	for r in 1:rows
		dbl_part1 = 0.0
		for c in 1:cols
			dbl_part1 = dbl_part1 + x[r,c] * beta[c]		
		end
		nll = nll -1.0 * y[r] * dbl_part1 + log(1+exp(dbl_part1))
	end
	return nll
end

sqrt(diag(inv(  ForwardDiff.hessian(mynllf,getvalue(B))  )))
getvalue(B)


# ----------------------- POISON --------------------------------
using JuMP,Ipopt,DataFrames, GLM
m=nothing
m = Model(solver=IpoptSolver()) 
dfd = readtable("/home/iriadmin/g/JLib/src/StatOpt/ZD/poison_gamma/y.csv")
names!(dfd,[:y,:x1, :x2, :x3,:x4,:x5])
x = Array(dfd[:,2:end])
y = Array(dfd[:,1])
(rows, cols) = size(x)
g = glm(y ~ -1 + x1 + x2 + x3 + x4 +x5 , dfd,  Poisson(), LogLink()) 
@variable(m, B[z=1:cols], start=0)
@NLobjective( m, Min, sum{ exp(sum{x[r,cc]*B[cc],cc in 1:cols}) - (y[r]*sum{x[r,c]*B[c],c in 1:cols}) ,r in 1:rows} )
statusM = solve( m )
df = DataFrame(beta=getvalue(B))
mynllf(beta::Vector) = sum(exp(x*beta) .- ((x*beta).*y))
sqrt(diag(inv(ForwardDiff.hessian(mynllf,getvalue(B)))))



# ----------------------- GAMMA --------------------------------
using JuMP,Ipopt,DataFrames, GLM
m=nothing
m = Model(solver=IpoptSolver()) 
dfd = readtable("/home/iriadmin/g/JLib/src/StatOpt/ZD/poison_gamma/g.csv")
names!(dfd,[:y,:x1, :x2, :x3,:x4,:x5])
x = Array(dfd[:,2:end])
y = Array(dfd[:,1])
(rows, cols) = size(x)
g = glm(y ~ -1 + x1 + x2 + x3 + x4 +x5 , dfd,  Gamma(), LogLink()) 
@variable(m, B[z=1:cols], start=0)
@NLobjective( m, Min, sum{ exp(-1.0*sum{x[r,cc]*B[cc],cc in 1:cols}) + (sum{x[r,c]*B[c],c in 1:cols}) ,r in 1:rows} )
statusM = solve( m )
df = DataFrame(beta=getvalue(B))
mynllf(beta::Vector) = sum(exp(-1.0*(x*beta)) .+ ((x*beta)))
sqrt(diag(inv(ForwardDiff.hessian(mynllf,getvalue(B)))))






str_path <- "C:\\Users\\mdzhd\\Documents\\Bayesian_Approach_with_Lampros\\"

myd <- read.csv( paste0( str_path, "gamma_regression_ex2.csv" ) )

myglm <- glm( y ~ x1 + x2 + x3, data = myd, family = Gamma( link = log ) )

x <- as.matrix( myd[, -c(1)] )
y <- as.matrix( myd[, 1] )

NegLogLike <- function( beta )
{
	xbeta <- x %*% beta
	Li <- y * exp( -1.0 * xbeta ) +  xbeta
	sum( Li ) 
}

beta_initial <- as.matrix( rep( 0, ncol( x ) ), ncol(x), 1 )

myoptim <- optim( beta_initial, f = NegLogLike, g = NULL, method = "BFGS", hessian = TRUE)

u <- exp( x %*% as.matrix( myoptim$par , ncol( x ), 1 ) )

v <- ( y - u )/u 

v2 <- v * v

phi <- sum( v2 ) / ( nrow( x ) - ncol(x ) )

myhessian <- myoptim$hessian / phi

var_cov_m <- solve( myhessian )

beta_se <- sqrt( diag( var_cov_m ) )

zvalue <- myoptim$par / beta_se

pvalue <- 2 * ( 1 - pnorm(abs( zvalue ) ) )

coeff_outputs <- cbind( myoptim$par, beta_se, zvalue, round( pvalue, 6 ) )

colnames( coeff_outputs ) <- c( "Estimates", "Std. Error", "z value", "Pr(>|z|)" )
rownames( coeff_outputs ) <- c()

print( coeff_outputs )
summary( myglm )$coefficients






