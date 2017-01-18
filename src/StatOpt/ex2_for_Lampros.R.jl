library( Matrix )

str_path <- "C:\\CBIG\\ZD\\"

myd <- read.csv( paste0( str_path, "logistic_regression_ex2.csv" ) )

x <- as.matrix( myd[, -c( 1 )] )
y <- myd[, c(1)]

## make sure its rank = # columns
ncol( x )
rankMatrix( x )

myglm <- glm( y ~ -1 + x, family = binomial( link = logit ) )

myglm_coeff <- summary( myglm )$coefficients

rownames( myglm_coeff ) <- c( paste0( "beta", seq( 1, 3 ) ), paste0( "creative_break_", seq( 1, 4 ) ) )

print( myglm_coeff )

## using optim to reproduce the estimates

Loglikelihood <- function( beta )
{
	xbeta <- x %*% beta
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






