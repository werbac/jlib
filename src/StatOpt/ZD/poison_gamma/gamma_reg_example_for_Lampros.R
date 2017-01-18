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






