x <- rnorm( 100, mean = 0, sd = 1 )
y <- rnorm( 120, mean = 0, sd = 1 )

ks.test( x, y )

q_grid <- seq( 0, 1, by = 0.01 )

q0 <- quantile( x, q_grid )
q1 <- quantile( y, q_grid )


plot( q_grid, q0, type = 'l', col ='red' )
lines( q_grid, q1, col = 'blue' )
