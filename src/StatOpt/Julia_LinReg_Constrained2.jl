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
using JuMP;
using Ipopt;
using DataFrames;
my_table = readtable("/home/iriadmin/g/jlib/src/StatOpt/logreg.csv")
( N, K1 ) = size( my_table )
m = Model(solver=IpoptSolver())  # sum( -1 * y .* (x*beta) + log(1+exp((x*beta))) )
@defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 ) 
@NLobjective( m, Min, sum{ ( my_table[i,1] - sum{ my_table[i,j] * b[ j - 1] , j in 2:K1} )^2 , i in 1:N} )
statusM = solve( m )
df = DataFrame( beta = getvalue( b ) )

# ----------------------- LOGISTIC --------------------------------


using JuMP, Ipopt, DataFrames
my_table = readtable( "/home/iriadmin/g/JStack/src/tst/newdata.csv");
( N, K1 ) = size( my_table );
m = Model(solver=IpoptSolver());
@defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 )
@NLobjective( m, Min, sum{ ( my_table[i,1] - sum{ my_table[i,j] * b[ j - 1] , j in 2:K1} )^2 , i in 1:N} )
statusM = solve( m )
df = DataFrame( beta = getvalue( b ) )

dfd = readtable( "/home/iriadmin/g/JStack/src/tst/logreg.csv")
dfd=dfd[setdiff(names(dfd),[:x0])]
mm = ModelMatrix(ModelFrame(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd))
g = glm(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd,  Bernoulli())   #m = glm(y ~ 1 + x1 + x2 + x3 + x4 + x6 + x7, dfd, Poisson())
a = convert(Array,dfd) #array(dfd[:, 1:end])
B = zeros(length(mm.m[1,2:end]))
x = a[:,2:end]
y = a[:,1]
# sum( -1 * y .* (x*beta) + log(1+exp((x*beta))) )

m = Model(solver=IpoptSolver());
@variable(m, B[i=1:6], start=0)

@defVar( m, 0.0 <= b[1:(K1-1)] <= 100.0 )
@NLobjective( m, Min, sum{ ( my_table[i,1] - sum{ my_table[i,j] * b[ j - 1] , j in 2:K1} )^2 , i in 1:N} )
statusM = solve( m )
df = DataFrame( beta = getvalue( b ) )


