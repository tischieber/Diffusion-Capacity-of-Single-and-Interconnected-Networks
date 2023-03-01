#creating a network

G=graph(c(1,2,1,3,2,3,3,4))

E(G)$weight=c(0.5,1,2,2)

G=as.undirected(G)

elos=get.edgelist(G)

L=-graph.laplacian(G)

auto=eigen(L)

inversaP=solve(auto$vectors)

N=length(V(G))

#temperature initial condition (in this case vertice 4 has 10 degrees and others zero)

X0=matrix(rep(0,N),ncol=1,nrow=N)

X0[4]=10

#

C=inversaP%*%X0

#choose a time

t=5

X=auto[["vectors"]]%*%diag(exp(auto$values*t))%*%C

#X gives the temperatures at t=5