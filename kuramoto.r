
require("pracma")

alpha=1



  
w=runif(25,min=0,max=5)*0 #random initial values for oscillators
  

  
  
  
  G <- make_lattice( c(5,5) )
  

  G=graph.adjacency(A,mode=c("directed"))  
  
  elos=get.edgelist(G)
  

  N=length(V(G))
  
  x=matrix(runif(N),ncol=1,nrow=N)*5
  
  
  t=0
  
  resposta_final=c()
  
  h=0.001
  
  resp=c()
  
  H=G    
  
  A=get.adjacency(H)
  
  
  
  lambda=0.01
  
  time=10
  
  
  for(interv in seq(lambda,time,0.01)){
    
    
    
    acop=c()
    for(i in (1:nrow(A))){
      
      acop=c(acop,(A[i,]*(sin(x-x[i]))))
      
    }
    
    acop=matrix(acop,nrow=25,byrow=TRUE)
    
    acop=acop*lambda/N
    
    k1=w+colSums(t(acop)) 
    
    acop=c()
    for(i in (1:nrow(A))){
      
      acop=c(acop,(A[i,]*(sin(x+k1*h/2-x[i]-k1[i]*h/2))))
      
    }
    
    acop=matrix(acop,nrow=25,byrow=TRUE)
    
    acop=acop*lambda/N
    
    k2=w+colSums(t(acop))
    
    
    acop=c()
    for(i in (1:nrow(A))){
      
      acop=c(acop,(A[i,]*(sin(x+k2*h/2-x[i]-k2[i]*h/2))))
      
    }
    
    acop=matrix(acop,nrow=25,byrow=TRUE)
    
    acop=acop*lambda/N
    
    k3=w+colSums(t(acop))
    
    
    acop=c()
    for(i in (1:nrow(A))){
      
      acop=c(acop,(A[i,]*(sin(x+k3*h-x[i]-k3[i]*h))))
      
    }
    
    acop=matrix(acop,nrow=25,byrow=TRUE)
    
    acop=acop*lambda/N
    
    k4=w+colSums(t(acop))
    

   x=x+(h/6)*(k1+2*k2+2*k3+k4)
    
    resp=c(resp,x)
    
  }
  
  x 
  
  #x is the final distribution

  