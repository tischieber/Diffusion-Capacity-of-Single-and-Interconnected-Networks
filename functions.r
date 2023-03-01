library(igraph)
library('Matrix')
require('MASS')


ent_mat<-function(s){
  
  quem=which(s>0)
  
  s[quem]=-s[quem]*log(s[quem])
  
  return(s)
  
}


J_mat<-function(p,q){
  
  valor=colSums(t(ent_mat((p+q)/2)-(ent_mat(p)+ent_mat(q))/2))
  
  valor[which(valor<=0)]=0
  
  return(sqrt(valor/log(2)))
  
  
}



Dist_cum<-function(p,q){
  
  n=max(c(max(which(colSums(p)>0)),max(which(colSums(q)>0))))
  
  p=p[,1:n]
  
  q=q[,1:n]
  
  r=J_mat(p,q)
  
  i=1
  
  while(i<n){
    
    p[,2]=p[,1]+p[,2]
    
    q[,2]=q[,1]+q[,2]
    
    k=length(p[1,])
    
    p=p[,2:k]
    
    q=q[,2:k]
    
    i=i+1
    
    r=r+J_mat(p,q)
    
    #print(i)
    
  }
  
  return(r)
  
}


DgDw<-function(dg,dw){
  
  N=length(dg[,1])
  
  #dw=shortest.paths(g,mode=c("out"))
  #options(warn=-1)#turn off warnings
  #dg=shortest.paths(g,,mode=c("out"),algorithm=c("unweighted"))#geodesic distance
  #options(warn=0)#turn on warnings
  
  quem=c(which(dg==Inf),which(dw==Inf))
  
  quem=intersect(quem,quem)
  
  dg[quem]=N
  
  dw[quem]=N
  
  v=dg/dw
  
  v[is.nan(v)]=0
  
  
  #ganhos
  
  ganho=which(as.numeric(v<1)*v>0)
  
  
  if(length(ganho)>0){
    
    coluna=dg[ganho]
    
    linha=N-(ceiling(ganho/N)*N-ganho)
    
    gganho=graph(c(t(matrix(c(linha,coluna),ncol=2))),n=N)
    
    E(gganho)$weight=1-v[ganho]
    
    
    
  }
  
  if(length(ganho)==0){
    
    gganho=graph(c(),n=N)
    
    E(gganho)$weight=1
    
  }
  
  
  #perdas
  
  ganho=which(as.numeric(v>1)*v>0)
  
  
  if(length(ganho)>0){
    
    coluna=dg[ganho]
    
    linha=N-(ceiling(ganho/N)*N-ganho)
    
    gperda=graph(c(t(matrix(c(linha,coluna),ncol=2))),n=N)
    
    E(gperda)$weight=1-1/v[ganho]
    
  }
  
  if(length(ganho)==0){
    
    gperda=graph(c(),n=N)
    
    E(gperda)$weight=1
    
  }
  
  
  
  #estavel
  
  v[which(v>1)]=1/v[which(v>1)]
  
  ganho=which(as.numeric(v>0)*v>0)
  
  coluna=dg[ganho]
  
  linha=N-(ceiling(ganho/N)*N-ganho)
  
  gest=graph(c(t(matrix(c(linha,coluna),ncol=2))),n=N)
  
  E(gest)$weight=v[ganho]
  
  return(list(gperda,gest,gganho))
  
  
}


reorganizapdf<-function(a){
  
  n=length(V(a[[1]]))
  
  resp=Matrix(0,ncol=3*n,nrow=n,sparse=TRUE)
  
  seq=seq(1,3*n-2,3)
  
  resp[1:n,seq]=resp[1:n,seq]+get.adjacency(a[[1]],attr="weight")
  
  seq=seq(2,3*n-1,3)
  
  resp[1:n,seq]=resp[1:n,seq]+get.adjacency(a[[2]],attr="weight")
  
  seq=seq(3,3*n,3)
  
  resp[1:n,seq]=resp[1:n,seq]+get.adjacency(a[[3]],attr="weight")
  
  resp[1:n,3*n]=resp[1:n,3*n]+resp[1:n,2*n]
  
  resp[1:n,2*n]=resp[1:n,2*n]*0
  
  return(resp/(n-1))
  
}


refresco<-function(arquivo){
  #retorna um refresh de um dado arquivo localizado em um dado diretorio guardando as variaveis
  
  
  makeActiveBinding("refresh", function() { shell("Rgui"); q("no") }, .GlobalEnv)
  
  s=sprintf("Rscript %s",arquivo)
  
  makeActiveBinding("refresh", function() { system(s); q("no") }, .GlobalEnv)
  
  refresh
  
  
  
  
}

dwinter<-function(dw1,dw2){
  
  n=matrix(length(dw1[,1]))
  
  r=matrix(0,ncol=n,nrow=n)
  
  for(i in (1:n)){
    
    r[i,]=sapply(data.frame(dw1[i,]+dw2),min)
    
  }
  
  return(r)
}


res_layer<-function(dg,dw,layers,i){
  
  tam=length(layers)
  
  
  
  pdf1=DgDw(dg[c(layers[i][[1]]),c(layers[i][[1]])],dw[c(layers[i][[1]]),c(layers[i][[1]])])
  
  pdf1=reorganizapdf(pdf1)
  
  ref=Matrix(0,ncol=ncol(pdf1),nrow=nrow(pdf1),sparse=TRUE)
  
  ref[,1]=ref[,1]+1
  
  valor=Dist_cum(pdf1,ref)
  
  for(k in setdiff((1:length(layers)),i)){
    
    dw12=dwinter(dw[layers[i][[1]],layers[k][[1]]],dw[layers[k][[1]],layers[i][[1]]])
    
    pdf2=DgDw(dg[layers[i][[1]],layers[i][[1]]],dw12)
    
    pdf2=reorganizapdf(pdf2)
    
    valor=valor+Dist_cum(pdf2,ref)/(tam-1)
    
  }
  
  valor=2/valor
  
  return(valor)
  
}



# g3 <- graph.graphdb("si6_r005_s100.B99") database de graph isomorphism 
#write(


entropia<-function(a){
  
  
  a[a>0]=(-a[a>0]*log(a[a>0]))
  
  return(a)
  
  
  
  
}

sldc<-function(g){
  
  
  dw=shortest.paths(g,mode=c("out"))
  options(warn=-1)#turn off warnings
  dg=shortest.paths(as.undirected(g),mode=c("out"),algorithm=c("unweighted"))#geodesic distance
  options(warn=0)#turn on warnings
  
  pdf1=DgDw(dg,dw)
  
  
  pdf1=reorganizapdf(pdf1)
  
  ref=Matrix(0,ncol=ncol(pdf1),nrow=nrow(pdf1),sparse=TRUE)
  
  ref[,1]=ref[,1]+1
  
  vl1=(1/((Dist_cum(pdf1,ref))))
  
  return(vl1)
  
}

mldc<-function(g,layers){
  
  
  dw=shortest.paths(g,mode=c("out"))
  options(warn=-1)#turn off warnings
  dg=shortest.paths(as.undirected(g),mode=c("out"),algorithm=c("unweighted"))#geodesic distance
  options(warn=0)#turn on warnings
  v=c()
  
  DC=c()
  
  for(i in 1:length(layers)){
    
    dc=res_layer(dg,dw,layers,i)
    
    v=c(v,rep(i,length(dc)))
    
    DC=c(DC,dc)
    
  }
  
  return(data.frame(node=1:length(V(g)),layer=v,diffusion_capacity=DC))
  
}


