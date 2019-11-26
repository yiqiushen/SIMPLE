
library(rARPACK)
library(pracma)
library(MASS)
library(foreach)
library(doParallel)

getG=function(i,j){
  Sigma_2=L[,i,]%*%diag(sigma_hat[i,])%*%t(L[,i,])+L[,j,]%*%diag(sigma_hat[j,])%*%t(L[,j,])-sigma_hat[i,j]*L[,i,j]%*%t(L[,j,i])
  Gij=t(Y[i,]-Y[j,])%*%pinv(Sigma_2)%*%(Y[i,]-Y[j,])
  return(Gij)
}

n=5000
adjmat=matrix(rnorm(n^2)*10+15,n,n) #dummy adjmat
adjmat=0.5*(adjmat+t(adjmat))
# set threshold and get adjmat X
n=as.numeric(nrow(adjmat))
threshold=10 #threshold
X=matrix(0,n,n)
X[adjmat>threshold]=1
rs=rowSums(X)
X=X[rs!=0,rs!=0]
n=as.numeric(nrow(X))

K=300 #number of clustering
eig_hat=eigs_sym(X,K,which='LM')
d_hat=eig_hat$values
V_hat=eig_hat$vectors
plot(sort(abs(d_hat),decreasing =TRUE))

#estimation of W_hat
d_diag=diag(d_hat)
W_0_hat=X-V_hat%*%d_diag%*%t(V_hat)
d_tilda=d_hat^3/diag((diag(d_hat^2)+diag(diag(t(V_hat)%*%diag(diag(W_0_hat^2))%*%V_hat))))
d_tilda_diag=diag(d_tilda)
W_hat=X-V_hat%*%d_tilda_diag%*%t(V_hat)
sigma_hat=W_hat

#t_k can be estimated by d_k very well
t=d_hat

Y=V_hat[,2:K]/V_hat[,1]
Y[is.nan(Y)]=1

L=array(0,c(K-1,n,n))
for (a in 1:(K-1)){
  L[a,,]=(t[1]*V_hat[,1]%*%t(V_hat[,(a+1)])-t[a+1]*V_hat[,(a+1)]%*%t(V_hat[,1]))/(t[a+1]*V_hat[,1])
}

G=matrix(0,n,n)

#foreach(i=2:n) %dopar%
# for (i in 2:n)
#   for(j in 1:(i-1))
#   {
#     print(c(i,j))
#     G[i,j]=getG(i,j)
#   }

registerDoParallel(cores=64)
for (i in 2:n){
  print(c(i,i-1))
  tempG=foreach(j=1:(i-1),.combine='cbind') %dopar% {
    Sigma_2=L[,i,]%*%diag(sigma_hat[i,])%*%t(L[,i,])+L[,j,]%*%diag(sigma_hat[j,])%*%t(L[,j,])-sigma_hat[i,j]*L[,i,j]%*%t(L[,j,i])
    t(Y[i,]-Y[j,])%*%pinv(Sigma_2)%*%(Y[i,]-Y[j,])
  }
  G[i,1:(i-1)]=tempG
}
G=G+t(G)
p_value=1-matrix(sapply(G,pchisq,df=K-1),n,n)
