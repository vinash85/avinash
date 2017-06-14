eeQTL.run <- function(obj)
{
  x1 <- obj$x
  x1[is.na(obj$x)] <- 0
  numF <- dim(obj$feature)[2]+1
  numSNP <- dim(obj$x)[2]
  R1.out=eeSNP.regress(x=x1, y=obj$y, feature=obj$feature, alpHa=numeric(numF), beTa=rnorm(numSNP), estimate.alpha=T, estimate.beta=T, B.inv.alpHa=rep(1/100,  numF), itermax = 1000, thin=40, burnIn=20, Debug=T)
  par(mfrow=c(2,2))
  boxplot(R1.out$beTa.mat,main = "eeSNP model with estimated alpha")
  boxplot(R1.out$gamMa.mat,main = "eeSNP model with estimated alpha")
  boxplot(R1.out$gamMa.prob.mat,main = "eeSNP model with estimated alpha")
  #points(obj$beTa,col="red" , pch=18)
  boxplot(R1.out$alpHa.mat)
  #points(obj$alpHa,col="red" , pch=18)
  R1.out
}
epi.eQTL.run <- function(obj)
{
  x1 <- obj$x
  x1[is.na(obj$x)] <- 0
  pairf = ifelse( is.null(obj$pairFeature),0,  ncol(obj$pairFeature[[(which(lapply(obj$mask, length) >0)[1])]]))
  numF <- dim(obj$feature)[2]+ pairf + 1
  numSNP <- dim(obj$x)[2]
  R1.out=epi.eQTL(x=x1, y=obj$y, feature=obj$feature, pairFeature=obj$pairFeature, mask = obj$mask, alpHa=NULL, gamMa = NULL, estimate.alpha=T, estimate.beta=T, B.inv.alpHa=rep(1/100,  numF), itermax = 500, thin=20, burnIn=100, Debug=T)
  R1.out
}
epi.eQTL_run <- function(obj, gamMa=NULL, alpHa=NULL, itermax=1000, thin=1, burnIn=0, cores=2, ...)
{
  x1 <- obj$x
  x1[is.na(obj$x)] <- 0
  pairf = ncol(obj$pairFeature.mat)
  numF <- dim(obj$feature)[2]+ pairf + 1
  numSNP <- dim(obj$x)[2]
  if(is.null(gamMa)) gamMa <- rep(0, nrow(obj$mask2))
  if(is.null(alpHa)) alpHa <- rep(0,numF)
  R1.out=epi_eQTL(x=x1, y=obj$y, feature=obj$feature, pairFeature=obj$pairFeature.mat, mask = obj$mask2, alpHa=alpHa, 
                  gamMa = gamMa,  estimate_alpha=T, estimate_beta=T, B_inv_alpHa=rep(1/10,  numF), itermax = itermax, thin= thin, burnIn= burnIn, threads=cores, ...)
  R1.out
}

eeQTL_run <- function(obj)
{
  x1 <- obj$x
  x1[is.na(obj$x)] <- 0
  numF <- dim(obj$feature)[2]+1
  numSNP <- dim(obj$x)[2]
  R1.out=eeSNP_regress(x=x1, y=obj$y, feature=obj$feature, alpHa=numeric(numF), beTa=rnorm(numSNP), estimate_alpha=T, estimate_beta=T, B_inv_alpHa=rep(1/100,  numF), itermax = 5000, thin=40, burnIn=1000)
  par(mfrow=c(2,1))
  boxplot(R1.out$beTa.mat,main = "eeSNP model with estimated alpha")
  points(obj$beTa,col="red" , pch=18)
  boxplot(R1.out$alpHa.mat)
  points(obj$alpHa,col="red" , pch=18)
  R1.out
}
