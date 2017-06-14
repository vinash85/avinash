eqtnminer.mle <- function(alpHa, feature, gamMa, mask, ... ){
  ll <- 0
  pi0  <- alpHa[1]
  alpHa1 <- alpHa[-1]
 pi1 <- exp(feature %*% alpHa1)
  for(maskCurr in mask){
    #featureCurr <- feature[maskCurr,]
    #piCurr <- logisitic(featureCurr %*% alpHa) 
    #gamMaCurr <- 
    piCurr <- pi1[maskCurr]
    piCurr <- piCurr/sum(piCurr)
    ll  <- ll + log(pi0 + (1 - pi0) * sum(gamMa[maskCurr] *piCurr) )
  }
  return(-ll)
}

eqtnminer.MLE <- function (alpHa, feature, bf, map, ...) 
{
    ll <- 0
  #pi0 have to be in (0,1)
    pi0 <- logistic(alpHa[1])
    alpHa1 <- alpHa[-1]
    map.unique <- unique(map)
    pi1 <- exp(feature %*% alpHa1)
    for (map.curr in map.unique) {
      maskCurr <- which(map == map.curr)
        piCurr <- pi1[maskCurr]
        piCurr <- piCurr/sum(piCurr)
        ll <- ll + log(pi0 + (1 - pi0) * sum(bf[maskCurr] * 
            piCurr))
    }
    return(-ll)
}

eqtnminer.MLE.ls <- function (alpHa,  bf.feature.list, nthreads=1,  ...) 
{
    ll <- 0
    pi0 <- logistic(alpHa[1])
    alpHa1 <- alpHa[-1]
    ll <- foreach (bf.feature= bf.feature.list, .combine="+")%dopar% {
        piCurr <- exp(bf.feature[,-1] %*% alpHa1)
        piCurr <- piCurr/sum(piCurr)
        log(pi0 + (1 - pi0) * sum(bf.feature[,1] * 
            piCurr))
    }
    return(-ll)
}

cal.eqtnProb <- function (bf, alpHa, map, feature) 
{
    ll <- 0
  #pi0 have to be in (0,1)
    pi0 <- logistic(alpHa[1])
    alpHa1 <- alpHa[-1]
    map.unique <- unique(map)
    pi1 <- exp(feature %*% alpHa1)
    post.p <- numeric(length(bf))
    for (map.curr in map.unique) {
      maskCurr <- which(map == map.curr)
        piCurr <- pi1[maskCurr]
        piCurr <- piCurr/sum(piCurr)
	num <- (1 - pi0) * bf[maskCurr] * piCurr
	deno <- pi0 + sum(num )
	post.p[maskCurr] <- num/deno
    }
    return(post.p)
}



eqtnminer <- function(obj.file, num.eeSNPg=1)
{

  library(eeSNP)
  #library(avinash)
  if(is.list(obj.file)) 
    obj <- obj.file
  else load(obj.file)

  if(is.null(obj$mask)) mask = create.mask(obj$mask2,ncol(obj$y))
  else mask = obj$mask
  #obj.fliping <- function(obj, mask)
  #{
    #for(mm in mask){
      #obj$x[,mm]= obj$x[,sort(mm,decreasing=T)]
      #obj$alpHa[mm]= obj$alpHa[sort(mm,decreasing=T)]
      #obj$beTa[mm]= obj$beTa[sort(mm,decreasing=T)]
      #obj$gamMa[mm]= obj$gamMa[sort(mm,decreasing=T)]
      #obj$feature[mm,]= obj$feature[sort(mm,decreasing=T),]
      ##obj$pairFeature.mat[mm,]= obj$pairFeature.mat[sort(mm,decreasing=T),]
    #}
    #obj
  #}
  #flip.back <- function(in1, mask){
    #for(mm in mask){
      #in1[mm] <- in1[seq(max(mm),min(mm))]
    #}
    #in1
  #}
  x1 <- obj$x
  x1[is.na(obj$x)] <- 0
  obj$dhs <- NULL
  obj$x <- NULL

  pairf = ncol(obj$pairFeature.mat)
  if(is.null(pairf)) pairf=0;
  numF <- dim(obj$feature)[2]+ pairf + 1
  gamMa <- as.matrix(rep(0, nrow(obj$mask2)) )
  alpHa <- as.matrix(rep(0,numF))
  temp <- rep(1: ceiling(nrow(obj$mask2)/2), each=2)
  #mask3 = cbind( obj$mask2, temp)  
numGenes <- ncol(obj$y)
numSNP <- ncol(obj$x)
prior <- num.eeSNPg * numGenes/numSNP
  regulator_prior = 1/100
  B_inv_alpHa = rep(1/1000, numF)
  B_inv_alpHa[1] = 1/10000
  #prior=.2
  alpHa.new <- alpHa
  feature = cbind(1, obj$feature, obj$pairFeature.mat)

  out <- list()
  maxiter <- 15
  out$alpHa.mat <- matrix(0, nrow=maxiter, ncol=length(alpHa))
  out$gamMa.mat <- matrix(0, nrow=maxiter, ncol=ncol(x1))
  out$gamMa_rate.mat <- matrix(0, nrow=maxiter, ncol=ncol(x1))
  for(iter in 1:maxiter){

    out.bf <- epi_eQTL_dissecting(x=x1, y=obj$y, feature=obj$feature, pairFeature=obj$pairFeature.mat, mask = obj$mask2, alpHa=alpHa.new, ratio=1, 
				  gamMa = gamMa,  estimate_alpha=F, estimate_beta=T, B_inv_alpHa=B_inv_alpHa, itermax =1, thin= 1, burnIn=0, threads=64, gamMa_thres = 0,
				  beTa_thres=0.0, balance=T, use_raoblackwell =F, logistic_variable_selection=F, num_logit_train=5000, regulator_prior=regulator_prior, accIter=1, prior=prior)

    save.image("temp.RData")
    gamMa_rate <- apply(out.bf$gamMa_rate.mat, 2, mean)
    alpHaPI <- c(0, alpHa.new)
    mle = optim(par=alpHaPI, fn=eqtnminer.mle, gamMa=gamMa_rate, mask=mask, feature=feature)
    alpHa.new  <- mle$par[-1]
    out$alpHa.mat[iter,] <- alpHa.new
    out$gamMa.mat[iter,]  <- apply(out.bf$gamMa.mat, 2, mean)
    out$gamMa_rate.mat[iter,]  <- gamMa_rate
  }


  #save(file="eqtnminer1.RData", out)
out
}
