logistic <- function(x)
  1/(1 + exp(-x))

bf.single <- function(x, y, gene.map, cc=100, log=T)
{
  x.norm = scale(x, center=T, scale=apply(x, 2, sd))
  y.norm = scale(y, center=T, scale=apply(y, 2, sd))
  xTy = t(x.norm)%*%y.norm
  y2 = apply(y.norm, 2, function(tt) sum(tt^2))
  n <- nrow(x)
  #snp  <- eeSNP.real[303] 
  bf <- numeric(ncol(x))
  #bf1 <- numeric(ncol(x))
  genes <- unique(gene.map)
  for(gene in genes){
    curr.inx <- which(gene.map == gene)
    bf[curr.inx] =  (n - 2.0)/2.0 * log(1.0 + cc) - (n -1.0 )/2.0 * log(1+cc*(1.0- ((xTy[curr.inx, gene])^2/(colSums((x.norm[,curr.inx])^2) * y2[gene])) ));
    #beta  <- xTy[curr.inx, gene]/colSums((x.norm[,curr.inx])^2)
    #x.curr  <-  x.norm[,curr.inx]
   #r2 <- 1- sapply(1:length(beta),  function(tt)  sum((y.norm[,gene] - x.curr[, tt] * beta[tt])^2) )/y2[gene]
   #bf1[curr.inx]  <- (n - 2.0)/2.0 * log(1.0 + cc) - (n -1.0 )/2.0 * log(1+cc*(1-r2))
  }
  if(!log) bf= exp(bf)
  #cbind(bf,bf1)
  bf
}

bf.single.guan <- function(x, y, gene.map, sigmaA=.4, log=T, nthreads=1)
{
  require(doMC)
  require(foreach)
   registerDoMC(cores=nthreads)
  x.norm = scale(x, center=T, scale=apply(x, 2, sd))
  y.norm = scale(y, center=T, scale=apply(y, 2, sd))
  sigmaA2 = sigmaA^2
  cc = .1
  #xTy = t(x.norm)%*%y.norm
  y2 = apply(y.norm, 2, function(tt) sum(tt^2))
  n <- nrow(x)
  #snp  <- eeSNP.real[303] 
  genes <- unique(gene.map)
  bf2 <- foreach(gene = genes, .inorder=T ) %dopar%{
    curr.inx <- which(gene.map == gene)
  xTy = t(x.norm[, curr.inx])%*%y.norm[,gene]
    bf.temp=  (n - 2.0)/2.0 * log(1.0 + cc) - (n -1.0 )/2.0 * log(1+cc*(1.0- ((xTy)^2/(colSums((x.norm[,curr.inx])^2) * y2[gene])) ));
    x.curr  <-  x.norm[,curr.inx]
    omega  <-  1/(1/sigmaA2 + colSums((x.norm[,curr.inx])^2)) 
    beta  <- xTy*omega
   #r2 <- 1- sapply(1:length(beta),  function(tt)  sum((y.norm[,gene] - x.curr[, tt] * beta[tt])^2) )/y2[gene]
   rxy <- 1 - beta * xTy/y2[gene] 
   bf1.temp <- 1/2 * log(n * omega/sigmaA2)  - n/2 * log(rxy)
   cbind(bf.temp, bf1.temp)
  }
  bf2 = do.call( rbind, bf2)
  #bf2 <- cbind(bf,bf1)
  if(!log) bf2= exp(bf2)
  bf2
}
bf.single.split <- function(yxmat,sigmaA=.4, log=T, nthreads=1)
{
  require(doMC)
  require(foreach)
   registerDoMC(cores=nthreads)
  sigmaA2 = sigmaA^2
  cc = .1
  n <- nrow(yxmat[[1]])
  bf2 <- foreach(yx = yxmat, .inorder=T ) %dopar%{
    x.norm <- yx[,-1]
    x.norm = scale(x.norm, center=T, scale=apply(x.norm, 2, sd))
    y.norm <- yx[,1]
    y.norm = scale(y.norm, center=T, scale=sd(y.norm)) 
    xTy = t(x.norm)%*%y.norm
    y2  <-  sum(y.norm^2)
    bf.temp=  (n - 2.0)/2.0 * log(1.0 + cc) - (n -1.0 )/2.0 * log(1+cc*(1.0- ((xTy)^2/(colSums((x.norm)^2) * (y2)) )));
    x.curr  <-  x.norm
    omega  <-  1/(1/sigmaA2 + colSums((x.norm)^2)) 
    beta  <- xTy*omega
   #r2 <- 1- sapply(1:length(beta),  function(tt)  sum((y.norm[,gene] - x.curr[, tt] * beta[tt])^2) )/y2[gene]
   rxy <- 1 - beta * xTy/y2 
   bf1.temp <- 1/2 * log(n * omega/sigmaA2)  - n/2 * log(rxy)
   cbind(bf.temp, bf1.temp)
  }
  bf2 = do.call( rbind, bf2)
  if(!log) bf2= exp(bf2)
  bf2
}
