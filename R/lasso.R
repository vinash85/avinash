lasso.eqtl<- function(x, y, gene.map, nthreads=1)
{
require(glmnet)
require(doMC)
require(foreach)

#registerDoMC(cores=64)
   registerDoMC(cores=nthreads)
  #x.norm = scale(x, center=T, scale=apply(x, 2, sd))
  #y.norm = scale(y, center=T, scale=apply(y, 2, sd))
  n <- nrow(x)
  coefs <- numeric(ncol(x))
  genes <- unique(gene.map)
  coefs <- foreach(gene = genes, .inorder=T ) %dopar%{
    curr.inx <- which(gene.map == gene)
    fit <- cv.glmnet(x=x[,curr.inx], y=y[, gene],  alpha=1,  parallel=T)
    coef = coef( fit, s="lambda.min")
    coef <- coef[-1]
    coef
  }
   coef = do.call( c, coefs)
  coefs
}


