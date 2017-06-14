precision.recall <- function(threshold, lab, pred)
{
    pos <- lab > 0.5; neg <- lab <=0.5
    pred.pos <- pred > threshold 
    pred.neg <-  pred <= threshold
    tp <- sum(pos & pred.pos); fp  <- sum(neg & pred.pos)
    fn <- sum(pos & pred.neg)
     c(precision=tp/(tp +fp), recall=tp/(tp +fn))
}


plot.ev <- function(obj, out1)
{
    gamMa.nz = apply(out1$gamMa.mat,2 , mean)
    numGenes = ncol(obj$y)
    create.mask = function(mask2, gene.range){
	mask <- list()
	for(gene in seq(gene.range)){
	    mask[[gene]] <- mask2[mask2[,1]==gene,2]
	}
	mask
    }
    obj$mask = create.mask(obj$mask2, seq(numGenes))
    regulators  = sapply(seq(numGenes), function(x) sum(gamMa.nz[obj$mask[[x]]] >= 0.5) )
    #regulators  = sapply(seq(numGenes), function(x) floor(sum(gamMa.nz[obj$mask[[x]]])) )
    x1 <- obj$x
    x1[is.na(obj$x)] <- 0
    y.back = obj$y 
    obj$y = t((t(obj$y)- colMeans(obj$y))/apply(obj$y,2,sd))
    numGenes =  ncol(obj$y)
    #x1 = t( (t(x1) - colMeans(x1)) / sd(x1)) 
    ## explained variance analysis
    ev = function(gene){
	mask = obj$mask[[gene]]
	gamMa = apply(out1$gamMa.mat[, mask, drop=F],2, median)
	nz = which(gamMa > 0.5)
	if(length(nz) ==0)
	    return(0) 
	beTa = apply(out1$beTa.mat[, mask] ,2 , mean)
	beTa.nz = matrix( beTa[nz, drop=F], ncol=1)
	y.var = var(obj$y[,gene])
	y.residue = obj$y[, gene] - x1[,mask[nz], drop=F] %*% beTa.nz
	residue.var = var(y.residue)
	out = (1 - residue.var/y.var)
	return(out)
    }

    #ev.100 = sapply(seq(numGenes), function(x) ev(x))
    if(is.null(out1$ev)) out1$ev = 1 - out1$sigma2.mat 
    ev.100 = (apply(out1$ev, 2, mean))
    #gene.analyzed = which(ev.100 > .15)
    #ev.ga = ev.100[gene.analyzed]
    #regulators.ga = c(unlist(regulators[gene.analyzed]))
    beTa.null = function(gene){
	require(MASS)
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	#x = center(x1[,mask])
	x = center(x1[,mask]) 
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	V.w = ginv((t(x)%*%x) + small * diag(p))
	m.w = V.w %*% t(x) %*% y
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out$gamMa = mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }
    #debug(beTa.null)
   null.ga = numeric(numGenes) 
   null.gamMa  <- list()
   for(gene in seq(numGenes)){
     temp  <-  beTa.null(gene)
     null.ga[gene]  <-  temp$ev
     null.gamMa[[gene]] <- temp$gamMa
    }
   null.gamMa <- c(unlist(null.gamMa))
    #null.ga = sapply(seq(numGenes), beTa.null)
    trad.eQTL = function(gene){
	require(MASS)
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	m.w = apply(x,2, function(xx) cor(xx, y))
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out <- list()
	out$gamMa <- mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }

    #eQTL.ev = sapply(seq(numGenes), trad.eQTL)
   eQTL.ev <-  numeric(numGenes)
   eQTL.gamMa  <- list()
   for(gene in seq(numGenes)){
     temp  <-  trad.eQTL(gene)
     eQTL.ev[gene]  <-  temp$ev
     eQTL.gamMa[[gene]] <- temp$gamMa
    }
   eQTL.gamMa <- c(unlist(eQTL.gamMa))
    ev.combined = data.frame(ev = c(null.ga, eQTL.ev,  ev.100)*100, 
			     label = as.factor(c( 
						 rep("multivariate", length(null.ga)), 
						 rep("univariate", length(null.ga)),
						 rep("epigenetic_eQTL",length(ev.100))
						 ) )
			     )
    #ev.combined = data.frame(ev = c(null.ga ,eQTL.ev)*100, 
			     #label = as.factor( c(rep("regression",length(ev.100)),
						  #rep("eQTL", length(null.ga))) ))
    library(ggplot2)
    #qplot(ev, data=ev.combined, geom="density", fill=label, alpha=I(.5),
	  #main="", xlab="Precentage of variance explained",
	  #ylab="Density")+theme(legend.position = c(.8, .8))
    #browser()
    p <- ggplot(ev.combined, aes(x=ev, fill=label)) + xlab("Precentage of variance explained") 
    p <-  p+ geom_density( aes( y=..scaled..), alpha=I(.5)) +theme(legend.position = c(.8, .8))
    
    epi.sd = apply(out1$ev, 2, sd)
    df1 = data.frame( epi=ev.100 , epi.sd=epi.sd,  eQTL = eQTL.ev, regression=null.ga)
    df1.order = df1[ order(ev.100,decreasing=T),]
    df1.order$inx = seq(numGenes)
    h1 <- ggplot(df1.order, aes(x=inx))
    h1 <- h1+ geom_ribbon(aes(ymin=epi - (epi.sd/2), ymax = epi+(epi.sd/2)), alpha=I(.15), fill="blue") +
   	 geom_line(aes( y=epi, color="epi" ))+ geom_line(aes( y=eQTL, color="eQTL")) + 
	geom_line(aes( y=regression, color="regression")) 
	p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())
	p = p+theme(strip.background=element_blank(), strip.text.y = element_text())   
    list(ev.combined=ev.combined, ev.100=ev.100, eQTL.ev=eQTL.ev, null.ga=null.ga, p=p, h1=h1, mask = obj$mask, null.gamMa=null.gamMa, eQTL.gamMa=eQTL.gamMa )
}




plot.EV <- function(obj, out1, threads=1, use_sum=F)
{
  require(eeSNP);
    gamMa.nz = apply(out1$gamMa.mat,2 , mean)
    numGenes = ncol(obj$y)
    create.mask = function(mask2, gene.range){
	mask <- list()
	for(gene in seq(gene.range)){
	    mask[[gene]] <- mask2[mask2[,1]==gene,2]
	}
	mask
    }
    obj$mask = create.mask(obj$mask2, seq(numGenes))
    if(!use_sum) regulators  = sapply(seq(numGenes), function(x) sum(gamMa.nz[obj$mask[[x]]] >= 0.5) )
    else regulators  = sapply(seq(numGenes), function(x) round(sum(gamMa.nz[obj$mask[[x]]])) )
    x1 <- obj$x
    x1[is.na(obj$x)] <- 0
    y.back = obj$y 
    obj$y = t((t(obj$y)- colMeans(obj$y))/apply(obj$y,2,sd))
    numGenes =  ncol(obj$y)
    #x1 = t( (t(x1) - colMeans(x1)) / sd(x1)) 
    ## explained variance analysis
    ev = function(gene){
	mask = obj$mask[[gene]]
	gamMa = apply(out1$gamMa.mat[, mask, drop=F],2, median)
	nz = which(gamMa > 0.5)
	if(length(nz) ==0)
	    return(0) 
	beTa = apply(out1$beTa.mat[, mask] ,2 , mean)
	beTa.nz = matrix( beTa[nz, drop=F], ncol=1)
	y.var = var(obj$y[,gene])
	y.residue = obj$y[, gene] - x1[,mask[nz], drop=F] %*% beTa.nz
	residue.var = var(y.residue)
	out = (1 - residue.var/y.var)
	return(out)
    }
    

    #ev.100 = sapply(seq(numGenes), function(x) ev(x))
    if(is.null(out1$ev)) out1$ev = 1 - out1$sigma2.mat 
    ev.100 = (apply(out1$ev, 2, mean))
    print("ev done")
    #gene.analyzed = which(ev.100 > .15)
    #ev.ga = ev.100[gene.analyzed]
    #regulators.ga = c(unlist(regulators[gene.analyzed]))
    beTa.null = function(gene){
	require(MASS)
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	V.w = ginv((t(x)%*%x) + small * diag(p))
	m.w = V.w %*% t(x) %*% y
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out$gamMa = mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }
    #debug(beTa.null)
   null.ga = numeric(numGenes) 
   null.gamMa  <- list()
   require(eeSNP)
   null.ga <- c(multivariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads))  
   #for(gene in seq(numGenes)){
     #temp  <-  beTa.null(gene)
     #null.ga[gene]  <-  temp$ev
     #null.gamMa[[gene]] <- temp$gamMa
    #}
   null.gamMa <- c(unlist(null.gamMa))
    print("NULL done")
    #null.ga = sapply(seq(numGenes), beTa.null)
    trad.eQTL = function(gene){
      require(parallel)
      require(multicore)
	require(MASS)
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	m.w = unlist(mclapply(seq(numGenes), function(xx) cor(x[,xx], y)))
	#browser()
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out <- list()
	out$gamMa <- mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }

    #eQTL.ev = sapply(seq(numGenes), trad.eQTL)
   eQTL.ev <-  numeric(numGenes)
   eQTL.gamMa  <- list()
   eQTL.ev <- c(univariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads))  
   #for(gene in seq(numGenes)){
     #temp  <-  trad.eQTL(gene)
     #eQTL.ev[gene]  <-  temp$ev
     #eQTL.gamMa[[gene]] <- temp$gamMa
    #}
   #eQTL.gamMa <- c(unlist(eQTL.gamMa))
    ev.combined = data.frame(ev = c(
				    null.ga, 
				    eQTL.ev,  ev.100)*100, 
			     label = as.factor(c( 
						 rep("multivariate", length(null.ga)), 
						 rep("univariate", length(null.ga)),
						 rep("epigenetic_eQTL",length(ev.100))
						 ) ))
    #ev.combined = data.frame(ev = c(null.ga ,eQTL.ev)*100, 
			     #label = as.factor( c(rep("regression",length(ev.100)),
						  #rep("eQTL", length(null.ga))) ))
    library(ggplot2)
    #qplot(ev, data=ev.combined, geom="density", fill=label, alpha=I(.5),
	  #main="", xlab="Precentage of variance explained",
	  #ylab="Density")+theme(legend.position = c(.8, .8))
    #browser()
    p <- ggplot(ev.combined, aes(x=ev, fill=label)) + xlab("Precentage of variance explained") 
    p <-  p+ geom_density( aes( y=..scaled..), alpha=I(.5)) +theme(legend.position = c(.8, .8))
    
    epi.sd = apply(out1$ev, 2, sd)
    df1 = data.frame( epi=ev.100 , epi.sd=epi.sd,  eQTL = eQTL.ev)
    df1.order = df1[ order(ev.100,decreasing=T),]
    df1.order$inx = seq(numGenes)
    h1 <- ggplot(df1.order, aes(x=inx))
    h1 <- h1+ geom_ribbon(aes(ymin=epi - (epi.sd/2), ymax = epi+(epi.sd/2)), alpha=I(.15), fill="blue") +
   	 geom_line(aes( y=epi, color="epi" ))+ geom_line(aes( y=eQTL, color="eQTL"))
	 #+ 
	#geom_line(aes( y=regression, color="regression")) 
	p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())
	p = p+theme(strip.background=element_blank(), strip.text.y = element_text())   
    
    list(ev.combined=ev.combined, ev.100=ev.100, eQTL.ev=eQTL.ev, null.ga=null.ga, p=p, h1=h1, mask = obj$mask, null.gamMa=null.gamMa, eQTL.gamMa=eQTL.gamMa )
}

plot.gamma.EV <- function(obj, out1, threads=1, gamMa.nz=NULL, threshold=0.5, multivariate=T)
{
  require(eeSNP); require(MASS)
      require(parallel)
      require(multicore)

    if(is.null(gamMa.nz)) gamMa.nz = apply(out1$gamMa_rate.mat,2 , median)
    numGenes = ncol(obj$y)
    create.mask = function(mask2, gene.range){
	mask <- list()
	for(gene in seq(gene.range)){
	    mask[[gene]] <- mask2[mask2[,1]==gene,2]
	}
	mask
    }
    obj$mask = create.mask(obj$mask2, seq(numGenes))
    regulators  = sapply(seq(numGenes), function(x) sum(gamMa.nz[obj$mask[[x]]] >= threshold) )
    x1 <- obj$x
    x1[is.na(obj$x)] <- 0
    x1 <- scale(x1, center=T, scale=apply(x1,2,sd))
    obj$y = scale( obj$y, center=T, apply(obj$y,2,sd))
    numGenes =  ncol(obj$y)
    
    ## explained variance analysis
    ev = function(gene){
	mask = obj$mask[[gene]]
	gamMa = gamMa.nz[mask] 
	nz = which(gamMa > threshold)
	if(length(nz) ==0)
	    return(0) 
	x = x1[,mask[nz]]
	y = obj$y[,gene] 
	p = length(nz)
	V.w = ginv((t(x)%*%x))
	m.w = V.w %*% t(x) %*% y
	y.var = var(obj$y[,gene])
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out = (1 - residue.var/y.var)
	return(out)
    }
    

    #ev.100 = sapply(seq(numGenes), function(x) ev(x))
    ev.100 = unlist(mclapply(seq(numGenes), function(xx) ev(xx)))
    #if(is.null(out1$ev)) out1$ev = 1 - out1$sigma2.mat 
    #ev.100 = 
    print("ev done")
    #gene.analyzed = which(ev.100 > .15)
    #ev.ga = ev.100[gene.analyzed]
    #regulators.ga = c(unlist(regulators[gene.analyzed]))
    beTa.null = function(gene){
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	V.w = ginv((t(x)%*%x) + small * diag(p))
	m.w = V.w %*% t(x) %*% y
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out$gamMa = mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }
    #debug(beTa.null)
   null.ga = numeric(numGenes) 
   null.gamMa  <- list()
   require(eeSNP)
   if(multivariate){ 
     null.ga <- c(multivariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads)) 
   }else{
     null.ga <- NULL
   }
   #for(gene in seq(numGenes)){
     #temp  <-  beTa.null(gene)
     #null.ga[gene]  <-  temp$ev
     #null.gamMa[[gene]] <- temp$gamMa
    #}
   null.gamMa <- c(unlist(null.gamMa))
    print("NULL done")
    #null.ga = sapply(seq(numGenes), beTa.null)
    trad.eQTL = function(gene){
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	m.w = unlist(mclapply(seq(numGenes), function(xx) cor(x[,xx], y)))
	#browser()
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out <- list()
	out$gamMa <- mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }

    #eQTL.ev = sapply(seq(numGenes), trad.eQTL)
   eQTL.ev <-  numeric(numGenes)
   eQTL.gamMa  <- list()
   eQTL.ev <- c(univariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads))  
   #for(gene in seq(numGenes)){
     #temp  <-  trad.eQTL(gene)
     #eQTL.ev[gene]  <-  temp$ev
     #eQTL.gamMa[[gene]] <- temp$gamMa
    #}
   #eQTL.gamMa <- c(unlist(eQTL.gamMa))
    ev.combined = data.frame(ev = c(
				    null.ga, 
				    eQTL.ev,  ev.100)*100, 
			     label = as.factor(c( 
						 rep("multivariate", length(null.ga)), 
						 rep("univariate", length(null.ga)),
						 rep("epigenetic_eQTL",length(ev.100))
						 ) ))
    #ev.combined = data.frame(ev = c(null.ga ,eQTL.ev)*100, 
			     #label = as.factor( c(rep("regression",length(ev.100)),
						  #rep("eQTL", length(null.ga))) ))
    library(ggplot2)
    #qplot(ev, data=ev.combined, geom="density", fill=label, alpha=I(.5),
	  #main="", xlab="Precentage of variance explained",
	  #ylab="Density")+theme(legend.position = c(.8, .8))
    #browser()
    p <- ggplot(ev.combined, aes(x=ev, fill=label)) + xlab("Precentage of variance explained") 
    p <-  p+ geom_density( aes( y=..scaled..), alpha=I(.5)) +theme(legend.position = c(.8, .8))
    
    #epi.sd = apply(out1$ev, 2, sd)
    #df1 = data.frame( epi=ev.100 , epi.sd=epi.sd,  eQTL = eQTL.ev)
    #df1.order = df1[ order(ev.100,decreasing=T),]
    #df1.order$inx = seq(numGenes)
    #h1 <- ggplot(df1.order, aes(x=inx))
    #h1 <- h1+ geom_ribbon(aes(ymin=epi - (epi.sd/2), ymax = epi+(epi.sd/2)), alpha=I(.15), fill="blue") +
   	 #geom_line(aes( y=epi, color="epi" ))+ geom_line(aes( y=eQTL, color="eQTL"))
	 #+ 
	#geom_line(aes( y=regression, color="regression")) 
	p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())
	p = p+theme(strip.background=element_blank(), strip.text.y = element_text())   
    
    list(ev.combined=ev.combined, ev.100=ev.100, eQTL.ev=eQTL.ev, null.ga=null.ga, p=p, mask = obj$mask, null.gamMa=null.gamMa, eQTL.gamMa=eQTL.gamMa )
}




plot.snp.EV <- function(obj, out1, threads=1, gamMa.nz=NULL, threshold=3)
{
  require(eeSNP); require(MASS)

    if(is.null(gamMa.nz)) gamMa.nz = apply(out1$gamMa_rate.mat,2 , median)
    numGenes = ncol(obj$y)
    create.mask = function(mask2, gene.range){
	mask <- list()
	for(gene in seq(gene.range)){
	    mask[[gene]] <- mask2[mask2[,1]==gene,2]
	}
	mask
    }
    numGenes =  ncol(obj$y)
    obj$mask = create.mask(obj$mask2, seq(numGenes))
    #regulators  = sapply(seq(numGenes), function(x) sum(gamMa.nz[obj$mask[[x]]] >= threshold) )
    regulators = rep(threshold, ncol(obj$x)) 
    x1 <- obj$x
    x1[is.na(obj$x)] <- 0
    x1 <- scale(x1, center=T, scale=apply(x1,2,sd))
    obj$y = scale( obj$y, center=T, apply(obj$y,2,sd))
    
    ## explained variance analysis
    ev = function(gene){
	mask = obj$mask[[gene]]
	gamMa = gamMa.nz[mask]
	nz = order(gamMa,decreasing=T)[1:threshold]
	#nz = which(gamMa > threshold)
	x = x1[,mask[nz]]
	y = obj$y[,gene] 
	p = length(nz)
	V.w = ginv((t(x)%*%x))
	m.w = V.w %*% t(x) %*% y
	y.var = var(obj$y[,gene])
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out = (1 - residue.var/y.var)
	return(out)
    }
    

    ev.100 = sapply(seq(numGenes), function(x) ev(x))
    #if(is.null(out1$ev)) out1$ev = 1 - out1$sigma2.mat 
    #ev.100 = 
    print("ev done")
    #gene.analyzed = which(ev.100 > .15)
    #ev.ga = ev.100[gene.analyzed]
    #regulators.ga = c(unlist(regulators[gene.analyzed]))
    beTa.null = function(gene){
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	V.w = ginv((t(x)%*%x) + small * diag(p))
	m.w = V.w %*% t(x) %*% y
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out$gamMa = mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }
    #debug(beTa.null)
   null.ga = numeric(numGenes) 
   null.gamMa  <- list()
   require(eeSNP)
   null.ga <- c(multivariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads))  
   #for(gene in seq(numGenes)){
     #temp  <-  beTa.null(gene)
     #null.ga[gene]  <-  temp$ev
     #null.gamMa[[gene]] <- temp$gamMa
    #}
   null.gamMa <- c(unlist(null.gamMa))
    print("NULL done")
    #null.ga = sapply(seq(numGenes), beTa.null)
    trad.eQTL = function(gene){
      require(parallel)
      require(multicore)
	small = .001
	nz.num = regulators[gene]
	out = list(gamMa=NULL, ev=0)
	if(nz.num == 0) return(out)
	mask = obj$mask[[gene]]
	x = center(x1[,mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	p = length(mask)
	m.w = unlist(mclapply(seq(numGenes), function(xx) cor(x[,xx], y)))
	#browser()
	inx = order(abs(m.w), decreasing=T)
	nz = inx[1:nz.num];
	out <- list()
	out$gamMa <- mask[nz]
	x = x[,nz, drop=F]
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out$ev = (1 - residue.var/y.var)
	return(out)
    }

    #eQTL.ev = sapply(seq(numGenes), trad.eQTL)
   eQTL.ev <-  numeric(numGenes)
   eQTL.gamMa  <- list()
   eQTL.ev <- c(univariate_eQTL(obj$mask2, numRegulators=regulators, x=x1, y =obj$y, threads=threads))  
   #for(gene in seq(numGenes)){
     #temp  <-  trad.eQTL(gene)
     #eQTL.ev[gene]  <-  temp$ev
     #eQTL.gamMa[[gene]] <- temp$gamMa
    #}
   #eQTL.gamMa <- c(unlist(eQTL.gamMa))
    ev.combined = data.frame(ev = c(
				    null.ga, 
				    eQTL.ev,  ev.100)*100, 
			     label = as.factor(c( 
						 rep("multivariate", length(null.ga)), 
						 rep("univariate", length(null.ga)),
						 rep("epigenetic_eQTL",length(ev.100))
						 ) ))
    #ev.combined = data.frame(ev = c(null.ga ,eQTL.ev)*100, 
			     #label = as.factor( c(rep("regression",length(ev.100)),
						  #rep("eQTL", length(null.ga))) ))
    library(ggplot2)
    #qplot(ev, data=ev.combined, geom="density", fill=label, alpha=I(.5),
	  #main="", xlab="Precentage of variance explained",
	  #ylab="Density")+theme(legend.position = c(.8, .8))
    #browser()
    p <- ggplot(ev.combined, aes(x=ev, fill=label)) + xlab("Precentage of variance explained") 
    p <-  p+ geom_density( aes( y=..scaled..), alpha=I(.5)) +theme(legend.position = c(.8, .8))
    
    epi.sd = apply(out1$ev, 2, sd)
    df1 = data.frame( epi=ev.100 , epi.sd=epi.sd,  eQTL = eQTL.ev)
    df1.order = df1[ order(ev.100,decreasing=T),]
    df1.order$inx = seq(numGenes)
    h1 <- ggplot(df1.order, aes(x=inx))
    h1 <- h1+ geom_ribbon(aes(ymin=epi - (epi.sd/2), ymax = epi+(epi.sd/2)), alpha=I(.15), fill="blue") +
   	 geom_line(aes( y=epi, color="epi" ))+ geom_line(aes( y=eQTL, color="eQTL"))
	 #+ 
	#geom_line(aes( y=regression, color="regression")) 
	p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())
	p = p+theme(strip.background=element_blank(), strip.text.y = element_text())   
    
    list(ev.combined=ev.combined, ev.100=ev.100, eQTL.ev=eQTL.ev, null.ga=null.ga, p=p, h1=h1, mask = obj$mask, null.gamMa=null.gamMa, eQTL.gamMa=eQTL.gamMa )
}


# matched ev
plot.matched.EV <- function(obj,  gamMa.nz, gamMa.target, hreads=1, threshold=0.5)
{
  require(eeSNP); require(MASS)

    numGenes = ncol(obj$y)
    create.mask = function(mask2, gene.range){
	mask <- list()
	for(gene in seq(gene.range)){
	    mask[[gene]] <- mask2[mask2[,1]==gene,2]
	}
	mask
    }
    obj$mask = create.mask(obj$mask2, seq(numGenes))
    x1 <- obj$x
    x1[is.na(obj$x)] <- 0
    x1 <- scale(x1, center=T, scale=apply(x1,2,sd))
    obj$y = scale( obj$y, center=T, apply(obj$y,2,sd))
    numGenes =  ncol(obj$y)
    small <- .001
    
    ## explained variance analysis
    ev = function(gene){
	mask = obj$mask[[gene]]
	#gamMa = gamMa.nz[mask]
	nz.num  <- sum(gamMa.nz[mask] > threshold)
	out <- 0
	if(nz.num == 0) return(out)
	inx <- order(gamMa.target[mask] , decreasing = T)
	nz = inx[1:nz.num];
	x = x1[,mask[nz], drop=F]
	y = obj$y[,gene] - mean( obj$y[,gene])
	V.w = ginv((t(x)%*%x) + small * diag(length(nz)))
	m.w = V.w %*% t(x) %*% y
	y.var = var(y)
	y.residue = y - (x %*% m.w)
	residue.var = var(y.residue)
	out = (1 - residue.var/y.var)
	return(out)
    }
      require(parallel)
      require(multicore)
    
      #browser()
    ev = unlist( mclapply(seq(numGenes), function(x) ev(x)))
    #ev = sapply(seq(numGenes), function(x) ev(x))
	#m.w = unlist(mclapply(seq(numGenes), function(xx) cor(x[,xx], y)))
    return(ev)
}

