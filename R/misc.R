.ls.objects <- function (pos = 1, pattern, order.by,
			 decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
				       fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
			   capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
		      as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
confusion.matrix  <- function(thres, gamMa, gamMa.prob)
{
  t.pos <- which(gamMa > 0.5)
  t.neg <-  which( gamMa <=0.5)
  pos  <- which(gamMa.prob > thres)
  neg  <- which( ! (seq(length(gamMa.prob) )%in% pos))
  out  <-  matrix(0,2,2)
  dimnames(out) = list(list(0,1), list(0,1))
  out[1,1] <- length(which(t.neg %in%  neg))
  out[2,1] <- length(which(t.neg %in%  pos))
  out[1,2] <-length(which( t.pos %in%  neg))
  out[2,2] <-length(which( t.pos %in%  pos))
  out
}
find.corr  <-  function(inx) {
  temp = cor.test(obj$y[,obj$mask2[inx,1]], obj$x[,inx])
  return(c(temp$p.value, temp$estimate))
}


multivariate.eQTL.inx = function( obj){
    require(MASS)
    small = .001
    multi <- numeric(ncol(obj$x))
    for(gene in seq(ncol(obj$y))){
	mask = obj$mask2[ obj$mask[,1] == gene,2]
	x  <- center(obj$x1[, mask])
	y = obj$y[,gene] - mean( obj$y[,gene])
	V.w = ginv((t(x)%*%x) + small * diag(length(mask)))
	multi[mask] = V.w %*% t(x) %*% y
    }
    multi
}
univariate.eQTL.inx = function( obj1, cores=1){
  require(MASS)
  library(doMC)
  registerDoMC(cores=cores)
  small = .001
  eQTL <- numeric(ncol(obj1$x))
  eQTL <- foreach(gene = seq(ncol(obj1$y)),  .inorder=T)%dopar%{
    #eQTL <- for(gene  in seq(ncol(obj1$y))){
      mask <- which(obj$mask2[1,] == gene)
    y = obj1$y[,gene]
    apply(obj1$x[, mask, drop=F], 2, function(xx) cor(xx, y))
  }
  eQTL
}
univariate.eQTL.eeSNP = function( xx, yy, eeSNP, mask2, cores=1){
  require(MASS)
  library(doMC)
  registerDoMC(cores=cores)
  small = .001
  eQTL <- numeric(ncol(xx))
  mask.gene <- mask2[,1]
  eQTL <- foreach(snp = seq(ncol(xx)),  .inorder=T)%dopar%{
     gene <- mask.gene[eeSNP[snp]] 
    y = yy[,gene]
    aa= cor.test(xx[,snp], y)
    c(aa$p.value, aa$estimate)
  }
  do.call( rbind, eQTL)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

featureExtract <- function(input, dir)
{

  library(data.table)

  files <- list.files(path = dir)
  feature <- data.table(input) 
  for(fid in files){
    name=sub("(*).500", "\\1", fid)
    name <- gsub(name, pattern="-" , replacement=".")
    name <- gsub(name, pattern="\\+" , replacement="p")
    print(name)
    a <- fread(paste(dir, fid, sep="/"))
    setkey(a, V5)
    nCol <- ncol(a) 
    if(nrow(a)!=nrow(feature)) stop(paste("number of rows not equal in: ", fid))
    expr <- parse(text = paste0(name, paste(':=a[,V', nCol, ']',sep="") ))
    #expr <- parse(text = paste0(name, ':=a[,V9]' ))
    #print(expr)
    feature[, eval(expr)]
  }
  feature
}
set.display <- function(){
  a =readLines("~/.last-display")
  Sys.setenv("DISPLAY"=a)
}
elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

sample.control <- function(data, size, target,  prob=seq(0,1,.05), replace=F)
{
  target.q <- quantile(target, prob)
  big <- 1e8
  small <- (max(target) - min(target))/(length(prob) *big)
  # matiching zeros
  inx <- numeric(size)
  last.ii = 0
  iter  <-  length(prob)
  prob <- c(0, prob)
  prob <- prob/max(prob)
  data.last  <-  min(data) -.01
  for(ii in seq(iter)){
    size.curr = floor(size*prob[ii+1]) - last.ii
    #if(ii == iter) sizer.curr = size - last.ii
    curr.data = which((data > (data.last-small) ) & (data < (target.q[ii] + small)))
    if(size.curr > 0) inx[(last.ii + 1) : (last.ii + size.curr)] = sample(curr.data, size=size.curr, replace=replace)
    last.ii = last.ii + size.curr
    data.last = target.q[ii]
  }
  inx
}

sample.distribution <- function(target, size, prob=seq(0,1,.05))
{
  # sampling arbitary continuous/discrete distribution and output Continuous data (inx is continuous) 
  target.q <- quantile(target, prob)
  # matiching zeros
  inx <- numeric(size)
  last.ii = min(target)
  iter  <-  length(prob)
  prob <- c(0, prob)
  prob <- prob/max(prob)
  data.last  <-  min(target)
  for(ii in seq(iter)){
    size.curr = floor(size*prob[ii+1]) - last.ii
    #if(ii == iter) sizer.curr = size - last.ii
    if(size.curr > 0) inx[(last.ii + 1) : (last.ii + size.curr)]  <- data.last + runif(size.curr) * (target.q[ii] -data.last) 

    #curr.data = which((data >= data.last) & (data <= (target.q[ii])))
    #if(size.curr > 0) inx[(last.ii + 1) : (last.ii + size.curr)] = sample(curr.data, size=size.curr, replace=replace)
    last.ii = last.ii + size.curr
    data.last = target.q[ii]
  }
  inx
}

ensembl2hgnc = function(){
# library(biomaRt)
ensembl <- biomaRt::useMart("ensembl", dataset= "hsapiens_gene_ensembl")
dt1 = biomaRt::getBM(attributes = c( "hgnc_symbol", "ensembl_gene_id"), mart=ensembl)
dt1= data.table(dt1)
setkey(dt1,ensembl_gene_id)
dt1
}
# 
# enterz2hgnc = function(){
# library(org.Hs.eg.db)
# x <- org.Hs.egREFSEQ2EG
# mapped_genes <- mappedkeys(x)
# xx <- as.list(x[mapped_genes])
# xx 
# }
detachAllPackages <- function() {
    basic.packages.blank <- c(    
        "stats",    
        "graphics",    
        "grDevices",    
        "utils",   
        "datasets",  
        "methods",    
        "base"    
    )    
    basic.packages <- paste("package:", basic.packages.blank, sep = "")   
    package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]   
    package.list <- setdiff(package.list, basic.packages)   
    if (length(package.list) > 0) {   
        for (package in package.list) {   
            detach(package, character.only = TRUE)   
        }   
    }    
}


get.synonym = function(gene){
  geneSynonym::humanSyno(gene,  caseSensitive = F) %>% 
  unlist %>% toupper 
}

convert.cd2hugo = function(genes, cd2hugo){
  out.all = list()
  for (gene in genes) {
    out  = NULL
    for (tt in cd2hugo) {
      if(length(out)== 0){
        # strsplit(tt,split = ",") %>% 
        out = grep(pattern=paste("^",gene, "$", sep=""), tt, ignore.case = T)
      }
    }
    hugo.id =cd2hugo$`Approved symbol`[out]
    hugo.id2 = geneSynonym::humanSyno(gene,  caseSensitive = F) %>% 
      unlist %>% toupper %>%
      intersect(cd2hugo$`Approved symbol`)
    out.all[[gene]] = c(hugo.id, hugo.id2) %>% unique
  }
  out.all
  
}