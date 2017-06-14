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
sample.control <- function(data, size, target,  prob=seq(0,1,.05), replace=F)
{
  target.q <- quantile(target, prob)
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
    curr.data = which((data >= data.last) & (data <= (target.q[ii])))
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
