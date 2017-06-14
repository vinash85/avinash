
writeSVMFeature <- function(outfile, pos=NULL, neg=NULL)
{
  outCon = file(outfile, 'w')
  dimLen = dim(pos)[2]
  if(!is.null(pos) && (dimLen > 0) ){
    for(i in seq(1:dim(pos)[1])){
      pos.text = paste("+1",  paste(1:dimLen, pos[i,], sep=":", collapse=" "))
      writeLines(pos.text, con = outCon)
    }

  }
  dimLen = dim(neg)[2]
  if(!is.null(neg) && (dimLen > 0)){
    for(i in seq(1:dim(neg)[1])){
      neg.text = paste("-1",  paste(1:dimLen, neg[i,], sep=":", collapse=" "))
      writeLines(neg.text, con = outCon)
    }
  }
  close(outCon)
}
extractFeatureGff <- function(regionFid, featureMask=NULL, useClust=F)
{

  numMotif = 981
  region = read.table(file=regionFid, header = F, sep="\t", stringsAsFactors =F,  strip.white =T)
  numRegions = dim(region)[1]/numMotif
  ######binaryization#############
  regionMotif =  rep(0, numRegions * numMotif) 
  regionMotif[region$V6 > 0.95 ] = 1 
  stopifnot(length(regionMotif)  ==  numRegions * numMotif)
  dim(regionMotif) = c(numRegions, numMotif)
  if(useClust){
    regionClust = t(clust2matMap %*% t(regionMotif))  
    regionClust[regionClust > 0] = 1
    regionMotif = regionClust
  }
  if(!is.null(featureMask)) regionMotif = regionMotif[,featureMask]
  regionMotif
}
# bar plot of primal variable
findPrimal <- function(modelFile, len=8)
{
  aa = read.table(file=modelFile, header = F, sep=" ", stringsAsFactors =F,  strip.white =T, skip=10)
  svCoef = aa[,1, drop=F]
  svCoef <- apply(svCoef,2,as.numeric)
  #mode(svCoef) = numeric
  mat = aa[,2:len] 
  #SVs = outer(seq(nrow(mat)), seq(ncol(mat)) , FUN=function(r,c) as.numeric(strsplit(mat[r,c],split=":")[[1]][2]))
  SVs = apply(mat, c(1,2), FUN=function(x) as.numeric(unlist(strsplit(x,split=":"))[2]))
  w = t(SVs) %*% svCoef
  return(w)

}
readSVMModel <- function(model)
{
  model.tab = read.table(file=model, header = F, sep=" ", skip=10, stringsAsFactors =F,  strip.white =T)
  out$label = model.tab[,1] 
  out$SV = apply(model.tab[,-1], c(1,2), FUN= function(x){
		 unlist(strsplit(x, split=":"))[2]
			 })
  out
}
