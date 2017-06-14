tfFisher <- function(eeSNP.file, control.file, eeSNPnum, controlnum, threshold=-8.5)
{
  require(data.table)
  require(fdrtool)
  tf = fread(eeSNP.file)
  tf=tf[V6 < threshold]
  #tf1 = tf[, list(value=max(V6)), by=list(V1,V2,V10)]
  tf.control = fread(control.file)
  tf.control=tf.control[V6 < threshold]
  #tf1.control = tf.control[, list(value=V6, V1,V2,V10)]
  tf.grp =  tf[, list( eeSNP=length(V6)) , by=V2]
  tf.control.grp =  tf.control[, list( motif=(strsplit(V10, split=" ")[[1]][4]),control=length(V6)), by=V2]
  tf.merge = merge(tf.grp, tf.control.grp, by="V2", all=T)
  tf.merge[is.na(tf.merge)] <- 0
  tf.merge$eeSNPneg  <- eeSNPnum - tf.merge$eeSNP
  tf.merge$controlneg  <- controlnum - tf.merge$control
  eeSNP.fisher =tf.merge[, list( motif=motif, eeSNP=eeSNP, eeSNPneg=eeSNPneg, conrol = control, controlneg = controlneg, 
				greater = fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="greater")$p.value, 
				less= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="less")$p.value,
				two.sided= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="two.sided")$p.value),
				by=V2]
  eeSNP.fisher$qval = fdrtool(eeSNP.fisher$two.sided, statistic="pvalue", plot=F)$qval
  setkey(eeSNP.fisher, greater)
  eeSNP.fisher
}

tfFisher.grp <- function(eeSNP.file, control.file, eeSNPnum, controlnum)
{
  require(data.table)
  require(fdrtool)
  tf = fread(eeSNP.file)
  tf1 = tf[, list(value=max(V6)), by=list(V1,V2,V10)]
  tf.control = fread(control.file)
  tf1.control = tf.control[, list(value=max(V6)), by=list(V1,V2,V10)]
  tf.grp =  tf1[, list( eeSNP=length(value)) , by=V2]
  tf.control.grp =  tf1.control[, list( motif=(strsplit(V10, split=" ")[[1]][4]),control=length(value)), by=V2]
  tf.merge = merge(tf.grp, tf.control.grp, by="V2", all=T)
  tf.merge[is.na(tf.merge)] <- 0
  tf.merge$eeSNPneg  <- eeSNPnum - tf.merge$eeSNP
  tf.merge$controlneg  <- controlnum - tf.merge$control
  eeSNP.fisher =tf.merge[, list( motif=motif, eeSNP=eeSNP, eeSNPneg=eeSNPneg, conrol = control, controlneg = controlneg, 
				greater = fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="greater")$p.value, 
				less= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="less")$p.value,
				two.sided= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="two.sided")$p.value),
				by=V2]
  eeSNP.fisher$qval = fdrtool(eeSNP.fisher$two.sided, statistic="pvalue", plot=F)$qval
  setkey(eeSNP.fisher, greater)
  eeSNP.fisher
}
tfFisher.majorminor <- function(eeSNP, control, eeSNPnum, controlnum, threshold=-8.5)
{
  eeSNP.major=paste(eeSNP, ".major.gff", sep="")
  eeSNP.minor=paste(eeSNP, ".minor.gff", sep="")
  control.major=paste(control, ".major.gff", sep="")
  control.minor=paste(control, ".minor.gff", sep="")

  require(data.table)
  require(fdrtool)
  tf = fread(eeSNP.major)
  setnames(tf, "V6", "major")
  tf$minor <- fread(eeSNP.minor)$V6
  tf[,change:=abs(major - minor)]
  tf[,sig:=xor((major<threshold), (minor<threshold))]
  tf[,disrupts:= ifelse(sig & (change > 1), 1, 0)  ]
  
  tf.control = fread(control.major)
  setnames(tf.control, "V6", "major")
  tf.control$minor <- fread(control.minor)$V6
  tf.control[,change:=abs(major - minor)]
  tf.control[,sig:=xor((major<threshold), (minor<threshold))]
  tf.control[,disrupts:= ifelse(sig & (change > 1), 1, 0)  ]
  
  #tf.control=tf.control[V6 < threshold]
  tf.grp =  tf[, list( motif=(strsplit(V10, split=" ")[[1]][4]), eeSNP=sum(disrupts)) , by=V2]
  tf.control.grp =  tf.control[, list( control=length(V6)), by=V2]
  tf.merge = merge(tf.grp, tf.control.grp, by="V2", all=T)
  tf.merge[is.na(tf.merge)] <- 0
  tf.merge$eeSNPneg  <- eeSNPnum - tf.merge$eeSNP
  tf.merge$controlneg  <- controlnum - tf.merge$control
  eeSNP.fisher =tf.merge[, list( motif=motif, eeSNP=eeSNP, eeSNPneg=eeSNPneg, conrol = control, controlneg = controlneg, 
				greater = fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="greater")$p.value, 
				less= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="less")$p.value,
				two.sided= fisher.test(matrix(c(eeSNP , eeSNPneg, control, controlneg),2,2),alternative="two.sided")$p.value),
				by=V2]
  eeSNP.fisher$qval = fdrtool(eeSNP.fisher$two.sided, statistic="pvalue", plot=F)$qval
  setkey(eeSNP.fisher, greater)
  eeSNP.fisher
}

tfFisher.majorminor.wilcox <- function(eeSNP, control,  mat2motif,threshold=-8.5)
{
  #write.table(file="mat2motif",x = temp1 , row.names = F, 
  #col.names =T,  sep="\t", quote=F )
  require(data.table)
  require(fdrtool)
    mat2motif <- fread(mat2motif)

  nMat <- nrow(mat2motif)
  if(is.character(eeSNP)) 
  	eeSNP <- scan(eeSNP,  sep="\n")
  if(is.character(control)) 
  	control <- scan(control,  sep="\n")
  tf <- matrix(abs(eeSNP), ncol=nMat)
  tf.control <- matrix(abs(control), ncol=nMat)
  mat2motif$delta <- colMeans(tf, na.rm=T)
  mat2motif$delta.control <- colMeans(tf.control, na.rm=T)
  mat2motif$greater <- 1
  mat2motif$less <- 1
  mat2motif$two.sided <- 1
  for(ii in seq(nMat)){
    a1 <- tf[,ii]
    a2 <- tf.control[,ii]
    if(! (is.nan(mat2motif$delta[ii]) |   is.nan(mat2motif$delta.control[ii]))){ 
      mat2motif$greater[ii]<- (wilcox.test(a1, a2, alternative="greater", exact=T))$p.value 
      mat2motif$less[ii]<- (wilcox.test(a1, a2, alternative="less", exact=T))$p.value 
      mat2motif$two.sided[ii]<- (wilcox.test(a1, a2, alternative="two.sided", exact=T))$p.value
    }

  }
  #mat2motif$qval = fdrtool(mat2motif$two.sided, statistic="pvalue", plot=F)$qval
  mat2motif
}

wilcox.mat<- function(foreground, background)
{
  require(data.table)
  matout <- data.table(name=colnames(foreground))
  matout$greater <- 1
  matout$less <- 1
  matout$two.sided <- 1
  nMat <- nrow(matout)
  for(ii in seq(nMat)){
    a1 <- foreground[,ii]
    a2 <- background[,ii]
    #if(! (is.nan(matout$delta[ii]) |   is.nan(matout$delta.control[ii]))){ 
      matout$greater[ii]<- (wilcox.test(a1, a2, alternative="greater", exact=T))$p.value 
      matout$less[ii]<- (wilcox.test(a1, a2, alternative="less", exact=T))$p.value 
      matout$two.sided[ii]<- (wilcox.test(a1, a2, alternative="two.sided", exact=T))$p.value
    #}

  }
  #matout$qval = fdrtool(matout$two.sided, statistic="pvalue", plot=F)$qval
  matout
}

 if(F) {
   mat2motif <- fread(mat2motif)
nMat <- nrow(mat2motif)
   eeSNP <- SNP.major
  threshold  <-  -8.5
#inx  <- grep(mat2motif$motif, pattern="NF-AT|Nkx|GATA|STAT|FOX|Evi|MEF-2|hox|crx|myc", ignore.case=T)
inx  <- grep(mat2motif$motif, pattern="Nkx6", ignore.case=T)
  tf <- matrix(eeSNP, ncol=nMat)
  tf.bin <- matrix(0, ncol=nMat, nrow=nrow(tf))
  tf.bin[tf< threshold]  <- 1
  aa <- rowSums(tf.bin[,inx])
  eeSNP.tf.frac <- sum(aa > 0)/length(aa)

  tf.control <- matrix(control.enhancercontroled.major , ncol=nMat)
  tf.control.bin <- matrix(0, ncol=nMat, nrow=nrow(tf.control))
  tf.control.bin[tf.control< threshold]  <- 1
  bb <- rowSums(tf.control.bin[, inx])
  control.tf.frac <- sum(bb > 0)/length(bb)
 }
