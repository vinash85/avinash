### 
#count.tpm
#
#
readcount2tpm <- function(count, gene.length){
	gene.length = gene.length/1E3
	count.per.kb = count/gene.length
	tot = sum(count.per.kb, na.rm=T)
	tpm.fac =tot/1E6
	count.per.kb/tpm.fac
}