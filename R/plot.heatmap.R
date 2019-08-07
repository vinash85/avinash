#' Plot heatmap
#'
#' @param dat  matrix 
#' @param filename  (optional) 
#' @param height of pdf file
#' @param width  of pdf file
#'
#' @return
#' @export
#'
#' @examples
plot.heatmap = function(dat, filename = NULL, height = 10, width =7, scale="none", showRowDendro=F,   showColDendro=F, cexRow=0.9, cexCol=0.9){
  hc = hclust(as.dist(1-cor(dat, method="spearman", use= "pairwise.complete.obs")), method="complete")
  hr = hclust(as.dist(1-cor(t(dat), method="spearman",  use= "pairwise.complete.obs")), method="complete")
  require(heatmap3)
  p = heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale=scale, balanceColor=T, showRowDendro=showRowDendro,   showColDendro=showColDendro, cexRow = cexRow, cexCol = cexCol)
  if (!is.null(filename)) {
    pdf( filename, height = height, width =width)
    print(p)
    dev.off()
  } 
}
