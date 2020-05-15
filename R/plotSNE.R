plotSNE <- function(data= data_tsne.merge, col, color.col = "condition",
                    title="t-SNE",size=0.25,do.discrete=T, filename=NULL, perplexity=30, theta=0.5, pca = FALSE, max_iter=5000, num_threads=32, tsne =NULL){
  library(Rtsne)
  library(ggplot2)
  set.seed(9)
  require(ggplot2)
  require(ggthemes)
  
  if(is.null(tsne)) tsne <- Rtsne(as.matrix(data[,col]), check_duplicates = FALSE, 
                                  pca = pca, perplexity=perplexity, theta=theta, dims=2, max_iter = max_iter, num_threads = num_threads)
  
  d_tsne_1 = as.data.frame(tsne$Y)
  d_tsne_1=cbind(d_tsne_1,data[,color.col,drop=F])
  
  ## plotting the results without clustering
  p=ggplot(d_tsne_1, aes(x=V1, y=V2)) +
    geom_point(size=size,aes_string(color=color.col), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("tSNE_1") + ylab("tSNE_2") +
    ggtitle(label = title) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  if (do.discrete) {
    p<- p+ scale_colour_brewer(palette = "Set2")
  }
  ##theme(legend.position = "none")
  if(!is.null(filename)) ggsave(file=filename, p)
  list(tsne, p)
}



plotUMAP <- function(data, umap.model = NULL, col = NULL, color.col = NULL,
                     title="UMAP",size=0.25, do.discrete=T, filename=NULL, n_neighbors =  15, learning_rate = 1, init = "spectral", min_dist = .01, pca = NULL,  n_threads=32, ...){
  require(uwot)
  require(ggthemes)
  require(ggplot2)
  if(is.null(umap.model)){
    set.seed(9)
    if(is.data.table(data)) data = as.data.frame(data)
    if(!is.null(col)) {
      umap.mat = as.matrix(data[,col])
    }else{
      umap.mat = as.matrix(data)
    }
    umap.model <- umap(umap.mat, 
                       pca = pca, n_neighbors =  n_neighbors, learning_rate = learning_rate, init = init, min_dist = min_dist, n_threads = n_threads, ret_model=T, ...)
  }
  
  d_umap_1 = as.data.frame(umap.model$embedding)
  d_umap_1$col = as.factor(0)
  if(!is.null(color.col)) {
    if(length(color.col) == nrow(data)){
      d_umap_1$col=color.col
    }else{
      d_umap_1$col=data[[color.col]]
    }
      
    # if(is.character(color.col)){
    #   d_umap_1$col=data[[color.col]]
    # }else if(length(color.col) == nrow(data)){
    #   d_umap_1$col=color.col
    # }
  }
  
  do.discrete = ifelse(is.factor(d_umap_1$col), T, F)
  ## plotting the results without clustering
  p=ggplot(d_umap_1, aes(x=V1, y=V2)) +
    geom_point(size=size,aes(color=col), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("umap_1") + ylab("umap_2") +
    ggtitle(label = title) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  if (do.discrete) {
    p<- p+ scale_colour_brewer(palette = "Set2")
  }else{
    p <- p+ scale_color_gradient2_tableau(palette = "Orange-Blue Diverging")
    # p <- p+ scale_color_gradient2_tableau(palette = "Orange-Blue Diverging", limits=c(-.3, .3))
  }
  
  ##theme(legend.position = "none")
  if(!is.null(filename)) ggsave(file=filename, p)
  list(umap.model, p)
}


plot.sequence.and.clustering <- function(data.clust, text, color.col = NULL, num.plot.seq =100, text.size =2) {
  require(viridis)
  require(ggthemes)
  require(ggrepel)
  # data.tcr.clust = data.table(umap.all.p[[1]]$embedding)
  
  d_umap_1 = cbind(data.clust, text = text)
  if(is.null(color.col)) {
    d_umap_1$col = as.factor(0)
  }else{
    d_umap_1$col=color.col
  }
  do.discrete = ifelse(is.factor(d_umap_1$col), T, F)
  num.plot.seq = min(num.plot.seq, nrow(d_umap_1))
  
  p=ggplot(d_umap_1, aes(x=V1, y=V2)) +
    geom_point(aes(color=col), alpha=0.8) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("umap_1") + ylab("umap_2") +
    # ggtitle(label = title) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) 
  
  if (do.discrete) {
    p<- p+ scale_colour_brewer(palette = "Set2")
  }else{
    p <- p+ scale_color_tableau()
  }
  
  p = p + geom_text_repel(
    # data.tcr = d_umap_1[sample(nrow(d_umap_1), size=100)],
    data = d_umap_1[sample(nrow(d_umap_1), size=num.plot.seq)],
    aes(x=V1, y=V2, label = text),
    size = text.size,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
    theme(legend.position = "none")
  p
}



plotSNE.array <- function(data= data_tsne.merge,col = c(2:114), color.cols = "condition",
                          title="t-SNE",size=0.25,do.discrete=T, filename=NULL, perplexity=30, theta=0.5, pca = FALSE, max_iter=5000, normalize=TRUE, num_threads=32){
  set.seed(9)
  
  tsne <- Rtsne(as.matrix(data[,col]), check_duplicates = FALSE, 
                pca = pca, perplexity=perplexity, theta=theta, dims=2, max_iter = max_iter, num_threads = num_threads)
  dt1 = as.data.frame(tsne$Y)
  ps = list()
  for (color.col in color.cols) {
    if(normalize) data[[color.col]] = znorm(data[[color.col]])
    d_tsne_1=cbind(dt1,col=data[[color.col]], shape=data$shape)
    title.curr = sprintf("%s_%s", title, color.col)
    
    ## plotting the results without clustering
    p=ggplot(d_tsne_1, aes(x=V1, y=V2)) +
      geom_point(size=size,aes(color=col, shape=as.factor(shape)), alpha=0.8) +
      scale_color_gradient2(low = "blue", mid = "white",
                            high = "red", space = "Lab" ) + 
      guides(colour=guide_legend(override.aes=list(size=2))) +
      xlab("tSNE_1") + ylab("tSNE_2") +
      ggtitle(label = title.curr) +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank()) 
    if (do.discrete) {
      p<- p+ scale_colour_brewer(palette = "Set2")
    }
    ##theme(legend.position = "none")
    if(!is.null(filename)) {
      filename.curr = sprintf("%s_%s.pdf", filename, color.col)
      
      ggsave(file=filename.curr, p)
      ps[[color.col]]  = p
    }
  }
  list(d_tsne_1, ps)
}


color.clusters.features <- function(data, cluster,  color.cols = "condition",
                                    title="t-SNE",size=0.25, do.discrete=T, filename=NULL, normalize=TRUE, shape = 1, remove.guide=T){
  require(viridis)
  require(ggthemes)
  dt1 = as.data.frame(cluster)
  colnames(dt1) = c("V1", "V2")
  ps = list()
  if(!is.null(filename)) dir.create(filename)
  for (color.col in color.cols) {
    tryCatch({
      color.col = gsub(color.col, pattern = "-", replacement = ".")
      if(normalize) data[[color.col]] = znorm(data[[color.col]])
      if(is.null(data$shape)) {
        d_cluster_1=cbind(dt1,col=data[[color.col]], shape=shape)
      }else{
        d_cluster_1=cbind(dt1,col=data[[color.col]], shape=data$shape)
      }
      title.curr = color.col
      if(title !="") title.curr = sprintf("%s_%s", title, color.col)
      
      do.discrete = ifelse(is.factor(d_cluster_1$col), T, F)
      p=ggplot(d_cluster_1, aes(x=V1, y=V2)) +
        geom_point(size=size,aes(color=col, shape=as.factor(shape)), alpha=0.7) +
        xlab("Dim1") + ylab("Dim2") +
        ggtitle(label = title.curr) +
        theme_classic() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      
      if (do.discrete) {
        colorCount = length(unique(d_cluster_1$col))
        if(colorCount > 8){
          getPalette = colorRampPalette(brewer.pal(9, "Set1"))
          p =  p + scale_fill_manual(values = getPalette(colorCount))
        }else{
          p<- p+ scale_colour_brewer(palette = "Set2")
        }
      }else{
        p <- p+ scale_color_gradientn(colours = heat.colors(20, alpha=0.7, rev=T))
      }
      if (remove.guide)  p<- p + theme(legend.position = "none")
      if(!is.null(filename)) {
        filename.curr = sprintf("%s_%s.pdf", filename, gsub(color.col, pattern="-", replacement = "_"))
        
        ggsave(file=filename.curr, p)
      }
      ps[[color.col]]  = p
    },
    error = function(e) ps[[color.col]]  = NA
    )
    
  }
  ps
}

plot.biauc <- function(ps, dir, indexes, cell.types, mat, ncol=3) {
  require(ggrepel)
  names.ps = names(ps)
  biauc.ps = dts = list() 
  for (tt in seq_along(ps) ) {
    p = ps[[tt]]$data
    name = p$name = names.ps[tt]
    title = gsub(ps[[tt]]$label$title, pattern = " ", replacement =  ".")
    p$title = title 
    xx = title
    curr.marker = p$marker
    responder.index.curr = intersect(indexes, which(cell.types==name & response=="Responder"))
    nonresponder.index.curr = intersect(indexes, which(cell.types==name & response!="Responder"))
    auc.pos = colMeans(mat[responder.index.curr,curr.marker], na.rm = T) > colMeans(mat[nonresponder.index.curr,curr.marker], na.rm = T)
    p[,biaucs:=ifelse(auc.pos, aucs-0.5, -(aucs-0.5))]
    dts[[name]] = p 
    
    di.aucs.dt = p[label=="signature"]
    aucs.dt = p[label=="gene"]
    
    m1 = di.aucs.dt[which(aucs > 0.7)]
    if(nrow(m1) > 20) m1 = di.aucs.dt[order(aucs,decreasing=T)][1:25]
    if(nrow(m1) < 2) m1 = di.aucs.dt[order(aucs,decreasing=T)[1:5]]
    m2 = aucs.dt[which(aucs > 0.7)]
    if(nrow(m2) > 20) m2 = aucs.dt[order(aucs,decreasing=T)[1:20]]
    if(nrow(m2) < 2) m2 = aucs.dt[order(aucs,decreasing=T)[1:5]]
    
    p[,sub:=ifelse(marker %in% c(m1$marker, m2$marker), 1, 0)]
    p1 = ggplot(p, aes(x = biaucs, y = logP)) +
      geom_point(aes(color=as.factor(label), alpha = alpha)) +
      
      theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
      labs(x="AUC - 0.5", y="Significance", title=title)+
      geom_text_repel(
        data = subset(p, sub==1 ),
        aes(x = biaucs, y = logP, label = marker),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )  
    
    biauc.ps[[name]] = p1
    filename = sprintf("%s/bi_aucs_%s.pdf",dir, title)
    ggsave(p1, file=filename, width =7, height = 7)
    
  }
  
  dts = do.call(rbind, dts)
  
  title1 = strsplit(title, split = "\\.")[[1]][1]
  auc.p = ggplot(dts, aes(x = aucs, y = logP)) +
    geom_point(aes(color=as.factor(label), alpha = alpha)) +
    
    theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
    labs(x="AUC", y="Significance", title=title1)+
    geom_text_repel(
      data = subset(dts, sub==1 ),
      aes(x = aucs, y = logP, label = marker),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )   + facet_wrap(~name, ncol=ncol, scales = "free")
  filename = sprintf("%s/aucs_%s.pdf",dir, title1)
  ggsave(auc.p, file=filename, width =20, height = 15)
  
  
  biauc.p = ggplot(dts, aes(x = biaucs, y = logP)) +
    geom_point(aes(color=as.factor(label), alpha = alpha)) +
    
    theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
    labs(x="AUC - 0.5", y="Significance", title=title1)+
    geom_text_repel(
      data = subset(dts, sub==1 ),
      aes(x = biaucs, y = logP, label = marker),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )   + facet_wrap(~name, ncol=ncol, scales = "free")
  filename = sprintf("%s/biaucs_%s.pdf",dir, title1)
  ggsave(biauc.p, file=filename, width =20, height = 15)
  
  list(biauc.p =biauc.p, biauc.ps=biauc.ps,auc.p=auc.p,  dts=dts)
}
