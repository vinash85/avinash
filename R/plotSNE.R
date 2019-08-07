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



plotUMAP <- function(data, col, color.col = "condition",
                     title="UMAP",size=0.25, do.discrete=T, filename=NULL, n_neighbors =  15, learning_rate = 1, init = "spectral", min_dist = .01, pca = NULL,  n_threads=32){
    
    require(uwot)
    require(ggthemes)
    set.seed(9)
    
    umap.model <- umap(as.matrix(data[,col]), 
                       pca = pca, n_neighbors =  n_neighbors, learning_rate = learning_rate, init = init, min_dist = min_dist, n_threads = n_threads, ret_model=T)
    
    
    d_umap_1 = as.data.frame(umap.model$embedding)
    d_umap_1=cbind(d_umap_1,data[,color.col,drop=F])
    
    ## plotting the results without clustering
    p=ggplot(d_umap_1, aes(x=V1, y=V2)) +
        geom_point(size=size,aes_string(color=color.col), alpha=0.8) +
        guides(colour=guide_legend(override.aes=list(size=2))) +
        xlab("umap_1") + ylab("umap_2") +
        ggtitle(label = title) +
        # theme_light(base_size=20) +
  scale_color_tableau() + 
  theme_classic() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank()) 
    if (do.discrete) {
        p<- p+ scale_colour_brewer(palette = "Set2")
    }
    ##theme(legend.position = "none")
    if(!is.null(filename)) ggsave(file=filename, p)
    list(umap.model, p)
}
