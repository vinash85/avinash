## all source utlitity function for R processing 


normalize.std = function(tt){
    (tt - min(tt, na.rm=T))/(max(tt,na.rm=T) - min(tt, na.rm=T))
}

# standardise
# install.packages("data.table")
minmax <- function(x) (x - min(x))/(max(x) - min(x))

big.prcomp = function(data, center=TRUE, scale=FALSE){
	data1 = scale(data, center=center, scale=scale)
	require(bigpca)
	out = list()
	pca = big.PCA(data1,  pcs.to.keep = ncol(data1))
	out$center = center
	if(center) out$center = colMeans(data)
	out$scale = scale
	if(scale) out$scale = apply(data,2,sd)
	out$rotation = pca$PCs
	out$x = data1 %*% out$rotation
	# out$x = pca$loadings
	out
}

get_pca = function(data, pca_obj=NULL, center = T, scale = F, subsample=1, match_standard_deviation = F){
	require(gmodels)
	sds =  apply(data, 2,  sd)
	zeros.sds = which(sds == 0)
	sds = ifelse(sds==0, 1E-3, sds)
	pca_out = NULL
	# calculate pca 
	if(is.null(pca_obj)){

		if(subsample < 1){
			data.sub = data[sample(nrow(data),size = subsample*nrow(data)),]
			pca_obj <- big.prcomp(data.sub,center = center, scale=scale)

			}else{
				pca_obj <- big.prcomp(data,center = center, scale=scale)
				pca_out = pca_obj$x
			}
		pca_obj$x = NULL ## to save space
		pca_obj$sds = sds
	}
	if(is.null(pca_out)) {
		aa = t(apply(data, 1, function(tt) {
			## matching sds
			if(match_standard_deviation){ 
				tt = tt * pca_obj$sds/sds ## matching the distribution
			}
			if(pca_obj$center[1]) 
				tt = tt-pca_obj$center
			if(pca_obj$scale[1]) 
				tt = tt/pca_obj$scale
			tt
		}))
		pca_out = aa %*% pca_obj$rotation 
	}
	list(pca_out=pca_out, pca_obj=pca_obj)
}


get_sel_pca = function(data, top.genes, pca_sel_obj=NULL, center = T, scale = F, subsample=1){
	data = data[,top.genes]
	temp_out = get_pca(data, pca_sel_obj, center , scale, subsample)
	temp_out$genes = top.genes
	temp_out
}

range01.norm = function(mat){
	m1 = min(mat,na.rm=T)
	m2 = max(mat,na.rm=T)
	if(m1==m2) m2 = m2 +1
	(mat - m1)/(m2-m1)
}

	
qnorm.array <- function(mat)
{
	mat.back = mat
	mat = mat.back[!is.na(mat.back)]
    mat = rank(mat,  rank, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat
    mat.back
}

#' Title
#'
#' @param common.genes common genes to imputed 
#' @param incomplete_exp_mat incomplete expression matrix sample x genes
#' @param ref  reference expression matrix sample x genes
#'
#' @return
#' @export
#'
#' @examples
impute.closest.gene = function(common.genes, incomplete_exp_mat, ref= NULL, subsample=NULL){
	  genes.imputed = setdiff(common.genes, colnames(incomplete_exp_mat))
	  if(length(genes.imputed) > 0) {
	  gene1 = colnames(incomplete_exp_mat)
	  if(is.null(class(ref))){
	  	ref.dt = fread("/liulab/asahu/data/ssgsea/xiaoman/Genetech_expression_TPM.txt")
	  	ref = t(as.matrix(ref.dt[,seq(2,ncol(ref.dt)),with=F]))
	  	setnames(ref.dt, 1, "gene_name")
	  	colnames(ref ) = ref.dt$gene_name
	  }
	  ref.dt.gene_name = colnames(ref)
	  if (!is.null(subsample)) {
	      ref = ref[sample.int(nrow(ref), size=nrow(ref)*subsample),]
	  }
	  impute = ref[,genes.imputed, drop=F]
	  only.genes = intersect(gene1, ref.dt.gene_name)
	  incomplete_exp_mat = incomplete_exp_mat[,only.genes, drop=F]
	  exp.present = ref[,only.genes,drop=F]
	  cors = cor(impute, exp.present, use="pairwise.complete.obs")
	  # cors = pcor(impute, exp.present, use="pairwise.complete.obs")
	  genes.inx = apply(cors,1, 
	  	function(tt) ifelse(sum(!is.na(tt)), which.max(tt), NA)
	  	)

	  imputed = incomplete_exp_mat[,genes.inx, drop=F]
	  imputed[is.na(imputed)] = 0
	  colnames(imputed) = genes.imputed
	  merged = cbind(incomplete_exp_mat, imputed) 
	  incomplete_exp_mat = merged
	}
	incomplete_exp_mat[,common.genes]
}


create.km = function( times1, times2,  labels=list(Del="not-rescued", Norm="rescued"), file="temp.pdf",ppt=F){
    require(data.table)
    require(survival)
    require(ggplot2)
    dt = data.table( rbind(times1, times2)  )
    dt$label = labels$Del
    dt$label[1:nrow(times1)] = labels$Norm
    setnames(dt, 1:2, c("time", "status"))
    sr.surv <- survfit(Surv(time,status==0)~label, data = dt)
    dt$marker <- c(rep(0,nrow(times1)), rep(1,nrow(times2)))
    dt$status1 = ifelse(dt$status==0,1,0)
    max_x = 3*max(times1[,1])/4 
    outf3 = logRank(times1, times2)
    aa <- ggsurv(sr.surv)
    text = paste0(
        drug, "\n",
        "P=", formatC(outf3[1], format="E", digits=2), "\n",
        "AUC=", formatC(outf3[8] - outf3[7],digits=2)
        )
    if(ppt){
      aa <-aa +
      annotate( "text", x = max_x, y = .6, label = text) 
  }else{
      aa <-aa +
      annotate( "text", x = max_x, y = .6, label = text) +
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank())  
  }
  # pdf(file)
  # print(aa)
  # dev.off()
  aa 
}


normalize.expression = function(mat, num.samp.thr=0){ 
	require(edgeR)
	require(limma)
	mat = t(mat)
	filter = rowSums(mat > 0) >= num.samp.thr
	mat1 = mat[filter,]
	aa =  calcNormFactors(mat1, method = "TMM")
	inp.norm = mat1/aa
	lcpm <- cpm(inp.norm, log=TRUE)
	mat[filter,] = lcpm
	t(mat)
}

match.distribution = function(dat, ref){
	dens.obs = density(ref, adjust=0.1, n = length(dat))
	dens.obs$x[order(dat)]
}

match.distribution.zeros = function(dat, ref){
	## match dat with ref. 
	#The minimum values and value < ref.min are replaced by ref.min
	ref.min = min(ref)
	min.inx = (dat ==min(dat))
	dens.obs = density(ref, adjust=0.1, n = length(dat))
	out = dens.obs$x[rank(dat,ties.method="min")]
	out[out<ref.min| min.inx] = ref.min
	out
}

match.expression.distribution = function(exp, ref.exp){
	require(parallel)
	out = t(do.call(rbind, mclapply(seq(ncol(exp)), function(tt) match.distribution.zeros(exp[,tt], ref.exp[,tt]), mc.cores=32)))
	rownames(out) =rownames(exp)
	colnames(out) =colnames(exp)
	out 

}

match.expression.distribution.combat = function(exp, ref.exp, num.samp.thr=2){
	require(sva) 
	mat =rbind(exp, ref.exp)
	exp.sd = apply(exp, 2, sd) >= 0
	ref.exp.sd = apply(exp, 2, sd) >= 0
	inp.dt = t(mat)
	my.groups =factor(c(rep(1,nrow(exp)), rep(2,nrow(ref.exp))))
	nz = apply(inp.dt, 1, function(tt) length(tt[tt!=0])>=num.samp.thr)
	nz.all = nz & exp.sd & ref.exp.sd 
    batch = as.numeric(my.groups)
    nzed = inp.dt[nz.all,]
    mod1 = model.matrix(~1, data=my.groups)
    newdata = ComBat(nzed, batch=batch, mod=mod1, ref.batch =2, par.prior = TRUE, prior.plots = FALSE)
    inp.dt[nz.all,] = newdata
    new.exp = inp.dt[,seq(nrow(exp))]
    t(new.exp)
}

write.dataset = function(output.dir,dataset, ref.dir =NULL, samples.names = NULL) {
	dir.create(output.dir)
	if(!is.null(ref.dir)){
		list.of.files <- c( 
			list.files(ref.dir, "datasets_*list.txt$", full.names=T), 
			list.files(ref.dir, "datasets_tsne_list.txt$", full.names=T), 
			list.files(ref.dir, "params.json$",full.names=T) 
			)
		file.copy(list.of.files, output.dir)
	}
	if(!is.null(samples.names)){
		write.table(file=paste0(output.dir, "/samples.names.txt"), x = samples.names,
		row.names = F, col.names =T,  sep="\t", quote=F )

	}
	write.table(file=paste0(output.dir, "/dataset.txt"),x = dataset,
		row.names = F, col.names =T,  sep="\t", quote=F )
	rand_inx = sample(nrow(dataset))
	dataset_shuffle = dataset[rand_inx,]
	train.inx = 1:ceiling(.85 * nrow(dataset_shuffle))
	val.inx = ceiling(.85 * nrow(dataset_shuffle)):nrow(dataset_shuffle)
	write.table(file=paste0(output.dir, "/dataset_train.txt"),x = dataset_shuffle[train.inx,],
		row.names = F, col.names =T,  sep="\t", quote=F )
	write.table(file=paste0(output.dir, "/dataset_val.txt"),x = dataset_shuffle[val.inx,],
		row.names = F, col.names =T,  sep="\t", quote=F )
}



cohensD.na.sign = function(x,y, ...) {
    tryCatch(
        cohensD.sign(x, y, ...),
        error = function(e) NA
        )
}

cohensD.sign = function(x,y, ...){
    require(lsr)
    aa = cohensD(x, y, ...)
    if(mean(x,na.rm=T)<mean(y,na.rm=T))
        aa = -aa
    aa
}

cohensD.na = function(x,y, ...) {
  require(lsr)
    tryCatch(
        cohensD(x, y, ...),
        error = function(e) NA
        )
}
ttest.na = function(x,y, ...) {
  require(lsr)
    tryCatch(
        t.test(x, y, ...)$p.value,
        error = function(e) NA
        )
}




plot.heatmap = function(dat, filename, height = 10, width =7){
	hc = hclust(as.dist(1-cor(dat, method="spearman")), method="complete")
	hr = hclust(as.dist(1-cor(t(dat), method="spearman")), method="complete")

	require(heatmap3)
	pdf( filename, height = height, width =width)

	heatmap3(dat, Rowv=as.dendrogram(hr),  Colv=as.dendrogram(hc), scale="col", balanceColor=T, showRowDendro=T ,   showColDendro=F, cexRow = .5, cexCol = .5)

	dev.off()
}

# https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2
theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){

  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space


  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))

  #Add axis
  theme.list <- 
  list(
    theme_void(), #Empty the current theme
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
    annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
    annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
    xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
    ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
  )

  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)

  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){

    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))

    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))

    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                         data = xtick, size = size, 
                                         color = color)

    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))


    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)

    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                            x = xgeo - 2.5*epsilon, 
                                            y = ticks_y[k],
                                            size = textsize,
                                            label = paste(ticks_y[k]))
  }

  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}


match.matrix.col = function(dat, header, fill=NA) {
	"match column of dat to header and fill with NA"
	dat.out = matrix(fill, ncol=length(header), nrow=nrow(dat))
	colnames(dat.out) = header
	aa = intersect(colnames(dat),header)
	xx = match(header,colnames(dat))
	match1.inx = which(!is.na(xx))
	match2.inx = xx[match1.inx]
	# match.inx = match(colnames(dat),header)
	dat.out[,match1.inx] = as.matrix(dat[,match2.inx, with=F])
	dat.out
}
plot.complexheatmap <- function(df, name) {
    library(ComplexHeatmap)
    Heatmap(df, 
            name = name, #title of legend
            column_title = "Variables", row_title = "Samples",
            row_names_gp = gpar(fontsize = 7) # Text size for row names
    )
    
}


znorm = function(xx) (xx - mean(xx,na.rm=T))/sd(xx,na.rm=T)

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
                                    title="t-SNE",size=0.25,do.discrete=T, filename=NULL, normalize=TRUE, shape = 1){
    require(viridis)
    require(ggthemes)
    dt1 = as.data.frame(cluster)
    colnames(dt1) = c("V1", "V2")
    ps = list()
    for (color.col in color.cols) {
        color.col = gsub(color.col, pattern = "-", replacement = ".")
        if(normalize) data[[color.col]] = znorm(data[[color.col]])
        if(is.null(data$shape)) {
            d_cluster_1=cbind(dt1,col=data[[color.col]], shape=shape)
        }else{
            d_cluster_1=cbind(dt1,col=data[[color.col]], shape=data$shape)
        }
        title.curr = sprintf("%s_%s", title, color.col)
        ## plotting the results without clustering
        p=ggplot(d_cluster_1, aes(x=V1, y=V2)) +
            geom_point(size=size,aes(color=col, shape=as.factor(shape)), alpha=0.7) +
            # scale_color_gradient2(low = "blue", mid = "white",
            # high = "red", space = "Lab" ) +
            
            # guides(colour=guide_legend(override.aes=list(size=2))) +
            xlab("Dim1") + ylab("Dim2") +
            ggtitle(label = title.curr) +
            theme_light(base_size=20) +
            theme(axis.text.x=element_blank(),
                  axis.text.y=element_blank()) 
        if (do.discrete) {
            p<- p+ scale_colour_brewer(palette = "Set2")
        }else{
            # p <- p+  scale_color_viridis() 
            p <- p+ scale_color_gradientn(colours = heat.colors(20, alpha=0.7, rev=T))
            
            # p <- p +   
            #     scale_colour_gradient_tableau(palette="Classic Red") +
            #     theme_classic() + theme(legend.position = "none") 
            # 
        }
        ##theme(legend.position = "none")
        if(!is.null(filename)) {
            filename.curr = sprintf("%s_%s.pdf", filename, gsub(color.col, pattern="-", replacement = "_"))
            
            ggsave(file=filename.curr, p)
        }
        ps[[color.col]]  = p
        
    }
    ps
}

