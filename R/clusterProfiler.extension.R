compareCluster.unverselist = function (geneClusters, fun = "enrichGO", data = "", universe.dt = NULL, universe.var=NULL,  ...) 
{
    # require(plyr)
    fun_name <- fun
    fun <- eval(parse(text = fun))

    if (typeof(geneClusters) == "language") {
        if (!is.data.frame(data)) {
            stop("no data provided with formula for compareCluster")
        }
        else {
            genes.var = all.vars(geneClusters)[1]
            grouping.formula = gsub("^.*~", "~", as.character(as.expression(geneClusters)))
            if (is.data.frame(data) && grepl("+", grouping.formula)) {
                groupVarName <- strsplit(grouping.formula, split = "\\+") %>% 
                    unlist %>% gsub("~", "", .) %>% gsub("^\\s*", "", .) %>% gsub("\\s*$", "", .)
                for (gg in groupVarName) {
                    data[[gg]]  = data[[gg]] %>% gsub("\\.", ., replacement = "_")
                    universe.dt[[gg]]  = universe.dt[[gg]] %>% gsub("\\.", ., replacement = "_")
                }
                                                         
            }
            geneClusters = plyr::dlply(.data = data, formula(grouping.formula), 
                                 .fun = function(x) {
                                     as.character(x[[genes.var]])
                                 })
            universeList = plyr::dlply(.data = universe.dt, formula(grouping.formula), 
                                 .fun = function(x) {
                                     as.character(x[[universe.var]])
                                 })
        }
    }
    stopifnot(all(names(geneClusters) %in% names(universeList)))
    clProf <- plyr::llply(names(geneClusters), .fun = function(i) {
        geneCluster = geneClusters[[i]]
        universe = universeList[[i]]
        x = suppressMessages(fun(geneCluster, universe=universe, ...))
        if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
            x= as.data.frame(x)
        }
        # browser()
        if(nrow(x) > 0) 
            x$Cluster = i
        x
    })
    clusters.levels = names(geneClusters)
    clProf.df <- plyr::ldply(clProf, rbind)
    if (nrow(clProf.df) == 0) {
        stop("No enrichment found in any of gene cluster, please check your input...")
    }
    # clProf.df <- plyr::rename(clProf.df, c(.id = "Cluster"))
    # clProf.df$Cluster = cluster.levels
    clProf.df$Cluster = factor(clProf.df$Cluster, levels = clusters.levels)
    if (is.data.frame(data) && grepl("+", grouping.formula)) {
        # groupVarName <- strsplit(grouping.formula, split = "\\+") %>% 
        #     unlist %>% gsub("~", "", .) %>% gsub("^\\s*", "", 
        #                                          .) %>% gsub("\\s*$", "", .)
        # browser()
        groupVars <- sapply(as.character(clProf.df$Cluster), 
                            strsplit, split = "\\.") %>% do.call(rbind, .)
        for (i in seq_along(groupVarName)) {
            clProf.df[, groupVarName[i]] <- groupVars[, i]
        }
        i <- which(colnames(clProf.df) %in% groupVarName)
        j <- (1:ncol(clProf.df))[-c(1, i)]
        clProf.df <- clProf.df[, c(1, i, j)]
    }
    new("compareClusterResult", compareClusterResult = clProf.df, 
        geneClusters = geneClusters, fun = fun_name, .call = match.call(expand.dots = TRUE))
}
