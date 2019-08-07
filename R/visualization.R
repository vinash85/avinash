## adopted from Seurat.v3 visualization.R
#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression.
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param group.bar Add a color bar showing group status for cells
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @param draw.lines Include white lines to separate the groups
#' @param lines.width Integer number to adjust the width of the separating white lines.
#' Corresponds to the number of "cells" between each group.
#' @param group.bar.height Scale the height of the color bar
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work
#' when plotting multiple dimensions
#'
#' @return A ggplot object
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build aes_string
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
    object,
    features = NULL,
    cells = NULL,
    group.by = 'ident',
    group.bar = TRUE,
    disp.min = -2.5,
    disp.max = NULL,
    slot = 'scale.data',
    assay = NULL,
    label = TRUE,
    size = 5.5,
    hjust = 0,
    angle = 45,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.02,
    combine = TRUE
) {
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = object)[cells]
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    features <- features %||% VariableFeatures(object = object)
    features <- rev(x = unique(x = features))
    disp.max <- disp.max %||% ifelse(
        test = slot == 'scale.data',
        yes = 2.5,
        no = 6
    )
    # make sure features are present
    possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if(length(x = features) == 0) {
            stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
        }
        warning("The following features were omitted as they were not found in the ", slot,
                " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
    }
    data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
        object = object,
        slot = slot)[features, cells, drop = FALSE])))
    
    object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
    group.by <- group.by %||% 'ident'
    groups.use <- object[[group.by]][cells, , drop = FALSE]
    # group.use <- switch(
    #   EXPR = group.by,
    #   'ident' = Idents(object = object),
    #   object[[group.by, drop = TRUE]]
    # )
    # group.use <- factor(x = group.use[cells])
    plots <- vector(mode = 'list', length = ncol(x = groups.use))
    for (i in 1:ncol(x = groups.use)) {
        data.group <- data
        group.use <- groups.use[, i, drop = TRUE]
        if (!is.factor(x = group.use)) {
            group.use <- factor(x = group.use)
        }
        names(x = group.use) <- cells
        if (draw.lines) {
            # create fake cells to serve as the white lines, fill with NAs
            lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
            placeholder.cells <- sapply(
                X = 1:(length(x = levels(x = group.use)) * lines.width),
                FUN = function(x) {
                    return(RandomName(length = 20))
                }
            )
            placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
            group.levels <- levels(x = group.use)
            names(x = placeholder.groups) <- placeholder.cells
            group.use <- as.vector(x = group.use)
            names(x = group.use) <- cells
            group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
            na.data.group <- matrix(
                data = NA,
                nrow = length(x = placeholder.cells),
                ncol = ncol(x = data.group),
                dimnames = list(placeholder.cells, colnames(x = data.group))
            )
            data.group <- rbind(data.group, na.data.group)
        }
        plot <- SingleRasterMap(
            data = data.group,
            raster = raster,
            disp.min = disp.min,
            disp.max = disp.max,
            feature.order = features,
            cell.order = names(x = sort(x = group.use)),
            group.by = group.use
        )
        if (group.bar) {
            # TODO: Change group.bar to annotation.bar
            group.use2 <- sort(x = group.use)
            if (draw.lines) {
                na.group <- RandomName(length = 20)
                levels(x = group.use2) <- c(levels(x = group.use2), na.group)
                group.use2[placeholder.cells] <- na.group
                cols <- c(hue_pal()(length(x = levels(x = group.use))), "#FFFFFF")
            } else {
                cols <- c(hue_pal()(length(x = levels(x = group.use))))
            }
            pbuild <- ggplot_build(plot = plot)
            names(x = cols) <- levels(x = group.use2)
            # scale the height of the bar
            y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
            y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
            y.max <- y.pos + group.bar.height * y.range
            
            plot <- plot + annotation_raster(
                raster = t(x = cols[group.use2]),,
                xmin = -Inf,
                xmax = Inf,
                ymin = y.pos,
                ymax = y.max
            ) +
                coord_cartesian(ylim = c(0, y.max), clip = 'off')
            if (label) {
                x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                x.divs <- pbuild$layout$panel_params[[1]]$x.major
                x <- data.frame(group = sort(x = group.use), x = x.divs)
                label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
                label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
                plot <- plot + geom_text(
                    stat = "identity",
                    data = label.x.pos,
                    aes_string(label = 'group', x = 'label.x.pos'),
                    y = y.max + y.max * 0.03 * 0.5,
                    angle = angle,
                    hjust = hjust,
                    size = size
                )
                plot <- suppressMessages(plot + coord_cartesian(
                    ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
                    clip = 'off')
                )
            }
        }
        plot <- plot + theme(line = element_blank())
        plots[[i]] <- plot
    }
    if (combine) {
        plots <- CombinePlots(plots = plots)
    }
    return(plots)
}
