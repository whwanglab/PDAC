#' ROI Plot
#'
#' Takes data matrix, including expression, signatures, pathways and segment annotations, and creates graphs based on the ROI ID
#'
#' @param data_df a data frame of targets / pathways / signatures scores / cell types vs AOIs to be plotted
#' @param annot_df a data frame of AOI annotations necessary to plot the data
#' @param ROI_ID name of ROI column or columns, list if multiple
#' @param AOI_ID name of AOI column
#' @param segment_ann segment annotation column
#' @param segment_names labels for segment, must be labels within the segment_ann column
#' @param segment_colors colors for segments used in joy plots and heatmap, not bar
#' @param title title to be added after ROI ID to each plot
#' @param plot_type type of plot, values include 'joy', 'bar', 'heatmap'
#' @param skip_ROIs convenience to just plot legend & grab gene sets
#' @param targets subset list or logic, either vector of booleans, character vector, or function for subseting including 'CV', 'IQR', 'sd', and other standard functions
#' @param target_thresh if a function is provided, use this threshold value to pick samples with value > thresh
#' @param target_quantile whether the target_thresh is a value or a quantile to filter to
#' @param plot_clusters whether clusters should be plotted in the static graph as well
#' @param cluster_pal cluster palette to use when plotting
#' @param cluster_method clustering type, 'average', 'complete', or 'ward' recommended
#' @param cluster_dist clustering distance, 'pearson', 'spearman', or values used by dist
#' @param cluster_n number of clusters to idenitfy
#' @param cluster_min the minimum size of a cluster for it to be considered distinct, otherwise breaks will be plotted with a neighboring cluster
#' @param transf 'log2','logit','fraction','scale'
#' @param cluster_transf 'log2','logit','fraction','scale'
#' @param shape type of plot to make based on original AOI shape - 'circle' or 'rect', rectangular is still in development
#' @param hole_diameter relative center diameter for window into ROI
#' @param width relative plotting area width
#' @param draw_bgd whether to draw the background behind plot if 'joy' or 'bar' used
#' @param split_groups whether to add a buffer between clusters and draw backgrounds seperately
#' @param bgd_color color for the border around backgrounds
#' @param bgd_fill color for the background fill
#' @param joy_scale scaling factor for joy plot
#' @param joy_alpha alpha for the fill of the joy plot
#' @param joy_lwd lwd for the joyplot
#' @param heat_theme Dark' = blacks with alpha, 'Light' = whites with alpha
#' @param heat_alpha range of alpha to used with the heatmap. 0 is pure color and 1 is a total overlay of the color
#' @param heat_color color for the border around the heatmap

#' @return a list of ROI plots and legends to be subsequently saved to pdf or plotted with ggplotly
#'
#' @export plot_ROI2D
#'
#' @export legend_dendro
#' @export filter_by
#' @export CV
#' @export make_ID
#' @export pcdist
#' @export scdist
#' @export data_trans
#' @export group_small


# pdf call - put in main hydra
# counter <- 0
# pb <- txtProgressBar(min = 0, max = length(ROI_list), style = 3)
# pdf(file = '~/Projects/Hydra/v0.3/TestPDF_median_heat.pdf')
# plot_list[[ROI_list[1]]]$interactive
# plot_list$legend
# for(roi in ROI_list) {
#   print(plot_list[[roi]]$plot)
#   counter <- counter + 1
#   setTxtProgressBar(pb = pb, value = counter)
# }
# dev.off()

plot_ROI2D <- function(data_df = NULL,               # a data frame of targets / pathways / signatures scores / cell types vs AOIs to be plotted
                       annot_df = NULL,              # a data frame of AOI annotations necessary to plot the data
                       title = NULL,                 # title to be added after ROI ID to each plot
                       plot_type = 'joy',            # 'joy', 'bar', 'heatmap'
                       skip_ROIs = TRUE,             # convenience to just plot legend & grab gene sets
                       targets = 'median',           # subset list or logic, either vector of booleans, character vector, or one of 'CV', 'IQR', 'sd
                       target_thresh = 0.5,          # if CV or IQR, use this threshold value to pick samples with value > thresh
                       target_quantile = TRUE,       # whether the target_thresh is a value or a quantile to filter to
                       # not implemented      grouped = FALSE,              # whether groups of tagets are provided or not
                       plot_clusters = TRUE,         # whether clusters should be plotted in the static graph as well
                       cluster_pal = 'Set1',        # cluster palette to use when plotting
                       cluster_method = 'complete',  # clustering type, 'average', 'complete', or 'ward' recommended
                       cluster_dist = 'pearson',     # clustering distance, 'pearson', 'spearman', or values used by dist
                       cluster_n = 6,                # number of clusters to idenitfy
                       cluster_min = 20,             # if small clusters are identified (n < 20), group them and keep cutting (in testing)
                       ROI_ID = 'ROI_ID',            # name of ROI column(s, list if plural)
                       AOI_ID = 'Sample_ID',         # name of AOI column
                       segment_ann = 'Segment',      # segment annotation column
                       segment_names = c('Tumor','Stroma'),     # labels for segment, must be labels
                       segment_colors = c('green3','magenta3'), # colors for segments (used in joy, heatmap, not bar)
                       transf = NULL,                # 'log2','logit','fraction','scale'
                       cluster_transf = NULL,        # 'log2','logit','fraction','scale'
                       scale = TRUE,                 # scale from 0-1 for plotting
                       shape = 'circle',             # 'circle' or 'rect'
                       hole_diameter = 5,            # expected center diameter, relative units
                       width = 1,                    # expected donut width, relative units
                       draw_bgd = TRUE,              # whether to draw the background behind plot if 'joy' or 'bar' used
                       split_groups = TRUE,          # whether to add a buffer between groups and draw backgrounds seperately
                       bgd_color = 'white',          # color for the border around backgrounds
                       bgd_fill = alpha('black', 0.5),     # color for the background fill
                       joy_scale = 5,                # scaling factor for joy plot
                       joy_alpha = 0.6,              # alpha for the fill of the joy plot
                       joy_lwd = 0.5,                # lwd for the joyplot
                       heat_theme = 'Dark',          # 'Dark' = blacks with alpha, 'Light' = whites with alpha
                       heat_alpha = c(0,0.9),        # range of alpha to used with the heatmap (0, mean pure color shown, 1 means all color gets overlaid)
                       heat_color = alpha('white', 0.5),   # color for the border around the heatmap
                       ...) {
  cat('   Preprocessing: Checking and Pairing AOIs with ROIs\n')

  data_df <- data_df[, !colnames(data_df) %in% c("Gene","TargetName")]

  #0.0 sanity check - names of AOIs in both datasets for subsetting and in same order
  if(ncol(data_df) != nrow(annot_df)) {
    if(all(annot_df[[AOI_ID]] %in% colnames(data_df))) {
      cat('      Subsetting columns to annotated AOIs\n')
      if(length(unique(data_df[,1])) == nrow(data_df)) {
        rownames(data_df) <- data_df[,1]
      }
      data_df <- data_df[, annot_df[[AOI_ID]]]
    } else {
      stop('Error: Some AOI_IDs from the annot_df were not found in your data_df\n  Please update so they match so plots can be properly subsetted')
    }
  }

  if(segment_ann == FALSE) {
    if(length(segment_names) > 1) {
      segment_names <- segment_names[1]
      cat('       Warning: segment names longer than expected, using only the first entry\n')
    }
    if(length(segment_colors) > 1) {
      segment_colors <- segment_colors[1]
      cat('       Warning: segment names longer than expected, using only the first entry\n')
    }
  }

  #0.1 Build ROI ID from multiple columns if list provided
  if(length(ROI_ID) > 1) {
    cat('      Concatenating ROI labels together to create unique ID\n')
    annot_df$ROI_ID <- make_ID(IDs = ROI_ID, df = annot_df)
    ROI_ID <- 'ROI_ID'
  }

  #0.2 sanity check # of AOIs / sample - remove AOIs not found in list, set up colors
  ROI_list <- unique(annot_df[[ROI_ID]])
  if(segment_ann != FALSE) {
    segments <- unique(annot_df[[segment_ann]])
    if(!all(segment_names %in% segments)) {
      stop('Error: Checking segment names for ROI plots\n\nConfirm that the segment_ann column in annot_df contains the values specified in segment_names')
    }
  } else {
    segments <- segment_names
  }

  #0.3 set up color palettes to be used throughout graphs
  names(segment_colors) <- segment_names
  if(cluster_n > 8) {
    cluster_pal <- colorRampPalette(brewer.pal(8, cluster_pal))(cluster_n)
  } else {
    cluster_pal <- brewer.pal(cluster_n, cluster_pal)
  }
  names(cluster_pal) <- as.character(1:cluster_n)
  full_pal <- c(segment_colors, cluster_pal)            # use this for all colors being plotted

  # print a note to the log that ## ROIs were found to plot, with ## ROIs having all segments, and ## having a subset of segments
  # length(ROI_list), length(segments), using segment_names & colors

  #0.4 filter data
  cat('   Step 1: Filtering data to selected Targets\n')
  if(length(targets) == 1) {
    targets <- filter_by(data = data_df,
                         method = targets,
                         value = target_thresh,
                         quantile = target_quantile)
  }
  data_df <- data_df[targets, ]

  #0.5 cluster data
  # transform data if necessary haead of clustering, otherwise cluster the same way
  cat('   Step 2: Clustering selected Targets\n')

  if(!is.null(cluster_transf)) {
    clust_df <- data_trans(data_df = data_df, method = cluster_transf)
  } else {
    clust_df <- data_df
  }

  if(cluster_dist == 'pearson') {
    data_clust <- hclust(pcdist(clust_df), method = cluster_method)
  } else if(cluster_dist == 'spearman') {
    data_clust <- hclust(scdist(clust_df), method = cluster_method)
  } else {
    data_clust <- hclust(dist(x = clust_df, method = cluster_dist), method = cluster_method)
  }
  clustering <- cutree(data_clust, cluster_n)

  #2. transform data
  cat('   Step 3: Data Transformation (if selected)\n')
  if(!is.null(transf)) {
    data_df <- data_trans(data_df = data_df, method = transf)
  }

  #3. make layout for dendrogram
  cat('   Step 4: Creating plot layouts\n')

  clust_layout <- create_layout(data_clust,
                                layout = 'dendrogram', #height = data_clust$height,
                                circular = TRUE)
  if(".ggraph.orig_index" %in% colnames(clust_layout)){
    clust_layout <- create_layout(data_clust,
                                  layout = 'dendrogram', height = height,
                                  circular = TRUE)
  }
  clust_layout$cluster <- NA
  clust_layout$cluster[match(data_clust$labels,
                             clust_layout$label)] <- as.character(clustering)

  if(!".ggraph.orig_index" %in% colnames(clust_layout)){
    clust_layout$'.ggraph.orig_index' <- clust_layout$ggraph.index
  }

  # #4. rotate / align layout (position 1 at 12 noon)
  # clust_layout <- rotate_graph(clust_layout,
  #                              rotate_to = ,
  #                              rotate_which = 1,
  #                              column_which = '.ggraph.orig_index')   #### UPDATE THIS!

  #5. Plot legend
  cat('   Step 5: Creating plot legends\n')

  plot_list <- list()
  plot_list$legend <- legend_dendro(clust_layout,         # a dendrogram object from hclust
                                    #                              fan = FALSE,          # whether the dendrogram should be made into a fan
                                    #                              fan_angle = 180       # arc of fan
                                    plot_start = TRUE,    # off plot indicator of starting position for plot rotation to align with ROI plots
                                    palette = full_pal)
  plot_list$clusters <- data.frame(Gene = names(clustering),
                                   Cluster = clustering,
                                   order = data_clust$order)
  plot_list$clusters <- plot_list$clusters[plot_list$clusters$order, ]
  plot_list$clusters$row <- 1:nrow(plot_list$clusters)

  cat('   Step 6: Calculating Break Locations\n')
  breaks <- data.frame(min = sapply(1:cluster_n, function(x) {min(subset(plot_list$clusters, Cluster == x)$row)}),
                       max = sapply(1:cluster_n, function(x) {max(subset(plot_list$clusters, Cluster == x)$row)}))
  breaks$count <- breaks$max - breaks$min + 1
  breaks_group <- group_small(breaks = breaks, min_g = cluster_min, cluster_n = cluster_n)
  plot_list$clusters_grouped <- breaks_group

  #5. plot ROIs by iterating through the list of all ROIs
  if(!skip_ROIs) {
    cat('   Step 7: Creating ROI Plots\n')
    counter <- 0
    pb <- txtProgressBar(min = 0, max = length(ROI_list), style = 3)
    break_start <- ifelse(plot_type == 'joy', 2, 1)
    breaks <- data.frame(x = sapply(break_start:max(breaks_group$group), function(x) {min(subset(breaks_group, group == x)$min)}))
    for(roi in ROI_list) {
      counter <- counter + 1
      setTxtProgressBar(pb = pb, value = counter)
      # sanity check that there are no AOI annotation duplications for multiple segments
      skip <- FALSE
      ROI_df <- NULL
      if(segment_ann != FALSE) {
        segs <- annot_df[annot_df[[ROI_ID]] == roi, segment_ann]
        if(length(segs) > length(segment_colors)) {
          cat(paste('       Warning: Skipping graphs - More segments than expected for ROI:', roi, '\n'))
          skip <- TRUE
        } else if((length(unique(segs)) < length(segment_colors)) &
                  (length(segs) == length(segment_colors))) {
          cat(paste('       Warning: Skipping graphs - duplicated annotations provided for ROI:', roi, '\n'))
          skip <- TRUE
        }
      }
      if(!skip) {
        # set up plot data frames by subsetting AOIs from each ROI
        for(seg in segment_names) {
          # add logic for ROI lookup based on annotation table
          if(segment_ann != FALSE) {
            ind <- annot_df[[segment_ann]] == seg & annot_df[[ROI_ID]] == roi
          } else {
            ind <- annot_df[[ROI_ID]] == roi
          }
          if(sum(ind) != 1) {
            cat(paste0('      Warning: Skipping ',seg,' Segment from ',roi,' due to no or multiple entries\n'))
            # create flog skipping annotation
          } else {
            # create dataset for ROI plotting
            tmp_df <- data.frame(exp = data_df[data_clust$labels, ind],
                                 gene = data_clust$labels,
                                 order = clust_layout[match(data_clust$labels,
                                                            clust_layout$label),
                                                      '.ggraph.orig_index'],
                                 cluster = as.character(cutree(data_clust, cluster_n)))
            tmp_df <- tmp_df[order(tmp_df$order), ]
            tmp_df$row <- row(tmp_df)[,1]
            tmp_df$segment <- seg
            tmp_df$segment_psn <- which(segment_names == seg)
            ROI_df <- rbind(ROI_df, tmp_df)
          }
        }

        # setup ROI plots
        #### Main Plots
        # static plots for overlays
        # 5.1 - Joyplot
        if(plot_type == 'joy') {
          static_plot <- ggplot(ROI_df)
          if(draw_bgd & !(shape == 'rect')) {
            static_plot <- static_plot +
              annotate(geom = 'rect', xmin = 1, xmax = max(ROI_df$row),
                       ymin = 0.75, ymax = width * length(segment_colors) + joy_scale + 0.25,
                       fill = bgd_fill, color = bgd_color, lwd = 1)
          }
          static_plot <- static_plot +
            geom_ridgeline(aes(x = row,
                               y = segment_psn,
                               height = (exp-min(exp))/max(exp),
                               group = segment_psn,
                               fill = segment,
                               color = segment),
                           alpha = joy_alpha, min_height = -20, scale = joy_scale, lwd = 0.5)
          if(split_groups) {
            static_plot <- static_plot +
              geom_segment(data = breaks, aes(x = x-0.5, xend = x-0.5),
                           y = 0.75, yend = width * length(segment_colors) + joy_scale + 0.25,
                           color = bgd_color)
          }
          if(plot_clusters) {
            static_plot <- static_plot +
              geom_tile(aes(x = row, fill = cluster, y = 0.25), height = 0.5)
          }
          if(shape == 'circle') {
            static_plot <- static_plot + coord_polar()
          }

          #else if(shape == 'rect') {
          #
          # to be implemented: square / rectangular plots
          #
          #  static_plot <- make_square(static_plot)
          #}

          static_plot <- static_plot +
            labs(title = ifelse(!is.null(title), paste0(roi,' : ', title), roi)) +
            ylim(-2 * hole_diameter, length(segment_colors) * width + joy_scale + 0.5) +
            scale_fill_manual(values = full_pal) +
            scale_color_manual(values = full_pal) +
            xlim(1, max(ROI_df$row)) +
            theme_void() +
            theme(legend.position = 'off')

          #5.2 - heatmap
        } else if(plot_type == 'heatmap') {
          static_plot <- ggplot(ROI_df, aes(x = row)) +
            geom_tile(aes(fill = segment,
                          y = segment_psn * width,
                          height = width)) +
            geom_tile(aes(alpha = 1 - (sqrt((exp-min(exp))/max(exp))),
                          y = segment_psn * width,
                          height = width),
                      fill = ifelse(heat_theme == 'Dark', 'black','white')) +
            geom_hline(yintercept = 1 - (width/2), lwd = 1, color = heat_color) +
            geom_hline(yintercept = max(ROI_df$segment_psn)*width + (width/2), lwd = 1, color = heat_color)
          if(split_groups) {
            static_plot <- static_plot +
              geom_segment(data = breaks, aes(x = x-0.5, xend = x-0.5),
                           y = 1 - (width/2), yend = max(ROI_df$segment_psn)*width + (width/2),
                           color = heat_color, lwd = 1)
          }
          if(plot_clusters) {
            #geom_tile(aes(fill = as.character(group+2), y = 0.4), height = 0.2) +
            static_plot <- static_plot +
              geom_tile(aes(fill = cluster), y = 1 - (1.5*width), height = width/2)
          }
          if(shape == 'circle') {
            static_plot <- static_plot + coord_polar()
          }
          static_plot <- static_plot +
            ylim(-1 * hole_diameter * 1.5, width * length(segment_colors) + 1) +
            labs(title = ifelse(!is.null(title), paste0(roi,' : ', title), roi)) +
            scale_fill_manual(values = full_pal) +
            scale_alpha_continuous(range = heat_alpha) +
            theme_void() +
            theme(legend.position = "none")
          # 5.3 - barplots
        } else {
          static_plot <- ggplot(ROI_df)
          if(draw_bgd & !(shape == 'rect')) {
            static_plot <- static_plot +
              annotate(geom = 'rect', xmin = 0.5, xmax = max(ROI_df$row)+.5,
                       ymin = 0, ymax = 1.1,
                       fill = bgd_fill, color = bgd_color, lwd = 1)
          }
          static_plot <- static_plot +
            geom_col(data = ROI_df, aes(y = exp/max(exp),
                                        x = row, fill = cluster), color = 'black') +
            coord_polar() +
            ylim(hole_diameter*-.5, 1.1) +
            theme_void() +
            scale_fill_manual(values = full_pal)
          if(split_groups) {
            static_plot <- static_plot +
              geom_segment(data = breaks, aes(x = x-0.5, xend = x-0.5),
                           y = 0, yend = 1.1,
                           color = bgd_color)
          }
        }

        # due to time constraints unlikely to be able to add inderactive versions of all plots:
        # ridgeline plot for interactive exploration
        int_plot <- suppressWarnings(                    # for labels which are otherwise unused by ggplot
          ggplot(ROI_df) +
            geom_line(aes(x = row,                       # position
                          y = exp/max(exp) + segment_psn,     # to be consistently a value between 0-1
                          label = gene,                  # gene name - causes warnings if you don't use suppressWarnings
                          color = segment), # palette, since it's not a fill it's different than the tile palette
                      alpha = 0.6, lwd = 0.5) +
            geom_tile(aes(x = row, y = 0.85, color = cluster, fill = cluster, label = gene), height = 0.2) +
            theme_classic() +
            ylim(0.7, length(segment_colors) + 1) +
            scale_color_manual(values = full_pal) +
            scale_fill_manual(values = full_pal) +
            xlab('Order around circle from top, clockwise')) + ylab('Expression')
        plot_list[[roi]] <- list(plot = static_plot, interactive = int_plot, df = ROI_df)
      }
    }
  }
  #6. return list of ROI data frames for passing to ggplot/plotly
  return(plot_list)
}

filter_by <- function(data = NULL,         # df to filter
                      method = 'IQR',      # function call 'IQR', 'CV', 'sd'
                      value = 0.5,         # cut value
                      quantile = TRUE,     # whether this is a real value or a quantile
                      sign = 1) {          # sign in case you want to select negatives
  stats <- apply(data, 1, eval(method))
  if(quantile) {
    filter_to <- stats > quantile(stats, value)
  } else {
    filter_to <- stats > value * sign
  }
  return(filter_to)
}

CV <- function(data) {
  sd(data) / mean(data)
}

# Create a legend plot (1 per all ROIs as clustering should be done once)
legend_dendro <- function(layout,               # a layout from create_layout
                          fan = FALSE,          # whether the dendrogram should be made into a fan
                          fan_angle = 180,     # arc of fan
                          plot_start = TRUE,    # off plot indicator of starting position for plot rotation to align with ROI plots
# not implemented         rotate = TRUE,        # apply rotation such that 1st target from 1st cluster shows up at 12o'clock
                          palette = NULL) {     # a palette for labels, either palette name (qual list), or a palette with sequential group numbers
  # grab the # of labels
  n_targets <- nrow(subset(layout, members == 1))

  # if FAN
  if(fan) {

  }

  # start graph
  dendro <-  ggraph(layout) +
    geom_node_point(aes(filter = leaf, x = 1.65*x, y = 1.65*y),  # expand the graph to make room for labels
                    color = 'white', alpha = 0) +
    geom_edge_elbow() +theme_void() + coord_fixed()
  # if n > 100 drop text and path plotting
  if(n_targets >= 100) {
    dendro <- dendro +
      geom_node_point(aes(filter = leaf, color = cluster), size = 2, shape = 'circle')
    if(plot_start) {
      dendro <- dendro +
        geom_node_point(data = subset(layout, .ggraph.orig_index == 1),
                        aes(x = 1.09*x,
                            y = 1.09*y, filter = leaf),
                        color = 'black', size = 1, shape = 'circle') +
        geom_node_text(data = subset(layout, .ggraph.orig_index == 1),
                       aes(filter = leaf, x = x*1.1, y = y*1.1, hjust = ifelse(x>0,0,1),
                           label = 'START', angle = ifelse(x>0, round(y*8)*10, round(y*-8)*10)))
    }
  # if < n use overlay & add text
  } else {
    g_group <- subset(layout, !is.na(cluster))
    g_group <- g_group[order(g_group$.ggraph.orig_index), ]
    dendro <- dendro + geom_path(aes(x = x, y = y, color = cluster),
                                   data = g_group,
                                   lwd = 4, lineend = 'square') +
      geom_node_text(aes(label = label, filter = leaf, x = x*1.1, y = y*1.1, hjust = ifelse(x>0,0,1),
                         angle = ifelse(x>0, round(y*8)*10, round(y*-8)*10)))
  }
  dendro <- dendro + scale_color_manual(values = palette)
  return(dendro)
}

#### helper functions:
# make_ID - concatenate IDs into a new ID
make_ID <- function(IDs, df) {
  unlist(apply(df, 1, function(x) {paste(trimws(x[IDs]), collapse = '_')}))
}

# correlation distance functions
pcdist <- function(y) {as.dist(1-cor(t(y), method="pearson", use = 'complete'))}
scdist <- function(y) {as.dist(1-cor(t(y), method="spearman", use = 'complete'))}

# data_trans - just a bunch of options
data_trans <- function(data_df = NULL,
                       method = NULL,         # 'logit', 'log2', 'fraction', 'scale'
                       floor = TRUE,          # remove values less than this quanitle, helps remove artificially low values
                       floor_quant = 0.01) {  # quantile value for floor
  if(is.null(method)) {
    stop('Method for data transformation not specified, please add method\n')
  }
  if(is.null(data_df)) {
    stop('No dataset provided for data transformation, please add data\n')
  }

  if(floor) {
    floor_val <- quantile(as.matrix(data_df), floor_quant)
    data_df <- data_df - floor_val
    data_df[data_df < 0] <- 0
  }

  #transformation cases
  if(method == 'logit') {                        # logit transform if data provided as 0-1 values
    logit <- function(x) {log(x / (1-x))}
    data_df <- logit(data_df)
  #} else if(trans == 'norm') {                  # figure out !!! what to do if they pass in NES/zscore/FC values that are -3, 3

  } else if(method == 'log2') {                  # log2 transform, flooring noise values
    data_df <- log2(data_df + 1)
  } else if(method == 'fraction') {              # to fraction of max-min
    data_df <- apply(data_df, 1, function(x){
      (x-min(x)) / (max(x)-min(x))
    })
  } else {                                      # add generic
    rnm <- rownames(data_df)
    cnm <- colnames(data_df)
    data_df <- t(apply(data_df, 1, eval(method)))
    colnames(data_df) <- cnm
    rownames(data_df) <- rnm
  }
  return(data_df)
}

# group_small - right join small groups to larger groups... may not always work well
group_small <- function(breaks, min_g = 15, cluster_n = 2) {
  breaks$orig <- 1:cluster_n
  breaks <- breaks[order(breaks$min),]
  breaks$group <- NA
  group <- 1
  for(i in 1:cluster_n) {
    if(breaks[i, 'count'] > min_g) {
      breaks[i, 'group'] <- group
      group <- group+1
    } else {
      breaks[i, 'group'] <- group
    }
  }
  return(breaks)
}

# functions to be updated:
#
# make_rect <- function(plot_df,         # data frame containing segments, values, etc
#                       value = 'exp',   # column containing value to be plotted
#                       segment = 'Segment',   # column containing the segment annotations
#                       center = TRUE,   # make sure center is (0,0)
#                       hieght = 1,      # scale multiplier
#                       width = 1,       # generalize for rectangles
#                       scale = TRUE,    # scale expression based on heighest observed value... may produce odd plots if set FALSE
#                       height = 1) {    # generalize for rectangles
#   # define where the corners change
#   if(width == height) {
#     count_split <- ceil(nrow(plot_df)/4)
#   } else {
#     stop("haven't implemented rectangles yet\n")
#   }
#
#   # capture x and y,s for transmutation
#   plot_df[, c('x','y')] <- (plot_df[[value]] / max(plot_df[[value]])) * height
#   segs <- levels(plot_df[[segment]])
#   for(i in 1:length(segs)) {
#     ind <- plot_df[[segment]] == segs[[i]]
#     tmp_df <- plot_df[ind, ]
#     # top
#     tmp_df[1:count_split, 'x'] <- (1:count_split)/count_split
#     # right side
#     tmp_df[(count_split+1):(2*count_split), 'x'] <- 1 + tmp_df[(count_split+1):(2*count_split), 'x']
#     tmp_df[(count_split+1):(2*count_split), 'y'] <- -1 * (1:count_split)/count_split
#     # bottom side
#     tmp_df[(2*count_split+1):(3*count_split), 'x'] <- (count_split:1)/count_split
#     tmp_df[(2*count_split+1):(3*count_split), 'y'] <- -1 * tmp_df[(2*count_split+1):(3*count_split), 'y']
#     # left side
#     tmp_df[(3*count_split+1):(nrow(tmp_df)), 'x'] <- -1 * tmp_df[(3*count_split+1):(nrow(tmp_df)), 'x']
#     tmp_df[(3*count_split+1):(nrow(tmp_df)), 'y'] <- ((nrow(tmp_df)-(3*count_split+1)):1)/(nrow(tmp_df)-(3*count_split+1))
#     # add corner values
#     crn_df <- tmp_df[1:4,]
#     crn_df[1:4, 1:ncol(crn_df)] <- NA
#     crn_df[, segment] <- i
#     crn_df[1:4, 'x'] <- c(0,1,1,0)
#     crn_df[1:4, 'y'] <- c(0,0,-1,-1)
#     # center
#     if(center) {
#       tmp_df[, 'x'] <- tmp_df[, 'x'] - 0.5
#       tmp_df[, 'y'] <- tmp_df[, 'y'] + 0.5
#     }
#     if(i == 1) {
#       out_df <- tmp_df
#     } else {
#       out_df <- rbind(out_df, tmp_df)
#     }
#   }
# }

# rotate graph:
# rotate_graph <- function(clust_layout,
#                          rotate_to = 0,
#                          rotate_which = 1,    # which sample should be where
#                          column_which = '.ggraph.orig_index') {
#
# }  #### UPDATE THIS!
# rotate_coords - rotate x, y values by a given rotation value (in pi units)
# rotate_coords <- function(x = 1, y = 0, rotation = 0*pi) {
#
# }
