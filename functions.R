# install.packages("tidyverse")
# install.packages("tidymodels")
# install.packages("sp")
# install.packages("scales")
# install.packages("cowplot")
devtools::install_github("sharlagelfand/dmc")
library(imager)
library(tidyverse)
library(tidymodels)
library(sp)
library(scales)
library(cowplot)
library(dmc)

change_resolution <- function(image_df, x_size) {
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}

############################################### Part (a) ###############################################
process_image <- function(image_file_name, k_list) {
  ## Runs clustering over an image.
  ## This function takes an image (PNG or JPG) and a list of the different number of possible
  ## clusters as input. Then it converts the image to a data frame with the RGB values of the 
  ## pixels and runs clustering on the image with the different values of 'k' i.e. the number of
  ## clusters as the values from k_list.
  ##
  ## Input:
  ##  – image_file_name - a PNG or JPEG image.
  ##  – out - the number of centres in the clustering
  ## Output:
  ##  – cluster_info: A list or tibble of information derived from the k_means that 
  ##                  will be sufficient to be the input to any other function you 
  ##                  write. This function computes a clustering. 
  ##                  This includes at least:
  ##                  * the original output of the kclust calls,
  ##                  * the tidied clusters, their associated RGB values and their 
  ##                    nearest DMC thread colour information.
  ## Example:
  ##   library(imager)
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(dmc)
  ##   cluster_info <- process_image("~/Downloads/d5496755a.jpg", c(2,4,6,8,10))
  
  # Load image and convert it to a data frame
  im <- imager::load.image(image_file_name)
  tidy_dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  dat <- select(tidy_dat,c(-x,-y))
  # Run k-means clustering on our image with 'k' as the values from 'k_list'
  kclusts <- tibble(k = k_list) %>% mutate(kclust = map(k, ~kmeans(x = dat , centers = .x, nstart=25)),
                                           tidied = map(kclust, tidy),
                                           glanced = map(kclust, glance),
                                           augmented = map(kclust, augment, tidy_dat))
  # The data separated by clusterings
  clusterings <- kclusts %>% unnest(cols = c(augmented))
  clusterings <- select(clusterings, c(-tidied, -glanced))
  # The cluster centers
  tidied_clusters <- kclusts %>% unnest(cols = c(tidied))
  tidied_clusters <- select(tidied_clusters, c(-augmented, -glanced))
  tidied_clusters <- tidied_clusters %>% mutate(col = rgb(R,G,B))
  tidied_clusters <- tidied_clusters %>% mutate(DMC = map(col, ~dmc(.x, visualize = FALSE)))
  tidied_clusters <- tidied_clusters %>% unnest(cols = c(DMC))
  tidied_clusters <- select(tidied_clusters, c(-red, -green, -blue))
  tidied_clusters <- tidied_clusters %>% rename("dmc_hex_value" = hex)
  # The single-row summary for all clusters
  cluster_summaries <- kclusts %>% unnest(cols = c(glanced))
  cluster_summaries <- select(cluster_summaries, c(-augmented, -tidied))
  # Combine the original kclust calls output, clusterings of the data, the tidied clusters
  # and the one line summary of each cluster into a list
  cluster_info <- list(kclusts, tidied_clusters, cluster_summaries, clusterings)
  return(cluster_info)
}

############################################### Part (b) ############################################### 
scree_plot <- function(cluster_info) {
  ## Returns a scree plot of the clustering performed on an image
  ##
  ## Input:
  ##  – cluster_info: A list of information derived process_image function, which gives
  ##                  us information about the k-means clustering for different possible 'k' with
  ##                  the input image, the sum of squares summaries for each clustering, and the
  ##                  cluster centers. 
  ##                  The list that cluster_info is, is shaped as follows:
  ##                  * list("original kclust calls output", "tidied clusters", "glanced data for
  ##                         clustering", "augmented clustering of data")
  ## Output:
  ##  – scree_plot: A scree plot of the total within-cluster sum of squares, i.e. tot.withinss.
  ## Example:
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   cluster_info <- process_image("~/Downloads/d5496755a.jpg", c(2,4,6,8,10))
  ##   scree_plot(cluster_info)
  
  # Get the clustering summary for all 'k' from cluster_info which was obtained using "glance()"
  cluster_summaries <- cluster_info[[3]]
  # Plot the scree plot for the cluster summaries
  scree_plot <- ggplot(cluster_summaries, aes(k, tot.withinss)) + 
                  xlab("Number of clusters - k") +
                  ylab("Total within-cluster sum of squares") +
                  geom_line() +
                  geom_point()
  return(scree_plot)
}

############################################### Part (c) ###############################################
colour_strips <- function(cluster_info) {
  ## Returns colour strips for the cluster centers for different 'k'
  ##
  ## Input:
  ##  – cluster_info: A list of information derived process_image function, which gives
  ##                  us information about the k-means clustering for different possible 'k' with
  ##                  the input image, the sum of squares summaries for each clustering, and the
  ##                  cluster centers. 
  ##                  The list that cluster_info is, is shaped as follows:
  ##                  * list("original kclust calls output", "tidied clusters", "glanced data for
  ##                          clustering", "augmented clustering of data")
  ## Output:
  ##  – colour_strips: colour strips for the cluster centers of clusterings for different 'k'
  ## Example:
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(scales)
  ##   cluster_info <- process_image("~/Downloads/d5496755a.jpg", c(2,4,6,8,10))
  ##   colour_strips(cluster_info)
  
  tidied_clusters <- cluster_info[[2]]
  k_list <- unique(tidied_clusters$k)
  k_list_len <- length(k_list)
  col_strip_list <- replicate(k_list_len, data.frame())
  i <- 1
  while (i < k_list_len+1) {
    new <- subset(tidied_clusters, k == k_list[i])
    new <- select(new, c(k, dmc, name, dmc_hex_value))
    col_strip_list[[i]] <- new
    i <- i + 1
  }
  par(mfrow = c(ceiling(k_list_len/2), 2))
  for (df in col_strip_list) {
   show_col(df$dmc_hex_value)
  }
}

############################################### Part (d) ###############################################
make_pattern <- function(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL) {
  ## Plots the cross-stitch pattern of the image that was clustered by cluster_info.
  ##
  ## Input:
  ##  – cluster_info - The output of process_image which is a list of information about the k-means 
  ##                   clustering for different possible 'k' with the input image, the sum of squares
  ##                   summaries for each clustering, and the cluster centers. 
  ##                   The list that cluster_info is, is shaped as follows:
  ##                   * list("original kclust calls output", "tidied clusters", "glanced data for
  ##                           clustering", "augmented clustering of data")
  ##  - k - The chosen cluster size
  ##  - x_size - The (approximate) total number of possible stitches in the horizontal direction
  ##  - black_white - (logical) Print the pattern in black and white (TRUE) or colour (FALSE, default)
  ##  - background_colour - The colour of the background, which should not be stitched in the
  ##                        pattern. (Default is to not have a colour)
  ## Output:
  ##  – This function is the only function that uses change_resolution(image_df, x_size). It should 
  ##    produce a cross-stitch pattern that can be followed, complete with a legend that has thread 
  ##    colour, and a guide grid.
  ## Example:
  ##   library(tidyverse)
  ##   library(tidymodels)
  ##   library(sp)
  ##   cluster_info <- process_image("~/Downloads/d5496755a.jpg", c(2,4,6,8,10))
  ##   make_pattern(cluster_info, 6, 50, FALSE, "#CF09E8")
  
  chosen_k <- k
  # Get the augmented data from our clustering info we got from process_image
  # and change the resolution of the image to a lower resolution using the 
  # change_resolution function.
  augmented_data <- cluster_info[[4]]
  augmented_data_for_chosen_k <- subset(augmented_data, k == chosen_k)
  augmented_data_for_chosen_k <- select(augmented_data_for_chosen_k, c(k, x, y, R, G, B, .cluster))
  # Get the tidied clusters for the chosen k, i.e. get the cluster centers
  centers_for_k_clusters <- subset(cluster_info[[2]], k == chosen_k)
  centers_for_k_clusters <- centers_for_k_clusters %>% mutate(name = paste(name, '(', dmc, ')'))
  # Check if background color argument is not NULL and if not NULL remove all points
  # of background color in data
  if(!is.null(background_colour)) {
    background_colour_cluster <- subset(centers_for_k_clusters, dmc_hex_value==background_colour)$cluster[1]
    centers_for_k_clusters <- subset(centers_for_k_clusters, cluster != background_colour_cluster)
    augmented_data_for_chosen_k <- subset(augmented_data_for_chosen_k, .cluster != background_colour_cluster)
  }
  # Pass augmented data for 'k' clusters to change res to get a low res image.
  low_res_img_data <- change_resolution(augmented_data_for_chosen_k, x_size)
  low_res_img_data <- low_res_img_data %>% rename("cluster" = .cluster)
  low_res_cluster_nums <- unique(low_res_img_data$cluster)
  # Create a cross-stitch for the image
  cluster_frame <- tibble(clust = centers_for_k_clusters$cluster, name = centers_for_k_clusters$name, 
                          col = centers_for_k_clusters$dmc_hex_value)
  cluster_frame <- subset(cluster_frame, clust %in% low_res_cluster_nums)
  if(black_white==TRUE) {
    cross_stitch <- low_res_img_data %>% ggplot(aes(x, y)) + geom_point(aes(shape = factor(cluster), col = NULL)) +
                                          scale_colour_manual(name = "",
                                                              values = NULL,
                                                              label =  cluster_frame %>% select(clust, name) %>% deframe) +
                                          scale_shape_manual(name = "",
                                                             values = low_res_cluster_nums,
                                                             labels  = cluster_frame %>% select(clust, name) %>% deframe) + 
                                          theme(legend.position="bottom", legend.direction = "horizontal") +
                                          scale_y_reverse() + theme_void()
  } else {
    cross_stitch <- low_res_img_data %>% ggplot(aes(x, y)) + geom_point(aes(shape = factor(cluster), col = factor(cluster))) +
                                          scale_colour_manual(name = "",
                                                              values = cluster_frame %>% select(clust, col) %>% deframe,
                                                              label =  cluster_frame %>% select(clust, name) %>% deframe) +
                                          scale_shape_manual(name = "",
                                                             values = low_res_cluster_nums,
                                                             labels  = cluster_frame %>% select(clust, name) %>% deframe) + 
                                          theme(legend.position="bottom", legend.direction = "horizontal") +
                                          scale_y_reverse() + theme_void()
  }
  return(cross_stitch)
}



