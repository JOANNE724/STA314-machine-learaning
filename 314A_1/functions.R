library(imager)
library(tidyverse) 
library(tidymodels) 
library(sp) 
library(scales)
library(devtools)
library(cowplot)
install_github("sharlagelfand/dmc")
library(dmc)
library(dplyr)
library(ggplot2)


change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ##  a lower resolution image. Any non-coordinate columns in the data
  ##  frame are summarized with their most common value in the larger
  ##  grid cell.
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

process_image <- function(image_file_name, k_list){

  ##  process_image(image_file_name, k_list) shows a list of cluster information
  ##   that we will use next for the image---image_file_name we uploaded. Here, 
  ##   k_list means the number of centres in the clustering.
  ##
  ## Input:
  ##   image_file_name - a PNG or JPEG image.
  ##   k_list: - the number of centres in the clustering 
  ## Output:
  ##   cluster_info: A list or tibble of information derived from the k_means 
  ##   that will be sufficient to be the input to any other function you write. 
  ##   This includes: 
  ##      The original output of the kclust calls,
  ##      The tidied clusters, their associated RGB values and their nearest DMC 
  ##      thread colour information.
  ##   
  ## Example: 
  ##      image_file <- "img.jpg" 
  ##      my_clusters <- process_image(image_file_name = image_file,k_list=2:10)
  ##
  im <- load.image("img.jpg")
  tidy_dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  
  dat <- select(tidy_dat, c(-x,-y))
  
  
  #   Let’s say we want to explore the effect of different choices of k, from 
  # 2 to 10, on this clustering. First cluster the data 9 times, each using a 
  # different value of k.
  kclusts <-
    tibble(k=c(k_list)) %>%
    mutate(
      kclust = map(k,~kmeans(x=dat, centers = .x, nstart = 4)),
      glanced = map(kclust, glance),
    )
  

  clusterings <- 
    kclusts %>% 
    unnest(cols = c(glanced)) %>%
    mutate(centres = map(kclust, tidy),
           tidy_dat = map(kclust, ~augment(.x, tidy_dat) %>% 
                               rename(cluster = .cluster)),
           )
  
  len <- length(k_list)
  for (i in 1:len){
    clusterings$centres[[i]] <- clusterings$centres[[i]] %>% 
      mutate(col = rgb(R,G,B), dmc = map(col, ~dmc(.x))) %>% select(-col)
  }
  
  return(clusterings)

}







scree_plot <- function(cluster_info)
{
  ## scree_plot(cluster_info) produces and plots a scree plot derived from
  ##   cluster_info. A scree plot plots the k-means objective function as a 
  ##   function of k.The production of this function will help us judge 
  ##   and specify a exact number k we should choose.
  ##
  ## Input:
  ## - cluster_info: Output of the process_image. A list or tibble of 
  ##   information derived from the k_means that will be sufficient to be the 
  ##   input to any other function you write. 
  ##   This includes: 
  ##    - The original output of the kclust calls,
  ##    - The tidied clusters, their associated RGB values and their nearest
  ##      DMC thread colour information.
  ## Output:
  ## - scree_plot: A ratio plot derived from cluster_info, which shows the 
  ##   k-means objective function as a function of k. 
  ##
  ## Example: 
  ##      scree_plot(cluster_info = my_clusters)
  ##
  nclust = length(cluster_info$k)
  ratio = rep(NA, nclust-1)
  for (kk in 2:nclust) {
    ratio[kk-1] = cluster_info$tot.withinss[kk]/cluster_info$tot.withinss[kk-1]
  }
  plot_data <- data.frame(k = cluster_info$k[2:nclust],ratio)
  ratio_plot <- ggplot(plot_data, aes(x=k, y = ratio)) + geom_line()
  
  return(ratio_plot)
}





colour_strips <- function(cluster_info)
{
  ## colour_strips(cluster_info) produces colour strips with the DMC
  ##  colour closest to the cluster centre colour, wihch is derived from
  ##  cluster_info. 
  ##
  ## Input:
  ## - cluster_info: Output of the process_image. A list or tibble of 
  ##   information derived from the k_means that will be sufficient to be the 
  ##   input to any other function you write. 
  ##   This includes: 
  ##    - The original output of the kclust calls,
  ##    - The tidied clusters, their associated RGB values and their nearest
  ##      DMC thread colour information.
  ## Output:
  ## - colour_strips: A plot shows the colour strips with the DMC colour
  ##   closest to the cluster centre colour.
  ##
  ## Example: 
  ##      colour_strips(cluster_info = my_clusters)
  ##
  
  square <- function(x, label_size) { 
    ggplot()  + 
      coord_fixed(xlim=c(0,1), ylim = c(0,1)) + theme_void() + 
      theme(plot.background = element_rect(fill = x)) + 
      geom_text(aes(0.5,0.5),label = x , size = label_size)
  }
  
  plot <- list()
  each_code <- vector()
  len1 <- length(cluster_info$k)
  for (n1 in 1:len1){  
    len2 <- nrow(cluster_info$centres[[n1]])
    for (n2 in 1:len2){
      each_code[n2] <- cluster_info$centres[[n1]]$dmc[[n2]]$hex
    }
    t <-tibble(colours=each_code,squares=purrr::map(colours,~square(.x,24/length(colours))))
    
    n_col <- length(t$colours)
    rect_dat <- tibble(x1 = c(0:(n_col-1)), x2 = c(1:n_col), y1 = rep(0,n_col),
                       y2 =rep(1,n_col), colour = t$colours)
    plot[[n1]] <- rect_dat %>% ggplot() + coord_fixed() + 
      geom_rect(aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=colour), color="black") +
      geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=colour), size=24/n_col) + 
      scale_fill_manual(values = rect_dat$colour)+ theme_void() + theme(legend.position = "none") 
  }
  q = paste0("k = ",c(cluster_info$k))
  plot_grid(plotlist = plot,
            labels = paste0(q), ncol = 1)
}








make_pattern <- function(cluster_info, k, x_size, black_white = FALSE, background_colour)
{
  ## make_pattern function is to plots the cross-stitch pattern, with input information
  ## to use the k cluster, create size with x_size, with in colour or just black and 
  ## white and remove the background_colour.
  ##
  ## Input:
  ## - cluster_info - Output of the process_image. A list or tibble of 
  ##   information derived from the k_means that will be sufficient to be the 
  ##   input to any other function you write. 
  ##   This includes: 
  ##    - The original output of the kclust calls,
  ##    - The tidied clusters, their associated RGB values and their nearest
  ##      DMC thread colour information.
  ##  - k - The chosen cluster size based on second and third function.
  ##  - x_size - The (approximate) total number of possible stitches in the 
  ## horizontal direction
  ##  - black_white - (logical) Print the pattern in black and white (TRUE) or 
  ## colour (FALSE, default)
  ##  - background_colour - The colour of the background, which should not be 
  ## stitched in the pattern. (Default is to not have a colour) 
  ## Output:
  ##  - This function should be the only function to use 
  ## change_resolution(image_df, x_size). It should produce a cross-stitch 
  ## pattern that can be followed, complete with a legend that has thread 
  ## colour, and a guide grid.
  ##
  ## Example:
  ##  make_pattern(cluster_info = my_clusters, k=7, x_size = 50, 
  ##               black_white = FALSE, background_color = "B5200")
  ##
 
  k_col <- cluster_info[cluster_info$k==k,]
  im_df <- k_col$tidy_dat[[1]]
  # Subsamples the im_df to produce a lower resolution image.
  low_im <- change_resolution(im_df, x_size = x_size)
  
  with_dmc <- k_col$centres[[1]]
  
  DMC_num <- vector()
  len <- nrow(with_dmc)
  for (n in 1:len){
    DMC_num[n] <- with_dmc$dmc[[n]]$dmc
    }
  with_dmc <- with_dmc %>% mutate(DMC=DMC_num)
  
  
  # Here I assigned DMC value to each cluster
  low_im$DMC<-NA

  for (i in 1:nrow(low_im)){
    with_num <- low_im[i,]$cluster
    low_im[i, ]$DMC <- with_dmc[with_dmc$cluster == with_num,]$DMC
  }
  
  dmc_info <- cluster_info[cluster_info$k== k,]$centres[[1]]$dmc
  dmc1 <- sapply(dmc_info, "[[", 1)
  name <- sapply(dmc_info, "[[", 2)
  col <- sapply(dmc_info, "[[", 3)
  combine <- paste(name, "(", dmc1, ")")
  img_frame <- tibble(dmc1, name, col,combine)
  
  if(black_white == TRUE){
    draw <- low_im%>%ggplot(aes(x,y))+
      geom_point(aes(shape=factor(DMC)))+
      scale_colour_manual(name="Stitch",values=img_frame%>%select(dmc1,col)%>%deframe,
                          label=img_frame%>%select(dmc1,combine)%>%deframe)+
      scale_shape_manual(name="Stitch",values=c(0:10),labels=img_frame%>%
                           select(dmc1,combine) %>%deframe)+
      scale_y_reverse()+theme_void()
    
  }
  else{
    draw <- low_im%>%ggplot(aes(x,y))+
      geom_point(aes(col=factor(DMC),
                     shape=factor(DMC)))+
      scale_colour_manual(name="Stitch",values=img_frame%>%select(dmc1,col)%>%deframe,
                          label=img_frame%>%select(dmc1,combine)%>%deframe)+
      scale_shape_manual(name="Stitch",values=c(0:10),labels=img_frame%>%
                           select(dmc1,combine) %>%deframe)+
      scale_y_reverse()+theme_void()
  }
  #Here we can find that the color "White" is the background color, with dmc code
  # "B5200". 
  background_colour <- "B5200"
  
    img_frame_no_back <- img_frame[img_frame[[1]] != background_colour, ]
    low_im_no_back <- low_im[low_im$DMC != background_colour,]
    if(black_white == TRUE){
      draw <- low_im_no_back%>%ggplot(aes(x,y))+
        geom_point(aes(shape=factor(DMC)))+
        scale_colour_manual(name="Stitch",values=img_frame_no_back%>%select(dmc1,col)%>%deframe,
                            label=img_frame_no_back%>%select(dmc1,combine)%>%deframe)+
        scale_shape_manual(name="Stitch",values=c(0:10),labels=img_frame_no_back%>%
                             select(dmc1,combine) %>%deframe)+
        scale_y_reverse()+theme(
          panel.grid = element_line(size = 0.5, linetype = 'solid',
                                    colour = "black"))+theme_bw()

      
    }
    else{
      draw <- low_im_no_back%>%ggplot(aes(x,y))+
        geom_point(aes(col=factor(DMC),
                       shape=factor(DMC)))+
        scale_colour_manual(name="Stitch",values=img_frame_no_back%>%select(dmc1,col)%>%deframe,
                            label=img_frame_no_back%>%select(dmc1,combine)%>%deframe)+
        scale_shape_manual(name="Stitch",values=c(0:10),labels=img_frame_no_back%>%
                             select(dmc1,combine) %>%deframe) +theme(
          panel.grid = element_line(size = 0.5, linetype = 'solid',
                                          colour = "black"))+ scale_y_reverse()+theme_bw()
        
    }


  draw
}

