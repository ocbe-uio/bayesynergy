# remodel the code for plotting functions 


#' Plot function for a bayesynergy object
#'
#' @description A function for plotting synergy surfaces and summary statistics of a \code{bayesynergy} object.
#'
#' @param x An object of class \code{bayesynergy}, the result of \code{\link{bayesynergy}}.
#' @param plot3D logical; if TRUE, interactive 3D surface plots of dose response function are displayed.
#' @param save_plots logical; if TRUE plots are saved locally.
#' @param path string; path for saving plots, if NULL defaults to work directory.
#' @param plotdevice string; device for saving plots locally, must be 'pdf' or 'png'
#' @param ... further arguments passed on to device for plotting, useful for setting width and height of saved plots.
#'
#' @details
#' This function extends \code{plot} to draw response and interaction surfaces of the fitted model. Both three-dimensional interactive plots, and two-dimensional contour plots can be displayed.
#'
#'
#' @examples
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' plot(fit)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point ylim xlim labs scale_x_continuous stat_contour_filled geom_contour_filled scale_fill_manual scale_y_continuous scale_y_discrete xlab ylab scale_fill_viridis_c guides guide_legend geom_hline annotate coord_cartesian
#' @importFrom ggridges geom_density_ridges_gradient theme_ridges
#' @importFrom plotly plot_ly add_surface add_trace %>% plotly_build as_widget add_paths
#' @importFrom scales math_format
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off
#' @importFrom Cairo CairoPNG CairoPDF
#' @importFrom stats quantile
#' @importFrom htmlwidgets saveWidget
#'
#' @export


# data to plot with ----

# data("mathews_DLBCL")
# y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
# x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
# fit <- bayesynergy(y_mat,x_mat)
fit <- bayesynergy(y_mat, x_mat, method = 'vb') # this is faster
# ?bayesynergy::bayesynergy()
# ?bayesynergy




plot.bayesynergy_chiversion <- function(x, 
                                        plot3D = T, 
                                        save_plots = FALSE, 
                                        path = NULL, 
                                        plotdevice = "pdf", ...){
  
  # x <- fit
  
  # Check for valid plotdevice
  if (!(plotdevice %in% c("pdf","png"))){
    stop("plotdevice must be one of {'pdf','png'}")
  }
  # Setting default path if note given
  if (is.null(path)){
    path = getwd()
  }
  # Check that path exists, if not create it
  if (!dir.exists(path)){
    dir.create(path)
  }
  # Put dirmark on file path
  path <- Sys.glob(path,dirmark = T)
  
  
  # processing ----
  # Creating some stuff needed for plots
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  # We add lower-asymptote parameters if these are not estimated
  if (!x$model$lower_asymptotes){
    posterior$la_1 <- rep(0,n.save)
    posterior$la_2 <- rep(0,n.save)
  }
  
  unqX1 <- log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 <- log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  dx1 <- mean(diff(unqX1))
  dx2 <- mean(diff(unqX2))
  nrep <- ncol(as.matrix(x$data$y))
  
  # Need to find coordinates for the observed variables in this new coordinate system
  Xgrid <- expand.grid(unqX1,unqX2)
  Xgrid <- Xgrid[order(Xgrid[,"Var1"],Xgrid[,"Var2"]),]
  
  mono1 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  # Remove NA here
  mono1 <- mono1[!is.na(mono1$y),]
  
  mono2 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  # Remove NA here
  mono2 <- mono2[!is.na(mono2$y),]
  
  # Pull out indices we want
  ii <- x$data$indices[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0))]
  # Also define residuals here
  
  combination <- data.frame(
    x1 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),1]),nrep),
    x2 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0)),]),
    f = 0,
    p0 = 0,
    Delta = 0,
    residuals = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0)),]) - as.vector(x$posterior_mean$p0)[ii]
  )
  # Remove NA here
  combination <- combination[!is.na(combination$y),]
  
  
  
  # Monotherapies ----
  # setting grid
  
  grid.size <- 100
  x.seq1 <- seq(min(unqX1)-dx1,max(unqX1)+dx1,length.out = grid.size)
  x.seq2 <- seq(min(unqX2)-dx2,max(unqX2)+dx2,length.out = grid.size)
  
  y.seq1 <- matrix(NA,nrow=grid.size,ncol=n.save)
  y.seq2 <- matrix(NA,nrow=grid.size,ncol=n.save)
  for (i in 1:grid.size){
    y.seq1[i,] <- as.vector(posterior$la_1)+as.vector((1-posterior$la_1))/(1+10^(as.vector(posterior$slope_1)*(x.seq1[i]-as.vector(posterior$log10_ec50_1))))
    y.seq2[i,] <- as.vector(posterior$la_2)+as.vector((1-posterior$la_2))/(1+10^(as.vector(posterior$slope_2)*(x.seq2[i]-as.vector(posterior$log10_ec50_2))))
  }
  df1 <- data.frame(
    x = x.seq1,
    mean = apply(y.seq1,1,mean),
    median = apply(y.seq1,1,median),
    lower = apply(y.seq1,1,quantile, probs=0.025),
    upper = apply(y.seq1,1,quantile, probs=0.975)
    # clower = apply(y.seq1,1,mean)-1.96*sqrt(mean(posterior$s)^2*(apply(y.seq1,1,mean)+x$model$lambda)),
    # cupper = apply(y.seq1,1,mean)+1.96*sqrt(mean(posterior$s)^2*(apply(y.seq1,1,mean)+x$model$lambda))
  )
  df2 <- data.frame(
    x = x.seq2,
    mean = apply(y.seq2,1,mean),
    median = apply(y.seq2,1,median),
    lower = apply(y.seq2,1,quantile, probs=0.025),
    upper = apply(y.seq2,1,quantile, probs=0.975)
    # clower = apply(y.seq2,1,mean)-1.96*sqrt(mean(posterior$s)^2*(apply(y.seq2,1,mean)+x$model$lambda)),
    # cupper = apply(y.seq2,1,mean)+1.96*sqrt(mean(posterior$s)^2*(apply(y.seq2,1,mean)+x$model$lambda))
  )
  
  
  
  
  # 3d (plotly) (TAKE OUT) ----
  #  if(plot3D){
  #   
  #   # Response
  #   z_response = x$posterior_mean$f[-1,-1]
  #   fig = plot_ly(x = unqX1, y = unqX2, z = z_response)
  #   fig = fig %>% add_surface(cmin=0,cmax=1)
  #   fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
  #                           type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="black",symbol=104),name = "Observed")
  #   fig = fig %>% plotly::layout(scene = list(zaxis = list(range=c(min(min(0,c(mono1$y,mono2$y,combination$y))),max(max(1,c(mono1$y,mono2$y,combination$y)))),
  #                                                          title="Viability",titlefont = list(size = 12)),
  #                                             xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
  #                                             yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
  #                                title = paste("Response surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
  #   fig = fig %>% add_paths(x = df1$x, y = (min(unqX2)-mean(diff(unqX2))), z = df1$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F)
  #   fig = fig %>% add_paths(x = (min(unqX1)-mean(diff(unqX1))), y = df2$x, z = df2$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F)
  #   fig = fig %>% add_trace(x = mono1$x, y = (min(unqX2)-mean(diff(unqX2))), z = mono1$y, type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="grey",symbol=104), showlegend = F)
  #   fig = fig %>% add_trace(x = (min(unqX1)-mean(diff(unqX1))), y = mono2$x, z = mono2$y, type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="grey",symbol=104), showlegend = F)
  #   for (i in 1:length(unqX1)){
  #     fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_response[,i]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   for (i in 1:length(unqX2)){
  #     fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_response[i,]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   
  #   response_3d = fig
  #   
  #   
  #   # Noninteraction ----#
  #   z_p0 = x$posterior_mean$p0[-1,-1]
  #   fig = plot_ly(x = unqX1, y = unqX2, z = z_p0)
  #   fig = fig %>% add_surface(cmin=0,cmax=1)
  #   fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
  #                           type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="black",symbol=104), name = "Observed")
  #   fig = fig %>% plotly::layout(scene = list(zaxis = list(range=c(min(min(0,c(mono1$y,mono2$y,combination$y))),max(max(1,c(mono1$y,mono2$y,combination$y)))),
  #                                                          title="Viability",titlefont = list(size = 12)),
  #                                             xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
  #                                             yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
  #                                title = paste("Non-interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
  #   fig = fig %>% add_paths(x = df1$x, y = (min(unqX2)-mean(diff(unqX2))), z = df1$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F)
  #   fig = fig %>% add_paths(x = (min(unqX1)-mean(diff(unqX1))), y = df2$x, z = df2$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F)
  #   fig = fig %>% add_trace(x = mono1$x, y = (min(unqX2)-mean(diff(unqX2))), z = mono1$y, type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="grey",symbol=104), showlegend = F)
  #   fig = fig %>% add_trace(x = (min(unqX1)-mean(diff(unqX1))), y = mono2$x, z = mono2$y, type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="grey",symbol=104), showlegend = F)
  #   for (i in 1:length(unqX1)){
  #     fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_p0[,i]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   for (i in 1:length(unqX2)){
  #     fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_p0[i,]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   noninter_3d = fig
  #   
  #   # Interaction ----#
  #   z_Delta = x$posterior_mean$Delta[-1,-1]
  #   fig = plot_ly(type = "mesh3d")
  #   fig = fig %>% add_trace(x = unqX1, y = unqX2, z = z_Delta,type = "surface",
  #                           colorscale = list(c(0,0.5, 1), c("2166AC","EAECCC", "B2182B")),cmin=-1,cmax=1)
  #   fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$residuals,
  #                           type = "scatter3d", mode = "markers",
  #                           marker = list(size=3,color="black",symbol=104), name = "y - p<sub>0</sub>", showlegend = T)
  #   fig = fig %>% plotly::layout(scene = list(zaxis= list(range=c(min(-1,min(combination$residuals)),max(1,max(combination$residuals))),
  #                                                         title="Interaction",titlefont = list(size = 12)),
  #                                             xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
  #                                             yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
  #                                title = paste("Interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
  #   for (i in 1:length(unqX1)){
  #     fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_Delta[,i]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   for (i in 1:length(unqX2)){
  #     fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_Delta[i,]+0.003, type="scatter3d", mode="lines",
  #                             showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
  #   }
  #   inter_3d = fig
  #   
  # }
  # 
  # 
  
  
  
  
  
  # p1: df1
  # p2: df2
  
  if (x$model$robust){
    # robust ----
    ## p1: drug1 ----

    p1 <- ggplot(data = df1, mapping = aes(x, median)) +
      geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
      geom_point(data = mono1, mapping = aes(x,y), na.rm = T) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(min(-0.5,min(mono1$y)),max(1.5,max(mono1$y))) +
      labs(x = bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")"), y = "Viability",
           title = paste("Monotherapy:",x$data$drug_names[1]),
           subtitle = paste0(x$data$experiment_ID)) +
      scale_x_continuous(labels = math_format(10^.x))
    
    ## p2: drug2 ----
    p2 <- ggplot(data = df2, mapping = aes(x, median)) +
      geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
      geom_point(data = mono2, mapping = aes(x,y), na.rm = T) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(min(-0.5,min(mono2$y)),max(1.5,max(mono2$y))) +
      labs(x = bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")"), y = "Viability",
           title = paste("Monotherapy:",x$data$drug_names[2]),
           subtitle = paste0(x$data$experiment_ID)) +
      scale_x_continuous(labels = math_format(10^.x))
    
    
  }else{
    ## not robust ----
    # p1: drug 1----
    p1 <- ggplot(data = df1, mapping = aes(x, mean)) +
      geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
      geom_point(data = mono1, mapping = aes(x,y), na.rm = T) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(min(-0.5,min(mono1$y)),max(1.5,max(mono1$y))) +
      labs(x = bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")"), y = "Viability",
           title = paste("Monotherapy:",x$data$drug_names[1]),
           subtitle = paste0(x$data$experiment_ID)) +
      scale_x_continuous(labels = math_format(10^.x))
    
    # p2: drug 2 ----
    p2 <- ggplot(data = df2, mapping = aes(x, mean)) +
      geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
      geom_point(data = mono2, mapping = aes(x,y), na.rm = T) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ylim(min(-0.5,min(mono2$y)),max(1.5,max(mono2$y))) +
      labs(x = bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")"), y = "Viability",
           title = paste("Monotherapy:",x$data$drug_names[2]),
           subtitle = paste0(x$data$experiment_ID)) +
      scale_x_continuous(labels = math_format(10^.x))
  
    
    
    }
  
  
  
  
  # Combinations ----
  # this df is for p3, p4
  # is Xgrid used in p1 p2?
  df <- data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    f = as.vector(x$posterior_mean$f[-1,-1]),
    p0 = as.vector(x$posterior_mean$p0[-1,-1]),
    Delta = as.vector(x$posterior_mean$Delta[-1,-1])
  )
  
  # Breaks
  breaks <- seq(0,1,by=0.1)
  
  # Function to create labels for legend
  breaklabel <- function(x, breaks){
    labels = paste0("[", 100*breaks[1],",",100*breaks[2],")")
    labels = c(labels,paste0("(",100*breaks[2:10], ",", 100*breaks[3:11],"]"))
    labels[1:x]
  }
  
  # p3: response surf ----
  p3 <- ggplot(data = df, aes(x = x1, y = x2, z = f)) +
    geom_contour_filled(breaks = breaks, show.legend = T,
                        colour="dark grey", linewidth=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(values=viridis::viridis(10), labels=breaklabel(10,breaks),
                      name="% Viability", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 1, x = -Inf, y = mono2$x, color = 'black') +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 1, x = mono1$x, y = -Inf, color = 'black') +
    coord_cartesian(clip = 'off') +
    labs(title = "Response surface",
         subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2]))
  

  # p4: non-interaction -----
  p4 <- ggplot(data = df, aes(x = x1, y = x2, z = p0)) +
    geom_contour_filled(breaks = breaks, show.legend = T,
                        colour="dark grey", linewidth=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(values=viridis::viridis(10), labels = breaklabel(10,breaks),
                      name="% Viability", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, x = -Inf, y = mono2$x, color = 'black') +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, x = mono1$x, y = -Inf, color = 'black') +
    coord_cartesian(clip = 'off') +
    labs(title = "Non-interaction surface",
         subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2]))
  
  
  # p5: interaction ----
  #Interaction surface has a different color scale
  # Delta_col_palette <- inlmisc::GetColors(scheme = "sunset")
  Delta_col_palette <- c("#354B99", "#4062A8", "#4A7AB7", "#5D90C2", "#6DA5CC", "#83B8D7", "#98CAE0",
                         "#ACD6E7", "#C2E3EE", "#D7E7DD", "#EAEBCC", "#F5E3AB", "#FED98B", "#FEC678",
                         "#FDB366", "#FA9958", "#F57D4A", "#EA603B", "#DC3D2D", "#C0242A", "#A50026")
  # Breaks
  eps <- 0.05
  breaks <- c(seq(-1,-0.1,by=0.1),0-eps,0+eps,seq(0.1,1,by=0.1))
  
  #Function to create labels for legend
  breaklabel <- function(x, breaks){
    labels = paste0("[",100*breaks[1],",",breaks[2],"]")
    labels = c(labels,paste0("(",100*breaks[2:10], ",", 100*breaks[3:11],"]"))
    labels = c(labels,paste0("(",100*breaks[11:21],",",100*breaks[12:22],"]"))
    labels[1:x]
  }
  

  p5 <- ggplot(data = df, aes(x = x1, y = x2, z = Delta)) +
    geom_contour_filled(breaks = breaks, show.legend = T,
                        colour="dark grey", size=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(values=Delta_col_palette, labels = breaklabel(21,breaks),
                      name="% Interaction", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, x = -Inf, y = mono2$x, color = 'black') +
    annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, x = mono1$x, y = -Inf, color = 'black') +
    coord_cartesian(clip = 'off') +
    labs(title = "Interaction surface",
         subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2])) +
    guides(fill=guide_legend(ncol=1))
  

  # Summary statistics

  # p6: dss ----
  df <- data.frame(
    dss = c(posterior$dss_1,posterior$dss_2),
    idx = factor(c(rep(x$data$drug_names[1],n.save),rep(x$data$drug_names[2],n.save)),
                 levels = c(x$data$drug_names[1],x$data$drug_names[2]))
  )
  p6 <- ggplot(df, aes(x = dss, y = idx, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, gradient_lwd = 1., from = 0, to = 100) +
    scale_fill_viridis_c(name = "DSS", option = "C",limits=c(0,100)) +
    labs(title = 'Estimated drug sensitivity scores',
         subtitle = x$data$experiment_ID,
         y = "") +
    xlim(0,100) +
    scale_y_discrete(limits = unique(rev(df$idx))) +
    theme_ridges(font_size = 13, grid = FALSE)
  
  
  # p7: dcs ----
  df <- data.frame(
    rVUS = c(posterior$rVUS_f, posterior$rVUS_p0,
             abs(posterior$VUS_syn), posterior$VUS_ant),
    idx = factor(c(rep('rVUS_f',n.save),rep('rVUS_p0',n.save),
                   rep('VUS_syn',n.save),rep('VUS_ant',n.save)),
                 levels = c('rVUS_f','rVUS_p0','VUS_syn','VUS_ant'))
  )
  
  p7 <- ggplot(df, aes(x = rVUS, y = idx, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 2, gradient_lwd = 1., from = 0, to = 100) +
    scale_fill_viridis_c(name = "", option = "C",limits=c(0,100)) +
    labs(title = 'Estimated drug combination scores',
         subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2]),
         y = "") +
    xlim(0,100) +
    xlab("") +
    scale_y_discrete(limits = unique(rev(df$idx)),
                     labels=c("Antagonism\nVUS(\u0394\u207a)",
                              "Synergy\n|VUS(\u0394\u207b)|",
                              "Bliss\nefficacy\nrVUS(p\u2080)",
                              "Overall\nefficacy\nrVUS(f)")
                     #expression(paste("Overall efficacy (",rVUS(f),")"))),
                     #expand = expansion(add = c(0.01))
    ) +
    # labels=c("Antagonism","Synergy","Overall Interaction","Overall Efficacy")) +
    theme_ridges(font_size = 12, grid = FALSE)
  
  
  
  # ________ -----
  # save functions -------
  if (!save_plots){
    # Displaying the plots
    grid.arrange(p1,p2,nrow=2) # Monotherapies
    readline("Press key for next plot")
    print(p6) # DSS scores
    readline("Press key for next plot")
    if (plot3D){
      # fig <- subplot(response_3d, noninter_3d,inter_3d)
      # plotly_build(fig)
      print(response_3d) # Response 3d
      readline("Press key for next plot")
      print(noninter_3d) # Noninteraction 3d
      readline("Press key for next plot")
      print(inter_3d) # interaction 3d
      readline("Press key for next plot")
    }
    print(p3) # Response
    readline("Press key for next plot")
    print(p4) # p0
    readline("Press key for next plot")
    print(p5) # Interaction
    readline("Press key for next plot")
    print(p7) # rVUS
  } else { # Save plots locally
    # Monotherapies
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Monotherapies",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    grid.arrange(p1,p2,nrow=2) # Monotherapies
    dev.off()
    
    # 3D plots ----
    if (plot3D){
      # Response
      file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Response3D.html",sep = "_"))
      htmlwidgets::saveWidget(as_widget(response_3d),file = file.name, libdir = "lib", selfcontained = F)
      # Non-interaction
      file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"NonInteraction3D.html",sep = "_"))
      htmlwidgets::saveWidget(as_widget(noninter_3d),file = file.name, libdir = "lib", selfcontained = F)
      # Interaction
      file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Interaction3D.html",sep = "_"))
      htmlwidgets::saveWidget(as_widget(inter_3d),file = file.name, libdir = "lib", selfcontained = F)
    }
    # DSS scores
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"DSS",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    print(p6) # DSS scores
    dev.off()
    # Response surface
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Response",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    print(p3) # Response surface
    dev.off()
    # Non-interaction
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"NonInteraction",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    print(p4) # Noninteraction surface
    dev.off()
    # Interaction
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Interaction",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    print(p5) # Interaction surface
    dev.off()
    # rVUS
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"rVUS",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      CairoPDF(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      CairoPNG(filename = file.name.png, ...)
    }
    print(p7) # rVUS scores
    dev.off()
  }
}






# _____________ ----
# _____________ ----

# multiple ----

#' data("mathews_DLBCL")
experiment1 <- list(y = mathews_DLBCL[[1]][[1]],
                    x = mathews_DLBCL[[1]][[2]],
                    drug_names = c("ispinesib","ibrutinib"))
experiment2 <- list(y = mathews_DLBCL[[2]][[1]],
                  x = mathews_DLBCL[[2]][[2]],
                   drug_names = c("canertinib","ibrutinib"))
experiments <- list(experiment1,experiment2)

# fit <- synergyscreen(experiments)
# plot(fit)

x <- fit

plot.synergyscreen_chiversion <- function(x, 
                                          groupbyExperimentID = T, 
                                          save_plots = FALSE, 
                                          path = NULL, 
                                          plotdevice = "pdf", ...){
  
  # Check for valid plotdevice
  if (!(plotdevice %in% c("pdf","png"))){
    stop("plotdevice must be one of {'pdf','png'}")
  }
  # Setting default path if note given
  if (is.null(path)){
    path = getwd()
  }
  # Check that path exists, if not create it
  if (!dir.exists(path)){
    dir.create(path)
  }
  # Put dirmark on file path
  path = Sys.glob(path,dirmark = T)
  
  
  # processing ----
  # Extracting the dataframe needed from the output
  synscores <- x$screenSummary
  
  
  # Do some preprocessing to make sure we are able to average across the pairs correctly
  # Particularly, if (drugA + drugB) and (drugB + drugA) are equal for two experimental IDs,
  # they should be averaged, as they represent the same pair
  # We do this alphabetically such that drug A comes earlier in the alphabet than drug B
  
  # add 2 col, replicates drug A and drug B columns
  synscores <- synscores %>% mutate(`drugATMP` = `Drug A`, `drugBTMP` = `Drug B`)
  synscores <- synscores %>% mutate(`Drug A` = ifelse(`Drug A` < `Drug B`,`Drug A`,`Drug B`))
  synscores <- synscores %>% mutate(`Drug B` = ifelse(`Drug A` < `Drug B`,`Drug B`,`drugATMP`))
  synscores <- select(synscores,-c(`drugATMP`,`drugBTMP`))
  synscores <- synscores %>% mutate(`Drug Pair` = paste(`Drug A`,"+",`Drug B`))
  
  # Create a list per experiment ID
  synscoreList <- synscores %>% dplyr::group_by(`Experiment ID`,.add=T) %>% group_split()
  
  # Average across all experimentID
  synscoresAverage <- synscores %>% 
    dplyr::group_by(`Drug A`,`Drug B`) %>% 
    dplyr::summarise(`mean_syn` = mean(`Synergy Score`,na.rm=T),
                     `max_syn` = min(`Synergy (mean)`,na.rm = T),
                     `mad_syn` = mad(`Synergy (mean)`,na.rm=T),
                     `mean_ant` = mean(`Antagonism Score`,na.rm=T),
                     `max_ant` = max(`Antagonism (mean)`,na.rm=T),
                     `mad_ant` = mad(`Antagonism (mean)`,na.rm=T),
                     `mean_int` = mean(`Interaction (mean)`,na.rm=T),
                     `abs_mean_int` = abs(`mean_int`),
                     `mad_int` = mad(`Interaction (mean)`,na.rm=T),
                     `n` = n(), .groups = "drop")
  
  
  # For p3, we also create the mirror-image so that we can display both antagonism and synergy
  synscores3AB <- synscoresAverage %>% mutate(`mean_int` = `mean_syn`)
  synscores3BA <- synscoresAverage %>% mutate(`DrugATMP` = `Drug A`,
                                              `Drug A` = `Drug B`,
                                              `Drug B` = `DrugATMP`,
                                              `mean_int` = `mean_ant`) %>% select(-c(`DrugATMP`))
  synscores3 <- rbind(synscores3AB,synscores3BA)
  
  # Breaks set by quantiles
  lower_breaks <- quantile(synscores3$mean_syn,probs=c(0,0.01,0.1,0.25))
  upper_breaks <- sort(quantile(synscores3$mean_ant,probs=1-c(0,0.01,0.1,0.25)))
  l <- max(abs(c(min(lower_breaks)-0.1,max(upper_breaks)+0.1)))
  lower_breaks[1] <- -l
  upper_breaks[4] <- l
  breaks <- as.numeric(c(lower_breaks,upper_breaks))
  breaks <- sort(jitter(breaks))
  # Delta_col_palette <- inlmisc::GetColors(scheme = "BuRd")
  Delta_col_palette <- c("#2066AC", "#60A3CC", "#BCDAEA", "#F7F7F7", "#FBC9AF", "#E0775E", "#B2182B")
  
  
  #Switching back some of Drug A and Drug B, as it makes p3 nicer
  # synscoresAverage <- synscoresAverage %>% mutate(`DrugATMP` = `Drug A`)
  # switchIDX = sample(1:nrow(synscoresAverage),floor(nrow(synscoresAverage)/2))
  # synscoresAverage$`Drug A`[switchIDX] = synscoresAverage$`Drug B`[switchIDX]
  # synscoresAverage$`Drug B`[switchIDX] = synscoresAverage$DrugATMP[switchIDX]
  # synscoresAverage <- select(synscoresAverage,-c(`DrugATMP`))
  
  
  # We need to give all experiments an "sd" value for the size in p3a, even the ones with only one observation
  synscoresAverage$mad_syn[is.na(synscoresAverage$mad_syn)] <- 0
  synscoresAverage$mad_ant[is.na(synscoresAverage$mad_ant)] <- 0
  synscoresAverage$mad_int[is.na(synscoresAverage$mad_int)] <- 0
  
  
  # Creating a function to label top synergy scores
  in_quantile <- function(x, n = 5) {
    # Pick out to get top and bottom 5
    vec = rep(F,length(x))
    vec[tail(order(x),n)] = T
    return(vec)
  }
  
  # Labeling top hits
  synscoresAverage <- synscoresAverage %>% 
    mutate(`Quantile` = ifelse(in_quantile(as.numeric(-`mean_syn`))|in_quantile(as.numeric(`mean_ant`)),
                               paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
    # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
    mutate(`Quantile` = ifelse(!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
  
  
  # Duplicate experimentID and drug pairs are also removed (there shouldn't be any)
  synscoresAverage <- synscoresAverage[!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))),]
  # Adding a colour variable for p1
  synscoresAverage <- synscoresAverage %>% mutate(`colp1` = ifelse(`mean_ant` > -`mean_syn`,`max_ant`,`max_syn`))
  
  
  
  # p1: synergy v antagonism ----
  # PLOT 1: Synergy vs. Antagonism across all experiments
  
  
  p1 <- ggplot(synscoresAverage, aes(x = `mean_ant`, y = `mean_syn`,fill=`colp1`)) +
    geom_point(color = "gray", shape = 21, aes(size = `mad_int`)) +
    # geom_abline(intercept=0,slope=1, linetype = "dashed", alpha=0.3) +
    scale_size_continuous(range = c(1.5, 4), name = "MAD") +
    scale_fill_gradientn(colours = c("#2066AC","#F7F7F7","#B2182B"),
                         values = scales::rescale(c(min(synscoresAverage$colp1,na.rm=T),max(synscoresAverage$colp1,na.rm=T))),
                         name = "Max\nInteraction",
                         limits = max(abs(synscoresAverage$colp1))*c(-1,1)) +
    # limits = c(min(0,min(synscoresAverage$colp1,na.rm=T)),max(0,max(synscoresAverage$colp1,na.rm=T)))) +
    geom_text_repel(aes(label = Quantile), color="black",  na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"),
                    segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
    coord_cartesian(xlim = c(), clip = "off") +
    
    xlab("Antagonism") +
    ylab("Synergy") +
    ggtitle("Average interaction across all experiments") +
    guides() +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
          legend.justification="center", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, vjust = 1, hjust = 1))
  
  
  # p2: synergy per drug -----
  #PLOT 2: Synergy per drug across all experiments

  
  # Label top hits for each experiment ID
  synscores <- synscores %>% dplyr::group_by(`Experiment ID`) %>% 
    mutate(`Quantile` = ifelse(in_quantile(-as.numeric(`Synergy (mean)`),n=1), paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA)))
  
  
  # Plotting only synergy for each experiment ID
  p2 <- ggplot(synscores, aes(x = factor(`Experiment ID`), y = `Synergy (mean)`, fill = `Synergy (mean)`)) +
    geom_point(aes(size = `Synergy (sd)`), colour="gray",pch = 21, alpha = 0.8) +
    scale_size_continuous(trans = 'reverse', range = c(1, 4), name = "Std. dev") +
    scale_fill_gradientn(colours = c("#2066AC","#F7F7F7"),
                         values = scales::rescale(c(min(synscores$`Synergy (mean)`,na.rm = T),max(synscores$`Synergy (mean)`,na.rm=T))),
                         limits = c(min(synscores$`Synergy (mean)`,na.rm = T),max(synscores$`Synergy (mean)`,na.rm=T))) +
    scale_y_continuous(name = "Synergy") +
    geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"),
                    segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
    xlab("") +
    ylab("Synergy") +
    ggtitle("Synergy per experiment ID") +
    guides()  +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),
          legend.position = "bottom", legend.title = element_text(size=10), legend.title.align = 0.5, legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(.1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
  
  
  # p3: synergy antagonism pairwise -----
  # PLOT 3: Synergy & Antagonism pairwise

  
  # Function for labels in p3
  label_Function <- function(x) {
    tmp <- c()
    # print(x)
    if (length(x) == 2){
      tmp = c("","Top 1%\nantagonistic")
    }
    else if (length(x) == 6){
      tmp = c("Top 1%\nsynergistic","Top 10%\nsynergistic","Top 25%\nsynergistic","\nNo effect","Top 25%\nantagonistic","Top 10%\nantagonistic")
    }
    return(tmp)
  }
  
  p3 <- ggplot(data = synscores3, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`mean_int`")) +
    geom_point(color = "gray", shape = 21, aes_string(size = "`mad_int`")) +
    scale_size_continuous(range = c(2, 5), name = "MAD") +
    scale_fill_stepsn(colors = Delta_col_palette,
                      breaks = breaks,
                      values = scales::rescale(breaks),
                      limits = c(breaks[1]*0.99,breaks[8]*0.99),
                      labels = label_Function,
                      guide = guide_coloursteps(even.steps = T,
                                                show.limits = T,
                                                title = NULL,
                                                barheight = unit(3.3, "in"),
                                                label.vjust = 1.5,order=1)) +
    ggtitle("Average pairwise drug interaction across all experiments") +
    coord_cartesian(clip = 'off') +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
          legend.justification="center", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, vjust = 0.3, hjust = 1))
  
  if (length(synscoreList) > 1){
    if (!save_plots){
      print(p1)
      readline("Press key for next plot")
      print(p2)
      readline("Press key for next plot")
      print(p3)
      readline("Press key for next plot")
    } else { # Save plots locally
      # First plot
      file.name = paste0(path,paste("Average_interaction",sep = "_"))
      if (plotdevice == "pdf"){
        file.name.pdf = paste0(file.name,".pdf")
        cairo_pdf(file=file.name.pdf, ...)
      } else {
        file.name.png = paste0(file.name,".png")
        png(filename = file.name.png, ...)
      }
      print(p1)
      dev.off()
      # Second plot
      file.name = paste0(path,paste("Average_synergy",sep = "_"))
      if (plotdevice == "pdf"){
        file.name.pdf = paste0(file.name,".pdf")
        cairo_pdf(file=file.name.pdf, ...)
      } else {
        file.name.png = paste0(file.name,".png")
        png(filename = file.name.png, ...)
      }
      print(p2)
      dev.off()
      # Third plot (a)
      file.name = paste0(path,paste("Average_pairwise_interaction",sep = "_"))
      if (plotdevice == "pdf"){
        file.name.pdf = paste0(file.name,".pdf")
        cairo_pdf(file=file.name.pdf, ...)
      } else {
        file.name.png = paste0(file.name,".png")
        png(filename = file.name.png, ...)
      }
      print(p3)
      dev.off()
    }
    
  }

  
  
    
  # _____ ----
  # PLOT 4-6: Synergy vs. Antagonism per experiment ID

  if (groupbyExperimentID){
    # this is specified in the plotting function
    # true by default
    
    # if (length(synscoreList) > 1){
    for (i in 1:length(synscoreList)){
      # i <- 1
      synSlice <- synscoreList[[i]]
      
      # Labeling top hits
      synSlice <- synSlice %>% mutate(`Quantile` = ifelse(in_quantile(-as.numeric(`Synergy (mean)`))|in_quantile(as.numeric(`Antagonism (mean)`)),
                                                          paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
        # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
        mutate(`Quantile` = ifelse(!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
      # Duplicate drug pairs are also removed.
      
      synSlice <- synSlice[!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))),]
      
      # Adding a colour variable for p1
      synSlice <- synSlice %>% mutate(`colp1` = ifelse(`Antagonism (mean)` > -`Synergy (mean)`,`Antagonism (mean)`, `Synergy (mean)`))
      
      
      # For p6, we also create the mirror-image so that we can display both antagonism and synergy
      synscores3AB <- synSlice %>% mutate(`mean_int` = `Synergy (mean)`)
      synscores3BA <- synSlice %>% mutate(`DrugATMP` = `Drug A`,
                                          `Drug A` = `Drug B`,
                                          `Drug B` = `DrugATMP`,
                                          `mean_int` = `Antagonism (mean)`) %>% select(-c(`DrugATMP`))
      synscores3 = rbind(synscores3AB,synscores3BA)
      
      # synscores3 <- synscores3 %>% mutate(`mean_int` = `mean_int`-mean(`mean_int`))
      
      # Breaks set by quantiles
      lower_breaks = quantile(synscores3$`Synergy (mean)`,probs=c(0,0.01,0.1,0.25))
      upper_breaks = sort(quantile(synscores3$`Antagonism (mean)`,probs=1-c(0,0.01,0.1,0.25)))
      l = max(abs(c(min(lower_breaks)-0.1,max(upper_breaks)+0.1)))
      lower_breaks[1] = -l
      upper_breaks[4] = l
      breaks = as.numeric(c(lower_breaks,upper_breaks))
      breaks = sort(jitter(breaks))
      
      
      # Center
      
      
      #Switching back some of Drug A and Drug B, as it makes p6 nicer
      # synSlice <- synSlice %>% mutate(`DrugATMP` = `Drug A`)
      # switchIDX = sample(1:nrow(synSlice),floor(nrow(synSlice)/2))
      # synSlice$`Drug A`[switchIDX] = synSlice$`Drug B`[switchIDX]
      # synSlice$`Drug B`[switchIDX] = synSlice$DrugATMP[switchIDX]
      # synSlice <- select(synSlice,-c(`DrugATMP`))
      
      ###########################################################################################
      ####### PLOT 4: Synergy vs. Antagonism per experiment ID
      ###########################################################################################
      
      p4 = ggplot(synSlice, aes(x = `Antagonism (mean)`, y = `Synergy (mean)`,fill=`colp1`)) +
        geom_point(color = "gray", shape = 21, aes(size = `Interaction (sd)`)) +
        scale_size_continuous(trans = "reverse",range = c(1.5, 4), name = "Std. dev.") +
        scale_fill_gradientn(colours = c("#2066AC","#F7F7F7","#B2182B"),
                             values = scales::rescale(c(min(synSlice$colp1,na.rm=T),max(synSlice$colp1,na.rm=T))),
                             name = "Interaction",
                             limits = max(abs(synSlice$colp1))*c(-1,1)) +
        geom_text_repel(aes(label = Quantile), color="black",  na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"),
                        segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
        coord_cartesian(xlim = c(), clip = "off") +
        facet_grid(`Experiment ID`~.) +
        xlab("Antagonism") +
        ylab("Synergy") +
        ggtitle(paste("Synergy vs. Antagonism: ",synSlice$`Experiment ID`[1])) +
        guides() +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
              legend.justification="center", legend.box.margin=margin(1,0,0,0),
              text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
              plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
              axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
              axis.text.y = element_text(size = 8, vjust = 1, hjust = 1))
      
      ###########################################################################################
      ####### PLOT 5: Synergy per drug per experiment ID
      ###########################################################################################
      
      
      
      ###########################################################################################
      ####### PLOT 6: Synergy & Antagonism pairwise per Experiment ID
      ###########################################################################################
      
      
      p6 = ggplot(data = synscores3, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`mean_int`")) +
        geom_point(color = "gray", shape = 21, aes_string(size = "`Interaction (sd)`")) +
        scale_size(trans = "reverse",range = c(1, 4), name = "Std. dev.") +
        scale_fill_stepsn(colors = Delta_col_palette,
                          breaks = breaks,
                          values = scales::rescale(breaks),
                          limits = c(breaks[1]*0.9999,breaks[8]*0.9999),
                          labels = label_Function,
                          guide = guide_coloursteps(even.steps = T,
                                                    show.limits = T,
                                                    title = NULL,
                                                    barheight = unit(3.3, "in"),
                                                    label.vjust = 1.5,
                                                    order=1)) +
        theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
        ggtitle(paste("Pairwise drug interaction: ",synSlice$`Experiment ID`[1])) +
        coord_cartesian(clip = 'off') +
        facet_grid(`Experiment ID`~.) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
              legend.justification="center", legend.box.margin=margin(1,0,0,0),
              text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
              plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
              axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
              axis.text.y = element_text(size = 8, vjust = 0.3, hjust = 1))
      
      
      
      if (!save_plots){
        print(p4)
        readline("Press key for next plot")
        print(p6)
        readline("Press key for next plot")
      } else { # Save plots locally
        # First plot
        file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Interaction",sep = "_"))
        if (plotdevice == "pdf"){
          file.name.pdf = paste0(file.name,".pdf")
          cairo_pdf(file=file.name.pdf, ...)
        } else {
          file.name.png = paste0(file.name,".png")
          png(filename = file.name.png, ...)
        }
        print(p4)
        dev.off()
        # Third plot (b)
        file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Pairwise_interaction",sep = "_"))
        if (plotdevice == "pdf"){
          file.name.pdf = paste0(file.name,".pdf")
          cairo_pdf(file=file.name.pdf, ...)
        } else {
          file.name.png = paste0(file.name,".png")
          png(filename = file.name.png, ...)
        }
        print(p6)
        dev.off()
      }
    }
    # }
    
  }
  
  
  
  
  
}













