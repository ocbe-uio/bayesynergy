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
#' @importFrom inlmisc GetColors
#' @importFrom grDevices dev.off
#' @importFrom Cairo CairoPNG CairoPDF
#' @importFrom stats quantile
#' @importFrom htmlwidgets saveWidget 
#' 
#' @export 


plot.bayesynergy <- function(x, plot3D = T, save_plots = FALSE, path = NULL, plotdevice = "pdf", ...){
  
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

  # Creating some stuff needed for plots
  posterior = rstan::extract(x$stanfit)
  n.save = length(posterior$lp__)
  # We add lower-asymptote parameters if these are not estimated
  if (!x$model$lower_asymptotes){
    posterior$la_1 = rep(0,n.save)
    posterior$la_2 = rep(0,n.save)
  }
  
  unqX1 = log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 = log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  dx1 = mean(diff(unqX1))
  dx2 = mean(diff(unqX2))
  nrep = ncol(as.matrix(x$data$y))
  # Need to find coordinates for the observed variables in this new coordinate system
  Xgrid = expand.grid(unqX1,unqX2)
  Xgrid = Xgrid[order(Xgrid["Var1"],Xgrid["Var2"]),]
  
  mono1 = data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  # Remove NA here
  idx = !is.na(mono1$y)
  mono1 = mono1[idx,]
  
  mono2 = data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  # Remove NA here
  idx = !is.na(mono2$y)
  mono2 = mono2[idx,]
  
  # Pull out indices we want
  ii = x$data$indices[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0))]
  # Also define residuals here
  
  combination = data.frame(
    x1 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),1]),nrep),
    x2 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0)),]),
    f = 0,
    p0 = 0,
    Delta = 0,
    residuals = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0)),]) - as.vector(x$posterior_mean$p0)[ii]
  ) 
  # Remove NA here
  idx = !is.na(combination$y)
  combination = combination[idx,]
  
  ####################################################################################
  # Monotherapies
  ####################################################################################
  grid.size = 100
  x.seq1 = seq(min(unqX1)-dx1,max(unqX1)+dx1,length.out = grid.size)
  x.seq2 = seq(min(unqX2)-dx2,max(unqX2)+dx2,length.out = grid.size)
  
  y.seq1 = matrix(NA,nrow=grid.size,ncol=n.save)
  y.seq2 = matrix(NA,nrow=grid.size,ncol=n.save)
  for (i in 1:grid.size){
    y.seq1[i,] = as.vector(posterior$la_1)+as.vector((1-posterior$la_1))/(1+10^(as.vector(posterior$slope_1)*(x.seq1[i]-as.vector(posterior$log10_ec50_1))))
    y.seq2[i,] = as.vector(posterior$la_2)+as.vector((1-posterior$la_2))/(1+10^(as.vector(posterior$slope_2)*(x.seq2[i]-as.vector(posterior$log10_ec50_2))))
  }
  df1 = data.frame(
    x = x.seq1,
    mean = apply(y.seq1,1,mean),
    median = apply(y.seq1,1,median),
    lower = apply(y.seq1,1,quantile, probs=0.025),
    upper = apply(y.seq1,1,quantile, probs=0.975)
  )
  df2 = data.frame(
    x = x.seq2,
    mean = apply(y.seq2,1,mean),
    median = apply(y.seq2,1,median),
    lower = apply(y.seq2,1,quantile, probs=0.025),
    upper = apply(y.seq2,1,quantile, probs=0.975)
  )
  
  # #################################
  # ## We include a 3d-plot option ##
  # #################################
  if(plot3D){
    
    # Response
    z_response = x$posterior_mean$f[-1,-1]
    fig = plot_ly(x = unqX1, y = unqX2, z = z_response)
    fig = fig %>% add_surface(cmin=0,cmax=1)
    fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
                            type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="black",symbol=104),name = "Observed")
    fig = fig %>% plotly::layout(scene = list(zaxis = list(range=c(min(min(0,c(mono1$y,mono2$y,combination$y))),max(max(1,c(mono1$y,mono2$y,combination$y)))),
                                                   title="% Viability",titlefont = list(size = 12)),
                                              xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
                                              yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
                         title = paste("Response surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    fig = fig %>% add_paths(x = df1$x, y = (min(unqX2)-mean(diff(unqX2))), z = df1$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F) 
    fig = fig %>% add_paths(x = (min(unqX1)-mean(diff(unqX1))), y = df2$x, z = df2$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F) 
    fig = fig %>% add_trace(x = mono1$x, y = (min(unqX2)-mean(diff(unqX2))), z = mono1$y, type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="grey",symbol=104), showlegend = F)
    fig = fig %>% add_trace(x = (min(unqX1)-mean(diff(unqX1))), y = mono2$x, z = mono2$y, type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="grey",symbol=104), showlegend = F)
    for (i in 1:length(unqX1)){
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_response[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_response[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    
    response_3d = fig
    
    
    # Noninteraction
    z_p0 = x$posterior_mean$p0[-1,-1]
    fig = plot_ly(x = unqX1, y = unqX2, z = z_p0)
    fig = fig %>% add_surface(cmin=0,cmax=1)
    fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
                            type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="black",symbol=104), name = "Observed")
    fig = fig %>% plotly::layout(scene = list(zaxis = list(range=c(min(min(0,c(mono1$y,mono2$y,combination$y))),max(max(1,c(mono1$y,mono2$y,combination$y)))),
                                                  title="% Viability",titlefont = list(size = 12)),
                                              xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
                                              yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
                         title = paste("Non-interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    fig = fig %>% add_paths(x = df1$x, y = (min(unqX2)-mean(diff(unqX2))), z = df1$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F) 
    fig = fig %>% add_paths(x = (min(unqX1)-mean(diff(unqX1))), y = df2$x, z = df2$mean, line = list(color = "grey", dash = "dash",width=4), showlegend = F) 
    fig = fig %>% add_trace(x = mono1$x, y = (min(unqX2)-mean(diff(unqX2))), z = mono1$y, type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="grey",symbol=104), showlegend = F)
    fig = fig %>% add_trace(x = (min(unqX1)-mean(diff(unqX1))), y = mono2$x, z = mono2$y, type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="grey",symbol=104), showlegend = F)
    for (i in 1:length(unqX1)){
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_p0[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_p0[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    noninter_3d = fig
    
    # Interaction
    z_Delta = x$posterior_mean$Delta[-1,-1]
    fig = plot_ly(type = "mesh3d")
    fig = fig %>% add_trace(x = unqX1, y = unqX2, z = z_Delta,type = "surface",
                            colorscale = list(c(0,0.5, 1), c("2166AC","EAECCC", "B2182B")),cmin=-1,cmax=1)
    fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$residuals,
                            type = "scatter3d", mode = "markers",
                            marker = list(size=3,color="black",symbol=104), name = "y - p<sub>0</sub>", showlegend = T)
    fig = fig %>% plotly::layout(scene = list(zaxis= list(range=c(-1,1),
                                                  title="Interaction",titlefont = list(size = 12)),
                                              xaxis = list(title=paste(x$data$units[1],x$data$drug_names[1]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>"),
                                              yaxis = list(title=paste(x$data$units[2],x$data$drug_names[2]),titlefont = list(size = 12),tickprefix="10<sup>",tickfont=list(size=10),ticksuffix="</sup>")),
                         title = paste("Interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    for (i in 1:length(unqX1)){
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX2)), y = unqX2, z = z_Delta[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX1)), z = z_Delta[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    inter_3d = fig

  }
  
  # Plotting monotherapies
  
  p1 = ggplot(data = df1, mapping = aes(x, mean)) +
    geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
    geom_point(data = mono1, mapping = aes(x,y), na.rm = T) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(min(-0.5,min(mono1$y)),max(1.5,max(mono1$y))) +
    labs(x = bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")"), y = "% Viability",
         title = paste("Monotherapy:",x$data$drug_names[1]),
         subtitle = paste0(x$data$experiment_ID)) +
    scale_x_continuous(labels = math_format(10^.x))
  
  p2 = ggplot(data = df2, mapping = aes(x, mean)) + 
    geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
    geom_point(data = mono2, mapping = aes(x,y), na.rm = T) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylim(min(-0.5,min(mono2$y)),max(1.5,max(mono2$y))) +
    labs(x = bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")"), y = "% Viability",
         title = paste("Monotherapy:",x$data$drug_names[2]),
         subtitle = paste0(x$data$experiment_ID)) +
    scale_x_continuous(labels = math_format(10^.x))
  
  
  ####################################################################################
  # Combinations (only plot these on the inside)
  ####################################################################################
  df = data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    f = as.vector(x$posterior_mean$f[-1,-1]),
    p0 = as.vector(x$posterior_mean$p0[-1,-1]),
    Delta = as.vector(x$posterior_mean$Delta[-1,-1])
  )
  
  # Breaks
  breaks = seq(0,1,by=0.1)

  #Function to create labels for legend
  breaklabel <- function(x, breaks){
    labels = paste0(100*breaks[1:10], "-", 100*breaks[2:11])
    labels[1:x]
  }
  
  # Response
  p3 = ggplot(data = df, aes(x = x1, y = x2, z = f)) +
    geom_contour_filled(breaks = breaks, show.legend = T, 
                        colour="dark grey", size=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(palette=viridis::viridis, values=breaklabel(10,breaks),
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
  
  # Noninteraction
  p4 = ggplot(data = df, aes(x = x1, y = x2, z = p0)) +
    geom_contour_filled(breaks = breaks, show.legend = T, 
                        colour="dark grey", size=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(palette=viridis::viridis, values=breaklabel(10,breaks),
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
  
  
  #Interaction surface has a different color scale
  Delta_col_palette <- inlmisc::GetColors(scheme = "sunset")
  # Breaks
  eps = 0.05
  breaks = c(seq(-1,-0.1,by=0.1),0-eps,0+eps,seq(0.1,1,by=0.1))
  
  #Function to create labels for legend
  breaklabel <- function(x, breaks){
    labels = paste0(100*breaks[1:10], "-", 100*breaks[2:11])
    labels = c(labels,paste0(100*breaks[11:21],"-",100*breaks[12:22]))
    labels[1:x]
  }
  
  p5 = ggplot(data = df, aes(x = x1, y = x2, z = Delta)) +
    geom_contour_filled(breaks = breaks, show.legend = T, 
                        colour="dark grey", size=.5, linetype = "dashed") +
    geom_point(aes(x = x1, y = x2),combination, color = "black", alpha = 0.1, size = 2,
               shape = 4) +
    scale_fill_manual(palette=Delta_col_palette, values=breaklabel(21,breaks),
                      name="Interaction %", drop = FALSE) +
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
  
  ####################################################################################
  # Summary statistics
  ####################################################################################
  # DSS scores
  df = data.frame(
    dss = c(posterior$dss_1,posterior$dss_2),
    idx = factor(c(rep(x$data$drug_names[1],n.save),rep(x$data$drug_names[2],n.save)),
                 levels = c(x$data$drug_names[1],x$data$drug_names[2]))
  )
  p6 = ggplot(df, aes(x = dss, y = idx, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, gradient_lwd = 1., from = 0, to = 100) +
    scale_fill_viridis_c(name = "DSS", option = "C",limits=c(0,100)) +
    labs(title = 'Estimated drug sensitivity scores',
         subtitle = x$data$experiment_ID,
         y = "") +
    xlim(0,100) +
    scale_y_discrete(limits = unique(rev(df$idx))) +
    theme_ridges(font_size = 13, grid = FALSE)
  
  
  # Combination scores
  df = data.frame(
    rVUS = c(posterior$rVUS_f, posterior$rVUS_p0,
             abs(posterior$VUS_syn), posterior$VUS_ant),
    idx = factor(c(rep('rVUS_f',n.save),rep('rVUS_p0',n.save),
                      rep('VUS_syn',n.save),rep('VUS_ant',n.save)),
                 levels = c('rVUS_f','rVUS_p0','VUS_syn','VUS_ant'))
  )
  
  p7 = ggplot(df, aes(x = rVUS, y = idx, fill = stat(x))) +
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
    
    # 3D plots
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
