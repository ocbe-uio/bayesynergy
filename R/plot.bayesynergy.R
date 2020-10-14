#' Plot function for a bayesynergy object
#' 
#' @description A function for plotting synergy surfaces and summary statistics of a \code{bayesynergy} object.
#' 
#' @param x bayesynergy object, the result of \code{\link{bayesynergy}}.
#' @param plot3D logical; if TRUE interactive 3d surface plots of dose response function are displayed.
#' @param save_plot logical; if TRUE plots are saved locally.
#' @param path string; path for saving plots, if NULL defaults to work directory.
#' @param plotdevice string; device for saving plots locally, must be 'pdf' or 'png'
#' @param ... further arguments passed on to device for plotting, useful for setting width and height of saved plots.
#' 
#' @details 
#' This function extends \code{plot} to draw response and interaction surfaces of the fitted model. Both three-dimensional interactive plots, and two-dimensional contour plots can be displayed.
#' 
#' 
#' @examples 
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
#' x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
#' fit <- bayesynergy(y_mat,x_mat)
#' plot(fit)
#' 
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point ylim xlim labs scale_x_continuous stat_contour_filled geom_contour_filled scale_fill_manual scale_y_continuous scale_y_discrete xlab ylab scale_fill_viridis_c guides guide_legend
#' @importFrom ggridges geom_density_ridges_gradient theme_ridges
#' @importFrom plotly plot_ly add_surface add_trace %>% plotly_build
#' @importFrom scales math_format
#' @importFrom gridExtra grid.arrange
#' @importFrom inlmisc GetColors
#' 
#' @export 


plot.bayesynergy <- function(x, plot3D = T, save_plot = FALSE, path = NULL, plotdevice = "pdf", ...){
  
  # Check for valid plotdevice
  if (!(plotdevice %in% c("pdf","png"))){
    stop("plotdevice must be one of {'pdf','png'}")
  }
  # Setting default path if note given
  if (is.null(path)){
    path = getwd()
  }
  # Put dirmark on file path
  path = Sys.glob(path,dirmark = T)

  
  #Close all graphic tools open
  # graphics.off()
  # 
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
  
  combination = data.frame(
    x1 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),1]),nrep),
    x2 = rep(log10(x$data$x[which((x$data$x[,1]!=0) & (x$data$x[,2]!=0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]!=0) & (x$data$x[,2] != 0)),])
  ) 
  # Remove NA here
  idx = !is.na(combination$y)
  combination = combination[idx,]
  
  # #################################
  # ## We include a 3d-plot option ##
  # #################################
  if(plot3D){
    
    # Response
    z_response = x$posterior_mean$p_ij[-1,-1]
    fig = plot_ly(x = unqX1, y = unqX2, z = z_response)
    fig = fig %>% add_surface(cmin=0,cmax=1)
    fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
                            type = "scatter3d", mode = "markers",
                            marker = list(size=5,color="black",symbol=104))
    fig = fig %>% plotly::layout(scene = list(zaxis = list(range=c(min(0,min(combination$y)),max(1,max(combination$y))),
                                                   title="% Viability"),
                                      xaxis = list(title=x$data$drug_names[1]),
                                      yaxis = list(title=x$data$drug_names[2])),
                         title = paste("Response surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    response_3d = fig
    # 
    
    # Noninteraction
    z_p0 = x$posterior_mean$p_0[-1,-1]
    fig = plot_ly(x = unqX1, y = unqX2, z = z_p0)
    fig = fig %>% add_surface(cmin=0,cmax=1)
    fig = fig %>% add_trace(x = combination$x1, y = combination$x2, z = combination$y,
                            type = "scatter3d", mode = "markers",
                            marker = list(size=5,color="black",symbol=104))
    fig = fig %>% plotly::layout(scene = list(zaxis= list(range=c(min(0,min(combination$y)),max(1,max(combination$y))),
                                                  title="% Viability"),
                                      xaxis = list(title=x$data$drug_names[1]),
                                      yaxis = list(title=x$data$drug_names[2])),
                         title = paste("Non-interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    
    noninter_3d = fig
    
    # Interaction
    z_Delta = x$posterior_mean$Delta[-1,-1]
    fig = plot_ly()
    fig = fig %>% add_trace(x = unqX1, y = unqX2, z = z_Delta,type = "surface",
                            colorscale = list(c(0,0.5, 1), c("2166AC","EAECCC", "B2182B")),cmin=-1,cmax=1)
    fig = fig %>% plotly::layout(scene = list(zaxis= list(range=c(-1,1),
                                                  title="Interaction"),
                                      xaxis = list(title=x$data$drug_names[1]),
                                      yaxis = list(title=x$data$drug_names[2])),
                         title = paste("Interaction surface:",x$data$experiment_ID,":",x$data$drug_names[1],"+",x$data$drug_names[2]))
    inter_3d = fig
    # 
  }
  
  # Plotting monotherapies
  
  ####################################################################################
  # Monotherapies
  ####################################################################################
  grid.size = 100
  x.seq1 = seq(min(unqX1)-dx1,max(unqX1)+dx1,length.out = grid.size)
  x.seq2 = seq(min(unqX2)-dx2,max(unqX2)+dx2,length.out = grid.size)
  
  y.seq1 = matrix(NA,nrow=grid.size,ncol=n.save)
  y.seq2 = matrix(NA,nrow=grid.size,ncol=n.save)
  for (i in 1:grid.size){
    y.seq1[i,] = as.vector(posterior$la_1)+as.vector((1-posterior$la_1))/(1+10^(posterior$slope_1*(x.seq1[i]-posterior$ec50_1)))
    y.seq2[i,] = as.vector(posterior$la_2)+as.vector((1-posterior$la_2))/(1+10^(posterior$slope_2*(x.seq2[i]-posterior$ec50_2)))
  }
  df1 = data.frame(
    x = x.seq1,
    mean = apply(y.seq1,1,mean),
    lower = apply(y.seq1,1,quantile, probs=0.025),
    upper = apply(y.seq1,1,quantile, probs=0.975)
  )
  df2 = data.frame(
    x = x.seq2,
    mean = apply(y.seq2,1,mean),
    lower = apply(y.seq2,1,quantile, probs=0.025),
    upper = apply(y.seq2,1,quantile, probs=0.975)
  )
  
  p1 = ggplot(data = df1, mapping = aes(x, mean)) +
    geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
    geom_point(data = mono1, mapping = aes(x,y), na.rm = T) +
    ylim(min(-0.5,min(mono1$y)),max(1.5,max(mono1$y))) +
    labs(x = bquote(log[10](x)~" ("~.(x$data$drug_names[1])~")"), y = "% Viability",
         title = "Monotherapy drug 1",
         subtitle = paste0(x$data$experiment_ID,": ", x$data$drug_names[1]))+
    scale_x_continuous(labels = math_format(10^.x))
  
  
  p2 = ggplot(data = df2, mapping = aes(x, mean)) + 
    geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) +
    geom_point(data = mono2, mapping = aes(x,y), na.rm = T) +
    ylim(min(-0.5,min(mono2$y)),max(1.5,max(mono2$y))) +
    labs(x = bquote(log[10](x)~" ("~.(x$data$drug_names[2])~")"), y = "% Viability",
         title = "Monotherapy drug 2",
         subtitle = paste0(x$data$experiment_ID,": ", x$data$drug_names[2]))+
    scale_x_continuous(labels = math_format(10^.x))
  
  
  ####################################################################################
  # Combinations (only plot these on the inside)
  ####################################################################################
  df = data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    pij = as.vector(x$posterior_mean$p_ij[-1,-1]),
    p0 = as.vector(x$posterior_mean$p_0[-1,-1]),
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
  p3 = ggplot(data = df, aes(x = x1, y = x2, z = pij)) +
    geom_contour_filled(breaks = breaks, show.legend = T, 
                        colour="dark grey", size=.5, linetype = "dashed") +
    scale_fill_manual(palette=viridis::viridis, values=breaklabel(10,breaks),
                      name="% Viability", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(log[10](x)~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(log[10](x)~" ("~.(x$data$drug_names[2])~")")) +
    labs(title = "Response surface",
        subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2]))
    
    
  
  p4 = ggplot(data = df, aes(x = x1, y = x2, z = p0)) +
    geom_contour_filled(breaks = breaks, show.legend = T, 
                        colour="dark grey", size=.5, linetype = "dashed") +
    scale_fill_manual(palette=viridis::viridis, values=breaklabel(10,breaks),
                      name="% Viability", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(log[10](x)~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(log[10](x)~" ("~.(x$data$drug_names[2])~")")) +
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
    scale_fill_manual(palette=Delta_col_palette, values=breaklabel(21,breaks),
                      name="Interaction %", drop = FALSE) +
    scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) +
    xlab(bquote(log[10](x)~" ("~.(x$data$drug_names[1])~")")) +
    ylab(bquote(log[10](x)~" ("~.(x$data$drug_names[2])~")")) +
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
    rVUS = c(posterior$rVUS_p, posterior$rVUS_Delta,
             posterior$rVUS_syn, posterior$rVUS_ant),
    idx = factor(c(rep('rVUS_p',n.save),rep('rVUS_Delta',n.save),
                      rep('rVUS_syn',n.save),rep('rVUS_ant',n.save)),
                 levels = c('rVUS_p','rVUS_Delta','rVUS_syn','rVUS_ant'))
  )
  
  p7 = ggplot(df, aes(x = rVUS, y = idx, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, gradient_lwd = 1., from = 0, to = 100) +
    scale_fill_viridis_c(name = "rVUS", option = "C",limits=c(0,100)) +
    labs(title = 'Estimated drug combination scores',
         subtitle = paste0(x$data$experiment_ID,": ",x$data$drug_names[1]," + ", x$data$drug_names[2]),
         y = "") +
    xlim(0,100) +
    scale_y_discrete(limits = unique(rev(df$idx)), 
                     labels=c(expression(rVUS(Delta^"+")),expression(rVUS(Delta^"-")),
                              expression(rVUS(paste("|", Delta, "|"))),expression(rVUS))) +
    theme_ridges(font_size = 13, grid = FALSE)
  
  
  if (!save_plot){
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
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    grid.arrange(p1,p2,nrow=2) # Monotherapies
    dev.off()
    # DSS scores
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"DSS",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p6) # DSS scores 
    dev.off()
    # Response surface
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Response",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p3) # Response surface
    dev.off()
    # Non-interaction
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"NonInteraction",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p4) # Noninteraction surface
    dev.off()
    # Interaction
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"Interaction",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p5) # Interaction surface
    dev.off()
    # rVUS
    file.name = paste0(path,paste(x$data$experiment_ID,x$data$drug_names[1],x$data$drug_names[2],"rVUS",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p7) # rVUS scores
    dev.off()
  }
}
