## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_example,cache=T-----------------------------------------------------
library(bayesynergy)
data("mathews_DLBCL")
y = mathews_DLBCL[[1]][[1]]
x = mathews_DLBCL[[1]][[2]]
head(cbind(y,x))

## ----regular_fit, cache = T, warning = FALSE----------------------------------
fit = bayesynergy(y,x, drug_names = c("ibrutinib", "ispinesib"),
                  units = c("nM","nM"),bayes_factor = T)

## ----summary_fit,cache=T------------------------------------------------------
summary(fit)

## ----regular_plots, cache = T, warning = FALSE, message = FALSE, fig.dim = c(7,7), fig.show="hold", results="hide", fig.keep ="all"----
plot(fit, plot3D = F)

## ----3dresponse, echo = F, cache = T, warning = FALSE, message = FALSE, fig.dim = c(7,7)----
x = fit
# Creating some stuff needed for plots
library(plotly)

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
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX1)), y = unqX2, z = z_response[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX2)), z = z_response[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    
    response_3d = fig
    
    response_3d
  
  

## ----3dnoninteraction, echo = F, cache = T, warning = FALSE, message = FALSE, fig.dim = c(7,7)----

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
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX1)), y = unqX2, z = z_p0[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX2)), z = z_p0[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    noninter_3d = fig
    noninter_3d


## ----3dinteraction, echo = F, cache = T, warning = FALSE, message = FALSE, fig.dim = c(7,7)----

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
      fig = fig %>% add_trace(x = rep(unqX1[i],length(unqX1)), y = unqX2, z = z_Delta[,i]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    for (i in 1:length(unqX2)){
      fig = fig %>% add_trace(x = unqX1, y = rep(unqX2[i],length(unqX2)), z = z_Delta[i,]+0.003, type="scatter3d", mode="lines",
                              showlegend = F, line = list(color="grey", width = 1, dash = "dot"))
    }
    inter_3d = fig
    inter_3d


