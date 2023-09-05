# plotting functions one for each 
# data("mathews_DLBCL")
# y_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[1]]
# x_mat <- mathews_DLBCL$`ispinesib + ibrutinib`[[2]]
# fit <- bayesynergy(y_mat,x_mat)
fit <- bayesynergy(y_mat, x_mat, method = 'vb') # this is faster


# test run ----

p1_p2 <- plot_monotherapy(fitobj = fit)
p1_p2$drug1
p1_p2$drug2

p3 <- plot_synergy_response(fitobj = fit)
p3

p4 <- plot_synergy_noninteract(fitobj = fit)
p4


p5 <- plot_synergy_interact(fitobj = fit)
p5

p6 <- plot_drug_sensitivity(fitobj = fit)
p6

p7 <- plot_drug_combscore(fitobj = fit)
p7


# p1, p2: mono ----

plot_monotherapy <- function(fitobj){
  
  x <- fitobj
  # x <- fit
  
  # processing ----# 
  # extract posterior object
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  if (!x$model$lower_asymptotes){
    posterior$la_1 <- rep(0,n.save)
    posterior$la_2 <- rep(0,n.save)
  }
  
  # unique doses (log) for two drugs
  unqX1 <- log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 <- log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  
  dx1 <- mean(diff(unqX1))
  dx2 <- mean(diff(unqX2))
  nrep <- ncol(as.matrix(x$data$y))
  
  # monotherapy for 2 drugs
  mono1 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  mono1 <- mono1[!is.na(mono1$y),] # rm NA
  
  mono2 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  mono2 <- mono2[!is.na(mono2$y),] # rm NA
  
  
  # set grid ----# 
  grid.size <- 100
  # dose
  x.seq1 <- seq(min(unqX1)-dx1,max(unqX1)+dx1,length.out = grid.size)
  x.seq2 <- seq(min(unqX2)-dx2,max(unqX2)+dx2,length.out = grid.size)
  # response
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
  
  if (x$model$robust){
    # robust: use median
    
    p1 <- ggplot(data = df1, mapping = aes(x, median)) 
    p2 <- ggplot(data = df2, mapping = aes(x, median)) 
    
    
  }else{
    # not robust: use mean
    p1 <- ggplot(data = df1, mapping = aes(x, mean)) 
    p2 <- ggplot(data = df2, mapping = aes(x, mean)) 
    
  }
  
  # p1
  p1 <- p1 + geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) 
  p1 <- p1 + geom_point(data = mono1, mapping = aes(x,y), na.rm = T) 
  p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed", color = "red") 
  p1 <- p1 + ylim(min(-0.5,min(mono1$y)),max(1.5,max(mono1$y))) 
  p1 <- p1 + labs(x = bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")"), 
                  y = "Viability",
                  title = paste("Monotherapy:",x$data$drug_names[1]),
                  subtitle = paste0(x$data$experiment_ID)) 
  p1 <- p1 + scale_x_continuous(labels = math_format(10^.x))
  
  
  # p2
  p2 <- p2 + geom_smooth(stat = "identity", aes(ymin = lower, ymax = upper)) 
  p2 <- p2 + geom_point(data = mono2, mapping = aes(x,y), na.rm = T) 
  p2 <- p2 + geom_hline(yintercept = 0, linetype = "dashed", color = "red") 
  p2 <- p2 + ylim(min(-0.5,min(mono2$y)),max(1.5,max(mono2$y))) 
  p2 <- p2 + labs(x = bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")"),
                  y = "Viability",
                  title = paste("Monotherapy:",x$data$drug_names[2]),
                  subtitle = paste0(x$data$experiment_ID)) 
  p2 <- p2 + scale_x_continuous(labels = math_format(10^.x))
  
  
  return(plots = list(
    drug1 = p1, 
    drug2 = p2))
  
}




# p3: response surface ----


plot_synergy_response <- function(fitobj){
  x <- fitobj
  # x <- fit
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  # We add lower-asymptote parameters if these are not estimated
  if (!x$model$lower_asymptotes){
    posterior$la_1 <- rep(0,n.save)
    posterior$la_2 <- rep(0,n.save)
  }
  
  unqX1 <- log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 <- log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  nrep <- ncol(as.matrix(x$data$y))
  
  # Need to find coordinates for the observed variables in this new coordinate system
  Xgrid <- expand.grid(unqX1,unqX2)
  Xgrid <- Xgrid[order(Xgrid[,"Var1"],Xgrid[,"Var2"]),]
  
  
  # monotherapy for 2 drugs
  mono1 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  mono1 <- mono1[!is.na(mono1$y),] # rm NA
  
  mono2 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  mono2 <- mono2[!is.na(mono2$y),] # rm NA
  
  df <- data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    f = as.vector(x$posterior_mean$f[-1,-1])
  )
  
  # Breaks
  breaks <- seq(0,1,by=0.1)
  
  # p3
  p <- ggplot(data = df, aes(x = x1, y = x2, z = f)) 
  p <- p + geom_contour_filled(breaks = breaks, 
                               show.legend = T,
                               colour="dark grey", 
                               linewidth=.5, 
                               linetype = "dashed") 
  p <- p + geom_point(data = df, aes(x = x1, y = x2),
                      color = "black", alpha = 0.1, size = 2, shape = 4) 
  p <- p + scale_fill_manual(values=viridis::viridis(10), 
                             labels= f_breaklabel(x=10,breaks),
                             name="% Viability", drop = FALSE) 
  p <- p + scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) 
  p <- p + ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 1, 
                    x = -Inf, y = mono2$x, color = 'black') 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 1, 
                    x = mono1$x, y = -Inf, color = 'black') 
  p <- p + coord_cartesian(clip = 'off') 
  p <- p + labs(title = "Response surface",
                subtitle = paste0(x$data$experiment_ID,": ",
                                  x$data$drug_names[1]," + ", 
                                  x$data$drug_names[2]))
  
  return(p)
  
}






# p4: response non interact ----

plot_synergy_noninteract <- function(fitobj){
  
  x <- fitobj
  # x <- fit
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  # We add lower-asymptote parameters if these are not estimated
  if (!x$model$lower_asymptotes){
    posterior$la_1 <- rep(0,n.save)
    posterior$la_2 <- rep(0,n.save)
  }
  
  unqX1 <- log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 <- log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  nrep <- ncol(as.matrix(x$data$y))
  
  # Need to find coordinates for the observed variables in this new coordinate system
  Xgrid <- expand.grid(unqX1,unqX2)
  Xgrid <- Xgrid[order(Xgrid[,"Var1"],Xgrid[,"Var2"]),]
  
  # monotherapy for 2 drugs
  mono1 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  mono1 <- mono1[!is.na(mono1$y),] # rm NA
  
  mono2 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  mono2 <- mono2[!is.na(mono2$y),] # rm NA
  
  
  df <- data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    p0 = as.vector(x$posterior_mean$p0[-1,-1])
  )
  
  
  # Breaks
  breaks <- seq(0,1,by=0.1)
  
  # p4
  p <- ggplot(data = df, aes(x = x1, y = x2, z = p0)) 
  p <- p + geom_contour_filled(breaks = breaks, 
                               show.legend = T,
                               colour="dark grey", 
                               linewidth=.5, 
                               linetype = "dashed") 
  p <- p + geom_point(data = df, aes(x = x1, y = x2),
                      color = "black", alpha = 0.1, size = 2, shape = 4) 
  p <- p + scale_fill_manual(values=viridis::viridis(10), 
                             labels = f_breaklabel(x=10,breaks),
                             name="% Viability", drop = FALSE) 
  p <- p + scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) 
  p <- p + ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, 
                    x = -Inf, y = mono2$x, color = 'black') 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, 
                    x = mono1$x, y = -Inf, color = 'black') 
  p <- p + coord_cartesian(clip = 'off') 
  p <- p + labs(title = "Non-interaction surface",
                subtitle = paste0(x$data$experiment_ID,": ",
                                  x$data$drug_names[1]," + ", 
                                  x$data$drug_names[2]))
  
  return(p)
  
}




# p5: interaction ----

plot_synergy_interact <- function(fitobj){
  
  x <- fitobj
  # x <- fit
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  # We add lower-asymptote parameters if these are not estimated
  if (!x$model$lower_asymptotes){
    posterior$la_1 <- rep(0,n.save)
    posterior$la_2 <- rep(0,n.save)
  }
  
  unqX1 <- log10(sort(unique(x$data$x[,1])))[-1] # Removing -Inf here
  unqX2 <- log10(sort(unique(x$data$x[,2])))[-1] # Removing -Inf here
  nrep <- ncol(as.matrix(x$data$y))
  
  
  # Need to find coordinates for the observed variables in this new coordinate system
  Xgrid <- expand.grid(unqX1,unqX2)
  Xgrid <- Xgrid[order(Xgrid[,"Var1"],Xgrid[,"Var2"]),]
  
  
  df <- data.frame(
    x1 = Xgrid[,1],
    x2 = Xgrid[,2],
    Delta = as.vector(x$posterior_mean$Delta[-1,-1])
  )
  
  
  # monotherapy for 2 drugs
  mono1 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),1]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,2]==0) & (x$data$x[,1] != 0)),])
  )
  mono1 <- mono1[!is.na(mono1$y),] # rm NA
  
  mono2 <- data.frame(
    x = rep(log10(x$data$x[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),2]),nrep),
    y = as.vector(as.matrix(x$data$y)[which((x$data$x[,1]==0) & (x$data$x[,2] != 0)),])
  )
  mono2 <- mono2[!is.na(mono2$y),] # rm NA
  
  
  # Interaction surface has a different color scale
  # Delta_col_palette <- inlmisc::GetColors(scheme = "sunset")
  Delta_col_palette <- c("#354B99", "#4062A8", "#4A7AB7", "#5D90C2", "#6DA5CC", "#83B8D7", "#98CAE0",
                         "#ACD6E7", "#C2E3EE", "#D7E7DD", "#EAEBCC", "#F5E3AB", "#FED98B", "#FEC678",
                         "#FDB366", "#FA9958", "#F57D4A", "#EA603B", "#DC3D2D", "#C0242A", "#A50026")
  # Breaks
  eps <- 0.05
  breaks <- c(seq(-1,-0.1,by=0.1),0-eps,0+eps,seq(0.1,1,by=0.1))
  
  p <- ggplot(data = df, aes(x = x1, y = x2, z = Delta)) 
  p <- p + geom_contour_filled(breaks = breaks, 
                               show.legend = T,
                               colour="dark grey", 
                               size=.5, 
                               linetype = "dashed") 
  p <- p + geom_point(data = df, aes(x = x1, y = x2), 
                      color = "black", alpha = 0.1, size = 2, shape = 4) 
  p <- p + scale_fill_manual(values = Delta_col_palette, 
                             labels = f_breaklabel(x=21, breaks = breaks), 
                             name="% Interaction", drop = FALSE) 
  p <- p + scale_x_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + scale_y_continuous(expand=c(0,0),labels = math_format(10^.x)) 
  p <- p + xlab(bquote(.(x$data$units[1])~" ("~.(x$data$drug_names[1])~")")) 
  p <- p + ylab(bquote(.(x$data$units[2])~" ("~.(x$data$drug_names[2])~")")) 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, 
                    x = -Inf, y = mono2$x, color = 'black') 
  p <- p + annotate(geom = "point", shape = 4, alpha = 0.1, size = 2, 
                    x = mono1$x, y = -Inf, color = 'black') 
  p <- p + coord_cartesian(clip = 'off') 
  p <- p + labs(title = "Interaction surface",
                subtitle = paste0(x$data$experiment_ID,": ",
                                  x$data$drug_names[1]," + ", 
                                  x$data$drug_names[2])) 
  p <- p + guides(fill=guide_legend(ncol=1))
  
  return(p)
  
}



# p6: dss ----

plot_drug_sensitivity <- function(fitobj){
  
  x <- fitobj
  # x <- fit
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  
  df <- data.frame(
    dss = c(posterior$dss_1,posterior$dss_2),
    idx = factor(c(rep(x$data$drug_names[1],n.save),
                   rep(x$data$drug_names[2],n.save)),
                 levels = c(x$data$drug_names[1],
                            x$data$drug_names[2]))
  )
  
  p <- ggplot(df, aes(x = dss, y = idx, fill = stat(x))) 
  p <- p + geom_density_ridges_gradient(scale = 3, gradient_lwd = 1., from = 0, to = 100) 
  p <- p + scale_fill_viridis_c(name = "DSS", option = "C",limits=c(0,100)) 
  p <- p + labs(title = 'Estimated drug sensitivity scores',
                subtitle = x$data$experiment_ID,
                y = "") 
  p <- p + xlim(0,100) 
  p <- p + scale_y_discrete(limits = unique(rev(df$idx))) 
  p <- p + theme_ridges(font_size = 13, grid = FALSE)
  
  
  return(p)
}


# p7: dcs ----
plot_drug_combscore <- function(fitobj){
  
  x <- fitobj
  # x <- fit
  posterior <- rstan::extract(x$stanfit)
  n.save <- length(posterior$lp__)
  
  df <- data.frame(
    rVUS = c(posterior$rVUS_f, posterior$rVUS_p0,
             abs(posterior$VUS_syn), posterior$VUS_ant),
    idx = factor(c(rep('rVUS_f',n.save),rep('rVUS_p0',n.save),
                   rep('VUS_syn',n.save),rep('VUS_ant',n.save)),
                 levels = c('rVUS_f','rVUS_p0','VUS_syn','VUS_ant'))
  )
  
  p <- ggplot(df, aes(x = rVUS, y = idx, fill = stat(x))) 
  p <- p + geom_density_ridges_gradient(scale = 2, gradient_lwd = 1., from = 0, to = 100) 
  p <- p + scale_fill_viridis_c(name = "", option = "C",limits=c(0,100))
  p <- p + labs(title = 'Estimated drug combination scores',
                subtitle = paste0(x$data$experiment_ID,": ",
                                  x$data$drug_names[1]," + ", 
                                  x$data$drug_names[2]),
                y = "", x = "") 
  p <- p + xlim(0,100)
  # p <- p + xlab("") 
  p <- p + scale_y_discrete(limits = unique(rev(df$idx)),
                            labels=c("Antagonism\nVUS(\u0394\u207a)",
                                     "Synergy\n|VUS(\u0394\u207b)|",
                                     "Bliss\nefficacy\nrVUS(p\u2080)",
                                     "Overall\nefficacy\nrVUS(f)"))
  #expression(paste("Overall efficacy (",rVUS(f),")"))),
  #expand = expansion(add = c(0.01))) 
  # labels=c("Antagonism","Synergy","Overall Interaction","Overall Efficacy")) +
  p <- p +theme_ridges(font_size = 12, grid = FALSE)
  
  return(p)
  
}





# ________ ----
# util ----

# function to make lables
# used for p3, p4, p5
f_breaklabel <- function(x, breaks){
  n <- length(breaks)
  label_1 <- paste0("[",100*breaks[1],",",100*breaks[2],"]")  
  label_x <- c(paste0("(",100*breaks[2:(n-1)], ",", 100*breaks[3:n],"]"))
  labels <- c(label_1, label_x)
  labels <- labels[1:x]
  return(labels)
}

# f_breaklabel(x=21, breaks = breaks)
# f_breaklabel(x=10, breaks = breaks)






