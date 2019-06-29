#' Predict function for BayeSyneRgy object
#' 
#' @description A function for predicting response and interaction at an unobserved location using interpolation from the Bayesian semi-parametric model.
#' 
#' @param object The fitted model from \code{\link{BayeSyneRgy}}.
#' @param x_mat_new Matrix of concentrations to predict.
#' @param log10_conc logical; if true, concentrations must be given on the log10 scale, and monotherapies must have value -Inf for the drug with concentration zero.
#' @param cred_level The credible level for equal-tailed credible intervals.
#' @param plot_pred logical; if true, posterior predictive densities are plotted.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details 
#' This function extends \code{predict} to the \code{BayeSyneRgy} class. The function provides point estimates and uncertainty quantifiaction for new points on the surface, as specified in the \code{x_mat_new}. In addition, the function can provide a visualisation of the predictions in terms of the full posterior predictive distribution.
#' 
#' 
#' @export

predict.BayeSyneRgy <- function(object, x_mat_new, log10_conc = FALSE, cred_level = 0.95, plot_pred = FALSE,...){
  
  data <- object$data
  MCMC_Output <- object$OUTPUT$MCMC_Output
  
  # Removed attach here, so need to declare some variables properly
  x_mat <- data$x_mat
  y_mat <- data$y_mat
  type <- object$type
  drug_names <- object$drug_names
  SLOPE_1 <- MCMC_Output$SLOPE_1
  SLOPE_2 <- MCMC_Output$SLOPE_2
  EC50_1 <- MCMC_Output$EC50_1
  EC50_2 <- MCMC_Output$EC50_2
  GAMMA_0 <- MCMC_Output$GAMMA_0
  GAMMA_1 <- MCMC_Output$GAMMA_1
  GAMMA_2 <- MCMC_Output$GAMMA_2
  B <- MCMC_Output$B
  C <- MCMC_Output$C
  B_K1 <- MCMC_Output$B_K1
  B_K2 <- MCMC_Output$B_K2
  B_K3 <- MCMC_Output$B_K3
  B1 <- MCMC_Output$B1
  B2 <- MCMC_Output$B2
  ELL <- MCMC_Output$ELL
  SIGMA2_F <- MCMC_Output$SIGMA2_F
  NU <- MCMC_Output$NU
  ALPHA <- MCMC_Output$ALPHA
  S2_EC50_1 <- MCMC_Output$S2_EC50_1
  S2_EC50_2 <- MCMC_Output$S2_EC50_2
  S2_GAMMA_0 <- MCMC_Output$S2_GAMMA_0
  S2_GAMMA_1 <- MCMC_Output$S2_GAMMA_0
  S2_GAMMA_2 <- MCMC_Output$S2_GAMMA_0
  S2_EPS <- MCMC_Output$S2_EPS
  LPML <- MCMC_Output$LPML
  
  n_save <- length(MCMC_Output[[1]])
  x1 <- unique(x_mat[,1])
  x2 <- unique(x_mat[,2])
  
  #Predict response and interaction at given new covariate points
  
  # Check concentration scale (we want log10)
  if(!log10_conc){
    if(min(x_mat_new) == 0){
      x_mat_new[x_mat_new[,1] == 0,1] <- min(x_mat_new[x_mat_new[,1] > 0,1])/100
      x_mat_new[x_mat_new[,2] == 0,2] <- min(x_mat_new[x_mat_new[,2] > 0,2])/100
      
      x_mat_new <- log10(x_mat_new)
    }
  }
  
  #Only unique pairs
  x_mat_new <- x_mat_new[!duplicated(x_mat_new),]
  drug_names <- colnames(x_mat_new)
  n_pred <- dim(x_mat_new)[1]
  x1_new <- x_mat_new[,1]
  x2_new <- x_mat_new[,2]
  n1_new <- length(x1_new)
  n2_new <- length(x2_new)
  
  
  #We need the indicator function for the existence of interactions
  id_new <- rep(1, n_pred)
  id_new[x_mat_new[,1] <= min(x_mat[,1])] <- 0
  id_new[x_mat_new[,2] <= min(x_mat[,2])] <- 0
  
  p_ij_mat <- matrix(0,n_pred,n_save)
  Int_mat <- matrix(0,n_pred,n_save)
  p_0_mat <- matrix(0,n_pred,n_save)
  
  #For plots
  n_grid <- 100
  yij_new <- seq(0, 1, length = n_grid)
  Int_new <- seq(-1, 1, length = 2*n_grid)
  yij_pred_mat <- matrix(0, n_pred, n_grid)
  Int_pred_mat <- matrix(0, n_pred, 2*n_grid)
  
  if(type == 1){#Splines model
    
    #We need to create the new spline matrices from old knots
    #Dim 1
    deg <- 3 #degree of the spline
    ndx <- 5
    xl <- min(unique(x1_new)) - 0.1
    xr <- max(unique(x1_new)) + 0.1
    # Construct a B-spline basis of degree 'deg'
    dx <- (xr - xl) / ndx
    knots_1 <- seq(xl - dx * deg, xr + dx * deg, dx)
    K1 <- length(knots_1)
    P <- matrix(0,n1_new,K1)
    for(i1 in 1:n1_new){
      # Truncated p-th power function
      P[i1,] <- (x1_new[i1] - knots_1)^ deg * (x1_new[i1] > knots_1)
    }
    D <- diff(diag(K1), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B_K1_new <- (-1) ^ (deg + 1) * P %*% t(D)
    K1 <- dim(B_K1_new)[2]
    # matplot(x1_new,B_K1_new, type = "l")
    
    #Dim
    deg <- 3 #degree of the spline
    ndx <- 5
    xl <- min(unique(x2_new)) - 0.1
    xr <- max(unique(x2_new)) + 0.1
    # Construct a B-spline basis of degree 'deg'
    dx <- (xr - xl) / ndx
    knots_2 <- seq(xl - dx * deg, xr + dx * deg, dx)
    K2 <- length(knots_2)
    P <- matrix(0,n2_new,K2)
    for(i2 in 1:n2_new){
      # Truncated p-th power function
      P[i2,] <- (x2_new[i2] - knots_2)^ deg * (x2_new[i2] > knots_2)
    }
    D <- diff(diag(K2), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B_K2_new <- (-1) ^ (deg + 1) * P %*% t(D)
    K2 <- dim(B_K2_new)[2]
    # matplot(x2_new,B_K2_new, type = "l")
    
    #Spline for interaction
    #This splines are product of the previous ones
    B_K3_new <- array(0, dim = c(K1,K2,n_pred))
    for(i1 in 1:K1){
      for(i2 in 1:K2){
        for(i in 1:n_pred){
          B_K3_new[i1,i2,i] <- B_K1_new[i,i1] * B_K2_new[i,i2]
        }
      }
    }
    
    for(g in 1:n_save){
      f_1_g <- (1 + 10^(SLOPE_1[g] * (x1_new - EC50_1[g])))^(-1)
      f_1_g[1] <- 1
      f_2_g <- (1 + 10^(SLOPE_2[g] * (x2_new - EC50_2[g])))^(-1)
      f_2_g[1] <- 1
      p_0_g <- f_1_g * f_2_g
      
      B_aux <- rep(0, n_pred)
      for(i in 1:n_pred){
        B_aux[i] <- sum(matrix(C[g,], K1, K2) * B_K3_new[,,i])
      }
      
      Delta <- GAMMA_0[g] + GAMMA_1[g] * x1_new + GAMMA_2[g] * x2_new + B_aux
      
      Delta_trans <- - p_0_g * (1 + exp(B1[g]*Delta))^(-1) + (1 - p_0_g) * (1 + exp(-B2[g]*Delta))^(-1)
      
      #No interaction with only one drug (based on fitted concentrations)
      Delta_trans <- Delta_trans * id_new
      
      Int_mat[,g] <- Delta_trans
      p_0_mat[,g] <- p_0_g
      p_ij_mat[,g] <- p_0_g + Delta_trans
      
      for(i in 1:n_pred){
        yij_pred_mat[i,] <- yij_pred_mat[i,] + dnorm(yij_new, p_ij_mat[i,g], sqrt(S2_EPS[g]))
        Int_pred_mat[i,] <- Int_pred_mat[i,] + dnorm(Int_new, p_ij_mat[i,g] - p_0_mat[i,g], sqrt(S2_EPS[g]))
      }
    }
    
  }else{#GP models (assuming hyperp random, otherwise they will just be saved as columns of ones)
    
    n3 <- dim(x_mat)[1]
    x_all <- rbind(x_mat, x_mat_new)
    for(g in 1:n_save){
      f_1_g <- (1 + 10^(SLOPE_1[g] * (x1_new - EC50_1[g])))^(-1)
      f_1_g[1] <- 1
      f_2_g <- (1 + 10^(SLOPE_2[g] * (x2_new - EC50_2[g])))^(-1)
      f_2_g[1] <- 1
      p_0_g <- f_1_g * f_2_g
      
      if(type == 2){#Double exponential
        model_spec <- list(type = type, ell = ELL[g], sigma2_f = SIGMA2_F[g])
      }
      if(type == 3){#Matern
        model_spec <- list(type = type, ell = ELL[g], nu = NU[g], sigma2_f = SIGMA2_F[g])
      }
      if(type == 4){#Rational Quadratic
        model_spec <- list(type = type, ell = ELL[g], alpha = ALPHA[g], sigma2_f = SIGMA2_F[g])
      }
      
      Kxx_all <- GPKernel(x_all, model_spec)
      
      Kxx_new <- Kxx_all[(n3+1):(n3+n_pred), (n3+1):(n3+n_pred)]
      Kxx_gg <- Kxx_all[1:n3, 1:n3]
      Kxx_newg <- Kxx_all[(n3+1):(n3+n_pred), 1:n3]
      
      aux_B <- Kxx_newg %*% solve(Kxx_gg)
      var_B <- Kxx_new - aux_B %*% t(Kxx_newg)
      var_B <- (var_B + t(var_B))/2
      mean_B <- aux_B %*% B[g,]
      B_vec_new <- rmvn(n = 1, c(mean_B), SIGMA2_F[g] * var_B)
      Delta <- GAMMA_0[g] + GAMMA_1[g] * x1_new + GAMMA_2[g] * x2_new + B_vec_new
      
      Delta_trans <- - p_0_g * (1 + exp(B1[g]*Delta))^(-1) + (1 - p_0_g) * (1 + exp(-B2[g]*Delta))^(-1)
      #No interaction with only one drug (based on fitted concentrations)
      Delta_trans <- Delta_trans * id_new
      
      Int_mat[,g] <- Delta_trans
      p_0_mat[,g] <- p_0_g
      p_ij_mat[,g] <- p_0_g + Delta_trans
      
      for(i in 1:n_pred){
        yij_pred_mat[i,] <- yij_pred_mat[i,] + dnorm(yij_new, p_ij_mat[i,g], sqrt(S2_EPS[g]))
        Int_pred_mat[i,] <- Int_pred_mat[i,] + dnorm(Int_new, p_ij_mat[i,g] - p_0_mat[i,g], sqrt(S2_EPS[g]))
      }
    }
    
  }
  
  #Average surfaces
  Int_pred <- round(t(apply(Int_mat,1,quantile, probs = c((1-cred_level)/2, 0.5, 1-(1-cred_level)/2))),4)
  p_ij_pred <- round(t(apply(p_ij_mat,1,quantile, probs = c((1-cred_level)/2, 0.5, 1-(1-cred_level)/2))),4)
  
  prediction_out <- cbind(p_ij_pred, Int_pred)
  
  colnames(prediction_out) <- paste(paste0(c(rep("Resp (", 3), rep("Int (", 3)), colnames(prediction_out)), ")", sep = "")
  rownames(prediction_out) <- paste("(x1, x2) = (", as.character(round(x1_new,2)), ", ", as.character(round(x2_new,2)), ")", sep = "")
  
  #Print out
  print(prediction_out)
  
  if(plot_pred){
    readline("Press key for next plot")
    
    ###
    #Plot predictive distirbutions of response and interactions
    graphics.off()
    
    ### Predictive distirbutions of response
    #Colors based on response median value
    pred_colors <- viridis(n = n_pred+1)[1 + round(n_pred*(p_ij_pred[,2] - min(p_ij_pred[,2]))/diff(range(p_ij_pred[,2])))]
    
    yij_pred_mat <- yij_pred_mat/n_save
    max_y <- max(yij_pred_mat)
    
    plot(0, 0, col = "white", xlim = c(0,1), ylim = c(0,max_y), xlab = "", ylab = "", main = "Predictive densities of responses", cex.main = 2)
    for(i in 1:n_pred){
      points(p_ij_pred[i,2], 0, col = pred_colors[i], pch = 20, cex = 3)
      lines(yij_new, yij_pred_mat[i,], col = pred_colors[i], lwd = 2.5)
    }
    legend("topleft", inset=c(0,0), legend = rownames(prediction_out), col = pred_colors, pch = 19, cex = 1, bty = "n")
    
    readline("Press key for next plot")
    
    ### Predictive distributions of interactions
    #Interaction surface has a different color scale
    Int_col_palette <- colorRampPalette(c("green", "yellow", "red"))
    pred_colors <- Int_col_palette(n = n_pred+1)[1 + round(n_pred*(Int_pred[,2] - min(Int_pred[,2]))/diff(range(Int_pred[,2])))]
    
    Int_pred_mat <- Int_pred_mat/n_save
    max_y <- max(Int_pred_mat)
    
    plot(0, 0, col = "white", xlim = c(-1,1), ylim = c(0,max_y), xlab = "", ylab = "", main = "Predictive densities of interactions", cex.main = 2)
    for(i in 1:n_pred){
      points(Int_pred[i,2], 0, col = pred_colors[i], pch = 20, cex = 3)
      lines(Int_new, Int_pred_mat[i,], col = pred_colors[i], lwd = 2.5)
    }
    legend("topleft", inset=c(0,0), legend = rownames(prediction_out), col = pred_colors, pch = 19, cex = 1, bty = "n")
  }
}