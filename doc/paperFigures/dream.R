### Script for processing the DREAM dataset in supplementary materials 

# Input: 
#         ch1_training_combinations folder (from DREAM challenge)
#         Available as supplementary material for [1] after applying for access


library(tidyverse)
library(bayesynergy)


extract_and_fit = function(mask){
  path = "path_to/ch1_training_combinations/"
  files = list.files(path)
  
  # This is a clear synergy case
  cell.line = "DMS-114"
  drug.A = "MAP2K_3"
  drug.B = "BRAF_M2"
  # Find a not so clear case
  # cell.line = "MDA-MB-436"
  # drug.A = "MAP2K_3"
  # drug.B = "BRAF_M2"
  
  files = files[stringr::str_detect(files,cell.line)]
  files = files[stringr::str_detect(files,drug.A)]
  files = files[stringr::str_detect(files,drug.B)]
  
  data = read.csv(paste0(path,files[1]),header=F)
  colnames(data) = as.character(data[1,])
  data$drug.A = data[9,2]
  data$drug.B = data[10,2]
  data$cell.Line = data[13,2]
  data$drug.A.conc = data[,1]
  
  data = as_tibble(data[-c(1,8:13),-c(1,8)])
  
  data = type_convert(data)
  
  # Apply a filter?
  data[2:6,2:6] = data[2:6,2:6]*mask
  
  data_long = pivot_longer(data,cols=c(1:6),names_to="drug.B.conc")
  data_long = type_convert(data_long)
  
  final = list(x = as.matrix(data_long[,c("drug.A.conc","drug.B.conc")]),
               y = as.matrix(data_long[,"value"])/100,
               experiment_ID = as.character(data_long[1,3]),
               drug_names= c(as.character(data_long[1,1]),as.character(data_long[1,2])),
               control=list(adapt_delta=0.99,max_treedepth=20))
  
  fit = do.call(bayesynergy,c(final,method="sampling"))
  return(fit)
}


# First we fit only monotherapies
mask_mono = matrix(NA,5,5)
monos = extract_and_fit(mask_mono)
# Then we fit full data
mask_full = matrix(1,5,5)
full = extract_and_fit(mask_full)

# We fit line and cross
# We find IC20s for the two compounds
monos_mean = monos$posterior_mean
IC20_1 = (1/monos_mean$slope_1)*log10((1-monos_mean$`la_1[1]`)/(0.8-monos_mean$`la_1[1]`)-1)+monos_mean$log10_ec50_1
IC20_2 = (1/monos_mean$slope_2)*log10((1-monos_mean$`la_2[1]`)/(0.8-monos_mean$`la_2[1]`)-1)+monos_mean$log10_ec50_2
# Which is closest to our values?
unqx1 = unique(monos$data$x[,1])[-1]
unqx2 = unique(monos$data$x[,2])[-1]
id1 = which.min(abs(10^IC20_1-unqx1))
id2 = which.min(abs(10^IC20_2-unqx2))
# Line fit
mask_line = matrix(NA,5,5)
mask_line[id1,] = 1
line = extract_and_fit(mask_line)
# Cross fit
mask_cross = matrix(NA,5,5)
mask_cross[id1,] = 1
mask_cross[,id2] = 1
cross = extract_and_fit(mask_cross)
# Diagonal
mask_diagonal = matrix(NA,5,5)
diag(mask_diagonal) = 1
diagonal = extract_and_fit(mask_diagonal)
# Single point
mask_single = matrix(NA,5,5)
mask_single[id1,id2] = 1
single = extract_and_fit(mask_single)

# Plots
pdf("DREAM_single.pdf")
plot(density(rstan::extract(single$stanfit)$VUS_Delta),xlim=c(-50,50),xlab="Interaction", main = "Single")
dev.off()
pdf("DREAM_line.pdf")
plot(density(rstan::extract(line$stanfit)$VUS_Delta),xlim=c(-50,50),xlab="Interaction", main="Line")
dev.off()
pdf("DREAM_cross.pdf")
plot(density(rstan::extract(cross$stanfit)$VUS_Delta),xlim=c(-50,50),xlab="Interacion",main="Cross")
dev.off()
pdf("DREAM_diagonal.pdf")
plot(density(rstan::extract(diagonal$stanfit)$VUS_Delta),xlim=c(-50,50),xlab="Interaction",main="Diagonal")
dev.off()
pdf("DREAM_full.pdf")
plot(density(rstan::extract(full$stanfit)$VUS_Delta),xlim=c(-50,50),xlab="Interaction",main="Full")
dev.off()

#
plot(single,save_plots=T, plotdevice="png", width=7,height=7,units="in",res=600, plot3D=F)
plot(line,save_plots=T, plotdevice="png", width=7,height=7,units="in",res=600, plot3D=F)
plot(cross,save_plots=T, plotdevice="png", width=7,height=7,units="in",res=600, plot3D=F)
plot(diagonal,save_plots=T, plotdevice="png", width=7,height=7,units="in",res=600, plot3D=F)
plot(full,save_plots=T, plotdevice="png", width=7,height=7,units="in",res=600, plot3D=F)

# References
# [1] Menden, M.P., Wang, D., Mason, M.J. et al. Community assessment to advance computational prediction of cancer drug combinations in a pharmacogenomic screen. Nat Commun 10, 2674 (2019). https://doi.org/10.1038/s41467-019-09799-2

