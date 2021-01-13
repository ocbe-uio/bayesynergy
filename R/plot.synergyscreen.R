#' Plot function for a synergyscreen object
#' 
#' @description A function for plotting summary statistics of a \code{synergyscreen} object.
#' 
#' @param x An object of class \code{synergyscreen}, the result of \code{\link{synergyscreen}}.
#' @param groupbyExperimentID logical; if TRUE, individual plots are produced per experiment ID.
#' @param save_plots logical; if TRUE plots are saved locally.
#' @param path string; path for saving plots, if NULL defaults to work directory.
#' @param plotdevice string; device for saving plots locally, must be 'pdf' or 'png'
#' @param ... further arguments passed on to device for plotting, useful for setting width and height of saved plots.
#' 
#' @details 
#' This function extends \code{plot} to summarize results from large screens.
#' 
#' 
#' @examples
#' \dontrun{
#' library(bayesynergy)
#' data("mathews_DLBCL")
#' experiment1 = list(y = mathews_DLBCL[[1]][[1]], 
#' x = mathews_DLBCL[[1]][[2]], 
#' drug_names = c("ispinesib","ibrutinib"))
#' experiment2 = list(y = mathews_DLBCL[[2]][[1]], 
#' x = mathews_DLBCL[[2]][[2]], 
#' drug_names = c("canertinib","ibrutinib"))
#' experiments = list(experiment1,experiment2)
#' fit <- synergyscreen(experiments)
#' plot(fit)
#' }
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous facet_grid xlab ylab ggtitle labs guides theme unit element_text element_blank margin aes_string guide_legend scale_size scale_fill_gradient2 scale_colour_gradientn scale_fill_gradientn 
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom dplyr mutate %>% group_split n select
#' @importFrom grDevices cairo_pdf
#' @importFrom stats mad median
#' @importFrom utils tail
#' 
#' @export 


plot.synergyscreen <- function(x, groupbyExperimentID = T, save_plots = FALSE, path = NULL, plotdevice = "pdf", ...){
  
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
  
  # Extracting the dataframe needed from the output
  synscores = x$screenSummary
  
  # Creating a function to label top synergy scores
  in_quantile <- function(x, n = 5) {
    # Pick out to get top and bottom 5
    vec = rep(F,length(x))
    vec[tail(order(x),n)] = T
    return(vec)
  }
  
  # Do some preprocessing to make sure we are able to average across the pairs correctly
  # Particularly, if (drugA + drugB) and (drugB + drugA) are equal for two experimental IDs,
  # they should be averaged, as they represent the same pair
  # We do this alphabetically such that drug A comes earlier in the alphabet than drug B
  
  synscores <- synscores %>% mutate(`drugATMP` = `Drug A`, `drugBTMP` = `Drug B`)
  synscores <- synscores %>% mutate(`Drug A` = ifelse(`Drug A` < `Drug B`,`Drug A`,`Drug B`))
  synscores <- synscores %>% mutate(`Drug B` = ifelse(`Drug A` < `Drug B`,`Drug B`,`drugATMP`))
  synscores <- select(synscores,-c(`drugATMP`,`drugBTMP`))
  synscores <- synscores %>% mutate(`Drug Pair` = paste(`Drug A`,"+",`Drug B`))
  
  # Create a list per experiment ID
  synscoreList = synscores %>% dplyr::group_by(`Experiment ID`,.add=T) %>% group_split()
  
  # Average across all experimentID
  synscoresAverage <- synscores %>% dplyr::group_by(`Drug A`,`Drug B`) %>% dplyr::summarise(`mean_syn` = mean(`Synergy (mean)`,na.rm=T),
                                                                                            `max_syn` = max(`Synergy (mean)`,na.rm = T),
                                                                                            `mad_syn` = mad(`Synergy (mean)`,na.rm=T),
                                                                                            `mean_ant` = mean(`Antagonism (mean)`,na.rm=T),
                                                                                            `max_ant` = max(`Antagonism (mean)`,na.rm=T),
                                                                                            `mad_ant` = mad(`Antagonism (mean)`,na.rm=T),
                                                                                            `mean_int` = mean(`Interaction (mean)`,na.rm=T),
                                                                                            `abs_mean_int` = abs(`mean_int`),
                                                                                            `mad_int` = mad(`Interaction (mean)`,na.rm=T),
                                                                                            `n` = n(), .groups = "drop")
  
    
  
  
  
  #Switching back some of Drug A and Drug B, as it makes p3 nicer
  synscoresAverage <- synscoresAverage %>% mutate(`DrugATMP` = `Drug A`)
  switchIDX = sample(1:nrow(synscoresAverage),floor(nrow(synscoresAverage)/2))
  synscoresAverage$`Drug A`[switchIDX] = synscoresAverage$`Drug B`[switchIDX]
  synscoresAverage$`Drug B`[switchIDX] = synscoresAverage$DrugATMP[switchIDX]
  synscoresAverage <- select(synscoresAverage,-c(`DrugATMP`))
  
  
  # We need to give all experiments an "sd" value for the size in p3a, even the ones with only one observation
  synscoresAverage$mad_syn[is.na(synscoresAverage$mad_syn)] = 0
  synscoresAverage$mad_ant[is.na(synscoresAverage$mad_ant)] = 0
  synscoresAverage$mad_int[is.na(synscoresAverage$mad_int)] = 0
  
  # Labeling top hits
  synscoresAverage <- synscoresAverage %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`mean_syn`))|in_quantile(as.numeric(`mean_ant`)),
                                                        paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
    # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
    mutate(`Quantile` = ifelse(!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
  # Duplicate experimentID and drug pairs are also removed (there shouldn't be any)
  synscoresAverage <- synscoresAverage[!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))),]
  # Adding a colour variable for p1
  synscoresAverage <- synscoresAverage %>% mutate(`colp1` = ifelse(`mean_ant` > `mean_syn`,`max_ant`,-`max_syn`))
  
  
  
  ###########################################################################################
  ####### PLOT 1: Synergy vs. Antagonism across all experiments
  ###########################################################################################
  p1 = ggplot(synscoresAverage, aes(x = `mean_ant`, y = `mean_syn`,fill=`colp1`)) +
    geom_point(color = "gray", shape = 21, aes(size = `mad_int`)) +
    # geom_abline(intercept=0,slope=1, linetype = "dashed", alpha=0.3) +
    scale_size_continuous(range = c(1.5, 4), name = "MAD(\u0394)") +
    scale_fill_gradientn(colours = c("#2166AC","#FFFFFF","#B2182B"),
                           values = scales::rescale(c(min(synscoresAverage$colp1,na.rm=T),max(synscoresAverage$colp1,na.rm=T))),
                           name = expression(paste("max(",Delta^'+',",",Delta^'-',")")),
                           limits = max(abs(synscoresAverage$colp1))*c(-1,1)) +
                           # limits = c(min(0,min(synscoresAverage$colp1,na.rm=T)),max(0,max(synscoresAverage$colp1,na.rm=T)))) +
    geom_text_repel(aes(label = Quantile), color="black",  na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"),
                    segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
    coord_cartesian(xlim = c(), clip = "off") +
    
    xlab("Antagonism") +
    ylab("Synergy") +
    ggtitle("Average interaction across all experiments") +
    guides() +
    theme(plot.title = element_text(hjust = 0.5, size = 12), legend.text = element_text(size = 8),
          legend.title = element_text(size=12),
          legend.position = "right", legend.justification="center", legend.direction = "vertical", legend.box.margin=margin(0,0,0,0),
          text = element_text(size = 8), plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 12), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  
  ###########################################################################################
  ####### PLOT 2: Synergy per drug across all experiments
  ###########################################################################################
  
  
  # Label top hits for each experiment ID
  synscores <- synscores %>% dplyr::group_by(`Experiment ID`) %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy (mean)`),n=1), paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA)))
  
  
  # Plotting only synergy for each experiment ID
  p2 = ggplot(synscores, aes(x = factor(`Experiment ID`), y = `Synergy (mean)`, fill = `Synergy (mean)`)) +
    geom_point(aes(size = `Synergy (sd)`), colour="gray",pch = 21, alpha = 0.8) +
    scale_size_continuous(trans = 'reverse', range = c(.5, 4), name = "Std. dev") +
    scale_fill_gradientn(colours = c("#FFFFFF","#2166AC"),
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
  
  ###########################################################################################
  ####### PLOT 3: Synergy & Antagonism pairwise 
  ###########################################################################################
  
  # Some minor changes here
  
  p3 = ggplot(data = synscoresAverage, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`mean_int`")) +
    geom_point(color = "gray", shape = 21, aes_string(size = "`mad_int`")) +
    scale_size(range = c(2, 5), name = "MAD(\u0394)") +
    guides() +
    # scale_fill_gradient2(low = "#2166AC", mid = "#EAECCC" , high = "#B2182B", space = "Lab", name = "Interaction", limits = c(min(synscoresAverage$mean_int),max(synscoresAverage$mean_int)),midpoint = 0) +
    scale_fill_gradientn(colours = c("#2166AC","#FFFFFF","#B2182B"),
                         values = scales::rescale(c(min(synscoresAverage$mean_int,na.rm=T),max(synscoresAverage$mean_int,na.rm=T))),
                         name = "Interaction",
                         limits = max(abs(synscoresAverage$mean_int))*c(-1,1)) +
    theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
    ggtitle("Average pairwise drug interaction across all experiments") +
    coord_cartesian(clip = 'off') +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
          legend.justification="center", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, vjust = 1, hjust = 1))
  
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
  
  ###########################################################################################
  ####### PLOT 4-6: Synergy vs. Antagonism per experiment ID
  ###########################################################################################
  
  if (groupbyExperimentID){
    # if (length(synscoreList) > 1){
      for (i in 1:length(synscoreList)){
        synSlice = synscoreList[[i]]
        # Labeling top hits
        synSlice <- synSlice %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy (mean)`))|in_quantile(as.numeric(`Antagonism (mean)`)),
                                                            paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
          # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
          mutate(`Quantile` = ifelse(!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
        # Duplicate drug pairs are also removed.
        synSlice <- synSlice[!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))),]
        
        
        # Adding a colour variable for p1
        synSlice <- synSlice %>% mutate(`colp1` = ifelse(`Antagonism (mean)` > `Synergy (mean)`,`Antagonism (mean)`, -`Synergy (mean)`))
        
        
        #Switching back some of Drug A and Drug B, as it makes p6 nicer
        synSlice <- synSlice %>% mutate(`DrugATMP` = `Drug A`)
        switchIDX = sample(1:nrow(synSlice),floor(nrow(synSlice)/2))
        synSlice$`Drug A`[switchIDX] = synSlice$`Drug B`[switchIDX]
        synSlice$`Drug B`[switchIDX] = synSlice$DrugATMP[switchIDX]
        synSlice <- select(synSlice,-c(`DrugATMP`))
        
        ###########################################################################################
        ####### PLOT 4: Synergy vs. Antagonism per experiment ID
        ###########################################################################################
        
        p4 = ggplot(synSlice, aes(x = `Antagonism (mean)`, y = `Synergy (mean)`,fill=`colp1`)) +
          geom_point(color = "gray", shape = 21, aes(size = `Interaction (sd)`)) +
          scale_size_continuous(trans = "reverse",range = c(1.5, 4), name = "Std. dev. (\u0394)") +
          scale_fill_gradientn(colours = c("#2166AC","#FFFFFF","#B2182B"),
                                 values = scales::rescale(c(min(synSlice$colp1,na.rm=T),max(synSlice$colp1,na.rm=T))),
                                 name = expression(paste("max(",Delta^'+',",",Delta^'-',")")),
                                 limits = max(abs(synSlice$colp1))*c(-1,1)) +
          geom_text_repel(aes(label = Quantile), color="black",  na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"),
                          segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
          coord_cartesian(xlim = c(), clip = "off") +
          facet_grid(`Experiment ID`~.) +
          xlab("Antagonism") +
          ylab("Synergy") +
          ggtitle(paste("Synergy vs. Antagonism: ",synSlice$`Experiment ID`[1])) +
          guides() +
          theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),
                legend.title = element_text(size=10),
                legend.position = "right", legend.justification="center", legend.direction = "vertical", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12), plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
        ###########################################################################################
        ####### PLOT 5: Synergy per drug per experiment ID
        ###########################################################################################
        
        
        
        ###########################################################################################
        ####### PLOT 6: Synergy & Antagonism pairwise oer Experiment ID
        ###########################################################################################
        
        p6 = ggplot(data = synSlice, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`colp1`")) +
          geom_point(color = "gray", shape = 21, aes_string(size = "`Interaction (sd)`")) +
          scale_size(trans = "reverse",range = c(2, 5), name = "Std. dev. (\u0394)") +
          guides() +
          scale_fill_gradientn(colours = c("#2166AC","#FFFFFF","#B2182B"),
                               values = scales::rescale(c(min(synSlice$colp1,na.rm=T),max(synSlice$colp1,na.rm=T))),
                               name = expression(paste("max(",Delta^'+',",",Delta^'-',")")),
                               limits = max(abs(synSlice$colp1))*c(-1,1)) +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
          ggtitle(paste("Pairwise drug interaction: ",synSlice$`Experiment ID`[1])) +
          coord_cartesian(clip = 'off') +
          facet_grid(`Experiment ID`~.) +
          theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
                legend.justification="center", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 8, vjust = 1, hjust = 1))
        
        
        
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
