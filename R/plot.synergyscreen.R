#' Plot function for a synergyscreen object
#' 
#' @description A function for plotting summary statistics of a \code{synergyscreen} object.
#' 
#' @param x An object of class \code{synergyscreen}, the result of \code{\link{synergyscreen}}.
#' @param groupbyExperimentID logical; if TRUE, individual plots are produced per experiment ID.
#' @param save_plot logical; if TRUE plots are saved locally.
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
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous facet_grid xlab ylab ggtitle labs guides theme unit element_text element_blank margin aes_string guide_legend scale_size scale_fill_gradient2 
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom dplyr mutate %>% group_split n
#' 
#' @export 


plot.synergyscreen <- function(x, groupbyExperimentID = T, save_plot = FALSE, path = NULL, plotdevice = "pdf", ...){
  
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
  
  # Get # of unique drugs
  unqDrugs = length(unique(synscores$`Drug A`))
  
  # Creating a function to label top synergy scores
  in_quantile <- function(x, n = 5) {
    # Pick out to get top and bottom 5
    vec = rep(F,length(x))
    vec[tail(order(x),n)] = T
    return(vec)
  }
  
  # Average across all experimentID
  synscoresAverage <- synscores %>% dplyr::group_by(`Drug A`,`Drug B`) %>% dplyr::summarise(`Synergy Score` = mean(`Synergy Score`),
                                                                              `Antagonism Score` = mean(`Antagonism Score`),
                                                                              `Mean (syn)` = mean(`Mean (syn)`),
                                                                              n = n(), .groups = "drop")
  
  # Labeling top hits
  synscoresAverage <- synscoresAverage %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy Score`))|in_quantile(as.numeric(`Antagonism Score`)),
                                                        paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
    # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
    mutate(`Quantile` = ifelse(!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
  # Duplicate experimentID and drug pairs are also removed.
  synscoresAverage <- synscoresAverage[!duplicated(t(apply(synscoresAverage[c("Drug A", "Drug B")], 1, sort))),]
  
  ###########################################################################################
  ####### PLOT 1: Synergy vs. Antagonism across all experiments
  ###########################################################################################
  p1 = ggplot(synscoresAverage, aes(x = `Antagonism Score`, y = `Synergy Score`)) +
    geom_point(aes(colour = `Drug A`, size = abs(ifelse(abs(`Antagonism Score`) > abs(`Synergy Score`), `Antagonism Score`, `Synergy Score`))), alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), guide = FALSE) +
    geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"), 
                    segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
    xlab("Antagonism") +
    ylab("Synergy") +
    ggtitle("Average interaction across all experiments") +
    guides(colour = guide_legend(title="Drugs", nrow = 4), size = FALSE) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 8),
          legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(0,0,0,0),
          text = element_text(size = 12), plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 12), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  
  ###########################################################################################
  ####### PLOT 2: Synergy per drug across all experiments
  ###########################################################################################
  
  
  # Relabel top hits for all drug pairs again and keep this time both transposed notations
  synscoresAverage <- synscoresAverage %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy Score`),n=10), paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA)))
  
  
  # Plotting only synergy for each drug
  p2 = ggplot(synscoresAverage, aes(x = factor(`Drug A`), y = `Synergy Score`)) +
    geom_point(aes(colour = `Drug B`, size = abs(ifelse(abs(`Antagonism Score`) > abs(`Synergy Score`), `Antagonism Score`, `Synergy Score`))), alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), guide = FALSE) +
    scale_y_continuous(name = "Synergy") +
    geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"), 
                    segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
    xlab("") +
    ylab("Synergy") +
    ggtitle("Average synergy across all experiments") +
    guides(colour = guide_legend(title="Drugs", nrow = 4), size = FALSE) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 8),
          legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
  
  ###########################################################################################
  ####### PLOT 3: Synergy & Antagonism pairwise 
  ###########################################################################################
  
  # Some minor changes here
  
  p3a = ggplot(data = synscoresAverage, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`Synergy Score`")) +
    geom_point(color = "gray", shape = 21, aes_string(size = "`Synergy Score`")) +
    scale_size(range = c(1, 6)) +
    guides(size = F) +
    scale_fill_gradient2(low = "#EAECCC", mid = "#7EB3D4" , high = "#354B99", space = "Lab", name = "Synergy Score", midpoint = 1) +
    theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
    ggtitle("Average pairwise drug synergy across all experiments") +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
          legend.justification="center", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
 
  
  p3b = ggplot(data = synscoresAverage, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`Antagonism Score`")) +
    geom_point(color = "gray", shape = 21, aes_string(size = "`Antagonism Score`")) +
    scale_size(range = c(1, 6)) +
    guides(size = F) +
    scale_fill_gradient2(low = "#EAECCC", mid = "#FBA35D" , high = "#A50026", space = "Lab", name = "Antagonism Score", midpoint = 1) +
    theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
    ggtitle("Average pairwise drug antagonism across all experiments") +
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8), legend.title = element_text(size=10),
          legend.justification="center", legend.box.margin=margin(1,0,0,0),
          text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
          plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))

  
  
  if (!save_plot){
    print(p1)
    readline("Press key for next plot")
    print(p2)
    readline("Press key for next plot")
    print(p3a)
    readline("Press key for next plot")
    print(p3b)
    readline("Press key for next plot")
  } else { # Save plots locally
    # First plot
    file.name = paste0(path,paste("Average_interaction",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
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
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p2) 
    dev.off()
    # Third plot (a)
    file.name = paste0(path,paste("Average_pairwise_synergy",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p3a) 
    dev.off()
    # Third plot (b)
    file.name = paste0(path,paste("Average_pairwise_antagonism",sep = "_"))
    if (plotdevice == "pdf"){
      file.name.pdf = paste0(file.name,".pdf")
      pdf(file=file.name.pdf, ...)
    } else {
      file.name.png = paste0(file.name,".png")
      png(filename = file.name.png, ...)
    }
    print(p3b) 
    dev.off()
  }
  
  
  ###########################################################################################
  ####### PLOT 4-6: Synergy vs. Antagonism per experiment ID
  ###########################################################################################
  
  if (groupbyExperimentID){
    synscoreList = group_split(synscores,`Experiment ID`)
    if (length(synscoreList) > 1){
      for (i in 1:length(synscoreList)){
        synSlice = synscoreList[[i]]
        # Labeling top hits
        synSlice <- synSlice %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy Score`))|in_quantile(as.numeric(`Antagonism Score`)),
                                                            paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA))) %>%
          # Duplicate top hit annotations are removed. All annotations for transposed drug treatments are removed, only the annotation for one drug pair is kept.
          mutate(`Quantile` = ifelse(!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))) == TRUE, `Quantile`, as.character(NA)))
        # Duplicate drug pairs are also removed.
        synSlice <- synSlice[!duplicated(t(apply(synSlice[c("Drug A", "Drug B")], 1, sort))),]
        
        
        ###########################################################################################
        ####### PLOT 4: Synergy vs. Antagonism per experiment ID
        ###########################################################################################
        p4 = ggplot(synSlice, aes(x = `Antagonism Score`, y = `Synergy Score`)) +
          geom_point(aes(colour = `Drug A`, size = abs(ifelse(abs(`Antagonism Score`) > abs(`Synergy Score`), `Antagonism Score`, `Synergy Score`))), alpha = 0.8) +
          scale_size_continuous(range = c(0.2, 5), guide = FALSE) +
          geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"), 
                          segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
          facet_grid(`Experiment ID`~.) +
          xlab("Antagonism") +
          ylab("Synergy") +
          ggtitle(paste("Interaction Scores:",synSlice$`Experiment ID`[1])) +
          guides(colour = guide_legend(title="Drugs", nrow = 4), size = FALSE) +
          theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 10),
                legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12), plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 12), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        
        ###########################################################################################
        ####### PLOT 5: Synergy per drug per experiment ID
        ###########################################################################################
        
        
        # Relabel top hits for all drug pairs again and keep this time both transposed notations
        synSlice <- synSlice %>% mutate(`Quantile` = ifelse(in_quantile(as.numeric(`Synergy Score`),n=10), paste(`Drug A`, "+", `Drug B`, sep = " "), as.character(NA)))
        
        
        # Plotting only synergy for each drug
        p5 = ggplot(synSlice, aes(x = factor(`Drug A`), y = `Synergy Score`)) +
          geom_point(aes(colour = `Drug B`, size = abs(ifelse(abs(`Antagonism Score`) > abs(`Synergy Score`), `Antagonism Score`, `Synergy Score`))), alpha = 0.8) +
          scale_size_continuous(range = c(0.2, 5), guide = FALSE) +
          scale_y_continuous(name = "Synergy") +
          geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.2, "lines"), 
                          segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 3) +
          facet_grid(`Experiment ID`~.) +
          xlab("") +
          ylab("Synergy") +
          ggtitle(paste("Synergy Scores:",synSlice$`Experiment ID`[1])) +
          guides(colour = guide_legend(title="Drugs", nrow = 4), size = FALSE) +
          theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 10),
                legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
        
        
        ###########################################################################################
        ####### PLOT 6: Synergy & Antagonism pairwise oer Experiment ID
        ###########################################################################################
        p6a = ggplot(data = synSlice, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`Synergy Score`")) +
          geom_point(color = "gray", shape = 21, aes_string(size = "`Synergy Score`")) +
          scale_size(range = c(1, 6)) +
          guides(size = F) +
          facet_grid(`Experiment ID`~.) +
          scale_fill_gradient2(low = "#EAECCC", mid = "#7EB3D4" , high = "#354B99", space = "Lab", name = "Synergy Score", midpoint = 1) +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
          ggtitle(paste("Pairwise synergies:",synSlice$`Experiment ID`[1])) +
          theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8),legend.title = element_text(size=10),
                legend.justification="center", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
        
        
        p6b = ggplot(data = synSlice, mapping = aes_string(x = "`Drug B`", y = "`Drug A`", fill = "`Antagonism Score`")) +
          geom_point(color = "gray", shape = 21, aes_string(size = "`Antagonism Score`")) +
          scale_size(range = c(1, 6)) +
          guides(size = F) +
          facet_grid(`Experiment ID`~.) +
          scale_fill_gradient2(low = "#EAECCC", mid = "#FBA35D" , high = "#A50026", space = "Lab", name = "Antagonism Score", midpoint = 1) +
          theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) +
          ggtitle(paste("Pairwise antagonism:",synSlice$`Experiment ID`[1])) +
          theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8), legend.title = element_text(size=10),
                legend.justification="center", legend.box.margin=margin(1,0,0,0),
                text = element_text(size = 12),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))
        
        
        
        if (!save_plot){
          print(p4)
          readline("Press key for next plot")
          print(p5)
          readline("Press key for next plot")
          print(p6a)
          readline("Press key for next plot")
          print(p6b)
          readline("Press key for next plot")
        } else { # Save plots locally
          # First plot
          file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Interaction",sep = "_"))
          if (plotdevice == "pdf"){
            file.name.pdf = paste0(file.name,".pdf")
            pdf(file=file.name.pdf, ...)
          } else {
            file.name.png = paste0(file.name,".png")
            png(filename = file.name.png, ...)
          }
          print(p4) 
          dev.off()
          # Second plot
          file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Synergy",sep = "_"))
          if (plotdevice == "pdf"){
            file.name.pdf = paste0(file.name,".pdf")
            pdf(file=file.name.pdf, ...)
          } else {
            file.name.png = paste0(file.name,".png")
            png(filename = file.name.png, ...)
          }
          print(p5) 
          dev.off()
          # Third plot (a)
          file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Pairwise_synergy",sep = "_"))
          if (plotdevice == "pdf"){
            file.name.pdf = paste0(file.name,".pdf")
            pdf(file=file.name.pdf, ...)
          } else {
            file.name.png = paste0(file.name,".png")
            png(filename = file.name.png, ...)
          }
          print(p6a) 
          dev.off()
          # Third plot (b)
          file.name = paste0(path,paste(synSlice$`Experiment ID`[1],"Pairwise_antagonism",sep = "_"))
          if (plotdevice == "pdf"){
            file.name.pdf = paste0(file.name,".pdf")
            pdf(file=file.name.pdf, ...)
          } else {
            file.name.png = paste0(file.name,".png")
            png(filename = file.name.png, ...)
          }
          print(p6b) 
          dev.off()
        }
      }
    }
    
  }
  

  
  
  
}
