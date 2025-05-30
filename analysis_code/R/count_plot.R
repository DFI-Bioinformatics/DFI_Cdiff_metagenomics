count_plot <- function(count, meta, rank.samps,
                          countcolor= "#264653"){
  # required tax df: has Kingdom:Species level information
  # required meta df: stores the samplenames (MUST have group_var)
  # required rank.samps: sample order determined outside of function
  # required TAXID: numeric that will be found as a ncbi taxid within
  # optional taxcolor: string for hexcode to use as fill of (otherwise drawn from NCBIpal) otherwise set at "#264653"
  
  # df to prepped for the plot
  count.sub <- count %>% 
    ungroup() %>% 
    # filter for meta genoMatchIDs
    filter(genoMatchID %in% meta$genoMatchID) %>%
    select(genoMatchID, count)
  # add taxcolor 
  
  
  if (nrow(count.sub) < 1 | sum(count.sub$count) == 0) {
    count.gg <- count.sub %>% 
      # add metadata even if missing
      full_join(meta,
                by = "genoMatchID",
                relationship = "many-to-many") %>%
      # refactor by ranking variable
      mutate(genoMatchID = factor(genoMatchID,
                                  levels = rank.samps)) %>%
      ggplot() +
      geom_bar(aes(x = genoMatchID,
                   y = 0),
               stat = "identity") +
      theme_prism(border=T)  +
      # label titles
      ylab("Reads") +
      xlab("") +      # facet by meta
      facet_grid(cols = vars(group_var),
                 scales = "free_x",
                 space = "free_x") +
      # theme
      theme_prism(#base_size = 5, 
        border = TRUE )  +
      # adjust plot themes
      theme(
        # keep minimal to adjust later
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        # strip.text.y = element_text(size=, angle=270),
        strip.text.x = element_blank(),
        # axis.text.y = element_text(size=5),
        legend.position = "none"
      )
  } 
  else {
    count.gg <- count.sub %>% 
      # add metadata even if missing
      full_join(meta,
                by = "genoMatchID",
                relationship = "many-to-many") %>%
      #if count == NA then replace with 0
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      # refactor by ranking variable
      mutate(genoMatchID = factor(genoMatchID,
                                  levels = rank.samps)) %>%
      ggplot() +
      # make barplots
      geom_bar(aes(x = genoMatchID,
                   y = count),
               fill = countcolor,
               color = "black",
               stat = "identity",
               width = 0.9) +
      #set y limits on log scale
      scale_y_log10("Reads",
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x)),
                    expand = c(0.00005, 0.00005)) +
      #remove random white space above min and max
      scale_x_discrete(expand = c(0.001, 0.001))+
      # label titles
      labs(x="") +
      # facet by meta
      facet_grid(cols = vars(group_var),
                 scales = "free_x",
                 space = "free_x") +
      # theme
      theme_prism(#base_size = 5, 
        border = TRUE )  +
      # adjust plot themes
      theme(
        # keep minimal to adjust later
        axis.text.x = element_text(angle = 90),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.y = element_text(size=, angle=270),
        strip.text.x = element_blank(),
        axis.text.y = element_text(size=5),
        legend.position = "none"
      )
  }
  return(count.gg)
  
}

