## taxid_RA_plot.R
## returns ggplot (bar plot) of taxid relative abundance

require(scales)

taxid_RA_plot <- function(tax, meta, rank.samps, TAXID, 
                          taxcolor= "#264653"){
  # required tax df: has Kingdom:Species level information
  # required meta df: stores the samplenames (MUST have group_var)
  # required rank.samps: sample order determined outside of function
  # required TAXID: numeric that will be found as a ncbi taxid within
  # optional taxcolor: string for hexcode to use as fill of (otherwise drawn from NCBIpal) otherwise set at "#264653"
  
  # df to prepped for the plot
  ra.sub <- tax %>% 
    ungroup() %>% 
    # filter for meta genoMatchIDs
    filter(genoMatchID %in% meta$genoMatchID) %>%
    # filter for when the min_taxid matches TAXID
    filter(min_taxid == TAXID) %>% 
    group_by(genoMatchID, min_taxid) %>%
    mutate(tax.pct = sum(pctseqs)) %>% 
    ungroup()
    # add taxcolor 
    
    
  if (nrow(ra.sub) < 1 | sum(ra.sub$tax.pct) == 0) {
    ra.gg <- ra.sub %>% 
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
      ylab("taxidRA") +
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
        # axis.ticks.x = element_blank(),
        # strip.text.y = element_text(size=, angle=270),
        # strip.text.x = element_blank(),
        # axis.text.y = element_text(size=5),
        legend.position = "none"
      )
  } 
    else {
    ra.gg <- ra.sub %>% 
      # add metadata even if missing
      full_join(meta,
                by = "genoMatchID",
                relationship = "many-to-many") %>%
      #if gen.rel.abund == NA then replace with 0
      mutate(tax.pct = ifelse(is.na(tax.pct), 0, tax.pct)) %>%
      # refactor by ranking variable
      mutate(genoMatchID = factor(genoMatchID,
                                  levels = rank.samps)) %>%
      ggplot() +
      # make barplots
      geom_bar(aes(x = genoMatchID,
                   y = tax.pct),
               fill = taxcolor,
               color = "black",
               stat = "identity",
               width = 0.9) +
      #set y limits 0-100%
      # coord_cartesian(ylim = c(0, 1)) +
      #remove random white space above min and max
      scale_y_continuous(
        # make percent (%) out of 100
        limits = c(0, max(ra.sub$tax.pct) + 0.00005),
        breaks = seq(0, max(ra.sub$tax.pct)+ 0.00005, length.out = 4),
        labels = label_percent(),
        expand = c(0.00005, 0.00005)) +
      #remove random white space above min and max
      scale_x_discrete(expand = c(0.001, 0.001))+
      # label titles
      ylab("taxidRA") +
      xlab("") +
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
        # axis.ticks.x = element_blank(),
        # strip.text.y = element_text(size=, angle=270),
        # strip.text.x = element_blank(),
        # axis.text.y = element_text(size=5),
        legend.position = "none"
      )
  }
  return(ra.gg)
  
}

# test------
# data for plotting (MP)
# tax = read.table("~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/syncom_abundance_results/SS.TAX.011-clinsamp/SS.TAX.011-mp_default_tax_df_2024-06-11.tsv")
# 
# 
# # genus level colors
# 
# meta = read.csv("~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/filter_samp_list/samp.list.full_2024-05-09.csv") %>%
#   select(-X) %>%
#   filter(grepl("HD_DFI", genoMatchID)) %>% 
#   mutate(group_var = "Healthy Donors")
# 
# 
# rank.samps = HD_meta %>% 
#   select(genoMatchID) %>% 
#   mutate(HD_num = as.numeric(str_remove(genoMatchID, "HD_DFI_"))) %>% 
#   arrange(HD_num) %>% 
#   pull(genoMatchID)

# taxid_RA_plot(tax, meta, rank.samps, TAXID = 562)
