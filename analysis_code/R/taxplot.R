## taxplot function

# colors for taxonomy 
source("~/Library/CloudStorage/Box-Box/SS-RR_CDiff-Clinical-Proj/code-share/getNCBIPal.R")



taxplot  <- function(tax, meta, rank.samps){
  # required tax df: has Kingdom:Species level information
  # required meta df: stores the samplenames (MUST have group_var)
  ## required rank.samps: sample order determined outside of function
  
    # df prepped for the plot
    tax.sub <- tax %>%
      ungroup() %>% 
      # filter for meta genoMatchIDs
      filter(genoMatchID %in% meta$genoMatchID) %>%
      # make sure the the min tax level is S (will only be bacteria)
      filter(min_tax_lvl == "S") %>%
      # things that are not bacteria will end up in the other category
      mutate(genLab = Genus,
             #if there isnt a genLab column, there will be even if this is redundant
             Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus)) %>% 
      group_by(genoMatchID, Kingdom, Phylum, Class, Order, Family, Genus, genLab) %>%
      #gathers all instances of a specific Genus in each sample
      summarize(gen.rel.adbund = sum(pctseqs)) %>% 
      ungroup() %>%
      arrange(Kingdom, Phylum, Class, Order, Family, Genus, genLab) %>%
      #mutate(Genus = factor(genLab, levels = unique(genLab))) %>%
      mutate(genLab = factor(genLab, levels = unique(genLab))) %>%
      group_by(genoMatchID) %>%
      #arrange(Genus) %>%
      arrange(genLab) %>% 
      mutate(cum.pct = cumsum(gen.rel.adbund),
             y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% #sets the plotting y factor for labels
      ungroup() %>%
      dplyr::select(-cum.pct) 
    
    
    colorpal = getNCBIPal(tax.sub)
    
    # color rank for phyla: 
    tax.order = tax.sub %>% 
      # go small to large since case_when is sequential
      mutate(
        phylo_rank = case_when(
          # genus
          grepl("Enterococcus", Genus) ~ "A",  
          grepl("Streptococcus", Genus) ~ "B", 
          grepl("Staphylococcus", Genus) ~ "C",
          #family
          grepl("Lactobacillaceae", Family) ~ "D", 
          grepl("Leuconostocaceae", Family) ~ "D", 
          grepl("Oscillospiraceae", Family) ~ "E", 
          grepl("Hungateiclostridiaceae", Family) ~ "E", 
          grepl("Ruminococcaceae", Family) ~ "E", 
          grepl("Erysipelotrichaceae", Family) ~ "F", 
          grepl("Lachnospiraceae", Family) ~ "F", 
          # order
          grepl("Eubacteriales", Order) ~ "H", 
          grepl("Clostridiales", Order) ~ "H",
          
          # phylum
          grepl("Actinomycetota", Phylum) ~ "I",
          grepl("Actinobacteria", Phylum) ~ "I",
          grepl("Actinobacteriota", Phylum) ~ "I",
          grepl("Actinobacteraeota", Phylum) ~ "I",
          
          grepl("Bacteroidota", Phylum) ~ "J",
          grepl("Bacteroidetes", Phylum) ~ "J",
          grepl("Sphingobacteria", Phylum) ~ "J",
          grepl("Bacteroidaeota", Phylum) ~ "J",
          
          grepl("Fusobacteriota", Phylum) ~ "K",
          grepl("Fusobacteraeota", Phylum) ~ "K",
          grepl("Fusobacteria", Phylum) ~ "K",
          
          grepl("Pseudomonadota", Phylum) ~ "L", 
          grepl( "Proteobacteria", Phylum) ~ "L", 
          grepl("Alphaproteobacteraeota", Phylum) ~ "L", 
          grepl("Proteobacteriota", Phylum) ~ "L", 
          
          grepl("Bacillota", Phylum) ~ "M",
          grepl("Firmicutes", Phylum) ~ "M",
          grepl("Firmicuteota", Phylum) ~ "M",
          grepl("Bacillaeota", Phylum) ~ "M",
          
          grepl("Verrucomicrobiota", Phylum) ~ "N",
          grepl("Verrucomicrobia", Phylum) ~ "N",
          grepl("Verrucomicrobaeota", Phylum) ~ "N",
          
          # Other
          .default = "O"
        )) %>% 
      mutate(phylo_rank = factor(phylo_rank,
                                 levels = c("A", "B", "C", "D", "E", "F",
                                            "G", "H", "I", "J",
                                            "K", "L","N", "M", 
                                            "O"))) %>% 
      arrange( phylo_rank #,Genus
              ) %>% 
      distinct(Genus) %>% 
      pull()
    
    
    tax.gg <- tax.sub %>%
      # add metadata even if missing
      full_join(meta,
                by = "genoMatchID",
                relationship = "many-to-many") %>%
      #if gen.rel.abund == NA then replace with 0
      mutate(gen.rel.adbund = ifelse(is.na(gen.rel.adbund), 0, gen.rel.adbund)) %>%
      # refactor by ranking variable
      mutate(genoMatchID = factor(genoMatchID, levels = rank.samps)) %>%
      # rerank the Genus colors
      mutate(Genus = factor(Genus, levels = tax.order)) %>% 
      ggplot() +
      # make barplots
      geom_bar(aes(x = genoMatchID,
                   y = gen.rel.adbund,
                   fill = Genus),
               stat = "identity",
               position = "fill", 
               width = 0.9) +
      #set y limits 0-100%
      # coord_cartesian(ylim = c(0, 1)) +
      #remove random white space above min and max
      scale_y_continuous(
        limits = c(0 , 1),
        expand = c(0.005, 0.005)) +
      #remove random white space above min and max
      scale_x_discrete(expand = c(0.001, 0.001)) +
      # colors from appropriate taxpal
      scale_fill_manual(values = colorpal) +
      # label titles
      ylab("Rel. Abundance") +
      xlab("") +
      # facet by meta
      facet_grid(cols = vars(group_var),
                 scales = "free_x",
                 space = "free_x") +
      # theme
      theme_prism(#base_size = 5, 
                  border = TRUE)  +
      # adjust plot themes
      theme(
        # keep minimal to adjust later
        axis.text.x=element_text(angle = 90),
        # axis.ticks.x = element_blank(),
        # strip.text.y=element_text(size=, angle=270),
        # strip.text.x = element_blank(),
        # axis.text.y=element_text(size=5),
        legend.position="none"
        )
    
    return(tax.gg)
}

# test------
# data for plotting (MP)
# tax = read.table("~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/syncom_abundance_results/SS.TAX.011-clinsamp/SS.TAX.011-mp_default_tax_df_2024-06-11.tsv")
# 
# 
# # genus level colors
# 
# meta = read.csv("~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/filter_samp_list/samp.list.full_2024-05-09.csv") %>% select(-X)
# 
# 

# rank.samps = HD_meta %>% 
#   select(genoMatchID) %>% 
#   mutate(HD_num = as.numeric(str_remove(genoMatchID, "HD_DFI_"))) %>% 
#   arrange(HD_num) %>% 
#   pull(genoMatchID)

# taxplot(tax, meta, rank.samps)

#meta = HD_meta
#tax = krak_default_tax_df
#rank.samps = HD_rank