
#rm(list = ls())

{
  pacman::p_load(tidyverse, 
                #holds: dplyr, tidyr, readr, tibble, stringr, purrr, forcats
                tidyr,
                janitor, #added tidying 
                lubridate, # for timeline
                yingtools2, #clinical and microbiome package #https://github.com/ying14/yingtools2
                ggplot2, #necessary plotting
                rstatix, #stats tools
                svglite, # to save svg files
                ggprism, #preferred theme
                patchwork, #add plots together
                cowplot, #add plots together
                RPostgreSQL, ## to access postgres db
                vegan, ##for diversity
                ggpmisc, ##additional plots
                ggview, ## helps preview huge plots
                ggtext, ## ggtext::element_markdown for species italics
                ComplexHeatmap, # heatmap
                gridtext, # adds gt_render and other text tools to grid plots
                ape,# to import and export phylogenetic files
                ggtree,          # to visualize phylogenetic files
                ggtreeExtra,     # extra tools for tree
                tidytree,        # tidy based tree tools
                RPostgreSQL,     # to access postgres db
                treeio,          # to visualize phylogenetic files
                svglite,         # to save as svg
                ggnewscale,     # to add additional layers of color schemes
                data.table, #for last()
                install = FALSE)
  
  library(here) # file management
  library(fs) # file management
  library(phyloseq) # for phyloseq tools
  library(ggplotify) # make plot or grid into ggplot}
  
  
  # paths -------------------
  file_PATH=here()
  fig_PATH=paste0(file_PATH, "/figs")
  fun_PATH=paste0(file_PATH, "/R")
  data_PATH=paste0(file_PATH,"/data")
  
  # load data -----------------
  ## sample redcap and pcr d ata
  samp_data <- list.files(data_PATH, pattern = "_samp_data.RData", full.names = TRUE)
  load(last(samp_data))
  
  ## SS73_samplist -----
  SS73_samplist_raw = fread(paste0(data_PATH,"/SS.TAX.007-SS73_MLST_final.csv")) %>% 
    select(-V1)
  
  SS73_samplist = SS73_samplist_raw %>% 
    mutate(strain_num = str_extract(seqid, "[0-9]\\.fna$")) %>% 
    mutate(strain_num = str_remove(strain_num, "\\.fna$")) %>% 
    mutate(strain = paste0(genoMatchID, "_", strain_num))
}

{
  # need: collection_date, env_broad_scale, env_local_scale, env_medium, geo_loc_name, iso_growth_condt, lat_lon, num_replicons, ref_biomaterial
  # see: https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/
  # see: https://github.com/EnvironmentOntology/envo/wiki/Using-ENVO-with-MIxS
  
  culture_id_ord = c("UCM_01 (ST14)", "UCM_02 (ST83)", "UCM_03 (ST475)", "UCM_04 (ST58)", 
                          "UCM_05 (ST14)", "UCM_06 (ST3)", "UCM_07 (ST3)", "UCM_08 (ST100)", 
                          "UCM_09 (ST2)", "UCM_10 (ST2)", "UCM_11 (ST4)", "UCM_12 (ST28)", 
                          "UCM_13 (ST-)", "UCM_14 (ST15)", "UCM_15 (ST35)", "UCM_16 (ST11)", 
                          "UCM_17 (ST1)", "UCM_18 (ST2)", "UCM_19 (ST42)", "UCM_20 (ST83)", 
                          "UCM_21 (ST-)", "UCM_22 (ST15)", "UCM_23 (ST58)", "UCM_24 (ST37)", 
                          "UCM_25 (ST34)", "UCM_26 (ST3)", "UCM_27 (ST8)", "UCM_28 (ST14)", 
                          "UCM_29 (ST2)", "UCM_30 (ST2)", "UCM_31 (ST1)", "UCM_32 (ST2)", 
                          "UCM_33 (ST26)", "UCM_34 (ST3)", "UCM_35 (ST1)", "UCM_36 (ST58)", 
                          "UCM_37 (ST26)", "UCM_38 (ST58)", "UCM_39 (ST3)", "UCM_40 (ST1)", 
                          "UCM_41 (ST1)", "UCM_42 (ST3)", "UCM_43 (ST15)", "UCM_44 (ST26)", 
                          "UCM_45 (ST1)", "UCM_46 (ST100)", "UCM_47 (ST15)", "UCM_48 (ST39)", 
                          "UCM_49 (ST42)", "UCM_50 (ST205)", "UCM_51 (ST43)", "UCM_52 (ST188)", 
                          "UCM_53 (ST45)", "UCM_54 (ST-)", "UCM_55 (ST15)", "UCM_56 (ST110)", 
                          "UCM_57 (ST2)", "UCM_58 (ST37)", "UCM_59 (ST1)", "UCM_60 (ST54)", 
                          "UCM_61 (ST83)", "UCM_62 (ST1)", "UCM_63 (ST10)", "UCM_64 (ST14)", 
                          "UCM_65 (ST43)", "UCM_66 (ST26)", "UCM_67 (ST34)", "UCM_68 (ST109)", 
                          "UCM_69 (ST42)", "UCM_70 (ST231)", "UCM_71 (ST53)", "UCM_72 (ST42)", 
                          "UCM_73 (ST58)")
  
  ncbi_df <- SS73_samplist %>% 
    left_join(study.redcap.TOTAL %>% select(genoMatchID, date_collect, db)) %>% 
    mutate(culture_id = factor(culture_id,
                               levels = culture_id_ord)) %>% 
    mutate(env_broad_scale = "clinical patient assessment facility [ENVO:03501136]", #clinical patient assessment facility
           env_local_scale = "fecal environment [ENVO:01001029]",
           env_medium = "feces [UBERON:0001988]",
           geo_loc_name = "USA: Chicago, IL",
           iso_growth_c ="Not Applicable",
           lat_lon="41.788200 N 87.604492 W",
           num_replicons = 1,
           ref_biomaterial="Not Applicable")
  
  date_col = ncbi_df %>% select(date_collect)
}

  