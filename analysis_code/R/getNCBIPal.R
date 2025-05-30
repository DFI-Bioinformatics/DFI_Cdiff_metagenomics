# 2024-03-27
## Sophie Son
# standardize the palette generator for phylo using taxids


## original by Eric Littman
getShades <- function(spdf){
  
  hexcol <- unique(spdf$color)
  
  if (nrow(spdf) > 3){
    resdf <- spdf %>%
      mutate(cols = rep(shades(hexcol, variation = 0.25),
                        length.out = nrow(spdf)
      )
      )
  } else {
    resdf <- spdf %>%
      mutate(cols = shades(hexcol, variation = 0.25, ncolor = nrow(spdf)))
  }
  
  return(resdf)
}


## original from Eric Littman
getRdpPal <- function(tax) {
  
  require(tidyverse)
  require(yingtools2)
  
  tax <- tax %>%
    ungroup()
  
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax.obj)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family",
             "Genus")
  
  if (!all(ranks %in% names(tax))) {
    stop("Error: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  
  tax.dict <- tax %>%
    dplyr::select(all_of(ranks)) %>%
    distinct()
  
  # set all color to gray as base
  tax.dict <- tax.dict %>%
    mutate(color = rep(shades("gray", variation = 0.25),
                       length.out = nrow(tax.dict)))
  
  # color for each level ------------------------------------------------------------
  phypal <- tibble(Phylum = c("Proteobacteria","Actinobacteria","Bacteroidetes"),
                   phycol = c("red","#A77097","#51AB9B"))
  ordpal <- tibble(Order = c("Clostridiales"),
                   ordcol = c("#9C854E"))
  fampal <- tibble(Family = c("Lachnospiraceae","Ruminococcaceae","Oscillospiraceae","Erysipelotrichaceae",
                              "Lactobacillaceae"),
                   famcol = c("#EC9B96","#9AAE73","#9AAE73","orange","#3b51a3"))
  genpal <- tibble(Genus = c("Enterococcus","Streptococcus","Staphylococcus",
                             "Lactobacillus"),
                   gencol = c("#129246","#9FB846","#f1eb25", "#3b51a3"))
  
  tax.split <- tax.dict %>%
    left_join(phypal) %>%
    left_join(ordpal) %>%
    left_join(fampal) %>%
    # ambiguous genus match
    mutate(gencol = case_when(
      grepl("Enterococcus$", Genus) ~ "#129246",
      grepl("Streptococcus$", Genus) ~ "#9FB846",
      grepl("Staphylococcus$", Genus) ~ "#f1eb25",
      TRUE ~ NA_character_
    )) %>%
    # left_join(genpal) %>%
    mutate(color = case_when(
      !is.na(gencol) ~ gencol,
      !is.na(famcol) ~ famcol,
      !is.na(ordcol) ~ ordcol,
      !is.na(phycol) ~ phycol,
      TRUE ~ color)
    ) %>%
    dplyr::select(Kingdom:Genus, color) %>%
    group_split(color)
  
  tax.color <- bind_rows(lapply(tax.split,getShades))
  tax.palette <- structure(tax.color$cols, names = as.character(tax.color$Genus))
  return(tax.palette)
  
}
    
    
    
## Sophie Fix for synonyms in NCBI
getNCBIPal <- function(tax) {
  require(tidyverse)
  require(yingtools2)
  
  tax <- tax %>%
    ungroup()
  
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax.obj)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family",
             "Genus")
  
  if (!all(ranks %in% names(tax))) {
    stop("Error: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  
  # just extract the taxonomy ranks
  tax.dict <- tax %>%
    dplyr::select(all_of(ranks)) %>%
    distinct()
  
  # set all color to gray as base
  tax.dict <- tax.dict %>%
    mutate(color = rep(shades("gray", variation = 0.25),
                       length.out = nrow(tax.dict)))
  
  # color for the top phylum in the gut microbiota ------------------------------------------------------------
  phypal = tibble(
    Phylum = 
      c("Actinomycetota", 
        "Pseudomonadota",
        "Bacteroidota", 
        "Bacillota", 
        "Fusobacteriota", 
        "Verrucomicrobiota",
        # synonyms from NCBI taxonomy browser
        ## "Actinomycetota"
        "Actinobacteria","Actinobacteriota", "Actinobacteraeota", 
        ## "Pseudomonadota"
        "Proteobacteria", "Alphaproteobacteraeota",
        "Proteobacteriota",
        ##"Bacteroidota"
        "Bacteroidetes", "Sphingobacteria", "Bacteroidaeota",
        ## "Bacillota"
        "Firmicutes", "Firmicuteota", "Bacillaeota",
        ##"Fusobacteriota"
        "Fusobacteraeota", "Fusobacteria",
        ## "Verrucomicrobiota"
        "Verrucomicrobia", "Verrucomicrobaeota"),
    phycol = c("#A77097",
               "red", 
               "#51AB9B",
               "#B56B45",
               "#4E8098", 
               "#4A1942",
               ## "Actinomycetota"
               "#A77097","#A77097", "#A77097",
               ## "Pseudomonadota"
               "red","red","red",
               ## "Bacteroidota"
               "#51AB9B","#51AB9B","#51AB9B",
               ## Bacillota
               "#B56B45","#B56B45","#B56B45",
               ## "Fusobacteriota"
               "#4E8098","#4E8098",
               ##"Verrucomicrobiota"
               "#4A1942","#4A1942"))
  
  fampal <- tibble(Family = c("Lachnospiraceae","Ruminococcaceae","Oscillospiraceae","Erysipelotrichaceae",
                              "Lactobacillaceae"),
                   famcol = c("#EC9B96","#9AAE73","#9AAE73","orange","#3b51a3"))
  
  tax.split <- tax.dict %>%
    left_join(phypal) %>%
    #left_join(fampal) %>%
    # ambiguous genus match
    mutate(ordcol = case_when(
      # in the p__Firmicute c__clostridia
      ## these are synonyms
      grepl("Clostridiales", Order) ~ "#9C854E",
      grepl("Eubacteriales", Order) ~ "#9C854E",
      TRUE ~ NA_character_
    )) %>%
    mutate(famcol = case_when(
      grepl("Lachnospiraceae", Family) ~ "#EC9B96",
      grepl("Ruminococcaceae", Family) ~ "#9AAE73",
      grepl("Oscillospiraceae", Family) ~ "#9AAE73",
      grepl("Hungateiclostridiaceae", Family) ~ "#9AAE73",
      grepl("Erysipelotrichaceae", Family) ~ "orange",
      grepl("Lactobacillaceae", Family) ~ "#3b51a3",
      grepl("Leuconostocaceae", Family) ~ "#3b51a3",
      TRUE ~ NA_character_
    )) %>%
    # ambiguous genus match
    mutate(gencol = case_when(
      grepl("Enterococcus$", Genus) ~ "#129246",
      grepl("Streptococcus$", Genus) ~ "#9FB846",
      grepl("Staphylococcus$", Genus) ~ "#f1eb25",
      TRUE ~ NA_character_
    )) %>%
    mutate(color = case_when(
      # if gencol !=NA then color == gencol
      !is.na(gencol) ~ gencol,
      # else if famcol !=NA then color == famcol
      !is.na(famcol) ~ famcol,
      !is.na(ordcol) ~ ordcol,
      !is.na(phycol) ~ phycol,
      TRUE ~ color)
    ) %>%
    dplyr::select(Kingdom:Genus, color) %>%
    group_split(color)
  
  tax.color <- bind_rows(lapply(tax.split,getShades))
  # the Genus col determines the color
  tax.palette <- structure(tax.color$cols, 
                           names = as.character(tax.color$Genus))
  
  return(tax.palette)
  
}



## add print legend for Sophie NCBI


getNCBIPal_colors <- function(tax) {
  
  require(tidyverse)
  require(yingtools2)
  
  tax <- tax %>%
    ungroup()
  
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax.obj)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family",
             "Genus")
  
  if (!all(ranks %in% names(tax))) {
    stop("Error: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  
  tax.dict <- tax %>%
    dplyr::select(all_of(ranks)) %>%
    distinct()
  
  # set all color to gray as base
  tax.dict <- tax.dict %>%
    mutate(color = rep(shades("gray", variation = 0.25),
                       length.out = nrow(tax.dict)))
  
  # color for each level ------------------------------------------------------------
  phypal <- tibble(Phylum = c("Proteobacteria","Actinobacteria","Bacteroidetes"),
                   phycol = c("red","#A77097","#51AB9B"))
  ordpal <- tibble(Order = c("Clostridiales"),
                   ordcol = c("#9C854E"))
  fampal <- tibble(Family = c("Lachnospiraceae","Ruminococcaceae","Oscillospiraceae","Erysipelotrichaceae",
                              "Lactobacillaceae"),
                   famcol = c("#EC9B96","#9AAE73","#9AAE73","orange","#3b51a3"))
  genpal <- tibble(Genus = c("Enterococcus","Streptococcus","Staphylococcus",
                             "Lactobacillus"),
                   gencol = c("#129246","#9FB846","#f1eb25", "#3b51a3"))
  
  tax.split <- tax.dict %>%
    left_join(phypal) %>%
    left_join(ordpal) %>%
    left_join(fampal) %>%
    # ambiguous genus match
    mutate(gencol = case_when(
      grepl("Enterococcus$", Genus) ~ "#129246",
      grepl("Streptococcus$", Genus) ~ "#9FB846",
      grepl("Staphylococcus$", Genus) ~ "#f1eb25",
      TRUE ~ NA_character_
    )) %>%
    # left_join(genpal) %>%
    mutate(color = case_when(
      !is.na(gencol) ~ gencol,
      !is.na(famcol) ~ famcol,
      !is.na(ordcol) ~ ordcol,
      !is.na(phycol) ~ phycol,
      TRUE ~ color)
    ) %>%
    dplyr::select(Kingdom:Genus, color) %>%
    group_split(color)
  
  tax.color <- bind_rows(lapply(tax.split,getShades))
  
  return(tax.color)
}

