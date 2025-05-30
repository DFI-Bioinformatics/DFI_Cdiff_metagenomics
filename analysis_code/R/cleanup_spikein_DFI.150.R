
## cleanup_spikein_DFI.150
## 2024-12-08
## adapt SS.TAX.011 code to grab all of the spikeins 

rm(list = ls())

{
  #load ------------------------------------------------------------
  
  ## packages -------
  {
    pacman::p_load(
      janitor, #added tidying 
      lubridate, # for timeline
      yingtools2, #clinical and microbiome package #https://github.com/ying14/yingtools2
      ggplot2, #necessary plotting
      rstatix, #stats tools
      svglite, # to save svg files
      ggprism, #preferred theme
      patchwork, #add plots together
      RPostgreSQL, ## to access postgres db
      vegan, ##for diversity
      ggpmisc, ##additional plots
      ComplexHeatmap, # heatmaps
      phyloseq, #phyloseq stuff
      data.table, #for last()
      tidyr,
      tidyverse, #holds: dplyr, tidyr, readr, tibble, stringr, purrr, forcats
      install = FALSE)
    
    library(here) # file management
    library(fs) # file management
    
    ## Paths -------
    file_PATH=here()
    fig_PATH=paste0(file_PATH, "/figs")
    fun_PATH=paste0(file_PATH, "/R")
    data_PATH=paste0(file_PATH, "/data")
    tool_PATH="~/Library/CloudStorage/Box-Box/SS-RR_CDiff-Clinical-Proj/tax_classification_results"
    
   setwd("~")
  }
  
  ## meta ------
  {
    ## cdc30 metadata
    cdc_meta_PATH <- list.files(data_PATH, pattern="CDC30_assembly_metadata", full.names = TRUE)
    cdc_meta <- readxl::read_xlsx(cdc_meta_PATH)
}
 
  # Functions -------
  # set dir to the location of the ncbi tax logs
  {
    source(paste0(fun_PATH, "/set_ncbi_taxid.R"))
    
    ## from Ram 
    readin_metaphlan = function(path)
    {
      files = list.files(path, pattern = ".*.txt", recursive = T, full.names = T)
      df_mp_raw = lapply(files, function(file){
        df = fread(file)  %>% 
          mutate(filename = file)
        return(df)
      })
      
      df_mp_raw = do.call("rbind", df_mp_raw)
      return(df_mp_raw)
    }
    
    ## read n clean functions 
    ## just for the default setting of stat_q
      ## metaphlan ----------
    clean_spikein_mp<- function(path)
      {
          files = list.files(path, pattern = ".*.txt", recursive = T, full.names = T)
          files <- files[grepl("cdc_isolates_150", files, ignore.case = TRUE)|
                           grepl("syn-control_statq_0.20_150", files, ignore.case = TRUE)]
          
          ### load mp -----
          mp_raw = lapply(files, function(file){
            df = fread(file)  %>% 
              mutate(filename = file)
          })
          mp_raw = do.call("rbind", mp_raw)
          
          
          ### filter and clean -----
        mp_spikein_df <- mp_raw %>% 
          mutate(seq_id = basename(filename),
                 exp_run = basename(dirname(dirname(filename)))) %>% 
          mutate(seq_id = str_remove(seq_id, "_R1\\.txt"),
                 split_id = seq_id) %>% 
          separate(split_id,
                   sep = "-",
                   into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
          mutate(
            bg_community = str_replace(bg_community, "syn", "DFI.150"),
            n_reads = as.numeric(n_reads),
             rep = as.numeric(rep)) %>%      
          filter(grepl("s__|UNCLASSIFIED", `#clade_name`),
                 !grepl("t__", `#clade_name`)) %>%  # filter if species level is unclassified due to duplicates
          mutate(pctseqs = relative_abundance/100) %>% #pctseqs is from 0-1 and relative_abundance is from 0-100
          mutate(Species = gsub(".*s__", "", `#clade_name`)) %>%
          mutate(`#clade_name` = gsub("[kpcofgst]__","",`#clade_name`)) %>% 
          separate(`#clade_name`,
                   into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","SGB"),
                   sep = "\\|") %>%
          # UNCLASSIFIED APPLIED ACROSS tax levels
          mutate(across(
            Phylum:Species,
            ~ ifelse(NCBI_tax_id == "-1", "UNCLASSIFIED", .))
          ) %>% 
          # if the taxid == -1 then the Species is UNCLASSIFIED
          mutate(Species = ifelse(NCBI_tax_id == "-1", "UNCLASSIFIED", Species))%>%
          # split taxid into pieces
          separate(col = `NCBI_tax_id`, 
                   sep = "\\|",
                   into = c("taxid.K", "taxid.P", "taxid.C", "taxid.O", "taxid.F", "taxid.G", "taxid.S", "taxid.SGB")) %>% 
          # if the value is NA then the taxid is -1
          mutate(across(
            "taxid.P":"taxid.S", 
            ~ ifelse(is.na(.), "-1", .))
          )  %>% 
          # get the smallest (S < G < F < O < C < P < K)
          mutate(min_taxid = ifelse(!is.na(taxid.S) & taxid.S != "", taxid.S,
                                    ifelse(!is.na(taxid.G) & taxid.G != "", taxid.G,
                                           ifelse(!is.na(taxid.F) & taxid.F != "", taxid.F,
                                                  ifelse(!is.na(taxid.O) & taxid.O != "", taxid.O,
                                                         ifelse(!is.na(taxid.C) & taxid.C != "", taxid.C,
                                                                ifelse(!is.na(taxid.P) & taxid.P != "", taxid.P, taxid.K)
                                                         ))))),
                 min_tax_lvl = ifelse(!is.na(taxid.S) & taxid.S != "", "S",
                                      ifelse(!is.na(taxid.G) & taxid.G != "", "G",
                                             ifelse(!is.na(taxid.F) & taxid.F != "", "F",
                                                    ifelse(!is.na(taxid.O) & taxid.O != "", "O",
                                                           ifelse(!is.na(taxid.C) & taxid.C != "", "C",
                                                                  ifelse(!is.na(taxid.P) & taxid.P != "", "P", "K")
                                                           )))))
          ) %>% 
          mutate(
            #taxid = as.numeric(min_taxid)
            taxid = as.numeric(taxid.S), ## version for species level taxid,
            min_taxid = as.numeric(min_taxid)
          ) %>% 
          dplyr::select(seq_id,  bg_community,assembly_accession, n_reads, rep, exp_run,
                        taxid, 
                        min_taxid, min_tax_lvl, 
                        Kingdom, Phylum, Class, Order, Family, Genus, species_lab = Species, pctseqs, relative_abundance,
                        filename) %>%  
          replace_na(list(n_reads=0)) %>% 
            mutate(
              bg_community = str_replace(bg_community, "dfi|syn", "DFI.150"),
              assembly_accession = str_replace(assembly_accession, "syn", "control"),
              seq_id = paste2(
                bg_community, assembly_accession, n_reads, rep,
                sep="-"
              )) %>% 
            # pivot to fill for n_reads such that the syn_control has a value for each
            # pivot_wider(
            #   names_from = n_reads,
            #   values_from = taxid
            # ) %>% 
            # pivot_longer(
            #   c("1000", "5000", "10000", "50000"),
            #   names_to = "n_reads",
            #   values_to = "taxid"
            # ) %>% 
            # mutate(taxid = ifelse(is.na(taxid), `NA`, taxid),
            #        n_reads = as.numeric(n_reads),
            #        bg_community = str_replace(bg_community, "dfi", "DFI.150"),
            #        assembly_accession = str_replace(assembly_accession, "syn", "control"),
            #        seq_id = paste2(
            #          bg_community, assembly_accession, n_reads, rep,
            #          sep="-")) %>% 
            # select(-"NA") %>% 
            # filter(!is.na(taxid)) %>% 
          filter(relative_abundance > 0 ) %>% 
          ungroup() 
          
          return(mp_spikein_df)
        }
        
    
      ## motu3 -------
    clean_spikein_motu <-function(path, ref_path)
      {
          files = list.files(path, pattern = ".*.txt", recursive = T, full.names = T)
          
          files <- files[grepl("cdc_isolates_150", files, ignore.case = TRUE)|
                                grepl("dfi-syn-control", files, ignore.case = TRUE)]

          ### load motu ----

          # uses gtdb for ref 
          {
            motu_ref <- fread(ref_path)  %>%
              # separate columns into taxid col and phylo col at first space only
              separate(col = kingdom, into = c("taxid.K", "Kingdom"), sep =" ", extra = "merge") %>%
              separate(col = phylum, into = c("taxid.P", "Phylum"), sep =" ", extra = "merge") %>%
              separate(col = class, into = c("taxid.C", "Class"), sep =" ", extra = "merge") %>%
              separate(col = order, into = c("taxid.O", "Order"), sep =" ", extra = "merge") %>%
              separate(col = family, into = c("taxid.F", "Family"), sep =" ", extra = "merge") %>%
              separate(col = genus, into = c("taxid.G", "Genus"), sep =" ", extra = "merge") %>%
              separate(col = mOTU, into = c("taxid.S", "Species"), sep =" ", extra = "merge") %>%
              #filter(taxid.K == "2") %>%  ##only bacteria
              # get the smallest (S < G < F < O < C < P < K)
              mutate(min_taxid = ifelse(taxid.S != "NA", taxid.S,
                                        ifelse(taxid.G != "NA", taxid.G,
                                               ifelse(taxid.F != "NA", taxid.F,
                                                      ifelse(taxid.O != "NA", taxid.O,
                                                             ifelse(taxid.C != "NA", taxid.C,
                                                                    ifelse(taxid.P != "NA", taxid.P, taxid.K)
                                                             ))))),
                     min_tax_lvl = ifelse(taxid.S != "NA", "S",
                                          ifelse(taxid.G != "NA", "G",
                                                 ifelse(taxid.F != "NA", "F",
                                                        ifelse(taxid.O != "NA", "O",
                                                               ifelse(taxid.C != "NA", "C",
                                                                      ifelse(taxid.P != "NA", "P", "K")
                                                               )))))
              ) %>%
              mutate(min_taxid = as.numeric(min_taxid))
            
            # grab the min_taxid
            motu_ref_min = motu_ref %>%
              select("ref-mOTU_v2_ID",
                     taxid = taxid.S, min_taxid, min_tax_lvl,
                     # taxid = taxid.S, ## version for Species level taxid only
                     "Kingdom","Phylum","Class","Order","Family","Genus","Species") %>%
              dplyr::rename("ref_mOTU" = "ref-mOTU_v2_ID") %>%
              mutate(taxid = as.numeric(taxid)) %>%
              # add row for UNCLASSIFIED
              add_row(
                taxid = -1,
                min_taxid = -1,
                min_tax_lvl = "S",
                Kingdom = "UNCLASSIFIED",
                Phylum = "UNCLASSIFIED",
                Class = "UNCLASSIFIED",
                Order = "UNCLASSIFIED",
                Family = "UNCLASSIFIED",
                Genus = "UNCLASSIFIED",
                Species = "UNCLASSIFIED",
              )
            
            remove(motu_ref)
          }
          
          motu_raw = lapply(files, function(file){
            df = fread(file)  %>% 
              mutate(filename = file)
          })
          motu_raw = do.call("rbind", motu_raw)
          
          ### filter and clean ------
          motu_spikein_df <- motu_raw %>% 
            mutate(seq_id = basename(filename),
                   exp_run = basename(dirname(dirname(filename)))) %>% 
            mutate(seq_id = str_remove(seq_id, "_motus.txt"),
                   split_id = seq_id) %>% 
            separate(split_id,
                     sep = "-",
                     into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
            mutate(
              bg_community = str_replace(bg_community, "syn", "DFI.150"),
              assembly_accession = str_replace(assembly_accession, "syn", "control"),
              n_reads = as.numeric(n_reads),
              rep = as.numeric(rep)) %>%      
            separate(`#consensus_taxonomy`,
                     sep = " \\[",
                     into = c("species_lab", "ref_mOTU")
            ) %>% 
            # remove remaining 
            mutate(ref_mOTU = str_remove(ref_mOTU, "\\]"),
                   species_lab = str_replace(species_lab, "unassigned", "UNCLASSIFIED")) %>%
            left_join(motu_ref_min, by = "ref_mOTU") %>% 
            # if taxid is NA then replace with -1
            replace_na(list(taxid = -1)) %>%
            # if min_taxid is NA then replace with -1 and 
            dplyr::select(
              seq_id, bg_community,assembly_accession, n_reads, rep, exp_run,
              taxid,
              min_taxid,
              min_tax_lvl, ref_mOTU, species_lab,
              pctseqs = `unnamed sample`,
              filename) %>%  
            replace_na(list(n_reads=0)) %>% 
            mutate(relative_abundance = pctseqs * 100) %>% 
            filter(relative_abundance > 0 ,
                   seq_id = paste2(
                     bg_community, assembly_accession, n_reads, rep,
                     sep="-")) %>%             
            mutate(
                       bg_community = str_replace(bg_community, "dfi|syn", "DFI.150"),
                       assembly_accession = str_replace(assembly_accession, "syn", "control"),
                       seq_id = paste2(
                         bg_community, assembly_accession, n_reads, rep,
                         sep="-"
                       )) %>% 
            ungroup()

          return(motu_spikein_df)
        }
  
    
    ## kraken2/bracken -------
    clean_spikein_krak <- function(path)
      {
      files = list.files(path, pattern = ".*.txt", recursive = T, full.names = T)
      
      files = files[grepl("cdc_isolates_150", files, ignore.case = TRUE)|
                       grepl("dfi-syn-control", files, ignore.case = TRUE)]
      
      ### load kraken ----
      krak_raw = lapply(files, function(file){
        df = fread(file)  %>% 
          mutate(filename = file)
      })
      
      krak_raw = do.call("rbind", krak_raw)
     
      krak_df =  krak_raw  %>% 
        mutate(seq_id = basename(filename),
               exp_run = basename(dirname(dirname(filename)))) %>%
        # only the species level
        filter(grepl("_bracken_species",seq_id)) %>%
        mutate(seq_id = str_remove(seq_id, "_bracken_species\\.txt"),
               split_id = seq_id)%>% 
        separate(split_id,
                 sep = "-",
                 into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
        mutate(
          bg_community = str_replace(bg_community, "syn", "DFI.150"),
          n_reads = as.numeric(n_reads),
          rep = as.numeric(rep)) %>%      
        filter(grepl("cdc_isolates_150", exp_run, ignore.case = TRUE)|
                 grepl("dfi-syn-control", seq_id, ignore.case = TRUE))%>% 
        replace_na(list(n_reads=0)) %>% 
        # pivot to fill for n_reads such that the syn_control has a value for each
        # pivot_wider(
        #   names_from = n_reads,
        #   values_from = taxonomy_id
        # ) %>% 
        # pivot_longer(
        #   c("1000", "5000", "10000", "50000"),
        #   names_to = "n_reads",
        #   values_to = "taxonomy_id"
        # ) %>% 
        # mutate(taxonomy_id = ifelse(is.na(taxonomy_id), `NA`, taxonomy_id),
        #        n_reads = as.numeric(n_reads)) %>% 
        # select(-"NA") %>% 
        #filter(!is.na(taxonomy_id)) %>% 
        # total reads of DFI.150=150*5e4
        # add the seq_id to allow for joining
        mutate(total_reads = (150*5e4)+n_reads,
               min_taxid = taxonomy_id,
               min_tax_lvl = "S"
        ) %>% 
        #select cols
        dplyr::select(seq_id, bg_community, assembly_accession, n_reads, rep, exp_run,
                      species_lab = name, taxid = taxonomy_id, min_taxid, min_tax_lvl,
                      reads = new_est_reads, total_reads,
                      filename) %>% 
        group_by(seq_id, exp_run, n_reads) %>%
        # add up reads from output
        mutate(sum_reads = sum(reads)) %>%
        # make sure taxid is numeric
        mutate(taxid = as.numeric(taxid),
               reads = as.numeric(reads),
               total_reads = as.numeric(total_reads),
               sum_reads = as.numeric(sum_reads)) %>% 
        ungroup()
      
      # check number for krak_df (missing 11-12 samples per n_read)
      # u = krak_spikein_tax_df %>%
      #   distinct(seq_id,bg_community,assembly_accession,n_reads,rep,exp_run,filename) %>%
      #   group_by( exp_run, n_reads, rep) %>%
      #   summarize(n = n_distinct(assembly_accession))
      
      ## calculate the pctseqs and relative abundance for unclassified
      grab_ids = krak_df %>% 
        select(seq_id,bg_community,assembly_accession,n_reads,rep,exp_run,
               total_reads, sum_reads, filename) %>% 
        mutate(calc_unclassed = total_reads - sum_reads) %>% 
        distinct()
      
      unclass_df = grab_ids %>% 
        mutate(
          taxid = -1,
          min_taxid = -1,
          min_tax_lvl = "S",
          species_lab = "UNCLASSIFIED",
          reads = calc_unclassed) %>% 
        dplyr::select(seq_id, bg_community, assembly_accession, n_reads, rep, exp_run,
                      taxid, min_taxid, min_tax_lvl, species_lab , reads, total_reads,
                      filename)
      
      krak_spikein_df = bind_rows(krak_df, unclass_df) %>% 
        # grouped by seq_id
        group_by(seq_id, exp_run,n_reads)%>%
        # calculate pctseqs 
        mutate(pctseqs = reads/total_reads) %>%
        # calculate relative_abundance and check if totalpctseqs==1
        mutate(relative_abundance = pctseqs*100,
               # check that this adds to 1
               tot_pctseqs = sum(pctseqs),
               bg_community = str_replace(bg_community, "dfi", "DFI.150"),
               assembly_accession = str_replace(assembly_accession, "syn", "control"),
               seq_id = paste2(
                 bg_community, assembly_accession, n_reads, rep,
                 sep="-"
               )) %>%
        #remove sum_reads = total classified reads
        select(seq_id, bg_community, assembly_accession, n_reads, rep, exp_run,
               taxid, min_taxid, min_tax_lvl, species_lab , pctseqs, relative_abundance, filename) %>%
        # remove if RA is 0
        filter(relative_abundance > 0 ) %>% 
        ungroup()%>% 
        distinct()
      # the length should be original length + number of ids
      
      return(krak_spikein_df)
    }
    
  }



  # read and clean -----
  {
  ## run paths ------
  krak_PATH=paste0(tool_PATH, "/bracken")
  
  mp_PATH=paste0(tool_PATH,"/metaphlan")
  
  motu_PATH=paste0(tool_PATH,"/motus")
  motu_ref_PATH=paste0(data_PATH,"/db_mOTU3_taxonomy_all_mOTU.txt")
  motu_gtdb_PATH=paste0(data_PATH,"/mOTUs_3.0.0_GTDB_tax.tsv")
  
  ## run -----
  krak_spikein_df <- clean_spikein_krak(krak_PATH)
  mp_spikein_df <-clean_spikein_mp(mp_PATH)
  motu_spikein_df <- clean_spikein_motu(motu_PATH, motu_ref_PATH)
  ## motu takes the longest by far
}
  
  # standardize taxid -----
  {
  ### standardize to ncbi taxid dump -----
  ## the taxids will be from the most recent taxdump since downloading locally (2024-06)
  krak_spikein_tax_df  = set_ncbi_taxid(krak_spikein_df) %>% 
    mutate(split_id = seq_id)%>% 
    separate(split_id,
             sep = "-",
             into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
    mutate(n_reads=as.numeric(n_reads),
           rep = as.numeric(rep)) %>% 
    ungroup()
  #saveRDS(krak_spikein_tax_df , "_DFI.150_CDC.30_krak_spikein_tax_df.rds")
  
  krak_spikein_sum = krak_spikein_tax_df %>%
    distinct(seq_id,bg_community,assembly_accession,n_reads,rep,exp_run) %>%
    group_by( exp_run, n_reads, rep) %>%
    summarize(n = n_distinct(assembly_accession))
  #saveRDS(krak_spikein_sum, "_DFI.150_CDC.30_krak_spikein_sum.rds")

  # the taxids will be from the most recent taxdump since downloading locally (2024-06)
  motu_spikein_tax_df  = set_ncbi_taxid(motu_spikein_df) %>% 
    mutate(split_id = seq_id)%>% 
    separate(split_id,
             sep = "-",
             into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
    mutate(n_reads=as.numeric(n_reads),
           rep = as.numeric(rep)) %>% 
    ungroup()
  #saveRDS(motu_spikein_tax_df, "_DFI.150_CDC.30_motu_spikein_tax_df.rds")
  
  motu_spikein_sum <- motu_spikein_tax_df %>%
    distinct(seq_id,bg_community,assembly_accession,n_reads,rep,exp_run) %>%
    group_by( exp_run, n_reads, rep) %>%
    summarize(n = n_distinct(assembly_accession)) 
  #saveRDS(motu_spikein_sum, "_DFI.150_CDC.30_motu_spikein_sum.rds")
  
  ## the taxids will be from the most recent taxdump since downloading locally (2024-06)
  mp_spikein_tax_df  = set_ncbi_taxid(mp_spikein_df) %>% 
    mutate(split_id = seq_id)%>% 
    separate(split_id,
             sep = "-",
             into = c("bg_community","assembly_accession", "n_reads", "rep")) %>% 
    mutate(n_reads=as.numeric(n_reads),
           rep = as.numeric(rep)) %>% 
    ungroup() 
  #saveRDS(mp_spikein_tax_df, "_DFI.150_CDC.30_mp_spikein_tax_df.rds")
  
  mp_spikein_sum <- mp_spikein_tax_df %>%
    distinct(seq_id,bg_community,assembly_accession,n_reads,rep,exp_run) %>%
    group_by( exp_run, n_reads, rep) %>%
    summarize(n = n_distinct(assembly_accession)) 
  #saveRDS(mp_spikein_sum, "_DFI.150_CDC.30_mp_spikein_sum.rds")
  
  
  }

# save -----------
  {
  ## LISTS ------
  syn_spikein_LIST = list(
    # mp4 cdc30 spikeins into synthetic communities (includes statq iterations)
    mp_spikein_tax_df,
    # motu cdc30 spikeins into synthetic communities
    motu_spikein_tax_df,
    # bracken cdc30 spikeins into synthetic communities
    krak_spikein_tax_df)
  names(syn_spikein_LIST)=c("mp_spikein_tax_df","motu_spikein_tax_df","krak_spikein_tax_df")
  
  syn_spikein_SUMMARY = list(
    # mp4 cdc30 spikeins into synthetic communities (includes statq iterations)
    mp_spikein_sum,
    # motu cdc30 spikeins into synthetic communities
    motu_spikein_sum,
    # bracken cdc30 spikeins into synthetic communities
    krak_spikein_sum
  )
  names(syn_spikein_SUMMARY) = c("mp_spikein_sum", "motu_spikein_sum","krak_spikein_sum")

  

  setwd(data_PATH)
  saveVars <- c(
    "syn_spikein_LIST",
    "syn_spikein_SUMMARY",
    "cdc_meta"
  )

  curVars <- ls() # remove list

  rm(list = curVars[!curVars %in% saveVars])

  # Save
  save.image(paste0(as.character(lubridate::today()),
                    "_DFI.150_CDC.30_spikeins.RData"))
  }
}  
