##set_ncbi_taxid.R



# Function -------
set_ncbi_taxid <- function(df){
  ## standardize the ncbi taxids to the most recent taxdump 
  # df has these columns seq_id, taxid, min_taxid, min_tax_lvl, species_lab , pctseqs, relative_abundance
  
  ## the taxids will be from the most recent taxdump since downloading locally (2024-06)
  fix_ncbi_key = tibble(
    taxid = c(
      2304686, #motu missing "Oscillospiraceae"
      3043409, #kraken "Pseudomonas sp. FLM 004-28"
      3043410, #kraken"Pseudomonas sp. FLM 11"
      1740163, #kraken"Burkholderia sp. Bp7605"
      2761539 #kraken"Vibrio virus 2019VC1"
    ),
    fixed_id = c(
      216572,  #motu missing "Oscillospiraceae"
      3043408, #kraken "Pseudomonas sp. FLM 004-28"
      3040598, #kraken"Pseudomonas sp. FLM 11"
      1503053, #kraken"Burkholderia sp. Bp7605"
      1480731 #kraken"Vibrio virus 2019VC1"
    )
  )
  
    # grab the necessary data from df
    vars = c("seq_id", "exp_run", "ref_mOTU",
             "taxid", "min_taxid", "min_tax_lvl", "species_lab" , "pctseqs", "relative_abundance")
    a = df  %>% 
      select(any_of(vars)) %>% 
      #select(seq_id, taxid, min_taxid, min_tax_lvl, species_lab , pctseqs, relative_abundance) %>% 
      # note: motu's references are much older so replace 2304686 Hungateiclostridiaceae with 216572 Oscillospiraceae
      left_join(fix_ncbi_key, by = c("min_taxid" = "taxid")) %>% 
      mutate(min_taxid = coalesce(fixed_id, min_taxid)) %>% 
      select(-fixed_id)
    
    # length(unique(a$taxid)) 
    # length(unique(a$min_taxid)) 
    # length(unique(a$seq_id)) 
    
    
    # use taxid for species level only or min_taxid for the LCA 
    b = taxonomizr::getTaxonomy(unique(a$min_taxid),
                                "accessionTaxa.sql",
                                desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus","species")) 
    
    # dim(b) 
    
    # convert the matrix into a data.frame
    c = data.frame(b) %>% 
      mutate(min_taxid = as.numeric(rownames(b))) %>% 
      select(min_taxid, 
             Kingdom = superkingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species) %>% 
      filter(!is.na(min_taxid))
    
    # add the taxonomy info back to the original data
    d = a %>% 
      left_join(c, by = "min_taxid")
    # dim(d) # same as a
    
    # there is a duplicated unclassified row that doesnt mean anything, remove
    e = d %>% 
      # remove the duplicate unclassified row
      mutate(across(
        Kingdom:Species,
        ~ ifelse((taxid == -1 & species_lab == "UNCLASSIFIED"), "UNCLASSIFIED", . )
      )) %>% 
      # if taxid is NA, replace with -1
      mutate(taxid = ifelse(is.na(taxid), -1, taxid)) %>%
      # # if the value is NA or empty for the tax level, replace with UNCLASSIFIED
      mutate(across(
        Kingdom:Species,
        ~ ifelse((is.na(.) | . == "" | . == " "), "UNCLASSIFIED", . )
      )) %>%
      distinct()
    
    
    f = d %>%
      filter(is.na(min_taxid))

  return(e)
} #end of function
