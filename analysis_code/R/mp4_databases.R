# mp4_databases.R
# 2025-02-25
# parse through the mp4 database lists

db_PATH = "~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/cdiff_tax_proj/data/metaphlan/mp4 databases/"

june2024_mg = fread(paste0(db_PATH, "mpa_vJun23_CHOCOPhlAnSGB_202403_marker_info.txt"),
                      header= FALSE,
                    fill=TRUE,
                    col.names = c(
                      "mg", "V2", "meta", 
                      "Phylum", "Class", "Order", "Family", "Genus", "Species",
                      "sgb"
                    )) %>% 
  mutate(
    sgb = str_remove(sgb,"}$"),
    SGB = str_extract(sgb, "SGB[0-9]+"))

june2024_spec = fread(paste0(db_PATH, "mpa_vJun23_CHOCOPhlAnSGB_202403_species.txt"),
                      header= FALSE,
                      col.names = c("SGB", "tax")) 


SGB_key = c("Clostridioides spp. 1"= "SGB109606",
            "Clostridioides spp. 2" ="SGB109607",
            "Clostridioides difficile" = "SGB6136")

SGB_key = june2024_spec %>% 
  filter(grepl("g__Clostridioides",tax)) %>% 
  pull(SGB)
  

c_spec = june2024_spec %>% 
  filter(SGB %in% SGB_key) %>% 
  separate_wider_delim(tax, delim = ",", 
                       names_sep = "_",
                       too_few = "align_end") %>% 
  pivot_longer(starts_with("tax"),
              names_to = "column",
              values_to = "species_entry") %>% 
  filter(!is.na(species_entry)) %>% 
  select(-column) %>% 
  separate(
    species_entry,
    sep = "\\|",
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ) 


#write.csv(c_spec, paste0(db_PATH, "mpa_vJun23_clostridioides.csv"))

c_spec <- read.csv(paste0(db_PATH, "mpa_vJun23_clostridioides.csv"))

c_spec %>% 
  group_by(SGB) %>% 
  summarize(n = n_distinct(Species))

c_spec %>% 
  group_by(SGB) %>% 
  summarize(n = n_distinct(Species))


c_mg <- june2024_mg %>% 
  filter(SG %in% SGB_key)


#write.csv(c_mg, paste0(db_PATH, "mpa_vJun23_clostridioides_markergenes.csv"))

c_mg %>% 
  group_by(SGB) %>% 
  summarize(n_markergenes = n_distinct(mg))


strp_c_genomes = fread("~/Library/CloudStorage/Box-Box/SS-tax-classifier-stuff/cdiff_tax_proj/data/strainphlan/cdiff_like_genomes.tsv")
