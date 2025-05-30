## Sophie Style Guide
## style guide for Sophie based projects


## preferred plotting theme
#library(ggprism)

## example
# ggplot(data, aes(x=x, y=y))+
#   geom_bar(stat="identity")+
#   #remove random white space above min and max
#   scale_y_discrete(expand = c(0.1, 0))+
#   #remove random white space above min and max
#   scale_x_discrete(expand = c(0.1, 0))+
#   theme_prism(
#     # borders
#     border = TRUE
#   )+
#   theme(
#     #aspect ratio to make square plots
#     aspect.ratio = 1,
#     # edit axis title text size
#     axis.title = element_text(size = 10),
#     # edit axis text size
#     axis.text = element_text(size = 8)
#   )

## functions ------
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

## taxplot colors -------
## for taxplots: see getNCBIPal.R
#source("~/Library/CloudStorage/Box-Box/SS-RR_CDiff-Clinical-Proj/code-share/getNCBIPal.R")

# ranking phyla: 
phylo_colors = c(
    # Genus
    "Enterococcus"="#129246", 
    "Streptococcus" ="#9FB846", 
    "Staphylococcus"="#f1eb25",
    #family
    "Lactobacillaceae"="#3b51a3",    
    "Oscillospiraceae"="#9AAE73", 
    "Erysipelotrichaceae"="orange", 
    "Lachnospiraceae"="#EC9B96",
    # order
    "Eubacteriales"="#9C854E",
    "Actinomycetota"="#A77097",
    "Bacteroidota"="#51AB9B",
    "Fusobacteriota"="#4E8098",
    "Pseudomonadota"="red",
    "Bacillota"="#B56B45",
    "Verrucomicrobiota"="#4A1942",

    "Other"= "gray"
)



## tax classifier order -----------
# plot in order: metaphlan4, motu, kraken


## tax classifiers --------
tax_class_cols = c(
  # colors checked for color blindness
  "Kraken2" = "#DC7474",
  "mOTU3" = "#779BD5",
  "Metaphlan4_default" = "#417447",
  
  # pcr
  "tcdB" = "#c86d8e",
  "Cdiff_16S"="#053359",
  
  # used for statQ
  "Metaphlan4_statq_0.1" = "#7BCC97",
  "Metaphlan4_statq_0.05" = "#8BADAD", #"#A7DDBA",
  "Metaphlan4_local_statq_0.05" = "#A7C786",
  # used for db Vers
  "Metaphlan4_V1" = "#ffd58e",
  "Metaphlan4_V2" = "#7352A9",
  "Metaphlan4_V3" = "#52b3b6",
  
  ## usually not needed
  "Kraken2_undetected" = "#FAC8C8",
  "Metaphlan4_default_undetected" = "#D9EAD3",
  "mOTU3_undetected" = "#CFE2F3"
)


## Cdiff MLST --------
cdiff_clade_cols = c(
  "Clade 1" = "#8ad3db",
  "Clade 2" = "#2bac13",
  "Clade 3" = "#bfb6e2",
  "Clade 4" = "#ff54d8",
  "Clade 5" =  "#f24949",
  "-" = "gray"
)

library(yingtools2)
library(here)
data_PATH=paste0(here(),"/data")

MLST_lookup <- read.csv(paste0(data_PATH, "/typing/pubMLST_clade_ST_lookup.csv")) 

MLST_col_assign <- MLST_lookup %>%
  mutate(ST = as.character(ST)) %>% 
  
  arrange(mlst_clade) %>% 
  mutate(color = case_when(
    mlst_clade == 1 ~ cdiff_clade_cols[[1]],
    mlst_clade == 2 ~ cdiff_clade_cols[[2]],
    mlst_clade == 3 ~ cdiff_clade_cols[[3]],
    mlst_clade == 4 ~ cdiff_clade_cols[[4]],
    mlst_clade == 5 ~ cdiff_clade_cols[[5]],
    is.na(mlst_clade) ~ cdiff_clade_cols[[6]]
  )) %>% 
  select(mlst_clade, ST, color) %>%
  group_split(color)
  

ST_col <- bind_rows(lapply(MLST_col_assign,getShades)) %>% 
  mutate(MLST = paste0("ST", ST))
ST_pal <- structure(ST_col$cols, names = ST_col$MLST)

## Classifiation results colors ------
class_result_cols = c(
  "TP" = "#006d77",
  "TN" = "#5e548e",
  "FP" = "#dc5b71",
  "FN" = "#f3722c",
  "target_val" = "#FF0000",
  "control_group" = "#234760"
)


## genome source colors --------

genome_cols = c(
  CDC30 = "#136cfb", # cdc color
  SS73 = "#800000", # uchicago color
  SS108 = "#caaaad", # likely not used
  DFI150 = "#245d38", # likey not used
  Knight260 = "#8c1b53"  # likely not used
)

## study cohort colors --------
db_cols = c(
  HT = "#e25272",
  LT = "#9bbb59",
  LD = 	"#ff9c59",
  MICU = "#80c5ca",
  `Multi-study` =  "#9d85a9",
  HD = "#343d46"
)

## plot_settings -----
nature <- c(
  "font" = "Helvetica",
  "fig_text_size" = 5, #pts
  "fig_lab_size" = 8,
  "table_text_size" = 7,
  "line_size" = 0.5, #pts
  # Figure elements
  "res" = 300,  # for dpi
  "pg_width" = 180, # mm
  "pg_height" = 247, # mm
  "pg_units" = "mm",
  "output_tabletype" = ".csv",
  "output_figtype" = ".svg",
  "label_type" = "auto"
)

# typical 16:9 ratio
ppt <- c(
  "font" = "Helvetica",
  "fig_text_size" = 15, #pts
  "fig_lab_size" = 18,
  "table_text_size" = 18,
  "line_size" = 0.5, #pts
  # Figure elements
  "res" = 300,  # for dpi
  "pg_width" = 13.33, 
  "pg_height" = 7.5, 
  "pg_units" = "in",
  "output_tabletype" = ".csv",
  "output_figtype" = ".png",
  "label_type" = "auto"
)

## backup colors ------
## not used in other groups OR used for the least used pallete option
extra_cols = c(
  # gray
  "#8a959b", "#a7adba", "#4f5b66", "#65737e", "#4f5b66",
  # blue
  "#9dafb8","#8BADAD", "#8aa5bc", "#6391b6", "#6aa2ab", "#38596e", "#278aa7", "#6a99ab", 
  # teal
   "#A7BEAE", "#7dceb3", "#20948B", "#4ebdae", 
  # greens
  "#a4a68e", "#737759", "#2f7e4b", "#6aab97",  "#535231", "#1a472a",
  # yellow/orange
  "#feebd4", "#f0d088", "#ffe599", "#febfa0","#ffa472",  "#c65333", 
  # reds/pinks
  "#caaaad", "#f0c5b4",	"#f49a9c", "#d08288",	"#9f5f5d", "#e1706d", "#f26249",	"#712727", "#ec2f3b", "#990011",
  # purples
  "#c08daa", "#bfb6e2", "#8c1b53"
)
