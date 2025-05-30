## genoMatchID_functions.R
## 2024-01-31


## Due to the messiness of all the sampleIDs, make a function that can convert an input column of characters into a genoMatchID

## add_genoMatchID() --------
## add genoMatchID and count infomation
## genoMatchID: 10 character total (XX_###_###) and healthy donors will be (HD_DFI_###)
#### used to prevent changing ID codes shared between other datasets

## output = vector that is pipe friendly


add_genoMatchID <- function(df, #pipefriendly
                            ID.col = NULL #character naming the ID column
                            ){
    df = df

    ## store the old ids in df
    old.ID.df = df %>%
      select(any_of(ID.col)) %>%
      dplyr::rename(old.ID = any_of(ID.col))
    
    
    new.id.df = old.ID.df %>%
      #genoMatchID is typically drawn from the sampleid
      mutate(genoMatchID = old.ID) %>%
      mutate(genoMatchID = gsub("\\.", "_", genoMatchID),
             genoMatchID = gsub("\\-", "_", genoMatchID)) %>%
      # Make MICU_ into MI_ samples
      mutate(genoMatchID = gsub("MICU_", "MI_", genoMatchID)) %>%
      ## one sample is listed as MIUC instead of MICU
      mutate(genoMatchID = gsub("MIUC_", "MI_", genoMatchID)) %>%
      # Make Healthy donors match this notation: 
      ## HD_DFI_### 
      mutate(genoMatchID = gsub("DFI0", "HD_DFI_", genoMatchID)) %>%
      # remove any weird stuff after the last number in sampleid
      mutate(genoMatchID = str_remove(genoMatchID, "[^0-9]*$")) %>%
      # remove any weird stuff BEFORE the first letter of the ID
      mutate(genoMatchID = sub("^_", "", genoMatchID)) %>%
      # make sure the number of charac in genoMatchID is 10 to find outliers
      mutate(num_char = nchar(genoMatchID)) %>%
      # if num_char == 8, then replace "_" with "_0" (for other samps) if HD then make sure its HD_DFI_00
      mutate(genoMatchID = ifelse(num_char==8 & grepl("HD_DFI", genoMatchID), 
                                  gsub("HD_DFI_", "HD_DFI_00", genoMatchID), 
                                  ifelse(num_char==8 & !grepl("HD_DFI", genoMatchID), gsub("_", "_0", genoMatchID), genoMatchID))) %>%
      # if num_char == 9, then replace the last "_" with "_0"
      mutate(genoMatchID = ifelse(num_char==9, 
                                  stringi::stri_replace_last_fixed(genoMatchID, '_', '_0'), 
                                  genoMatchID)) %>%
      mutate(genoMatchID = str_remove_all(genoMatchID, " ")) %>%
      mutate(num_char = nchar(genoMatchID))
    
    # vector to store new ids
    new.ids = new.id.df %>%
      pull(genoMatchID)
  
 return(new.ids) 
}
## END add_genoMatchID()


## is_genoMatchID() -------
## check if the ID col is a genoMatchID format

## returns a df that lists the tested ID as test.ID and if it was in genoMatchID notation

is_genoMatchID <- function(df,
                           ID.col){
  # if ID.col is missing, set ID.col as "genoMatchID
  if(missing(ID.col)){
    ID.col = "genoMatchID"
  }
  
  ## check if genoMatchID format is achieved
  id.df = df %>%
    select(any_of(ID.col))%>%
    dplyr::rename(test.ID = any_of(ID.col))
  
  test.df = id.df %>%
    mutate(is_geno = grepl("^[A-Z]{2}_[A-Z0-9]{3}_[0-9]{3}$",
                           id.df$test.ID)) %>%
    distinct(test.ID, is_geno)

  return(test.df)
  
}
## END is_genoMatchID()


## healthy_mrns() -------
## generate fake mrns for healthy donors for easier 
### must be unique and must be reproducible for each person
## extract the numeric component of the DFI name and add 1000000


healthy_mrns<- function(df,
                        ID.col){
  df = df
  
  # if ID.col is missing, set ID.col as "genoMatchID
  if(missing(ID.col)){
    ID.col = "genoMatchID"
  }
  
  ## check if genoMatchID format is achieved
  old.df = df %>%
    dplyr::rename(ID.samp = any_of(ID.col))
  
  
  fake.mrn = old.df %>%
    mutate(mrn.new = ifelse(
    is.na(mrn) & grepl("DFI", ID.samp),
    as.character(1000000 + as.numeric(substring(ID.samp, nchar(ID.samp)-2, nchar(ID.samp)))),
      NA)) %>%
    distinct(ID.samp, mrn.new, mrn)
  
  ## return the same 
  mrn.df = old.df %>%
    left_join(fake.mrn %>%
                select(ID.samp, mrn.new)) %>%
    mutate(mrn = ifelse(!is.na(mrn.new), mrn.new, mrn)) %>%
    select(-mrn.new)
  
  return(mrn.df)
  
}





# ## test zone -------

#a = "HT_51_01"
#b = "LD_61_01"
#c = "DFI005"

# temp = data.frame(test = c(a,b,c))
# 
# test_df = temp %>%
#   mutate(genoMatchID = add_genoMatchID(temp, ID.col = "test"))

# test.red = postgres_redcap %>%
#   head()
# ID.col = "sampleid"
# 
# test.dfi = postgres_redcap %>%
#   filter(grepl("DFI", ID)) %>%
#   distinct(ID)
# ID.col = "ID"
# 
# 
# v = add_genoMatchID(test.dfi, ID.col = "ID")
# 
# test.b = test.red %>%
#   mutate(genoMatchID = add_genoMatchID(test.red, ID.col = "sampleid"))
# 
# test.c = test.dfi %>%
#   mutate(genoMatchID = add_genoMatchID(test.dfi, ID.col = "ID"))
# 
# 
# pattern = "^[A-Z]{2,}_[A-Z0-9]{3,}_[0-9]{3,}$"
# 
# grepl("^[A-Z]{2,}_[A-Z0-9]{3,}_[0-9]{3,}$", 
#       test.b$ID)
# 
# is_genoMatchID(test.b)
# is_genoMatchID(test.c, 
#                ID.col = "ID")
# 
# test.dfi.red = postgres_redcap %>%
#   filter(grepl("DFI", ID)) %>%
#   distinct(ID, .keep_all = T)
# 
# z = healthy_mrns(test.dfi.red, 
#                  ID.col = "ID")
# 
# h = healthy_mrns(test.red, 
#                      ID.col = "sampleid")




