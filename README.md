# Optimized MetaPhlAn4 and StrainPhlAn

## Run the analysis
To run the analysis as described in the publication, download the analysis_code folder and data folder from here: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15557191.svg)](https://doi.org/10.5281/zenodo.15557191)

Make sure you have compatible versions of RStudio and R installed.
Open DFI_Cdiff_metagenomics.Rproj using RStudio and execute code as required or knit to generate the html with all figures.

The analysis code was compiled in RStudio (2025.05.0+496) and R(4.5.0 (2025-04-11))
The versions of packages used are in R_package_versions.txt.


## MetaPhlAn4 - Local StatQ Implementation

1. Install this version of MetaPhlAn:
(https://github.com/biobakery/MetaPhlAn/tree/local_stat_q)

2. Modify `stat_q.cfg` to your choice of species and associated `stat_q` value (default 0.20, lower value corresponds to higher sensitivity)
e.g. `s__Clostridioides_difficile = 0.05`


## MetaPhlAn4 - Optimize markers database for your strains

1. Use extract_markers.py to extract markers for a given SGB of your choice from the metaphlan database pkl file.
```
mkdir -p db_markers
extract_markers.py -c t__SGB16955 -o db_markers/ -d /path/to/metaphlan_db/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl`
```

2. Use the extracted markers to check alignment/presence with your set of strains (using BLAST, Bowtie2 etc.)

3. Create a text file with markers to exclude (that have low mapping score to your set of strains from your analysis)
Say exclude_markers.txt, with the names of markers to be excluded:
```
UniRef90_E6DA28|7__10|SGB16955
UniRef90_U7M4A5|1__8|SGB16955
UniRef90_D4HFG6|3__8|SGB16955
UniRef90_Q6A6M1|5__8|SGB16955
UniRef90_A0A1N4YMG9|1__6|SGB16955
UniRef90_A0A2B7JWE2|1__4|SGB16955
```

4. Run `metaphlan` with `--ignore_markers exclude_markers.txt`


## StrainPhlAn Optimized for C. diff

1. Download latest version of MetaPhlAn4

2. Follow the instructions given in strainphlan_tuned.sh.


# Citation
Coming soon.
