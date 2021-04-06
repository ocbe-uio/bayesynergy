##### Script for pre-processing the dataset from [1]
# Used for generating data for Figure 5 and Table 1

# Input: 
#         156849_1_supp_0_w2lh45.xlsx (single agent viabilities)
#         156849_1_supp_1_w2lrww.xls  (combination viabilities) 
#         N.B. This file needs to be opened and resaved as xlsx
#         Both of these available as supplementary material for [1]

# Output:
#         OCUBM : Gemcitabine + MK-8776.csv
#         .csv file containing the single example experiment used for Figure 5 and Table 1

# set path to the location of these files

#path = "path to raw input files"

# Use dplyr and tidyr for preprocessing, readxl for input
library(tidyr)
library(dplyr)
library(readxl)

# Read in files
mono = read_excel(paste0(path,"/156849_1_supp_0_w2lh45.xlsx"))
combo = read_excel(paste0(path,"/156849_1_supp_1_w2lrww.xlsx"))

# Select only OCUBM cell line and the drugs c("Gemcitabine","MK-8776")
drugs =  c("Gemcitabine","MK-8776")
cellLine = "OCUBM"

mono = mono %>% filter(`cell_line` == cellLine) %>% filter(`drug_name` %in% drugs)
combo = combo %>% filter(`cell_line` == cellLine) %>% filter(`drugA_name` == drugs[1] & `drugB_name` == drugs[2] |
                                                             `drugA_name` == drugs[2] & `drugB_name` == drugs[1])

# Now work everything into the correct format, viabilities are numeric
mono = mono %>% mutate(across(.cols=starts_with("viability"),as.numeric))
combo = combo %>% mutate(across(.cols=starts_with("viability"),as.numeric))

# Convert to long format and remove some columns
mono  = pivot_longer(mono,cols=starts_with("viability"),values_to = "viability") %>% 
  select(cell_line, drug_name, `Drug_concentration (µM)`,viability)
combo = pivot_longer(combo, cols=starts_with("viability"),values_to="viability") %>% 
  select("cell_line", "drugA_name", "drugA Conc (µM)", "drugB_name","drugB Conc (µM)","viability")

# Now get mono into the same format as the combo
mono = mono %>% mutate("drugA_name" = combo$drugA_name[1],
                "drugB_name" = combo$drugB_name[1]) %>%
  mutate(`drugA Conc (µM)` = ifelse(drug_name == drugA_name,`Drug_concentration (µM)`,0)) %>%
  mutate(`drugB Conc (µM)` = ifelse(drug_name == drugB_name,`Drug_concentration (µM)`,0)) %>%
  select("cell_line", "drugA_name", "drugA Conc (µM)", "drugB_name","drugB Conc (µM)","viability")
# Combine the two sets
final = bind_rows(mono,combo)


# Write csv
final %>% write.csv(paste0(path,"/OCUBM : Gemcitabine + MK-8776.csv"))


# References:
# [1]
# Jennifer O'Neil, Yair Benita, Igor Feldman, Melissa Chenard, Brian Roberts, Yaping Liu, Jing Li, Astrid Kral, Serguei Lejnine, Andrey Loboda, William Arthur, Razvan Cristescu, Brian B. Haines, Christopher Winter, Theresa Zhang, Andrew Bloecher and Stuart D. Shumway
# Mol Cancer Ther June 1 2016 (15) (6) 1155-1162; DOI: 10.1158/1535-7163.MCT-15-0843