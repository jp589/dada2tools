#'@title Data set from Theis 2019 "Does the Human Placenta..."
#'
#'@description 16s rRNA sequencing data processed using Dada2.
#'Placental samples are included with blank and control samples.
#'Samples normalized using `phyloseq::rarefy_even_depth()`.
#'
#'@format A dataframe with 201 rows and 99 variables:
#'\describe{
#'\item{SRR number}{NCBI Sample number}
#'}
#'@source Files from NCBI website SRA Run Selector
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_rare_ASV"

#'@title Data set from Theis 2019 "Does the Human Placenta..."
#'
#'@description 16s rRNA sequencing data processed using Dada2.
#'Placental samples are included with blank and control samples.
#'
#'@format A dataframe with 1925 rows and 250 variables:
#'\describe{
#'\item{X}{Sequence from 16s rRNA region of Taxa genome}
#'\item{Kingdom}{Kingdom Classification of ASV}
#'\item{Phylum}{Phylum Classification of ASV}
#'\item{Class}{Class Classification of ASV}
#'\item{Order}{Order Classification of ASV}
#'\item{Family}{Family Classification of ASV}
#'\item{Genus}{Genus Classification of ASV}
#'\item{Species}{Species Classification of ASV}
#'\item{SRR number}{NCBI Sample number}
#'}
#'@source Files from NCBI website SRA Run Selector
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_merged"

#'@title Cleaned up Data Set from Theis 2019 "Does the Human Placenta..."
#'
#'@description 16s rRNA sequencing data processed using Dada2.
#'Placental samples (AC and V samples) are included with control samples.
#'
#'@format A dataframe with 967 rows and 106 variables:
#'\describe{
#'\item{X}{Sequence from 16s rRNA region of Taxa genome}
#'\item{Kingdom}{Kingdom Classification of ASV}
#'\item{Phylum}{Phylum Classification of ASV}
#'\item{Class}{Class Classification of ASV}
#'\item{Order}{Order Classification of ASV}
#'\item{Family}{Family Classification of ASV}
#'\item{Genus}{Genus Classification of ASV}
#'\item{Species}{Species Classification of ASV}
#'\item{Sample Name}{Descriptive sample name.}
#'}
#'@source Filtered and cleaned merged data set.
#'
"Theis_merged_clean"

#'@title Sample Metadata for Theis 2019 "Does the Human Placenta..."
#'
#'@description All columns from the sample metadata.
#'
#'@format A dataframe with 242 rows and 4 variables:
#'\describe{
#'\item{Run}{SRR number given to each sample}
#'\item{Host}{Sample origin}
#'\item{Sample.Name}{Entered sample name}
#'}
#'@source NCBI website SRA Run Selector Metadata
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_META"

#'@title Sample Metadata for Theis 2019 "Does the Human Placenta..."
#'
#'@description Select rows of sample metadata for samples run with 30 cycles of PCR after normalization.
#'
#'@format A dataframe with 99 rows and 5 variables:
#'\describe{
#'\item{Run}{SRR number given to each sample}
#'\item{Sample}{Entered sample name}
#'\item{Type}{Sample type}
#'\item{GA}{Gestational age}
#'\item{Delivery}{Placental sample delivery method}
#'}
#'@source NCBI website SRA Run Selector Metadata
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_META_order"

#'@title Sample Metadata for Theis 2019 "Does the Human Placenta..."
#'
#'@description Select rows of sample metadata for samples run with 30 cycles of PCR after normalization.
#'
#'@format A dataframe with 99 rows and 7 variables:
#'\describe{
#'\item{Run}{SRR number given to each sample}
#'\item{Sample}{Entered sample name}
#'\item{Type}{Sample type}
#'\item{GA}{Gestational age}
#'\item{Delivery}{Placental sample delivery method}
#'\item{Del_GA}{Delivery method followed by gestational age}
#'\item{Level}{Placental region sampled}
#'}
#'@source NCBI website SRA Run Selector Metadata
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_rare_META"

#'@title Data Set from Theis 2019 "Does the Human Placenta..." Prepared for Heatmapping
#'
#'@description 16s rRNA sequencing data processed using DADA2tools functions.
#'Genera are on columns and Samples are on rows.
#'
#'@format A data frame with 98 rows and 15 variables:
#'\describe{
#'\item{Achromobacter}{Taxa Genus Classification 1}
#'\item{Delftia}{Taxa Genus Classification 2}
#'\item{Phyllobacterium}{Taxa Genus Classification 3}
#'\item{Stenotrophomonas}{Taxa Genus Classification 4}
#'\item{Blastomonas}{Taxa Genus Classification 5}
#'\item{Clostridium_sensu_stricto_5}{Taxa Genus Classification 6}
#'\item{Methylobacterium-Methylorubrum}{Taxa Genus Classification 7}
#'\item{Methylobacterium-Methylorubrum.1}{Taxa Genus Classification 8}
#'\item{0319-6G20 Unclassified}{Taxa Genus Classification 9}
#'\item{Ralstonia}{Taxa Genus Classification 10}
#'\item{Leucobacter}{Taxa Genus Classification 11}
#'\item{Bacillus}{Taxa Genus Classification 12}
#'\item{Acinetobacter.1}{Taxa Genus Classification 13}
#'\item{Leucobacter.1}{Taxa Genus Classification 14}
#'\item{Cutibacterium}{Taxa Genus Classification 15}
#'}
#'@source Data set prepared for heatmapping.
#'
"Theis_preheat_ab"

#'@title Taxonomy from Theis 2019 "Does the Human Placenta..."
#'
#'@description Taxonomy of 16s rRNA sequencing data processed using Dada2
#'after filtering and normalization.
#'
#'@format A dataframe with 201 rows and 8 variables:
#'\describe{
#'\item{X}{Sequence from 16s rRNA region of Taxa genome}
#'\item{Kingdom}{Kingdom Classification of ASV}
#'\item{Phylum}{Phylum Classification of ASV}
#'\item{Class}{Class Classification of ASV}
#'\item{Order}{Order Classification of ASV}
#'\item{Family}{Family Classification of ASV}
#'\item{Genus}{Genus Classification of ASV}
#'\item{Species}{Species Classification of ASV}
#'}
#'@source Files from NCBI website SRA Run Selector
#'\url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA397876&o=acc_s\%3Aa}
#'
"Theis_rare_TAX"
