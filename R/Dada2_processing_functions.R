
#' Preparation for Barplotting
#'
#' Given a dataframe with genera as columns and samples as rows, this will produce a taxa by sample type dataframe.
#'
#' @param df dataframe with genera as columns and samples as rows.
#' @param sample_types named list where names are sample types in order, and values are number of samples per type.
#'
#' @return Produces a taxa by sample type dataframe.
#' @export
#'
#' @examples
#'
#' types_list <- list("Amnion-Chorion" = 28, "Control" = 41, "Villous Tree" = 29)
#' prebar <- barplot_prep(preheat_ab, sample_types = types_list)
barplot_prep <- function(df, sample_types) {
  #data frame input is transposed and converted back into a dataframe.
  df1 <- as.data.frame(t(df))

  #loops for every sample type
  for (i in 1:length(sample_types)){
    #adds a new column to the dataframe which represents the row sum for the columns of each sample type.
    df1[,1 +NCOL(df1)] <- rowSums(df1[,1:sample_types[[i]]])
    #subsets the dataframe to only those columns newly added.
    df1 <- df1[,(sample_types[[i]]+1):NCOL(df1)]
  }
  #uses the names of the sample types for column names.
  colnames(df1) <- names(sample_types)
  return(df1)
}


#' Barplotting
#'
#' Function used to create stacked barplots of absolute and relative reads of sample types.
#'
#' @param df dataframe prepared with barplot_prep()
#' @param study name of the study or origin of the data
#' @param rt_mar adjustable right margin parameter. For when taxa names exceed space. Default is 10.
#'
#' @return returns a list of barplots.
#' @export
#'
#' @examples
#'
#' types_list <- list("Amnion-Chorion" = 28, "Control" = 41, "Villous Tree" = 29)
#' prebar <- barplot_prep(preheat_ab, sample_types = types_list)
#' barplotting(df = prebar, study = "Example", rt_mar = 18)
barplotting <- function(df, study, rt_mar = 10) {
  #ensures dataset is a dataframe if not already.
  df <- as.data.frame(df)

  #creates the local function getPalette which uses the colorRampPalette function of the RColorBrewer package and the Spectral palette which is 11 colors.
  getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

  #getPalette() is used to extrapolate the Spectral palette to any number of steps beyond the original 11.
  colr <- getPalette(NROW(df))

  #sets the graphical parameters xpd to NA which allows plotting outside the normal bounds, and increases the right margin by 10 by default.

  graphics::par(xpd = NA, mar = c(5.1,4.1,4.1, rt_mar))

  #plots the absolute reads per sample type as a stacked barplot.
  plot_abs <- graphics::barplot(as.matrix(df), col = colr, border = "white", ylab = "Absolute Reads", main = paste(study, "Samples by Type"))

  #includes a legend to the right of the plot. x and y coordinates can be adjusted as necessary to fit properly.
  graphics::legend(x = (NCOL(df) + 1), y = max(colSums(df)), legend = rownames(df[order(nrow(df):1),]), fill = rev(colr), bty = "n")

  #transforms average absolute abundance read counts to percentages of sample types.
  df.prop <- as.data.frame(prop.table(as.matrix(df),2)); df.prop <- df.prop*100

  #plots the taxa percentages of sample types
  plot_prop <- graphics::barplot(as.matrix(df.prop), col = colr, border = "white", ylab = "Avg Percentage of Sample", main = paste(study, "Samples by Type"))

  #creates a legend of taxa and their colors and places it to the right of the plot.
  graphics::legend(x = (NCOL(df.prop) + 1), y = 100, legend = rownames(df.prop[order(nrow(df.prop):1),]), fill = rev(colr), bty = "n")

  #resets the plotting margins back to original settings.
  graphics::par(mar = c(5.1,4.1,4.1,4.1))

  #returns a list of the two plots.
  plots <- list(Abs_Abundance = plot_abs, Proportional = plot_prop)
  return(plots)
}


#' PCoA plotting
#'
#' Create a PCoA plot based on Bray-Curtis dissimilarity to map samples and visualize similarity.
#'
#' @param df a dataset that has been processed by dadaset_clean() and col_reorder()
#' @param sample_types a named list of sample types with counts per sample types as values
#' @param study name of study from which dataset originated
#' @param inset alters the placement of the legend. -0.4 is the default.
#' @param rt_margin adjusts the size of the right margin if needed to fit the legend.
#' @param invert_x flips plot over the y-axis.
#' @param invert_y flips plot over the x-axis.
#'
#' @return returns a Bray-Curtis PCoA plot of sample types with a legend on the right.
#' @export
#'
#' @examples
#' types_list <- list("Amnion-Chorion" = 28, "Control" = 41, "Villous Tree" = 29)
#' beta_div(merged_clean, sample_types = types_list, study = "Example Study")
beta_div <- function(df, sample_types, study, inset = -0.4, rt_margin = 9.1, invert_x = FALSE, invert_y = FALSE) {

  #converts dataframe to datatable.
  dt <- data.table::as.data.table(df)

  #subsets data table to only columns that are numeric.
  dt_numeric <- dt[,sapply(dt, is.numeric), with = FALSE]

  #transposes the dataset so that species are on columns and samples are on rows.
  df_numeric_t <- as.data.frame(t(dt_numeric))

  #instantiates X
  X <- c()

  #creates a character vector X which contains defined length repeats of each sample type name based on the value the sample type.
  for (i in 1:length(sample_types)) {
    X <- append(X, values = rep(names(sample_types)[i],sample_types[[i]]))
  }
  #creates a Type column and fills it with values from X
  df_numeric_t$Type <- X

  #uses the vegdist() function to calculate bray-curtis dissimilarity indices
  df.bray <- vegan::vegdist(df_numeric_t[1:(NCOL(df_numeric_t)-1)])

  #transforms B-C indices to PCoA coordinates
  df.pcoa <- ecodist::pco(df.bray)

  #options for inverting over the X or Y axis.
  ifelse(invert_x == TRUE, invert_x <- -1, invert_x <- 1)
  ifelse(invert_y == TRUE, invert_y <- -1, invert_y <- 1)

  #sets the plot margins with adjustable right margin.
  graphics::par(xpd = NA, mar=c(5.1, 4.1, 4.1, rt_margin))

  #plots the PCoA. pch parameter is used to change the shape of point used(i.e. star instead of dot)
  graphics::plot(x = invert_x*(df.pcoa$vectors[,1]), y = invert_y*(df.pcoa$vectors[,2]), col= as.numeric(as.factor(df_numeric_t$Type)), pch=as.numeric(as.factor(df_numeric_t$Type)), main = paste(study, "by Sample Type"), xlab = "PCO1", ylab = "PCO2")
  #adds a legend in the topright corner of the plot with adjustable inset.
  graphics::legend("topright", inset = c(inset,0), legend = sort(unique(df_numeric_t$Type)), col = 1:length(unique(df_numeric_t$Type)), pch = 19, bty = "n")
  #resets the plotting settings to original.
  graphics::par(xpd = FALSE, mar=c(5.1, 4.1, 4.1, 4.1))
}


#' Column Reorder
#'
#' Reorders numeric columns alphabetically based on sample names.
#'
#' @param dt input which is a dada2 taxa table fused to the ASV table.
#'
#' @return Reordered data table
#' @export
#'
#' @examples
#' merged_clean_reordered <- col_reorder(merged_clean)
col_reorder <- function(dt){
  #converts input to data table if it wasn't already.
  dt <- data.table::as.data.table(dt)

  #subsets data table to only columns that are numeric.
  mcn <- dt[,sapply(dt, is.numeric), with = FALSE]

  #Stores the correct order of samples based on their column names alphabetically.
  mco <- order(names(mcn))

  #performs the reordering of "mcn" based on order in "mco"
  col_reordered <- data.table::setcolorder(mcn, mco)

  #merges original taxa with reordered ASV table.
  all <- cbind(dt[,1:8], col_reordered)

  all
}

#' Dadaset Clean
#'
#' Clean a merged dataset output from dada2 (include species) (taxa:ASV table)
#'
#' @param df merged dataset
#' @param read_thresh Total ASV reads for samples to pass (>= read_thresh)
#'
#' @return returns dataset with only samples that passed.
#' @importFrom dplyr %>%
#' @importFrom purrr pmap
#' @export
#'
#' @examples
#' merged_clean <- dadaset_clean(merged, read_thresh = 200)
dadaset_clean <- function(df, read_thresh = 200) {
  Kingdom <- Phylum <- Family <- Order <- NULL

  #records how many taxa were in the original dataset before processing.
  nrow_taxa_ori <- NROW(df)

  #Dplyr way of subsetting dataset to taxa from the Bacterial kingdom and those which are classified at Phylum level.
  df <- df %>% dplyr::filter(Kingdom == "Bacteria" & !is.na(Phylum))

  #records how many taxa are left after subsetting.
  nrow_taxa_bac <- NROW(df)

  #Dplyr way of eliminating Mitochondria if they exist at the Family level
  if (!is.na(table(df$Family == "Mitochondria")[2])) {
    df1 <- df %>% dplyr::filter(Family != "Mitochondria" | is.na(Family))
  }
  else {df1 <- df}

  #records how many taxa are left after Mitochondria elimination
  nrow_taxa_Mit <- NROW(df1)

  #Checks for taxa classified as Chloroplast at Order level and then removes them.
  if (!is.na(table(df1$Order == "Chloroplast")[2])) {
    df2 <- df1 %>% dplyr::filter(Order != "Chloroplast" | is.na(Order))
  }
  else {df2 <- df1}

  #records how many taxa are left after Chloroplast elimination.
  nrow_taxa_Chl <- NROW(df2)

  #converts data frame to data table
  dt <- data.table::as.data.table(df2)

  #subsets data table to columns that are numeric.
  dt_numeric <- dt[,sapply(dt, is.numeric), with = FALSE]

  #dplyr way of removing samples with less than read_thresh reads.
  #first calculate sums of all columns
  dt_num_pass <- dt_numeric %>% dplyr::select_if(function(col) sum(col) >= read_thresh)
  elim <- dt_numeric %>% dplyr::select_if(function(col) sum(col) < read_thresh)
  elim_names <- names(elim)
  cleaned <- cbind(dt[,sapply(dt, is.character), with = FALSE], dt_num_pass)

  #Replacement of NA's in Genus column with "{Family} unclassified". Calls function unclassified_replace.
  cleaned$Genus <- cleaned %>% pmap(unclassified_replace)

  #Ensures Genus column is of type character.
  cleaned$Genus <- as.character(cleaned$Genus)

  #Outputs changes made to dataset.
  cat("Originally there were", nrow_taxa_ori, "taxa.", nrow_taxa_ori - nrow_taxa_bac, "were not bacterial or classified at phylum level.", nrow_taxa_bac - nrow_taxa_Mit, "were mitochondrial.", nrow_taxa_Mit - nrow_taxa_Chl, "were chloroplast.", nrow_taxa_Chl, "taxa remain. ")
  cat(length(elim_names), ifelse(length(elim_names) > 1 | length(elim_names) < 1, "samples", "sample"), "eliminated with less than", read_thresh, "reads:", elim_names)

  #at this point the dataset is cleaned and ready to move to preparation for other analyses.
  #sample type counts should be made in a list with list_types or manually for heatmap_prep

  cleaned
}

#' Preparation for Decontam
#'
#' Prepares a merged_clean dataset that has been run through dadaset_clean() and col_reorder() for running Decontam.
#'
#' @param df merged_clean dataset
#' @param type a named list where each sample is either "sample" or "control" (as in negative control).
#'
#' @return returns a phyloseq object which can then be used in decontam_histo_prev_plots()
#' @export
#'
#' @examples
#' decontam_types <- list("sample" = 28, "control" = 41, "sample" = 29)
#' pre_decontam <- decontam_prep(df = merged_clean, type = decontam_types)
decontam_prep <- function(df, type) {

  #instantiates type1 vector
  type1 <- c()
  #creates character vector where named types ("sample" or "control") are repeated as many times as denoted in list.
  for(i in 1:length(type)){
    type1 <- append(type1, values = rep(names(type)[i], type[i]))
  }

  #ensures data is of type dataframe
  df <- as.data.frame(df)
  #selects for numeric columns only
  df_num <- dplyr::select_if(df, is.numeric)

  #creates an "OTU Table" from the numeric columns of the merged_clean dataset
  OTU = phyloseq::otu_table(df_num, taxa_are_rows = TRUE)

  #taxa table created from taxa classifications Kingdom through Genus
  TAX = phyloseq::tax_table(as.matrix(df[,2:7]))

  #sample data phyloseq object SAM created with Type column populated by type1 chr vector.
  SAM = phyloseq::sample_data(data.frame(
    Type = type1,
    row.names = names(df_num),
    stringsAsFactors = FALSE
  ))
  #combine all phyloseq tables into one phyloseq object
  physeq <- phyloseq::phyloseq(OTU, TAX, SAM)
  physeq
}


#' Generate Decontam Histogram and Prevalence Plots
#'
#' To be used on a Phyloseq object to generate a histogram of Decontam scores and a prevalence plot.
#'
#' @param physeq phyloseq object with OTU, TAX, and SAM components
#' @param thresh Decontam score threshold. Default is 0.5
#' @param study_name name of the study where the dataset is from.
#'
#' @return returns histogram and prevalence plots.
#' @export
#'
#' @examples
#' decontam_types <- list("sample" = 28, "control" = 41, "sample" = 29)
#' pre_decontam <- decontam_prep(df = merged_clean, type = decontam_types)
#' decontam_histo_prev_plots(pre_decontam, thresh = 0.5, study_name = "Study Name")
decontam_histo_prev_plots <- function(physeq, thresh = 0.5, study_name){
  p <- prev <- physeq.neg <- physeq.pos <- contaminant <- NULL

  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$Type == "control"

  #runs isNotContaminant() which splits taxa up into contaminants and true taxa.
  #isNotContaminant() is chosen here over isContaminant() since samples are presumed low biomass and the majority of taxa are assumed to be contaminants.
  #Keep in mind that taxa with true are not contaminants.
  contamdf.prev <- decontam::isNotContaminant(physeq, neg = "is.neg", detailed = TRUE, threshold = thresh)
  #counts how many true taxa and contaminants there are and prints it to the console.
  print("True represents true taxa:")
  print(table(contamdf.prev$not.contaminant))

  #code for histogram plotting to evaluate for an appropriate threshold
  #Green tones used for the prevalence palette
  prevalencePalette <- c("2" = "#edf8e9", "3-5" = "#bae4b3", "6-10" = "#74c476", "11+" = "#238b45")
  #getPalette local function calls colorRamp Palette on the prevalence palette
  getPalette = grDevices::colorRampPalette(prevalencePalette)
  #determines the number of steps that the histogram should have (how many unique counts were there at each decontam score level)
  steps <- length(table(contamdf.prev$prev))
  #splits the palette to number of steps determined.
  colr <- getPalette(steps)
  #creates the title for histogram
  histo_title <- paste("Decontam Histgram for", study_name)
  #fill has to be in factor form and then use scale fill manual() to select palatte to use.
  decontam_histogram <- ggplot2::ggplot(contamdf.prev, ggplot2::aes(x = p, fill = factor(prev))) + ggplot2::geom_histogram(bins = 100) + ggplot2::labs(x = 'Decontam score', y = 'Number of ASVs') + ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05)) + ggplot2::scale_fill_manual(values = colr) + ggplot2::theme_dark() + ggplot2::ggtitle(histo_title)

  #code for prevalence plotting
  #converts sample counts into prevalence counts by turning values greater than 0 into 1's.
  physeq.pa <- phyloseq::transform_sample_counts(physeq, function(abund) 1*(abund>0))

  #subsets samples which are controls into phyloseq object physeq.pa.neg
  physeq.pa.neg <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$Type == "control", physeq.pa)

  #subsets samples which are actual samples into phyloseq object physeq.pa.pos
  physeq.pa.pos <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$Type == "sample", physeq.pa)

  #creates a data frame where taxa prevalence is summed up across control and true samples. Contaminants are indicated in !contamdf.prev$not.contaminant
  df.pa <- data.frame(physeq.pos= phyloseq::taxa_sums(physeq.pa.pos), physeq.neg = phyloseq::taxa_sums(physeq.pa.neg),
                      contaminant=!contamdf.prev$not.contaminant)

  #creates title for prevalence plot
  prev_title <- paste("Contaminant Prevalence at Threshold:", thresh, "for", study_name, "Study")

  #plots prevalence with taxa negative control prevalence on the x-axis and taxa sample prevalence on the y-axis.
  prev_plot <- ggplot2::ggplot(data=df.pa, ggplot2::aes(x=physeq.neg, y=physeq.pos, color=contaminant)) + ggplot2::geom_point() +
    ggplot2::xlab("Prevalence in Negative Controls") + ggplot2::ylab("Prevalence in Samples") + ggplot2::ggtitle(prev_title)

  #returns plot list
  list(Histogram = decontam_histogram, Prevalence_Plot = prev_plot)
}

#' Decontaminate Dataset
#'
#' Final function for using Decontam.
#'
#' @param df dataframe from which the phyloseq object was generated with decontam_prep()
#' @param physeq phyloseq object created from running decontam_prep()
#' @param thresh threshold chosen from looking at histogram and prevalence plots
#'
#' @return gives a list of dataframes split into true taxa, contaminants, and the chosen threshold.
#' @export
#'
#' @examples
#' decontam_types <- list("sample" = 28, "control" = 41, "sample" = 29)
#' pre_decontam <- decontam_prep(df = merged_clean, type = decontam_types)
#' decontaminated <- decontaminate(merged_clean, pre_decontam, thresh = 0.1)
#' true_taxa_0.1 <- decontaminated$TrueTaxa; contaminants_0.1 <- decontaminated$Contaminants
decontaminate <- function(df, physeq, thresh) {

  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$Type == "control"

  #isNotContaminant() is chosen here over isContaminant() since samples are presumed low biomass and the majority of taxa are assumed to be contaminants.
  #Keep in mind that Trues are not contaminants
  contamdf.prev <- decontam::isNotContaminant(physeq, neg = "is.neg", detailed = TRUE, threshold = thresh)

  #subsetting contaminants and 'true taxa'
  true_taxa <- subset(df, contamdf.prev$not.contaminant)
  contaminants <- subset(df, !contamdf.prev$not.contaminant)

  #returns list of dataframes true taxa and contaminants along with the threshold used.
  list(TrueTaxa = true_taxa, Contaminants = contaminants, Threshold = thresh)
}

#' Create Sample Type List
#'
#' Crude way to create Sample Type Lists for later PCoA plotting, heatmapping, or barplotting analyses.
#'
#' @param df dataset with sample names on columns or rows.
#' @param types_on_col if sample names are on columns this should be true, if on rows this should be false
#' @param split_by character to separate sample type name from unique identifier. Default is underscore.
#'
#' @return Returns list of sample types and their counts based on the first letter of every sample.
#' @export
#'
#' @examples
#' Sample_Types <- list_types(merged_clean, types_on_col = TRUE, split_by = "_")
list_types <- function(df, types_on_col = FALSE, split_by = "_"){
  X <- NULL

  #Determines whether column names or row names should be split into sample types.
  if (types_on_col == TRUE) {
    #Takes first letter of each column name and puts it into a dataframe.
    df2 <- as.data.frame((substr(colnames(df)[9:NCOL(df)], 0, 1)))

    #Names column "X"
    colnames(df2) <- "X"

    #makes tally of each unique first letter.
    type_counts <- df2 %>% dplyr::count(X)
    #makes note of unique column names using the first part of the sample name and splits by user input character.
    type_names <- unique(sapply(strsplit(colnames(df)[9:NCOL(df)], split = split_by), `[`, 1))
  }
  else {
    #counts of rownames.
    type_counts <- df %>% dplyr::count(substr(rownames(df), 0, 1))

    #names of unique row names using the first part of the sample name and splits by user input character.
    type_names <- unique(sapply(strsplit(rownames(df)[1:NROW(df)], split = split_by), `[`, 1))
  }

  #Replaces underscores in sample type names to spaces and capitalizes first letter.
  type_names <- stringr::str_to_title(stringr::str_replace(type_names, pattern = "_", replacement = " "))

  #initiates sample type list with sample type counts
  type_list <- as.list(as.vector(type_counts$n))

  #gives names to each sample type
  names(type_list) <- type_names

  type_list
}

#' Preparation for Heatmapping
#'
#' Takes a cleaned merged dataset from dadaset_clean() and prepares the dataset for heatmapping. Only one cutoff metric should be used at a time.
#'
#' @param df dataset in merged format
#' @param taxa_cutoff how many taxa to keep. Default is zero
#' @param mean_ab_cutoff whether or not to apply a 1% mean abundance cutoff.
#'
#' @return returns a dataset ready to input into a heatmapping function.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' preheat_top_10 <- heatmap_prep(merged_clean, taxa_cutoff = 10, mean_ab_cutoff = FALSE)
#' preheat_mean_ab <- heatmap_prep(merged_clean, taxa_cutoff = 0, mean_ab_cutoff = TRUE)
heatmap_prep <- function(df, taxa_cutoff = 0, mean_ab_cutoff = TRUE){
  #does rank in order of abundance.

  #Removes all taxa columns except for "Genus". Checks if "Species" column is included in taxa classification.
  if("Species" %in% colnames(df)){
    df_preheat <- df %>% dplyr::select(!c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Species"))
  }
  else {
    df_preheat <- df %>% dplyr::select(!c("X", "Kingdom", "Phylum", "Class", "Order", "Family"))
  }

  #transposes dataset
  dt_preheat <- t(df_preheat)

  #Converts values in Genus Row to column names
  colnames(dt_preheat) <- unlist(dt_preheat[row.names(dt_preheat)=='Genus',])

  #Removes Genus row and converts dataset to dataframe
  dt_preheat <- as.data.frame(dt_preheat[!row.names(dt_preheat)=='Genus',])

  #Makes sure that values are of type numeric and not factors. Then converts dataset to dataframe
  dt_preheat1 <- as.data.frame(sapply(dt_preheat, function(x) as.numeric(as.character(x))))

  #row names lost in the process are added back from the original dataset.
  rownames(dt_preheat1) <- rownames(dt_preheat)

  #column names lost in the process are added back from the original dataset but also made unique.
  colnames(dt_preheat1) <- make.unique(colnames(dt_preheat1))

  #Now, with a dataset with unique genera as columns and samples on rows abundance cutoffs or taxa cutoffs are applied.
  if (mean_ab_cutoff == TRUE) {
    #read counts are divided by total read counts per sample (per row sums).
    data.prop <- dt_preheat1/rowSums(dt_preheat1)

    #calculates the average fraction of sample that each taxa, which are on columns, make up.
    mean_abundance <- apply(data.prop, 2, mean)

    #stores the column names(taxa) for which mean abundance is less than 1%.
    n2 <- names(which(mean_abundance < 0.01))

    #removes taxa less than 1% average relative abundance
    dt_preheat1.ab <- dt_preheat1[, -which(names(data.prop) %in% n2)]

    #Samples are eliminated which no longer have taxa abundance.
    dt_preheat1.ab <- row_sample_elim(dt_preheat1.ab, thresh = 1)

    return(dt_preheat1.ab)
  }

  #if a non-zero taxa cutoff is given, proceeds to rank taxa based on total read count across all samples.
  if (taxa_cutoff != 0) {

    #stores ranks of column sums and splits ties based on which count came up first.
    dt_rank <- rank(-colSums(dt_preheat1), ties.method = "first")

    #subsets taxa based on their rank.
    dt_preheat1 <- as.data.frame(t(subset(t(dt_preheat1), dt_rank <= taxa_cutoff)))

    #checks to make sure all samples still have a read count of at least one. Employs function "row_sample_elim()"
    dt_preheat1 <- row_sample_elim(dt_preheat1, thresh = 1)

    return(dt_preheat1)
  }
}


#' Heatmapping
#'
#' Generates a heatmap of a merged dataset after processing with heatmap_prep()
#'
#' @param df dataset after running heatmap_prep()
#' @param scalecolor color theme to be used for the heatmap. Options are "red" or "blue"
#' @param types sample types as a named list.
#' @param title title for heatmap
#' @param by_percent logical. If true, taxa in heatmap are shown as percentage of samples. If false, taxa are shown as absolute abundances.
#'
#' @return returns a heatmap generated from library(ComplexHeatmap)
#' @export
#'
#' @examples
#' types_list <- list("Amnion-Chorion" = 28, "Control" = 41, "Villous Tree" = 29)
#' preheat_top_10 <- heatmap_prep(merged_clean, taxa_cutoff = 10, mean_ab_cutoff = FALSE)
#' heatmapping(df = preheat_top_10, scalecolor = "blue", types = types_list, title = "Example Heatmap")
heatmapping <- function(df, scalecolor, types, title, by_percent = TRUE) {

  #converts dataset to proportionalized dataset (taxa percentage of sample).
  data.prop <- (df/rowSums(df))*100
  #create color palette with options for red or blue.
  if(scalecolor == "red"){
    finalscale <- c("lightyellow", "red")
  }
  else if(scalecolor == "blue"){
    #alternate black and blue color palette
    finalscale <- c("#000033", "#66CCFF")
  }
  #instantiate factortypes and repeat type names based on values in named list types.
  factortypes <- c()
  for(i in 1:length(types)){
    factortypes <- as.factor(append(factortypes, values = rep(names(types)[i], types[[i]])))
  }

  if(by_percent == TRUE) {
    #for percentage heatmap
    ht = ComplexHeatmap::Heatmap(matrix = as.matrix(t(data.prop)), column_title = names(types), show_column_dend = FALSE, col = finalscale, name = "Percentage of Sample", cluster_columns = FALSE, row_names_side = "left", row_dend_side = "right", row_dend_width = grid::unit(1.5, "cm"), column_split = factortypes, border = TRUE)
    ComplexHeatmap::draw(ht, column_title = title)
  }
  else{
    #for absolute abundance
    ht = ComplexHeatmap::Heatmap(matrix = as.matrix(t(df)), column_title = names(types), show_column_dend = FALSE, col = finalscale, name = "Absolute Reads", cluster_columns = FALSE, row_names_side = "left", row_dend_side = "right", row_dend_width = grid::unit(1.5, "cm"), column_split = factortypes, border = TRUE)
    ComplexHeatmap::draw(ht, column_title = title)
  }
}

#' Remove Control Taxa
#'
#' Removes taxa found in negative control samples.
#'
#' @param df merged and cleaned dataset to be worked on
#' @param control_sample_columns numeric indices of negative control columns.
#'
#' @return returns a dataset without taxa found in negative controls.
#' @export
#'
#' @examples
#' merged_no_control_taxa <- neg_taxa_remove(merged_clean, control_sample_columns = 37:77)
neg_taxa_remove <- function(df, control_sample_columns) {

  #selects for numeric columns only
  df_num <- dplyr::select_if(df, is.numeric)
  #selects control sample columns
  dt <- data.table::as.data.table(df)
  dt_controls <- dt[,control_sample_columns, with = FALSE]
  df_no_control_taxa <- subset(df, !(rowSums(dt_controls) >= 1))
  df_control_taxa <- subset(df, (rowSums(dt_controls) >= 1))
  list("No Control Taxa" = df_no_control_taxa, "Control Taxa" = df_control_taxa)
}

#' Taxa Prevalence Determination
#'
#' Determines taxa prevalence for a subset of a merged dataset.
#'
#' @param df merged_clean dataset or other dataset where taxa are on rows and samples are on columns.
#' @param subset Logical. If true, must include lt_col and rt_col parameters.
#' @param lt_col leftmost column of subset. Can be numeric or name
#' @param rt_col rightmost column of subset. Can be numeric or name
#'
#' @return returns prevalence of taxa in descending order.
#' @export
#'
#' @examples
#' Taxa_Prev <- prevalence(merged_clean, subset = TRUE, lt_col = 9, rt_col = NCOL(merged_clean))
prevalence <- function(df, subset = FALSE, lt_col, rt_col) {
  prev <-NULL
  #must have Genus column

  #
  ifelse(subset == TRUE, total <- (rt_col - lt_col), total <- sum(sapply(df, is.numeric)))

  #converts data frame to data table.
  dt <- data.table::as.data.table(df)
  #subsets datatable to only columns of type numeric. Or numeric columns within specified bounds.
  ifelse (subset == FALSE, dt_num <- dt[, sapply(dt, is.numeric), with = FALSE], dt_num <- dt %>% dplyr::select(lt_col:rt_col))

  #rownames for subsetted dataset are taken from Genus column of original dataset.
  rownames(dt_num) <- make.unique(dt$Genus)

  #converting dataset to prevalence.
  dt_prev <- sapply(dt_num, function(x) x/x)
  dt_prev[is.nan(dt_prev)] <- 0
  dt_prev <- as.data.frame(dt_prev)

  #calculating prevalence of taxa in new prev column.
  dt_prev$prev <- rowSums(dt_prev)
  rownames(dt_prev) <- make.unique(dt$Genus)

  #prevalence listed in descending order.
  dt_ordered <- dt_prev %>% dplyr::arrange(dplyr::desc(prev))

  #Stores reordered rownames in chr vector Prev_names
  Prev_names <- rownames(dt_ordered)

  #Takes prevalence calculated from row sums and stores it in data frame.
  Prev_nums <- as.data.frame(dt_ordered$prev)

  #Combines names with prevalence.
  Prev <- cbind(Prev_names, Prev_nums, stringsAsFactors = FALSE)

  #renames columns
  colnames(Prev)[1:2] <- c("Taxa", paste("Prevalence out of", total))
  Prev
}

#' Create Excel Prevalence Workbook
#'
#' Combines prevalence objects into different sheets of an excel workbook.
#'
#' @param objects Objects as a list i.e. list(prev_x, prev_y, prev_z)
#' @param sheetnames Names for sheets of workbook i.e. c("Prev X", "Prev Y", "Prev Z")
#' @param wb_name Workbook Name
#'
#' @return creates an excel workbook with all prevalence objects as their own sheet.
#' @export
#'
#' @examples
#' Prev_AC <- prevalence(merged_clean, subset = TRUE, lt_col = 9, rt_col = 36)
#' Prev_C <- prevalence(merged_clean, subset = TRUE, lt_col = 37, rt_col = 77)
#' Prev_V <- prevalence(merged_clean, subset = TRUE, lt_col = 78, rt_col = 106)
#' Prevs <- list(Prev_AC, Prev_C, Prev_V)
#' prev_write_xlsx(objects = Prevs, sheetnames = c("Prev AC", "Prev Ctrl", "Prev V"), wb_name = "Prev")
prev_write_xlsx <- function(objects, sheetnames, wb_name = "Prevalence") {

  #instantiates workbook
  wb <- openxlsx::createWorkbook()

  #for every Prevalence object adds a sheet and names it. Then the data is written to the newly created sheet with column names included.
  for (i in 1:length(objects)) {
    openxlsx::addWorksheet(wb = wb, sheetName = sheetnames[i])
    openxlsx::writeDataTable(wb, sheet = sheetnames[i], as.data.frame(objects[i]), colNames = TRUE)
  }
  #Checks to see if file already exists with same name. If so, workbook will not save.
  if(file.exists(paste0(wb_name, ".xlsx"))) {
    #Will fail to write if xlsx with indicated wb_name already exists.
    print(paste0("File: ", wb_name, ".xlsx ", "already exists. Please retry with a different wb_name."))
    Failed <- TRUE
  }
  else {
    Failed <- FALSE
    #saves workbook
    openxlsx::saveWorkbook(wb, paste0(wb_name, ".xlsx"), overwrite = FALSE)
  }

  if(file.exists(paste0(wb_name, ".xlsx")) && Failed == FALSE) {
    print(paste0("File: ", wb_name, ".xlsx ", "created!"))
  }

}

#' Eliminate Samples in Rows
#'
#' Checks samples to verify that they have at least a certain number of reads left. Those that don't pass are eliminated.
#'
#' @param df dataset of type data frame with samples on rows and taxa in columns.
#' @param thresh Read threshold. Default is 1 i.e. samples with no reads are eliminated.
#'
#' @return returns a dataframe in which all samples remaining have reads.
#' @export
#'
#' @examples
#' merged_preheat_pass <- row_sample_elim(preheat_ab, thresh = 1)
row_sample_elim <- function(df,thresh = 1) {
  #Calculates row sums for data frame assuming samples are on rows.
  df$rowsums <- rowSums(df)
  #subsets data frame to samples below threshold
  elim <- df[df$rowsums<thresh,]
  #subsets data frame to samples above or equal to threshold
  df <- df[df$rowsums>=thresh,]
  #removes the rowsums column added in the beginning.
  df <- df[,-ncol(df)]
  cat("Samples that had less than ", thresh, "reads left:", row.names(elim))
  df
}

#' Replace Unclassified Genus
#'
#' This function is used in dadaset_clean to replace NA's with "{Family} Unclassified"
#'
#' @param Family Family column of merged dataset
#' @param Genus Genus column of merged dataset
#' @param ... For additional arugments
#'
#' @return Returns Genus column with replaced NA's
#' @export
unclassified_replace <- function(Family, Genus, ...) {
  #This function is used in dadaset_clean to replace NA's with "{Family} Unclassified"
  #If cell in Genus is NA, makes replacement, otherwise it leaves the cell alone.
  ifelse(is.na(Genus) == TRUE,
         Genus <- glue::glue("{Family} Unclassified"),
         Genus <- Genus)
}
