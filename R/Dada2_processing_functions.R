if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Agglomeration
#'
#' Takes taxonomy and ASV/OTU table, and merges ASVs which share an indicated taxonomic classification.
#'
#' @param df numeric ASV/OTU table dataframe with ASVs on rows as row names and samples on columns as column names
#' @param df_tax taxonomy table with ASVs on rows as rownames (matching `df`).
#' @param agg_class character vector of column name in `df_tax` by which to agglomerate.
#'
#' @return a dataframe with merged ASV counts and agglomerated taxonomy as rownames.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' Theis_ASV <- Theis_merged_clean[,9:ncol(Theis_merged_clean)]
#' Theis_TAX <- Theis_merged_clean[,1:8]
#' Theis_agg <- agglomerate(df = Theis_ASV, df_tax = Theis_TAX, agg_class = "Genus")
agglomerate <- function(df, df_tax, agg_class){
  df_final <- c()
  class_list <- unique(df_tax[[agg_class]])
  for(i in 1:length(class_list)){
    df_sub <- df %>% dplyr::filter(df_tax[[agg_class]] == class_list[i])
    df_sums <- as.data.frame(t(data.frame(CS = colSums(df_sub))))
    rownames(df_sums)[1] <- unique(df_tax[[agg_class]])[i]
    df_final <- rbind(df_final, df_sums)
  }
  df_final
}

#' Barplot Preparation
#'
#' Prepares data for barplotting. Recommended to use result from `heatmap_prep()` for `df`.
#'
#' @param df numeric dataframe with ASVs on columns and samples on rows.
#' @param df_meta sample metadata with samples on rows.
#' @param type_col character vector of column name in `df_meta` to be used for sample grouping
#' @param study character vector of study name/title.
#' @param DECONTAM_dataset boolean vector indicating whether dataset has been processed through DECONTAM.
#'
#' @return dataframe ready for `barplotting()`
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' library(dplyr)
#' Theis_META_P <- Theis_rare_META %>% filter(Type == "Placenta")
#' Theis_agg <- agglomerate(df = Theis_rare_ASV, df_tax = Theis_rare_TAX, agg_class = "Genus")
#' Theis_agg_tax <- data.frame(Genus = rownames(Theis_agg), row.names = rownames(Theis_agg))
#' Theis_agg_P <- Theis_agg %>% select_if(colnames(.) %in% Theis_META_P$Run)
#' Theis_agg_tax_P <- Theis_agg_tax %>% filter(rownames(.) %in% rownames(Theis_agg_P))
#' Theis5 <- get_top_ASVs(Theis_agg_P, 5)
#' Theis_preheatbar <- heatmap_prep(df = Theis_agg_P,
#' df_tax = Theis_agg_tax_P,
#' class_col = "Genus",
#' df_meta = Theis_META_P,
#' select_ASVs = Theis5,
#' mean_ab_cutoff = FALSE)
#' Theis_prebar <- barplot_prep(df = Theis_preheatbar,
#' df_meta = Theis_META_P,
#' type_col = "Del_GA",
#' study = "Theis",
#' DECONTAM_dataset = FALSE)
barplot_prep <- function (df, df_meta, type_col, study, DECONTAM_dataset = FALSE) {
  prebar <- c()
  df1 <- as.data.frame(t(df))
  df_ordered <- df1[df1 %>% rowSums(.) %>% sort(., decreasing = TRUE) %>% names(.),]
  uniques <- sort(unique(df_meta[[type_col]]))
  if("Technical Control" %in% uniques){
    uniques <- uniques[c(match("Technical Control", uniques), which(uniques %ni% "Technical Control"))]
  }
  for(i in 1:length(uniques)){
    df_subset <- df_ordered %>% dplyr::select_if(df_meta[[type_col]] == uniques[i])
    prebar <- as.data.frame(rbind(prebar, data.frame(
      ASVs = factor(levels = rev(rownames(df_subset)), x = rownames(df_subset)),
      Type = factor(levels = uniques, x = rep(uniques[i], dim(df_subset)[1])),
      Genus = factor(levels = make.unique(rev(sub(pattern = "(ASV\\d+)(-)(.*)", x = rownames(df_subset), replacement = "\\3"))), x = make.unique(sub(pattern = "(ASV\\d+)(-)(.*)", x = rownames(df_subset), replacement = "\\3"))),
      Total_Reads = as.numeric(rowSums(df_subset)),
      Relative_Reads = as.numeric(100*prop.table(as.matrix(rowSums(df_subset)),2)),
      Study = rep(study, dim(df_subset)[1]),
      DECONTAM = rep(DECONTAM_dataset, dim(df_subset)[1])
    )))
  }
  return(prebar)
}

#' Barplotting
#'
#' Function used to create stacked barplots of absolute and relative reads of sample types. Plots will be arranged so that after DECONTAM dataset is below regular version.
#'
#' @param df dataframe prepared with barplot_prep()
#' @param study name of the study or origin of the data
#' @param deco_df optional dataframe prepared with barplot_prep() but processed with DECONTAM
#' @param legend_cols numeric vector of number of columns to split legend into.
#' @param leg_title character vector to name legend.
#'
#' @return returns a plot object.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' library(dplyr)
#' Theis_META_P <- Theis_rare_META %>% filter(Type == "Placenta")
#' Theis_agg <- agglomerate(df = Theis_rare_ASV, df_tax = Theis_rare_TAX, agg_class = "Genus")
#' Theis_agg_tax <- data.frame(Genus = rownames(Theis_agg), row.names = rownames(Theis_agg))
#' Theis_agg_P <- Theis_agg %>% select_if(colnames(.) %in% Theis_META_P$Run)
#' Theis_agg_tax_P <- Theis_agg_tax %>% filter(rownames(.) %in% rownames(Theis_agg_P))
#' Theis5 <- get_top_ASVs(Theis_agg_P, 5)
#' Theis_preheatbar <- heatmap_prep(df = Theis_agg_P,
#' df_tax = Theis_agg_tax_P,
#' class_col = "Genus",
#' df_meta = Theis_META_P,
#' select_ASVs = Theis5,
#' mean_ab_cutoff = FALSE)
#' Theis_prebar <- barplot_prep(df = Theis_preheatbar,
#' df_meta = Theis_META_P,
#' type_col = "Del_GA",
#' study = "Theis",
#' DECONTAM_dataset = FALSE)
#' Theis_barplot <- barplotting(df = Theis_prebar,
#' study = "Theis",
#' legend_cols = 2,
#' leg_title = "Genera")
barplotting <- function(df, deco_df, study, legend_cols = 1, leg_title) {
  #this will plot the barplot in one section of the screen and the legend in another.

  if(!missing(deco_df)){
    #ensures dataset is a dataframe if not already.
    df <- as.data.frame(df)
    deco_df <- as.data.frame(deco_df)
    totalsD <- deco_df %>% dplyr::group_by(.data$Type) %>% dplyr::summarise(sum = sum(.data$Total_Reads))
    totals <- df %>% dplyr::group_by(.data$Type) %>% dplyr::summarise(sum = sum(.data$Total_Reads))

    All_combined <- rbind(df, deco_df)
    abundance_totals <- All_combined %>% dplyr::group_by(.data$ASVs) %>% dplyr::summarise(sum = sum(.data$Total_Reads))
    abundance_totals <- abundance_totals[order(abundance_totals$sum, decreasing = TRUE),]
    abundance_totals$ASVs <- factor(x = abundance_totals$ASVs, levels = abundance_totals$ASVs)
    #creates the local function getPalette which uses the colorRampPalette function of the RColorBrewer package and the Spectral palette which is 11 colors.
    getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    colr <- getPalette(length(unique(All_combined$ASVs)))
    abundance_totals$Color <- colr
    abundance_totals$Color[match("Other", abundance_totals$ASVs)] <- "#808080"


    #getPalette() is used to extrapolate the Spectral palette to any number of steps beyond the original 11.
    Leg_vals <- abundance_totals %>% dplyr::filter(.data$ASVs %in% All_combined$ASVs) %>% dplyr::pull(.data$Color); names(Leg_vals) <- abundance_totals %>% dplyr::filter(.data$ASVs %in% All_combined$ASVs) %>% dplyr::pull(.data$ASVs)
    for_legend <- ggplot2::ggplot(data = All_combined, aes(x=.data$Type, y = .data$Relative_Reads, fill = .data$ASVs)) + ggplot2::geom_bar(stat="identity") + ggplot2::theme(legend.text = element_text(face = "italic")) + ggplot2::ylab("Percent Relative Abundance") + ggplot2::scale_fill_manual(name = leg_title, values = Leg_vals, guide = guide_legend(ncol = legend_cols)) + ggplot2::ggtitle(paste("For Legend", "by type"))
    Legend_alone <- cowplot::get_legend(for_legend)

    Rel_vals <- abundance_totals %>% dplyr::filter(.data$ASVs %in% df$ASVs) %>% dplyr::pull(.data$Color); names(Rel_vals) <- abundance_totals %>% dplyr::filter(.data$ASVs %in% df$ASVs) %>% dplyr::pull(.data$ASVs)
    Relative <- ggplot2::ggplot(data = df, aes(x= .data$Type, y = .data$Relative_Reads, fill = .data$ASVs)) + ggplot2::geom_bar(stat="identity") + ggplot2::ylab("Percent Relative Abundance") + ggplot2::labs(x= "") + ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_manual(values = Rel_vals)

    RelD_vals <- abundance_totals %>% dplyr::filter(.data$ASVs %in% deco_df$ASVs) %>% dplyr::pull(.data$Color); names(RelD_vals) <- abundance_totals %>% dplyr::filter(.data$ASVs %in% deco_df$ASVs) %>% dplyr::pull(.data$ASVs)
    RelativeD <- ggplot2::ggplot(data = deco_df, aes(x= .data$Type, y = .data$Relative_Reads, fill = .data$ASVs)) + ggplot2::geom_bar(stat="identity") + ggplot2::ylab("Percent Relative Abundance") + ggplot2::labs(x= "") + ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_manual(values = RelD_vals)
    plot <- cowplot::plot_grid(Relative, RelativeD, Legend_alone, labels = c("A", "B", ""), ncol = 1, nrow = 3, rel_heights = c(1, 1,.75), vjust = 1)
  } else {

    df <- as.data.frame(df)
    All_combined <- df
    totals <- df %>% dplyr::group_by(.data$Type) %>% dplyr::summarise(sum = sum(.data$Total_Reads))
    abundance_totals <- All_combined %>% dplyr::group_by(.data$ASVs) %>% dplyr::summarise(sum = sum(.data$Total_Reads))
    abundance_totals <- abundance_totals[order(abundance_totals$sum, decreasing = TRUE),]
    abundance_totals$ASVs <- factor(x = abundance_totals$ASVs, levels = abundance_totals$ASVs)
    #creates the local function getPalette which uses the colorRampPalette function of the RColorBrewer package and the Spectral palette which is 11 colors.
    getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    colr <- getPalette(length(unique(All_combined$ASVs)))
    abundance_totals$Color <- colr
    abundance_totals$Color[match("Other", abundance_totals$ASVs)] <- "#808080"


    #getPalette() is used to extrapolate the Spectral palette to any number of steps beyond the original 11.
    Leg_vals <- abundance_totals %>% dplyr::filter(.data$ASVs %in% All_combined$ASVs) %>% dplyr::pull(.data$Color); names(Leg_vals) <- abundance_totals %>% dplyr::filter(.data$ASVs %in% All_combined$ASVs) %>% dplyr::pull(.data$ASVs)
    for_legend <- ggplot2::ggplot(data = All_combined, aes(x=.data$Type, y = .data$Relative_Reads, fill = .data$ASVs)) + ggplot2::geom_bar(stat="identity") + ggplot2::theme(legend.text = element_text(face = "italic")) + ggplot2::ylab("Percent Relative Abundance") +  ggplot2::scale_fill_manual(name = leg_title, values = Leg_vals, guide = guide_legend(ncol = legend_cols)) + ggplot2::ggtitle(paste("For Legend", "by type"))
    Legend_alone <- cowplot::get_legend(for_legend)

    Rel_vals <- abundance_totals %>% dplyr::filter(.data$ASVs %in% df$ASVs) %>% dplyr::pull(.data$Color); names(Rel_vals) <- abundance_totals %>% dplyr::filter(.data$ASVs %in% df$ASVs) %>% dplyr::pull(.data$ASVs)
    Relative <- ggplot2::ggplot(data = df, aes(x= .data$Type, y = .data$Relative_Reads, fill = .data$ASVs)) + ggplot2::geom_bar(stat="identity") + ggplot2::ylab("Percent Relative Abundance") + ggplot2::labs(x= "") + ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_manual(values = Rel_vals)

    plot <- cowplot::plot_grid(Relative, Legend_alone, labels = c("A", ""), ncol = 1, nrow = 2, rel_heights = c(1, .75), vjust = 1)
  }

  return(plot)
}

#' Beta Diversity Plotting
#'
#' A beta diversity plotting function built in base R to make a PCoA plot of Bray-Curtis distances between samples.
#'
#' @param df merged data frame with taxonomy on rows.
#' @param df_meta metadata dataframe with samples as rows
#' @param type character vector indicating `df_meta` column name of column for primary sample grouping (by color).
#' @param study character vector indicating dataset study name or title for plot.
#' @param inset numeric vector designating horizontal position of legend on plot. Positive values move the legend left, negative values move it to the right. Default is 0
#' @param invert_x boolean vector indicating whether plot should be inverted over the y-axis.
#' @param invert_y boolean vector indicating whether plot should be inverted over the x-axis.
#' @param ASP numeric vector. Aspect ratio for the plot. Default is 1
#' @param Legend_Yspace numeric vector which modifies the spacing between legend entries. Default is 1, Range between 0 and 1.
#' @param rel boolean vector indicating whether data in `df` is already converted to percent relative abundances or not. Default is `FALSE`
#' @param red_blue boolean vector indicating whether red/blue color palette is desired. Default is `FALSE`
#' @param all_O boolean vector indicating whether all points should be plotted as open circles. Defualt is `FALSE`
#' @param taxa_on_rows boolean vector indicating whether ASVs/OTUs in `df` taxonomy are on rows or on columns. Default is `TRUE`
#' @param highlight optional character vector of all samples to be colored normally
#' @param nonhgrey optional boolean vector indicating whether non-highlighted samples should be grey (`TRUE`) or left blank (`FALSE`)
#' @param type2 optional character vector of column name from `df_meta` to perform secondary grouping on dataset (by shape)
#' @param legendyn boolean vector indicating whether to include legend in plot. Default is `TRUE`.
#' @param size optional numeric vector indicating size of all points plotted.
#' @param transp optional numeric vector indicating transparency of all points to be plotted. Range from 0 to 1.
#'
#' @return plots beta diversity PCoA plot and returns a dataframe containing eigenvalues and vectors for all possible plotting axes.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' Beta_plot <- beta_div(df = Theis_rare_ASV,
#' df_meta = Theis_rare_META,
#' type = "Type",
#' study = "Example Plot",
#' all_O = TRUE,
#' taxa_on_rows = TRUE,
#' size = 2,
#' legendyn = FALSE)
beta_div <- function(df, df_meta, type, study, inset = 0, invert_x = FALSE, invert_y = FALSE, ASP = 1, Legend_Yspace = 1, rel = FALSE, red_blue = FALSE, all_O = FALSE, taxa_on_rows = TRUE, highlight, nonhgrey = FALSE, type2, legendyn = TRUE, size, transp) {
  #records current palette
  pal <- grDevices::palette()

  #subsets data frame to only columns that are numeric.
  dt_numeric <- df %>% dplyr::select_if(., is.numeric)

  #converts total counts to relative abundances.
  if(rel == TRUE) {
    dt_numeric <- dt_numeric/colSums(dt_numeric)
  }

  #transposes the dataset so that species are on columns and samples are on rows.
  if(taxa_on_rows == TRUE){
    df_numeric_t <- as.data.frame(t(dt_numeric))
  }
  else {
    df_numeric_t <- as.data.frame(dt_numeric)
  }

  #gets factor data for sample type from metadata
  X <- df_meta[[type]]

  #uses the vegdist() function to calculate bray-curtis dissimilarity indices
  df.bray <- vegan::vegdist(df_numeric_t[1:NCOL(df_numeric_t)])

  #transforms Bray-Curtis indices to PCoA coordinates
  df.pcoa <- ecodist::pco(df.bray, negvals = "rm")

  #calculates Relative Eigenvalues for axes 1 and 2
  axis1_rel_eigen <- round(100*df.pcoa$values[1]/sum(df.pcoa$values), digits = 2)
  axis2_rel_eigen <- round(100*df.pcoa$values[2]/sum(df.pcoa$values), digits = 2)

  #options for inverting over the X or Y axis.
  ifelse(invert_x == TRUE, invert_x <- -1, invert_x <- 1)
  ifelse(invert_y == TRUE, invert_y <- -1, invert_y <- 1)

  #blue and red color scheme
  if (red_blue == TRUE){
    cols <- c("#353795", "#a41d24")
  } else {
    cols <- 0
  }

  if(!missing(highlight)){
    getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    colr <- getPalette(length(unique(highlight)))
    highlight_bool <- X %in% highlight
    col_vec <- color_rep(colr, X, highlight, highlight_bool, nonhgrey)
  }
  if(!missing(transp)){
    transparent_colr <- grDevices::palette.colors(n = length(unique(df_meta[[type]])), palette = "R4", alpha = transp)
    transp_vec <- color_rep(transparent_colr, X, nonhgrey = FALSE)
  }

  graphics::par(xpd = FALSE)
  #plots the PCoA. pch parameter is used to change the shape of point used(i.e. star instead of dot)
  if(!missing(type2)){
    graphics::plot(x = invert_x*(df.pcoa$vectors[,1]), y = invert_y*(df.pcoa$vectors[,2]), col= if(!missing(highlight)){col_vec}else{if(!missing(transp)){transp_vec}else{as.numeric(as.factor(df_meta[[type]]))}}, lwd = 3, cex =if(!missing(size)){size}else{1}, pch= as.numeric(as.factor(df_meta[[type2]])), main = study, xlab = paste0("PC1 (", axis1_rel_eigen, "%)"), ylab = paste0("PC2 (", axis2_rel_eigen, "%)"), asp = ASP)
  } else{
    graphics::plot(x = invert_x*(df.pcoa$vectors[,1]), y = invert_y*(df.pcoa$vectors[,2]), col= if(!missing(highlight)){col_vec}else{if(!missing(transp)){transp_vec}else{as.numeric(as.factor(df_meta[[type]]))}}, lwd = 3, cex =if(!missing(size)){size}else{1}, pch= if(all_O == TRUE){1}else{as.numeric(as.factor(df_meta[[type2]]))}, main = study, xlab = paste0("PC1 (", axis1_rel_eigen, "%)"), ylab = paste0("PC2 (", axis2_rel_eigen, "%)"), asp = ASP)
  }

  if(legendyn == TRUE){
    #adds a legend in the topright corner of the plot with adjustable inset.
    if(!missing(type2)){
      graphics::legend("bottomright", inset = c(inset, 0), legend = sort(unique(df_meta[[type2]])), col = "black", pch = as.numeric(as.factor(df_meta[[type2]])), bty = "n", cex = 1.5, y.intersp = Legend_Yspace)
      graphics::legend("topright", inset = c(inset,0), legend = sort(unique(df_meta[[type]])), col = if(red_blue == TRUE){cols[as.numeric(as.factor(sort(unique(X))))]}else{1:length(unique(X))}, pch = 19, bty = "n", cex = 1.5, y.intersp = Legend_Yspace)
    } else{
      graphics::legend("topright", inset = c(inset,0), legend = sort(unique(df_meta[[type]])), col = if(red_blue == TRUE){cols[as.numeric(as.factor(sort(unique(X))))]}else{1:length(unique(X))}, pch = 19, bty = "n", cex = 1.5, y.intersp = Legend_Yspace)
    }
  }
  #resets palette to original
  grDevices::palette(pal)

  df.pcoa
}


#' Blaster
#'
#' Sequence queries are crosschecked to a local 16S rRNA gene database. Subject titles returned are the names of the database sequences which match the query sequences.
#' This function needs to run in the top level project directory i.e. /Analysis_Ready_Files.
#'
#' @param seqs a character vector of 16S rRNA gene sequences to be cross checked against a local 16S rRNA gene database
#' @param names a character vector equal in length to `seqs` which denotes the query name for each sequence.
#' @param alignments an integer value designating how many alignments should be returned per sequence query.
#'
#' @return A summary table listing the query names, Expect value, Bitscore, Subject TAX ID, Percent Identical, Query Coverage, and Subject Title.
#' @export
#'
#' @examples
#' Blast_results <- blaster(Theis_merged_clean$X[1:2], names = c("ASV1", "ASV2"), alignments = 5)
blaster <- function(seqs, names, alignments = 1) {
  #store sequence data in fasta formatted object
  seqinr::write.fasta(sequences = as.list(seqs), names = names, file.out = "blaster.fasta")
  #Take fasta object and run it against 16S rRNA gene database.
  system(paste0("blastn -query ", paste0(here::here(), "/blaster.fasta -num_alignments "), alignments, " -num_threads 3 -out ", paste0(here::here(), "/Blasted.tsv -outfmt \"6 qseqid evalue bitscore staxids pident qcovs stitle\" -db "), paste0(here::here(), "/dada2tools/extdata/16S_ribosomal_RNA")))
  results <- utils::read.table(file = "Blasted.tsv", sep = "\t", stringsAsFactors = FALSE)
  colnames(results) <- c("Query Seq ID", "Expect value", "Bitscore", "Subject TAX ID", "Percent Identical", "Query Coverage", "Subject Title")
  unlink(c("blaster.fasta", "Blasted.tsv"))
  results
}

#' Color Replication
#'
#' Auxillary function to `beta_div()` which allows for modification of plotted colors. In particular, it allows for certain points to be highlighted while others are colored grey or left blank.
#'
#'
#' @param colors color palette as character vector of colors to be used.
#' @param full_set factor vector indicating sample group for each sample.
#' @param highlight character vector of all samples to be colored normally
#' @param highlight_bool boolean vector indicating which samples are in `highlight`.
#' @param nonhgrey boolean vector indicating whether non-highlighted samples should be grey (`TRUE`) or left blank (`FALSE`)
#'
#' @return a new character vector containing the new colors by sample.
#'
color_rep <- function(colors, full_set, highlight, highlight_bool, nonhgrey){
  uniques <- unique(full_set)
  color_vector <- c()
  if (missing(highlight)){
    for(i in 1:length(full_set)){
      color_vector <- append(color_vector, colors[full_set[i] == uniques])
    }
  } else {
    for(i in 1:length(full_set)){
      color_vector <- append(color_vector, ifelse(highlight_bool[i] == TRUE, colors[match(full_set[i], highlight)], NA))
    }
  }
  if(nonhgrey == TRUE){
    color_vector[is.na(color_vector)] <- "#D3D3D3"
  }
  color_vector
}

#' Convert Merged Dataset
#'
#' Converts a merged dataset taxonomy columns to character and ASV/OTU counts to numeric to ensure that all columns are of the proper type.
#'
#' @param merged Column bound taxonomy with ASVs/OTUs on rows to ASV/OTU counts with samples on columns.
#' @param char_range numeric vector giving column range for taxonomy. Default is `1:8`
#' @param num_range numeric vector giving range for ASV/OTU count table. Default is all other columns to the right of taxonomy.
#'
#' @return merged dataset with taxonomy columns as character and ASV/OTU counts per sample as numeric.
#' @export
#'
#' @examples
#' Theis_converted <- convert_merged(merged = Theis_merged, char_range = 1:8)
convert_merged <- function(merged, char_range = 1:8, num_range){
  if(missing(num_range)){
    num_range <- (max(char_range)+1):ncol(merged)
  }
  for(i in 1:ncol(merged)){
    if(i %in% char_range & is.factor(merged[,i])){
      merged[,i] <- as.character(merged[,i])
    }
    if(i %in% num_range & is.factor(merged[,i])){
      merged[,i] <- as.numeric(as.character(merged[,i]))
    }
    if(i %in% num_range & is.character(merged[,i])){
      merged[,i] <- as.numeric(merged[,i])
    }
  }
  merged
}

#' Dadaset Clean
#'
#' Removes taxonomy which is not classified as bacteria, or down to the phylum level, mitochondria, and chloroplast. Also removes samples which are below a designated read threshold after filtering.
#'
#' @param df a merged dataframe i.e. taxonomy on rows and samples on columns.
#' @param read_thresh a numeric vector determining the number of reads per sample below which a sample is removed from the dataset. This cutoff is applied after filtering taxonomy. The default is 100.
#'
#' @return a merged dataframe which is free of unclassified taxonomy at the phylum level, and non-bacterial mitochondrial, and chloroplast sequences.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' Theis_merged_clean <- dadaset_clean(df = Theis_merged, read_thresh = 100)
dadaset_clean <- function(df, read_thresh = 100) {
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

  df2_numeric <- df2 %>% dplyr::select_if(., is.numeric)

  #dplyr way of removing samples with less than read_thresh reads.
  #first calculate sums of all columns
  df2_num_pass <- df2_numeric %>% dplyr::select_if(function(col) sum(col) >= read_thresh)
  elim <- df2_numeric %>% dplyr::select_if(function(col) sum(col) < read_thresh)
  elim_names <- names(elim)
  cleaned <- cbind(df2 %>% dplyr::select_if(., is.character), df2_num_pass)

  for(i in 2:6){
    cleaned[,(i + 1)] <- purrr::map2_chr(cleaned[,i], cleaned[,(i+1)], unclassified_replace)
  }

  #Outputs changes made to dataset.
  cat("Originally there were", nrow_taxa_ori, "taxa.", nrow_taxa_ori - nrow_taxa_bac, "were not bacterial or classified at phylum level.", nrow_taxa_bac - nrow_taxa_Mit, "were mitochondrial.", nrow_taxa_Mit - nrow_taxa_Chl, "were chloroplast.", nrow_taxa_Chl, "taxa remain. ")
  cat(length(elim_names), ifelse(length(elim_names) > 1 | length(elim_names) < 1, "samples", "sample"), "eliminated with less than", read_thresh, "reads:", elim_names)

  #at this point the dataset is cleaned and ready to move to preparation for other analyses.
  #sample type counts should be made in a list with list_types or manually for heatmap_prep

  cleaned <- as.data.frame(cleaned)
  cleaned
}

#' Preparation for DECONTAM
#'
#' Prepares a merged_clean dataset that has been run through dadaset_clean() for processing with DECONTAM
#'
#' @param df merged_clean dataset with taxonomy and ASV table combined. ASVs should be on rows.
#' @param type the name of a column in `meta` where each sample is either "sample" or "control" (as in Technical control).
#' @param meta a metadata data frame with column `type`
#' @param sample_col character vector of column name in `meta` which indicates sample names. Also should match column names in merged dataset.
#'
#' @return returns a phyloseq object which can then be used in decontam_histo_prev_plots()
#' @export
#'
#' @examples
#' Theis_META_order$Deco_type <- dplyr::case_when(Theis_META_order$Type == "Placenta" ~ "sample",
#' Theis_META_order$Type == "Technical Control" ~ "control")
#' DECO_prep <- decontam_prep(df = Theis_merged_clean,
#' meta = Theis_META_order,
#' type = "Deco_type",
#' sample_col = "Run")
decontam_prep <- function(df, meta, type, sample_col){

  #ensures data is of type dataframe
  df <- as.data.frame(df)
  #selects for numeric columns only
  df_num <- dplyr::select_if(df, is.numeric)

  #creates an "OTU Table" from the numeric columns of the merged_clean dataset
  OTU = phyloseq::otu_table(object = df_num, taxa_are_rows = TRUE)

  #taxa table created from taxa classifications Kingdom through Genus
  TAX = phyloseq::tax_table(as.matrix(df[,2:7]))

  #sample data phyloseq object SAM created with Type column populated by type1 chr vector.
  SAM = phyloseq::sample_data(data.frame(data = meta[[type]], row.names = meta[[sample_col]]))

  #combine all phyloseq tables into one phyloseq object
  physeq <- phyloseq::phyloseq(OTU, TAX, SAM)
  physeq
}

#' Generate DECONTAM Histogram and Prevalence Plots
#'
#' To be used on a Phyloseq object to generate a histogram of DECONTAM scores and a prevalence plot.
#'
#' @param physeq phyloseq object with OTU, TAX, and SAM components
#' @param thresh DECONTAM score threshold. Default is 0.5
#' @param study_name character vector name of the study/title for plot.
#'
#' @return returns a list object of the histogram and prevalence plots.
#' @export
#'
#' @examples
#' Theis_META_order$Deco_type <- dplyr::case_when(Theis_META_order$Type == "Placenta" ~ "sample",
#' Theis_META_order$Type == "Technical Control" ~ "control")
#' DECO_prep <- decontam_prep(df = Theis_merged_clean,
#' meta = Theis_META_order,
#' type = "Deco_type",
#' sample_col = "Run")
#' decontam_plots <- decontam_histo_prev_plots(physeq = DECO_prep,
#' thresh = 0.1,
#' study_name = "Theis")
decontam_histo_prev_plots <- function(physeq, thresh = 0.5, study_name){

  p <- prev <- physeq.neg <- physeq.pos <- contaminant <- NULL
  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$data == "control"

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
  physeq.pa.neg <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$data == "control", physeq.pa)

  #subsets samples which are actual samples into phyloseq object physeq.pa.pos
  physeq.pa.pos <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$data == "sample", physeq.pa)

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
#' Final function for using DECONTAM.
#'
#' @param df dataframe from which the phyloseq object was generated with decontam_prep()
#' @param physeq phyloseq object created from running decontam_prep()
#' @param thresh numeric vector of threshold chosen based on histogram and prevalence plots.
#'
#' @return gives a list of dataframes split into true taxa, contaminants, and the chosen threshold.
#' @export
#'
#' @examples
#' Theis_META_order$Deco_type <- dplyr::case_when(Theis_META_order$Type == "Placenta" ~ "sample",
#' Theis_META_order$Type == "Technical Control" ~ "control")
#' DECO_prep <- decontam_prep(df = Theis_merged_clean,
#' meta = Theis_META_order, type = "Deco_type", sample_col = "Run")
#' Theis_decontaminated <- decontaminate(df = Theis_merged_clean, physeq = DECO_prep, thresh = 0.1)
decontaminate <- function(df, physeq, thresh) {

  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$data == "control"

  #isNotContaminant() is chosen here over isContaminant() since samples are presumed low biomass and the majority of taxa are assumed to be contaminants.
  #Keep in mind that Trues are not contaminants
  contamdf.prev <- decontam::isNotContaminant(physeq, neg = "is.neg", detailed = TRUE, threshold = thresh)

  #subsetting contaminants and 'true taxa'
  true_taxa <- subset(df, contamdf.prev$not.contaminant)
  contaminants <- subset(df, !contamdf.prev$not.contaminant)

  #returns list of dataframes true taxa and contaminants along with the threshold used.
  list(TrueTaxa = true_taxa, Contaminants = contaminants, Threshold = thresh)
}

#' DePhyloseqize
#'
#' Takes a phyloseq object and converts it back into a merged taxonomy and ASV/OTU table.
#'
#' @param phyloseq_obj the phyloseq object to be converted
#'
#' @return a merged taxonomy and ASV/OTU table as a data frame.
#' @export
#'
#' @examples
#' phy <- phyloseqize(merged_df = Theis_merged_clean, taxa_as_rows = TRUE, keep_ASV_nums = TRUE)
#' merged_dephy <- dephy(phyloseq_obj = phy)
dephy <- function(phyloseq_obj){
  ASV <- phyloseq::otu_table(phyloseq_obj)
  ASV <- as.data.frame(ASV)
  TAX <- phyloseq::tax_table(phyloseq_obj)
  TAX <- as.data.frame(TAX)
  if(nrow(TAX) == nrow(ASV)){
    merged <- cbind(TAX, ASV)
  } else {
    merged <- cbind(TAX, as.data.frame(t(ASV)))
  }
  merged
}

#' Get Top ASVs
#'
#' A simple function to return the names of the top ASVs in a dataset.
#'
#' @param df numeric ASV/OTU table.
#' @param number numeric vector indicating how many top ASVs to return.
#' @param ASVs_on_Rows boolean vector indicating whether ASVs in `df` are on rows (`TRUE`) or columns (`FALSE`)
#'
#' @return character vector of the top ASVs/OTUs in a dataset.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' Theis_topASVs <- get_top_ASVs(df = Theis_rare_ASV, number = 10, ASVs_on_Rows = TRUE)
get_top_ASVs <- function(df, number, ASVs_on_Rows = TRUE){
  if(ASVs_on_Rows == FALSE){
    df <- as.data.frame(t(df))
  }
  if(dim(df)[2] == 1){
    df$dummy_col <- 0
  }
  df <- df[order(rowSums(df), decreasing = TRUE),]
  top_ASVs <- rownames(df) %>% utils::head(., number)
  top_ASVs
}

#' Heatmap Preparation
#'
#' This function takes an ASV/OTU table, taxonomy table, and metadata file, and prepares the ASVs for heatmapping.
#' In particular it groups all other ASVs not selected into an Other column.
#' There are multiple options for ASV selection: `mean_ab_by_type`, `mean_ab_cutoff`, `select_ASVs`, and `taxa_cutoff`.
#' Only one option can be used at a time. They will be evaluated in the order listed so if `mean_ab_by_type` is passed a character vector, no other selection methods will be used.
#'
#' @param df an ASV/OTU table dataframe with ASVs/OTUs as rows and samples on columns. Row/Column names should be ASV number/Sample name respectively.
#' @param df_tax an associated taxonomy table dataframe. Needs tomatch ASVs on rows in `df`
#' @param df_meta an associated metadata dataframe with optional column for ASV selection based on type column
#' @param class_col character vector of column in `df_tax` to designate which level of classification to name ASVs
#' @param taxa_cutoff integer number of top ASVs to select. Default is 0. Only used if `!= 0`.
#' @param mean_ab_cutoff boolean vector to select ASVs greater than 1% relative abundance across the entire dataset. Default is `TRUE`
#' @param select_ASVs character vector of ASVs to select by ASV number, which should be rownames in `df_tax`.
#' @param mean_ab_by_type no default. If passed a character vector of a column in `df_meta`, will be used to split samples into types before calculating mean percent relative abundance.
#'
#' @return A dataframe with columns named by ASV numbers and selected taxonomic classification and samples on rows.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' library(dplyr)
#' Theis_META_P <- Theis_rare_META %>% filter(Type == "Placenta")
#' Theis_agg <- agglomerate(df = Theis_rare_ASV, df_tax = Theis_rare_TAX, agg_class = "Genus")
#' Theis_agg_tax <- data.frame(Genus = rownames(Theis_agg), row.names = rownames(Theis_agg))
#' Theis_agg_P <- Theis_agg %>% select_if(colnames(.) %in% Theis_META_P$Run)
#' Theis_agg_tax_P <- Theis_agg_tax %>% filter(rownames(.) %in% rownames(Theis_agg_P))
#' Theis5 <- get_top_ASVs(Theis_agg_P, 5)
#' Theis_preheatbar <- heatmap_prep(df = Theis_agg_P,
#' df_tax = Theis_agg_tax_P,
#' class_col = "Genus",
#' df_meta = Theis_META_P,
#' select_ASVs = Theis5,
#' mean_ab_cutoff = FALSE)
heatmap_prep <- function (df, df_tax, df_meta, class_col, taxa_cutoff = 0, mean_ab_cutoff = TRUE, select_ASVs, mean_ab_by_type){
  n2 <- c()
  dt_preheat <- t(df)
  dt_preheat <- as.data.frame(dt_preheat)
  colnames(dt_preheat) <- paste(colnames(dt_preheat),  unlist(df_tax[[class_col]]), sep = "-")
  dt_preheat1 <- as.data.frame(sapply(dt_preheat, function(x) as.numeric(as.character(x))))
  if(dim(dt_preheat1)[2] == 1){
    dt_preheat1 <- as.data.frame(t(dt_preheat1))
  }
  rownames(dt_preheat1) <- rownames(dt_preheat)
  colnames(dt_preheat1) <- make.unique(colnames(dt_preheat1))
  if(!missing(mean_ab_by_type)){
    for(i in 1:length(unique(df_meta[[mean_ab_by_type]]))){
      data.prop <- dt_preheat1/rowSums(dt_preheat1)
      data.prop_sub <- data.prop %>% dplyr::filter(df_meta[[mean_ab_by_type]] == unique(df_meta[[mean_ab_by_type]])[i])
      means <- apply(data.prop_sub, 2, mean)
      n2 <- append(n2, names(which(means > 0.01)))
    }
    dt_preheat1.ab <- dt_preheat1[, which(names(data.prop) %in% n2)]
    preheat_other <- dt_preheat1[, -which(names(data.prop) %in% n2)]
    dt_preheat1.ab$Other <- rowSums(preheat_other)
    dt_preheat1.ab <- row_sample_elim(dt_preheat1.ab, thresh = 1)
    return(dt_preheat1.ab)
  }
  if (mean_ab_cutoff == TRUE) {
    data.prop <- dt_preheat1/rowSums(dt_preheat1)
    mean_abundance <- apply(data.prop, 2, mean)
    n2 <- names(which(mean_abundance < 0.01))
    dt_preheat1.ab <- dt_preheat1[, -which(names(data.prop) %in% n2)]
    preheat_other <- dt_preheat1[, which(names(data.prop) %in% n2)]
    dt_preheat1.ab$Other <- rowSums(preheat_other)
    dt_preheat1.ab <- row_sample_elim(dt_preheat1.ab, thresh = 1)
    return(dt_preheat1.ab)
  }
  if (!missing(select_ASVs)){
    dt_preheat1.ab <- dt_preheat1 %>% dplyr::select_if(rownames(df_tax) %in% select_ASVs)
    others <- dt_preheat1 %>% dplyr::select_if(rownames(df_tax) %ni% select_ASVs)
    dt_preheat1.ab$Other <- rowSums(others)
    return(dt_preheat1.ab)
  }
  if (taxa_cutoff != 0) {
    dt_rank <- rank(-colSums(dt_preheat1), ties.method = "first")
    dt_preheat1 <- as.data.frame(t(subset(t(dt_preheat1),
                                          dt_rank <= taxa_cutoff)))
    dt_preheat1 <- row_sample_elim(dt_preheat1, thresh = 1)
    return(dt_preheat1)
  }
}

#' Heatmapping
#'
#' Produces a heatmap which is grouped by a column within the associated metadata dataframe. Heatmap takes into account ASVs not plotted in 'Other' category.
#' ASV/OTUs are grouped on rows by k-means clustering.
#'
#' @param df a dataframe with samples on rows and ASV/OTUs on columns. Both should be row/column names respectively.
#' @param df_meta an associated metadata dataframe with samples on rows and a column to group samples.
#' @param types_col character vector of column name within metadata dataframe to group samples by.
#' @param scalecolor currently only one color scale. "red-blue" is the default.
#' @param title character vector for the title of the heatmap if desired.
#' @param by_percent boolean vector to determine if heatmap should be visualized based on percent relative abundance.
#' @param show_other boolean vector. If only a subset of ASVs/OTUs are to be plotted, should ASVs/OTUs grouped into "other" catogory by `heatmap_prep()` be plotted as well?
#'
#' @return a heatmap clustered by sample type designated by column within metadata dataframe.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' library(dplyr)
#' Theis_META_P <- Theis_rare_META %>% filter(Type == "Placenta")
#' Theis_agg <- agglomerate(df = Theis_rare_ASV, df_tax = Theis_rare_TAX, agg_class = "Genus")
#' Theis_agg_tax <- data.frame(Genus = rownames(Theis_agg), row.names = rownames(Theis_agg))
#' Theis_agg_P <- Theis_agg %>% select_if(colnames(.) %in% Theis_META_P$Run)
#' Theis_agg_tax_P <- Theis_agg_tax %>% filter(rownames(.) %in% rownames(Theis_agg_P))
#' Theis5 <- get_top_ASVs(Theis_agg_P, 5)
#' Theis_preheatbar <- heatmap_prep(df = Theis_agg_P,
#' df_tax = Theis_agg_tax_P,
#' class_col = "Genus",
#' df_meta = Theis_META_P,
#' select_ASVs = Theis5,
#' mean_ab_cutoff = FALSE)
heatmapping <- function(df, df_meta, types_col, scalecolor = "red-blue", title, by_percent = TRUE, show_other = FALSE) {
  #converts dataset to proportionalized dataset (taxa percentage of sample).
  data.prop <- (df/rowSums(df))*100
  if(show_other == FALSE){
    data.prop <- data.prop %>% dplyr::select_if(colnames(.) != "Other")
    df <- df %>% dplyr::select_if(colnames(.) != "Other")
  }

  #instantiate factortypes and repeat type names based on values in named list types.
  factortypes <- as.factor(df_meta[[types_col]])

  if(by_percent == TRUE) {
    if(scalecolor == "red-blue"){
      col_fun <- circlize::colorRamp2(c(0, 25, 75, 100), c("white", "blue", "yellow", "red"))
    }

    #for plotting based on percent relative abundance
    ht = ComplexHeatmap::Heatmap(matrix = as.matrix(t(data.prop)),
                                 column_title = sort(unique(df_meta[[types_col]])),
                                 show_column_names = FALSE,
                                 show_column_dend = FALSE,
                                 col = col_fun,
                                 name = "Percentage of Sample",
                                 cluster_columns = TRUE,
                                 cluster_column_slices = FALSE,
                                 row_names_gp = grid::gpar(fontface = 3),
                                 row_names_side = "left",
                                 row_dend_side = "right",
                                 row_dend_width = grid::unit(1.5, "cm"),
                                 heatmap_legend_param = list(at =c(0, 50, 100), labels = c(0, 50, 100)),
                                 column_split = factortypes,
                                 border = TRUE)

    ComplexHeatmap::draw(ht, column_title = title)
  }
  else{
    if(scalecolor == "red-blue"){
      mx <- max(df)
      col_fun <- circlize::colorRamp2(c(0, (mx*1/4), (mx*3/4), mx), c("white", "blue", "yellow", "red"))
    }
    #for plotting based on total reads
    ht = ComplexHeatmap::Heatmap(matrix = as.matrix(t(df)),
                                 column_title = sort(unique(df_meta[[types_col]])),
                                 show_column_names = FALSE,
                                 show_column_dend = FALSE,
                                 col = col_fun,
                                 name = "Total Reads",
                                 cluster_columns = TRUE,
                                 cluster_column_slices = FALSE,
                                 row_names_gp = grid::gpar(fontface = 3),
                                 row_names_side = "left",
                                 row_dend_side = "right",
                                 row_dend_width = grid::unit(1.5, "cm"),
                                 column_split = factortypes,
                                 border = TRUE)
    ComplexHeatmap::draw(ht, column_title = title)
  }
}

#' Jensen-Shannon Sample Clustering by Group
#'
#' Calculates the average Jensen-Shannon Divergence within sample groups. Can be used to determine degree of similarity.
#'
#' @param df an ASV/OTU table with samples as column names and ASV/OTUs as row names. The entire dataframe should be ASV/OTU counts, not taxonomy.
#' @param meta an associated metadata file which is ordered by the order of samples in the ASV/OTU table, and which has a column to group samples by.
#' @param Grouping_Col character vector of the column name by which to group by.
#'
#' @return A dataframe with group names, and mean Jensen-Shannon divergences.
#' @export
#'
#' @examples
#' JS_clusters <- JS_group_div(df = t(Theis_rare_ASV),
#' meta = Theis_rare_META,
#' Grouping_Col = "Type")
JS_group_div <- function(df, meta, Grouping_Col){
  JSD_means <- c()
  unique_dups <- unique(meta[[Grouping_Col]][duplicated(meta[[Grouping_Col]])])
  #needs to be iterative.
  grepl_result <- grepl(pattern = paste0("^", unique_dups[[1]], "$"), meta[[Grouping_Col]])
  #be able to isolate each group.
  for(i in 2:length(unique_dups)){
    grepl_result <- append(grepl_result, grepl(pattern = paste0("^", unique_dups[[i]], "$"), meta[[Grouping_Col]]))
  }
  for (i in 1:length(unique_dups)){
    df_group <- df[grepl_result[(1+(dim(df)[1]*(i-1))):(dim(df)[1]+(dim(df)[1]*(i-1)))],1:ncol(df)]
    JSD_mat <- philentropy::JSD(as.matrix(df_group), test.na = FALSE, unit = "log", est.prob = "empirical")
    if(length(dim(JSD_mat)[2]) == 1){
      col_len <- dim(JSD_mat)[2]
      skip_TF <- c(rep(c(FALSE, rep(TRUE, col_len)), col_len-1), FALSE)
      JSD_vec <- as.vector(JSD_mat)
      JSD_vec <- JSD_vec[skip_TF]
    } else{
      JSD_vec <- JSD_mat
    }
    JSD_means <- append(JSD_means, mean(JSD_vec))
  }
  JSD_out <- data.frame(Group = unique_dups, Means = JSD_means)
  JSD_out
}

#allows "not in" syntax
'%ni%' <- Negate('%in%')

#' Phyloseqize
#'
#' Convert data in R to a phyloseq object
#'
#' @param merged_df a numeric dataframe of ASV sample counts. Can include taxonomy as well. `merged_df` and `tax_df` row names should be identical.
#' @param tax_df optional, if taxonomy is not included in `merged_df` include taxonomy table. For best results, create a taxonomy data frame with ASV row names and taxonomic levels as column names. Then convert to a matrix to use as `tax_df`.
#' @param taxa_as_rows boolean vector indicating whether ASVs are on rows or on columns. Default is `TRUE`
#' @param keep_ASV_nums boolean vector indicating whether to keep row names of ASV/OTU table.
#'
#' @return a phyloseq object with ASV/OTU count data and associated taxonomy.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' phy <- phyloseqize(merged_df = Theis_merged_clean, taxa_as_rows = TRUE, keep_ASV_nums = TRUE)
phyloseqize <- function(merged_df, tax_df, taxa_as_rows = TRUE, keep_ASV_nums = TRUE) {
  if(sum(sapply(merged_df, is.character)) == 8){
    pre_otu <- as.data.frame(merged_df[,9:ncol(merged_df)])
  }
  if(!missing(tax_df)){
    pre_otu <- merged_df
  }
  if(keep_ASV_nums == FALSE){
    rownames(pre_otu) <- paste0("ASV", 1:nrow(pre_otu))
  }

  OTU <- phyloseq::otu_table(pre_otu, taxa_are_rows = taxa_as_rows)

  if(sum(sapply(merged_df, is.character)) == 8){
    pre_tax <- as.data.frame(merged_df[,1:8])
  } else {
    pre_tax <- tax_df
  }
  if(keep_ASV_nums == FALSE){
    rownames(pre_tax) <- paste0("ASV", 1:nrow(pre_tax))
  }
  pre_tax <- as.matrix(pre_tax)
  TAX <- phyloseq::tax_table(pre_tax)

  phy_object <- phyloseq::phyloseq(OTU, TAX)
  phy_object
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
#' Theis_preheat_pass <- row_sample_elim(Theis_preheat_ab, thresh = 1)
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

#' Top ASVs Above Cutoff
#'
#' Determines ASVs above a certain cutoff and returns taxonomy and abundance data.
#'
#' @param ASV numeric dataframe with ASVs on rows and samples on columns.
#' @param TAX taxonomy matching ASV/OTU dataset with taxa on rows.
#' @param META a metadata dataframe associated with samples on rows
#' @param cutoff percent relative abundance cutoff in decimal form above which ASVs will be listed.
#' @param top integer of ASVs above cutoff to return
#' @param ASVs_on_Rows boolean vector indicating if ASVs are on rows in `ASV`
#' @param study character vector name of study/title.
#' @param subject_col character vector of column name in `META` to group samples by.
#' @param subjects character vector of groups within `subject_col` to evaluate.
#'
#' @return a dataframe with ASV name, abundance in dataset, genus and species classifications, and exact sequence.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' Theis_top5 <- top_ASVs_above_cutoff(ASV = Theis_rare_ASV,
#' TAX = Theis_rare_TAX,
#' cutoff = 0.01,
#' top = 5,
#' ASVs_on_Rows = TRUE,
#' study = "Theis")
top_ASVs_above_cutoff <- function(ASV, TAX, META, cutoff, top, ASVs_on_Rows = TRUE, study, subject_col, subjects){
  ASVs_Greater_than <- c()
  if(ASVs_on_Rows == TRUE){
    ASV <- as.data.frame(t(ASV))
  }
  if(!missing(subject_col)){
    for (i in 1:length(subjects)){
      ASV_subject <- as.data.frame(ASV) %>% dplyr::filter(META[[subject_col]] == subjects[i])
      ASV_colsums <- colSums(ASV_subject)
      ASV_colsums_pct <- ASV_colsums/sum(ASV_colsums)
      ASV_cutoff <- ASV_colsums_pct[ASV_colsums_pct > cutoff]
      ASVs_Greater_than <- append(ASVs_Greater_than, ASV_cutoff)
      if(missing(top)){
        top <- length(ASVs_Greater_than)
      }
    }
  }
  if(missing(subject_col)){
    ASV_colsums <- colSums(ASV)
    ASV_colsums_pct <- ASV_colsums/sum(ASV_colsums)
    ASV_cutoff <- ASV_colsums_pct[ASV_colsums_pct > cutoff]
    ASVs_Greater_than <- append(ASVs_Greater_than, ASV_cutoff)
    if(missing(top)){
      top <- length(ASVs_Greater_than)
    }
    ASVs_Greater_than <- ASVs_Greater_than %>% sort(., decreasing = TRUE) %>% utils::head(., top)
  }

  df <- data.frame(ASV = names(ASVs_Greater_than), Rel_Abund = unlist(ASVs_Greater_than))
  TAX_filt <- TAX %>% dplyr::filter(rownames(.) %in% names(ASVs_Greater_than)) %>% dplyr::select(c("X", "Genus", "Species"))
  TAX_filt$ASV <- rownames(TAX_filt)
  df_merged <- merge(df, TAX_filt, by = "ASV")

  if(dim(df_merged)[1] != top){
    current_rows <- dim(df_merged)[1]
    missing_rows <- (top - current_rows)
    for(i in 1:missing_rows){
      df_merged <- rbind(df_merged, data.frame(ASV = NA, Rel_Abund = 0, X = NA, Genus = NA, Species = NA))
    }
  }

  df_merged <- cbind(df_merged, data.frame(Study = rep(study, nrow(df_merged))))
  df_merged <- df_merged[order(df_merged$Rel_Abund, decreasing = TRUE),]
  df_merged
}

#' Weighted Average Genus Labels
#'
#' Plots weighted average genus labels over a beta diversity plot generated by `beta_div()`.
#'
#' @param Plot_object object returned from `beta_div()`
#' @param ASV_table numeric ASV_table dataframe used to produce beta diversity plot. ASVs can be on rows or columns but should be indicated with `ASVs_on_rows`
#' @param TAX_table dataframe of taxonomy matching ASVs in `ASV_table`. Taxa should be on rows.
#' @param ASVs_to_plot a character vector of ASV numbers to plot
#' @param Tax_col chracter vector of column name in `TAX_table` to use for labels if `ASV_with_Genus` is set to `FALSE`
#' @param ASV_with_Genus a boolean vector to determine if ASV numbers should also be plotted along with genus classifications of ASVs.
#' @param ASVs_on_rows a boolean vector indicating if ASVs are on rows in `ASV_table`
#' @param color a character vector of color to plot points.
#' @param Greater_than_1per optional boolean vector. If `TRUE`, calculates ASVs to plot based on a cutoff of greater than 1% mean relative abundance.
#' @param invert_x optional boolean vector indicating whether to invert points over y-axis. Default is FALSE.
#' @param invert_y optional boolean vector indicating whether to invert points over x-axis. Default is FALSE
#'
#' @return Does not return an object, only adds points and labels to a current beta diversity plot.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' Theis_topASVs <- get_top_ASVs(Theis_rare_ASV, number = 10)
#' Theis_Beta_plot_PTCWA <- beta_div(df = Theis_rare_ASV,
#' df_meta = Theis_rare_META,
#' type = "Type",
#' study = "Example Plot",
#' all_O = TRUE,
#' taxa_on_rows = TRUE,
#' size = 2,
#' legendyn = FALSE)
#' WA_labels(Plot_object = Theis_Beta_plot_PTCWA,
#' ASV_table = Theis_rare_ASV,
#' TAX_table = Theis_rare_TAX,
#' ASVs_on_rows = TRUE,
#' ASV_with_Genus = FALSE,
#' Tax_col = "Genus",
#' color = "grey",
#' Greater_than_1per = FALSE,
#' ASVs_to_plot = Theis_topASVs)
WA_labels <- function(Plot_object, ASV_table, TAX_table, ASVs_to_plot, Tax_col, ASV_with_Genus = TRUE, ASVs_on_rows = TRUE, color, Greater_than_1per = FALSE, invert_x = FALSE, invert_y = FALSE){
  if(ASVs_on_rows == TRUE){
    ASV_table <- t(ASV_table)
  }

  if(Greater_than_1per == TRUE){
    ASV_df <- t(ASV_table)
    ASV_rel <- ASV_df/colSums(ASV_df)
    mean_abundance <- rowMeans(ASV_rel)
    ASV_names <- names(which(mean_abundance > 0.01))
    ASVs_to_plot <- ASV_names
  }


  WAscores <- vegan::wascores(x = Plot_object$vectors, w = ASV_table)
  WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
  if(invert_x == TRUE){
    WAscores_df[,1] <- WAscores_df[,1]*-1
  }
  if(invert_y == TRUE){
    WAscores_df[,2] <- WAscores_df[,2]*-1
  }
  WAscores_select <- WAscores_df %>% dplyr::filter(rownames(.) %in% ASVs_to_plot)
  TAX_select <- TAX_table %>% dplyr::filter(rownames(.) %in% ASVs_to_plot)
  if(ASV_with_Genus == TRUE){
    WAscores_select$labels <- paste(rownames(TAX_select), TAX_select$Genus, sep = "-")
  } else {
    WAscores_select$labels <- make.unique(TAX_select[[Tax_col]])
  }
  graphics::text(WAscores_select[,1:2], labels = WAscores_select$labels, col = "grey50", pos = 4, font = 3)
  graphics::points(WAscores_select[,1:2], pch = 18, col = color)
}

#' Unclassified Replace
#'
#' Auxillary function to `dadaset_clean()`. Used in `dadaset_clean()` to replace NA's with next highest taxonomy.
#'
#' @param Col_upper dataframe column corresponding to higher taxonomic rank.
#' @param Col_lower dataframe column corresponding to lower taxonomic rank.
#' @param ... left blank.
#'
#' @return modified lower taxonomic rank
#'
#' @examples
#' Theis_merged$Genus <- dada2tools:::unclassified_replace(Col_upper = Theis_merged$Family,
#' Col_lower = Theis_merged$Genus)
unclassified_replace <- function(Col_upper, Col_lower, ...) {
  #If cell in lower level is NA, makes replacement, otherwise it leaves the cell alone.
  ifelse(is.na(Col_lower) == TRUE,
         Col_lower <- Col_upper,
         Col_lower <- Col_lower)
}
