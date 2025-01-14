# Libraries
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggrepel)
library(DT)

#' Retrieve a file from a linked Google Drive profile
#'
#' @param filename Name of the file to be retrieved with string and extension.
#'
#' @returns
#' @export
#'
#' @examples
retrieve_from_google <- function(filename, extension = "tsv") {
    asm_sheet <- googledrive::drive_find(pattern = "Assembly Stats")
    googledrive::drive_download(file = asm_sheet, path = filename, type = extension, overwrite = TRUE)   #authenticate as kcollier@alaska.edu
  return(filename)
}

#' Convert genome assembly sheet to data frame
#'
#' @param filename Name of the file to be wrangled.
#'
#' @returns
#' @export
#'
#' @examples
wrangle_asm_sheet <- function(filename){
  character_cols <- c("assembly_name", "genus", "species", "taxon", "common_name", "project")
  integer_cols <- c("num_contigs", "n50", "l50", "n90", "l90")
  bool_cols <- c("has_ont", "has_pacbio", "has_illumina")
  
  # Some very large columns ("raw_input_data_bp", "trimmed_input_data_bp", "total_length") are treated as percents.
  # This is because they exceed R's maximum allowed integer range, resulting in NA values if coerced to integer
  percent_cols <- c("raw_input_data_bp", "trimmed_input_data_bp", "reads_over_50000bp_bp", "total_length",
                    "coverage", "percent_gc_content", "percent_genome_in_fragments_1_000_000_bp", "percent_complete_busco")
  asm_tib <- janitor::clean_names(readr::read_tsv(filename, show_col_types = FALSE))
  asm_tib <- asm_tib %>%
    select(all_of(c(character_cols, integer_cols, percent_cols, bool_cols)))
  
  asm_tib <- asm_tib %>% na.omit() %>%
    mutate(across(all_of(character_cols), as.character)) %>%
    mutate(across(all_of(integer_cols), as.integer)) %>% 
    mutate(across(all_of(percent_cols), as.numeric)) %>%
    mutate(across(all_of(bool_cols), as.logical)) %>%
    mutate(percent_longrds = round((reads_over_50000bp_bp/raw_input_data_bp)*100,2))
  
  # replace long project name with shorter abbreviation
  asm_tib$project <- str_replace_all(asm_tib$project, "Abu Dhabi Environment Agency", "ADEA")
  
  return(asm_tib)
}

#GET SHEET
asm_sheet <- "Assembly Stats - Assembly stats.tsv"
if(!file.exists(asm_sheet)) {
  asm_sheet <- retrieve_from_google(asm_sheet) # this retrieves the file if not present
}
asm_tib <- wrangle_asm_sheet(asm_sheet)


#' Display a tibble or data frame as an html table.
#'
#' @param tib A tibble to display as a data frame.  
#'
#' @returns
#' @export
#'
#' @examples
create_dt <- function(tib) {
  # Copied from "https://martinctc.github.io/blog/vignette-downloadable-tables-in-rmarkdown-with-the-dt-package/".
  DT::datatable(tib,
                extensions = 'Buttons',
                options = list(dom = 'lfrtipB',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")))) %>% 
    DT::formatStyle(columns = colnames(asm_tib_sumstats), backgroundColor = "white", background = "white")
}


#' Title
#'
#' @param tib 
#'
#' @returns
#' @export
#'
#' @examples
get_asm_sumstats <- function(tib){
  sumstats <- tib %>% group_by(project, taxon) %>%
    summarize(num = n(), 
              mean_cov = stringr::str_c(round(mean(coverage),2), "x"), 
              mean_percent_n50 = round(mean(n50/total_length)*100, 2), 
              mean_percent_n90 = round(mean(n90/total_length)*100, 2), 
              mean_percent_large_ctgs = round(mean(percent_genome_in_fragments_1_000_000_bp),2), 
              mean_BUSCO = round(mean(percent_complete_busco),2))
  return(sumstats)
}



# CONTIGUITY_PLOT FUNCTIONS
#' Create attractive scatterplots with N50 or N90 as the Y-axis
#'
#' @param in_data 
#' @param dependent_metric 
#' @param contiguity_metric 
#' @param fill 
#'
#' @returns
#' @export
#'
#' @examples
plot_contiguity <- function(in_data = NULL, dependent_metric, contiguity_metric, fill) {
  
  # Select dependent var
  if (dependent_metric == "coverage") {
    dependent_var = in_data$coverage
  } else if (dependent_metric == "BUSCO") {
    dependent_var = in_data$percent_complete_busco
  } else if (dependent_metric == "longrds") {
    dependent_var=in_data$percent_longrds
  } else {
    stop("No dependent metric provided")
  }
  
  # Select contiguity metric (independent var)
  if (contiguity_metric == "N50") {
    contiguity_value = (in_data$n50 / in_data$total_length)
    point_shape = 21
  } else if (contiguity_metric == "N90") {
    contiguity_value = (in_data$n90 / in_data$total_length)
    point_shape = 24
  } else {
    stop("No contiguity metric provided")
  }
  
  if (fill == "project") {
    colors = factor(in_data$project); project_legend = "Project"
  } else if (fill == "highly_repetitive") {
    colors = factor(in_data$highly_repetitive); project_legend = "Repetitive genome?"
  } else if (fill == "taxon") {
    colors = factor(in_data$taxon); project_legend = "Taxon"
  } else {
    stop("No fill given")
  }
  
  # Format species labels
  species_labels = str_c(str_c(substr(in_data$genus, 1, 1), ".", sep=""), in_data$species, sep = " ")
  
  # Plot
  plot <- in_data %>% ggplot(aes(x = dependent_var, y = contiguity_value)) +
    geom_point(shape = point_shape, color = "black", cex = 2, aes(fill = colors)) +
    geom_label_repel(label = species_labels, fontface = "italic") +
    stat_poly_line() +
    stat_poly_eq(use_label("R2")) +
    xlab("") + ylab(contiguity_metric) +
    scale_y_continuous(labels = scales::label_percent()) +
    guides(fill=guide_legend(title = project_legend)) +
    theme_bw()
  
  if (dependent_metric == "coverage") {
    plot<-plot + scale_x_continuous(labels = scales::label_currency(prefix = "", suffix = "x"))
  } else if (dependent_metric =="longrds" | dependent_metric =="BUSCO") {
    plot<-plot + scale_x_continuous(labels = scales::label_currency(prefix = "", suffix = "%"))
    
  }
  return(plot)
}
compose_contiguity_plots <- function(plot1, plot2, title = "", x_title = "") {
  plot <- ggarrange(plot1, plot2, nrow=2, common.legend = TRUE, legend = "right")
  plot <- annotate_figure(plot, top = text_grob(title, color = "black", face = "bold", size = 14),
                          bottom = text_grob(x_title),
                          left = text_grob("Contiguity metric as percent of genome length", rot = 90))
  return(plot)
}
