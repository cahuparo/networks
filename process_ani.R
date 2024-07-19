#!/usr/bin/env Rscript
#' Process skani clustering matrix to extract network and detailled tables
#'
#' Description: This script processes a skani clustering matrix to extract network graphs and detailed tables. 
#' It reads the ANI (Average Nucleotide Identity) matrix, processes it to create a network of high similarity clusters, 
#' and generates visualizations and tables summarizing the clustering results. The script identifies species groups 
#' based on a given threshold, visualizes the network of relationships, and extracts detailed information about the groups.
#'
#' Author: Camilo H. Parada Rojas
#' Contact: paradar@oregonstate.edu
#' Institution: Oregon State University - BPP
#' Date: 2024-07-19
#'
#' License: This script is licensed under the MIT License.
#' 
#' Citation:
#' Camilo H. Parada Rojas (2024). ANI network analysis. Version 1.0. https://github.com/cahuparo/networks/edit/main/process_ani.R
#'
#' Acknowledgments:
#' This work was supported by [Brady grant]. Special thanks to [Riley B., Jeff C.].
#'

# Load necessary libraries
libraries <- c("tidyr", "dplyr", "igraph", "ComplexHeatmap", "circlize", 
               "tidygraph", "ggraph", "ggplot2", "RColorBrewer", "rentrez")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

lapply(libraries, install_if_missing)

# Function to read and process the ANI matrix
read_ani_matrix <- function(file) {
  data <- readLines(file)
  items <- as.integer(data[1])
  labels <- c()
  matrix <- matrix(0, nrow=items, ncol=items)
  
  for (i in 2:length(data)) {
    line <- unlist(strsplit(data[i], "\t"))
    labels <- c(labels, gsub(".*/|\\.fna$", "", line[1]))
    if (length(line) > 2) {
      matrix[i-1, 1:(i-2)] <- as.numeric(line[2:length(line)])
    }
  }
  
  rownames(matrix) <- labels
  colnames(matrix) <- labels
  matrix <- t(matrix)
  matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(matrix)]
  
  return(matrix)
}

# Function to process ANI data for a given threshold
process_ani_data <- function(file, threshold) {
  ani_matrix <- read_ani_matrix(file)
  diag(ani_matrix) <- 100
  
  # Convert the matrix to a data frame
  ani_df <- as.data.frame(ani_matrix)
  ani_df$strain1 <- rownames(ani_df)
  
  # Convert to long format
  ani_long <- gather(ani_df, strain2, ani, -strain1)
  
  # Filter for high similarity using the threshold
  high_similarity <- ani_long %>%
    filter(ani >= threshold) %>%
    select(strain1, strain2)
  
  # Create an undirected graph
  species_graph <- graph_from_data_frame(high_similarity, directed=FALSE)
  
  # Identify connected components (species groups)
  components <- components(species_graph)
  
  # Assign each strain to a species group based on connected components
  ani_df$species_group <- components$membership[match(ani_df$strain1, names(components$membership))]
  
  # Calculate the number of nodes in each species group
  group_counts <- data.frame(
    species_group = components$membership,
    count = table(components$membership)[components$membership]
  )
  group_counts <- unique(group_counts)
  
  # Filter out groups with 5 or fewer nodes
  large_groups <- group_counts %>%
    filter(count.Freq > 5)
  head(large_groups)
  # Get the list of strains that belong to large groups
  large_group_strains <- ani_df %>%
    filter(species_group %in% large_groups$species_group)
  
  # Filter the graph to include only nodes in large groups
  large_group_graph <- induced_subgraph(species_graph, V(species_graph)$name %in% large_group_strains$strain1)
  
  # Recalculate components for the filtered graph
  components_large <- components(large_group_graph)
  
  # Create a tidygraph object with species group membership
  tg_large <- as_tbl_graph(large_group_graph) %>%
    activate(nodes) %>%
    mutate(species_group = factor(components_large$membership))
  
  # Ensure node names are properly assigned
  V(tg_large)$name <- V(large_group_graph)$name
  
  # Mark nodes as GCF or not
  tg_large <- tg_large %>%
    activate(nodes) %>%
    mutate(is_gcf = ifelse(grepl("^GCF", name), TRUE, FALSE))
  
  # Create a color palette for non-GCF genomes using viridis
  species_groups <- unique(V(tg_large)$species_group)
  
  n <- length(species_groups)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col = sample(col_vector, n)
  
  palette <- setNames(col, species_groups)
  
  # Ensure species_group is a factor
  tg_large <- tg_large %>%
    activate(nodes) %>%
    mutate(species_group = factor(species_group),
           is_gcf = ifelse(is_gcf, "GCF", as.character(species_group)))
  
  # Update color palette to include GCF
  palette <- c(palette, "GCF" = "#262B4C")
  
  # Error handling for network visualization
  tryCatch({
    network_plot <- ggraph(tg_large, layout = "fr") + 
      geom_edge_link(aes(alpha = 0.5)) +
      geom_node_point(aes(color = is_gcf), size = 5) +
      scale_color_manual(values = palette) +
      theme_void() +
      labs(color = "Species Group") +
      ggtitle(paste("Network Graph of Species Groups Larger than 2 Nodes at", threshold, "ANI threshold. Dark blue circles correspond to GCF Genomes"))
    
    # Save network plot as PDF
    pdf(paste0("network_at_", threshold, "_ANI.pdf"), width = 16, height = 9)
    print(network_plot)
    dev.off()
  }, error = function(e) {
    message("An error occurred during network visualization: ", e$message)
  })
  
  # Ensure is_gcf is correctly populated
  tg_large <- tg_large %>%
    activate(nodes) %>%
    mutate(is_gcf = grepl("^GCF", name))
  
  # Identify clusters with the largest number of non-GCF strains
  non_gcf_strains <- tg_large %>%
    activate(nodes) %>%
    filter(!is_gcf) %>%
    as_tibble() %>%
    group_by(species_group) %>%
    tally(name = "count_non_gcf")
  
  # Get the largest clusters based on the count of non-GCF strains
  largest_clusters <- non_gcf_strains %>%
    filter(count_non_gcf >= 5)
  
  # Extract GCF genomes from these large clusters
  gcf_genomes_in_large_clusters <- V(tg_large)$name[V(tg_large)$species_group %in% largest_clusters$species_group & V(tg_large)$is_gcf]
  
  if (length(gcf_genomes_in_large_clusters) > 0) {
    # Function to retrieve taxonomic ID and scientific name using rentrez
    get_taxid_and_name <- function(accession) {
      accession <- sub("_genomic", "", accession) # Remove suffix
      accession <- sub("_", ".", accession) # Replace underscore with dot before version number
      search_result <- entrez_search(db = "assembly", term = accession)
      if (length(search_result$ids) > 0) {
        summary_result <- entrez_summary(db = "assembly", id = search_result$ids[1])
        taxid <- summary_result$taxid
        scientific_name <- summary_result$speciesname
        return(c(taxid, scientific_name))
      } else {
        return(c(NA, NA))
      }
    }
    
    # Retrieve taxonomic IDs and scientific names for GCF genomes in the largest clusters
    taxids_and_names <- t(sapply(gcf_genomes_in_large_clusters, get_taxid_and_name))
    
    # Create a data frame with GCF genomes and their taxonomic IDs and scientific names
    gcf_taxids_names <- data.frame(genome = gcf_genomes_in_large_clusters, taxid = taxids_and_names[, 1], scientific_name = taxids_and_names[, 2], stringsAsFactors = FALSE)
  } else {
    gcf_taxids_names <- data.frame(genome = character(), taxid = character(), scientific_name = character(), stringsAsFactors = FALSE)
  }
  
  # Combine all the information into a single data frame
  final_table <- largest_clusters %>%
    mutate(
      gcf_genomes = sapply(species_group, function(group) {
        paste(V(tg_large)$name[V(tg_large)$species_group == group & V(tg_large)$is_gcf], collapse = ", ")
      }),
      gcf_taxids = sapply(species_group, function(group) {
        paste(gcf_taxids_names$taxid[gcf_taxids_names$genome %in% V(tg_large)$name[V(tg_large)$species_group == group & V(tg_large)$is_gcf]], collapse = ", ")
      }),
      gcf_scientific_names = sapply(species_group, function(group) {
        paste(gcf_taxids_names$scientific_name[gcf_taxids_names$genome %in% V(tg_large)$name[V(tg_large)$species_group == group & V(tg_large)$is_gcf]], collapse = ", ")
      }),
      non_gcf_strains = sapply(species_group, function(group) {
        paste(V(tg_large)$name[V(tg_large)$species_group == group & !V(tg_large)$is_gcf], collapse = ", ")
      }),
      count_gcf = sapply(species_group, function(group) {
        sum(V(tg_large)$is_gcf[V(tg_large)$species_group == group])
      })
    ) %>%
    select(
      species_group,
      count_non_gcf,
      count_gcf,
      gcf_genomes,
      gcf_taxids,
      gcf_scientific_names,
      non_gcf_strains
    )
  
  # Ensure non-GCF strains are included if no GCF strains are found
  if (nrow(final_table) == 0) {
    final_table <- non_gcf_strains %>%
      mutate(
        gcf_genomes = NA,
        gcf_taxids = NA,
        gcf_scientific_names = NA,
        non_gcf_strains = sapply(species_group, function(group) {
          paste(V(tg_large)$name[V(tg_large)$species_group == group & !V(tg_large)$is_gcf], collapse = ", ")
        }),
        count_gcf = 0
      ) %>%
      select(
        species_group,
        count_non_gcf,
        count_gcf,
        gcf_genomes,
        gcf_taxids,
        gcf_scientific_names,
        non_gcf_strains
      )
  }
  
  # Print the final table
  print(final_table)
  
  # Write the table to a tab-delimited file
  write.table(final_table, paste0(threshold, "_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save the list of non-GCF genomes to a file
  non_gcf_genomes <- unique(unlist(strsplit(final_table$non_gcf_strains, ", ")))
  writeLines(non_gcf_genomes, paste0(threshold, "_non_gcf_genomes.txt"))
  
  # Save the list of GCF genomes to a file
  gcf_genomes <- unique(unlist(strsplit(final_table$gcf_genomes, ", ")))
  writeLines(gcf_genomes, paste0(threshold, "_gcf_genomes.txt"))
  
  # Combine genomes from all species groups in final_table
  all_genomes <- unique(c(non_gcf_genomes, gcf_genomes))
  
  # Subset the ANI matrix for all genomes in the final table
  ani_matrix_subset <- ani_matrix[all_genomes, all_genomes]
  
  # Create a data frame for row and column annotations
  annotation_df <- data.frame(
    genome = all_genomes,
    species_group = factor(tg_large %>% 
                             activate(nodes) %>% 
                             as_tibble() %>% 
                             filter(name %in% all_genomes) %>% 
                             pull(species_group))
  )
  
  # Create column and row annotations for the heatmap based on species group
  row_annotation <- rowAnnotation(
    species_group = annotation_df$species_group,
    col = list(species_group = palette)
  )
  
  col_annotation <- HeatmapAnnotation(
    species_group = annotation_df$species_group,
    col = list(species_group = palette)
  )
  
  # Improved color palette for heatmap
  color_palette <- colorRamp2(c(0, 50, 100), c("#313695", "#ffffbf", "#a50026"))
  
  # Visualize the ANI matrix as a heatmap
  heatmap <- ComplexHeatmap::Heatmap(
    ani_matrix_subset, 
    cluster_rows = TRUE, 
    cluster_columns = TRUE, 
    show_row_names = TRUE, 
    show_column_names = TRUE, 
    col = color_palette,
    heatmap_legend_param = list(title = "ANI (%)"),
    row_split = annotation_df$species_group,
    column_split = annotation_df$species_group,
    left_annotation = row_annotation,
    top_annotation = col_annotation
  )
  
  # Save heatmap as PDF
  pdf(paste0("clustered_at_", threshold, "_ANI_heatmap.pdf"), width = 16, height = 9)
  draw(heatmap, annotation_legend_side = "bottom", heatmap_legend_side = "right")
  dev.off()
  
  # Create individual files for each species group
  for (group in unique(tg_large %>% activate(nodes) %>% pull(species_group))) {
    non_gcf_genomes <- V(tg_large)$name[V(tg_large)$species_group == group & !V(tg_large)$is_gcf]
    gcf_genomes <- V(tg_large)$name[V(tg_large)$species_group == group & V(tg_large)$is_gcf]
    
    combined_genomes <- unique(c(non_gcf_genomes, gcf_genomes))
    
    file_name <- paste0("species_group_", group, "_genomes_at_", threshold, "ani.txt")
    writeLines(combined_genomes, file_name)
  }
}

# List of threshold values
thresholds <- c(90, 92, 93, 94, 95, 96, 97, 97.5, 98, 99, 100)

# Read input file from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide the ANI matrix file as a command-line argument.")
}
input_file <- args[1]

# Iterate over the threshold values
for (threshold in thresholds) {
  process_ani_data(input_file, threshold)
}

