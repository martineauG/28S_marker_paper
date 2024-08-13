##################################################################
######### A. Load the working environment and packages ###########
##################################################################


## A.1 Set working directory
# setwd("your_directory")


## A.2 download necessary librairies
necessarypackages = c("dplyr", "magrittr", "tidyr", "taxize", "purrr", "data.table", "knitr","tidyverse", "plyr", "ggplot2", "reshap2", "RColorBrewer", "vegan", "stringr", "ggpubr", "extrafont", "effects")
# Install packages
for (package in necessarypackages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}



##############################################################################
######### B. Load the data and general data wrangling: clean tables ##########
##############################################################################
# Function to read three files and prepare the data
read_and_prepare_data <- function(taxonomy_df, abundance_ASV_df, meta_data_df, phylum) {
  if (missing(phylum)) {
    taxonomy_DNA <- fread(taxonomy_df)
    ASV_table_DNA <- fread(abundance_ASV_df, header = TRUE)
    metadata_df_DNA <- as.data.frame(fread(meta_data_df))
    
    # Set row names
    rownames(ASV_table_DNA) <- ASV_table_DNA$ID_ASV
    rownames(taxonomy_DNA) <- taxonomy_DNA$ID_ASV
    rownames(metadata_df_DNA) <- metadata_df_DNA$Sample
    
    assign("metadata_df_DNA", metadata_df_DNA, envir = .GlobalEnv)
    
    # Replace NAs with 0 in ASV table
    ASV_table_DNA[is.na(ASV_table_DNA)] <- 0
    
    # Merge taxonomy and ASV tables
    combined_table <- inner_join(taxonomy_DNA, ASV_table_DNA, by = "ID_ASV")
    
    # Remove rows where taxa column is empty
    combined_df_working <- combined_table %>% 
      as.data.frame()
    
    rownames(combined_df_working) <- combined_df_working$ID_ASV
    combined_df_working$ID_ASV <- NULL
    
    return(combined_df_working)
  } else {
    taxonomy_DNA <- fread(taxonomy_df)
    ASV_table_DNA <- fread(abundance_ASV_df, header = TRUE)
    metadata_df_DNA <- as.data.frame(fread(meta_data_df))
    
    # Set row names
    rownames(ASV_table_DNA) <- ASV_table_DNA$ID_ASV
    rownames(taxonomy_DNA) <- taxonomy_DNA$ID_ASV
    rownames(metadata_df_DNA) <- metadata_df_DNA$Sample
    
    assign("metadata_df_DNA", metadata_df_DNA, envir = .GlobalEnv)
    
    # Replace NAs with 0 in ASV table
    ASV_table_DNA[is.na(ASV_table_DNA)] <- 0
    
    # Merge taxonomy and ASV tables
    combined_table <- inner_join(taxonomy_DNA, ASV_table_DNA, by = "ID_ASV")
    
    # Remove rows where taxa column is empty
    combined_df_working <- combined_table %>% filter(!is.na(phylum) & phylum != "")
    
    # Keep only Porifera results
    combined_df_working <- combined_df_working %>% filter(phylum == "Porifera") %>%
      as.data.frame()
    
    rownames(combined_df_working) <- combined_df_working$ID_ASV
    combined_df_working$ID_ASV <- NULL
    
    return(combined_df_working)
    
  }
}

# Function to subset numeric columns
subset_numeric_columns <- function(df) {
  return(df[sapply(df, is.numeric)])
}

# Function to subset character columns
subset_character_columns <- function(df) {
  return(df[sapply(df, is.character)])
}

# Function to merge flipped ASV table with metadata
merge_flipped_ASV_with_metadata <- function(flipped_ASV_df, metadata_df) {
  merged_df <- merge(flipped_ASV_df, metadata_df, by = "row.names", all.y = TRUE)
  rownames(merged_df) <- merged_df[, 1]
  merged_df$Row.names <- NULL
  return(merged_df)
}

# Function to make the data frame longer 
elongate_df <- function(combined_df, meta_data) {
  combined_full_df <- as.data.frame(t((combined_df)))
  combined_df_Meta_data <- merge (combined_full_df,  meta_data, by = "row.names",all.y = T)
  rownames(combined_df_Meta_data) = combined_df_Meta_data$Sample
  combined_df_Meta_data$Row.names <- NULL
  long_combined_ASV_Meta_data_df <- as.data.frame(reshape2::melt(combined_df_Meta_data, id.vars = names(meta_data)))
  long_combined_ASV_Meta_data_df$value = as.numeric( long_combined_ASV_Meta_data_df$value)
  #return(long_combined_ASV_Meta_data_df)
  return(long_combined_ASV_Meta_data_df)
  # assign("long_df", long_combined_ASV_Meta_data_df , envir = .GlobalEnv)
}

# Select only the species column and set ID_ASV as variable and add to the long table
extract_taxa_col <- function(working_taxonomy_table, taxa_col, long_df) {
  # Extract the specified column from working_taxonomy_table
  taxa_data <-working_taxonomy_table %>% 
    select ({{taxa_col}}) %>%
    as.data.frame()
  
  # Convert to data frame and merge with long_df
  merged_data <- as.data.frame(merge(long_df, data.frame(variable = rownames(taxa_data), value = taxa_data), by = "variable"))
  #assign("long_df_variable", merged_data , envir = .GlobalEnv)
  return(merged_data)
}

# Usage example
combined_df <- read_and_prepare_data("sponge_ntARMS17_abundance_table.csv",
                                     "sponge_ntARMS17_taxonomy_table.csv", 
                                     "Martineau_TestingPrimers_metadata_17.csv")

ASV_table_DNA <- subset_numeric_columns(combined_df)
taxonomy_table_DNA <- subset_character_columns(combined_df)
flipped_ASV_with_metadata <- merge_flipped_ASV_with_metadata(t(ASV_table_DNA), metadata_df_DNA)

long_df <- elongate_df(combined_df, metadata_df_DNA)
long_df_variable <- extract_taxa_col(taxonomy_table_DNA, "species", long_df)



#########################################################
######### C. Taxonomic barplots to combine both #########
#########################################################

## C.1 Select the 20 most abundant families/phylum (taxa) and label and combine the other ones as other labelled the "others" together
# Define the first 19 species
top_species <- long_df_variable %>%
  dplyr::group_by(species) %>%
  drop_na() %>%
  dplyr::summarize(sum_asv = sum(value)) %>%
  arrange(-sum_asv) %>%
  top_n(19) %>%
  pull(species)


# C.2 Replace the other species by "others"
Taxa_ASV_abundant <- long_df_variable %>%
  drop_na() %>%
  mutate(species = if_else(species %in% top_species, species, "Other")) %>%
  mutate(value = as.numeric(value))

## C.5 Summarise the number of species per sample
ASV_per_group_per_taxa <- Taxa_ASV_abundant %>%
  group_by(Sample, species) %>%
  mutate(mean_abundance_ASV_group = mean(value, na.rm = TRUE)) %>%
  mutate(species = str_replace_all(species, "JV-2020", ""))

## C.6 Make and extend color palette for graph
make_color_palette <- function(number_palette, number_taxa) {
  color_barplot <- brewer.pal(number_palette, "Paired")
  color_barplot <- colorRampPalette(color_barplot)(number_taxa)
  return(color_barplot)
}

color_barplot <- make_color_palette(12, 20)

color_barplot_2 <- c("#A6CEE3", "#7EBA98", "#F06C45", "#D9A295", "#F0EB99", "#F32B65", "#7D54A5", "#FE982C",
                     "#98D277", "#5B9EC9", "#6F9E4C", "#F16667", "#E42022", "#FDBB69", "#B15928", "#B294C7",
                     "#2D82AF", "#9E8099", "#52AF43", "#DBB466", "#B15928")

## C.7 Plot everything
taxonomic_barplot <- ggplot(ASV_per_group_per_taxa, aes(fill = species, y = value, x = Sample)) +
  geom_bar(position = "fill", stat = "identity", color = "black") + 
  theme_classic() + 
  ylab("Species composition (read relative abundance)") +
  xlab("Replicate") +
  facet_grid(. ~ cocktail, space = "free_x", scales = "free_x") +
  scale_x_discrete(breaks = c("ARMS17R1", "ARMS17R2", "ARMS17R3", "Equimolar17R1", "Equimolar17R2", "Equimolar17R3", "%cover17R1", "%cover17R2", "%cover17R3"),
                   labels = c("1", "2", "3", "1", "2", "3", "1", "2", "3")) +
  scale_fill_manual(values = color_barplot_2, breaks = str_replace_all(c(top_species, "Other"), "JV-2020", "")) +
  labs(fill = "Species") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank())
taxonomic_barplot
