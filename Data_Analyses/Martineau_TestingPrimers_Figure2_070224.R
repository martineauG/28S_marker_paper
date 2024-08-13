##################################################################
######### A. Load the working environment and packages ###########
##################################################################


## A.1 Set working directory
# setwd("your_directory")


# A.2 download necessary librairies
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
combined_df <- read_and_prepare_data("taxonomy_all_OTUs_070224.csv",
                                     "abundance_all_OTUs_070224.csv", 
                                     "Martineau_TestingPrimers_CO1AND28S_metadata_allOTUs.csv", "Porifera")

# Remove singletons from combined df
combined_df_nosingletons = as.data.frame(combined_df)
columns_keep <- colnames(combined_df_nosingletons)[colnames(combined_df_nosingletons) %in% c("AllOTUs", "ExpectedAllOTUs")]
columns_keep
for (i in columns_keep) {
  combined_df_nosingletons[, i] <- ifelse(combined_df_nosingletons[, i] > 1, combined_df_nosingletons[, i], 0)
}

ASV_table_DNA <- subset_numeric_columns(combined_df_nosingletons)
taxonomy_table_DNA <- subset_character_columns(combined_df_nosingletons)
flipped_ASV_with_metadata <- merge_flipped_ASV_with_metadata(t(ASV_table_DNA), metadata_df_DNA)

long_df <- elongate_df(combined_df_nosingletons, metadata_df_DNA)
long_df_variable <- extract_taxa_col(taxonomy_table_DNA, "species", long_df)

# Expected_species: make a vector with species that are in the expected cocktail
Expected_species <- long_df_variable %>%
  filter(Sample == "ExpectedAllOTUs") %>%
  filter(value > 0)


###############################################
######### C. Figure mock community 2 ##########
###############################################

# C.2 Select only the 3 ARMS units
# Filter data frame
long_taxa_ASV_df_filtered <- long_df_variable %>%
  filter(cocktail == "All_OTU(4)") %>%
  group_by(species) %>%
  mutate(sum_species = sum(value)) %>%
  filter(sum_species > 0) %>%
  # Only select the expected species
  filter(species %in% Expected_species$species)

# C.3 Add the other taxonomic groups to data frame: genus, class, order
other_taxonomic_levels <- taxonomy_table_DNA %>%
  filter(rownames(.) %in% long_taxa_ASV_df_filtered$variable) %>%
  select(-species) %>%
  mutate(variable = rownames(.))

# Merge elongated df with other taxonomic information
factor_boxplot <- left_join(long_taxa_ASV_df_filtered, other_taxonomic_levels, by = "variable")
factor_boxplot$ARMUNIT <- as.character(factor_boxplot$ARMUNIT)


# C.4 Find the number of detected OTUs per order
factor_boxplot_order <- factor_boxplot %>%
  filter(value > 1) %>%
  group_by(Sample, order) %>%
  dplyr::summarize(sum_OTUs_order = sum(n()))

# Make df wider
factor_boxplot_order_wider <- factor_boxplot_order %>%
  pivot_wider(names_from = Sample, values_from = sum_OTUs_order, values_fill = 0)

# Calculate % of observed species per order (observed/expected)
factor_boxplot_order_wider_perc <- factor_boxplot_order_wider %>%
  mutate(perc = round(AllOTUs / ExpectedAllOTUs * 100, 0),
         #  perc_corr = replace_inf(perc),
         not_detected = 100 - perc)


# C.8 Calculate the proportion of species that were detected vs. not-detected in the data set
proportional_barplot <- factor_boxplot_order_wider %>%
  mutate(not_detected = ExpectedAllOTUs - AllOTUs) %>%
  dplyr::rename(Detected = AllOTUs) %>%
  pivot_longer(cols = c("Detected", "not_detected"), names_to = "Status", values_to = "values") %>%
  arrange(ExpectedAllOTUs) 

# Add the percentage to the dataframe to code it on the plot
proportional_barplot_perc <- left_join(proportional_barplot, factor_boxplot_order_wider_perc[, c("order", "perc")], by = "order") %>%
  distinct()
proportional_barplot_perc_1 <- proportional_barplot_perc[!(duplicated(proportional_barplot_perc)),]

x.axis.order.2 = unique(proportional_barplot$order)
proportional_barplot_perc_1$Status <- factor(proportional_barplot_perc_1$Status, levels = c("not_detected", "Detected"))


# C.9 Barplot with bars proportional to values length
barplot_perc_prop <- ggplot(proportional_barplot_perc_1, aes(x = order, y = values)) +
  geom_bar (aes(fill = Status), stat= "identity",  width = 0.5 ,col = "black") +
  xlab("Porifera order") +
  ylab("Number of sponge species present in the mock community") +
  ylim(0, 25) +
  theme_classic()+
  scale_fill_manual (values = c("mediumpurple2","lightgoldenrod"), 
                     breaks = c("Detected", "not_detected"),
                     labels =c("Detected", "Not detected")) +
  scale_x_discrete (limits = x.axis.order.2)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  # Add V line for false positives
  geom_vline(aes(xintercept = order), 
             data = proportional_barplot %>% filter (values < 0),
             col = "red", alpha = 0.5, 
             size = 2, linetype= 2) +
  geom_text(aes(label = paste0(perc, "%"), y = ExpectedAllOTUs), 
            position = position_dodge(width = 0.5),
            hjust = -0.5,
            color = "red", size = 3) +
  coord_flip() +
  theme(text = element_text(family = "Open Sauce Sans")) 


barplot_perc_prop

## C.10 Calculate the number and the total % of species detected and not-detected
number_species_detected <- proportional_barplot_perc_1 %>%
  dplyr::filter(Status == "Detected") %>%
  dplyr::summarise(sum = sum(values))
number_species_detected

total_number_species <- sum(proportional_barplot_perc_1$ExpectedAllOTUs)/2
total_number_species


perc_species_detected <- (number_species_detected/total_number_species)*100
perc_species_detected
