##################################################################
######### A. Load the working environment and packages ###########
##################################################################


## A.1 Set working directory
setwd("your_directory")


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
combined_df <- read_and_prepare_data("combined28SCOI_taxonomy_table_071224.csv",
                                     "combined28SCOI_abundance_table_071224.csv", 
                                     "Martineau_TestingPrimers_CO1AND28S_metadata.csv", "Porifera")
ASV_table_DNA <- subset_numeric_columns(combined_df)
taxonomy_table_DNA <- subset_character_columns(combined_df)
flipped_ASV_with_metadata <- merge_flipped_ASV_with_metadata(t(ASV_table_DNA), metadata_df_DNA)

long_df <- elongate_df(combined_df, metadata_df_DNA)
long_df_variable <- extract_taxa_col(taxonomy_table_DNA, "species", long_df)



#################################################################################################
######### C. Figure 3: Primer biases lollipop chart and corresponding chi-square tests #########
#################################################################################################

## C.1  Function to index species supposed to be there
index_species <- function(species_list, combined_table, unit, cocktail) {
  combined_table$ID_ASV <- row.names(combined_table)
  observed_species <- combined_table[which(combined_table$species %in% species_list), "ID_ASV"]
  curated_data <- combined_table %>%
    filter(ID_ASV %in% observed_species) %>%
    select(c(ID_ASV, matches(paste0(cocktail, unit))))
  return(na.omit(curated_data))
}

## C.2 Make vectors for each expected species and extract ASV table for each unit
Species_17 = c( "Clathrinidae sp. 4 JV-2020",
                "Leucosolenida sp. 1 JV-2020",
                "Leucosolenida sp. 3 JV-2020",
                "Leucosolenida sp. 7 JV-2020",
                "Clathrinida sp. 1 JV-2020",
                "Haliclona sp. 1 JV-2020",
                "Haliclona sp. 2 JV-2020",
                "Haplosclerida sp. 10 JV-2020",
                "Ancorinidae sp. 3 JV-2020",
                "Tethya sp. 5 JV-2020",
                "Suberitidae sp. 1 JV-2020",
                "Hymeniacidon sp. 2 JV-2020",
                "Oscarella sp. 6 JV-2020")

Species_18 = c("Clathrinidae sp. 3 JV-2020",
               "Clathrinidae sp. 4 JV-2020",
               "Leucosolenida sp. 1 JV-2020",
               "Leucosolenida sp. 5 JV-2020",
               "Leucosolenida sp. 7 JV-2020",
               "Clathrinida sp. 1 JV-2020",
               "Haliclona sp. 1 JV-2020",
               "Haliclona sp. 2 JV-2020",
               "Ancorinidae sp. 3 JV-2020",
               "Dendroceratida sp. 3 JV-2020",
               "Tethya sp. 3 JV-2020",
               "Tethya sp. 5 JV-2020",
               "Poecilosclerida sp. 2 JV-2020",
               "Suberitidae sp. 1 JV-2020",
               "Halichondria sp. 1 JV-2020",
               "Ircinia sp. 2 JV-2020")

Species_20 = c("Leucosolenida sp. 11 JV-2020",
               "Leucosolenida sp. 13 JV-2020",
               "Clathrinidae sp. 3 JV-2020",
               "Clathrinidae sp. 5 JV-2020",
               "Clathrinidae sp. 4 JV-2020",
               "Leucosolenida sp. 1 JV-2020",
               "Leucosolenida sp. 5 JV-2020",
               "Leucosolenida sp. 7 JV-2020",
               "Haliclona sp. 1 JV-2020",
               "Ancorinidae sp. 3 JV-2020",
               "Aplysilla rosea",
               "Tethya sp. 3 JV-2020",
               "Dysidea sp. 5 JV-2020",
               "Poecilosclerida sp. 1 JV-2020",
               "Poecilosclerida sp. 7 JV-2020",
               "Suberitidae sp. 1 JV-2020",
               "Suberitida sp. 3 JV-2020",
               "Plakina sp. 1 JV-2020")

# Index species for each unit
OTU_17 <- index_species(Species_17, combined_df, 17, "Equimolar")
OTU_18 <- index_species(Species_18, combined_df, 18, "Equimolar")
OTU_20 <- index_species(Species_20, combined_df, 20, "Equimolar")


## C.3 Make an abundance table for all of the species and units
# Clean up the table
ASV_table_DNA$ID_ASV <- NULL
ASV_table_DNA[is.na(ASV_table_DNA)] <- 0
ASV_table_DNA_curated <- ASV_table_DNA %>%
  select(matches("Equimolar"))

# Merge all units together
ASV_table_curated_raw <- join_all(list(OTU_17, OTU_18, OTU_20),
                                  by = "ID_ASV", type = "full") %>%
  left_join(rownames_to_column(taxonomy_table_DNA, var = "ID_ASV"), by = "ID_ASV") %>%
  select(matches("Equimolar"), order, species) %>%
  dplyr::rename (UnitA = Equimolar17,
                 UnitB = Equimolar18,
                 UnitC = Equimolar20)


## C.4 Compute the average # of reads per species and the expected number of reads in each unit
# Compute the average per species and make df longer for plotting
lolipop_merged_df_long_raw <- ASV_table_curated_raw %>%
  mutate(average = rowMeans(select(., matches("Unit")), na.rm = TRUE)) %>%
  pivot_longer(cols = matches("Unit"), names_to = "Unit", values_to = "value")

# Rename the species without JV
lolipop_merged_df_long_raw$species_short_name <- str_replace_all(lolipop_merged_df_long_raw$species, "JV-2020", "")

# Determine the average number of reads
sum_of_reads <- colSums(ASV_table_DNA_curated, na.rm = TRUE)
expected_read_average <- sum_of_reads / c(18, 16, 13)
average_expected_average <- mean(expected_read_average)


## C.5 Calculate deviation fold from expected number of reads
# Deviation function 
deviation_calculation <- function(df, expected_read_average, unit_name) {
  deviation_col_name <- paste0("Equimolar", unit_name)
  deviation <- (df[[deviation_col_name]] - expected_read_average) / expected_read_average
  df[[deviation_col_name]] <- deviation
  return(df)
}

OTU_17_deviation <- deviation_calculation(OTU_17, expected_read_average["Equimolar17"], 17)
OTU_18_deviation <- deviation_calculation(OTU_18, expected_read_average["Equimolar18"], 18)
OTU_20_deviation <- deviation_calculation(OTU_20, expected_read_average["Equimolar20"], 20)

# Merge and clean up the deviation tables
ASV_table_deviation <- join_all(list(OTU_17_deviation, OTU_18_deviation, OTU_20_deviation), by = "ID_ASV", type = "full") %>%
  left_join(rownames_to_column(taxonomy_table_DNA, var = "ID_ASV"), by = "ID_ASV") %>%
  select(matches("Equimolar"), order, species) %>%
  dplyr::rename (UnitA = Equimolar17,
                 UnitB = Equimolar18,
                 UnitC = Equimolar20)
write.csv(ASV_table_deviation, "./tables/Martineau_TestingPrimers_folddeviation_071224.csv", row.names = FALSE)

# Average units together
lolipop_merged_df_deviation  <- ASV_table_deviation %>%
  mutate(average_deviation = rowMeans(select(., matches("Unit")), na.rm = TRUE)) %>%
  select(starts_with("Unit"), species, order, average_deviation)

# Make df longer for plotting
lolipop_merged_df_long_deviation  <- lolipop_merged_df_deviation %>%
  pivot_longer(cols = starts_with("Unit"), names_to = "Unit", values_to = "Deviation") %>%
  mutate(species_short_name = str_replace_all(species, "JV-2020", "")) %>%
  na.omit()

# Make order for deviation
lolipop_order <- c("Tethya sp. 5 ",
                   "Tethya sp. 3 " ,
                   "Hymeniacidon sp. 2 ",
                   "Suberitidae sp. 1 ",
                   "Suberitida sp. 3 ",
                   "Halichondria sp. 1 ",
                   "Poecilosclerida sp. 1 ",
                   "Poecilosclerida sp. 7 ",
                   "Poecilosclerida sp. 2 ",
                   "Ancorinidae sp. 3 ",
                   "Haplosclerida sp. 10 ",
                   "Haliclona sp. 2 ",
                   "Haliclona sp. 1 ",
                   "Dendroceratida sp. 3 ",
                   "Aplysilla rosea",
                   "Dysidea sp. 5 ",
                   "Ircinia sp. 2 ",
                   "Plakina sp. 1 ",
                   "Oscarella sp. 6 ",
                   "Clathrinidae sp. 4 ",
                   "Clathrinidae sp. 3 ",
                   "Clathrinida sp. 1 ", 
                   "Clathrinidae sp. 5 ",   
                   "Leucosolenida sp. 1 ",
                   "Leucosolenida sp. 11 ",
                   "Leucosolenida sp. 3 ",
                   "Leucosolenida sp. 5 ",
                   "Leucosolenida sp. 13 ",
                   "Leucosolenida sp. 7 ")


## C.6 Plot everything
lolipop_plot_deviation_white <- ggplot(lolipop_merged_df_long_deviation, aes(x = species_short_name, y = average_deviation, label = average_deviation)) + 
  geom_point(stat = 'identity', aes(fill = order), shape = 21, size = 5, colour = "black") +
  geom_segment(aes(y = 0, x = species_short_name, yend = average_deviation, xend = species_short_name), color = "black") +
  scale_fill_manual(values = c("darkcyan", "cadetblue1", "darkseagreen3", "palegoldenrod", "antiquewhite3", "mistyrose", "darksalmon", "coral", "sienna3", "firebrick3"), 
                    breaks = c("Clathrinida", "Haplosclerida", "Homosclerophorida", "Leucosolenida", "Poecilosclerida", "Suberitida", "Tethyida" ,"Tetractinellida", "Dendroceratida", "Dictyoceratida")) +
  scale_x_discrete(limits = rev(lolipop_order)) +
  geom_hline(yintercept = 0, col = "red", size = 1, alpha = 0.7) +
  coord_flip() +
  theme_bw() +
  # ylim(-100,100)+
  theme(axis.text.y = element_text(hjust = 0), text = element_text(family = "Open Sauce Sans")) +
  xlab("Porifera species") +
  ylab("Average deviation from expected number of reads") +
  guides(shape = guide_legend(override.aes = list(colour = "white")))

lolipop_plot_deviation_white


# Plot 
## C.5 Chi-Square tests 
# Make data frame for the Chi-Square test
chi_square_df <- ASV_table_curated_raw %>%
  pivot_longer(cols = starts_with("Unit"),
               values_to = "observed",
               names_to = "expected",
               values_drop_na = TRUE) %>%
  mutate(Unit = expected) %>%
  mutate(expected = case_when(
    expected == "UnitA" ~ expected_read_average["Equimolar17"],
    expected == "UnitB" ~ expected_read_average["Equimolar18"],
    expected == "UnitC" ~ expected_read_average["Equimolar20"]
  )) %>%
  mutate(expected = as.numeric(expected)) %>%
  mutate(Chisquare = ((observed - expected)^2) / expected)

# Manually compute Chi-Square statistics for each ARMS units
compute_chisquare_stats <- function(data, unit_name) {
  Chisquare_df <- data %>%
    filter(Unit == unit_name) %>%
    select(-Unit)
  
  Chisquare_stats <- sum(Chisquare_df$Chisquare)
  df <- nrow(Chisquare_df) - 1
  p_value <- 1 - pchisq(q = Chisquare_stats, df = df)
  
  return(list(Chisquare_stats = Chisquare_stats, df = df, p_value = p_value))
}

# Apply the function to each ARMS unit
unit_names <- unique(chi_square_df$Unit)
chisquare_stats_list <- map(unit_names, ~ compute_chisquare_stats(chi_square_df, .x))

# Export the Chi-Square calculated values for tables S4-S6
write.csv(chi_square_df, "./tables/Martineau_TestingPrimers_MFP_ChiSquarevalues_020724.csv", row.names = FALSE)



################################################################################
######### D. Figure 4: Scatter plot %cover and read abundance and GLM ##########
################################################################################

## D.1 Similar to C.1: make a data frame with observed and expected for % coverage
# Adjust species 20, because Suberitidae sp. 3 is removed from that one 
Species_20.2 = c("Leucosolenida sp. 11 JV-2020",
                 "Leucosolenida sp. 13 JV-2020",
                 "Clathrinidae sp. 3 JV-2020",
                 "Clathrinidae sp. 5 JV-2020",
                 "Clathrinidae sp. 4 JV-2020",
                 "Leucosolenida sp. 1 JV-2020",
                 "Leucosolenida sp. 5 JV-2020",
                 "Leucosolenida sp. 7 JV-2020",
                 "Haliclona sp. 1 JV-2020",
                 "Ancorinidae sp. 3 JV-2020",
                 "Aplysilla rosea",
                 "Tethya sp. 3 JV-2020",
                 "Dysidea sp. 5 JV-2020",
                 "Poecilosclerida sp. 1 JV-2020",
                 "Poecilosclerida sp. 7 JV-2020",
                 "Suberitidae sp. 1 JV-2020",
                 "Plakina sp. 1 JV-2020")

# Index species for each unit: make 3 small abundance tables
pc_17 <- index_species(Species_17, combined_df, 17, "%cover")
pc_18 <- index_species(Species_18, combined_df, 18, "%cover")
pc_20 <- index_species(Species_20.2, combined_df,20,"%cover")
ARMS17 <- index_species(Species_17,combined_df, 17, "ARMS")[,1:2]
ARMS18 <- index_species(Species_18,combined_df, 18, "ARMS")[,1:2]
ARMS20 <- index_species(Species_20.2,combined_df, 20, "ARMS")[,1:2]

# Merge all units together
observed_observations <- function (cocktail, data_frame_vector){
  
  observed_ASV_table_pc <-  purrr::reduce(data_frame_vector, full_join, by = "ID_ASV") %>%
    dplyr::rename(UnitA = matches("17"), 
                  UnitB = matches("18"),
                  UnitC = matches("20")) %>%
    pivot_longer(cols = starts_with("Unit"), names_to = "Unit", values_to = paste("observed",cocktail, sep = "_")) %>%
    na.omit()
}

# Make a data frame for PC and for the ARMS (for GLM)
observed_ASV_table_pc <- observed_observations(cocktail = "%cover", data_frame_vector = list(pc_17, pc_18, pc_20))
observed_ASV_table_ARMS <- observed_observations(cocktail = "ARMS", data_frame_vector = list(ARMS17, ARMS18, ARMS20))

# Make Expected table
pc_17_expected <- index_species(Species_17, combined_df, 17, "ExpectedPC")
pc_18_expected <- index_species(Species_18, combined_df, 18, "ExpectedPC")
pc_20_expected <- index_species(Species_20.2, combined_df, 20,"ExpectedPC")

# Merge the expected tables
expected_ASV_table_pc <- join_all(list(pc_17_expected, pc_18_expected, pc_20_expected),
                                  by = "ID_ASV", type = "full") %>%
  dplyr::rename(UnitA = ExpectedPC17, 
                UnitB = ExpectedPC18,
                UnitC = ExpectedPC20) %>%
  pivot_longer(cols = starts_with("Unit"), names_to = "Unit", values_to = "expected") %>%
  na.omit()

# Join observed and expected and add species
merge_observed_expected <- function (observed_df, cocktail) {
  observed_col <- paste0("observed_", cocktail)
  
  pc_observed_expected <- observed_df %>%
    full_join(expected_ASV_table_pc, by = c("ID_ASV", "Unit")) %>%
    left_join(rownames_to_column(taxonomy_table_DNA, var = "ID_ASV"), by = "ID_ASV") %>%
    select(colnames(observed_df), expected, order, species) %>%
    
    # Add log transformation
    dplyr::mutate(logExpected =  log10 (expected+1),
                  logObserved =  log10(.data[[observed_col]] + 1))
  return(pc_observed_expected)
}
pc_observed_expected <- merge_observed_expected(observed_ASV_table_pc, "%cover")
ARMS_observed_expected <- merge_observed_expected(observed_ASV_table_ARMS, "ARMS")


## D.2 Plot everything
# Define color and shape values
color_values <- c("darkcyan", "cadetblue1", "darkseagreen3", "lightgoldenrod2", "antiquewhite3",
                  "rosybrown3", "darksalmon", "sienna3", "firebrick3", "tomato4", "coral4", "white")
shape_values <- c(21, 22, 24)

# Create the scatter plot
corr_plot <- ggscatter(pc_observed_expected, x = "logExpected", y = "logObserved", fill = "order", shape = "Unit",
                       size = 4, xlab = "log(Percent cover)", ylab = "log(Read Abundance from percent cover-based DNA concentration)",
                       font.family = "Open Sauce Sans", legend = "right", cor.coef.size = 15) +
  # Customize the color and shape scales
  scale_fill_manual(values = color_values, breaks = c("Clathrinida", "Haplosclerida", "Homosclerophorida",
                                                      "Leucosolenida", "Poecilosclerida", "Suberitida",
                                                      "Tethyida", "Tetractinellida", "Dendroceratida",
                                                      "Dictyoceratida", "Axinellida", "unclassified")) +
  scale_shape_manual(values = shape_values) +
  # Add linear regression line
  geom_smooth(method = "lm", color = "black") +
  # Customize theme
  theme(text = element_text(family = "Open Sauce Sans"), legend.position = "top") +
  # Remove fill legend
  guides(fill = "none") +
  # Add correlation coefficient label
  stat_cor(label.x = 0, label.y = 4)

corr_plot


## D.3 Stats: GLM 
# for PC and Expected
poisson.model.pc <- glm(`observed_%cover` ~ expected  , family= quasipoisson, data = pc_observed_expected)
summary(poisson.model.pc)
plot(allEffects(poisson.model.pc))
confint(poisson.model.pc, level = 0.95)
plot(poisson.model.pc)
plot(poisson.model.pc$y ~ fitted(poisson.model.pc))

# for ARMS and Expected
poisson.model.ARMS <- glm(observed_ARMS ~ expected  , family= quasipoisson, data = ARMS_observed_expected)
summary(poisson.model.ARMS)
plot(allEffects(poisson.model.ARMS))
confint(poisson.model.ARMS, level = 0.95)
plot(poisson.model.pc)
plot(poisson.model.pc$y ~ fitted(poisson.model.pc))



##########################################################################
######### E. Figure 5A-B: Detection ARMS heat map and pie charts #########
##########################################################################

## E.1 Data wrangling: get a data frame of species for 1) detected and 2) not detected
# Select only the 3 ARMS and delete the species that have an abundance <0
long_taxa_ASV_df_filtered <- long_df_variable %>%
  dplyr::filter(test %in% c("ARMS CO1", "ARMS 28S", "Expected equimolar")) %>%
  group_by(species) %>%
  mutate(sum_species = sum(value)) %>%
  filter(sum_species > 0) %>%
  ungroup()

# Compute the RRA per sample for each species 
factor_boxplot <- long_taxa_ASV_df_filtered %>%
  group_by(Sample) %>%
  mutate(sum_ASV_per_sample = sum(value),
         RRA = value / sum_ASV_per_sample * 100,
         RRA_per_sample = mean(RRA)) %>%
  ungroup()

# Add other taxonomic information (not just species) 
other_taxonomic_levels <- taxonomy_table_DNA %>%
  filter(rownames(.) %in% factor_boxplot$variable) %>%
  select(-species) %>%
  mutate(variable = rownames(.))
factor_boxplot <- left_join(factor_boxplot, other_taxonomic_levels, by = "variable")
factor_boxplot$ARMUNIT <- as.character(factor_boxplot$ARMUNIT)

# Compute RRA per sample at other taxonomic ranks (genus, family, order, class)
factor_boxplot_other_ranks <- factor_boxplot %>%
  group_by(Sample, genus) %>%
  mutate(sum_ASV_per_sample_genus = sum(value),
         RRA_genus = sum_ASV_per_sample_genus / sum_ASV_per_sample * 100) %>%
  group_by(Sample, family) %>%
  mutate(sum_ASV_per_sample_family = sum(value),
         RRA_family = sum_ASV_per_sample_family / sum_ASV_per_sample * 100) %>%
  group_by(Sample, order) %>%
  mutate(sum_ASV_per_sample_order = sum(value),
         RRA_order = sum_ASV_per_sample_order / sum_ASV_per_sample * 100) %>%
  ungroup()

# Round the numerical values of the table to 3 decimals
factor_boxplot_other_ranks <- factor_boxplot_other_ranks %>%
  mutate_if(is.numeric, round, digits = 3)

# Make a different vector for ARMS 20: Suberitida sp. 3 = Suberitidae sp. 3
Species_20.3 = c("Leucosolenida sp. 11 JV-2020",
                 "Leucosolenida sp. 13 JV-2020",
                 "Clathrinidae sp. 3 JV-2020",
                 "Clathrinidae sp. 5 JV-2020",
                 "Clathrinidae sp. 4 JV-2020",
                 "Leucosolenida sp. 1 JV-2020",
                 "Leucosolenida sp. 5 JV-2020",
                 "Leucosolenida sp. 7 JV-2020",
                 "Haliclona sp. 1 JV-2020",
                 "Ancorinidae sp. 3 JV-2020",
                 "Aplysilla rosea",
                 "Tethya sp. 3 JV-2020",
                 "Dysidea sp. 5 JV-2020",
                 "Poecilosclerida sp. 1 JV-2020",
                 "Poecilosclerida sp. 7 JV-2020",
                 "Suberitidae sp. 1 JV-2020",
                 "Suberitidae sp. 3 JV-2020",
                 "Plakina sp. 1 JV-2020")

# Combine all species in ARMS into a single vector
all_correct_species = unique(c(Species_17, Species_18, Species_20.3))

# Function to generate correctly and incorrectly identified vectors
generate_identified_vectors <- function(correct_species, all_species, factor_boxplot_other_ranks, arm_unit) {
  incorrectly_identified <- factor_boxplot_other_ranks %>%
    filter(ARMUNIT == arm_unit & !species %in% correct_species) %>%
    pull(species) %>%
    unique()
  
  correctly_identified <- factor_boxplot_other_ranks %>%
    filter(ARMUNIT == arm_unit & species %in% correct_species) %>%
    filter(RRA > 0 | RRA == 0 & species %in% all_species) %>%
    pull(species) %>%
    unique()
  
  list(incorrectly_identified = incorrectly_identified, correctly_identified = correctly_identified)
}

# Generate correctly and incorrectly identified vectors for each ARM unit
identified_17 <- generate_identified_vectors(Species_17, all_correct_species, factor_boxplot_other_ranks, "17")
identified_18 <- generate_identified_vectors(Species_18, all_correct_species, factor_boxplot_other_ranks, "18")
identified_20 <- generate_identified_vectors(Species_20.3, all_correct_species, factor_boxplot_other_ranks, "20")

# Combine all identified vectors into a list
identified_list <- list(identified_17, identified_18, identified_20)

# Label the species in diagnostic_df_test
diagnostic_df_test <- factor_boxplot_other_ranks %>%
  mutate(status = case_when(
    ARMUNIT == "17" & species %in% identified_list[[1]]$incorrectly_identified & RRA > 0 ~ "false positive",
    ARMUNIT == "17" & species %in% identified_list[[1]]$correctly_identified & RRA > 0 ~ "true positive",
    ARMUNIT == "17" & species %in% identified_list[[1]]$incorrectly_identified & RRA == 0 ~ "true negative",
    ARMUNIT == "17" & species %in% identified_list[[1]]$correctly_identified & RRA == 0 ~ "false negative",
    
    ARMUNIT == "18" & species %in% identified_list[[2]]$incorrectly_identified & RRA > 0 ~ "false positive",
    ARMUNIT == "18" & species %in% identified_list[[2]]$correctly_identified & RRA > 0 ~ "true positive",
    ARMUNIT == "18" & species %in% identified_list[[2]]$incorrectly_identified & RRA == 0 ~ "true negative",
    ARMUNIT == "18" & species %in% identified_list[[2]]$correctly_identified & RRA == 0 ~ "false negative",
    
    ARMUNIT == "20" & species %in% identified_list[[3]]$incorrectly_identified & RRA > 0 ~ "false positive",
    ARMUNIT == "20" & species %in% identified_list[[3]]$correctly_identified & RRA > 0 ~ "true positive",
    ARMUNIT == "20" & species %in% identified_list[[3]]$incorrectly_identified & RRA == 0 ~ "true negative",
    ARMUNIT == "20" & species %in% identified_list[[3]]$correctly_identified & RRA == 0 ~ "false negative",
    
    TRUE ~ "NO"
  )) %>%
  arrange(factor(order))

# Plot only true positives, false negatives, and false positives
heatmap_df <- diagnostic_df_test %>%
  mutate(status = if_else(status == "true positive", "detected", status)) %>%
  filter(status %in% c("false negative", "detected", "false positive")) %>%
  mutate(RRA2 = if_else(status == "false negative", 10, RRA),
         value = if_else(value <= 1, 0, value),
         status = if_else(value == 0 & status == "detected", "false negative", status),
         species2 = if_else(species %in% c("dropped", "unclassified"), variable, species),
         status = str_replace_all(status, c("false negative" = "not detected")),
         Sample = str_replace_all(Sample, c("ARMS2028S" = "C", "ARMS1828S" = "B", "ARMS1728S" = "A", "ARMS20CO1" = "C", "ARMS18CO1" = "B", "ARMS17CO1" = "A")),
         test = str_replace(test, "ARMS CO1", "ARMS COI")) %>%
  mutate(status = as.factor(status)) %>%
  filter(cocktail == "ARMS unit (3)")

# Make species name cleaner without the JV-2020
heatmap_df$species_short_name <- str_replace_all(heatmap_df$species2, "JV-2020", "")

# Select only sponges and only the species present in the ARMS for the heatmap
Only_sponges_heatmap_df <- heatmap_df %>%
  filter(class %in% c("Calcarea", "Homoscleromorpha", "Demospongiae"),
         species %in% all_correct_species,
         status %in% c("detected", "not detected")) %>%
  dplyr::arrange(factor(class), species)


# Make order vector for y axis
y.axis.order <- rev(unique(Only_sponges_heatmap_df$species_short_name))
y.axis.order <- c("Oscarella sp. 6 ", "Plakina sp. 1 ",
                  "Ancorinidae sp. 3 ", "Tethya sp. 5 ",  "Tethya sp. 3 ",
                  "Suberitidae sp. 3 " ,   "Hymeniacidon sp. 2 "  , "Halichondria sp. 1 " , "Suberitidae sp. 1 ",
                  "Poecilosclerida sp. 7 ", "Poecilosclerida sp. 1 ","Poecilosclerida sp. 2 ",
                  "Haplosclerida sp. 10 " ,  "Haliclona sp. 2 "    ,  "Haliclona sp. 1 " ,
                  "Ircinia sp. 2 "  ,      "Dysidea sp. 5 "     ,   "Aplysilla rosea"   ,    "Dendroceratida sp. 3 ",
                  "Leucosolenida sp. 13 ", "Leucosolenida sp. 3 ", "Leucosolenida sp. 5 "  ,  "Leucosolenida sp. 11 " ,  "Leucosolenida sp. 7 ", "Leucosolenida sp. 1 ",
                  "Clathrinida sp. 1 "  ,  "Clathrinidae sp. 5 " , "Clathrinidae sp. 3 " , "Clathrinidae sp. 4 ")
y.axis.order


## E.2 Data wrangling: false positives for heat map
# Filter and preprocess data for Venn analysis
VENN_sponges <- diagnostic_df_test %>%
  dplyr::mutate(status = str_replace(status, "true positive" , "detected")) %>%
  filter(status %in% c("false negative", "detected", "false positive"))  %>%
  mutate(RRA2 = ifelse(status == "false negative", 10, RRA),
         value = ifelse(value <= 1, 0, value),
         status = ifelse(value == 0 & status == "detected", "false negative", status),
         species2 = ifelse(species %in% c("dropped", "unclassified"), variable, species)) %>%
  dplyr::filter(!(value == 0 & status == "false positive")) %>%
  dplyr::mutate(status = str_replace(status, "false negative" , "not detected")) %>%
  arrange(factor(class))  %>%
  filter(class %in% c("Calcarea", "Homoscleromorpha", "Demospongiae"),
         cocktail == "ARMS unit (3)", RRA2 > 0) %>%
  mutate(status = as.factor(status)) 

# Clean species names
VENN_sponges$species_short_name <- str_replace_all(VENN_sponges$species2, "JV-2020", "")

# Summarize values for each ARM
VENN_sponges_summary <- VENN_sponges %>%
  group_by(Sample, status, test) %>%
  dplyr::summarise(n = n()) %>%
  select(Sample, n, test, status)

# Calculate the % of false positives for 28S and COI
VENN_sponge_fp_perc <- VENN_sponges_summary %>%
  group_by(test) %>%
  mutate(sum_species = sum(n)) %>%
  group_by(test, status) %>%
  dplyr::reframe(sum = sum(n) /sum_species * 100) %>%
  distinct()

# Total number of species per sample: it should be 13, 16 and 18
number_n_sponges_fp <- VENN_sponges_summary %>%
  group_by(Sample) %>%
  dplyr::summarise(total = sum(n))
number_n_sponges_fp

# Join with total species count and calculate the percentage of success: value for average % of false positives
VENN_sponges_pie_fp <- inner_join(VENN_sponges_summary,  VENN_sponges) %>%
  select(n, Sample, test, status, ARMUNIT, species) %>%
  dplyr::mutate(nunber_species = case_when ((Sample == "ARMS1728S")~ 19,
                                            (Sample == "ARMS17CO1")~ 16,
                                            (Sample == "ARMS1828S")~ 17,
                                            (Sample == "ARMS18CO1")~ 22,
                                            (Sample == "ARMS2028S")~ 30,
                                            (Sample == "ARMS20CO1")~ 25,
                                            TRUE ~ NA))%>%
  mutate(per_succ = (n / nunber_species) * 100) %>%
  distinct(per_succ, test, status) %>%
  ungroup() %>%
  # Aggregate mean values for pie charts
  group_by(test, status) %>%
  dplyr::summarise(avg_status = mean(per_succ, na.rm = TRUE)) %>%
  mutate(avg_status = round(avg_status, 1))

# Summarize false positive values
FP_values <- VENN_sponges_summary %>%
  filter(status == "false positive") %>%
  mutate(Sample = str_replace_all(Sample, c("ARMS2028S" = "C", "ARMS1828S" = "B", "ARMS1728S" = "A",
                                            "ARMS20CO1" = "C", "ARMS18CO1" = "B", "ARMS17CO1" = "A"))) %>%
  dplyr::mutate(test = str_replace(test, "ARMS CO1", "ARMS COI")) %>%
  ungroup() %>%
  dplyr::select(c(Sample, test, n)) 

# Add the false postives to the heatmap dataset and adjust the axes
Only_sponges_df_FP <- right_join(FP_values, Only_sponges_heatmap_df) %>%
  mutate(n = as.numeric(n),
         Sample.2 = case_when((primer == "28S" & Sample == "C") ~ "ARMS2028S",
                              (primer == "28S" & Sample == "B") ~ "ARMS1828S",
                              (primer == "28S" & Sample == "A") ~ "ARMS1728S",
                              (primer == "CO1" & Sample == "C") ~ "ARMS20CO1",
                              (primer == "CO1" & Sample == "B") ~ "ARMS18CO1",
                              (primer == "CO1" & Sample == "A") ~ "ARMS17CO1",
                              TRUE ~ NA_character_)) %>%
  mutate(Sample = as.factor(Sample),
         n = as.factor(n))


## E.3 Plot Heatmap
heatmap_status_white_sponges_shading_fp <- ggplot(Only_sponges_df_FP, aes(x = Sample, y = species_short_name, fill = test, alpha = status)) + 
  geom_tile(color = "black") +
  scale_alpha_discrete(range = c(1, 0.25), guide = guide_legend(override.aes = list(fill = "black"))) +
  theme_bw() +
  scale_y_discrete(limits = rev(y.axis.order)) +
  theme(axis.text = element_text(size = 10, family = "Open Sauce Sans"),
        axis.ticks = element_line(size = 0.4),
        strip.background = element_rect(colour = "black", fill = "white")) +
  facet_grid(~ factor(test), drop = TRUE, space = "free_x", scales = "free_x", col = NULL) +
  scale_x_discrete(labels = c("A\n(n =    )",  "B\n(n =     )", "C\n(n =     )"),  breaks = c("A", "B", "C")) +
  labs(x = "Sample number", y = "Sponge Species") +
  scale_fill_manual(breaks = c("ARMS 28S", "ARMS COI"), values = c("violetred1", "steelblue1"), na.value = "white") + 
  theme(legend.position = "none") +
  guides(fill = FALSE)

heatmap_status_white_sponges_shading_fp



# Summarize the false positives to add to table S9
sponges_fp <- VENN_sponges %>%
  filter(status == "false positive") %>%
  dplyr::mutate(test = str_replace(test, "ARMS CO1", "ARMS COI")) %>% 
  dplyr::select(c(species, Sample, value, species_short_name)) %>%
  pivot_wider(names_from = Sample, values_from = value) %>%
  dplyr::select(species_short_name, ARMS1728S, ARMS1828S, ARMS2028S, ARMS17CO1, ARMS18CO1, ARMS20CO1)

# Save it in excel  
write.csv(sponges_fp, file = "./tables/FalsePositives_Porifera_ARMS_MFPrarefied_050424.csv", row.names = FALSE)


## E.4 Make pie charts
# Summarise the number of species of each status for each sample
VENN_sponges_summary_nofp <- VENN_sponges %>%
  filter(status != "false positive") %>%
  group_by(Sample, status) %>%
  dplyr::summarise(n = n())

# Total number of species per sample: it should be 13, 16 and 18
number_n_sponges <- VENN_sponges_summary_nofp %>%
  group_by(Sample) %>%
  dplyr::summarise(total = sum(n))
number_n_sponges

# Join with total species count and calculate the percentage of success
VENN_sponges_pie <- inner_join(VENN_sponges_summary_nofp, VENN_sponges) %>%
  select(n, Sample, test, status, ARMUNIT, species) %>%
  mutate(nunber_species = case_when(Sample %in% c("ARMS1728S", "ARMS17CO1") ~ 13,
                                    Sample %in% c("ARMS1828S", "ARMS18CO1") ~ 16,
                                    Sample %in% c("ARMS2028S", "ARMS20CO1") ~ 18,
                                    TRUE ~ NA_real_)) %>%
  mutate(per_succ = (n / nunber_species) * 100) %>%
  distinct(per_succ, test, status) %>%
  ungroup()

# Aggregate mean values for pie charts
VENN_sponges_pie_mean <- VENN_sponges_pie %>%
  group_by(test, status) %>%
  dplyr::summarise(avg_status = mean(per_succ, na.rm = TRUE)) %>%
  mutate(avg_status = round(avg_status, 1))

# Split data for 28S and COI pie charts
pie_28S_sponges <- VENN_sponges_pie_mean %>%
  filter(test == "ARMS 28S") %>%
  arrange(avg_status) %>%
  mutate(ypos = cumsum(avg_status) - 0.5 * avg_status)

pie_COI_sponges <- VENN_sponges_pie_mean %>%
  filter(test == "ARMS CO1") %>%
  arrange(-(avg_status))%>%
  mutate(ypos = cumsum(avg_status) - 0.5 * avg_status)

# Plot 28S pie chart
pie_28S_sponges_white_nofp <- ggplot(pie_28S_sponges, aes(x = "", y = avg_status, fill = test)) + 
  geom_bar(aes(alpha = status), width = 1, stat = "identity", color = "black") +
  scale_alpha_discrete(range = c(1, 0.25), guide = "none") +
  scale_fill_manual(values = "violetred1") +
  theme_void() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top",
        text = element_text(family = "Open Sauce Sans"),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  coord_polar("y", start=0, direction = -1) +
  geom_text(aes(label = avg_status), position = position_stack(vjust = 0.5), color = "black", size = 6) +
  guides(alpha = "none", fill = "none")
pie_28S_sponges_white_nofp

# Plot COI pie chart
pie_COI_sponges_white_nofp <- ggplot(pie_COI_sponges, aes(x = "", y = avg_status, fill = test)) + 
  geom_bar(aes(alpha = status), width = 1, stat = "identity", color = "black") +
  scale_alpha_discrete(range = c(1, 0.25), guide = "none") +
  scale_fill_manual(values = "steelblue1") +
  theme_void() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white",  color = "white"),
        legend.position = "top",
        text = element_text(family = "Open Sauce Sans"),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  coord_polar("y") +
  geom_text(aes(label = avg_status, y = ypos), color = "black", size = 6) +
  guides(alpha = "none", fill = "none")
pie_COI_sponges_white_nofp




#################################################################
######### F. Figure S2: bar chart of taxonomic coverage #########
#################################################################

## F.1 Find how many species for are detected at each taxononic level for 28S, COI and expected data sets
# Data wrangling from status data frame
Taxonomic_coverage_df <- diagnostic_df_test %>%
  mutate(status = if_else(status == "true positive", "detected", status)) %>%
  filter(status %in% c("false negative", "detected", "false positive")) %>%
  mutate(RRA2 = if_else(status == "false negative", 10, RRA),
         value = if_else(value <= 1, 0, value),
         status = if_else(value == 0 & status == "detected", "false negative", status),
         species2 = if_else(species %in% c("dropped", "unclassified"), variable, species),
         status = str_replace_all(status, c("false negative" = "not detected")),
         #Sample = str_replace_all(Sample, c("ARMS2028S" = "C", "ARMS1828S" = "B", "ARMS1728S" = "A", "ARMS20CO1" = "C", "ARMS18CO1" = "B", "ARMS17CO1" = "A")),
         test = str_replace(test, "ARMS CO1", "ARMS COI"))%>%
  filter(class %in% c("Calcarea", "Homoscleromorpha", "Demospongiae")) %>%
  mutate(status = as.factor(status)) %>%
  dplyr::filter(status == "detected")

# Make species name cleaner without the JV-2020
Taxonomic_coverage_df$species_short_name <- str_replace_all(Taxonomic_coverage_df$species2, "JV-2020", "")

# Function to summarize values for different taxonomic categories
summarize_taxonomic_coverage <- function(data, group_var, taxa) {
  taxa_name <- deparse(substitute(taxa))  # Get the name of the taxa argument as a string
  data %>%
    dplyr::group_by({{ taxa }}, {{ group_var }}) %>%
    dplyr::filter(!(duplicated({{ group_var }}, {{ taxa }}))) %>%
    dplyr::group_by({{ group_var }}) %>%
    dplyr::summarise(n = length(unique(species))) %>%
    dplyr::mutate(taxa = taxa_name)
}

# Apply the function to each taxa column
taxo_cov_genus <- summarize_taxonomic_coverage(Taxonomic_coverage_df, Sample, genus)
taxo_cov_fam <- summarize_taxonomic_coverage(Taxonomic_coverage_df, Sample, family)
taxo_cov_order <- summarize_taxonomic_coverage(Taxonomic_coverage_df, Sample, order)
taxo_cov_class <- summarize_taxonomic_coverage(Taxonomic_coverage_df, Sample, class)
taxo_cov_species <- Taxonomic_coverage_df %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(taxa = "species")

# Get all objects in the environment that match the pattern "taxo_coverage"
matching_objects <- mget(ls(pattern = "taxo_cov"))

# Bind rows from all matching objects
Taxonomic_coverage_histogram <- bind_rows(matching_objects)

# Add unit number
UnitNumber <- Taxonomic_coverage_df %>%
  select(Sample, ARMUNIT, primer) %>%
  distinct(Sample, ARMUNIT, .keep_all = TRUE)

# Merge the ARM units together
Taxonomic_coverage_histogram_plus <- left_join(Taxonomic_coverage_histogram, UnitNumber) %>%
  group_by(primer, taxa) %>%
  mutate(n_average_arms = mean(n))


## F.2 Make histogram
barplot_taxonomic_coverage <- ggplot(Taxonomic_coverage_histogram_plus, aes(x = taxa, y = n_average_arms, fill = primer)) +
  geom_bar(stat= "identity", position = "dodge", width = 0.5, col= "black") +
  xlab("Taxa") +
  ylab("Average number of taxa detected per ARMS unit") +
  labs(fill = "Communities") + 
  theme_bw() +
  theme(text = element_text(family = "Open Sauce Sans")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank()) +
  scale_fill_manual(breaks = c("Expected", "28S", "CO1"), 
                    values = c("goldenrod", "violetred1", "steelblue1"),
                    labels = c("Expected", "28S", "COI")) 
barplot_taxonomic_coverage



##############################################################
########  G. Figure 5C: Histogram of classes detection #######
################################################### ##########

# G.1 Find how many species have been detected per order for the 28S, COI and expected data sets
# Select necessary columns and remove duplicated rows
Taxonomic_coverage_histo <- Taxonomic_coverage_df %>%
  select(species, genus, family, order, class, primer, test, species_short_name) %>%
  distinct() %>%
  dplyr::group_by(order, primer) %>%
  dplyr::summarise(n = n()) 

# Extract class level information
Taxonomy <- Taxonomic_coverage_df %>%
  select(order, primer, class) %>%
  distinct()

# Merge the summarize df to the taxonomic info
Taxonomic_coverage_histo_2 <-  left_join(Taxonomic_coverage_histo, Taxonomy, keep = NULL, multiple = "first") %>%
  # Add rows for Calcarea CO1 and modify row for Suberitidae sp. 3
  ungroup() %>%
  add_row(primer = "CO1", n = 0, order = c("Clathrinida", "Leucosolenida"), class = "Calcarea") %>%
  add_row(primer = "CO1", n = 0, order = "Dictyoceratida", class = "Demospongiae") %>%
  mutate(n = if_else(primer == "Expected" & order == "Suberitida", 4, n))

# Verify sum of all species: should be 29
Sum_expected <- Taxonomic_coverage_histo_2 %>%
  filter(primer == "Expected") %>%
  summarise(sum(n)) %>%
  pull()

## G.2 Plot everything
# Plot
barplot_tax_orders <- ggplot(Taxonomic_coverage_histo_2, aes(x = order, y = n, fill = primer)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, col = "black") +
  xlab("Sponge order") +
  ylab("Number of sponge species detected") +
  theme_bw() +
  labs(fill = "Communities") +
  theme(text = element_text(family = "Open Sauce Sans")) +
  facet_grid(. ~ class, drop = TRUE, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5), 
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 9)) +
  scale_fill_manual(breaks = c("Expected", "28S", "CO1"), 
                    labels = c("Expected", "ARMS 28S", "ARMS COI"), 
                    values = c("goldenrod", "violetred1", "steelblue1"))
barplot_tax_orders


####################################################
########  H. Make figure 5: combine 3 plots  #######
####################################################

## H.1 Combine the piecharts into one plot
pie_chart_pannel <- ggarrange(pie_28S_sponges_white_nofp, pie_COI_sponges_white_nofp,
                              ncol = 1,
                              nrow = 2,
                              common.legend = TRUE)
pie_chart_pannel

## H.2 Combine the pie charts and the heatmap
heatmap_pannel <- ggarrange(heatmap_status_white_sponges_shading_fp, pie_chart_pannel,
                            labels = c("A", "B"),
                            ncol = 2,
                            nrow = 1,
                            heights = c(1, 0.5),
                            widths = c(1, 0.5))
heatmap_pannel

## H.3 Combine with the barchart
pannel_plot_cocktail_4 <- ggarrange(heatmap_pannel, barplot_tax_orders,
                                    labels = c(" ","C"), 
                                    ncol = 1, 
                                    nrow = 2)+
  bgcolor("white") +
  border ("white")

pannel_plot_cocktail_4



############################################################################################
######## I. Make figure S3: bubble plot on read relative abundance, perform t-tests ########
############################################################################################

## I.1 Reload the data: delete singletons and change to read relative abundance (RRA)
# Remove singletons from combined df
combined_df_nosingletons = as.data.frame(combined_df)
columns_keep <- colnames(combined_df_nosingletons)[colnames(combined_df_nosingletons) %in% c("ARMS1728S", "ARMS1828S", "ARMS2028S",
                                                                                             "ARMS17CO1", "ARMS18CO1", "ARMS20CO1",
                                                                                             "Equimolar17", "Equimolar18", "Equimolar20",
                                                                                             "%cover17", "%cover18", "%cover20", "AllOTUs")]
columns_keep
for (i in columns_keep) {
  combined_df_nosingletons[, i] <- ifelse(combined_df_nosingletons[, i] > 1, combined_df_nosingletons[, i], 0)
}

# Remove species that have abundance <0 and the labeling mistake one (Suberitidae/a sp.3)
sponges_combined_table <- combined_df_nosingletons %>%
  mutate(ID_ASV = rownames(.)) %>%
  group_by(species) %>%
  mutate(sum_species = sum(c_across(where(is.numeric)))) %>%
  filter(sum_species > 0 & !species %in% c("Suberitidae sp. 3 JV-2020", "Suberitida sp. 3 JV-2020")) %>%
  select(-sum_species) %>%
  as.data.frame()

rownames(sponges_combined_table) = sponges_combined_table$ID_ASV
sponges_combined_table$ID_ASV <- NULL

# Make everything RRA in the combine table
# Replace NA values with 0
sponges_combined_table[is.na(sponges_combined_table)] <- 0

# Find numeric columns and calculate column sums
numeric_columns <- sapply(sponges_combined_table, is.numeric)
sum_of_reads <- colSums(sponges_combined_table[, numeric_columns])

# Divide numeric columns by their column sums
sponges_combined_table[, numeric_columns] <- sweep(sponges_combined_table[, numeric_columns], 2, sum_of_reads, "/")

# Round numeric columns to 4 decimal places
sponges_combined_table[, numeric_columns] <- round(sponges_combined_table[, numeric_columns], 4)


# Double check if equals to 1
sum_of_reads_RRA <- colSums(sponges_combined_table[, numeric_columns])
sum_of_reads_RRA


# Run basic functions to get data frames
RRA_ASV_table_DNA <- subset_numeric_columns(sponges_combined_table) 
RRA_taxonomy_table_DNA <- subset_character_columns(sponges_combined_table)
RRA_flipped_ASV_with_metadata <- merge_flipped_ASV_with_metadata(t(RRA_ASV_table_DNA), metadata_df_DNA)

RRA_long_df <- elongate_df(sponges_combined_table, metadata_df_DNA)
RRA_extract_taxa_col<- extract_taxa_col(RRA_taxonomy_table_DNA, "species", RRA_long_df)


## I.2 Data wrangling to get the abundance of each species for each cocktail
# Select only the desired cocktails and delete species with abundance < 0 for all Samples
RRA_long_taxa_ASV_df_filtered <- RRA_extract_taxa_col %>%
  filter(!(Sample %in% c("AllOTUs", paste0("Equimolar", 17:20), paste0("ExpectedEq", 17:20),
                         "Negative17", "Negative1820", paste0("ARMS", 17:20, "CO1")))) %>%
  group_by(species) %>%
  mutate(sum_species = sum(value)) %>%
  filter(sum_species > 0) %>%
  ungroup() %>%
  select(-sum_species)

# Filter necessary metadata
Meta_data_heatmap_primer <- metadata_df_DNA%>%
  filter(!(row.names(metadata_df_DNA) %in% c("AllOTUs", paste0("Equimolar", 17:20),
                                             paste0("ExpectedEq", 17:20), "Negative17",
                                             "Negative1820", paste0("ARMS", 17:20, "CO1"))))
# Combine long dataframe with metadata
heatmap_df_primers <- RRA_long_taxa_ASV_df_filtered %>%
  left_join(rownames_to_column(RRA_taxonomy_table_DNA, var = "variable"), by = "variable") %>%
  select(-species.y) %>%
  filter(species.x %in% unique(c(Species_17, Species_18, Species_20.2))) %>%
  arrange(ARMUNIT) %>%
  mutate(species.x = str_replace_all(species.x, "JV-2020", ""),
         cocktail2 = case_when(Sample %in% paste0("ARMS", 17:20, "28S") ~ "ARMS Unit",
                               Sample %in% paste0("%cover", 17:20) ~ "% cover-based DNA concentration mock community", 
                               TRUE ~ "Traced coverage"),
         value2 = round(value * 100, 2),
         ARMUNIT = str_replace_all(ARMUNIT, c("17" = "A", "18" = "B", "20" = "C"))) %>%
  filter(class %in% c("Calcarea", "Homoscleromorpha", "Demospongiae")) %>%
  arrange(class, order)

y.axis.order <- unique(heatmap_df_primers$species.x)


## I.3 Plot corncob plot
# Define the formatting function
format_size_labels <- function(x) {
  return(sprintf("%.2f", x))
}
corncob <- ggplot(heatmap_df_primers, aes(x = Sample, y = species.x, fill = order, size = value2)) + 
  geom_point(shape = 21, alpha = 1) +
  scale_size_continuous(trans = "log", labels = format_size_labels) +
  xlab("Sample") +
  ylab("Sponge species") + 
  theme_bw(base_size = 10) +
  scale_alpha(range = c(0.1, 2)) +
  labs(fill = "Order", size = "Sponge read relative abundance (%)") +
  scale_y_discrete(limits = y.axis.order) +
  facet_grid(~ factor(ARMUNIT), drop = TRUE, space = "free", scales = "free", col = NULL) +
  theme(text = element_text(family = "Open Sauce Sans"),
        axis.ticks = element_line(size = 0.4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(colour = "black", fill = "white")) +
  scale_fill_manual(values = c("darkcyan", "cadetblue1", "darkseagreen3", "palegoldenrod", "antiquewhite3",
                               "mistyrose", "darksalmon", "coral", "sienna3", "firebrick3", "tomato4", "white"),
                    breaks = c("Clathrinida", "Haplosclerida", "Homosclerophorida", "Leucosolenida",  
                               "Poecilosclerida", "Suberitida", "Tethyida", "Tetractinellida", "Dendroceratida",
                               "Dictyoceratida", "Axinellida", "unclassified")) + 
  scale_x_discrete(labels = c("% cover mock community", "ARMS unit", "Traced coverage")) +
  guides(fill = guide_legend(override.aes = list(size = 5)))


#######################################################################################
######## J. Make a taxonomic table that includes all OTUs for Mock community 2 ########
#######################################################################################

## J.1 Make taxonomic swaps or updated taxonomy
# If present, make the following swaps:
manual_swaps <- function(data, swaps) { 
  all_OTUs_taxonomy <- data
  
  for (swap in swaps) {
    all_OTUs_taxonomy  <- all_OTUs_taxonomy %>%
      mutate(species = case_when(species == swap$from ~ swap$to,
                                 TRUE ~ species
      ))
  }
  
  return(all_OTUs_taxonomy)
}
# Define swaps for each cocktail
swaps_all_OTUs <- list(
  list(from = "Ancorinidae sp. 1 JV-2020", to = "Stelletta sp. 1"),
  list(from = "Ancorinidae sp. 2 JV-2020", to = "Stelletta sp. 1"),
  list(from = "Ancorinidae sp. 7 JV-2020", to = "Stelletta sp. 5"),
  list(from = "Ancorinidae sp. 5 JV-2020", to = "Stelletta sp. 5"),
  list(from = "Haliclona sp. 23 JV-2020", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Heteroscleromorpha sp. 3 JV-2020", to = "Haplosclerida sp. 21"),
  list(from = "Haplosclerida sp. 12 JV-2020", to = "Haliclona (Soestella) caerulea"),
  # Other swaps that were wrong from GenBank
  list(from = "Biemna sp. P60", to = "Biemna fistulosa"),
  list(from = "Clionaida sp. 1 JV-2020", to = "Clionaida sp. 2 JV-2020"),
  list(from = "Dictyonella sp. SI_121Y", to = "Bubarida sp. 1 JV-2020"),
  list(from = "Clathrinidae sp. 2 JV-2020", to = "Clathrinida sp. 2 JV-2020"),
  list(from = "Mycale sp. SI_110Y", to = "Mycale sp. 14 JV-2020")
)
taxonomy_table_DNA_swaps <- as.data.frame(manual_swaps(taxonomy_table_DNA, swaps_all_OTUs))


# J.2 Delete duplicates in swap table
duplicated_species <- function (abundance, taxonomy){
  # J.1 Delete duplicates
  abundance <- abundance
  taxonomy <- taxonomy
  
  taxonomy_colnames <- colnames(taxonomy)
  # Remove "ID_ASV" from the column names
  taxonomy_colnames <- setdiff(taxonomy_colnames, "ID_ASV")
  
  row.names(taxonomy) = taxonomy$ID_ASV
  
  ## Verify if there are duplicates
  nrow(taxonomy) == length(unique(row.names(taxonomy)))
  
  ## Join combined tax table and abundance_global tables together
  a.combined_table <- inner_join(abundance, taxonomy, by ="ID_ASV")
  
  ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
  a.final.abundance.table <-  a.combined_table %>%
    dplyr::select(colnames(abundance))
  row.names(a.final.abundance.table)<- a.final.abundance.table$ID_ASV
  
  unique_taxa <- a.combined_table %>%
    distinct(species)
  
  ## J.2 Delete duplicates in both the abundance and tax tables
  # Find the number of unique taxa
  rownames(a.combined_table) = a.combined_table$ID_ASV
  
  extensive_taxonomy <- taxonomy
  
  # when species are duplicated, pick the one that has the most extensive taxonomy
  # For each duplicated species, return the number of "unclassified" in each row, if not return NA 
  extensive_taxonomy$unclassified <- ifelse(duplicated(extensive_taxonomy$species) | duplicated(extensive_taxonomy$species, fromLast = TRUE),rowSums(extensive_taxonomy == "unclassified"), NA)
  
  extensive_taxonomy_verified <- subset(extensive_taxonomy, !(duplicated(species) & unclassified != min(unclassified)))
  extensive_taxonomy_verified$unclassified <- NULL
  
  # Verify if there are still duplicated species: should have zero columns in total
  duplicated_species_still <- extensive_taxonomy_verified[which(duplicated(extensive_taxonomy_verified$ID_ASV)),]
  
  # J.3 Map those taxa
  unique_taxa <- extensive_taxonomy_verified %>%
    distinct(species)
  unique_taxa$mapping_ASV <- rownames(unique_taxa)
  
  number_unique_taxa <- nrow(unique_taxa) 
  
  # Remove duplicates from taxonomy
  # Delete the ones that are not mapping taxa
  taxonomy_mapping <- extensive_taxonomy_verified %>%
    dplyr::filter(ID_ASV %in% unique_taxa$mapping_ASV)
  
  abundance_mapping <- a.combined_table %>%
    dplyr::select(c(colnames(abundance), species)) 
  
  abundance_mapping_join <- right_join(unique_taxa, abundance_mapping) %>%
    dplyr::select(-c("ID_ASV", "species"))
  
  # Merge abundance by mapping ASV
  abundance_no_duplicates_2 = aggregate (. ~ abundance_mapping_join$mapping_ASV, abundance_mapping_join[,!(colnames(abundance_mapping_join) %in% c("mapping_ASV"))], sum)
  colnames(abundance_no_duplicates_2)[1] <- "ID_ASV"
  rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$ID_ASV
  
  # Check if all the mapping taxa are present
  if (identical(sort(unique_taxa$mapping_ASV), sort(abundance_no_duplicates_2$ID_ASV)))  {
    print( "V3 good" )
  } else {
    print ("V3 is wrong")
  }
  
  # J.4 Saved finished tables
  # Save that df in the working directory
  return(list(taxonomy_nodupes = taxonomy_mapping, abundance_nodupes = abundance_no_duplicates_2))
  
}

# Add both ID_ASV to the data frames on the left
taxonomy_table_DNA_swaps$ID_ASV <- rownames(taxonomy_table_DNA_swaps)
abundance_swaps = ASV_table_DNA
abundance_swaps$ID_ASV = rownames(abundance_swaps)
deduplicated <- duplicated_species(taxonomy = taxonomy_table_DNA_swaps, abundance = abundance_swaps)
deduplicated_taxonomy_table_DNA_swaps = deduplicated$taxonomy_nodupes
dedupliated_abundance <- deduplicated$abundance_nodupes

## J.3 Add missing OTUs to the abundance and taxonomy tables
# Find species that need to be added to the table (were not detected by metabarcoding but are supposed to be there)
All_OTUs <- c("Leucosolenida sp. 11 JV-2020",
              "Leucosolenida sp. 13 JV-2020",
              "Clathrinidae sp. 2 JV-2020",
              "Clathrinidae sp. 3 JV-2020",
              "Clathrinidae sp. 4 JV-2020",
              "Clathrinidae sp. 5 JV-2020",
              "Leucosolenida sp. 1 JV-2020",
              "Leucosolenida sp. 3 JV-2020",
              "Leucosolenida sp. 5 JV-2020",
              "Leucosolenida sp. 7 JV-2020",
              "Haliclona sp. 24 JV-2020",
              "Haliclona sp. 1 JV-2020",
              "Haliclona sp. 2 JV-2020",
              "Haplosclerida sp. 10 JV-2020",
              "Stelletta sp. 1",
              "Ancorinidae sp. 3 JV-2020",
              "Stelletta sp. 5",
              "Aplysilla rosea",
              "Dendroceratida sp. 3 JV-2020",
              "Tethya sp. 3 JV-2020",
              "Tethya sp. 5 JV-2020",
              "Chelonaplysilla sp. 1 JV-2020",
              "Dysidea sp. 5 JV-2020",
              "Mycale sp. 9 JV-2020",
              "Poecilosclerida sp. 7 JV-2020",
              "Poecilosclerida sp. 1 JV-2020",
              "Poecilosclerida sp. 2 JV-2020",
              "Suberitidae sp. 1 JV-2020",
              "Hymeniacidon sp. 1 JV-2020",
              "Hymeniacidon sp. 2 JV-2020",
              "Halichondria sp. 1 JV-2020",
              "Oscarella sp. 3 JV-2020",
              "Oscarella sp. 4 JV-2020",
              "Oscarella sp. 5 JV-2020",
              "Plakina sp. 1 JV-2020",
              "Plakortis sp. 1 JV-2020",
              "Ircinia sp. 2 JV-2020",
              "Cladocroce sp. 1 JV-2020",
              "Oscarella sp. 6 JV-2020",
              "Lissodendoryx sp. 2 JV-2020",
              "Clathrinidae sp. 1 JV-2020",
              "Leucosolenida sp. 10 JV-2020",
              "Leucosolenida sp. 2 JV-2020",
              "Haplosclerida sp. 1 JV-2020",
              "Haplosclerida sp. 2 JV-2020",
              "Haliclona sp. 3 JV-2020",
              "Haliclona sp. 4 JV-2020",
              "Haplosclerida sp. 3 JV-2020",
              "Callyspongia sp. 1 JV-2020",
              "Haplosclerida sp. 6 JV-2020",
              "Haplosclerida sp. 7 JV-2020",
              "Haplosclerida sp. 8 JV-2020",
              "Gelliodes wilsoni",
              "Stelletta sp. 3",
              "Tethya sp. 2 JV-2020",
              "Tethya sp. 4 JV-2020",
              "Chondrillida sp. 1 JV-2020",
              "Mycale sp. 10 JV-2020",
              "Monanchora clathrata",
              "Tedania sp. 1 JV-2020",
              "Tedania cf. klausi JV-2020",
              "Suberitidae sp. 2 JV-2020",
              "Suberitidae sp. 4 JV-2020",
              "Suberitidae sp. 5 JV-2020",
              "Hymeniacidon sp. 3 JV-2020",
              "Hymeniacidon sp. 5 JV-2020",
              "Halichondria sp. 2 JV-2020",
              "Terpios sp. 1 JV-2020",
              "Oscarella sp. 1 JV-2020",
              "Oscarella sp. 2 JV-2020",
              "Verongimorpha sp. 1 JV-2020",
              "Dysidea sp. 3 JV-2020",
              "Leucosolenida sp. 9 JV-2020",
              "Tethyidae sp. 1 JV-2020",
              "Lissodendoryx sp. 1 JV-2020",
              "Halichondria coerulea",
              "Poecilosclerida sp. 4 JV-2020",
              "Haplosclerida sp. 21",
              "Haliclona (Soestella) caerulea",
              "Chalinidae sp. 2 JV-2020",
              "Suberites sp. 1 JV-2020",
              "Suberitida sp. 3 JV-2020",
              "Keratosa sp. 1 JV-2020",
              "Chelonaplysilla erecta",
              "Dysidea cf. arenaria cf. JV-2020",
              "Mycale parishii",
              "Chelonaplysilla sp. 3 JV-2020",
              "Suberitidae sp. 6 JV-2020",
              "Ancorinidae sp. 8 JV-2020",
              "Haplosclerida sp. 15 JV-2020",
              "Poecilosclerida sp. 10 JV-2020",
              "Poecilosclerida sp. 11 JV-2020",
              "Mycale sp. 14 JV-2020",
              "Keratosa sp. 2 JV-2020",
              "Suberitida sp. 1 JV-2020",
              "Suberitida sp. 2 JV-2020",
              "Halichondria sp. 4 JV-2020",
              "Halichondria sp. 5 JV-2020",
              "Mycale sp. 4 JV-2020",
              "Mycale sp. 15 JV-2020",
              "Lissodendoryx hawaiiana",
              "Poecilosclerida sp. 9 JV-2020",
              "Poecilosclerida sp. 3 JV-2020",
              "Iotrochota protea",
              "Poecilosclerida sp. 6 JV-2020",
              "Poecilosclerida sp. 12 JV-2020",
              "Echinodictyum sp. 1 JV-2020",
              "Raspailiidae sp. 2 JV-2020",
              "Calcinea sp. 1 JV-2020",
              "Tethya sp. 1 JV-2020",
              "Dictyoceratida sp. 5 JV-2020",
              "Dictyoceratida sp. 4 JV-2020",
              "Dictyoceratida sp. 6 JV-2020",
              "Dictyoceratida sp. 3 JV-2020",
              "Dictyoceratida sp. 2 JV-2020",
              "Dysidea cf. pallescens JV-2020",
              "Dysidea sp. 6 JV-2020",
              "Lamellodysidea cf. chlorea JV-2020",
              "Dysideidae sp. 1 JV-2020",
              "Dysidea sp. 9 JV-2020",
              "Chondrosiida sp. 1 JV-2020",
              "Spirastrella sp. 1 JV-2020",
              "Spheciospongia solida",
              "Axinellida sp. 1 JV-2020",
              "Pseudoceratina sp. 1 JV-2020",
              "Chondrillida sp. 3 JV-2020",
              "Biemna fistulosa",
              "Bubarida sp. 1 JV-2020",
              "Leucosolenida sp. 12 JV-2020",
              "Clathrinidae sp. 8 JV-2020",
              "Clathrinidae sp. 9 JV-2020",
              "Clathrinidae sp. 10 JV-2020",
              "Leucettidae sp. 3 JV-2020",
              "Clathrinidae sp. 11 JV-2020",
              "Leucosolenida sp. 19 JV-2020",
              "Leucosolenida sp. 20 JV-2020",
              "Leucosolenida sp. 21 JV-2020",
              "Haliclona sp. 21 JV-2020",
              "Clionaida sp. 2 JV-2020",
              "Callyspongia sp. 2 JV-2020",
              "Oscarella sp. 7 JV-2020",
              "Plakinastrella sp. 1 JV-2020",
              "Corticium sp. 1 JV-2020",
              "Haplosclerida sp. 14 JV-2020",
              "Suberitida sp. 4 JV-2020",
              "Tethya sp. 6 JV-2020")

species_not_present <- All_OTUs[which(!All_OTUs %in% deduplicated_taxonomy_table_DNA_swaps$species)]
view(species_not_present)

# Missing taxonomy
taxonomy_values <- list ("Clathrinidae sp. 2 JV-2020" = c("Clathrinida", "Clathrinidae","unclassified"),
                         "Oscarella sp. 4 JV-2020" = c("Homosclerophorida", "Oscarellidae", "Oscarella"),
                         "Oscarella sp. 7 JV-2020" = c("Homosclerophorida", "Oscarellidae", "Oscarella"),
                        "Lissodendoryx sp. 2 JV-2020" = c("Poecilosclerida","Coelosphaeridae", "Lissodendoryx"),
                        "Clathrinidae sp. 1 JV-2020" = c("Clathrinida", "Clathrinidae", "unclassified"),
                        "Leucosolenida sp. 2 JV-2020" = c("Leucosolenida", "unclassified", "unclassified"),
                        "Haplosclerida sp. 1 JV-2020" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Haliclona sp. 3 JV-2020" = c("Haplosclerida", "Chalinidae", "Haliclona"),
                        "Haplosclerida sp. 3 JV-2020" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Haplosclerida sp. 6 JV-2020" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Haplosclerida sp. 7 JV-2020" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Stelletta sp. 3" = c("Tetractinellida", "Ancorinidae", "unclassified"),
                        "Tethya sp. 2 JV-2020" = c("Tethyida", "Tethyidae", "Tethya"),
                        "Tethya sp. 4 JV-2020" = c("Tethyida", "Tethyidae", "Tethya"),
                        "Suberitidae sp. 4 JV-2020" = c("Suberitida", "Suberitidae", "unclassified"),
                        "Suberitidae sp. 5 JV-2020" = c("Suberitida", "Suberitidae", "unclassified"),
                        "Terpios sp. 1 JV-2020" = c ("Suberitida", "Suberitida incertae sedis", "Terpios"),
                        "Verongimorpha sp. 1 JV-2020" = c ("unclassified", "unclassified", "Verongimorpha"),
                        "Leucosolenida sp. 9 JV-2020" = c("Leucosolenida", "unclassified", "unclassified"),
                        "Haplosclerida sp. 21" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Haliclona (Soestella) caerulea" = c("Haplosclerida", "Chalinidae", "Haliclona"),
                        "Dysidea cf. arenaria cf. JV-2020" = c("Dictyoceratida", "Dysideidae", "Dysidea"),
                        "Mycale parishii" = c("Poecilosclerida", "Mycalidae", "Mycale"),
                        "Suberitidae sp. 6 JV-2020" = c("Suberitida", "Suberitidae", "unclassified"),
                        "Haplosclerida sp. 15 JV-2020" = c("Haplosclerida", "unclassified", "unclassified"),
                        "Poecilosclerida sp. 11 JV-2020" = c("Poecilosclerida", "unclassified", "unclassified"),
                        "Mycale sp. 15 JV-2020" = c("Poecilosclerida", "Mycalidae", "Mycale"),
                        "Iotrochota protea" = c("Poecilosclerida", "Iotrochotidae", "Iotrochota"),
                        "Poecilosclerida sp. 6 JV-2020" = c("Poecilosclerida", "unclassified", "unclassified"),
                        "Poecilosclerida sp. 12 JV-2020" = c("Poecilosclerida", "unclassified", "unclassified"),
                        "Raspailiidae sp. 2 JV-2020" = c ("Axinellida", "Raspailiidae", "unclassified"),
                        "Tethya sp. 1 JV-2020" = c("Tethyida", "Tethyidae", "Tethya"),
                        "Dictyoceratida sp. 5 JV-2020" = c("Dictyoceratida", "unclassified", "unclassified"),
                        "Dictyoceratida sp. 6 JV-2020" = c("Dictyoceratida", "unclassified", "unclassified"),
                        "Dictyoceratida sp. 2 JV-2020" = c("Dictyoceratida", "unclassified", "unclassified"),
                        "Dysidea sp. 6 JV-2020" = c("Dictyoceratida", "Dysideidae", "Dysidea"),
                        "Lamellodysidea cf. chlorea JV-2020" = c("Dictyoceratida", "Dysideidae", "Lamellodysidea"),
                        "Dysidea sp. 9 JV-2020" = c("Dictyoceratida", "Dysideidae", "Dysidea"),
                        "Chondrosiida sp. 1 JV-2020" = c("Chondrosiida", "unclassified", "unclassified"),
                        "Pseudoceratina sp. 1 JV-2020" = c("Verongiida", "Pseudoceratinidae", "Pseudoceratina"),
                        "Leucettidae sp. 3 JV-2020" = c("Clathrinida", "Lecuettidae", "unclassified"),
                        "Clathrinidae sp. 11 JV-2020" = c("Clathrinida", "Clathrinidae", "unclassified"),
                        "Leucosolenida sp. 20 JV-2020" = c("Leucosolenida", "unclassified", "unclassified"),
                        "Corticium sp. 1 JV-2020" = c("Homosclerophorida","Plakinidae", "Corticium"),
                        "Suberitida sp. 4 JV-2020" = c("Suberitida", "Suberitidae", "unclassified"))
                


# Check if all missing species are added to the table
species_added <- as.vector(names(taxonomy_values))
missing_ones <- species_added[which(!species_added %in% species_not_present)]
missing_ones
missing_ones_2 <- species_not_present[which(!species_not_present %in% species_added)]
missing_ones_2
length(species_added)
length(species_not_present)
sort(species_added) == sort(species_not_present)


# Add missing species to the taxonomic and abundance tables
add_missing_species <- function(taxonomy_table, missing_species, taxonomy_values, ASVnumber, abundance, All_OTUs) {
  rows <- nrow(taxonomy_table)
  num_missing <- length(missing_species)
  taxonomy_table[rows + num_missing, ] <- NA
  # Create a sequence of ASV numbers
  asv_numbers <- paste0("ASV", ASVnumber:(ASVnumber + num_missing - 1))
  taxonomy_table[(rows + 1):nrow(taxonomy_table), c("ID_ASV", "species")] <- c(asv_numbers, missing_species)
  
  for (i in 1:num_missing) {
    if (missing_species[i] %in% names(taxonomy_values)) {
      taxonomy_table[which(taxonomy_table$species == missing_species[i]), c("order", "family", "genus")] <- taxonomy_values[[missing_species[i]]]
    } else {
      warning(paste("No taxonomy values found for", missing_species[i]))
    }
    taxonomy_table[, c("kingdom")]<- c("Metazoa")
    taxonomy_table[, c("phylum")]<- c("Porifera")
  }
  # Update abundance table
  abundance$ID_ASV = rownames(abundance)
  abundance[rows + num_missing, ] <- NA
  abundance[(rows + 1):nrow(abundance), c("ID_ASV")] <- c(asv_numbers)
  abundance[is.na(abundance)] <- 0 
  ID_ASV_needed <- as.vector(taxonomy_table[which(taxonomy_table$species %in% All_OTUs), "ID_ASV"])
  abundance_2 <- abundance %>%
    dplyr::mutate(ExpectedAllOTUs = case_when((ID_ASV  %in% ID_ASV_needed) ~ "10",
                                              TRUE ~ "0"))
  row.names(abundance_2)= abundance_2$ID_ASV
  # Return the updated taxonomy table
  return(list(abundance_missing = abundance_2, taxonomy_table_missing = taxonomy_table))
}

All_OTUs_table <-  add_missing_species(deduplicated_taxonomy_table_DNA_swaps, species_not_present, taxonomy_values, 20000, dedupliated_abundance, All_OTUs)
abundance_all_OTUs <- All_OTUs_table$abundance_missing
taxonomy_all_OTUs <- All_OTUs_table$taxonomy_table_missing

# Save those tables
write.csv(abundance_all_OTUs, "./abundance_all_OTUs_070224.csv", row.names = FALSE)
write.csv(taxonomy_all_OTUs, "./taxonomy_all_OTUs_070224.csv", row.names = FALSE)



#####################################################################################
######## K. Make Figure 6: sunburst plot for COI and 28S: all metazoan phyla ########
#####################################################################################

## K.1 Data wrangling to obtain data frame with the % of species in each class and phylum
# Load Re-load data but include all OTUs: not only Porifera
Metazoans_combined_df <- read_and_prepare_data("combined28SCOI_taxonomy_table_071224.csv",
                                               "combined28SCOI_abundance_table_071224.csv", 
                                               "Martineau_TestingPrimers_CO1AND28S_metadata.csv")
Metazoans_ASV_table_DNA <- subset_numeric_columns(Metazoans_combined_df)
Metazoans_taxonomy_table_DNA <- subset_character_columns(Metazoans_combined_df)
Metazoans_flipped_ASV_with_metadata <- merge_flipped_ASV_with_metadata(t(Metazoans_ASV_table_DNA), metadata_df_DNA)
Metazoans_flipped_ASV_with_metadata$ARMUNIT = as.character(Metazoans_flipped_ASV_with_metadata$ARMUNIT)

Metazoans_long_df <- elongate_df(Metazoans_combined_df, metadata_df_DNA)
Metazoans_long_df_variable <- extract_taxa_col(Metazoans_taxonomy_table_DNA, "species", Metazoans_long_df)

# Find species for which phylum was unclassified (can do this earlier in a different script)
unclassified <- rownames(Metazoans_taxonomy_table_DNA[which(Metazoans_taxonomy_table_DNA$phylum == "unclassified"),])


# Select only the 3 ARMS units and delete species that have abundane < 0 for all 3 ARMS
long_taxa_ASV_df_filtered <- Metazoans_long_df_variable %>%
  filter(test %in% c("ARMS CO1", "ARMS 28S")) %>%
  # Remove the ASVs for which phylum is unclassified
  filter (!(variable %in% unclassified)) %>%
  dplyr::mutate(value = ifelse(value <= 1, 0, value)) %>%
  group_by(species) %>%
  mutate(sum_species = sum(value)) %>%
  filter(sum_species > 0) %>%
  ungroup() %>%
  # compute the RRA per Primer and add another column with it
  group_by(primer) %>%
  mutate(sum_ASV_per_sample = sum(value)) %>%
  group_by(primer, species) %>%
  mutate(sum_species_primer = sum(value)) %>%
  mutate(RRA = sum_species_primer / sum_ASV_per_sample * 100) %>%
  mutate(RRA_per_sample = mean(RRA)) %>%
  ungroup()

# Add the other taxonomic groups to data frame: genus, class, order, family, phylum
Metazoans_2_other_taxonomic_levels <- Metazoans_taxonomy_table_DNA %>%
  filter(rownames(.) %in% long_taxa_ASV_df_filtered$variable) %>% 
  #  mutate_at(vars(class:phylum), ~replace(., is.na(.), 0)) %>%
  mutate(variable = rownames(.)) %>%
  select(-(species))

# Merge the other taxonomic levels
Metazoans_factor_sunplot <- join(long_taxa_ASV_df_filtered, Metazoans_2_other_taxonomic_levels, by = "variable") %>%
  mutate(ARMUNIT = as.character(ARMUNIT))

## Make another data frame with RRA for other taxonomic ranks (class, phylum)
factor_sunplot_other_ranks <- Metazoans_factor_sunplot %>%
  group_by(class, phylum, primer) %>%
  mutate(sum_ASV_per_sample_class = sum(value)) %>%
  mutate(RRA_class = sum_ASV_per_sample_class / sum_ASV_per_sample * 100) %>%
  group_by(phylum, primer) %>%
  mutate(sum_ASV_per_sample_phylum = sum(value)) %>%
  mutate(RRA_phylum = sum_ASV_per_sample_phylum / sum_ASV_per_sample * 100) %>%
  ungroup()


## K.2 Sunplot 28S
# Add the % of species per class and phylum
factor_sunplot_28S <- factor_sunplot_other_ranks %>%
  dplyr::select(c(phylum, class, order, family, genus, species, primer, RRA_class, RRA))%>%
  dplyr::filter (primer == "28S") %>%
  dplyr::distinct() %>%
  dplyr::filter(RRA > 0) %>%
  dplyr::ungroup() %>%
  # Find the number of unique species (no species are duplicated, so it is equal to the number of rows)
  dplyr::mutate (sum_species_sample = n()) %>%
  # Find the % of each species that belong to each class and phylum
  dplyr::group_by (phylum) %>%
  dplyr::mutate(sum_species_phylum = n(),
                perc_species_phylum = sum_species_phylum/sum_species_sample *100) %>%
  # Do the same for classes 
  dplyr::group_by (class,phylum) %>%
  dplyr::mutate(sum_species_class = n(),
                perc_species_class = sum_species_class/sum_species_sample *100)%>%
  ungroup() %>%
  # Delete unecessary columns and duplicated rows
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(c(class, phylum,  perc_species_class, perc_species_phylum, RRA_class, sum_species_class, sum_species_phylum)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(desc(perc_species_phylum)) 


## Pick colour pallette for each Phylum
n_phyla <- length(unique(factor_sunplot_28S$phylum))

colors_phylum_28S <- 
  c ("006666",
     "plum1",
     "#A6CEE3",
     "lightgoldenrod",
     "lightgreen",
     "#B294C7" ,
     "bisque3",
     "lemonchiffon1",
     "mediumturquoise",
     "thistle3",
     "lightskyblue4",
     "mediumpurple4",
     "darksalmon",
     'rosybrown',
     "darkorange1",
     "yellow4") %>%
  set_names(unique(factor_sunplot_28S[["phylum"]]))
colors_phylum_28S
test28s <- unique(factor_sunplot_28S$phylum)
test28s
view(colors_phylum_28S)
# Outer ring (Phylum)
ring_outer_28S <- factor_sunplot_28S %>%
  # Add colours to df
  left_join (enframe(colors_phylum_28S, name = "phylum", value = "color_phylum"), by = "phylum") %>%
  # Arrange descending phylum by ascending class for compute title color
  dplyr::arrange(desc(perc_species_phylum), phylum, perc_species_class, class) %>%
  group_by (phylum) %>%
  # Number of class per phylum 
  dplyr::mutate (id_class = row_number(),
                 num_class = max(id_class)) %>%
  dplyr::ungroup() %>%
  # Degree of transparency
  dplyr::mutate (color_trans = id_class / num_class,
                 color_class = map2_chr(color_phylum, color_trans, ~ adjustcolor (.x, .y))) %>%
  # Order the figure clockwise
  dplyr::arrange(-perc_species_phylum, phylum, -perc_species_class, desc(class)) %>%
  dplyr::mutate(phylum = fct_inorder(phylum),
                class = fct_inorder(class)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(phylum) %>%
  dplyr:: mutate(cum_classpct = cumsum(perc_species_class))%>%
  dplyr::ungroup() %>%
  dplyr::mutate(cum_phylapct = cumsum(cum_classpct),
                cum_all = cumsum (perc_species_class)) %>%
  # Compute coordinates of the rectangles 
  # 0.3 just adds small amount of white space between rectangles
  dplyr::mutate(rect_x_max = cumsum(perc_species_class) - 0.3,
                rect_x_min = rect_x_max - perc_species_class + 0.3,
                rect_x_mid = (rect_x_min+ rect_x_max)/2,
                # Angles for adding text
                angle = 90 - 360 * rect_x_mid/100,
                hjust = ifelse (angle < -90, 1, 0),
                angle = ifelse(angle < -90, angle + 180, angle),
                # Labels for gggraph
                label_class = glue::glue( '{scales::percent(perc_species_phylum/100, accuracy = 0.001)} {phylum} {scales::percent(perc_species_class/100, accuracy = 0.001)} {class}\n '
                ),
                # add a column with log10 of class RRA
                log10classRRA = log(RRA_class + 1),
                sqrtclassRRA =  sqrt(RRA_class + 1)
  )

# Combine phylum and class color in a single vector
colors_class_28S <- ring_outer_28S[["color_class"]] %>% 
  set_names(ring_outer_28S[["class"]])
all_colors_28S <- c(colors_phylum_28S, colors_class_28S)

# Inner ring
ring_inner_28S <- ring_outer_28S %>%
  select(phylum, perc_species_phylum, rect_x_max, rect_x_min) %>%
  dplyr::arrange(phylum, desc(rect_x_max)) %>%
  dplyr::group_by(phylum) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(phylum, perc_species_phylum, stop = rect_x_max) %>%
  left_join(ring_outer_28S) %>%
  dplyr::select(phylum, perc_species_phylum, rect_x_max, rect_x_min, stop) %>%
  dplyr::arrange(phylum, rect_x_max) %>%
  dplyr::group_by(phylum) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  select(phylum, perc_species_phylum, start = rect_x_min, stop) %>%
  bind_rows(tibble(attr_category = NA, start = 99.7, stop = 100.1)) %>%
  mutate(mid = (stop - start) / 2 + start,
         angle = 90 - 360 * mid / 100,
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle),
         label_phylum = ifelse(is.na(phylum), NA, glue::glue('{scales::percent(perc_species_phylum / 100, accuracy = 0.001)} {phylum} ')))

# Add colour to all_colours_vector
test_colours_phylum28S <- c(colors_phylum_28S,  "goldenrod1","yellow") %>% 
  set_names(c(unique(factor_sunplot_28S[["phylum"]]), "Arthropoda", "Haptophyta"))

test_colours_phylum28S

all_colors_test_28S <- c(test_colours_phylum28S, colors_class_28S)
unique(factor_sunplot_28S[["phylum"]])

# Plot 28S sunburst
ARMS_28S_sunburst_prop <- ggplot(ring_outer_28S) + 
  geom_rect(
    aes(ymin = 11, ymax = 11 + sum_species_class, xmin = rect_x_min, xmax = rect_x_max, fill = class),
    show.legend = FALSE) +
  coord_polar() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = 'white', colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(data = ring_outer_28S, aes(rect_x_mid, y = sum_species_class + 30, label = class, angle = angle), size = 3) +
  theme(text = element_text(family = "Open Sauce Sans")) +
  ylim(-40, 58) +
  geom_rect(data = ring_inner_28S, aes(ymin = -35, ymax = 10, xmin = start, xmax = stop, fill = phylum), show.legend = TRUE) +
  scale_fill_manual(values = c(all_colors_test_28S), limits = c(levels(ring_outer_28S$phylum), levels(ring_outer_28S$class), "Arthropoda", "Haptophyta"), breaks = c(levels(ring_outer_28S$phylum), "Arthropoda", "Haptophyta"), drop = FALSE, guide = guide_legend(override.aes = list(linetype = c(1), color = "black"))) +
  labs(fill = "phylum") +
  geom_text(data = ring_inner_28S, aes(mid, -5, label = phylum, angle = angle), size = 2.5) 
ARMS_28S_sunburst_prop

# Extract some numbers for analyses
Number_phyla_detected_28S <- length(unique(factor_sunplot_28S$phylum))
Number_classes_detected_28S <- nrow(ring_outer_28S)  

excel_results_28S <- factor_sunplot_28S %>%
  dplyr::mutate(Number_phyla_detected = Number_phyla_detected_28S,
                Number_classes_detected = Number_classes_detected_28S)
write.csv(excel_results_28S,"./tables/Martineau_TestingPrimers_28S_sunplot_abundance_results_070224.csv", row.names =  FALSE)


## K.3 Sunplot COI
# Add the % of species per class and phylum
factor_sunplot_CO1 <- factor_sunplot_other_ranks %>%
  dplyr::select(c(phylum, class, order, family, genus, species, primer, RRA_class, RRA))%>%
  dplyr::filter (primer == "CO1") %>%
  dplyr::distinct() %>%
  dplyr::filter(RRA > 0) %>%
  dplyr::ungroup() %>%
  # Find the number of unique species (no species are duplicated, so it is equal to the number of rows)
  dplyr::mutate (sum_species_sample = n()) %>%
  # Find the % of each species that belong to each class and phylum
  dplyr::group_by (phylum) %>%
  dplyr::mutate(sum_species_phylum = n(),
                perc_species_phylum = sum_species_phylum/sum_species_sample *100) %>%
  # Do the same for classes 
  dplyr::group_by (class,phylum) %>%
  dplyr::mutate(sum_species_class = n(),
                perc_species_class = sum_species_class/sum_species_sample *100)%>%
  ungroup() %>%
  # Delete unecessary columns and duplicated rows
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(c(class, phylum,  perc_species_class, perc_species_phylum, RRA_class, sum_species_class, sum_species_phylum)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(desc(perc_species_phylum)) %>%
  dplyr::mutate(class_unique =  glue::glue("{phylum}/{class}"))

## Pick colour pallette for each Phylum
n_phyla <- length(unique(factor_sunplot_CO1$phylum))

colors_phylum_CO1 <- 
  c ("plum1",
     "006666",
     "#A6CEE3",
     "lightgreen",
     "lightgoldenrod",
     "#B294C7" ,
     "goldenrod1",
     "rosybrown",
     "thistle3", 
     "mediumpurple4",
     "mediumturquoise",
     "bisque3") %>%
  set_names(unique(factor_sunplot_CO1[["phylum"]]))
view(colors_phylum_CO1)
testCO1 <- unique(factor_sunplot_CO1$phylum)
testCO1
colors_phylum_CO1
# Outer ring (Phylum)
ring_outer_CO1 <- factor_sunplot_CO1 %>%
  # Add colours to df
  left_join (enframe(colors_phylum_CO1, name = "phylum", value = "color_phylum"), by = "phylum") %>%
  # Arrange descending phylum by ascending class for compute title color
  dplyr::arrange(desc(perc_species_phylum), phylum, perc_species_class, class_unique) %>%
  group_by (phylum) %>%
  # Number of class per phylum 
  dplyr::mutate (id_class = row_number(),
                 num_class = max(id_class)) %>%
  dplyr::ungroup() %>%
  # Degree of transparency
  dplyr::mutate (color_trans = id_class / num_class,
                 color_class = map2_chr(color_phylum, color_trans, ~ adjustcolor (.x, .y))) %>%
  # Order the figure clockwise
  dplyr::arrange(-perc_species_phylum, phylum, -perc_species_class, desc(class)) %>%
  dplyr::mutate(phylum = fct_inorder(phylum),
                class = fct_inorder(class)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(phylum) %>%
  dplyr:: mutate(cum_classpct = cumsum(perc_species_class))%>%
  dplyr::ungroup() %>%
  dplyr::mutate(cum_phylapct = cumsum(cum_classpct),
                cum_all = cumsum (perc_species_class)) %>%
  # Compute coordinates of the rectangles 
  # 0.3 just adds small amount of white space between rectangles
  dplyr::mutate(rect_x_max = cumsum(perc_species_class) - 0.3,
                rect_x_min = rect_x_max - perc_species_class + 0.3,
                rect_x_mid = (rect_x_min+ rect_x_max)/2,
                # Angles for adding text
                angle = 90 - 360 * rect_x_mid/100,
                hjust = ifelse (angle < -90, 1, 0),
                angle = ifelse(angle < -90, angle + 180, angle),
                # Labels for gggraph
                label_class = glue::glue( '{scales::percent(perc_species_phylum/100, accuracy = 0.001)} {phylum} {scales::percent(perc_species_class/100, accuracy = 0.001)} {class}\n '
                ),
                # add a column with log10 of class RRA
                log10classRRA = log(RRA_class + 1),
                sqrtclassRRA =  sqrt(RRA_class + 1)
  )

# Combine phylum and class color in a single vector
colors_class_CO1 <- ring_outer_CO1[["color_class"]] %>% 
  set_names(ring_outer_CO1[["class_unique"]])
all_colors_CO1 <- c(colors_phylum_CO1, colors_class_CO1)

# Inner ring
ring_inner_CO1 <- ring_outer_CO1 %>%
  select(phylum, perc_species_phylum, rect_x_max, rect_x_min) %>%
  dplyr::arrange(phylum, desc(rect_x_max)) %>%
  dplyr::group_by(phylum) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(phylum, perc_species_phylum, stop = rect_x_max) %>%
  left_join(ring_outer_CO1) %>%
  dplyr::select(phylum, perc_species_phylum, rect_x_max, rect_x_min, stop) %>%
  dplyr::arrange(phylum, rect_x_max) %>%
  dplyr::group_by(phylum) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  select(phylum, perc_species_phylum, start = rect_x_min, stop) %>%
  bind_rows(tibble(attr_category = NA, start = 99.7, stop = 100.1)) %>%
  mutate(mid = (stop - start) / 2 + start,
         angle = 90 - 360 * mid / 100,
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle),
         label_phylum = ifelse(is.na(phylum), NA, glue::glue('{scales::percent(perc_species_phylum / 100, accuracy = 0.001)} {phylum} ')))

# Plot CO1 sunburst
ARMS_CO1_sunburst_prop <- ggplot(ring_outer_CO1) + 
  geom_rect(
    aes(ymin = 11, ymax = 11 + sum_species_class, xmin = rect_x_min, xmax = rect_x_max, fill = class_unique),
    show.legend = FALSE) +
  coord_polar() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = 'white', colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(data = ring_outer_CO1, aes(rect_x_mid, y = sum_species_class + 30, label = class, angle = angle), size = 3) +
  theme(text = element_text(family = "Open Sauce Sans")) +
  ylim(-40, 58) +
  geom_rect(data = ring_inner_CO1, aes(ymin = -35, ymax = 10, xmin = start, xmax = stop, fill = phylum), show.legend = TRUE) +
  
  
  
  scale_fill_manual(values = c(all_colors_CO1), breaks = ring_outer_CO1$phylum, guide = guide_legend(override.aes = list(linetype = c(1), color = "black"))) +
  labs(fill = "phylum") +
  geom_text(data = ring_inner_CO1, aes(mid, -5, label = phylum, angle = angle), size = 2.5) 
ARMS_CO1_sunburst_prop


# Extract some numbers for analyses
Number_phyla_detected_CO1 <- length(unique(factor_sunplot_CO1$phylum))
Number_classes_detected_CO1 <- nrow(ring_outer_CO1)  

excel_results_CO1 <- factor_sunplot_CO1 %>%
  dplyr::mutate(Number_phyla_detected = Number_phyla_detected_CO1,
                Number_classes_detected = Number_classes_detected_CO1)
write.csv(excel_results_CO1,"./tables/Martineau_TestingPrimers_CO1_sunplot_abundance_results_070224.csv", row.names =  FALSE)


## K.4 Combine plots
sunplot_COI <- ARMS_CO1_sunburst_prop
sunplot_28S = ARMS_28S_sunburst_prop

pannel_sunplot <- ggarrange(sunplot_28S, sunplot_COI,
                            ncol=2,
                            labels = c("A ","B"),
                            common.legend = TRUE, 
                            legend = "right")+
  bgcolor("white") +
  rremove ("legend") +
  border ("white") 
pannel_sunplot


## K.5 Make data frame for t-test: species richness by unit
ttest_species_across_units <- Metazoans_long_df_variable %>%
  filter(test %in% c("ARMS CO1", "ARMS 28S")) %>%
  # Remove the ASVs for which phylum is unclassified
  filter (!(variable %in% unclassified)) %>%
  dplyr::mutate(value = ifelse(value <= 1, 0, value)) %>%
  filter(value > 0) 


# Number of species in total in for 28S and COI
COIvs28S <- ttest_species_across_units %>%
  group_by(primer)  %>%
  filter(!(duplicated(species))) %>% 
  dplyr::summarise(species_richness = n()) %>%
  as.data.frame()
rownames(COIvs28S) = COIvs28S$primer

# Number of species across units
COIvs28S_units <- ttest_species_across_units %>%
  group_by(Sample)  %>%
  filter(!(duplicated(species))) %>% 
  dplyr::summarise(species_richness = n()) %>%
  as.data.frame() %>%
  dplyr::mutate(primer = case_when(grepl("28", Sample) ~ "28",
                                   TRUE ~ "COI")) %>%
  dplyr::mutate(unit = case_when(grepl("17", Sample) ~ "17",
                                 grepl("18", Sample) ~ "18",
                                   TRUE ~ "20")) %>%
  select(-(Sample)) %>%
  pivot_wider(names_from = primer, values_from= species_richness)
  
# Shapiro and variance tests
shapiro.test (COIvs28S_units$`28`)
shapiro.test (COIvs28S_units$COI)

# variance not equal
var.test(COIvs28S_units$`28`, COIvs28S_units$COI)

# T.test
ttestrichnessacrossunits <- t.test(COIvs28S_units$`28`, COIvs28S_units$COI)
ttestrichnessacrossunits

