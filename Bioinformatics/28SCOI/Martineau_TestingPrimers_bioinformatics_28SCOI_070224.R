############################################
########  A. Load working directory ########
############################################

# setwd("change for your WD")



########################################################################
########  B. Load blast results and select only the ones needed ########
########################################################################

## WRITE IT FOR DIFFERENT FORMATS OF OUTPUTS
## B.1 Load blastn results, format in table and save to csv
load_blast <- function (blast_tab_file, db_name){
  
  # Install packages for that function
  if (!require("data.table")) install.packages("data.table")
  library("data.table")
  
  # Lod blastn file in the tab format
  blastn_file <- fread(blast_tab_file)
  
  # Rename the columns
  names(blastn_file) <- c("query.id", 
                          "identity",   
                          "subject.id",
                          "scientific.name",
                          "common.name",
                          "super.kingdom",
                          "percent.identity", 
                          "alignment.length",
                          "query.length",
                          "sequence.length",
                          "mismatches", 
                          "gap.opens",
                          "gaps",
                          "qstart", 
                          "sequence.end",
                          "start.of.alignment",
                          "qend",
                          "subject.title",
                          "evalue", 
                          "bitscore",
                          "query.coverage",
                          "query.coverage.HSP")
  # Save the raw output to csv in the same folder directory recursive directory: create it
  current_directory <- getwd()
  subdirectory_name <- "outputs"
  if (!file.exists(file.path(current_directory, subdirectory_name))) {
    # If it doesn't exist, create the subdirectory
    dir.create(file.path(current_directory, subdirectory_name))
  } else {
  }
  write.csv(blastn_file, paste0("./outputs/", db_name, "_raw_blastn_output.csv"), row.names = FALSE)
  # Save that table to global.env (rename)
  blastn_var_name <- paste0(db_name, "_blastn_file_raw")
  assign(blastn_var_name, blastn_file, envir = .GlobalEnv)
}

load_blast(blast_tab_file = "blast_results_spongesdb.tab", db_name = "sponge")



## B.2 Filter blastn outputs based on a percentage identity and query.id parameters, delete blastn with the same names and and save the others to a fasta file

filter_blast <- function (raw_blastn, min_perc_id, min_query_cov, sequences_abundance_path, abundance_path, db_name, return_fasta_and_remaining) {
  # Install necessary packages if needed
  necessarypackages = c("dplyr", "magrittr", "data.table","seqinr")
  # Install packages for that function
  for (package in necessarypackages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
    }
  }
  # Load necessary libraries
  lapply(necessarypackages, library, character.only = TRUE)
  
  # Condition based on return_fasta_remaining arguemnt)     
  if (return_fasta_and_remaining == TRUE){
    
    #Filter the blastn results based on query ID and % identity
    top_blastn <- raw_blastn %>%
      dplyr::filter(percent.identity  >= min_perc_id & query.coverage >= min_query_cov) %>%
      dplyr::group_by (query.id) %>%
      # Among them, keep the ones with top % identity
      filter(percent.identity == max (percent.identity, na.rm = TRUE)) %>%
      # Delete duplicates by scientific name
      dplyr::group_by (query.id) %>%
      filter(!(duplicated(scientific.name)))
    
    # Save object to environment
    assign(paste0(db_name, "_filtered_blast"), top_blastn ,envir = .GlobalEnv)
    
    # Import abundance 
    abundance_table <- as.data.frame(data.table::fread (abundance_path), header = TRUE) 
    assign(paste0(db_name, "_abundance_table"), abundance_table, envir = .GlobalEnv)
    
    # Find the ASVs that did not make the cut out of the total ASV from beginning
    all_ASV <- abundance_table[[1]]
    ASV_not_retained <- all_ASV[which(!(all_ASV %in% (unique(top_blastn$query.id))))]
    
    # Import sequences 
    sequence_table <- as.data.frame(data.table::fread (sequences_abundance_path), header = TRUE)
    sequences <- colnames(sequence_table)[2: ncol(sequence_table)]
    
    # Make a df that has ASV names and sequences
    sequences_and_ASVs <- as.data.frame(cbind(all_ASV, sequences)) 
    colnames(sequences_and_ASVs) <- c("all_ASV", "sequences")
    sequences_and_ASVs <- sequences_and_ASVs %>%
      dplyr::filter(!(all_ASV %in% (unique(top_blastn$query.id))))
    
    # Make output directory if not already made
    current_directory <- getwd()
    subdirectory_name <- "outputs"
    if (!file.exists(file.path(current_directory, subdirectory_name))) {
      # If it doesn't exist, create the subdirectory
      dir.create(file.path(current_directory, subdirectory_name))
    }
    
    # Save that df in the working directory and change name
    write.csv(sequences_and_ASVs, paste0("./outputs/", db_name, "_remaining_asvs.csv"), row.names = FALSE)
    
    # Return a fasta file with only unique 
    uniqueseqs <- as.list(sequences_and_ASVs$sequences)
    seqnum <- sequences_and_ASVs$all_ASV
    # Save fasta file
    write.fasta(uniqueseqs, seqnum, paste0("./outputs/", db_name, "_sequences_not_retained.fasta"))
    
    # Verify if the names of sequences in fasta files are the opposite of top blast results
    verify_test <- c(seqnum, unique(top_blastn$query.id))
    if (identical(sort(verify_test), sort(all_ASV))) {
      print("initial filtering successful")
    } else {
      print("initial filtering corrupted")
    }
    
    # Condition based on return_fasta_remaining arguemnt)          
  } else {
    
    #Filter the blastn results based on query ID and % identity
    top_blastn <- raw_blastn %>%
      dplyr::filter(percent.identity  >= min_perc_id & query.coverage >= min_query_cov) %>%
      dplyr::group_by (query.id) %>%
      # Among them, keep the ones with top % identity
      filter(percent.identity == max (percent.identity, na.rm = TRUE)) %>%
      # Delete duplicates by scientific name
      dplyr::group_by (query.id) %>%
      filter(!(duplicated(scientific.name)))
    
    # Save object to environment
    assign(paste0(db_name, "_filtered_blast"), top_blastn ,envir = .GlobalEnv)
    
    # Import abundance and assign to environment
    abundance_table <- as.data.frame(data.table::fread (abundance_path), header = TRUE) 
    assign(paste0(db_name, "_abundance_table"), abundance_table, envir = .GlobalEnv)
  }
}


# Run the filter_blastn function
filter_blast (raw_blastn = sponge_blastn_file_raw,
              min_perc_id = 95,
              min_query_cov = 95,
              sequences_abundance_path = "Martineau_Testing_primers_abundance_sequences_070224.rarefied.csv",
              abundance_path = "Martineau_Testing_primers_abundance_ASV_070224.rarefied.csv",
              db_name = "sponge",
              return_fasta_and_remaining = TRUE)



## B.3 Delete Duplicates function 
delete_duplicated_ID <- function (filtered_blastn, db_name){
  
  #Install necessary libraries for part C
  necessarypackages = c("dplyr", "magrittr")
  # Install packages for that function
  for (package in necessarypackages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
    }
  }
  # Delete duplicates
  no.duplicated.blast.searches <- filtered_blastn %>%
    group_by(query.id) %>%
    filter(!(duplicated(query.id)))
  
  # Verify if it worked
  if (length(unique(filtered_blastn$query.id)) == length(no.duplicated.blast.searches$query.id)) {
    print("remove duplicates successful")
  } else {
    print("remove duplicates failure")
  }
  # Save deduplicated df to environment
  assign(paste0(db_name, "_no_duplicates"), no.duplicated.blast.searches, envir = .GlobalEnv)
}


## B.5 Delete duplicated species
duplicated_species <- function (abundance, taxonomy, filename){
  # J.1 Delete duplicates
  abundance <- abundance
  taxonomy <- taxonomy
  
  taxonomy_colnames <- colnames(taxonomy)
  # Remove "ID" from the column names
  taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
  
  row.names(taxonomy) = taxonomy$ID
  
  ## Verify if there are duplicates
  nrow(taxonomy) == length(unique(row.names(taxonomy)))
  
  ## Join combined tax table and abundance_global tables together
  a.combined_table <- inner_join(abundance, taxonomy, by ="ID")
  
  ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
  a.final.abundance.table <-  a.combined_table %>%
    dplyr::select(colnames(abundance))
  row.names(a.final.abundance.table)<- a.final.abundance.table$ID
  
  unique_taxa <- a.combined_table %>%
    distinct(species)
  
  ## J.2 Delete duplicates in both the abundance and tax tables
  # Find the number of unique taxa
  rownames(a.combined_table) = a.combined_table$ID
  
  extensive_taxonomy <- taxonomy
  
  # when species are duplicated, pick the one that has the most extensive taxonomy
  # For each duplicated species, return the number of "unclassified" in each row, if not return NA 
  extensive_taxonomy$unclassified <- ifelse(duplicated(extensive_taxonomy$species) | duplicated(extensive_taxonomy$species, fromLast = TRUE),rowSums(extensive_taxonomy == "unclassified"), NA)
  
  extensive_taxonomy_verified <- subset(extensive_taxonomy, !(duplicated(species) & unclassified != min(unclassified)))
  extensive_taxonomy_verified$unclassified <- NULL
  
  # Verify if there are still duplicated species: should have zero columns in total
  duplicated_species_still <- extensive_taxonomy_verified[which(duplicated(extensive_taxonomy_verified$ID)),]
  
  # J.3 Map those taxa
  unique_taxa <- extensive_taxonomy_verified %>%
    distinct(species)
  unique_taxa$mapping_ASV <- rownames(unique_taxa)
  
  number_unique_taxa <- nrow(unique_taxa) 
  
  # Remove duplicates from taxonomy
  # Delete the ones that are not mapping taxa
  taxonomy_mapping <- extensive_taxonomy_verified %>%
    dplyr::filter(ID %in% unique_taxa$mapping_ASV)
  
  abundance_mapping <- a.combined_table %>%
    dplyr::select(c(colnames(abundance), species)) 
  
  abundance_mapping_join <- right_join(unique_taxa, abundance_mapping) %>%
    dplyr::select(-c("ID", "species"))
  
  # Merge abundance by mapping ASV
  abundance_no_duplicates_2 = aggregate (. ~ abundance_mapping_join$mapping_ASV, abundance_mapping_join[,!(colnames(abundance_mapping_join) %in% c("mapping_ASV"))], sum)
  colnames(abundance_no_duplicates_2)[1] <- "ID"
  rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$ID
  
  # Check if all the mapping taxa are present
  if (identical(sort(unique_taxa$mapping_ASV), sort(abundance_no_duplicates_2$ID)))  {
    print( "V3 good" )
  } else {
    print ("V3 is wrong")
  }
  
  # J.4 Saved finished tables
  # Save cleaned ASV table to environment
  assign(paste0(filename, "_taxonomy_table"),taxonomy_mapping, envir = .GlobalEnv)
  assign(paste0(filename, "_abundance_table"),abundance_no_duplicates_2, envir = .GlobalEnv)
  
  # Make output directory if not already made
  current_directory <- getwd()
  subdirectory_name <- "outputs"
  if (!file.exists(file.path(current_directory, subdirectory_name))) {
    # If it doesn't exist, create the subdirectory
    dir.create(file.path(current_directory, subdirectory_name))
  }
  
  # Save that df in the working directory
  write.csv(taxonomy_mapping, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
  write.csv(abundance_no_duplicates_2, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
}


## B.6 Merge tables function
merge_taxo_abundance <- function (taxonomyA, taxonomyB, abundanceA,abundanceB,abundance_global, filename){
  
  # Install packages for that function
  necessarypackages = c("dplyr", "magrittr", "tidyr", "taxize", "purrr", "data.table", "knitr")
  # Install packages for that function
  for (package in necessarypackages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
    }
  }
  
  if (missing(abundanceA)){
    
    taxonomyA <-  as.data.frame(fread (paste0("./outputs/", taxonomyA), header = TRUE))
    names(taxonomyA)[1] = "ID"
    taxonomyB <-  as.data.frame(fread (paste0("./outputs/", taxonomyB), header = TRUE))
    names(taxonomyB)[1] = "ID"
    abundance_global <-  as.data.frame(data.table::fread (abundance_global), header = TRUE)
    names(abundance_global)[1] <- "ID"
    
    # Get the column names of taxonomyA
    taxonomy_colnames <- colnames(taxonomyA)
    # Remove "ID" from the column names
    taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
    
    ## Join two taxonomy tables into total.taxonomy
    total.taxonomy <- rbind(taxonomyA, taxonomyB)
    row.names(total.taxonomy) = total.taxonomy$ID
    
    ## Verify if there are duplicates
    nrow(total.taxonomy) == length(unique(row.names(total.taxonomy)))
    
    ## Join combined tax table and abundance_global tables together
    a.combined_table <- inner_join(abundance_global, total.taxonomy, by ="ID")
    
    ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
    a.final.abundance.table <-  a.combined_table %>%
      dplyr::select(colnames(abundance_global))
    row.names(a.final.abundance.table)<- a.final.abundance.table$ID
    
    
    ## G.2 Delete duplicates in both the abundance and tax tables
    # Find the number of unique taxa
    rownames(a.combined_table) = a.combined_table$ID
    unique_taxa <- a.combined_table %>%
      distinct(across(all_of(taxonomy_colnames)))
    
    
    number_unique_taxa <- nrow(unique_taxa)
    # Find the unique taxa and their mapping ASV
    unique_taxa$mapping_ASV <- rownames(unique_taxa)
    mapping_taxa <- unique_taxa %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV)
    
    # Remove duplicates from total.taxonomy
    total.taxonomy$ID <- rownames(total.taxonomy)
    
    duplicated_tax_table <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n) 
    
    # Add mapping taxa and replace taxa name by mapping taxa ID in duplicated abundance table
    duplicated_tax_table_2 = right_join(mapping_taxa, duplicated_tax_table)
    rownames(duplicated_tax_table_2) <- duplicated_tax_table_2$mapping_ASV
    duplicated_tax_table_2 = as.data.frame(t(duplicated_tax_table_2[,!(colnames(duplicated_tax_table_2) %in% c("mapping_ASV", taxonomy_colnames))]))
    
    duplicated_test <- which(!is.na(duplicated_tax_table_2), arr.ind=TRUE)
    duplicated_tax_table_2[duplicated_test] <- rownames(duplicated_test)
    
    # Try to clean it up and remove NAs
    #define a function
    RemoveNAs <- function(x, size) {
      temp <- x[!is.na(x)]
      c(temp, rep("", size - length(temp)))
    }
    
    # Calculate longest column size in df (non-NA)
    Max <- max(colSums(!is.na(duplicated_tax_table_2)))
    
    # Remove NAs from duplicated_tax_table_2
    Indexes <- setDT(duplicated_tax_table_2)[, lapply(.SD, RemoveNAs, Max)]
    
    # Some verification:
    # Check if the most abundant unique species has the right number of columns
    max.count = duplicated_tax_table <- total.taxonomy%>%
      add_count(across(all_of(taxonomy_colnames)))
    max.count = max(max.count$n)
    
    if (max.count == Max & max.count == nrow(Indexes)){
      print ("V1 good")
    } else {
      print ("V1 is wrong")
    }
    
    
    # Check if the number of unique ASV is good
    if (ncol(Indexes) == number_unique_taxa){
      print( "V2 good" )
    } else {
      print ("V2 is wrong")
    }
    
    # Delete duplicated rows on tax table and abundace tables
    duplicated_tax_table_1 <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n)
    
    tax_table_no_duplicate = right_join (mapping_taxa, duplicated_tax_table_1) %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV) 
    rownames(tax_table_no_duplicate) <- tax_table_no_duplicate$mapping_ASV
    
    
    # Delete duplicates on abundance table
    # Add mapping ASV column to the original tax table
    index_taxonomy <- merge(total.taxonomy, mapping_taxa) 
    
    # Check if all the mapping taxa are present
    if (length(unique(mapping_taxa$mapping_ASV)) == length(unique(index_taxonomy$mapping_ASV)))  {
      print( "V3 good" )
    } else {
      print ("V3 is wrong")
    }
    
    abundance_no_duplicates <- inner_join(abundance_global, index_taxonomy)%>%
      dplyr::select(-c(all_of(taxonomy_colnames)))
    
    abundance_no_duplicates_1 <- abundance_no_duplicates %>%
      dplyr::select(-c(ID))
    
    abundance_no_duplicates_2 = aggregate (. ~ mapping_ASV, abundance_no_duplicates_1, sum)
    ## G.3 Verification steps: does the total equals to the sum of its part for different indexes?
    # Create data frames to store the results
    indexes_test_df <-  data.frame(matrix(nrow = nrow(Indexes), ncol = ncol(Indexes)))
    colnames(indexes_test_df) <- colnames(Indexes)
    indexes_test_df[] <- NA
    caculation_test_df <-  data.frame(matrix(nrow = ncol(abundance_no_duplicates_2) -1, ncol = ncol(Indexes)))
    colnames(caculation_test_df) <- colnames(Indexes)
    
    for (index in colnames(Indexes)){
      I <- Indexes[[index]]
      V <- abundance_no_duplicates %>%
        dplyr::filter(mapping_ASV %in% index)%>%
        as.data.frame()
      
      v.vector <-  V$ID %in% I
      replacement_length <- length(v.vector)
      indexes_test_df[, index] <- c(as.character(v.vector), rep(NA, nrow(indexes_test_df) - replacement_length))
      
      columns_to_exclude <- c("mapping_ASV", "ID")
      V.1 <- colSums(V[ , !(names(V) %in% columns_to_exclude)])
      caculation_test_df[,index] <- as.character(V.1 ==  abundance_no_duplicates_2[which(abundance_no_duplicates_2$mapping_ASV == index), !colnames(abundance_no_duplicates_2) %in% columns_to_exclude])
    }
    if (any(is.na(indexes_test_df)) || any(is.na(caculation_test_df))){
      print ("V4 is good")
    } else {
      print ("V4 wrong")
    }
    
    
    # Check if there are still duplicated species
    if (length(unique(abundance_no_duplicates$mapping_ASV)) == length(unique(index_taxonomy$species))) {
      print( "V5 good" )
      
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),test_taxonomy_table, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),test_abudance_no_duplicates_2, envir = .GlobalEnv)
      
      # Make output directory if not already made
      current_directory <- getwd()
      subdirectory_name <- "outputs"
      if (!file.exists(file.path(current_directory, subdirectory_name))) {
        # If it doesn't exist, create the subdirectory
        dir.create(file.path(current_directory, subdirectory_name))
      }
      
      # Save that df in the working directory
      write.csv(tax_table_no_duplicate, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
      write.csv(abundance_no_duplicates_2, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
      write.csv (Indexes, paste0("./outputs/", filename, "_Indexes.csv"), row.names = FALSE)
      
      
      
    } else {
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      duplicated_species(taxonomy = tax_table_no_duplicate, abundance = abundance_no_duplicates_2, filename = "test" )
      if   (length(unique(test_abundance_table$ID)) == length(unique(index_taxonomy$species))) {
        print( "V5 good" )
      } else {
        print ("V5 is wrong")
      }
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),test_taxonomy_table, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),test_abundance_table, envir = .GlobalEnv)
      
      # Make output directory if not already made
      current_directory <- getwd()
      subdirectory_name <- "outputs"
      if (!file.exists(file.path(current_directory, subdirectory_name))) {
        # If it doesn't exist, create the subdirectory
        dir.create(file.path(current_directory, subdirectory_name))
      }
      
      # Save that df in the working directory
      write.csv(test_taxonomy_table, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
      write.csv(test_abundance_table, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
      write.csv (Indexes, paste0("./outputs/", filename, "_Indexes.csv"), row.names = FALSE)
      
      
    }
    
    
    # Condition based on abundance A argument missing or not   
  } else {
    
    taxonomyA <-  as.data.frame(fread (paste0("./", taxonomyA), header = TRUE))
    names(taxonomyA)[1] = "ID"
    taxonomyB <-  as.data.frame(fread (paste0("./", taxonomyB), header = TRUE))
    names(taxonomyB)[1] = "ID"
    abundanceA <-  as.data.frame(data.table::fread (abundanceA), header = TRUE)
    names(abundanceA)[1] <- "ID"
    abundanceB <- as.data.frame(data.table::fread (abundanceB), header = TRUE)
    names(abundanceB)[1] <- "ID"
    
    # Change the names for abundance and taxonomy tables B (just to make sure they are not duplicated)
    combined_table_B <- merge(abundanceB, taxonomyB, by = "ID")
    
    # Rename the ASV
    rownames(combined_table_B)=   ASV.sequences = paste("ASV",  seq(from = 10000,to = (9999  + nrow(combined_table_B))), sep ="")
    combined_table_B$mapping = row.names(combined_table_B)
    
    # Save the mapping in another data frame
    Mapping_ASV <- as.data.frame(cbind(combined_table_B$ID, combined_table_B$mapping))
    
    # Resplit them
    combined_table_B$ID <- NULL
    combined_table_B$mapping <- NULL
    
    # Find out how many taxa are in taxonomy table
    taxonomy_colnames_B <- colnames(taxonomyB)
    # Remove "ID" from the column names
    taxonomy_colnames_B <- setdiff(taxonomy_colnames_B, "ID")
    
    taxonomy.B.rename <- combined_table_B[, colnames(combined_table_B) %in% taxonomy_colnames_B]
    taxonomy.B.rename$ID = row.names(taxonomy.B.rename)
    
    abundance.B.rename <- combined_table_B[, !(colnames(combined_table_B) %in% taxonomy_colnames_B)]
    abundance.B.rename$ID <- row.names(abundance.B.rename)
    
    # Get the common taxonomy colnames
    taxonomy.A.colnames <- colnames(taxonomyA)
    
    # Take the intersection of taxonomy A and B colnames
    taxonomy_colnames <- intersect(taxonomy.A.colnames, taxonomy_colnames_B)
    
    # Remove "ID" from the column names
    taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
    
    ## Join two taxonomy tables into total.taxonomy
    total.taxonomy <- rbind(taxonomyA, taxonomy.B.rename)
    row.names(total.taxonomy) = total.taxonomy$ID
    
    ## Verify if there are duplicates
    nrow(total.taxonomy) == length(unique(row.names(total.taxonomy)))
    
    # Combine abundance table
    abundance_combined <- full_join(abundanceA, abundance.B.rename)
    abundance_combined[is.na(abundance_combined)] <- 0
    
    
    ## Join combined tax table and  abundance_combined tables together
    a.combined_table <- inner_join( abundance_combined, total.taxonomy, by ="ID")
    
    ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
    a.final.abundance.table <-  a.combined_table %>%
      dplyr::select(colnames( abundance_combined))
    row.names(a.final.abundance.table)<- a.final.abundance.table$ID
    
    
    ## G.2 Delete duplicates in both the abundance and tax tables
    # Find the number of unique taxa
    rownames(a.combined_table) = a.combined_table$ID
    unique_taxa <- a.combined_table %>%
      distinct(across(all_of(taxonomy_colnames)))
    
    
    number_unique_taxa <- nrow(unique_taxa)
    # Find the unique taxa and their mapping ASV
    unique_taxa$mapping_ASV <- rownames(unique_taxa)
    mapping_taxa <- unique_taxa %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV)
    
    # Remove duplicates from total.taxonomy
    total.taxonomy$ID <- rownames(total.taxonomy)
    
    duplicated_tax_table <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n) 
    
    # Add mapping taxa and replace taxa name by mapping taxa ID in duplicated abundance table
    duplicated_tax_table_2 = right_join(mapping_taxa, duplicated_tax_table)
    rownames(duplicated_tax_table_2) <- duplicated_tax_table_2$mapping_ASV
    duplicated_tax_table_2 = as.data.frame(t(duplicated_tax_table_2[,!(colnames(duplicated_tax_table_2) %in% c("mapping_ASV", taxonomy_colnames))]))
    
    duplicated_test <- which(!is.na(duplicated_tax_table_2), arr.ind=TRUE)
    duplicated_tax_table_2[duplicated_test] <- rownames(duplicated_test)
    
    # Try to clean it up and remove NAs
    #define a function
    RemoveNAs <- function(x, size) {
      temp <- x[!is.na(x)]
      c(temp, rep("", size - length(temp)))
    }
    
    # Calculate longest column size in df (non-NA)
    Max <- max(colSums(!is.na(duplicated_tax_table_2)))
    
    # Remove NAs from duplicated_tax_table_2
    Indexes <- setDT(duplicated_tax_table_2)[, lapply(.SD, RemoveNAs, Max)]
    
    # Some verification:
    # Check if the most abundant unique species has the right number of columns
    max.count = duplicated_tax_table <- total.taxonomy%>%
      add_count(across(all_of(taxonomy_colnames)))
    max.count = max(max.count$n)
    
    if (max.count == Max & max.count == nrow(Indexes)){
      print ("V1 good")
    } else {
      print ("V1 is wrong")
    }
    
    
    # Check if the number of unique ASV is good
    if (ncol(Indexes) == number_unique_taxa){
      print( "V2 good" )
    } else {
      print ("V2 is wrong")
    }
    
    # Delete duplicated rows on tax table and abundace tables
    duplicated_tax_table_1 <- total.taxonomy %>%
      add_count(across(all_of(taxonomy_colnames))) %>%
      pivot_wider(names_from = ID, values_from = n)
    
    tax_table_no_duplicate = right_join (mapping_taxa, duplicated_tax_table_1) %>%
      dplyr::select(c(all_of(taxonomy_colnames)), mapping_ASV) 
    rownames(tax_table_no_duplicate) <- tax_table_no_duplicate$mapping_ASV
    
    # Delete duplicates on abundance table
    # Add mapping ASV column to the original tax table
    index_taxonomy <- merge(total.taxonomy, mapping_taxa) 
    
    # Check if all the mapping taxa are present
    if (length(unique(mapping_taxa$mapping_ASV)) == length(unique(index_taxonomy$mapping_ASV)))  {
      print( "V3 good" )
    } else {
      print ("V3 is wrong")
    }
    
    abundance_no_duplicates <- inner_join( abundance_combined, index_taxonomy)%>%
      dplyr::select(-c(all_of(taxonomy_colnames)))
    
    abundance_no_duplicates_1 <- abundance_no_duplicates %>%
      dplyr::select(-c(ID))
    
    abundance_no_duplicates_2 = aggregate (. ~ mapping_ASV, abundance_no_duplicates_1, sum)
    
    ## G.3 Verification steps: does the total equals to the sum of its part for different indexes?
    # Create data frames to store the results
    indexes_test_df <-  data.frame(matrix(nrow = nrow(Indexes), ncol = ncol(Indexes)))
    colnames(indexes_test_df) <- colnames(Indexes)
    indexes_test_df[] <- NA
    caculation_test_df <-  data.frame(matrix(nrow = ncol(abundance_no_duplicates_2) -1, ncol = ncol(Indexes)))
    colnames(caculation_test_df) <- colnames(Indexes)
    
    for (index in colnames(Indexes)){
      I <- Indexes[[index]]
      V <- abundance_no_duplicates %>%
        dplyr::filter(mapping_ASV %in% index)%>%
        as.data.frame()
      
      v.vector <-  V$ID %in% I
      replacement_length <- length(v.vector)
      indexes_test_df[, index] <- c(as.character(v.vector), rep(NA, nrow(indexes_test_df) - replacement_length))
      
      columns_to_exclude <- c("mapping_ASV", "ID")
      V.1 <- colSums(V[ , !(names(V) %in% columns_to_exclude)])
      caculation_test_df[,index] <- as.character(V.1 ==  abundance_no_duplicates_2[which(abundance_no_duplicates_2$mapping_ASV == index), !colnames(abundance_no_duplicates_2) %in% columns_to_exclude])
    }
    if (any(is.na(indexes_test_df)) || any(is.na(caculation_test_df))){
      print ("V4 is good")
    } else {
      print ("V4 wrong")
    }
    
    # Check if there are still duplicated species
    if (length(unique(abundance_no_duplicates$mapping_ASV)) == length(unique(index_taxonomy$species))) {
      print( "V5 good" )
      
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),tax_table_no_duplicate, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),abundance_no_duplicates_2, envir = .GlobalEnv)
      
      # Make output directory if not already made
      current_directory <- getwd()
      subdirectory_name <- "outputs"
      if (!file.exists(file.path(current_directory, subdirectory_name))) {
        # If it doesn't exist, create the subdirectory
        dir.create(file.path(current_directory, subdirectory_name))
      }
      
      # Save that df in the working directory
      write.csv(tax_table_no_duplicate, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
      write.csv(abundance_no_duplicates_2, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
      write.csv (Indexes, paste0("./outputs/", filename, "_Indexes.csv"), row.names = FALSE)
      
      
      
    } else {
      rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$mapping_ASV
      mapping_ASV <- "mapping_ASV"
      abundance_no_duplicates_2 <- abundance_no_duplicates_2[, c(mapping_ASV, setdiff(names(abundance_no_duplicates_2), mapping_ASV))]
      colnames(abundance_no_duplicates_2)[1] = "ID"
      
      rownames(tax_table_no_duplicate) = tax_table_no_duplicate$mapping_ASV
      tax_table_no_duplicate <- tax_table_no_duplicate[, c(mapping_ASV, setdiff(names(tax_table_no_duplicate), mapping_ASV))]
      colnames(tax_table_no_duplicate)[1] = "ID"
      
      duplicated_species(taxonomy = tax_table_no_duplicate, abundance = abundance_no_duplicates_2, filename = "test" )
      if   (length(unique(test_abundance_table$ID)) == length(unique(index_taxonomy$species))) {
        print( "V5 good" )
      } else {
        print ("V5 is wrong")
      }
      # Save cleaned ASV table to environment
      assign(paste0(filename, "_taxonomy_table"),test_taxonomy_table, envir = .GlobalEnv)
      assign(paste0(filename, "_abundance_table"),test_abundance_table, envir = .GlobalEnv)
      
      # Make output directory if not already made
      current_directory <- getwd()
      subdirectory_name <- "outputs"
      if (!file.exists(file.path(current_directory, subdirectory_name))) {
        # If it doesn't exist, create the subdirectory
        dir.create(file.path(current_directory, subdirectory_name))
      }
      
      # Save that df in the working directory
      write.csv(test_taxonomy_table, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
      write.csv(test_abundance_table, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
      write.csv (Indexes, paste0("./outputs/", filename, "_Indexes.csv"), row.names = FALSE)
      
      
    }
    
  }
}


###########################################################
########  C. Delete duplicates (manual filtration) ########
###########################################################

## C.1 Find duplicates and data wrangling
#Install necessary libraries for part C
necessarypackages = c("dplyr", "magrittr", "tidyr")
# Install packages for that function
for (package in necessarypackages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

## C.2 Just delete all duplicates 
# Identify duplicated ASVs
dupes <- sponge_filtered_blast %>%
  filter(duplicated(query.id))

duplicates.vector = unique(dupes$query.id)
length(duplicates.vector)

# Filter out non-duplicated ASVs
non_dupes <- sponge_filtered_blast %>%
  filter(!(query.id %in% unique(dupes$query.id)))

# Filter duplicated ASVs
dupes_filtration <- sponge_filtered_blast %>%
  filter(query.id %in% unique(dupes$query.id)) %>%
  select(ID_ASV = query.id, species = scientific.name, subject = subject.id)

# Merge abundance table and species
dupes_filtration_merged <- right_join(sponge_abundance_table %>%
                                        mutate(ID_ASV = V1) %>%
                                        select(-V1), dupes_filtration)

# Convert to long format
dupes_filtration_long <- dupes_filtration_merged %>%
  pivot_longer(cols = ARMS20:negative,
               names_to = "Sample",
               values_to = "values") %>%
  mutate(values = as.numeric(values))

##  EQ20
# Find duplicates that are present in Eq20 cocktail
species.in.cocktail <- c("Leucosolenida sp. 11 JV-2020",
                         "Leucosolenida sp. 13 JV-2020",
                         "Clathrinidae sp. 3 JV-2020",
                         "Clathrinidae sp. 4 JV-2020",
                         "Clathrinidae sp. 5 JV-2020",
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


# Among replicates, find the ones that are in cocktail eq20
duplicates_eq20 <- dupes_filtration_long %>%
  filter(Sample == "Equimolar20" & values > 0 & species %in% species.in.cocktail)

# Eq18

# Select the duplicates that are left (not in eq20) and remove them from the long data set
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_eq20$ID_ASV)), ]

# Species in cocktail 18
# Find false positives fro Equimolar 20
species.in.cocktail.2   <- c("Clathrinidae sp. 3 JV-2020",
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

# from those, inspect the ones present in eq18
duplicates_eq18 <- dupes_filtration_long %>%
  filter(Sample == "Equimolar18" & values > 0 & species %in% species.in.cocktail.2)


# PC 20
# Select the duplicates that are left (not in eq20, not in eq18) and remove them from the long data set
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_eq18$ID_ASV)), ]

# Among duplicates left, find the ones that are in cocktail pc20
duplicates_pc20 <- dupes_filtration_long %>%
  filter(Sample == "%cover20" & values > 0 & species %in% species.in.cocktail)


# PC 18
# Select the duplicates that are left (not in eq20, not in eq18, not in PC20) and remove them from the long data set
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_pc20$ID_ASV)), ]

# Among duplicates left, find the ones that are in cocktail pc18
duplicates_pc18 <- dupes_filtration_long %>%
  filter(Sample == "%cover18" & values > 0 & species %in% species.in.cocktail.2)

# ARMS 20
# Select the duplicates that are left (not in eq20, not in eq18, not in PC20, not PC18) and remove them from the long data set
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_pc18$ID_ASV)), ]

# Among duplicates left, find the ones that are in cocktail ARMS20
duplicates_ARMS20 <- dupes_filtration_long %>%
  filter(Sample == "ARMS20" & values > 0 & species %in% species.in.cocktail)

# ARMS 18
# Select the duplicates that are left (not in eq20, not in eq18, not in PC20, not PC18, not ARMS20) and remove them from the long data set
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_ARMS20$ID_ASV)), ]

# Among duplicates left, find the ones that are in cocktail ARMS18
duplicates_ARMS18 <- dupes_filtration_long %>%
  filter(Sample == "ARMS18" & values > 0 & species %in% species.in.cocktail.2)

# Inspect what is left: not in eq, not in pc and not in ARMS
dupes_filtration_long <- dupes_filtration_long[which(!(dupes_filtration_long$ID_ASV %in% duplicates_ARMS18$ID_ASV)), ]
length(unique(dupes_filtration_long$ID_ASV))


# Make the replacements 
## Selection of correct duplicates

# Make a data set to select off
dupes_filtration_long_clean <-  dupes_filtration_merged %>%
  pivot_longer(cols = ARMS20:negative,
               names_to = "Sample",
               values_to= "values") %>%
  dplyr::filter (values > 0)

dupes_filtration_long_clean$values <- as.numeric(dupes_filtration_long_clean$values)

# Select the duplicates identified in the cocktails 
# Label the inspected ones as "inspected"
clean_duplicated_blast_search <- dupes_filtration_long_clean %>%
  dplyr::mutate(species_new = case_when ((ID_ASV %in% duplicates_eq20$ID_ASV)& (Sample == "Equimolar20") & (species %in% species.in.cocktail) ~ species,
                                         (ID_ASV %in% duplicates_eq18$ID_ASV) & (Sample == "Equimolar18")& (species %in% species.in.cocktail.2) ~ species,
                                         (ID_ASV %in% duplicates_pc20$ID_ASV) & (Sample == "%cover20") &    (species %in% species.in.cocktail) ~ species,
                                         (ID_ASV %in% duplicates_pc18$ID_ASV) & (Sample == "%cover18")  &  (species %in% species.in.cocktail.2) ~ species,
                                         (ID_ASV %in% duplicates_ARMS20$ID_ASV) & (Sample == "ARMS20") &  (species %in% species.in.cocktail) ~ species,
                                         (ID_ASV %in% duplicates_ARMS18$ID_ASV) & (Sample == "ARMS18") & (species %in% species.in.cocktail.2) ~ species,
                                         (ID_ASV %in% dupes_filtration_long$ID_ASV) ~ "inspected", 
                                         TRUE ~ "NO" )) %>%
  # Deletes the duplicates that say "NO"
  dplyr::filter (!(species_new == "NO")) %>%
  # Delete the duplicated inspected
  dplyr::filter(!(duplicated(ID_ASV) & species_new == "inspected"))


## C.3 Verification
# Verify if there are still duplicates
length(unique(clean_duplicated_blast_search$ID_ASV)) == nrow(clean_duplicated_blast_search)

# Inspect the rest of the duplicates
duplicates_left <- clean_duplicated_blast_search[which(duplicated(clean_duplicated_blast_search$ID_ASV) == TRUE),]

# Those are still left because both species left are present in the cocktail (Leucosolenida sp. 13 and Leucosolenida sp. 7). We just have to pick one or the other
clean_duplicated_blast_search_2  <- clean_duplicated_blast_search %>%
  filter(!(duplicated(ID_ASV)))

length(unique(clean_duplicated_blast_search_2$ID_ASV)) == nrow(clean_duplicated_blast_search_2)

# Verify if all the duplicates are there
missing_duplicates <- duplicates.vector[which(!(duplicates.vector %in% clean_duplicated_blast_search_2$ID_ASV))] #should be empty
length(unique(duplicates.vector)) ==  length(unique(clean_duplicated_blast_search_2$ID_ASV))

# Replace at the end
duplicates_keep <- clean_duplicated_blast_search_2[, 1:2]
names(duplicates_keep) <- c("query.id","scientific.name")
duplicates_keep_df <- left_join (duplicates_keep, sponge_filtered_blast)

# Merge the corrected duplicates and the non-duplicates
no.duplicated.blast.searches_1 <- rbind (duplicates_keep_df, non_dupes)

# Last verification if same number of unique ASVs
length(unique(sponge_filtered_blast$query.id)) == length(unique(no.duplicated.blast.searches_1$query.id))



############################################
########  D. Create taxonomic table ########
############################################

# REQUIRES a blast search without duplicates that has the ASV/OTU name column labelled as query.id and tax.id as subject.id
# Change to if not specified, taxonomy vector is X

## D.1 Classification function:
classify_blast <- function (blastn_result, taxonomy_vector, db_name) {
  #Install necessary libraries for part D
  necessarypackages = c("dplyr", "magrittr", "tidyr", "taxize", "purrr", "data.table", "knitr")
  # Install packages for that function
  for (package in necessarypackages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
    }
  }
  
  # Change column name 'query.id' to 'ASV'
  blastn_result <-blastn_result %>%
    dplyr::rename(ASV = query.id)
  
  # Separate subject.id and remove junk columns (sometimes the taxonomic ID is duplicated in the column)
  NAs.IDnumber <-blastn_result %>%
    separate(subject.id, into = c("subject.id", "subject.id.junk"), sep = ';') %>%
    dplyr::rename("subject.id" = "subject.id") %>%
    dplyr::select(-subject.id.junk)
  
  # Extract ID numbers
  all_ID_numbers <- as.numeric(NAs.IDnumber$subject.id)
  
  # Split in intervals of 500 because sometimes the classify functions lags and store them into a list
  split_into_intervals <- function(char_vector, interval_size) {
    num_intervals <- ceiling(length(char_vector) / interval_size)
    interval_list <- split(char_vector, ceiling(seq_along(char_vector) / interval_size))
    return(interval_list)
  }
  
  all_ID_numbers_list <- split_into_intervals(all_ID_numbers,500)
  
  # Write classification function
  classify_vector <- function(id_vector) {
    classification(sci_id = id_vector, callopts = list(), db = "ncbi", return_id = TRUE)
  }
  
  classification_test <- lapply(all_ID_numbers_list, classify_vector)
  classified.asv <-  do.call("c", classification_test)
  
  # Quick check to verify if the right numbers of IDs have been researched
  taxonomic_assignment = classified.asv 
  length(taxonomic_assignment) == length(all_ID_numbers)
  
  #Store ASV names into vector (add name for ASV that are non-unique in case the data frame still has duplicates)
  NAs.IDnumber$ASV.unique <- make.unique (NAs.IDnumber$ASV)
  ASVs_taxonomy <- NAs.IDnumber$ASV.unique
  
  # If there are no duplicates, ASV and ASV_unique should be the same
  NAs.IDnumber$ASV.unique == NAs.IDnumber$ASV
  
  # Extract each database for each ASV and store them in a list
  df_list <- lapply(taxonomic_assignment, \(x) as.data.frame(t(x))) |>
    setNames(paste0('database', ASVs_taxonomy))
  
  # Initialize a list to store the processed dataframes
  df_list_2 <- vector("list", length(all_ID_numbers))
  
  # Perform data wrangling for each table to clean up the data frame
  for (i in seq_along(df_list)) {
    # Extract the dataframe and transpose
    df <- df_list[[i]]
    # Set colnames as rank
    colnames(df) = c(df[which(rownames(df) == "rank"), ])
    # Remove the rank row
    df <- df[-which(rownames(df) == "rank"), ]
    # Remove unnecessary columns
    drop <- c("no rank", "clade")
    df <- df[, !(names(df) %in% drop)]
    # Bind the ASV name
    df <- cbind(df, ASV = ASVs_taxonomy[i])
    # Rename the last column to "taxonomic.ID"
    names(df)[length(names(df))] <- "taxonomic.ID"
    # Remove the column that says ID
    df <- df[-which(rownames(df) == "id"), ]
    # Store the dataframe in the list with a unique name
    df_list_2[[i]] <- df
  }
  
  # Merge taxonomic information for each ASV in the taxonomy matrix
  amerged.taxo.df.2 <- rbindlist(df_list_2, fill = TRUE) 
  
  # Ensure that amerged.taxo.df.2 has one row for each element of df_list_2
  nrow(amerged.taxo.df.2) == length(df_list_2) 
  
  
  # Get the common columns between what I wanna keep and what is found in the data frame
  common_cols <- intersect(colnames(amerged.taxo.df.2), taxonomy_vector)
  
  # Keep only columns in vector taxa
  amerged.taxo.df.3 <- amerged.taxo.df.2 %>%
    dplyr::select(c("taxonomic.ID", common_cols))
  
  # Change the NA as unclassified
  amerged.taxo.df.3[is.na(amerged.taxo.df.3)] <- "unclassified"
  amerged.taxo.df.3 = as.data.frame(amerged.taxo.df.3)
  
  # Change the ASV names in case there were some duplicates
  amerged.taxo.df.4 <- amerged.taxo.df.3 %>%
    dplyr::rename(ASV.unique = taxonomic.ID)
  
  amerged.taxo.df.4 <- left_join(amerged.taxo.df.4, NAs.IDnumber, by = "ASV.unique")
  
  amerged.taxo.df.5 <- amerged.taxo.df.4 %>%
    dplyr::select(c(common_cols, ID = ASV.unique)) 
  
  #Final tax table
  rownames(amerged.taxo.df.5) <- amerged.taxo.df.5$ID
  
  # Verify if nothing in the order has been messed up
  input_classification <- as.data.frame(blastn_result) %>%
    dplyr::select(ID = ASV, species = scientific.name)%>%
    separate(species, into = c("species1", "species2"), sep = ' ') %>%
    dplyr::select(c(ID, species1, species2))
  
  output_classification <- as.data.frame(amerged.taxo.df.5)%>%
    dplyr::select(ID, species) %>%
    separate(species, into = c("species1", "species2"), sep = ' ') %>%
    dplyr::select(c(ID, species1, species2))
  
  # Check if columns in both data frames are identical
  all_successful <- TRUE
  for (col in colnames(input_classification)) {
    if (!identical(input_classification[[col]], output_classification[[col]])) {
      all_successful <- FALSE
      print(paste("Column", col, "comparison: corrupted"))
    }
  }
  
  if(all_successful) {
    print("Classification successful")
  }
  # Save cleaned ASV table to environment
  assign(paste0(db_name, "_taxonomy_table"),amerged.taxo.df.5, envir = .GlobalEnv)
  
  amerged.taxo.df.5 <- amerged.taxo.df.5 %>%
    select(ID, everything())
  # Make output directory if not already made
  current_directory <- getwd()
  subdirectory_name <- "outputs"
  if (!file.exists(file.path(current_directory, subdirectory_name))) {
    # If it doesn't exist, create the subdirectory
    dir.create(file.path(current_directory, subdirectory_name))
  } else {
  }
  
  # Save that df in the working directory
  write.csv(amerged.taxo.df.5, paste0("./outputs/", db_name, "_taxonomy_table.csv"), row.names = FALSE)
  
}

classify_blast(blastn_result = no.duplicated.blast.searches_1, taxonomy_vector = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), db_name = "sponge")



############################################
########  E. Taxonomic Replacements ########
############################################

## E.1 Data wrangling 
# Add taxonomic information to the abundance table
ASV_fp <- sponge_abundance_table %>%
  dplyr::rename(ID = V1)
rownames(ASV_fp) <- ASV_fp$ID

ASV_fp_long <- right_join(ASV_fp, sponge_taxonomy_table)%>%
  dplyr::select(c(species, colnames(ASV_fp))) %>%
  pivot_longer(cols = ARMS20:negative,
               names_to = "Sample",
               values_to= "values") %>%
  dplyr::mutate(values = as.numeric(values))

## E.2 Function: select but the expected species in each cocktail
# Select only the one cocktail
species_cocktail <- function (input_sample, species_list) {
  cocktails_df <- ASV_fp_long%>%
    dplyr::filter(Sample == input_sample) %>%
    dplyr::filter(values > 0) %>%
    dplyr::filter(!(species %in% species_list)) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(sum = sum(values))
  return(cocktails_df)
}

# Store the results in a data frame (Example for Equimolar 20) 
cocktails_df <- species_cocktail(input_sample =  "Equimolar20", species_list = species.in.cocktail)


# Do that for all samples and store the data frames in a list for each Sample
species_cocktails_list <- list()

# Make a list of cocktails and of species list
species.in.cocktail.3 <- c("Leucosolenida sp. 11 JV-2020",
                           "Leucosolenida sp. 13 JV-2020",
                           "Clathrinidae sp. 3 JV-2020",
                           "Clathrinidae sp. 4 JV-2020",
                           "Clathrinidae sp. 5 JV-2020",
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
species.in.all.otus <- c("Leucosolenida sp. 11 JV-2020",
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




cocktails <- c("Equimolar20", "Equimolar18", "%cover20", "%cover18", "ARMS20", "ARMS18")
cocktails_expected_species <- list(species.in.cocktail, species.in.cocktail.2, species.in.cocktail, species.in.cocktail.2, species.in.cocktail.3, species.in.cocktail.2)

# Loop over cocktails
for (i in seq_along(cocktails)) {
  sample_name <- cocktails[i]
  # Name each data frame based on the cocktail names 
  species_cocktails_list[[paste("species_cocktail", sample_name, sep = "_")]] <- 
    species_cocktail(sample_name, cocktails_expected_species[[i]])
}

# Check if Equimolar 20 df is the same as expected
nrow(species_cocktails_list[["species_cocktail_Equimolar20"]]) == nrow(cocktails_df)


## E.3 Manual swaps for each cocktail
# Define a function for manual swaps
manual_swaps <- function(data, cocktail_name, swaps) {
  cocktail_species <- species_cocktails_list[[paste("species_cocktail", cocktail_name, sep = "_")]]
  
  for (swap in swaps) {
    data <- data %>%
      mutate(species = case_when(
        ID %in% cocktail_species$ID & species == swap$from ~ swap$to,
        TRUE ~ species
      ))
  }
  
  return(data)
}
# Define swaps for each cocktail
swaps_eq20 <- list(
  list(from = "Suberitida sp. 4 JV-2020", to = "Suberitidae sp. 1 JV-2020"),
  list(from = "Suberitidae sp. 4 JV-2020", to = "Suberitidae sp. 1 JV-2020"),
  list(from = "Haplosclerida sp. 2 JV-2020", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Tethya sp. 5 JV-2020", to = "Tethya sp. 3 JV-2020"),
  list(from = "Haliclona sp. 23 JV-2020", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Haliclona laubenfelsi", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Haplosclerida sp. 2 JV-2020", to = "Haliclona sp. 1 JV-2020")
)
swaps_eq18 <- list(
  list(from = "Clathrinidae sp. 2 JV-2020", to = "Clathrinidae sp. 3 JV-2020"),
  list(from = "Haplosclerida sp. 15 JV-2020", to = "Haliclona sp. 2 JV-2020"),
  list(from = "Tethya sp. 4 JV-2020", to = "Tethya sp. 5 JV-2020")   
)
swaps_ARMS20 <- list(
  list(from = "Haplosclerida sp. 10 JV-2020", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Haplosclerida sp. 3 JV-2020", to = "Haliclona sp. 1 JV-2020"),
  list(from = "Tethyidae sp. 1 JV-2020", to = "Tethya sp. 3 JV-2020"),
  list(from = "Clathrinidae sp. 1 JV-2020", to = "Clathrinidae sp. 4 JV-2020"),
  list(from = "Leucosolenida sp. 10 JV-2020", to = "Leucosolenida sp. 7 JV-2020"),
  list(from = "Tethya sp. 1 JV-2020", to = "Tethya sp. 3 JV-2020")
)
swaps_pc20 <- list(
  list(from = "Tethya sp. 2 JV-2020", to = "Tethya sp. 3 JV-2020")
)


# Apply manual swaps for each cocktail
ASV_fp_long_modified <- manual_swaps(ASV_fp_long, "Equimolar20", swaps_eq20)
ASV_fp_long_modified <- manual_swaps(ASV_fp_long_modified, "Equimolar18", swaps_eq18)
ASV_fp_long_modified <- manual_swaps(ASV_fp_long_modified, "%cover20", swaps_pc20)
ASV_fp_long_modified <- manual_swaps(ASV_fp_long_modified, "ARMS20", swaps_ARMS20)


## E.4 Verifications
# Make same function but for verificaton
species_cocktail_verify <- function (input_sample, species_list) {
  cocktails_df_verify <- ASV_fp_long_modified%>%
    dplyr::filter(Sample == input_sample) %>%
    dplyr::filter(values > 0) %>%
    dplyr::filter(!(species %in% species_list)) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(sum = sum(values))
  return(cocktails_df_verify)
}
# Try a test for Equimolar 20
cocktails_df_verify  <- species_cocktail_verify(input_sample =  "Equimolar20", species_list = species.in.cocktail)

# Do that for all samples and store the data frames in a list for each Sample
species_cocktails_list_verify <- list()

# Loop over cocktails
for (i in seq_along(cocktails)) {
  sample_name <- cocktails[i]
  # Name each data frame based on the cocktail names 
  species_cocktails_list_verify[[paste("verify_species_cocktail", sample_name, sep = "_")]] <- 
    species_cocktail_verify(sample_name, cocktails_expected_species[[i]])
}

# Check if Equimolar 20 df is the same as expected
nrow(species_cocktails_list_verify[["verify_species_cocktail_Equimolar20"]]) == nrow(cocktails_df_verify)


## E.5 Recreate abundance and taxonomy tables
# Make the table wider again
new_abundance <- ASV_fp_long_modified %>%
  select(species, Sample, values, ID) %>%  
  pivot_wider(values_from = values,
              names_from = Sample) %>%
  as.data.frame()

nrow(new_abundance)== nrow(sponge_taxonomy_table)

# Extract taxonomy for all ASV: include everything but species
taxo_join <- sponge_taxonomy_table%>%
  dplyr::select(c(kingdom:genus, ID))%>%
  dplyr::filter(!(duplicated(.)))
rownames(taxo_join) <- NULL

# Assign the proper taxonomy to the new abundance table: join and then merge
new_taxo <- left_join (new_abundance, taxo_join, by =  "ID") %>%
  dplyr::select(c(ID,kingdom,phylum,class,order,family,genus,species))
rownames(new_taxo) <- new_taxo$ID
new_taxo$ID <- NULL


# Clean the abundance table
new_abundance$species <- NULL
rownames(new_abundance) <- new_abundance$ID
new_abundance$ID <- NULL 


# E.6 Save tables and working environment
# Save that df in the working directory
write.csv(new_taxo, "./outputs/sponge_taxonomy_table_swaps.csv", row.names = TRUE)
#Final cleaned table
write.csv(new_abundance,"./outputs/sponge_abundance_table_swaps.csv", row.names = TRUE)



############################################################################
########  F. Load, de-duplicate, filter and classify nt db searches ########
############################################################################

## F.1 Load blastn results
load_blast(blast_tab_file = "blast_results_ntdb.tab", db_name = "nt")

## F.2 Filter nt blast results
filter_blast (raw_blastn = nt_blastn_file_raw,
              min_perc_id = 95,
              min_query_cov = 95,
              sequences_abundance_path = "Martineau_Testing_primers_abundance_sequences_070224.rarefied.csv",
              abundance_path = "Martineau_Testing_primers_abundance_ASV_070224.rarefied.csv",
              db_name = "nt",
              return_fasta_and_remaining = FALSE)

## F.3 Delete duplicates
delete_duplicated_ID(filtered_blastn = nt_filtered_blast,
                     db_name = "nt")

## F.4 Classify the results
classify_blast(blastn_result = nt_no_duplicates, taxonomy_vector = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), db_name = "nt")



#######################################
########  G. Merge 2 data sets ########
#######################################

## G.1 Merge data sets
### THE data sets need unique ASV numbers... MAKE AN OPTION IF NON UNIQUE 
merge_taxo_abundance(taxonomyA = "sponge_taxonomy_table_swaps.csv",
                     taxonomyB = "nt_taxonomy_table.csv",
                     abundance_global = "Martineau_Testing_primers_abundance_ASV_070224.rarefied.csv",
                     filename = "sponge_nt")



################################################
########  H. Merge ARMS 18/20 to ARMS 17 #######
################################################

merge_taxo_abundance(taxonomyA = "outputs/sponge_nt_taxonomy_table.csv",
                     taxonomyB = "sponge_ntARMS17_taxonomy_table.csv",
                       #original path was: "~/U_of_H/SPONGES/Testing_primers/Data_analysis/Bioinformatics_040523/ARM17_clean_052723/ARM17_clean_Taxonomy_and_clean_tables/070224/outputs/sponge_ntARMS17_taxonomy_table.csv",
                     abundanceA = "outputs/sponge_nt_abundance_table.csv",
                     abundanceB = "sponge_ntARMS17_abundance_table.csv",
                       #original path was: "~/U_of_H/SPONGES/Testing_primers/Data_analysis/Bioinformatics_040523/ARM17_clean_052723/ARM17_clean_Taxonomy_and_clean_tables/070224/outputs/sponge_ntARMS17_abundance_table.csv",
                     filename = "combined171820")



################################################################################
########  I. Add expected values for equimolar and % coverage cocktails ########
################################################################################

## I.1 Add columns to abundance table
# Fill them with zeroes
fill_abundance_with_zeros <- function(abundance_table) {
  abundance_table[, paste0("ExpectedEq", c(17, 18, 20))] <- 0
  abundance_table
}
abundance_expected <- fill_abundance_with_zeros(combined171820_abundance_table)


## I.2 Verify if all species are already present in tjhe taxonomy table: if not add it and its taxonomy (manually)
species_17 <- c("Clathrinidae sp. 4 JV-2020",
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

species_18 <- c("Clathrinidae sp. 3 JV-2020",
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
species_20 <-  c("Leucosolenida sp. 11 JV-2020",
                 "Leucosolenida sp. 13 JV-2020",
                 "Clathrinidae sp. 3 JV-2020",
                 "Clathrinidae sp. 4 JV-2020",
                 "Clathrinidae sp. 5 JV-2020",
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
all_species <- c(species_17, species_18, species_20)
taxonomy_expected <- combined171820_taxonomy_table

# Return the species that are not present in the taxonomic table
species_not_present <- all_species[which(!all_species %in% taxonomy_expected$species)]
print(species_not_present)

# Add the taxonomy of it to the table
# Add row
add_missing_species <- function(taxonomy_table, missing_species, taxonomy_values) {
  rows <- nrow(taxonomy_table)
  taxonomy_table[rows + length(missing_species), ]  <- NA
  taxonomy_table[(rows + 1):nrow(taxonomy_table), c("ID", "species")] <- c(paste0("ASV", 8000), missing_species)
  if (length(taxonomy_values) == 6) {
    taxonomy_table[which(taxonomy_table$species %in% missing_species), c("kingdom", "phylum", "class", "order", "family", "genus")] <- taxonomy_values
  } else {
    warning("Taxonomy values provided do not match the expected length. Using default values.")
  }
  rownames(taxonomy_table) <- taxonomy_table$ID
  taxonomy_table
}
custom_taxonomy_values <- c("Metazoa", "Porifera", "Calcarea", "Leucosolenida", "unclassified", "unclassified")
taxonomy_expected <- add_missing_species(taxonomy_expected, species_not_present, custom_taxonomy_values)


## I.3 Map species to their ASV ID and add proper expected abundance to the table
# Merge abundance and taxonomy tables
abundance_expected_merge <- full_join(taxonomy_expected, abundance_expected) %>%
  dplyr::select(c(colnames(abundance_expected), species)) %>%
  as.data.frame()
abundance_expected_merge[is.na(abundance_expected_merge)]  <- 0

# Add values for equimolar
abundance_expected_equimolar <- abundance_expected_merge %>%
  mutate(ExpectedEq20 = ifelse(species %in% species_20, sum(Equimolar20) / length(species_20), 0), 
         ExpectedEq17 = ifelse(species %in% species_17, sum(Equimolar17R1) / length(species_17), 0),
         ExpectedEq18 = ifelse(species %in% species_18, sum(Equimolar18) / length(species_18), 0))

# percent cover, make a vector of the surface area in the same order as the species_20
surface_area_20 <- c(0.04090311, 0.030987205, 0.00600675, 0.90562092, 0.045352555, 
                     1.527971117, 0.017511743, 0.142922523, 0.178518081, 0.964448936, 
                     0.51715261, 0.055363806, 0.215607381, 0.408745066, 0.71238789, 
                     0.179121934, 0.314988908, 1.912530272)
surface_area_18 <- c(0.024916891, 0.237091843, 0.052662357, 0.017511743, 0.041951908, 
                     0.003051048, 0.663253307, 0.009439179, 0.460295062, 0.078405573, 
                     0.055363806, 0.023486712, 0.772900339, 0.109710595, 0.32738379, 
                     0.117083961)
surface_area_17 <- c(1.621123421, 0.018751232, 0.018433414, 0.061179866, 0.000921671,
                     0.011123612, 0.033370836, 0.0454479, 1.697685653, 0.00381381,
                     0.330720874, 6.030332501, 0.718903148)

expected_pc_samples <- c("ExpectedPC17", "ExpectedPC18", "ExpectedPC20")
abundance_expected_pc = abundance_expected_equimolar

# Create a forloop that replaces the expected abundance for each pccocktail
for (i in seq_along(species_17)) {
  abundance_expected_pc[which(abundance_expected_pc$species == species_17[i]), "ExpectedPC17"] <-  sum(abundance_expected_pc$`%cover17R1`)*(surface_area_17[i])/100
  
}

for (i in seq_along(species_18)) {
  abundance_expected_pc[which(abundance_expected_pc$species == species_18[i]), "ExpectedPC18"] <-  sum(abundance_expected_pc$`%cover18`)*(surface_area_18[i])/100
  
}

for (i in seq_along(species_20)) {
  abundance_expected_pc[which(abundance_expected_pc$species == species_20[i]), "ExpectedPC20"] <-  sum(abundance_expected_pc$`%cover20`)*(surface_area_20[i])/100
  
}

## I.4 Clean_up tables
expected_abundance <- abundance_expected_pc
rownames(expected_abundance) = expected_abundance$ID
expected_abundance$species <- NULL
expected_abundance[is.na(expected_abundance)]  <- 0

# Check if any ASV are missing
expected_verification <- taxonomy_expected$ID == expected_abundance$ID

if (all(expected_verification)){
  "expectation success"
} else {
  "expectation failure"
}
expected_taxonomy <- taxonomy_expected
length(unique(expected_taxonomy$species)) == nrow(expected_taxonomy)



##########################################################
########  J. Renmove duplicates in expected table ########
##########################################################

duplicated_species <- function (abundance, taxonomy, filename){
  # J.1 Delete duplicates
  abundance <- abundance
  taxonomy <- taxonomy
  
  taxonomy_colnames <- colnames(taxonomy)
  # Remove "ID" from the column names
  taxonomy_colnames <- setdiff(taxonomy_colnames, "ID")
  
  row.names(taxonomy) = taxonomy$ID
  
  ## Verify if there are duplicates
  nrow(taxonomy) == length(unique(row.names(taxonomy)))
  
  ## Join combined tax table and abundance_global tables together
  a.combined_table <- inner_join(abundance, taxonomy, by ="ID")
  
  ## Save abundance table: that one contains only the IDs that were retained in the pipeline (in abundance)
  a.final.abundance.table <-  a.combined_table %>%
    dplyr::select(colnames(abundance))
  row.names(a.final.abundance.table)<- a.final.abundance.table$ID
  
  unique_taxa <- a.combined_table %>%
    distinct(species)
  
  ## J.2 Delete duplicates in both the abundance and tax tables
  # Find the number of unique taxa
  rownames(a.combined_table) = a.combined_table$ID
  
  extensive_taxonomy <- taxonomy
  
  # when species are duplicated, pick the one that has the most extensive taxonomy
  # For each duplicated species, return the number of "unclassified" in each row, if not return NA 
  extensive_taxonomy$unclassified <- ifelse(duplicated(extensive_taxonomy$species) | duplicated(extensive_taxonomy$species, fromLast = TRUE),rowSums(extensive_taxonomy == "unclassified"), NA)
  
  extensive_taxonomy_verified <- subset(extensive_taxonomy, !(duplicated(species) & unclassified != min(unclassified)))
  extensive_taxonomy_verified$unclassified <- NULL
  
  # Verify if there are still duplicated species: should have zero columns in total
  duplicated_species_still <- extensive_taxonomy_verified[which(duplicated(extensive_taxonomy_verified$ID)),]
  
  # J.3 Map those taxa
  unique_taxa <- extensive_taxonomy_verified %>%
    distinct(species)
  unique_taxa$mapping_ASV <- rownames(unique_taxa)
  
  number_unique_taxa <- nrow(unique_taxa) 
  
  # Remove duplicates from taxonomy
  # Delete the ones that are not mapping taxa
  taxonomy_mapping <- extensive_taxonomy_verified %>%
    dplyr::filter(ID %in% unique_taxa$mapping_ASV)
  
  abundance_mapping <- a.combined_table %>%
    dplyr::select(c(colnames(abundance), species)) 
  
  abundance_mapping_join <- right_join(unique_taxa, abundance_mapping) %>%
    dplyr::select(-c("ID", "species"))
  
  # Merge abundance by mapping ASV
  abundance_no_duplicates_2 = aggregate (. ~ abundance_mapping_join$mapping_ASV, abundance_mapping_join[,!(colnames(abundance_mapping_join) %in% c("mapping_ASV"))], sum)
  colnames(abundance_no_duplicates_2)[1] <- "ID"
  rownames(abundance_no_duplicates_2) = abundance_no_duplicates_2$ID
  
  # Check if all the mapping taxa are present
  if (identical(sort(unique_taxa$mapping_ASV), sort(abundance_no_duplicates_2$ID)))  {
    print( "V3 good" )
  } else {
    print ("V3 is wrong")
  }
  
  # J.4 Saved finished tables
  # Save cleaned ASV table to environment
  assign(paste0(filename, "_taxonomy_table"),taxonomy_mapping, envir = .GlobalEnv)
  assign(paste0(filename, "_abundance_table"),abundance_no_duplicates_2, envir = .GlobalEnv)
  
  # Make output directory if not already made
  current_directory <- getwd()
  subdirectory_name <- "outputs"
  if (!file.exists(file.path(current_directory, subdirectory_name))) {
    # If it doesn't exist, create the subdirectory
    dir.create(file.path(current_directory, subdirectory_name))
  }
  
  # Save that df in the working directory
  write.csv(taxonomy_mapping, paste0("./outputs/", filename, "_taxonomy_table.csv"), row.names = FALSE)
  write.csv(abundance_no_duplicates_2, paste0("./outputs/", filename, "_abundance_table.csv"), row.names = FALSE)
}


duplicated_species(taxonomy = expected_taxonomy,
                   abundance = expected_abundance,
                   filename = "combined171820_expected_no_duplicates")



##################################################################
########  J. Reformat the 28s and COI data before merging ########
##################################################################

## J.1 Merge 3 ARMS for ARMS17-18-20
# Read meta data file (so that information can be changed in there too)
metadata_df_DNA <- as.data.frame(fread ("Martineau_TestingPrimers_28SANDCO1_Metadata.csv"))
rownames(metadata_df_DNA) <- metadata_df_DNA$Sample

# Merge 3 Taxonomic replicates for ARM 17 for abundance nad metadata
library("stringr")
ASV_table_DNA <- combined171820_expected_no_duplicates_abundance_table %>%
  mutate(Equimolar17 = rowMeans(select(., starts_with("Equimolar17R"))),
         `%cover17` = rowMeans(select(., starts_with("%cover17R"))),
         ARMS17 = rowMeans(select(., starts_with("ARMS17R")))) %>%
  select(-starts_with("Equimolar17R"), -starts_with("%cover17R"), -starts_with("ARMS17R"))

metadata_df_DNA = metadata_df_DNA[which(!(metadata_df_DNA$Sample %in% c("Equimolar17R2", "Equimolar17R3", 
                                                                        "%cover17R3", "%cover17R2",
                                                                        "ARMS17R2", "ARMS17R3"))), ]
metadata_df_DNA <- metadata_df_DNA %>%
  dplyr::mutate(Sample = str_replace_all(Sample, c("%cover17R1" = "%cover17",
                                                   "Equimolar17R1" = "Equimolar17R",
                                                   "ARMS17R1"=  "ARMS17")))
rownames(metadata_df_DNA) <- metadata_df_DNA$Sample


## J.2 Reformat the 28S table 
taxonomy.28s <- combined171820_expected_no_duplicates_taxonomy_table

abundance.28S <- ASV_table_DNA %>%
  dplyr::rename(ARMS1728S = ARMS17,
                ARMS1828S = ARMS18,
                ARMS2028S = ARMS20,
                Negative1820 = negative,
                Negative17 = Negative)
abundance.28S$ID <- row.names(abundance.28S)


## J.3 Reformat COI table: change name
#The COI table was found at this path: U_of_H/SPONGES/Testing_primers/Data_analysis/Bioinformatics_040523/CO1_Taxonomy/CO1_Taxonomy_and_clean_Tables/072024")
taxonomy.CO1.raw <- as.data.frame(fread ("sponge_ntCOI_taxonomy_table.csv",  header = TRUE)) #upload clean tax table
rownames(taxonomy.CO1.raw) = taxonomy.CO1.raw$ID

# Replace NAS with unclassified
abundance.CO1.raw <- as.data.frame(fread ("sponge_ntCOI_abundance_table.csv",header = TRUE)) #upload clean tax table
rownames(abundance.CO1.raw) = abundance.CO1.raw$ID

# Only select the 3 ARMS that I want and rename columns and delete OTUs not present in those 3 ARMS
abundance.CO1.raw.1 <- abundance.CO1.raw %>%
  mutate(across(starts_with("ARMS"), as.numeric)) %>%
  dplyr::rename(ARMS17CO1 = ARMS17,
                ARMS18CO1 = ARMS18,
                ARMS20CO1 = ARMS20) %>%
  filter(rowSums(select(., starts_with("ARMS"))) > 0)
abundance.CO1 = abundance.CO1.raw.1

# Map the taxonomy_table to delete those uncessary OTUs
COI_combine <- left_join(abundance.CO1.raw.1, taxonomy.CO1.raw)
taxonomy.CO1 <- COI_combine %>%
  dplyr::select(colnames(taxonomy.CO1.raw))
row.names(taxonomy.CO1) = taxonomy.CO1$ID
taxonomy.CO1[taxonomy.CO1 == ''] <- NA
taxonomy.CO1[is.na(taxonomy.CO1)] <- "unclassified"


# Save the clean tables
write.csv(abundance.28S, "./28SCOI/abundance_all_28S.csv", row.names = FALSE)
write.csv(taxonomy.28s, "./28SCOI/taxonomy_all_28S.csv", row.names = FALSE)
write.csv(abundance.CO1,"./28SCOI/abundance_all_COI.csv", row.names = FALSE)
write.csv(taxonomy.CO1, "./28SCOI/taxonomy_all_COI.csv", row.names = FALSE)


###########################################
########  J. Merge with COI tables ########
###########################################

merge_taxo_abundance(taxonomyA = "28SCOI/taxonomy_all_28S.csv",
                     taxonomyB = "28SCOI/taxonomy_all_COI.csv",
                     abundanceA = "28SCOI/abundance_all_28S.csv",
                     abundanceB = "28SCOI/abundance_all_COI.csv",
                     filename = "combined28SCOI")
