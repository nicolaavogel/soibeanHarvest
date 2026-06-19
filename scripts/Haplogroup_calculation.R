
library(tidyverse)

# simulation details file
sim <- read_delim("simulation_details.csv", delim=";")

# reference genome → species → haplogroup table
ref <- read_delim("ref_mitogenome_clade_info.csv", delim=";")

################################################################################

### For 2-source simulations:

################################################################################

#filter for available result files
sim_2s <- filter(sim, !is.na(sim$`Result file k1`))

################################################################################

#### for k = 1

# Function: read one result file and annotate with reference info
read_and_annotate <- function(filename, ref_table) {
  # Construct full path
  full_path <- file.path("Beluga_2_res", filename)
  
  read_lines(full_path) %>%
    tibble(genome = .) %>%
    left_join(ref_table, by = "genome")}

# Process ONLY the k1 result files
results <- sim_2s %>%
  mutate(row_id = row_number()) %>%
  mutate(
    k1_data = map(`Result file k1`, read_and_annotate, ref_table = ref))

target_species <- "beluga"

proportions <- results %>%
  transmute(
    row_id,
    Species = Species,
    No_of_sources = `No. of sources`,
    haplotype_proportions = map(
      k1_data,
      ~ .x %>%
        filter(species == target_species) %>%
        count(haplogroup, name = "n") %>%
        mutate(proportion = n / sum(n))))

final_table <- proportions %>% unnest(haplotype_proportions)
print(final_table)

#################################################################

#### for k = 2

read_k2 <- function(filename, ref_table) {
  full_path <- file.path("Beluga_2_res", filename)
  
  read_lines(full_path) %>%
    tibble(raw = .) %>%
    separate(raw, into = c("genome1", "genome2"), sep = "\t", remove = TRUE) %>%
    pivot_longer(cols = c(genome1, genome2),
                 values_to = "genome",
                 names_to = "source_index") %>%   # optional indicator: first vs second genome
    drop_na(genome) %>%
    left_join(ref_table, by = "genome")
}

results <- sim_2s %>%
  mutate(row_id = row_number()) %>%
  mutate(
    #k1_data = map(`Result file k1`, read_k1, ref_table = ref),
    k2_data = map(`Result file k2`, read_k2, ref_table = ref)
  )


compute_props <- function(df, species) {
  df %>%
    filter(species == !!species) %>%
    count(haplogroup, name = "n") %>%
    mutate(proportion = n / sum(n))
}

proportions <- results %>%
  transmute(
    row_id,
    Species,
    No_of_sources = `No. of sources`,
    #k1_haplo = map(k1_data, compute_props, species = target_species),
    k2_haplo = map(k2_data, compute_props, species = target_species)
  )


k2_table <- proportions %>%
  select(row_id, Species, No_of_sources, k2_haplo) %>%
  unnest(k2_haplo)


################################################################################
################################################################################


################################################################################

### For 1-source simulations:

################################################################################

#filter for available result files
sim_1s <- filter(sim, sim$`No. of sources` ==1)
sim_1s <- select(sim_1s, `Species`, `Source genome 1`, `Simulated sequences`, `Haplogroup source genome 1`)

# fix ref double __ for harp seal genome IDs
ref <- ref %>%
  mutate(genome = gsub("Phoca2cut", "harpSeal", genome)) %>%
  mutate(genome = gsub("_P", ".1_P", genome)) %>%
  mutate(genome = gsub("isolate.1_", "isolate_", genome))

# Base directory containing the species subfolders
base_path <- "Beluga_2_res"

# 1. Recursively list all Diagnostics10 files in *any* subdirectory
all_files <- list.files(
  path = base_path,
  pattern = "Diagnostics10\\.txt$",
  recursive = TRUE,
  full.names = TRUE)

# 2. Extract file and folder information
file_info <- tibble(
  full_path = all_files,
  file = basename(all_files),
  
  # Name of immediate parent folder (e.g., beluga, harp, narwhal)
  folder = basename(dirname(all_files)))

# 3. Parse metadata from filenames:
# Example: Beluga_Capture_50x_1HighDiagnostics10.txt
file_info <- file_info %>%
  mutate(
    species = str_extract(file, "(?i)beluga|harp|narwhal") %>% tolower(),
    seq_type = str_extract(file, "(?i)Capture|Shotgun") %>% tolower(),
    reads = str_extract(file, "\\d+(?=x)") %>% as.integer(),
    
    # FIXED: Extract bootstrap without look-behind
    bootstrap = str_extract(file, "_\\d+x_(\\d+)") %>%
      str_extract("\\d+$") %>%
      as.integer(),
    
    damage = str_extract(file, "(?i)low|high|none") %>% tolower())

# 4. remove files which were not part of the single source simulations:
file_info <- filter(file_info, !is.na(reads))
# 4b. keep only beluga for now
file_info <- filter(file_info, species == "beluga")

# 5. Read results and add to dataframe

combined_results <- file_info %>%
  arrange(species, seq_type, reads, damage, bootstrap) %>%
  group_by(species, seq_type, reads, damage) %>%
  group_modify(~{
    
    files <- .x$full_path
    
    # Read all files as raw lines
    lines_list <- map(files, readLines)
    
    # Keep header from first file only
    header <- lines_list[[1]][1]
    
    # Drop headers from all files
    data_lines <- map(lines_list, ~ .x[-1])
    
    # Combine everything
    all_lines <- c(header, unlist(data_lines))
    
    # Parse once, consistently
    read_tsv(
      I(all_lines),
      col_types = cols(),
      progress = FALSE)})

# 6. add haplogroup info

combined_results <- combined_results %>%
  left_join(select(ref, "genome", "haplogroup"), by = c("Source" = "genome"))

# 7. Adjust metadata

sim_meta <- sim_1s %>%
  rename(
    species = Species,
    source_genome = `Source genome 1`,
    simulated_reads = `Simulated sequences`,
    true_haplogroup = `Haplogroup source genome 1`
  )

# 8. Add metadata to results table
results_annotated <- combined_results %>%
  rename(
    predicted_source = Source,
    simulated_reads = reads,
    predicted_haplogroup = haplogroup
  ) %>%
  left_join(
    sim_meta,
    by = c("species", "simulated_reads")
  )

# 9. Calculate proportions

# 9. Calculate proportions

haplogroup_proportions <- results_annotated %>%
  group_by(
    species,
    seq_type,
    simulated_reads,
    damage,
    true_haplogroup
  ) %>%
  count(predicted_haplogroup, name = "n") %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

#haplogroup_proportions <- results_annotated %>%
#  group_by(species, seq_type, simulated_reads, damage, true_haplogroup) %>%
 # count(predicted_haplogroup, name = "n") %>%
 # mutate(proportion = n / 500000) %>%
 # ungroup()

# sanity check: sum of proportions = 1

summed_proportions <- haplogroup_proportions %>%
  group_by(
    species,
    seq_type,
    simulated_reads,
    damage
  ) %>%
  summarise(
    summed_proportion = sum(proportion, na.rm = TRUE),
    .groups = "drop")

# 10. Save file / Write table

#write.csv(haplogroup_proportions, "haplogroup_proportions_single_source_from_Diagnostics10_files.csv")

