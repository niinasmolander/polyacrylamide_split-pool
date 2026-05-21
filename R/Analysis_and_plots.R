library(blaster)
library(tidyverse)
library(Biostrings)
library(patchwork)

# load the functions
source("../input/BC_functions.R")

#### read FASTAs ####

#16S
fasta1 <- blaster::read_fasta("../input/S1_R1.fasta")

#A.t. specific
fasta_g1 <- blaster::read_fasta("../input/S2_R1.fasta")

#### BC sequences ####

bcs <- c("AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT", "ACATTGGC",
         "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA", "CTGTAGCC", "AGTACAAG",
         "AACAACCA", "AACCGAGA", "AACGCTTA", "AAGACGGA", "AAGGTACA", "ACACAGAA",
         "ACAGCAGA", "ACCTCCAA", "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA",
         "AGATCGCA", "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
         "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA", "CCAGTTCA",
         "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA", "CGACTGGA", "CGCATACA",
         "CTCAATGA", "CTGAGCCA", "CTGGCATA", "GAATCTGA", "CAAGACTA", "GAGCTGAA",
         "GATAGACA", "GCCACATA", "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA",
         "GGTGCGAA", "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
         "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA", "TGGCTTCA",
         "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC", "AAGGACAC", "AATCCGTC",
         "AATGTTGC", "ACACGACC", "ACAGATTC", "AGATGTAC", "AGCACCTC", "AGCCATGC",
         "AGGCTAAC", "ATAGCGAC", "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC",
         "CCATCCTC", "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
         "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC", "GCCAAGAC")

# reverse the bc sequences
bcs2 <- as.data.frame(Biostrings::reverseComplement(Biostrings::DNAStringSet(bcs))) %>% pull(x) 


#### Filter out bad BCs ####

# 16S
seq_filt <- fasta1 %>%
  filter(str_detect(Seq, "CGAATGCTCTGGCCCTCAAGCACGTGGAT")) %>% #non-variable region between BC3 and BC2
  filter(str_detect(Seq, "AGTCGTACGCCGATGCGAAACATCGGCCAC")) %>% #non-variable region between BC2 and BC1
  filter(str_detect(Seq, "ACCCTTGCTCAGAACACCACGCTCCAATTAAGCG")) %>% #non-variable region between BC1 and amplicon
  separate_wider_delim(Seq, "CGAATGCTCTGGCCCTCAAGCACGTGGAT", names = c("BC3", "rest")) %>%
  mutate(start = str_sub(BC3, 1, -9),
         BC3 = str_sub(BC3, -8, -1)) %>% #separate the BC3
  separate_wider_delim(rest, "AGTCGTACGCCGATGCGAAACATCGGCCAC", names = c("BC2", "rest")) %>%
  separate_wider_delim(rest, "ACCCTTGCTCAGAACACCACGCTCCAATTAAGCG", names = c("BC1", "bact")) %>%
  filter(nchar(BC1) == 18 & nchar(BC2) == 8 & nchar(BC3) == 8) %>% #BC1 includes a UMI
  mutate(UMI = str_sub(BC1, -10, -1),
         BC1 = str_sub(BC1, 1, 8)) %>% #separate the UMI and BC1
  filter(BC1 %in% bcs2 & BC2 %in% bcs2 & BC3 %in% bcs2) %>% #remove impossible BCs
  mutate(full_BC = paste0(BC3, BC2, BC1)) %>%
  relocate(Id, full_BC, start, BC3, BC2, BC1, UMI)

# A.t. specific
seq_filt_g <- fasta_g1 %>%
  filter(str_detect(Seq, "CGAATGCTCTGGCCCTCAAGCACGTGGAT")) %>% #non-variable region between BC3 and BC2
  filter(str_detect(Seq, "AGTCGTACGCCGATGCGAAACATCGGCCAC")) %>% #non-variable region between BC2 and BC1
  filter(str_detect(Seq, "ACCCTTGCTCAGAACACCACGCTCCAATTAAGCG")) %>% #non-variable region between BC1 and amplicon
  separate_wider_delim(Seq, "CGAATGCTCTGGCCCTCAAGCACGTGGAT", names = c("BC3", "rest")) %>%
  mutate(start = str_sub(BC3, 1, -9),
         BC3 = str_sub(BC3, -8, -1)) %>% #separate the BC3
  separate_wider_delim(rest, "AGTCGTACGCCGATGCGAAACATCGGCCAC", names = c("BC2", "rest")) %>%
  separate_wider_delim(rest, "ACCCTTGCTCAGAACACCACGCTCCAATTAAGCG", names = c("BC1", "bact")) %>%
  filter(nchar(BC1) == 18 & nchar(BC2) == 8 & nchar(BC3) == 8) %>% #BC1 includes a UMI
  mutate(UMI = str_sub(BC1, -10, -1),
         BC1 = str_sub(BC1, 1, 8)) %>% #separate the UMI and BC1
  filter(BC1 %in% bcs2 & BC2 %in% bcs2 & BC3 %in% bcs2) %>% #remove impossible BCs
  mutate(full_BC = paste0(BC3, BC2, BC1)) %>%
  relocate(Id, full_BC, start, BC3, BC2, BC1, UMI)

# How many BCs out of the possible 96?
length(unique(c(seq_filt$BC1, seq_filt_g$BC1))) #BC1
length(unique(c(seq_filt$BC2, seq_filt_g$BC2))) #BC2
length(unique(c(seq_filt$BC3, seq_filt_g$BC3))) #BC3

#### Assign the strain ####

# reference sequences
DB <- blaster::read_fasta("../input/ref_sequences.fasta")

# format the 16S query 
query <- seq_filt %>%
  select(Id, Seq = bact)

# format the A.t. specific query 
query_g <- seq_filt_g %>%
  select(Id, Seq = bact)

# 16S and A.t. specific analysed separately, hence different blast databases
DB_g <- DB[4,]
DB <- DB[1:3,]

# run blaster
blast_res <- blaster::blast(query, DB, maxAccepts = 3, minIdentity = 0.95)

blast_res_g <- blaster::blast(query_g, DB_g, maxAccepts = 3, minIdentity = 0.95)

# remove seqs assigned to multiple targets
blast_res <- blast_res %>%
  filter(!duplicated(QueryId))

# correct the names
annotated <- left_join(seq_filt, select(blast_res, QueryId, TargetId), by = join_by(Id == QueryId)) %>%
  filter(!is.na(TargetId)) %>%
  mutate(TargetId = str_split_i(TargetId, "_", 1),
         TargetId = str_replace_all(TargetId, "-", "_"),
         TargetId = dplyr::case_when(TargetId == "HAMBI_0403" ~ "C. testosteroni",
                                     TargetId == "HAMBI_0105" ~ "A. tumefaciens",
                                     TargetId == "HAMBI_0262" ~ "B. bullata"))

annotated_g <- left_join(seq_filt_g, select(blast_res_g, QueryId, TargetId), by = join_by(Id == QueryId)) %>%
  filter(!is.na(TargetId)) %>%
  mutate(TargetId = dplyr::case_when(TargetId == "HAMBI-0105-Specific" ~ "A. tumefaciens specific"))


# combine the 16S and A.t. specific
full_annotated <- bind_rows(annotated, annotated_g)

#### Some statistics ####

# how many full BCs (BC3,BC2 and BC1 concatenated)
length(unique(full_annotated$full_BC))

# how many times each BC sampled
BC_stats <- full_annotated %>% group_by(full_BC) %>%
  summarise(n = n()) %>% ungroup()

summary(BC_stats$n)

# how many UMIs
length(unique(full_annotated$UMI))

#### UMIs ####

# UMI distribution
UMI_stats <- full_annotated %>%
  group_by(UMI) %>%
  mutate(n_targets_per_UMI = n_distinct(TargetId),
         n_BCs_per_UMI = n_distinct(full_BC),
         n_sampled = n()) %>%
  ungroup() %>%
  select(UMI, n_targets_per_UMI, n_BCs_per_UMI, n_sampled) %>%
  distinct()

# How many UMIs associated with only one target and one BC
sum(UMI_stats$n_targets_per_UMI == 1 & UMI_stats$n_BCs_per_UMI == 1)

summary(UMI_stats$n_sampled)

# Prepare filtering column based on how frequently each UMI is associated with each 
# target and BC. Each UMI should have a strong association with only one BC / target.

umi_counts_target <-  full_annotated %>%
  group_by(UMI) %>%
  mutate(n_targets_per_UMI = n_distinct(TargetId)) %>%
  ungroup() %>%
  group_by(UMI, TargetId) %>%
  mutate(n_UMI = n()) %>% # how many times each UMI associated with each target
  ungroup() %>%
  select(UMI, n_targets_per_UMI, TargetId, n_UMI) %>%
  distinct() %>%
  group_by(UMI) %>%
  mutate(umifilter = n_UMI / sum(n_UMI)) %>% #how frequently the UMI is associated with a specific target as opposed to all targets
  ungroup()

umi_counts_BC <- full_annotated %>%
  group_by(UMI) %>%
  mutate(n_BCs_per_UMI = n_distinct(full_BC),
         n_BC1_per_UMI = n_distinct(BC1),
         n_BC2_per_UMI = n_distinct(BC2)) %>% # how many times each UMI associated with BCs
  ungroup() %>%
  group_by(UMI, full_BC) %>%
  mutate(n_UMI = n()) %>%
  ungroup() %>%
  select(UMI, n_BCs_per_UMI, n_BC1_per_UMI, n_BC2_per_UMI, full_BC, n_UMI) %>%
  distinct() %>%
  group_by(UMI) %>%
  mutate(umifilter = n_UMI / sum(n_UMI)) %>% #how frequently the UMI is associated with a specific BC as opposed to all BCs
  ungroup()
  
##### Filter based on UMIs #####
# This filtering removes those cases where UMI is clearly associated with multiple BCs / targets.
# It is based on how frequent the association is compared to all of the associations
# of the UMI. Hence, if a UMI is associated with one target once and with another target
# 20 times, the read with only one association is removed as it is likely not a genuine
# UMI-target or UMI-BC association. Conversely, if a UMI is associated 20 times with one target
# and 20 times with another target, all of those are filtered out.

# BC-UMI
umis_BC <- umi_counts_BC %>%
  filter(umifilter >= 0.7) %>%
  group_by(UMI) %>%
  mutate(new_BC_per_UMI = n_distinct(full_BC)) %>% #calculate new BCs per UMI after filtering
  ungroup()
# check that all new_BC_per_UMI == 1

umis_BC <- umis_BC %>%
  select(UMI, full_BC) %>% distinct()

# strain-UMI
umis_strain <- umi_counts_target %>%
  filter(umifilter >= 0.7) %>%
  group_by(UMI) %>%
  mutate(new_target_per_UMI = n_distinct(TargetId)) %>%
  ungroup()
# check that all new_target_per_UMI == 1

umis_strain <- umis_strain %>%
  select(UMI, TargetId) %>% distinct()

# How many unique UMIs
length(unique(c(umis_strain$UMI, umis_BC$UMI)))

# prepare the final df with the following columns: full_BC (sequence), TargetId,
# n_associated (number of times each BC associates with each target) and
# BC_shared_between (how many targets share the BC).

final_df <- full_annotated %>%
  inner_join(umis_BC, by = join_by(full_BC, UMI)) %>% # UMI-BC filtering
  inner_join(umis_strain, by = join_by(TargetId, UMI)) %>% # UMI-strain filtering
  group_by(full_BC, TargetId) %>%
  mutate(n_associated = n()) %>% # how many times each BC associated with each target
  ungroup() %>%
  select(full_BC, TargetId, n_associated) %>% distinct() %>%
  group_by(full_BC) %>%
  mutate(BC_shared_between = n_distinct(TargetId)) %>%
  ungroup()
  
# how many unique
length(unique(final_df$full_BC))

# how many BC associated with 3 targets?
final_df %>%
  summarise(n = n_distinct(full_BC[BC_shared_between == 3]))

# how many BC associated with 4 targets
final_df %>%
  summarise(n = n_distinct(full_BC[BC_shared_between == 4]))

#### BC distribution table ####

# There are 4 targets, hence setting the max number of target combinations as 4. 
BC_associations <- BC_tables(final_df,4)

# how many BCs associated with each target (combination) and how many BCs
# associated ≥ 10 times with the targets (solely).
BC_associations$n_table

# For each target pair, how many BCs associated with either of them (n_all column)
# and percentage of BC associated solely with both of the respective targets (prop)
# at least 10 times (prop_min_10)
BC_associations$p_table

# save
write_tsv(BC_associations$n_table, "../output/BC_associations.tsv")
write_tsv(BC_associations$p_table, "../output/BC_association_pairwise_proportions.tsv")

#### BC distribution plots ####

# Prepare the plots
BC_plots <- BC_plots(final_df)

# adjust the order of the plots to match the preprint
plot_order <- c(4,1,2,6,5,3)

# Blue points highlight the BCs associated with any three targets
# Orange points with all four targets

wrap_plots(BC_plots[plot_order], ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

# save
ggsave(filename = "../output/BC_association_plots.pdf",
       device = "pdf", dpi = 300, width = 10, height = 7, units = "in")

