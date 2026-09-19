###Libraries and Database####
library(phytools)
library(tidyverse)
library(DBI)
library(RPostgres)
library(dbplyr)
library(simona)
library(simplifyEnrichment)

###First restore the database available here - https://russelllab.org/funcevol/ using the pg_restore command
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_revision", bigint = "integer")

species_tree <- sptree_revised <- read.tree("data/species_tree_cleaned_final.nwk")
species_details <- read_tsv("data/TableS1.tsv")
major_clades <- read_tsv("data/TableS2.tsv")
# go_obo <- ontologyIndex::get_ontology(file = "data/go-basic.obo", propagate_relationships = c("is_a", "part_of"))

full_go_lineage <- function(term, go_obo) {
  xx <- unlist(go_obo$name[unlist(go_obo$ancestors[go_obo$id[match(term, go_obo$name)]])])
  xx
}


###Fig 4B-C####

euk_tree <- keep.tip(sptree_revised, species_details$tree_tip_label[species_details$superkingdom == "Eukaryota"])

####Get the terrestrial transition nodes####
load("data/20260722_habitat_anc_reconstruction.simmap")

habitat_anc_summary <- summary(habitat_anc)

ace_matrix <- habitat_anc_summary$ace

habitat_transitions <- euk_tree$edge %>% 
  as_tibble() %>% 
  mutate(V1 = as.character(V1), V2 = if_else(V2 <= Ntip(euk_tree), euk_tree$tip.label[V2], as.character(V2))) %>% 
  left_join(ace_matrix %>% as_tibble(rownames = "node"), by = c("V1" = "node"))

p_state <- NULL
for(x in 1:nrow(habitat_transitions)) {
  row_oi <- habitat_transitions[x,] %>% pivot_longer(cols = c(-V1, -V2))
  if(any(row_oi$value >= 0.5) | all(row_oi$V1 == 200)) {
    p_state <- rbind(
      p_state,
      tibble(
        V1 = unique(row_oi$V1),
        V2 = unique(row_oi$V2),
        p_state = row_oi$name[row_oi$value == max(row_oi$value)],
        p_pp = max(row_oi$value)
      )
    )
  } else {
    y <- which(p_state$V2 == unique(row_oi$V1))
    parent_row <- p_state[y,] 
    p_state <- rbind(
      p_state,
      tibble(
        V1 = unique(row_oi$V1),
        V2 = unique(row_oi$V2),
        p_state = parent_row$p_state,
        p_pp = row_oi$value[row_oi$name == parent_row$p_state]
      ) 
    )
  }
}

habitat_transitions <- habitat_transitions %>% 
  left_join(p_state) %>% 
  select(p_node = V1, V2, p_state, p_pp) %>% 
  left_join(ace_matrix %>% as_tibble(rownames = "node"), by = c("V2" = "node"))

d_state <- NULL
for(x in 1:nrow(habitat_transitions)) {
  row_oi <- habitat_transitions[x,] %>% pivot_longer(cols = -c(p_node:p_pp))
  if(any(row_oi$value >= 0.5) | all(row_oi$p_node == 200)) {
    d_state <- rbind(
      d_state,
      tibble(
        p_node = unique(row_oi$p_node),
        V2 = unique(row_oi$V2),
        d_state = row_oi$name[row_oi$value == max(row_oi$value)],
        d_pp = max(row_oi$value)
      )
    )
  } else {
    y <- which(d_state$V2 == unique(row_oi$p_node))
    parent_row <- d_state[y,] 
    d_state <- rbind(
      d_state,
      tibble(
        p_node = unique(row_oi$p_node),
        V2 = unique(row_oi$V2),
        d_state = parent_row$d_state,
        d_pp = row_oi$value[row_oi$name == parent_row$d_state]
      ) 
    )
  }
}

habitat_transitions <- habitat_transitions %>% 
  left_join(d_state) %>% 
  select(p_node, d_node = V2, p_state, p_pp, d_state, d_pp) %>% 
  mutate(change = if_else(d_state != p_state, 1, 0))

aq2ter_nodes <- habitat_transitions %>% 
  rowwise() %>% 
  mutate(n_desc = sum(!is.na(euk_tree$tip.label[getDescendants(euk_tree, d_node)]), na.rm = T)) %>% 
  ungroup() %>% 
  filter(n_desc > 2, change == 1, d_state == "Terrestrial") %>% 
  pull(d_node) %>% 
  as.numeric() %>% 
  na.omit() %>% 
  as.vector()

##change node ids to the wider species tree
aq2ter_nodes_orig <- gsub("Node", "", euk_tree$node.label[aq2ter_nodes - Ntip(euk_tree)])
aq2ter_nodes_parents <- gsub("Node", "", euk_tree$node.label[sapply(aq2ter_nodes, function(x) getParent(euk_tree, x)) - Ntip(euk_tree)])

#####Analyze GO-BP for terrestrial nodes####

aq2ter_all_parents <- list()
for(x in aq2ter_nodes_orig) {
  parent_oi <- getParent(species_tree, as.numeric(x))
  aq2ter_all_parents[[x]] <- c(aq2ter_all_parents[[x]], parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(species_tree, parent_oi)
    aq2ter_all_parents[[x]] <- c(aq2ter_all_parents[[x]], parent_oi)
  }
}


aq2ter_go_bp_term <- tbl(con, "og_gains_losses_parsimony") %>% 
  pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
  mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
  # filter(node %in% !!c(aq2ter_nodes_orig)) %>% 
  filter(node %in% !!c(aq2ter_nodes_orig, unique(unlist(aq2ter_all_parents)))) %>%
  left_join(tbl(con, "ptns_og_phylo") %>% select(HOG, org, uniprot_acc:gene_desc), by = join_by(HOG)) %>% 
  inner_join({
    tbl(con, "orthologs_gains_losses_parsimony") %>% 
      pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
      mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
      filter(node %in% !!c(aq2ter_nodes_orig, unique(unlist(aq2ter_all_parents)))) %>%
      # filter(node %in% !!c(aq2ter_nodes_orig)) %>% 
      distinct(uniprot_acc, node)
  }) %>% 
  inner_join(tbl(con, "ptns_go_bp") %>% select(uniprot_acc, bp, term = go_id) %>% dplyr::rename(value = bp)) %>% 
  # select(-level, -matches("l\\d+")) %>% 
  # dplyr::rename(value = bp) %>% 
  # filter(!is.na(value)) %>% 
  distinct() %>% 
  collect() %>% 
  group_by(node) %>% 
  filter(
    case_when(
      all(node %in% aq2ter_nodes_orig) ~ (
        org %in% (
          habitat_species %>% 
            filter(OSCODE %in% sptree_revised$tip.label[getDescendants(sptree_revised, node = unique(node))]) %>% 
            filter(habitat == "Terrestrial") %>% 
            pull(OSCODE)
        )
      ),
      T ~ T
    )
  ) %>% 
  ungroup()

write_tsv(aq2ter_go_bp_term, "data/20260723_aq2ter_bp_terms_by_group.tsv.gz")

source("scripts/fishers_test_fast.R")

infile = "data/20260723_aq2ter_bp_terms_by_group.tsv.gz"
outfile = "data/20260723_aq2ter_fisher_test_results_bp.csv.gz"
require(data.table)
data <- fread(infile)

# Step 2: Create mappings of nodes to orthogroups and GO terms to orthogroups
node_to_orthogroups <- data %>% distinct(HOG, node) %>% with(., split(HOG, node))
go_to_orthogroups <- data %>% distinct(HOG, value) %>% with(., split(HOG, value))
node_to_go <- data %>% distinct(value, node) %>% with(., split(value, node))
cat("Mappings created...\n")

# Get all unique orthogroups
all_orthogroups <- unique(data$HOG)
run_fast_fishers(infile = infile, outfile = outfile)

aq2ter_go_bp_term_test <- read_csv("data/20260723_aq2ter_fisher_test_results_bp.csv.gz", col_names = c(
  "node",
  "value",
  "has_term_at_node",
  "not_term_at_node",
  "has_term_not_at_node",
  "not_term_not_at_node",
  "odds_ratio",
  "p_value"
), col_types = c("node" = "c"), skip = 1)
aq2ter_go_bp_term_test <- aq2ter_go_bp_term_test %>% 
  mutate(padj = p.adjust(p_value, method = "BH"))

aq2ter_go_bp_term_sig <- aq2ter_go_bp_term_test %>% 
  filter(p_value < (0.05), has_term_not_at_node > 0, has_term_at_node + has_term_not_at_node > 4, odds_ratio > 1)

gained_func_term_tbl <- NULL
for(i in 1:length(aq2ter_nodes_orig)) {
  desc_func <- aq2ter_go_bp_term_sig %>% 
    filter(node == aq2ter_nodes_orig[i]) %>% 
    pull(value) %>% 
    unique()
  
  row_add <- tibble(node = aq2ter_nodes_orig[i], name = major_clades$name[major_clades$node_name == node], gained_func = desc_func)
  
  gained_func_term_tbl <- rbind(gained_func_term_tbl, row_add)
}

gained_func_term_tbl <- gained_func_term_tbl %>% 
  inner_join(aq2ter_go_bp_term_test, by = join_by(node, gained_func == value)) %>% 
  unnest()

annotation_tbl <- tbl(con, "ptns_go_bp") %>%
  select(uniprot_acc, go_id, term = bp) %>%
  filter(term %in% !!gained_func_term_tbl$gained_func) %>%
  inner_join(
    (tbl(con, "orthologs_gains_losses_parsimony") %>%
       pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>%
       mutate(node = sql("unnest(string_to_array(node, ', '))")) %>%
       filter(node %in% !!c(aq2ter_nodes_orig)) %>%
       select(uniprot_acc))
  ) %>%
  collect() %>%
  distinct()

annotation_list <- annotation_tbl %>%
  with(., split(uniprot_acc, factor(go_id, levels = unique(go_id))))
require(simona)
dag <- import_obo("data/go.obo", relation_type = c("is_a", "part_of"), annotation = annotation_list)
all_terms <- unique(annotation_tbl$go_id)
all_terms_go <- all_terms[!is.na(all_terms)]
all_terms_sim <- term_sim(dag = dag, terms = all_terms_go, method = "Sim_Lin_1998")
all_terms_IC <- term_IC(dag = dag, method = "IC_offspring")

require(simplifyEnrichment)
all_terms_sim_clustered <- cluster_terms(all_terms_sim, method = "louvain", control = list(resolution = 1))

all_terms_sim_clustered_dist <- annotation_tbl %>%
  filter(!is.na(go_id)) %>%
  mutate(go_cluster = all_terms_sim_clustered[match(go_id, rownames(all_terms_sim))]) %>%
  distinct(term, go_id, go_cluster) %>%
  drop_na(go_cluster) %>%
  group_by(go_cluster) %>%
  filter(n() > 5) %>% 
  do({
    sub_tbl <- .
    sub_tbl %>%
      rowwise() %>%
      do({
        sub_sub_tbl <- .
        tibble(term = sub_sub_tbl$term, go_id = sub_sub_tbl$go_id, go_cluster = sub_sub_tbl$go_cluster, avg_dist_to_others = mean(all_terms_sim[sub_sub_tbl$go_id, sub_tbl$go_id[sub_tbl$go_id != sub_sub_tbl$go_id]]))
      })
  })

gained_func_term_tbl2 <- gained_func_term_tbl %>% 
  inner_join(all_terms_sim_clustered_dist, join_by(gained_func == term)) %>% 
  group_by(go_cluster) %>% 
  ungroup()

parents_func_term_tbl <- NULL
for(i in 1:length(aq2ter_nodes_orig)) {
  parents1 <- getParent(species_tree, aq2ter_nodes_orig[i])
  parents_oi <- parents1
  while(!is.null(parents1)) {
    parents1 <- getParent(species_tree, parents_oi[length(parents_oi)])
    parents_oi <- c(parents_oi, parents1)
  }
  parent_func <- aq2ter_go_bp_term_test %>% 
    filter(node %in% parents_oi, value %in% gained_func_term_tbl2$gained_func) %>% 
    dplyr::rename(parent_node = node, parent_func = value) %>% 
    mutate(node = aq2ter_nodes_orig[i], .before = everything())
  
  row_add <- parent_func
  
  parents_func_term_tbl <- rbind(parents_func_term_tbl, row_add)
}

parents_func_term_tbl2 <- parents_func_term_tbl %>% 
  inner_join(all_terms_sim_clustered_dist, join_by(parent_func == term)) %>% 
  filter(p_value < 0.05)

col_palette <- setNames(rep(RColorBrewer::brewer.pal(n = length(aq2ter_nodes_orig), name = "Dark2"), 1) , c(aq2ter_nodes_orig))

aq2ter_go_bp_odds_tbl_kw_cands <- gained_func_term_tbl2 %>% 
  distinct(go_cluster, gained_func) %>% 
  rowwise() %>% 
  mutate(parent_terms = list(full_go_lineage(term = gained_func, go_obo = go_obo))) %>% 
  unnest() %>% 
  group_by(go_cluster) %>%
  mutate(n_go = n_distinct(gained_func)) %>%
  ungroup() %>% 
  group_by(go_cluster, parent_terms) %>% 
  summarize(n_pt = n(), n_go = unique(n_go), prop = n_pt/n_go, .groups = "drop") %>% 
  arrange(go_cluster, desc(prop)) %>% 
  mutate(IC = all_terms_IC[go_obo$id[match(parent_terms, go_obo$name)]]) %>% 
  group_by(go_cluster) %>%
  ungroup()

aq2ter_go_bp_kw_cands_final <- aq2ter_go_bp_odds_tbl_kw_cands %>% 
  filter(prop > 0.5) %>% 
  group_by(go_cluster) %>%
  mutate(rel_IC = IC/max(IC)) %>% 
  do({
    sub_tbl <- .
    
    print(sub_tbl)
    sub_tbl %>% 
      slice(as.numeric(readline("Choose the row to keep: ")))
  })

go_bp_odds_tbl <- gained_func_term_tbl2 %>% 
  group_by(go_cluster) %>% 
  filter(n_distinct(gained_func) > 2) %>% 
  ungroup() %>% 
  inner_join(aq2ter_go_bp_kw_cands_final %>% select(go_cluster, keywords = parent_terms)) %>% 
  left_join(
    parents_func_term_tbl2 %>% 
      group_by(go_cluster) %>% 
      filter(n_distinct(parent_func) > 2) %>% 
      ungroup() %>% 
      inner_join(aq2ter_go_bp_kw_cands_final %>% select(go_cluster, keywords = parent_terms)),
    join_by(node, go_id, go_cluster, avg_dist_to_others, keywords),
    suffix = c("_gained", "_parent")
  )

f4b1 <- gained_func_term_tbl2 %>% 
  inner_join(aq2ter_go_bp_kw_cands_final %>% select(go_cluster, keywords = parent_terms, prop)) %>% 
  group_by(node, name, keywords, prop) %>% 
  summarize(mean_log2_odds_ratio = mean(log2(odds_ratio), na.rm = T), .groups = "drop") %>% 
  group_by(keywords) %>% 
  mutate(sum_or = sum(mean_log2_odds_ratio, na.rm = T), label = str_wrap(paste0(keywords, "(", round(prop*100, 1), "%)"), 20)) %>% 
  ungroup() %>% 
  ggplot(aes(x = name, y = fct_reorder(label, sum_or), fill = mean_log2_odds_ratio)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(1, NA)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Transition node", y = "MICA (% terms parent to)", fill = "Mean Log2 Odds Ratio")

sf7a <- gained_func_term_tbl2 %>% 
  inner_join(uni2multi_go_bp_kw_cands_final %>% select(go_cluster, keywords = parent_terms)) %>% 
  group_by(gained_func) %>%
  filter(n_distinct(name) > 1) %>% 
  mutate(sum_or = sum(log2(odds_ratio), na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = name, y = fct_reorder(gained_func, sum_or), fill = log2(odds_ratio))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Transition node", y = "GO term", fill = "Log2 Odds Ratio")

##fig4b2
# ptns_gained_at_trans <- tbl(con, "orthologs_gains_losses_parsimony") %>% 
#   pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
#   mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
#   filter(node %in% !!c(aq2ter_nodes_orig)) %>% 
#   distinct(uniprot_acc) %>% 
#   pull()

# terms_oi <- go_bp_odds_tbl %>% 
#   distinct(node, name, gained_func, keywords)

aq2ter_go_bp_term <- read_tsv("data/20260723_aq2ter_bp_terms_by_group.tsv.gz")

# ptn_oi <- aq2ter_go_bp_term %>% 
#   filter(node %in% aq2ter_nodes_orig, uniprot_acc %in% ptns_gained_at_trans) %>% 
#   inner_join(terms_oi %>% mutate(node = as.numeric(node)), join_by(node, value == gained_func)) %>% 
#   distinct(node, name, uniprot_acc, HOG, value, keywords)
# 
# export_tbl <- tbl(con, "ptns_og_phylo") %>% 
#   inner_join(ptn_oi, join_by(uniprot_acc, HOG), copy = T) %>% 
#   filter(gene_desc != "Uncharacterized protein") %>%
#   collect() %>% 
#   arrange(keywords, node, org, value, gene_desc) %>% 
#   left_join(aq2ter_go_bp_kw_cands_final %>% select(keywords = parent_terms, prop_go_terms_w_cluster_name = prop)) %>%
#   select(lineage = node, lineage_name = name, cluster_name = keywords, gained_func = value, gene_name, uniprot_acc, org, gene_desc, HOG) %>% 
#   distinct()

# export_tbl %>% 
#   write_tsv("data/TableS6.tsv.gz")

library(treemap)
export_tbl <- read_tsv("data/TableS6.tsv.gz")

stuff <- export_tbl %>% 
  group_by(cluster_name, gained_func) %>% 
  summarize(n_hogs = n_distinct(HOG)) %>% 
  group_by(cluster_name) %>% 
  mutate(prop_hogs = n_hogs/max(n_hogs)) %>% 
  ungroup() %>% 
  filter(prop_hogs > 0.1)

treemap(dtf = stuff, index = c("cluster_name", "gained_func"), vSize = "prop_hogs", type = "categorical", vColor = "cluster_name", inflate.labels = F, lowerbound.cex.labels = 0, bg.labels = "#CCCCCCD9", fontcolor.labels = c("#34449c", "black"), border.lwds = c(1, 0.5), position.legend = "none", force.print.labels = T, title = "", aspRatio = 0.5, overlap.labels = 1, fontfamily.labels = "sans", palette = "Pastel2")

####Get the multicellularity transition nodes####
load("data/20260722_cellularity_anc_reconstruction.simmap")

cellularity_anc_summary <- summary(cell_anc)

ace_matrix <- cellularity_anc_summary$ace
cellularity_transitions <- euk_tree$edge %>% 
  as_tibble() %>% 
  mutate(V1 = as.character(V1), V2 = if_else(V2 <= Ntip(euk_tree), euk_tree$tip.label[V2], as.character(V2))) %>% 
  left_join(ace_matrix %>% as_tibble(rownames = "node"), by = c("V1" = "node"))

p_state <- NULL
for(x in 1:nrow(cellularity_transitions)) {
  row_oi <- cellularity_transitions[x,] %>% pivot_longer(cols = c(-V1, -V2))
  if(any(row_oi$value >= 0.75) | all(row_oi$V1 == 200)) {
    p_state <- rbind(
      p_state,
      tibble(
        V1 = unique(row_oi$V1),
        V2 = unique(row_oi$V2),
        p_state = row_oi$name[row_oi$value == max(row_oi$value)],
        p_pp = max(row_oi$value)
      )
    )
  } else {
    y <- which(p_state$V2 == unique(row_oi$V1))
    parent_row <- p_state[y,] 
    p_state <- rbind(
      p_state,
      tibble(
        V1 = unique(row_oi$V1),
        V2 = unique(row_oi$V2),
        p_state = parent_row$p_state,
        p_pp = row_oi$value[row_oi$name == parent_row$p_state]
      ) 
    )
  }
}

cellularity_transitions <- cellularity_transitions %>% 
  left_join(p_state) %>% 
  select(p_node = V1, V2, p_state, p_pp) %>% 
  left_join(ace_matrix %>% as_tibble(rownames = "node"), by = c("V2" = "node"))

d_state <- NULL
for(x in 1:nrow(cellularity_transitions)) {
  row_oi <- cellularity_transitions[x,] %>% pivot_longer(cols = -c(p_node:p_pp))
  if(any(row_oi$value >= 0.75) | all(row_oi$p_node == 200)) {
    d_state <- rbind(
      d_state,
      tibble(
        p_node = unique(row_oi$p_node),
        V2 = unique(row_oi$V2),
        d_state = row_oi$name[row_oi$value == max(row_oi$value)],
        d_pp = max(row_oi$value)
      )
    )
  } else {
    y <- which(d_state$V2 == unique(row_oi$p_node))
    parent_row <- d_state[y,] 
    d_state <- rbind(
      d_state,
      tibble(
        p_node = unique(row_oi$p_node),
        V2 = unique(row_oi$V2),
        d_state = parent_row$d_state,
        d_pp = row_oi$value[row_oi$name == parent_row$d_state]
      ) 
    )
  }
}

cellularity_transitions <- cellularity_transitions %>% 
  left_join(d_state) %>% 
  select(p_node, d_node = V2, p_state, p_pp, d_state, d_pp) %>% 
  mutate(change = if_else(d_state != p_state, 1, 0))

uni2multi_nodes <- cellularity_transitions %>% filter(change == 1, p_state == "Unicellular", d_state == "Multicellular") %>% rowwise() %>% filter(length(getDescendants(euk_tree, d_node)) > 2) %>% pull(d_node) %>% as.numeric() %>% na.omit() %>% as.vector()

uni2multi_nodes_orig <- gsub("Node", "", euk_tree$node.label[uni2multi_nodes - Ntip(euk_tree)])
uni2multi_nodes_parents <- gsub("Node", "", euk_tree$node.label[sapply(uni2multi_nodes, function(x) getParent(euk_tree, x)) - Ntip(euk_tree)])

#####Analyze GO-CC for multicellular nodes####

uni2multi_all_parents <- list()
for(x in uni2multi_nodes_orig) {
  parent_oi <- getParent(species_tree, as.numeric(x))
  uni2multi_all_parents[[x]] <- c(uni2multi_all_parents[[x]], parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(species_tree, parent_oi)
    uni2multi_all_parents[[x]] <- c(uni2multi_all_parents[[x]], parent_oi)
  }
}


uni2multi_go_cc_term <- tbl(con, "og_gains_losses_parsimony") %>% 
  pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
  mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
  filter(node %in% !!c(uni2multi_nodes_orig, unique(unlist(uni2multi_all_parents)))) %>% 
  left_join(tbl(con, "ptns_og_phylo") %>% select(HOG, org, uniprot_acc:gene_desc), by = join_by(HOG)) %>% 
  inner_join({
    tbl(con, "orthologs_gains_losses_parsimony") %>% 
      pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
      mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
      filter(node %in% !!c(uni2multi_nodes_orig, unique(unlist(uni2multi_all_parents)))) %>% 
      distinct(uniprot_acc, node)
  }) %>% 
  inner_join(tbl(con, "ptns_go_cc") %>% select(uniprot_acc, cc) %>% dplyr::rename(value = cc)) %>% 
  distinct() %>% 
  collect() %>% 
  group_by(node) %>% 
  filter(
    case_when(
      all(node %in% uni2multi_nodes_orig) ~ (
        org %in% (
          euk_cell %>% 
            filter(OSCODE %in% sptree_revised$tip.label[getDescendants(sptree_revised, node = unique(node))]) %>% 
            filter(cellularity == "Multicellular") %>% 
            pull(OSCODE)
        )
      ),
      T ~ T
    )
  ) %>% 
  ungroup()

write_tsv(uni2multi_go_cc_term, "data/20260824_uni2multi_cc_terms_by_group.tsv.gz")

source("scripts/fishers_test_fast.R")
infile = "data/20260824_uni2multi_cc_terms_by_group.tsv.gz"
outfile = "data/20260824_uni2multi_fisher_test_results_cc.csv.gz"

require(data.table)
data <- fread(infile)

# Step 2: Create mappings of nodes to orthogroups and GO terms to orthogroups
node_to_orthogroups <- data %>% distinct(HOG, node) %>% with(., split(HOG, node))
go_to_orthogroups <- data %>% distinct(HOG, value) %>% with(., split(HOG, value))
node_to_go <- data %>% distinct(value, node) %>% with(., split(value, node))
cat("Mappings created...\n")

# Get all unique orthogroups
all_orthogroups <- unique(data$HOG)
run_fast_fishers(infile = infile, outfile = outfile)

uni2multi_go_cc_term_test <- read_csv("data/20260824_uni2multi_fisher_test_results_cc.csv.gz", col_names = c(
  "node",
  "value",
  "has_term_at_node",
  "not_term_at_node",
  "has_term_not_at_node",
  "not_term_not_at_node",
  "odds_ratio",
  "p_value"
), col_types = c("node" = "c"), skip = 1)
uni2multi_go_cc_term_test <- uni2multi_go_cc_term_test %>% 
  mutate(padj = p.adjust(p_value, method = "BH"))

uni2multi_go_cc_term_sig <- uni2multi_go_cc_term_test %>% 
  filter(p_value < (0.05), has_term_not_at_node > 0, has_term_at_node + has_term_not_at_node > 4, odds_ratio > 1)

gained_func_term_tbl <- NULL
for(i in 1:length(uni2multi_nodes_orig)) {
  desc_func <- uni2multi_go_cc_term_sig %>% 
    filter(node == uni2multi_nodes_orig[i]) %>% 
    pull(value) %>% 
    unique()
  
  row_add <- tibble(node = uni2multi_nodes_orig[i], name = major_clades$name[major_clades$node_name == node], gained_func = desc_func)
  
  gained_func_term_tbl <- rbind(gained_func_term_tbl, row_add)
}

gained_func_term_tbl <- gained_func_term_tbl %>% 
  inner_join(uni2multi_go_cc_term_test, by = join_by(node, gained_func == value)) %>% 
  unnest()

annotation_tbl <- tbl(con, "orthologs_gains_losses_parsimony") %>%
  pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>%
  mutate(node = sql("unnest(string_to_array(node, ', '))")) %>%
  filter(node %in% !!c(uni2multi_nodes_orig)) %>%
  distinct(uniprot_acc) %>%
  inner_join(
    (tbl(con, "ptns_go_cc") %>%
       select(uniprot_acc, go_id, term = cc) %>%
       filter(term %in% !!gained_func_term_tbl$gained_func))
  ) %>%
  distinct(uniprot_acc, go_id, term) %>%
  collect()

annotation_list <- annotation_tbl %>%
  with(., split(uniprot_acc, factor(go_id, levels = unique(go_id))))
dag <- import_obo("/net/home.isilon/ds-russell/GO/go.obo", relation_type = c("is_a", "part_of"), annotation = annotation_list)
all_terms <- unique(annotation_tbl$go_id)
all_terms_go <- all_terms[!is.na(all_terms)]
all_terms_sim <- term_sim(dag = dag, terms = all_terms_go, method = "Sim_Lin_1998")
require(simplifyEnrichment)
all_terms_sim_clustered <- cluster_terms(all_terms_sim, method = "louvain", control = list(resolution = 1))
all_terms_IC <- term_IC(dag = dag, method = "IC_offspring")

all_terms_sim_clustered_dist <- annotation_tbl %>%
  filter(!is.na(go_id)) %>%
  mutate(go_cluster = all_terms_sim_clustered[match(go_id, rownames(all_terms_sim))]) %>%
  distinct(term, go_id, go_cluster) %>%
  drop_na(go_cluster) %>%
  group_by(go_cluster) %>%
  filter(n() > 5) %>% 
  do({
    sub_tbl <- .
    sub_tbl %>%
      rowwise() %>%
      do({
        sub_sub_tbl <- .
        tibble(term = sub_sub_tbl$term, go_id = sub_sub_tbl$go_id, go_cluster = sub_sub_tbl$go_cluster, avg_dist_to_others = mean(all_terms_sim[sub_sub_tbl$go_id, sub_tbl$go_id[sub_tbl$go_id != sub_sub_tbl$go_id]]))
      })
  })

gained_func_term_tbl2 <- gained_func_term_tbl %>% 
  inner_join(all_terms_sim_clustered_dist, join_by(gained_func == term)) %>% 
  group_by(go_cluster) %>% 
  ungroup()

parents_func_term_tbl <- NULL
for(i in 1:length(uni2multi_nodes_orig)) {
  parents1 <- getParent(species_tree, uni2multi_nodes_orig[i])
  parents_oi <- parents1
  while(!is.null(parents1)) {
    parents1 <- getParent(species_tree, parents_oi[length(parents_oi)])
    parents_oi <- c(parents_oi, parents1)
  }
  parent_func <- uni2multi_go_cc_term_test %>% 
    filter(node %in% parents_oi, value %in% gained_func_term_tbl2$gained_func) %>% 
    dplyr::rename(parent_node = node, parent_func = value) %>% 
    mutate(node = uni2multi_nodes_orig[i], .before = everything())
  
  row_add <- parent_func
  
  parents_func_term_tbl <- rbind(parents_func_term_tbl, row_add)
}

parents_func_term_tbl2 <- parents_func_term_tbl %>% 
  inner_join(all_terms_sim_clustered_dist, join_by(parent_func == term)) %>% 
  filter(p_value < 0.05)

col_palette <- setNames(rep(RColorBrewer::brewer.pal(n = length(uni2multi_nodes_orig), name = "Dark2"), 1) , c(uni2multi_nodes_orig))

uni2multi_go_cc_odds_tbl_kw_cands <- gained_func_term_tbl2 %>% 
  distinct(go_cluster, gained_func) %>% 
  rowwise() %>% 
  mutate(parent_terms = list(full_go_lineage(term = gained_func, go_obo = go_obo))) %>% 
  unnest() %>% 
  group_by(go_cluster) %>%
  mutate(n_go = n_distinct(gained_func)) %>%
  ungroup() %>% 
  group_by(go_cluster, parent_terms) %>% 
  summarize(n_pt = n(), n_go = unique(n_go), prop = n_pt/n_go, .groups = "drop") %>% 
  arrange(go_cluster, desc(prop)) %>% 
  mutate(IC = all_terms_IC[go_obo$id[match(parent_terms, go_obo$name)]])

uni2multi_go_cc_kw_cands_final <- uni2multi_go_cc_odds_tbl_kw_cands %>% 
  filter(prop > 0.5) %>% 
  group_by(go_cluster) %>%
  mutate(rel_IC = IC/max(IC)) %>% 
  do({
    sub_tbl <- .
    
    print(sub_tbl)
    sub_tbl %>% 
      slice(as.numeric(readline("Choose the row to keep: ")))
  })

go_cc_odds_tbl <- gained_func_term_tbl2 %>% 
  group_by(go_cluster) %>% 
  filter(n_distinct(gained_func) > 2) %>% 
  ungroup() %>% 
  inner_join(uni2multi_go_cc_kw_cands_final %>% select(go_cluster, keywords = parent_terms, prop)) %>% 
  left_join(
    parents_func_term_tbl2 %>% 
      group_by(go_cluster) %>% 
      filter(n_distinct(parent_func) > 2) %>% 
      ungroup() %>% 
      inner_join(uni2multi_go_cc_kw_cands_final %>% select(go_cluster, keywords = parent_terms)),
    join_by(node, go_id, go_cluster, avg_dist_to_others, keywords),
    suffix = c("_gained", "_parent")
  )

f4c1 <- gained_func_term_tbl2 %>% 
  inner_join(uni2multi_go_cc_kw_cands_final %>% select(go_cluster, keywords = parent_terms, prop)) %>% 
  group_by(node, name, keywords, prop) %>% 
  summarize(mean_log2_odds_ratio = mean(log2(odds_ratio), na.rm = T), .groups = "drop") %>% 
  group_by(keywords) %>% 
  mutate(sum_or = sum(mean_log2_odds_ratio, na.rm = T), label = str_wrap(paste0(keywords, "(", round(prop*100, 1), "%)"), 20)) %>% 
  ungroup() %>% 
  ggplot(aes(x = name, y = fct_reorder(label, sum_or), fill = mean_log2_odds_ratio)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(1, NA)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Transition node", y = "MICA (% of terms parent to)", fill = "Mean Log2 Odds Ratio")

sf7b <- gained_func_term_tbl2 %>% 
  inner_join(uni2multi_go_cc_kw_cands_final %>% select(go_cluster, keywords = parent_terms)) %>% 
  group_by(gained_func) %>%
  filter(n_distinct(name) > 1) %>% 
  mutate(sum_or = sum(log2(odds_ratio), na.rm = T)) %>% 
  ungroup() %>% 
  # summarize(mean_log2_odds_ratio = mean(log2(odds_ratio), na.rm = T), .groups = "drop") %>% 
  ggplot(aes(x = name, y = fct_reorder(gained_func, sum_or), fill = log2(odds_ratio))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Transition node", y = "GO term", fill = "Log2 Odds Ratio")

# ptns_gained_at_trans <- tbl(con, "orthologs_gains_losses_parsimony") %>% 
#   pivot_longer(cols = c(gain_nodes, gain_tips), names_to = "gain_type", values_to = "node") %>% 
#   mutate(node = sql("unnest(string_to_array(node, ', '))")) %>% 
#   filter(node %in% !!c(uni2multi_nodes_orig)) %>% 
#   distinct(uniprot_acc) %>% 
#   pull()
# 
# terms_oi <- go_cc_odds_tbl %>% 
#   distinct(node, name, gained_func, keywords)

uni2multi_go_cc_term <-  read_tsv("main_set/analysis/revised_tree/20260824_uni2multi_cc_terms_by_group.tsv.gz")
# ptn_oi <- uni2multi_go_cc_term %>% 
#   filter(node %in% uni2multi_nodes_orig, uniprot_acc %in% ptns_gained_at_trans) %>% 
#   inner_join(terms_oi %>% mutate(node = as.numeric(node)), join_by(node, value == gained_func)) %>% 
#   distinct(node, name, uniprot_acc, HOG, value, keywords)

# export_tbl <- tbl(con, "ptns_og_phylo") %>% 
#   inner_join(ptn_oi, join_by(uniprot_acc, HOG), copy = T) %>% 
#   filter(gene_desc != "Uncharacterized protein") %>% 
#   collect() %>% 
#   # left_join(all_terms_sim_clustered_dist %>% select(term, avg_dist_to_others), join_by(value == term)) %>% 
#   arrange(keywords, node, org, value, gene_desc) %>% 
#   left_join(uni2multi_go_cc_kw_cands_final %>% select(keywords = parent_terms, prop_go_terms_w_cluster_name = prop)) %>% 
#   select(lineage = node, lineage_name = name, cluster_name = keywords, prop_go_terms_w_cluster_name, gained_func = value, gene_name, uniprot_acc, org, gene_desc, HOG) %>% 
#   distinct()
# 
# export_tbl %>% 
#   write_tsv("data/TableS8.tsv.gz")

library(treemap)

export_tbl <- read_tsv("data/TableS8.tsv.gz")

stuff <- export_tbl %>% 
  group_by(cluster_name, gained_func) %>% 
  summarize(n_hogs = n_distinct(HOG)) %>% 
  group_by(cluster_name) %>% 
  mutate(prop_hogs = n_hogs/sum(n_hogs)) %>% 
  ungroup() %>% 
  filter(prop_hogs > 0.01)

treemap(dtf = stuff, index = c("cluster_name", "gained_func"), vSize = "prop_hogs", type = "categorical", vColor = "cluster_name", inflate.labels = F, lowerbound.cex.labels = 0, bg.labels = "#CCCCCCD9", fontcolor.labels = c("#34449c", "black"), border.lwds = c(1, 0.5), position.legend = "none", force.print.labels = T, title = "", aspRatio = 0.5, overlap.labels = 1, fontfamily.labels = "sans", palette = "Pastel2")
