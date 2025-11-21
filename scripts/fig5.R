###Libraries and Database####
library(phytools)
library(tidyverse)
library(DBI)
library(RPostgres)
library(dbplyr)

###First restore the database available here - https://russelllab.org/funcevol/ using the pg_restore command
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_pub", bigint = "integer")

species_tree <- read.tree("data/species_tree_cleaned_final.nwk")
tt_dendo <- chronos(species_tree)
species_details <- read_tsv("data/TableS1.tsv")
major_clades <- read_tsv("data/TableS2.tsv")
prop_orthologs_present <- read_tsv("data/genomic_expectation.tsv.gz")

euk_tree <- keep.tip(read.newick("data/species_tree_for_Count.nwk"), species_details$tree_tip_label[species_details$superkingdom == "Eukaryota"])

source("scripts/fig5_functions.R")

####Get the terrestrial transition nodes####
load("data/20250514_habitat_anc_reconstruction.simmap")

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

###Parents, node itself and descendants
aq2ter_all_parents <- list()
for(x in aq2ter_nodes_orig) {
  parent_oi <- getParent(species_tree, as.numeric(x))
  aq2ter_all_parents[[x]] <- c(aq2ter_all_parents[[x]], parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(species_tree, parent_oi)
    aq2ter_all_parents[[x]] <- c(aq2ter_all_parents[[x]], parent_oi)
  }
}

aq2ter_nodes_full <- sort(unique(c(
  aq2ter_nodes_orig,
  unlist(sapply(aq2ter_nodes_orig, function(x) {
    desc_oi <- getDescendants(species_tree, node = x)
    tips_oi <- species_tree$tip.label[desc_oi]
    desc_oi[!is.na(tips_oi)] <- tips_oi[!is.na(tips_oi)]
    out_nodes <- habitat_transitions_orig %>% filter(d_node %in% desc_oi, d_state == "Terrestrial") %>% pull(d_node)
    out_tips <- match(out_nodes, species_tree$tip.label)
    out_nodes[!is.na(out_tips)] <- out_tips[!is.na(out_tips)]
    out_nodes[out_nodes > 508]
  }))
)))

aq2ter_all_desc <- lapply(aq2ter_nodes_orig, function(x){
  xx <- getDescendants(species_tree, x)
  xx[(xx > Ntip(species_tree) & xx %in% as.numeric(aq2ter_nodes_full)) | (xx <= Ntip(species_tree) & xx %in% which(species_tree$tip.label %in% habitat_transitions$d_node[habitat_transitions$d_state == "Terrestrial"]))]
  
})
names(aq2ter_all_desc) <- aq2ter_nodes_orig

#####Calculate relative node positions for all terrestrial lineages####

aq2ter_nodes_posn_tbl <- rbind(
  tibble(
    lineage = rep(names(aq2ter_all_parents), sapply(aq2ter_all_parents, length)),
    nodes = unlist(lapply(aq2ter_all_parents, rev)),
    type = rep("Parent nodes", length(unlist(aq2ter_all_parents)))
  ),
  tibble(
    lineage = aq2ter_nodes_orig,
    nodes = aq2ter_nodes_orig,
    type = "Transition node"
  ),
  tibble(
    lineage = rep(names(aq2ter_all_desc), sapply(aq2ter_all_desc, length)),
    nodes = unlist(aq2ter_all_desc),
    type = rep("Descendant nodes", length(unlist(aq2ter_all_desc)))
  )
) %>% 
  rowwise() %>% 
  mutate(node_height = nodeheight(tt_dendo, nodes)) %>% 
  ungroup() %>% 
  group_by(lineage) %>% 
  arrange(node_height, .by_group = T) %>% 
  ungroup()

aq2ter_midpoints <- aq2ter_nodes_posn_tbl %>% 
  group_by(lineage) %>% 
  summarize(
    midpoint = node_height[type == "Transition node"],
    min_val = min(node_height),
    max_val = max(node_height),
    .groups = "drop"
  )

aq2ter_nodes_posn_tbl <- aq2ter_nodes_posn_tbl %>% 
  left_join(aq2ter_midpoints) %>% 
  rowwise() %>% 
  mutate(
    rel_posn = case_when(
      node_height <= midpoint ~ node_height/midpoint/2,
      node_height > midpoint ~ scales::rescale(node_height, to = c(0.5, 1), from = c(midpoint, max_val)))
  ) %>% 
  ungroup() %>% 
  select(-midpoint, -min_val, -max_val)

#####Terrestrial known genes evolutionary history####

genes_oi <- list(
  "KEAP1/NRF-2" = c("Q14145", "Q16236"), 
  "Hox Genes" = readLines("data/homeobox_ptn_acc"),
  "Aquaporins" = readLines("data/aquaporins.txt"),
  "Noggin" = "Q13253",
  "BMP-4" = "P12644",
  "Eyeless/Pax6" = c("Q05201", "P63015"),
  "Pax2" = "P32114",
  "GRAS TF" = c("A7U4T7", "Q9M384", "Q9SZF7", "Q9LRW3", "Q9ZWC5", "Q9SUF5", "Q9M000", "Q9LDL7", "Q9CAN3", "G7L166", "G7JMM0", "Q3EDH0")
)

aq2ter_genes_oi_evol_history <- NULL
for(i in 1:length(genes_oi)) {
  print(i)
  histories_oi <- list()
  for(x in genes_oi[[i]]) {
    print(x)
    histories_oi[[x]] <- orthologs_anc_reconst(gene = x, method = "mp") %>% mutate(gene = x)
  } 
  histories_oi <- histories_oi %>% 
    bind_rows() %>% 
    filter(change < 2) %>% 
    mutate(category = names(genes_oi)[i])
  aq2ter_genes_oi_evol_history <- rbind(aq2ter_genes_oi_evol_history, histories_oi)
}

aq2ter_genes_og <- tbl(con, "ptns_og_phylo") %>% 
  filter(uniprot_acc %in% !!unique(aq2ter_genes_oi_evol_history$gene)) %>% 
  distinct(uniprot_acc, HOG) %>% 
  collect()

aq2ter_genes_og_oi <- aq2ter_genes_oi_evol_history %>% 
  inner_join(aq2ter_genes_og, join_by(gene == uniprot_acc)) %>% 
  group_by(gene, category) %>% 
  summarize(HOG = toString(unique(sort(HOG))))


aq2ter_og_oi_evol_history <- aq2ter_genes_og_oi %>% 
  group_by(HOG) %>% 
  do({
    og_oi <- unique(.$HOG)
    print(og_oi)
    og_anc_reconst(og = og_oi)
  })

aq2ter_nodes_posn_category_props_og <- aq2ter_nodes_posn_tbl %>% 
  group_by(lineage, nodes) %>% 
  summarize(rel_posn = min(rel_posn), .groups = "drop") %>% 
  left_join((aq2ter_genes_oi_evol_history %>% inner_join(aq2ter_genes_og, join_by(gene == uniprot_acc)) %>% distinct(category, HOG) %>% left_join(aq2ter_og_oi_evol_history) %>% mutate(node = as.character(node))), join_by(nodes == node)) %>% 
  group_by(category, HOG) %>% 
  mutate(earliest_gain = if_else(gain == 1 & rel_posn == min(rel_posn[gain == 1]), 1, 0)) %>% 
  ungroup()

aq2ter_nodes_posn_category_props <- aq2ter_nodes_posn_tbl %>% 
  left_join(aq2ter_genes_oi_evol_history %>% mutate(node = as.character(node)), join_by(nodes == node)) %>% 
  left_join(major_clades %>% transmute(lineage = as.character(node), name)) %>% 
  group_by(lineage, name, nodes, rel_posn, category) %>% 
  summarize(prop_present = sum(State == "Present")/n_distinct(gene), prop_absent = sum(State == "Absent")/n_distinct(gene), n_genes = n_distinct(gene), .groups = "drop") %>% 
  mutate(non_zeros = n_distinct(lineage[prop_present > 0.25]), .by = c(category)) %>% 
  mutate(category = fct_reorder(category, non_zeros)) %>%
  left_join(prop_orthologs_present %>% mutate(nodes = as.character(node)) %>% select(nodes, prop_expectation = prop_present), join_by(nodes)) %>%
  mutate(
    # observed successes = how many genes in this category originate at this node
    k_obs = round(prop_present * n_genes),
    
    # expected successes under null
    k_exp = prop_expectation * n_genes,
    
    # enrichment ratio just for plotting / effect size
    enrichment_ratio = prop_present / prop_expectation
  ) %>% 
  rowwise() %>%
  mutate(
    # run one-sided binomial test: P(X >= k_obs | n_set, p = prop_expectation)
    binom_pval = {
      bt <- binom.test(
        x = k_obs,
        n = n_genes,
        p = prop_expectation,
        alternative = "greater"
      )
      bt$p.value
    },
    
    # optional: observed proportion CI from binomial test
    ci_low  = binom.test(k_obs, n_genes, p = prop_expectation, alternative = "greater")$conf.int[1],
    ci_high = binom.test(k_obs, n_genes, p = prop_expectation, alternative = "greater")$conf.int[2]
  ) %>%
  ungroup() %>%
  mutate(
    binom_padj = p.adjust(binom_pval, method = "BH")
  ) %>% 
  mutate(pt_shape = if_else(binom_padj < 0.05, "*", "NS")) %>% 
  mutate(cat_label = paste0(category, "\n", n_genes, " gene(s)\n", "XXXXX"))

col_palette <- setNames(rep(RColorBrewer::brewer.pal(n = length(aq2ter_nodes_orig), name = "Dark2"), 1) , c(aq2ter_nodes_orig))

aq2ter_nodes_posn_category_props$col_cat <- col_palette[aq2ter_nodes_posn_category_props$lineage]
aq2ter_nodes_posn_category_props_og$col_cat <- col_palette[aq2ter_nodes_posn_category_props_og$lineage]

aq2ter_nodes_posn_category_props_tbl <- aq2ter_nodes_posn_category_props %>% 
  nest_by(lineage, category, col_cat, rel_posn <= 0.25, rel_posn > 0.25 & rel_posn <= 0.5, rel_posn > 0.5 & rel_posn <= 0.75, rel_posn > 0.75 & rel_posn <= 1) %>% 
  mutate(mod = list(lm(prop_present~rel_posn, data = data))) %>% 
  reframe(broom::augment(mod)) %>% 
  mutate(.fitted = if_else(.fitted < 0, 0, if_else(.fitted > 1, 1, .fitted))) %>% 
  left_join(aq2ter_nodes_posn_category_props) %>% 
  mutate(non_zeros = n_distinct(lineage[prop_present > 0.25]), .by = c(category)) %>% 
  mutate(category = fct_reorder(category, non_zeros))

f5a <- aq2ter_nodes_posn_category_props %>% 
  ggplot(aes(x = rel_posn, y = prop_present*100, colour = col_cat)) +
  geom_vline(aes(xintercept = 0.5), colour = "darkblue", linetype = "dashed") +
  geom_line(data = prop_orthologs_present %>% left_join(aq2ter_nodes_posn_tbl %>% transmute(node = as.numeric(nodes), rel_posn)), mapping = aes(x = rel_posn, y = prop_present*100), colour = "darkred", linetype = "dotted", linewidth = 0.5) +
  geom_line(data = aq2ter_nodes_posn_category_props_tbl %>% filter(!category %in% c("Type I ORs", "Type II ORs")), mapping = aes(x = rel_posn, y = .fitted*100, colour = col_cat), linewidth = 1.25, inherit.aes = F) +
  geom_point(aes(shape = pt_shape), size = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = rel_posn), colour = "gray40", linewidth = 2.5, alpha = 0.5, data = (aq2ter_nodes_posn_category_props_og  %>% filter(!category %in% c("Type I ORs", "Type II ORs")) %>% filter(earliest_gain > 0) %>% distinct(rel_posn, category, earliest_gain)), inherit.aes = F) +
  facet_wrap(category~., ncol = 2, scales = "free_y", labeller = labeller(category = setNames(aq2ter_nodes_posn_category_props$cat_label, aq2ter_nodes_posn_category_props$category))) +
  scale_colour_identity(labels = (major_clades %>% dplyr::slice(match(aq2ter_nodes_orig, node)) %>% pull(name) %>% as.character()), breaks = as.character(col_palette[1:length(aq2ter_nodes_orig)]), guide = "legend") +
  scale_shape_manual(values = c("*" = 8, "NS" = 19)) +
  ylim(c(0, 100)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom") +
  labs(
    x = "Relative node position (Root --> Tips)",
    y = "% present",
    colour = "Taxonomy of transition node"
  )

f5a

####Get the multicellularity transition nodes####
load("data/20250428_cellularity_anc_reconstruction.simmap")

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

uni2multi_nodes_full <- sort(unique(c(
  uni2multi_nodes_orig,
  unlist(sapply(uni2multi_nodes_orig, function(x) {
    desc_oi <- getDescendants(species_tree, node = x)
    tips_oi <- species_tree$tip.label[desc_oi]
    desc_oi[!is.na(tips_oi)] <- tips_oi[!is.na(tips_oi)]
    out_nodes <- cellularity_transitions_orig %>% filter(d_node %in% desc_oi, d_state == "Multicellular") %>% pull(d_node)
    out_tips <- match(out_nodes, species_tree$tip.label)
    out_nodes[!is.na(out_tips)] <- out_tips[!is.na(out_tips)]
    out_nodes[out_nodes > Ntip(species_tree)]
  }))
)))

all_parents <- list()
for(x in uni2multi_nodes_orig) {
  parent_oi <- getParent(species_tree, as.numeric(x))
  all_parents[[x]] <- c(all_parents[[x]], parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(species_tree, parent_oi)
    all_parents[[x]] <- c(all_parents[[x]], parent_oi)
  }
}

all_desc <- lapply(uni2multi_nodes_orig, function(x){
  xx <- getDescendants(species_tree, x)
  xx[(xx > Ntip(species_tree) & xx %in% as.numeric(aq2ter_nodes_full)) | (xx <= Ntip(species_tree) & xx %in% which(species_tree$tip.label %in% cellularity_transitions$d_node[cellularity_transitions$d_state == "Multicellular"]))]
})
names(all_desc) <- uni2multi_nodes_orig


#####Calculate relative node positions for all multicellular lineages####
cellularity_transitions_orig <- cellularity_transitions %>% 
  mutate(p_node = gsub("Node", "", euk_tree$node.label[as.numeric(p_node) - Ntip(euk_tree)])) %>% 
  mutate(d_node = if_else(!is.na(as.numeric(d_node)), gsub("Node", "", euk_tree$node.label[as.numeric(d_node) - Ntip(euk_tree)]), d_node))

nodes_posn_tbl <- rbind(
  tibble(
    lineage = rep(names(all_parents), sapply(all_parents, length)),
    nodes = unlist(lapply(all_parents, rev)),
    type = rep("Parent nodes", length(unlist(all_parents)))
  ) %>% 
    mutate(posn = row_number(), .by = lineage),
  tibble(
    lineage = uni2multi_nodes_orig,
    nodes = uni2multi_nodes_orig,
    type = "Transition node",
    posn = sapply(all_parents, length) + 1
  ),
  tibble(
    lineage = rep(names(all_desc), sapply(all_desc, length)),
    nodes = unlist(all_desc),
    type = rep("Descendant nodes", length(unlist(all_desc)))
  ) %>% 
    mutate(rn = row_number(), .by = lineage) %>% 
    mutate(posn = sapply(all_parents, length)[lineage] + 1 + rn, .by = lineage) %>% 
    select(-rn)
) %>% 
  rowwise() %>% 
  mutate(node_height = nodeheight(tt_dendo, nodes)) %>% 
  ungroup() %>% 
  group_by(lineage) %>% 
  arrange(node_height, .by_group = T) %>% 
  ungroup()

midpoints <- nodes_posn_tbl %>% 
  group_by(lineage) %>% 
  summarize(
    midpoint = node_height[type == "Transition node"],
    min_val = min(node_height),
    max_val = max(node_height),
    .groups = "drop"
  )

nodes_posn_tbl <- nodes_posn_tbl %>% 
  left_join(midpoints) %>% 
  rowwise() %>% 
  mutate(
    rel_posn = case_when(
      node_height <= midpoint ~ node_height/midpoint/2,
      node_height > midpoint ~ scales::rescale(node_height, to = c(0.5, 1), from = c(midpoint, max_val)))
  ) %>% 
  ungroup() %>% 
  select(-midpoint, -min_val, -max_val)

#####Multicellularity known genes evolutionary history####

genes_oi <- list(
  "Cadherins" = tbl(con, "pfam_domains_gains_losses_parsimony") %>% filter(grepl("Cadherin$|[C]adherin_\\d{1}|^GPS$|7tm_2|Protocadherin", pfam_name)) %>% distinct(pfam_id) %>% pull(),
  "Animal LCA genes" = c(
    tbl(con, "ptns_og_phylo") %>% filter(grepl("^Protein Wnt|^Transforming growth factor beta$|^Transforming growth factor beta-\\d.proprotein$|^Frizzled$|^Tyrosine-protein kinase JAK|^Segment polarity protein dishevelled|^Mothers against decapentaplegic|^Catenin beta|^Catenin delta|^Espin$|^Dystroglycan|^Hemicentin|^Fermitin family", gene_desc)) %>% filter(!grepl(".protein|.isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull()
  ),
  "Integrin adhesion complex" = c(
    tbl(con, "ptns_og_phylo") %>% filter(grepl("[I]ntegrin.[Bb]eta|[I]ntegrin.[Aa]lpha|Focal adhesion kinase|^Paxillin|^Talin|^Alpha-actinin|^Vinculin|^Zyxin|^KN motif and ankyrin repeat domain|^Vasodilator-stimulated phosphoprotein|^Tensin", gene_desc)) %>% filter(!grepl(".protein|.isoform|.like|.containing|Fragment|^Putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull()
  ),
  "Flocculation proteins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("[Ff]locculation.protein", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Yeast clonal" = c("P29029", "P21192"),
  "Chlamydomonas cyclins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("[Cc]yclin|CYCLIN", gene_desc), org == "CHLRE") %>% distinct(uniprot_acc) %>% pull(),
  "NMD machinery and SR splicing factors" = tbl(con, "ptns_og_phylo") %>% distinct(uniprot_acc, org, gene_name, gene_desc) %>% filter(org == "HUMAN") %>% filter(gene_name %in% c("UPF1", "UPF2", "UPF3A", "UPF3B", "SMG1", "SMG5", "SMG6", "SMG7", "SMG8", "SMG9", "ETF1", "GSPT1", "GSPT2", "EIF4A3", "MAGOH", "MAGOHB", "PYM1", "CASC3")) %>% pull(uniprot_acc),
  "Fungal multicellularity genes" = tbl(con, "ptns_og_phylo") %>% distinct(uniprot_acc, org, gene_name, gene_desc) %>% filter(org == "YEAST") %>% filter(gene_name %in% !!readLines("data/multicellular_fungi_gene_list.txt")) %>% pull(uniprot_acc)
)

animal_lca_genes_split <- list(
  "Wnt" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Protein Wnt", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "TGF-beta" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Transforming growth factor beta$|^Transforming growth factor beta-\\d.proprotein$", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Frizzled" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Frizzled$", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "JAK" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Tyrosine-protein kinase JAK", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Dishevelled" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Segment polarity protein dishevelled", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "SMADs" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Mothers against decapentaplegic", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "beta-Catenins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Catenin beta", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "delta-Catenins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Catenin delta", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Espin" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Espin$", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Dystroglycans" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Dystroglycan", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Hemicentins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Hemicentin", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull(),
  "Fermitins" = tbl(con, "ptns_og_phylo") %>% filter(grepl("^Fermitin family", gene_desc)) %>% filter(!grepl(".isoform|.like|.containing|Fragment|^Putative|putative", gene_desc)) %>% distinct(uniprot_acc) %>% pull()
) %>% 
  enframe(name = "category_lca", value = "gene") %>% 
  unnest()

genes_oi_evol_history <- NULL
for(i in 1:length(genes_oi)) {
  print(i)
  histories_oi <- list()
  for(x in genes_oi[[i]]) {
    print(x)
    if(grepl("^PF", x)) {
      histories_oi[[x]] <- pfam_node_states(pfam_oi = x) %>% mutate(gene = x)
    } else if(grepl("^N0.", x)) {
      histories_oi[[x]] <- og_anc_reconst(og = x, method = "mp") %>% mutate(gene = x)
    } else {
      histories_oi[[x]] <- orthologs_anc_reconst(gene = x, method = "mp") %>% mutate(gene = x)
    }
  } 
  histories_oi <- histories_oi %>% 
    bind_rows() %>% 
    # filter(node %in% nodes_posn_tbl$nodes, change < 2) %>%
    filter(change < 2) %>% 
    mutate(category = names(genes_oi)[i])
  genes_oi_evol_history <- rbind(genes_oi_evol_history, histories_oi)
}

genes_og <- tbl(con, "ptns_og_phylo") %>% 
  filter(uniprot_acc %in% !!(genes_oi_evol_history %>% left_join(animal_lca_genes_split) %>% filter(is.na(category_lca) | category_lca == "beta-Catenins") %>% mutate(category = if_else(!is.na(category_lca), category_lca, category)) %>% pull(gene))) %>% 
  distinct(uniprot_acc, HOG) %>% 
  collect()

genes_og_oi <- genes_oi_evol_history %>% 
  inner_join(genes_og, join_by(gene == uniprot_acc)) %>% 
  group_by(gene, category) %>% 
  summarize(HOG = toString(unique(sort(HOG))))

genes_og_oi_lineage <- tbl(con, "ptns_og_phylo") %>% 
  filter(uniprot_acc %in% !!unique(genes_og_oi$gene)) %>% 
  collect() %>% 
  inner_join(lapply(all_desc, function(x) species_tree$tip.label[x]) %>% enframe(name = "lineage", value = "org") %>% unnest() %>% drop_na()) %>% 
  distinct(uniprot_acc, HOG, lineage)


og_oi_evol_history <- genes_og_oi %>% 
  group_by(HOG) %>% 
  do({
    og_oi <- unique(.$HOG)
    print(og_oi)
    og_anc_reconst(og = og_oi)
  })

nodes_posn_category_props_og <- nodes_posn_tbl %>% 
  group_by(nodes) %>% 
  summarize(rel_posn = min(rel_posn), .groups = "drop") %>% 
  left_join((genes_oi_evol_history %>% left_join(animal_lca_genes_split) %>% filter(is.na(category_lca) | category_lca == "beta-Catenins") %>% mutate(category = if_else(!is.na(category_lca), category_lca, category)) %>% inner_join(genes_og, join_by(gene == uniprot_acc)) %>% distinct(category, HOG) %>% left_join(og_oi_evol_history) %>% mutate(node = as.character(node))), join_by(nodes == node)) %>% 
  group_by(category, HOG) %>% 
  mutate(earliest_gain = if_else(gain == 1 & rel_posn == min(rel_posn[gain == 1]), 1, 0)) %>% 
  ungroup() %>% 
  left_join(genes_og_oi_lineage %>% distinct(HOG, lineage))

nodes_posn_category_props <- nodes_posn_tbl %>% 
  left_join(genes_oi_evol_history %>% mutate(node = as.character(node)), join_by(nodes == node)) %>% 
  left_join(animal_lca_genes_split) %>% 
  filter(is.na(category_lca) | category_lca == "beta-Catenins") %>% 
  mutate(category = if_else(!is.na(category_lca), category_lca, category)) %>% 
  left_join(major_clades %>% transmute(lineage = as.character(node), name)) %>% 
  group_by(lineage, name, nodes, rel_posn, category) %>% 
  summarize(prop_present = sum(State == "Present")/n(), prop_absent = sum(State == "Absent")/n(), n_genes = n(), .groups = "drop") %>% 
  mutate(non_zeros = n_distinct(lineage[prop_present > 0.2]), .by = c(category)) %>% 
  mutate(category = fct_reorder(category, non_zeros)) %>%
  left_join(prop_orthologs_present %>% mutate(nodes = as.character(node)) %>% select(nodes, prop_expectation = prop_present), join_by(nodes)) %>%
  mutate(
    # observed successes = how many genes in this category originate at this node
    k_obs = round(prop_present * n_genes),
    
    # expected successes under null
    k_exp = prop_expectation * n_genes,
    
    # enrichment ratio just for plotting / effect size
    enrichment_ratio = prop_present / prop_expectation
  ) %>% 
  rowwise() %>%
  mutate(
    # run one-sided binomial test: P(X >= k_obs | n_set, p = prop_expectation)
    binom_pval = {
      bt <- binom.test(
        x = k_obs,
        n = n_genes,
        p = prop_expectation,
        alternative = "greater"
      )
      bt$p.value
    },
    
    # optional: observed proportion CI from binomial test
    ci_low  = binom.test(k_obs, n_genes, p = prop_expectation, alternative = "greater")$conf.int[1],
    ci_high = binom.test(k_obs, n_genes, p = prop_expectation, alternative = "greater")$conf.int[2]
  ) %>%
  ungroup() %>%
  mutate(
    binom_padj = p.adjust(binom_pval, method = "BH")
  ) %>% 
  mutate(pt_shape = if_else(binom_padj < 0.05, "*", "NS")) %>% 
  mutate(cat_label = paste0(category, "\n", n_genes, " gene(s)\n", "XXXXX"))

col_palette <- setNames(rep(RColorBrewer::brewer.pal(n = length(uni2multi_nodes_orig), name = "Dark2"), 1) , c(uni2multi_nodes_orig))
nodes_posn_category_props$col_cat <- col_palette[nodes_posn_category_props$lineage]
nodes_posn_category_props_og$col_cat <- col_palette[nodes_posn_category_props_og$lineage]
nodes_posn_category_props_og <- nodes_posn_category_props_og %>% 
  mutate(col_cat = if_else(is.na(col_cat), "gray40", col_cat))

nodes_posn_category_props_tbl <- nodes_posn_category_props %>% 
  nest_by(lineage, category, col_cat, rel_posn <= 0.25, rel_posn > 0.25 & rel_posn <= 0.5, rel_posn > 0.5 & rel_posn <= 0.75, rel_posn > 0.75 & rel_posn <= 1) %>% 
  mutate(mod = list(lm(prop_present~rel_posn, data = data))) %>% 
  reframe(broom::augment(mod)) %>% 
  mutate(.fitted = if_else(.fitted < 0, 0, if_else(.fitted > 1, 1, .fitted))) %>% 
  left_join(nodes_posn_category_props) %>% 
  mutate(non_zeros = n_distinct(lineage[prop_present > 0]), .by = c(category)) %>% 
  mutate(category = fct_reorder(category, non_zeros))

f5b <- nodes_posn_category_props %>% 
  filter(category != "Chlamydomonas Histones") %>% 
  ggplot(aes(x = rel_posn, y = prop_present*100, colour = col_cat)) +
  geom_vline(aes(xintercept = 0.5), colour = "darkblue", linetype = "dashed") +
  geom_line(data = prop_orthologs_present %>% left_join(aq2ter_nodes_posn_tbl %>% transmute(node = as.numeric(nodes), rel_posn)), mapping = aes(x = rel_posn, y = prop_present*100), colour = "darkred", linetype = "dotted", linewidth = 0.5) +
  geom_line(data = nodes_posn_category_props_tbl, mapping = aes(x = rel_posn, y = .fitted*100, colour = col_cat), linewidth = 1.25, inherit.aes = F) +
  geom_point(aes(shape = pt_shape), size = 2, alpha = 0.3) +
  geom_vline(aes(xintercept = rel_posn), colour = "gray40", linewidth = 2.5, alpha = 0.5, data = (nodes_posn_category_props_og  %>% filter(category != "Chlamydomonas Histones") %>% filter(earliest_gain > 0) %>% distinct(rel_posn, category, col_cat, earliest_gain)), inherit.aes = F) +
  # geom_smooth(method = "gam", se = F) +
  facet_wrap(category~., labeller = labeller(category = setNames(nodes_posn_category_props$cat_label, nodes_posn_category_props$category)), ncol = 2, scales = "free_y") +
  scale_colour_identity(labels = (major_clades %>% dplyr::slice(match(uni2multi_nodes_orig, node)) %>% pull(name) %>% as.character()), breaks = col_palette[1:length(uni2multi_nodes_orig)], guide = "legend") +
  scale_shape_manual(values = c("*" = 8, "NS" = 19)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom") +
  labs(
    x = "Relative node position (Root --> Tips)",
    y = "% present",
    colour = "Taxonomy of transition node"
  )

f5b
