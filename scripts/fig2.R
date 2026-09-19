###Libraries and Database####
library(phytools)
library(tidyverse)
library(DBI)
library(RPostgres)
library(dbplyr)

###First restore the database available here - https://russelllab.org/funcevol/ using the pg_restore command
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_revision", bigint = "integer")

#### Fig 2A-D ####

nodes_list <- list(
  "Eukaryota" = c("548"),
  "Euteleostomi" = c("693"),
  "Deuterostomia" = c("689"),
  "Metazoa + Choanoflagellate" = c("648", "649"),
  "Embryophyta" = c("572"),
  "Magnoliopsida" = c("575"),
  "Halobacteriales" = c("545"),
  "Dikarya" = c("603"),
  "Dictyostelia" = c("595"),
  "Cyanobacteria" = c("950", "953"),
  "Pezizomycotina" = c("629")
)

nodes_sig_terms <- readRDS("data/20260806_go_enrichment_data.rds")

go_obo <- ontologyIndex::get_ontology(file = "data/go.obo", propagate_relationships = c("is_a", "part_of", "regulates", "occurs_in", "in_taxon"))

full_go_lineage <- function(term, go_obo) {
  xx <- unlist(go_obo$name[unlist(go_obo$ancestors[go_obo$id[match(term, go_obo$name)]])])
  xx[-length(xx)]
}

nodes_sig_terms_desc <- lapply(1:length(nodes_sig_terms), function(x) {
  nodes_sig_terms[[x]] %>% 
    mutate(clade = names(nodes_list)[x])
}) %>% 
  bind_rows() %>%
  arrange(-fishers_odds_ratio)

#####treemaps####

library(simona)
library(simplifyEnrichment)
bp_dag <- create_ontology_DAG_from_GO_db(
  namespace = "BP",
  org_db = NULL
)

stuff_list <- NULL
for(i in 1:length(nodes_sig_terms)) {
  print(i)
  go_ids_oi <- unique(go_obo$id[match(nodes_sig_terms[[i]]$bp[1:50], go_obo$name)])
  go_sim <- term_sim(
    bp_dag,
    terms = go_ids_oi,
    method = "Sim_Wang_2007"
  )
  
  go_clusters <- cluster_terms(
    go_sim,
    method = "louvain"
  )
  
  cluster_key <- tibble(
    go_id = rownames(go_sim),
    bp = go_obo$name[match(go_id, go_obo$id)],
    semantic_cluster = as.integer(go_clusters)
  )
  
  stuff <- nodes_sig_terms[[i]] %>% 
    inner_join(cluster_key) %>% 
    group_by(semantic_cluster) %>% 
    arrange(-term_at_node_totals, .by_group = T) %>% 
    mutate(keyword = paste(na.omit(keyword_enrichment_from_GO(go_id = go_id)$keyword[1:2]), collapse = "|\n"))
  
  stuff_list[[i]] <- stuff
  names(stuff_list)[i] <- names(nodes_list)[i]
}

library(treemap)
for(nm in names(stuff_list)) {
  pdf(file = paste0(nm, "_treemap.pdf"), width = 4, height = 4, family = "ArialMT")
  treemap(dtf = stuff_list[[nm]], index = c("keyword", "bp"), vSize = "term_at_node_totals", type = "categorical", vColor = "keyword", inflate.labels = F, lowerbound.cex.labels = 0, bg.labels = "#CCCCCCD9", position.legend = "none", force.print.labels = T, title = "", aspRatio = 1, overlap.labels = 1)
  dev.off()
}

####KEGG gains#####
kegg_hierarchy <- read_tsv("data/20230713_kegg_pathway_hierarchy.tsv.gz")

data <- tbl(con, "kegg_gains_parsimony") %>% 
  inner_join(enframe(nodes_list, name = "clade", value = "node") %>% unnest() %>% mutate(node = as.numeric(node)), copy = T) %>% 
  group_by(clade, pathway_name) %>% 
  summarize(prop_comp_gained = sum(prop_comp_gained), n_comp_gained = sum(n_comp_gained)) %>% 
  collect() %>% 
  left_join(kegg_hierarchy) %>% 
  filter(H1 != "Human Diseases") %>%
  group_by(clade, pathway_name) %>% 
  summarize(avg_prop_gained = mean(prop_comp_gained), total_comp_gained = sum(n_comp_gained)) %>% 
  group_by(clade) %>% 
  arrange(-avg_prop_gained) %>% 
  slice_max(order_by = avg_prop_gained, n = 10) %>% 
  mutate(pathway_name = tidytext::reorder_within(x = pathway_name, by = avg_prop_gained, within = clade))

ggplot(data, aes(y = pathway_name, x = avg_prop_gained*100, fill = total_comp_gained)) +
  geom_bar(stat = "identity") +
  # scale_x_continuous(expand = c(0, NA)) +
  scale_fill_viridis_c(option = "F", limits = c(0, max(data$total_comp_gained))) +
  theme_bw(base_size = 14) +
  facet_wrap(~factor(clade, levels = names(nodes_list)[c(1, 4, 3, 2, 7, 9, 8, 11, 10, 5, 6)]), scales = "free", ncol = 2) +
  tidytext::scale_y_reordered() +
  labs(x = "% Pathway gained", y = "KEGG Pathway", fill = "# of components\ngained")

####Figure 2E####
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_revision", bigint = "integer")

species_tree <- read.tree("data/species_tree_cleaned_final.nwk")
major_clades <- read_tsv("data/TableS2.tsv")

all_parents_finder <- function(x) {
  all_parents <- NULL
  parent_oi <- getParent(sptree_revised, as.numeric(x))
  all_parents <- c(all_parents, parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(sptree_revised, parent_oi)
    all_parents <- c(all_parents, parent_oi)
  }
  toString(unique(all_parents))
}

term_origins <- function(search_term, mode = c("gene", "hog")[1], drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T) {
  if(mode == "gene") {
    og_gains_losses <- tbl(con, "orthologs_gains_losses_parsimony")
  } else if(mode == "hog") {
    og_gains_losses <- tbl(con, "og_gains_losses_parsimony")
  }
  print(paste("Starting with search term:", search_term, "..."))
  t1 <- tbl(con, "ptns_go_bp") %>% 
    mutate(keyword = as.character(search_term)) 
  
  if(drop_iea) {
    t1 <- t1 %>% 
      filter(evidence_code != "IEA")
  }
  
  if(use_full_lineage) {
    t1 <- t1 %>% 
      select(uniprot_acc, bp, keyword, matches("l\\d+")) %>% 
      filter(if_any(.cols = -keyword, .fns = ~ grepl(pattern = keyword, x = .x, ignore.case = TRUE)))
  } else {
    t1 <- t1 %>% 
      filter(grepl(pattern = keyword, x = bp))
  }
  
  t2 <- t1 %>% 
    # distinct(uniprot_acc, bp) %>%
    inner_join(tbl(con, "ptns_og_phylo"))
  
  if(drop_model_orgs) {
    t2 <- t2 %>% 
      filter(!org %in% model_orgs)
  }
  
  sleep_og_hist <- t2 %>%
    inner_join(og_gains_losses) %>% 
    collect()
  
  if(drop_taxon_violations) {
    all_node_desc <- lapply(c((Ntip(sptree_revised) + 1) : (Ntip(sptree_revised) + Nnode(sptree_revised))), function(x) {
      xx <- sptree_revised$tip.label[getDescendants(sptree_revised, x)]
      xx[!is.na(xx)]
    })
    names(all_node_desc) <- c((Ntip(sptree_revised) + 1) : (Ntip(sptree_revised) + Nnode(sptree_revised)))
    all_node_desc <- all_node_desc %>% enframe(name = "node", value = "desc") %>% unnest()
    
    
    prop_candidates <- sleep_og_hist %>%
      select(uniprot_acc, bp, gain_nodes) %>%
      separate_longer_delim(gain_nodes, delim = ", ") %>%
      rename(node = gain_nodes) %>%
      filter(!is.na(node), node != "") %>%
      left_join(all_node_desc, by = "node") %>%
      rename(organism_id = desc) %>% 
      filter(!is.na(organism_id)) %>%
      distinct(uniprot_acc, bp, node, organism_id)
    
    message("Beaming candidate protein histories to the DB for testing...")
    copy_to(
      con,
      prop_candidates,
      name = "tmp_prop_candidates",
      temporary = TRUE,
      overwrite = TRUE,
      indexes = list(
        c("organism_id", "bp")
      )
    )
    
    prop_candidates_db <- tbl(con, "tmp_prop_candidates")
    
    taxon_db <- tbl(con, "go_bp_taxon_constraints") %>%
      select(
        organism_id,
        go_term,
        can_exist
      )
    
    taxon_test_db <- prop_candidates_db %>%
      left_join(
        taxon_db,
        join_by(
          organism_id,
          bp == go_term
        )
      ) %>%
      mutate(
        taxon_allowed = coalesce(can_exist, TRUE)
      )
    message("Testing if GO terms are allowed in the decendant species...")
    taxon_test <- taxon_test_db %>%
      collect()
    
    allowed_taxon_go <- taxon_test %>% filter(taxon_allowed) %>% distinct(uniprot_acc, bp, node)
    
    sleep_og_hist2 <- sleep_og_hist %>%
      separate_longer_delim(gain_nodes, delim = ", ") %>%
      rename(node = gain_nodes) %>%
      filter(!is.na(node), node != "") %>% 
      inner_join(allowed_taxon_go)
    
    dbExecute(con, "DROP TABLE IF EXISTS tmp_prop_candidates;")
    
  } else {
    sleep_og_hist2 <- sleep_og_hist %>%
      separate_longer_delim(gain_nodes, delim = ", ") %>%
      rename(node = gain_nodes) %>%
      filter(!is.na(node), node != "")
  }
  
  print("Counting gains...")
  counts_table <- sleep_og_hist2 %>% 
    # pivot_longer(cols = c(gain_nodes), names_to = "gain_type", values_to = "node") %>% 
    # mutate(node = strsplit(node, split = ", ")) %>% 
    # unnest() %>% 
    # filter(node != "") %>% 
    mutate(node = as.numeric(node)) %>% 
    group_by(node) %>%
    summarize(n = length(unique(HOG))) %>% 
    left_join(major_clades %>% dplyr::select(node = node_name, name, level, full_taxonomy, full_taxonomy_levels, other_names), by = c("node")) %>% 
    filter(as.numeric(node) > Ntip(sptree_revised)) %>%
    arrange(-n)
  
  
  if(nrow(counts_table) > 0) {
    # print(paste("Building taxonomy for", nrow(counts_table), "entries ..."))
    # 
    # counts_table <- counts_table %>%
    #   mutate(
    #     taxa_names = str_split(full_taxonomy, ", ", simplify = FALSE),
    #     taxa_levels = str_split(full_taxonomy_levels, ", ", simplify = FALSE),
    #     phylum_name = map2_chr(taxa_levels, taxa_names, ~if_else("order" %in% .x, .y[which(.x == "order")][1], NA_character_))
    #   ) %>% 
    #   mutate(
    #     phylum = case_when(
    #       level %in% c("kingdom", "superkingdom", "phylum", "subphylum", "class", "order") ~ name,
    #       !is.na(phylum_name) ~ phylum_name,
    #       TRUE ~ name
    #     )
    #   ) %>%
    #   select(-taxa_names, -taxa_levels, -phylum_name)
    
    counts_table %>% 
      mutate(phylum = paste(node, name, sep = "|")) %>% 
      group_by(phylum) %>%
      summarize(n = sum(n)) %>% 
      mutate(
        perc = round(n/sum(n)*100, digits = 2),
        relative_prop = n/max(n), 
        z_score = (n - mean(n))/sd(n), 
        p_value = 2 * pnorm(-abs(z_score)),
        term = as.character(search_term), 
        n_groups = if_else(mode == "gene", paste0(nrow(sleep_og_hist), " genes in ", length(unique(sleep_og_hist$org)), " organisms\n with ", sum(n), " independent gains"), paste0(nrow(sleep_og_hist), " Orthogroups"))
      ) %>% 
      arrange(-n)
  } else {
    tibble(
      name = NA,
      n = 0,
      perc = 0,
      relative_prop = 0,
      term = search_term,
      n_groups = "no hits found"
    )
  }
}

terms_oi_table_new <- rbind(
  term_origins(search_term = "photosynthesis", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "brain", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "nuclear pore", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "adaptive immun", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "heart development", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "flagella|flagellum", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T),
  term_origins(search_term = "CRISPR", mode = "gene", drop_taxon_violations = T, drop_iea = F, drop_model_orgs = F, use_full_lineage = T)
)

terms_wider <- terms_oi_table_new %>% 
  mutate(term = paste0(term, " (", n_groups, ")")) %>% 
  group_by(term) %>% 
  mutate(c_perc = cumsum(perc)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = term, names_from = phylum, values_from = perc) %>% 
  mutate(across(-term, ~ replace_na(., 0)))

col_order_temp <- colnames(terms_wider)[-1][match(unique(major_clades$name[as.numeric(major_clades$node_name)>508]), str_split_i(colnames(terms_wider)[-1], pattern = "[|]", i = 2), nomatch = F)]
col_order <- unique(col_order_temp[!is.na(col_order_temp)])

scale_high <- if(ceiling(max(terms_oi_table$perc)) %% 5 == 0) {
  ceiling(max(terms_oi_table$perc))
} else {
  xx <- ceiling(max(terms_oi_table$perc))
  while(xx %% 5 != 0) {
    xx <- xx + 1
  }
}
breaks_oi <- seq(0, xx, 1)
f3 <- terms_wider %>% 
  pivot_longer(-term) %>% 
  left_join(terms_oi_table_new %>% dplyr::rename(name = phylum) %>% mutate(term = paste0(term, " (", n_groups, ")"), Significant = if_else(p_value < 0.05, "*", "")), by = c("term", "name")) %>% 
  replace_na(list(n = 0)) %>% 
  mutate(perc_2.5 = value >= 2) %>% 
  filter(any(perc_2.5), .by = name) %>% 
  mutate(name = factor(name, levels = col_order)) %>%
  drop_na(name) %>% 
  ggplot() +
  geom_tile(aes(x = name, y = term, fill = value), colour = "black", linewidth = 0.1) + 
  geom_text(aes(x = name, y = term, label = Significant, colour = value > breaks_oi[ceiling(length(breaks_oi)/2)]), size = 8) +
  scale_fill_viridis_b(breaks = breaks_oi, labels = function(x) ifelse(x %% 2 == 0, as.character(x), ""), limits = c(min(breaks_oi), max(breaks_oi)), option = "G", direction = -1) +
  scale_colour_manual(guide = "none", values = c("black", "white")) +
  theme_minimal(base_size = 20, base_family = "ArialMT") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    legend.key.size = unit(2, 'cm')
  ) +
  labs(
    x = "",
    y = "Search Term (details of hits)",
    fill = "% of gains"
  )

x_axis_labels <- unique(f3$data$name)[order(match(unique(f3$data$name), levels(f3$data$name)))]

kngdm_oi <- major_clades %>% 
  filter(as.numeric(node_name) > 508) %>% 
  dplyr::slice(match(str_split_i(x_axis_labels, pattern = "[|]", i = 1), node_name)) %>% 
  rowwise() %>% 
  mutate(parents = all_parents_finder(x = node_name)) %>% 
  mutate(kng_name = case_when(
    grepl("746", parents) ~ "Bacteria",
    grepl("548", parents) ~ "Eukaryota",
    grepl("511", parents) ~ "Archaea",
    TRUE ~ name
  )) %>% 
  distinct(node_name, kng_name) %>% 
  deframe()

euk_line_x <- match(names(kngdm_oi)[kngdm_oi == "Eukaryota"][1], str_split_i(x_axis_labels, pattern = "[|]", i = 1))
bac_line_x <- match(names(kngdm_oi)[kngdm_oi == "Bacteria"][1], str_split_i(x_axis_labels, pattern = "[|]", i = 1))

f3 <- f3 + geom_vline(xintercept = euk_line_x - 0.5, linewidth = 1)
f3 <- f3 + geom_vline(xintercept = bac_line_x - 0.5, linewidth = 1)
f3

#####Figure S13####
breaks_oi <- seq(0, ceiling(max(terms_oi_table$perc)), 1)
sf3 <- terms_wider %>% 
  pivot_longer(-term) %>% 
  left_join(terms_oi_table_new %>% dplyr::rename(name = phylum) %>% mutate(term = paste0(term, " (", n_groups, ")")), by = c("term", "name")) %>% 
  replace_na(list(n = 0)) %>%
  mutate(name = factor(name, levels = col_order)) %>%
  drop_na(name) %>% 
  ggplot(aes(x = name, y = term, fill = value)) +
  geom_tile(colour = "black", linewidth = 0.1) + 
  geom_text(aes(label = n, colour = value > breaks_oi[ceiling(length(breaks_oi)/2)]), size = 3, angle = 90) +
  scale_fill_viridis_b(breaks = breaks_oi, labels = function(x) ifelse(x %% 2 == 0, as.character(x), ""), option = "G", direction = -1) +
  scale_colour_manual(guide = "none", values = c("black", "white")) +
  theme_minimal(base_size = 16, base_family = "ArialMT") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    legend.key.size = unit(2, 'cm')
  ) +
  labs(
    x = "",
    y = "Search Term (# of hits)",
    fill = "% of gains"
  )

x_axis_labels <- unique(sf3$data$name)[order(match(unique(sf3$data$name), levels(sf3$data$name)))]

kngdm_oi <- major_clades %>% 
  filter(as.numeric(node_name) > 508) %>% 
  dplyr::slice(match(str_split_i(x_axis_labels, pattern = "[|]", i = 1), node_name)) %>% 
  rowwise() %>% 
  mutate(parents = all_parents_finder(x = node_name)) %>% 
  mutate(kng_name = case_when(
    grepl("746", parents) ~ "Bacteria",
    grepl("548", parents) ~ "Eukaryota",
    grepl("511", parents) ~ "Archaea",
    TRUE ~ name
  )) %>% 
  distinct(node_name, kng_name) %>% 
  deframe()

euk_line_x <- match(names(kngdm_oi)[kngdm_oi == "Eukaryota"][1], str_split_i(x_axis_labels, pattern = "[|]", i = 1))
bac_line_x <- match(names(kngdm_oi)[kngdm_oi == "Bacteria"][1], str_split_i(x_axis_labels, pattern = "[|]", i = 1))

sf3 <- sf3 + geom_vline(xintercept = euk_line_x - 0.5, size = 1.5)
sf3 <- sf3 + geom_vline(xintercept = bac_line_x - 0.5, size = 1.5)

sf3
