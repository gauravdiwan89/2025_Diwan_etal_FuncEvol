library(phytools)
library(tidyverse)
library(DBI)
library(RPostgres)
library(dbplyr)

###First restore the database available here - xxx.xxx.xxx using the pg_restore command
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_pub", bigint = "integer")

nodes_list <- list(
  "Eukaryota" = c("781"),
  "Craniata" = c("876", "879"),
  "Metazoa" = c("834", "835"),
  "Magnoliopsida" = c("797", "798"),
  "Arthropoda" = c("867", "868"),
  "Cyanobacteria" = c("614", "617"),
  "Evosea" = c("823")
)
nodes_sig_terms <- readRDS("data/20250428_go_enrichment_data.rds")

go_obo <- ontologyIndex::get_ontology(file = "data/go-basic.obo", propagate_relationships = c("is_a", "part_of"))

full_go_lineage <- function(term, go_obo) {
  xx <- unlist(go_obo$name[unlist(go_obo$ancestors[go_obo$id[match(term, go_obo$name)]])])
  xx[-length(xx)]
}

nodes_sig_terms_desc <- lapply(1:length(nodes_sig_terms), function(x) {
  
  number <- 50
  only_desc_terms <- unique(unlist(sapply(nodes_sig_terms[[x]]$BP[1:number], function(x) full_go_lineage(term = x, go_obo = go_obo))))
  
  nodes_sig_terms[[x]] %>% 
    filter(BP %in% only_desc_terms) %>% 
    mutate(clade = names(nodes_list)[x]) #, node = list(nodes_list[[x]]))
}) %>% 
  bind_rows()

#order by log odds ratio
unique_items <- nodes_sig_terms_desc %>%
  dplyr::slice(1:10, .by = clade) %>%
  group_by(BP) %>%
  mutate(n = n_distinct(clade)) %>% 
  mutate(unique_in_group = n == 1) %>%
  distinct(BP, clade, unique_in_group) %>%
  filter(unique_in_group == TRUE) %>%
  ungroup()

opts <- tibble(
  clade = c("Eukaryota", "Metazoa", "Craniata", "Magnoliopsida", "Evosea", "Cyanobacteria"),
  empty_bar = c(3, 4, 4, 4, 3, 2),
  data_prep_cmd = c(
    "rbind(data[1:8,], to_add[1,], data[9:10,], to_add[2:3,])", #Eukaryota
    "rbind(data[1:2,], to_add[1,], data[3:8,], to_add[2,], data[9:10,], to_add[3:4,])", #Metazoa
    "rbind(data[1:2,], to_add[1,], data[3:8,], to_add[2,], data[9:10,], to_add[3:4,])", #Craniata
    "rbind(data[1:2,], to_add[1,], data[3:8,], to_add[2,], data[9:10,], to_add[3:4,])", #Angiosperms
    "rbind(data[1:2,], to_add[1:2,], data[3:8,], data[9:10,], to_add[3,])", #Evosea
    "rbind(data, to_add)" #Cyanobacteria
  )
)

go_plot_list <- list()
for(i in 1:nrow(opts)) {
  clade_oi <- opts$clade[i]
  data <- split(x = nodes_sig_terms_desc, f = nodes_sig_terms_desc$clade)[[clade_oi]] %>% 
    filter(term_not_at_node_totals > 0) %>%
    dplyr::slice(1:10, .by = clade)
  
  empty_bar <- opts$empty_bar[i]
  
  to_add <- matrix(NA, empty_bar, ncol(data))
  colnames(to_add) <- colnames(data)
  
  data <- eval(parse(text = opts$data_prep_cmd[i]))
  
  data$id <- seq(1, nrow(data))
  
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- (90 - ((1)*15)) - 360 * (label_data$id-0.5) / number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  max_y <- max(data$term_at_node_totals, na.rm = T)
  
  go_plot_list[[clade_oi]] <- ggplot(data, aes(x = factor(id), y = term_at_node_totals, fill = log2(fishers_odds_ratio))) +
    geom_bar(stat = "identity") +
    coord_polar(start = pi/(12)) +
    scale_y_continuous(
      limits = c(-25, NA),
      expand = c(0, 0)
    ) +
    scale_fill_viridis_c(limits = c(0.9, 6)) +
    theme_minimal() +
    theme(
      # Remove axis ticks and text
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "bottom",
      # Ensure labels are not clipped
      plot.margin = unit(rep(0, 4), "cm"),
      panel.grid.major.y = element_line(color = "gray80", size = 0.5)
    ) + 
    geom_text(data = label_data, aes(x = id, y = term_at_node_totals + 5, label = BP, hjust = hjust), color = "black",alpha = 1, size = 3.5, angle =  label_data$angle, inherit.aes = FALSE )  +
    # Add y-axis labels manually
    annotate("text", x = 0.3, y = seq(0, max_y+25, by = 25), label = seq(0, max_y+25, by = 25), angle = 0, vjust = -1, hjust = 1, color = "black", size = 3.5) +
    labs(
      x = "",
      y = "Number of Genes",
      fill = "Log2 Odds Ratio"
    )
}

go_plot_list$Eukaryota
go_plot_list$Metazoa
go_plot_list$Craniata
go_plot_list$Magnoliopsida
go_plot_list$Evosea
go_plot_list$Cyanobacteria

###KEGG gains####
kegg_hierarchy <- read_tsv("data/20230713_kegg_pathway_hierarchy.tsv.gz")

kegg_plot_list <- list()
for(i in 1:nrow(opts)) {
  clade_oi <- opts$clade[i]
  
  data <- tbl(con, "kegg_gains_parsimony") %>% 
    filter(node %in% !!nodes_list[[clade_oi]]) %>% 
    group_by(pathway_name) %>% 
    summarize(prop_comp_gained = sum(prop_comp_gained), n_comp_gained = sum(n_comp_gained)) %>% 
    collect() %>% 
    left_join(kegg_hierarchy) %>% 
    filter(H1 != "Human Diseases") %>%
    {if(i <= 3) group_by(., H2) else group_by(., pathway_name)} %>% 
    summarize(avg_prop_gained = mean(prop_comp_gained), total_comp_gained = sum(n_comp_gained)) %>% 
    arrange(-avg_prop_gained) %>% 
    slice_max(order_by = avg_prop_gained, n = 10)
  
  empty_bar <- 2
  
  to_add <- matrix(NA, empty_bar, ncol(data))
  colnames(to_add) <- colnames(data)
  
  data <- rbind(data, to_add)
  
  data$id <- seq(1, nrow(data))
  
  label_data <- data
  colnames(label_data)[colnames(label_data) == "pathway_name"] <- "H2"
  number_of_bar <- nrow(label_data)
  angle <- (90 - 30) - 360 * (label_data$id-0.5) / number_of_bar  # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  max_y <- max(data$avg_prop_gained*100, na.rm = T)
  
  kegg_plot_list[[clade_oi]] <- ggplot(data, aes(x = factor(id), y = avg_prop_gained*100, fill = total_comp_gained)) +
    geom_bar(stat = "identity") +
    coord_polar(start = pi/7) +
    scale_y_continuous(
      limits = c(-20, NA),
      expand = c(0, 0)
    ) +
    scale_fill_viridis_c(option = "F", limits = c(0, if(i <= 3) 650 else max(data$total_comp_gained))) +
    theme_minimal() +
    theme(
      # Remove axis ticks and text
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "bottom",
      # Ensure labels are not clipped
      plot.margin = unit(rep(0, 4), "cm"),
      panel.grid.major.y = element_line(color = "gray80", size = 0.5)
    ) + 
    geom_text(data = label_data, aes(x = id, y = (avg_prop_gained*100) + 5, label = H2, hjust = hjust), color = "black",alpha = 1, size = 3.5, angle =  label_data$angle, inherit.aes = FALSE )  +
    # Add y-axis labels manually
    annotate("text", x = 0.4, y = seq(0, max_y + 5, by = 20), label = seq(0, max_y + 5, by = 20), angle = 0, vjust = -1, hjust = 1, color = "black", size = 3.5) +
    labs(
      x = "",
      y = "",
      fill = "# components gained"
    )
}

kegg_plot_list$Eukaryota
kegg_plot_list$Metazoa
kegg_plot_list$Craniata
kegg_plot_list$Magnoliopsida
kegg_plot_list$Evosea
kegg_plot_list$Cyanobacteria

###Figure 2E####

species_tree <- read.tree("data/species_tree_cleaned_final.nwk")
major_clades <- read_tsv("data/TableS2.tsv")

all_parents_finder <- function(x) {
  all_parents <- NULL
  parent_oi <- getParent(species_tree, as.numeric(x))
  all_parents <- c(all_parents, parent_oi)
  while(!is.null(parent_oi)) {
    parent_oi <- getParent(species_tree, parent_oi)
    all_parents <- c(all_parents, parent_oi)
  }
  toString(unique(all_parents))
}

term_origins <- function(search_term, mode = "gene") {
  if(mode == "gene") {
    og_gains_losses <- tbl(con, "orthologs_gains_losses_parsimony")
    print(paste("Starting with search term:", search_term, "..."))
    sleep_og_hist <- tbl(con, "ptns_go_bp") %>% 
      mutate(keyword = as.character(search_term)) %>% 
      filter(if_any(.cols = -keyword, .fns = ~ grepl(pattern = keyword, x = .x, ignore.case = TRUE))) %>% 
      distinct(uniprot_acc) %>%
      inner_join(tbl(con, "ptns_og_phylo")) %>%
      inner_join(og_gains_losses) %>% 
      collect()
    
    print("Counting gains...")
    counts_table <- sleep_og_hist %>% 
      pivot_longer(cols = c(gain_nodes), names_to = "gain_type", values_to = "node") %>% 
      mutate(node = strsplit(node, split = ", ")) %>% 
      unnest() %>% 
      filter(node != "") %>% 
      mutate(node = as.numeric(node)) %>% 
      group_by(node) %>%
      summarize(n = length(unique(HOG))) %>%
      left_join(major_clades %>% dplyr::select(node, name, level = taxa_level, full_taxonomy, full_taxonomy_levels, other_names), by = c("node")) %>% 
      filter(as.numeric(node) > Ntip(species_tree)) %>%
      arrange(-n)
  }
  
  if(nrow(counts_table) > 0) {
    print(paste("Building taxonomy for", nrow(counts_table), "entries ..."))
    
    counts_table <- counts_table %>%
      mutate(
        taxa_names = str_split(full_taxonomy, ", ", simplify = FALSE),
        taxa_levels = str_split(full_taxonomy_levels, ", ", simplify = FALSE),
        phylum_name = map2_chr(taxa_levels, taxa_names, ~if_else("order" %in% .x, .y[which(.x == "order")][1], NA_character_))
      ) %>% 
      mutate(
        phylum = case_when(
          level %in% c("kingdom", "superkingdom", "phylum", "subphylum", "class", "order") ~ name,
          !is.na(phylum_name) ~ phylum_name,
          TRUE ~ name
        )
      ) %>%
      select(-taxa_names, -taxa_levels, -phylum_name)
    
    counts_table %>% 
      group_by(phylum) %>%
      summarize(n = sum(n)) %>% 
      mutate(
        perc = round(n/sum(n)*100, digits = 1),
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

terms_oi_table <- rbind(
  term_origins(search_term = "photosynthesis", mode = "gene"),
  term_origins(search_term = "brain", mode = "gene"),
  term_origins(search_term = "nuclear pore", mode = "gene"),
  term_origins(search_term = "adaptive immun", mode = "gene"),
  term_origins(search_term = "heart development", mode = "gene"),
  term_origins(search_term = "flagella|flagellum", mode = "gene"),
  term_origins(search_term = "CRISPR", mode = "gene")
)

terms_wider <- terms_oi_table %>% 
  mutate(term = paste0(term, " (", n_groups, ")")) %>% 
  group_by(term) %>% 
  mutate(c_perc = cumsum(perc)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = term, names_from = phylum, values_from = perc) %>% 
  mutate(across(-term, ~ replace_na(., 0)))

col_order_temp <- colnames(terms_wider)[-1][match(unique(major_clades$name[as.numeric(major_clades$node)>508]), colnames(terms_wider)[-1], nomatch = F)]
col_order <- unique(col_order_temp[!is.na(col_order_temp)])

scale_high <- if(ceiling(max(terms_oi_table$perc)) %% 5 == 0) {
  ceiling(max(terms_oi_table$perc))
} else {
  xx <- ceiling(max(terms_oi_table$perc))
  while(xx %% 5 != 0) {
    xx <- xx + 1
  }
}
breaks_oi <- seq(0, xx, 2)
f3 <- terms_wider %>% 
  pivot_longer(-term) %>% 
  left_join(terms_oi_table %>% dplyr::rename(name = phylum) %>% mutate(term = paste0(term, " (", n_groups, ")"), Significant = if_else(p_value < 0.05, "*", "")), by = c("term", "name")) %>% 
  replace_na(list(n = 0)) %>% 
  mutate(perc_2.5 = value > 2) %>% 
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
  filter(as.numeric(node) > 508) %>% 
  dplyr::slice(match(unique(f3$data$name), name)) %>% 
  rowwise() %>% 
  mutate(parents = all_parents_finder(x = node)) %>% 
  mutate(kng_name = case_when(
    grepl("510", parents) ~ "Bacteria",
    grepl("781", parents) ~ "Eukaryota",
    grepl("979", parents) ~ "Archaea",
    TRUE ~ name
  )) %>% 
  distinct(name, kng_name) %>% 
  deframe()

euk_line_x <- match(names(kngdm_oi)[kngdm_oi == "Eukaryota"][1], x_axis_labels)
arc_line_x <- match(names(kngdm_oi)[kngdm_oi == "Archaea"][1], x_axis_labels)

f3 <- f3 + geom_vline(xintercept = euk_line_x - 0.5, linewidth = 1)
f3 <- f3 + geom_vline(xintercept = arc_line_x - 0.5, linewidth = 1)

f3

##Figure S10####
breaks_oi <- seq(0, ceiling(max(terms_oi_table$perc)), 1)
sf3 <- terms_wider %>% 
  pivot_longer(-term) %>% 
  left_join(terms_oi_table %>% dplyr::rename(name = phylum) %>% mutate(term = paste0(term, " (", n_groups, ")")), by = c("term", "name")) %>% 
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
  filter(as.numeric(node) > 508) %>% 
  dplyr::slice(match(unique(sf3$data$name), name)) %>% 
  rowwise() %>% 
  mutate(parents = all_parents_finder(x = node)) %>% 
  mutate(kng_name = case_when(
    grepl("510", parents) ~ "Bacteria",
    grepl("781", parents) ~ "Eukaryota",
    grepl("979", parents) ~ "Archaea",
    TRUE ~ name
  )) %>% 
  distinct(name, kng_name) %>% 
  deframe()

euk_line_x <- match(names(kngdm_oi)[kngdm_oi == "Eukaryota"][1], x_axis_labels)
arc_line_x <- match(names(kngdm_oi)[kngdm_oi == "Archaea"][1], x_axis_labels)

sf3 <- sf3 + geom_vline(xintercept = euk_line_x - 0.5, size = 1.5)
sf3 <- sf3 + geom_vline(xintercept = arc_line_x - 0.5, size = 1.5)

sf3
