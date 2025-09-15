library(phytools)
library(tidyverse)
library(ggtree)

###Load Data####
species_details <- read_tsv("data/TableS1.tsv")
tables2 <- read_tsv("data/TableS2.tsv")
species_tree <- read.tree("data/species_tree_cleaned_final.nwk")
phylum_tree_final <- read.tree("data/20250508_Phylum_Tree_for_Figure_1.nwk")

###Change tip labels for phylum tree in fig1####
orig_ph_tree_tips <- phylum_tree_final$tip.label

ph_tree_tips <- species_details$phylum[match(phylum_tree_final$tip.label, species_details$tree_tip_label)]
ph_tree_tips[is.na(ph_tree_tips)] <- species_details$family[match(phylum_tree_final$tip.label[is.na(ph_tree_tips)], species_details$tree_tip_label)]

species_details_longer <- species_details %>% 
  pivot_longer(cols = phylum:family, names_to = "taxon")

tips_split <- split(orig_ph_tree_tips, ph_tree_tips)

ph_tree_clades <- sapply(tips_split[sapply(tips_split, length) > 1], function(x) getMRCA(phylum_tree_final, x))

ph_tree_tips[ph_tree_tips %in% names(ph_tree_clades)] <- ""

ph_tree_clade_labs <- ph_tree_clades %>% enframe() %>% dplyr::rename(clade = name, node = value)

ph_tree_tips[ph_tree_tips != ""] <- paste0(ph_tree_tips[ph_tree_tips != ""], " (", sapply(ph_tree_tips[ph_tree_tips != ""], function(x) {
  
  spp_oi <- species_details_longer %>% 
    filter(value == x) %>% 
    pull(tree_tip_label) %>% 
    unique()
  
  species_details_lng_oi <- species_details %>% 
    filter(tree_tip_label %in% spp_oi) %>% 
    pivot_longer(cols = phylum:family, names_to = "taxon") %>% 
    mutate(rn = row_number(), .by = tree_tip_label)
  
  rn_oi <- max(unique(na.omit(species_details_lng_oi$rn[species_details_lng_oi$value == x])))
  
  xx <- sum(species_details_lng_oi$value == x & species_details_lng_oi$rn == rn_oi, na.rm = T) - (
    species_details_lng_oi %>% 
      filter(rn > rn_oi) %>% 
      filter(value %in% c(ph_tree_tips, ph_tree_clade_labs$clade)) %>% 
      nrow()
  )
  xx
}), ")")

ph_tree_clade_labs$clade <- paste0(ph_tree_clade_labs$clade, " (", sapply(ph_tree_clade_labs$clade, function(x) {
  spp_oi <- species_details_longer %>% 
    filter(value == x) %>% 
    pull(tree_tip_label) %>% 
    unique()
  
  species_details_lng_oi <- species_details %>% 
    filter(tree_tip_label %in% spp_oi) %>% 
    pivot_longer(cols = phylum:family, names_to = "taxon") %>% 
    mutate(rn = row_number(), .by = tree_tip_label)
  
  rn_oi <- max(unique(na.omit(species_details_lng_oi$rn[species_details_lng_oi$value == x])))
  
  xx <- sum(species_details_lng_oi$value == x & species_details_lng_oi$rn == rn_oi, na.rm = T) - (
    species_details_lng_oi %>% 
      filter(rn > rn_oi) %>% 
      filter(value %in% c(gsub(" [(]\\d+[)]", "", ph_tree_tips), ph_tree_clade_labs$clade)) %>% 
      nrow()
  )
  
  xx
}), ")")

ph_tree_clade_labs[1, 1] <- "Actinobacteria (41)"

phylum_tree_final$tip.label <- ph_tree_tips

###Combined origins of functions, pathways and domains - Fig.1####
top_events_all <- tables2 %>% 
  select(node, name, n_descendants, ends_with("_gained")) %>% 
  pivot_longer(cols = ends_with("_gained"), names_to = "Type", values_to = "n") %>% 
  mutate(Type = str_replace(gsub("_", " ", gsub("_gained", "", Type)), "^[a-z]", ~str_to_upper(.x))) %>% 
  filter(n_descendants > 2)

set.seed(2398473)
pt <- ggtree::ggtree(phylum_tree_final, size = 0.3, colour = gray(level = 0.3), layout = "circular", branch.length = "none")
ph_tree_alt <- pt %>% 
  rotate(201) %>% 
  rotate(214) %>% 
  rotate(200) +
  geom_rootedge(rootedge = .5) +
  geom_tiplab(align = T, offset = 0.5, linesize = 0, size = 5) +
  geom_cladelab(
    data = ph_tree_clade_labs,
    mapping = aes(
      node = node,
      label = clade
    ), 
    geom = "text", 
    lineheight = 0, 
    fontsize = 5,
    angle = "auto",
    align = T,
    show.legend = F
  )

final_top_events <- top_events_all %>% 
  filter(node %in% as.numeric(str_match(phylum_tree_final$node.label, "\\d+")[,1])) %>% 
  select(node, n, Type)

t1_final <- nest(final_top_events, n = n, Type = Type)

t1_final <- t1_final %>% 
  mutate(node = if_else(paste0("Node", node) %in% ph_tree_alt$data$label, ph_tree_alt$data$node[match(paste0("Node", node), ph_tree_alt$data$label)], node))

ph_tree_alt$data <- ph_tree_alt$data %>% left_join(t1_final, by = "node")
ph_tree_alt$data <- ph_tree_alt$data %>% 
  mutate(n = if_else(as.numeric(str_match(label, "\\d+")[,1]) %in% c(1, 2), list(NULL), n), Type = if_else(as.numeric(str_match(label, "\\d+")[,1]) %in% c(1, 2), list(NULL), Type))

ph5_alt <- ph_tree_alt +
  geom_nodepoint(mapping = aes(node = node, size = n, colour = Type), shape = 21, alpha = 1, stroke = 1, data = td_unnest(c(n, Type))) + 
  scale_size(range = c(0.5, 15), breaks = c(ceiling(min(final_top_events$n)), 10, 25, 50, 100, 250, 500, 1000, 1500, 2500, ceiling(max(final_top_events$n))), labels = as.character(c(ceiling(min(final_top_events$n)), 10, 25, 50, 100, 250, 500, 1000, 1500, 2500, ceiling(max(final_top_events$n)))), limits = c(floor(min(final_top_events$n)), ceiling(max(final_top_events$n)))) +
  scale_color_manual(values = c(
    "Genes" = "#4daf4a",
    "Domains" = "#e41a1c",
    "KEGG Pathways" = "#377eb8",
    "Reactome Pathways" = "#984ea3"
  )) +
  guides(colour = guide_legend(override.aes = list(size = 4)), size = guide_legend(ncol = 2)) +
  theme(legend.position = "left", text = element_text(size = 14), legend.text = element_text(size = 14), legend.spacing = unit(300, "pt")) +
  labs(
    size = "Number of\nFounder Events\n(size)",
    colour = "Biological Entity (colour)"
  )

ph5_alt
