###Libraries and Database####
library(phytools)
library(tidyverse)
library(DBI)
library(RPostgres)
library(dbplyr)

###First restore the database available here - xxx.xxx.xxx using the pg_restore command
con <- dbConnect(drv = RPostgres::Postgres(), dbname = "orthologs_pub", bigint = "integer")

species_tree <- read.tree("data/species_tree_cleaned_final.nwk")
species_details <- read_tsv("data/TableS1.tsv")

###Fig 3####

euk_tree <- keep.tip(read.newick("data/species_tree_for_Count.nwk"), species_details$tree_tip_label[species_details$superkingdom == "Eukaryota"])

####Habitat transitions####
##ancestral reconstruction of habitat
habitat_species <- species_details %>% 
  select(OSCODE = tree_tip_label, species, habitat) %>% 
  filter(!is.na(habitat)) #drops the bacteria and archaea

# habitat_anc <- make.simmap(tree = euk_tree, x = habitat_species %>% select(OSCODE, habitat) %>% deframe(), model = "ER", nsim = 100, pi = c("Aquatic" = 1, "Terrestrial" = 0, "Host-associated"= 0, "Amphibian" = 0), Q = "mcmc")
# save(habitat_anc, file = "data/20250514_habitat_anc_reconstruction.simmap")

##analyze the output
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

####Cellularity transitions####
##ancestral reconstruction of cellularity
euk_cell <- species_details %>% 
  select(OSCODE = tree_tip_label, species, cellularity) %>% 
  filter(!is.na(cellularity)) #drops the bacteria and archaea

# cell_anc <- make.simmap(tree = euk_tree, x = euk_cell %>% select(OSCODE, cellularity) %>% deframe(), model = "ER", nsim = 100, pi = c(0, 1), Q = "mcmc")
# save(cell_anc, file = "data/20250428_cellularity_anc_reconstruction.simmap")

##analyze the output
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
###Fig 3####

par(mfcol = c(1, 2))

plot(habitat_anc_summary, ftype = "off", fsize = 0.5, offset = 2, lwd = 1, cex = 0.4, align.tip.label= T)
tiplabels(text = species_details$species[match(euk_tree$tip.label, species_details$tree_tip_label)], frame = "none", cex = 0.5, offset = 0.3)
pp.2 <- get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:Ntip(euk_tree)) lines(c(pp.2$xx[i] + 0.02, max(pp.2$xx)), rep(pp.2$yy[i], 2), lty = "dotted")
tiplabels(tip = habitat_transitions %>% filter(change == 1, d_state == "Terrestrial") %>% pull(d_node) %>% as.numeric(), pch = 1, col = "#eba2c7", cex = 3, lwd = 5)

plot(cellularity_anc_summary, ftype = "off", offset = 2, cex = 0.4, direction = "leftwards")
tiplabels(tip = cellularity_transitions %>% filter(change == 1, p_state == "Unicellular", d_state == "Multicellular") %>% rowwise() %>% filter(length(getDescendants(euk_tree, d_node)) > 2) %>% pull(d_node) %>% as.numeric(), pch = 1, col = "#00aad4", cex = 3, lwd = 5)

##this plot was exported and then formatted in Adobe Illustrator to produce figure 3