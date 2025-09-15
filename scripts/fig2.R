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
  "Evosea" = c("823"),
  "Chytridiomycetes" = "423"
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
  
  plot_list[[clade_oi]] <- ggplot(data, aes(x = factor(id), y = term_at_node_totals, fill = log2(fishers_odds_ratio))) +
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