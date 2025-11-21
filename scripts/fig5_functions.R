og_anc_reconst <- function(og, method = "mp") {
  og_gains_losses <- tbl(con, "og_gains_losses_parsimony")
  
  og_oi_table <- og_gains_losses %>% 
    filter(HOG == og) %>% 
    collect()
  
  anc_state_oi <- og_oi_table$anc_state[1]
  
  all_states <- NULL
  all_states[as.character(getDescendants(species_tree, node = Ntip(species_tree) + 1))] <- anc_state_oi
  
  if(anc_state_oi == 1) {
    
    all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
    loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
    if(!all(is.na(loss_nodes))) {
      all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
      all_states[c(loss_nodes, all_loss_desc)] <- 0
    }
    
    gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
    if(!all(is.na(gain_nodes))) {
      all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
      if(any(all_gain_desc %in% loss_nodes)) {
        all_drop_desc <- as.character(unique(unlist(lapply(all_gain_desc[all_gain_desc %in% loss_nodes], function(x) getDescendants(species_tree, as.numeric(x))))))
        all_gain_desc <- all_gain_desc[!all_gain_desc %in% all_drop_desc]
      }
      all_states[c(gain_nodes, all_gain_desc)] <- 1
    }
    
    loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
    if(!all(is.na(loss_tips))) all_states[loss_tips] <- 0
    
    gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
    if(!all(is.na(gain_tips))) all_states[gain_tips] <- 1
    
  } else if(anc_state_oi == 0) {
    
    all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
    
    gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
    if(!all(is.na(gain_nodes))) {
      all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
      all_states[c(gain_nodes, all_gain_desc)] <- 1
    }
    
    loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
    if(!all(is.na(loss_nodes))) {
      all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
      if(any(all_loss_desc %in% gain_nodes)) {
        all_drop_desc <- as.character(unique(unlist(lapply(all_loss_desc[all_loss_desc %in% gain_nodes], function(x) getDescendants(species_tree, as.numeric(x))))))
        all_loss_desc <- all_loss_desc[!all_loss_desc %in% all_drop_desc]
      }
      all_states[c(loss_nodes, all_loss_desc)] <- 0
    }
    
    gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
    if(!all(is.na(gain_tips))) all_states[gain_tips] <- 1
    
    loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
    if(!all(is.na(loss_tips))) all_states[loss_tips] <- 0
    
  }
  
  node_states <- tibble(
    node = as.numeric(names(all_states)),
    State = if_else(all_states == 0, "Absent", "Present"),
    change = if_else(names(all_states) %in% c(gain_nodes, loss_nodes, gain_tips, loss_tips), 1, 0),
    gain = if_else(names(all_states) %in% c(gain_nodes, gain_tips), 1, 0),
    loss = if_else(names(all_states) %in% c(loss_nodes, loss_tips), 1, 0)
  )
  
  node_states
}

orthologs_anc_reconst <- function(gene, method = "mp") {
  
  if(method == "mp") {
    og_gains_losses <- tbl(con, "orthologs_gains_losses_parsimony")
    
    og_oi_table <- og_gains_losses %>% 
      filter(uniprot_acc == gene) %>% 
      collect() %>% 
      mutate(across(gain_tips:loss_nodes, ~ replace_na(.x, "")))
    
    if(nrow(og_oi_table) > 0) {
      anc_state_oi <- as.numeric(og_oi_table$anc_state[1])
      
      all_states <- NULL
      all_states[as.character(getDescendants(species_tree, node = Ntip(species_tree) + 1))] <- anc_state_oi
      
      if(anc_state_oi == 1) {
        
        all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
        loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
        if(!all(is.na(loss_nodes))) {
          all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
          all_states[c(loss_nodes, all_loss_desc)] <- 0
        }
        
        gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
        if(!all(is.na(gain_nodes))) {
          all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
          if(any(all_gain_desc %in% loss_nodes)) {
            all_drop_desc <- as.character(unique(unlist(lapply(all_gain_desc[all_gain_desc %in% loss_nodes], function(x) getDescendants(species_tree, as.numeric(x))))))
            all_gain_desc <- all_gain_desc[!all_gain_desc %in% all_drop_desc]
          }
          all_states[c(gain_nodes, all_gain_desc)] <- 1
        }
        
        loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
        if(!all(is.na(loss_tips))) all_states[loss_tips] <- 0
        
        gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
        if(!all(is.na(gain_tips))) all_states[gain_tips] <- 1
        
      } else if(anc_state_oi == 0) {
        
        all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
        
        gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
        if(!all(is.na(gain_nodes))) {
          all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
          all_states[c(gain_nodes, all_gain_desc)] <- 1
        }
        
        loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
        if(!all(is.na(loss_nodes))) {
          all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
          if(any(all_loss_desc %in% gain_nodes)) {
            all_drop_desc <- as.character(unique(unlist(lapply(all_loss_desc[all_loss_desc %in% gain_nodes], function(x) getDescendants(species_tree, as.numeric(x))))))
            all_loss_desc <- all_loss_desc[!all_loss_desc %in% all_drop_desc]
          }
          all_states[c(loss_nodes, all_loss_desc)] <- 0
        }
        
        gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
        if(!all(is.na(gain_tips))) all_states[gain_tips] <- 1
        
        loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
        if(!all(is.na(loss_tips))) all_states[loss_tips] <- 0
        
      }
      
      org_oi <- tbl(con, "ptns_og_phylo") %>% filter(uniprot_acc == gene) %>% pull(org)
      all_states[as.character(match(org_oi, species_tree$tip.label))] <- 1
      
      node_states <- tibble(
        node = as.numeric(names(all_states)),
        State = if_else(all_states == 0, "Absent", "Present"),
        change = if_else(names(all_states) %in% c(gain_nodes, loss_nodes, gain_tips, loss_tips), 1, 0),
        gain = if_else(names(all_states) %in% c(gain_nodes, gain_tips), 1, 0),
        loss = if_else(names(all_states) %in% c(loss_nodes, loss_tips), 1, 0)
      )
    } else {
      all_states <- NULL
      all_states[as.character(getDescendants(species_tree, node = Ntip(species_tree) + 1))] <- 0
      org_oi <- tbl(con, "ptns_og_phylo") %>% filter(uniprot_acc == gene) %>% pull(org)
      all_states[as.character(match(org_oi, species_tree$tip.label))] <- 1
      node_states <- tibble(
        node = as.numeric(names(all_states)),
        State = if_else(all_states == 0, "Absent", "Present"),
        change = 2,
        gain = 0,
        loss = 0
      )
    }
    node_states
    
  }
}

pfam_node_states <- function(pfam_oi) {
  og_gains_losses <- tbl(con, "pfam_domains_gains_losses_parsimony")
  
  og_oi_table <- og_gains_losses %>% 
    filter(pfam_id == pfam_oi) %>% 
    collect() %>% 
    mutate(across(gain_tips:loss_nodes, ~ replace_na(.x, "")))
  
  anc_state_oi <- as.numeric(og_oi_table$anc_state[1])
  
  all_states <- NULL
  all_states[as.character(getDescendants(species_tree, node = Ntip(species_tree) + 1))] <- anc_state_oi
  
  if(anc_state_oi == 1) {
    
    all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
    loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
    all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
    all_states[c(loss_nodes, all_loss_desc)] <- 0
    
    gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
    all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
    all_states[c(gain_nodes, all_gain_desc)] <- 1
    
    loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
    all_states[loss_tips] <- 0
    
    gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
    all_states[gain_tips] <- 1
    
  } else if(anc_state_oi == 0) {
    
    all_states[as.character(Ntip(species_tree) + 1)] <- anc_state_oi
    
    gain_nodes <- strsplit(og_oi_table$gain_nodes[1], split = ", ")[[1]]
    all_gain_desc <- as.character(unique(unlist(lapply(gain_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
    all_states[c(gain_nodes, all_gain_desc)] <- 1
    
    loss_nodes <- strsplit(og_oi_table$loss_nodes[1], split = ", ")[[1]]
    all_loss_desc <- as.character(unique(unlist(lapply(loss_nodes, function(x) getDescendants(species_tree, as.numeric(x))))))
    all_states[c(loss_nodes, all_loss_desc)] <- 0
    
    gain_tips <- strsplit(og_oi_table$gain_tips[1], split = ", ")[[1]]
    all_states[gain_tips] <- 1
    
    loss_tips <- strsplit(og_oi_table$loss_tips[1], split = ", ")[[1]]
    all_states[loss_tips] <- 0
    
  }
  
  node_states <- tibble(
    node = as.numeric(names(all_states)),
    State = if_else(all_states == 0, "Absent", "Present"),
    change = if_else(names(all_states) %in% c(gain_nodes, loss_nodes, gain_tips, loss_tips), 1, 0),
    gain = if_else(names(all_states) %in% c(gain_nodes, gain_tips), 1, 0),
    loss = if_else(names(all_states) %in% c(loss_nodes, loss_tips), 1, 0)
  )
  
  node_states
}