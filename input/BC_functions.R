
#### BC distribution table ####
# x = BC association table with columns full_BC, TargetId, n_associated (number
# of times each BC associates with each target) and BC_shared_between (how many
# targets share the BC).
# n = the max size of different combinations of targets.

BC_tables <- function(x, n){
  
  combinations <- list()
  p_table <- data.frame()
  
  n_table <- x %>%
    filter(BC_shared_between == 1) %>%
    group_by(TargetId) %>%
    summarise(n = n_distinct(full_BC),
              n_min_10 = n_distinct(full_BC[n_associated >= 10])) %>%
    ungroup()
  
  for(m in 2:n){
    c_temp <- combn(unique(x$TargetId), m = m, simplify = F)
    combinations <- c(combinations, c_temp)
  }
  
  for(i in 1:length(combinations)){
    
    temp_n <- x %>%
      filter(TargetId %in% combinations[[i]] & BC_shared_between == length(combinations[[i]])) %>%
      group_by(full_BC) %>% filter(n() == length(combinations[[i]])) %>% ungroup() %>%
      summarise(n = n_distinct(full_BC)) %>%
      mutate(TargetId = paste((combinations[[i]]), collapse = " & "))
    
    temp_n_10 <- x %>%
      filter(TargetId %in% combinations[[i]] & BC_shared_between == length(combinations[[i]])) %>%
      group_by(full_BC) %>% 
      filter(n() == length(combinations[[i]]) & all(n_associated >= 10)) %>% 
      ungroup() %>%
      summarise(n_min_10 = n_distinct(full_BC)) %>%
      mutate(TargetId = paste((combinations[[i]]), collapse = " & "))
    
    temp <- full_join(temp_n, temp_n_10)
    
    if(length(combinations[[i]]) == 2){
      
      temp_f <- x %>%
        filter(TargetId %in% combinations[[i]]) %>%
        summarise(n_all = n_distinct(full_BC)) %>%
        mutate(TargetId = paste((combinations[[i]]), collapse = " & ")) %>%
        full_join(temp)
      
      p_table <- rbind(p_table, temp_f)
    }
    
    n_table <- rbind(n_table, temp)
    
  }
  
  p_table <- p_table %>%
    mutate(prop = round(n/n_all, 4)*100,
           prop_min_10 = round(n_min_10/n_all, 4)*100) %>%
    relocate(TargetId)
  
  list(p_table = p_table, n_table = n_table)
}


#### BC distribution plots ####
# x = BC association table with columns full_BC, TargetId, n_associated (number
# of times each BC associates with each target) and BC_shared_between (how many
# targets share the BC).
BC_plots <- function(x){
  
  BC_plots <- list()
  
  combinations <- combn(unique(x$TargetId), m = 2, simplify = F)
  
  input <- x %>%
    mutate(colour_col = ifelse(BC_shared_between == 1 | BC_shared_between == 2,
                           "black",
                           ifelse(BC_shared_between == 3,
                                  "skyblue2",
                                  ifelse(BC_shared_between >= 4,
                                         "orange",
                                         "red"))))
  
  for(i in 1:length(combinations)){
    x_val <- combinations[[i]][1]
    y_val <- combinations[[i]][2]
    
    name <- paste(combinations[[i]], collapse = " & ")
    
    BC_plots[[name]] <- input %>%
      filter(TargetId %in% combinations[[i]]) %>%
      pivot_wider(names_from = TargetId, values_from = n_associated, values_fill = 0) %>%
      ggplot(aes(x = .data[[x_val]], y = .data[[y_val]], colour = colour_col)) +
      geom_point(size = 1) +
      scale_color_identity() +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(0, 10, 100, 1000, 10000), minor_breaks = NULL) +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(0, 10, 100, 1000, 10000), minor_breaks = NULL) +
      theme_minimal()
      
  }
  
  return(BC_plots)
  
}
  
