aggregate_local_metrics <- function(position_vector, target_vector, half_window_size) {
  map(position_vector, function(x) target_vector[position_vector %in% (x - half_window_size):(x + half_window_size)])
}

local_mean_func <- function(scores_vector, window_size) {
  sum(scores_vector)/window_size 
}

local_sd_func <- function(scores_vector, window_size) {
  scores_vector_2 <- c(scores_vector, rep(0, window_size - length(scores_vector)))
  local_sd = sd(scores_vector_2)
  return(local_sd)
}

normalize <- function(x) {
  (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
}

get_threshold_value <- function(y) {
  y <- sort(y)
  y_smooth <- rollmean(y, k = length(y)*0.2, fill = NA) %>% log2() %>% normalize()
  x = seq(0, 1, length.out = length(y_smooth))
  spline_function <- splinefun(x,y_smooth)
  threshold = which(spline_function(x, deriv = 1) > 1 & spline_function(x) > median(y_smooth, na.rm = T)) %>% min()
  threshold = 2^y[threshold]
  return(threshold)
}

find_tsr <- function(input_bed, half_window_size, technique = c("center_to_mean_ratio", "maxTSS"), pseudo_count = 1) {
  
  # Save fragment_lengths
  input_bed$frag_length <- input_bed %>% width()
  
  # Get position of 5' ends 
  five_prime_ends <- input_bed %>% resize(width = 1, fix = "start")
  five_prime_ends.tibble <- five_prime_ends %>% as_tibble() 
  
  # Summarize reads starting at each potential TSR 
  df <- five_prime_ends.tibble %>%
    dplyr::select(seqnames, strand, start, score, frag_length, frag_name = name) %>%
    chop(c(score, frag_length, frag_name)) %>%
    mutate(score = map_dbl(.$score, sum)) %>%
    group_by(strand) %>%
    mutate(local_scores = aggregate_local_metrics(start, score, half_window_size)) %>%
    ungroup()
  
  if("maxTSS" %in% technique) {
    df <- df %>% 
      mutate(max_score = map_dbl(.$local_scores, max))
    
    threshold_value <- get_threshold_value(df$max_score)
    
    df <- df %>% 
      mutate(maxTSS = score > threshold_value & score == max_score)
    
  } 
  
  if("center_to_mean_ratio" %in% technique) {
    df <- df %>% 
      mutate(local_mean = map_dbl(.$local_scores, local_mean_func, window_size = half_window_size*2)) %>%
      mutate(local_sd = map_dbl(.$local_scores, local_sd_func, window_size = half_window_size*2)) %>%
      mutate(score_to_mean_ratio = (score+pseudo_count)/(local_mean+pseudo_count)) %>%
      mutate(z_score = (score - local_mean)/local_sd)
    
    threshold_value <- get_threshold_value(df$score_to_mean_ratio)
    
    df <- df %>% 
      mutate(center_to_mean_ratio = score_to_mean_ratio > threshold_value)
  }
  
  return(df)
}





