library(MatchIt)
library(dplyr)
library(MatchIt)

compare_hungtoprop <- function(n_treat, mu_treat, sd_treat, mu_control, sd_control, percentile) {   


    library(dplyr)
    
    t_1 <- rnorm(n_treat ,mean = mu_treat, sd = sd_treat) #Higher
    t_2 <- rnorm(n_treat ,mean = mu_treat, sd = sd_treat) #Lower 
    t_3 <- rnorm(n_treat ,mean = mu_treat, sd = sd_treat)
    t_4 <- rnorm(n_treat ,mean = mu_treat, sd = sd_treat)
    t_5 <- rnorm(n_treat ,mean = mu_treat, sd = sd_treat)
    t_assignment<- rep(1, n_treat)
    
    c_1 <- rnorm(n_treat ,mean = mu_control, sd = sd_control)
    c_2 <- rnorm(n_treat,mean = mu_control, sd = sd_control)
    c_3 <- rnorm(n_treat ,mean = mu_control, sd = sd_control)
    c_4 <- rnorm(n_treat ,mean = mu_control, sd = sd_control)
    c_5 <- rnorm(n_treat ,mean = mu_control, sd = sd_control)
    c_assignment <- rep(0,n_treat)
    
    
    data_t <- data.frame("X1" = t_1, "X2" = t_2, "X3" = t_3, "X4" = t_4, "X5" = t_5, "Assignment" = t_assignment )
    data_c <- data.frame("X1" = c_1, "X2" = c_2, "X3" = c_3, "X4" = c_4, "X5" = c_5, "Assignment"= c_assignment)
    
    data_combined <- bind_rows(data_t, data_c)
    
    
    
    
    #---------------------------------------------------------Propensity Score------------------------
    
    # Create Model
    
    propensity_model_1 <- glm(Assignment ~ X1 + X2 + X3 + X4 + X5, family = binomial() , data = data_combined)
    
    # Predict on the data:
    
    predictions_1 <- predict.glm(propensity_model_1, data_combined, type = 'response')
    
    # Add propensity score column to the dataset:
    
    data_combined$Propensity_score <- predictions_1
    
    
    
   
    mod_match <- matchit(Assignment ~ X1 + X2 + X3 + X4 + X5,
                         method = "nearest",discard = "none", data = data_combined)
    prop_pairs <- mod_match$match.matrix
    
    # Vector of the difference in proportion for each treatment - control pair.
    prop_difference <- c()
    for(i in 1:n_treat){
      prop_difference <- c(abs((data_combined[i, "Propensity_score"]  - data_combined[prop_pairs[i], "Propensity_score"])), prop_difference)
    }
    
    treatment_number <- rev(seq(1, n_treat))
    
    df_propdiff <- data.frame("Treatment" = treatment_number, "Control"= rev(prop_pairs), "Propensity_Difference" = prop_difference)
    
    sorted_propdiff <-  df_propdiff[order(Propensity_Difference),] 
    
    percentile_cutoff <- percentile * n_treat
    
    propensity_pairs <- sorted_propdiff[1:percentile_cutoff, 1:2]
    
    # NOW CALCULATE STANDARDIZED MEAN DIFFERENCE FOR FIRST VARIABLE
    
    x1_treat <- data_combined[propensity_pairs$Treatment, "X1"]
    x1_control <- data_combined[propensity_pairs$Control, "X1"]
    
    x2_treat <- data_combined[propensity_pairs$Treatment, "X2"]
    x2_control <- data_combined[propensity_pairs$Control, "X2"]
    
    x3_treat <- data_combined[propensity_pairs$Treatment, "X3"]
    x3_control <- data_combined[propensity_pairs$Control, "X3"]
    
    x4_treat <- data_combined[propensity_pairs$Treatment, "X4"]
    x4_control <-data_combined[propensity_pairs$Control, "X4"]
    
    x5_treat <- data_combined[propensity_pairs$Treatment, "X5"]
    x5_control <- data_combined[propensity_pairs$Control, "X5"]
    
    smd_x1 <- 100 * abs((mean(x1_treat) - mean(x1_control)) / (sqrt((var(x1_treat) + var(x1_control))/2)))
    smd_x2 <- 100 * abs((mean(x2_treat) - mean(x2_control)) / (sqrt((var(x2_treat) + var(x2_control))/2)))
    smd_x3 <- 100 * abs((mean(x3_treat) - mean(x3_control)) / (sqrt((var(x3_treat) + var(x3_control))/2)))
    smd_x4 <- 100 * abs((mean(x4_treat) - mean(x4_control)) / (sqrt((var(x4_treat) + var(x4_control))/2)))
    smd_x5 <- 100 * abs((mean(x5_treat) - mean(x5_control)) / (sqrt((var(x5_treat) + var(x5_control))/2)))
    
    propensity_asmd <- mean(smd_x1, smd_x2 ,smd_x3 , smd_x4 , smd_x5)
    
    
    
    
    #-----------------------------------------------------Hungarian Method-------------------------------
    library(dplyr)
    library(cluster)
    library(clue)
    
    n <-nrow(data_combined)
    
    data_diss <- data_combined %>%
      select(X1, X2, X3, X4, X5)
    
    dissimilarity_matrix <- as.matrix(daisy(data_diss, stand = TRUE))
    
    hungarian_matrix <- dissimilarity_matrix[1:(n/2),((n/2)+ 1): n]
    
    
    f <-solve_LSAP(hungarian_matrix, maximum = FALSE)
    
    hungarian_control <- as.vector(f) + n_treat
    hungarian_treat <- seq(1, n_treat)
    
    # For each treatment pair calcualte the distance.
    pair_distance <- c()
    for (i in 1:n_treat){
      pair_distance <- c(pair_distance, hungarian_matrix[i, f[i]])
    }
    
    hungarian_df <- data.frame("Treatment" = hungarian_treat, "Control" = hungarian_control, "Distance" = pair_distance)
    
    
    #--- Now remove the pairs that don't make the cutoff:
    
    
    sorted_hungarian_pairs <-  hungarian_df[order(hungarian_df$Distance),] # Sort matrix by distance
    
    hungarian_pairs <- sorted_hungarian_pairs[1:percentile_cutoff, 1:2] # Select only the pairs that make the cutoff
    
    # --- Now compute smd for the remaining pairs
    
    
    hung_x1_treat <- data_combined[hungarian_pairs$Treatment, "X1"]
    hung_x1_control <- data_combined[hungarian_pairs$Control, "X1"]
    
    hung_x2_treat <- data_combined[hungarian_pairs$Treatment, "X2"]
    hung_x2_control <- data_combined[hungarian_pairs$Control, "X2"]
    
    hung_x3_treat <- data_combined[hungarian_pairs$Treatment, "X3"]
    hung_x3_control <- data_combined[hungarian_pairs$Control, "X3"]
    
    hung_x4_treat <- data_combined[hungarian_pairs$Treatment, "X4"]
    hung_x4_control <-data_combined[hungarian_pairs$Control, "X4"]
    
    hung_x5_treat <- data_combined[hungarian_pairs$Treatment, "X5"]
    hung_x5_control <- data_combined[hungarian_pairs$Control, "X5"]
    
    hung_smd_x1 <- 100 * abs((mean(hung_x1_treat) - mean(hung_x1_control)) / (sqrt((var(hung_x1_treat) + var(hung_x1_control))/2)))
    hung_smd_x2 <- 100 * abs((mean(hung_x2_treat) - mean(hung_x2_control)) / (sqrt((var(hung_x2_treat) + var(hung_x2_control))/2)))
    hung_smd_x3 <- 100 * abs((mean(hung_x3_treat) - mean(hung_x3_control)) / (sqrt((var(hung_x3_treat) + var(hung_x3_control))/2)))
    hung_smd_x4 <- 100 * abs((mean(hung_x4_treat) - mean(hung_x4_control)) / (sqrt((var(hung_x4_treat) + var(hung_x4_control))/2)))
    hung_smd_x5 <- 100 * abs((mean(hung_x5_treat) - mean(hung_x5_control)) / (sqrt((var(hung_x5_treat) + var(hung_x5_control))/2)))
    
    hung_asmd <- mean(hung_smd_x1, hung_smd_x2 ,hung_smd_x3 , hung_smd_x4 , hung_smd_x5)
    
    
    
    #----------------------------------------Calculate original ASMD---------------------------------------------
    orig_smd_x1 <- 100 * abs((mean(t_1) - mean(c_1)) / (sqrt((var(t_1) + var(c_1))/2)))
    orig_smd_x2 <- 100 * abs((mean(t_2) - mean(c_2)) / (sqrt((var(t_2) + var(c_2))/2)))
    orig_smd_x3 <- 100 * abs((mean(t_3) - mean(c_3)) / (sqrt((var(t_3) + var(c_3))/2)))
    orig_smd_x4 <- 100 * abs((mean(t_4) - mean(c_4)) / (sqrt((var(t_4) + var(c_4))/2)))
    orig_smd_x5 <- 100 * abs((mean(t_5) - mean(c_5)) / (sqrt((var(t_5) + var(c_5))/2)))
    
    orig_asmd <- mean(orig_smd_x1, orig_smd_x2, orig_smd_x3, orig_smd_x4, orig_smd_x5)
    
    #--------------------------------------Compare hungarian to propensity
    if(hung_asmd < propensity_asmd){
      return(1)
    } else{
      return(0)
    }
    
}


compare_hungtoprop(20, 100, 10, 100, 10, .8)
