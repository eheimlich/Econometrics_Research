library(dplyr)

t_1 <- rnorm(100 ,mean = 100, sd = 10) #Higher
t_2 <- rnorm(100 ,mean = 100, sd = 10) #Lower 
t_3 <- rnorm(100 ,mean = 100, sd = 10)
t_4 <- rnorm(100 ,mean = 100, sd = 10)
t_5 <- rnorm(100 ,mean = 100, sd = 10)
t_assignment<- rep(1, 100)

c_1 <- rnorm(100 ,mean = 95, sd = 10)
c_2 <- rnorm(100 ,mean = 95, sd = 10)
c_3 <- rnorm(100 ,mean = 95, sd = 10)
c_4 <- rnorm(100 ,mean = 95, sd = 10)
c_5 <- rnorm(100 ,mean = 95, sd = 10)
c_assignment <- rep(0,100)


data_t <- data.frame("X1" = t_1, "X2" = t_2, "X3" = t_3, "X4" = t_4, "X5" = t_5, "Assignment" = t_assignment )
data_c <- data.frame("X1" = c_1, "X2" = c_2, "X3" = c_3, "X4" = c_4, "X5" = c_5, "Assignment"= c_assignment)

data_combined <- bind_rows(data_t, data_c)




#---------------------------------------------------------Propensity Score------------------------

# Create Model

propensity_model_1 <- glm(Assignment ~ X1 + X2 + X3 + X4 + X5, family = binomial() , data = data_combined)

# Predict on the data:

predictions_1 <- predict.glm(propensity_model_1, data_combined, type = 'response')

# Add propensity score column to the dataset:

data_combined$Propensity_score <- signif(predictions_1, 4)



#library(Matching)

#Perform Matching
#x_1 <- data_combined$Propensity_score
#t_1 <- data_combined$Assignment
#rr_1 <- Match(Tr = t_1, X = x_1, ties = FALSE)

#rr_1$MatchLoopC


# Check covariate balance:
#matching_1 <- MatchBalance(Assignment ~ X1 + X2 + X3 + X4 + X5 ,match.out = rr_1 ,data = data_combined, nboots = 500 )

library(MatchIt)
mod_match <- matchit(Assignment ~ X1 + X2 + X3 + X4 + X5,
                     method = "nearest",discard = "none", data = data_combined)
prop_pairs <- mod_match$match.matrix

# Vector of the difference in proportion for each treatment - control pair.
prop_difference <- c()
end <- nrow(data_combined)/2
for(i in 1:end){
  prop_difference <- c(abs((data_combined[i, "Propensity_score"]  - data_combined[prop_pairs[i], "Propensity_score"])), prop_difference)
}

treatment_number <- rev(seq(1, end))

df_propdiff <- data.frame("Treatment" = treatment_number, "Control"= rev(prop_pairs), "Propensity_Difference" = prop_difference)

sorted_propdiff <-  df_propdiff[order(Propensity_Difference),] 

percentile_cutoff <- .8 * end

resulting_pairs <- sorted_propdiff[1:percentile_cutoff, 1:2]

# NOW CALCULATE STANDARDIZED MEAN DIFFERENCE FOR FIRST VARIABLE

x1_treat <- data_combined[resulting_pairs$Treatment, "X1"]
x1_control <- data_combined[resulting_pairs$Control, "X1"]

x2_treat <- data_combined[resulting_pairs$Treatment, "X2"]
x2_control <- data_combined[resulting_pairs$Treatment, "X2"]

x3_treat <- data_combined[resulting_pairs$Treatment, "X3"]
x3_control <- data_combined[resulting_pairs$Treatment, "X3"]

x4_treat <- data_combined[resulting_pairs$Treatment, "X4"]
x4_control <-data_combined[resulting_pairs$Treatment, "X4"]

x5_treat <- data_combined[resulting_pairs$Treatment, "X5"]
x5_control <- data_combined[resulting_pairs$Treatment, "X5"]

smd_x1 <- abs()

# diff_prop <- c()
# for(i in 1:5){
#   diff_prop <- c((matching_1$AfterMatching[[i]]$sdiff), diff_prop)
# }
# average_diff_prop <- mean(diff_prop)
# average_diff_prop

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

vec <- c()
for (i in 1:(n/2)){
  vec <- c(vec, hungarian_matrix[i, f[i]])
  vec <- sort(vec, decreasing = FALSE)
}

#return(vec[vec <= quantile(vec, percentile)])
#return (hungarian_matrix)
vec

library(dplyr)
data_combined %>%
  filter(Assignment == 1) %>%
  select(X1) %>%
  summarise(sd = var(X1))



