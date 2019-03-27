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

head(data_combined)



# Create Model
propensity_model_1 <- glm(Assignment ~ X1 + X2 + X3 + X4 + X5, family = binomial() , data = data_combined)

# Predict on the data:
predictions_1 <- predict.glm(propensity_model_1, data_combined, type = 'response')

# Add propensity score column to the dataset:
data_combined$Propensity_score <- signif(predictions_1, 4)



library(Matching)

#Perform Matching
x_1 <- data_combined$Propensity_score
t_1 <- data_combined$Assignment
rr_1 <- Match(Tr = t_1, X = x_1, M = 1)

# Check covariate balance:
matching_1 <- MatchBalance(Assignment ~ X1 + X2 + X3 + X4 + X5 ,match.out = rr_1 ,data = data_combined, nboots = 500)

for(i in 1:5){
  print(abs(matching_1$AfterMatching[[i]]$sdiff))
}
