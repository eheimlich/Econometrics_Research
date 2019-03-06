library(ggplot2)
library(dplyr)

treatment_age <- smoking_data %>% 
  filter(qsmk == 1) %>%
  summarise(mean(age))

control_age <- smoking_data %>% 
  filter(qsmk == 0) %>%
  summarise(mean(age))

treatment_sex <- smoking_data %>% 
  filter(qsmk == 1) %>%
  summarise(mean(sex))

control_sex <- smoking_data %>% 
  filter(qsmk == 0) %>%
  summarise(mean(sex))



covariate_inbalance <- data.frame("Group" = c("Treatment", "Treatment", "Control", "Control"),
                                  "Variable" = c("Average Age", "Proportion Men", "Average Age", "Proportion Men"), 
                                  "Value" = c(as.numeric(treatment_age), as.numeric(treatment_sex) * 100, 
                                              as.numeric(control_age), as.numeric(control_sex) *100))

ggplot(covariate_inbalance, aes(x = Group, y = Value)) + geom_bar(stat = "identity", aes(fill = Group)) +scale_fill_brewer(palette="Set1")+facet_grid(. ~Variable ) + coord_cartesian(ylim=c(40,55))

aasmd_prop <- data.frame("Type" = c( "After Matching", "Before Matching"), "Value" = c( prop_log_total_aasmd, before_aasmd))
library(ggplot2)
aasmd_prop_plot <- ggplot(aasmd_prop, aes(x = reorder(Type, -Value), Value, fill = Type)) + geom_bar(stat = "identity") + xlab("Before vs After")+ ylab("Average Absolute Standardized Mean Difference")  


aasmd_prop_plot
