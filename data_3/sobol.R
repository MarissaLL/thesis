library(tidyverse)
library(sensitivity)


#############
# FUNCTIONS #
#############


make_unif_dist <- function(num, dataset, model_runs){
  lower_bound <- pull(dataset[num,2])
  upper_bound <- pull(dataset[num,3])
  type <- pull(dataset[num,4])
  if (type == "discrete"){
    samples <- matrix(sample(x = seq(from = lower_bound, to = upper_bound, by = 1),
                             size = model_runs,
                             replace=T))
  } 
  if (type == "continuous"){
    samples <- matrix(runif(model_runs, min = lower_bound, max = upper_bound))
  }
  return(samples)
}



param_vals <- tibble(name = c("max_allowed_pairings_m", 
                  "age_first_breeding_f", "age_first_breeding_m", "maximum_breeding_interval", "prop_multimating_m",
                  "survival_rate", "loc_prop_females", "scale_prop_females", "obs_mean_bi", "obs_sd_bi"),
       lower_bound = c(1,3,7,2,0,0.89, 20,5,2.2,0.1),
       upper_bound = c(6,7,13,4,1,0.99,80,35,5,2),
       type = c("discrete","discrete","discrete","discrete", "continuous", "continuous","continuous", "continuous","continuous", "continuous"))


 # Number of samples
n_model_runs = 100





param_samples_list_m1 <- lapply(1:nrow(param_vals), make_unif_dist, dataset = param_vals, model_runs = n_model_runs)
param_samples_list_m2 <- lapply(1:nrow(param_vals), make_unif_dist, dataset = param_vals, model_runs = n_model_runs)


m1 <- matrix(unlist(param_samples_list_m1), ncol=length(param_samples_list_m1), byrow=FALSE)
m2 <- matrix(unlist(param_samples_list_m2), ncol=length(param_samples_list_m2), byrow=FALSE) 


# Perform the cross-sampling
sobol_seq = sobolmartinez(model = NULL, X1 = data.frame(m1), X2 = data.frame(m2), nboot = 100)

to_be_run <- sobol_seq$X 
names(to_be_run) <- param_vals$name
parameter_settings_file <- to_be_run %>%
  rownames_to_column(var = "iter")



## DON'T uncomment unless intending to re-run simulations
write_delim(parameter_settings_file, "sobol_params_to_run.csv", delim = ",")
save(sobol_seq, file = "sobol_seq.RData")




# Run the model and then load the results here
result <- read_delim("time_to_carrying_capacity.txt", delim = ",") %>% 
  mutate(X2 = case_when(year_passed == Inf ~ 2300,
                        TRUE ~ year_passed)) %>% 
  pull(X2)

result <- read_delim("lifetime_offspring.txt", delim = ",") %>% 
  pull(mean)

load("sobol_seq.RData")

tell(sobol_seq, result)
first_order <- tibble(sobol_seq$S) %>% rownames_to_column(var = "parameter") %>% mutate(index = "First order")
total_indices <- tibble(sobol_seq$T) %>% rownames_to_column(var = "parameter") %>% mutate(index = "Total")
 
# ggplot(sobol_seq, ylim = c(0, 1))
# 
# 
sobol_to_plot <- bind_rows(first_order, total_indices) %>% 
  rename(std_err = `std. error`,
         min_ci = `min. c.i.`,
         max_ci = `max. c.i.`)
  
sobol_to_plot$parameter <- factor(sobol_to_plot$parameter, levels = c(1:10), labels = c("max_allowed_pairings_m", "age_first_breeding_f", "age_first_breeding_m", 
                                                                                        "maximum_breeding_interval", "prop_multimating_m", "survival_rate", "loc_prop_females", 
                                                                                        "scale_prop_females", "obs_mean_bi", "obs_sd_bi"))




# There is a known effect where zero value sobol indicies can appear to be negative, 
# but they should be treated as zero, especially if the CI includes 0.
# Tidy this in the data to plot
sobol_to_plot <- sobol_to_plot %>% 
  mutate(original = case_when(original < 0 ~ 0,
                              TRUE ~ original)) %>% 
  filter(parameter != "max_allowed_pairings_m") # remove defunct parameter that doesn't change the model


ggplot(sobol_to_plot, aes(x = reorder(parameter, 
                                      as.numeric(original)), y = original, colour = index ))+
  geom_point(position = position_dodge(0.4), size = 2)+
  geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width = 0.5, position = position_dodge(0.4), 
                size = 1.2)+
  labs(x = "Parameter",  y = "Sobol index", colour = "Effect") +
  theme_classic()+ 
  theme(legend.position = "top",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))+
  scale_colour_manual(values = c("#440154FF","#27AD81FF"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # scale_x_discrete(labels = function(labels) {
  #   fixedLabels <- c()
  #   for (l in 1:length(labels)) {
  #     fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
  #   }
  #   return(fixedLabels)
  # })





#The sobol2002 function and others use an approach that takes two input matrices 𝐴 and 𝐵, both of with 𝑁 rows and 𝑘 columns. Here, 𝑘 is the number of model parameters (factors) and 𝑁 is the number of model evaluations.
# The two matrices are combined to obtain random parameter sets that differ only in one parameter. This is achieved by replacing the 𝑖-th column in 𝐴 with the 𝑖-th column in 𝐵. This gives 𝑘 matrices 𝐴(𝑖)𝐵.
# The model is then evaluated for the 𝑘𝑁 rows from the matrices 𝐴(1)𝐵…𝐴(𝑘)𝐵. Intuitively, this allows to assess the variation of the model measurements when varying each of the factors while keeping the others fixed.
# 
# n <- 1000
# X1 <- data.frame(matrix(runif(8 * n), nrow = n))
# X2 <- data.frame(matrix(runif(8 * n), nrow = n))
# 
# #Sensitivity analysis
# x <- sobol2002(model = NULL, X1, X2, nboot = 100)
# y <- ishigami.fun(x$X)
# tell(x,y)
# print(x)
# plot(x)


# ## Doesn't handle discrete parameter values
# # Multiply by the range and then add the minimum value - I'm guessing this is so that each matrix contains uniformly-distributed samples within the bounds of the parameter range
# if (includes_discrete == FALSE){
#   m1 = matrix(runif(n_params*n_model_runs), nrow=n_model_runs)
#   m1 = sweep(m1, MARGIN=2, param_vals$upper_bound - param_vals$lower_bound, '*') 
#   m1 = sweep(m1, MARGIN=2, param_vals$lower_bound, '+');
#   
#   m2 = matrix(runif(n_params*n_model_runs), nrow=n_model_runs)
#   m2 = sweep(m2, MARGIN=2,  param_vals$upper_bound - param_vals$lower_bound, '*')
#   m2 = sweep(m2, MARGIN=2, param_vals$lower_bound, '+')
# }
# elif(includes_discrete == TRUE){
#   
# }