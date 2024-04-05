
rm(list = ls())

# libraries
xfun::pkg_attach2(c("tidyverse",
                    "dplyr", 
                    "here",      
                    "MatchIt",  
                    "gridExtra",
                    "ggplot2", 
                    "optmatch",  
                    "cobalt"))   

options(scipen = 999)

# Load ypsps data
ypsps <- read_csv('/Users/brenda/github/Computational-Social-Science-Projects/Project 6/ypsps.csv')
head(ypsps)
nrow(ypsps)

##########################
## Simulation I: NEAREST #
##########################

# Remove post-treatment covariates
filtered_vars <- names(ypsps)[!grepl("1973|1982", names(ypsps))]
baseline_df <- ypsps %>% select(all_of(filtered_vars)) %>% 
  filter(!is.na(parent_GPHighSchoolPlacebo))
baseline_df <- baseline_df[complete.cases(baseline_df), ] 
## generating df with complete cases for function used below 

set.seed(3141)
simulate_and_get_att_nearest <- function(num_simulations) {
  
  ## Initializing empty lists
  att_list <- numeric(num_simulations)
  prop_balanced_covariates <- numeric(num_simulations) 
  summary_df <- list() 
  num_covariates <- list()
  mean_pct_improvement <- numeric(num_simulations)
  match_nearest <- vector("list", length = num_simulations)  
  
  # Simulate samples, perform propensity score matching, and calculate ATT
  for (i in 1:num_simulations) {
    
    num_covariates[[i]] <- sample(5:20, 1) # Randomly choose the number of covariates
    available_cols <- setdiff(names(baseline_df), c('college', 'student_ppnscal'))
    sample_size <- min(num_covariates[[i]] + 2, length(available_cols))
    sampled_cols <- sample(available_cols, sample_size)
    sampled_data <- baseline_df[, c(sampled_cols, "college", "student_ppnscal"), drop = FALSE] 
    # sampled df with sampled cols and treatment  
    
    # Create a matching object with matching using the nearest model 
    match_nearest[[i]] <- matchit(formula = college ~ . - college - student_ppnscal, 
                                  data = sampled_data,             
                                  method = "nearest",
                                  distance = "glm",           
                                  link = "logit",             
                                  discard = "control",
                                  replace = FALSE,
                                  ratio = 2)
    
    ## Summary to store models
    summary_df[[i]] <- summary(match_nearest[[i]]) 
    
    matched_data <- match.data(match_nearest[[i]]) ## creating matched data
    
    # Fit linear model to calculate ATT
    lm_nearest <- lm(student_ppnscal ~ . - student_ppnscal, data = sampled_data)
    ATT_nearest <- coef(lm_nearest)["college"] # Extract ATT
    att_list[i] <- ATT_nearest # Store ATT in the list
    
    ## SMD
    df_all <- as.data.frame(summary_df[[i]]$sum.all)
    df_matched <- as.data.frame(summary_df[[i]]$sum.matched)
    
    # Calculate mean percent improvement in SMD
    if (i > 1) {
      mean_smd_all <- mean(df_all$`Std. Mean Diff.`)
      mean_smd_matched <- mean(df_matched$`Std. Mean Diff.`)
      mean_pct_improvement[i] <- ((mean_smd_matched - mean_smd_all) / mean_smd_matched) * 100
      prop_balanced_covariates[i] <- nrow(df_matched %>% filter(`Std. Mean Diff.` <= 0.1)) / nrow(df_matched)
    }
    
  }
  
  # Create a dataframe to store simulation results
  results_df <- data.frame(
    simulation = 1:num_simulations,
    num_covariates = unlist(num_covariates),
    prop_balanced_covariates = prop_balanced_covariates,
    att = att_list,
    mean_pct_improvement = mean_pct_improvement
    
  )
  
  return(list(match_nearest = match_nearest, results_df = results_df))  
}

# calling function
result_nearest <- simulate_and_get_att_nearest(num_simulations = 1000)


##########################
## Simulation I: OPTIMAL #
##########################

set.seed(3141)

simulate_and_get_att_optimal <- function(num_simulations) {
  
  ## Initializing empty lists
  att_list <- numeric(num_simulations)
  prop_balanced_covariates <- numeric(num_simulations) 
  summary_df <- list() 
  num_covariates <- list()
  mean_pct_improvement <- numeric(num_simulations)
  
  # Simulate samples, perform propensity score matching, and calculate ATT
  for (i in 1:num_simulations) {
    
    num_covariates[[i]] <- sample(1:123, 1) # Randomly choose the number of covariates
    available_cols <- setdiff(names(baseline_df), c('college', 'student_ppnscal'))
    sample_size <- min(num_covariates[[i]] + 2, length(available_cols))
    sampled_cols <- sample(available_cols, sample_size)
    sampled_data <- baseline_df[, c(sampled_cols, "college", "student_ppnscal"), drop = FALSE] # sampled df with sampled cols, treatment and outcome 
    
    # Create a matching object with matching using the optimal model 
    match_optimal <- matchit(formula = college ~ . - college - student_ppnscal, 
                             data = sampled_data,             
                             method = "optimal",
                             distance = "glm",           
                             link = "logit",             
                             discard = "control",
                             replace = FALSE,
                             ratio = 2)
    
    ## Summary to store models
    summary_df[[i]] <- summary(match_optimal) 
    
    matched_data <- match.data(match_optimal) ## creating matched data
    
    # Fit linear model to calculate ATT
    lm_optimal <- lm(student_ppnscal ~ . - student_ppnscal, data = sampled_data)
    ATT_optimal <- coef(lm_optimal)["college"] # Extract ATT
    att_list[i] <- ATT_optimal # Store ATT in the list
    
    ## SMD
    df_all <- as.data.frame(summary_df[[i]]$sum.all)
    df_matched <- as.data.frame(summary_df[[i]]$sum.matched)
    
    # Calculate mean percent improvement in SMD
    if (i > 1) {
      mean_smd_all <- mean(df_all$`Std. Mean Diff.`)
      mean_smd_matched <- mean(df_matched$`Std. Mean Diff.`)
      mean_pct_improvement[i] <- ((mean_smd_matched - mean_smd_all) / mean_smd_matched) * 100
      prop_balanced_covariates[i] <- nrow(df_matched %>% filter(`Std. Mean Diff.` <= 0.1)) / nrow(df_matched)
    }
    
  }
  
  # Create a dataframe to store simulation results
  results_df <- data.frame(
    simulation = 1:num_simulations,
    num_covariates = unlist(num_covariates),
    prop_balanced_covariates = prop_balanced_covariates,
    att = att_list,
    mean_pct_improvement = mean_pct_improvement
  )
  
  return(list(match_optimal = match_optimal, results_df = results_df))  
}

# calling function
result_optimal <- simulate_and_get_att_optimal(num_simulations = 1000)

