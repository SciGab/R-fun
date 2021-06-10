### --------------------------------------------------------------------- ###
### --------------- Helper functions for correlations  ------------------ ###
### ----------------------- Version: 10.06.2021 ------------------------- ###
### -------------------------- Gabriela Hofer --------------------------- ###
### --------------------------------------------------------------------- ###


#### 1. significance boundaries for intercorrelations  -------------------------
# based on a MATLAB script of Mathias Benedek

get_rsig <- function(p, n) {
  
  # compute t from desired significance level p
  t <- qnorm((1 - p / 2), 0, 1)
  
  # Solve formula t = r*(sqrt(n-2))/sqrt(1-r^2) for r
  rsig <- t / sqrt(n - 2 + t ^ 2)
  
  # return correlation that would become sig. at desired p-level
  return(rsig)

}

#### 2. function for bootstrapping of (a lot of) correlations ------------------

if (!require(wBoot)) {
  install.packages("wBoot", repos = "http://cran.us.r-project.org")
  library(wBoot)
}

#' A function for computing correlations with bootstrap confidence intervals
#' @xy_names_df A dataframe with two columns where each line contains the names of
#' a variable x to be correlated with a variable y
#' @data A dataframe containing the data the correlation are to be based on
#' @bootsamples The number of bootstrap samples
#' @seed A seed for random number generation to obtain reproducible results    

boot_cor_df <- function(xynames_df, data, bootsamples = 2000, seed = 214035) {
  
  # set seed
  set.seed(seed)
  
  # add 5 empty columns for results
  corrs  <- data.frame(cbind(xynames_df, 
                             rep(NA, nrow(xynames_df)),
                             rep(NA, nrow(xynames_df)),
                             rep(NA, nrow(xynames_df)),
                             rep(NA, nrow(xynames_df)),
                             rep(NA, nrow(xynames_df))))
  
  # name columns of results dataframe
  names(corrs) <- c("x", "y", "Observed_r", "Mean_boot_r", "p", "Perc_CI_lower", "Perc_CI_upper")
  
  corrs <- corrs %>%
    mutate(across(Observed_r:Perc_CI_upper, ~as.numeric(.))) 
  
  # for each desired combination of variables (x and y)...
  for (i in 1:nrow(corrs)){
    x                          <- corrs$x[i]
    y                          <- corrs$y[i]
    
    # ... compute percentile bootstrapping with desired number of bootsamples
    bootinfo                   <- boot.cor.per(data[[x]], data[[y]], R = bootsamples, null.hyp = 0, 
                                               alternative = "two.sided", type = "two-sided")
    
    # save relevant statistics into boot dataframe
    corrs$Observed_r[i]        <- bootinfo$Observed
    corrs$Mean_boot_r[i]       <- bootinfo$Mean
    corrs$p[i]                 <- bootinfo$p.value
    corrs$Perc_CI_lower[i]     <- bootinfo$Confidence.limits[1]
    corrs$Perc_CI_upper[i]     <- bootinfo$Confidence.limits[2]
  }

# return table with results    
return(corrs)
  
}

#### 3. function to compute several "regular" correlations ---------------------

cor_df <- function(xynames_df, data) {
  
  names(xynames_df)   <- c("x", "y", "r")
  sample_corrs        <- xynames_df
  sample_corrs$r      <- as.numeric(sample_corrs$r)

  # calculate correlations and include in df
  for (i in 1:nrow(sample_corrs)){
    x                 <- sample_corrs$x[i]
    y                 <- sample_corrs$y[i]
    data              <- data_whole %>% select(all_of(x), all_of(y)) # all_of to tell dplyr that x and y don't come from data_whole
    sample_corrs$r[i] <- cor(data)[2, 1]
  }

# return table with results    
return(sample_corrs)  
  
}


