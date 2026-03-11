# Libs --------------------------------------------------
library(tidyverse)
library(poLCA)
library(gt)
library(tidyLPA)

df <- df %>% dplyr::select(-c('year', 'age', 'reason','sui')) %>% 
                    mutate(occupation = ifelse(occupation == '12', '1', '0') %>% as.factor())


# relevel the data for poLCA, all levels have to start with 1
df[] <- lapply(df, function(x) {
  if (is.factor(x) && min(as.integer(levels(x))) != 1){
    levels(x) <- as.integer(levels(x)) - min(as.integer(levels(x))) + 1} 
  x })


# LCA --------------------------------------------------

f <- cbind(gender, age_cat, religion, race, occupation, marital_stat, education, lifetime_psyk_hos, past_sui_attempt,
           sui_thoughts, self_injury, psyk_dis, past_illness, alcohol_drug_cons,
           anger, sleep_prob, social_iso, sad_weary, humiliated) ~ 1 



classes <- 2:10  

set.seed(2025) 

results <- lapply(classes, function(c) {
  poLCA::poLCA(f, df,
               nclass = c,
               verbose = FALSE,
               nrep = 100,
               maxiter = 3000)})

names(results) <- paste0("class_", classes)


## Metrics ---------------------------------------------

### Entropy --------------------------------------------
machine_tolerance <- sqrt(.Machine$double.eps)
entropy.R2 <- function(fit) {
  entropy <- function(p) {
    p <- p[p > machine_tolerance] # since Lim_{p->0} p log(p) = 0
    sum(-p * log(p))
  }
  error_prior <- entropy(fit$P) # Class proportions
  error_post <- mean(apply(fit$posterior, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  R2_entropy
}


### LMR test (Lo–Mendell–Rubin) -----------------------
LMR_values <- c("-")  # blank for first row (no comparison)

for (i in 2:length(results)) {
  M1 <- results[[i - 1]]  
  M2 <- results[[i]]      
  
  Nobs <- if (!is.null(M1$Nobs)) M1$Nobs else M1$N
  
  LMR <- calc_lrt(
    Nobs,
    M1$llik,
    M1$npar,
    length(M1$P),
    M2$llik,
    M2$npar,
    length(M2$P))
  
  LMR_values <- c(LMR_values, round(LMR[[4]],2))}

metrics <- data.frame(
  classes = classes,
  BIC = sapply(results, function(m) round(m$bic)),
  cAIC = sapply(results, function(m) round((-2*m$llik)+ m$npar*(1+log(m$N)))),  
  entropy = sapply(results, function(m) (entropy.R2(m) %>% round(2))),
  smallest_class_expected = sapply(results, function(m) round(min(m$P) * m$N)))


metrics$LMR_p <- LMR_values


### predicted probability of membership per class --------------------
for (k in 1:max(classes)) {
  metrics[[paste0("avg_post_C", k)]] <- sapply(results, function(m) {
    if (ncol(m$posterior) >= k) {
      round(mean(m$posterior[, k]), 2)
    } else {
      ""
    }
  })
}



metric_table <- metrics %>% gt() %>% 
  cols_label(
    classes = "Class",
    entropy = "Entropy",
    smallest_class_expected = "Smallest Class Count",
    avg_post_C1 = "1",
    avg_post_C2 = "2",
    avg_post_C3 = "3",
    avg_post_C4 = "4",
    avg_post_C5 = "5",
    avg_post_C6 = "6",
    avg_post_C7 = "7",
    avg_post_C8 = "8",
    avg_post_C9 = "9",
    avg_post_C10 = "10") %>% 
  tab_spanner(label = "Average Posterior Probability per Class", columns = starts_with("avg_post")) %>% 
  tab_footnote(footnote = "BIC: Bayesian Information Criterion",
               locations = cells_column_labels(columns = BIC)) %>% 
  tab_footnote(footnote = "cAIC: Consistent Akaike Information Criterion",
               locations = cells_column_labels(columns = cAIC)) %>% 
  tab_footnote(footnote = "LMR_p: Lo-Mendell-Rubin likelihood ratio test p value",
               locations = cells_column_labels(columns = LMR_p)) %>% 
  cols_align(align = "center",
             columns = everything())