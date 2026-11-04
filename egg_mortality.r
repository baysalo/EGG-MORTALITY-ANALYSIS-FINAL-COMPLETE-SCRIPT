# ============================================================
# EGG MORTALITY ANALYSIS
# FINAL COMPLETE SCRIPT - 17 GROUPS
#
# 16 treatment groups + 1 untreated control
# 8 independent Petri-dish replicates per group
# 10 eggs per replicate
#
# Analysis:
# 1. Clean start
# 2. Data import and validation
# 3. Grouped-binomial response
# 4. Descriptive statistics
# 5. Mean bias-reduced binomial logistic regression
# 6. Replicate-level omnibus permutation test
# 7. Exact pairwise permutation post-hoc tests
# 8. Holm adjustment
# 9. Compact letter display
# ============================================================


# ============================================================
# 0. CLEAN START
# ============================================================

rm(list = ls(all.names = TRUE))

while (dev.cur() > 1) {
  dev.off()
}

invisible(gc())

options(stringsAsFactors = FALSE)

cat(
  "\n====================================================\n",
  " EGG MORTALITY ANALYSIS - 17 GROUPS\n",
  " CLEAN R SESSION INITIALIZED\n",
  "====================================================\n",
  sep = ""
)


# ============================================================
# 1. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "brglm2",
  "readxl",
  "dplyr",
  "multcompView"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  
  cat("\nInstalling missing packages:\n")
  print(missing_packages)
  
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

still_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing) > 0) {
  
  stop(
    paste(
      "The following packages could not be loaded:",
      paste(still_missing, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages(library(brglm2))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(multcompView))

cat("\nAll required packages are available.\n")


# ============================================================
# 2. IMPORT DATA
# ============================================================

cat(
  "\nSelect the Excel file containing the 'yumurta' sheet.\n"
)

dat <- readxl::read_excel(
  file.choose(),
  sheet = "yumurta"
)

cat("\nExcel file imported successfully.\n")


# ============================================================
# 3. CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "Treat",
  "Replication",
  "%mortality"
)

if (!all(required_columns %in% names(dat))) {
  
  stop(
    paste(
      "Required columns not found.",
      "Expected:",
      paste(required_columns, collapse = ", "),
      "Found:",
      paste(names(dat), collapse = ", ")
    )
  )
}

if (anyNA(dat[, required_columns])) {
  
  stop(
    "Missing values detected in required columns."
  )
}


# ============================================================
# 4. CONVERT MORTALITY TO DEAD / ALIVE COUNTS
# ============================================================

dead_unrounded <-
  (dat$`%mortality` / 100) * 10

if (
  any(
    abs(
      dead_unrounded -
      round(dead_unrounded)
    ) > 1e-8
  )
) {
  
  stop(
    "At least one mortality percentage cannot be converted ",
    "to an integer number of deaths out of 10 eggs."
  )
}

dat <- dat %>%
  dplyr::mutate(
    Dead = as.integer(
      round((`%mortality` / 100) * 10)
    ),
    Alive = as.integer(10 - Dead),
    Total = 10L
  )


# ============================================================
# 5. VALIDATE COUNTS
# ============================================================

if (any(dat$Dead < 0)) {
  stop("Negative death count detected.")
}

if (any(dat$Dead > 10)) {
  stop("Death count greater than 10 detected.")
}

if (any(dat$Alive < 0)) {
  stop("Negative survival count detected.")
}

if (any(dat$Dead + dat$Alive != dat$Total)) {
  stop("Dead + Alive does not equal Total.")
}


# ============================================================
# 6. DEFINE TREATMENT ORDER
# ============================================================

if (is.numeric(dat$Treat)) {
  
  treatment_levels <- as.character(
    sort(unique(dat$Treat))
  )
  
} else {
  
  treatment_levels <- unique(
    as.character(dat$Treat)
  )
}

dat$Treat <- factor(
  dat$Treat,
  levels = treatment_levels
)


# ============================================================
# 7. CHECK EXPERIMENTAL DESIGN
# ============================================================

cat(
  "\n================ EXPERIMENTAL DESIGN ================\n"
)

cat(
  "Number of groups =",
  length(treatment_levels),
  "\n"
)

if (length(treatment_levels) != 17) {
  
  stop(
    "Exactly 17 groups are required: ",
    "16 treatments + untreated control."
  )
}

rep_check <- table(dat$Treat)

cat("\nReplicates per group:\n")
print(rep_check)

if (any(rep_check != 8)) {
  
  stop(
    "Each group must contain exactly ",
    "8 independent replicates."
  )
}

if (any(dat$Total != 10)) {
  
  stop(
    "Each replicate must contain exactly 10 eggs."
  )
}


# ============================================================
# 8. CHECK UNIQUE REPLICATE IDENTIFIERS
# ============================================================

duplicate_rep <- dat %>%
  dplyr::count(
    Treat,
    Replication,
    name = "n"
  ) %>%
  dplyr::filter(n != 1)

if (nrow(duplicate_rep) > 0) {
  
  stop(
    "Duplicate Treatment ?? Replication identifiers detected."
  )
}


# ============================================================
# 9. DESCRIPTIVE RESULTS
# ============================================================

observed <- dat %>%
  dplyr::group_by(Treat) %>%
  dplyr::summarise(
    Replicates = dplyr::n(),
    Dead = sum(Dead),
    Alive = sum(Alive),
    Total = sum(Total),
    Mortality_percent = 100 * Dead / Total,
    .groups = "drop"
  )

cat(
  "\n================ OBSERVED MORTALITY ================\n"
)

print(observed)

cat(
  "\nTotal dead eggs =",
  sum(dat$Dead),
  "\n"
)

cat(
  "Total surviving eggs =",
  sum(dat$Alive),
  "\n"
)

cat(
  "Total eggs =",
  sum(dat$Total),
  "\n"
)


# ============================================================
# 10. CHECK UNTREATED CONTROL
# ============================================================

# In this dataset Treatment 17 is the untreated control.

control_data <- dat %>%
  dplyr::filter(Treat == "17")

cat(
  "\n================ UNTREATED CONTROL ================\n"
)

cat(
  "Control group = Treatment 17\n"
)

cat(
  "Control replicates =",
  nrow(control_data),
  "\n"
)

cat(
  "Control deaths =",
  sum(control_data$Dead),
  "/",
  sum(control_data$Total),
  "\n"
)

cat(
  "Observed control mortality (%) =",
  100 *
    sum(control_data$Dead) /
    sum(control_data$Total),
  "\n"
)


# ============================================================
# 11. BIAS-REDUCED BINOMIAL LOGISTIC REGRESSION
# ============================================================

fit_br <- glm(
  cbind(Dead, Alive) ~ Treat,
  family = binomial(link = "logit"),
  data = dat,
  method = brglm2::brglmFit,
  type = "AS_mean"
)

cat(
  "\n================ BIAS-REDUCED MODEL ================\n"
)

print(summary(fit_br))


# ============================================================
# 12. MODEL-BASED PREDICTIONS
# ============================================================

pred <- data.frame(
  Treat = factor(
    treatment_levels,
    levels = treatment_levels
  )
)

pred$Predicted_probability <- predict(
  fit_br,
  newdata = pred,
  type = "response"
)

pred$Predicted_percent <-
  100 * pred$Predicted_probability

results <- observed %>%
  dplyr::left_join(
    pred,
    by = "Treat"
  )

cat(
  "\n================ OBSERVED AND PREDICTED ================\n"
)

print(
  results[
    ,
    c(
      "Treat",
      "Dead",
      "Total",
      "Mortality_percent",
      "Predicted_percent"
    )
  ]
)


# ============================================================
# 13. LIKELIHOOD FUNCTION
# ============================================================

# Binomial log-likelihood excluding combinatorial constants.
# The constants cancel in likelihood-ratio comparisons.

binom_ll <- function(y, n) {
  
  if (y == 0 || y == n) {
    return(0)
  }
  
  p <- y / n
  
  y * log(p) +
    (n - y) * log1p(-p)
}


# ============================================================
# 14. OMNIBUS LR STATISTIC FUNCTION
# ============================================================

omnibus_lr <- function(data) {
  
  group_summary <- data %>%
    dplyr::group_by(Treat) %>%
    dplyr::summarise(
      Dead = sum(Dead),
      Total = sum(Total),
      .groups = "drop"
    )
  
  total_dead <- sum(group_summary$Dead)
  total_n <- sum(group_summary$Total)
  
  ll_null <- binom_ll(
    total_dead,
    total_n
  )
  
  ll_full <- sum(
    mapply(
      binom_ll,
      group_summary$Dead,
      group_summary$Total
    )
  )
  
  2 * (ll_full - ll_null)
}


# ============================================================
# 15. OBSERVED OMNIBUS LR
# ============================================================

LR_obs <- omnibus_lr(dat)

cat(
  "\n================ OBSERVED OMNIBUS LR ================\n"
)

cat(
  "Observed LR/deviance =",
  LR_obs,
  "\n"
)


# ============================================================
# 16. OMNIBUS REPLICATE-LEVEL PERMUTATION TEST
# ============================================================

# Complete Petri-dish outcomes are retained.
# Only treatment labels are permuted.
#
# 17 groups ?? 8 replicates = 136 independent units.

set.seed(20261029)

B <- 100000L

LR_perm <- numeric(B)

cat(
  "\n================ OMNIBUS PERMUTATION TEST ================\n"
)

cat(
  "Running",
  B,
  "permutations...\n"
)

pb <- txtProgressBar(
  min = 0,
  max = B,
  style = 3
)

for (b in seq_len(B)) {
  
  dat_perm <- dat
  
  dat_perm$Treat <- sample(
    dat$Treat,
    size = nrow(dat),
    replace = FALSE
  )
  
  LR_perm[b] <- omnibus_lr(
    dat_perm
  )
  
  if (b %% 1000L == 0L) {
    
    setTxtProgressBar(
      pb,
      b
    )
  }
}

close(pb)


# ============================================================
# 17. OMNIBUS PERMUTATION P VALUE
# ============================================================

exceed <- sum(
  LR_perm >= LR_obs - 1e-12
)

p_perm <-
  (exceed + 1) /
  (B + 1)

mc_se <- sqrt(
  p_perm *
    (1 - p_perm) /
    (B + 1)
)

alpha <- 0.05

cat(
  "\n================ FINAL OMNIBUS RESULT ================\n"
)

cat(
  "Observed LR/deviance =",
  round(LR_obs, 5),
  "\n"
)

cat(
  "Number of permutations =",
  B,
  "\n"
)

cat(
  "Exceedances =",
  exceed,
  "\n"
)

cat(
  "Permutation P =",
  format.pval(
    p_perm,
    digits = 8
  ),
  "\n"
)

cat(
  "Monte Carlo SE =",
  round(mc_se, 6),
  "\n"
)


omnibus_significant <-
  p_perm < alpha

cat(
  "Omnibus significant =",
  omnibus_significant,
  "\n"
)


# ============================================================
# 18. EXACT PAIRWISE LR FUNCTION
# ============================================================

pair_lr <- function(
    y1,
    n1,
    y2,
    n2
) {
  
  total_y <- y1 + y2
  total_n <- n1 + n2
  
  ll_null <- binom_ll(
    total_y,
    total_n
  )
  
  ll_alt <-
    binom_ll(y1, n1) +
    binom_ll(y2, n2)
  
  LR <- 2 * (ll_alt - ll_null)
  
  if (LR < 0 && abs(LR) < 1e-10) {
    LR <- 0
  }
  
  LR
}


# ============================================================
# 19. EXACT PAIRWISE PERMUTATION FUNCTION
# ============================================================

exact_pair_test <- function(
    treatment1,
    treatment2,
    data
) {
  
  d <- data[
    data$Treat %in%
      c(treatment1, treatment2),
  ]
  
  d$Treat <- droplevels(d$Treat)
  
  n_by_group <- table(d$Treat)
  
  if (
    length(n_by_group) != 2 ||
    any(n_by_group != 8)
  ) {
    
    stop(
      "Each pair must contain 8 + 8 replicates."
    )
  }
  
  d1 <- d[d$Treat == treatment1, ]
  d2 <- d[d$Treat == treatment2, ]
  
  dead1 <- sum(d1$Dead)
  total1 <- sum(d1$Total)
  
  dead2 <- sum(d2$Dead)
  total2 <- sum(d2$Total)
  
  LR_observed <- pair_lr(
    dead1,
    total1,
    dead2,
    total2
  )
  
  # All possible allocations:
  # C(16,8) = 12,870
  
  combinations <- combn(
    nrow(d),
    8
  )
  
  n_allocations <- ncol(
    combinations
  )
  
  if (
    n_allocations != choose(16, 8)
  ) {
    
    stop(
      "Unexpected number of exact allocations."
    )
  }
  
  all_dead <- sum(d$Dead)
  all_total <- sum(d$Total)
  
  dead_matrix <- matrix(
    d$Dead[combinations],
    nrow = 8
  )
  
  total_matrix <- matrix(
    d$Total[combinations],
    nrow = 8
  )
  
  perm_dead1 <- colSums(
    dead_matrix
  )
  
  perm_total1 <- colSums(
    total_matrix
  )
  
  perm_dead2 <-
    all_dead -
    perm_dead1
  
  perm_total2 <-
    all_total -
    perm_total1
  
  LR_exact <- numeric(
    n_allocations
  )
  
  for (k in seq_len(n_allocations)) {
    
    LR_exact[k] <- pair_lr(
      perm_dead1[k],
      perm_total1[k],
      perm_dead2[k],
      perm_total2[k]
    )
  }
  
  extreme <- sum(
    LR_exact >=
      LR_observed - 1e-12
  )
  
  exact_p <-
    extreme /
    n_allocations
  
  data.frame(
    Treatment1 = as.character(treatment1),
    Treatment2 = as.character(treatment2),
    
    Dead1 = dead1,
    Total1 = total1,
    Mortality1_percent =
      100 * dead1 / total1,
    
    Dead2 = dead2,
    Total2 = total2,
    Mortality2_percent =
      100 * dead2 / total2,
    
    LR = LR_observed,
    
    Extreme_allocations =
      extreme,
    
    Total_allocations =
      n_allocations,
    
    Exact_P =
      exact_p,
    
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 20. GENERATE ALL 136 PAIRS
# ============================================================

treatment_pairs <- combn(
  treatment_levels,
  2,
  simplify = FALSE
)

number_pairwise_tests <-
  length(treatment_pairs)

cat(
  "\n================ POST-HOC STRUCTURE ================\n"
)

cat(
  "Number of groups = 17\n"
)

cat(
  "Number of pairwise comparisons =",
  number_pairwise_tests,
  "\n"
)

cat(
  "Exact allocations per comparison =",
  choose(16, 8),
  "\n"
)

if (
  number_pairwise_tests != 136
) {
  
  stop(
    "Expected exactly 136 pairwise comparisons."
  )
}


# ============================================================
# 21. RUN ALL EXACT PAIRWISE TESTS
# ============================================================

if (!omnibus_significant) {
  
  cat(
    "\nOmnibus test was not significant.\n",
    "Post-hoc analysis will not be interpreted.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nRunning 136 exact pairwise tests...\n"
  )
  
  pairwise_results_list <- vector(
    "list",
    number_pairwise_tests
  )
  
  pb_exact <- txtProgressBar(
    min = 0,
    max = number_pairwise_tests,
    style = 3
  )
  
  for (
    i in seq_along(treatment_pairs)
  ) {
    
    pair <- treatment_pairs[[i]]
    
    pairwise_results_list[[i]] <-
      exact_pair_test(
        treatment1 = pair[1],
        treatment2 = pair[2],
        data = dat
      )
    
    setTxtProgressBar(
      pb_exact,
      i
    )
  }
  
  close(pb_exact)
  
  
  # ============================================================
  # 22. COMBINE PAIRWISE RESULTS
  # ============================================================
  
  pairwise_results <- do.call(
    rbind,
    pairwise_results_list
  )
  
  rownames(pairwise_results) <- NULL
  
  
  # ============================================================
  # 23. HOLM ADJUSTMENT
  # ============================================================
  
  pairwise_results$Holm_P <- p.adjust(
    pairwise_results$Exact_P,
    method = "holm"
  )
  
  pairwise_results$Significant <- ifelse(
    pairwise_results$Holm_P < alpha,
    "Yes",
    "No"
  )
  
  pairwise_results_sorted <-
    pairwise_results[
      order(
        pairwise_results$Holm_P,
        pairwise_results$Exact_P,
        pairwise_results$Treatment1,
        pairwise_results$Treatment2
      ),
    ]
  
  rownames(
    pairwise_results_sorted
  ) <- NULL
  
  
  # ============================================================
  # 24. SIGNIFICANT PAIRS
  # ============================================================
  
  significant_pairs <-
    pairwise_results_sorted[
      pairwise_results_sorted$Holm_P <
        alpha,
    ]
  
  n_significant <-
    nrow(significant_pairs)
  
  cat(
    "\n================ POST-HOC RESULT ================\n"
  )
  
  cat(
    "Pairwise comparisons =",
    nrow(pairwise_results),
    "\n"
  )
  
  cat(
    "Significant comparisons after Holm =",
    n_significant,
    "\n"
  )
  
  cat(
    "Minimum unadjusted exact P =",
    min(pairwise_results$Exact_P),
    "\n"
  )
  
  cat(
    "Minimum Holm-adjusted P =",
    min(pairwise_results$Holm_P),
    "\n"
  )
  
  
  # ============================================================
  # 25. COMPACT LETTER DISPLAY
  # ============================================================
  
  p_for_letters <-
    pairwise_results$Holm_P
  
  names(p_for_letters) <- paste(
    pairwise_results$Treatment1,
    pairwise_results$Treatment2,
    sep = "-"
  )
  
  letter_object <-
    multcompView::multcompLetters(
      p_for_letters,
      threshold = alpha
    )
  
  letters_df <- data.frame(
    Treat =
      names(letter_object$Letters),
    
    Letter =
      unname(letter_object$Letters),
    
    stringsAsFactors = FALSE
  )
  
  letters_df$Treat <- factor(
    letters_df$Treat,
    levels = treatment_levels
  )
  
  letters_df <- letters_df[
    order(letters_df$Treat),
  ]
  
  letters_df$Treat <-
    as.character(letters_df$Treat)
  
  
  # ============================================================
  # 26. FINAL LETTER TABLE
  # ============================================================
  
  observed_for_letters <- observed %>%
    dplyr::mutate(
      Treat =
        as.character(Treat)
    ) %>%
    dplyr::select(
      Treat,
      Replicates,
      Dead,
      Alive,
      Total,
      Mortality_percent
    )
  
  final_letters <- dplyr::left_join(
    observed_for_letters,
    letters_df,
    by = "Treat"
  )
  
  final_letters$Treat <- factor(
    final_letters$Treat,
    levels = treatment_levels
  )
  
  final_letters <- final_letters[
    order(final_letters$Treat),
  ]
  
  final_letters$Treat <-
    as.character(final_letters$Treat)
  
  rownames(final_letters) <- NULL
  
  
  # ============================================================
  # 27. LETTER CONSISTENCY CHECK
  # ============================================================
  
  share_letter <- function(x, y) {
    
    if (is.na(x) || is.na(y)) {
      return(FALSE)
    }
    
    lx <- strsplit(x, "")[[1]]
    ly <- strsplit(y, "")[[1]]
    
    length(
      intersect(lx, ly)
    ) > 0
  }
  
  letter_lookup <- setNames(
    final_letters$Letter,
    final_letters$Treat
  )
  
  letter_check <- pairwise_results
  
  letter_check$Shared_letter <- mapply(
    function(t1, t2) {
      
      share_letter(
        letter_lookup[
          as.character(t1)
        ],
        letter_lookup[
          as.character(t2)
        ]
      )
    },
    
    letter_check$Treatment1,
    letter_check$Treatment2
  )
  
  letter_check$Correct <- ifelse(
    letter_check$Holm_P < alpha,
    !letter_check$Shared_letter,
    letter_check$Shared_letter
  )
  
  if (all(letter_check$Correct)) {
    
    cat(
      "\nLetter-display consistency check: PASSED\n"
    )
    
  } else {
    
    warning(
      "Letter-display consistency check failed."
    )
  }
  
  
  # ============================================================
  # 28. SAVE OUTPUT FILES
  # ============================================================
  
  write.csv(
    pairwise_results_sorted,
    "Egg_17groups_exact_pairwise_Holm.csv",
    row.names = FALSE
  )
  
  write.csv(
    significant_pairs,
    "Egg_17groups_significant_pairs.csv",
    row.names = FALSE
  )
  
  write.csv(
    final_letters,
    "Egg_17groups_Holm_letters.csv",
    row.names = FALSE
  )
  
  write.csv(
    letter_check,
    "Egg_17groups_letter_check.csv",
    row.names = FALSE
  )
  
  
  # ============================================================
  # 29. FINAL SUMMARY
  # ============================================================
  
  cat(
    "\n====================================================\n",
    " FINAL ANALYSIS SUMMARY\n",
    "====================================================\n",
    sep = ""
  )
  
  cat(
    "\nEXPERIMENTAL DESIGN\n"
  )
  
  cat(
    "Groups = 17\n"
  )
  
  cat(
    "Replicates per group = 8\n"
  )
  
  cat(
    "Eggs per replicate = 10\n"
  )
  
  cat(
    "Total independent replicates =",
    nrow(dat),
    "\n"
  )
  
  cat(
    "Total eggs =",
    sum(dat$Total),
    "\n"
  )
  
  cat(
    "Untreated control = Treatment 17\n"
  )
  
  cat(
    "Control mortality =",
    sum(control_data$Dead),
    "/",
    sum(control_data$Total),
    "\n"
  )
  
  cat(
    "\nOMNIBUS PERMUTATION TEST\n"
  )
  
  cat(
    "Observed LR/deviance =",
    round(LR_obs, 5),
    "\n"
  )
  
  cat(
    "Number of permutations =",
    B,
    "\n"
  )
  
  cat(
    "Exceedances =",
    exceed,
    "\n"
  )
  
  cat(
    "Permutation P =",
    format.pval(
      p_perm,
      digits = 8
    ),
    "\n"
  )
  
  cat(
    "Monte Carlo SE =",
    round(mc_se, 6),
    "\n"
  )
  
  cat(
    "\nEXACT POST-HOC TESTS\n"
  )
  
  cat(
    "Pairwise comparisons =",
    nrow(pairwise_results),
    "\n"
  )
  
  cat(
    "Exact allocations per comparison =",
    choose(16, 8),
    "\n"
  )
  
  cat(
    "Total allocations evaluated =",
    nrow(pairwise_results) *
      choose(16, 8),
    "\n"
  )
  
  cat(
    "Multiplicity adjustment = Holm\n"
  )
  
  cat(
    "Family-wise alpha =",
    alpha,
    "\n"
  )
  
  cat(
    "Significant pairwise comparisons =",
    n_significant,
    "\n"
  )
  
  cat(
    "\nFINAL TREATMENT GROUPS\n"
  )
  
  print(
    final_letters[
      ,
      c(
        "Treat",
        "Mortality_percent",
        "Letter"
      )
    ],
    row.names = FALSE
  )
}


# ============================================================
# 30. SAVE DESCRIPTIVE AND OMNIBUS RESULTS
# ============================================================

write.csv(
  results,
  "Egg_17groups_descriptive_results.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    LR_observed = LR_obs,
    Permutations = B,
    Exceedances = exceed,
    Permutation_P = p_perm,
    Monte_Carlo_SE = mc_se,
    Significant = omnibus_significant
  ),
  "Egg_17groups_omnibus_permutation.csv",
  row.names = FALSE
)


# ============================================================
# 31. SESSION INFORMATION
# ============================================================

cat(
  "\n================ SESSION INFORMATION ================\n"
)

print(
  sessionInfo()
)


# ============================================================
# 32. COMPLETED
# ============================================================

cat(
  "\n====================================================\n",
  " ANALYSIS COMPLETED SUCCESSFULLY\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nOutput directory:\n",
  normalizePath(getwd()),
  "\n",
  sep = ""
)
