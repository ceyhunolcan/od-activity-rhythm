# Primary regression analysis: M1..M4 progressive adjustment with NHANES survey
# weights and PSU-cluster-robust SEs, plus BH-FDR control within outcome families.
#
# Inputs:
#   analytic_full.csv         from stage1
#   fragmentation_features.csv from stage8
#   activity_summary.csv      per-person summaries (mean_mims, mvpa_min, IS, IV,
#                             ASTP, RA, M10, L5, total_sleep_min, WASO,
#                             sleep_efficiency) -- produced by stage6 (not in
#                             this deposit; build from PAXMIN_H per-day file)
#
# Outputs:
#   table2_Aim1.csv                       -- volume / intensity
#   table3_Aim2_3.csv                     -- rhythm / sleep / fragmentation
#   tableS1_M4_full_coefficients.csv      -- M4 with all covariates
#   tableS2_M4_with_BHFDR.csv             -- M4 primary outcomes, q-values
#   tableS5_MICE_pooled.csv               -- MICE-pooled M4 (mitools::MIcombine)

suppressPackageStartupMessages({
    library(survey)
    library(mice)
    library(mitools)
    library(dplyr)
    library(readr)
})

set.seed(42)


# --- load and merge ---------------------------------------------------------

df_main <- read_csv("analytic_full.csv",          show_col_types = FALSE)
df_frag <- read_csv("fragmentation_features.csv", show_col_types = FALSE)
df <- df_main %>% left_join(df_frag, by = "SEQN")

# activity_summary.csv supplies the volume / rhythm / sleep outcomes
if (file.exists("activity_summary.csv")) {
    df_act <- read_csv("activity_summary.csv", show_col_types = FALSE)
    df <- df %>% left_join(df_act, by = "SEQN")
} else {
    stop("activity_summary.csv not found -- this file holds mean_mims, mvpa_min, ",
         "IS, IV, ASTP, total_sleep_min, WASO, sleep_efficiency. Build it from ",
         "the PAXMIN_H per-day file before running stage30.")
}


# --- outcome families ------------------------------------------------------

primary_outcomes   <- c("mean_mims", "mvpa_min", "IS", "IV", "ASTP")
secondary_outcomes <- c("total_sleep_min", "WASO", "sleep_efficiency",
                        "sed_bout_mean", "sed_bout_p90", "act_bout_max",
                        "frac_sed_in_long_bouts")


# --- progressive adjustment formulae ---------------------------------------

m1_rhs <- "od_binary"
m2_rhs <- "od_binary + age + female + race_eth + education + pir"
m3_rhs <- paste(m2_rhs, "+ bmi + smoker_status + diabetes + comorbidity_count")
m4_rhs <- paste(m3_rhs, "+ phq9 + sinus + head_injury + nmedications")
models <- list(M1 = m1_rhs, M2 = m2_rhs, M3 = m3_rhs, M4 = m4_rhs)


# --- complete-case primary inference ---------------------------------------

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                 data = df, nest = TRUE)

fit_one <- function(outcome, rhs, design) {
    f <- as.formula(paste(outcome, "~", rhs))
    mod <- tryCatch(svyglm(f, design = design), error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    co <- summary(mod)$coefficients
    if (!"od_binary" %in% rownames(co)) return(NULL)
    r <- co["od_binary", ]
    n <- sum(!is.na(model.frame(mod)[[outcome]]))
    list(outcome = outcome, beta = r["Estimate"], se = r["Std. Error"],
         t = r["t value"], p = r["Pr(>|t|)"],
         lo95 = r["Estimate"] - 1.96 * r["Std. Error"],
         hi95 = r["Estimate"] + 1.96 * r["Std. Error"],
         n = n)
}

results <- list()
for (mod_name in names(models)) {
    rhs <- models[[mod_name]]
    for (out in c(primary_outcomes, secondary_outcomes)) {
        if (!out %in% colnames(df)) next
        r <- fit_one(out, rhs, des)
        if (!is.null(r)) {
            results[[length(results) + 1]] <- c(model = mod_name, r)
        }
    }
}
res_df <- as.data.frame(do.call(rbind, lapply(results, function(x) {
    data.frame(model = x$model, outcome = x$outcome,
               beta = x$beta, se = x$se, t = x$t, p = x$p,
               lo95 = x$lo95, hi95 = x$hi95, n = x$n,
               stringsAsFactors = FALSE)
})))

# Cohen's d standardized by survey-weighted SD (consistent with the
# survey-weighted regression coefficient in the numerator).
sds <- sapply(c(primary_outcomes, secondary_outcomes), function(v) {
    if (!v %in% colnames(df)) return(NA_real_)
    tryCatch({
        # svyvar may fail if outcome is fully NA; return unweighted SD as fallback
        sqrt(as.numeric(svyvar(as.formula(paste("~", v)), design = des, na.rm = TRUE)))
    }, error = function(e) sd(df[[v]], na.rm = TRUE))
})
res_df$cohen_d <- res_df$beta / sds[res_df$outcome]

# BH-FDR within the primary outcome family at M4
m4_primary <- res_df %>% filter(model == "M4", outcome %in% primary_outcomes)
m4_primary$q_bh <- p.adjust(m4_primary$p, method = "BH")

write_csv(res_df,     "tableS1_M4_full_coefficients.csv")
write_csv(m4_primary, "tableS2_M4_with_BHFDR.csv")

# table 2: M1..M4 for volume / intensity
t2 <- res_df %>% filter(outcome %in% c("mean_mims", "mvpa_min"))
write_csv(t2, "table2_Aim1.csv")

# table 3: M1..M4 for rhythm / sleep
t3 <- res_df %>% filter(outcome %in% c("IS", "IV", "ASTP",
                                       "total_sleep_min", "WASO", "sleep_efficiency"))
write_csv(t3, "table3_Aim2_3.csv")


# --- MICE pooled M4 (mitools::MIcombine for survey-design variance) --------
#
# pool() from mice cannot propagate complex-survey variance correctly.
# The right approach is: per imputation, build the survey design and fit
# svyglm; then combine with mitools::MIcombine, which applies Rubin's rules
# to the per-imputation survey-corrected coefficient vectors and variance
# matrices.

mice_vars <- unique(c("age", "female", "race_eth", "education", "pir",
                      "bmi", "smoker_status", "diabetes", "comorbidity_count",
                      "phq9", "sinus", "head_injury", "nmedications",
                      "od_binary",
                      intersect(c(primary_outcomes, secondary_outcomes),
                                colnames(df))))

mids <- mice(df[, mice_vars], m = 10, seed = 42, printFlag = FALSE)

mice_rows <- list()
for (out in primary_outcomes) {
    if (!out %in% colnames(df)) next
    f <- as.formula(paste(out, "~", m4_rhs))

    # fit svyglm on each imputed copy
    fits <- lapply(seq_len(mids$m), function(i) {
        di <- mice::complete(mids, i)
        di$SDMVSTRA <- df$SDMVSTRA
        di$SDMVPSU  <- df$SDMVPSU
        di$WTMEC2YR <- df$WTMEC2YR
        dsg <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = ~WTMEC2YR, data = di, nest = TRUE)
        svyglm(f, design = dsg)
    })

    combined <- MIcombine(fits)
    s <- summary(combined)
    # mitools::summary.MIcombine returns columns: results, se, (lower, upper),
    # missInfo. The 'df' column is in the print method, not the summary frame,
    # so we don't reference it here.
    if ("od_binary" %in% rownames(s)) {
        r <- s["od_binary", ]
        mice_rows[[out]] <- data.frame(
            outcome  = out,
            beta     = r[["results"]],
            se       = r[["se"]],
            lo95     = r[["(lower"]],
            hi95     = r[["upper)"]],
            missInfo = r[["missInfo"]],
            stringsAsFactors = FALSE
        )
    }
}
if (length(mice_rows)) {
    mice_df <- do.call(rbind, mice_rows)
    write_csv(mice_df, "tableS5_MICE_pooled.csv")
}

cat("Done. Wrote:\n",
    "  table2_Aim1.csv\n",
    "  table3_Aim2_3.csv\n",
    "  tableS1_M4_full_coefficients.csv\n",
    "  tableS2_M4_with_BHFDR.csv\n",
    "  tableS5_MICE_pooled.csv\n", sep = "")
