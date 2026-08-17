################################################################################
# Figure3.R
# Trait importance manuscript
#
# Generates:
#   Figure 3            - species whose test R2 (NSE) significantly improves with
#                         each variable category
#   Figure S3           - the same comparison using nRMSE
#   Figure 3 (delta)    - paired delta-NSE effect sizes, all retained species
#   CategoryImprovementCounts - summary bar chart
#
# Replaces the "Figure 3" and "Figure 3 alt" sections of full_protocol.R
# (approximately lines 726-1010). Two changes from that version:
#
#   1. PAIRED tests. In xgboost.R the train/test split is seeded by iteration
#      (set.seed(123 + iter)) and does not depend on the variable set, so within
#      an iteration all seven variable subsets see the identical split. The
#      comparison is therefore paired and a paired test is appropriate. The
#      point estimate is unchanged (the mean of the paired differences equals
#      the difference of the means), but the standard deviation of the paired
#      differences is roughly 2.4x smaller than the unpaired equivalent, so the
#      paired test is considerably more powerful.
#
#   2. BENJAMINI-HOCHBERG FDR CORRECTION. There are 6 categories x ~86 species
#      tests. At alpha = 0.05 uncorrected, roughly 26 false positives would be
#      expected by chance alone.
#
# Note on nRMSE: the response is normalised by NFPD(), which sets each species'
# maximum to 1 and mean to 0.5, so every species is already on a common [0, 1]
# scale and RMSE_test is already range-normalised. nRMSE is therefore equal to
# RMSE_test. See NRMSE_METHOD below to change this.
#
# Written 7 August 2026
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readr)

# ---------------------------------------------------------------- settings ----

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

MODEL_PARAMS <- "data/xgb_modelParameters_20260622.rds"
OUTDIR       <- "figures"
DATESTAMP    <- format(Sys.Date(), "%Y%m%d")

ALPHA <- 0.05

# Scope of the Benjamini-Hochberg correction:
#   "all"      - one family of all species x category tests (recommended)
#   "category" - correct within each category separately (more permissive)
FDR_SCOPE <- "all"

# nRMSE definition:
#   "none"  - RMSE_test as-is. Correct here because NFPD already bounds the
#             response to [0, 1] with a maximum of exactly 1, so RMSE_test is
#             already range-normalised.
#   "range" - divide by (max - min) of the observed test values
#   "sd"    - divide by the SD of the observed test values. NOTE: this yields
#             RSR = sqrt(1 - NSE), a monotone transform of the test R2, so
#             Figure S3 would carry no information beyond Figure 3.
NRMSE_METHOD <- "none"

dir.create(OUTDIR, showWarnings = FALSE)

# ------------------------------------------------------------------- data ----

xgb_modelParameters <- readRDS(MODEL_PARAMS)

# category lookup: varSet name, printed label, outline colour, fill colour.
# Colours match the labelColor column of TableS1_Biocube_var_description.xlsx.
cats <- data.frame(
  varSet = c("notTraits", "notStr", "notClim", "notPheno", "notTerr", "notDist"),
  label  = c("Foliar Traits", "Forest Structure", "Climate",
             "Phenology", "Terrain", "Disturbance"),
  line   = c("#842B3B", "#FDC71B", "#70A4AF", "#D97C55", "#A8BBA3", "#7C4584"),
  fill   = c("#aa384c", "#FFE090", "#91b9c1", "#e29c7f", "#c4d1c0", "#AC86B0"),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------- species retention ----
# Exclude species whose full model cannot beat the mean predictor.
#
# NOTE: full_protocol.R averaged R2_test across all seven variable sets when
# applying this filter. Here the filter uses the full (spatVars) model only,
# which is the more defensible reading of the rule. Both give the same 8
# exclusions for the 20260622 run, so results are unaffected.

species_all <- unique(xgb_modelParameters$species)

keep <- xgb_modelParameters %>%
  filter(varSet == "spatVars") %>%
  group_by(species) %>%
  summarise(meanNSE = mean(R2_test), .groups = "drop") %>%
  filter(meanNSE > 0) %>%
  pull(species)

message(sprintf("Retained %d of %d species (mean full-model NSE > 0). Excluded: %s",
                length(keep), length(species_all),
                paste(sort(setdiff(species_all, keep)), collapse = ", ")))

# ------------------------------------------------------ paired differences ----

wide <- xgb_modelParameters %>%
  select(species, iteration, varSet, R2_test, RMSE_test) %>%
  pivot_wider(names_from = varSet, values_from = c(R2_test, RMSE_test))

# apply the chosen nRMSE scaling (a per-species constant, so it cancels out of
# the paired difference under "none" and "range")
rmse_scale <- switch(
  NRMSE_METHOD,
  none  = setNames(rep(1, length(species_all)), species_all),
  range = setNames(rep(1, length(species_all)), species_all),  # NFPD range = 1
  sd    = {
    # sd(y_test) backed out from RMSE = sd * sqrt(1 - NSE)
    xgb_modelParameters %>%
      filter(varSet == "spatVars", R2_test < 1) %>%
      group_by(species) %>%
      summarise(s = median(RMSE_test / sqrt(pmax(1 - R2_test, 1e-9))),
                .groups = "drop") %>%
      { setNames(.$s, .$species) }
  },
  stop("NRMSE_METHOD must be one of 'none', 'range', 'sd'")
)

# long table of paired differences, one row per species x iteration x category.
# Both differences are signed so that POSITIVE = the category improved the model.
paired <- lapply(cats$varSet, function(vs) {
  wide %>%
    transmute(
      species,
      iteration,
      varSet  = vs,
      dNSE    = .data[["R2_test_spatVars"]] - .data[[paste0("R2_test_", vs)]],
      dnRMSE  = (.data[[paste0("RMSE_test_", vs)]] -
                   .data[["RMSE_test_spatVars"]]) / rmse_scale[species],
      errRed  = 1 - .data[["RMSE_test_spatVars"]] / .data[[paste0("RMSE_test_", vs)]]
    )
}) %>%
  bind_rows() %>%
  filter(species %in% keep)

# ------------------------------------------------------------------ stats ----

# guard against a degenerate (zero-variance) vector
safe_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || length(unique(x)) < 2) return(NA_real_)
  stats::t.test(x, alternative = "greater")$p.value
}

stats_tbl <- paired %>%
  group_by(varSet, species) %>%
  summarise(
    n             = dplyr::n(),
    mean_dNSE     = mean(dNSE),
    median_dNSE   = stats::median(dNSE),
    lwr_dNSE      = as.numeric(stats::quantile(dNSE, 0.025)),
    upr_dNSE      = as.numeric(stats::quantile(dNSE, 0.975)),
    mean_dnRMSE   = mean(dnRMSE),
    median_errRed = stats::median(errRed),
    p_NSE         = safe_p(dNSE),
    p_nRMSE       = safe_p(dnRMSE),
    .groups       = "drop"
  )

if (FDR_SCOPE == "all") {
  stats_tbl <- stats_tbl %>%
    mutate(q_NSE   = p.adjust(p_NSE,   method = "BH"),
           q_nRMSE = p.adjust(p_nRMSE, method = "BH"))
} else {
  stats_tbl <- stats_tbl %>%
    group_by(varSet) %>%
    mutate(q_NSE   = p.adjust(p_NSE,   method = "BH"),
           q_nRMSE = p.adjust(p_nRMSE, method = "BH")) %>%
    ungroup()
}

stats_tbl <- stats_tbl %>%
  mutate(sig_NSE   = !is.na(q_NSE)   & q_NSE   < ALPHA,
         sig_nRMSE = !is.na(q_nRMSE) & q_nRMSE < ALPHA) %>%
  left_join(cats, by = "varSet")

write_csv(stats_tbl,
          file.path("data", paste0("categoryImprovement_pairedBH_", DATESTAMP, ".csv")))

# console summary -- these are the numbers reported in Results 3.2
summary_tbl <- stats_tbl %>%
  group_by(label) %>%
  summarise(
    n_sig_NSE      = sum(sig_NSE),
    n_sig_nRMSE    = sum(sig_nRMSE),
    n_sig_uncorr   = sum(!is.na(p_NSE) & p_NSE < ALPHA),
    median_dNSE    = round(stats::median(mean_dNSE[sig_NSE]), 4),
    max_dNSE       = round(max(mean_dNSE), 4),
    median_errRed  = round(100 * stats::median(median_errRed[sig_NSE]), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_NSE))
print(as.data.frame(summary_tbl))

# =============================================================== FIGURE 3 =====
# With / without boxplots, restricted to species passing the FDR threshold.

panel_withwithout <- function(vs, metric = c("NSE", "nRMSE")) {

  metric <- match.arg(metric)
  info   <- cats[cats$varSet == vs, ]

  sigcol <- if (metric == "NSE") "sig_NSE" else "sig_nRMSE"
  ssp    <- stats_tbl$species[stats_tbl$varSet == vs & stats_tbl[[sigcol]]]

  valcol <- if (metric == "NSE") "R2_test" else "RMSE_test"
  xlab_t <- if (metric == "NSE") "Best Model test R2" else "Best Model RMSE"

  # empty panel if nothing survives correction
  if (length(ssp) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, size = 3, colour = "grey40",
                 label = paste0("No species improved\nwith ", info$label)) +
        theme_void()
    )
  }

  dat <- xgb_modelParameters %>%
    filter(species %in% ssp, varSet %in% c("spatVars", vs)) %>%
    mutate(
      value = .data[[valcol]],
      model = ifelse(varSet == "spatVars",
                     paste("With", info$label),
                     paste("Without", info$label))
    )
  dat$model <- factor(dat$model,
                      levels = c(paste("With", info$label),
                                 paste("Without", info$label)))

  # order species by full-model performance; best models at the top
  ord <- dat %>%
    filter(varSet == "spatVars") %>%
    group_by(species) %>%
    summarise(m = mean(value), .groups = "drop")
  ord <- if (metric == "NSE") ord[order(ord$m), ] else ord[order(-ord$m), ]
  dat$species <- factor(dat$species, levels = ord$species)

  ggplot(dat, aes(x = value, y = species, colour = model, fill = model)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_colour_manual(values = setNames(c(info$line, "black"), levels(dat$model))) +
    scale_fill_manual(values   = setNames(c(info$fill, "grey30"), levels(dat$model))) +
    theme_minimal() +
    theme(legend.title    = element_blank(),
          legend.position = "bottom",
          axis.text.y     = element_text(size = 7)) +
    xlab(xlab_t) +
    ylab("")
}

# panel-height weights so rows scale with the number of species shown
row_weight <- function(vsets, sigcol) {
  max(1, max(sapply(vsets, function(v)
    sum(stats_tbl$varSet == v & stats_tbl[[sigcol]]))))
}

fig3 <- lapply(cats$varSet, panel_withwithout, metric = "NSE")
hts  <- c(row_weight(cats$varSet[1:3], "sig_NSE"),
          row_weight(cats$varSet[4:6], "sig_NSE"))

annotate_figure(
  ggarrange(plotlist = fig3, ncol = 3, nrow = 2, heights = hts, align = "hv"),
  left = text_grob("Species with significant improvement (BH FDR q < 0.05)",
                   rot = 90, size = 14, vjust = 1)
)
ggsave(file.path(OUTDIR, paste0("Figure3_NSE_", DATESTAMP, ".pdf")),
       height = 10, width = 15)

# ============================================================== FIGURE S3 =====
# Identical layout, nRMSE instead of NSE.

figS3 <- lapply(cats$varSet, panel_withwithout, metric = "nRMSE")
htsS3 <- c(row_weight(cats$varSet[1:3], "sig_nRMSE"),
           row_weight(cats$varSet[4:6], "sig_nRMSE"))

annotate_figure(
  ggarrange(plotlist = figS3, ncol = 3, nrow = 2, heights = htsS3, align = "hv"),
  left = text_grob("Species with significant improvement (BH FDR q < 0.05)",
                   rot = 90, size = 14, vjust = 1)
)
ggsave(file.path(OUTDIR, paste0("FigureS3_nRMSE_", DATESTAMP, ".pdf")),
       height = 10, width = 15)

# ========================================================= DELTA-NSE FIGURE ===
# One box per species showing the distribution of the 100 paired differences.
# ALL retained species are shown, not only the significant ones, so that the
# null distribution is visible and the figure is not selectively reported.

panel_delta <- function(vs) {

  info <- cats[cats$varSet == vs, ]
  st   <- stats_tbl %>% filter(varSet == vs)

  dat <- paired %>%
    filter(varSet == vs) %>%
    left_join(select(st, species, sig_NSE), by = "species")
  dat$sig     <- factor(ifelse(dat$sig_NSE, "TRUE", "FALSE"),
                        levels = c("FALSE", "TRUE"))
  dat$species <- factor(dat$species, levels = st$species[order(st$median_dNSE)])

  ggplot(dat, aes(x = dNSE, y = species, fill = sig, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.15, linewidth = 0.2) +
    scale_fill_manual(values   = c("FALSE" = "grey88", "TRUE" = info$fill),
                      guide = "none") +
    scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = info$line),
                        guide = "none") +
    ggtitle(info$label) +
    xlab(expression(Delta * " test R"^2 * " (full model - category removed)")) +
    ylab("") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 4),
          plot.title  = element_text(size = 11, face = "bold"))
}

figD <- lapply(cats$varSet, panel_delta)

ggarrange(plotlist = figD, ncol = 3, nrow = 2, heights = c(1, 1), align = "hv")
ggsave(file.path(OUTDIR, paste0("Figure3_deltaNSE_", DATESTAMP, ".pdf")),
       height = 16, width = 15)

# ==================================================== IMPROVEMENT COUNTS ======

cat_improvement <- stats_tbl %>%
  group_by(label, line) %>%
  summarise(count = sum(sig_NSE), .groups = "drop")

ggplot(cat_improvement, aes(x = reorder(label, -count), y = count, fill = label)) +
  geom_col() +
  scale_fill_manual(values = setNames(cat_improvement$line, cat_improvement$label)) +
  geom_text(aes(label = count), vjust = -0.5) +
  theme_minimal() +
  ylab("Species models with significant improvement in accuracy") +
  xlab("Variable category") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, paste0("CategoryImprovementCounts_", DATESTAMP, ".pdf")),
       width = 6, height = 5)

################################################################################
# end
################################################################################

