# =============================================================================
#
#   RNAi Screen Analysis — C. elegans Progeny Counts
#
#   Experiment: 384-well RNAi screen, 4 biological repeats, multiple strains.
#   Baseline:   empty_vector (EV) wells define the expected relationship
#               between control and treated strains.
#
#   This script produces six figure types:
#     Fig 1 — Mean progeny scatter (linear scale)
#     Fig 2 — Mean progeny scatter (log2-transformed data)
#     Fig 3 — Mean progeny scatter (log2-scaled axes, raw tick labels)
#     Fig 4 — Control vs treated: log2(gene / EV)   [replicate-split designs]
#     Fig 5 — Control vs treated: log2(mean progeny) [replicate-split designs]
#     Fig 6 — Effect-size correlation between TIR strains (wBT182 vs wBT186)
#
# =============================================================================


# =============================================================================
# 1.  SETUP — libraries, paths, global parameters
# =============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(stats)
library(writexl)
library(ggrepel)
library(svglite)

# ---- file paths -------------------------------------------------------------

data_files <- c(
  "Y:/shared/agraf/20251118_Kinetix_4x_RNAi_1000uM/report/progeny_number_results.csv",
  "Y:/shared/agraf/20251202_Kinetix_4x_RNAi_repeat_2/report/progeny_number_results.csv",
  "Y:/shared/agraf/20251209_Kinetix_4x_RNAi_repeat_3/report/progeny_number_results.csv",
  "Y:/shared/agraf/20251216_Kinetix_4x_RNAi_repeat_4/report/progeny_number_results.csv"
)

layout_AG1     <- "Y:/shared/agraf/20251202_Kinetix_4x_RNAi_repeat_2/report/AG1.csv"
layout_control <- "Y:/shared/agraf/20251202_Kinetix_4x_RNAi_repeat_2/report/control_plate.csv"

out_dir <- "Y:/shared/agraf/Pooled_data_analysis_RNAi/log2_axis/"

# ---- analysis parameters ----------------------------------------------------

P_CUTOFF       <- 0.01    # binomial test significance threshold
FC_CUTOFF      <- 0.5     # maximum fold-change to call a hit (suppressor screen)
MIN_CTRL_PROG  <- 50      # minimum mean control progeny to consider a gene
MAX_PROGENY    <- 1000     # upper progeny filter (remove outlier wells)
HIGHLIGHT_GENE <-  ("")# gene to highlight on every plot
N_TOP_GENES    <- 0      # number of hit genes to label (split plots)

# ---- strain comparisons to run ----------------------------------------------
#   Each row: control strain  →  treated strain

comparisons <- list(
  c("wBT160", "wBT182"),
  c("wBT160", "wBT186"),
  c("wBT356", "wBT369")
)


# =============================================================================
# 2.  PUBLICATION THEME
# =============================================================================

theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      # axes
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      axis.text        = element_text(colour = "black", size = rel(0.85)),
      axis.title       = element_text(size = rel(1)),
      
      # legend
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_text(size = rel(0.85)),
      legend.text       = element_text(size = rel(0.8)),
      legend.position   = "bottom",
      
      # panel
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      
      # strip (for faceted plots)
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(size = rel(0.8), face = "bold"),
      
      # title / subtitle
      plot.title    = element_text(size = rel(1.1), face = "bold", hjust = 0),
      plot.subtitle = element_text(size = rel(0.85), hjust = 0, colour = "grey30"),
      plot.margin   = margin(8, 10, 8, 8)
    )
}

# colour constants (kept consistent throughout all figures)
COL_EV          <- "black"
COL_NONSIG      <- "grey65"
COL_SIG         <- "#2CA02C"       # green
COL_HIGHLIGHT   <- "#D62728"       # red / deeppink
COL_GENE        <- "#1F77B4"       # steel-blue (split plots)
COL_TOP         <- "#FF7F0E"       # orange (top-N labels)
COL_REFLINE     <- "grey30"


# =============================================================================
# 3.  HELPER FUNCTIONS — well mapping
# =============================================================================

#' Normalize 96-well IDs: "A1" → "A01"
normalize_96_well <- function(well) {
  row <- substr(well, 1, 1)
  col <- as.numeric(substr(well, 2, nchar(well)))
  paste0(row, sprintf("%02d", col))
}

#' Map a 384-well position to the corresponding 96-well position.
#' Adjacent 2×2 blocks in the 384 plate map to one 96-well.
map_384_to_96_well <- function(well) {
  row_letter <- substr(well, 1, 1)
  col_num    <- as.numeric(substr(gsub("^0+", "", substr(well, 2, nchar(well))), 1, 2))
  row_pos    <- match(row_letter, LETTERS)
  mapped_row <- LETTERS[(row_pos - 1) %/% 2 + 1]
  mapped_col <- (col_num - 1) %/% 2 + 1
  paste0(mapped_row, sprintf("%02d", mapped_col))
}


# =============================================================================
# 4.  DATA LOADING
# =============================================================================

message("Loading progeny data...")

df_all <- lapply(seq_along(data_files), function(i) {
  f  <- data_files[i]
  df <- read_csv(f, col_types = cols(.default = col_character())) %>%
    mutate(Progeny  = as.numeric(Progeny),
           RepeatID = i)
  # Repeat 1 contains extra strains — keep only the three used across all repeats
  if (grepl("20251118_Kinetix_4x_RNAi_1000uM", f)) {
    df <- df %>% filter(grepl("160|182|186", wBTNumber))
  }
  df
}) %>% bind_rows()

# Standardize Strain column
if (!"Strain" %in% names(df_all) && "wBTNumber" %in% names(df_all)) {
  df_all <- df_all %>% mutate(Strain = paste0("wBT", wBTNumber))
}

# Plate layouts (96-well gene map)
plate_layout <- bind_rows(
  read.csv(layout_AG1, stringsAsFactors = FALSE) %>%
    mutate(Plate = "AG1", Well96 = normalize_96_well(Well)) %>%
    dplyr::select(Plate, Well96, Gene),
  read.csv(layout_control, stringsAsFactors = FALSE) %>%
    mutate(Plate = "Control", Well96 = normalize_96_well(Well)) %>%
    dplyr::select(Plate, Well96, Gene)
)

# Map 384 → 96 and attach gene identity
df_all <- df_all %>%
  mutate(
    Plate   = as.character(Plate),
    Well384 = Well,
    Well96  = map_384_to_96_well(Well)
  ) %>%
  left_join(plate_layout, by = c("Plate", "Well96")) %>%
  filter(Gene != "-", Gene != "", Progeny <= MAX_PROGENY) %>%
  mutate(Group96 = paste0(Plate, "_", Well96))

message("  Loaded ", nrow(df_all), " wells across ",
        n_distinct(df_all$Strain), " strains and ",
        n_distinct(df_all$Gene), " gene annotations (incl. empty_vector).")


# =============================================================================
# 5.  ANALYSIS FUNCTIONS — EV model, binomial test, pairwise comparison
# =============================================================================

# ---- 5a. EV slope: mean(treated EV) / mean(control EV) ----------------------

compute_ev_slope <- function(df_ev, strainA, strainB) {
  mean_A <- df_ev %>% filter(Strain == strainA) %>% pull(Progeny) %>% mean(na.rm = TRUE)
  mean_B <- df_ev %>% filter(Strain == strainB) %>% pull(Progeny) %>% mean(na.rm = TRUE)
  
  if (!is.finite(mean_A) || mean_A <= 0) {
    warning("EV mean for ", strainA, " is non-positive or NA.")
    return(NA_real_)
  }
  mean_B / mean_A
}


# ---- 5b. Day-stratified binomial test ---------------------------------------
#   Within each repeat, every (strainA well, strainB well) pair is a trial.
#   Success = observed strainB progeny < predicted (slope × strainA progeny).
#   Under H0 (gene behaves like EV), P(success) = 0.5.

binom_test_gene <- function(gene_name, df_pair, slope, strainA, strainB) {
  df_gene <- df_pair %>% filter(Gene == gene_name)
  repeats <- unique(df_gene$RepeatID)
  
  count_below <- 0L
  trials      <- 0L
  
  for (r in repeats) {
    vals_A <- df_gene %>% filter(Strain == strainA, RepeatID == r) %>% pull(Progeny)
    vals_B <- df_gene %>% filter(Strain == strainB, RepeatID == r) %>% pull(Progeny)
    if (length(vals_A) == 0 || length(vals_B) == 0) next
    
    count_below <- count_below + sum(vapply(vals_A, function(x) {
      sum(vals_B < slope * x)
    }, integer(1)))
    trials <- trials + length(vals_A) * length(vals_B)
  }
  
  if (trials == 0) return(NA_real_)
  1 - pbinom(count_below - 1, trials, 0.5)
}


# ---- 5c. Full pairwise comparison -------------------------------------------
#   Returns a list with: $results (per-gene table), $slope, $ev_points.

compare_strains <- function(df, strainA, strainB,
                            p_cutoff  = P_CUTOFF,
                            fc_cutoff = FC_CUTOFF,
                            min_ctrl  = MIN_CTRL_PROG) {
  
  message("Comparing ", strainA, " (control) → ", strainB, " (treated)")
  
  df_pair    <- df %>% filter(Strain %in% c(strainA, strainB))
  df_ev_pair <- df_pair %>% filter(Gene == "empty_vector")
  df_genes   <- df_pair %>% filter(Gene != "empty_vector")
  
  slope <- compute_ev_slope(df_ev_pair, strainA, strainB)
  if (is.na(slope)) {
    warning("  EV slope is NA — skipping.")
    return(list(results = tibble(), slope = NA_real_, ev_points = tibble()))
  }
  message("  EV slope = ", signif(slope, 3))
  
  # EV points for scatter overlays (mean per 96-well group, per strain)
  ev_points <- inner_join(
    df_ev_pair %>% filter(Strain == strainA) %>%
      group_by(Group96) %>%
      summarize(mean_control = mean(Progeny, na.rm = TRUE), .groups = "drop"),
    df_ev_pair %>% filter(Strain == strainB) %>%
      group_by(Group96) %>%
      summarize(mean_treated = mean(Progeny, na.rm = TRUE), .groups = "drop"),
    by = "Group96"
  )
  
  # Per-gene statistics
  gene_list <- unique(df_genes$Gene)
  
  results <- bind_rows(lapply(gene_list, function(g) {
    pval   <- binom_test_gene(g, df_pair, slope, strainA, strainB)
    mean_A <- df_genes %>% filter(Gene == g, Strain == strainA) %>%
      summarize(m = mean(Progeny, na.rm = TRUE)) %>% pull(m)
    mean_B <- df_genes %>% filter(Gene == g, Strain == strainB) %>%
      summarize(m = mean(Progeny, na.rm = TRUE)) %>% pull(m)
    pred_B <- slope * mean_A
    fc     <- mean_B / pred_B
    
    tibble(
      Gene              = g,
      control_strain    = strainA,
      treated_strain    = strainB,
      mean_control      = mean_A,
      mean_treated      = mean_B,
      predicted_treated = pred_B,
      fold_change       = fc,
      p_value           = pval,
      significant       = !is.na(pval) & pval < p_cutoff &
        fc <= fc_cutoff & mean_A >= min_ctrl
    )
  }))
  
  list(results = results, slope = slope, ev_points = ev_points)
}


# =============================================================================
# 6.  PLOTTING FUNCTIONS — Figures 1–3 (pairwise strain comparisons)
# =============================================================================
get_top_hits <- function(gene_df,
                         n_top = N_TOP_GENES,
                         highlight_gene = HIGHLIGHT_GENE) {
  gene_df %>%
    filter(significant, Gene != highlight_gene) %>%
    arrange(p_value, fold_change) %>%
    slice_head(n = n_top)
}

# Shared scatter logic used by Figs 1–3
.base_scatter <- function(gene_df, ev_df,
                          highlight_gene = HIGHLIGHT_GENE,
                          n_top = N_TOP_GENES) {
  top_hits <- get_top_hits(gene_df, n_top = n_top, highlight_gene = highlight_gene)
  
  list(
    # EV wells
    geom_point(data = ev_df,
               aes(x = x, y = y), colour = COL_EV, size = 4, alpha = 0.6,
               shape = 16),
    
    # gene wells (coloured by significance)
    geom_point(data = gene_df,
               aes(x = x, y = y, colour = significant),
               size = 4, alpha = 0.6, shape = 16),
    
    # top hit genes in orange
    geom_point(data = top_hits, aes(x = x, y = y),
               colour = COL_TOP, size = 4, shape = 16),
    geom_text_repel(data = top_hits, aes(x = x, y = y, label = Gene),
                    colour = COL_TOP, size = 3, fontface = "bold.italic",
                    max.overlaps = Inf, seed = 42),
    
    # significance colour scale
    scale_colour_manual(
      name   = NULL,
      values = c(`FALSE` = COL_NONSIG, `TRUE` = COL_SIG),
      labels = c(`FALSE` = "n.s.", `TRUE` = "significant")
    )
  )
}

# ---- Fig 1: linear scale ----------------------------------------------------

plot_scatter_linear <- function(res, strainA, strainB) {
  df  <- res$results; ev <- res$ev_points; s <- res$slope
  if (nrow(df) == 0) return(invisible(NULL))
  
  gene_df <- df %>% mutate(x = mean_control, y = mean_treated)
  ev_df   <- ev %>% mutate(x = mean_control, y = mean_treated)
  n_sig   <- sum(df$significant)
  
  p <- ggplot() +
    .base_scatter(gene_df, ev_df) +
    geom_abline(intercept = 0, slope = s,
                linetype = "dashed", colour = COL_REFLINE, linewidth = 0.5) +
    labs(
      x        = bquote("Mean progeny (" * .(strainA) * ")"),
      y        = bquote("Mean progeny (" * .(strainB) * ")"),
      subtitle = paste0(nrow(df), " genes | ", n_sig, " hits (p < ",
                        P_CUTOFF, ", FC \u2264 ", FC_CUTOFF, ")")
    ) +
    theme_pub()
  
  fname <- file.path(out_dir, paste0("fig1_scatter_linear_",
                                     strainA, "_vs_", strainB, ".svg"))
  ggsave(fname, p, width = 11, height = 10, device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# ---- Fig 2: log2-transformed data -------------------------------------------

plot_scatter_log2_data <- function(res, strainA, strainB) {
  df  <- res$results; ev <- res$ev_points; s <- res$slope
  if (nrow(df) == 0) return(invisible(NULL))
  
  gene_df <- df %>%
    filter(mean_control > 0, mean_treated > 0) %>%
    mutate(x = log2(mean_control), y = log2(mean_treated))
  ev_df <- ev %>%
    filter(mean_control > 0, mean_treated > 0) %>%
    mutate(x = log2(mean_control), y = log2(mean_treated))
  n_sig <- sum(gene_df$significant)
  
  p <- ggplot() +
    .base_scatter(gene_df, ev_df) +
    # In log2 space: log2(y) = log2(slope) + log2(x)  →  slope=1, intercept=log2(s)
    geom_abline(intercept = log2(s), slope = 1,
                linetype = "dashed", colour = COL_REFLINE, linewidth = 0.5) +
    labs(
      title = paste("RNAi gene effects vs EV model:", strainA, "→", strainB),
      x        = bquote(log[2] ~ "mean progeny (" * .(strainA) * ")"),
      y        = bquote(log[2] ~ "mean progeny (" * .(strainB) * ")"),
      subtitle = paste0(nrow(df), " genes | ", n_sig, " hits (p < ",
                        P_CUTOFF, ", FC \u2264 ", FC_CUTOFF, ")")
    ) +
    theme_pub()
  
  fname <- file.path(out_dir, paste0("fig2_scatter_log2data_",
                                     strainA, "_vs_", strainB, ".svg"))
  ggsave(fname, p, width = 11, height = 10, device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# ---- Fig 3: log2-scaled axes (tick labels = raw counts) ---------------------

plot_scatter_log2_axis <- function(res, strainA, strainB) {
  df  <- res$results; ev <- res$ev_points; s <- res$slope
  if (nrow(df) == 0) return(invisible(NULL))
  
  gene_df <- df %>% filter(mean_control > 0, mean_treated > 0) %>%
    mutate(x = mean_control, y = mean_treated)
  ev_df <- ev %>% filter(mean_control > 0, mean_treated > 0) %>%
    mutate(x = mean_control, y = mean_treated)
  n_sig <- sum(gene_df$significant)
  
  # EV model line (y = slope * x) as an explicit set of points
  x_range <- range(c(gene_df$x, ev_df$x), na.rm = TRUE)
  ev_line <- data.frame(
    x = seq(max(x_range[1], 0.5), x_range[2], length.out = 300)
  ) %>% mutate(y = s * x)
  
  brks <- c(1, 2, 5, 10, 25, 50, 100, 200, 400, 800)
  
  p <- ggplot() +
    geom_line(data = ev_line, aes(x = x, y = y),
              linetype = "dashed", colour = COL_REFLINE, linewidth = 0.5) +
    .base_scatter(gene_df, ev_df) +
    scale_x_continuous(trans = "log2", breaks = brks, labels = brks) +
    scale_y_continuous(trans = "log2", breaks = brks, labels = brks) +
    labs(
      x        = bquote("Mean progeny " * .(strainA) * " (log2 scale)"),
      y        = bquote("Mean progeny " * .(strainB) * " (log2 scale)"),
      subtitle = paste0(nrow(gene_df), " genes | ", n_sig, " hits | log2-scaled axes")
    ) +
    theme_pub()
  
  fname <- file.path(out_dir, paste0("fig3_scatter_log2axis_",
                                     strainA, "_vs_", strainB, ".svg"))
  ggsave(fname, p, width = 11, height = 10, device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# =============================================================================
# 7.  ANALYSIS FUNCTIONS — replicate-split designs (Figs 4–6)
# =============================================================================

# ---- 7a. log2(gene mean / EV mean) within a strain × rep subset -------------

compute_log2_gene_over_ev <- function(df, strain, reps) {
  sub <- df %>% filter(Strain == strain, RepeatID %in% reps)
  ev_mean <- sub %>% filter(Gene == "empty_vector") %>%
    pull(Progeny) %>% mean(na.rm = TRUE)
  if (!is.finite(ev_mean) || ev_mean <= 0) return(tibble())
  
  genes <- df %>% filter(!Gene %in% c("empty_vector", "-", "")) %>%
    pull(Gene) %>% unique()
  
  bind_rows(lapply(genes, function(g) {
    vals <- sub %>% filter(Gene == g) %>% pull(Progeny)
    if (length(vals) == 0) return(NULL)
    gm <- mean(vals, na.rm = TRUE)
    if (!is.finite(gm) || gm <= 0) return(NULL)
    tibble(Gene = g, value = log2(gm / ev_mean))
  })) %>% filter(is.finite(value))
}


# ---- 7b. log2(mean progeny) per gene ----------------------------------------

compute_log2_mean_progeny <- function(df, strain, reps) {
  sub <- df %>% filter(Strain == strain, RepeatID %in% reps)
  
  genes <- df %>% filter(!Gene %in% c("empty_vector", "-", "")) %>%
    pull(Gene) %>% unique()
  
  bind_rows(lapply(genes, function(g) {
    vals <- sub %>% filter(Gene == g) %>% pull(Progeny)
    if (length(vals) == 0) return(NULL)
    gm <- mean(vals, na.rm = TRUE)
    if (!is.finite(gm) || gm <= 0) return(NULL)
    tibble(Gene = g, value = log2(gm))
  })) %>% filter(is.finite(value))
}


# ---- 7c. EV well values per Group96 (for plotting as sanity check) ----------

compute_ev_values <- function(df, strain, reps, metric = c("ratio", "raw")) {
  metric <- match.arg(metric)
  sub <- df %>% filter(Strain == strain, RepeatID %in% reps, Gene == "empty_vector")
  ev_mean <- mean(sub$Progeny, na.rm = TRUE)
  if (!is.finite(ev_mean) || ev_mean <= 0) return(tibble())
  
  sub %>%
    group_by(Group96) %>%
    summarize(grp_mean = mean(Progeny, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(grp_mean), grp_mean > 0) %>%
    mutate(value = if (metric == "ratio") log2(grp_mean / ev_mean) else log2(grp_mean)) %>%
    dplyr::select(Group96, value)
}


# ---- 7d. log2(fold-change vs EV) per gene -----------------------------------
#   FC = mean_treated / (slope × mean_control)

compute_log2fc <- function(df, ctrl_strain, trt_strain, ctrl_reps, trt_reps) {
  dc <- df %>% filter(Strain == ctrl_strain, RepeatID %in% ctrl_reps)
  dt <- df %>% filter(Strain == trt_strain,  RepeatID %in% trt_reps)
  
  ev_ctrl <- dc %>% filter(Gene == "empty_vector") %>% pull(Progeny) %>% mean(na.rm = TRUE)
  ev_trt  <- dt %>% filter(Gene == "empty_vector") %>% pull(Progeny) %>% mean(na.rm = TRUE)
  
  if (!is.finite(ev_ctrl) || ev_ctrl <= 0 || !is.finite(ev_trt) || ev_trt <= 0) {
    return(tibble(Gene = character(), log2_fc = numeric()))
  }
  
  slope <- ev_trt / ev_ctrl
  
  genes <- df %>%
    filter(!Gene %in% c("empty_vector", "-", "")) %>%
    pull(Gene) %>%
    unique()
  
  out <- bind_rows(lapply(genes, function(g) {
    cv <- dc %>% filter(Gene == g) %>% pull(Progeny)
    tv <- dt %>% filter(Gene == g) %>% pull(Progeny)
    if (length(cv) == 0 || length(tv) == 0) return(NULL)
    
    mc <- mean(cv, na.rm = TRUE)
    mt <- mean(tv, na.rm = TRUE)
    pred <- slope * mc
    
    if (!is.finite(pred) || pred <= 0 || !is.finite(mt) || mt <= 0) return(NULL)
    
    tibble(Gene = g, log2_fc = log2(mt / pred))
  }))
  
  if (nrow(out) == 0) {
    return(tibble(Gene = character(), log2_fc = numeric()))
  }
  
  out %>% filter(is.finite(log2_fc))
}


# ---- 7e. EV fold-change per Group96 (sanity, should cluster near 0) ---------

compute_ev_log2fc <- function(df, ctrl_strain, trt_strain, ctrl_reps, trt_reps) {
  dc <- df %>% filter(Strain == ctrl_strain, RepeatID %in% ctrl_reps, Gene == "empty_vector")
  dt <- df %>% filter(Strain == trt_strain,  RepeatID %in% trt_reps,  Gene == "empty_vector")
  
  mc_all <- mean(dc$Progeny, na.rm = TRUE)
  mt_all <- mean(dt$Progeny, na.rm = TRUE)
  
  if (!is.finite(mc_all) || mc_all <= 0 || !is.finite(mt_all) || mt_all <= 0) {
    return(tibble(Group96 = character(), log2_fc = numeric()))
  }
  
  slope <- mt_all / mc_all
  
  out <- inner_join(
    dc %>% group_by(Group96) %>% summarize(mc = mean(Progeny, na.rm = TRUE), .groups = "drop"),
    dt %>% group_by(Group96) %>% summarize(mt = mean(Progeny, na.rm = TRUE), .groups = "drop"),
    by = "Group96"
  ) %>%
    mutate(log2_fc = log2(mt / (slope * mc))) %>%
    filter(is.finite(log2_fc)) %>%
    dplyr::select(Group96, log2_fc)
  
  if (nrow(out) == 0) {
    return(tibble(Group96 = character(), log2_fc = numeric()))
  }
  
  out
}


# =============================================================================
# 8.  REPLICATE-SPLIT DATA BUILDERS
# =============================================================================

# The 4 repeats are split between wBT182 and wBT186 control pools to ensure
# statistical independence.  Three designs:
#   Version A — 2 control reps + 2 TIR reps  (exclusive, most conservative)
#   Version B — 2 control reps + 4 TIR reps  (maximizes TIR power)
#   Version C — 4 control reps + 4 TIR reps  (shared control, maximal power)

all_reps     <- 1:4
rep_splits   <- combn(all_reps, 2, simplify = FALSE)   # 6 ways to choose 2 of 4
tir_halves   <- list(c(1, 2), c(3, 4))                 # TIR 2-rep splits
tir_all      <- 1:4


# ---- Builder for Fig 4 (EV-normalised) and Fig 5 (raw log2 mean) -----------

build_ctrl_vs_treated <- function(df, treated_strain, ctrl_reps, trt_reps,
                                  label, metric = c("ratio", "raw")) {
  metric <- match.arg(metric)
  compute_fn <- if (metric == "ratio") compute_log2_gene_over_ev else compute_log2_mean_progeny
  
  ctrl_genes <- compute_fn(df, "wBT160", ctrl_reps)
  trt_genes  <- compute_fn(df, treated_strain, trt_reps)
  
  ctrl_ev_mean <- df %>%
    filter(Strain == "wBT160", RepeatID %in% ctrl_reps, Gene == "empty_vector") %>%
    pull(Progeny) %>%
    mean(na.rm = TRUE)
  
  trt_ev_mean <- df %>%
    filter(Strain == treated_strain, RepeatID %in% trt_reps, Gene == "empty_vector") %>%
    pull(Progeny) %>%
    mean(na.rm = TRUE)
  
  ev_slope <- trt_ev_mean / ctrl_ev_mean
  
  genes <- inner_join(
    ctrl_genes %>% dplyr::select(Gene, log2_ctrl = value),
    trt_genes  %>% dplyr::select(Gene, log2_trt  = value),
    by = "Gene"
  ) %>%
    mutate(
      split = label,
      type = "gene",
      ev_slope = ev_slope,
      ev_ctrl_center = if (metric == "raw") log2(ctrl_ev_mean) else 0,
      ev_trt_center  = if (metric == "raw") log2(trt_ev_mean)  else 0
    )
  
  ctrl_ev <- compute_ev_values(df, "wBT160", ctrl_reps, metric)
  trt_ev  <- compute_ev_values(df, treated_strain, trt_reps, metric)
  
  ev <- inner_join(
    ctrl_ev %>% dplyr::select(Group96, log2_ctrl = value),
    trt_ev  %>% dplyr::select(Group96, log2_trt  = value),
    by = "Group96"
  ) %>%
    mutate(
      split = label,
      type = "empty_vector",
      Gene = "EV",
      ev_slope = ev_slope,
      ev_ctrl_center = if (metric == "raw") log2(ctrl_ev_mean) else 0,
      ev_trt_center  = if (metric == "raw") log2(trt_ev_mean)  else 0
    )
  
  bind_rows(genes, ev)
}


# ---- Builder for Fig 6 (effect-size correlation 182 vs 186) -----------------

build_effect_corr <- function(df, reps_182, reps_186, tir_reps, label) {
  fc_182 <- compute_log2fc(df, "wBT160", "wBT182", reps_182, tir_reps)
  fc_186 <- compute_log2fc(df, "wBT160", "wBT186", reps_186, tir_reps)
  ev_182 <- compute_ev_log2fc(df, "wBT160", "wBT182", reps_182, tir_reps)
  ev_186 <- compute_ev_log2fc(df, "wBT160", "wBT186", reps_186, tir_reps)
  
  genes <- inner_join(
    fc_182 %>% dplyr::select(Gene, log2_fc_182 = log2_fc),
    fc_186 %>% dplyr::select(Gene, log2_fc_186 = log2_fc),
    by = "Gene"
  ) %>% mutate(split = label, type = "gene")
  
  ev <- inner_join(
    ev_182 %>% dplyr::select(Group96, log2_fc_182 = log2_fc),
    ev_186 %>% dplyr::select(Group96, log2_fc_186 = log2_fc),
    by = "Group96"
  ) %>% mutate(split = label, type = "empty_vector", Gene = "EV")
  
  bind_rows(genes, ev)
}


# =============================================================================
# 9.  PLOTTING FUNCTIONS — Figures 4–6 (split-based analyses)
# =============================================================================

# ---- Shared helper: top-N outlier genes + correlation annotation ------------

.annotate_split_plot <- function(p, plot_data, xcol, ycol, n_top = N_TOP_GENES) {
  
  
  top_genes <- plot_data %>%
    filter(type == "gene", Gene != HIGHLIGHT_GENE) %>%
    group_by(Gene) %>%
    summarize(
      mean_dist = mean(sqrt(.data[[xcol]]^2 + .data[[ycol]]^2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    slice_max(mean_dist, n = n_top) %>% pull(Gene)
  
  top_df <- plot_data %>% filter(Gene %in% top_genes)
  
  cor_labs <- plot_data %>%
    filter(type == "gene") %>%
    group_by(split) %>%
    summarize(
      r = cor(.data[[xcol]], .data[[ycol]], use = "complete.obs"),
      n = n(), .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", round(r, 3), "  (n = ", n, ")"))
  
  p +
    # top-N outliers
    geom_point(data = top_df, aes(x = .data[[xcol]], y = .data[[ycol]]),
               colour = COL_TOP, size = 2.5, shape = 16) +
    geom_text_repel(data = top_df,
                    aes(x = .data[[xcol]], y = .data[[ycol]], label = Gene),
                    colour = COL_TOP, size = 2.3, fontface = "bold.italic",
                    max.overlaps = Inf, seed = 42) +
   
    # correlation annotation
    geom_text(data = cor_labs,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.3, size = 2.8, colour = "grey20")
}


# ---- Fig 4: Control vs Treated — log2(gene / EV) ----------------------------

plot_ctrl_vs_treated_ev <- function(plot_data, treated_label, version_label,
                                    faceted = TRUE) {
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = COL_REFLINE, linewidth = 0.4) +
    geom_point(data = plot_data %>% filter(type == "empty_vector"),
               aes(x = log2_ctrl, y = log2_trt),
               colour = COL_EV, size = 1.6, alpha = 0.5, shape = 16) +
    geom_point(data = plot_data %>% filter(type == "gene"),
               aes(x = log2_ctrl, y = log2_trt),
               colour = COL_GENE, size = 1.6, alpha = 0.6, shape = 16) +
    labs(
      x = bquote(log[2] * "(gene / EV) in wBT160"),
      y = bquote(log[2] * "(gene / EV) in " * .(treated_label))
    ) +
    theme_pub()
  
  p <- .annotate_split_plot(p, plot_data, "log2_ctrl", "log2_trt")
  
  if (faceted) p <- p + facet_wrap(~ split, ncol = 3)
  dims <- .facet_dims(plot_data, faceted)
  
  fname <- file.path(out_dir, paste0("fig4_ev_normalised_",
                                     gsub(" ", "_", treated_label), "_",
                                     gsub("[^a-zA-Z0-9]", "_", version_label), ".svg"))
  ggsave(fname, p, width = dims[1], height = dims[2], device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# ---- Fig 5: Control vs Treated — log2(mean progeny) -------------------------

plot_ctrl_vs_treated_raw <- function(plot_data, treated_label, version_label,
                                     faceted = TRUE) {
  
  ref_lines <- plot_data %>%
    group_by(split) %>%
    summarize(
      ev_slope       = first(ev_slope),
      ev_ctrl_center = first(ev_ctrl_center),
      ev_trt_center  = first(ev_trt_center),
      .groups = "drop"
    ) %>%
    mutate(log2_intercept = log2(ev_slope))
  
  p <- ggplot() +
    geom_vline(data = ref_lines,
               aes(xintercept = ev_ctrl_center),
               linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_hline(data = ref_lines,
               aes(yintercept = ev_trt_center),
               linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_abline(data = ref_lines,
                aes(intercept = log2_intercept, slope = 1),
                linetype = "dashed", colour = COL_REFLINE, linewidth = 0.4) +
    geom_point(data = plot_data %>% filter(type == "empty_vector"),
               aes(x = log2_ctrl, y = log2_trt),
               colour = COL_EV, size = 1.6, alpha = 0.5, shape = 16) +
    geom_point(data = plot_data %>% filter(type == "gene"),
               aes(x = log2_ctrl, y = log2_trt),
               colour = COL_GENE, size = 1.6, alpha = 0.6, shape = 16) +
    labs(
      x = bquote(log[2] * "(mean progeny) in wBT160"),
      y = bquote(log[2] * "(mean progeny) in " * .(treated_label))
    ) +
    theme_pub()
  
  p <- .annotate_split_plot(p, plot_data, "log2_ctrl", "log2_trt")
  
  if (faceted) p <- p + facet_wrap(~ split, ncol = 3)
  dims <- .facet_dims(plot_data, faceted)
  
  fname <- file.path(out_dir, paste0("fig5_raw_progeny_",
                                     gsub(" ", "_", treated_label), "_",
                                     gsub("[^a-zA-Z0-9]", "_", version_label), ".svg"))
  ggsave(fname, p, width = dims[1], height = dims[2], device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# ---- Fig 6: Effect-size correlation between TIR strains --------------------

plot_effect_correlation <- function(plot_data, version_label, faceted = TRUE) {
  
  # For top-gene selection, use distance from origin (large effect in either TIR)
  top_genes <- plot_data %>%
    filter(type == "gene", Gene != HIGHLIGHT_GENE) %>%
    group_by(Gene) %>%
    summarize(
      mean_dist = mean(sqrt(log2_fc_182^2 + log2_fc_186^2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    slice_max(mean_dist, n = N_TOP_GENES) %>% pull(Gene)
  top_df <- plot_data %>% filter(Gene %in% top_genes)
  
  cor_labs <- plot_data %>%
    filter(type == "gene") %>%
    group_by(split) %>%
    summarize(
      r = cor(log2_fc_182, log2_fc_186, use = "complete.obs"),
      n = n(), .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", round(r, 3), "  (n = ", n, ")"))
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = COL_REFLINE, linewidth = 0.4) +
    geom_point(data = plot_data %>% filter(type == "empty_vector"),
               aes(x = log2_fc_182, y = log2_fc_186),
               colour = COL_EV, size = 1.6, alpha = 0.5, shape = 16) +
    geom_point(data = plot_data %>% filter(type == "gene"),
               aes(x = log2_fc_182, y = log2_fc_186),
               colour = COL_GENE, size = 1.6, alpha = 0.6, shape = 16) +
    geom_point(data = top_df, aes(x = log2_fc_182, y = log2_fc_186),
               colour = COL_TOP, size = 2.5, shape = 16) +
    geom_text_repel(data = top_df,
                    aes(x = log2_fc_182, y = log2_fc_186, label = Gene),
                    colour = COL_TOP, size = 2.3, fontface = "bold.italic",
                    max.overlaps = Inf, seed = 42) +
    geom_text(data = cor_labs,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.3, size = 2.8, colour = "grey20") +
    labs(
      x = bquote(log[2] * "(fold change) — wBT182"),
      y = bquote(log[2] * "(fold change) — wBT186")
    ) +
    theme_pub()
  
  if (faceted) p <- p + facet_wrap(~ split, ncol = 3)
  dims <- .facet_dims(plot_data, faceted)
  
  fname <- file.path(out_dir, paste0("fig6_effect_correlation_",
                                     gsub("[^a-zA-Z0-9]", "_", version_label), ".svg"))
  ggsave(fname, p, width = dims[1], height = dims[2], device = "svg")
  message("  Saved: ", fname)
  invisible(p)
}


# ---- Helper: figure dimensions for faceted/single panels --------------------

.facet_dims <- function(plot_data, faceted) {
  if (faceted) {
    n_panels <- length(unique(plot_data$split))
    n_rows   <- ceiling(n_panels / 3)
    c(14, 4 * n_rows + 1)
  } else {
    c(5.5, 5)
  }
}


# =============================================================================
# 10.  RUN — Part A: pairwise strain comparisons (Figs 1–3)
# =============================================================================

message("\n", strrep("=", 60))
message("PART A: Pairwise strain comparisons")
message(strrep("=", 60))

all_results <- list()

for (comp in comparisons) {
  sA <- comp[1]; sB <- comp[2]
  res <- compare_strains(df_all, sA, sB)
  
  plot_scatter_linear(res, sA, sB)
  plot_scatter_log2_data(res, sA, sB)
  plot_scatter_log2_axis(res, sA, sB)
  
  all_results[[paste(sA, sB, sep = "_vs_")]] <- res$results
}

# Export combined table
results_table <- bind_rows(all_results)
write_xlsx(results_table,
           file.path(out_dir, "RNAi_screen_all_pairwise_results.xlsx"))
message("  Exported results table: RNAi_screen_all_pairwise_results.xlsx")


# =============================================================================
# 11.  RUN — Part B: replicate-split analyses (Figs 4–6)
#
#   These analyses apply only to the wBT160 → wBT182 / wBT186 comparisons.
#   Repeats are split between the two TIR strains to ensure independence.
# =============================================================================

message("\n", strrep("=", 60))
message("PART B: Replicate-split analyses (Versions A / B / C)")
message(strrep("=", 60))

# Storage lists
fig4_A_182 <- list(); fig4_A_186 <- list()
fig5_A_182 <- list(); fig5_A_186 <- list()
fig6_A      <- list()

fig4_B_182 <- list(); fig4_B_186 <- list()
fig5_B_182 <- list(); fig5_B_186 <- list()
fig6_B      <- list()

counter_A <- 1

for (i in seq_along(rep_splits)) {
  reps_182 <- rep_splits[[i]]
  reps_186 <- setdiff(all_reps, reps_182)
  base_lbl <- paste0("160[", paste(reps_182, collapse = ","),
                     "]\u2192182 / 160[", paste(reps_186, collapse = ","), "]\u2192186")
  
  # --- Version A: 2 ctrl reps × 2 TIR reps (exclusive) ---
  for (tir2 in tir_halves) {
    lbl_A <- paste0(base_lbl, " | TIR[", paste(tir2, collapse = ","), "]")
    
    fig4_A_182[[counter_A]] <- build_ctrl_vs_treated(df_all, "wBT182", reps_182, tir2, lbl_A, "ratio")
    fig4_A_186[[counter_A]] <- build_ctrl_vs_treated(df_all, "wBT186", reps_186, tir2, lbl_A, "ratio")
    fig5_A_182[[counter_A]] <- build_ctrl_vs_treated(df_all, "wBT182", reps_182, tir2, lbl_A, "raw")
    fig5_A_186[[counter_A]] <- build_ctrl_vs_treated(df_all, "wBT186", reps_186, tir2, lbl_A, "raw")
    fig6_A[[counter_A]]     <- build_effect_corr(df_all, reps_182, reps_186, tir2, lbl_A)
    
    counter_A <- counter_A + 1
  }
  
  # --- Version B: 2 ctrl reps × 4 TIR reps ---
  lbl_B <- paste0(base_lbl, " | TIR[1,2,3,4]")
  
  fig4_B_182[[i]] <- build_ctrl_vs_treated(df_all, "wBT182", reps_182, tir_all, lbl_B, "ratio")
  fig4_B_186[[i]] <- build_ctrl_vs_treated(df_all, "wBT186", reps_186, tir_all, lbl_B, "ratio")
  fig5_B_182[[i]] <- build_ctrl_vs_treated(df_all, "wBT182", reps_182, tir_all, lbl_B, "raw")
  fig5_B_186[[i]] <- build_ctrl_vs_treated(df_all, "wBT186", reps_186, tir_all, lbl_B, "raw")
  fig6_B[[i]]     <- build_effect_corr(df_all, reps_182, reps_186, tir_all, lbl_B)
}

# --- Version C: 4 ctrl × 4 TIR (shared control) ---
lbl_C <- "4\u00d74 shared control"

fig4_C_182 <- build_ctrl_vs_treated(df_all, "wBT182", all_reps, tir_all, lbl_C, "ratio")
fig4_C_186 <- build_ctrl_vs_treated(df_all, "wBT186", all_reps, tir_all, lbl_C, "ratio")
fig5_C_182 <- build_ctrl_vs_treated(df_all, "wBT182", all_reps, tir_all, lbl_C, "raw")
fig5_C_186 <- build_ctrl_vs_treated(df_all, "wBT186", all_reps, tir_all, lbl_C, "raw")
fig6_C <- build_effect_corr(df_all, all_reps, all_reps, tir_all, lbl_C)


# =============================================================================
# 12.  GENERATE SPLIT PLOTS
# =============================================================================

message("\nGenerating Fig 4 (EV-normalised ctrl vs treated)...")
plot_ctrl_vs_treated_ev(bind_rows(fig4_A_182), "wBT182", "Version_A_2x2_exclusive")
plot_ctrl_vs_treated_ev(bind_rows(fig4_A_186), "wBT186", "Version_A_2x2_exclusive")
plot_ctrl_vs_treated_ev(bind_rows(fig4_B_182), "wBT182", "Version_B_2x4")
plot_ctrl_vs_treated_ev(bind_rows(fig4_B_186), "wBT186", "Version_B_2x4")
plot_ctrl_vs_treated_ev(fig4_C_182, "wBT182", "Version_C_4x4_shared", faceted = FALSE)
plot_ctrl_vs_treated_ev(fig4_C_186, "wBT186", "Version_C_4x4_shared", faceted = FALSE)

message("\nGenerating Fig 5 (raw log2 mean progeny, ctrl vs treated)...")
plot_ctrl_vs_treated_raw(bind_rows(fig5_A_182), "wBT182", "Version_A_2x2_exclusive")
plot_ctrl_vs_treated_raw(bind_rows(fig5_A_186), "wBT186", "Version_A_2x2_exclusive")
plot_ctrl_vs_treated_raw(bind_rows(fig5_B_182), "wBT182", "Version_B_2x4")
plot_ctrl_vs_treated_raw(bind_rows(fig5_B_186), "wBT186", "Version_B_2x4")
plot_ctrl_vs_treated_raw(fig5_C_182, "wBT182", "Version_C_4x4_shared", faceted = FALSE)
plot_ctrl_vs_treated_raw(fig5_C_186, "wBT186", "Version_C_4x4_shared", faceted = FALSE)

message("\nGenerating Fig 6 (effect-size correlation, wBT182 vs wBT186)...")
plot_effect_correlation(bind_rows(fig6_A), "Version_A_2x2_exclusive")
plot_effect_correlation(bind_rows(fig6_B), "Version_B_2x4")
plot_effect_correlation(fig6_C, "Version_C_4x4_shared", faceted = FALSE)


# =============================================================================
# 13.  DONE
# =============================================================================

message("\n", strrep("=", 60))
message("All figures and tables saved to:\n  ", out_dir)
message(strrep("=", 60))

