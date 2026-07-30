# ==============================================================
# fig_phmrc_paper1_verification.R  --  extra PHMRC verification figures from
# results_phmrc_batch_paper1.csv. Three sections:
#   (1) LCVA multi vs LCVA pooled {with target / no target / sources only},
#       across no_shift / mild / severe, ALL panels on a shared y-axis so the
#       three shifts are directly comparable.
#   (2) YZ variants vs normal BFL variants (base/domain/partial/mix), for each
#       shift x metric (csmf / top1 / balanced), dodged Normal vs YZ.
#   (3) OpenVA phi vs LCVA phi (no_shift), dodged, per metric.
# MetBrewer "Cross" palette throughout.
#
#   cd "BFL_results/Figures/scripts" && Rscript fig_phmrc_paper1_verification.R
# Needs: dplyr, ggplot2, tibble, tidyr, grid
# ==============================================================
suppressMessages({ library(dplyr); library(ggplot2); library(tibble); library(tidyr); library(grid) })

RESULTS_DIR <- "../CSV/"
CSV     <- file.path("../CSV/results_phmrc_batch_paper1.csv")
OUT_DIR <- file.path("../Figures", "PHMRC")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_DPI <- 300

if (!file.exists(CSV)) stop("Missing CSV: ", CSV, " -> run collect_phmrc_paper1.R")

regions      <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
shift_levels <- c("no_shift", "mild", "severe")
shift_labs   <- c(no_shift = "No shift", mild = "Mild shift", severe = "Severe shift")

# MetBrewer "Cross" palette (swatch order)
.cross <- tryCatch(MetBrewer::met.brewer("Cross", n = 9, override.order = TRUE),
  error = function(e) c("#c969a1","#ce4441","#ee8577","#eb7926",
                        "#ffbb44","#859b6c","#62929a","#004f63","#122451"))

res_all <- read.csv(CSV, stringsAsFactors = FALSE)

metric_cfg <- list(
  csmf     = list(col = "csmf_acc",     lab = "CSMF Accuracy"),
  top1     = list(col = "top1_acc",     lab = "Top-1 Accuracy"),
  balanced = list(col = "balanced_acc", lab = "Balanced Accuracy"))

big_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_line(colour = "gray85", linewidth = 0.25),
      panel.grid.minor.y = element_line(colour = "gray92", linewidth = 0.15),
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      strip.text       = element_text(face = "bold", color = "black", size = 10),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.25, "cm"),
      legend.position = "none",
      plot.title = element_blank())
}

# ==========================================================================
# Section 1: LCVA multi vs pooled (with/no target, sources only) across shifts
#   facet_grid(target ~ shift) with a SHARED y-axis (scales = "fixed"), so the
#   three shift columns are directly comparable.  (lcva_multi is sparse: n~86.)
# ==========================================================================
lcva_meta <- tribble(
  ~method,                   ~family,
  "lcva_multi",              "LCVA multi",
  "lcva_pool_with_target",   "Pool (w/ target)",
  "lcva_pool_no_target",     "Pool (no target)",
  "lcva_pool_sources_only",  "Pool (sources only)")
lcva_fill <- c("LCVA multi"          = .cross[2],   # red
               "Pool (w/ target)"    = .cross[4],   # orange
               "Pool (no target)"    = .cross[7],   # teal
               "Pool (sources only)" = .cross[1])   # magenta

lcva_section <- function(metric) {
  cfg <- metric_cfg[[metric]]
  d <- res_all %>%
    filter(phi_source == "LCVA", method %in% lcva_meta$method) %>%
    left_join(lcva_meta, by = "method") %>%
    transmute(region = factor(target, levels = regions),
              shift  = factor(shift, levels = shift_levels, labels = shift_labs[shift_levels]),
              family = factor(family, levels = lcva_meta$family),
              value  = .data[[cfg$col]]) %>%
    filter(!is.na(value))
  if (!nrow(d)) { message("no rows (lcva): ", metric); return(invisible()) }

  p <- ggplot(d, aes(x = family, y = value, fill = family)) +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.35, outlier.shape = NA) +
    scale_fill_manual(values = lcva_fill, name = NULL) +
    facet_grid(region ~ shift, scales = "fixed") +   # shared y across ALL panels
    labs(title = paste0("PHMRC: LCVA multi vs pooled variants  |  ", cfg$lab),
         x = NULL, y = "") +
    big_theme()
  out <- file.path(OUT_DIR, sprintf("verify_lcva_pool_%s.pdf", metric))
  ggsave(out, p, width = 11, height = 13, dpi = FIG_DPI); cat("wrote", out, "\n")
}

# ==========================================================================
# Section 2: YZ variants vs normal BFL variants, per shift x metric (dodged)
# ==========================================================================
yz_map <- tribble(
  ~method,          ~variant,  ~kind,
  "bfl_base",       "Base",    "Standard: Ensemble P(X|Y)",
  "bfl_domain",     "BFL-Domain",  "Standard: Ensemble P(X|Y)",
  "bfl_partial",    "BFL-Partial", "Standard: Ensemble P(X|Y)",
  "bfl_mix",        "BFL-Mix",     "Standard: Ensemble P(X|Y)",
  "bfl_yz_base",    "Base",    "Expanded: Ensemble P(X|Y,Z)",
  "bfl_yz_domain",  "BFL-Domain",  "Expanded: Ensemble P(X|Y,Z)",
  "bfl_yz_partial", "BFL-Partial", "Expanded: Ensemble P(X|Y,Z)",
  "bfl_yz_mix",     "BFL-Mix",     "Expanded: Ensemble P(X|Y,Z)")
yz_fill <- c("Standard: Ensemble P(X|Y)" = .cross[5], "Expanded: Ensemble P(X|Y,Z)" = .cross[8])   # yellow vs dark teal

yz_section <- function(metric, shift) {
  cfg <- metric_cfg[[metric]]
  # boxes: Domain/Partial/Mix only (Base is drawn as a dashed reference line)
  d <- res_all %>%
    filter(phi_source == "LCVA", shift == !!shift, method %in% yz_map$method) %>%
    left_join(yz_map, by = "method") %>%
    filter(variant != "Base") %>%
    transmute(region  = factor(target, levels = regions),
              variant = factor(variant, levels = c("BFL-Domain", "BFL-Partial", "BFL-Mix")),
              kind    = factor(kind, levels = c("Standard: Ensemble P(X|Y)", "Expanded: Ensemble P(X|Y,Z)")),
              value   = .data[[cfg$col]]) %>%
    filter(!is.na(value)) %>% droplevels()
  if (!nrow(d)) { message("no rows (yz): ", metric, " ", shift); return(invisible()) }

  # base as dashed reference line: mean bfl_base (Normal) / bfl_yz_base (YZ) per region
  base_ln <- res_all %>%
    filter(phi_source == "LCVA", shift == !!shift, method %in% c("bfl_base", "bfl_yz_base")) %>%
    left_join(yz_map, by = "method") %>%
    group_by(region = factor(target, levels = regions),
             kind = factor(kind, levels = c("Standard: Ensemble P(X|Y)", "Expanded: Ensemble P(X|Y,Z)"))) %>%
    summarise(y = mean(.data[[cfg$col]], na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(y))

  p <- ggplot(d, aes(x = variant, y = value, fill = kind)) +
    geom_boxplot(width = 0.7, position = position_dodge(preserve = "single"),
                 colour = "black", linewidth = 0.35, outlier.shape = NA) +
    { if (nrow(base_ln))
        geom_hline(data = base_ln, aes(yintercept = y, colour = kind),
                   linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) } +
    scale_fill_manual(values = yz_fill, name = NULL) +
    scale_colour_manual(values = yz_fill, guide = "none") +
    facet_wrap(~region, ncol = 3, scales = "free_y") +
    labs(title = paste0("PHMRC ", shift_labs[[shift]],
                        ": YZ vs normal BFL  |  ", cfg$lab,
                        "   (dashed = base)"), x = NULL, y = "") +
    big_theme() + theme(legend.position = "bottom")
  out <- file.path(OUT_DIR, sprintf("verify_yz_%s_%s.pdf", metric, shift))
  ggsave(out, p, width = 12, height = 8, dpi = FIG_DPI); cat("wrote", out, "\n")
}

# ==========================================================================
# Section 3: OpenVA phi vs LCVA phi (no_shift), dodged, per metric
# ==========================================================================
# Each x position pairs an LCVA-phi method (red) with its OpenVA-phi counterpart
# (green). openva_self/openva_avg are OpenVA-only baselines; their LCVA
# counterparts are local_self/local_avg (same estimator, LCVA phi).
ova_pairs <- tribble(
  ~label,        ~lcva_method,  ~ova_method,
  "Local self",  "local_self",  "openva_self",
  "Local avg",   "local_avg",   "openva_avg",
  "BFL-base",    "bfl_base",    "bfl_base",
  "BFL-domain",  "bfl_domain",  "bfl_domain",
  "BFL-partial", "bfl_partial", "bfl_partial",
  "BFL-mix",     "bfl_mix",     "bfl_mix")
ova_fill <- c(LCVA = .cross[2], InSilicoVA = .cross[7])   # red vs blue

ova_section <- function(metric) {
  cfg <- metric_cfg[[metric]]
  box_pairs <- ova_pairs %>% filter(label != "BFL-base")   # base -> dashed line, not a box
  lc <- res_all %>% filter(shift == "no_shift", phi_source == "LCVA") %>%
    inner_join(box_pairs, by = c("method" = "lcva_method")) %>%
    transmute(region = factor(target, levels = regions),
              label  = factor(label, levels = box_pairs$label),
              phi = "LCVA", value = .data[[cfg$col]])
  ov <- res_all %>% filter(shift == "no_shift", phi_source == "OpenVA") %>%
    inner_join(box_pairs, by = c("method" = "ova_method")) %>%
    transmute(region = factor(target, levels = regions),
              label  = factor(label, levels = box_pairs$label),
              phi = "InSilicoVA", value = .data[[cfg$col]])
  d <- bind_rows(lc, ov) %>%
    mutate(phi = factor(phi, levels = c("LCVA", "InSilicoVA"))) %>%
    filter(!is.na(value)) %>% droplevels()
  if (!nrow(d)) { message("no rows (openva): ", metric); return(invisible()) }

  # bfl_base as dashed reference line: LCVA (red) + OpenVA (green) mean per region
  base_ln <- res_all %>%
    filter(shift == "no_shift", method == "bfl_base", phi_source %in% c("LCVA", "OpenVA")) %>%
    group_by(region = factor(target, levels = regions),
             phi = factor(phi_source, levels = c("LCVA", "OpenVA"), labels = c("LCVA", "InSilicoVA"))) %>%
    summarise(y = mean(.data[[cfg$col]], na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(y))

  p <- ggplot(d, aes(x = label, y = value, fill = phi)) +
    geom_boxplot(width = 0.7, position = position_dodge(preserve = "single"),
                 colour = "black", linewidth = 0.35, outlier.shape = NA) +
    { if (nrow(base_ln))
        geom_hline(data = base_ln, aes(yintercept = y, colour = phi),
                   linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) } +
    scale_fill_manual(values = ova_fill, name = NULL) +
    scale_colour_manual(values = ova_fill, guide = "none") +
    facet_wrap(~region, ncol = 3, scales = "free_y") +
    labs(title = paste0("PHMRC no_shift: LCVA phi vs InSilicoVA phi  |  ", cfg$lab,
                        "   (dashed = BFL-base)"), x = NULL, y = "") +
    big_theme() + theme(legend.position = "bottom")
  out <- file.path(OUT_DIR, sprintf("verify_openva_%s.pdf", metric))
  ggsave(out, p, width = 12, height = 8, dpi = FIG_DPI); cat("wrote", out, "\n")
}

# ---- run all ----
for (metric in c("csmf", "top1", "balanced")) lcva_section(metric)                 # 3
for (metric in c("csmf", "top1", "balanced"))
  for (shift in shift_levels)                 yz_section(metric, shift)            # 9
for (metric in c("csmf", "top1", "balanced")) ova_section(metric)                 # 3
 cat("Done ->", OUT_DIR, "\n")
