# ==============================================================
# fig_champs_paper1.R  --  imitation of Supplement figure_08/09 (CHAMPS), but
# CSV-driven from results_champs_batch8_500.csv (nchain=3 paper run). Boxes colored by
# method family, faceted by partial-label fraction (20% / 40%), dashed = BFL-base
# @frac0. One PNG per metric (csmf / balanced / top1). BFL methods use gibbs_dir.
#
#   cd "BFL_results/Figures/scripts" && Rscript fig_champs_paper1.R
# Needs: dplyr, ggplot2, tibble, tidyr, grid, patchwork
# ==============================================================
suppressMessages({ library(dplyr); library(ggplot2); library(tibble); library(tidyr); library(grid); library(patchwork) })

RESULTS_DIR <- "../CSV/"
CSV     <- file.path("../CSV/results_champs_500.csv")
OUT_DIR <- file.path("../Figures", "CHAMPS")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_W <- 10; FIG_H <- 6; FIG_DPI <- 300

SAMPLER <- Sys.getenv("BFL_FIG_SAMPLER", "gibbs_dir")

# MetBrewer "Cross" palette -> BFL = yellow, LCVA = red, GBQL = green
.cross <- tryCatch(MetBrewer::met.brewer("Cross", n = 9, override.order = TRUE),
  error = function(e) c("#c969a1","#ce4441","#ee8577","#eb7926",
                        "#ffbb44","#859b6c","#62929a","#004f63","#122451"))
fill_colors <- c(local = .cross[1], BFL = .cross[5], LCVA = .cross[2], GBQL = .cross[6])


meta <- tribble(
  ~method,       ~family,       ~fill_group,
  "local_self",  "Local self",  "local",
  "local_avg",   "Local avg",   "local",
  "bfl_domain",  "BFL-domain",  "BFL",
  "bfl_partial", "BFL-partial", "BFL",
  "bfl_mix",     "BFL-mix",     "BFL",
  "lcva_multi",  "LCVA",        "LCVA",
  "gbql_0.5",    "GBQL 0.5",    "GBQL",
  "gbql_50",     "GBQL 50",     "GBQL"
)
is_bfl <- function(m) m %in% c("bfl_base", "bfl_domain", "bfl_partial", "bfl_mix")

res <- read.csv(CSV, stringsAsFactors = FALSE)
res <- res %>% filter(!is_bfl(method) | sampler == SAMPLER)   # one sampler for BFL

frac_labels <- c("0.2" = "20% partial labels", "0.4" = "40% partial labels")
make_group_bands <- function() tibble(xmin = c(0.5, 5.5), xmax = c(2.5, 6.5))

big_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(colour = "gray85", linewidth = 0.25),
      panel.grid.minor.y = element_line(colour = "gray92", linewidth = 0.15),
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.25, "cm"),
      legend.position = "none",
      plot.title = element_blank())
}

metric_cfg <- list(
  csmf     = list(col = "csmf_acc",     methods = meta$method,      lab = "CSMF Accuracy"),
  balanced = list(col = "balanced_acc", methods = meta$method[1:6], lab = "Balanced Accuracy"),
  top1     = list(col = "top1_acc",     methods = meta$method[1:6], lab = "Top-1 Accuracy"))

make_fig <- function(metric) {
  cfg  <- metric_cfg[[metric]]
  fams <- meta$family[match(cfg$methods, meta$method)]
  d <- res %>%
    filter(frac %in% c(0.2, 0.4), method %in% cfg$methods) %>%
    left_join(meta, by = "method") %>%
    transmute(frac = factor(as.character(frac), levels = c("0.2", "0.4"), labels = frac_labels),
              family = factor(family, levels = fams),
              fill_group = factor(fill_group, levels = c("local", "BFL", "LCVA", "GBQL")),
              value = .data[[cfg$col]]) %>%
    filter(!is.na(value))
  if (!nrow(d)) { message("no rows: ", metric); return(invisible()) }

  base_y <- res %>% filter(method == "bfl_base", frac == 0) %>%
    summarise(y = mean(.data[[cfg$col]], na.rm = TRUE)) %>% pull(y)

  band_df <- make_group_bands()
  p <- ggplot(d, aes(x = family, y = value, fill = fill_group)) +
    geom_rect(data = band_df, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "gray45", alpha = 0.10) +
    { if (is.finite(base_y))
        geom_hline(yintercept = base_y, linetype = "dashed", linewidth = 0.5, colour = "black") } +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.35, outlier.shape = NA) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    facet_wrap(~frac, ncol = 2, scales = "fixed") +   # shared y-axis so base line + boxes are comparable across facets
    labs(title = paste0("CHAMPS champs batch8 500-run (nchain=3, ", SAMPLER,
                        ")  |  ", cfg$lab, "   (dashed = BFL-base @frac0)"),
         x = NULL, y = "") +
    big_theme()
  out <- file.path(OUT_DIR, sprintf("champs_%s.pdf", metric))
  ggsave(out, p, width = 12, height = 5.5, dpi = 300); cat("wrote", out, "\n")
}

for (metric in c("csmf", "balanced", "top1")) make_fig(metric)

# --------------------------------------------------------------------------
# Stacked figures (figure_09 style): CSMF row on top (all 8 methods, incl GBQL),
# second-metric row below (6 methods; GBQL has no top-1 / balanced score). The two
# rows are built as SEPARATE plots and stacked with patchwork, so the bottom row
# carries NO empty GBQL slots. Group bands (Local, LCVA) fall at the same x in both.
# --------------------------------------------------------------------------
make_stacked <- function(metric2) {
  cfg2     <- metric_cfg[[metric2]]
  fams8    <- meta$family          # top (CSMF) row: all 8 methods
  fams6    <- meta$family[1:6]     # bottom row: 6 methods (GBQL dropped, no empty slots)
  methods6 <- meta$method[1:6]

  base_d <- res %>%
    filter(frac %in% c(0.2, 0.4)) %>%
    left_join(meta, by = "method") %>%
    mutate(frac       = factor(as.character(frac), levels = c("0.2", "0.4"), labels = frac_labels),
           fill_group = factor(fill_group, levels = c("local", "BFL", "LCVA", "GBQL")))
  d_top <- base_d %>%
    mutate(family = factor(family, levels = fams8)) %>%
    transmute(frac, family, fill_group, value = csmf_acc) %>%
    filter(!is.na(value))
  d_bot <- base_d %>%
    filter(method %in% methods6) %>%
    mutate(family = factor(family, levels = fams8)) %>%
    transmute(frac, family, fill_group, value = .data[[cfg2$col]]) %>%
    filter(!is.na(value))
  if (!nrow(d_top) && !nrow(d_bot)) { message("no rows (stacked): ", metric2); return(invisible()) }

  base_line <- res %>%
    filter(method == "bfl_base", frac == 0) %>%
    summarise(CSMF = mean(csmf_acc, na.rm = TRUE),
              M2   = mean(.data[[cfg2$col]], na.rm = TRUE))
  band_df <- make_group_bands()

  row_plot <- function(d, y_base, ylab, show_x, show_strip) {
    ggplot(d, aes(x = family, y = value, fill = fill_group)) +
      geom_rect(data = band_df, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                fill = "gray45", alpha = 0.10) +
      { if (is.finite(y_base))
          geom_hline(yintercept = y_base, linetype = "dashed", linewidth = 0.45, colour = "black") } +
      geom_boxplot(width = 0.6, colour = "black", linewidth = 0.35, outlier.shape = NA) +
      scale_x_discrete(drop = FALSE) + 
      scale_fill_manual(values = fill_colors, name = NULL, drop = FALSE) +
      facet_grid(. ~ frac) +
      labs(x = NULL, y = ylab) +
      big_theme() +
      theme(axis.title.y = element_text(face = "bold", size = 12, angle = 90)) +
      { if (!show_x)     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) } +
      { if (!show_strip) theme(strip.text.x = element_blank(),
                               strip.background.x = element_blank()) }
  }

  p_top <- row_plot(d_top, base_line$CSMF, "CSMF Accuracy", show_x = FALSE, show_strip = TRUE)
  p_bot <- row_plot(d_bot, base_line$M2,   cfg2$lab,        show_x = TRUE,  show_strip = FALSE)

  p <- (p_top / p_bot) + plot_layout(axes = "collect")
    # plot_annotation(
    #   title = paste0("CHAMPS champs_paper1 (nchain=3, ", SAMPLER,
    #                  ")  |  CSMF + ", cfg2$lab, "   (dashed = BFL-base @frac0)"),
    #   theme = theme(plot.title = element_blank()))
  out <- file.path(OUT_DIR, sprintf("champs_csmf_%s_stacked.pdf", metric2))
  ggsave(out, p, width = FIG_W, height = FIG_H, dpi = FIG_DPI); cat("wrote", out, "\n")
}

for (metric2 in c("balanced", "top1")) make_stacked(metric2)
cat("Done ->", OUT_DIR, "\n")
