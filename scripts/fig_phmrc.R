# ==============================================================
# fig_phmrc_paper1.R  --  imitation of Supplement figure_02-07 (CSMF/Top-1) and
# S01-03 (Balanced), but CSV-driven from results_phmrc_batch_paper1.csv (LCVA phi).
# One PNG per (metric x shift): 2x3 facets by target, boxes colored by method
# family, gray group bands, per-site dashed bfl_base reference (no_shift only).
#
# Run on the Mac:
#   cd "BFL_results/Figures/scripts" && Rscript fig_phmrc_paper1.R
#   # override CSV / out dir via BFL_RESULTS_DIR
# Needs: dplyr, ggplot2, tibble, grid
# ==============================================================
suppressMessages({ library(dplyr); library(ggplot2); library(tibble); library(grid) })

RESULTS_DIR <- "../CSV/"
CSV     <- file.path("../CSV/results_phmrc_batch_paper1.csv")
OUT_DIR <- file.path("../Figures", "PHMRC")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_W <- 10; FIG_H <- 5; FIG_DPI <- 300

if (!file.exists(CSV)) stop("Missing CSV: ", CSV, " -> run collect_phmrc_paper1.R")

regions <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
# MetBrewer "Cross" palette -> BFL = yellow, LCVA = red, GBQL = green
.cross <- tryCatch(MetBrewer::met.brewer("Cross", n = 9, override.order = TRUE),
  error = function(e) c("#c969a1","#ce4441","#ee8577","#eb7926",
                        "#ffbb44","#859b6c","#62929a","#004f63","#122451"))
fill_colors <- c(local = .cross[1], BFL = .cross[5], LCVA = .cross[2], GBQL = .cross[6])

# method id -> (family label, fill_group)   [same families/colors as figure_02]
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

res <- read.csv(CSV, stringsAsFactors = FALSE) %>% filter(phi_source == "LCVA")

# gray95 group bands: local (x1-2) and LCVA (x6), exactly as figure_02
make_group_bands <- function() tibble(xmin = c(0.5, 5.5), xmax = c(2.5, 6.5))

big_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(colour = "gray85", linewidth = 0.25),
      panel.grid.minor.y = element_line(colour = "gray92", linewidth = 0.15),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text       = element_text(face = "bold", color = "black", size = 12),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.2, "cm"),
      legend.position = "none", plot.title = element_blank())
}

metric_cfg <- list(
  csmf     = list(col = "csmf_acc",     methods = meta$method),        # 8 (incl GBQL)
  top1     = list(col = "top1_acc",     methods = meta$method[1:6]),   # 6 (drop GBQL)
  balanced = list(col = "balanced_acc", methods = meta$method[1:6]))

make_fig <- function(metric, shift) {
  cfg  <- metric_cfg[[metric]]
  fams <- meta$family[match(cfg$methods, meta$method)]

  d <- res %>%
    filter(shift == !!shift, method %in% cfg$methods) %>%
    left_join(meta, by = "method") %>%
    transmute(region = factor(target, levels = regions),
              family = factor(family, levels = fams),
              fill_group = factor(fill_group, levels = c("local", "BFL", "LCVA", "GBQL")),
              value = .data[[cfg$col]]) %>%
    filter(!is.na(value))
  if (!nrow(d)) { message("  no rows: ", metric, " ", shift); return(invisible()) }

  base_line <- res %>%
    filter(shift == !!shift, method == "bfl_base") %>%
    transmute(region = factor(target, levels = regions), value = .data[[cfg$col]]) %>%
    filter(!is.na(value)) %>%
    group_by(region) %>% summarise(y_base = mean(value, na.rm = TRUE), .groups = "drop")

  band_df <- make_group_bands()
  p <- ggplot(d, aes(x = family, y = value, fill = fill_group)) +
    geom_rect(data = band_df, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "gray45", alpha = 0.10) +
    { if (nrow(base_line))
        geom_hline(data = base_line, aes(yintercept = y_base), inherit.aes = FALSE,
                   linetype = "dashed", linewidth = 0.45, colour = "black") } +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.35, outlier.shape = NA) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    facet_wrap(~region, ncol = 3, scales = "free_y") +
    labs(x = NULL, y = "") + big_theme()

  out <- file.path(OUT_DIR, sprintf("phmrc_%s_%s.pdf", metric, shift))
  ggsave(out, p, width = FIG_W, height = FIG_H, dpi = FIG_DPI)
  cat("wrote", out, "\n")
}

for (shift in c("no_shift", "mild", "severe"))
  for (metric in c("csmf", "top1", "balanced"))
    make_fig(metric, shift)
cat("Done ->", OUT_DIR, "\n")
