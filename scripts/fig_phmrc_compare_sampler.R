# ==============================================================
# fig_phmrc_sampler_similarity.R --- BFL methods only (domain / partial / mix),
# the three samplers (stan, gibbs_dir, gibbs_ln) side-by-side per site, to show
# the samplers give near-identical accuracy. CSV-driven from batch1_all_scores.csv
# (the only PHMRC CSV that carries all three samplers). No title, PDF output.
#   cd "<project>/Results/scripts" && Rscript fig_phmrc_sampler_similarity.R
# Needs: dplyr, ggplot2
# ==============================================================
suppressMessages({ library(dplyr); library(ggplot2) })


RESULTS_DIR <- "../CSV/"
CSV     <- file.path("../CSV/results_comparison.csv")
OUT_DIR <- file.path("../Figures", "PHMRC")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_W <- 10; FIG_H <- 5; FIG_DPI <- 300


regions     <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
fam_of      <- c(bfl_domain = "BFL-domain", bfl_partial = "BFL-partial", bfl_mix = "BFL-mix")
samp_lv     <- c("stan", "gibbs_dir", "gibbs_ln")
fill_colors <- c(stan = "#4472C4", gibbs_dir = "#E69F00", gibbs_ln = "#CC79A7")

res <- read.csv(CSV, stringsAsFactors = FALSE) %>%
  filter(shift == "no_shift", method %in% names(fam_of), sampler %in% samp_lv)

metric_cfg <- list(
  csmf     = list(col = "csmf_acc",     lab = "CSMF Accuracy"),
  top1     = list(col = "top1_acc",     lab = "Top-1 Accuracy"),
  balanced = list(col = "balanced_acc", lab = "Balanced Accuracy"))

big_theme <- function() theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.spacing = unit(0.25, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_blank())

make_fig <- function(metric) {
  cfg <- metric_cfg[[metric]]
  d <- res %>%
    transmute(region  = factor(target, levels = regions),
              family  = factor(fam_of[method], levels = unname(fam_of)),
              sampler = factor(sampler, levels = samp_lv),
              value   = .data[[cfg$col]]) %>%
    filter(!is.na(value))
  if (!nrow(d)) { message("no rows: ", metric); return(invisible()) }
  p <- ggplot(d, aes(x = family, y = value, fill = sampler)) +
    geom_boxplot(width = 0.7, colour = "black", linewidth = 0.3, outlier.shape = NA,
                 position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    facet_wrap(~region, ncol = 3, scales = "free_y") +
    labs(x = NULL, y = cfg$lab) +
    big_theme()
  out <- file.path(OUT_DIR, sprintf("phmrc_sampler_%s_no_shift.pdf", metric))
  ggsave(out, p, width = 12, height = 7, dpi = 300); cat("wrote", out, "\n")
}

for (m in c("csmf", "top1", "balanced")) make_fig(m)
cat("Done ->", OUT_DIR, "\n")
