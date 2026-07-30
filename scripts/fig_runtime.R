# ==============================================================
# fig_runtime.R --- runtime / computational-cost figures for the supplement,
# styled to match the paper (theme_minimal + black panel border + lightgray
# strips + MetBrewer "Renoir" palette).
#
# Reads the two timing CSVs and writes two PDFs into
#   Results/Figures/Supplemnt Figures/Runtime/ :
#     runtime_vs_iterations.pdf  (Gibbs samplers vs iterations, Stan reference)
#     runtime_breakdown.pdf      (per-step end-to-end breakdown)
#
# Naming: "Gibbs-Dir" = Gibbs sampler with the Dirichlet prior on the model
# weights; "Gibbs-LN" = Gibbs sampler with the logistic-normal prior.
# "(4k)/(20k)" = number of MCMC iterations (4,000 / 20,000). A single
# "Stan/NUTS" bar is shown for the Hamiltonian-Monte-Carlo aggregation.
#
#   cd "Results/Figures/scripts" && Rscript fig_runtime.R
# Needs: dplyr, ggplot2  (MetBrewer optional -- falls back to hard-coded Renoir)
# ==============================================================
suppressMessages({ library(dplyr); library(ggplot2) })

RESULTS_DIR <- Sys.getenv("BFL_RESULTS_DIR",
  "/Users/toastymac/Desktop/BFL/BFL Reports/BFL Project/Results")
DIR <- file.path(RESULTS_DIR, "Figures", "Supplemnt Figures", "Runtime")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

# ---- paper palette + theme -------------------------------------------------
renoir <- tryCatch(as.character(MetBrewer::met.brewer("Renoir", n = 12, override.order = FALSE)),
  error = function(e) c("#17154f","#2f357c","#6c5d9e","#9d9cd5","#b0799a","#f6b3b0",
                        "#e48171","#bf3729","#e69b00","#f5bb50","#ada43b","#355828"))
COL_DIR <- renoir[3]; COL_LN <- renoir[5]; COL_STAN <- renoir[8]; COL_LCVA <- renoir[1]

theme_paper <- function(base = 14) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray88", linewidth = 0.3),
        panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "lightgray", colour = "black"),
        strip.text = element_text(face = "bold", colour = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_blank())

# ---- Figure 1: runtime vs MCMC iterations ----------------------------------
sweep <- read.csv(file.path(DIR, "champs_time_vs_iter_results.csv"), stringsAsFactors = FALSE)

g <- sweep %>%
  filter(method %in% c("gibbs_dir", "gibbs_ln")) %>%
  group_by(method, iter) %>%
  summarise(time = mean(time_sec), sd = sd(time_sec), .groups = "drop") %>%
  mutate(sampler = recode(method, gibbs_dir = "Gibbs-Dir", gibbs_ln = "Gibbs-LN"))

stan_mean <- mean(sweep$time_sec[sweep$method == "prototype_stan_nuts"])

p1 <- ggplot(g, aes(iter, time, colour = sampler, shape = sampler)) +
  geom_hline(yintercept = stan_mean, linetype = "dashed", colour = COL_STAN, linewidth = 0.9) +
  annotate("text", x = max(g$iter), y = stan_mean + 2.5,
           label = sprintf("Stan / NUTS reference (4k iter, %.0fs)", stan_mean),
           hjust = 1, vjust = 0, colour = COL_STAN, size = 4, fontface = "bold") +
  geom_errorbar(aes(ymin = time - sd, ymax = time + sd), width = 250, linewidth = 0.6) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = c("Gibbs-Dir" = COL_DIR, "Gibbs-LN" = COL_LN)) +
  scale_shape_manual(values = c("Gibbs-Dir" = 16, "Gibbs-LN" = 15)) +
  scale_x_continuous(breaks = c(4000, 6000, 8000, 10000)) +
  coord_cartesian(ylim = c(0, stan_mean + 16)) +
  labs(x = "MCMC iterations", y = "Wall-clock time (seconds)") +
  theme_paper()

out1 <- file.path(DIR, "runtime_vs_iterations.pdf")
ggsave(out1, p1, width = 8, height = 5.2, dpi = 300); cat("wrote", out1, "\n")

# ---- Figure 2: per-step end-to-end breakdown -------------------------------
# Relabel to short, intuitive names; drop the redundant reference-Stan rows
# ("Prototype Stan ..." was the standalone validation model = same algorithm
# as the package's Stan option, so a single Stan/NUTS bar is shown).
squish  <- function(x) gsub("\\s+", " ", trimws(x))
relabel <- c(
  "LCVA train [KE]" = "LCVA train (KE)",  "LCVA pred [KE]" = "LCVA predict (KE)",
  "LCVA train [MZ]" = "LCVA train (MZ)",  "LCVA pred [MZ]" = "LCVA predict (MZ)",
  "LCVA train [ZA]" = "LCVA train (ZA)",  "LCVA pred [ZA]" = "LCVA predict (ZA)",
  "BFL stan (4k)" = "Stan/NUTS (4k)",
  "BFL gibbs Dirichlet (4k)"      = "Gibbs-Dir (4k)",
  "BFL gibbs Dirichlet (20k)"     = "Gibbs-Dir (20k)",
  "BFL gibbs LogisticNormal (4k)" = "Gibbs-LN (4k)",
  "BFL gibbs LogisticNormal (20k)"= "Gibbs-LN (20k)")
drop_steps <- squish(c("Prototype Stan compile", "Prototype Stan sampling", "TOTAL (timed steps)"))

ord <- c("LCVA train (KE)","LCVA predict (KE)","LCVA train (MZ)","LCVA predict (MZ)",
         "LCVA train (ZA)","LCVA predict (ZA)","Stan/NUTS (4k)",
         "Gibbs-Dir (4k)","Gibbs-Dir (20k)","Gibbs-LN (4k)","Gibbs-LN (20k)")

brk <- read.csv(file.path(DIR, "example_champs_timings.csv"), stringsAsFactors = FALSE)
brk$step <- squish(brk$step)
brk <- brk %>%
  filter(!step %in% drop_steps) %>%
  mutate(label = unname(relabel[step]),
         cat = case_when(grepl("LCVA", label)      ~ "LCVA (source)",
                         grepl("Stan", label)      ~ "Stan / NUTS",
                         grepl("Gibbs-Dir", label) ~ "Gibbs-Dir",
                         grepl("Gibbs-LN", label)  ~ "Gibbs-LN"),
         cat = factor(cat, levels = c("LCVA (source)","Stan / NUTS","Gibbs-Dir","Gibbs-LN")),
         label = factor(label, levels = rev(ord)))

fill_cols <- c("LCVA (source)" = COL_LCVA, "Stan / NUTS" = COL_STAN,
               "Gibbs-Dir" = COL_DIR, "Gibbs-LN" = COL_LN)

p2 <- ggplot(brk, aes(elapsed_sec, label, fill = cat)) +
  geom_col(colour = "black", linewidth = 0.4, width = 0.75) +
  geom_text(aes(label = sprintf("%.1fs", elapsed_sec)), hjust = -0.15, size = 4) +
  scale_fill_manual(values = fill_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Elapsed time (seconds)", y = NULL) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank())

out2 <- file.path(DIR, "runtime_breakdown.pdf")
ggsave(out2, p2, width = 9, height = 6, dpi = 300); cat("wrote", out2, "\n")
cat("Done ->", DIR, "\n")
