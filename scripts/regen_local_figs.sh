#!/bin/bash
###############################################################################
# regen_local_figs.sh
# Regenerates all CSV-driven (local) figures at 12pt, straight into their final
# folders (PHMRC/, CHAMPS/, Supplemnt Figures/...). Run on your Mac where R +
# the packages (dplyr, ggplot2, tidyr, tibble, grid, patchwork) are installed.
#
#   cd ".../BFL Project/Results/Figures/scripts" && bash regen_local_figs.sh
#
# Does NOT touch figure_08 / figure_10 (those recompute from HB caches).
###############################################################################
set -e
cd "$(dirname "$0")"

SCRIPTS=(
  fig_phmrc.R                 # PHMRC/  phmrc_<metric>_<shift>.pdf
  fig_champs.R                # CHAMPS/ champs_<metric>.pdf + stacked
  fig_phmrc_compare_method.R    # Supplemnt Figures/{LCVA Pooled,flattened BFL,InsilicoVA}/
  fig_phmrc_compare_sampler.R                # CHAMPS/champs_batch8/champs_<metric>.pdf
)

for s in "${SCRIPTS[@]}"; do
  echo "==================== $s ===================="
  Rscript "$s"
done
echo "All local (CSV-driven) figures regenerated."
