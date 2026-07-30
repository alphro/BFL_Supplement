# ==============================================================
# helpers.R --- data + LCVA + BFL glue for the PHMRC / AP reproducible example.
#
# Faithful port of the paper's four BFL variants (base, domain, partial, mix) for
# the NO within-target label-shift setting, distilled from the tested pipeline
# (Scripts/phmrc). Each SOURCE site trains a local LCVA model and shares only its
# posterior_phi; DOMAIN and MIX additionally train a target "self" model on the
# labeled target rows. The package (run_BFL / score_BFL) does the aggregation.
#
# Dependencies: LCVA (remotes::install_github("richardli/LCVA")) and BFL.
# The PHMRC data (phmrc_clean.csv) is bundled here; run from this folder.
# ==============================================================
suppressMessages({ library(LCVA); library(BFL) })

# ---- LCVA settings ----------------------------------------------------------
# Paper recipe is train Nitr = 2000, predict Nitr = 4000 (Burn_in 500, drop 2000).
# LCVA.pred(return_likelihood = TRUE) allocates an [Nitr x N x C] array, which can
# blow past R's vector-memory limit on a laptop (esp. with the 5 source fits kept
# in memory for local-avg). So the PREDICT draws are halved here (2000, drop 1000);
# the phi is still a mean over 1000 draws. Bump PRED Nitr/DROP back to 4000/2000 if
# you have the RAM (or raise the limit: Sys.setenv(R_MAX_VSIZE = ...)/mem.maxVSize).
LCVA_K          <- 5L
LCVA_TRAIN_ARGS <- list(model = "S", Nitr = 2000L, thin = 2L, nchain = 3L)
LCVA_PRED_ARGS  <- list(model = "C", Burn_in = 500L, Nitr = 2000L)
LCVA_DROP_DRAWS <- 1000L
MISS_PROP       <- 0.8   # unlabeled fraction (=> 20% labeled)

# =============================================================================
# 1. PHMRC data access (bare-minimum helpers, from BFL-VA-Sub)
# =============================================================================
get_sim_data_from_phmrc_by_location <- function(site, phmrc_csv = "phmrc_clean.csv") {
  phmrc  <- read.csv(phmrc_csv)
  causes <- unique(phmrc$cause)
  s      <- phmrc[phmrc$site == site, ]
  s$cause <- factor(s$cause, levels = causes)
  s$Y     <- as.integer(s$cause)
  list(data       = list(cause_level = causes, cause = s$cause, Y.t = s$Y,
                        X = s[, 1:168], G = rep(1, nrow(s))),
       data.truth = list(Y.t = s$Y), C = 34, NG = 1)
}

filter_sparse_causes <- function(sim_data, threshold = 0) {
  freq   <- table(sim_data$data.truth$Y.t)
  sparse <- as.numeric(names(freq[freq <= threshold]))
  keep   <- which(!sim_data$data.truth$Y.t %in% sparse)
  sim_data$data.truth$Y.t <- sim_data$data.truth$Y.t[keep]
  sim_data$data$Y.t       <- sim_data$data$Y.t[keep]
  sim_data$data$X         <- sim_data$data$X[keep, ]
  sim_data$data$G         <- sim_data$data$G[keep]
  uniq <- unique(sim_data$data.truth$Y.t)
  o2n  <- setNames(seq_along(uniq), uniq)
  n2o  <- setNames(uniq, seq_along(uniq))
  sim_data$data.truth$Y.t <- as.numeric(o2n[as.character(sim_data$data.truth$Y.t)])
  sim_data$data$Y.t       <- as.numeric(o2n[as.character(sim_data$data$Y.t)])
  sim_data$C <- length(uniq)
  list(filtered_data = sim_data, mapping_origin_to_new = o2n, mapping_new_to_origin = n2o)
}

# labeled/unlabeled split for the no-shift setting (20% labeled)
check_valid_sample <- function(Y, missing_idx) {
  all(unique(Y[missing_idx]) %in% unique(Y[-missing_idx]))
}
generate_missing_Yt <- function(Y.t, miss_prop) {
  num_samples <- round(miss_prop * length(Y.t))
  counts <- table(Y.t)
  small  <- as.numeric(names(counts[counts <= 10]))
  missing_indices <- integer(0)
  for (cod in small) {
    idx <- which(Y.t == cod); n <- length(idx)
    if (n > 1) {
      nm  <- min(round(miss_prop * n), n - 1)
      missing_indices <- c(missing_indices, idx[sample.int(n, nm)])
    }
  }
  repeat {
    gen  <- which(!Y.t %in% small)
    gmis <- sample(gen, num_samples - length(missing_indices), replace = FALSE)
    fin  <- c(missing_indices, gmis)
    if (check_valid_sample(Y.t, fin)) break
  }
  fin
}

# --- label-shift split generators (Cases II & III) ---------------------------
# MILD: draw labeled and unlabeled class prevalences from independent Dirichlet(1)
# vectors, then sample rows WITH replacement -> a STACKED frame with prevalence
# shift between the two halves. (rdirichlet(1, rep(1, K)) is just K iid Exp(1)
# draws normalised to sum 1; inlined here to avoid a gtools dependency.)
.rdirichlet1 <- function(K) { g <- rgamma(K, 1); g / sum(g) }
generate_missing_Yt_unbalanced <- function(Y.t, train_prop) {
  cats <- unique(Y.t); K <- length(cats)
  pi1  <- .rdirichlet1(K); pi2 <- .rdirichlet1(K)
  draw <- function(counts) unlist(mapply(function(cat, n) {
    idx <- which(Y.t == cat); if (length(idx) == 0) integer(0) else sample(idx, n, replace = TRUE)
  }, cats, counts))
  tr <- draw(rmultinom(1, round(train_prop * length(Y.t)),       pi1))
  te <- draw(rmultinom(1, round((1 - train_prop) * length(Y.t)), pi2))
  miss_cat <- setdiff(cats, unique(Y.t[tr]))          # guarantee every cause is labeled
  if (length(miss_cat)) tr <- c(tr, unlist(lapply(miss_cat, function(cat) which(Y.t == cat))))
  list(training_indices = tr, testing_indices = te)
}
# SEVERE: per-class Beta(0.2, 0.2) unlabeled fraction (bimodal -> some causes go
# almost fully unlabeled, others fully labeled), rescaled to the target unlabeled
# size. Disjoint, no replacement; samples POSITIONS so singleton causes are safe.
generate_missing_Yt_severe <- function(Y.t, miss_prop) {
  cats <- sort(unique(Y.t))
  target_u <- round(miss_prop * length(Y.t))
  raw    <- vapply(cats, function(c) round(sum(Y.t == c) * rbeta(1, 0.2, 0.2)), numeric(1))
  scaled <- round(raw * target_u / sum(raw))
  is.lab <- no.lab <- integer(0)
  for (i in seq_along(cats)) {
    idx <- which(Y.t == cats[i]); n_u <- min(length(idx), scaled[i])
    u   <- idx[sample.int(length(idx), n_u)]
    is.lab <- c(is.lab, setdiff(idx, u)); no.lab <- c(no.lab, u)
  }
  list(training_indices = is.lab, testing_indices = no.lab)
}

# make_split: single entry point for all three cases. Returns labeled/missing (as
# POSITIONS in the ctx frame) + row_map (original target row per frame row).
#   no_shift / severe : identity frame, disjoint partition of the original rows.
#   mild              : STACKED with-replacement frame; labeled/missing index INTO
#                       the stacked frame, row_map = c(labeled_rows, unlabeled_rows).
make_split <- function(shift, Y_new, seed, miss_prop = MISS_PROP) {
  set.seed(seed)
  if (shift == "no_shift") {
    missing <- sort(as.integer(generate_missing_Yt(Y_new, miss_prop)))
    labeled <- sort(setdiff(seq_along(Y_new), missing))
    return(list(labeled = labeled, missing = missing, row_map = seq_along(Y_new)))
  }
  if (shift == "mild") {
    res <- generate_missing_Yt_unbalanced(Y_new, 1 - miss_prop)   # arg = LABELED fraction
    tr  <- as.integer(res$training_indices); te <- as.integer(res$testing_indices)
    return(list(labeled = seq_along(tr), missing = length(tr) + seq_along(te),
                row_map = c(tr, te)))
  }
  res     <- generate_missing_Yt_severe(Y_new, miss_prop)         # severe
  labeled <- sort(as.integer(res$training_indices))
  missing <- sort(setdiff(seq_along(Y_new), labeled))
  list(labeled = labeled, missing = missing, row_map = seq_along(Y_new))
}

# labeled rows only, causes relabeled to a compact 1..C_obs space
get_observed_sim_data <- function(sim_data) {
  non_na <- which(!is.na(sim_data$data$Y.t))
  raw_Y  <- sim_data$data$Y.t[non_na]
  obs_causes <- sort(unique(raw_Y))
  remap  <- setNames(seq_along(obs_causes), obs_causes)
  Yc     <- as.integer(remap[as.character(raw_Y)])
  list(obs_data = list(
    data.truth = list(Y.t = Yc),
    data       = list(Y.t = Yc, X = sim_data$data$X[non_na, , drop = FALSE],
                      G = rep(1L, length(Yc))),
    NG = 1L, C = length(obs_causes), obs_causes = obs_causes))
}

# =============================================================================
# 2. PHMRC site accessors
# =============================================================================
load_site    <- function(site) filter_sparse_causes(get_sim_data_from_phmrc_by_location(site), 0)
site_X       <- function(s) { X <- as.matrix(s$filtered_data$data$X); dimnames(X) <- NULL; X }
site_Y_new   <- function(s) as.integer(s$filtered_data$data.truth$Y.t)
site_map_n2o <- function(s) s$mapping_new_to_origin
site_map_o2n <- function(s) s$mapping_origin_to_new
site_origin  <- function(s) as.numeric(s$mapping_new_to_origin)

# =============================================================================
# 3. LCVA fit + posterior phi + local-summary assembly
# =============================================================================
lcva_fit <- function(X, Y, seed) {
  a <- LCVA_TRAIN_ARGS
  LCVA::LCVA.train(X = as.matrix(X), Y = as.integer(Y), Domain = rep(1L, length(Y)),
                   K = LCVA_K, model = a$model, Nitr = a$Nitr, thin = a$thin,
                   nchain = a$nchain, seed = seed, verbose = FALSE)
}
lcva_phi <- function(fit, X_pred, seed, drop = LCVA_DROP_DRAWS) {
  a <- LCVA_PRED_ARGS; set.seed(seed)
  out <- LCVA::LCVA.pred(fit = fit, X_test = as.matrix(X_pred), model = a$model,
                         Burn_in = a$Burn_in, Nitr = a$Nitr,
                         return_likelihood = TRUE, verbose = FALSE)
  apply(out$x_given_y_prob[-seq_len(drop), , , drop = FALSE], c(2, 3), mean)   # N x C
}
make_local_summaries <- function(models, row_hash, P) {
  lapply(models, function(m) list(
    posterior_phi = m$phi, cause_ids = as.character(m$origin),
    target_info   = list(row_hash = row_hash, N = nrow(m$phi), P = P)))
}

# target self-model: LCVA trained on the labeled target rows, scored on predict_rows
target_self_model <- function(simf, labeled, missing, phi_seed, predict_rows, X, Y_new) {
  Y_part <- Y_new; Y_part[missing] <- NA
  obs <- get_observed_sim_data(list(data = list(X = X, Y.t = Y_part),
                                    data.truth = list(Y.t = Y_new)))$obs_data
  fit      <- lcva_fit(obs$data$X, obs$data.truth$Y.t, phi_seed)
  phi_full <- lcva_phi(fit, X, phi_seed, drop = LCVA_DROP_DRAWS)
  list(phi      = phi_full[predict_rows, , drop = FALSE],
       origin   = as.numeric(site_map_n2o(simf)[as.character(obs$obs_causes)]),
       newspace = as.numeric(obs$obs_causes))
}

# pad the target self-phi to the full target cause list (zero-fill absent causes)
expand_target_phi_full <- function(phi, origin, Y_truth_origin) {
  all_o <- sort(unique(as.numeric(Y_truth_origin)))
  out   <- matrix(0, nrow(phi), length(all_o), dimnames = list(NULL, as.character(all_o)))
  idx   <- match(as.numeric(origin), all_o)
  if (anyNA(idx)) stop("expand_target_phi_full: target cause(s) not in truth: ",
                       paste(origin[is.na(idx)], collapse = ", "))
  out[, idx] <- phi
  list(phi = out, origin = all_o)
}

# SEVERE placement: like expand_target_phi_full but keyed on the self-model's
# NEW-SPACE codes (the "self-generated" cause placement the paper's severe case
# uses), and a self-cause absent from the target truth is dropped, not an error.
expand_target_phi_selfgen_placement <- function(phi, newspace, Y_truth_origin) {
  all_o <- sort(unique(as.numeric(Y_truth_origin)))
  out   <- matrix(0, nrow(phi), length(all_o), dimnames = list(NULL, as.character(all_o)))
  idx   <- match(as.numeric(newspace), all_o); keep <- !is.na(idx)
  out[, idx[keep]] <- phi[, keep, drop = FALSE]
  list(phi = out, origin = all_o)
}

# split the labeled rows ~50/50 per cause (singletons -> domain) for MIX
.split_labeled_mix <- function(labeled, Y_new, seed) {
  set.seed(seed); dom <- par <- integer(0)
  for (c in unique(Y_new[labeled])) {
    idx <- labeled[Y_new[labeled] == c]
    if (length(idx) <= 1L) { dom <- c(dom, idx) }
    else { sd <- sample(idx, floor(length(idx) / 2)); dom <- c(dom, sd); par <- c(par, setdiff(idx, sd)) }
  }
  list(domain = sort(as.integer(dom)), partial = sort(as.integer(par)))
}

# =============================================================================
# 4. Prepare the AP context (loads sites, fits source LCVA, splits labels)
# =============================================================================
prepare_ap <- function(target  = "AP",
                       sources = c("Mexico", "Bohol", "Dar", "Pemba", "UP"),
                       seed = 1, shift = c("no_shift", "mild", "severe"),
                       miss_prop = MISS_PROP) {
  shift <- match.arg(shift)
  ap    <- load_site(target)
  Xf    <- site_X(ap); Yf <- site_Y_new(ap)
  Yof   <- as.numeric(site_map_n2o(ap)[as.character(Yf)])

  # split -> labeled/missing (POSITIONS in the frame) + row_map (original row per
  # frame row). Identity frame for no_shift/severe; stacked (duplicated) for mild.
  sp      <- make_split(shift, Yf, seed, miss_prop)
  row_map <- sp$row_map
  X              <- Xf[row_map, , drop = FALSE]
  Y_new          <- Yf[row_map]
  Y_truth_origin <- Yof[row_map]

  message("Fitting source LCVA models (this is the slow part)...")
  # for-loop (not lapply) so we can gc() the big [Nitr x N x C] prediction array
  # after each source, keeping peak memory to ~one source at a time.
  src <- setNames(vector("list", length(sources)), sources)
  for (s in sources) {
    sm  <- load_site(s)
    fit <- lcva_fit(site_X(sm), site_Y_new(sm), seed)
    phi <- lcva_phi(fit, Xf, seed)                                   # scored on ALL original AP rows
    src[[s]] <- list(phi = phi[row_map, , drop = FALSE], origin = site_origin(sm),  # sliced to frame
                     fit = fit, o2n = site_map_o2n(sm), n2o = site_map_n2o(sm))      # kept for local-avg
    rm(sm, fit, phi); gc(FALSE)
  }

  list(simf = ap, X = X, Y_new = Y_new, Y_truth_origin = Y_truth_origin,
       labeled = sp$labeled, missing = sp$missing, row_map = row_map,
       src_sites = sources, source_phi = src, shift = shift,
       target = target, phi_seed = seed, seed = seed,
       label_shift = shift != "no_shift",
       run_sampler = "gibbs", mcmc_args = list(iter = 2000, chains = 2, seed = seed),
       gibbs_args = list())
}

# =============================================================================
# 5. Per-variant INPUT builders. Each returns the exact arguments for a
#    run_BFL() + score_BFL() call, which the vignette/Rmd makes explicitly. This
#    keeps the LCVA/phi assembly here and the *package usage* visible in the Rmd.
# =============================================================================

# base: sources only, no target labels, all rows.
inputs_base <- function(ctx) {
  ls <- make_local_summaries(ctx$source_phi, BFL::compute_row_hashes(ctx$X), ncol(ctx$X))
  list(local_summaries = ls, X_target = ctx$X, Y_target = NULL, Y_add = NULL,
       Y_eval = as.character(ctx$Y_truth_origin), eval_idx = NULL)
}

# domain: sources + target self-model, on the unlabeled rows; labeled truth -> Y_add.
inputs_domain <- function(ctx) {
  rows <- ctx$missing
  tsm  <- target_self_model(ctx$simf, ctx$labeled, ctx$missing, ctx$phi_seed,
                            predict_rows = rows, X = ctx$X, Y_new = ctx$Y_new)
  tgt  <- expand_target_phi(ctx, tsm)
  models <- setNames(lapply(ctx$src_sites, function(s)
    list(phi = ctx$source_phi[[s]]$phi[rows, , drop = FALSE],
         origin = ctx$source_phi[[s]]$origin)), ctx$src_sites)
  models[[ctx$target]] <- list(phi = tgt$phi, origin = tgt$origin)
  ls <- make_local_summaries(models, BFL::compute_row_hashes(ctx$X[rows, , drop = FALSE]), ncol(ctx$X))
  list(local_summaries = ls, X_target = ctx$X[rows, , drop = FALSE],
       Y_target = NULL, Y_add = as.character(ctx$Y_truth_origin[ctx$labeled]),
       Y_eval = as.character(ctx$Y_truth_origin[rows]),
       eval_idx = if (isTRUE(ctx$label_shift)) seq_along(rows) else NULL)
}

# target self-phi placement: origin slots normally; NEW-SPACE (self-gen) placement
# under severe. Shared by domain and mix.
expand_target_phi <- function(ctx, tsm) {
  if (identical(ctx$shift, "severe"))
    expand_target_phi_selfgen_placement(tsm$phi, tsm$newspace, ctx$Y_truth_origin)
  else
    expand_target_phi_full(tsm$phi, tsm$origin, ctx$Y_truth_origin)
}

# partial: sources on all rows; observed labels fed into the model (NA on unlabeled).
inputs_partial <- function(ctx) {
  ls <- make_local_summaries(ctx$source_phi, BFL::compute_row_hashes(ctx$X), ncol(ctx$X))
  Y_target <- as.character(ctx$Y_truth_origin); Y_target[ctx$missing] <- NA
  list(local_summaries = ls, X_target = ctx$X, Y_target = Y_target, Y_add = NULL,
       Y_eval = as.character(ctx$Y_truth_origin), eval_idx = ctx$missing)
}

# mix: labeled rows split per cause; half -> self-model, half -> partial labels.
inputs_mix <- function(ctx) {
  sp  <- .split_labeled_mix(ctx$labeled, ctx$Y_new, ctx$seed)
  dom <- sp$domain; par <- sp$partial; mis <- ctx$missing
  stan_rows <- c(par, mis); n_par <- length(par); n_mis <- length(mis)
  tsm <- target_self_model(ctx$simf, labeled = dom, missing = setdiff(seq_along(ctx$Y_new), dom),
                           phi_seed = ctx$phi_seed, predict_rows = stan_rows, X = ctx$X, Y_new = ctx$Y_new)
  tgt <- expand_target_phi(ctx, tsm)
  models <- setNames(lapply(ctx$src_sites, function(s)
    list(phi = ctx$source_phi[[s]]$phi[stan_rows, , drop = FALSE],
         origin = ctx$source_phi[[s]]$origin)), ctx$src_sites)
  models[[ctx$target]] <- list(phi = tgt$phi, origin = tgt$origin)
  ls <- make_local_summaries(models, BFL::compute_row_hashes(ctx$X[stan_rows, , drop = FALSE]), ncol(ctx$X))
  Y_target <- rep(NA_character_, length(stan_rows))
  if (n_par > 0) Y_target[seq_len(n_par)] <- as.character(ctx$Y_truth_origin[par])
  list(local_summaries = ls, X_target = ctx$X[stan_rows, , drop = FALSE],
       Y_target = Y_target, Y_add = as.character(ctx$Y_truth_origin[dom]),
       Y_eval = as.character(ctx$Y_truth_origin[stan_rows]), eval_idx = n_par + seq_len(n_mis))
}

# =============================================================================
# 6. LCVA baselines (local-self, local-avg). These do NOT use run_BFL; they score
#    LCVA predictions directly, and are shown for comparison (as in the paper).
# =============================================================================
lcva_pred_full <- function(fit, X_pred, Y_test = NULL, seed) {
  a <- LCVA_PRED_ARGS; set.seed(seed)
  LCVA::LCVA.pred(fit = fit, X_test = as.matrix(X_pred), Y_test = Y_test,
                  model = a$model, Burn_in = a$Burn_in, Nitr = a$Nitr, verbose = FALSE)
}
# modal predicted class per observation across the post-burn-in draws
get_assignment <- function(Y_test, drop = LCVA_DROP_DRAWS) {
  d <- Y_test[-seq_len(drop), , drop = FALSE]
  apply(d, 2, function(col) as.integer(names(which.max(table(col)))))
}
# macro-average recall (balanced accuracy) over causes present in the truth
balanced_acc <- function(ytrue, ypred) {
  cs   <- sort(unique(ytrue))
  conf <- table(factor(ytrue, cs), factor(ypred, cs))
  rec  <- diag(conf) / rowSums(conf); rec[rowSums(conf) == 0] <- NA
  mean(rec, na.rm = TRUE)
}
# CSMF accuracy over the union of estimated/true causes (via the package metric)
csmf_acc_origin <- function(pi_hat, pi_true) {
  all_o <- union(names(pi_hat), names(pi_true))
  a <- setNames(rep(0, length(all_o)), all_o); a[names(pi_hat)]  <- pi_hat
  b <- setNames(rep(0, length(all_o)), all_o); b[names(pi_true)] <- pi_true
  BFL::CSMF_acc(a / sum(a), b / sum(b))
}
# true CSMF over the eval subset (no_shift: whole target; label_shift: unlabeled)
csmf_true_origin <- function(ctx, all_o) {
  idx <- if (isTRUE(ctx$label_shift)) ctx$missing else seq_along(ctx$Y_new)
  tab <- table(factor(ctx$Y_truth_origin[idx], levels = as.numeric(all_o)))
  setNames(as.numeric(tab / sum(tab)), as.character(all_o))
}

# local-self: the target's OWN LCVA (trained on its labeled rows), predicting the
# unlabeled rows. Returns c(top1, top1bal, csmf).
local_self <- function(ctx) {
  Y_part <- ctx$Y_new; Y_part[ctx$missing] <- NA
  obs <- get_observed_sim_data(list(data = list(X = ctx$X, Y.t = Y_part),
                                    data.truth = list(Y.t = ctx$Y_new)))$obs_data
  fit <- lcva_fit(obs$data$X, obs$data.truth$Y.t, ctx$phi_seed)
  out <- lcva_pred_full(fit, ctx$X[ctx$missing, , drop = FALSE], seed = ctx$phi_seed)
  obs_origin <- as.numeric(site_map_n2o(ctx$simf)[as.character(obs$obs_causes)])
  ypred <- obs_origin[get_assignment(out$Y_test)]
  ytrue <- ctx$Y_truth_origin[ctx$missing]
  pi_hat <- setNames(colMeans(out$pi_test[-seq_len(LCVA_DROP_DRAWS), , drop = FALSE]),
                     as.character(obs_origin))
  all_o  <- union(names(pi_hat), as.character(sort(unique(ytrue))))
  c(top1 = mean(ypred == ytrue), top1bal = balanced_acc(ytrue, ypred),
    csmf = csmf_acc_origin(pi_hat, csmf_true_origin(ctx, all_o)))
}

# local-avg: each SOURCE LCVA predicts the target (with the target's partial
# labels); the per-source metrics are averaged. Returns c(top1, top1bal, csmf).
local_avg <- function(ctx) {
  mis <- ctx$missing; ytrue <- ctx$Y_truth_origin[mis]
  per <- lapply(ctx$src_sites, function(s) {
    S  <- ctx$source_phi[[s]]
    Yp <- as.numeric(S$o2n[as.character(ctx$Y_truth_origin)]); Yp[mis] <- NA
    out  <- lcva_pred_full(S$fit, ctx$X, Y_test = Yp, seed = ctx$phi_seed)
    pred <- as.numeric(S$n2o[as.character(get_assignment(out$Y_test))])[mis]
    pi_c <- colMeans(out$pi_test[-seq_len(LCVA_DROP_DRAWS), , drop = FALSE])
    pi_hat <- setNames(pi_c, as.character(as.numeric(S$n2o[as.character(seq_along(pi_c))])))
    all_o  <- union(names(pi_hat), as.character(sort(unique(ytrue))))
    c(top1 = mean(pred == ytrue), top1bal = balanced_acc(ytrue, pred),
      csmf = csmf_acc_origin(pi_hat, csmf_true_origin(ctx, all_o)))
  })
  colMeans(do.call(rbind, per))
}
