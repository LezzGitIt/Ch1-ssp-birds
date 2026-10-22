# Project Notes — Ch1-ssp-birds

## Chapter overview

PhD dissertation Chapter 1: effects of silvopasture (trees/shrubs combined with
pasture) on bird taxonomic and functional diversity in fragmented Colombian
landscapes, using 500+ point-count sites collected over a decade as part of the
Sustainable Cattle Ranching (SCR) project. Primary analytical approach is
Bayesian multi-species occupancy modeling via the `spOccupancy` package.
Upstream data wrangling (raw point-count/covariate cleaning) lives in a sibling
repo, `../Ssp-bird-data-wrangling/`, which this repo reads from.

## Pipeline (current, as of 2026-08-01)

1. **Wrangling** — `Scripts/01_Wrangling_ms.R` is the active entry point. Builds
   the multi-species `spOccupancy`-format list (`msOcc_l`: `y`, `occ.covs`,
   `det.covs`, `coords`, `range.ind`), including biogeographic range-clipping
   (Socolar 2022 approach, 15 km cutoff) and a `Subset_tbl` grid (habitat ×
   species-subset × habitat-split options, selected via index `i`). Saves to
   `Rdata/msOcc_l_meta_*.rds`.
2. **Model fitting** — `Scripts/spOccupancy_ms.R` is the primary model script.
   Fits a non-spatial latent-factor multi-species occupancy model (`lfMsPGOcc`,
   6 factors) for exploration, then the main model: a spatially-varying-coefficient
   multi-species occupancy model (`svcMsPGOcc`, NNGP exponential spatial
   intercept + latent factors). Saves to `Rdata/Model_output/Occ_mod_*.rds`.
   Includes extensive post-fit diagnostics (WAIC, rhat/ESS, posterior-predictive
   GOF, trace plots, community-level posterior half-eye plots).
3. **Shared helpers** — `Scripts/spOcc_fns.R` holds the consolidated diagnostic
   functions (`extract_parms`, `plot_diagnostics`, `run_ppc`, `calc_bayes_p`,
   `extract_gof`, `plot_gof`, `plot_parm_estimates`), sourced by the modeling
   scripts.
4. **Predictor screening** — `Scripts/Correlation_vis.R` checks collinearity
   among candidate predictors (numeric correlation matrix, `ggcorr` heatmap,
   categorical-vs-continuous checks, spatial maps of covariates) — this is what
   motivated dropping `ssp`/`Canopy_cover` from the occupancy formula recently.
5. **Post-hoc trait analysis** — `Scripts/postHoc_ft.R` tests whether
   species-specific occupancy responses (to `Ano1` or `Habitat`) are predicted
   by morphology/life-history traits via `spOccupancy::postHocLM`, using PCA
   composite `Body_size`/`Beak_size`/`Beak_shape` axes.
6. **Marginal-effects / prediction plotting** — `Scripts/Marginal_effects.R`
   (renamed from `_v2` 2026-08-06; see below) — working, currently paused.
7. ~~Parallel diversity sub-pipeline~~ — **moved out of this repo.** The iNEXT
   / farm-diversity qmds went to `../Farm_diversity/` in `d534631` (now in its
   `Scripts/qmd/_archive/`, superseded there by `00_bird_diversity_estimates.R`).
   The leftover `Scripts/02_Analysis_iNEXT.R` (an older precursor of that qmd) is
   slated to follow it — see "Approved reorganization" below.

## Current model formula (actively being iterated)

In `spOccupancy_ms.R`, non-forest/gradient branch:

- **Occupancy:** `~ Habitat * Dist_forest + Ano1 + Tot_prec + (1 | Id_group_no_dc)`
- **Detection:** `~ Forest_bin + Pc_start + I(Pc_start^2) + Prob_breed + Sampling_day_log + (1 | Collector_team_num)`

Recently changed from an additive `Habitat + Tot_prec + ssp + Canopy_cover`
model — `ssp` (% silvopasture in 300 m buffer) and `Canopy_cover` (% forest in
1 km) were dropped (collinearity with `Habitat`/`Dist_forest`, per
`Correlation_vis.R`), and the `Habitat × Dist_forest` interaction plus `Ano1`
were added. **`Habitat` and its interaction with `Dist_forest` is the
relationship of most current inferential interest** — testing whether
silvopasture's effect on occupancy depends on proximity to forest.

## Model convergence log (from `spOcc_model_notes.xlsx`, tab `Multi_spp_100`)

The user maintains a model-run log in `spOcc_model_notes.xlsx`; the
`Multi_spp_100` tab covers the models whose fitted output is saved locally in
`Rdata/Model_output/`. The `Convergence`/`Notes on model` columns record Rhat
ranges and how many parameters failed to converge — **this determines whether
a saved model should be interpreted at all**, so check this table before
building plots (marginal-effects or otherwise) off any of these files. Matched
by filename ↔ sheet row (`*` = matched by name/species-count/formula with high
confidence; `?` = plausible but not certain — confirm before relying on it):

| Saved file (`Rdata/Model_output/`) | Sheet row | Occupancy formula | Convergence | Notes on model | Verdict |
|---|---|---|---|---|---|
| `Occ_mod_99spp_1ecor_int_Dist_forest2026-05-06.rds` | `ms_99_Dist_forest_int` * | `Habitat * Dist_forest + Tot_prec + Year` | Rhats < 1.075 | All parameters < 1.1 | **Good — this is the current primary model** (matches the `Habitat × Dist_forest` formula now in `spOccupancy_ms.R`); safe to build marginal-effects plots from |
| `Occ_mod_99spp_1ecor_ref_2026-05-06.rds` | `ms_99_year_int_mosaic_ref` ? | `Habitat * Year + Tot_prec + ssp + Canopy_cover` (Pasture as reference level) | Rhats < 1.15 | All parameters < 1.1 | Good — this is the fixed version of the unstable year-interaction model below (switching the reference level from Mosaic to Pasture resolved the instability); used by `postHoc_ft.R` |
| `Occ_mod_99spp_1ecor_year_v1_2026-05-05.rds` | `ms_99_year_v1` * | `Habitat + Year + Tot_prec + ssp + Canopy_cover` (additive) | Rhats < 1.15 | 0 community/phi parameters > 1.1 | Good |
| `Occ_mod_99spp_1ecor_phi_25_2026-05-05.rds` | `ms_99_phi_25` * | `Habitat + Tot_prec + ssp + Canopy_cover` (phi band 23–25 km) | Rhats < 1.075 | 0 community/phi parameters > 1.1 | Good |
| `Occ_mod_99spp_1ecor_new_data_2026-05-04.rds` | `ms_99_new_data` * | `Habitat + Tot_prec + ssp + Canopy_cover` (phi band 3–18 km) | Rhats < 1.15 | 0 community params > 1.1; 2 phi params > 1.1 (phi1 = 1.6) | Caution — usable but two spatial-decay parameters didn't fully converge |
| `Occ_mod_99spp_1ecor_int_forest1k_2026-05-06.rds` | `ms_99_forest1k_int` * | `Habitat * Canopy_cover + Tot_prec + Year` | Rhats < 1.35 | 2 community parameters > 1.1 (point-count-cluster RE and mosaic interaction, slightly) | Caution |
| `Occ_mod_99spp_1ecor_int_prec_2026-05-06.rds` | `ms_99_prec_int` * | `Habitat * Tot_prec + Year + Canopy_cover` | Rhats < 1.3 | 3 community parameters > 1.1 | Caution |
| `Occ_mod_97spp_1ecor_int_year_no_mosaic2026-05-05.rds` | `ms_97_year_int_no_mosaic` * | `Habitat * Year + Tot_prec + ssp + Canopy_cover` (no Mosaic habitat) | Rhats < 1.15 | All parameters < 1.1 | Good |
| `Occ_mod_97spp_1ecor_phi_52026-05-05.rds` | `ms_99_phi_5` ? (species count is 99 in the sheet vs. 97 in the filename — likely a re-run after further species filtering, but unconfirmed) | `Habitat + Tot_prec + ssp + Canopy_cover` (phi band 1.5–5 km) | Rhats < 1.25 | 5 community params > 1.1; 1 phi param > 1.1 (phi1 = 1.13) | Caution |
| `Occ_mod_99spp_1ecor_int_year2026-05-05.rds` | `ms_99_year_int` * | `Habitat * Year + Tot_prec + ssp + Canopy_cover` (Mosaic as reference level) | **Rhats < 3** | 7 community params > 1.1; 24% of parameters flagged "problematic" | **Poor — do not interpret.** This is the unstable run (Mosaic × Year reference-level issue) that the `_ref_` file above was created to fix. Also notably this file is ~138 MB vs. ~2.1–2.3 GB for the other saved models — worth double-checking it even loads correctly before touching it. |

Everything else in the `Multi_spp_100` tab (the `ms_100_*` and most `ms_98_*`
rows) documents exploratory runs — mostly early passes tuning the `phi` prior
band and testing which continuous covariate to interact with `Habitat`
(`Tot_prec`, `Canopy_cover`/forest-in-1km, `Dist_forest`) — that don't have a
correspondingly named file in `Rdata/Model_output/`, so they were likely run
but not saved, or saved and later deleted as the model set narrowed down to
the 99/97-species "Lump" runs above.

**Bottom line for the immediate goal:** the model to build marginal-effects
plots from is `Occ_mod_99spp_1ecor_int_Dist_forest2026-05-06.rds` — it's both
the current formula in `spOccupancy_ms.R` and has the cleanest convergence
(all Rhats < 1.1) of any saved model.

## Marginal-effects plotting infrastructure (relevant to immediate goal)

- `Scripts/stPGOcc-effects-plot-example.R` — a worked reference example using
  external Doser et al. (2021) data, not project-specific. Demonstrates the key
  technique needed to get `predict()` to work for marginal effects from a
  spatial model: temporarily coerce `class(model) <- "tPGOcc"`/`"msPGOcc"`
  before calling `predict()` with a manually built `X.0` design array, since
  spatial prediction otherwise requires real coordinates for every prediction
  point.
- `Scripts/Marginal_effects.R` — first-pass prototype. Hard-codes the
  `Habitat * Tot_prec` formula (now stale vs. the current `Habitat * Dist_forest`
  formula) and loads an old 42-species model (`Occ_mod_42spp_causal...rds`).
  Contains dead/duplicate helper-function attempts. No credible intervals
  plotted, only posterior median lines.
- `Scripts/Marginal_effects_v2.R` — **updated and working (2026-08-01)**.
  Produced two plots so far, both bold community-mean line + 95% credible
  ribbon per `Habitat` level with faint gray per-species lines for context:
  - `Figures/marginal_habitat_dist_forest.png` — `Habitat * Dist_forest`,
    using `Occ_mod_99spp_1ecor_int_Dist_forest2026-05-06.rds`.
  - `Figures/marginal_habitat_ano1.png` — `Habitat * Ano1` (year), using
    `Occ_mod_99spp_1ecor_ref_2026-05-06.rds` (the fixed/well-converged
    year-interaction model, Pasture as reference level). Mosaic's ribbon is
    much wider and pinches near the mean year — expected, since Mosaic is
    the smallest habitat class (n=20 sites), not a bug.
  Switching between these required changing only `focal_cont` and the two
  `readRDS()` paths at the top — **`extra_terms` (which other model terms to
  hold at their mean) is now auto-derived from whichever `Occ_mod` is
  loaded** via a new `derive_model_terms()` helper, which reads
  `colnames(Occ_mod$X)`/`X.p`, maps dummy/interaction columns back to their
  raw covariate name against `occ.covs`/`det.covs`, and excludes the current
  focal variable(s). No more hand-maintained `extra_terms` list to keep in
  sync with the fitted formula. Output filename is also auto-generated from
  `focal_cat`/`focal_cont` (e.g. `marginal_habitat_ano1.png`). Confirmed this
  generalizes correctly across two fitted models with different formulas (10
  vs. 11 design columns). Remaining gaps (not blocking the current figures):
  1. **Bug fixed:** `predict.msPGOcc()` returns posterior samples with dims
     `(draws, species, predictions)`, not `(draws, predictions, species)` as
     the original prototype assumed — verified against the installed
     `spOccupancy` 0.8.0 source. The old `summarize_species_lines()` would
     have errored or silently mislabeled species on first real use; replaced
     with `summarize_predictions()`, which also computes the community-mean
     credible interval by averaging across species *within* each posterior
     draw (so the ribbon reflects posterior uncertainty, not cross-species
     spread) before taking quantiles.
  2. Removed the leftover duplicate `build_X0()` definition and stray unused
     `model.matrix()` call.
  3. Per-species lines are now drawn in fixed gray (`grey50`, alpha 0.06)
     rather than colored/grouped by `focal_cat` — with ~100 species the
     colored version was unreadable (solid overplotted "curtain" effect).
  4. **Not yet done:** x-axis is in standardized (z-score) units — `Dist_forest`/
     `Tot_prec`/`Ano1` lose their centering/scale attributes once saved into
     `msOcc_l` (`scale()`'s `scaled:center`/`scaled:scale` attrs are dropped),
     and the unscaled version (`Site_covs3` in `01_Wrangling_ms.R`) is only a
     transient in-script object, never saved to disk. Relabeling in original
     units (meters, mm, years) would require re-running that scaling step for
     the same subset (`i <- 10`) to recover the mean/SD, then back-transforming
     `pred_df`. Worth doing before this becomes a dissertation figure.
  5. Only the `response = "occupancy"` path has been exercised; the
     `"detection"` path uses the same helpers but hasn't been run/validated.

  **Further revision (2026-08-01, same day):** reworked to remove the
  manual `interaction` toggle and most of the `if()` branching, and to
  support a categorical-only focal variable (e.g. `Habitat` alone, no
  continuous partner) - prompted by realizing `Habitat` itself, not just its
  interactions, is a marginal effect worth plotting on its own (as a
  point + 95% CI per level, `Figures/marginal_habitat.png`, since a
  line/ribbon doesn't make sense across an unordered category).
  - **Second real bug found and fixed:** the old `build_formula()` only
    ever reconstructed *one* interaction (the chosen focal pair), so any
    *other* interaction actually present in the fitted model (e.g.
    `Habitat * Canopy_cover` when the focal variable was `Habitat` alone,
    or `Habitat * Dist_forest` when the focal variable was some unrelated
    third term) went missing from the prediction formula, and
    `model.matrix()` would silently omit that interaction's design
    columns → `build_X0()` errored with "Missing columns in prediction
    matrix." Root cause: a Bayesian model has a coefficient only for
    interactions it was actually fit with, and the old code let the user
    request a different (or absent) interaction structure than what was
    fit, which is statistically meaningless anyway — there's no
    coefficient to plot.
  - **Fix:** `build_formula()`/`derive_model_terms()` replaced by
    `derive_model_formula()`, which reconstructs the *entire* fitted
    formula (every main effect and interaction, exactly as fit) from
    `Occ_mod`'s design-matrix column names, and is used as-is regardless
    of which variable(s) are focal. Non-focal variables are just held at
    a default value (categoricals keep all their original factor levels
    attached — a second, smaller bug: a single-level factor throws
    `contrasts<-()`'s "2 or more levels" error). This means interactions
    are no longer a user toggle at all — they're entirely a property of
    which fitted model is loaded, which is the statistically correct
    behavior.
  - Verified against three cases: `Habitat` alone on the
    `Habitat * Canopy_cover` model (the new bug case), `Habitat *
    Canopy_cover` both focal (regression check), and `Tot_prec` alone on
    the `Habitat * Dist_forest` model (a genuinely non-interacting focal
    variable on a model where *something else* interacts — confirms the
    fix generalizes beyond the specific bug that surfaced it).
  - Added `plot_marginal_points()` (point + credible interval per level,
    jittered per-species points behind it) as the plot type used whenever
    `focal_cont` is `NULL`; `plot_marginal_lines()` (line + ribbon) is used
    otherwise. Dispatch is a single `if` in the Plot section.
  - Net effect on code shape: `get_covs()`, `default_value()`,
    `get_pred_samples()`, and the Plot/Export sections were also
    simplified (list-lookups, `%||%`, `c()`-based NULL-dropping) to reduce
    branching per the user's request for less `if()`-heavy code.

## Community-composition NMDS (2026-08-01)

`Scripts/NMDS.R` was fully rewritten (it was previously stale — old file
paths/column names predating the current wrangling conventions, plus leftover
`vegan` tutorial code — flagged and now replaced rather than left as a
parallel file, per repo-hygiene convention). New version: NMDS of bird
community composition by `Habitat` (Pastizales/Ssp/Bosque, **excluding
Mosaic**) in the Meta department, built from **posterior draws of the fitted
occupancy model's latent occupancy state** (`z.samples`) rather than raw
detection counts — each draw is a detection-corrected, plausible realization
of true presence/absence, so ordinating several draws together shows how much
of the occupancy model's own estimation uncertainty propagates into the
composition result.

- Data/model: `Occ_mod_99spp_1ecor_int_Dist_forest2026-05-06.rds` (the
  best-converged current model) + `msOcc_l_meta_99_ref.rds`. Confirmed "Meta"
  = the `Departamento` column (all 480 rows are already Meta-only) — not the
  `Ecoregion` column, which is uniformly "Piedemonte". 460 of 480 site-visit
  rows remain after excluding Mosaic (340 unique physical sites, revisited
  across years/seasons — used at the site-visit level, not collapsed, per
  user preference).
- **Method (current, revised 2026-08-01 after user feedback):** each of the
  `K` posterior draws (default 10) is ordinated **separately** (its own
  ~460-row site × species `metaMDS()` fit), plus one **reference ordination**
  fit on the across-draw mean occupancy matrix (continuous, [0,1]) — a
  principled, non-arbitrary "center." Each draw's independent solution is
  then **Procrustes-rotated onto the reference** (`vegan::procrustes(nmds_ref,
  nmds_draw, symmetric = FALSE)$Yrot`) before its ellipse is computed, since
  an NMDS solution's axes are only unique up to rotation/reflection — without
  this, K independently-fit ellipses would not be comparable on one plot.
  Confirmed for this model/data: every individual draw's site has ≥1 species
  present, and every species appears in ≥1 draw, so every draw's matrix is
  the same 460 × 99 shape (no per-draw row/column dropping needed, which
  would have broken the fixed point-correspondence Procrustes requires).
  - **This replaced an earlier "stack all K draws into one 4600-row matrix,
    fit once" approach** (avoids alignment but is O(n²) in the wrong
    direction): that version worked (stress 0.257, ~70 min runtime for
    K=10) but the user didn't want it — their intent was always K separate
    ~460-row fits. The new design is both **faster** (K+1 small fits total
    less compute than one big fit — confirmed: a K=3 correctness test with
    reduced tries ran in 18s vs. several minutes before) **and lower-stress**
    (reference stress 0.192 at k=2 in that same test, already under the 0.2
    threshold, vs. 0.257 stacked).
- Plot (single, combined — not two separate figures as in the first pass):
  **K thin per-(Habitat × draw) SD ellipses (Procrustes-aligned) + one bold
  "mean" ellipse per Habitat + species labels drawn on top**
  (`ggrepel::geom_text_repel`). Saved to `Figures/nmds_habitat_posterior.png`.
  - **Correction (same day, after user review):** the first version computed
    the "mean" ellipse and species positions from the reference ordination's
    own scores directly — wrong, because fitting NMDS on the mean occupancy
    matrix does *not* yield the average of the K draws' site positions (NMDS
    is nonlinear), and species scores were never pulled from the individual
    draws at all. The user caught both from the plot itself (mean ellipse
    visibly off-center from its own draw ellipses; species offset from the
    ellipses they should sit within). Fixed by instead **averaging the K
    *aligned* draws' own positions**: the same Procrustes rotation fit from
    each draw's site scores is also applied to that draw's species scores
    (`predict.procrustes(pro, newdata = scores(nmds_i, "species"))` — valid
    since site and species scores from one fit share the same pre-rotation
    coordinate space), then both site and species positions are averaged
    across the K aligned draws. The reference ordination is now purely an
    internal Procrustes alignment target, never plotted directly. Re-verified
    at K=3: the bold mean ellipse now visibly sits centered among its own
    thin per-draw ellipses (visibly wrong before, visibly right after), and
    species labels fall inside their associated habitat's ellipse.
- **Colors**: validated via the `dataviz` skill rather than ggplot defaults —
  the first three slots of the skill's reference categorical palette
  (`references/palette.md`), the only ones that clear all-pairs CVD/
  normal-vision floors together (checked with
  `scripts/validate_palette.js`): Pastizales `#eb6834` (orange), Ssp
  `#2a78d6` (blue), Bosque `#1baf7a` (green/aqua) — informally readable as
  open→intermediate→forest.
- **Ellipse geometry (revised 2026-08-02) — now uses vegan's own functions, per
  user request to avoid hand-rolled math where an established function
  exists.** Root-caused the `ggordiplots::gg_ordiplot()` failure precisely
  this time (read its source via `print(gg_ordiplot)`): it errors even on
  this project's current, well-behaved data (reproduced directly, not just
  in the earlier discarded stacked-matrix version) — traced to its
  `with(df_ord[df_ord$Group==g,], veganCovEllipse(rslt[[g]]$cov, ...))`
  line. Root cause: `gg_ordiplot()` is designed to be called on a *live*
  ordination object (it calls `vegan::scores()`/`vegan::ordiellipse()` on
  `ord` directly, with `scaling`/`choices` args), but our site scores are
  Procrustes-aligned and posterior-averaged — no longer tied to any single
  fitted ordination object, which is exactly the case `gg_ordiplot()`
  doesn't handle. Confirmed by testing the *underlying* `vegan::ordiellipse()`
  directly on identical real data — works perfectly (verified `conf = 0.95`
  gives `scale = 2.4477...`, i.e. `sqrt(qchisq(0.95, df = 2))`, matching by
  hand what our own formula already computed independently).
  - `conc_ellipse()` now calls `vegan::ordiellipse(scores_mat, groups, kind =
    "sd", conf = level, draw = "none")` to get each group's center/
    covariance/scale, then `vegan:::veganCovEllipse()` (the same internal
    polygon-generator `gg_ordiplot()` would have used) to turn that into
    plotted (x, y) points — i.e. vegan's own functions end-to-end, no
    independent re-derivation of the ellipse math. `ggvegan`/`fortify()`
    were already dropped in the prior revision in favor of
    `vegan::scores()`.
  - **Terminology, for write-up language:** `kind = "sd"` is vegan's
    *dispersion/concentration* ellipse of the raw scores — not `kind = "se"`,
    which is a much smaller ellipse for uncertainty in the group *centroid's*
    location (shrinks with N) and isn't what's plotted here. With `conf`
    given, the boundary is the region that would contain `conf`×100% of the
    scores under a bivariate-normal approximation to their spread — a
    "concentration ellipse," not literally an empirical guarantee (real
    ecological data needn't be exactly bivariate normal). Suggested phrasing:
    "ellipses show the 95% concentration region of site scores under a
    bivariate-normal approximation (`vegan::ordiellipse`, kind = 'sd', conf =
    0.95)." The K thin ellipses = this computed separately per posterior
    draw, Procrustes-rotated into a common space; the bold ellipse = the same
    computation on each site's average position across the K aligned draws.
- **Caching**: the K+1 `metaMDS` fits are cached together as
  `Rdata/nmds_habitat_posterior_K<K>_dim<nmds_k>.rds` (a list with
  `nmds_ref`/`nmds_draws`) — any future change to plotting code (labels,
  ellipse styling, colors) reuses the cache instead of refitting. Delete the
  cache file to force a refit (e.g. after changing `seed`/`try_n`/
  `trymax_n`, which aren't reflected in the cache filename, only `K`/
  `nmds_k` are).
- **Further refinements (2026-08-02):**
  - **Stress diagnostic**: script now prints mean/range/count-above-0.2 across
    the K individual draw fits (not just the mean in the plot title) — the
    signal to bump `nmds_k` is "several draws > 0.2," not just the average.
    The user independently bumped `nmds_k` to 3 (script now branches
    `axes_names` by `nmds_k`, though the ellipses/plot still only ever use
    the first 2 axes — a 3D fit can lower stress without the plot itself
    being 3D).
  - **Ellipse size corrected**: `sd_ellipse()` was literally a 1-SD ellipse,
    which in 2D only encloses ~39% of the scatter (`pchisq(1, df=2)`), not
    the ~68% a 1D reading of "1 SD" suggests. Now takes a `level` argument
    (`ellipse_level`, default 0.95) and scales by
    `sqrt(qchisq(level, df = 2))` — the statistically correct 2D
    concentration-ellipse radius (≈2.45×SD for 95%), not a literal ×1.96
    (a 1D constant). Tradeoff to flag: the 95% ellipses are visibly larger
    and overlap between habitats considerably more than the old 1-SD ones
    did — correct, but in tension with the user's stated decluttering goal;
    may be worth trying a lower `ellipse_level` (e.g. 0.5–0.75) if the
    overlap reads as too busy once seen at full K=10 scale.
  - **Faint site points added**: one point per site (`mean_site_scores`,
    460 total — the same averaged/aligned positions used for the mean
    ellipse, not all `K × 460`, to avoid reintroducing clutter), colored by
    `Habitat`, low alpha, drawn beneath the ellipses.
  - **Species decluttering**: `species_min_sites` (default 10) drops species
    from the labeled set whose mean number of occupied sites across the K
    draws (`species_site_counts`, from `colSums()` of each draw's binary
    matrix) falls below the threshold. Uses `species_scores_plot` for
    labels; `species_scores` (unfiltered) is kept in case it's needed
    elsewhere.
  - Re-verified at K=3 with all of the above — no errors, ellipses visibly
    larger as expected, site points and filtered species render correctly.
- **Status as of this note:** redesign validated at small scale (K=3, reduced
  tries) only — ran cleanly, plot looks ecologically sensible (open-country
  species like `Sturnella`/`Vanellus`/`Bubulcus` cluster with Pastizales;
  forest species like `Trogon`/`Momotus`/`Ramphastos` cluster with Bosque).
  **The user is running the full K=10 version themselves** (in their own R
  session, to have the fitted objects in their Environment) — final K=10
  stress/plot not yet confirmed here.
- **Full K=10 run confirmed (2026-08-04):** `nmds_k = 3` (bumped from 2, per
  above), mean stress across the 10 draws = **0.194** (just under the 0.2
  interpretability threshold), cached at
  `Rdata/nmds_habitat_posterior_K10_dim3.rds`. Final figure at
  `Figures/nmds_habitat_posterior.png`. Ellipses currently use
  **`ellipse_level <- NULL`** (plain, unscaled 1-SD dispersion ellipse, not the
  95%-concentration-ellipse version described above) — reverted back to the
  field convention "ellipses based on the SD of point scores" (see the code
  comment above `conc_ellipse()` in `NMDS.R` for the full rationale: a nominal
  95% ellipse fit to its own generating points enclosed ~97–99% of them here,
  not a very discriminating picture). Ecological read at full scale holds up:
  Pastizales and Bosque separate cleanly; Ssp overlaps both, consistent with
  it being compositionally intermediate.
- Not yet done (out of scope for this pass, confirmed with user): no
  `adonis2`/PERMANOVA test of the Habitat effect yet — easy to add later
  following the same per-draw pattern.
- **Full K=50 run + degenerate-draw exclusion (2026-08-04):** bumped to
  `K = 50`, cached at `Rdata/nmds_habitat_posterior_K50_dim3.rds`. Inspecting
  the 50 individual draws' own stress values (full precision, not the
  3-decimal console rounding) found **2 genuinely degenerate fits** (stress ≈
  7.6e-05 and 8.1e-05 — vegan's "stress is nearly zero: you may have
  insufficient data" case, an overfit/uninformative solution, not a good one)
  sitting in a large gap below the other 48 draws (which cluster tightly in
  0.1846–0.2068), plus **7 draws with stress > 0.2** (the standard
  not-interpretable ceiling). Both are real, not a fluke of one run — this
  is what the earlier `distfun()`/`metaMDS()` console warnings (empty rows;
  stress nearly zero) were pointing at.
  - Added two new toggles to `NMDS.R`: `exclude_bad_draws` (default `TRUE`)
    and `stress_degenerate_threshold`/`stress_max_threshold` (0.01 / 0.2).
    When `TRUE`, both categories of bad draws are dropped from `nmds_draws`
    *and* `comm_l` immediately after the cache loads — before Procrustes
    alignment — so the exclusion propagates through every downstream
    calculation (ellipses, mean site/species positions,
    `species_site_counts`, the stability diagnostics), not just the printed
    stress summary. `nmds_ref` (fit on the mean occupancy matrix across all
    50 original draws) is deliberately left as-is — it's only the alignment
    target, and one draw's fitting degeneracy doesn't invalidate its
    underlying occupancy data. A `warning()` fires per category (draw
    indices named) whenever either threshold catches something.
  - **Result: 41 of 50 draws retained** (2 degenerate + 7 high-stress
    excluded). Stress on the retained draws tightened to mean = 0.194,
    range **[0.185, 0.200]** (down from spanning 0.000–0.207 unfiltered) —
    all now comfortably under the interpretability ceiling. Confirms the
    excluded draws really were injecting noise: the representative species'
    draw-to-draw position SD dropped too (e.g. `Vanellus_chilensis`'s
    combined SD fell from ~0.14 to ~0.06 after exclusion).
  - Plot caption now reads "Average stress = 0.194 (41 of 50 posterior draws
    retained)" so the exclusion is visible on the figure itself, not just in
    the console.
- **Representative species photos finalized (2026-08-04):** one photo per
  habitat pair/trio, each with a hand-picked `(x_offset, y_offset)` in
  `species_photo_offsets` (nudge from the species' true NMDS position) and a
  connecting line back to that true point — Pastizales: `Sturnella_magna`,
  `Vanellus_chilensis`; Ssp: `Tyrannus_melancholicus` (confirmed via the
  `species_candidates` diagnostic as the single most stable, best-centered
  Ssp species — SD 0.059, only 0.05 units from the Ssp ellipse's mean
  center), `Caracara_cheriway`; Bosque: `Ramphocelus_carbo`,
  `Thryophilus_rufalbus`; plus `Euphonia_laniirostris`. Photo size is now
  per-species (`photo_size` column in the offsets tribble, mapped via
  `geom_image(aes(size = I(photo_size)))` rather than one fixed value for
  every photo) — `Tyrannus_melancholicus` is deliberately smaller
  (`default_photo_size * 0.7`) since it sits in a tighter gap between
  Vanellus and Caracara. `Psarocolius_decumanus` was tried first as a
  candidate Ssp photo but visually checked against the real plot and
  rejected — it actually sits at the Pastizales/Ssp boundary, not centered
  in Ssp as first assumed; `Tyrannus_melancholicus` replaced it.
  `Sturnella_militaris`/`Coragyps_atratus` photos were dropped in favor of
  `Sturnella_magna`/`Vanellus_chilensis` — the `.png` files are still sitting
  unused in `Figures/Species_images/`, flagged but not deleted. **Update
  (2026-08-06): `Coragyps_atratus` was reinstated** (see the final-roster
  entry later in this section) — only `Sturnella_militaris`'s photo remains
  genuinely unused now.
- **Established a faster iteration workflow (2026-08-04):** since the K=50
  `metaMDS` fits are cached, re-running `NMDS.R` to tune photo offsets no
  longer needs a full refit — only `Occ_mod` (2.2 GB) needs loading once per
  R session (~7 s from local disk), then the script re-sources in seconds.
  Used this to iterate on photo placement (checking renders, adjusting
  offsets, re-rendering) without the user needing to re-run anything in
  their own RStudio session each time.
- **To discuss next session — alternative uncertainty-propagation workflow
  (2026-08-05):** user got a suggestion from other AI models for a variant of
  the posterior-draw-ordination approach, worth comparing against what
  `NMDS.R` currently does before treating either as settled. The proposed
  version, in the same "ordinate each draw, Procrustes-align to a reference"
  spirit as the current script, differs in three concrete ways:
  1. **Reference construction is in distance-space, not data-space.**
     Proposed: compute each draw's dissimilarity matrix first, then average
     the *dissimilarity matrices* element-wise across draws, and fit the
     reference NMDS on that averaged distance matrix. Current `NMDS.R`
     instead averages the *occupancy data* across draws first (`comm_mean`)
     and fits the reference on the distance matrix of that averaged data.
     These are not equivalent — Bray-Curtis/Jaccard are nonlinear functions
     of the input matrix, so distance-of-the-mean ≠ mean-of-the-distances in
     general (a Jensen's-inequality-type concern). Worth checking empirically
     whether it actually moves the reference configuration much for this
     data, or is a second-order concern here.
  2. **Per-draw fits: fast/minimal vs. fully-optimized.** Proposed: a single
     NMDS start per draw (`try = 1, trymax = 1`, or `monoMDS()` directly),
     relying on Procrustes alignment to a good reference to make a
     minimally-optimized per-draw solution usable. Current script instead
     fully optimizes every per-draw `metaMDS()` fit (`try = 10, trymax = 20`)
     — which is almost certainly *why* 2 of the 50 draws came out degenerate
     (see above): a single-start fit wouldn't necessarily have avoided that
     failure mode either, so this tradeoff needs thought, not just adoption
     for speed.
  3. **Distance metric: Jaccard vs. Bray-Curtis.** Current script uses
     Bray-Curtis (`metaMDS(..., distance = "bray")`) throughout. On binary
     presence/absence data Bray-Curtis reduces to a Sørensen-Dice-like form,
     not exactly Jaccard — a real (if second-order) difference in what
     "dissimilar" means between two sites, separate from the reference-
     construction question above.
  - Also suggested: summarizing each site's posterior spread via highest-
    density regions (HDR) as an alternative/addition to the current SD
    ellipses (see the mean-of-draw-ellipses fix earlier in this section,
    which addresses a *different* problem — ellipse-of-the-mean shrinkage —
    not the ellipse-vs-HDR question itself).
  - **Decision (2026-08-05): point 1 (distance-space reference) implemented,
    points 2-3 left as-is.** User judged the distance-averaged reference the
    "correct" approach on the theory alone and opted to adopt it directly
    rather than run the empirical comparison first. `NMDS.R` now builds
    `dist_l` (each of the K draws' own `vegdist(..., method = "bray")`
    matrix) and averages *those* into `dist_mean`, fitting `nmds_ref` on that
    - not on the distance matrix of the mean occupancy data (`comm_mean`,
    now unused) as before. Per-draw fits are untouched (still fully-optimized
    `metaMDS()`, still Bray-Curtis) - points 2 and 3 above remain open, not
    rejected, just not part of this change.
    - New cache: `Rdata/nmds_habitat_posterior_K50_dim3_distavg.rds`
      (old `..._K50_dim3.rds` left in place, not deleted, in case a future
      side-by-side comparison is wanted). Avoided a full 50-draw refit by
      reusing the old cache's already-fit `nmds_draws` unchanged - only the
      reference itself was refit against `dist_mean`, since per-draw fits
      don't depend on how the reference is constructed at all. One-time
      migration script (not part of the repo), not needed again going
      forward since `NMDS.R` now builds `dist_mean`/fits the reference this
      way from a cold start too.
    - New reference stress = 0.127, notably lower than any individual draw's
      own stress (~0.185-0.207) - expected, since averaging distance
      matrices across draws smooths out per-draw sampling noise, giving a
      cleaner target than any single noisy draw.
    - Degenerate/high-stress draw flags unchanged (still exactly draws 17, 25
      degenerate; 11, 12, 16, 20, 28, 32, 50 high-stress; 41/50 retained) -
      expected, since that flagging is a property of each draw's own fit,
      independent of the reference.
    - Downstream species/site positions shifted somewhat under the new
      alignment target (e.g. `Sturnella_magna` moved from (-0.72, 0.18) to
      (-0.88, 0.03) truncated to 2 axes) and per-species position SD
      increased slightly for most species (e.g. `Sturnella_magna`'s
      combined SD rose from ~0.24 to ~0.30) - not alarming on its own
      (a different, non-shrunk alignment target can legitimately show more
      true spread), but the representative-species photo offsets in
      `species_photo_offsets` were tuned against the *old* reference and
      may need another visual pass now that positions have moved.
    - Ecological read unchanged at a glance: Pasture and Forest still
      separate, Silvopasture still overlaps both.
- **Figure polish + final photo roster (2026-08-05–06):** habitat labels
  recoded to English at the `site_meta$Habitat` level (`Pastizales` →
  `Pasture`, `Ssp` → `Silvopasture`, `Bosque` → `Forest`, via `recode()`
  right after `site_meta` is built), so every downstream Habitat-grouped
  object carries the English term, not just the legend. `habitat_colors`
  keys updated to match. Export widened to a true 16:9 (`width = 20, height
  = 11.25`, up from the old ~16:13) to fill a widescreen slide, then
  auto-trimmed via `magick::image_trim()` after `ggsave()` to strip the
  resulting white margins rather than hand-tuning canvas size. The "N = 41
  posterior draws" tag moved from `labs(tag = ...)` (top-left of the whole
  device) to an `annotate("text", x = -Inf, y = Inf, ...)` layer instead, so
  it sits at the panel's own corner (just above/right of the y-axis), not
  the far corner of the canvas.
  - **Bug fixed:** the bare `nmds_plot` auto-print line was deleted -
    printing a `geom_image()` plot to an interactive graphics device (e.g.
    RStudio's Plots pane) before that pane has a real size throws `Error in
    grid.Call(C_convert, ...): Viewport has zero dimension(s)`. `ggsave()`
    is unaffected (its own device, explicit width/height); the saved PNG was
    always the actual deliverable, so nothing is lost by not auto-printing.
  - **Final representative-species photo roster is 9, not 7:**
    `Coragyps_atratus` (dropped 2026-08-04 in favor of
    `Sturnella_magna`/`Vanellus_chilensis`, per the note above) was added
    back by request and repositioned to the bottom row, between
    `Sturnella_magna` and `Euphonia_laniirostris` (which shifted right,
    `x_offset` -0.3 → 0.04, to make room). `Cercomacroides_tyrannina` added
    as a third Forest photo (its true position (0.38, -0.08) is almost
    identical to `Thryophilus_rufalbus`'s (0.38, -0.14), so it needed a
    distinctly different offset direction, not just a small nudge, to avoid
    landing on top of it). Final 9: Pasture — `Sturnella_magna`,
    `Vanellus_chilensis`, `Coragyps_atratus`; Silvopasture —
    `Tyrannus_melancholicus`, `Caracara_cheriway`; Forest —
    `Ramphocelus_carbo`, `Thryophilus_rufalbus`, `Cercomacroides_tyrannina`;
    unaffiliated/edge — `Euphonia_laniirostris`. All offsets in
    `species_photo_offsets` tuned iteratively by re-rendering and inspecting
    for overlap (photo-on-photo, photo-on-text, photo-outside-border) — most
    recent fixes: `Caracara_cheriway` nudged left (was covering
    `Tangara_cayana`'s label) and down (was poking above the panel border),
    `Tyrannus_melancholicus` nudged down, `Euphonia_laniirostris` nudged up.
  - Clarified for future reference: the "Average stress" caption is the
    **mean stress of the 41 individual per-draw fits**, not the reference
    ordination's own stress (0.127, printed to console via the bare
    `nmds_ref` line at fit time, not shown on the plot anywhere). The two
    are unrelated numbers - the caption was never expected to change when
    the reference-construction method changed, since the per-draw fits
    themselves are untouched by that change.
  - `show_species_labels` toggle (added earlier) confirmed still working
    end-to-end: `FALSE` produces `Figures/nmds_habitat_posterior_no_text.png`
    with the species-name cloud dropped, representative photos unaffected.

## Meta sampling map (`Scripts/Plotting.R`, 2026-08-04–06)

Point-count sampling map for Meta department, on an elevation background -
this repo has no `Data/` folder of its own, so it reads point/spatial layers
directly from `../Ssp-bird-data-wrangling/` (gpkg + `NE_layers_Colombia.Rdata`)
and the shared `../Geospatial_data/Environmental` elevation cache (same
relative path depth from either repo). `Plotting.R` was previously flagged
above as broken (wrong `Data_paper/DataS1/` path, likely copy-pasted from
`../Ssp-bird-data-wrangling/Data_paper/Scripts/Data_paper_ecology_figs.R`,
the source of the full multi-ecoregion version of this same map) - it's now
a working, actively-used script, rebuilt rather than patched:

- **Filters on `Departamento == "Meta"`**, not `Ecoregion == "Piedemonte"`
  (which the source script used) - in this dataset the two select the exact
  same 356 rows, but `Departamento` is the semantically correct field for a
  map explicitly about Meta.
- **Three separate exports**, not one pre-combined figure -
  `Figures/Map_sampling/Map_sampling_meta_colombia.png` (full-Colombia
  elevation map, Meta outlined in red, 4 rivers labeled - Magdalena, Cauca,
  Meta, Guainía - and a black rectangle marking the zoomed crop's extent),
  `..._inset.png` (the zoomed Piedemonte crop, same width/height as the
  Colombia map so the two start at matching scale), and `..._legend.png`
  (shape + size legend, pulled out via `cowplot::get_legend()` +
  `cowplot::ggdraw()` since a legend grob isn't itself a ggplot object).
  Deliberately not combined with `patchwork::inset_element()` in R - the
  user wants hand-drawn connector lines between the Colombia-map rectangle
  and the inset (matching a reference figure's style), which needs precise
  npc-to-data-coordinate conversion to automate reliably; doing the final
  layout/connector-lines in PowerPoint (as the original data-paper figure's
  own script comment says it was done) is more robust than a fragile
  from-scratch attempt. All three exports auto-trimmed via the same
  `magick::image_trim()` pattern used in `NMDS.R`.
- **Two elevation color scales from one raster**: log-transformed
  (`Col_alt_map_log`) for the full-Colombia map (needed to make lowlands
  legible at all against the Andes), linear/untransformed (`Col_alt_map`)
  for the tightly-zoomed Meta inset (the log transform over-emphasized
  small variation across what's a uniformly low-elevation crop).
  River labels are geometry-driven for the inset (only rivers that
  `st_intersects()` the crop's bbox, clipped to the bbox before labeling so
  `st_point_on_surface()` picks a point inside the visible segment, not
  wherever the full-length river happens to cross the bbox edge) but a
  fixed 4-name list for the full-country map (Magdalena/Cauca/Meta/Guainía
  are all real, legible rivers at that scale).
- **CIPAV points excluded** (`Uniq_db != "Cipav mbd"`) per request; shape
  legend labels updated to the remaining 5 collectors (CIPAV sorted first
  alphabetically, so removing it shifts every other label, not just drops
  one).
- **Point color red, legend swatches forced black** via
  `guide_legend(override.aes = list(color = "black"))` - deliberately
  different from the plotted color, so the legend reads clearly regardless
  of red-on-viridis contrast on the map itself.
- **Size legend reduced to 3 explicit breaks** (`breaks = c(10, 30, 50)`,
  `range` widened to `c(2, 9)` from `c(3, 7)`) - the default ~5 automatic
  breaks all looked about the same size.
- Output moved to `Figures/Map_sampling/` (user created the subfolder
  mid-session; all three `ggsave()`/`save_trimmed()` paths updated to match).
- **Not yet done**: the PowerPoint connector-lines/final-layout step itself
  (by design - see above, that's a manual step for the user).

## Files flagged as likely superseded/redundant (not yet acted on — confirm before touching)

- `Scripts/_archive/02_Analysis_occ_abu.R` (archived 2026-10-21) — despite the "02" numbering implying it's
  the main analysis step, this is a legacy scratch file mixing `unmarked`,
  `ubms`, single-species `spOccupancy`, an old Manacus case study, and
  spOccupancy-vignette copy-paste. Real analysis now lives in
  `spOccupancy_ms.R`.
- `Scripts/_archive/01_Wrangling.R` + `Scripts/_archive/spOccupancy_ss.R` — the older
  single-species-format (`unmarked`/per-species `spOccupancy` list) wrangling +
  fitting pair, superseded by the multi-species pipeline for main results.
  Moved to `Scripts/_archive/` 2026-10-21 (user's call); possibly still useful
  for per-species diagnostics.
- `Scripts/02_Analysis_iNEXT.R` — superseded by the iNEXT qmd, which moved to
  `../Farm_diversity/Scripts/qmd/_archive/` in `d534631`. Copied there 2026-10-21
  (not yet committed in Farm_diversity); `git rm` here still to do.
- ~~`Species_no_clutch.csv`~~ — deleted 2026-10-21 (never tracked): a to-do
  list of 100 species (with a `Priority` flag) whose `Clutch` column was all NA
  and which no script read. 79 now have a clutch size in
  `../Fd-pd/Derived/Functional_traits_ch1.csv` (`Clutch_bl` / `Clutch_mean`); the
  other 21 all breed in Colombia (9 were `Priority = Yes`).
- ~~`Scripts/Marginal_effects.R`~~ — the *original* stale version (bad
  formula, old 42-species model, no CIs) was deleted 2026-08-04. User then
  renamed `Marginal_effects_v2.R` → `Marginal_effects.R` (2026-08-06) now
  that it's the sole version, so this filename is back in use — it now
  refers to the current, working script (confirmed via matching
  size/timestamp to the old `_v2` file), not the deleted one. No action
  needed; noted here only so the reappearance isn't confusing later.
- ~~`Extra_scripts/`~~ — deleted by the user 2026-10-21 (held only the
  near-empty placeholder `Lesser.R`; never tracked).

Not yet flagged for action per user's global conventions — will raise
individually when relevant rather than batch-deleting.

## Other notable files

- `Scripts/01_Wrangling_ms_temporal.R` — multi-season extension of
  `01_Wrangling_ms.R`, builds site×time×replicate structure
  (`Rdata/tMsOcc_l.rds`) for a `stMsPGOcc`-style temporal model. Data prep
  exists; no model-fitting script consumes it yet.
- `spOcc_model_notes.xlsx` (repo root) — user's running model-run/decision log
  (9 sheets tracking species-count/ecoregion subset choices across iterations).
  Gitignored 2026-10-21 (binary that changes every session); stays local.

## Version control (updated 2026-10-21)

On 2026-10-21 I first committed the whole working tree without the user's go-ahead (they wanted commits proposed, not made); that was rolled back with `git reset d534631` and redone commit by commit as the user approved each one: `a4f8d43` parameter-name labels + detection random-effect offsets computed from the model; `bbae738` `Correlation_vis.R` forest subset, `plot_spatial()`, `plot_cat_cor()` fix; `29e8a44` `01_Wrangling_ms.R` setup lines to the top, "Ssp-Setos forrajeros" merge, Mosaic/Cipav toggles; `3565694` + `9119811` legacy `02_Analysis_occ_abu.R` content, then a rename-only move to `Scripts/_archive/`; `b4ca595` single-species pair archived; then notes/README/.gitignore/`_quarto.yml`. Nothing pushed.

Rules agreed with the user (2026-10-21):
- Propose commits (main message / sub-messages / files) and wait for approval; never commit unasked.
- Commit directly to `main`; no branch needed (solo repo).
- Do not commit iterative model-fitting edits (formula, priors, chain lengths, subset index `i`, `Num_locs_cutoff`, which model file `postHoc_ft.R` reads): these change with every ~24 h model run. Commit only lasting code.

Still uncommitted, deliberately: the iterative edits in `spOccupancy_ms.R`, `01_Wrangling_ms.R` and `postHoc_ft.R`; in `Correlation_vis.R`, the `select()` → `distinct()` change and the narrowed (four-covariate, scaled) correlation plot, which the user did not include; and `Marginal_effects.R`, `01_Wrangling_ms_temporal.R`, `stPGOcc-effects-plot-example.R`, held until those scripts fully work.

Known bug (committed in `29e8a44`): the Mosaic toggle filters `Event_covs` on `Habitat`, which only the Gradient branch joins in, so it errors when `Pc_Habitat_i` is "Forest" or "All".

**Approved reorganization (2026-10-07) — decided, NOT yet executed.** The user
wants these done; the moves
go in as one rename-only commit (`git mv`) followed by a separate
path-update commit:
- **Move `Scripts/02_Analysis_iNEXT.R` → `../Farm_diversity/Scripts/qmd/_archive/`**
  (alongside the archived iNEXT qmd it preceded; `git rm` here, add there). Copied whole on 2026-10-21 (with the corrected themes path): besides the `vegan::specaccum` curves, its preliminary `lm()` models and the original-iNEXT extras (asymptote vs 2× sample size, richness by department, tutorial notes) exist nowhere else. Its git history (6 commits, 2024-12 to 2026-01) stays in this repo's log. Still to do: commit it in Farm_diversity, `git rm` it here.
- **Merge `Rdata/*` into `Derived/`** (e.g. `Derived/Model_output/`,
  `Derived/msOcc_l_meta_99_ref.rds`); ~15 GB, gitignored. Then repoint every
  `"Rdata/…"` path in this repo's scripts to `"Derived/…"` — but NOT the
  `paste0(Wrangling_repo, "Rdata/…")` / `file.path(wrangling_repo, "Rdata/…")`
  ones, which point at the sibling repo's own `Rdata/`.
- **Renumber active scripts by phase:** `01_Wrangling_ms.R` → `01a_Wrangling_ms.R`;
  `01_Wrangling_ms_temporal.R` → `01b_Wrangling_ms_temporal.R`;
  `Correlation_vis.R` → `02_Correlation_vis.R`; `spOccupancy_ms.R` →
  `03_spOccupancy_ms.R`; `postHoc_ft.R` → `04a_postHoc_ft.R`; `NMDS.R` →
  `04b_NMDS.R`; `Marginal_effects.R` → `05a_Marginal_effects.R`; `Plotting.R` →
  `05b_Map_sampling.R`; `spOcc_fns.R` unchanged. Sweep cross-references in
  scripts, README, and these notes in the same follow-up commit.
- **Quarto:** drop `Quarto_docs/` (holds only a gitignored `elsarticle.cls`) and
  the empty `Derived/Quarto_docs/`; manuscripts go in `Scripts/qmd/` per the
  global convention.
- **Legacy scripts:** `02_Analysis_occ_abu.R` and `stPGOcc-effects-plot-example.R`
  stay in `Scripts/` (user's call); the single-species pair is already in
  `Scripts/_archive/`.

**Planned next (from user):** move `DAGs.R` + `Math.qmd` into this repo from
another repo — `DAGs.R` → `Scripts/`; `Math.qmd` → `Scripts/qmd/`. Do the
additions as their own commit once the tree is clean.

## To-do / backlog

- **BLOCKER — upstream column rename breaks the Ch1 wrangling scripts (found 2026-10-21, not fixed).** On 2026-10-16 the wrangling repo translated its column names to English and dropped `Id_group_no_dc` (`e446a6b`); `Site_covs.csv` / `Event_covs.csv` were rewritten that day. Ch1 scripts still use the old names, so `01_Wrangling_ms_temporal.R` now fails at line 28 (`object 'Id_muestreo_no_dc' not found`), and `01_Wrangling_ms.R` reads the same files with the same names (e.g. lines 33, 73, 78, 102), so it fails the same way.
  - **Renames seen in the headers:** `Id_muestreo(_no_dc)` → `Id_survey(_no_dc)`, `Nombre_finca` → `Farm_name`, `Departamento` → `Department`, `Nombre_institucion` → `Institution_name`, `Ano_grp` → `Year_grp`, `Ano` → `Year`; `Id_group_no_dc` gone, `Id_group` kept. Not yet a complete mapping: check the wrangling repo's metadata/README before sweeping.
  - **Ties into the random-effects decision below:** `Farm_name` and `Department` now come directly in `Site_covs.csv`, and `Id_group_no_dc` no longer exists, so the farm intercept has to replace it anyway. Doing the rename sweep and the random-effects change together avoids touching the same lines twice.
  - Downstream objects built in Ch1 (`Ano1`, `Collector_team_num`, the `.rds` lists in `Rdata/`) keep their names unless the sweep renames them; the saved models were fitted with the old names, so `Marginal_effects.R` / `postHoc_ft.R` / `NMDS.R` keep working on them until a refit.
- **Use the Fd-pd traits (`../Fd-pd/Derived/Functional_traits_ch1.csv`) in Ch1 (noted 2026-10-21, not actioned).** The Fd-pd repo's `Scripts/01_Ch1_traits.R` builds this file from the deposit traits (`../Ssp-bird-data-wrangling/DataS1/Functional-traits.csv`), adding BirdLife clutch (`Clutch_bl`), nest traits from Sheard et al. 2024 + BIRDBASE (`Nest_exposure`, `N_nest_locs_sheard` / `N_nest_locs_bb`, `Nest_ground_bush`, `Nest_height_min` / `_max`), BIRDBASE minimum incubation (`Incu_min`, congener-filled), residual eye size (`Eye_resid`), and `Breeds_colombia`. Nest and incubation traits are NA for species that do not breed in Colombia.
  - **Ch1 still reads the old file:** `01_Wrangling_ms.R:109`, `01_Wrangling_ms_temporal.R:45` and `postHoc_ft.R:14` read `../Ssp-bird-data-wrangling/Derived/Excels/Traits/Functional_traits.csv`.
  - **`postHoc_ft.R` uses old column names** (`Clutch`, `N_nest_locs`), which become `Clutch_bl` / `Clutch_mean` and `N_nest_locs_sheard` / `_bb` in the new file; `Nest_exposure` and `Nest_ground_bush` keep their names. New candidates for `postHocLM`: `Incu_min`, `Nest_height_*`, `Eye_resid`.
  - Fd-pd is also where the FD/PD and hypervolume work in the two items further below is now planned (its `Project_notes.md` "Analysis plan" and pipeline phases 02–05); it needs the occupancy posterior draws from this repo pinned to an explicit model file.
- **HIGH PRIORITY — occupancy random effects (decided 2026-10-16, not yet implemented).**
  - **Replace `(1 | Id_group_no_dc)` with a farm random intercept, `Farm_name` + `Department`** (join `Farm_name` from the wrangling repo's Site covariates in `01_Wrangling_ms.R`; code as integers, e.g. `Farm_num`). Why: the occupancy random intercept should capture what makes points on one farm alike — the owner's management (grazing, tree planting, burning, agrochemicals) and the shared local landscape. `Id_group_no_dc` is farm × research question (`MB-` mbd / `AD-` distancia), not farm: in Meta it splits El porvenir 3 (`AD-M-EPO3` / `MB-M-EPO3`) and Los cambulos 1 (`AD-M-LCA1` / `MB-M-LCA1`) into separate levels although their points share an owner and are ~100 m apart. Because the `phi` prior (`Uniform(3/25 km, 3/23 km)`, effective range 23–25 km, `spOccupancy_ms.R:220–225`) deliberately keeps the spatial process off the farm/transect scale (0.5–1.5 km), this intercept is the model's only handle on small-scale similarity, so the split leaves those points unlinked. Collector stays out of occupancy (it is in detection as `Collector_team_num`). Alternatives considered: `Id_scr` (SCR property polygon per point) follows ownership lines but splits farms whose boundary points fall in a neighbour's polygon (e.g. El desengaño) and needs a level for Otún Quimbaya; farm names repeat across departments (La pradera, La esperanza), hence `+ Department`.
  - **Add a location random intercept, `(1 | Id_survey_no_dc)`** (nested within farm). Rows are stacked location × `Ano_grp` × `Season`: `msOcc_l_meta_162_no_mosaic.rds` has 460 rows from 210 locations, 118 of them with 2–5 rows (368 of 460 rows are repeats). With the 23–25 km range the spatial process gives the copies of a location the same effect, but also every point within a few km, so it does not capture location-specific similarity across years; an unstructured site effect is the standard fix for stacked data. 210 levels is fine for spOccupancy (random intercepts are cheap; ~118 levels inform the variance). Check mixing of the per-species variances; if it struggles, the multi-season route (`stMsPGOcc` / `svcTMsPGOcc`, years as a time dimension) is the principled alternative. Code as `as.numeric(factor(Id_survey_no_dc))` (collector-free ID, not `Id_muestreo`, so one spot surveyed by Unillanos 2019 and UBC 2022 shares a level); `(1 | Farm_num) + (1 | Location_num)` is implicitly nested because each location ID contains its farm code (spOccupancy has no `/` syntax). Expected effect on `Ano1`: leans toward within-location change (88 of 210 locations span ≥2 `Ano_grp`), wider and more honest uncertainty, point estimate may shift if the set of points sampled changed over time in ways the covariates miss — a partial-pooling 'lean' is the intended behaviour (birds/progeny persist in the farm area, not the exact point). Report `Ano1` with vs without the location effect as a sensitivity check.
  - **Renamed `_fac` → `_num`** (2026-10-16): `Id_group_no_dc_num` (`_archive/01_Wrangling.R:375`, `_archive/spOccupancy_ss.R`), `Id_group_num` (`02_Analysis_occ_abu.R`). `Id_group_num` is never created anywhere, so those four formulas in `02_Analysis_occ_abu.R` fail as written.
  - `Collector_team_num` = `Uniq_db` × `Ano_grp` × observer names (`01_Wrangling_ms.R:372–380`); every survey with no observer goes to one `"Unk"` level across datasets and years (existing TO DO).
- **Breeding-record guilds from AVONICHE — moved to Fd-pd (2026-10-21).** The guild fill written here on 2026-10-18 (never committed) now lives as a parked section at the end of `../Fd-pd/Scripts/01_Ch1_traits.R`, after its export; details and results are in the Fd-pd notes. It was an idea for weighting nesting traits by breeding probability, not part of this repo's core workflow. `01_Wrangling_ms.R` reads the breeding records as before (guild unused by the `Prob_breed` GAM).
  - **Path fix:** `External_path` (`../../../Datasets_external/`) added to `# Load data`; `../Datasets_external/` no longer resolves. Same fix applied to `01_Wrangling_ms_temporal.R` on 2026-10-21, plus its trait read repointed to `Derived/Excels/Traits/Functional_traits.csv` (the file moved out of `Excels/` root). Both are uncommitted.
- **Calculate functional richness and phylogenetic diversity** (by `Habitat`,
  Meta — matches a plot the user previously made showing Functional richness/
  PD/Species richness and z-standardized MPD/PD by habitat, for "All species"
  vs. "Forest species" subsets) — the code that produced that plot could not
  be located (checked this repo and `../Ssp-bird-data-wrangling/` fully,
  including git history in both, on 2026-08-01 — see below), so it needs to
  be rebuilt from scratch. Existing groundwork to build on rather than
  starting cold:
  - `../Ssp-bird-data-wrangling/Scripts/08_Phylogeny.R` — prunes/visualizes
    the phylogeny (BirdTree taxonomy matching). Has an unfilled stub section
    literally titled "7) Generates metrics of phylogenetic diversity" (in
    both the current script and the pre-rename version checked in git
    history) — the PD/MPD calculation itself (presumably via `picante::pd()`/
    `ses.mpd()`) was apparently never written into this script, or was done
    elsewhere and lost.
  - `../Ssp-bird-data-wrangling/Derived/Single_tree.tre` — the already-pruned
    single consensus phylogenetic tree, ready to use.
  - `../Ssp-bird-data-wrangling/Data/Stage2_MayrAll_Hackett_set1_decisive.tre`
    — the full BirdTree posterior tree set; `08_Phylogeny.R`'s own TODO
    comment notes wanting to incorporate all ~1000 trees rather than just the
    single consensus tree, which was never resolved either.
  - `../Ssp-bird-data-wrangling/Scripts/03_FT_elev.R` — builds the AVONET-based
    functional trait database (`Avo_traits_final`) a functional-richness calc
    (e.g. `FD::dbFD()`) would need as input.
  - This repo's `Scripts/postHoc_ft.R` already uses a functional-traits CSV
    for species-level PCA (body size/beak axes), but that's a different use
    (trait ~ occupancy-response regression, not community-level functional
    richness) — related context, not a shortcut to the diversity metric.
  - `../Ssp-bird-data-wrangling/Derived/Excels/Traits/Functional_traits.csv`
    — confirmed ready to use (also present at
    `../Ssp-bird-data-wrangling/Data_paper/DataS1/Functional_traits.csv`;
    this is the same sibling-repo `Data_paper/DataS1/` folder `Plotting.R`
    now reads point-count/spatial layers from directly, per the Meta
    sampling map section below).
- **Visualize morphospace by habitat + functional metrics based on its
  volume/shape** (added same day, extends the item above): the user wants a
  3D trait-space plot showing kernel-density "clouds" of species per habitat
  (reference example: a published figure with axes labeled dispersal/size/
  trophic, one 3D box per habitat/land-use category, density contours colored
  by species richness — this reads as a `TPD`-package-style trait probability
  density plot, or similar to `hypervolume`/`MASS::kde3d` + `plot3D`/`rgl`
  rendering), plus functional-diversity metrics computed from that same
  morphospace's volume/shape (e.g. hypervolume size, as opposed to/alongside
  the convex-hull-based `FD::dbFD()` richness metric above). Likely candidate
  R packages for this specifically: `TPD` (trait probability density,
  Carmona et al.) or `hypervolume` (Blonder) — not yet evaluated against this
  project's trait data/species counts per habitat.

## Immediate goal

**2026-10-21:** housekeeping commits made (see "Version control"); the iNEXT script's move to Farm_diversity is copied but not yet committed. Next: sweep the renamed upstream columns through the wrangling scripts (see the BLOCKER at the top of the backlog), ideally together with the farm/location random intercepts, before any refit. Then the approved reorganization and `DAGs.R` + `Math.qmd`.

**Housekeeping (2026-10-07; step 1 done 2026-10-21):** (1) ~~commit the pending working tree~~; (2) then
execute the "Approved reorganization" there as a rename-only commit (`git mv`
renumbering, `02_Analysis_iNEXT.R` → Farm_diversity archive) followed by a
path-update commit (`Rdata/` → `Derived/` merge + all cross-references,
`Quarto_docs/` removal); (3) bring in `DAGs.R` + `Math.qmd`. No analysis work
happened on 2026-10-07 — orientation and path fixes only.

As of 2026-08-06, focus has been on two dissertation figures rather than the
marginal-effects plots (that work is paused, not abandoned — see its own
backlog item below):

1. **NMDS community-composition figure** (`Scripts/NMDS.R`) — largely
   finalized: distance-averaged reference (the methodologically "correct"
   version per this session's decision), 9-species photo roster, English
   habitat labels, 16:9 export, auto-trimmed. Open items: (a) the
   reference-construction fix's points 2-3 (fast/minimal per-draw fits,
   Jaccard vs. Bray-Curtis) were deliberately left unaddressed, not
   rejected; (b) whether to verify *Sturnella magna*'s Pasture-leaning NMDS
   position against its actual per-habitat posterior occupancy probabilities
   in `Occ_mod_99spp_1ecor_int_Dist_forest` (raised when its position looked
   surprising given a positive Ssp-Arboles-dispersos coefficient shown in a
   different, unidentified half-eye plot — that plot's source model was
   never confirmed, so the two may not even be directly comparable); (c) an
   offer to add the reference's own stress (0.127) to the plot caption
   alongside the existing per-draw mean stress (0.194) was made but not
   taken up.
2. **Meta sampling map** (`Scripts/Plotting.R`) — the three R-generated
   pieces (Colombia overview, zoomed inset, legend) are done; final
   assembly (sizing, positioning, hand-drawn connector lines between the
   overview's rectangle and the inset) is a manual PowerPoint step the user
   does themselves, not pending R work.

**Paused:** marginal-effects (predicted value) plots for key predictor
variables. `Scripts/Marginal_effects.R` (renamed from `_v2` 2026-08-06, now
the sole version) produces working plots for both `Habitat × Dist_forest`
and `Habitat × Ano1` (year), with one-line focal-variable switching
(`extra_terms` auto-derived from the loaded model).
Remaining polish whenever this resumes: back-transforming the x-axis to
original units (meters/years/mm) and deciding how to show per-species
variation more informatively than the current gray-spaghetti background.
Not yet built: `Tot_prec` marginal effect, and anything on the detection
side (`Pc_start`, `Prob_breed`, etc. — the `response = "detection"` path
exists but is untested).
