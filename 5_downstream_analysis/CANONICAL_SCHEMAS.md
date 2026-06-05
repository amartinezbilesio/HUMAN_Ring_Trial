# Canonical annotation schemas (Phase 1 contract)

The two CSVs below are the **single source of truth** for downstream
annotation analyses. Producers write them; consumers read them. No
analysis code in a producer, no extraction/SIRIUS/peak-picking in a
consumer.

```
sirius_annotation_all_ms2.qmd ─┐
                               ├─> annotation_identity.csv      ─> 9_per_compound_diagnostics.qmd
                               └─> annotation_chrom_metrics.csv ─> 8_data_level_progression.qmd
                                                                  (P4 reproducibility)
```

Both files live next to the SIRIUS producer at
`3_annotation_manual/HE/3_Sirius_curation/results_all_ms2/`.

---

## `annotation_identity.csv`

One row per **(mixture, lab, compound, polarity)**. MS2-based identity
of expected compounds.

| Column                  | Type      | Source                                  | Notes |
|-------------------------|-----------|-----------------------------------------|-------|
| `mixture`               | character | producer §4                             | e.g. `"1.1"` |
| `lab`                   | character | producer §3                             | one of `afekta, cembio, hmgu, icl` |
| `compound`              | character | `expected_all.compound` (ChEBI name)    | the matched expected compound |
| `inchikey_2d`           | character | first 14 chars of SIRIUS `inchiKey`     | join key against expected |
| `polarity`              | character | producer §3                             | `"pos"` / `"neg"` |
| `best_formula_rank`     | numeric   | producer §4 `min(rank.x)`               | NA if not in candidate list |
| `best_structure_rank`   | numeric   | producer §4 `min(rank.y)`               | NA if not in candidate list |
| `adducts`               | character | producer §4                             | `;`-separated list (e.g. `"[M+H]+;[M+Na]+"`) |
| `molecularFormula`      | character | producer §4                             | first candidate's formula |
| `expected_formula`      | character | `expected_all`                          | from `standards.xlsx` |
| `expected_inchikey`     | character | `expected_all`                          | full InChIKey |
| `mz_observed`           | numeric   | `all_summaries$ionMass` (best-csiScore candidate per row) | precursor m/z of the candidate that supplied the COSMIC scores; basis for inter-lab m/z deviation |
| `confidenceExactMatch`  | numeric   | `all_summaries$confidenceExactMatch`    | COSMIC exact-match score (best-csiScore candidate per row) |
| `confidenceApproxMatch` | numeric   | `all_summaries$confidenceApproxMatch`   | COSMIC approximate-match score (best-csiScore candidate per row) |
| `confidence_tier`       | character | P1 placeholder = `"unfiltered"`         | P5 overwrites with `"A"` / `"B"` / `"C"`; consumers filter on this column (default keeps all) |

**Key:** `(mixture, lab, compound, polarity)`.

**Why COSMIC scores are aggregated per row:** SIRIUS reports COSMIC at
the candidate level; identity is one-per-compound, so the producer picks
the candidate with the best `csiScore` per `(mixture, lab, compound,
polarity)` and carries its COSMIC scores. This matches §4's existing
`slice_max(csiScore, n = 1)` pattern.

---

## `annotation_chrom_metrics.csv`

One row per **(mixture, lab, polarity, compound, adduct, source)**. MS1
peak metrics from the SNR-aware single-pass EIC pipeline (producer §7).

| Column             | Type      | Source                    | Notes |
|--------------------|-----------|---------------------------|-------|
| `mixture`          | character | producer §8               | |
| `lab`              | character | producer §8               | |
| `polarity`         | character | producer §8               | `"pos"` / `"neg"` |
| `compound`         | character | producer §8               | matched expected name |
| `adduct`           | character | producer §8               | single adduct (one row per adduct) |
| `source`           | character | producer §8               | `"own"` = lab's MS2 annotation; `"consensus"` = MS1-only borrow of median (m/z, RT) |
| `areaUnderTic`     | numeric   | producer §7.6 picker      | trapezoidal area between valleys; **NA when `detection != "detected"`** |
| `xicFwhm`          | numeric   | producer §7.6 picker      | full-width-half-maximum (s); **NA when `detection != "detected"`** |
| `maxIntensity`     | numeric   | producer §7.6 picker      | raw apex intensity; **primary abundance metric**; **NA when `detection != "detected"`** |
| `apex_rt`          | numeric   | producer §7.6 picker      | apex_rt of the (real or noise) maximum inside the RT prior window; **populated even when below_detection** so consumers can inspect noise-grab locations |
| `snr`              | numeric   | producer §7.6 picker      | `(apex_int − baseline) / noise`; `noise = max(MAD, IQR/2, 1)` of tail at `|rt − apex| > 25 s` |
| `baseline`         | numeric   | producer §7.6 picker      | median of tail intensities |
| `noise`            | numeric   | producer §7.6 picker      | floor-1 noise estimate (see §7.4 guard #1) |
| `n_scans_in_shape` | integer   | producer §7.6 picker      | count of scans within apex_idx ± 2 with intensity ≥ 30 % of apex (§7.4 guard #2) |
| `prior_kind`       | character | producer §7.5 priors      | `"own"` (lab's SIRIUS RT, ±8 s search) or `"consensus"` (cross-lab median, ±15 s); see §7.5 for the demotion guards |
| `detection`        | character | producer §7.6 picker      | `"detected"` / `"below_detection"` / `"no_data"`; see §7.4 for the three-guard logic |

**Key:** `(mixture, lab, polarity, compound, adduct, source)`.

**Why `maxIntensity` is primary:** scan rates differ ~15× across labs
(cembio 1.5 s vs icl 0.1 s); area is scan-count- and boundary-sensitive,
so it is kept as a secondary metric only.

**Why both `own` and `consensus` rows exist for the same compound:** the
"consensus" row is a propagation of the median annotation (m/z, RT) from
labs that did annotate, into labs that did not. Lets P4 see cross-lab
abundance for the 275 propagated compounds without requiring MS2-in-all-4.

**Why `detection` is on the row, not via filtering:** consumers that
filter `!is.na(maxIntensity)` (legacy contract — keep "real signal only")
keep working without modification. Consumers that want the explicit
classification (e.g. per-lab sensitivity, presence/absence matrices)
read the `detection` column directly. Three classes:

- `"detected"` — picker passed all three guards (SNR ≥ 5, shape ≥ 2/5,
  absolute floor when noise=0). Trust `maxIntensity` / `areaUnderTic`.
- `"below_detection"` — the EIC was extracted and a tentative apex was
  located, but it failed at least one guard. `apex_rt` is populated so
  the noise location can be inspected; metric columns are NA. Counts
  toward cross-lab presence-absence as "not seen here".
- `"no_data"` — the EIC was empty, the prior had < 5 scans of context,
  or the prior itself was missing. Everything NA.

**`prior_kind` interpretation:** controls extraction confidence, not
identity confidence. A `prior_kind = "consensus"` row is not a weaker
identification; it just means we located the lab's MS1 signal using the
cross-lab consensus RT rather than its own SIRIUS hit (either because
SIRIUS missed it, or because the SIRIUS hit failed the COSMIC ≥ 0.1 /
|Δrt| ≤ 30 s sanity check that demotes spurious matches).

---

## `confidence_tier` lifecycle

- **P1** (this contract): producer writes the column with constant value
  `"unfiltered"`. P1–P4 consumers can run end-to-end without P5 existing.
- **P5** overwrites the column with `"A"` / `"B"` / `"C"`:
  - **A**: high COSMIC × MS1-RT-concordant across detecting labs
  - **B**: one axis strong, one axis marginal
  - **C**: low on both axes
- **Consumer rule:** identity-only analyses (P6 confidence panel, P5
  sensitivity re-run of P4) filter on `confidence_tier %in% c("A", "B")`
  or similar. Default for everything else: keep all rows.

## Gotcha: reading `mixture` as character

The `mixture` column uses dotted decimal labels (`"1.1"`, `"2.10"`, `"3.3"`).
Default `read.csv()` does type-inference on quoted values and coerces
to numeric, which silently collapses `"2.10" → 2.1` — colliding with
real mixture `"2.1"`. Every consumer that reads these CSVs **must**
force it:

```r
read.csv(path, stringsAsFactors = FALSE,
         colClasses = c(mixture = "character"))
```

(Or use `readr::read_csv()` with explicit `col_types`, which doesn't
have this bug.)

## Producer / consumer rules

- **Producers** (`sirius_annotation_all_ms2.qmd`) write the canonical
  CSVs and render no comparison plots. No identity-vs-RT, no cross-lab
  presence, no abundance heatmaps — those move to consumers.
- **Consumers** (`9_per_compound_diagnostics.qmd`, `8_data_level_progression.qmd`)
  start from `read.csv(...)` of these two files. They do not call
  `chromExtract`, run SIRIUS, or re-pick peaks.
- **`1_full_detected_objects.qmd`** is no longer responsible for the annotated
  EIC table — the producer's chrom_metrics covers it. Object creation
  keeps Full and Detected; the annotated *bulk* signal (P4 input) is
  derived from `annotation_chrom_metrics.csv` by summing area per
  sample.
- **Schema changes** require an entry below.

## Changelog

- **P10 (this commit) — Cluster-based consensus RT + proximity
  grouping for cembio.** Producer §7.1 now derives the consensus RT
  for each `(mixture, compound)` from a **cross-lab + cross-polarity +
  cross-feature clustering** instead of collapsing features with
  `median(rt_sec)`. Greedy 1-D clustering (`rt_cluster_window = 30 s`)
  groups the candidate RTs; clusters are scored as `n_labs × 10 +
  n_polarities × 5 + sum_cosmic`; the winning cluster's mean RT becomes
  the consensus prior for all labs in both polarities. Labs whose own
  MS2 fell in the winning cluster get `source = "own"`; labs whose
  MS2 fell in a different cluster (DDA triggered on the wrong peak)
  get `source = "consensus"`. Dry-run on the existing SIRIUS caches:
  77 % of (mixture, compound) tuples have a single cluster (clean);
  23 % have multiple clusters and the algorithm picks the supported
  one. 531 wrong-peak DDA features get demoted from "own" to
  "consensus".

  In parallel, §1's cembio acquisition handler replaces
  `rle(precursorMz)` with an RT-proximity grouping (`cembio_block_gap_s
  = 30 s`). The legacy `rle` was shattering one chromatographic peak
  into many features whenever the inclusion-list cycled through other
  precursors between two MS2 of the same compound — e.g. ritalinic
  acid mix 2.9 pos: 11 features at 198–258 s collapse to **1**
  feature under the new rule. **This requires re-running SIRIUS for
  cembio** (other labs are unaffected because they use
  `fragmentGroupIndex`).
- **P9 — All adducts + cross-polarity sanity.** Producer
  §7.1 emits one row per `(mixture, lab, polarity, compound,
  adduct)` instead of collapsing adducts (legacy `adduct =
  first(adduct)` dropped every adduct beyond the first). The
  cross-polarity demotion that was added here is now subsumed by
  P10's clustering algorithm (which uses cross-polarity as one of the
  cluster-scoring axes).
- **P6 — SNR-aware extraction.** Replaced `peakBoundary()`
  + 2-pass refinement with a single-pass SNR + shape + absolute-floor
  picker (producer §7.4). Added columns to `annotation_chrom_metrics.csv`:
  `snr`, `baseline`, `noise`, `n_scans_in_shape`, `prior_kind`,
  `detection`. Existing metric columns now go NA for `below_detection`
  rows (legacy `!is.na(maxIntensity)` filter still works as a "real
  signal only" guard). Motivation: the previous picker reported noise
  spikes at random RTs as "the peak" whenever the compound was below
  detection at a given lab, producing ~70 % spurious hmgu rows and
  corrupting cross-lab correlation. See `MANUSCRIPT_NOTES.md` § hmgu.
- *(initial)* — P1.1 schema freeze.
