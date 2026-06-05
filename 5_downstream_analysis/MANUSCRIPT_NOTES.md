# Manuscript notes — current state after P6 extraction overhaul

Working draft. Most recent update: **the entire ring-trial story** was
rewritten after auditing the chromExtract producer and finding two
extraction bugs that were corrupting headline numbers. See the §P6
findings section below for the corrected picture. The pre-overhaul
sections (§ "Findings (numbers from §9 ICC analysis)" and below) are
preserved for diff-history but their numbers are stale.

## P6 extraction overhaul — what changed and what we now know

### Extraction bugs identified and fixed (Jun 2026)

Two distinct failure modes in the legacy chromExtract pipeline
(`peakBoundary()` + two-pass refinement) were producing systematically
wrong abundance numbers. Both are fixed in the producer (§7.4-§7.6).

**Bug 1 — noise grabs at sub-floor compounds.** For compounds below
a given lab's per-compound detection floor, the EIC is near-empty;
`peakBoundary(threshold = 0.1)` returns the tallest local maximum —
i.e. a noise spike at a random RT — and the pipeline reports its
intensity (~10-100 counts) as the "measurement". Impact: ~70 % of
hmgu's annotated rows were noise grabs. Mechanism: replaced with
SNR + shape + absolute-floor guards in §7.6. Verified on a 6-compound
test panel and on the actual production output.

**Bug 2 — noise estimate contaminated by other peaks at the same m/z.**
The noise estimator used the full EIC tail (|rt − apex| > 25 s) which
frequently contained other chromatographic peaks at the same 20 ppm
m/z window (isomers, late-eluting contaminants). Their intensity was
interpreted as "noise" and crushed SNR for real peaks. Cembio
L-phenylalanine pos: noise = 894 k, SNR = 0.01 → flagged
below_detection despite COSMIC = 0.992. **Fix**: noise estimated from
the **bottom 50 % of the tail** (intensities below the tail median),
which is the true baseline. Verified by recovering 77 high-COSMIC
own/own false negatives.

**Bug 3 — same `peakBoundary` failure for NAPS extraction.** The
NAPS producer `3_naps_extraction.qmd` used the same legacy
peakBoundary pattern with a ±10 s search. icl pos NAPS C10 reported
median 386 counts (noise) when the real peak was at 502 k–538 k
counts, just 11 s before the cross-lab consensus RT. The legacy
two-pass retry never fired because peakBoundary returned a noise spike
inside the ±10 s window. Fix: same SNR-aware picker with a ±20 s
apex search and ±60 s wide EIC.

### Producer outputs after both re-renders

`annotation_chrom_metrics.csv` schema gains columns: `snr`,
`baseline`, `noise`, `n_scans_in_shape`, `prior_kind`, `detection`.
Metric columns are NA for `below_detection` rows so legacy consumer
filters `!is.na(maxIntensity)` keep working.

**Detection rates per (lab, polarity) — first sanity number:**

|  | neg detected | neg below_det | pos detected | pos below_det |
|---|---|---|---|---|
| afekta | 196 / 231 (84.8 %) | 32 | 203 / 242 (83.9 %) | 29 |
| cembio | 192 / 231 (83.1 %) | 31 | 176 / 242 (72.7 %) | 38 |
| hmgu   | 134 / 231 (58.0 %) | 25 | 109 / 242 (45.0 %) | 43 |
| icl    | 213 / 231 (92.2 %) | 17 | 194 / 242 (80.2 %) | 46 |

Three of four labs detect ~75-92 % of the annotated panel. **hmgu
sees only ~45-58 %** — about half the panel falls under hmgu's
per-compound MS1 floor. Sciex Zeno TOF 7600's MS1 sensitivity is
not low globally (its NAPS responses are highest of the four labs)
but per-compound ionization efficiency differs from Agilent Jet
Stream and Waters Z-spray, and the spike-in panel was tuned for
Agilent sensitivity.

NAPS chrom_metrics: icl pos within-lab CV dropped from 57 % to 7 %.
Variance decomposition: icl pos went from 14 % anchor / 34 %
replicate / 52 % residual (was "instrument chaos") to 76 % / 3 %
/ 21 % (clean anchor-explained signal).

### Headline cross-lab ICC ladder (clean data)

ICC computed as `var_between_compound / (var_between + var_within)`
on log10(maxIntensity + 1) via `lme4::lmer(~ 1 + (1|compound))`,
restricted to compounds detected in all labs of the scope.

| Scope                         | n pos | n neg | Raw pos | Raw neg | Global pos | Global neg |
|-------------------------------|------:|------:|--------:|--------:|-----------:|-----------:|
| Agilent pair (afekta+cembio)  |   166 |   175 |   0.396 |   0.589 |  **0.504** |  **0.563** |
| afekta + icl (cross-vendor)   |   169 |   187 |   0.170 |   0.248 |  **0.520** |  **0.481** |
| 3-lab (afekta + cembio + icl) |   143 |   171 |   0.220 |   0.333 |  **0.429** |  **0.440** |
| 4-lab (all)                   |    82 |   105 |   0.089 |   0.170 |   0.287    |   0.274    |

**Cross-vendor (afekta + icl, Agilent ↔ Waters) reaches ICC ≈ 0.50**
under Global normalization — nearly matching the same-platform
Agilent pair ceiling (0.50-0.56). This is the strongest single
finding of the paper.

**3-lab (no hmgu) reaches ICC ≈ 0.43-0.44** — defensible
cross-platform reproducibility.

The 4-lab number (0.27-0.29) is capped by hmgu's detection rate, not
by reproducibility per se: for compounds hmgu detects, its pairwise
Spearman ρ with the other three matches icl's (ρ = 0.15-0.36).

### Cross-lab Spearman pairs (Raw, since Spearman ρ is invariant to monotonic scaling)

|  | pos ρ (n) | neg ρ (n) |
|---|---|---|
| afekta ↔ cembio | 0.479 (167) | 0.619 (177) |
| afekta ↔ icl    | 0.563 (170) | 0.512 (189) |
| afekta ↔ hmgu   | 0.271 (102) | 0.359 (115) |
| cembio ↔ icl    | 0.331 (152) | 0.351 (184) |
| cembio ↔ hmgu   | 0.157 (092) | 0.225 (116) |
| hmgu ↔ icl      | 0.344 (099) | 0.150 (127) |

The afekta ↔ icl cross-vendor ρ of 0.51-0.56 sits between the two
same-Agilent pair (0.48-0.62), reinforcing that cross-vendor LC-MS1
reproducibility is achievable with harmonized chromatography.

### NAPS-anchored intensity normalization — does NOT outperform Global

We built a NAPS-anchored normalization: per (lab, polarity), fit a
LOESS curve through `log2(cross-lab median anchor intensity / lab
anchor intensity)` as a function of anchor RT. For any annotated peak
at (lab, polarity, apex_rt), multiply its intensity by the LOESS
prediction. Compared to Global (lab-median scaling) on the all-4-lab
detected subset:

| Strategy           | 4-lab pos ICC | 4-lab neg ICC | 3-lab pos | 3-lab neg |
|--------------------|--------------:|--------------:|----------:|----------:|
| Raw                | 0.089         | 0.170         | 0.220     | 0.333     |
| **Global**         | **0.287**     | **0.274**     | **0.429** | **0.440** |
| NAPS               | 0.070         | 0.100         | 0.285     | 0.194     |
| NAPS + Global      | 0.252         | 0.202         | (similar) | (similar) |

NAPS norm underperforms simple Global at every scope and polarity.
**Mechanism**: NAPS standards are designed to span the chromatographic
range (RT calibration), but their MS1 response profile differs from
the annotated panel chemistry. icl Waters neg: NAPS standards
(phenolic acids) are weak in early RT → NAPS-derived factor for
icl neg early RT = log2 5.86 = 58 ×. Multiplying icl's early-RT
amino acids by 58 over-corrects — they don't share Waters Z-spray
ionization efficiency with the NAPS phenolic acids.

**Paper framing**: NAPS-anchored normalization is theoretically
attractive (uses external reference, captures per-RT lab-specific
response) but **empirically does not transfer across compound classes**.
Simple lab-median (Global) scaling is the recommended primary
normalization for this dataset.

### Within-lab RT alignment from preprocessing — not used in current pipeline

All 4 labs run `adjustRtime(mse, PeakGroupsParam(rt_matrix))` in their
preprocessing. The aligned RTs are stored in the saved mse object,
NOT in the raw mzML. Both producers (annotation + NAPS) read raw
mzML via `Spectra(file_path, MsBackendMzR)` and therefore work with
**un-aligned** RTs. Effect on the analysis: small. Within-lab drift
is typically <2 s; cross-lab drift is much larger (10-20 s) and is
what the §7 RT alignment in `4_rt_alignment.qmd` corrects via
per-lab Hyman monotonic splines. A future cleanup could route
chromExtract through the aligned mse Spectra object instead of raw
mzML, but is not on the critical path for current findings.

### Hardware setup (verified `lab_setup_phase1.csv`)

| Lab | LC | MS instrument |
|---|---|---|
| afekta | Agilent 1290 Infinity II UHPLC | Agilent 6546 Q-ToF |
| cembio | Agilent 1290 Infinity II UHPLC | Agilent 6545 Q-ToF |
| icl    | Waters Acquity UHPLC            | Waters Synapt G2S Q-ToF |
| hmgu   | Sciex ExionLC AD                | Sciex Zeno TOF 7600 |

All four use ESI sources; harmonized LC gradient, mobile phases,
column. Lab identity remains confounded with MS platform (n = 1
instrument per lab).

---

## P10–P13 updates (Jun 2026, post-cluster-algorithm)

### Cluster-based consensus RT and multi-adduct retention (P10)

The producer §7.1 was rewritten so that each candidate annotation
contributes one row per FEATURE (not collapsed to dominant adduct),
and per-(mixture, compound) RT priors come from a greedy 1-D
agglomeration over labs (30 s window, score = n_labs × 10 +
n_polarities × 5 + Σ COSMIC). The winning cluster's median becomes the
RT prior used for chromExtract, and `apex_pad` widens with the
cluster's observed `rt_spread`.

**Effect on detection rates (vs pre-cluster):**

| Lab    | neg (was) | neg (now) | pos (was) | pos (now) |
|--------|----------:|----------:|----------:|----------:|
| afekta |    84.0 % |  **89.2** |    64.0 % |  **87.8** |
| cembio |    83.0 % |  **86.6** |    73.0 % |  **76.2** |
| hmgu   |    58.0 % |    58.4   |    45.0 % |    49.3   |
| icl    |    92.0 % |  **91.8** |    80.0 % |  **85.8** |

**Effect on cross-polarity RT agreement** (the biggest manuscript
talking point): zero compounds with ≥30 s pos/neg disagreement (was
~15 %); maximum disagreement now ≤ 30 s vs ~100 s+ outliers before.
**Effect on MS1 rescue**: 1205 consensus-prior detected rows (was 636),
mostly recovered through multi-lab MS2 agreement at compounds the
local lab missed.

The producer §7.1 also includes a RT-proximity grouping replacement
for the legacy `rle(precursorMz)` on cembio (block-MS2 — see below).
Effect not yet realized because the cembio SIRIUS caches were not
deleted, but is in place for any future cembio re-run.

### cembio block-MS2 explainer (P11) → supplementary figure

`9_per_compound_diagnostics.qmd` includes a reproducible figure
showing why cembio's targeted-inclusion-list MS2 strategy (MS1 and
MS2 in temporally-disjoint files) cannot be apex-anchored from MS2
trigger times alone — ritalinic acid mix 2.9 as the worked example.
The figure pulls live from `D_1_pos_021_*.mzML` and
`D_1_pos_MS2_021_*.mzML` so it stays current with re-renders.
**Recommended for the supplementary figure slot** explaining the
"acquisition strategy" axis of the ring trial.

### hmgu DDA explainer (P12) → supplementary figure / discussion

`9_per_compound_diagnostics.qmd` now also includes a parallel section
on hmgu's open-DDA acquisition. For mix 2.3 pos:

- **hmgu fires 8 147 MS2** across **5 238 unique precursor m/z**
  (DDA across the whole m/z range, single CE = 35 eV)
- **Only 0.8 %** of those MS2 fire within ±50 ppm of an expected
  spike-in compound m/z. The other 99.2 % target matrix.
- By contrast, cembio's 23-m/z inclusion list spends 100 % of its
  duty cycle on known targets; afekta uses 202 precursors with
  stepped CE (10/20/40); icl uses 176 precursors at low CE.

**Implication for framing**: hmgu's pairwise ICC gap (0.12-0.36 vs
afekta+cembio 0.87, afekta+icl 0.81, cembio+icl 0.79) is dominated
by **acquisition policy mismatch**, not Sciex Zeno-TOF sensitivity per
se. Open-bore DDA is the wrong instrument-configuration mismatch for a
targeted spike-in panel; it would likely excel at untargeted biomarker
discovery.

**Additional structural data deficit**: 162 of hmgu's 167 `no_data`
compound tuples have a consensus RT from other labs, confirming the
compounds *are* present in the chromatographic run — they're just
below hmgu's per-compound MS1 detection floor. For 229 compounds
detected in all 3 non-hmgu labs AND in hmgu, hmgu's median
`maxIntensity` is **8.8 %** of the median of the other three
(i.e. ~12× lower). 52 % of hmgu detections are < 10 % of the
other-three intensity.

### Headline ICC numbers — current

3-lab Strict ICC (afekta + cembio + icl) under best strategy:

| metric        | polarity | Lab-median ICC | Aligned ICC | Decile ICC |
|---------------|----------|---------------:|------------:|-----------:|
| areaUnderTic  | Negative |      **0.658** |       0.509 |      0.596 |
| areaUnderTic  | Positive |      **0.596** |       0.577 |      0.537 |
| maxIntensity  | Negative |      **0.569** |       0.385 |      0.412 |
| maxIntensity  | Positive |      **0.467** |       0.294 |      0.297 |

Pairwise Spearman on `areaUnderTic` (Strict, all 6 lab pairs):

| pair                | neg ρ | pos ρ |
|---------------------|------:|------:|
| afekta + cembio     | 0.868 | 0.796 |
| afekta + icl        | 0.815 | 0.619 |
| cembio + icl        | 0.792 | 0.634 |
| afekta + hmgu       | 0.287 | 0.279 |
| cembio + hmgu       | 0.119 | 0.193 |
| hmgu   + icl        | 0.357 | 0.482 |

**hmgu framing decision**: keep all four labs in the main figure;
report the 4-lab ICC table with hmgu as a labeled outlier in the
pairwise heatmap; quote **3-lab ICC 0.66** (areaUnderTic neg
Lab-median) as the headline reproducibility number when the
acquisition strategy is uniform (scheduled or inclusion-list MS2,
stepped CE); attribute hmgu's pairwise gap to the DDA section in
`9_per_compound_diagnostics.qmd`.

### RT-alignment threshold fix (P13.1)

`standards_anchor_max_aligned_sd` raised 3 → 7 s in
`4_rt_alignment.qmd`. At 3 s only 2 standards qualified as splines
anchors, leaving the per-lab Hyman spline floating freely between
NAPS anchors → median cross-lab aligned SD ended up 5.86 → 7.80 s
**WORSE** than raw. With threshold 7 s, 45 standards qualify and
the spline tightens:

| threshold | n_std | median aligned SD | p95 aligned SD | pct worse |
|----------:|------:|------------------:|---------------:|----------:|
|       3 s |     2 |            7.10 s |        20.63 s |    60.8 % |
|   **7 s** |**45** |       **5.26 s** |   **18.62 s** |  **46.1 %** |
|      10 s |   109 |           10.03 s |        37.54 s |    60.8 % |

Beyond 10 s, low-quality standards dominate and degrade the fit.
Rank-inverted NAPS anchor detection (icl_neg C2/C3) was a separate
investigation; the PAV step already absorbs them so dropping them
explicitly changes median SD by < 0.1 s.

### Scope clarification — `2_consensus_peaks` re-render

`2_consensus_peaks.qmd` does NOT need re-rendering after the cluster
algorithm change. Its inputs (`detected_peaks_*_HE.csv`) are xcms
peak-picker outputs from Feb 6 and are independent of the SIRIUS MS2
grouping. Phase 5 of `8_data_level_progression.qmd` correctly
reflects xcms-detected cross-lab matched peaks; this is a different
question from the Annotated row's SIRIUS-confirmed spike-in panel.
Phase 5 noted explicitly in the QMD as a scope reminder.

---

## Legacy notes (pre-P6 numbers; superseded by the above)

## Title under consideration

> Defining and overcoming the limits of multi-center metabolomics:
> a computational blueprint for harmonized-LC ring trials

(Earlier framing emphasized "platform identity dominates". Current
findings show that's only partly true — see Phase 3.)

## Core question

What are the computational steps required to make multi-center
metabolomics data integration possible under harmonized LC, and what
are the irreducible limits of cross-platform measurement of annotated
abundance?

## Hardware setup (verified from `lab_setup_phase1.csv`, 2024-12)

| Lab | LC | MS instrument |
|---|---|---|
| afekta | Agilent 1290 Infinity II UHPLC | **Agilent 6546 Q-ToF** |
| cembio | Agilent 1290 Infinity II UHPLC | **Agilent 6545 Q-ToF** |
| icl    | Waters Acquity UHPLC            | **Waters Synapt G2S Q-ToF** |
| hmgu   | Sciex ExionLC AD                | **Sciex Zeno TOF 7600** |

All four use ESI sources; harmonized LC gradient, mobile phases, and
column. Lab identity remains confounded with instrument platform —
n=1 instrument per lab.

## Findings (numbers from §9 ICC analysis)

### Identity (P1, working)

- MS2 yield varies by inclusion-list strategy: cembio (full-MS2) 92%
  vs icl (targeted inclusion list) 2% of expected compounds annotated.
- 275-set of well-detected compounds = compounds with raw maxIntensity
  > 0 in all 4 labs at the matching theoretical m/z. This is the
  annotated subset used throughout downstream analysis.

### Retention (P2, working)

- Per-`(lab, polarity)` Hyman monotonic spline mapping local apex_rt
  → consensus apex_rt, built from 20 APS anchors + spiked standards
  where ≥3 labs co-detect.
- Cross-lab apex-RT residual: median 7-8 s after alignment
  (down from ~10 s raw); p95 ≈ 50 s after clamping to anchor range.
- Within-lab NAPS replicate SD ≈ 0.2 s for afekta+cembio+icl; cembio-pos
  had wrong-peak outliers fixed by chromExtract producer; **icl
  early-RT (< 60 s) is genuinely unreliable** — that's the chromatography,
  not the alignment.

### Abundance (P3 — most of the new findings)

**Cross-platform Annotated bulk ICC under best strategy (LOESS, Strict,
n=253 well-detected and aligned compounds):**

|  | maxIntensity | areaUnderTic |
|---|---|---|
| Negative | 0.252 | 0.223 |
| Positive | 0.177 | 0.149 |

**Same-platform (afekta + cembio Agilent pair) ICC under best strategy:**

|  | maxIntensity | areaUnderTic |
|---|---|---|
| Negative | 0.688 | 0.608 |
| Positive | 0.807 | 0.769 |

**3-lab ICC (excluding hmgu — afekta + cembio + icl) under best strategy:**

|  | maxIntensity | areaUnderTic |
|---|---|---|
| Negative | 0.667 | 0.667 |
| Positive | 0.512 | 0.561 |

**Implication:** the 4-lab cross-platform Annotated ICC of ~0.2 is
largely driven by hmgu's outlier behavior. The 3-lab subset
(afekta + cembio + icl, spanning Agilent Q-TOF and Waters Synapt
platforms) reaches **ICC ≈ 0.5–0.67 under best normalization**.
That's a defensible cross-platform reproducibility number.

The Agilent same-platform pair tops out at ~0.7–0.9; the gap to
3-lab cross-platform is the platform-identity cost among these 3
instruments (~0.1–0.3 ICC).

### hmgu (Sciex Zeno TOF 7600) — root cause identified

hmgu's per-(mixture, compound) maxIntensity has near-zero bulk
correlation with the other 3 labs (Spearman ρ ≈ 0.03–0.05). The
preceding draft of this section hypothesized DDA-interleave and/or
generic Sciex sensitivity. Raw-file inspection (mzML metadata + per-lab
NAPS) rules both out and identifies the real mechanism.

**Raw-file diagnostics on a paired (Human Endosome 2.1, pos) injection:**

| Property | afekta | cembio | hmgu | icl |
|---|---|---|---|---|
| MS1 scans / file | 1499 | 1326 | 1107 | 3180 |
| MS1 cycle gap median (s) | 0.60 | 0.68 | **0.67** | 0.23 |
| MS1 cycle gap max (s) | 0.60 | 0.68 | 0.81 | 0.25 |
| Gaps > 2s | 0 | 0 | **0** | 0 |
| Centroided | TRUE | TRUE | TRUE | TRUE |

hmgu's interleaved Top-12 DDA does **not** starve MS1 — the Zeno cycles
fast enough that MS1 spacing is regular (median 0.67 s, no gaps > 2 s).

**NAPS QC standards (same RT-ladder spike-in at every lab) intensity:**

| Lab | NAPS pos median | NAPS pos geo-mean |
|---|---|---|
| afekta | 315,508 | 185,226 |
| cembio | 178,499 | 75,065 |
| **hmgu** | **540,889** | **354,473** |
| icl | 112,725 | 30,769 |

**hmgu's NAPS responses are the HIGHEST of all four labs.** The
instrument is fully sensitive in MS1. So the floor is not the
instrument — it is **per-compound** for the annotated spike-in
set.

**Annotated-compound chromExtract diagnostics (`apex_rt` vs cross-lab
consensus, stratified by hmgu maxIntensity tertile):**

| hmgu intensity tertile | geo-mean intensity | median \|RT dev\| (s) | p95 \|RT dev\| (s) |
|---|---|---|---|
| low | 11 | 17–20 | up to 184 |
| mid | 122–534 | 11 | 67–91 |
| high | 10,953–113,334 | **8.5–9.2** | 29–41 |

And restricting hmgu's correlation by RT-placement quality:

- Well-placed hmgu peaks (\|rt_dev\| ≤ 5 s): ρ_pos = **0.34**
  vs other-lab median (n = 65). Defensible cross-lab signal.
- Poorly-placed hmgu peaks: ρ ≈ −0.05 to 0.02 (n = 144–153). Pure noise.

**Mechanism (concrete):**

The chromExtract producer (`sirius_annotation_all_ms2.qmd` §7) computes
maxIntensity by `Chromatograms::peakBoundary()` — a local-maximum search
on each EIC. When the compound is **below hmgu's per-compound detection
floor at this dilution**, the EIC contains only noise; `peakBoundary()`
returns the tallest noise spike, and that noise apex_rt + tiny
maxIntensity gets recorded as a "measurement." The result is that
~70 % of hmgu's annotated rows are noise grabs scattered across the
chromatogram, not real peaks. They dominate the bulk ρ and crush
cross-lab correlation.

The instrument is fine. The pipeline is doing exactly what its EIC-based
design specifies. The mismatch is that **the spike-in concentrations
are tuned for the Agilent labs' sensitivity profile** and many specific
compounds fall below hmgu's per-compound floor (likely a Sciex ESI
vs Agilent Jet Stream ionization-efficiency difference). The wiff→mzML
conversion (pwiz_Reader_ABI 3.0.25071) is standard and not the issue;
m/z calibration is fine (median \|dev\| 0.85 ppm neg, 4.0 ppm pos).

**Proposed honest-reporting fix (downstream of producer):**

In the consumer (7_icc_reproducibility.qmd, or as a post-hoc step), classify
a row as "below detection" rather than "low signal" when both:

- \|apex_rt − cross-lab median apex_rt\| > 5 s, AND
- `signalToNoiseRatio < SNR_floor` (e.g. < 3)

…and set its maxIntensity to NA. This stops the noise grabs from
contaminating ρ and lets hmgu's real measurements show their true
correlation with the others (ρ ≈ 0.3–0.45 by polarity).

**Implication for paper framing:**

"hmgu has lower cross-lab reproducibility" is replaced by **"hmgu has
strong MS1 sensitivity (highest NAPS responses) but ~70 % of the
spike-in compound set falls below its per-compound ionization-efficiency
threshold at the dilution used; for the 30 % above threshold,
cross-lab correlation is ρ ≈ 0.3–0.45, comparable to icl."**

We **do not** drop hmgu. The 4-lab number stays the headline. A
detection-floor-aware re-analysis is the methodologically honest
companion.

## Strategy ranking (P3)

In order of best paper-default given current data:

1. **Lab-median** (one factor per `(lab, polarity)`, derived from the
   bulk distribution; matches `8_data_level_progression.qmd`'s `tic_norm`
   formulation). Most robust across lab pairs — never hurts
   substantially, usually helps. Recommended primary.
2. **Decile / LOESS** (per-`(lab, injection, RT decile)` or smooth).
   Helps a small additional amount when there's RT-localized
   non-linearity. Can over-fit and reduce same-platform ICC.
3. **Global** (per-`(lab, injection)`). Less stable than Lab-median
   because each injection's factor is estimated on smaller data.

## "Platform identity dominates" framing — needs correction

The earlier reading was: "same-platform Agilent → 0.84+; cross-platform
4-lab → 0.2; the 0.6 gap is the platform cost." Current data shows
this isn't quite right:

- afekta + icl (Agilent + Waters) pairwise ICC = 0.67–0.87 in negative
  polarity. Cross-vendor reproducibility **can** be high.
- afekta + cembio (Agilent same-platform) ICC = 0.61–0.94 — only
  modestly higher than the best cross-platform pair in neg.
- The 4-lab number is low primarily because of hmgu, not because of
  platform-identity per se.

Better framing for the paper: **"harmonized LC + careful processing
delivers cross-platform reproducibility for most platforms; one of
four labs in this study showed outlier behavior whose specific source
we cannot identify from the available data."**

## Phase 1 (Baseline problem) — original wording stands

- Lab variance dominates raw data (~97–100% of TIC variance is
  lab-attributable in Full Signal). Variance decomposition is in
  `8_data_level_progression.qmd` §1.2.
- PCA / UMAP / presence-absence visualizations in
  `8_data_level_progression.qmd` §1.3.
- Consensus peak matrices are sparse and block-structured — each
  consensus peak appears in only one mixture across labs.

## Phase 2 (RT alignment) — see `4_rt_alignment.qmd` for detail

Per-lab Hyman monotonic spline against cross-lab consensus, NAPS-based
anchors + standards extension. See file for math + diagnostics.

## Phase 3 (Normalization + ICC ladder) — see `8` + `9`

`5_normalization.qmd` builds the factors and writes the
normalized annotated table. `7_icc_reproducibility.qmd` reads it and computes
the ICC ladders + pairwise comparisons + per-lab diagnostics + per-compound CV.

## Limitations

- n = 4 labs, no technical replicates per lab. ICC point estimates are
  noisy — bootstrap or LMM would be appropriate but not currently in
  the pipeline.
- Lab identity fully confounded with MS platform (n=1 instrument per
  lab).
- Annotated set (~275 well-detected compounds) is spiked standards —
  may not represent the full feature space. The Full/Detected backdrop
  from `2_consensus_peaks.qmd` is the broader
  context.
- Same-platform pair (afekta + cembio, both Agilent) is two different
  Agilent Q-TOF models (6546 vs 6545). The "same-platform ceiling" is
  an Agilent-line ceiling, not strictly identical-instrument.

## Open items / decisions for the paper

- **Should the ICC table report Lab-median + LOESS, or just Lab-median?**
  LOESS helps a small additional amount but introduces a more involved
  methodology. Single-strategy reporting (Lab-median) is cleaner;
  multi-strategy is more transparent.
- **3-lab analysis as primary or sensitivity?** Reporting 4-lab as
  primary with a 3-lab sensitivity panel is the standard ring-trial
  approach.
- **Confidence tiering (P5) for compounds with ambiguous
  identification.** The `confidence_tier` column is currently
  `"unfiltered"` placeholder; P5 should populate from COSMIC scores ×
  cross-lab apex_rt concordance.

## Notes for future me

- The `Consensus Signal` data level mentioned in older notes (consensus
  peaks across labs from `2_consensus_peaks.qmd`)
  is preserved as PCA panel in `8_data_level_progression.qmd` §1.3.2 but not
  in the Annotated bulk ICC ladder. Could add if helpful.
- The strategy color palette is centralized in
  `8_data_level_progression.qmd` setup and `7_icc_reproducibility.qmd` setup. Keep
  them in sync.
- hmgu instrument: Sciex Zeno TOF 7600 (NOT Bruker — earlier draft of
  this file had the wrong identification).
