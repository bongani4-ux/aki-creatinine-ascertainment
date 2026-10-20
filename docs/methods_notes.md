# Methods notes

## First episode per patient

Repeat episodes create overlapping windows: an earlier admission can supply the "baseline"
for a later one. The analysis keeps one episode per patient, the earliest in time.

`first_episode_per_patient()` in `analysis/pipeline.R` orders episodes by t0 where t0 carries real time
(a datetime, or an offset that varies between a patient's episodes). Where every episode
has t0 = 0, time order cannot be recovered from the export and `episode_id` order is used.
The rule applied is recorded per source in `results/measurement_and_provenance.csv`
(`first_episode_rule`). In this dataset that fallback applied to eICU-CRD and Zigong.

An earlier version of the analysis ordered all sources by `episode_id`. Because episode IDs
are not chronological in several sources, it sometimes selected a later admission, whose
earlier admission then sat inside the baseline window. This inflated hierarchy-baseline
availability in MIMIC-IV (44.5% vs 40.3%), NWICU, INSPIRE, SICdb and MOVER. The published
figures use time ordering. No conclusion changed.

## Data quality notes

- **NWICU** `died_inhosp` is mis-mapped in the export (flags 118 patients, none of whom died
  before discharge). Death timestamps are correct (12.2% in-hospital mortality). The analysis
  uses timestamps throughout.
- **ZD** contains 2,350 exact duplicate result rows (1.1%). They do not affect evaluability.
- **PIC** An initial export returned a filtered subset of 7,488 ICU stays rather than 13,941,
  inflating apparent ascertainment by 30 percentage points. The full cohort is used. See
  Supplementary Item S4 of the paper.

## Definitions as implemented

- **Windows** (lead = days before t0). Hierarchy baseline 7 <= lead <= 365 (tier 1 is
  [7, 30]); proximal tier 0 is 0 < lead < 7. The same boundaries are used for provenance
  (Table 3), the cascade (Table 5) and misclassification (Table 4). An earlier version put a
  value exactly 7 days before t0 in tier 0 for Table 3 but in the hierarchy elsewhere.
- **KDIGO staging** compares against thresholds with a 1e-9 tolerance: in floating point
  1.2 - 0.9 = 0.29999999999999993, so a rise from 0.9 to 1.2 mg/dL was not counted.
  The 0.3 mg/dL criterion uses the peak within 48 h after t0 against the chosen baseline
  (a fixed-baseline approximation of KDIGO's 48 h window). If there is no value within 48 h
  the criterion is not assessed; it is no longer replaced by the 7-day peak. Stage 3 by an
  absolute value >= 4.0 mg/dL requires a stage-1 criterion to be met first, so a stable
  patient with a baseline above 4.0 (e.g. on dialysis) is not classified as AKI.
- **Duplicates.** Rows with the same episode, t_days and value are one draw recorded more
  than once and are removed before counting (`DEDUPLICATE` in `config.R`). In the exports
  this removed 60% of HiRID rows, 8.5% of INSPIRE, 7% of SICdb and PIC, and <0.1% of MIMIC-IV.
- **Undefined vs unmeasured** (Tables S2-S3). Patients with no pre-index value cannot be
  evaluated whatever happens later, so they are a separate category; the death / discharge
  split is made among patients with a pre-index value, and the "% of measurable" denominator
  is patients with a pre-index value alive at day 7 (AKI) or day 90 (AKD). The tables cover
  every source that records death timing (previously: sources with >5% deaths, which dropped
  INSPIRE and MOVER).
- **Medians** of results per patient are reported as computed (e.g. 2.5), not truncated.

## Language and validation

The analysis is implemented in R (`analysis/*.R`). An earlier implementation in Python
was used to develop and verify the method; the R version was checked against it source by
source. All reported quantities agreed to machine precision, except for:

1. **Numeric time offsets.** eICU-CRD records discharge and death as minutes from t0, and
   SICdb and VitalDB as seconds. The earlier implementation read these integers as
   datetimes, collapsing every death to day 0. The R pipeline converts them using the
   per-source `offset_unit` in `config.R`. This changed only the undefined-versus-unmeasured
   decomposition (Supplementary Tables S2-S3): for SICdb the death-attributable share of
   AKD non-evaluability is 12.0% (previously 17.9%), because most SICdb deaths occur after
   day 90. eICU-CRD totals were unchanged; only the split between early and late deaths moved.
2. **One MIMIC-IV patient** with two creatinine results at an identical timestamp, where the
   two languages break the tie differently (0.005 percentage points; invisible at reported
   precision).

## A note on date parsing in R

`as.POSIXct(x, tryFormats = ...)` chooses a format only if it fits every non-missing value.
A single empty string makes it return NA for the whole vector, silently. `parse_time()` in
`pipeline.R` therefore parses element by element. This matters for any export in which some
discharge or death times are blank.

## Corrections after the first submission version (October 2026)

The analysis was ported to Python and the port, run with the original logic, reproduced the
committed `results/` exactly (every column, including the timing histograms) for MIMIC-IV,
NWICU, HiRID, Zigong, PIC, INSPIRE, VitalDB and SICdb. The fixes above were then applied one
at a time to measure their effect. `results/` has not yet been regenerated: re-run
`run_all.R` on all fourteen sources. Three extraction errors were also found:

1. **Zigong times are hours, not days** (see `extraction/README.md`).
2. **INSPIRE t0** was midnight of the operation date; it is now the operation start.
3. **HiRID** rows were repeated 2-5 times per draw.
