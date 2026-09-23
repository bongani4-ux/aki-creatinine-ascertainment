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
