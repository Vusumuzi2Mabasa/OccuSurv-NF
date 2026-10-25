# Executed validation

Validation date: 25 October 2026. All scientific inputs are synthetic.

| Check | Observed result |
|---|---|
| Unit suite | 14 tests passed; `unit_tests.txt` |
| Full Nextflow execution | Five processes completed with exit code 0; 500 generated workers |
| Full analysis | 499 eligible workers, 476 incident-cohort workers, 137 events, 2793.1064 person-years |
| Cox model | 461 complete cases, 130 events, three parameters, maximum absolute score 9.28e-14, no warnings |
| Low-event smoke execution | Five processes completed; 40 generated workers, 39 retained, 37 incident-cohort workers, 13 events |
| Model-skip behaviour | 12 complete-case events; Cox skipped with explicit reason |
| Resume | All five processes CACHED with exit code 0; `resume_trace.tsv` |
| Numerical checks | Explicit Cox partial likelihood, finite-difference score, independent Kaplan–Meier comparison, interval boundary cases |
| Input/output integrity | Raw SHA-256 hashes, cohort eligibility, counts, person-time, rates, model gate, figures and report verified |
| Word review | Final eight-page report exported through Microsoft Word to PDF, rasterised using bundled Poppler, every page visually inspected |
| HPC and Docker | Not tested; no remote jobs submitted and no container validation claimed |

The standard `render_docx.py` was attempted and failed because the runtime had no LibreOffice executable. Word COM export provided the rendering fallback; the scientific Word deliverable remains a normal DOCX file. QA PDFs and page images remain local and are excluded from Git.

## Version 0.2 extension validation — 25 October 2026

The earlier five-process and eight-page records above describe version 0.1. The extended workflow completed all seven processes successfully (`extension_trace.tsv`). It generated a separate follow-up scenario for 477 workers, with 2743 start-stop intervals, 93 diagnoses, 83 competing deaths and 2603.3552 person-years. Both cause-specific models converged without warnings; maximum absolute scores were below 1e-12.

The 40-worker smoke input produced 38 longitudinal subjects, 10 diagnoses and 5 deaths. Both regressions were explicitly skipped because their complete-case cause-specific event counts were below 30. All seven smoke processes completed (`extension_smoke_trace.tsv`). All seven full-run processes were then CACHED on resume (`extension_resume_trace.tsv`).

The expanded suite passed all 26 tests (`extended_tests.txt`), including independent cumulative-incidence comparisons, probability conservation, start-stop partial likelihood and finite-difference score checks, and malformed interval fixtures. `verify_extension.py` passed for full and smoke outputs, checking hashes, interval QC, counts, person-time, probability mass, monotonicity, confidence bounds, regression gates and deliverable existence. Four new figures and the final 12-page Word report were rendered; every report page was visually inspected using the Word COM/Poppler fallback described above.

This is one synthetic realisation and a small-data execution test, not a power, bias or coverage study. The 200-replicate worker bootstrap provides pointwise intervals under independent-worker assumptions. Neither Fine–Gray modelling nor remote HPC execution is claimed.

## Version 0.3 R and Quarto conversion — 25 October 2026

The production statistical implementation is now `bin/analyse.R`, using R 4.3.3, survival 3.5-8 and ggplot2 3.4.4. Quarto 1.10.19 executes `report/OccuSurv.qmd` and writes the Word report. Python simulation and QC remain unchanged; `validation/python-v0.2/` preserves the earlier numerical reference. R packages and session details are recorded in `R-packages.csv` and `results/analysis/R-sessionInfo.txt`.

The full six-process run completed with exit code zero. All 23 independent numerical comparisons passed. Baseline Cox coefficients, standard errors, hazard ratios, bounds and P values differed by at most 7.8e-15; cause-specific comparisons differed by at most 4.7e-13. Aalen-Johansen estimates agreed within 5.4e-15 and shared-resample bootstrap bounds within 5.3e-16. The exact observed differences and tolerances are in `results/crosscheck/comparison.csv`.

The 40-worker smoke run passed all 19 applicable comparisons. Baseline Cox and both cause-specific models were explicitly skipped below their event thresholds. The original 26 Python tests passed, and R checks passed for hand-computed AJ values, probability conservation, all-censored data and a start-stop likelihood with events exactly at interval boundaries. Full and smoke integrity checks passed (`r_integrity.txt`).

During development, the smoke comparison detected that the historical rounded normal critical value 1.96 produces larger KM interval differences in small samples. Rather than relax tolerance, the cross-check now independently recomputes Python bounds with the exact normal quantile used by R and compares them at 1e-12. The original Python output is retained. Bootstrap samples are recorded explicitly so R and Python recompute exactly the same resamples; the R intervals can differ slightly from the archived v0.2 intervals because the historical RNG stream was different.

Initial report execution preceded completion of the rmarkdown dependency installation and failed; it was rerun successfully after installation. The fs package required libuv development headers. Quarto was installed on the Linux filesystem after Toshiba extraction could not preserve executable metadata. No statistical failures were suppressed. Development and successful run logs are retained separately. The bundled document renderer lacked LibreOffice, so final visual verification used Microsoft Word export and bundled Poppler, as in earlier versions.

The final Quarto report contains ten pages, nine figures and native Word equations. Every page was visually inspected after the final layout revision; title borders, table formatting and a detached caption were corrected before delivery. The full six-process graph subsequently resumed entirely from cache (`r_resume_trace.tsv`); the final smoke report also completed successfully. The README was rewritten as a research report, with mathematical methods and interpretation before each figure. All results remain explicitly labelled as synthetic.

## Research manuscript reporting revision — 25 October 2026

The Quarto source was reorganised to match the README's research-paper structure: Abstract, Introduction, Methods with mathematical specifications, Results, Discussion and Conclusion, followed by data availability and references. Vusumuzi Mabasa is named in the title block and Word author metadata, as confirmed by the author. Interpretations precede all nine figures. The title, opening statement and results explicitly identify the data as simulated; no affiliation or institutional endorsement was added.

Report values and tables remain linked to the staged R outputs. The full and 40-worker workflows retained passing numerical comparisons, including explicit regression skips in the small cohort. Scientific analysis code and parameter values were unchanged. The document renderer again lacked LibreOffice, so Word COM export and Poppler were used for page-by-page visual verification. Execution logs and traces for this revision use the `manuscript_` prefix.

The final manuscript has 15 pages, nine figures, 43 native Word equations and 9 tables. All 15 pages were visually inspected after the final render. The full run passed 23 comparisons and the smoke run passed 19; both final report tasks completed with exit code zero. The delivered Word file is byte-identical to the final full-run output.

## Runtime details for the historical Python report

Nextflow 25.10.2 and OpenJDK 21.0.12.1 executed under Ubuntu 24.04 in WSL. A project-specific Python 3.13 environment at `/home/vusi/.venvs/occusurv-nf` contains the scientific dependencies listed in `python-lock.txt`. The complete workflow source, task work and outputs are on Toshiba at `D:/occusurv-nf` (`/mnt/toshiba/occusurv-nf` in the tested WSL session).

The report process ran the bundled Windows Python 3.12.14 with python-docx 1.2.0 through WSL interoperability and copy staging. The local config is excluded from Git; the exact settings are documented in `../RUN_WORKFLOW.txt`. The default portable configuration uses Python from PATH. This default all-Linux report-runtime combination has not been separately executed; the report script itself was tested as a Nextflow process with the bundled runtime.

Initial integration caught a floating-point boundary problem in employment duration after CSV rounding; the QC now allows a 1e-8-year tolerance at that boundary. Initial Nextflow config interpolation and smoke report paths were corrected. These were development failures, not ignored failures in the final run.

The full trace, smoke trace, cache trace, integration checks and numerical test output are included. `full_trace.tsv` records the first successful execution; subsequent report-only changes improved layout without changing the statistical analysis. All raw and clean data, statistical tables and plotting source are available for inspection.

Text sources were normalised to LF for portable Git checkout. The source manifest records these canonical bytes. Generated Nextflow HTML reports are retained locally but excluded from Git; TSV execution traces are committed.
