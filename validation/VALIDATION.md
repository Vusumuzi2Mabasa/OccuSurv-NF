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

## Runtime details

Nextflow 25.10.2 and OpenJDK 21.0.12.1 executed under Ubuntu 24.04 in WSL. A project-specific Python 3.13 environment at `/home/vusi/.venvs/occusurv-nf` contains the scientific dependencies listed in `python-lock.txt`. The complete workflow source, task work and outputs are on Toshiba at `D:/occusurv-nf` (`/mnt/toshiba/occusurv-nf` in the tested WSL session).

The report process ran the bundled Windows Python 3.12.14 with python-docx 1.2.0 through WSL interoperability and copy staging. The local config is excluded from Git; the exact settings are documented in `../RUN_WORKFLOW.txt`. The default portable configuration uses Python from PATH. This default all-Linux report-runtime combination has not been separately executed; the report script itself was tested as a Nextflow process with the bundled runtime.

Initial integration caught a floating-point boundary problem in employment duration after CSV rounding; the QC now allows a 1e-8-year tolerance at that boundary. Initial Nextflow config interpolation and smoke report paths were corrected. These were development failures, not ignored failures in the final run.

The full trace, smoke trace, cache trace, integration checks and numerical test output are included. `full_trace.tsv` records the first successful execution; subsequent report-only changes improved layout without changing the statistical analysis. All raw and clean data, statistical tables and plotting source are available for inspection.

Text sources were normalised to LF for portable Git checkout. The source manifest records these canonical bytes. Generated Nextflow HTML reports are retained locally but excluded from Git; TSV execution traces are committed.
