nextflow.enable.dsl=2
params.n_workers=500
params.seed=2026
params.outdir='results'
params.python='python3'
params.report_python='python3'
params.report_script="${projectDir}/bin/report.py"
process SIMULATE {
 publishDir "${params.outdir}", mode:'copy'
 input: path script
 output: path 'raw'
 script:
 """
 ${params.python} ${script} --n ${params.n_workers} --seed ${params.seed} --out raw
 """
}
process QC {
 publishDir "${params.outdir}", mode:'copy'
 input: path raw; path script
 output: path 'clean'
 script:
 """
 ${params.python} ${script} ${raw} --out clean
 """
}
process ANALYSE {
 publishDir "${params.outdir}", mode:'copy'
 input: path clean; path raw; path script
 output: path 'analysis'
 script:
 """
 ${params.python} ${script} ${clean} ${raw}/truth.json --out analysis
 """
}
process FIGURES {
 publishDir "${params.outdir}", mode:'copy'
 input: path analysis; path script
 output: path 'figures'
 script:
 """
 ${params.python} ${script} ${analysis} --out figures
 """
}
process REPORT {
 publishDir "${params.outdir}", mode:'copy'
 input: path analysis; path clean; path figures; path methods; path equations; path report_source
 output: path 'report'
 script:
 """
 ${params.report_python} ${params.report_script} --analysis analysis --clean clean --figures figures --methods methods.json --equations equations --out report
 """
}
workflow {
 raw=SIMULATE(file("${projectDir}/bin/simulate.py"))
 clean=QC(raw,file("${projectDir}/bin/qc.py"))
 analysis=ANALYSE(clean,raw,file("${projectDir}/bin/analyse.py"))
 figures=FIGURES(analysis,file("${projectDir}/bin/figures.py"))
 REPORT(analysis,clean,figures,file("${projectDir}/docs/methods.json"),file("${projectDir}/assets/equations"),file("${projectDir}/bin/report.py"))
}
