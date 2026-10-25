nextflow.enable.dsl=2
params.n_workers=500
params.seed=2026
params.outdir='results'
params.python='python3'
params.rscript='Rscript'
params.quarto='quarto'
process SIMULATE {
 publishDir "${params.outdir}", mode:'copy'
 input: path source
 output: path 'raw'
 script:
 """
 ${params.python} ${source} --n ${params.n_workers} --seed ${params.seed} --out raw
 """
}
process QC {
 publishDir "${params.outdir}", mode:'copy'
 input: path raw; path source
 output: path 'clean'
 script:
 """
 ${params.python} ${source} ${raw} --out clean
 """
}
process LONGITUDINAL_INPUT {
 publishDir "${params.outdir}", mode:'copy'
 input: path clean; path sources
 output: path 'longitudinal_input'
 script:
 """
 ${params.python} prepare_longitudinal.py ${clean}/workers.csv ${params.seed + 1}
 """
}
process ANALYSE_R {
 publishDir "${params.outdir}", mode:'copy'
 input: path clean; path raw; path longitudinal_input; path source
 output: path 'analysis', emit: analysis
 path 'longitudinal', emit: longitudinal
 path 'figures', emit: figures
 path 'longitudinal_figures', emit: longitudinal_figures
 script:
 """
 ${params.rscript} ${source} ${clean} ${raw} ${longitudinal_input}
 """
}
process CROSSCHECK {
 publishDir "${params.outdir}", mode:'copy'
 input: path clean; path raw; path analysis; path longitudinal; path sources
 output: path 'crosscheck'
 script:
 """
 ${params.python} crosscheck.py ${clean} ${raw} ${analysis} ${longitudinal}
 """
}
process REPORT {
 publishDir "${params.outdir}", mode:'copy'
 input: path analysis; path clean; path figures; path longitudinal; path longitudinal_figures; path crosscheck; path source; path reference
 output: path 'report'
 script:
 """
 ${params.quarto} render ${source} --to docx --output OccuSurv_report.docx
 mkdir report
 mv OccuSurv_report.docx report/
 """
}
workflow {
 raw=SIMULATE(file("${projectDir}/bin/simulate.py"))
 clean=QC(raw,file("${projectDir}/bin/qc.py"))
 inputs=LONGITUDINAL_INPUT(clean,[file("${projectDir}/bin/prepare_longitudinal.py"),file("${projectDir}/bin/longitudinal.py")])
 r=ANALYSE_R(clean,raw,inputs,file("${projectDir}/bin/analyse.R"))
 check=CROSSCHECK(clean,raw,r.analysis,r.longitudinal,[file("${projectDir}/bin/crosscheck.py"),file("${projectDir}/bin/analyse.py"),file("${projectDir}/bin/longitudinal.py")])
 REPORT(r.analysis,clean,r.figures,r.longitudinal,r.longitudinal_figures,check,file("${projectDir}/report/OccuSurv.qmd"),file("${projectDir}/report/reference.docx"))
}
