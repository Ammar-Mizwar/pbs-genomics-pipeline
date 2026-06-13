nextflow.enable.dsl=2

params.outdir = "results"

process POPULATION_SCAN {
    publishDir "${params.outdir}/01_scans", mode: 'copy'

    input:
    path bed
    path bim
    path fam

    output:
    path "Master_Genomic_Scan_FST_DAF.csv", emit: daf
    path "Master_PBS_Results.csv", emit: pbs
    path "results_fst_chr*.weir.fst"

    script:
    """
    python ${projectDir}/scripts/fst_daf_admixture.py
    python ${projectDir}/scripts/pbs.py
    """
}

process FILTER_AND_PLOT {
    publishDir "${params.outdir}/02_filtered", mode: 'copy'

    input:
    path daf_csv
    path pbs_csv

    output:
    path "Target_Regions_Global_FDR_Significant.csv", emit: sig_csv
    path "Manhattan_FDR.png"

    script:
    """
    python ${projectDir}/scripts/filtering_fdr.py
    python ${projectDir}/scripts/plotting.py
    """
}

// ==========================================
// ANNOTATION PHASE (Now Parallelized!)
// ==========================================

process RUN_VEP {
    publishDir "${params.outdir}/03_annotations", mode: 'copy'

    input:
    path sig_csv

    output:
    path "VEP_Annotations_Results.csv", emit: vep_csv

    script:
    """
    python ${projectDir}/scripts/VEP.py
    """
}

process RUN_ENCODE {
    publishDir "${params.outdir}/03_annotations", mode: 'copy'

    input:
    path vep_csv // Waits for RUN_VEP to finish

    output:
    path "Intergenic_ENCODE_Annotations.csv"

    script:
    """
    python ${projectDir}/scripts/ENCODE.py
    """
}

process RUN_POLYPHEN {
    publishDir "${params.outdir}/03_annotations", mode: 'copy'

    input:
    path vep_csv  // Waits for RUN_VEP to finish
    path bim_file // Also pulls in the original .bim file

    output:
    path "Missense_PolyPhen_SIFT_Results.csv"

    script:
    """
    python ${projectDir}/scripts/PolyPhen.py
    """
}

// ==========================================
// WORKFLOW EXECUTION
// ==========================================

workflow {
    bed_file = file("${projectDir}/Plink_Output.cleaned.top.snpqc.bed")
    bim_file = file("${projectDir}/Plink_Output.cleaned.top.snpqc.bim")
    fam_file = file("${projectDir}/Plink_Output.cleaned.top.snpqc.fam")

    // 1. Run FST, DAF, and PBS 
    scans = POPULATION_SCAN(bed_file, bim_file, fam_file) 

    // 2. Filter results and generate Manhattan plot
    filtered = FILTER_AND_PLOT(scans.daf, scans.pbs)

    // 3. Run VEP Annotation first (depends on filtered output)
    vep_results = RUN_VEP(filtered.sig_csv)

    // 4. Run ENCODE and PolyPhen simultaneously!
    RUN_ENCODE(vep_results.vep_csv)
    RUN_POLYPHEN(vep_results.vep_csv, bim_file)
}