nextflow.enable.dsl=2

params.samples     = 'config/samples.tsv'
params.results     = 'results'
params.threads     = 10
params.target_gene = 'blaCTX-M-15'

// Database path-গুলোকে projectDir দিয়ে Absolute করা হলো
params.card_json    = "${projectDir}/resources/card_data/card.json"
params.resfinder_db = "${projectDir}/resources/resfinder/db_resfinder"

process FETCH_SRA {
    tag "$sample"
    conda "${projectDir}/workflow/envs/qc.yaml"
    publishDir "${params.results}/00_raw_reads", mode: 'copy'
    
    input:
    tuple val(sample), val(accession), val(layout)
    
    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), emit: reads
    
    script:
    """
    python ${projectDir}/workflow/scripts/fetch_genomes.py \
        --mode sra --layout ${layout} --accession ${accession} \
        --sample-id ${sample} --outdir . --threads ${task.cpus}
    """
    
    stub:
    """
    touch ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
    """
}

process FASTQC_RAW {
    tag "$sample"
    conda "${projectDir}/workflow/envs/qc.yaml"
    publishDir "${params.results}/01_qc/fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample), path(r1), path(r2)
    
    output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip"), emit: reports
    
    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
    
    stub:
    """
    touch ${sample}_R1.fastq_fastqc.html ${sample}_R2.fastq_fastqc.html
    touch ${sample}_R1.fastq_fastqc.zip ${sample}_R2.fastq_fastqc.zip
    """
}

process TRIMMOMATIC {
    tag "$sample"
    conda "${projectDir}/workflow/envs/qc.yaml"
    publishDir "${params.results}/01_qc/trimmed", mode: 'copy'
    
    input:
    tuple val(sample), path(r1), path(r2)
    
    output:
    tuple val(sample), path("${sample}_R1.paired.fastq.gz"), path("${sample}_R2.paired.fastq.gz"), emit: trimmed
    
    script:
    """
    trimmomatic -Xmx8g PE -threads ${task.cpus} -phred33 ${r1} ${r2} \
        ${sample}_R1.paired.fastq.gz /dev/null \
        ${sample}_R2.paired.fastq.gz /dev/null \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    """
    
    stub:
    """
    touch ${sample}_R1.paired.fastq.gz ${sample}_R2.paired.fastq.gz
    """
}

process ASSEMBLY {
    tag "$sample"
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/02_assembly", mode: 'copy'
    
    input:
    tuple val(sample), path(r1), path(r2)
    
    output:
    tuple val(sample), path("${sample}.assembly.fasta"), emit: assembly
    
    script:
    """
    shovill --assembler skesa --R1 ${r1} --R2 ${r2} \
        --outdir assembly --cpus ${task.cpus} --force
    cp assembly/contigs.fa ${sample}.assembly.fasta
    """
    
    stub:
    """
    printf '>stub\\nACGT\\n' > ${sample}.assembly.fasta
    """
}

process BUSCO {
    tag "$sample"
    conda "${projectDir}/workflow/envs/busco.yaml"
    publishDir "${params.results}/03_busco", mode: 'copy'
    
    input:
    tuple val(sample), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.busco.txt"), emit: summary
    
    script:
    """
    busco -i ${fasta} -l ${projectDir}/resources/busco_db \
        -o ${sample} --out_path busco --mode genome -c ${task.cpus} -f --offline
    cp busco/${sample}/short_summary*.txt ${sample}.busco.txt
    """
    
    stub:
    """
    printf 'C:100%%[S:100%%,D:0%%],F:0%%,M:0%%,n:1\\n' > ${sample}.busco.txt
    """
}

process COVERAGE {
    tag "$sample"
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/02_assembly", mode: 'copy'
    
    input:
    tuple val(sample), path(r1), path(r2), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.coverage.txt"), emit: coverage
    
    script:
    """
    python ${projectDir}/workflow/scripts/estimate_coverage.py \
        --r1 ${r1} --r2 ${r2} --assembly ${fasta} --output ${sample}.coverage.txt
    """
    
    stub:
    """
    printf '100.00\\t100\\t1\\n' > ${sample}.coverage.txt
    """
}

process FILTER_ASSEMBLY {
    tag "$sample"
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/04_filtered_assemblies", mode: 'copy'
    
    input:
    tuple val(sample), path(fasta), path(busco), path(coverage)
    
    output:
    tuple val(sample), path("${sample}.filter_status.json"), emit: status
    
    script:
    """
    python ${projectDir}/workflow/scripts/validate_assemblies.py filter-one \
        --sample-id ${sample} --fasta ${fasta} --busco-summary ${busco} \
        --coverage-file ${coverage} --busco-min 95 --coverage-min 50 \
        --max-contigs 200 --min-contig-len 500 --outdir . \
        --status-json ${sample}.filter_status.json
    """
    
    stub:
    """
    printf '{"sample_id":"%s","pass":true}\\n' '${sample}' > ${sample}.filter_status.json
    """
}

process PROMOTE {
    tag "$sample"
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/04_filtered_assemblies", mode: 'copy'
    
    input:
    tuple val(sample), path(status), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.pass.fasta"), emit: passing
    
    script:
    """
    python ${projectDir}/workflow/scripts/validate_assemblies.py promote \
        --status-json ${status} --fasta ${fasta} --output ${sample}.pass.fasta
    """
    
    stub:
    """
    cp ${fasta} ${sample}.pass.fasta
    """
}

process RGI {
    tag "$sample"
    conda "${projectDir}/workflow/envs/rgi.yaml"
    publishDir "${params.results}/05_amr/rgi", mode: 'copy'
    
    input:
    tuple val(sample), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.rgi.txt"), path("${sample}.rgi.json"), emit: results
    
    script:
    """
    python ${projectDir}/workflow/scripts/amr_screening.py rgi \
        --fasta ${fasta} --card-json ${params.card_json} \
        --aligner BLAST --out-prefix ${sample}.rgi --threads ${task.cpus}
    """
    
    stub:
    """
    touch ${sample}.rgi.txt ${sample}.rgi.json
    """
}

process RESFINDER {
    tag "$sample"
    conda "${projectDir}/workflow/envs/resfinder.yaml"
    publishDir "${params.results}/05_amr/resfinder", mode: 'copy'
    
    input:
    tuple val(sample), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.ResFinder_results.txt"), emit: results
    
    script:
    """
    python ${projectDir}/workflow/scripts/amr_screening.py resfinder \
        --fasta ${fasta} --db-path ${params.resfinder_db} \
        --species Klebsiella --outdir .

    if [[ -f ResFinder_results_tab.txt ]]; then
        cp ResFinder_results_tab.txt ${sample}.ResFinder_results.txt
    elif [[ -f ResFinder_results.txt ]]; then
        cp ResFinder_results.txt ${sample}.ResFinder_results.txt
    else
        echo "ResFinder completed without a recognized tabular output" >&2
        exit 1
    fi
    """
    
    stub:
    """
    touch ${sample}.ResFinder_results.txt
    """
}

process SUMMARY {
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/06_metadata", mode: 'copy'
    
    input:
    path samples_tsv
    path statuses
    path rgi
    path resfinder
    
    output:
    path 'cohort_metadata_summary.tsv'
    path 'cohort_metadata_summary.xlsx'
    
    script:
    """
    shopt -s nullglob
    
    # RGI ফাইলগুলো সাজানো
    mkdir -p rgi_dir
    for result in *.rgi.txt; do
        mv "\${result}" rgi_dir/
    done

    # ResFinder ফাইলগুলো সাজানো
    mkdir -p resfinder_dir
    for result in *.ResFinder_results.txt; do
        sample=\${result%.ResFinder_results.txt}
        mkdir -p resfinder_dir/\${sample}
        mv "\${result}" resfinder_dir/\${sample}/ResFinder_results.txt
    done

    python ${projectDir}/workflow/scripts/validate_assemblies.py build-summary \
        --samples-tsv ${samples_tsv} \
        --filter-dir . \
        --rgi-dir rgi_dir \
        --resfinder-dir resfinder_dir \
        --target-gene ${params.target_gene} \
        --out-tsv cohort_metadata_summary.tsv \
        --out-xlsx cohort_metadata_summary.xlsx
    """
    
    stub:
    """
    printf 'sample_id\tqc_pass\n' > cohort_metadata_summary.tsv
    touch cohort_metadata_summary.xlsx
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample = row.sample_id?.toString()?.trim()
            def layout = row.source_type?.toString()?.trim()
            def accession = row.sra_accession?.toString()?.trim()
            if (!sample || !accession || !['SRA_PE', 'PE'].contains(layout)) {
                error "Unsupported or incomplete sample row: sample_id=${sample}, source_type=${layout}, sra_accession=${accession}"
            }
            tuple(sample, accession, layout)
        }

    FETCH_SRA(samples_ch)
    FASTQC_RAW(FETCH_SRA.out.reads)
    TRIMMOMATIC(FETCH_SRA.out.reads)
    ASSEMBLY(TRIMMOMATIC.out.trimmed)
    BUSCO(ASSEMBLY.out.assembly)
    
    COVERAGE(TRIMMOMATIC.out.trimmed
        .join(ASSEMBLY.out.assembly, by: 0)
        .map { sample, r1, r2, fasta -> tuple(sample, r1, r2, fasta) })
        
    FILTER_ASSEMBLY(ASSEMBLY.out.assembly
        .join(BUSCO.out.summary, by: 0)
        .join(COVERAGE.out.coverage, by: 0)
        .map { sample, fasta, busco_summary, coverage_file ->
            tuple(sample, fasta, busco_summary, coverage_file)
        })
        
    passing_status = FILTER_ASSEMBLY.out.status
        .filter { sample, status -> status.text =~ /"pass"\s*:\s*true/ }

    PROMOTE(passing_status.join(ASSEMBLY.out.assembly, by: 0))
    passing_assemblies = PROMOTE.out.passing
        .filter { sample, fasta -> fasta.size() > 0 }
    RGI(passing_assemblies)
    RESFINDER(passing_assemblies)
    
    SUMMARY(
        Channel.fromPath(params.samples),
        FILTER_ASSEMBLY.out.status.map { it[1] }.collect(),
        RGI.out.results.map { it[1] }.collect(),
        RESFINDER.out.results.map { it[1] }.collect()
    )
}