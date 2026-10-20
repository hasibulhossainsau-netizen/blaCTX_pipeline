nextflow.enable.dsl=2


process FETCH_SRA {
    tag "$sample"
    conda "${projectDir}/workflow/envs/qc.yaml"
    // raw reads are intermediate only (not published); FASTP deletes them once trimmed

    input:
    tuple val(sample), val(accession), val(layout)
    
    output:
    tuple val(sample), val(layout), path("${sample}_R*.fastq.gz"), emit: reads
    
    script:
    """
    python ${projectDir}/workflow/scripts/fetch_genomes.py \
        --mode sra --layout ${layout} --accession ${accession} \
        --sample-id ${sample} --outdir . --threads ${task.cpus}
    """
    
    stub:
    """
    touch ${sample}_R1.fastq.gz
    if [[ "${layout}" == *PE ]]; then touch ${sample}_R2.fastq.gz; fi
    """
}

process FASTP {
    tag "$sample"
    conda "${projectDir}/workflow/envs/qc.yaml"
    publishDir "${params.results}/01_qc/fastp", mode: 'copy', pattern: '*.fastp.{html,json}'

    input:
    tuple val(sample), val(layout), path(reads)

    output:
    tuple val(sample), val(layout), path("${sample}_R*.{paired,trimmed}.fastq.gz"), emit: trimmed
    tuple val(sample), path("${sample}.fastp.html"), path("${sample}.fastp.json"), emit: reports

    script:
    def common = "--html ${sample}.fastp.html --json ${sample}.fastp.json --report_title ${sample} " +
        "--thread ${task.cpus} --cut_front --cut_tail --cut_mean_quality ${params.fastp_min_quality} " +
        "--cut_window_size ${params.fastp_window_size} --length_required ${params.fastp_min_length}"
    // storage: raw reads are no longer needed once trimmed, so delete them from the upstream work dir
    def cleanup = params.delete_intermediates ? "rm -f -- \$(readlink -f ${reads})" : ''
    if (layout == 'PE')
        """
        fastp --in1 ${reads[0]} --in2 ${reads[1]} \
            --out1 ${sample}_R1.paired.fastq.gz --out2 ${sample}_R2.paired.fastq.gz \
            --detect_adapter_for_pe ${common}
        ${cleanup}
        """
    else
        """
        fastp --in1 ${reads} --out1 ${sample}_R1.trimmed.fastq.gz ${common}
        ${cleanup}
        """

    stub:
    def suffix = layout == 'PE' ? 'paired' : 'trimmed'
    """
    touch ${sample}_R1.${suffix}.fastq.gz
    if [[ "${layout}" == PE ]]; then touch ${sample}_R2.${suffix}.fastq.gz; fi
    printf '{}\\n' > ${sample}.fastp.json
    touch ${sample}.fastp.html
    """
}

process ASSEMBLY {
    tag "$sample"
    conda "${projectDir}/workflow/envs/assembly.yaml"
    publishDir "${params.results}/02_assembly", mode: 'copy'
    
    input:
    tuple val(sample), val(layout), path(reads)
    
    output:
    tuple val(sample), path("${sample}.assembly.fasta"), emit: assembly
    
    script:
    // Paired-end: Shovill + SKESA. Single-end: Shovill does not support it, so use SPAdes directly.
    if (layout == 'PE')
        """
        shovill --assembler skesa --R1 ${reads[0]} --R2 ${reads[1]} --gsize ${params.genome_size} \
            --outdir assembly --cpus ${task.cpus} --force
        cp assembly/contigs.fa ${sample}.assembly.fasta
        ${params.delete_intermediates ? 'rm -rf assembly' : ''}
        """
    else
        """
        spades.py --isolate -s ${reads} -o assembly -t ${task.cpus} -m ${params.spades_memory}
        cp assembly/contigs.fasta ${sample}.assembly.fasta
        ${params.delete_intermediates ? 'rm -rf assembly' : ''}
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
    busco -i ${fasta} -l ${params.busco_db} \
        -o ${sample} --out_path busco --mode genome -c ${task.cpus} -f --offline
    cp busco/${sample}/short_summary*.txt ${sample}.busco.txt
    ${params.delete_intermediates ? 'rm -rf busco' : ''}
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
    tuple val(sample), val(layout), path(reads), path(fasta)
    
    output:
    tuple val(sample), path("${sample}.coverage.txt"), emit: coverage
    
    script:
    """
    python ${projectDir}/workflow/scripts/estimate_coverage.py \
        --r1 ${reads[0]} ${layout == 'PE' ? "--r2 ${reads[1]}" : ''} --assembly ${fasta} --output ${sample}.coverage.txt
    ${params.delete_intermediates ? "rm -f -- \$(readlink -f ${reads})" : ''}
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
        --coverage-file ${coverage} --busco-min ${params.busco_min} --coverage-min ${params.coverage_min} \
        --max-contigs ${params.max_contigs} --min-contig-len ${params.min_contig_len} --outdir . \
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
        --aligner ${params.rgi_aligner} --out-prefix ${sample}.rgi --threads ${task.cpus}
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
        --species ${params.resfinder_species} --outdir .

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
    
    mkdir -p rgi_dir
    for result in *.rgi.txt; do
        mv "\${result}" rgi_dir/
    done

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

process MULTIQC {
    conda "${projectDir}/workflow/envs/qc.yaml"
    publishDir "${params.results}/01_qc/multiqc", mode: 'copy'

    input:
    path fastp_json
    path busco
    path status
    path rgi
    path resfinder

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    shopt -s nullglob
    mkdir -p mqc_input

    # BUSCO's MultiQC module only recognises short_summary.* file names
    for f in *.busco.txt; do
        cp "\${f}" "mqc_input/short_summary.specific.${params.busco_lineage}.\${f%.busco.txt}.txt"
    done
    cp *.fastp.json mqc_input/

    python ${projectDir}/workflow/scripts/make_multiqc_custom.py \
        --input-dir . --target-gene ${params.target_gene} --outdir mqc_input

    multiqc mqc_input --filename multiqc_report.html --force \
        --title "blaCTX pipeline QC"
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data
    """
}

process BAKTA {
    tag "$sample"
    conda "${projectDir}/workflow/envs/bakta.yaml"
    publishDir "${params.results}/07_phylogenomics/bakta", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}.gff3"), emit: gff
    tuple val(sample), path("${sample}.gbk"), path("${sample}.faa"), emit: extras

    script:
    """
    bakta --db ${params.bakta_db} --threads ${task.cpus} \
        --prefix ${sample} --locus-tag ${sample} --output bakta_out --force ${fasta}
    cp bakta_out/${sample}.gff3 ${sample}.gff3
    cp bakta_out/${sample}.gbff ${sample}.gbk
    cp bakta_out/${sample}.faa ${sample}.faa
    ${params.delete_intermediates ? 'rm -rf bakta_out' : ''}
    """

    stub:
    """
    touch ${sample}.gff3 ${sample}.gbk ${sample}.faa
    """
}

process PANAROO {
    conda "${projectDir}/workflow/envs/panaroo.yaml"
    publishDir "${params.results}/07_phylogenomics/panaroo", mode: 'copy'

    input:
    path gffs

    output:
    path "panaroo_out/*", emit: all
    path "panaroo_out/core_gene_alignment.aln", emit: core_aln

    script:
    """
    panaroo -i ${gffs} -o panaroo_out -t ${task.cpus} \
        --clean-mode ${params.panaroo_mode} --remove-invalid-genes -a core --aligner mafft \
        --core_threshold ${params.panaroo_core_threshold}
    """

    stub:
    """
    mkdir panaroo_out
    touch panaroo_out/core_gene_alignment.aln
    """
}

process GUBBINS {
    conda "${projectDir}/workflow/envs/gubbins.yaml"
    publishDir "${params.results}/07_phylogenomics/gubbins", mode: 'copy'

    input:
    path core_aln

    output:
    path "gubbins.*", emit: all
    path "gubbins.filtered_polymorphic_sites.fasta", emit: filtered_aln

    script:
    """
    run_gubbins.py --prefix gubbins --threads ${task.cpus} \
        --tree-builder ${params.gubbins_tree_builder} --iterations ${params.gubbins_iterations} \
        ${core_aln}
    """

    stub:
    """
    touch gubbins.filtered_polymorphic_sites.fasta
    """
}

process IQTREE {
    conda "${projectDir}/workflow/envs/iqtree.yaml"
    publishDir "${params.results}/07_phylogenomics/iqtree", mode: 'copy'

    input:
    path aln

    output:
    path "iqtree.*", emit: all
    path "iqtree.treefile", emit: tree

    script:
    """
    iqtree2 -s ${aln} --prefix iqtree -m ${params.iqtree_model} -bb ${params.iqtree_bootstraps} \
        -nt ${task.cpus} -seed ${params.iqtree_seed}
    """

    stub:
    """
    touch iqtree.treefile iqtree.iqtree
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            // sample id comes from 'sample_id', falling back to the SRA 'Run' column
            def sample = row.sample_id?.toString()?.trim() ?: row.Run?.toString()?.trim()
            def layout = row.source_type?.toString()?.trim()
            def accession = row.sra_accession?.toString()?.trim() ?: row.Run?.toString()?.trim()
            
            if (!sample || !accession || !['SRA_PE', 'PE', 'SRA_SE', 'SE'].contains(layout)) {
                error "Unsupported or incomplete sample row: sample_id=${sample}, source_type=${layout}, sra_accession=${accession}"
            }
            tuple(sample, accession, layout)
        }

    FETCH_SRA(samples_ch)
    // downstream steps only need PE or SE
    reads_ch = FETCH_SRA.out.reads.map { sample, layout, reads -> tuple(sample, layout.endsWith('SE') ? 'SE' : 'PE', reads) }
    FASTP(reads_ch)
    ASSEMBLY(FASTP.out.trimmed)
    BUSCO(ASSEMBLY.out.assembly)
    
    COVERAGE(FASTP.out.trimmed
        .join(ASSEMBLY.out.assembly, by: 0))
        
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

    // Phylogenomics: each step needs a minimum number of genomes
    BAKTA(passing_assemblies)
    PANAROO(BAKTA.out.gff.map { it[1] }.collect()
        .filter { it.size() >= params.min_phylo_samples })
    BAKTA.out.gff.count().subscribe { n ->
        if (n < params.min_phylo_samples)
            log.warn "Only ${n} genome(s) passed QC (< min_phylo_samples=${params.min_phylo_samples}); skipping PANAROO/GUBBINS/IQTREE"
    }
    GUBBINS(PANAROO.out.core_aln)
    IQTREE(GUBBINS.out.filtered_aln)
    
    SUMMARY(
        Channel.fromPath(params.samples),
        FILTER_ASSEMBLY.out.status.map { it[1] }.collect(),
        RGI.out.results.map { it[1] }.collect().ifEmpty([]),
        RESFINDER.out.results.map { it[1] }.collect().ifEmpty([])
    )

    MULTIQC(
        FASTP.out.reports.map { it[2] }.collect(),
        BUSCO.out.summary.map { it[1] }.collect(),
        FILTER_ASSEMBLY.out.status.map { it[1] }.collect(),
        RGI.out.results.map { it[1] }.collect().ifEmpty([]),
        RESFINDER.out.results.map { it[1] }.collect().ifEmpty([])
    )
}
