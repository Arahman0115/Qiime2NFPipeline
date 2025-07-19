#!/usr/bin/env nextflow

params.manifest = "manifest.csv"
params.metadata = "metadata.tsv"
params.sampling_depth = 10000

process import_fastq {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path manifest_file
    output:
    path "paired-end-demux.qza"
    script:
    """
    qiime tools import \\
      --type 'SampleData[PairedEndSequencesWithQuality]' \\
      --input-path ${manifest_file} \\
      --output-path paired-end-demux.qza \\
      --input-format PairedEndFastqManifestPhred33V2
    """
}

process demux_summary {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path imported_qza
    output:
    path "demux.qzv"
    script:
    """
    qiime demux summarize \\
    --i-data ${imported_qza} \\
    --o-visualization demux.qzv
    """
}

process dada2_denoise {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path demux_qza
    output:
    path "table.qza"
    path "rep-seqs.qza"
    path "dada2-stats.qza"
    script:
    """
    qiime dada2 denoise-paired \\
      --i-demultiplexed-seqs ${demux_qza} \\
      --p-trim-left-f 0 \\
      --p-trim-left-r 0 \\
      --p-trunc-len-f 240 \\
      --p-trunc-len-r 240 \\
      --o-table table.qza \\
      --o-representative-sequences rep-seqs.qza \\
      --o-denoising-stats dada2-stats.qza
    """
}

process dada2_visuals {
    container 'quay.io/qiime2/core:2023.9'
    publishDir "visuals", mode: 'copy'
    input:
    path table_qza
    path repseqs_qza
    path stats_qza
    output:
    path "dada2-stats.qzv"
    path "feature-table.qzv"
    path "rep-seqs.qzv"
    script:
    """
    qiime metadata tabulate \\
      --m-input-file ${stats_qza} \\
      --o-visualization dada2-stats.qzv

    qiime feature-table summarize \\
      --i-table ${table_qza} \\
      --o-visualization feature-table.qzv

    qiime feature-table tabulate-seqs \\
      --i-data ${repseqs_qza} \\
      --o-visualization rep-seqs.qzv
    """
}

process get_silva_data {
    container 'quay.io/qiime2/metagenome-workshop:2024.10'
    output:
    path "silva-138.2-ssu-nr99-rna-seqs.qza"
    path "silva-138.2-ssu-nr99-tax.qza"
    script:
    """
    qiime rescript get-silva-data \\
      --p-version '138.2' \\
      --p-target 'SSURef_NR99' \\
      --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \\
      --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza
    """
}

process assign_taxonomy {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path repseqs_qza
    path silva_seqs
    path silva_tax
    output:
    path "taxonomy.qza"
    path "taxonomy.qzv"
    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \\
      --i-reference-reads ${silva_seqs} \\
      --i-reference-taxonomy ${silva_tax} \\
      --o-classifier silva-classifier.qza

    qiime feature-classifier classify-sklearn \\
      --i-classifier silva-classifier.qza \\
      --i-reads ${repseqs_qza} \\
      --o-classification taxonomy.qza

    qiime metadata tabulate \\
      --m-input-file taxonomy.qza \\
      --o-visualization taxonomy.qzv
    """
}

process generate_phylogeny {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path repseqs_qza
    output:
    path "rooted-tree.qza"
    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \\
      --i-sequences ${repseqs_qza} \\
      --o-alignment aligned-rep-seqs.qza \\
      --o-masked-alignment masked-aligned-rep-seqs.qza \\
      --o-tree unrooted-tree.qza \\
      --o-rooted-tree rooted-tree.qza
    """
}

process diversity_analysis {
    container 'quay.io/qiime2/core:2023.9'
    input:
    path table_qza
    path rooted_tree_qza
    path metadata_file
    output:
    path "core-metrics-results"
    script:
    """
    qiime diversity core-metrics-phylogenetic \\
      --i-phylogeny ${rooted_tree_qza} \\
      --i-table ${table_qza} \\
      --p-sampling-depth ${params.sampling_depth} \\
      --m-metadata-file ${metadata_file} \\
      --output-dir core-metrics-results
    """
}

workflow {
    manifest_file = file(params.manifest)
    metadata_file = file(params.metadata)

    import_fastq(manifest_file)
    demux_summary(import_fastq.out)
    dada2_denoise(import_fastq.out)
    dada2_visuals(
        dada2_denoise.out[0], // table.qza
        dada2_denoise.out[1], // rep-seqs.qza
        dada2_denoise.out[2]  // dada2-stats.qza
    )
    get_silva_data()
    assign_taxonomy(
        dada2_denoise.out[1], // rep-seqs.qza
        get_silva_data.out[0], // silva-seqs
        get_silva_data.out[1]  // silva-tax
    )
    generate_phylogeny(dada2_denoise.out[1])
    diversity_analysis(
        dada2_denoise.out[0],    // table.qza
        generate_phylogeny.out,  // rooted-tree.qza
        metadata_file
    )
}
