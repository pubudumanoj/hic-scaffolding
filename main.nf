#!/usr/bin/env nextflow

nextflow.enable.dsl=2

println \
"""
=================================
 H I - C   S C A F F O L D I N G   P I P E L I N E
 BY PUBUDU NAWARATHNA
=================================
Workflow Information:
---------------------
  Project Directory:  ${workflow.projectDir}
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Container Engine:   ${workflow.containerEngine}
  Profiles:           ${workflow.profile}

"""

// Acknowledgement: Rob Syme


workflow {

    channel.fromPath(params.REF) | index
    ch_raw_fastq = channel.fromFilePairs( params.in_dir + params.fastq )  // [ id, [fwd, rev] ]
    align_R1 ( ch_raw_fastq, index.out.first() )  // | view\\
    // align_R1.out.verbiage.view()
    align_R2 ( ch_raw_fastq, index.out.first() )  // | view\\
    filter_R2(align_R2.out)
    filter_R1(align_R1.out)
    // filter_R2.out.verbiage.view()

    filter_R1.out
    | join(filter_R2.out) // Joins two channels that share a common first element in the tuple
    | set { ch_pairs }

    pair_reads( ch_pairs, index.out.first() )

    // This gives a channel with form [id, read1, read2]
    // pair_reads(filter_R1.out, filter_R2.out, index.out.first())
    // pair_reads.out.verbiage.view()
    add_read_group(channel.of(params.LABEL).first(),  pair_reads.out)
    
    merge_technical_replicates(add_read_group.out.toList(), channel.of(params.LABEL).first())
    // merge_technical_replicates.out.verbiage.view()

    mark_duplicates(merge_technical_replicates.out, channel.of(params.LABEL).first())
    // mark_duplicates.out.verbiage.view()

    sam_index(mark_duplicates.out)

    stats(mark_duplicates.out)

    bam_to_bed(mark_duplicates.out)
    
    sort_bed(bam_to_bed.out)

    fasta_faidx_index(channel.fromPath(params.REF))

    salsa2_scaffolding(channel.fromPath(params.REF).first(), fasta_faidx_index.out.first(), sort_bed.out, channel.fromList(params.iterations))

    index_scaffolded_fasta(salsa2_scaffolding.out)

    make_chromosome_sizes(index_scaffolded_fasta.out) 

    sorted_iterated_alignment(salsa2_scaffolding.out)

    make_chromosome_sizes.out | join(sorted_iterated_alignment.out) | set {ch_hic}

    creating_hic_file(ch_hic)


}

process align_R1 {
    module 'mugqic/bwa/0.7.17'
    cpus 1
    publishDir "bam_out"
    label 'align_R1'

    input:
    tuple val(sample_id), path(reads)
    path(index)

    
    output:
    tuple val(sample_id), path("*.bam")
    // stdout emit: verbiage

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    echo ${reads[0]} && \
    bwa mem -t ${task.cpus} -5SP \$INDEX ${reads[0]} | samtools view --threads ${task.cpus} -Sb --output ${sample_id}_R1.bam
    """
}

process align_R2 {
    module 'mugqic/bwa/0.7.17'
    cpus 1
    publishDir "bam_out"
    label 'align_R2'

    input:
    tuple val(sample_id), path(reads)
    path(index)

    
    output:
    tuple val(sample_id), path("*.bam")
    // stdout emit: verbiage

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    echo ${reads[1]} && \
    bwa mem -t ${task.cpus} -5SP \$INDEX ${reads[1]} | samtools view --threads ${task.cpus} -Sb --output ${sample_id}_R2.bam
    """
}

process filter_R1 {
    module 'mugqic/samtools/1.14'
    cpus 1
    publishDir "."
    label 'filter'


    input:
    tuple val(sample_id), path(input)

    // path FILTER

    
    output:
    tuple val(sample_id), path("filtered/*.bam")
    // stdout emit: verbiage

    script:
    def oo = 10
    """
    echo ${input.getSimpleName()}
    mkdir filtered
    samtools view -h $input | filter_five_end.pl | samtools view --threads ${task.cpus} -Sb --output filtered/${sample_id}_R1.bam 
    """
}


process filter_R2 {
    module 'mugqic/samtools/1.14'
    cpus 1
    publishDir "."
    label 'filter'

    input:
    tuple val(sample_id), path(input)
    // path FILTER

    
    output:
    tuple val(sample_id), path("filtered/*.bam")
    // stdout emit: verbiage

    """
    echo ${input.getSimpleName()}
    mkdir filtered
    samtools view -h $input | filter_five_end.pl | samtools view --threads ${task.cpus} -Sb --output filtered/${sample_id}_R2.bam 
    """
}

process index {
    module 'mugqic/bwa/0.7.17'
    // publishDir "out"
    label 'index'

    input:
    path(REF)

    output:
    path("*")
    
    "bwa index -a bwtsw $REF"
}

process pair_reads{

    module 'mugqic/samtools/1.14'
    cpus 1
    publishDir "."
    label 'pair_reads'

    input:
    // tuple val(sample_id_R1), path(input1)
    // tuple val(sample_id_R2), path(input2)
    // Rob's edit Start
    tuple val(sample_id), path(r1), path(r2)
    // Rob's edit End
    path(index)
    

    output:
    // stdout emit: verbiage
    tuple val(sample_id), path("*.bam")

    //sample_id =`printf '%s\\n' "$sample_id_R1" "$sample_id_R2" | sed -e "N;s/^\\(.*\\).*\\n\\1.*$/\\1/"`
    // sample_id =`printf '%s\n' "$sample_id_R1" "$sample_id_R2" | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d'`
    """

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    two_read_bam_combiner.pl ${r1} ${r2} samtools ${params.MAPQ_FILTER} | samtools view -bS -t \$INDEX | samtools sort --threads ${task.cpus} -o ${sample_id}.bam 
    """

}

process add_read_group{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."
    label 'add_read_group'

    input:
    val LABEL
    tuple val(sample_id), path(input)

    output:
    path("paired/*.bam")
    // stdout emit: verbiage


    """
    mkdir temp
    mkdir paired
    
    java -Xmx${params.java_memory} -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT ${input} -OUTPUT paired/${input} -ID ${sample_id} -LB ${sample_id} -SM ${LABEL} -PL ILLUMINA -PU none
    """
}

process merge_technical_replicates{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."
    label 'merge_technical_replicates'

    input:
    path input
    val LABEL
    

    output:
    path("merged/*.bam")
    // stdout emit: verbiage

    """

    mkdir temp
    mkdir merged
    inputs=`echo ${input} | sed 's\\ \\ -I \\g'`
    echo \$inputs
    java -Xmx${params.java_memory} -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar MergeSamFiles -I \$inputs -OUTPUT merged/${LABEL}.bam -USE_THREADING TRUE -ASSUME_SORTED TRUE -VALIDATION_STRINGENCY LENIENT
    """

}

process mark_duplicates{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."
    label 'mark_duplicates'

    input:
    path input
    // stdin
    val LABEL
    

    output:
    path("deduplicated/*.bam")
    // stdout emit: verbiage

    """
    echo $input
    mkdir temp
    mkdir deduplicated
    java -Xmx${params.java_memory} -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar MarkDuplicates INPUT=${input} OUTPUT=deduplicated/${LABEL}.bam METRICS_FILE=deduplicated/metrics.${LABEL}.txt TMP_DIR=temp ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

    """

}

process sam_index{

    module 'mugqic/samtools/1.14'
    cpus 1
    publishDir "."
    label 'sam_index'

    input:
    path input
    

    output:
    path("deduplicated/*.bai")
    stdout emit: verbiage

    """
    mkdir deduplicated
    samtools index ${input} deduplicated/${input}.bai
    """

}

process stats{

    cpus 1
    publishDir "."
    label 'stats'

    input:
    path input
    

    output:
    path("deduplicated/*.stats")
    stdout emit: verbiage

    """
    mkdir deduplicated
    get_stats.pl ${input} > deduplicated/${input}.stats
    """

}

///////
// Starting sals2 processing

process bam_to_bed{
    
    module 'mugqic/bedtools/2.30.0'
    cpus 1
    label 'bam_to_bed'

    input:
    path input
  
    output:
    path("*.bed")
  
    """
    bamToBed -i $input > ${input}.bed

    """
}

process sort_bed{

    cpus 1
    label 'sort_bed'

    input:
    path input
    
    output:
    path("*.sort.bed")

    
    """
    temp=$input
    bamfile=\${temp%.*}
    echo \$bamfile
    sort -k 4 $input -S ${params.java_memory} > \$bamfile.sort.bed
    
    """
}

process fasta_faidx_index {

    module 'mugqic/samtools/1.14'
    publishDir "."
    cpus 1

    input:
    path fasta
        
    output:
    path("*.fai")

    """
    samtools faidx $fasta > ${fasta}.fai
    """
}
process salsa2_scaffolding {

    module 'mugqic/python/2.7.14'
    publishDir "scaffolding/iteration_$iteration"
    cpus 1
    label 'salsa2_scaffolding'

    input:
    path fasta
    path fasta_index
    path input_bed
    val iteration
        
    output:
    // path("scaffolds_i"+ iteration + "/*")
    tuple val("scaffolds_i$iteration"), path("*")

    """
    mkdir scaffolds_i${iteration}
    run_pipeline.py -a $fasta -l ${fasta_index} -b ${input_bed} -e ${params.hic_enzyme} -o scaffolds_i${iteration} -i $iteration -m yes
    """


}


////////////////
// steps for creating .hic file

process index_scaffolded_fasta {
    
    module 'mugqic/samtools/1.14'
    cpus 1
    label 'index_scaffolded_fasta'

    input:
    tuple val(iteration), path(salsa2_output)

    output:
    tuple val(iteration), path("*.fai")
    
    """
    FASTA=`find -L ./ -name "*FINAL.fasta"`
    samtools faidx $iteration/scaffolds_FINAL.fasta -o scaffolds_FINAL.fasta.fai
    """
}

process make_chromosome_sizes {
    
    module 'mugqic/samtools/1.14'
    cpus 1
    label 'make_chromosome_sizes'

    input:
    tuple val(iteration), path(index)

    output:
    tuple val(iteration), path("chromosome_sizes.tsv")
    //stdout

    """
    echo $index $iteration
    cut -f 1,2 $index | head -n ${params.scaffolding_count} > chromosome_sizes.tsv


    """
}

process sorted_iterated_alignment{

    module 'mugqic/python/2.7.14'
    cpus 1
    label 'sorted_iterated_alignment'

    input:
    tuple val(iteration), path(salsa2_output)

    output:
    tuple val(iteration), path("alignments_sorted.txt")
    // stdout

    """
    alignments2txt.py -b $iteration/alignment_iteration_1.bed  -a $iteration/scaffolds_FINAL.agp -l $iteration/scaffold_length_iteration_1 | \\
    awk -v OFS="\\t" '{if (\$2 > \$6) {print \$1","\$6","\$7","\$8","\$5","\$2","\$3","\$4} else {print }}'  > alignments_sorted.txt
    """

}


process creating_hic_file{

    module 'mugqic/java/openjdk-jdk1.8.0_72'
    cpus 1
    memory '2GB'
    label 'creating_hic_file'
    publishDir "scaffolding/$iteration"

    input:
    tuple val(iteration), path(chromosome_sizes), path(alignments_sorted)

    output:
    tuple val(iteration), path("*.hic")
    // stdout

    """
    mkdir ${iteration}
    java -Xmx${params.java_memory} -jar ${params.juicer} pre -j ${task.cpus} $alignments_sorted ${iteration}/salsa_${iteration}.hic $chromosome_sizes
    """

}