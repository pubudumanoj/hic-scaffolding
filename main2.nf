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

workflow {

    // fastq = channel.fromFilePairs(params.in_dir + params.fastq)
    // ref = channel.fromPath(params.REF)
    // prefix= channel.of(params.LABEL+'_1', params.LABEL+'_2')


    // dir_names = channel.fromPath([params.RAW_DIR, params.RAW_DIR, params.PAIR_DIR, params.REP_DIR, params.MERGE_DIR])
    // make_dirs(dir_names)
    // index(prefix.combine(ref))
    // align(fastq, prefix.combine(channel.fromPath(params.RAW_DIR)))
    // align.out.verbiage.view()

    channel.fromPath(params.REF) | index
    ch_raw_fastq = channel.fromFilePairs( params.in_dir + params.fastq )  // [ id, [fwd, rev] ]
    align_R1 ( ch_raw_fastq, index.out.first() )  // | view\\
    // align_R1.out.verbiage.view()
    align_R2 ( ch_raw_fastq, index.out.first() )  // | view\\
    filter_R2(align_R2.out)
    filter_R1(align_R1.out)
    // filter_R2.out.verbiage.view()

    // Rob's Edits Start
    filter_R1.out
    | join(filter_R2.out) // Joins two channels that share a common first element in the tuple
    | set { ch_pairs }

    pair_reads( ch_pairs, index.out.first() )

    // This gives a channel with form [id, read1, read2]
    // Rob's Edits End
   // pair_reads(filter_R1.out, filter_R2.out, index.out.first())
    // pair_reads.out.verbiage.view()
    add_read_group(channel.of(params.LABEL).first(),  pair_reads.out)
    
    merge_techincal_replicates(add_read_group.out.toList(), channel.of(params.LABEL).first())
    // merge_techincal_replicates.out.verbiage.view()

    mark_duplicates(merge_techincal_replicates.out, channel.of(params.LABEL).first())
    // mark_duplicates.out.verbiage.view()

    sam_index(mark_duplicates.out)

    stats(mark_duplicates.out)

    bam_to_bed(mark_duplicates.out)
    
    sort_bed(bam_to_bed.out)

    fasta_faidx_index(channel.fromPath(params.REF))

    salsa2_scaffolding(channel.fromPath(params.REF).first(), fasta_faidx_index.out.first(), sort_bed.out, channel.of(3,5,10))

    index_scaffolded_fasta(salsa2_scaffolding.out)

    make_chromosome_sizes(index_scaffolded_fasta.out) 


}



process align_R1 {
    module 'mugqic/bwa/0.7.17'
    cpus 1
    publishDir "bam_out"

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

process make_dirs{

input:
    path dir_names

output:
    path "${dir_names}"

script:
"""
    mkdir -p $projectDir/${dir_names} \

"""

}


process index {
    module 'mugqic/bwa/0.7.17'
    // publishDir "out"

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
    two_read_bam_combiner.pl ${r1} ${r2} samtools 1 | samtools view -bS -t \$INDEX | samtools sort --threads ${task.cpus} -o ${sample_id}.bam 
    """

}

process add_read_group{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."

    input:
    val LABEL
    tuple val(sample_id), path(input)

    output:
    path("paired/*.bam")
    // stdout emit: verbiage


    """
    mkdir temp
    mkdir paired
    
    java -Xmx20G -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT ${input} -OUTPUT paired/${input} -ID ${sample_id} -LB ${sample_id} -SM ${LABEL} -PL ILLUMINA -PU none
    """
}

process merge_techincal_replicates{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."

    input:
    path input
    val LABEL
    

    output:
    path("merged/*.bam")
    //stdout emit: verbiage

    """
    echo $input
    mkdir temp
    mkdir merged
    inputs=`printf I=${input} | sed 's\\ \\ I=\\'`
    java -Xmx20G -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar MergeSamFiles \$inputs OUTPUT=merged/${LABEL}.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
    """

}

process mark_duplicates{

    module 'mugqic/java/openjdk-jdk1.8.0_72:mugqic/picard/2.26.6'
    cpus 1
    publishDir "."

    input:
    path input
    val LABEL
    

    output:
    path("deduplicated/*.bam")
    // stdout emit: verbiage

    """
    echo $input
    mkdir temp
    mkdir deduplicated
    java -Xmx20G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar \$PICARD_HOME/picard.jar MarkDuplicates INPUT=${input} OUTPUT=deduplicated/${LABEL}.bam METRICS_FILE=deduplicated/metrics.${LABEL}.txt TMP_DIR=temp ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

    """

}

process sam_index{

    module 'mugqic/samtools/1.14'
    cpus 1
    publishDir "."

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

    input:
    path input
  
    output:
    path("*.bed")
  
    """
    bamToBed -i $input > ${input}.bed

    """
}

process sort_bed{

    input:
    path input
    
    output:
    path("*.sort.bed")

    
    """
    temp=$input
    bamfile=\${temp%.*}
    echo \$bamfile
    sort -k 4 $input -S 2G > \$bamfile.sort.bed
    
    """
}

process fasta_faidx_index {

    module 'mugqic/samtools/1.14'
    publishDir "."

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

    input:
    path fasta
    path fasta_index
    path input_bed
    val iteration
        
    output:
    // path("scaffolds_i"+ iteration + "/*")
    tuple path("*"), val("scaffolds_i$iteration")

    """
    mkdir scaffolds_i${iteration}
    run_pipeline.py -a $fasta -l ${fasta_index} -b ${input_bed} -e GATC,GANTC -o scaffolds_i${iteration} -i $iteration -m yes
    """


}


////////////////
// steps for creating .hic file

process index_scaffolded_fasta {
    
    module 'mugqic/samtools/1.14'

    input:
    tuple path(index), val(iteration)

    output:
    tuple path("*/*.fai"), val(iteration)
    
    """
    FASTA=`find -L ./ -name "*FINAL.fasta"`
    samtools faidx \$FASTA -o \$FASTA.fai
    """
}

process make_chromosome_sizes {
    
    module 'mugqic/samtools/1.14'

    input:
    tuple path(index), val(iteration)

    output:
    tuple path("chromosome_sizes.tsv"), val(iteration)
    stdout

    """
    echo $index $iteration
    cut -f 1,2 $index | head -n 500 > chromosome_sizes.tsv


    """
}

