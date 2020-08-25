#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2019 Sergej Nowoshilow (sergej.nowoshilow@imp.ac.at)
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

params.genome = ''
params.outputdir = ''
params.reads = ''
params.pairmask = '/*_{left,right}*.fastq.gz'
params.chunkSize = 10_000_000
params.packChunks = false
params.adapters = ''
params.minQual5 = 10
params.minQual3 = 10
params.avgQual = 15
params.minLen = 50

def pairmask = (params.pairmask[0] == '/') ? params.pairmask : '/'+params.pairmask

def helpMessage() {
    log.info"""
    ===========================================================================================================================

    nf-merqury (v${workflow.manifest.version}) is a NextFlow-based genome analysis pipeline.
    ===========================================================================================================================

    Usage:

    Typical command for running the pipeline are as follows:

    nextflow run Pipeline.nf --genome <genome.fa> --outputroot <output directory> --reads <reads directory> --pairs-mask <pairs mask> 
                             --chunkSize <chunk size> --packChunks 
                             -resume

    Mandatory parameters:
      --genome            path to the genome FASTA file to correct
      --outputroot        path to the output directory
      --reads             path to the directory containing the read pairs

    Optional parameters:
      --pairmask          mask describing how to identify pairs in the reads directory. Default: ${params.pairmask}.
      --chunkSize         number of reads in one processing unit. Default is ${params.chunkSize}.
      --packChunks        if set, the chunks are gzipp'ed before mapping. Default is unset
      --adapters          FASTA file containing Illumina adapters to remove
      --minQual5          minimum required quality at the beginning of the read. Default is ${params.minQual5}.
      --minQual3          minimum required quality at the end of the read. Default is ${params.minQual3}.
      --avgQual           minimum required quality throughout the read. Default is ${params.avgQual}.
      --minLen            minimum required read length after trimming. Default is ${params.minLen}.

    """.stripIndent()
   return ""
}

// Check the parameters
if ( !params.genome ){
   exit 1, helpMessage() + "--genome not specified"
}

if ( !params.outputroot ){
   exit 2, helpMessage() + "--outputroot not specified"
}

if ( !params.reads ){
   exit 3, helpMessage() + "--reads not specified"
}

log.info"""
    ===========================================================================================================================

    nf-merqury (v${workflow.manifest.version}) is a NextFlow-based genome analysis pipeline.
    ===========================================================================================================================
""".stripIndent()

log.info "Running the pipeline with the following parameters: "
log.info "====================================================="
log.info "Genome path                : ${params.genome}"
log.info "Output directory           : ${params.outputroot}"
log.info "DNA-seq reads directory    : ${params.reads}"
log.info "Pair mask                  : ${params.pairmask}"
log.info "Chunk size                 : ${params.chunkSize}"
log.info "Pack chunks                : " + (params.packChunks ? 'Yes' : 'No')
log.info "Adapters                   : ${params.adapters}"
log.info "Leading quality            : ${params.minQual5}"
log.info "Trailing quality           : ${params.minQual3}"
log.info "Average quality            : ${params.avgQual}"
log.info "Minimum read length        : ${params.minLen}"

/********************************************************************************************************************************************
  STEP 1
  Identify read file pairs to use
********************************************************************************************************************************************/
Channel
 .fromFilePairs(params.reads + pairmask, flat: true)
 .into { ch_pairs_for_qc; ch_pairs_for_trimming }


/********************************************************************************************************************************************
  STEP 2
  Run FastQC for each file separately
********************************************************************************************************************************************/
ch_pairs_for_qc
  .flatMap {pair_id, fwd_reads, rev_reads -> [[pair_id, fwd_reads], [pair_id, rev_reads]]}
  .set { ch_fastqc_in }

process runFastQC {

  tag "${pair_id}"

  input:
    tuple val(pair_id), file(reads_file) from ch_fastqc_in

  output:
    file("*_fastqc.zip") into ch_fastqc_result

  script:
    """
    fastqc --outdir . --noextract --threads ${task.cpus} ${reads_file}
    """
}


/********************************************************************************************************************************************
  STEP 3
  Trim the reads
********************************************************************************************************************************************/
def iFragment = 1

ch_pairs_for_trimming
 .splitFastq(by: params.chunkSize, pe: true, compress: params.packChunks, file: true)
 .map {pair_id, fwd_reads, rev_reads -> [pair_id, iFragment++, fwd_reads, rev_reads] }
 .set { ch_read_chunks_trim }

process trimReads {

  tag "${pair_id} - ${fragment_id}"

  input:
    tuple val(pair_id), val(fragment_id), file(fwd_reads), file(rev_reads) from ch_read_chunks_trim

  output:
    tuple val(pair_id), 
          val(fragment_id), 
          file("${pair_id}.${fragment_id}_fwd.trimmed.fastq"), 
          file("${pair_id}.${fragment_id}_rev.trimmed.fastq") into ch_paired_reads
    tuple val(pair_id), 
          val(fragment_id),
          file("${pair_id}.${fragment_id}_fwd.unpaired.trimmed.fastq"),
          file("${pair_id}.${fragment_id}_rev.unpaired.trimmed.fastq") into ch_unpaired_reads
    tuple val(pair_id), file("${pair_id}.${fragment_id}.trimmingstats") into ch_trimming_stats
  
  script:
    """
    ILLUMINACLIP=
    if [ ! -z "${params.adapters}" ]
    then
      ILLUMINACLIP=ILLUMINACLIP:${params.adapters}:2:30:10
    fi
    
    trimmomatic \
      PE \
      -threads ${task.cpus} \
      ${fwd_reads} ${rev_reads} \
      -baseout ${pair_id}.fastq \
      \${ILLUMINACLIP} \
      LEADING:${params.minQual5} \
      TRAILING:${params.minQual3} \
      SLIDINGWINDOW:4:${params.avgQual} \
      MINLEN:${params.minLen} 2>${pair_id}.${fragment_id}.trimmingstats

    mv ${pair_id}_1P.fastq ${pair_id}.${fragment_id}_fwd.trimmed.fastq
    mv ${pair_id}_2P.fastq ${pair_id}.${fragment_id}_rev.trimmed.fastq
    mv ${pair_id}_1U.fastq ${pair_id}.${fragment_id}_fwd.unpaired.trimmed.fastq
    mv ${pair_id}_2U.fastq ${pair_id}.${fragment_id}_rev.unpaired.trimmed.fastq
    """
}


/********************************************************************************************************************************************
  STEP 4
  Count k-mers
********************************************************************************************************************************************/
ch_paired_reads
 .map {pair_id, frg_id, fwd_reads, rev_reads -> [pair_id, frg_id, fwd_reads, rev_reads, 'PE']}
 .set { ch_tmp_paired_trimmed }

ch_unpaired_reads
 .map {pair_id, frg_id, fwd_reads, rev_reads -> [pair_id, frg_id, fwd_reads, rev_reads, 'SE']}
 .set { ch_tmp_unpaired_trimmed }

ch_tmp_paired_trimmed
 .mix(ch_tmp_unpaired_trimmed)
 .set { ch_files_for_count }

process countKmers {
  tag "${pair_id} - ${fragment_id}"

  input:
    tuple val(pair_id), val(fragment_id), file(fwd_reads), file(rev_reads), val(mode) from ch_files_for_count

  output:
    tuple val(pair_id), val(fragment_id), val(mode), file("meryl.${pair_id}.${fragment_id}.${mode}.meryl") into ch_count_result

  script:
    """
    cat ${fwd_reads} ${rev_reads} > tmp.fastq
    MEM=\$(echo ${task.memory} | grep -Po "[0-9]+")
    meryl k=23 memory=\${MEM} threads=${task.cpus} count output meryl.${pair_id}.${fragment_id}.${mode}.meryl tmp.fastq
    rm tmp.fastq
    """
}


/********************************************************************************************************************************************
  STEP 5
  Merge the results
********************************************************************************************************************************************/
ch_count_result
 .collectFile { pair_id, fragment_id, mode, dir -> ["meryl.list", dir.collect{ '/' + it.toString() }.join('') + '\n' ] }
 .set { ch_merge_in }

process mergeMerylResults {

  input:
    file(list) from ch_merge_in

  output:
    file('merged.meryl') into ch_merge_out

  script:
  """
  MEM=\$(echo ${task.memory} | grep -Po "[0-9]+")

  # Meryl cannot handle too many files at once. Therefore, iteratively process 100 chunks at a time and merge the results afterwards
  ITERATION=1
  cp ${list} input_\${ITERATION}
  N_INPUT=\$(cat input_\${ITERATION} | wc -l)
  while [ \$N_INPUT -gt 100 ]
  do
    BLOCKDIR=iter_\${ITERATION}
    mkdir -p \${BLOCKDIR}

    split --numeric-suffixes=1 --lines=100 input_\${ITERATION} \${BLOCKDIR}/block
    BLOCK=1
    for FILE in \$(find \${BLOCKDIR} -name "block*")
    do
      meryl memory=\${MEM} threads=${task.cpus} union-sum output \${BLOCKDIR}/merged_\${BLOCK}.meryl \$(cat \$FILE)
      echo \${BLOCKDIR}/merged_\${BLOCK}.meryl >> tmp
      BLOCK=\$((BLOCK + 1))
    done
    ITERATION=\$((ITERATION + 1))
    mv tmp input_\${ITERATION}
    N_INPUT=\$(cat input_\${ITERATION} | wc -l)
  done

  meryl memory=\${MEM} threads=${task.cpus} union-sum output merged.meryl \$(cat input_\${ITERATION})
  """
}


/********************************************************************************************************************************************
  STEP 6
  Run merqury
********************************************************************************************************************************************/
process runQV {

  input:
    file(meryl) from ch_merge_out

  output:
    file('qv.out.*') into ch_qv_out

  script:
  """
  qv.sh ${meryl} ${params.genome} qv.out
  """
}

ch_qv_out.println()



/********************************************************************************************************************************************
  STEP 6
  Collect the FastQC result files and generate a single MultiQC report.
********************************************************************************************************************************************
ch_trimming_stats
  .groupTuple(by: 0)
  .collectFile { pair_id, files -> ["${pair_id}.trimmingstatslist", files.collect{ it.toString() }.join('\n') ] }
  .map { f -> [f.simpleName, f]}
  .set { ch_trimmingstats_list }

process mergeTrimmomaticLogs {

  tag "${pair_id}"

  input:
    tuple val(pair_id), file(statslist) from ch_trimmingstats_list

  output:
    file("${pair_id}.merged.trimmingstats") into ch_merged_trimmingstats

  script:
    """
    mergeTrimmomaticStats.py --logs ${statslist} > ${pair_id}.merged.trimmingstats
    """
}
*/