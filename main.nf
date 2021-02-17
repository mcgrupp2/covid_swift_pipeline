#!/usr/bin/env nextflow
/*
=======================================================================
                        TAYLOR
=======================================================================
 Trimming Amplications You LOve Rapidly
 #### Homepage / Documentation
https://github.com/greninger-lab/covid_swift_pipeline
-----------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2

// Print help message
def helpMessage() {
    log.info"""
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/covid_swift_pipeline -resume -with-docker ubuntu:18.04 --INPUT example/ --OUTDIR example/output/ --PRIMERS v2
    
    Parameters:
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --PRIMERS       Primer masterfile to run. By default, this pipeline uses the original Swift V1 primers. --PRIMERS V2 can be specified for the Swift V2 primer set.
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.
        --NO_CLIPPING   Skip primerclip option.
        --SGRNA_COUNT   Add extra step to count sgRNAs. 

        -with-docker ubuntu:18.04   [REQUIRED]
        -resume [RECOMMENDED]
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Initializing flags
params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false
params.PRIMERS = false
params.SGRNA_COUNT = false
params.NO_CLIPPING = false

// Checking for argument validity
// Throw error if --INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
    println("Make sure your input directory ends with trailing slash.")
   exit(1)
}
// Throw error if --OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   println("Make sure your output directory ends with trailing slash.")
   exit(1)
}
// Use Swift V2 masterfile if --PRIMERS v2 indicated
if (params.PRIMERS == "V2" | params.PRIMERS == "v2") {
    MASTERFILE = file("${baseDir}/sarscov2_v2_masterfile.txt")
    println("Using Swift V2 primerset...")
}
else {
    MASTERFILE = file("${baseDir}/sarscov2_masterfile.txt")
    println("Using Swift V1 primerset [default]...")
}
// Print when using SGRNA_COUNT
if (params.SGRNA_COUNT == false) {
    println("--SGRNA_COUNT not specified, skipping counting sgRNAs...")
} else {
    println("--SGRNA_COUNT specified. Will count sgRNAs...")
}

// Setting up files 
REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
REFERENCE_FASTA_FAI = file("${baseDir}/NC_045512.2.fasta.fai")
TRIM_ENDS=file("${baseDir}/trim_ends.py")
VCFUTILS=file("${baseDir}/vcfutils.pl")
SPLITCHR=file("${baseDir}/splitchr.txt")
ADAPTERS = file("${baseDir}/All_adapters.fa")
FIX_COVERAGE = file("${baseDir}/fix_coverage.py")
PROTEINS = file("${baseDir}/NC_045512_proteins.txt")
AT_REFGENE = file("${baseDir}/annotation/AT_refGene.txt")
AT_REFGENE_MRNA = file("${baseDir}/annotation/AT_refGeneMrna.fa")
LAVA_GFF = file("${baseDir}/annotation/lava_ref.gff")
MAT_PEPTIDES = file("${baseDir}/annotation/mat_peptides.txt")
MAT_PEPTIDE_ADDITION = file("${baseDir}/annotation/mat_peptide_addition.py")
RIBOSOMAL_START = file("${baseDir}/annotation/ribosomal_start.txt")
RIBOSOMAL_SLIPPAGE = file("${baseDir}/annotation/ribosomal_slippage.py")
PROTEINS = file("${baseDir}/annotation/proteins.csv")
CORRECT_AF = file("${baseDir}/annotation/correct_AF.py")
CORRECT_AF_BCFTOOLS = file("${baseDir}/annotation/correct_AF_bcftools.py")
SGRNAS = file("${baseDir}/sgRNAs_120bp.fasta")

// Import processes 
include Trimming from './Modules.nf'
include Fastqc from './Modules.nf'
include Aligning from './Modules.nf'
include Trimming_SE from './Modules.nf' 
include Fastqc_SE from './Modules.nf'
include CountSubgenomicRNAs from './Modules.nf'
include NameSorting from './Modules.nf'
include Clipping from './Modules.nf'
include BamSorting from './Modules.nf'
include GenerateConsensus from './Modules.nf'
include AnnotateVariants from './Modules.nf'

// Import reads depending on single end vs. paired end
if(params.SINGLE_END == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.INPUT}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    // Paired end first few steps
    if(params.SINGLE_END == false) {
        Trimming (
            input_read_ch, 
            ADAPTERS
        )
        Fastqc (
            Trimming.out[1]
        )
        Aligning (
            Trimming.out[0],
            REFERENCE_FASTA
        )
        // Optional step for counting sgRNAs 
        if (params.SGRNA_COUNT != false) {
            CountSubgenomicRNAs (
                Trimming.out[2],
                SGRNAS
            )
        }
    } else {
    // Single end first few steps
        Trimming_SE (
            input_read_ch,
            ADAPTERS
        )
        Fastqc_SE (
            Trimming_SE.out[1]
        )
        Aligning (
            Trimming_SE.out[0],
            REFERENCE_FASTA
        )
        // Optional step for counting sgRNAs 
        if (params.SGRNA_COUNT != false) {
            CountSubgenomicRNAs (
                Trimming_SE.out[2],
                SGRNAS
            )
        }
    }

    // Primerclip options for Swift runs
    if(params.NO_CLIPPING == false) {
        NameSorting (
            Aligning.out[0]
        )
        Clipping (
            NameSorting.out[0],
            MASTERFILE
        )
        GenerateConsensus (
            Clipping.out[0],
            REFERENCE_FASTA,
            TRIM_ENDS,
            FIX_COVERAGE,
            VCFUTILS,
            REFERENCE_FASTA_FAI,
            SPLITCHR
        )
    } else {
    // Skip primerclip for non-Swift runs
        BamSorting (
            Aligning.out[0]
        )
        GenerateConsensus (
            BamSorting.out[0],
            REFERENCE_FASTA,
            TRIM_ENDS,
            FIX_COVERAGE,
            VCFUTILS,
            REFERENCE_FASTA_FAI,
            SPLITCHR
        )
    }

    AnnotateVariants (
        GenerateConsensus.out[4],
        MAT_PEPTIDES,
        MAT_PEPTIDE_ADDITION,
        RIBOSOMAL_SLIPPAGE,
        RIBOSOMAL_START,
        PROTEINS,
        AT_REFGENE,
        AT_REFGENE_MRNA,
        CORRECT_AF_BCFTOOLS
    )
}