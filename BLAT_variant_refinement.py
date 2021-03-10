#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import datetime
import os,sys,re
import ntpath
import json
import glob
import urllib.request
import logging
import subprocess
import multiprocessing as mp
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)


def blat_variant_refinement(vcf_infile, outdir, genome_fasta_file, threads, bamFile, logger):

    logger.info("\n Running BLAT for realignment to remap reads containing variants.")


    ##############
    # Constants #
    ##############
    
    Phredscore_endcoding_offset = 33
    minimum_base_quality = 25
    minimum_mismatch = 1
    score_limit = 0.95
    
    # File Paths input and output 
    outfile_passed_path = "{}/passed_variants.vcf".format(outdir)
    outfile_failed_path = "{}/failed_variants.vcf".format(outdir)

    outprefix = "{}/BlatVarRefine".format(outdir)
    fa_file_path = "{}.fa".format(outprefix)
    psl_file_path = "{}.psl".format(outprefix)
    TEMP = outprefix + '_tmp'

    if os.path.exists(fa_file_path):
        logger.warn("-reusing existing reads file: {}".format(fa_file_path))
    else:
        write_reads_fasta(vcf_infile, fa_file_path, TEMP, logger)
        

    if os.path.exists(psl_file_path):
        logger.warn("-reusing existing blat psl file: {}".format(psl_file_path))
    else:
        run_pblat(fa_file_path, genome_fasta_file, psl_file_path, threads, logger)
        

    
    #----------------------------------------
    # Process the PBLAT output
    #----------------------------------------
    message = "Process the PBLAT output"
    logger.info(message)

    # separate the file by variant-supporting read 
    psl_dict = {}
    psl = open(psl_file_path, "r")
    for i in psl:
        line = i.strip().split("\t")

        if len(line) < 21:
            logger.warn("-psl line has fewer fields than expected: {}. Skipping.".format(i))
            continue
        
        name = line[9]

        line_list = [line[0], # alignment score
                     line[13], # chromosome
                     line[17], # 
                     line[18],
                     line[20]]
        # add to dictionary 
        psl_dict.setdefault(name, []).append(line_list)
    psl.close()



    logger.info("Examining multimappings")
    
    keep_dict = {}
    filter_dict = {}

    # Key format example 
    #   'chr1-14542-1-readname'

    # track the reads that fail this check.
    failed_reads_filename = outprefix + ".failed_reads"
    failed_reads_ofh = open(failed_reads_filename, 'w')

    passed_reads_filename = outprefix + ".passed_reads"
    passed_reads_ofh = open(passed_reads_filename, 'w')

    
    for i in psl_dict:
        tmp = i.split("^")
        current_chromosome = tmp[0]
        current_position = int(tmp[1])
        psl_id = tmp[0] + "_" + tmp[1]  #chrpos

        read_name = tmp[3]
        
        # set the max line to the first line 
        max_scored_line = psl_dict[i][0]
        max_score = max_scored_line[0]
        list_scores = []
        # loop over the multiple mappings for read i
        for j in psl_dict[i]:
            if int(j[0]) > int(max_score):
                max_scored_line = j
                max_score = int(j[0])
            list_scores.append(int(j[0]))

        #-------------------------------------------
        # Check to see if the reads are overlapping 
        # make sure the second best hit (based on score) has a score thats is < 95% of the best hit 
        #-------------------------------------------
        overlaping = 0
        # if more than one read was found
        if len(list_scores) == 1:
            list_scores.append(0) # handles the unique mapping scenario

        list_scores.sort(reverse = True) 
        # check chromosomes to make sure match 
        if max_scored_line[1] == current_chromosome:
            
            # check second best hit score 
            if list_scores[1] < (list_scores[0]*score_limit):

                # have multi-mapping read.
                # ensure that the best match corresponds to our variant location in the genome:
                block_count, block_sizes, block_starts = max_scored_line[2], max_scored_line[3].split(","), max_scored_line[4].split(",")
                for k in range(int(block_count)):
                    start_position = int(block_starts[k])+1
                    end_position = int(block_starts[k]) + int(block_sizes[k])
                    # check if positions are overlapping 
                    if (current_position >= start_position) and (current_position <= end_position):
                        overlaping = 1
                if overlaping:
                    # track the number of reads with best matches hitting the expected variant position according to blat. 
                    if (psl_id in keep_dict): 
                        keep_dict[psl_id] += 1
                    else:
                        keep_dict[psl_id] = 1
                    passed_reads_ofh.write(read_name + "\n")
                    
        if overlaping == 0:
            # multiple high-scoring hits
            # track count in filter_dict
            if (psl_id in filter_dict): 
                filter_dict[psl_id] += 1
            else:
                filter_dict[psl_id] = 1
            failed_reads_ofh.write(read_name + "\n")

    failed_reads_ofh.close()
    
    ######################
    ## re-process the VCF, partitioning into pass vs. fail sites.

    logger.info("-reprocessing VCF, partitioning pass vs. fail sites")
    
    outfile_passed = open(outfile_passed_path, "w")
    outfile_failed = open(outfile_failed_path, "w")

    
    # counters to see how many variants pass and fail 
    passed, failed = 0,0
        
    infile = open(vcf_infile, "r")
    for i in infile:
        if i[0] == "#" : continue
        line = i.strip().split("\t")
        
        psl_id = line[0] + "_" + line[1]

        if psl_id in keep_dict:
            new_alter = int(keep_dict[psl_id])

            if psl_id in filter_dict:
                discard_counter = filter_dict[psl_id]
            else:
                discard_counter = 0

            
            if (new_alter >= minimum_mismatch) and (new_alter > discard_counter):
                # have at least one pass mismatch read and number passed exceeds number discarded.
                
                new_output_line =line[0] + "\t" + line[1] + "\t" + str(new_alter) + "\t" + line[3] + "\t" + line[4] + "\n"
                outfile_passed.write(new_output_line)
                passed += 1
            else:
                failed_output_line = i + "\tfailed_freq #mismatches: " + str(new_alter) + "\tminimumMismatchesNecessary: " + str(minimum_mismatch) + "\tdiscarded mismatch reads: " + str(discard_counter) + "\n"
                outfile_failed.write(failed_output_line)
                failed += 1
        else:
            failed_output_line = i[:-2] + "\ttotalcover_failed\n"
            outfile_failed.write(failed_output_line)
            failed += 1

    outfile_passed.close()
    outfile_failed.close()
    
    print("Variants kept:", passed)
    print("Variants filtered:", failed)

    return

        

def write_reads_fasta(vcf_infile, fa_file_path, TEMP, logger):

    infile = open(vcf_infile, "r")
    fa_file = open(fa_file_path, "w")

    ########
    # 1) Prep file for BLAT
    #       use samtools view to get the alignments and get reads containing variants 
    # 2) Use BLAT
    #       use pblat (parallel BLAT) for realignment of the reads that contain the variants
    ########

    ####################################
    # Prep file for BLAT
    ####################################
    # get reads containing variants 
    # Loop over each variant 

    for i in infile:
        #-----------------------------------------
        # Get the reads that contain the variants 
        # samtools view: prints alignments that are in the given input alignment file
        #-----------------------------------------
        if i[0] == "#":
            continue
        
        line = i.split("\t")
        chrom, position = line[0], line[1]
        bamposition = chrom + ':' + position+'-'+position

        cmd = "samtools view {} {} > {}".format(bamFile, bamposition, TEMP)
        #subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()
        subprocess.check_call(cmd, shell=True)
        
        editnuc = line[4]
        newmismatch = 0
        mismatchreadcount = 0
        newcov, newmismatch = 0, 0 

        sam_alignments = open(TEMP,"r")

        #----------------------------------------
        # read through the input file 
        #----------------------------------------
        for sam_align in sam_alignments:
            bam_fields = sam_align.strip().split("\t")
            (read_name,
             alignment,
             readstart,
             cigar,
             sequence_bases,
             quality_scores) = (bam_fields[0],
                                bam_fields[1],
                                bam_fields[3],
                                bam_fields[5],
                                bam_fields[9],
                                bam_fields[10])
            

            m = re.search("NH:i:\d+", sam_align)
            if not m:
                print("Error, sam line {} is missing NH:i:number hit count".format(sam_align), file=sys.stderr)
                sys.exit(1)
            NH_tag = m.group(0)

            current_pos, readpos = int(readstart), 1
            base_readpos = False
            
            # leters 
            cigarletters = re.findall(r'[MIDNSHP=X]',cigar)
            # numbers
            cigarnums = re.split('[MIDNSHP]', cigar)


            for k in range(len(cigarletters)):

                #### it is now faster ####
                position = int(position)
                if current_pos > position:
                    break
                ##### error corrected ####
                if cigarletters[k] == "S" or cigarletters[k] == "I":
                # if cigarletters[i] =~ m/[SI]/) {
                    readpos = readpos + int(cigarnums[k])

                elif cigarletters[k] == "D" or cigarletters[k] == "N": 
                    current_pos = current_pos + int(cigarnums[k])
                
                elif cigarletters[k] == "M":
                    for j in range(int(cigarnums[k])):
                        if current_pos == position:
                            ## at variant site position in genome and read
                            if sequence_bases[readpos-1] == editnuc:
                                ## have a read containing the variant base
                                if ord(quality_scores[readpos-1]) >= minimum_base_quality + Phredscore_endcoding_offset:
                                    ## consider a quality base call at variant position.
                                    base_readpos = True
                                    break
                        current_pos += 1
                        readpos += 1

                        
            if base_readpos:
                # increment miss matched read count by one 
                mismatchreadcount+=1
                # write the line to the FA file 
                fa_line = "^".join([">" + chrom,
                                    str(position),
                                    str(mismatchreadcount),
                                    read_name + "," + NH_tag]) + "\n" + sequence_bases + "\n"
                # print(fa_line)
                fa_file.write(fa_line)

    fa_file.close() # fa file containing mismatch-containing reads
    infile.close() # input vcf file
    
    return



def run_pblat(fa_file_path, genome_fasta_file, psl_file_path, threads, logger):
        
    
    #----------------------------------------
    # RUN PBLAT: 
    #       use pblat (parallel BLAT) for realignment of the reads that contain the variants 
    #----------------------------------------

    message = "Runnning PBLAT to remap variant reads. Threads: {}".format(threads)
    logger.info(message)

    # Path to blat/pblat
    cmd_pblat = "pblat \
                 -threads={} \
                 -stepSize=5 \
                 -repMatch=2253 \
                 -minScore=20 \
                 -minIdentity=0 \
                 -noHead {} {} {}".format(threads, genome_fasta_file, fa_file_path, psl_file_path)
    print(cmd_pblat)
    subprocess.Popen(cmd_pblat, shell=True).communicate()

    return


def make_menu():
    ###########################################
    # Gather the arguments passed to the SNPiR script in comand line 
    ###########################################

    ## Input Arguments
    # Description 
    args_parser = argparse.ArgumentParser(
        description = "Run SNPiR steps."
        )
    #---------------------
    ## Required arguments
    #---------------------
    required = args_parser.add_argument_group('required arguments')

    required.add_argument( "-I", metavar = "Input_File", dest = "str_input_file",
                           required = True, help = "The input VCF file to run SNPiR on.")

    required.add_argument( "-B", metavar = "Input_BAM_File", dest = "str_BAM_file",
                           required = True, help = "The input BAM file.")

    required.add_argument("-G", metavar="genome_fasta_file", dest="str_genome_fasta_file",
                          required=True, help="genome fasta file")
    
    required.add_argument( "-O", metavar = "Output_Directory", dest = "str_out_dir",
                           required = False, help = "Where to put the results from SNPiR.",
                           default = ".")
    #---------------------
    ## Optional arguments
    #---------------------
    optional = args_parser.add_argument_group('optional arguments')

    optional.add_argument("--threads", metavar = "Process_threads", dest = "i_number_threads",
                          type = int, default = 1, help = "The number of threads to use for multi-threaded steps.")
    optional.add_argument("--debug", action="store_true", help="sets debug mode for logger")

    return args_parser

if __name__ == "__main__":

    # Parse the arguments given to SNPiR
    args_parser = make_menu()
    args_parsed = args_parser.parse_args()
    
    logger = logging.getLogger(__name__)
    if args_parsed.debug:
        logger.setLevel(logging.DEBUG)
    
    ###################
    # Set Constants
    ###################

    infile = args_parsed.str_input_file
    outdir = args_parsed.str_out_dir
    refined_bam = args_parsed.str_BAM_file
    genome_fasta_file = args_parsed.str_genome_fasta_file
    
    # Get the output directory and create it if it does not exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)


    blat_variant_refinement(infile, outdir, genome_fasta_file, threads=args_parsed.i_number_threads, bamFile = refined_bam, logger=logger)
    

    sys.exit(0)
