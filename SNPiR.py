#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, re
import datetime
import subprocess
import argparse
import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger('ctat_mutations')
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../PyLib"]))
from Pipeliner import Pipeliner, Command, run_cmd, ParallelCommandList


SNPIR_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "util"])
CTAT_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../src"])


def main():
    parser = argparse.ArgumentParser(description="SNPiR filtering applied to annotated vcf", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument("--work_dir", type=str, required=True, help="working directory name for intermediates")

    parser.add_argument("--output_filename", required=True, help="name of output file containing final list of variants")

    args = parser.parse_args()
     
    
    ## prep for run
    output_dir = args.work_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpts_dir = os.path.join(output_dir, "__chckpts")

    pipeliner = Pipeliner(checkpts_dir)

    attribute_listing = ",".join(['RS', 'QUAL', 'TPR', 'RPT', 'DJ', 'Homopolymer', 'ED', 'RNAEDIT', 'LEN'])
    
    # extract feature matrix
    matrix_file = os.path.join(output_dir, "features.matrix")
    
    cmd = " ".join([ os.path.join(CTAT_UTILDIR, "annotated_vcf_to_feature_matrix.py"),
                    "--vcf", args.input_vcf,
                    "--features", attribute_listing,
                    "--output", matrix_file])
    pipeliner.add_commands([Command(cmd, "feature_extraction_to_matrix.ok")])



    snpir_selected_matrix_file = matrix_file + ".SNPiR_selected"
    
    cmd = " ".join([ os.path.join(SNPIR_UTILDIR, "SNPiR_select_variants.Rscript"),
                    "--input_matrix", matrix_file,
                    "--output_matrix", snpir_selected_matrix_file])

    pipeliner.add_commands([Command(cmd, "snpir_selected_matrix.ok")])


    cmd = " ".join([ os.path.join(CTAT_UTILDIR, "extract_boosted_vcf.py"),
                    "--vcf_in", args.input_vcf,
                    "--boosted_variants_matrix", snpir_selected_matrix_file,
                    "--vcf_out", args.output_filename])

    pipeliner.add_commands([Command(cmd, "snpir_filt_vcf.ok")])
    
    pipeliner.run()


if __name__=='__main__':
    main()

