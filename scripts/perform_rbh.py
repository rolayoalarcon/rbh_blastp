import os
import time
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline


def do_blastp(file1, file2, out1, out2, outdir):
    
    fwd_blastp = NcbiblastpCommandline(query=file1, subject=file2, out=os.path.join(outdir, out1), outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue", max_target_seqs=1)
    rev_blastp = NcbiblastpCommandline(query=file2, subject=file1, out=os.path.join(outdir, out2), outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue", max_target_seqs=1)
    
    print("FORWARD: %s" % fwd_blastp)
    fwd_stdout, fwd_stderr = fwd_blastp()
    
    print("REVERSE: %s" % rev_blastp)
    rev_stdout, rev_stderr = rev_blastp()

    while not os.path.exists(os.path.join(outdir, out1)) or not os.path.exists(os.path.join(outdir, out2)):
        time.sleep(1) 
    
    return (os.path.join(outdir, out1), os.path.join(outdir, out2))

def combine_rbh(file1, file2, outname, outdir):
    fwd_results = pd.read_csv(file1, sep="\t", header=None, 
                              names=["query", "subject", "identity", "coverage",
                              "qlength", "slength", "alength",
                              "bitscore", "E-value"])
    
    rev_results = pd.read_csv(file2, sep="\t", header=None, 
                              names=["query", "subject", "identity", "coverage", 
                              "qlength", "slength", "alength", 
                              "bitscore", "E-value"])
    
    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                    left_on='subject', right_on='query',
                    how='outer')

    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]
    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()

    rbbh.to_csv(os.path.join(outdir, outname), sep='\t')

def main():

    cjejuni_pathogenex = "../../cjejuni_pathogenex/data/genome_info/GCF_000015525.1_ASM1552v1/translated_cds_loctag_headers.faa"
    cjejuni_nfcore = "../../cjejuni_nextflow/data/genome_info/translated_cds_loctag_headers.faa"
    salmonella_pathogenex = "../../salmonella_pathogenex/data/genome_info/GCF_000210855.2_ASM21085v2/translated_cds_loctag_headers.faa"
    ecoli_k12 = "../fastas/ecolik12_translated_cds_headers.faa"

    cpathogenex_vs_cnfcore = do_blastp(cjejuni_pathogenex, cjejuni_nfcore,
                                       "cpath_vs_cnfcore.tab", "cnfcore_vs_cpath,tab", "../blast_results/")
    
    combine_rbh(cpathogenex_vs_cnfcore[0], cpathogenex_vs_cnfcore[1], "cpath_vs_cnfcore.tsv", "../rbh_dataframes/")
    
    cpathogenex_vs_eck12 = do_blastp(cjejuni_pathogenex, ecoli_k12,
                                       "cpath_vs_eck12.tab", "eck12_vs_cpath,tab", "../blast_results/")
    combine_rbh(cpathogenex_vs_eck12[0], cpathogenex_vs_eck12[1], "cpath_vs_eck12.tsv", "../rbh_dataframes/")

    cpathogenex_vs_spathogenex = do_blastp(cjejuni_pathogenex, salmonella_pathogenex,
                                       "cpath_vs_spath.tab", "spath_vs_cpath,tab", "../blast_results/")
    combine_rbh(cpathogenex_vs_spathogenex[0], cpathogenex_vs_spathogenex[1], "cpath_vs_spath.tsv", "../rbh_dataframes/")
    
    spathogenex_vs_eck12 = do_blastp(salmonella_pathogenex, ecoli_k12, 
                                     "spath_vs_eck12.tab", "eck12_vs_spath.tab", "../blast_results/")
    combine_rbh(spathogenex_vs_eck12[0], spathogenex_vs_eck12[1], "spath_vs_eck12.tsv", "../rbh_dataframes/")

if __name__ == "__main__":
    main()