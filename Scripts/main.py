#! /usr/bin/python3

import gws_fetch, compare, annotate

def main(args)
    ###########################
    ###Filter and Group SNPs###
    ###########################
    args.out_fname="fetch.out"
    args.annotate_fpath=args.out_fname
    args.compare_fname=args.out_fname
    gws_fetch.fetch_gws(args)
    
    ###########################
    ##########Finemap##########
    ###########################

    ###########################
    #######Annotate SNPs#######
    ###########################
    annotate.main(args)
    ###########################
    ######Compare results######
    ###########################
    compare.compare(args)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="FINNGEN automatic hit reporting tool")
    #gws_fetch
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    #parser.add_argument("--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    #finemap
    
    #annotate
    #parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("--gnomad-path",dest="gnomad_path",type=str,help="Gnomad annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default=None,help="Finngen annotation file filepath")
    #parser.add_argument("--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    
    #compare results
    #parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use 'database' or 'file'")
    parser.add_argument("--summary-fpath",type=str,help="comparison summary filepath")
    
    
    args=parser.parse_args()
    main(args)