This is the Github repository for the scrips and input files used to perform the analysis described in the following project: "PhD thesis chapter 3- Ether lipid diversity in marine sponges"

## Requirements:

Please keep in mind that these scripts have only been tested in python3.

### Software:

- antiSMASH v5.0
- CORASON
- BiG-MAP
- dREP
- GTDB-Tk
- Pyhton 3+

### Pyhton Packages:

- argparse
- pandas
- numpy
- collections
- json
- SeqIO
- subprocess


## Guidelines:

### Extract and trim BGCs
vepe_trimming.xls: xls file with BGC, ORF_edge1 (starting ORF), ORF_edge2 (ending ORF)
```
python3 trim_genbank_orfnr.py -i vepe_trimming.xls -g /input_gbks/ -o /trimmed_gbks/
```


### Extract houskeeping genes from MAGs and calculate metatranscriptome expression metrics

extract coordinates of the HKG pfam hits and make bedfiles for each bin
```
python3 extract_hkgcds_coords.py -d /MAGs_hmmscan_out/ -p hkg_pfams.txt -o /hkg_bedfiles/ -n /MAGs_prodigal_out/
```

calculate metrics for hkgs and bgcs (avg_cov, min_cov, max_cov, percent_cov, tpm, rpkm)
```
python3 hkg_metat_metrics.py -b /hkg_sorted_count_files/ -l MAG_names.txt -p /hkg_bedgraphfiles/ -o hkg_metrics.tsv
```
```
python3 bgc_metat_metrics.py -c /bgc_sorted_count_files/ -b /bgc_bedgraphfiles/ -l BGC_headers.txt -o BGC_metrics.tsv
```

### Process BiG-MAP, dRep & GTDK-tk outputs to generate BGG/GCF annotations & iTOL datasets

merge BiG-MAP json files into a single all encompassing BGC -> GCF relationship
```
python3 merge_bigmap_jsons.py -m /BiG-MAP/BiG-MAP.GCs.json -b /BiG-MAP/BiG-MAP.GCF.json -j BiG-MAP.GC_GCF.json -t BiG-MAP.GC_GCF.tsv
```

add contig information
```
python3 gbk_extract_contigname.py -g /bgcs/ -t BiG-MAP.GC_GCF.tsv
# outfile: BiG-MAP.GC_GCF_contigs.tsv
```

add bin information
```
python3 from_contig_get_bin.py -i BiG-MAP.GC_GCF_contigs_sample.tsv -o BiG-MAP.GC_GCF_contigs_sample_bin.tsv -b /all_refined_bins/
```

compile bin clusters of dRep representative and member bins
```
python3 compile_bin_cluster_info.py -r dRep_representative_binnames.txt -c /dRep/data_tables/Cdb.csv -o rep_bins_samples.tsv
```

add representative bin and taxonomic classification info to  BiG-MAP.GC_GCF_contigs_sample_bin.tsv
```
python3 add_binrep_domexpinfo.py -b rep_bins_samples.tsv -t BiG-MAP.GC_GCF_contigs_sample_bin.tsv -o BiG-MAP.GC_GCF_contigs_sample_bin_binrep.tsv
```
```
python3 add_binrepclass_domexpinfo.py -c /gtdbtk/classify/gtdbtk.bac120.summary.tsv -t BiG-MAP.GC_GCF_contigs_sample_bin_binrep.tsv -o BiG-MAP.GC_GCF_contigs_sample_bin_binrep_class.tsv
```

build GCF oriented table with sample and bin info
```
python3 build_vepedomexp_repinfo_gtdbtk.py -i BiG-MAP.GC_GCF_contigs_sample_bin_binrep_class.tsv -j BiG-MAP.GC_GCF.json -o BiG-MAP.GC_GCF_reps_samples_binreps_class.tsv
```

make dRep representative bin abundance table, scaling values to meet a smaller range and improve visual aspect of iTOL dataset
```
python3 EL_binrep_quant_itol.py -q all_bin_gtdbtkquant.tsv -d rep_bins_samples.tsv -s rep_bin_quant_scaled.tsv -i rep_bin_quant.tsv
```

extract A-domain amino acid sequence for muscle alignment
```
python3 gbk_extract_adomain.py -i BGC_names.txt -g /bgcs_gbks/ -o extracted_Adomain.faa
```

