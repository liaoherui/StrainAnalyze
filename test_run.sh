# Mode-1 : Run with only reference genomes (give the target strain)
python StrainAnalyze.py -i /home/heruiliao2/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq -t GCF_016496025.1_ASM1649602v1_genomic.fna -r Ref_genome_select_Sau -o Out_mode1_give_tstrain

# Mode-1 Run with only reference genomes (no target strain)
python StrainAnalyze.py -i /home/heruiliao2/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq -r Ref_genome_select_Sau -o Out_mode1_no_tstrain

# Mode-2: Run with reference genomes and cds sequences (give the target strain)
python StrainAnalyze.py -i /home/heruiliao2/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq -t GCF_016496025.1_ASM1649602v1_genomic.fna -r Ref_genome_select_Sau -g Strain_cds -o Out_mode2_give_tstrain

# Mode-2: Run with reference genomes and cds sequences (no target strain)
python StrainAnalyze.py -i /home/heruiliao2/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq -r Ref_genome_select_Sau -g Strain_cds -o Out_mode2_no_tstrain

# Mode-3: Run with reference genomes and StrainScan output
python StrainAnalyze.py -i /home/heruiliao2/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq  -r /home/heruiliao2/Liangjun_test/Sau_ref_New_2022 -s ../StrainScan/Sau_test/z211_exRegion/final_report.txt -o Out_mode3

# Mode-4: Run with reference genomes, cds sequences, and StrainScan output
python StrainAnalyze.py -i /home/heruiliao3/Liangjun_test/Bac_HQ_data/MgPL201911674/z211/2_HQData/z211.fq  -r /home/heruiliao2/Liangjun_test/Sau_ref_New_2022 -g Strain_cds -s ../StrainScan/Sau_test/z211_exRegion/final_report.txt -o Out_mode4


