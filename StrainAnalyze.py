import re
import os
import sys
import argparse

__author__="Liao Herui - PhD of City University of Hong Kong"
usage="StrainAnalyze - A user-friendly bioinformatics pipeline used to analyze given identified bacterial strains in short reads."

def choose_mode(genes_seq,sc_out):
    mode=0
    if genes_seq=='' and sc_out=='':
        mode=1 # Only reference genomes are given
    if not genes_seq=='' and sc_out=='':
        mode=2
    if not sc_out=='' and genes_seq=='':
        mode=3
    if not sc_out=='' and not genes_seq=='':
        mode=4

    return mode

def select_target(ref_seq):
    d={} # filename -> genome seq
    for filename in os.listdir(ref_seq):
        f=open(ref_seq+'/'+filename,'r')
        d[filename]=0
        while True:
            line=f.readline().strip()
            if not line:break
            if re.search('>',line):continue
            d[filename]+=len(line)
    res=sorted(d.items(), key = lambda kv:(kv[1], kv[0]), reverse = True)
    return res_seq+'/'+res[0][0]

def build_dir(indir):
    if not os.path.exists(indir):
        os.makedirs(indir)

def pred_anno_genes(ref_seq,out_dir):
    print('StrainAnalyze:: Run Prodigal and Prokka for gene prediction and annotation...')
    # Predict genes using prodigal
    pred_genes=out_dir+'/Pred_genes'
    build_dir(pred_genes)
    for filename in os.listdir(ref_seq):
        pre=re.split('\.',filename)[0]
        print('prodigal -i '+ref_seq+'/'+filename+' -o '+pred_genes+'/'+pre+'.txt -d '+pred_genes+'/'+pre+'.fasta')
        os.system('prodigal -i '+ref_seq+'/'+filename+' -o '+pred_genes+'/'+pre+'.txt -d '+pred_genes+'/'+pre+'.fasta')
    # Annotate genes using prokka
    gene_list=[]
    for filename in os.listdir(pred_genes):
        if not re.search('fasta',filename):continue
        pre=re.split('\.',filename)[0]
        anno_genes=out_dir+'/Anno_genes/'+pre
        #build_dir(anno_genes)
        print('prokka --kingdom Bacteria --outdir '+anno_genes+' --prefix '+pre+'  --locustag '+pre+' '+pred_genes+'/'+filename)
        os.system('prokka --kingdom Bacteria --outdir '+anno_genes+' --prefix '+pre+' --locustag '+pre+' '+pred_genes+'/'+filename)
        gene_list.append(anno_genes+'/'+pre+'.fna')
    print('StrainAnalyze:: Gene prediction and annotation - Done.')
    return gene_list

def cluster_genes(gene_list,out_dir):
    print('StrainAnalyze:: Use CD-Hit to cluster genes from all input strains')
    cls_dir=out_dir+'/Cls_genes'
    build_dir(cls_dir)
    merge_cls=cls_dir+'/all_strains_genes.fasta'
    o=open(merge_cls,'w+')
    for g in gene_list:
        f=open(g,'r')
        pre=re.split('/',g)[-1]
        pre=re.split('\.',pre)[0]
        while True:
            line=f.readline().strip()
            if not line:break
            if re.search('>',line):
                line=re.sub('>','',line)
                o.write('>'+pre+'|'+line+'\n')
            else:
                o.write(line+'\n')
    o.close()
    cls_genes=cls_dir+'/reduced_strains_genes.fasta'
    print('cd-hit-est -g 1 -c 0.95 -T 0 -M 0 -d 1000 -i '+merge_cls+' -o '+cls_genes)
    os.system('cd-hit-est -g 1 -c 0.95 -T 0 -M 0 -d 1000 -i '+merge_cls+' -o '+cls_genes)
    f2=open(cls_genes+'.clstr','r')
    d={} # unique genes head name
    c=0
    while True:
        line=f2.readline()
        if not line:break
        
        if re.search('>Cluster',line):
            if not c==0:
                if c==1:
                    #print(name)
                    d[name]=''
                c=0
        else:
            c+=1
            name=line.split()[2]
            name=re.sub('\.\.\.','',name)
    if c==1:
        d[name]=''
    #exit()
    f3=open(cls_genes,'r') 
    # Search unique genes
    cls_genes_new=cls_dir+'/reduced_strains_genes_anno.fasta'
    o2=open(cls_genes_new,'w+')
    while True:
        line=f3.readline().strip()
        if not line:break
        if re.search('>',line):
            ele=line.split()
            if ele[0] in d:
                name=re.sub('>','',line)
                o2.write('>(Unique)|'+name+'\n')
            else:
                o2.write(line+'\n')
        else:
            o2.write(line+'\n')
    o2.close()
    print('StrainAnalyze:: Gene clustering - Done.')
    return cls_genes_new

def mapping_reads_analyze(fq_dir,tag_strain,gene_cls,out_dir):
    print('StrainAnalyze:: Map reads to ref with BWA and analyze using StrainGR...')
    map_dir=out_dir+'/Map_file'
    build_dir(map_dir)
    # BWA index
    os.system('bwa index '+gene_cls)
    os.system('bwa index '+tag_strain)
    # BWA map
    gbam=map_dir+'/gene_map.bam'
    sbam=map_dir+'/tagstrain_map.bam'
    os.system('bwa mem -I 300 -t 8 '+gene_cls+' '+fq_dir+' | samtools sort -@ 8 -O BAM -o '+map_dir+'/gene_map.bam -')
    os.system('bwa mem -I 300 -t 8 '+tag_strain+' '+fq_dir+' | samtools sort -@ 8 -O BAM -o '+map_dir+'/tagstrain_map.bam -')
    # Samtools index
    os.system('samtools index '+gbam)
    os.system('samtools index '+sbam)
    # Analyze alignment using StrainGR
    analyze_dir=out_dir+'/Analyze_res'
    build_dir(analyze_dir)

    os.system('straingr call '+gene_cls+' '+gbam+' --hdf5-out '+analyze_dir+'/gene_map.hdf5 --summary '+analyze_dir+'/gene_map.tsv --tracks all')
    os.system('straingr call '+tag_strain+' '+sbam+' --hdf5-out '+analyze_dir+'/tagstrain_map.hdf5 --summary '+analyze_dir+'/tagstrain_map.tsv --tracks all')

    os.system('straingr view '+analyze_dir+'/gene_map.hdf5 -V '+analyze_dir+'/gene_map.vcf')
    os.system('straingr view '+analyze_dir+'/tagstrain_map.hdf5 -V '+analyze_dir+'/tagstrain_map.vcf')
    print('StrainAnalyze:: All analysis are finished!')

def search_strainscan_res(sc_out,ref_seq,out_dir):
    f=open(sc_out,'r')
    strain_pre_list=[]
    line=f.readline()
    new_ref_seq=out_dir+'/Ref_seq'
    build_dir(new_ref_seq)
    while True:
        line=f.readline().strip()
        if not line:break
        ele=line.split('\t')
        st=ele[1].split()[0]
        strain_pre_list.append(st)
    bfix=['fna','fa','fasta']
    tag_strain=''
    
    for filename in os.listdir(ref_seq):
        pre=re.split('\.',filename)[0]
        bre=re.split('\.',filename)[-1]
        if bre not in bfix:continue
        '''
        if pre == strain_pre_list[0]:
            tag_strain=ref_seq+'/'+filename
        '''
        if pre in strain_pre_list:
            os.system('cp '+ref_seq+'/'+filename+' '+new_ref_seq)
        if pre == strain_pre_list[0]:
            tag_strain=new_ref_seq+'/'+filename

    return tag_strain,strain_pre_list,new_ref_seq

def choose_genes(plist,genes_seq):
    gene_list=[]
    for filename in os.listdir(genes_seq):
        pre=re.split('\.',filename)[0]
        if pre in plist:
            gene_list.append(genes_seq+'/'+filename)
    return gene_list



def analyze(fq_dir,ref_seq,tag_strain,genes_seq,sc_out,out_dir,mode):
    if mode==1:
        if tag_strain=='':
            tag_strain=select_target(ref_seq)
        if not re.search('/',tag_strain):
            tag_strain=ref_seq+'/'+tag_strain

        gene_list=pred_anno_genes(ref_seq,out_dir)
        if len(gene_list)==1:
            gene_cls=gene_list[0]
        else:
            gene_cls=cluster_genes(gene_list,out_dir)
        mapping_reads_analyze(fq_dir,tag_strain,gene_cls,out_dir)    

    if mode==2:
        if tag_strain=='':
            tag_strain=select_target(ref_seq)
        if not re.search('/',tag_strain):
            tag_strain=ref_seq+'/'+tag_strain

        gene_list=[]
        for filename in os.listdir(genes_seq):
            gene_list.append(genes_seq+'/'+filename)
        if len(gene_list)==1:
            gene_cls=gene_list[0]
        else:
            gene_cls=cluster_genes(gene_list,out_dir)
        mapping_reads_analyze(fq_dir,tag_strain,gene_cls,out_dir)

    if mode==3:
        if tag_strain=='':
            tag_strain,plist,new_ref_seq=search_strainscan_res(sc_out,ref_seq,out_dir)
        gene_list=pred_anno_genes(new_ref_seq,out_dir)
        if len(gene_list)==1:
            gene_cls=gene_list[0]
        else:
            gene_cls=cluster_genes(gene_list,out_dir)
        mapping_reads_analyze(fq_dir,tag_strain,gene_cls,out_dir)



    if mode==4:
        if tag_strain=='':
            tag_strain,plist,new_ref_seq=search_strainscan_res(sc_out,ref_seq,out_dir)
        gene_list=choose_genes(plist,genes_seq)
        #print(gene_list)
        if len(gene_list)==1:
            gene_cls=gene_list[0]
        else:
            gene_cls=cluster_genes(gene_list,out_dir)
        #exit()
        mapping_reads_analyze(fq_dir,tag_strain,gene_cls,out_dir)
         


def main():
    parser=argparse.ArgumentParser(prog='StrainAnalyze.py',description=usage)
    parser.add_argument('-i','--input_fastq',dest='input_fq',type=str,required=True,help="The dir of input fastq data --- Required")
    parser.add_argument('-r','--identified_strains',dest='iden_strains',type=str,required=True,help="The dir of reference genomes of identified strains --- Required")
    parser.add_argument('-t','--target_strain',dest='tag_strain',type=str,help="The filename of the targeted strain. If neither this parameter nor \"-s\" parameters are set, then will choose the longest genome as the targeted strain by default and map reads to the targeted strain to call SNVs.")
    parser.add_argument('-g','--strains_genes',dest='st_genes',type=str,help="The dir of CDS sequences of identified strains. If this parameter is not provided, then will predict genes from reference genomes using prodigal and annotate predicted genes using prokka.")
    parser.add_argument('-s','--strainscan_output',dest='sc_output',type=str,help="The dir of the StrainScan output file. If this parameter is given, then will take the dominant strain as the targeted strain. All possible strains and strains with extraRegion covered will be analyzed.")
    parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir (default: current dir/StrainAnalyze_Result)')
    

    args=parser.parse_args()
    
    pwd=os.getcwd()
    fq_dir=args.input_fq
    ref_seq=args.iden_strains
    tag_strain=args.tag_strain
    genes_seq=args.st_genes
    sc_out=args.sc_output
    out_dir=args.out_dir

    if not genes_seq:
        genes_seq=''
    if not sc_out:
        sc_out=''
    if not tag_strain:
        tag_strain=''
    if not out_dir:
        out_dir=pwd+'/StrainAnalyze_Result'
    
    if not os.path.exists(out_dir) :
        os.makedirs(out_dir)

    #if genes_seq=='' and sc_out=='':

    mode=choose_mode(genes_seq,sc_out)
    #print(mode,genes_seq,sc_out)

    analyze(fq_dir,ref_seq,tag_strain,genes_seq,sc_out,out_dir,mode)
    

if __name__=="__main__":
    sys.exit(main())
