#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 16:58:20 2012

@author: Marc Tormo
"""

import argparse
import csv
import datetime
import glob
import multiprocessing
import os
import subprocess
import sys
import vcf


#canviar rutes!!!!!!!!!!!!!!!!!
home='/home/sequentia1'
home_soft=os.path.join(home,'opt')
# pindel_t=os.path.join(home_soft,'pindel024t/pindel')#old version (024t)
# pindel2vcf_t=os.path.join(home_soft,'pindel024t/pindel2vcf')#old version (024t)
# pindel_v=os.path.join(home_soft,'pindel/trunk/pindel')#new version (024v)
# pindel2vcf_v=os.path.join(home_soft,'pindel/trunk/pindel2vcf')#new version (024v)
bwa=os.path.join(home_soft,'bwa/bwa-0.7.5a/bwa')

pindel_v='pindel'
pindel_t='pindel'
pindel2vcf_v='pindel2vcf'
pindel2vcf_t='pindel2vcf'

bgzip='bgzip'
bcftools='bcftools'
picard='java -jar /usr/local/bin/picard-tools-1.112/AddOrReplaceReadGroups.jar'
samtools='samtools'
tabix='tabix'
vcf_annotate='vcf-annotate'
gatk_ref='java -jar /usr/local/bin/picard-tools-1.112/CreateSequenceDictionary.jar'
gatk_target='java -Xmx4g -jar /home/sequentia1/opt/GenomeAnalysisTK.jar'


##local!!!!!!!!!!!!!
#samtools='samtools'
#bwa=os.path.join(home_soft,'bwa-0.7.5a/bwa')
#bcftools='bcftools'
#picard=' '.join(['java','-jar',os.path.join(home_soft,'picard-tools-1.56/AddOrReplaceReadGroups.jar')])
#pindel_t=os.path.join(home_soft,'pindel024t/pindel')#old version (024t)
#pindel2vcf_t=os.path.join(home_soft,'pindel024t/pindel2vcf')#old version (024t)
#pindel_v=os.path.join(home_soft,'pindel/trunk/pindel')#new version (024v)
#pindel2vcf_v=os.path.join(home_soft,'pindel/trunk/pindel2vcf')#new version (024v)
#bgzip='bgzip'
#tabix='tabix'
#vcf_annotate='vcf-annotate'


def InputPar():
    ####Introducing arguments
    parser = argparse.ArgumentParser(description='S.U.P.E.R. (Simply Unified Paired-End Resequencing Pipeline)')
    parser.add_argument('-p','--path', metavar='PATH', required=True, 
                        help='Path with files to analyse (Required)')
    parser.add_argument('-g','--genome',  metavar='GENOME', 
                        help='Genome file with path (Required for fastq and bam files)')
    parser.add_argument('-b','--bwa_dir',  metavar='./PATH', default='BWA_output',
                        help='Name of Bwa output directory without entire path (default BWA_output/)')
    parser.add_argument('-E','--ends', metavar='s/p', default='p', choices='sp' , 
                        help='Single or Paired-Ends?[p]')
    parser.add_argument('-c','--pe_distance', metavar='INT', default='500',
                        help='Paired-end distance [500]')
    parser.add_argument('-s','--analysis', metavar='INT', default='1' , choices='123',
                        help='Analyse 1.SNPs/INDELs or 2.SVs or 3.Both [1]')
    parser.add_argument('-f','--flank', metavar='INT', default=1,type=int, 
                        help='Flanking region for SNPs [1]')
   # parser.add_argument('-a','--bwa_aln', metavar='-a"OPTS BWA ALN"', default=' ', const=str, nargs="?", 
   #                    help='Options for alignment in bwa aln, in quotes, without spaces -a" " [None]')
   #parser.add_argument('-e','--bwa_sam', metavar='-e"OPTS BWA SAMPE/SAMSE"', default=' ', const=str, nargs="?", 
   #                    help='Options for alignment in bwa sampe/samse, in quotes, without spaces -e" " [None]')
    parser.add_argument('-q','--qsam', metavar='INT', default='10', 
                        help='Quality for samtools mpileup [10]')
    parser.add_argument('-R','--rem_dup', metavar='s/S', default=' ', choices='sS' , 
                        help='Option for samtools rmdup (s for single ends) [None]')
    parser.add_argument('-d','--pindel_config', metavar='-d"OPTS PINDEL CONFIG"', default=' ', const=str, nargs="?", 
                        help='Options for configuration in Pindel, in quotes, without spaces -d" " [None]')
    parser.add_argument('-D','--pindel2vcf_config', metavar='-D"OPTS PINDEL2VCF CONFIG"', default=' ', const=str, nargs="?", 
                        help='Options for pindel2vcf, in quotes, without spaces -D" " [None]')
    parser.add_argument('-n','--ann_file', metavar='ANNOTATION', 
                        help='Annotation File Decompressed (CHR-ANN-FROM-TO)')
    parser.add_argument('-T','--type_geno', metavar='c/s', default='c', choices='cs', 
                        help='Type of genome file information: c (chromosome) or s (scaffold) [c]')
    parser.add_argument('-N','--threads', metavar='INT', default='8', 
                        help='Number of threads to run [8]')

    
    return parser.parse_args()

def CheckInput(args):
    ###check extension of fastq files
    if glob.glob(args.path+'/*.f*q'):
        files,n_file=FilesExt(args.path,'.f*q')
        fastq=os.path.splitext(files[0])[1]
    elif glob.glob(args.path+'/*.f*q.gz'):
        files,n_file=FilesExt(args.path,'.f*q.gz')
        root1,ext = os.path.splitext(files[0])
        root2, first_ext = os.path.splitext(root1)
        ext = first_ext + ext
        fastq=ext
    else:
        fastq='.fastq'

    if not glob.glob(args.path+'/*.vcf') and not args.genome:
        sys.exit(' '.join(['!!!You need to add a genome file!!!','\nExiting process']))
    if glob.glob(args.path+'/*.vcf') and int(args.analysis)==3:
        sys.exit(' '.join(['!!! Can not analyse SNPs and SVs with the same vcf !!!','\nExiting process']))
    if args.ends=='p' and args.rem_dup=='s':
        sys.exit(' '.join(['!!! Error in parameter --rem_dup (s is only for single ends) !!!','\nExiting process']))

    if args.ends=='s':
        args.rem_dup='-s'#rmdup for samtools
    elif args.ends=='p' and args.rem_dup==' ':
        args.rem_dup=' '#rmdup for samtools
    elif args.ends=='p' and args.rem_dup=='S':
        args.rem_dup='-S' #rmdup for samtools

    return fastq,args.rem_dup

def PrintLog(command):
    ###Print all commands in a log file
    cl = open('Commands.log','a')
    cl.write(command)
    cl.write('\n')
    cl.close()

###Wait until all parallel processes are finished (python)
##########################################################
###list_script=list of processes
def parallel(list_script):
    ps = {}
    for script in list_script:
        p = subprocess.Popen(script,stdout=subprocess.PIPE,stdin=subprocess.PIPE,shell=True)
        ps[p.pid] = p
        PrintLog(script)

    print "Waiting for %d processes..." % len(ps)
    while ps:
        pid, status = os.wait()
        if pid in ps:
            del ps[pid]
            print "Waiting for %d processes..." % len(ps)        

def FilesExt(path,ext):#from an extension, returns all files in a list and nº of files
    filename=''.join(['*',ext])
    files = glob.glob(os.path.join(path, filename))

    return files,len(files)

def CreateIndex(genome):
    ###Create bwa & samtool index if they don't exist
    genome_bwt = '.'.join([genome,'bwt'])
    genome_fai = '.'.join([genome,'fai'])
    genome_dict = '.'.join([genome,'dict'])
    #list_for_dict_index=[]
    if not os.path.exists(genome_bwt) and not os.path.exists(genome_fai) and not os.path.exists(genome_dict) :
        ##bwa
        os.system(' '.join([bwa ,'index -a bwtsw', genome]))
        PrintLog(' '.join([bwa ,'index -a bwtsw', genome]))
        ##samtool
        os.system(' '.join([samtools, 'faidx', genome]))
        PrintLog(' '.join([samtools, 'faidx', genome]))
        ##index for GATK
        os.system(' '.join([gatk_ref, '='.join(['R',genome]), '='.join(['O',os.path.join([genome,'dict'])])]))
        

def BWAMap(param,name_trim,bam_dir,no_of_files,fastq):#align and map with BWA
    #### Run BWA for single/paired-ends
    ###################################
    ###join all commands in a list
    sam_list=[]
    sam_file=[]
    ###Single-Ends
    if param.ends == 's':

        for in_trim in name_trim:
            sam_list.append(" ".join([bwa, 'mem', '-t', param.threads, param.bwa_sam, param.genome, os.path.join(param.path,''.join([in_trim,fastq])), '>', os.path.join(bam_dir,'.'.join([in_trim,'sam']))]))
            sam_file.append(in_trim)
    
        ###run & wait until all parallel processes are finished
        parallel(sam_list)
    ###Paired-Ends
    elif param.ends == 'p':
        for k in xrange(no_of_files):
            ###take only pairs of reads: 0-1,2-3,4-5,...
            if not k % 2 == 0:
                continue
            sam_list.append(" ".join([bwa, 'mem', '-t', param.threads, param.bwa_sam, param.genome, os.path.join(param.path,''.join([name_trim[k],fastq])),os.path.join(param.path,''.join([name_trim[k+1],fastq])), '>', os.path.join(bam_dir,'.'.join([name_trim[k],'sam']))]))
            sam_file.append(name_trim[k])
        ###run & wait until all parallel processes are finished
        parallel(sam_list)

    ###extract an identifier
    iden=[]        
    for item in sam_file:
        iden.append(item.partition('.')[0])
    
    return sam_file,iden

def SamToBam(bam_dir,sam_file,rem_dup,typ,iden,genome):#convert sam to bam, in snps (uniques reads) or svs
    map_list2=[]
    map_list3=[]
    map_list4=[]
    map_list5=[]
    map_list6=[]
    map_list7=[]
    map_list8=[]
    map_list9=[]
    map_list10=[]
    for n in xrange(len(sam_file)):
        ###convert sam to bam
        map_list2.append(' '.join([samtools ,'view', '-Sb', os.path.join(bam_dir,'.'.join([sam_file[n],'sam'])), '>', os.path.join(bam_dir, '.'.join([sam_file[n],'bam']))]))
        ###sort bam
        map_list3.append(' '.join([samtools, 'sort', os.path.join(bam_dir, '.'.join([sam_file[n],'bam'])), os.path.join(bam_dir, '.'.join([sam_file[n],'sort']))]))
        ###do index (bai)
        map_list4.append(' '.join([samtools, 'index', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.bam']))]))
        ###remove duplicates
        map_list5.append(' '.join([samtools, 'rmdup', rem_dup, os.path.join(bam_dir, '.'.join([sam_file[n],'sort.bam'])), os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.bam']))]))
        ###Add RG
        map_list6.append(' '.join([picard, '='.join(['I', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.bam']))]), '='.join(['O',os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG',typ,'bam']))]), '='.join(['ID',iden[n]]), 'LB=1', 'PL=illumina', 'PU=1', '='.join(['SM',iden[n]]), 'VALIDATION_STRINGENCY=SILENT']))
        ###do index (bai)
        map_list7.append(' '.join([samtools, 'index', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG',typ,'bam']))]))
        ###extract_region_to_realign_with_GATK
        map_list8.append(' '.join([gatk_target, '-T', 'RealignerTargetCreator', '-R', genome, '-I', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG',typ,'bam'])), '-o', 'RTC.intervals']))
        ###realign_with_GATK
        map_list9.append(' '.join([gatk_target, '-T', 'IndelRealigner', '-R', genome, '-I', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG',typ,'bam'])),'-targetIntervals', 'RTC.intervals', '-o', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG.realigned',typ,'bam']))]))
       ###do index with RG (bai)
        map_list10.append(' '.join([samtools, 'index', os.path.join(bam_dir, '.'.join([sam_file[n],'sort.pcr_rem.RG.realigned',typ,'bam']))]))
    ###run & wait until all parallel processes are finished
    parallel(map_list2)
    parallel(map_list3)
    parallel(map_list4)
    parallel(map_list5)
    parallel(map_list6)
    parallel(map_list7)
    parallel(map_list8)
    parallel(map_list9)
    parallel(map_list10)
    ###remove intermediate bam snps
    print '### Removing intermediate files ###'
   # for rem in glob.glob(bam_dir+'/*.uniq_reads.bam')+glob.glob(bam_dir+'/*.sort.bam')+glob.glob(bam_dir+'/*.sort.pcr_rem.bam')+glob.glob(bam_dir+'/*.sort.bam.bai')+glob.glob(bam_dir+'/*.sort.pcr_rem.RG.SNP.bam')+glob.glob(bam_dir+'/*.sort.pcr_rem.RG.SNP.bam.bai'):
    for rem in glob.glob(bam_dir+'/*.uniq_reads.bam')+glob.glob(bam_dir+'/*.sort.bam')+glob.glob(bam_dir+'/*.sort.pcr_rem.bam')+glob.glob(bam_dir+'/*.sort.pcr_rem.RG.SNP.bam')+glob.glob(bam_dir+'/*SNP.bai')+glob.glob(bam_dir+'/*.sort.pcr_rem.RG.SV.bam')+glob.glob(bam_dir+'/*SV.bai')+glob.glob(bam_dir+'/*.sort.pcr_rem.RG.SV.bam.bai'):
        os.remove(rem)
    

def TakeBam(bam_dir,typ):
    ###take all bam files
    print '### Analysing bam files... ###'
    ###Create vcf files directory for SNPs
    spec_dir = os.path.join(bam_dir,typ)
    if not os.path.exists(spec_dir):
        os.makedirs(spec_dir)

    files_bam,no_of_bam=FilesExt(bam_dir,'.bam')
    files_bam = sorted(files_bam)

    ###name without extension
    name = []
    index_list = []
    for each_file in files_bam:
        name.append((os.path.split(each_file)[1]).partition('.bam')[0])
        index_list.append(' '.join([samtools, 'index', each_file]))
    ###run & wait until all parallel processes are finished
    parallel(index_list)
    name = sorted(name)
    
    return spec_dir,files_bam,no_of_bam,name

def CheckEOF(files_bam):
    bad_eof=[]    
    for filename in files_bam:
        eof = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC" + \
              "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
        size = os.path.getsize(filename)
        h = open(filename, "rb") #read only for now
        #Check if it has the EOF already
        h.seek(size - 28)
        data = h.read(28)
        h.close()
        if data != eof:
            bad_eof.append(filename)
#            sys.stderr.write("EOF NOT present in %s\n" % filename)
    if len(bad_eof)>0:
        sys.exit("\nEOF marker absent in %s\nTry to create bam file again" % ','.join(map(str, bad_eof)))
        
def FlagBam(no_of_bam,files_bam):
    flag_list= []
    for l in xrange(no_of_bam):
         flag_list.append(' '.join([samtools, 'flagstat',files_bam[l], '>','.'.join([files_bam[l],'flag'])]))
	    ###run & wait until all parallel processes are finished
    parallel(flag_list)

###Create vcf files
###################    
def VcfForSNPs(param,no_of_bam_snp,files_bam_snp,name_snp,snp_dir):
    ##create vcf files
    print '### Calling SNP/INDEL... ###'
    ###all bam files in uniq vcf
    all_bam=""
    for n in xrange(no_of_bam_snp):
        all_bam = all_bam+' '+files_bam_snp[n]

    ###create a vcf file with annotation...
    if param.ann_file:
        line_ann=open(param.ann_file,'r')
        annotation=line_ann.readline().split('\t')[1]
        line_ann.close()

        print '### Creating vcf files with annotation... ####'
        all_bam_uniq = os.path.join(snp_dir,'All_mpileup.ann.vcf')
        os.system(' '.join([samtools, 'mpileup', '-Q', param.qsam,'-uDf', param.genome, all_bam, '|', bcftools, 'view -Ncvg - >', os.path.join(snp_dir,'All_mpileup.vcf')]))
        PrintLog(' '.join([samtools, 'mpileup', '-Q', param.qsam,'-uDf', param.genome, all_bam, '|', bcftools, 'view -Ncvg - >', os.path.join(snp_dir,'All_mpileup.vcf')]))
        os.system(' '.join([bgzip,'-c', param.ann_file,'>', '.'.join([param.ann_file,'gz'])]))
        PrintLog(' '.join([bgzip,'-c', param.ann_file,'>', '.'.join([param.ann_file,'gz'])]))
        os.system(' '.join([tabix, '-f', '-b', '3', '-e', '4', '.'.join([param.ann_file,'gz'])]))
        PrintLog(' '.join([tabix, '-f', '-b', '3', '-e', '4', '.'.join([param.ann_file,'gz'])]))
    
        if os.path.exists('.'.join([param.ann_file,'gz.tbi'])):
            os.system(" ".join(["cat", os.path.join(snp_dir,"All_mpileup.vcf"), "|", vcf_annotate, "-a", ".".join([param.ann_file,"gz"]), "-c CHROM,INFO/ANN,FROM,TO -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' >", os.path.join(snp_dir,'All_mpileup_ann_dirty.vcf')]))
            PrintLog(" ".join(["cat", os.path.join(snp_dir,"All_mpileup.vcf"), "|", vcf_annotate, "-a", ".".join([param.ann_file,"gz"]), "-c CHROM,INFO/ANN,FROM,TO -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' >", os.path.join(snp_dir,'All_mpileup_ann_dirty.vcf')]))
            sed_ann="'"+"/".join(["s",annotation,"1","g"])+"'"
            os.system(" ".join(["sed", sed_ann, os.path.join(snp_dir,'All_mpileup_ann_dirty.vcf'), ">", all_bam_uniq]))
            PrintLog(" ".join(["sed", sed_ann, os.path.join(snp_dir,'All_mpileup_ann_dirty.vcf'), ">", all_bam_uniq]))
            os.remove(os.path.join(snp_dir,'All_mpileup.vcf'))
            os.remove(os.path.join(snp_dir,'All_mpileup_ann_dirty.vcf'))
        else:
            print '### Error in .tbi file, analysing without annotation... ####'
            all_bam_uniq = os.path.join(snp_dir,'All_mpileup.vcf')

    ### ...or without
    else:
        print '### Creating vcf files... ####'
        all_bam_uniq = os.path.join(snp_dir,'All_mpileup.vcf')
        os.system(' '.join([samtools, 'mpileup', '-Q', param.qsam,'-uDf', param.genome, all_bam, '|', bcftools, 'view -Ncvg - >', all_bam_uniq]))
        PrintLog(' '.join([samtools, 'mpileup', '-Q', param.qsam,'-uDf', param.genome, all_bam, '|', bcftools, 'view -Ncvg - >', all_bam_uniq]))

    ###create single vcf directory
    single_dir = os.path.join(snp_dir,'single')
    if not os.path.exists(single_dir):
        os.makedirs(single_dir)

    ###single vcf for each bam file
    single_bam = []
    single_list= []
    for l in xrange(no_of_bam_snp):
        single_bam.append(os.path.join(single_dir,'.'.join([name_snp[l],'vcf'])))
        single_list.append(' '.join([samtools, 'mpileup', '-Q', '1', '-uDf', param.genome, files_bam_snp[l], '|', bcftools, 'view -Ncvg - >', single_bam[l]]))
    ###run & wait until all parallel processes are finished
    parallel(single_list)

    ###name of vcf uniq file without extension
    uniq_vcf = all_bam_uniq.partition('.vcf')[0]
    print '### vcf files created, calling SNP/INDEL... ###'
    
    return uniq_vcf,single_dir

###Extract AF of single files for SNPs
######################################
def extractAF(files_vcf):
    z=0
    snp_dict=dict()
    lista_di_AF = []
    linea_precedente=None
    linea=None
    for infile_name in sorted(files_vcf):
        ###single files
        vcf_reader = vcf.Reader(open(infile_name, 'r'))
        for record in vcf_reader:
            
            if not linea_precedente and not linea:
                linea_precedente=record
    
            elif linea_precedente and not linea:
                linea=record
                
            if linea_precedente and linea:
                snp_position='_'.join([linea_precedente.CHROM, str(linea_precedente.POS)])
        
                ref_F = float(linea_precedente.INFO['DP4'][0])
                ref_R = float(linea_precedente.INFO['DP4'][1])
                alt_F = float(linea_precedente.INFO['DP4'][2])
                alt_R = float(linea_precedente.INFO['DP4'][3])
            
                AF = (alt_F+alt_R)/(alt_F+alt_R+ref_F+ref_R)
                if not snp_position in snp_dict:
                    snp_dict[snp_position]=list(('0') for _ in range(len(files_vcf)))
    
                snp_dict[snp_position][z] = "%.3f" %AF #record.INFO['DP4']
                        
            linea_precedente=linea
            linea=record
        z+=1
    for snp_position in sorted(snp_dict):
        AF_freq = snp_dict[snp_position]
#        print snp_dict[snp_position]
        lista_di_AF.append(AF_freq)
       
    return snp_dict,lista_di_AF

###Create a genotype.name file if not exist
def GenotypeName(geno_name,result):
    vcf_reader_sample = vcf.Reader(open('.'.join([result,'vcf']), 'r')) ### vcf file
    gn = open(geno_name,'w')
    for record in vcf_reader_sample:
        for call in record.samples:
            gn = open(geno_name,'a')
            gn.write(call.sample)
            gn.write('\n')
        break
    gn.close()

def CheckVcf_snp(ind_geno,allele_freq,ref,alt):
    if ind_geno==ref+ref and float(allele_freq)>=0.251:#invariable alleles,variable frequency
        ind_geno=ind_geno.lower()
    elif ind_geno!=ref+ref and float(allele_freq)<=0.250:#variable alleles, invariable frequency
        ind_geno=ind_geno.lower()
    elif ind_geno==alt+ref and (float(allele_freq)<=0.250 or float(allele_freq)>=0.750):#heterozygous alleles, not heteroz frequency
        ind_geno=ind_geno.lower()
    elif ind_geno!=alt+ref and (float(allele_freq)>0.250 and float(allele_freq)<0.750):#homozygous alleles, heteroz frequency
        ind_geno=ind_geno.lower()
    return ind_geno

def CheckVcf_indel(ind_geno,allele_freq,ref,alt):
    if ind_geno==ref and float(allele_freq)>=0.251:#invariable alleles,variable frequency
        ind_geno=ind_geno.lower()
    elif ind_geno!=ref and float(allele_freq)<=0.250:#variable alleles, invariable frequency
        ind_geno=ind_geno.lower()
    elif ind_geno==alt+'/'+ref and (float(allele_freq)<=0.250 or float(allele_freq)>=0.750):#heterozygous alleles, not heteroz frequency
        ind_geno=ind_geno.lower()
    elif ind_geno!=alt+'/'+ref and (float(allele_freq)>0.250 and float(allele_freq)<0.750):#homozygous alleles, heteroz frequency
        ind_geno=ind_geno.lower()
    return ind_geno
   
###Call SNPs/DIPs
#################
def SNPAnalysis(param,uniq_vcf,name_snp,snp_dir,single_dir):
    print '### Calling SNP/INDEL... ###'
    
    ###create scaffold.list file
    scaffold = os.path.join(param.path,'scaffold.list')
    if not os.path.exists(scaffold):
        sc = open(scaffold,'w')
        sc.close()

    ###create genotype.name file
    geno_name = os.path.join(param.path,'genotype.name')
    if not os.path.exists(geno_name):
        GenotypeName(geno_name,uniq_vcf)

    ###Definition of parameters
    vcf_reader2 = vcf.Reader(open('.'.join([uniq_vcf,'vcf']), 'r')) ### uniq vcf file
    vcf_reader = vcf.Reader(open('.'.join([uniq_vcf,'vcf']), 'r')) ### uniq vcf file
    tab_out = csv.writer(open(os.path.join(snp_dir,'Results_snp_dip.tab'), 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE) ###name of results.tab
    bed_out = csv.writer(open(os.path.join(snp_dir,'Results_snp_dip.bed'), 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE) ###name of results.bed
    genotype_name = open(geno_name, 'rb') ###file with genotype names
    scaffold_list = open(scaffold, 'rb') ###file with scaffold names (optional)

    ###take single vcf files
    filename_vcf='*.vcf'
    files_vcf = glob.glob(os.path.join(single_dir, filename_vcf))
    files_vcf.sort()
    
    ###impute dictionary of function to a new dictionary
    all_vcf,lista_di_AF = extractAF(files_vcf)
    
    ###fill list with all single variations
    lista_coord_single= []
    for key in all_vcf:
        lista_coord_single.append(key)
    lista_coord_single=sorted(lista_coord_single)
    
    titles=['#ID','Chr/scaffold','pos','SNP_quality']
    ###extract names of scaffold and genotype files
    call = []
    scaffold_name = list()
    for nome_genotipo in sorted(genotype_name):
        nome_genotipo = nome_genotipo.rstrip('\n')
        call.append(nome_genotipo)
        titles.append(nome_genotipo)
        titles.extend(('gt_quality','gt_depth','AF'))
    for nome_scaffold in scaffold_list:
        nome_scaffold = nome_scaffold.rstrip('\n')
        scaffold_name.append(nome_scaffold)

    titles.extend(('type','subtype','n_of_invariable_genotypes','ref'))
    if param.ann_file:
        titles.append('annotation')

    tab_out.writerow((titles))

    linea_di_AF = []
    allTheLists = [[] for x in range(len(call))]
    lista_coord_global=[]
    linea_precedente=None
    linea=None
    
    lista_join_coord =[]
    lista_join_final =[]
    ###fill list with all uniq variations
    for record in vcf_reader2:
        coord_global='_'.join([record.CHROM,str(record.POS)])
        lista_coord_global.append(coord_global)
    lista_coord_global=sorted(lista_coord_global)

    ###join variations of both files
    lista_join_coord=set(lista_coord_single) & set(lista_coord_global)

    lista_join_final=lista_join_coord

    ###read all vcf
    ###start analysis in third line: line1=linea_precedente, line2=linea, line3=record
    for record in vcf_reader:
        ###impute reading lines
        ###this condition is only for the first line of the file
        if not linea_precedente and not linea:
            linea_precedente=record
        ###this condition is only for the second line of the file
        elif linea_precedente and not linea:
            linea=record

        if linea:
            snp_pos='_'.join([linea.CHROM,str(linea.POS)])
            ###read only coincident lines of vcf 
            if snp_pos in lista_join_final:
                linea_di_AF=all_vcf[snp_pos]

        ###this condition is for third line and so on
        if linea_precedente and linea:
            tmp2=list()
            geno_qual=list()
            geno_deep=list()
            all_GT=list()
            ###impute values of each genotype
            for sample in call:
                genotype = linea.genotype(sample)
#                print genotype['GT']
                tmp2.append(genotype['GT'])###genotype
                geno_qual.append(genotype['GQ'])###quality
                geno_deep.append(genotype['DP'])###deep
    
            ###checking snps in flanking region
            if (
                snp_pos in lista_join_final
                and linea.POS < (record.POS-param.flank)
                and linea.POS > (linea_precedente.POS+param.flank)
                and (linea.CHROM in scaffold_name or  not scaffold_name)
                ):
                
                final_line=[snp_pos,linea.CHROM,linea.POS,linea.QUAL]     
                final_bed_line=[linea.CHROM,linea.POS-param.flank,linea.POS+param.flank,linea.POS,linea.QUAL, '+']
                ###reference or alternative bases
                for GT in tmp2:
                    alternative=str(linea.ALT[0])
                    if linea.var_type=='snp':
                        if GT == '0/0':
                            GT=linea.REF+linea.REF
                        elif GT == '1/1':
                            GT=alternative+alternative
                        elif GT == '0/1':
                            GT=alternative+linea.REF
                    else:
                        if GT == '0/0':
                            GT=linea.REF
                        elif GT == '1/1':
                            GT=alternative
                        elif GT == '0/1':
                            GT=alternative+'/'+linea.REF
                        
                    all_GT.append(GT)

                ###append values of each genotype
                gt_counter=0
                for i in range(len(geno_qual)):
                    allTheLists[i].append((linea_di_AF[i]))

                    if linea.var_type=='snp': 
                        ind_geno=CheckVcf_snp(all_GT[i],linea_di_AF[i],str(linea.REF),str(linea.ALT[0]))
                        final_line.append(ind_geno)
                    else:
                        ind_geno=CheckVcf_indel(all_GT[i],linea_di_AF[i],str(linea.REF),str(linea.ALT[0]))
                        final_line.append(ind_geno)

                    final_line.append(geno_qual[i])
                    final_line.append(geno_deep[i])
                    final_line.append((linea_di_AF[i]))
                    if linea.var_type == 'snp' and ind_geno == linea.REF+linea.REF:
                        gt_counter+=1#count number of invariable genotypes
                    elif linea.var_type == 'indel' and ind_geno == linea.REF:
                        gt_counter+=1#count number of invariable genotypes

                final_line.append(linea.var_type)
                final_line.append(linea.var_subtype)
                final_line.append(gt_counter)
                final_line.append(linea.REF)
                for key in linea.INFO.iteritems() :
                    if 'ANN' in key:
                        final_line.append(linea.INFO['ANN']) 
                tab_out.writerow(final_line)
                bed_out.writerow(final_bed_line)

            linea_precedente=linea
            linea=record
                
    #Statistics
    invariable = []
    homozigous = []
    ethero= []
    for h in range(len(call)):
        invariabili=0
        omozigoti=0
        eterozigoti=0
        
        for j in range(len(allTheLists[0])):
            if float(allTheLists[h][j]) <= 0.250:
                invariabili+=1
            elif float(allTheLists[h][j]) >= 0.750:
                omozigoti+=1
            else:
                eterozigoti+=1    
                
        invariable.append(invariabili)
        homozigous.append(omozigoti)
        ethero.append(eterozigoti)

    ###create statistics file
    p_invariable=str(invariable)        
    p_homozigous=str(homozigous)
    p_ethero=str(ethero)
    stats = os.path.join(snp_dir,'SNPs.stats')
    st = open(stats,'w')
    st.write(' '.join(['No of invariable site', p_invariable, '\n']))
    st.write(' '.join(['No of homozygous site', p_homozigous, '\n']))
    st.write(' '.join(['No of heterozygous site', p_ethero, '\n']))
    st.close()


###Call SVs
###########
###Only for Paired-Ends
def VcfForSVs(param,no_of_bam_sv,files_bam_sv,sv_dir):
    ###create pindel configuration file with bam files
    conf_pindel = os.path.join(param.path,'configuration.pindel')
    cp = open(conf_pindel,'w')
    for c in xrange(no_of_bam_sv):
        cp = open(conf_pindel,'a')
        cp.write('\t'.join([files_bam_sv[c], param.pe_distance, os.path.split(files_bam_sv[c])[1].partition('.')[0]]))
        cp.write('\n')
    cp.close()

    ###extract genome name
    genome_name=(os.path.split(param.genome)[1]).partition('.f')[0]
    ###name of output file of results
    result_pindel=os.path.join(sv_dir,'res_pindel')

    ###run pindel 024t version for genomes with scaffolds
    if param.type_geno=='s':
        ###run pindel (old version)
        os.system(' '.join([pindel_t,'-f',param.genome,'-i', conf_pindel, '-o', result_pindel, param.pindel_config]))
        PrintLog(' '.join([pindel_t,'-f',param.genome,'-i', conf_pindel, '-o', result_pindel, param.pindel_config]))
        ###convert results in vcf SV 
        os.system(' '.join([pindel2vcf_t, '-G','-P', result_pindel, '-r ', param.genome, '-R', genome_name, '-d',str(datetime.date.today()),param.pindel2vcf_config]))
        PrintLog(' '.join([pindel2vcf_t, '-G','-P', result_pindel, '-r ', param.genome, '-R', genome_name, '-d',str(datetime.date.today()),param.pindel2vcf_config]))

    ###run pindel 024v version for genomes with chromosomes
    elif param.type_geno=='c':
        ###take all chromosomes of the genome in a list
        os.system(" ".join(["grep '>'", param.genome, "| sed 's/>//g' >", os.path.join(sv_dir, 'genome_keys.txt')]))
        lines_geno=[]
        g_keys = open(os.path.join(sv_dir, 'genome_keys.txt'),'r')
        lines_geno=g_keys.readlines()
        g_keys.close()
        result_pindel_v=os.path.split(result_pindel)[1]

	    ###variables for join all results in one pindel file
        BP=''
        CEM=''
        D=''
        INV=''
        LI=''
        SI=''
        TD=''

	    ###create a directory for long insertions
        li_dir=os.path.join(sv_dir,'LI')
        if not os.path.exists(li_dir):
            os.makedirs(li_dir)
	    
	    ###analyse each chromosome
        for chrs in lines_geno:
            chrs=chrs.partition(' ')[0].rstrip('\n')
            chr_dir=os.path.join(sv_dir,chrs)                
            if not os.path.exists(chr_dir):
                os.makedirs(chr_dir)
		###without parallelization
            os.system(' '.join([pindel_v, '-f', param.genome, '-i', conf_pindel, '-o', os.path.join(chr_dir,'.'.join([result_pindel_v,chrs])), '-T', param.threads, '-c', chrs,param.pindel_config]))
            PrintLog(' '.join([pindel_v, '-f', param.genome, '-i', conf_pindel, '-o', os.path.join(chr_dir,'.'.join([result_pindel_v,chrs])), '-T', param.threads, '-c', chrs,param.pindel_config]))
            ###join all results of all chromosomes
            BP=BP+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'BP'])]))
            CEM=CEM+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'CloseEndMapped'])]))
            D=D+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'D'])]))
            INV=INV+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'INV'])]))
            LI=LI+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'LI'])]))
            SI=SI+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'SI'])]))
            TD=TD+' '+os.path.join(chr_dir,'.'.join([result_pindel_v,'_'.join([chrs,'TD'])]))


	    ###create an unique file for each result
        os.system(' '.join(['cat', BP, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'BP']))]))
        os.system(' '.join(['cat', CEM, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'CloseEndMapped']))]))
        os.system(' '.join(['cat', D, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'D']))]))
        os.system(' '.join(['cat', INV, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'INV']))]))
        os.system(' '.join(['cat', LI, '>', os.path.join(li_dir,'_'.join([result_pindel_v,'LI']))]))###LI directory
        os.system(' '.join(['cat', SI, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'SI']))]))
        os.system(' '.join(['cat', TD, '>', os.path.join(sv_dir,'_'.join([result_pindel_v,'TD']))]))

        ###convert results in vcf SV 
    os.system(' '.join([pindel2vcf_v, '-G','-P', result_pindel, '-r ', param.genome, '-R', genome_name, '-d',str(datetime.date.today()),param.pindel2vcf_config]))
    PrintLog(' '.join([pindel2vcf_v, '-G','-P', result_pindel, '-r ', param.genome, '-R', genome_name, '-d',str(datetime.date.today()),param.pindel2vcf_config]))
    
    return result_pindel

def SVAnalysis(param,sv_dir,result_pindel):
    print '### Calling SVs... ###'
    ###Analyse pindel results
    vcf_reader_p = vcf.Reader(open('.'.join([result_pindel,'vcf']), 'r')) ### vcf file
    tab_r_out = csv.writer(open(os.path.join(sv_dir,'Results_svs.tab'), 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE) ###name of results.tab

    ###create genotype.name file
    geno_name = os.path.join(param.path,'genotype.name')
    if not os.path.exists(geno_name):
        GenotypeName(geno_name,result_pindel)

    genotype_name = open(geno_name, 'rb') ###file with genotype names

    titles=['#ID','Chr/scaffold','from','to','diff_length']
    ###extract names of genotype file
    for nome_genotipo in sorted(genotype_name):
        nome_genotipo = nome_genotipo.rstrip('\n')
        titles.append(nome_genotipo)
        titles.extend(('ref_depth','allele_depth'))

    titles.extend(('reference','alternative','length_of_base_pairs','seq_of_base_pairs_identical','n_of_bases_inserted','type'))

    tab_r_out.writerow((titles))

    for record in vcf_reader_p:
        sample_GT=[]    
        sample_RD=[]    
        sample_AD=[]
        sv_pos='_'.join([record.CHROM,str(record.POS)])

        for genos in record.samples:
            sample_GT.append(genos['GT'])
            if param.type_geno=='c' :
                sample_RD.append(int(genos['RD']))
            else:
                sample_RD.append('-')
            sample_AD.append(genos['AD'])
        if int(record.INFO['HOMLEN'])==0:
            hom_seq='-'
        else:
            hom_seq=record.INFO['HOMSEQ']
        if 'NTLEN' in record.INFO:
            nt_len=record.INFO['NTLEN'][0]
        else:
            nt_len='-'

        final_sv=[sv_pos, record.CHROM, record.POS, record.INFO['END'], record.INFO['SVLEN']]

        for n in range(len(sample_GT)):
            if sample_GT[n]:
                final_sv.append(sample_GT[n])
            else:
                final_sv.append('.')
            final_sv.append(sample_RD[n])
            final_sv.append(sample_AD[n])
        
        final_sv.append(record.REF)
        final_sv.append(record.ALT[0])
        final_sv.append(record.INFO['HOMLEN'])
        final_sv.append(hom_seq)
        final_sv.append(nt_len)
        final_sv.append(record.INFO['SVTYPE'])

        tab_r_out.writerow(final_sv)

###Search reads with unique alignments
def FilterUniqSam(sam_file):
    bad=set()
    with open(sam_file,'rb') as file_new:
        for line_new in file_new:
            fields=line_new.rstrip('\n').split('\t')
            if fields[0].find('@')==0:
                continue
            if int(fields[1])>=256 or int(fields[4])==0:#searching for double alignments and bad mapq
                bad.add(fields[0])
    
    file_out=open('.'.join([sam_file.partition('.sam')[0],'uniq_reads','sam']),'wb')
    with open(sam_file,'rb') as file_new:
        for line_new in file_new:
            fields=line_new.rstrip('\n').split('\t')
            if fields[0] in bad or (not fields[0].find('@')==0 and bool(int(fields[1]) & 0x2)==False) :#searching for properly paired
                continue

            file_out.write('\t'.join(fields))
            file_out.write('\n')


def main():
    #take arguments
    args=InputPar()

    print_args=vars(args)#take args in a dict and print in a log
    for k,v in print_args.iteritems():
        strategy_out='='.join([str(k),str(v)])
        PrintLog(strategy_out)
        
    #checking input args
    fastq,args.rem_dup=CheckInput(args)#extract fastq extension

    ####Process with fastq files
    if glob.glob(args.path+'/*'+fastq):
        ###take all fastq or fastq.gz files
        files,no_of_files=FilesExt(args.path,fastq)

        ###take the same name as we have created (without extension)
        name_ok=[]
        for infile_name in files:
            name_ok.append((os.path.split(infile_name)[1]).partition(fastq)[0])
        
        name_trim = sorted(name_ok)
    
        print '### Alignment in course... ###'
        CreateIndex(args.genome)
            
        ###Create bam files directory
        bam_dir = os.path.join(args.path,args.bwa_dir)
        if not os.path.exists(bam_dir):
            os.makedirs(bam_dir)
    
        ###Perform the BWA alghoritm to map pair-end reads on the reference genome
        ##########################################################################
        sam_file,iden=BWAMap(args,name_trim,bam_dir,no_of_files,fastq)
    
        ######convert sam to bam and purify for SNP calling
        if int(args.analysis) == 1 or int(args.analysis) == 3:
            print '### Creating bam files for SNPs... ###'
            ###Create vcf files directory for SNPs
            snp_dir = os.path.join(bam_dir,"SNP")
            if not os.path.exists(snp_dir):
                os.makedirs(snp_dir)

            if int(args.analysis) == 1:#only snps. create dirty bam
                map_list2=[]
                for n in xrange(len(sam_file)):
                    map_list2.append(' '.join([samtools ,'view', '-Sb', os.path.join(bam_dir,'.'.join([sam_file[n],'sam'])), '>', os.path.join(bam_dir, '.'.join([sam_file[n],'bam']))]))
                ###run & wait until all parallel processes are finished
                parallel(map_list2)

            ###select uniq reads
            sam_file_snp=[]
            full_sam_file=[]
            for each_sam in sam_file:
                full_sam_file.append(os.path.join(bam_dir,'.'.join([each_sam,'sam'])))
                sam_file_snp.append('.'.join([each_sam,'uniq_reads']))
            #parallelize processes
            pool = multiprocessing.Pool()
            pool.map(FilterUniqSam, full_sam_file)

            typ='SNP'
            SamToBam(bam_dir,sam_file_snp,args.rem_dup,typ,iden,args.genome)

        ######convert sam to bam and purify for extract SVs (without filtering uniq reads) ONLY Paired-Ends
        if (int(args.analysis) == 2 or int(args.analysis) == 3) and args.ends == 'p':
            print '### Creating bam files for SVs... ###'
            ###Create vcf files directory for SVs
            sv_dir = os.path.join(bam_dir,"SV")
            if not os.path.exists(sv_dir):
                os.makedirs(sv_dir)
            typ='SV'
            SamToBam(bam_dir,sam_file,args.rem_dup,typ,iden,args.genome)
        
        for rem in glob.glob(bam_dir+'/*.sai')+glob.glob(bam_dir+'/*.sam'):
            os.remove(rem)

        ###take all bam files
        if int(args.analysis) == 1 or int(args.analysis) == 3:
            files_bam_snp,no_of_bam_snp=FilesExt(bam_dir,'.SNP.bam')
            files_bam_snp = sorted(files_bam_snp)
            ###name without extension
            name_snp = []
            for each_file in files_bam_snp:
                name_snp.append((os.path.split(each_file)[1]).partition('.bam')[0])
            name_snp = sorted(name_snp)
    
        if (int(args.analysis) == 2 or int(args.analysis) == 3) :
            files_bam_sv,no_of_bam_sv=FilesExt(bam_dir,'.SV.bam')
            files_bam_sv = sorted(files_bam_sv)
            ###name without extension
            name_sv = []
            for each_file in files_bam_sv:
                name_sv.append((os.path.split(each_file)[1]).partition('.bam')[0])
            name_sv = sorted(name_sv)
    
        print '### Alignment finished, analysing variation... ###'
        
    ###Check extension of input files
    elif glob.glob(args.path+'/*.bam'):
        bam_dir = args.path
        if int(args.analysis) == 1 or int(args.analysis) == 3:
            typ='SNP'
            snp_dir,files_bam_snp,no_of_bam_snp,name_snp=TakeBam(bam_dir,typ)#take bam files form path

        if (int(args.analysis) == 2 or int(args.analysis) == 3) :
            typ='SV'
            sv_dir,files_bam_sv,no_of_bam_sv,name_sv=TakeBam(bam_dir,typ)#take bam files form path

    ###Check differences between bam files
    if not glob.glob(args.path+'/*.vcf'):
        if int(args.analysis) == 1 or int(args.analysis) == 3:
            CheckEOF(files_bam_snp)
            FlagBam(no_of_bam_snp,files_bam_snp)
        if (int(args.analysis) == 2 or int(args.analysis) == 3) :
            CheckEOF(files_bam_sv)
            FlagBam(no_of_bam_sv,files_bam_sv)

    if int(args.analysis) == 1 or int(args.analysis) == 3:
        ###Check extension of input files
        if glob.glob(args.path+'/*.vcf') and glob.glob(os.path.join(args.path,'single')+'/*.vcf'):
            snp_dir=args.path
            name_snp = []
            single_dir = os.path.join(snp_dir,'single')
            ###check existence and size of single vcf files
            for v_size in glob.glob(single_dir+'/*.vcf'):
                sizes = os.path.getsize(v_size)
                if sizes == 0:
                    sys.exit(' '.join(['!!!Empty',v_size,'!!!','\nExiting process']))
                name_snp.append((os.path.split(v_size)[1]).partition('.vcf')[0])
            name_snp = sorted(name_snp)
                
            ###name of vcf uniq file without extension
            filename_vcf='*.vcf'
            uniq_vcf = glob.glob(os.path.join(snp_dir, filename_vcf))[0].partition('.vcf')[0]
        else:
            uniq_vcf,single_dir=VcfForSNPs(args,no_of_bam_snp,files_bam_snp,name_snp,snp_dir)
    
        SNPAnalysis(args,uniq_vcf,name_snp,snp_dir,single_dir)


    if (int(args.analysis) == 2 or int(args.analysis) == 3) and args.ends == 'p':
        ###taking vcf files from path
        if glob.glob(args.path+'/*.vcf'):
            sv_dir=args.path
            name_sv = []
            ###name of vcf uniq file without extension
            filename_vcf='*.vcf'
            result_pindel = glob.glob(os.path.join(sv_dir, filename_vcf))[0].partition('.vcf')[0]

        else:
            result_pindel=VcfForSVs(args,no_of_bam_sv,files_bam_sv,sv_dir)

        SVAnalysis(args,sv_dir,result_pindel)

    print '### Process finished ###'


if __name__ == "__main__":
    main()

