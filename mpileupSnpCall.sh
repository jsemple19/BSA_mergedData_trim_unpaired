#! /bin/bash
#last modified 20160818

GENOME_VER="PRJNA13758.WS250"
BASEDIR=/users/blehner/jsemple/seqResults/mergedData
GENOMEDIR=/users/blehner/jsemple/seqResults/GenomeBuilds/$GENOME_VER

# to create env variable for samtools:
export SAMTOOLS_HOME="/software/bl/el6.3/samtools-1.3.1"

#mkdir -p $BASEDIR/bamFiles/$GENOME_VER
#mkdir -p $BASEDIR/vcfFiles/$GENOME_VER


genomefile=`ls $GENOMEDIR/*.fa`

#index genome
#$SAMTOOLS_HOME/samtools faidx $genomefile

# create arrays of file names
samfiles=(`ls $BASEDIR/samFiles/$GENOME_VER/*.sam`)

# ordinal number of samFile input from command line arg
i=$1
let i=i-1
r1=`echo ${samfiles[$i]%_*_aln_pe.sam} | cut -f9 -d'/'`
outbam=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".bam"`

##convert samfiles to bamfiles (slow):
#$SAMTOOLS_HOME/samtools view -b -S -o $outbam ${samfiles[$i]}

##sort the bam file:
outsort=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".sorted.bam"`	
#$SAMTOOLS_HOME/samtools sort $outbam -T $outsort -o $outsort 

##count variants at SNP sites from .bed file
#bedfile=`ls $GENOMEDIR/AnnSNPs*.pos`
bedfile=`ls $GENOMEDIR/SNVs_N2-CB4856*.pos`

#names of sorted files now end in .bam
#insorted=`echo $outsort".bam"`

#create variant count file
vcffile=`echo $BASEDIR/vcfFiles/$GENOME_VER/$r1"_raw.vcf"`

export BCFTOOLS_HOME="/software/bl/el6.3/bcftools-1.3.1"

#run samtools mpileup:qq and bcftools
#$SAMTOOLS_HOME/samtools mpileup -f $genomefile -l $bedfile -Q 20 -uBg $outsort | $BCFTOOLS_HOME/bcftools view -> $vcffile

#$SAMTOOLS_HOME/samtools mpileup -f $genomefile -uBg $outsort | $BCFTOOLS_HOME/bcftools call -r $bedfile -O u -m | $BCFTOOLS_HOME/bcftools view -> $vcffile

#$SAMTOOLS_HOME/samtools mpileup -f $genomefile -l $bedfile -Q 20 -uBg $outsort | $BCFTOOLS_HOME/bcftools view | $BCFTOOLS_HOME/bcftools norm -m-both -f $genomefile -> $vcffile

mpfile=`echo $BASEDIR/vcfFiles/$GENOME_VER/$r1"_mp.bcf"`
#$SAMTOOLS_HOME/samtools mpileup -f $genomefile -Q 20 -l $bedfile -BIg $outsort -o $mpfile 

# ver6
$BCFTOOLS_HOME/bcftools index $mpfile
$BCFTOOLS_HOME/bcftools call -R $bedfile -m -P0 -Ov $mpfile -o $vcffile

#ver7 
#$BCFTOOLS_HOME/bcftools call -R $bedfile -c -Ov $mpfile -o $vcffile

# ver8
#$BCFTOOLS_HOME/bcftools call -R $bedfile -c -P0 -Ov $mpfile -o $vcffile

