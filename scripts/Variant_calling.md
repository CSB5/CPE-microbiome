# Variant Calling
## Get Alignment
```bash
NB_MISMATCH=3 #2,3 or 4

# Using samtools mismatch definition
# bwa mem -t 12 $REF $FASTQ1 $FASTQ2 | awk '"'{if($6 !~ /S/){print $0}}'"' | samtools view  -Sub -F 4 -e '[NM]<=${NB_MISMATCH}' - | samtools sort - -o $OUTPUTPREFIX.bam && samtools index $OUTPUTPREFIX.bam

# Using custom perl script to remove reads with n mismatches based on CIGAR
bwa mem -t 12 $REF $FASTQ1 $FASTQ2 | awk '"'{if($6 !~ /S/){print $0}}'"' | perl utilities/filtermm.pl ${NB_MISMATCH} | samtools view  -Sub -F 4 - | samtools sort - -o $OUTPUTPREFIX.bam && samtools index $OUTPUTPREFIX.bam
```

## SNV calling
```bash
# Input - Metagenome reads
lofreq viterbi -f $REF $OUTPUTPREFIX.bam | samtools sort - | lofreq indelqual --dindel -f $REF -o $INDEL_OUT - && lofreq call --call-indels -f $REF -o $VCF_OUT --verbose $INDEL_OUT

# Input - Assembled contigs
nucmer --maxmatch --prefix=$OUTPUTPREFIX $REF $ASSEMBLY
show-snps -Tlr $OUTPUTPREFIX.delta | cut --complement -f7,8 > $OUTPUTPREFIX.snps
```

## Translate coordinates using flo + liftover or crossmap
```bash
# Generate chain file
rake -f $FLO_CODE_DIR/Rakefile; cp ./run/liftover.chn $CHAIN_OUTPUT
liftOver -bedPlus=3 temp.bed $CHAIN_OUTPUT output.bed unlifted.bed
#CrossMap.py bed $CHAIN $OUTPUTPREFIX.bam ${OUTPUTPREFIX}_indel.bam.vcf
```

## Filter Ambiguous Region with Rescue
```bash
AMBREGION='../data/ec_amb_rescued.bed'# Regions where kp reads mapped to ec reference with high coverage. Variants with low frequency in the alignment are rescued.
#AMBREGION='../data/kp_amb_rescued.bed'
bedtools intersect -v -a $VCF_IN -b $AMBREGION -wa > $VCF_OUT"
```
