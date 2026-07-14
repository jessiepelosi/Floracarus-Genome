# Assembly and Annotation of the <i>Floracarus perrepae</i> genome

## Contaminant Filtering 

To identify and remove reads from contaminating sources, raw reads were mapped to the human genome (GRCh38.p14) and <i>Lygodium microphyllum</i> reference genome and plastome with minimap2 v 2.24-r1122.  

```
minimap2 -ax map-hifi GCF_000001405.40_GRCh38.p14_genomic.fna \
m84082_260417_030759_s2.hifi_reads.16_UDI_14_F02_F--16_UDI_14_F02_R.hifi_reads.fastq.gz | samtools fastq -n -f 4 - > 9534_no_human.fastq.gz

minimap2 -ax map-hifi Lygodium_microphyllum_v1.1.9.fa 9534_no_human.fastq.gz > 9534_no_human_lygodium.fastq.gz
```

## k-mer Profile 

Unmapped reads were used for all downstream assembly steps. From these reads, we generated a k-mer profile with KMC v3.2.4 (Kokot et al. 2017) with k=21, and estimated genome size with GenomeScope2.0 (Ranallo-Benavidez et al. 2020). 

```
$kmc -k21 9534_no_human_plastome_lygo.fastq flpe_21kmers kmc_tmp
```

## Mitochondrial Genome Assembly 

The mitochondrial genome was assembled using MitoHiFi v3.2.1 with the <i>Nothopoda</i> sp. reference mitogenome (NCBI Accession OQ934080) as the reference. We verified the assembled sequenced belong to Floracarus perrepae using BLASTn of the COX 1 sequence against the nr database on NCBI (e-value = 9e-179, 99.7-100% identity match to <i>F. perrepae</i>). 

```
mitohifi='/groups/kdlugosch/jessiepelosi/local_prgms/mitohifi.sif' 
singularity exec $mitohifi mitohifi.py -r 9534_no_human_lygodium.fastq.gz -f OQ934080.1.fasta -g OQ934080.1.gb -t 16 -a animal
```

## Metagenome Assembly and Characterization

Given that we had a pooled sample, we used myloasm v0.5.1 (Shaw et al. 2026), a metagenome assembler, to generate a preliminary assembly of the filtered reads. The resulting assembly graph was inspected with Bandage v0.8.1 (Wick et al. 2015) and contigs were taxonomically categorized with Blobtools2 v4.5.3 (Challis et al. 2020) after identifying BLAST hits to the non-redundant nucleotide database (Camacho et al. 2009) and mapping reads to the assembly with minimap2 v2.24-r1122 (Li 2018). 

```
myloasm 9534_no_human_lygodium.fastq.gz -t 64 --hifi -o FLPE_myloasm
minimap2 -ax map-hifi assembly_primary.fa 9534_no_human_lygodium.fastq.gz | samtools view -S -b | samtools sort -o aligned_sorted.bam
blastn -db $blastdb -query assembly_primary.fa -outfmt "6 qseqid staxids bitscore std" -evalue 1e-25 -max_hsps 1 -max_target_seqs 10 -num_threads 128 -out blast_hits_may3.out
```

## Hifiasm Assembly 

Given the high coverage of our target species in the filtered data, we generated five random subsets with seqkit v2.8.2 (Shen et al. 2016) each totalling 5% of the data (~120x coverage). Hifiasm v0.25.0-r726 Cheng et al. 2021) was run on each subset with the following parameter modifications: k-mer size of 21 (-k 21), aggressive purging of alternate haplotypes (-l3), targeted homozygous genome size of 40 Mb (--hg-size 40m), and specified the number of haplotypes (--n-hap 492). 

```
seqkit sample 
hifiasm -o 9534_no_human_plastome_lygo_5Perc_rep1 -t 96 $reads --hg-size 40m --n-hap 492 -k 21 -l3
```

## Assembly Merging 

We used mummer v4.0.1 (Marçais et al. 2018) and quickmerge vXX (Chakraborty et al. 2016) to merge the best hifiasm assembly (based on BUSCO score) and the 21.3Mb component from myloasm and retained all contigs that were not merged from both assemblies. 


```
nucmer -l 100 myloasm_component2.fa 9534_5Perc_rep3.bp.p_ctg.fa
delta-filter -r -q -l 10000 out.delta > out.rq.delta
quickmerge -d out.rq.delta -q 9534_5Perc_rep1.bp.p_ctg.fa -r myloasm_component2.fa -hco 5.0 -c 1.5 -l 1000 -ml 5000 -p rep3_myloasm

grep ">" myloasm_component2.fa | seqkit grep -f merged_rep3_myloasm.fasta > missing.fa
cat missing.fa merged_rep3_myloasm.fasta > merged_rep3_myloasm+missing.fa
```

## Competitive Mapping

```
minimap2 -ax map-hifi $assembly $reads --secondary=no | samtools view -S -b | samtools sort -o aligned_sorted.bam
samtools index aligned_sorted.bam
samtools idxstats aligned_sorted.bam

samtools idxstats aligned_sorted.bam | grep "^Ref2_" | awk '{print $1"\t0\t"$2}' > ref2_regions.bed
samtools view -b -L ref2_regions.bed aligned_sorted.bam > component2_reads.bam

samtools idxstats aligned_sorted.bam | grep "^Ref1_" | awk '{print $1"\t0\t"$2}' > ref1_regions.bed
samtools view -b -L ref1_regions.bed aligned_sorted.bam > component1-3-4_reads.bam

```

## Annotation 

Custom repeat libraries were generated for the <i>Floracarus perrepae</i> and bycatch genomes with RepeatModeler2 with LTRStruct enabled. 
```
asm="Floracarus_perrepae_v1.fa"
dfam='/groups/kdlugosch/jessiepelosi/local_prgms/dfam-tetools-latest.sif'

singularity exec $dfam BuildDatabase -name Floracarus $asm
singularity exec $dfam RepeatModeler -database Floracarus -threads 96 -LTRStruct
singularity exec $dfam RepeatMasker -pa 96 -norna -lib Floracarus-families.fa -no_is -gff -a -xsmall $asm
```

The resulting soft-masked genomes were used as input for the BRAKER4 pipeline. We used the predicted proteomes of 18 mite genomes and the Arthopoda ODB12 database as external evidence. 
```
snakemake --executor slurm --default-resources slurm_partition=standard mem_mb=120000 \
slurm_account=kdlugosch --cores 48 --jobs 48 --use-singularity \
--singularity-args "-B /xdisk -B /groups" --singularity-prefix .singularity_cache \
--snakefile /groups/kdlugosch/jessiepelosi/local_prgms/BRAKER4/Snakefile \
--latency-wait 120 --restart-times 3
```
