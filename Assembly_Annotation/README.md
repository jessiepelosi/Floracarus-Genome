# Assembly and Annotation of the <i>Floracarus perrepae</i> genome

## Contaminant Filtering 

To identify and remove reads from contaminating sources, raw reads were mapped to the human genome (GRCh38.p14) and <i>Lygodium microphyllum</i> reference genome and plastome with minimap2 v 2.24-r1122.  

```
minimap2 -ax map-hifi GCF_000001405.40_GRCh38.p14_genomic.fna \
m84082_260417_030759_s2.hifi_reads.16_UDI_14_F02_F--16_UDI_14_F02_R.hifi_reads.fastq.gz | samtools fastq -n -f 4 - > 9534_no_human.fastq.gz

minimap2 -ax map-hifi Lygodium_microphyllum_v1.1.9.fa 9534_no_human.fastq.gz > 9534_no_human_lygodium.fastq.gz
```


## Mitochondrial Genome Assembly 

The mitochondrial genome was assembled using MitoHiFi v3.2.1 with the <i>Nothopoda</i> sp. reference mitogenome (NCBI Accession OQ934080) as the reference. We verified the assembled sequenced belong to Floracarus perrepae using BLASTn of the COX 1 sequence against the nr database on NCBI (e-value = 9e-179, 99.7-100% identity match to <i>F. perrepae</i>). 

```
mitohifi='/groups/kdlugosch/jessiepelosi/local_prgms/mitohifi.sif' 
singularity exec $mitohifi mitohifi.py -r 9534_no_human_lygodium.fastq.gz -f OQ934080.1.fasta -g OQ934080.1.gb -t 16 -a animal
```

