
# Work flow and scripts for *Typha laxmannii* genome assembly and annotation, and population genomics analysis in Tianshan mountains region

  Overview: we *de novo* assembled the genome of *Typha laxmannii* using a combination of Oxford Nanopore long reads, DNBSEQ-T7 short reads, and Hi-C data. We then annotated the genome for genes and repetitive elements. Finally, we performed population genomics analyses using resequencing data from multiple individuals collected across the Tianshan mountains region.


## 1. Genome Assembly
  The DNA sequence data are available in the National Genomics Data Center, BioProject no. PRJCA045581, reference genome is available under GWH No. GWHGPYE00000000.1. https://ngdc.cncb.ac.cn/gsub/submit/bioproject/list

### 1.1 Sequencing Data 
  - Oxford Nanopore long reads: ~126x coverage; ONT_pass.fq.gz
  - DNBSEQ-T7 short reads: ~140x coverage; NGS_1.fq.gz, NGS_2.fq.gz
  - Hi-C data: ~116x coverage; HiC_1.fq.gz, HiC_2.fq.gz
  - RNA-seq data:
    ROOT: root1_1.fq.gz, root1_2.fq.gz, root2_1.fq.gz, root2_2.fq.gz
    LEAF: leaf1_1.fq.gz, leaf1_2.fq.gz, leaf2_1.fq.gz, leaf2_2.fq.gz
    STEM: stem1_1.fq.gz, stem1_2.fq.gz, stem2_1.fq.gz, stem2_2.fq.gz
    Shoot: shoot1_1.fq.gz, shoot1_2.fq.gz, shoot2_1.fq.gz, shoot2_2.fq.gz

### 1.2 Quality Control of the raw data
    NanoFilt: https://github.com/wdecoster/nanofilt
    NanoFilt -q 10 -l 500 -t 40 -o ONT_pass_filter.fq.gz ONT_pass.fq.gz

    Fastp: https://github.com/OpenGene/fastp
    fastp -i NGS_1.fq.gz -I NGS_2.fq.gz -o NGS_1_filter.fq.gz -O NGS_2_filter.fq.gz
    fastp -i HiC_1.fq.gz -I HiC_2.fq.gz -o HiC_1_filter.fq.gz -O HiC_2_filter.fq.gz
    fastp -i root1_1.fq.gz -I root1_2.fq.gz -o root1_1_filter.fq.gz -O root1_2_filter.fq.gz
    fastp -i root2_1.fq.gz -I root2_2.fq.gz -o root2_1_filter.fq.gz -O root2_2_filter.fq.gz
    fastp -i leaf1_1.fq.gz -I leaf1_2.fq.gz -o leaf1_1_filter.fq.gz -O leaf1_2_filter.fq.gz
    fastp -i leaf2_1.fq.gz -I leaf2_2.fq.gz -o leaf2_1_filter.fq.gz -O leaf2_2_filter.fq.gz
    fastp -i stem1_1.fq.gz -I stem1_2.fq.gz -o stem1_1_filter.fq.gz -O stem1_2_filter.fq.gz
    fastp -i stem2_1.fq.gz -I stem2_2.fq.gz -o stem2_1_filter.fq.gz -O stem2_2_filter.fq.gz
    fastp -i shoot1_1.fq.gz -I shoot1_2.fq.gz -o shoot1_1_filter.fq.gz -O shoot1_2_filter.fq.gz
    fastp -i shoot2_1.fq.gz -I shoot2_2.fq.gz -o shoot2_1_filter.fq.gz -O shoot2_2_filter.fq.gz

### 1.3 Genome Size Estimation
    GenomeScope: http://qb.cshl.edu/genomescope/
    jellyfish count -m 19 -s 1000000000 -t 40 -C NGS_1_filter.fq.gz NGS_2_filter.fq.gz -o typha19.jf
    jellyfish histo -t 40 typha19.jf > typha19.histo
    Rscript genomescope.R typha19.histo 19 150 typha_genomescope

### 1.4 Genome Assembly
    NextDenovo: https://github.com/Nextomics/NextDenovo
    nextDenovo Nextdenovo.cfg

### 1.5 Genome Polishing
#### Polish the assembly both by NGS and nanopore reads twice.
    NextPolish: https://github.com/Nextomics/NextPolish
    nextPolish NextPolish.cfg

#### estimate the assembly quality with BUSCO 
    BUSCO: https://busco.ezlab.org/
    busco -i typha_polish.fasta -l embryophyta_odb10 -o busco_polish -m genome -c 40

### 1.6 Hi-C Scaffolding

    Bowtie2: https://github.com/BenLangmead/bowtie2
    HiCUP: https://github.com/StevenWingett/HiCUP
    ALLHIC: https://github.com/tangerzhang/ALLHiC

    bowtie2-build typha_polish.fasta typha_polish
    hicup_digester --genome typha_polish --re DpnII typha_polish.fa
    hicup --config hicup.config
    ALLHiC_partition -b HiCUP.bam -r typha_polish.fasta -e GATC -k 15
    ALLHiC_rescue -b HiCUP.bam -r typha_polish.fasta -c clusters.txt -i counts_RE.txt
    allhic extract -b HiCUP.bam -r typha_polish.fasta -RE GATC
    allhic optimize group1.txt sample.clean.clm
    ...
    allhic optimize group15.txt sample.clean.clm
    ALLHiC build typha_polish.fasta
  #### here we got the groups.asm.fasta
  #### next we used juicer and 3D-DNA to further scaffold the genome

    3d-dna: https://github.com/aidenlab/3d-dna
    juicer: https://github.com/aidenlab/juicer
    BWA: https://github.com/lh3/bwa

    bwa index groups.asm.fasta
    python path/juicer/misc/generate_site_positions.py DpnII typha_allhic groups.asm.fasta
    awk 'BEGIN{OFS="\t"}{print $1, $NF}' typha_allhic_DpnII.txt > typha_allhic_DpnII.bed
    juicer.sh -g typha_allhic -d ./ -s DpnII -p typha_allhic_DpnII.bed -z groups.asm.fasta -D path/juicer -t 16
    run-asm-pipeline.sh groups.asm.fasta aligned/merged_nodups.txt

  #### after the 3D-DNA, we used Juicebox to manually correct the assembly errors
    Juicebox:https://github.com/aidenlab/Juicebox
  #### rerun the 3d-dna with the chromosome level draft assembly
    run-asm-pipeline-review.sh -r groups.asm.fasta aligned/merged_nodups.txt
#### the chrosome level assembly : final.assembly.fasta (Typha_laxmannii_genome_assembly.fasta)
#### we got the final chromosome level assembly *Typha laxmannii* genome and polish with the NGS and TGS reads
    NextPolish: https://github.com/Nextomics/NextPolish
    SAMtools: https://github.com/lh3/samtools
    BWA: https://github.com/lh3/bwa
        round=2
        threads=32
        read1=NGS_1.fq.gz
        read2=NGS_2.fq.gz
        input= final.assembly.fasta
        for ((i=1; i<=${round};i++)); do
        #step 1:
        #index the genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
        input=genome.polishtemp.fa;
        #step2:
        #index genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
        input=genome.nextpolish.fa;
        done;
        #Finally polished genome file: genome.nextpolish.fa

#### rename the final assembly
    mv genome.nextpolish.fa Typha_laxmannii.fasta

### 1.7 Assembly Quality Assessment
#### realign the NGS reads to the final assembly
    bwa index Typha_laxmannii.fasta
    bwa mem -t 32 Typha_laxmannii.fasta NGS_1.fq.gz NGS_2.fq.gz | samtools view -bS - | samtools sort -o Typha_NGS.bam
#### count the mapping rate
    samtools flagstat Typha_NGS.bam > Typha_NGS.flagstat
    # the mapping rate is 99.57%

#### Estimate the assembly quality with BUSCO
    busco -i Typha_laxmannii.fasta -l embryophyta_odb10 -o busco_final -m genome -c 40

## 2. Genome Annotation
### materials
  - Final genome assembly: Typha_laxmannii.fasta
  - RNA-seq data:
    ROOT: root1_1.fq.gz, root1_2.fq.gz, root2_1.fq.gz, root2_2.fq.gz
    LEAF: leaf1_1.fq.gz, leaf1_2.fq.gz, leaf2_1.fq.gz, leaf2_2.fq.gz
    STEM: stem1_1.fq.gz, stem1_2.fq.gz, stem2_1.fq.gz, stem2_2.fq.gz
    Shoot: shoot1_1.fq.gz, shoot1_2.fq.gz, shoot2_1.fq.gz, shoot2_2.fq.gz
  - Protein sequences from related species for homology-based prediction:
    Oryza sativa (accession no. GCF_001433935.1), Zea mays (no. GCF_902167145.1), Brachypodium distachyon (no. GCF_000005505.3), and Ananas comosus (no. GCF_001540865.1).
### 2.1 Repeat Annotation
    RepeatModeler: https://github.com/Dfam-consortium/RepeatModeler
    RepeatMasker: https://github.com/rmhubley/RepeatMasker
    LTR_FINDER_parallel: https://github.com/oushujun/LTR_FINDER_parallel
    LTR_retriever: https://github.com/oushujun/LTR_retriever
    CD-HIT https://github.com/weizhongli/cdhit-web-server
    BuildDatabase -name typha -engine ncbi typha_draft.fa
    RepeatModeler -database typha -engine ncbi -pa 40
    RepeatMasker -pa 40 -s -lib path/typha-families.fa typha_draft.fa
    LTR_FINDER_parallel -seq typha_draft.fa -size 1000000 -time 300 -threads 40 > ltr_finder.out
    LTR_retriever -genome typha_draft.fa -inharvest ltr_finder.out -threads 40
    cat consensi.fa LTRlib.fa > repeats_combined.fa
    cd-hit-est -i repeats_combined.fa -o repeats_nr.fa -c 0.95 -n 10 -T 40 -M 0
    RepeatMasker -pa 40 -s -lib repeats_nr.fa typha_draft.fa

### 2.2 Gene Structure Prediction
### prepared files
    Typha_laxmannii.fasta.masked (from RepeatMasker)
    RNA-seq data (from fastp)
    Protein sequences from related species (from NCBI)
### 2.2.1 RNA-seq based prediction
    StringTie: https://github.com/gpertea/stringtie
    HISAT2: https://github.com/DaehwanKimLab/hisat2
    Samtools: https://github.com/samtools/samtools
    TransDecoder: https://github.com/TransDecoder/TransDecoder
    hisat2-build Typha_laxmannii.fasta.masked Typha_index
    hisat2 -p 40 -x Typha_index -1 root1_1.fq.gz -2 root1_2.fq.gz | samtools view -bS - | samtools sort -o root1.bam
    hisat2 -p 40 -x Typha_index -1 root2_1.fq.gz -2 root2_2.fq.gz | samtools view -bS - | samtools sort -o root2.bam
    hisat2 -p 40 -x Typha_index -1 leaf1_1.fq.gz -2 leaf1_2.fq.gz | samtools view -bS - | samtools sort -o leaf1.bam
    hisat2 -p 40 -x Typha_index -1 leaf2_1.fq.gz -2 leaf2_2.fq.gz | samtools view -bS - | samtools sort -o leaf2.bam
    hisat2 -p 40 -x Typha_index -1 stem1_1.fq.gz -2 stem1_2.fq.gz | samtools view -bS - | samtools sort -o stem1.bam
    hisat2 -p 40 -x Typha_index -1 stem2_1.fq.gz -2 stem2_2.fq.gz | samtools view -bS - | samtools sort -o stem2.bam
    hisat2 -p 40 -x Typha_index -1 shoot1_1.fq.gz -2 shoot1_2.fq.gz | samtools view -bS - | samtools sort -o shoot1.bam
    hisat2 -p 40 -x Typha_index -1 shoot2_1.fq.gz -2 shoot2_2.fq.gz | samtools view -bS - | samtools sort -o shoot2.bam
    samtools merge -@ 40 root1.bam root2.bam leaf1.bam leaf2.bam stem1.bam stem2.bam shoot1.bam shoot2.bam > rna_seq.bam
    stringtie -p 40 -o rna_seq.gtf -l Typha rna_seq.bam
    gtf_genome_to_cdna_fasta.pl rna_seq.gtf Typha_laxmannii.fasta.masked > rna_seq_transcripts.fa
    gtf_to_alignment_gff3.pl rna_seq.gtf > rna_seq.gff3
    TransDecoder.LongOrfs -t rna_seq_transcripts.fa
### 2.2.2 Homology-based prediction
    BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
    Exonrate:
    cat Oryza_sativa.fasta Zea_mays.fasta Brachypodium_distachyon.fasta Ananas_comosus.fasta > homology_proteins.fa
    makeblastdb -in homology_proteins.fa -dbtype prot -out homology_proteins_db
    blastx -query Typha_laxmannii.fasta.masked -db homology_proteins_db -outfmt 6 -evalue 1e-5 -num_threads 40 -out blastx.out
    exonerate --model protein2genome --softmasktarget -showtargetgff yes --percent 70 --geneseed 200 --minintron 50 --maxintron 2000 --fsmmemory 1000 homology_proteins.fa Typha_laxmannii.fasta.masked > Typha_exonerate.gff3
#### 2.2.3 Ab initio prediction

    GeneMark-ES: http://exon.gatech.edu/GeneMark/
    AUGUSTUS: https://github.com/Gaius-Augustus/Augustus

    gmes_petap.pl --sequence Typha_laxmannii.fasta.masked --ES --cores 40

    # Output: genemark.gtf (gene predictions), output directory with trained model

    # ==============================
    # AUGUSTUS (train and predict)
    # ==============================
    # Create new species profile
    new_species.pl --species=typha_laxmannii

    # Train using GeneMark predictions (optional: high-confidence set)
    etraining --species=typha_laxmannii genemark.gtf

    # Predict with trained model
    augustus --species=typha_laxmannii Typha_laxmannii.fasta.masked > augustus.gff3

#### 2.2.4 Integrate all the evidences to get the final gene models
    Maker: http://www.yandell-lab.org/software/maker.html  



### 2.3 Functional Annotation
    SiwssProt: https://www.uniprot.org/
    NR: https://ftp.ncbi.nlm.nih.gov/blast/db/
    KEGG: https://www.genome.jp/kegg/
    GO: http://geneontology.org/
    Pfam: http://pfam.xfam.org/
    InterPro: https://www.ebi.ac.uk/interpro/
    Diamond: https://github.com/bbuchfink/diamond
    KOBAS: http://kobas.cbi.pku.edu.cn/
    InterProScan: https://github.com/ebi-pf-team/interproscan
    hmmscan: http://hmmer.org/

    diamond makedb --in swissprot.fasta -d swissprot
    diamond blastp -d swissprot -q typha_proteins.fa -o diamond_swissprot.out -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue 1e-5 --threads 40

    diamond makedb --in nr.fasta -d nr
    diamond blastp -d nr -q typha_proteins.fa -o diamond_nr.out -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue 1e-5 --threads 40

    diamond makedb --in go.fasta -d go
    diamond blastp -d go -q typha_proteins.fa -o diamond_go.out -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue 1e-5 --threads 40

    diamond makedb --in kegg.fasta -d kegg
    diamond blastp -d kegg -q typha_proteins.fa -o diamond_kegg.out -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --evalue 1e-5 --threads 40

    # annotate the kegg pathways with kobas web server

    interproscan.sh -i typha_proteins.fa -f tsv,xml -o interproscan.out -dp -T ./temp -cpu 40

    hmmscan --cpu 40 --domtblout pfam.out Pfam-A.hmm typha_proteins.fa > pfam.log

    # integrate all the functional annotation results

### 2.4 Non-coding RNA Annotation
    tRNAscan-SE: http://lowelab.ucsc.edu/tRNAscan-SE/  
    RNAmmer: http://rna.informatik.uni-freiburg.de/RNA.html
    Infernal: http://eddylab.org/infernal/
    Rfam: http://rfam.xfam.org/
    tRNAscan-SE -B -o typha_tRNAs.out Typha_laxmannii.fasta.masked
    rnammer -S euk -m tsu,ssu,lsu -gff typha_rRNAs.gff Typha_laxmannii.fasta.masked
    cmsearch --cpu 40 --tblout typha_ncRNA.out Rfam.cm Typha_laxmannii.fasta.masked > typha_ncRNA.log

### 2.5 Genome analysis

#### 2.5.1 species phylogenetic tree
#### Species
| Assembly  | Species   | Source |
| :---: | :---: | :---: |
|GCF_004353265.1    |*Vutis ripapria* |https://www.ncbi.nlm.nih.gov/assembly/GCF_004353265.1/|
|GCA_902729315.2	|*Spirodela intermedia*   |https://www.ncbi.nlm.nih.gov/assembly/GCA_902729315.2/|
|GCF_000001605.2	|*Sorghun bicolor*    |https://www.ncbi.nlm.nih.gov/assembly/GCF_000003195.3/|
|GCF_001263595.1	|*Phalaenopsis equestris*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001263595.1/|
|GCF_001433935.1	|*Oryza sativa*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001433935.1/|
|GCF_000313855.2	|*Musa acuminata*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000313855.2/|
|GCF_000442705.1	|*Elaeis guineensis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000442705.1/|
|GCF_001876935.1	|*Asparagus officinalis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_001876935.1/|
|GCF_001540865.1	|*Ananas comosus* |https://www.ncbi.nlm.nih.gov/assembly/GCF_001540865.1/|
|GCF_000130695.1	|*Amborella trichopoda*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_000130695.1/|
|GWHCBIL00000000	|*Typha latifolia*   |https://download.cncb.ac.cn/gwh/Plants/Typha_latifolia_Typha_latifolia_jzg2022_GWHCBIL00000000/GWHCBIL00000000.genome.fasta.gz|
|GCF_048772165.1  |*Typha angustiflia*   |https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048772165.1/|
#### build the tree with OrthoFinder
    OrthoFinder:https://github.com/davidemms/OrthoFinder
    mafft: https://mafft.cbrc.jp/alignment/software/
    BEAST: https://beast.community/
#### prepare files
    typha_laxmannii.pep
    typha_angustifolia.pep
    typha_latifolia.pep
    Vutis.pep
    Spirodela.pep
    Sorghun.pep
    Phalaenopsis.pep
    Oryza.pep
    Musa.pep
    Elaeis.pep
    Asparagus.pep
    Ananas.pep
    Amborella.pep
    
    orthofinder -f ./ -t 40 -S diamond

    #### tree file is in the OrthoFinder/Results_data/WorkingDirectory/Species_Tree/ directory
    cp OrthoFinder/Results_data/Singe_Copy_Orthologue_Sequences/* ./data/
    for file in $(ls data/*.fa);do
        mafft --auto $file > $file.mafft.fas
    done
    python3 merge_species.py data merged.fas
    mafft --auto merged.fas > merged.fas.mafft.fas
    #  we use the beast to estimate the divergence time between species.

### 2.6 WGD analysis
    MCScanX:   https://github.com/wyp1125/MCScanX
    ParaAT: https://github.com/wonaya/ParaAT
    Kaks_Calculator: https://sourceforge.net/projects/kakscalculator2/
## prepare files
    typha.pep
    typha.gff3
    typha_final.gff3
    typha_mcscanx.gff
    axt_cal.py 
    
    makeblastdb -in typha.pep -dbtype prot -out typha.pep
    blastp -query typha.pep -db typha.pep -outfmt 6 -evalue 1e-10 -num_threads 40 -out typha.pep
    cp typha_mcscanx.gff typha.gff
    MCScanX typha.gff
    cat typha.gff.collinearity | rep "Chr" | awk '{print $3"\t"$4}'> typha.homolog
    cat 40 > threads
    perl path/ParaAT.pl -h typha.homolog -n typha.pep -a typha.pep -p threads -m clustalw -f axt -o typha_out
    for i in `ls *.axt`;do KaKs_Calculator -i $i -o ${i}.kaks -m YN;done
    for i in `ls *.axt`;do python axt_cal.py $i ${i}.one-line;done
    for i in `ls *.kaks`;do awk 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' $i >>all-kaks.txt;done
    sort all-kaks.txt|uniq >all-kaks.results
#### here we got the all-kaks.results file, by this method, we can get the WGD events in *Oryza sativa*, *Sorghun bicolor*, and *Ananas comosus*, *Typha latifolia*, *Typha angustifolia* and *Typha laxmannii* genomes.


## 3. Population Genomics Analysis
### 3.1 Materials
  - Resequencing data of 126 individuals from 31 populations across the Tianshan mountain region
  - Reference genome: Typha_laxmannii.fasta
  - Outgroup genome: Typha angustifolia (GCF_048772165.1)

### 3.2 Build the index of the reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    GATK4: https://github.com/broadgsa/gatk
    bwa index -p typha Typha_laxmannii.fasta
    samtools faidx Typha_laxmannii.fasta
    java -Djava.io.tmpdir=/tmp -Xmx32g -jar path/gatk-package-4.2.4.1-local.jar CreateSequenceDictionary -R typha_laxmannii.fa -O typha.dict 
  
#### 3.3 quality control of the raw data
    fastp: https://github.com/OpenGene/fastp
    for sample in $(cat sample_list.txt);
    do 
    fastp -i ${sample}_1.fq.gz -I ${sample}_2.fq.gz -o ${sample}_1.fq.gz -O ${sample}_2.fq.gz -h ${sample}.html -j ${sample}.json -w 8 -q 20 -u 20 -n 5 -l 50 -y 20 -x --umi --umi_loc read1 --umi_len 10;
    done 

### 3.4 Mapping the reads to the reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    for sample in $(cat sample_list.txt);
    do bwa mem -t 8 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" typha ${sample}_1.fq.gz ${sample}_2.fq.gz | samtools view -bS | samtools sort -o ${sample}.bam;
    samtools index ${sample}.bam;
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}_dedup.metrics.txt;
    samtools index ${sample}_dedup.bam;
    done

### 3.5 Variant Calling
    GATK4:  https://github.com/broadgsa/gatk
    for sample in $(cat sample_list.txt);
    do java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar HaplotypeCaller --native-pair-hmm-threads 16 -R typha_laxmannii.fa -I ${sample}_dedup.bam -ERC GVCF -O ${sample}.g.vcf.gz;
    done

### 3.6 Joint Genotyping
    GATK4:  https://github.com/broadgsa/gatk
    for file in $(cat list);do echo -en --variant ${file}.g.vcf' '; done > mergevcf
    echo  >> mergevcf
    cat mergevcf | while read line;do
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar CombineGVCFs -R typha_laxmannii.fa $line -O merge.g.vcf.gz
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar GenotypeGVCFs -R typha_laxmannii.fa -V merge.g.vcf.gz -O typha_raw.vcf.gz
    done
#### here we go the raw variants vcf file, and we need to do the filter of the variants
    typha_raw.vcf.gz

### 3.7 Variant Filtering
    GATK4: https://github.com/broadgsa/gatk
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantsToTable -R typha_laxmannii.fa -V typha_raw.vcf.gz -O typha_raw_variants.csv -F CHROM -F POS -F DP -F AF -F QUAL -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantFiltration -R typha_laxmannii.fa -V typha_raw.vcf.gz -O typha_raw_recode.vcf.gz --filter-expression "DP < 30||DP>8000 || QD < 15.0 || MQ < 58.0 || FS > 2.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -4.0" --filter-name "my_snp_filter"
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar SelectVariants -R typha_laxmannii.fa -V typha_raw_recode.vcf.gz  -O typha_raw_filter_variants.vcf.gz --exclude-filtered true
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar SelectVariants -R typha_laxmannii.fa -V typha_raw_filter_variants.vcf.gz --select-type-to-include SNP -O typha_raw_filter_variants_snp.vcf.gz
    # filter the SNPs with missing rate > 10% and MAF < 0.01 and MAF > 0.99
    vcftools --gzvcf typha_raw_filter_variants_snp.vcf.gz --max-missing 0.9 --min-alleles 2 --max-alleles 2 --maf 0.01 --max-maf 0.99 --recode --recode-INFO-all --out typha_final_filter_variants_snp

#### here we got the final SNPs file typha_final_filter_variants_snp.recode.vcf

### 3.8 group classification with pca
    #### here we use the Rscript to do the PCA analysis, please see 3_8_pca.r in the scripts directory.


### 3.9 build the phylogenetic tree of whole genome
    vcf2phylip: https://github.com/edgardomortiz/vcf2phylip
    IQ-TREE: https://github.com/Cibiv/IQ-TREE
    python /path/vcf2phylip.py -i typha_final_filter_variants_snp.recode.vcf -o typha_filter_variants_snp.phy
    iqtree -s typha_filter_variants_snp.phy -m TVMe+R6 -bb 1000 -nt 40

#### here we got the phylogenetic tree file typha_filter_variants_snp.phy.treefile, and we visualize the tree with itol (https://itol.embl.de/)

### 3.10 admixture analysis
    admixturePipeline: https://github.com/stevemussmann/admixturePipeline
    admixturePipeline.py -m popmap.txt -v typha_final_filter_variants_snp.recode.vcf -k 2 -K 6 -R 8 -n 48
    ## here we use the CLUMPAK pipeline to basicly visualize the admixture results, and than use python script to plot the result

### 3.11 Infer the effective population size (Ne) history with PSMC
    PSMC: https://github.com/lh3/psmc
    for i in $(cat list)
    do
    samtools mpileup -C50 -uf typha_laxmannnii.fa ${i}_markdup.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 200 | gzip > ${i}.fq.gz
    fq2psmcfa -q20 ${i}.fq.gz > ${i}.psmcfa
    psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o ${i}.psmc ${i}.psmcfa
    python script/3_11_plot_psmc.py -namelist popmap.txt -psmcdir psmc_data -o psmc_result.svg

### 3.12 demographic history with momi2
    # we first construct the Fst-based NJ tree with the 126 individuals in five regional groups
    # we have the pop classification in five regions in the popmap_five_group.txt
    # generate the pairwise Fst matrix, and buid the NJ tree in Rscript (see 3_12_fst_nj_tree.r in the scripts directory)
    # the final Fst matrix was visualized with python script (see 3_12_fst_heatmap.py in the scripts directory)
    # the NJ tree was visualized with itol (https://itol.embl.de/) without the length information

    # based this Fst supported NJ tree as the topology structure, we used the momi2 to infer the demographic history of the five regional groups
    momi2: https://github.com/popgenmethods/momi2
    momi.ipynb
    # and we plot the demographic history with python script (see 3_12_momi2_plot.py in the scripts directory)

### 3.13 Estimated effective migration surface (EEMS)
    Ref：dipetkov/eems: Estimating Effective Migration Surfaces
    runeems_snp requires three main input files:
    *.diffs genetic distance matrix
    *.coord sampling point coordinates
    *.outer habitat coordinates
    
    #params-chain1.ini
    datapath = ./xiangpu
    mcmcpath = ./xiangpu-eems_result 
    nIndiv = 126                
    nSites = 1735925
    nDemes = 800
    diploid = false
    numMCMCIter = 20000000
    numBurnIter = 1000000
    numThinIter = 9999

    # run
    ../syz_eems_src/runeems_snps --params ./params-chain1.ini



### 3.14 Niche analysis
    Maxent: https://biodiversityinformatics.amnh.org/open_source/maxent/
    wordclim bio data: http://www.worldclim.com/past; https://www.worldclim.com/paleo-climate1;
    Last inter-glacial (LIG; ~120,000 - 140,000 years BP) : http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/pst/lig/lig_30s_bio.zip
    Last glacial maximum (LGM; ~21,000 years BP):http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/pst/21k/wc_2_5m_CCSM_21k_bio.zip
    Mid-Holocene (~6000 BP):http://biogeo.ucdavis.edu/data/climate/cmip5/mid/cnmidbi_2-5m.zip
    Current climate: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip
    # we preformance the niche build with the maxent in the arcgis software.
### 3.15 Environmental factor analysis

    # we first extract the enviromental date of sites, and performance the t-test between North and South groups with python scripts (see 3_14_t_test4envs.py in the scripts directory)
    # then we do the pca analysis of the enviromental factors with python scripts (see 3_14_envs_pca.py in the scripts directory)


### 3.16 Adaptation-related genes and regions identification
    # we first perform the genome-wide scan, to find identify the similarly genomic regions, and differently genomic regions between North and South groups
    vcftools --vcf North_clean.vcf --out North_pi --window-pi 50000
    vcftools --vcf South_clean.vcf --out South_pi --window-pi 50000
    vcftools --vcf North_clean.vcf --out North_tajimaD --TajimaD 50000
    vcftools --vcf South_clean.vcf --out South_tajimaD --TajimaD 50000
    vcftools --vcf typha_final_filter_variants_snp.recode.vcf --weir-fst-pop sample_S.txt --weir-fst-pop sample_N.txt --out fst_50kb --fst-window-size 50000 --keep sample_S.txt --keep sample_N.txt

    # we visualize the pi, tajimaD and Fst results with the python scripts (see 3_16_plot_pi_tajimaD_fst.py in the scripts directory)

    # Next, we performance the positive selection and selective sweep analysis with RAiSd
    RAiSD: https://github.com/alachins/raisd
    RAiSD -n North_clean -I North_clean.vcf -O North_clean -R
    RAiSD -n South_clean -I South_clean.vcf -O South_clean -R
    
    # we visualize the RAiSD results with the python scripts (see 3_16_plot_raid.py in the scripts directory)

    # here we find two large blocks for parallel sweep on chromosome 5, and genome-wide differentiation regions between North and South groups on chromosome 6.
    # and identitied the core sweep region, genes both in parallel sweep regions and differentiation regions with python scripts
    # see 3_16_sweep_region_chr5.py and 3_16_diff_region_chr6.py in the scripts directory


### 3.17 GWAS-based North-South differentiation related SNPs identification
    vcf2gwas:https://github.com/frankvogt/vcf2gwas
    vcf2gwas -v typha_final_filter_variants_snp.recode.vcf - -pf North_South.csv -ap -lmm -T 6 -nl
    # We visualize the result and anotated the genes with the python scripts (see 3_17_plot_GWAS.py)

