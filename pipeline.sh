#Install SAW pipeline via Singularity
singularity build SAW_v4.1.0.sif docker://stomics/saw:04.1.0

#Prepare input files for building reference
referenceDir=/data1/dansonloi/STOmics/reference
refName=MusMusculus
# make a new directory called reference specie name
mkdir -p $referenceDir/$refName
cd $referenceDir/$refName
mkdir genome genes

# download the genome fasta file, like "Genome sequence, primary assembly (GRCh38)", from Ensembl
cd genome
wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #wget https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# download the annotation file that correspond to your reference genome (.gtf/.gff) from Ensembl
cd ../genes
wget http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz #wget https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz 
gunzip Mus_musculus.GRCm39.107.gtf.gz

#build the genome index file
referenceDir=/data1/dansonloi/STOmics/reference
export SINGULARITY_BIND=$referenceDir
mkdir $referenceDir/STAR_SJ100

singularity exec ~/STOmics/SAW_v4.1.0.sif mapping \
    --runMode genomeGenerate \
    --genomeDir $referenceDir/STAR_SJ100 \
    --genomeFastaFiles $referenceDir/MusMusculus/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \ #changed to GRCm38 build
    --sjdbGTFfile $referenceDir/MusMusculus/genes/Mus_musculus.GRCm38.102.gtf \ #changed to GRCm38 build
    --sjdbOverhang 99 \
    --runThreadN 12

#GRCm38 build
wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz

#Prepare mapping --stParaFile input file {lane}.bcPara
#This mapping pipeline only supports one lane each time, change the "lane" variable below
workspace=/data1/dansonloi/STOmics
mkdir $workspace/00.mapping
lane=V350080546_L01
SN=SS200000954BR_D4
vim $workspace/00.mapping/$lane.bcPara
#Copy and paste the codes below into the .bcPara file
in=workspace/mask/SN.barcodeToPos.h5
in1=workspace/20220726_V350080546_raw/lane_read_1.fq.gz
in2=workspace/20220726_V350080546_raw/lane_read_2.fq.gz
encodeRule=ACTG
out=lane
barcodeReadsCount=workspace/00.mapping/lane.barcodeReadsCount.txt
action=4
platform=T10
barcodeStart=0
barcodeLen=25
umiStart=25
umiLen=10
umiRead=1
mismatch=1
#press esc & enter ":x" to save and exit vim 
#specify lane/chip numbers 
sed -i "s/lane/$lane/g" $workspace/00.mapping/$lane.bcPara
sed -i "s/SN/$SN/g" $workspace/00.mapping/$lane.bcPara
sed -i "s|workspace|$workspace|g" $workspace/00.mapping/$lane.bcPara #| is used instead of / to avoid confusion with the delimiter 
cat $workspace/00.mapping/$lane.bcPara #optional: check to see if lane/chip numbers are specified 

#vimout=/data1/dansonloi/STOmics/00.mapping
#fastq=/data1/dansonloi/STOmics/20220726_V350080546_raw
#cd $fastq
#for i in `ls *_read_1.fq.gz`
#do
#	vim $vimout/${i%_read_1.fq.gz}.bcPara
#done

#Run mapping pipeline
referenceDir=/data1/dansonloi/STOmics/reference
cd $workspace
singularity exec SAW_v4.1.0.sif mapping \
                --outSAMattributes spatial \
                --outSAMtype BAM SortedByCoordinate \
                --genomeDir $referenceDir/STAR_SJ100 \
                --runThreadN 8 \
                --outFileNamePrefix $workspace/00.mapping/$lane. \
                --stParaFile $workspace/00.mapping/$lane.bcPara \
                --readNameSeparator \" \" \
                --limitBAMsortRAM 38582880124 \
                --limitOutSJcollapsed 10000000 \
                --limitIObufferSize=280000000 \
                --outBAMsortingBinsN 50 \
                > $workspace/00.mapping/${lane}_barcodeMap.stat

#merge barcodeReadsCount.txt files of all 4 lanes
mkdir $workspace/01.merge
lane=V350080546_L01
singularity exec SAW_v4.1.0.sif merge \
              -i $workspace/00.mapping/${lane}.barcodeReadsCount.txt,$workspace/00.mapping/${lane%1}2.barcodeReadsCount.txt,$workspace/00.mapping/${lane%1}3.barcodeReadsCount.txt,$workspace/00.mapping/${lane%1}4.barcodeReadsCount.txt \
              --out $workspace/01.merge/${SN}.barcodeReadsCount.txt \
              --action 2

#Run count pipeline
mkdir -p $workspace/02.count
geneExp=$workspace/02.count/${SN}.raw.gef
saturationSamplingFile=$workspace/02.count/$SN_raw_barcode_gene_exp.txt
gtf=/data1/dansonloi/STOmics/reference/MusMusculus/genes/Mus_musculus.GRCm39.107.gtf
singularity exec SAW_v4.1.0.sif count \
        -i $workspace/00.mapping/${lane}.Aligned.sortedByCoord.out.bam,$workspace/00.mapping/${lane%1}2.Aligned.sortedByCoord.out.bam,$workspace/00.mapping/${lane%1}3.Aligned.sortedByCoord.out.bam,$workspace/00.mapping/${lane%1}4.Aligned.sortedByCoord.out.bam \
        -o $workspace/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
        -a $gtf \
        -s $workspace/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -e ${geneExp} \
        --umi_len 10 \
        --sat_file ${saturationSamplingFile} \
        --sn ${SN} \
        --umi_on \
        --save_lq \
        --save_dup \
        -c 8 \
        -m 128

#Run register pipeline
image=$workspace/ImageQC
image4register=$(find ${image} -maxdepth 1 -name *.tar.gz | head -1)
imageQC=$(find ${image} -maxdepth 1 -name *.json | head -1)
mkdir -p $workspace/03.register
singularity exec SAW_v4.1.0.sif register \
                      -i ${image4register} \
                      -c ${imageQC} \
                      -v $workspace/02.count/${SN}.raw.gef \
                      -o $workspace/03.register

#Run tissueCut pipeline
mkdir -p $workspace/04.tissuecut
singularity exec SAW_v4.1.0.sif tissueCut \
        --dnbfile $workspace/01.merge/${SN}.barcodeReadsCount.txt \
        -i $workspace/02.count/${SN}.raw.gef \
        -o $workspace/04.tissuecut \
        -s $workspace/03.register/7_result \
        -t tissue \
        --snId ${SN} \
        --platform T10

#spatialCluster
mkdir -p $workspace/05.spatialcluster
singularity exec SAW_v4.1.0.sif spatialCluster \
    -i $workspace/04.tissuecut/${SN}.tissue.gef \
    -o $workspace/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -s 200

#saturation
mkdir -p $workspace/06.saturation
lane=V350080546_L01
singularity exec SAW_v4.1.0.sif saturation \
	-i $workspace/02.count/raw_barcode_gene_exp.txt \
	--tissue $workspace/04.tissuecut/${SN}.tissue.gef \
	-o $workspace/06.saturation \
	--bcstat $workspace/00.mapping/${lane}_barcodeMap.stat,$workspace/00.mapping/${lane%1}2_barcodeMap.stat,$workspace/00.mapping/${lane%1}3_barcodeMap.stat,$workspace/00.mapping/${lane%1}4_barcodeMap.stat \
	--summary $workspace/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat

#report
mkdir -p $workspace/07.report
singularity exec SAW_v4.1.0.sif report \
    -m $workspace/00.mapping/${lane}_barcodeMap.stat,$workspace/00.mapping/${lane%1}2_barcodeMap.stat,$workspace/00.mapping/${lane%1}3_barcodeMap.stat,$workspace/00.mapping/${lane%1}4_barcodeMap.stat \
    -a $workspace/00.mapping/${lane}.Log.final.out,$workspace/00.mapping/${lane%1}2.Log.final.out,$workspace/00.mapping/${lane%1}3.Log.final.out,$workspace/00.mapping/${lane%1}4.Log.final.out \
    -g $workspace/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -l $workspace/04.tissuecut/tissuecut.stat \
    -n $workspace/04.tissuecut/${SN}.gef \
    -d $workspace/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -t $workspace/06.saturation/plot_200x200_saturation.png \
    -b $workspace/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
    -v $workspace/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
    -r standard_version \
    -i $workspace/04.tissuecut/tissue_fig/${SN}.ssDNA.rpi \
    -s ${SN} \
    --pipelineVersion SAW_v4.1.0 \
    -o $workspace/07.report

