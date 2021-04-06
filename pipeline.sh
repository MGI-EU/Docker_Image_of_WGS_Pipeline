#!/bin/env sh

#*****************************************************************************************
# FileName: pipeline.sh
# Creator: <chengbochen@genomics.cn>
# Date: 2020-12-11

# Description: Whole-genome sequencing data analysis 
#              reference: GRCh37/hg19
#              call genetic variants by DeepVariant
# CopyRight: 
# version: 0.1
#*****************************************************************************************

function usage {
    echo -e "    USAGE:sh $0 -d file_database -o path2outdir -n name_outdir -p image_pipeline -v image_deepvariant"
    echo -e "================================================"
    echo -e "    Example:sh $0 -d DATABASE_DIR -o OUTPUT_DIR -n RESULT_DIR_NAME -m PIPELINE_IMAGE_NAME -v DEEPVARIANT_IMAGE_NAME\n"
    echo -e "    option:"
    echo -e "           -d path of reference and config file"
    echo -e "           -o path to output directory"
    echo -e "           -n name of output directory"
    echo -e "           -p name of pipeline image"
    echo -e "           -v name of deepvariant image\n"
}

if [ $# != 10 ]
then
    usage
    exit
fi

while getopts ":d:o:n:m:" opt
do
    case $opt in
        d)
        DATABASE=$OPTARG
        ;;
        o)
        OUTPUT_DIR=$OPTARG
        ;;
        n)
        DIR_NAME=$OPTARG
        ;;
        p)
        PIPELINE_IMAGE_NAME=$OPTARG
        ;;
        v)
        DEEPVARIANT_IMAGE_NAME=$OPTARG
        ;;
        ?)
        echo "Unknown parameter"
        exit 1;;
    esac
done

#if [ $# != 8 ] ; then
#    echo "USAGE: $0 TABNAME"
#    echo " e.g.: sh $0 -d DATABASE_DIR -o OUTPUT_DIR -n RESULT_DIR_NAME -m PIPELINE_PIPELINE_IMAGE_NAME"
#    exit 1;
#fi

#Generate WGS pipeline script
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} perl /opt/pipeline/Pipeline_WGS_BGISEQ.bash2.pl /opt/pipeline/database/sample.info /opt/pipeline/database/docker.config.txt /opt/pipeline/result/${DIR_NAME}
#sudo docker run -it -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME}

## cleaning raw fq files
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/clean.sh

## aligning using BWA
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/bwa.sh

## checking the results of BWA
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/bwaCheck.sh 

## sorting bam
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/sortBam.sh

## checking sorted BAM
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/sortCheck.sh 

## merging bam
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/mergeBam.sh

## marking duplicates in BAMs
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/mkdup.sh

## Bam coverage statistics
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/covdep.sh 

## checking duplicates-marked BAM
sudo docker run -v ${DATABASE}:/opt/pipeline/database -v ${OUTPUT_DIR}:/opt/pipeline/result ${PIPELINE_IMAGE_NAME} sh /opt/pipeline/result/${DIR_NAME}/shell_run/mkdupCheck.sh 

## SNP calling by Deepvariant
for i in `ls ${OUTPUT_DIR}/${DIR_NAME}/*/mkdup-realn/*.bam`
do
BAM_FILE=`basename $i`;
SAMPLE=${BAM_FILE%%.*}
#echo ${SAMPLE}
#echo "${i%/*}"
#sudo /usr/bin/docker pull google/deepvariant:"0.10.0-gpu"
sudo /usr/bin/docker run -m 20480M -c 4  --privileged=true -v "${i%/*}":"/input":rw -v "${i%/*}/../SNP":"/output":rw -v "${DATABASE}/reference/hg19":"/reference":rw  ${DEEPVARIANT_IMAGE_NAME} /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/reference/hg19.fa --reads=/input/${BAM_FILE}  --output_vcf=/output/${SAMPLE}.vcf.gz --output_gvcf=/output/${SAMPLE}.g.vcf.gz --num_shards=4;
done
