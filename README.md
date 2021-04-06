# Docker_Image_of_WGS_Pipeline

### Introduction

This is used to build the WGS pipeline docker image.  

It can also be used to build local images, especially if you've made changes to the code.  

Ensure that docker is installed locally.

The download links reference genome (database file) is:  
https://bgitech-my.sharepoint.com/:u:/g/personal/chengbochen_genomics_cn/EaLtQNezQFJLojuoP9opw78BYFDUpVMftT-XIeFo6RsM5A?e=IKyMbI

The download links reference genome (software) is: https://pan.genomics.cn/ucdisk/s/Urmq6b


Example command:
```
unzip WGS_pipeline.zip
cd WGS_pipeline
sudo docker build -t wgspipeline
```

After finishing the installation, you can get the wgspipeline image locally

Use the command ```sudo docker images``` to view

#### 1.1 Required software

(1) Docker

(2) DeepVariant


The following software has been installed in the pipeline image：

SOAPnuke (v1.5.6)

bwa (v0.7.17)

gnuplot (v4.6.7)

java (v1.8.0_271)

perl (v5.16.3)

picard (v2.10.9)

samtools (v1.11) (using htslib 1.11)

R (v3.6.0)

#### 1.2 installation

Use docker to import this image file:

```sudo /usr/bin/docker load -i wgspipeline.tar```

install deepvariant :

```sudo /usr/bin/docker pull google/deepvariant```


#### 1.3 Configuration file


**pipeline.sh** : The script to run the pipeline.

Files in the database folder：

**reference** : Reference genome

**TMP** : Fq files of each sample

**docker.config.txt** : docker configuration file

**singularity.config.txt** : singularity configuration file

**sample.info** : information of samples


#### 1.4 run

Put the docker.config.txt, singularity.config.txt, sample.info in the database folder

Put the TMP folder in the database folder

```shell
 USAGE:sh pipeline.sh -d file_database -o path2outdir -n name_outdir -p image_pipeline -v image_deepvariant
==================================================================================================
    Example:sh pipeline.sh -d DATABASE_DIR -o OUTPUT_DIR -n RESULT_DIR_NAME -m PIPELINE_IMAGE_NAME -v DEEPVARIANT_IMAGE_NAME

    option:
           -d path of reference and config file
           -o path to output directory
           -n name of output directory
           -p name of pipeline image
           -v name of deepvariant image
```

**1.5 Result path**

OUTPUT_DIR/RESULT_DIR_NAME
