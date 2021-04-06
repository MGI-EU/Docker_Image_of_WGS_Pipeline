#This is used to build the WGS pipeline docker image
#Auother: chengbochen@genomics.cn
#Example command:
#sudo docker build -t wgspipeline .

#  Centos 7 base image
FROM centos:7

RUN echo "bash;" >> /docker-entrypoint.sh
ENTRYPOINT ["sh", "/docker-entrypoint.sh"]
#ENV HOME /opt
WORKDIR /opt

# Copy the Dockerfile to root for reference
ADD Dockerfile /opt/.Dockerfile

# Update and software install
RUN yum update -y && \
 yum install -qy wget perl epel-release bzip2 ncurses-libs ncurses-devel ImageMagick && \
 yum install -qy R && \
 mkdir /opt/pipeline && \
 mkdir /opt/pipeline/bin && \
 mkdir /opt/pipeline/database && \
 mkdir /opt/pipeline/lib && \
 mkdir /usr/local/software && \
 mkdir /usr/local/software/SOAPnuke && \
 cd /usr/local/software/SOAPnuke && \
# SOAPnuke install
 wget https://github.com/BGI-flexlab/SOAPnuke/archive/1.5.6-linux.zip && \
 unzip 1.5.6-linux.zip && \
 #cp -r /usr/local/software/SOAPnuke/SOAPnuke-1.5.6-linux/bin/SOAPnuke /opt/pipeline/bin/SOAPnuke1.5.6 && \
 ln -s /usr/local/software/SOAPnuke/SOAPnuke-1.5.6-linux/bin/SOAPnuke /opt/pipeline/bin/SOAPnuke1.5.6 && \
 rm -rf 1.5.6-linux.zip && \
 cd /usr/local/software/ && \
# bwa install
 wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
 tar -jxvf bwa-0.7.17.tar.bz2 && \
 rm -rf bwa-0.7.17.tar.bz2 && \
 mv bwa-0.7.17 bwa && \
 cd /usr/local/software/bwa && \
 make && \
 mkdir /usr/local/software/picard && \
 cd /usr/local/software/picard && \
# picard install
 wget https://github.com/broadinstitute/picard/releases/download/2.10.9/picard.jar && \
 mkdir /usr/local/software/samtools && \
 cd /usr/local/software/samtools && \
# samtools install
 wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
 tar -jxvf samtools-1.11.tar.bz2 && \
 rm -rf samtools-1.11.tar.bz2 && \
 cd /usr/local/software/samtools/samtools-1.11 && \
 ./configure && \
 make 
# gnuplot install
ADD software/gnuplot-4.6.7.tar.gz /opt/pipeline/bin/
RUN cd /opt/pipeline/bin/gnuplot-4.6.7 && \
 ./configure --prefix=/usr/local/gnuplot && \
 make && \
 make install && \
 mkdir /usr/local/java
#java install
ADD software/jdk-8u271-linux-x64.tar.gz /usr/local/java/
RUN ln -s /usr/local/java/jdk1.8.0_271 /usr/local/java/jdk && \
 cd /opt/pipeline/bin && \
 ln -s /usr/local/software/bwa/bwa bwa && \
 ln -s /usr/local/gnuplot/bin/gnuplot gnuplot && \
 ln -s /usr/local/java/jdk1.8.0_271/bin/java java && \
 ln -s /usr/bin/perl perl && \
 ln -s /usr/local/software/samtools/samtools-1.11/samtools samtools && \
 ln -s /usr/local/software/picard/picard.jar /opt/pipeline/bin/picard-2.10.9.jar && \
 yum clean all && rm -rf /tmp/* /var/tmp/* && rm -rf /var/cache/yum/* 

#files in /bin
COPY software/bam2depth \
 software/fqcheck \
 software/fqcheck33 \
 scripts/SLURM.pl \
 scripts/SLURM.v2.pl \
 scripts/TEST.Phred64_33.pl \
 scripts/bam2Check.pl \
 scripts/bamCheck_v1.pl \
 scripts/depthV2.0.pl \
 scripts/fqcheck_distribute.pl \
 scripts/picard.pl \
 scripts/soapnuke_stat.pl /opt/pipeline/bin/
#config files
COPY config/docker.config.txt \
 config/Pipeline_WGS_BGISEQ.bash2.pl \
 config/Pipeline_WGS_BGISEQ.gpu.pl \
 config/sample.info \
 config/singularity.config.txt \
 config/test.sh /opt/pipeline/
#files in /lib
COPY lib/PLOT.pm /opt/pipeline/lib/
#env set
ENV JAVA_HOME=/usr/local/java/jdk 
ENV JRE_HOME=${JAVA_HOME}/jre 
ENV CLASSPATH=.:${JAVA_HOME}/lib:${JRE_HOME}/lib 
ENV PATH=${JAVA_HOME}/bin:$PATH 
ENV PATH="/usr/local/software/samtools/samtools-1.11/:${PATH}"

CMD ["bash"]
