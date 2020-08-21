FROM    centos:8

RUN     yum install -y wget which git java-1.8.0-openjdk perl python3-pip

RUN     mkdir /apps && \
        cd /apps && \
        git clone https://github.com/marbl/merqury.git && \
        wget -q -O meryl.tar.gz 'https://github.com/marbl/meryl/releases/download/v1.0/meryl-1.0.Linux-amd64.tar.xz' && \
        tar xf meryl.tar.gz && \
        rm meryl.tar.gz

RUN     wget -q -O fastqc_v0.11.9.zip 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip' && \
	    unzip fastqc_v0.11.9.zip && \
        rm fastqc_v0.11.9.zip

RUN     wget -q -O Trimmomatic-0.39.zip 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip' && \
        unzip Trimmomatic-0.39.zip && \
        rm Trimmomatic-0.39.zip

RUN     yum install -y java-1.8.0-openjdk libgomp perl python3-pip tbb && \
        pip3 install biopython && \
        pip3 install multiqc

RUN	    useradd -ms /bin/bash localuser

RUN     chmod +x /apps/meryl-1.0/Linux-amd64/bin/* && \
        chmod +x /apps/merqury/*.sh && \
        chmod +x /apps/merqury/eval/*.sh && \
        chmod +x /apps/trimmomatic && \
	    chmod +x /apps/FastQC/fastqc && \
	    chown -R localuser /apps

ENV     PATH=${PATH}:/apps/FastQC:/apps/Trimmomatic-0.39:/apps/meryl-1.0/Linux-amd64/bin:/apps/merqury/eval

USER 	localuser
WORKDIR /home/localuser

LABEL   author="Sergej Nowoshilow"
LABEL   maintainer="Sergej Nowoshilow"
LABEL   email="Sergej.Nowoshilow@imp.ac.at"
LABEL   description="This container image contains tools required to run nf-merqury pipeline"
LABEL   version="1.0"