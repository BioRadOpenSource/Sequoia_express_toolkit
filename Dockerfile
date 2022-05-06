FROM condaforge/mambaforge:4.12.0-0

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

ARG CONDA_ENV="SequoiaExpress"

LABEL Bio-Rad Support <support@bio-rad.com>

RUN apt-get --allow-releaseinfo-change update && \
#	apt-get install -y procps=2:3.3.15-2 && \
	apt-get clean -y

COPY $CONDA_ENV.yaml /opt/biorad/env/
RUN conda env create -f /opt/biorad/env/$CONDA_ENV.yaml

RUN conda clean -afy

#RUN conda update conda -y
#RUN conda update --all
#RUN apt-get update && apt-get install -y \
#    curl \
#    unzip \
#    perl \
#    parallel \
#    pigz \
#    wget \
#    samtools \
#    zlib1g-dev \
#    libbz2-dev \
#    liblzma-dev

######### Fix time and date interaction
######################################

######### FastQC Setup ###############
#RUN apt-get install -y \
#    openjdk-8-jre-headless

#ENV FASTQC_URL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#ENV FASTQC_VERSION 0.11.7
#ENV FASTQC_RELEASE fastqc_v${FASTQC_VERSION}.zip
#ENV DEST_DIR /opt/
# Make destination directory
#RUN mkdir -p $DEST_DIR
#RUN echo "This is a test"
# Do this in one command to avoid caching the zip file and its removal in separate layers
#RUN curl -SLO ${FASTQC_URL}/${FASTQC_RELEASE} && unzip ${FASTQC_RELEASE} -d ${DEST_DIR} && rm ${FASTQC_RELEASE}
# Make the wrapper script executable
#RUN chmod a+x ${DEST_DIR}/FastQC/fastqc
# Include it in PATH
#ENV PATH ${DEST_DIR}/FastQC:$PATH
######### End FastQC Setup ###########

######### Cutadapt Setup #############
#RUN apt-get install -y \
#    python3-pip
#RUN pip3 install --upgrade pip setuptools

#RUN pip3 install --upgrade pip setuptools
#RUN pip3 install --user --upgrade 'cutadapt==2.7'
#RUN ln -s ~/.local/bin/cutadapt /usr/bin/
######### End Cutadapt Setup #########

######### UMI Tools Setup ############
#RUN pip3 install --user --upgrade umi_tools
#RUN ln -s ~/.local/bin/umi_tools /usr/bin/umi_tools
######### END UMI Tools Setup ########

######### STAR Setup #################
#ENV STAR_VERSION 2.7.0f
# Same deal as above with FASTQC; using precompiled executable
#RUN curl -SLO https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && tar -zxvf ${STAR_VERSION}.tar.gz --directory /opt/ && rm ${STAR_VERSION}.tar.gz
#ENV PATH /opt/STAR-${STAR_VERSION}/bin/Linux_x86_64:$PATH
######### End STAR Setup #############

######### PICARD Setup ###############
#ENV PICARD_VERSION 2.20.0
#RUN mkdir -p ${DEST_DIR}/picard
#RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -P ${DEST_DIR}/picard/
#ENV PATH ${DEST_DIR}/picard:$PATH
######### End PICARD Setup ##########

######### Subread Setup #############
#ENV SUBREAD_VERSION 1.6.4 
#RUN mkdir -p ${DEST_DIR}/subread
#RUN curl -SLO https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz \
#    && tar -zxvf subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz --directory /opt/ \
#    && rm subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz
#ENV PATH ${DEST_DIR}/subread-${SUBREAD_VERSION}-Linux-x86_64/bin/:$PATH
######### End Subread Setup ##########

######### Sambamba Setup #############
#ENV SAMBAMBA_VERSION 0.6.9
#RUN mkdir -p ${DEST_DIR}/sambamba
#RUN curl -SLO https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-static.gz \
#    && unpigz sambamba-${SAMBAMBA_VERSION}-linux-static.gz && mv sambamba-${SAMBAMBA_VERSION}-linux-static ${DEST_DIR}/sambamba/sambamba
#RUN chmod a+x ${DEST_DIR}/sambamba/sambamba
#ENV PATH ${DEST_DIR}/sambamba/:$PATH
######### End Sambamba Setup #########

##########Install Pandoc #############
#ENV PANDOC_VERSION 1.16.0.2
#RUN mkdir -p ${DEST_DIR}/pandoc
#RUN curl -SLO https://hackage.haskell.org/package/pandoc-${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}.tar.gz \
#	&& tar xvzf pandoc-${PANDOC_VERSION}.tar.gz && mv pandoc-${PANDOC_VERSION} ${DEST_DIR}/pandoc/pandoc
#RUN chmod a+x ${DEST_DIR}/pandoc/pandoc
#ENV PATH ${DEST_DIR}/pandoc/:$PATH
#RUN rm pandoc-${PANDOC_VERSION}.tar.gz
#####################################

######### Pysam Setup ################
#RUN pip3 install pysam
######### End Pysam Setup ############

######### RPKM/TPM Setup #############
#RUN pip3 install pandas
######### End RPKM/TPM Setup #########

######### pythong excel info ########
#RUN pip3 install openpyxl

####################################

######### R Setup ###############
#RUN apt-get update && apt-get install -y \
#    apt-transport-https \
#    libcurl4-openssl-dev \
#    pandoc \
#    libssl-dev \
#    libxml2-dev \
#    fonts-freefont-ttf \
#    libfontconfig1-dev
#ENV R_VER=3.6.3
#RUN curl -SLO https://cran.rstudio.com/src/base/R-3/R-${R_VER}.tar.gz
#RUN tar -xzvf R-${R_VER}.tar.gz
#RUN mkdir -p ${DEST_DIR}/R
#RUN mv R-${R_VER}/ ${DEST_DIR}/R/
#RUN cd ${DEST_DIR}/R/R-${R_VER}/
#RUN ./configure --with-x=no
#RUN make
#RUN make install
#RUN cd 
#RUN ln -s ${DEST_DIR}/R/R-${R_VER}/R /usr/local/bin/R
#RUN ln -s ${DEST_DIR}/R/R-${R_VER}/Rscript /usr/local/bin/Rscript
#RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list 
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 
#RUN apt update

#RUN apt-get install -y r-base -q
#RUN Rscript -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
#RUN Rscript -e 'install.packages(c("dplyr", "knitr", "rmarkdown", "kableExtra", "ggplot2", "plotly", "fastqcr", "data.table", "tibble", "rlist", "tinytex", "webshot", "DT"), repos = "http://cran.r-project.org")'
#RUN Rscript -e 'tinytex::install_tinytex()'
#RUN Rscript -e 'webshot::install_phantomjs()'
#RUN Rscript -e 'tinytex::tlmgr_update()'
#RUN Rscript -e 'tinytex::tlmgr_install(pkgs = c("xcolor", "colortbl", "multirow", "wrapfig", "float", "tabu", "varwidth", "threeparttable", "threeparttablex", "environ", "trimspaces", "ulem", "makecell", "titling","mathspec","fancyhdr"))'
######### End R Setup ###########

#Integrate RUST for DEAD and rumi#########
#RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
#RUN apt-get update
#RUN apt-get install -f -y clang
	
#RUN export LIBCLANG_PATH=/usr/lib32/

RUN rm root/.bashrc
RUN echo "source /etc/container.bashrc" >> /etc/bash.bashrc && \
	echo "set +u" > /etc/container.bashrc && \
	echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/container.bashrc && \
	echo "conda activate $CONDA_ENV" >> /etc/container.bashrc

# Activating environment when using non-login, non-interactive shell
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

# Adding Bio-Rad bin to path
ENV PATH /opt/biorad/bin/:$PATH

WORKDIR /opt/biorad

COPY . .


#install rumi and DEAD
#RUN ["/bin/bash", "-c", "source ~/.cargo/env"]
#RUN apt-get update -y
#RUN apt-get install git -y
#RUN cargo install src/rumi/
ENV PATH=$PATH:/opt/biorad/src
ENV TZ=US
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get install texlive-xetex texlive-fonts-recommended texlive-plain-generic texlive-fonts-extra -y

WORKDIR /opt/conda/envs/$CONDA_ENV/lib

#fix for issues / conflict between samtools and pandoc / R 
#samtools needs libcrypto.so.1.0.0 and wont accept other versions as of May22
RUN cp libcrypto.so.1.1 libcrypto.so.1.0.0

#pdf report wants freesans

WORKDIR /opt/biorad 
# Pull in some ARGS for defining container name
ARG IMAGE_NAME
ARG SOURCE_BRANCH
ARG SOURCE_COMMIT
RUN printf "Container Name: ${IMAGE_NAME:-local}\n" > imageInfo.txt
RUN printf "Source Branch: ${SOURCE_BRANCH:-local}\n" >> imageInfo.txt
RUN printf "Source Commit: ${SOURCE_COMMIT:-local}\n" >> imageInfo.txt

