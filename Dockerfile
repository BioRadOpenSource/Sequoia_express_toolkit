FROM condaforge/mambaforge:24.3.0-0

# Building image using bash shell
SHELL ["/bin/bash", "-c"]

ARG CONDA_ENV="SequoiaExpress"

LABEL Bio-Rad Support <support@bio-rad.com>

RUN apt-get --allow-releaseinfo-change update && \
	apt-get clean -y

COPY $CONDA_ENV.yaml /opt/biorad/env/
RUN conda env create -f /opt/biorad/env/$CONDA_ENV.yaml

RUN conda clean -afy
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

#set install for the pandoc and PDF needs as well as time data 
ENV PATH=$PATH:/opt/biorad/src
ENV TZ=US
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# add fonts needs for the pdf report
RUN apt-get update
RUN apt-get install texlive-xetex texlive-fonts-recommended texlive-plain-generic texlive-fonts-extra -y

WORKDIR /opt/conda/envs/$CONDA_ENV/lib

#fix for issues / conflict between samtools and pandoc / R 
#samtools needs libcrypto.so.1.0.0 and wont accept other versions as of May22
RUN cp libcrypto.so.1.1 libcrypto.so.1.0.0

WORKDIR /opt/biorad 
# Pull in some ARGS for defining container name
ARG IMAGE_NAME
ARG SOURCE_BRANCH
ARG SOURCE_COMMIT
RUN printf "Container Name: ${IMAGE_NAME:-local}\n" > imageInfo.txt
RUN printf "Source Branch: ${SOURCE_BRANCH:-local}\n" >> imageInfo.txt
RUN printf "Source Commit: ${SOURCE_COMMIT:-local}\n" >> imageInfo.txt

