FROM bioraddbg/sequoia-express:latest
LABEL Bio-Rad Support <support@bio-rad.com>
######## Nextflow Setup ########
RUN apt-get install -y coreutils
#ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow
RUN curl -fsSL https://get.nextflow.io | bash && cp nextflow /usr/local/bin/
RUN nextflow info

COPY ./src/entrypoint.sh /opt/biorad/src/entrypoint.sh

ENTRYPOINT ["/bin/bash", "/opt/biorad/src/entrypoint.sh"]
