# Sequoia_express_toolkit
Analysis toolkit for sequoia express RNAseq kits 

## Docker enviorment
To use the toolkit a virtual enviorement is required to run the software, prepared here as a docker container. To use please ensure docker is both installed and running. Either generate (via docker build) or retieve (docker pull) the container to continue. Future state as the container is made available publically nextflow will be able pull the container automatically.

## Analysis via Nextflow
Nextlfow is the primary software the runs and coodinates the pipeline (groovy / Java language base) so you will need Java 8 or higher with nextflow installed to run.

### Running the pipeline 
For the majority of users there are only some basic commands that will need to be done but for a full list of options please see the nextflow.config file, using `nextflow run main.nf --help` will only lis the basic options at the moment enough to get a basic run started. 

### Basic updates
This pipeline has been set up with bulk runs in mind, meaning that the predecessor Sequoia Complete took one file at time while this pipeline takes a whole directory of files at the same time. 
With this however your fastq files must have at a minium R1 / R2 in the file name to specify that they are reads

## Support
If you encounter and error / bug / issue please contact support@bio-rad.com to let us know or let us know via github so that we can address it. 
