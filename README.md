# Sequoia_express_toolkit
Analysis toolkit for sequoia express RNAseq kits 

## Docker enviorment
To use the toolkit a virtual enviorement is required to run the software, prepared here as a docker container. To use please ensure docker is both installed and running. Either generate (via docker build) or retieve (docker pull) the container to continue. Future state as the container is made available publically nextflow will be able pull the container automatically.

## Analysis via Nextflow
Nextlfow is the primary software the runs and coodinates the pipeline (groovy / Java language base) so you will need Java 8 or higher with nextflow installed to run.

### Running the pipeline 
For the majority of users there are only some basic commands that will need to be done but for a full list of options please see the nextflow.config file, using `nextflow run main.nf --help` will only lis the basic options at the moment enough to get a basic run started. 

#### Generate the docker image needed for the virtual enviorment
This pipeline uses a docker conainer as a virtual enviorment to run the software as OS agnostic as possible. So outside of intsalling docker and nextflow no other software is required. To use this docker container one can simply build it from the included docker file. 

```
Docker build -t bioraddbg/sequoia-express [path to Dockerfile]

```
Alernatively this docker file will also be created and pushed after finalization to dockerhub, where it can be pulled directly with no extra fuss.

```
Docker pull -t bioraddbg/sequoia-express:latest
```

#### Running the pipeline for analysis with nextflow 
For basic options
```
nextflow run repos/Sequoia_express_toolkit/main.nf --help
```

For a full run something like this will be your system call
```
nextflow run repos/Sequoia_express_toolkit/main.nf  --outDir ./output/express --reads ./data/ --genome hg38 --genomes_base ./genome/ --max_cpus 16 --max_memory 60 -with-docker bioraddbg/sequoia-express:latest -resume --seqType="PE"

```

### Basic updates
This pipeline has been set up with bulk runs in mind, meaning that the predecessor Sequoia Complete took one file at time while this pipeline takes a whole directory of files at the same time. 
With this however your fastq files must have at a minium R1 / R2 in the file name to specify that they are reads

### Expected Outputs
This pipeline creates output similar to those used for Sequoia Complete, each individual sample will have a report in csv, html and pdf format for ease of viewing. Additionally each batch that is run will have its own high level report that is created to have a side by side comparision of metrics as well.

## Support
If you encounter and error / bug / issue please contact support@bio-rad.com to let us know or let us know via github so that we can address it. 
