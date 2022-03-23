![Bio-Rad Laboratories](src/vendor-logo.png?raw=true "Title")

# SEQuoia Express Toolkit
Analysis toolkit for SEQuoia Express Stranded RNA-Seq kits 

## Docker Environment
To use the toolkit a virtual environment is required to run the software, prepared here as a docker container. To use please ensure docker is both installed and running. Either generate (via docker build) or retrieve (docker pull) the container to continue. Please note that that nextflow can call docker directly and will be able to pull the container automatically.

## Analysis via Nextflow
Nextflow is the primary software the runs and coordinates the pipeline (groovy / Java language base) so you will need Java 8 or higher with nextflow installed to run.

### Installing Nextflow 
```
wget -qO- https://get.nextflow.io | bash
#or
curl -s https://get.nextflow.io | bash
```
If you are more comfortable with conda it can also be done there. 
```
conda install -c bioconda nextflow
```

### Cleanup the Work Directory
When running the toolkit, nextflow will produce intermediate files required to complete the processes. To do this please follow the instructions from nextflow. There are options to keep the logs or to do as part of a run after complete. 
One example would be:
```
nextflow clean -f ./work
```
for the full options:
```
nextflow clean -h 
```

### Downloading Reference Genomes
It is suggested that you copy the tar of the reference that you want locally. These commands will take a while to run.
For full list of options see: [Sequoia Genomes](https://www.dropbox.com/sh/kqy6kt9qewqsmbl/AABSjlIs87-cWMLdLPd8eDOja?dl=0) 

When downloading from Dropbox, it will add `?dl=0` to the end of each link. This needs to be removed. 
```
mkdir ./ref_data/genome-annotations
cd ./ref_data/genome-annotations
wget -O hg38.tar.gz https://www.dropbox.com/s/hm6kyp70dtbqovr/hg38.tar.gz?dl=0
tar xvzf hg38.tar.gz
```

### Running the pipeline 
For most users there are only some basic commands that will need to be done to run the pipeline. For a full list of options, please see the `nextflow.config` file. Using `nextflow run main.nf --help` will only list the basic options. 

#### Generate the Docker Image
This pipeline uses a docker container as a virtual environment to run the software. Outside of installing docker and nextflow, no other software is required. The recommendation for running this analysis is to pull the docker container from dockerhub. However, should the user choose to modify the docker container for a customized analysis, the Dockerfile is provided in this repository. 

(Recommended) the container will be pulled automatically by nextflow, however this is the command required if needed: 

```
docker pull -t bioraddbg/sequoia-express:latest
```

For the custom analysis described above, the following command will generate a docker container for the analysis:

```
docker build -t bioraddbg/sequoia-express [path to Dockerfile]
```

#### Running the pipeline for analysis with nextflow 

#### Typical Usage:

```
nextflow run Sequoia_express_toolkit/main.nf  --outDir ./output/ --reads '~/read/express/' --genome hg38 --genomes_base ./genomes/
```
#### Help
```
$ nextflow run main.nf --help

/-----------------------------------------------------------\ 
| __________.__                __________             .___  |
|  \_____   \__|____           \______   \____      __| _/  |
|   |  |  _/  |/  _ \   ______   |     _/\__  \   / __ |    |
|   |  |   \  (  <_> ) /_____/   |  |   \ / __ \_/ /_/ |    |
|   |____  /__|\____/            |__|_  /(____  /\____ |    |
|        \/                           \/      \/      \/    |
\___________________________________________________________/



Usage:

The typical command for running the pipeline is as follows:
nextflow run Sequoia_express_toolkit/main.nf  --outDir ./output/ --reads '~/read/express/' --genome hg38 --genomes_base ./genomes/

Args:

REQUIRED:
    genome               (string )           Genome to align to and annotate against                                                                 [hg38, mm10, rnor6]       
    genomes_base         (string )           Bio-Rad formatted refence genomes and annotations                                                                                 
    reads                (string )           The path to the fastq files must be wrapped in single quotes.                           

OPTIONAL:
    fivePrimeQualCutoff  (integer)           The read quality below which bases will be trimmed on the 5' end        		                     [0, 42]                   
    max_cpus             (integer)           The max number of cpus the pipeline may use. Defaults provided by -profile.                                                       
    max_memory           (integer)           The max memory in GB that the pipeline may use. Defaults provided by -profile.                                                    
    minBp                (intger ) 15        Reads with fewer base pairs will be rejected                                                            [0, 500]                  
    minGeneCutoff        (double )           Provide double value to cutt off how many reads are minimum                                             [0, 9E+7]                 
    minGeneType          (string )           Provide metric to be used                                                                               [none, reads, RPKM, TPM]  
    minMapqToCount       (integer) 1         The minimum MapQ socre for an aligned read to count toward a feature count                              [0, 255]                  
    noTrim               (boolean)           Indicates whether or not trimming skipped on the reads                                                                            
    outDir               (string ) ./results Indicate the output directory to write to                                                                                         
    reverseStrand        (boolean)           Indicate if your library is reverse stranded                            
                                                          
    seqType              (string )           Provide sequencing method used, if SE provided deduplication will not occur                             [SE, PE]                  
    skipUmi              (boolean)           Indicate no UMI processing is required                                                        
    spikeType            (string ) NONE      The type of spike in used, if any                                                                       [NONE, ercc]              
    threePrimeQualCutoff (integer)           The read quality below which bases will be trimmed on the 3' end                                        [0, 42]                   
    validateInputs       (boolean) true      Ensure input meets standards and is below 500 million reads
```


### Basic updates
This pipeline has been set up with multiple-sample bulk runs in mind, meaning that the predecessor Sequoia Complete took one file at time while this pipeline takes a whole directory of files at the same time. 
With this however your fastq files must have at a minimum R1 / R2 in the file name to specify that they are paired reads.

### Expected Outputs
This pipeline creates output like those used for Sequoia Complete, each individual sample will have a report in csv, html, and pdf formats. Additionally, each batch that is run will have its own high-level report that is created to have a side by side comparison of metrics as well.

## Support
If you encounter an error / bug / issue, please contact support@bio-rad.com or submit and issue to this repository so that we can address it.

### Tips and Tricks
If you find that you are getting an error where nextflow cannot find your files, check your path, and if needed use an absolute path, or check the formatting on your relative path. Also check your reads have R1 / R2 (in caps) and end with .fastq or fastq.gz

The pipeline runs with paired end as default (assumes you have both R1 and R2) if this is not the case you can run `--seqType=SE` to use just the R1 reads 
