![Bio-Rad Laboratories](src/vendor-logo.png?raw=true "Title")

# Sequoia Express Toolkit
Analysis toolkit for Sequoia Express RNAseq kits 

## Docker Enviorment
To use the toolkit a virtual enviorement is required to run the software, prepared here as a docker container. To use please ensure docker is both installed and running. Either generate (via docker build) or retieve (docker pull) the container to continue. Please note that that nextflow can call docker directly and will be able to pull the container automatically.

## Analysis via Nextflow
Nextflow is the primary software the runs and coodinates the pipeline (groovy / Java language base) so you will need Java 8 or higher with nextflow installed to run.

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
When running the toolkit, nextflow will produce intermediate files required to complete the processes. To do this please follow the instuctions from nextflow. There are options to keep the logs or to do as part of a run after complete. 
One example would be:
```
nextflow clean -f ./work

```
for the full options:
```
nextflow clean -h 
```


### Downloading Reference Genomes
The reference genomes are stored in S3 for convenience. As of this writing. The reference genome can be found here: s3://dbg-cloudpipeline-data-us-west-2-prod/ref\_data/sequoia\_analysis/latest/.

It is suggested that you copy the tar of the reference that you want locally. These commands will take a while to run.
```
mkdir /mnt/ref_data/genome-annotations
cd /mnt/ref_data/genome-annotations
aws s3 cp s3://dbg-cloudpipeline-data-us-west-2-prod/ref_data/sequoia_analysis/latest/hg38.tar.gz ./
aws s3 cp s3://dbg-cloudpipeline-data-us-west-2-prod/ref_data/sequoia_analysis/latest/mm10.tar.gz ./
aws s3 cp s3://dbg-cloudpipeline-data-us-west-2-prod/ref_data/sequoia_analysis/latest/rnor6.tar.gz ./
tar xvzf hg38.tar.gz
tar xvzf rnor6.tar.gz
tar xvzf mm10.tar.gz
md5sum -c ./*/*.chk
```

### Running the pipeline 
For the majority of users there are only some basic commands that will need to be done. For a full list of options, please see the nextflow.config file. Using `nextflow run main.nf --help` will only list the basic options. 

#### Generate the Docker Image
This pipeline uses a docker conainer as a virtual enviorment to run the software. Outside of installing docker and nextflow, no other software is required. The recommendation for running this analysis is to pull the docker container from dockerhub. However, should the user choose to modify the docker container for a customized analysis, the Dockerfile is provided in this repository. 

(Recommended) the container should be pulled automatically by nextflow, however this is the command required if needed: 

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
    reads                (string )           Tese must be wrapped in single quotes. If R{1,2} is specified, UMI dedu$lication processes will be run.                           

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
                                                          
    seqType              (string )           Provide sequecing method used                                                                           [SE, PE]                  
    skipUmi              (boolean)           Indicate that only R1 has been passed in and no UMI processing is required                                                        
    spikeType            (string ) NONE      The type of spike in used, if any                                                                       [NONE, ercc]              
    threePrimeQualCutoff (integer)           The read quality below which bases will be trimmed on the 3' end                                        [0, 42]                   
    validateInputs       (boolean) true      Ensure input meets standards and is below 500 million reads
```


### Basic updates
This pipeline has been set up with multiple-sample bulk runs in mind, meaning that the predecessor Sequoia Complete took one file at time while this pipeline takes a whole directory of files at the same time. 
With this however your fastq files must have at a minium R1 / R2 in the file name to specify that they are paired reads.

### Expected Outputs
This pipeline creates output similar to those used for Sequoia Complete, each individual sample will have a report in csv, html, and pdf formats. Additionally, each batch that is run will have its own high level report that is created to have a side by side comparision of metrics as well.

## Support
If you encounter an error / bug / issue please contact support@bio-rad.com or submit and issue to this repository so that we can address it.

### Tips and Tricks
If you find that you are getting an error where nextflow can not find you files check your path, and if needed use an absolute path, or check the formatting on your relative path. Also check your reads have R1 / R2 (in caps) and end with .fastq or fastq.gz

The pipeline runs with paired end as default (assumes you have both R1 and R2) if this is not the case you can run --seqType=SE to use just the R1 reads 
