# Express Tookit Containerization

This is a Docker wrapper around the SEQuoia Express Toolkit.

## Pull

```bash
docker pull bioraddbg/sequoia-express
```

## Build Container 
```
docker build -t bioreaddbg/sequoia-express:docker ./
```

## To run 
This will require docker commands and mapping out directories to run and will run with docker run instead of nextflow run 

```
docker run -rm -t -v /local_path/work_dir:/work \
	-v /local_path/ref_dir:/ref_data \
	-v /local_path/output_dir:/output \
	-v /local_path/fastq_dir:/data \
	bioreaddbg/sequoia-express:docker \
	--reads '/data/' \
	--outDir /output \
	--genomes_base /ref_data \
	-w /work \
	--profile indocker \
	--genome hg38	
```	
