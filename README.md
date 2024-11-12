# Viral Expression Quantification Pipeline
## Data preparation steps
Run the research_viral_annotations.ipynb script to prepare all the necessary files.

## If you just want to run the quantification of the viral expression for several cell lines then you can follow these steps. 

We need a csv list of the cell lines of interest. The csv needs a column with an identifier for the cell line and the location of bam file of the google bucket.
Then we need to set up a normal VM (we will name it here “my-vm”) in the Google Cloud Services and a Google Bucket (“virus_expression_results”) were the results shall then be stored, as well as some files for our quantification.
On the VM you need to install sparklespray (https://github.com/broadinstitute/sparklespray) as this will create several VMs for us and run our script.
The script will first run a start up process to install all needed software.
Then the script will make a bam to fastq file conversion (takes 40-60 minutes), is needed as salmon (the software for our quantifications) needs the data in fastq format, but the data is stored in a compressed bam format to save memory.
The script will then run the quantification with salmon and give out the results to the google bucket.

Synchronize needed data with the your VM 
rsync -avz ~/experiment_data/virus_genes/ my-vm:/home/ubuntu/experiments/virus_genes


For eventual updates to your script you can synchronize it like this
###
    rsync -avz ~/experiment_data/virus_genes/scripts/ my-vm:/home/ubuntu/experiments/virus_genes/scripts

Upload the indexing files from your VM for the quantification to your google bucket
###
    gsutil -m cp -r ~/experiments/virus_genes/salmon_index/ gs://virus_expression_results/salmon_data/salmon_index

Run the script with sparkles
###
    conda activate sparkles
    sparkles sub -n run_viral_annotation_decoy_all -u vm-startup-script_decoy.sh -u run_viral_annotation.py --params cell_line_google_bucket_index.csv --nodes 250 bash vm-startup-script_decoy.sh '{entity:sample_id}' '{hg38_rna_bam}'

Get the result to your machine and analyse with Python script
###
    gsutil -m cp -r -c gs://virus_expression_results/transcript_quant_results/* ~/experiment_data/virus_genes/transcripts_quant

Google bucket with all the results:
gs://virus_expression_results/transcript_quant_results

## If you want to understand the full process then read the details below.
Have your VMs and Google Bucket be in the same region as the storage buckets, as this reduces the costs.
The bam files are accessed directly from the google bucket and converted to fastq files with samtools and bedtools.
I you want to make the salmon index yourself you also need to install salmon to the VM (https://salmon.readthedocs.io/en/latest/index.html) and use a VM with enough CPU and RAM, n2-highmem-4 works (but anything less can make problems, mostly the task will get killed from the OS as it uses to much RAM), this then can still take around 1.5 hours.

###
    gsutil -u {your-project-name} cat gs://cclebams/rnasq_hg38/CDS-37go6g.Aligned.sortedByCoord.out.bam | samtools sort -n | bedtools bamtofastq -i - -fq reads_1.fq -fq2 reads_2.fq

We then also need an index for the salmon quantification. We can either do this without decoy, is easier to set up, but with decoy is supposed to be even more exact.

Download the human genome and viral genes with the Python script.

### Without the decoy
Make a reference.fasta file with Python script
Run this terminal script
###
    conda activate salmon
    salmon index -t reference_genes.fasta -i transcripts_index -k 31

### With the decoy
Make the decoy.txt file 
Make a gentrome fasta file and a decoy.txt with Python script
Run this terminal script
###
    conda activate salmon
    salmon index -t gentrome.fa -d decoys.txt -p 12 -i salmon_index --gencode -k 31

Now you can run the salmon quantification either on individual cell lines, like
###
    conda activate salmon
    salmon quant -i transcripts_index -l A -1 cell_line_data/ACH-000556/reads_1.fq -2 cell_line_data/ACH-000556/reads_2.fq --validateMappings -o transcripts_quant_no_decoy/ACH-000556

Or run lots of VMs at the same time with the sparkles command
###
    conda activate sparkles
    sparkles sub -n run_viral_annotation -u vm-startup-script.sh -u run_viral_annotation.py --params cell_line_google_bucket_index.csv --nodes 2 bash vm-startup-script.sh '{entity:sample_id}' '{hg38_rna_bam}'


Sparkle needs a config file, set it up like this
Use your project name and docker image, use a VM with high memory and enough disk space
###
    nano ~/.sparkles

This config file works for most parts, but some files are so large that they need >300GB disk space. Also the calculation needs lots of RAM and CPU, so a bigger one works better

```[config]
default_url_prefix=gs://virus_expression_results
project=philipp-trollmann
default_image=us-central1-docker.pkg.dev/philipp-trollmann/hermit-dev/hermit-ubuntu:v09112024
machine_type=n2-highmem-2
zones=us-central1-a
region=us-central1
boot_volume_in_gb=100```

config file for sparkles
(a smaller VM and less memory space work for most files, but some are too big and need more memory)

```[config]
default_url_prefix=gs://virus_expression_results
project=philipp-trollmann
default_image=us-central1-docker.pkg.dev/philipp-trollmann/hermit-dev/hermit-ubuntu:v09112024
machine_type=n2-highmem-8
zones=us-central1-a
region=us-central1
boot_volume_in_gb=400
mount_1_size_in_gb=50```


The output are quant.sf files for every cell line that contains the TPM values for every gene


## Data analysis
Run the analysis Python script provided to get correlations and plots.

