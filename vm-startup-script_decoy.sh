#!/bin/bash
# This is the startup script to configure the VM with necessary tools and download the transcript index.

df -h

# Install required packages and update the package list if needed
sudo apt-get install -y libtbb2 libtbb-dev # needed for salmon

# Check if apt-get update is needed
if ! dpkg -s samtools &> /dev/null || ! dpkg -s bedtools &> /dev/null || ! dpkg -s python3 &> /dev/null || ! dpkg -s default-jre &> /dev/null; then
    sudo apt-get update
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "samtools not found, installing..."
    sudo apt-get install -y samtools
else
    echo "samtools is already installed."
fi

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "bedtools not found, installing..."
    sudo apt-get install -y bedtools
else
    echo "bedtools is already installed."
fi

df -h

# Check if python3 is installed
if ! command -v python3 &> /dev/null; then
    echo "python3 not found, installing..."
    sudo apt-get install -y python3 python3-pip
else
    echo "python3 is already installed."
fi

# Check if default-jre is installed
if ! command -v java &> /dev/null; then
    echo "default-jre not found, installing..."
    sudo apt-get install -y default-jre
else
    echo "default-jre is already installed."
fi

# Check if Salmon is installed
if ! command -v salmon &> /dev/null; then
    echo "Salmon not found, installing..."
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz || { echo "Salmon download failed"; exit 1; }
    tar -xvzf salmon-1.10.0_linux_x86_64.tar.gz
    sudo mv salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/
else
    echo "Salmon is already installed."
fi

df -h

# Set environment variables or paths (if not already set)
if ! [[ "$PATH" =~ "/usr/local/bin" ]]; then
    export PATH=$PATH:/usr/local/bin
    echo "Path updated."
else
    echo "Path already includes /usr/local/bin."
fi

# Define the path to google bucket
GOOGLE_BUCKET=gs://virus_expression_results

# Define the path to the transcript index
TRANSCRIPT_INDEX=~/transcripts_index/

# Check if the transcript index already exists
if [ ! -d "$TRANSCRIPT_INDEX" ]; then
    echo "Transcript index not found. Copying from Google Cloud Storage..."

    # Verify access to the Google Cloud bucket by listing its contents
    gsutil ls "$GOOGLE_BUCKET"/salmon_data/ || { echo "Access to bucket failed"; exit 1; }

    mkdir -p $TRANSCRIPT_INDEX

    # Copy the transcript index to the VM
    gsutil -m cp -r "$GOOGLE_BUCKET"/salmon_data/salmon_index/* "$TRANSCRIPT_INDEX" || { echo "Failed to copy transcript index"; exit 1; }

    echo "Transcript index successfully copied."
else
    echo "Transcript index already exists."
fi

# Create directories needed for python script
TMP_FQ_DATA=~/tmp/cell_line_fq_data/"$1"/
TMP_RESULTS=~/tmp/transcript_quant_results/"$1"/
mkdir -p $TMP_FQ_DATA
mkdir -p $TMP_RESULTS

# Define the path to the transcript index
GCLOUD_INIT=~/.config/gcloud/
GCLOUD_MY_ACC=~/.config/gcloud/cache/ptrollma@broadinstitute.org/

# Check if the correct account already exists
if [ ! -d "$GCLOUD_MY_ACC" ]; then
    echo "gcloud init files not found. Copying from Google Cloud Storage..."

    # Verify access to the Google Cloud bucket by listing its contents
    gsutil ls "$GOOGLE_BUCKET"/gsutil_config/ || { echo "Access to bucket failed"; exit 1; }

    mkdir -p $GCLOUD_INIT

    # Copy the transcript index to the VM
    gsutil -m cp -r "$GOOGLE_BUCKET"/gsutil_config/gcloud/* "$GCLOUD_INIT" || { echo "Failed to copy transcript index"; exit 1; }

    echo "gcloud credentials successfully copied."
else
    echo "gcloud credentials already exists."
fi

df -h

# Call the Python script with the parameters passed from Sparkle
# The parameters are: $1 = cell_line, $2 = bam_location
python3 run_viral_annotation.py "$1" "$2" "$GOOGLE_BUCKET" "$TRANSCRIPT_INDEX" "$TMP_FQ_DATA" "$TMP_RESULTS"

df -h

# Cleanup temporary FASTQ files
rm -rf $TMP_FQ_DATA

df -h

# Calculation completion
echo "Calculations finished" >> /var/log/calculation-script.log
