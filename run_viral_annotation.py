import argparse
import subprocess
import sys

def create_gcs_placeholder(output_folder):
    # Create a placeholder to simulate folder structure in GCS
    placeholder_file = f"{output_folder}placeholder"
    subprocess.run(f"gsutil cp /dev/null {placeholder_file}", shell=True, check=True)


def process_cell_line(cell_line, bam_location, transcript_index, output_bucket, local_output_folder, local_results_folder):
    # Convert BAM to FASTQ
    read1_fq = f"{local_output_folder}reads_1.fq"
    read2_fq = f"{local_output_folder}reads_2.fq"

    try:
        print(f'Bam conversion for {cell_line} starting')
        bam_conversion_command = f"gsutil -u philipp-trollmann cat {bam_location} | samtools sort -n | bedtools bamtofastq -i - -fq {read1_fq} -fq2 {read2_fq}"
        subprocess.run(bam_conversion_command, shell=True, check=True)
        print(f'Bam conversion for {cell_line} finished')

        # Run Salmon quantification
        print(f'Salmon quantification for {cell_line} starting')
        salmon_command = f"salmon quant -i {transcript_index} -l A -1 {read1_fq} -2 {read2_fq} --validateMappings -o {local_results_folder}"
        subprocess.run(salmon_command, shell=True, check=True)
        print(f'Salmon quantification for {cell_line} finished')

        # # Compress FASTQ files
        # print(f'Compressing FASTQ files for {cell_line}')
        # subprocess.run(f"gzip {read1_fq} {read2_fq}", shell=True, check=True)
        # print(f'FASTQ files for {cell_line} compressed')
        
        # # Upload compressed FASTQ files to Google Cloud Storage
        # output_folder = f"{output_bucket}/cell_line_fq_data/{cell_line}/"
        # create_gcs_placeholder(output_folder)  # Create placeholder for FASTQ output folder
        # subprocess.run(f"gsutil -m cp -r {local_output_folder}*.gz {output_folder}", shell=True, check=True)
        # print(f'FASTQ files for {cell_line} uploaded to GCS')

        # Upload Salmon quantification results to Google Cloud Storage
        results_folder = f"{output_bucket}/transcript_quant_results/{cell_line}/"
        create_gcs_placeholder(results_folder)  # Create placeholder for results folder
        subprocess.run(f"gsutil -m cp -r {local_results_folder}* {results_folder}", shell=True, check=True)
        print(f'Salmon quantification results for {cell_line} uploaded to GCS')

    except subprocess.CalledProcessError as e:
        print(f"Error occurred during processing of {cell_line}: {e}")
        sys.exit(1)


def main():
    # Argument parser to accept cell line name and BAM location
    parser = argparse.ArgumentParser(description="Process cell line and BAM location.")
    parser.add_argument("cell_line", type=str, help="The cell line ID")
    parser.add_argument("bam_location", type=str, help="The location of the BAM file in Google Cloud Storage")
    parser.add_argument("output_bucket", type=str, help="The location of the google bucket to use")
    parser.add_argument("transcript_index", type=str, help="The location of the local transcript index for the salmnon quantification")
    parser.add_argument("local_output_folder", type=str, help="The location of the local tmp folder for the fq files")
    parser.add_argument("local_results_folder", type=str, help="The location of the local tmp folder for the restuls files")

    args = parser.parse_args()

    print(args)

    # Process the cell line and BAM location passed as arguments
    process_cell_line(args.cell_line, args.bam_location, args.transcript_index, args.output_bucket, args.local_output_folder, args.local_results_folder)   


if __name__ == "__main__":
    main()
