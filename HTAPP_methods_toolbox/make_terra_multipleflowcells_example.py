#!/usr/bin/env python
# coding: utf-8

# # Make mkfastq, cellranger_workflow, cumulus (pegasus), and cellbender remove-background files for running in Terra
#
# Last updated: 2021/10/15
#
# Author: Orr Ashenberg & Caroline Porter
#
# This script creates the `mkfastq` script for running Cell Ranger `mkfastq` locally on compute cluster to generate `fastq` files from `bcl` files. It then uploads the `fastq` files to your Google Cloud Storage Bucket. Next, the script creates the files for running `cellranger_workflow` and `cumulus` pipelines on Terra and uploads them to your workspace. The user must specify all the input settings in the first code cell of this notebook, including local  directories and Google Cloud Buckets, and the user must create the sample tracking file that records information on all their samples. No other code cells need to be modified, unless you want to change settings for running any of the pipelines (for example, running mkfastq on Terra using cellranger_workflow). When running this script as a jupyter notebook, terminal commands to copy files to and from the Google Cloud Storage buckets are printed as output from the code cells. When running this script as a `python` script from your terminal, these commands are written to the standard out of your terminal. This code was originally written to process Human Tumor Atlas Pilot Project (HTAPP) samples, but has been modified to process other samples with different organizational conventions.
#
# This code is written to process all the samples listed in a sample sheet containing information on all your samples. The code can process multiple biological samples coming from a sample tracking file. Individual count matrices are made for each 10x channel by `cellranger_workflow`.
# **Sample tracking file**
#
# The sample tracking file, in csv format, is a useful way to track the important information for each sample, and is needed to run this script. Each sample requires the following text fields.
# - date: The date your samples are processed in yyyy_mm_dd format.
# - run_pipeline: Boolean (True or False) that determines what samples are processed. Set this to TRUE for all samples you want to processs. All other samples must be set to FALSE. How this works in operation is that as you add your new samples, set them to run_pipeline = TRUE and set the previously run samples to run_pipeline = FALSE.
# - Channel Name: This is the sample name that is used by the experimental team to name a 10x channel.
# - sampleid: This is the sample id that will be used to name all of the output files.
# - condition: Any biological or technical condition used in collecting this sample. This could be a buffer used or a flow cytometry sorting gate, or just single cell (sc) or single nucleus (sn). If there is no special condition, label this as none.
# - tissue: The tissue of origin. This column is optional.
# - replicate: This is used to designate which 10x channel (up to 8 different channels on a Chromium chip) this sample was run on, and is useful when multiple 10x channels are run for the same biological sample. The channel number must be an integer but it does not matter what integers you choose, as long as different channels use different integers. We suggest 1, 2, 3... If there is only a single channel, label this as channel1.
# - Lane: The lane that the sample was sequenced on within the flowcell. It can be a single lane (ex: 5), several lanes (ex: 5-6), or * (all lanes on the flow cell).
# - Index: The 10x index for the sample.
# - project: The name of the project you'd like to see attached to your directories.
# - reference: The genome reference to use when Cell Ranger `count` is creating the counts matrices. Please choose from one of references listed in Cumulus read the docs. https://cumulus.readthedocs.io/en/stable/
# - chemistry: The sequencing chemistry used.
# - flowcell: The flowcell id from your sequencing run.
# - seq_dir: The directory of your sequencing results.
# - min_umis: the min number of UMIs you'd like to use for filtering when you run cumulus pegasus.
# - min_genes: the min number of genes you'd like to use for filtering when you run cumulus pegasus.
# - percent_mito: the max percentage of expression coming from mito genes that you'd like to set for filtering when you run cumulus pegasus.
# - cellbender_expected_cells: the expected number of cells to input to cellbender remove-background. This value cannot be selected until after you see the cellranger count output because it's based on the barcode rank plot. See the cellbender documentation for information on how to choose this value: https://cellbender.readthedocs.io/en/latest/getting_started/remove_background/index.html
# - cellbender_total_droplets_included: the expected number of total droplets to input to cellbender remove-background. This value cannot be selected until after you see the cellranger count output because it's based on the barcode rank plot. See the cellbender documentation for information on how to choose this value: https://cellbender.readthedocs.io/en/latest/getting_started/remove_background/index.html
# - calc_signature_scores: a path to a .gmt file you can use for testing gene signatures using cumulus pegasus.
#
# The input cell also requires the user to specify the location of a Cell Ranger `mkfastq` sample sheet template, that the script will modify with the appropriate settings. The template is provided in github with this python script.
#
# **Setting up Google Cloud authorization, Terra workspaces, and cellranger_workflow and cumulus pipelines**
#
# Follow the [instructions here](https://cumulus.readthedocs.io/en/latest/) to set up Google Cloud authorization. This only needs to be done a single time. You need to set up a Terra Workspace and have both writer permission and computes permission to write data and to run the computational pipelines.
#
# The directions for the KCO `cellranger_workflow` pipeline are located in https://cumulus.readthedocs.io/en/latest/. When you run this pipeline, at that point Cell Ranger `mkfastq` has already been run to generate `fastq` files, although it is possible to use this pipeline to also run `mkfastq`.
#
# After the count matrices are generated, we prepare the files to run the `cumulus` pipeline. These contain sample metadata and paths to the count matrices. The directions for the `cumulus` pipeline are located in https://cumulus.readthedocs.io/en/latest/)
#
# After running cumulus on the cellranger counts matrices, we prepare files for running the 'cellbender remove-background' pipeline to remove ambient RNA and empty droplets.
#
# After the new count matrices are produced by cellbender, we prepare files for again running the `cumulus` pipeline.
#
# Importing the pipelines only needs to be done once in the processed data workspace. To import the pipelines, do the following:
# - Go to Method Configurations.
# - Import Configuration.
# - Import from Method Repository.
# - Search for `cellranger_workflow`, `cumulus`, or `cellbender remove-background` and Select Configuration.
# - Use Blank Configuration.
# - Select Name to be `cellranger_workflow`, `cumulus`, or 'remove-background' and Root Entity Type to be participant.
# - Import Method.
#
# This python script will automatically generate commands to run all workflows from the command line using alto cumulus.
# To instead run the pipelines in the Terra workspace directly (without using alto cumulus), do the following:
# - Select Method Configurations, select `cellranger_workflow`, `cumulus`, or `remove-background`.
# - Check the box for "Run workflow with inputs defined by file paths."
# - Unceck "use call caching"
# - Click Populate with a .json file, upload the `json` file from your local computer, and verify the text field inputs are filled as expected. Ensure all previous text field inputs (from any earlier runs) have been cleared after uploading the `json` file.
# - Select Save and then Launch Analysis.
#
# Outputs of this script and steps for running the pipeline in Terra
# All outputs are organized by flow cell, but all commands except the download bash script only need to be run once (you can select the command from the output for any of the flow cells).
#
# The first few lines that are printed to the terminal are visual checks that the sample and its information, such as lane and index, are being parsed correctly by the script.
#
# To run mkfastq on your computer cluster, you'll submit a bash script `run_mkfastq_uploadfastq.sh` that runs Cell Ranger `mkfastq`, and when that is done, it uploads the `fastq` files to the Google Cloud Storage Bucket.
#
# If you prefer to run mkfastq on Terra, you'll need to upload your bcl files to Terra and adjust this script so that cellranger_workflow runs mkfastq and stores the fastq files in the right location for the cellranger count files in this script to find them.
#
# You can launch cellranger count in Terra from the command line using the script: run_alto_cellranger_workflow.sh
#
# You can monitor the cellranger workflow jobs in Terra under "Job History". When this workflow finishes, run cumulus pegasus from the command line using the bash script run_alto_cumulus.sh.
#
# Alternatively, you could use the .json files and sample sheets that were uploaded automatically to Terra to run these workflows from the user interface.
#
# Download the count matrices and cumulus results using the bash download script `run_download_FLOWCELL-ID.sh`. You'll have to run this script once per flowcell.
#
# To run cellbender remove-background, you'll first need to select the total number of cells and total number of droplet parameters for each sample based on the barcode rank plot in each cellranger count websummary. Then, you'll need to rerun this python script to generate new cellbender remove-background input files that will automatically upload to Google Cloud.
#
# After adding the cellbender remove-background parameters to your sample tracking sheet and rerunning this python script, you can launch cellbender remove-background from the command line using the run_alto_cellbender.sh script. Like the other workflows, this can alternatively be run from the user interface using the files automatically uploaded to the goolge bucket.
#
# When cellbender remove-background finishes, you can rerun cumulus pegasus on the cellbender counts matrices using the bash script run_alto_cellbender_cumulus.sh.
#
# # Download the cellbender and clustering results using the bash download script `cellbenderV2_cumulus/run_download_FLOWCELL-ID.sh`. You'll have to run this script once per flowcell.

import os
import pandas as pd
import subprocess
import re

# Local Broad UGER paths. fastqdir, countsdir, resultsdir should already exist, otherwise this program creates them.
basedir = "/path/to/project/directory/on/compute/cluster"
fastqdir = basedir + "/fastq"  # directory in which folder with Cell Ranger mkfastq sample sheet file will be created
if not os.path.exists(fastqdir):
    os.makedirs(fastqdir)
countsdir = basedir + "/counts"  # directory in which folders with count matrices will be created
if not os.path.exists(countsdir):
    os.makedirs(countsdir)
resultsdir = basedir + "/cumulus"  # directory in which folders with clustering results will be created
if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)
cellbenderdir = basedir + "/cellbender"  # directory in which folders with cellbender results will be created
if not os.path.exists(cellbenderdir):
    os.makedirs(cellbenderdir)
cellbender_resultsdir = basedir + "/cellbender_cumulus"  # directory in which folders with clustering results will be created
if not os.path.exists(cellbender_resultsdir):
    os.makedirs(cellbender_resultsdir)

# The sample tracking file contains all sample specific information. Only samples with run_pipeline=True are processed.
# sampletrackingfile is set by the user.
sampletrackingfile = basedir + "/fastq/sampletracking.csv"
sampletracking_alldata = pd.read_csv(sampletrackingfile)
sampletracking = sampletracking_alldata[sampletracking_alldata.run_pipeline == True]
project = sampletracking['project'].tolist()[0]  # tumor type in this workspace
seqdirs = set(sampletracking['seq_dir'])  # get the sequencing directories to prepare data from.

# Google Cloud Storage bucket paths set by the user. The paths can be different across fastq, counts, and results folders.
fastqbucket = "gs://fc-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/fastq_" + project  # raw data
countsbucket = "gs://fc-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/counts_" + project  # processed data
resultsbucket = "gs://fc-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/cumulus_" + project  # processed data
cellbenderbucket = "gs://fc-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/cellbender_" + project  # processed data
cellbender_resultsbucket = "gs://fc-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/cellbender_cumulus_" + project  # processed data

# Additional user settings.
email = "name@broadinstitute.org"  # email address for your Google Cloud Platform account

# User may specify Google Cloud workspace for running alto (optional, not needed to run this script).
# alto allows the user to run the pipelines on Terra (cellranger_workflow, cumulus) from the command line.
alto_workspace = "'kco/htapp'"  # single quotes are required to parse this as a single command-line argument
alto_counts_folder = re.sub(r'^gs://.*/', "", countsbucket)
alto_results_folder = re.sub(r'^gs://.*/', "", resultsbucket)
alto_cellbender_folder = re.sub(r'^gs://.*/', "", cellbenderbucket)
alto_cellbender_results_folder = re.sub(r'^gs://.*/', "", cellbender_resultsbucket)
run_alto_file = "%s/run_alto_cellranger_workflow.sh" % (countsdir)
open(run_alto_file, "w").close()
run_alto_file = "%s/run_alto_cumulus.sh" % (resultsdir)
open(run_alto_file, "w").close()
run_alto_file = "%s/run_alto_cellbender.sh" % (cellbenderdir)
open(run_alto_file, "w").close()
run_alto_file = "%s/run_alto_cellbender_cumulus.sh" % (cellbender_resultsdir)
open(run_alto_file, "w").close()
# User must specify location of template file for cellranger mkfastq for Cell Ranger 6.
# Cell Ranger 6 is used. User may specify version used by cellranger_workflow WDL.
cellranger_version = "6.0.1"
count_matrix_name = "raw_feature_bc_matrix.h5"
cellbender_matrix_name = "out_FPR_0.01_filtered.h5"
mkfastq_template_file = "/path/to/mkfast/template/cellranger_mfastq_local_6.0.2.sh"

# Loop over the sequencing directory path for each flow cell, and prepare files to process the data.
for seqdir in seqdirs:
    # Get sample data for this sequencing directory.
    sampletracking = sampletracking_alldata[(sampletracking_alldata.run_pipeline == True) &
                                            (sampletracking_alldata.seq_dir == seqdir)]

    flowcellid = sampletracking['flowcell'].tolist()[0]  # found in the sequencing directory, and used for fastqs folder in cellranger count
    date = sampletracking['date'].tolist()[0]  # placed in directory name for a sample in alignreads yyyy_mm_dd

    # ## Parse sample tracking sheet file and make Cell Ranger mkfastq samplesheet
    # Sample naming is determined by the sample names used in the cellranger mkfastq samplesheet. The sample naming uses the following conventions. Sampleid is the sampleid from the tumor tracking worksheet and sample is the name of a 10x channel processed from that sampleid. As an example, the sampleid could be NB335 and the samples could be NB335_TST_channel1 and NB335_TST_channel1.

    # Assign Sample name & collect info from sample tracking sheet
    sampletracking['Sample'] = sampletracking['sampleid']
    sampletracking = sampletracking[['date', 'run_pipeline', 'Channel Name', 'Sample', 'sampleid', 'condition', 'replicate', 'tissue', 'Lane', 'Index', 'project', 'reference', 'introns', 'chemistry', 'flowcell', 'seq_dir', 'min_umis', 'min_genes', 'percent_mito', 'cellbender_expected_cells', 'cellbender_total_droplets_included']]
    print(sampletracking)

    # Write Cell Ranger mkfastq sample sheet file.
    fastqsampledir = fastqdir + "/%s_%s" % (date, flowcellid)
    if not os.path.isdir(fastqsampledir):
        os.mkdir(fastqsampledir)
    samplesheetfile = '%s/samplesheet_mkfastq.csv' % fastqsampledir
    samplesheet = sampletracking[['Lane', 'Sample', 'Index']]
    samplesheet.to_csv(samplesheetfile, index=False)
    print(samplesheet)

    # Get sampleids from sample tracking sheet.
    sampleids = sampletracking['sampleid'].tolist()

    # Initialize sampledict where key is sampleid and value is a list of sample names associated with that sampleid.
    # e.g. {'NB335': ['NB335_TST_channel1', 'NB335_TST_channel2']}
    # Initialize mkfastqdict where key is sample name and value is flowcell lane, 10x index, reference, sequencing chemistry.
    # e.g. {'NB200_TST_channel1': ['7', 'SI-GA-C5', 'GRCh38_premrna', 'threeprime'], 'NB200_TST_channel2': ['7', 'SI-GA-D5', 'GRCh38_premrna', 'threeprime']}
    # Initialize sampledict where key is sampleid and value is cumulus settings for this sampleid.
    # e.g. {'NB335': [400, 200, 0.2]}
    sampledict = dict([(sample, []) for sample in sampleids])
    mkfastqdict = dict()
    cumulusdict = dict()
    cellbenderdict = dict()
    cellrangerdict = dict()
    for index, row in sampletracking.iterrows():
        sampledict[row['sampleid']].append(row['Sample'])
        mkfastqdict[row['Sample']] = [row['Lane'], row['Index'], row['reference'], row['chemistry']]
        cumulusdict[row['sampleid']] = [row['min_umis'], row['min_genes'], row['percent_mito']]
        cellbenderdict[row['sampleid']] = [row['cellbender_expected_cells'], row['cellbender_total_droplets_included']]
        cellrangerdict[row['sampleid']] = [row['introns']]


    print(sampledict)
    print(mkfastqdict)
    print(cumulusdict)
    print(cellbenderdict)

    # File to store the commands to run Cell Ranger mkfastq and to upload the fastq files for this flowcell.
    run_mkfastq_uploadfastq_file = "%s/run_mkfastq_uploadfastq.sh" % fastqsampledir

    # ## Cellranger mkfastq
    # Generate Cell Ranger `mkfastq` file.

    # Read cellranger mkfastq template file.
    with open(mkfastq_template_file, "r") as f:
        template = f.readlines()
    template[18] = "cellranger mkfastq --run=%s --samplesheet=%s --output-dir=%s" % (seqdir, samplesheetfile, fastqsampledir)

    # Write new mkfastq file with updated options.
    mkfastq_file = "%s/%s" % (fastqsampledir, os.path.basename(mkfastq_template_file))
    with open(mkfastq_file, "w") as f:
        f.writelines(template)

    # Write to bash file that will run Cell Ranger mkfastq and upload the fastq files.
    with open(run_mkfastq_uploadfastq_file, "w") as f:
        f.write("bash %s\n" % mkfastq_file)

    # ## Uploading fastq files to Google Cloud Storage bucket in Terra workspace
    # Write bash script to upload fastq files to Google Cloud Storage bucket in Terra workspace.

    uploadfastq_file = "%s/uploadfastq.sh" % fastqsampledir
    with open(uploadfastq_file, "w") as f:
        f.write("source /broad/software/scripts/useuse\n")
        f.write("use Google-Cloud-SDK\n")
        f.write("gcloud auth login %s\n" % email)
        f.write("cd %s/%s\n" % (fastqsampledir, flowcellid))
        for sample in mkfastqdict:
            f.write("gsutil -m cp -r %s*fastq*gz %s/%s/\n" % (sample, fastqbucket, sample))

    # Write to bash file that will run Cell Ranger mkfastq and upload the fastq files.
    with open(run_mkfastq_uploadfastq_file, "a") as f:
        f.write("bash %s\n" % uploadfastq_file)

    # Terminal commands to run cellranger mkfastq (40 minutes - 2 hours) and to upload fastq files.
    print("\n## Copy and run terminal commands below to run cellranger mkfastq (40 minutes - 2 hours), and to run bash script that uploads fastq files to Google Cloud Storage Bucket. ##")
    print("use UGER")
    print("ish -l h_vmem=3g -l os=RedHat7 -pe smp 12 -R y -binding linear:12")
    print("bash %s" % run_mkfastq_uploadfastq_file)

    # ## Uploading `cellranger_workflow` files to Google Cloud Storage bucket in Terra

    # Make samplesheet for cellranger_workflow.
    for sampleid in sampledict.keys():
        if not os.path.isdir("%s/%s" % (countsdir, sampleid)):
            os.mkdir("%s/%s" % (countsdir, sampleid))
        samplesheet_cellranger_file = "%s/%s/samplesheet_cellranger.csv" % (countsdir, sampleid)

        with open(samplesheet_cellranger_file, "w") as f:
            f.write("Sample,Reference,Flowcell,Chemistry\n")
            for sample in sampledict[sampleid]:
                lane = mkfastqdict[sample][0]
                index = mkfastqdict[sample][1]
                reference = mkfastqdict[sample][2]
                chemistry = mkfastqdict[sample][3]
                f.write("%s,%s,%s,%s\n" % (sample, reference, fastqbucket, chemistry))

    # Make input_cellranger file for cellranger_workflow.
    for sampleid in sampledict.keys():
        input_cellranger_file = "%s/%s/input_cellranger.json" % (countsdir, sampleid)

        with open(input_cellranger_file, "w") as f:
            f.write("{\n")
            f.write("\t\"cellranger_workflow.input_csv_file\" : \"%s/%s/samplesheet_cellranger.csv\",\n" % (countsbucket, sampleid))
            f.write("\t\"cellranger_workflow.output_directory\" : \"%s\",\n" % (countsbucket))
            f.write("\t\"cellranger_workflow.cellranger_version\" : \"%s\",\n" % (cellranger_version))
            f.write("\t\"cellranger_workflow.run_mkfastq\" : false,\n")
            f.write("\t\"cellranger_workflow.include_introns\" : %s\n" % str(cellrangerdict[sampleid][0]).lower())
            f.write("}\n")

    # Running bash script below to upload cellranger samplesheet and input file to Google Cloud Storage Bucket.
    print("\n## Bash script below is AUTOMATICALLY RUNNING to upload cellranger samplesheet and input file to Google Cloud Storage Bucket. ##")
    uploadcellranger_file = "%s/uploadcellranger.sh" % fastqsampledir
    with open(uploadcellranger_file, "w") as f:
        f.write("source /broad/software/scripts/useuse\n")
        f.write("use Google-Cloud-SDK\n")
        f.write("gcloud auth login %s\n" % email)
        for sampleid in sampledict.keys():
            samplesheet_cellranger_file = "%s/%s/samplesheet_cellranger.csv" % (countsdir, sampleid)
            input_cellranger_file = "%s/%s/input_cellranger.json" % (countsdir, sampleid)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cellranger_file, countsbucket, sampleid))
            f.write("gsutil cp %s %s/%s/\n" % (input_cellranger_file, countsbucket, sampleid))
    command = "bash %s" % uploadcellranger_file
    subprocess.call(command, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    print(command)

    print("\n## Run cellranger_workflow in Terra (or run alto bash script below) and wait for jobs to complete. ##")


    # ## Running `cellranger_workflow` WDL in Terra from the command line.

    # Write bash script below to run alto to kick off cellranger_workflow jobs from command line to Terra.
    run_alto_file = "%s/run_alto_cellranger_workflow.sh" % (countsdir)
    alto_method = "cumulus/cellranger_workflow/26"

    bash_alto = open(run_alto_file, "a")
    bash_alto.write("source /broad/software/scripts/useuse\n")
    bash_alto.write("reuse Anaconda3\n")
    bash_alto.write("use Google-Cloud-SDK\n")
    bash_alto.write("gcloud auth login %s\n" % email)
    bash_alto.write("source activate /seq/regev_genome_portal/conda_env/cumulus\n")
    for sampleid in sampledict.keys():
        input_cellranger_file = "%s/%s/input_cellranger.json" % (countsdir, sampleid)
        bash_alto.write("alto run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (alto_method, input_cellranger_file, alto_workspace, alto_counts_folder, sampleid))
    bash_alto.close()

    # Terminal commands to run alto cellranger_workflow bash script.
    print("\n## Run terminal commands below to run alto, which runs Terra cellranger_workflow pipeline via command line. ##")
    print("bash %s" % run_alto_file)


    # ## Uploading `cumulus` files to Google Cloud Storage bucket in Terra workspace

    # Make samplesheet for cumulus.
    for sampleid in sampledict.keys():
        if not os.path.isdir("%s/%s" % (resultsdir, sampleid)):
            os.mkdir("%s/%s" % (resultsdir, sampleid))
        samplesheet_cumulus_file = "%s/%s/samplesheet_cumulus.csv" % (resultsdir, sampleid)

        with open(samplesheet_cumulus_file, "w") as f:
            f.write("Sample,Location\n")
            for sample in sampledict[sampleid]:
                f.write("%s,%s/%s/%s\n" % (sample, countsbucket, sample, count_matrix_name))

    # Make input_cumulus file for cumulus.
    for sampleid in sampledict.keys():
        input_cumulus_file = "%s/%s/input_cumulus.json" % (resultsdir, sampleid)

        with open(input_cumulus_file, "w") as f:
            f.write("{\n")
            f.write("\t\"cumulus.input_file\" : \"%s/%s/samplesheet_cumulus.csv\",\n" % (resultsbucket, sampleid))
            f.write("\t\"cumulus.output_directory\" : \"%s/\",\n" % (resultsbucket))
            f.write("\t\"cumulus.output_name\" : \"%s\",\n" % (sampleid))
            f.write("\t\"cumulus.min_umis\" : %s,\n" % cumulusdict[sampleid][0])
            f.write("\t\"cumulus.min_genes\" : %s,\n" % cumulusdict[sampleid][1])
            f.write("\t\"cumulus.percent_mito\" : %s,\n" % cumulusdict[sampleid][2])
            f.write("\t\"cumulus.infer_doublets\" : true,\n")

            f.write("\t\"cumulus.run_louvain\" : false,\n")
            f.write("\t\"cumulus.run_leiden\" : true,\n")
            f.write("\t\"cumulus.leiden_resolution\" : 1,\n")
            f.write("\t\"cumulus.run_diffmap\" : false,\n")
            f.write("\t\"cumulus.perform_de_analysis\" : true,\n")
            f.write("\t\"cumulus.cluster_labels\" : \"leiden_labels\",\n")
            f.write("\t\"cumulus.annotate_cluster\" : false,\n")
            f.write("\t\"cumulus.fisher\" : true,\n")
            f.write("\t\"cumulus.t_test\" : true,\n")
            f.write("\t\"cumulus.find_markers_lightgbm\" : false,\n")

            f.write("\t\"cumulus.run_tsne\" : false,\n")
            f.write("\t\"cumulus.run_umap\" : true,\n")
            f.write("\t\"cumulus.umap_K\" : 15,\n")
            f.write("\t\"cumulus.umap_min_dist\" : 0.5,\n")
            f.write("\t\"cumulus.umap_spread\" : 1,\n")
            f.write("\t\"cumulus.plot_umap\" : \"leiden_labels\",\n")
            f.write("\t\"cumulus.output_h5ad\" : true\n")
            f.write("}\n")

    # Running bash script below to upload cumulus samplesheet and input file to Google Cloud Storage Bucket.
    print("\n## Bash script below is AUTOMATICALLY RUNNING to upload cumulus samplesheet and input file to Google Cloud Storage Bucket. ##")
    uploadcumulus_file = "%s/uploadcumulus_%s.sh" % (resultsdir, sampletracking['flowcell'].iloc[0])
    with open(uploadcumulus_file, "w") as f:
        f.write("source /broad/software/scripts/useuse\n")
        f.write("use Google-Cloud-SDK\n")
        f.write("gcloud auth login %s\n" % email)
        for sampleid in sampledict.keys():
            samplesheet_cumulus_file = "%s/%s/samplesheet_cumulus.csv" % (resultsdir, sampleid)
            input_cumulus_file = "%s/%s/input_cumulus.json" % (resultsdir, sampleid)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cumulus_file, resultsbucket, sampleid))
            f.write("gsutil cp %s %s/%s/\n" % (input_cumulus_file, resultsbucket, sampleid))
    command = "bash %s" % uploadcumulus_file
    subprocess.call(command, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    print(command)

    print("\n## Run cumulus in Terra (or run alto bash script below) and wait for jobs to complete. ##")


    # ## Running `cumulus` WDL in Terra from the command line.

    # Write bash script below to run alto to kick off cumulus jobs from command line to Terra.
    run_alto_file = "%s/run_alto_cumulus.sh" % (resultsdir)
    alto_method = "cumulus/cumulus/41" # UPDATE TO MOST CURRENT SNAPSHOT

    bash_alto = open(run_alto_file, "a")
    bash_alto.write("source /broad/software/scripts/useuse\n")
    bash_alto.write("reuse Anaconda3\n")
    bash_alto.write("use Google-Cloud-SDK\n")
    bash_alto.write("gcloud auth login %s\n" % email)
    bash_alto.write("source activate /seq/regev_genome_portal/conda_env/cumulus\n")
    for sampleid in sampledict.keys():
        input_cumulus_file = "%s/%s/input_cumulus.json" % (resultsdir, sampleid)
        bash_alto.write("alto run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (alto_method, input_cumulus_file, alto_workspace, alto_results_folder, sampleid))
    bash_alto.close()

    # Terminal commands to run alto cumulus bash script.
    print("\n## Run terminal commands below to run alto, which runs Terra cumulus pipeline via command line. ##")
    print("bash %s" % run_alto_file)


    # ## Uploading `remove-background` files to Google Cloud Storage bucket in Terra
    # For each 10x channel processed for each sample, we need to remove the ambient background.

    # Make input_cellbender file for cellbender.
    for sampleid in sampledict.keys():
        if not os.path.isdir("%s/%s" % (cellbenderdir, sampleid)):
            os.mkdir("%s/%s" % (cellbenderdir, sampleid))

        for sample in sampledict[sampleid]:
            if not os.path.isdir("%s/%s" % (cellbenderdir, sampleid)):
                os.mkdir("%s/%s" % (cellbenderdir, sampleid))
            input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbenderdir, sampleid)

            with open(input_cellbender_file, "w") as f:
                f.write("{\n")
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.total_droplets_included\" : %s,\n" % cellbenderdict[sampleid][1]),
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.z_dim\" : 100,\n"),
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.epochs\" : 150,\n"),
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.expected_cells\" : %s,\n" % cellbenderdict[sampleid][0]),
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.sample_name\" : \"%s\",\n" % sampleid)
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.hardware_preemptible_tries\" : 0,\n"),
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.input_10x_h5_file_or_mtx_directory\" : \"%s/%s/%s\",\n" % (countsbucket, sampleid, count_matrix_name))
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.output_directory\" : \"%s/%s\",\n" % (cellbenderbucket, sampleid))
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.fpr\" : \"0.01 0.05 0.1\",\n")
                f.write("\t\"cellbender_remove_background.run_cellbender_remove_background_gpu.learning_rate\" : %f\n" % float("0.00005"))
                f.write("}\n")

    # Running bash script below to upload cellbender input file to Google Cloud Storage Bucket.
    print("\n## Bash script below is AUTOMATICALLY RUNNING to upload cellbender input file to Google Cloud Storage Bucket. ##")
    uploadcellbender_file = "%s/uploadcellbender_%s.sh" % (cellbenderdir, sampletracking['flowcell'].iloc[0])
    with open(uploadcellbender_file, "w") as f:
        f.write("source /broad/software/scripts/useuse\n")
        f.write("use Google-Cloud-SDK\n")
        f.write("gcloud auth login %s\n" % email)
        for sampleid in sampledict.keys():
            for sample in sampledict[sampleid]:
                input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbenderdir, sampleid)
                f.write("gsutil cp %s %s/%s/\n" % (input_cellbender_file, cellbenderbucket, sampleid))
    command = "bash %s" % uploadcellbender_file
    subprocess.call(command, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    print(command)

    print("\n## Run cellbender in Terra (or run alto bash script below) and wait for jobs to complete. ##")


    # ## Running `remove-background-` WDL in Terra from the command line.

    # Write bash script below to run alto to kick off cumulus jobs from command line to Terra.
    run_alto_file = "%s/run_alto_cellbender.sh" % (cellbenderdir)
    alto_method = "cellbender/remove-background/11"

    bash_alto = open(run_alto_file, "a")
    bash_alto.write("source /broad/software/scripts/useuse\n")
    bash_alto.write("reuse Anaconda3\n")
    bash_alto.write("use Google-Cloud-SDK\n")
    bash_alto.write("gcloud auth login %s\n" % email)
    bash_alto.write("source activate /seq/regev_genome_portal/conda_env/cumulus\n")
    for sampleid in sampledict.keys():
        for sample in sampledict[sampleid]:
            input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbenderdir, sampleid)
            bash_alto.write("alto run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (alto_method, input_cellbender_file, alto_workspace, alto_cellbender_folder, sampleid))
    bash_alto.close()

    # Terminal commands to run alto cumulus bash script.
    print("\n## Run terminal commands below to run alto, which runs Terra remove-background pipeline via command line. ##")
    print("bash %s" % run_alto_file)


    # ## Uploading `cumulus` files (post remove-backgrounda) to Google Cloud Storage bucket in Terra workspace

    # Make samplesheet for cumulus FOLLOWING CELLBENDER.
    for sampleid in sampledict.keys():
        if not os.path.isdir("%s/%s" % (cellbender_resultsdir, sampleid)):
            os.mkdir("%s/%s" % (cellbender_resultsdir, sampleid))
        samplesheet_cellbender_cumulus_file = "%s/%s/samplesheet_cellbender_cumulus.csv" % (cellbender_resultsdir, sampleid)

        with open(samplesheet_cellbender_cumulus_file, "w") as f:
            f.write("Sample,Location\n")
            for sample in sampledict[sampleid]:
                f.write("%s,%s/%s/%s/%s_%s\n" % (sample, cellbenderbucket, sampleid, sampleid, sampleid, cellbender_matrix_name))


    # Make input_cumulus file for cumulus.
    for sampleid in sampledict.keys():
        input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_resultsdir, sampleid)

        with open(input_cellbender_cumulus_file, "w") as f:
            f.write("{\n")
            f.write("\t\"cumulus.input_file\" : \"%s/%s/samplesheet_cellbender_cumulus.csv\",\n" % (cellbender_resultsbucket, sampleid))
            f.write("\t\"cumulus.output_directory\" : \"%s/\",\n" % (cellbender_resultsbucket))
            f.write("\t\"cumulus.output_name\" : \"%s\",\n" % (sampleid))
            f.write("\t\"cumulus.min_umis\" : %s,\n" % cumulusdict[sampleid][0])
            f.write("\t\"cumulus.min_genes\" : %s,\n" % cumulusdict[sampleid][1])
            f.write("\t\"cumulus.percent_mito\" : %s,\n" % cumulusdict[sampleid][2])
            f.write("\t\"cumulus.infer_doublets\" : true,\n")

            f.write("\t\"cumulus.run_louvain\" : false,\n")
            f.write("\t\"cumulus.run_leiden\" : true,\n")
            f.write("\t\"cumulus.leiden_resolution\" : 1,\n")
            f.write("\t\"cumulus.run_diffmap\" : false,\n")
            f.write("\t\"cumulus.perform_de_analysis\" : true,\n")
            f.write("\t\"cumulus.cluster_labels\" : \"leiden_labels\",\n")
            f.write("\t\"cumulus.annotate_cluster\" : false,\n")
            f.write("\t\"cumulus.fisher\" : true,\n")
            f.write("\t\"cumulus.t_test\" : true,\n")
\            f.write("\t\"cumulus.find_markers_lightgbm\" : false,\n")

            f.write("\t\"cumulus.run_tsne\" : false,\n")
            f.write("\t\"cumulus.run_umap\" : true,\n")
            f.write("\t\"cumulus.umap_K\" : 15,\n")
            f.write("\t\"cumulus.umap_min_dist\" : 0.5,\n")
            f.write("\t\"cumulus.umap_spread\" : 1,\n")
            f.write("\t\"cumulus.plot_umap\" : \"leiden_labels\",\n")
            f.write("\t\"cumulus.output_h5ad\" : true\n")
            f.write("}\n")

    # Running bash script below to upload cumulus samplesheet and input file to Google Cloud Storage Bucket.
    print("\n## Bash script below is AUTOMATICALLY RUNNING to upload post-cellbender cumulus samplesheet and input file to Google Cloud Storage Bucket. ##")
    uploadcellbendercumulus_file = "%s/uploadcellbendercumulus_%s.sh" % (cellbender_resultsdir, sampletracking['flowcell'].iloc[0])
    with open(uploadcellbendercumulus_file, "w") as f:
        f.write("source /broad/software/scripts/useuse\n")
        f.write("use Google-Cloud-SDK\n")
        f.write("gcloud auth login %s\n" % email)
        for sampleid in sampledict.keys():
            samplesheet_cellbender_cumulus_file = "%s/%s/samplesheet_cellbender_cumulus.csv" % (cellbender_resultsdir, sampleid)
            input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_resultsdir, sampleid)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cellbender_cumulus_file, cellbender_resultsbucket, sampleid))
            f.write("gsutil cp %s %s/%s/\n" % (input_cellbender_cumulus_file, cellbender_resultsbucket, sampleid))
    command = "bash %s" % uploadcellbendercumulus_file
    subprocess.call(command, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    print(command)

    print("\n## Run post-cellbender cumulus in Terra (or run alto bash script below) and wait for jobs to complete. ##")


    # ## Running `cumulus` WDL in Terra from the command line.

    # Write bash script below to run alto to kick off cumulus jobs from command line to Terra.
    run_alto_file = "%s/run_alto_cellbender_cumulus.sh" % (cellbender_resultsdir)
    alto_method = "cumulus/cumulus/41"

    bash_alto = open(run_alto_file, "a")
    bash_alto.write("source /broad/software/scripts/useuse\n")
    bash_alto.write("reuse Anaconda3\n")
    bash_alto.write("use Google-Cloud-SDK\n")
    bash_alto.write("gcloud auth login %s\n" % email)
    bash_alto.write("source activate /seq/regev_genome_portal/conda_env/cumulus\n")
    for sampleid in sampledict.keys():
        input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_resultsdir, sampleid)
        bash_alto.write("alto run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (alto_method, input_cellbender_cumulus_file, alto_workspace, alto_results_folder, sampleid))
    bash_alto.close()

    # Terminal commands to run alto cumulus bash script.
    print("\n## Run terminal commands below to run alto, which runs Terra cumulus pipeline via command line. ##")
    print("bash %s" % run_alto_file)


    # ## Downloading count matrices from `cellranger_workflow` and clustering results from `cumulus` in the Terra workspace

    for sampleid in sampledict.keys():
        # Bash script to download counts folders and results folders.
        download_file = "%s/%s/download.sh" % (resultsdir, sampleid)

        with open(download_file, "w") as f:
            f.write("source /broad/software/scripts/useuse\n")
            f.write("use Google-Cloud-SDK\n")

            # Download cellranger_workflow molecule files for each sample.
            for sample in sampledict[sampleid]:
                if not os.path.isdir("%s/%s" % (countsdir, sample)):
                    os.mkdir("%s/%s/%s" % (countsdir, sample))
                f.write("gsutil cp %s/%s/%s %s/%s\n" % (countsbucket, sample, count_matrix_name, countsdir, sample))
                f.write("gsutil cp %s/%s/molecule_info.h5 %s/%s\n" % (countsbucket, sample, countsdir, sample))
                f.write("gsutil cp %s/%s/metrics_summary.csv %s/%s\n" % (countsbucket, sample, countsdir, sample))
                f.write("gsutil cp %s/%s/web_summary.html %s/%s\n" % (countsbucket, sample, countsdir, sample))

            # Download cumulus results folder.
            f.write("\ngsutil -m cp -r %s/%s/* %s/%s/\n" % (resultsbucket, sample, resultsdir, sample))

    # Write to bash file that will download the counts matrices and results folders once Cloud pipelines finish running.
    run_download_file = "%s/run_download_%s.sh" % (resultsdir, sampletracking['flowcell'].iloc[0])
    bash_download = open(run_download_file, "w")

    # Terminal commands to run bash script.
    print("\n## Copy and run terminal commands below to run bash script that downloads file from Google Cloud Storage Bucket. ##")
    print("use UGER")
    print("ish -l h_vmem=4g -pe smp 3 -R y -binding linear:3")
    for sampleid in sampledict.keys():
        download_file = "%s/%s/download.sh" % (resultsdir, sampleid)
        # print("bash %s" % download_file)
        bash_download.write("bash %s\n" % download_file)
    print("bash %s\n" % run_download_file)

    bash_download.close()

    # DOWNLOAD CELLBENDER RESULTS
    for sampleid in sampledict.keys():
        if not os.path.isdir("%s/%s" % (cellbenderdir, sampleid)):
            os.mkdir("%s/%s" % (cellbenderdir, sampleid))

        # Bash script to download counts folders and results folders.
        download_cellbender_file = "%s/%s/download.sh" % (cellbender_resultsdir, sampleid)

        with open(download_cellbender_file, "w") as f:
            f.write("source /broad/software/scripts/useuse\n")
            f.write("use Google-Cloud-SDK\n")

            # Download cellbender results
            f.write("\ngsutil -m cp -r %s/%s/* %s/%s/\n" % (cellbenderbucket, sampleid, cellbenderdir, sampleid))

            # Download cumulus results folder.
            f.write("\ngsutil -m cp -r %s/%s/* %s/%s/\n" % (cellbender_resultsbucket, sampleid, cellbender_resultsdir, sampleid))

    # Write to bash file that will download the counts matrices and results folders once Cloud pipelines finish running.
    run_cellbender_download_file = "%s/run_download_%s.sh" % (cellbender_resultsdir, sampletracking['flowcell'].iloc[0])
    bash_cellbender_download = open(run_cellbender_download_file, "w")

    # Terminal commands to run bash script.
    print("\n## Copy and run terminal commands below to run bash script that downloads post-cellbender files from Google Cloud Storage Bucket. ##")
    print("use UGER")
    print("ish -l h_vmem=4g -pe smp 3 -R y -binding linear:3")
    for sampleid in sampledict.keys():
        download_cellbender_file = "%s/%s/download.sh" % (cellbender_resultsdir, sampleid)
        # print("bash %s" % download_file)
        bash_cellbender_download.write("bash %s\n" % download_cellbender_file)
    print("bash %s\n" % run_cellbender_download_file)

    bash_cellbender_download.close()
