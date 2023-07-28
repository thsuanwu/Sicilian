# Wrapper script for running SICILIAN

import glob
import os
import subprocess
import sys
import argparse
import time
import argparse

def detect_10x_fastq_files(data_path, r1_suffix="_R1_001.fastq.gz", r2_suffix="_R2_001.fastq.gz"):
    """Automatically detect R1 and R2 fastq files generated from 10x demultiplexing.

    Args:
        data_path (str): Path to the directory containing the fastq files.
        r1_suffix (str, optional): Suffix for R1 files. Defaults to "_R1_001.fastq.gz".
        r2_suffix (str, optional): Suffix for R2 files. Defaults to "_R2_001.fastq.gz".

    Returns:
        tuple: A tuple containing the paths to R1 and R2 fastq files if both are found, or None otherwise.
    """
    r1_files = glob.glob(os.path.join(data_path, f"*{r1_suffix}"))
    r2_files = glob.glob(os.path.join(data_path, f"*{r2_suffix}"))

    if len(r1_files) == 1 and len(r2_files) == 1:
        return r1_files[0], r2_files[0]
    else:
        print("Error: R1 or R2 files not found or multiple files detected.")
        return None

def sbatch_file(file_name,out_path, name, job_name, time, mem, command, dep="", dep_type = "afterok"):
  """Write sbatch script given parameters"""
  job_file = open(file_name, "w")
  job_file.write("#!/bin/bash\n#\n")
  job_file.write("#SBATCH --job-name=" + job_name + "\n")
  job_file.write("#SBATCH --output={}{}/log_files/{}.%j.out\n".format(out_path, name,job_name))
  job_file.write("#SBATCH --error={}{}/log_files/{}.%j.err\n".format(out_path, name,job_name))
  job_file.write("#SBATCH --time={}\n".format(time))
  #job_file.write("#SBATCH --qos=high_p\n")
  job_file.write("#SBATCH --partition=batch\n")
  job_file.write("#SBATCH --account=mnicolls\n")
#  job_file.write("#SBATCH --account=horence\n")
#  job_file.write("#SBATCH --partition=nih_s10\n")
  job_file.write("#SBATCH --mail-user=thsuanwu@stanford.edu")
  job_file.write("#SBATCH --mail-type=ALL")
  job_file.write("#SBATCH --nodes=1\n")
  job_file.write("#SBATCH --mem={}\n".format(mem)) 
  if dep != "":
    job_file.write("#SBATCH --dependency={}:{}\n".format(dep_type,dep))
    job_file.write("#SBATCH --kill-on-invalid-dep=yes\n")
  job_file.write("date\n")
  job_file.write(command + "\n")
  job_file.write("date\n")
  job_file.close()


def GLM(out_path, name, gtf_file, single, tenX, stranded_library, domain_file, exon_pickle_file, splice_pickle_file, dep = ""):
  """Run the GLM script to compute the statistical scores for junctions in the class input file"""
  command = "Rscript scripts/GLM_script_light.R {}{}/ {} ".format(out_path, name, gtf_file)
  if single:
    command += " 1 "
  else:
    command += " 0 "
  if tenX:
    command += " 1 "
  else:
    command += " 0 "
  if stranded_library:
    command += " 1 "
  else:
    command += " 0 "
  command += "{} {} {} ".format(domain_file, exon_pickle_file, splice_pickle_file)
  sbatch_file("run_GLM.sh", out_path, name,"GLM_{}".format(name), "48:00:00", "150Gb", command, dep=dep)  # used 200Gb for CML 80Gb for others and 300 for 10x blood3 
  return submit_job("run_GLM.sh")

def whitelist(data_path,out_path, name, bc_pattern, r_ends):
  command = "mkdir -p {}{}\n".format(out_path, name)
  command += "umi_tools whitelist "
  command += "--stdin {}{}{} ".format(data_path, name, r_ends[0])
  command += "--bc-pattern={} ".format(bc_pattern)
  command += "--log2stderr > {}{}_whitelist.txt ".format(data_path,name)
  command += "--plot-prefix={}{} ".format(data_path, name)
 # command += "--knee-method=density "
  sbatch_file("run_whitelist.sh",out_path, name, "whitelist_{}".format(name), "24:00:00", "20Gb", command)
  return submit_job("run_whitelist.sh")

def extract(out_path, data_path, name, bc_pattern, r_ends, dep = ""):
  command = "umi_tools extract "
  command += "--bc-pattern={} ".format(bc_pattern)
  command += "--stdin {}{}{} ".format(data_path, name, r_ends[0])
  command += "--stdout {}{}_extracted{} ".format(data_path, name, r_ends[0])
  command += "--read2-in {}{}{} ".format(data_path, name, r_ends[1])
  command += "--read2-out={}{}_extracted{} ".format(data_path, name, r_ends[1])
#  command += "--read2-stdout "
#  command += "--filter-cell-barcode "
  command += "--whitelist={}{}_whitelist.txt ".format(data_path, name)
  command += "--error-correct-cell "
  sbatch_file("run_extract.sh", out_path, name,"extract_{}".format(name), "24:00:00", "20Gb", command, dep = dep)
  return submit_job("run_extract.sh")


def class_input(out_path, name, gtf_file, annotator_file, tenX, single, stranded_library, dep=""):
  """Run script to create class input file"""
  command = "python3 scripts/light_class_input.py --outpath {}{}/ --gtf {} --annotator {} --bams ".format(out_path, name, gtf_file,annotator_file) 
  if single:
    command += "{}{}/2Aligned.out.bam ".format(out_path,name)
  else:
    command += "{}{}/1Aligned.out.bam ".format(out_path,name)
    command += "{}{}/2Aligned.out.bam ".format(out_path,name)
  if tenX:
    command += "--UMI_bar "
#  if stranded_library:
  command += "--stranded_library "
  if not single:
    command += "--paired "
  sbatch_file("run_class_input.sh", out_path, name,"class_input_{}".format(name), "48:00:00", "200Gb", command, dep=dep)  # 96:00:00, and 210 Gb for Lu, 100 for others
  return submit_job("run_class_input.sh")

def STAR_map(out_path, data_path, name, r_ends, gzip, single, gtf_file, tenX, star_path, star_ref_path, dep = ""):
  """Run script to perform mapping job for STAR"""
  print(single)
  # Pattern matching to find R1 and R2 files
  r1_files = glob.glob(os.path.join(data_path, f"{name}_R1*fastq*"))
  r2_files = glob.glob(os.path.join(data_path, f"{name}_R2*fastq*"))

  print("Detected R1 files:")
  print(r1_files)
  print("Detected R2 files:")
  print(r2_files)


  # Automatically detect single or paired-end data if not specified
  if not single:
      single = not (r1_files and r2_files)
      print(f"Automatically detected single-end: {single}")

  if single:
      # Single-end data
      if r1_files:
          read_files_in = [r1_files[0]]
      else:
          raise ValueError("R1 file not found for single-end data.")
  else:
      # Paired-end data
      if r1_files and r2_files:
          read_files_in = [r1_files[0], r2_files[0]]
      elif r1_files or r2_files:
          raise ValueError("Both R1 and R2 files should be present for paired-end data.")
      else:
          raise ValueError("No matching R1 and R2 files found.")

  # Construct the call to STAR
  command = "{} --runThreadN 4 ".format(star_path)
  command += "--genomeDir {} ".format(star_ref_path)
  command += "--readFilesIn {} ".format(" ".join(read_files_in))

  if gzip:
      command += "--readFilesCommand zcat "

  command += "--twopassMode Basic "
  command += "--alignIntronMax 1000000 "
  command += "--outFileNamePrefix {}{}/ ".format(out_path, name)
  command += "--outSAMtype BAM Unsorted "
  command += "--outSAMattributes All "
  command += "--chimOutType WithinBAM SoftClip Junctions "
  command += "--chimJunctionOverhangMin 10 "
  command += "--chimSegmentReadGapMax 0 "
  command += "--chimOutJunctionFormat 1 "
  command += "--chimSegmentMin 12 "
  command += "--quantMode GeneCounts "
  command += "--sjdbGTFfile {} ".format(gtf_file)
  command += "--outReadsUnmapped Fastx \n\n"

  # Submit the job
  sbatch_file_name = "run_map.sh"
  sbatch_file(sbatch_file_name, out_path, name, "map_{}".format(name), "24:00:00", "60Gb", command, dep=dep)
  return submit_job(sbatch_file_name)
    
def STAR_map_depracated(out_path, data_path, name, r_ends, gzip, single, gtf_file, tenX, star_path, star_ref_path, dep = ""):
  """Run script to perform mapping job for STAR"""

  command = "mkdir -p {}{}\n".format(out_path, name)
  command += "{} --version\n".format(star_path)
  if single:
    l = 1
  else:
    l = 0
  for i in range(l,2):
    command += "{} --runThreadN 4 ".format(star_path)
    command += "--genomeDir {} ".format(star_ref_path)
    if tenX:
      command += "--readFilesIn {}{}_extracted{} ".format(data_path, name, r_ends[i])
    else:
      command += "--readFilesIn {}{}{} ".format(data_path, name, r_ends[i])
    if gzip:
      command += "--readFilesCommand zcat "
    command += "--twopassMode Basic "
    command += "--alignIntronMax 1000000 "
    command += "--outFileNamePrefix {}{}/{} ".format(out_path, name, i + 1)
    command += "--outSAMtype BAM Unsorted "
    command += "--outSAMattributes All "
    command += "--chimOutType WithinBAM SoftClip Junctions "
    command += "--chimJunctionOverhangMin 10 "
    command += "--chimSegmentReadGapMax 0 "
    command += "--chimOutJunctionFormat 1 "
    command += "--chimSegmentMin 12 "
  #  command += "--chimMultimapNmax 20 "
    command += "--chimScoreJunctionNonGTAG -4 "
    command += "--chimNonchimScoreDropMin 10 "
  #  command += "--chimMultimapScoreRange 3 "
    command += "--quantMode GeneCounts "
    command += "--sjdbGTFfile {} ".format(gtf_file)
    command += "--outReadsUnmapped Fastx \n\n"
  sbatch_file("run_map.sh", out_path, name,"map_{}".format(name), "24:00:00", "60Gb", command, dep = dep)
  return submit_job("run_map.sh")


def submit_job(file_name):
  """Submit sbatch job to cluster"""
  status, job_num = subprocess.getstatusoutput("sbatch {}".format(file_name))
  if status == 0:
    print("{} ({})".format(job_num, file_name))
    return job_num.split()[-1]
  else:
    print("Error submitting job {} {} {}".format(status, job_num, file_name))

def main():

##########################################################################################
################## Input arguments that should be set by the user  ########################
###########################################################################################
  data_path = "/scratch/PI/horence/Roozbeh/single_cell_project/data/TSP2_SS2/RUN2/"
  out_dir = "/scratch/PI/horence/Roozbeh/single_cell_project/output"
  run_name = "test"
  r_ends = ["_R1_001.fastq.gz", "_R2_001.fastq.gz"]
  names = ["TSP2_SI_distal_SS2_B114584_B133323_Epithelial_D18_S42"]
  star_path = "/oak/stanford/groups/horence/Roozbeh/software/STAR-2.7.5a/bin/Linux_x86_64/STAR"
  star_ref_path = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/hg38_ERCC_STAR_2.7.5.a"
  gtf_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/ucsc_known_genes/grch38_known_genes.gtf"
  annotator_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/hg38_refseq.pkl"
  exon_pickle_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/hg38_refseq_exon_bounds.pkl"
  splice_pickle_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/SICILIAN_references/human/hg38_refseq_splices.pkl"
  domain_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/ucscGenePfam.txt"
  single = False
  tenX = False
  stranded_library = False
  bc_pattern = "C"*16 + "N"*12
#########################################################################################
#########################################################################################
#########################################################################################

  parser = argparse.ArgumentParser(description="Wrapper script for running SICILIAN")
  parser.add_argument("--runs", nargs="*", help="List of runs to process (optional)")
  parser.add_argument("--data_path", type=str, required=True, help="Path to the directory that contains the fastq files for the input RNA-Seq data.")
  parser.add_argument("--out_dir", type=str, required=True, help="Path to the directory that will contain the folder specified by project_name for SICILIAN output files.")
  parser.add_argument("--project_name", type=str, required=True, help="Folder name for the SICILIAN output files.")
  parser.add_argument("--r_ends", nargs='+', type=str, required=True, help="List of unique endings for the file names of R1 and R2 fastq files. Example: --r_ends _1.fastq.gz _2.fastq.gz")
  #parser.add_argument("--names", nargs='+', type=str, required=False, help="List of sample names.")
  parser.add_argument("--star_path", type=str, required=True, help="Path to the STAR executable file.")
  parser.add_argument("--star_ref_path", type=str, required=True, help="Path to the STAR index files.")
  parser.add_argument("--gtf_file", type=str, required=True, help="Path to the GTF file used as the reference annotation file for the genome assembly.")
  parser.add_argument("--annotator_file", type=str, required=True, help="Path to the annotation_name.pkl file.")
  parser.add_argument("--exon_pickle_file", type=str, required=False, help="Path to the annotation_name_exon_bounds.pkl file.")
  parser.add_argument("--splice_pickle_file", type=str, required=False, help="Path to the annotation_name_splices.pkl file.")
  parser.add_argument("--domain_file", type=str, required=False, help="Path to the reference file for annotated protein domains downloaded from UCSC used for finding missing and inserted protein domains in the splice junction.")
  parser.add_argument("--single", action="store_true", help="Set this flag if the data is single-end.")
  parser.add_argument("--tenX", action="store_true", help="Set this flag if the input RNA-Seq data is 10x.")
  parser.add_argument("--stranded_library", action="store_true", help="Set this flag if input RNA-Seq data is based on a stranded library.")
# parser.add_argument("--bc_pattern", type=str, default="C"*16 + "N"*12, help="Barcode pattern for 10x data.")

  args = parser.parse_args()

  # Use the arguments from command line
  data_path = args.data_path
  out_dir = args.out_dir
  project_name = args.project_name
  r_ends = args.r_ends
  #names = args.names
  star_path = args.star_path
  star_ref_path = args.star_ref_path
  gtf_file = args.gtf_file
  annotator_file = args.annotator_file
  exon_pickle_file = args.exon_pickle_file
  splice_pickle_file = args.splice_pickle_file
  domain_file = args.domain_file
  single = args.single
  tenX = args.tenX
  stranded_library = args.stranded_library
# bc_pattern = args.bc_pattern

## Toggles for deciding which steps in SICILIAN should be run #####
  run_whitelist = False
  run_extract = False
  run_map = True
  run_class = True
  run_GLM = True
###################################################################


  project_path = out_dir + "/{}/".format(project_name) 

# If the --runs argument is provided, use those runs; otherwise, detect all runs in the data_path
  if args.runs:
    runs_to_process = args.runs
  else:
    runs_to_process = [os.path.basename(run_path) for run_path in glob.glob(os.path.join(data_path, "*")) if os.path.isdir(run_path)]

  print("Detected runs:")
  print(runs_to_process)
  
  total_jobs = []
  total_job_names = []

  for run in runs_to_process:
    print(f"Processing run: {run}")
    name = run


    data_path = os.path.join(data_path, name, '')
    print(f"Data path: {data_path}")


    # Create a directory for the current run
    #out_parent_dir = os.path.join(project_path)
    #print(out_parent_dir)
    os.makedirs(os.path.join(project_path, name), exist_ok=True)
    
    if not single:
      run_whitelist = False
      run_extract = False

    if r_ends[0].split(".")[-1] == "gz":
      gzip = True
    else:
      gzip = False


    jobs = []
    job_nums = []
    #for name in names:
    #  jobs = []
    #  job_nums = []

    if not os.path.exists("{}{}/log_files".format(project_path, name)):
      os.makedirs("{}{}/log_files".format(project_path, name))

    if run_whitelist:
      whitelist_jobid = whitelist(data_path,project_path, name, bc_pattern, r_ends)
      jobs.append("whitelist_{}.{}".format(name, whitelist_jobid))
      job_nums.append(whitelist_jobid)
    else:
      whitelist_jobid = ""

    if run_extract:
      extract_jobid = extract(project_path, data_path, name, bc_pattern, r_ends, dep = ":".join(job_nums))
      jobs.append("extract_{}.{}".format(name, extract_jobid))
      job_nums.append(extract_jobid)
    else:
      extract_jobid = ""
    if run_map:
      map_jobid = STAR_map(project_path, data_path, name, r_ends, gzip, single, gtf_file, tenX, star_path, star_ref_path, dep = ":".join(job_nums))
      jobs.append("map_{}.{}".format(name,map_jobid))
      job_nums.append(map_jobid)

    if run_class:
      class_input_jobid = class_input(project_path, name, gtf_file, annotator_file, tenX, single, stranded_library, dep=":".join(job_nums))
      jobs.append("class_input_{}.{}".format(name,class_input_jobid))
      job_nums.append(class_input_jobid)
    else:
      class_input_jobid = ""

    if run_GLM:
      GLM_jobid = GLM(project_path, name, gtf_file, single, tenX, stranded_library, domain_file, exon_pickle_file, splice_pickle_file, dep=":".join(job_nums))
      jobs.append("GLM_{}.{}".format(name,GLM_jobid))
      job_nums.append(GLM_jobid)
    else:
      GLM_jobid =  ""
    break

main()
