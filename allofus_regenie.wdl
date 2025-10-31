

version 1.0

## Version 01-29-2025
##
## This WDL workflow runs Regenie.
## This workflow assumes users have thoroughly read the Regenie docs for caveats and details.
## Regenie's documentation: https://rgcgithub.github.io/regenie/options/
## This is addapted from a script originally developed by Brian Sharber to work in the All of Us Cromwell server

struct Step2 {
    File bed_step2
    File bim_step2
    File fam_step2
    File extract
    String chr
    String range
}

workflow Regenie {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bed_step1
        File bim_step1
        File fam_step1
        File? step1_extract
        Array[Step2] chroms
        String fit_bin_out_name = "fit_bin_out" # File prefix for the list of predictions produced in Step1.
        File? keep
        File? phenoFile
        File? covarFile
        String? covarColList
        String? catCovarList
        String phenoColList
        Boolean? apply_rint
        Boolean? bt
        Array[String] phenotype_names


        
        String regenie_docker = "briansha/regenie:v4.0"
        String r_base_docker = "briansha/regenie_r_base:4.1.0"  # Ubuntu 18.04, R 4.1.0, and a few Ubuntu and R packages.
    }

    call RegenieStep1WholeGenomeModel {
        input:
            bed_step1 = bed_step1,
            bim_step1 = bim_step1,
            fam_step1 = fam_step1,
            fit_bin_out_name = fit_bin_out_name,
            phenoColList = phenoColList,
            docker = regenie_docker
    }

    scatter (chrom in chroms) {
      call RegenieStep2AssociationTesting {
          input:
              bed_step2 = chrom.bed_step2,
              bim_step2 = chrom.bim_step2,
              fam_step2 = chrom.fam_step2,
              chr = chrom.chr,
              extract=chrom.extract,
              range=chrom.range,
              phenoColList = phenoColList,
              pred = RegenieStep1WholeGenomeModel.fit_bin_out,
              docker = regenie_docker,
              output_locos = RegenieStep1WholeGenomeModel.output_locos
      }
    }

    call join_Output {
      input:
        output_files = RegenieStep2AssociationTesting.test_bin_out_firth, # Refers implicitly to the entire array of files that were scattered.
        phenotype_names = phenotype_names,
        docker = r_base_docker
    }


    output {
         Array[File] output_regenie = join_Output.outputs
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow runs Regenie - see the README on the github for more information - https://github.com/briansha/Regenie_WDL."
    }
}

# In Step 1, the whole genome regression model is fit to the traits
# and a set of genomic predictions are produced as output.
task RegenieStep1WholeGenomeModel {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bed_step1
        File bim_step1
        File fam_step1
        File? phenoFile
        File? keep
        File? step1_extract
        File? covarFile
        String phenoColList
        String? covarColList
        String? catCovarList
        String? eventColList
        String fit_bin_out_name
        Boolean apply_rint = false
        Boolean bt = false
        Boolean t2e = false 
        Int bsize = 1000

        # Runtime
        String docker
        Float memory = 10
        Int? disk_size_override
        Int cpu = 8
        Int preemptible = 1
        Int maxRetries = 1
    }
    Float genotype_size = size(bed_step1, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * genotype_size)])
    String bfile = basename(bed_step1, ".bed")

    command <<<
        cp ~{bed_step1} .
        cp ~{bim_step1} .
        cp ~{fam_step1} .
        set -euo pipefail
        regenie \
        --step 1 \
        --threads=~{cpu} \
        --bed=~{bfile} \
        --phenoFile=~{phenoFile} \
        --phenoColList=~{phenoColList} \
        --lowmem \
        --lowmem-prefix tmp_rg \
        ~{if defined(keep) then "--keep=~{keep} " else " "} \
        ~{if defined(step1_extract) then "--extract=~{step1_extract} " else " "} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(covarColList) then "--covarCol=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        ~{if defined(eventColList) then "--eventColList=~{eventColList} " else " "} \
        ~{if bt then "--bt " else " "} \
        ~{if apply_rint then "--apply-rint" else " "} \
        ~{if t2e then "--t2e" else " "} \
        --bsize=~{bsize} \
        --out fit_bin_out
    >>>

    output {
        File fit_bin_out = "${fit_bin_out_name}_pred.list" # Refers to the list of loco files written. Lists the files as being located in the current working directory.
        Array[File] output_locos = glob("*.loco") # Writes n loco files for n phenotypes.
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

task RegenieStep2AssociationTesting {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bed_step2
        File bim_step2
        File fam_step2
        File extract
        String chr
        String range
        Array[File] output_locos
        File? phenoFile
        String phenoColList
        File? covarFile
        String? covarColList
        String? catCovarList
        String? eventColList
        String? interaction
        Boolean apply_rint = false
        Boolean bt = false
        Boolean t2e = false
        Int bsize = 400
        File? keep
        File? pred
        Boolean firth = false
        Boolean  approx = false

        # Runtime
        String docker
        Float memory = 50
        Int? disk_size_override
        Int cpu = 16
        Int preemptible = 1
        Int maxRetries = 1
    }
    Float genotype_size = size(bed_step2, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 3.0 * genotype_size)])
    String bfile = basename(bed_step2, ".bed")

    # Loco files are moved to the current working directory due to the list of predictions (pred) expecting them to be there.
    command <<<
        set -euo pipefail
        for file in ~{sep=' ' output_locos}; do \
          mv $file .; \
        done
    
        cp ~{bed_step2} .
        cp ~{bim_step2} .
        cp ~{fam_step2} .
    
        regenie \
        --step 2 \
        --threads=~{cpu} \
        --bed=~{bfile} \
        --phenoFile=~{phenoFile} \
        ~{if defined(keep) then "--keep=~{keep} " else " "} \
        --extract=~{extract} \
        --phenoFile=~{phenoFile} \
        --phenoColList=~{phenoColList} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        ~{if defined(eventColList) then "--eventColList=~{eventColList} " else " "} \
        ~{if t2e then "--t2e" else " "} \
        --pred=~{pred} \
        ~{if bt then "--bt " else " "} \
        --bsize=~{bsize} \
        ~{if firth then "--firth " else " "} \
        ~{if approx then "--approx " else " "} \
        ~{if firth then "--firth-se " else " "} \
        ~{if apply_rint then "--apply-rint" else " "} \
        --chr=~{chr} \
        --range=~{range} \
        ~{if defined(interaction) then "--interaction=~{interaction} " else " "} \
        --out ~{range}
    >>>

    output {
        Array[File] test_bin_out_firth = glob("*.regenie")
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Join together all output relating to the scattered tasks in Step 2.
# For each phenotype, one file will be made containing all analysis from Step 1 on each chromosome used in testing.
task join_Output {
  input {
      Array[Array[File]] output_files # All output from Step 2.
      Array[String] phenotype_names
      Array[String] range_list
      String docker
      Float memory = 10
      Int? disk_size_override
      Int cpu = 1
      Int preemptible = 1
      Int maxRetries = 1
  }
  Array[File] all_output_files = flatten(output_files)
  Float regenie_files_size = size(all_output_files, "GiB")
  Int disk = select_first([disk_size_override, ceil(30.0 + 2.0 * regenie_files_size)])

  command <<<
        set -euo pipefail
        for array in ~{sep=' ' all_output_files}; do \
          for file in $array; do \
            mv $file .; \
        done done
        for pheno in ~{sep=' ' phenotype_names}; do \
          for chr in ~{sep= ' ' range_list}; do \
            cat ${chr}_${pheno}.regenie >> $pheno.unfil_regenie; \
            rm ${chr}_${pheno}.regenie; \
        done done
        for pheno in ~{sep=' ' phenotype_names}; do \
            awk '!x[$0]++' $pheno.unfil_regenie > $pheno.regenie
        done
  >>>

  output {
        Array[File] outputs = glob("*.regenie")
  }

  runtime {
      docker: docker
      memory: memory + " GiB"
      disks: "local-disk " + disk + " HDD"
      cpu: cpu
      preemptible: preemptible
      maxRetries: maxRetries
  }
}
