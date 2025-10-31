
version 1.0

## Version 02-25-2025
##
## This WDL workflow runs Regenie.
## This workflow assumes users have thoroughly read the Regenie docs for caveats and details.
## Regenie's documentation: https://rgcgithub.github.io/regenie/options/
## This is addapted from a script originally developed by Brian Sharber to work in the All of Us Cromwell server
##
## Cromwell version support - Successfully tested on v77
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

struct Step2 {
    File bed_step2
    File bim_step2
    File fam_step2
    File anno_file
    File set_file
    String chr
}

workflow Regenie {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bed_step1
        File bim_step1
        File fam_step1
        Array[Step2] chroms
        Array[String] chr_list
        String fit_bin_out_name = "fit_bin_out" # File prefix for the list of predictions produced in Step1.
        File? keep
        File? phenoFile
        File? mask_file
        String? phenoCol
        String? phenoColList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList
        Boolean apply_rint
        String trait
        File? pred
        Array[String] phenotype_names # Phenotypes you want to analyze. (Column names).

        String regenie_docker = "briansha/regenie:v4.0" # Compiled with Boost IOSTREAM: https://github.com/rgcgithub/regenie/wiki/Using-docker
        String r_base_docker = "briansha/regenie_r_base:4.1.0"  # Ubuntu 18.04, R 4.1.0, and a few Ubuntu and R packages.
    }

    call RegenieStep1WholeGenomeModel {
        input:
            bed_step1 = bed_step1,
            bim_step1 = bim_step1,
            fam_step1 = fam_step1,
            keep = keep,
            fit_bin_out_name = fit_bin_out_name,
            covarFile = covarFile,
            phenoFile = phenoFile,
            phenoColList = phenoColList,
            covarColList = covarColList,
            catCovarList = catCovarList,
            apply_rint = apply_rint,
            trait = trait,
            docker = regenie_docker
    }

    scatter (chrom in chroms) {
      call RegenieStep2AssociationTesting {
          input:
              bed_step2 = chrom.bed_step2,
              bim_step2 = chrom.bim_step2,
              fam_step2 = chrom.fam_step2,
              anno_file = chrom.anno_file,
              set_file = chrom.set_file,
              keep = keep,
              chr = chrom.chr,
              covarFile = covarFile,
              mask_file = mask_file,
              phenoFile = phenoFile,
              phenoColList = phenoColList,
              covarColList = covarColList,
              catCovarList = catCovarList,
              apply_rint = apply_rint,
              trait = trait,
              pred = RegenieStep1WholeGenomeModel.fit_bin_out,
              docker = regenie_docker,
              output_locos = RegenieStep1WholeGenomeModel.output_locos
      }
    }

    call join_Output {
      input:
        output_files = RegenieStep2AssociationTesting.test_bin_out_firth, # Refers implicitly to the entire array of files that were scattered.
        chr_list = chr_list,
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
        String fit_bin_out_name # File prefix for the list of predictions produced.

        # Basic options
        File bed_step1
        File bim_step1
        File fam_step1
        File? keep
        File? extract
        File? phenoFile
        String? phenoCol
        String? phenoColList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList

        # Options
        Int bsize = 1000  # 1000 is recommended by regenie's documentation
        Boolean apply_rint = true
        Boolean firth = false
        Boolean approx = false
        Boolean firth_se = false
        Boolean write_null_firth = false
        File? use_null_firth
        Boolean spa = false
        Boolean debug = false
        Boolean verbose = true
        String trait

      
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
        --~{trait} \
        ~{if defined(keep) then "--keep=~{keep} " else " "} \
        --phenoFile=~{phenoFile} \
        ~{if defined(phenoCol) then "--phenoCol=~{phenoCol} " else " "} \
        ~{if defined(phenoColList) then "--phenoColList=~{phenoColList} " else " "} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(covarCol) then "--covarCol=~{covarCol} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        --bsize=~{bsize} \
        ~{if firth then "--firth " else " "} \
        ~{if approx then "--approx " else " "} \
        ~{if firth_se then "--firth-se " else " "} \
        ~{if write_null_firth then "--write-null-firth " else " "} \
        ~{if defined(use_null_firth) then "--use-null-firth=~{use_null_firth} " else " "} \
        ~{if spa then "--spa " else " "} \
        ~{if apply_rint then "--apply-rint" else " "} \
        ~{if debug then "--debug " else " "} \
        ~{if verbose then "--verbose " else " "} \
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

# In Step 2, a set of imputed SNPs are tested for association using a Firth logistic regression model.
task RegenieStep2AssociationTesting {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bed_step2
        File bim_step2
        File fam_step2
        File anno_file
        File? mask_file
        File set_file
        String chr
        Array[File] output_locos

        # Basic options
        File? keep
        File? phenoFile
        String? phenoCol
        String? phenoColList
        String? phenoExcludeList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList
        String? covarExcludeList
        File? pred     # File containing predictions from Step 1
        String trait

        # Options
        Int bsize = 400      # 400 is recommended by regenie's github
        Boolean apply_rint = true
        Boolean firth = false
        Boolean approx = false
        Boolean firth_se = false
        Boolean write_null_firth = false
        File? use_null_firth
        Boolean spa = false
        Boolean debug = false
        Boolean verbose = false

  
        # Options
        Float? aaf_bins
        String? build_mask
        Boolean singleton_carrier = false
        Boolean write_mask = false
        Float? vc_maxAAF
        Float? skat_params
        Float? skat_rho
        Float? vc_MACthr
        String? joint
        Boolean skip_test = false
        String? mask_lovo
        Boolean mask_lodo = false
        Boolean write_mask_snplist = false
        Boolean check_burden_files = false
        Boolean strict_check_burden = false

        # Runtime
        String docker
        Float memory = 10
        Int? disk_size_override
        Int cpu = 8
        Int preemptible = 1
        Int maxRetries = 1
    }
    Float genotype_size = size(bed_step2, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 4.0 * genotype_size)])
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
        --phenoFile=~{phenoFile} \
        ~{if defined(phenoCol) then "--phenoCol=~{phenoCol} " else " "} \
        ~{if defined(phenoColList) then "--phenoColList=~{phenoColList} " else " "} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(covarCol) then "--covarCol=~{covarCol} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        ~{if defined(pred) then "--pred=~{pred} " else " "} \
        --bsize=~{bsize} \
        ~{if firth then "--firth " else " "} \
        ~{if approx then "--approx " else " "} \
        ~{if firth_se then "--firth-se " else " "} \
        ~{if apply_rint then "--apply-rint" else " "} \
        ~{if write_null_firth then "--write-null-firth " else " "} \
        ~{if defined(use_null_firth) then "--use-null-firth=~{use_null_firth} " else " "} \
        ~{if spa then "--spa " else " "} \
        ~{if debug then "--debug " else " "} \
        ~{if verbose then "--verbose " else " "} \
        ~{if defined(anno_file) then "--anno-file=~{anno_file} " else " "} \
        ~{if defined(set_file) then "--set-list=~{set_file} " else " "} \
        ~{if defined(mask_file) then "--mask-def=~{mask_file} " else " "} \
        ~{if defined(build_mask) then "--build-mask=~{build_mask} " else " "} \
        ~{if write_mask then "--write-mask " else " "} \
        --vc-tests skato \
        --~{trait} \
        ~{if defined(vc_maxAAF) then "--vc-maxAAF=~{vc_maxAAF} " else " "} \
        ~{if defined(skat_params) then "--skat-params=~{skat_params} " else " "} \
        ~{if defined(skat_rho) then "--skat-rho=~{skat_rho} " else " "} \
        ~{if defined(vc_MACthr) then "--vc-MACthr=~{vc_MACthr} " else " "} \
        ~{if defined(mask_lovo) then "--mask-lovo=~{mask_lovo} " else " "} \
        ~{if mask_lodo then "--mask-lodo " else " "} \
        ~{if write_mask_snplist then "--write-mask-snplist " else " "} \
        ~{if check_burden_files then "--check-burden-files " else " "} \
        ~{if strict_check_burden then "--strict-check-burden " else " "} \
        --out ~{chr}
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
      Array[String] chr_list
      String docker
      Float memory = 10
      Int? disk_size_override
      Int cpu = 1
      Int preemptible = 1
      Int maxRetries = 1
  }
  Array[File] all_output_files = flatten(output_files)
  Float regenie_files_size = size(all_output_files, "GiB")
  Int disk = select_first([disk_size_override, ceil(30.0 + 3.0 * regenie_files_size)])

  command <<<
        set -euo pipefail
        for array in ~{sep=' ' all_output_files}; do \
          for file in $array; do \
            mv $file .; \
        done done
        for pheno in ~{sep=' ' phenotype_names}; do \
          for chr in ~{sep= ' ' chr_list}; do \
            cat ${chr}_${pheno}.regenie >> $pheno.regenie; \
            rm ${chr}_${pheno}.regenie; \
        done done
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
