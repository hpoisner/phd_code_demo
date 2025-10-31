version 1.0

workflow agd_regenie_qc_multi_ancestry {
    input {
        Array[File] chr_pgen
        Array[File] chr_pvar
        Array[File] chr_psam
        File pgenx
        File pvarx
        File psamx
        File keep_all
        File keep_eur
        File keep_afr
        File keep_sas
        File keep_eas
        File keep_amr
        File keep_ad_maj
        Array[String] chr_names
        String output_prefix
  }


    scatter (i in range(length(chr_pgen))) {
        call Step1_and_Step2_filter {
            input:
                pgen = chr_pgen[i],
                pvar = chr_pvar[i],
                psam = chr_psam[i],
                chr_name = chr_names[i],
                keep = keep_all,
                keep_afr = keep_afr,
                keep_eur= keep_eur,
                keep_sas = keep_sas,
                keep_eas = keep_eas,
                keep_ad_maj = keep_ad_maj,
                keep_amr = keep_amr,
                output_prefix = output_prefix
            }
  }
    call Step1_and_Step2_filter_chrX {
        input: 
            pgen = pgenx,
            pvar = pvarx,
            psam = psamx,
            keep = keep_all,
            keep_afr = keep_afr,
            keep_eur = keep_eur,
            keep_sas = keep_sas,
            keep_eas = keep_eas,
            keep_ad_maj = keep_ad_maj,
            keep_amr = keep_amr,
            output_prefix = output_prefix 
    }

    call MergeStep1 {
        input:
            bed_files = Step1_and_Step2_filter.bed_file,
            bim_files = Step1_and_Step2_filter.bim_file,
            fam_files = Step1_and_Step2_filter.fam_file,
            chrx_bed = Step1_and_Step2_filter_chrX.bed_file,
            chrx_bim = Step1_and_Step2_filter_chrX.bim_file,
            chrx_fam = Step1_and_Step2_filter_chrX.fam_file
    }

    output {
        Array[File] snplists = Step1_and_Step2_filter.snplist
        File snplist = Step1_and_Step2_filter_chrX.snplist
        File merged_bed = MergeStep1.bed_file
        File merged_bim = MergeStep1.bim_file
        File merged_fam = MergeStep1.fam_file
    }
}

    task Step1_and_Step2_filter {
        input {
            File pgen
            File pvar
            File psam
            File keep
            File keep_ad_maj
            File keep_afr
            File keep_amr
            File keep_eur
            File keep_sas
            File keep_eas
            String chr_name
            String output_prefix 
            Int cpu = 8
            Int preemptible = 1
            Int maxRetries = 1
        }

        Float pgen_gb = size(pgen, "GB")
        Int disk_gb = ceil(pgen_gb * 4.0) + 10

        command <<<
        ##Step 1 basic qc
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep} \
            --maf 0.25 \
            --mac 25000 \
            --hwe 1e-15 \
            --geno 0.1 \
            --make-bed \
            --out ~{output_prefix}_~{chr_name}

        ##Step 2 basic qc AFR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_afr} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_afr_filtered


        ##Step 2 basic qc EUR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_eur} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_eur_filtered


        ##Step 2 basic qc SAS
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_sas} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_sas_filtered


        ##Step 2 basic qc EAS
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_eas} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_eas_filtered


        ##Step 2 basic qc AMR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_amr} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_amr_filtered

        ##Step 2 basic qc No Madj ADM
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_ad_maj} \
            --hwe 1e-15 \
            --write-snplist \
            --out ~{chr_name}_ad_maj_filtered

        cat ~{chr_name}_*_filtered.snplist | sort | uniq -d > HWE_filtered_~{chr_name}.txt

        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep} \
            --extract HWE_filtered_~{chr_name}.txt \
            --maf 0.01 \
            --geno 0.1 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --out ~{chr_name}_filtered_maf0.01_geno0.01_hwe1e15
        
        >>>

        output {
            File snplist = "~{chr_name}_filtered_maf0.01_geno0.01_hwe1e15.snplist"
            File bed_file = "~{output_prefix}_~{chr_name}.bed"
            File bim_file = "~{output_prefix}_~{chr_name}.bim"
            File fam_file = "~{output_prefix}_~{chr_name}.fam"
        }

        runtime {
            docker: "briansha/plink2:terra"
            memory: "50G"
            cpu: cpu
            disks: "local-disk ~{disk_gb} HDD"
            preemptible: preemptible
            maxRetries: maxRetries
        }
    }

    task MergeStep1 {
        input {
            Array[File] bed_files
            Array[File] bim_files
            Array[File] fam_files
            File chrx_bed
            File chrx_bim
            File chrx_fam
            String output_prefix
            Int disk = 100
            Float memory = 10
            Int cpu = 8
            Int preemptible = 1
            Int maxRetries = 1
        }


        command <<<
            # Create the merge list file for PLINK
            rm -f merge_list.txt

            # Initialize arrays from WDL inputs
            bed_files_array=(~{sep=' ' bed_files})
            bim_files_array=(~{sep=' ' bim_files})
            fam_files_array=(~{sep=' ' fam_files})

                    # Loop over the arrays starting from index 1
            for ((i=1; i<${#bed_files_array[@]}; i++)); do
                echo "${bed_files_array[i]} ${bim_files_array[i]} ${fam_files_array[i]}" >> merge_list.txt
            done

            echo ~{chrx_bed} ~{chrx_bim} ~{chrx_fam} >> merge_list.txt
            # Perform the merge using the merge list
            plink \
            --bfile ${bed_files_array[0]} \
            --merge-list merge_list.txt \
            --thin-count 500000 \
            --make-bed \
            --out ~{output_prefix}_merged_step1

        >>>

        output {
            File bed_file = "~{output_prefix}_merged_step1.bed"
            File bim_file = "~{output_prefix}_merged_step1.bim"
            File fam_file = "~{output_prefix}_merged_step1.fam"

        }

        runtime {
            docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.12"
            memory: memory + " GiB"
            disks: "local-disk " + disk + " HDD"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
        }

    }

    task Step1_and_Step2_filter_chrX {
        input {
            File pgen
            File pvar
            File psam
            File keep
            File keep_ad_maj
            File keep_afr
            File keep_amr
            File keep_eur
            File keep_sas
            File keep_eas
            String output_prefix 
            Int cpu = 8
            Int preemptible = 1
            Int maxRetries = 1
        }

        Float pgen_gb = size(pgen, "GB")
        Int disk_gb = ceil(pgen_gb * 4.0) + 10

        command <<<
        ##Step 1 basic qc
        ##Step 1 basic qc
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep} \
            --split-par hg38 \
            --maf 0.25 \
            --mac 25000 \
            --hwe 1e-15 \
            --geno 0.1 \
            --make-bed \
            --out ~{output_prefix}_chrX

        ##Step 2 basic qc AFR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_afr} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            -no-id-header \
            --split-par hg38 \
            --out chrX_afr_filtered


        ##Step 2 basic qc EUR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_eur} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_eur_filtered


        ##Step 2 basic qc SAS
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_sas} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_sas_filtered


        ##Step 2 basic qc EAS
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_eas} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_eas_filtered


        ##Step 2 basic qc AMR
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_amr} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_amr_filtered

        ##Step 2 basic qc No Madj ADM
        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep_ad_maj} \
            --hwe 1e-15 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_ad_maj_filtered

        cat chrX_*_filtered.snplist | sort | uniq -d > HWE_filtered_chrX.txt

        plink2 \
            --pgen ~{pgen} \
            --pvar ~{pvar} \
            --psam ~{psam} \
            --keep ~{keep} \
            --extract HWE_filtered_chrX.txt \
            --maf 0.01 \
            --geno 0.1 \
            --write-snplist \
            --write-samples \
            --no-id-header \
            --split-par hg38 \
            --out chrX_filtered_maf0.01_geno0.01_hwe1e15        

        >>>

        output {
            File snplist = "chrX_filtered_maf0.01_geno0.01_hwe1e15.snplist"
            File bed_file = "~{output_prefix}_chrX.bed"
            File bim_file = "~{output_prefix}_chrX.bim"
            File fam_file = "~{output_prefix}_chrX.fam"
        }

        runtime {
            docker: "briansha/plink2:terra"
            memory: "50G"
            cpu: cpu
            disks: "local-disk ~{disk_gb} HDD"
            preemptible: preemptible
            maxRetries: maxRetries
        }
    }