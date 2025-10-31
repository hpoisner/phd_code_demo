version 1.0

## Make Cohort Specific PCs
## This will use custom variant counts per chromosome

struct chrom_inputs {
    File bed_file
    File bim_file
    File fam_file
    String chr
    Int var_count
}


workflow CustomPCs {
    input {
        Array[chrom_inputs] chroms
        File keep
        String cohort_string
        Int pc_count

        String plink_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.12"
    }


    scatter (chrom in chroms) {
        call ExtractVariants {
            input:
            bed_file = chrom.bed_file,
            bim_file = chrom.bim_file,
            fam_file = chrom.fam_file,
            chr = chrom.chr,
            var_count = chrom.var_count,
            keep = keep,
            docker = plink_docker
        }
    }

    call MergeBfiles_RunPCA {
        input:
        bed_files = ExtractVariants.bed_out,
        bim_files = ExtractVariants.bim_out,
        fam_files = ExtractVariants.fam_out,
        cohort_string = cohort_string,
        pc_count = pc_count,
        docker=plink_docker
    }

    output {
         File eigenvec = MergeBfiles_RunPCA.eigenvec
         File eigenval  = MergeBfiles_RunPCA.eigenval
    }
}

task ExtractVariants {
    input {
        File bed_file
        File bim_file
        File fam_file
        String chr
        Int var_count
        File keep
        String docker

        Int? disk_size_override
        Float memory = 10
        Int cpu = 8
        Int preemptible = 1
        Int maxRetries = 1
    }    
    Float genotype_size = size(bed_file, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * genotype_size)])
    String bfile = basename(bed_file, ".bed")

    command <<<
    
    plink2 \
    --bed ~{bed_file} \
    --bim ~{bim_file} \
    --fam ~{fam_file} \
    --keep ~{keep} \
    --maf 0.05 --geno 0.01 \
    --autosome \
    --snps-only "just-acgt" \
    --make-bed --out filtered_data_~{chr} \
    --threads 8 

    plink2 \
    --bfile filtered_data_~{chr} \
    --indep-pairwise 50 5 0.1 \
    --thin-count ~{var_count} \
    --out pruned_data_~{chr} \
    --threads 8 

    plink2 \
    --bfile filtered_data_~{chr} \
    --extract pruned_data_~{chr}.prune.in \
    --make-bed \
    --threads 8 \
    -out final_data_~{chr} 
    >>>

    output {
        File bed_out = "final_data_${chr}.bed"
        File bim_out = "final_data_${chr}.bim"
        File fam_out = "final_data_${chr}.fam"
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

task MergeBfiles_RunPCA {
    input {
        Array[File] bed_files
        Array[File] bim_files
        Array[File] fam_files
        String cohort_string
        Int pc_count
        String docker

        Int disk = 100
        Float memory = 300
        Int cpu = 32
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

        # Loop over the arrays starting from index 0
        for ((i=0; i<${#bed_files_array[@]}; i++)); do
            echo "${bed_files_array[i]} ${bim_files_array[i]} ${fam_files_array[i]}" >> merge_list.txt
        done

        # Perform the merge using the merge list
        plink \
        --merge-list merge_list.txt \
        --make-bed \
        --threads 8 \
        --out merged_bfiles

        plink2 --bfile merged_bfiles --pca approx ~{pc_count} --out ~{cohort_string} --threads 32 --memory 200000

    >>>

    output {
        File eigenvec = "${cohort_string}.eigenvec"
        File eigenval  = "${cohort_string}.eigenval"
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