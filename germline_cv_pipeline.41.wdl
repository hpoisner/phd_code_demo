version 1.0

## Version 12-2022

## This workflow runs 4 steps
##DONE
## 1  Combines vcfs & Removes individual sample columns 
    ## Provide all unzipped vcfs as set and cohort name
## 2. Annotates merged vcf
## 3. Convert to txt 
## 4. Format file as desired

## Distributed under terms of the MIT License
## Copyright (c) 2022 Hannah Poisner
## Contact hannah.m.poisner@vanderbilt.edu

workflow germline_cv {
    input {
    Array[File] germline_vcfs
    String cohort_string
    File clinvar
    File clinvar_index
    File cad_snp
    File cad_snp_index
    File arr_snp
    File arr_snp_index
    File twist_panel
    File twist_panel_index
    File caddIndel
    File caddIndelTbi
    File caddSnv
    File caddSnvTbi
    File cache
    File python_script
    String species = "homo_sapiens"
    String assembly = "GRCh38"
}

call zip_vcfs {
        input:
        germline_vcfs = germline_vcfs,
        cohort_string = cohort_string
}

call annotate_vcf {
    input:
    input_vcf = zip_vcfs.merged_vcf_gz,
    input_vcf_index = zip_vcfs.merged_vcf_index,
    cohort_string = cohort_string,
    species = species,
    assembly = assembly,
    clinvar = clinvar,
    clinvar_index = clinvar_index,
    cad_snp = cad_snp,
    cad_snp_index = cad_snp_index,
    arr_snp = arr_snp,
    arr_snp_index = arr_snp_index,
    twist_panel = twist_panel,
    twist_panel_index = twist_panel_index,
    caddIndel =caddIndel,
    caddIndelTbi = caddIndelTbi,
    caddSnv = caddSnv,
    caddSnvTbi = caddSnvTbi,
    cache = cache
}

call variant_to_table {
    input:
    input_vcf = annotate_vcf.anno_vcf,
    cohort_string = cohort_string
}

call convert_to_excel {
    input:
    input_table = variant_to_table.output_txt,
    cohort_string = cohort_string,
    python_script = python_script
}

output {
    File anno_txt = variant_to_table.output_txt
    File anno_xlsx = convert_to_excel.output_xlsx
}
}


task zip_vcfs {
    input {
        Array[File] germline_vcfs
        String cohort_string
        Float memory = 3.5
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    File f = germline_vcfs[0]
    Float vcf_size = size(f, "GiB")
    Int samples = length(germline_vcfs)
    Int disk_size = ceil(15 + samples * vcf_size)
    command <<<
    set -eux -o pipefail
        for x in ~{sep=' ' germline_vcfs}
        do
            ## bgzips files
            bgzip -cf $x > $(basename "${x}").gz
            ## indexes files
            tabix -p vcf $(basename "${x}").gz
        done;
        ls *vcf.gz >> vcf_list.txt
        bcftools merge -l vcf_list.txt -Oz -o ~{cohort_string}.merged.vcf.gz
        tabix -p vcf ~{cohort_string}.merged.vcf.gz 
        bcftools sort ~{cohort_string}.merged.vcf.gz -Oz -o ~{cohort_string}.merged.sorted.vcf.gz
        tabix -p vcf ~{cohort_string}.merged.sorted.vcf.gz
    >>>

    output {
        File merged_vcf_gz = "~{cohort_string}.merged.sorted.vcf.gz"
        File merged_vcf_index = "~{cohort_string}.merged.sorted.vcf.gz.tbi"
    }

    runtime {
    	docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
		cpu: cpu
        preemptible: 1
        maxRetries: 0
    }

}


task annotate_vcf {
    input {
        File input_vcf
        File input_vcf_index
        String cohort_string
        String assembly
        File ref_fasta
        File ref_fasta_index
        File clinvar
        File clinvar_index
        File cad_snp
        File cad_snp_index
        File arr_snp
        File arr_snp_index
        File twist_panel
        File twist_panel_index
        File caddSnv
        File caddSnvTbi
        File caddIndel
        File caddIndelTbi
        File cache
        String species
        String assembly
        Float memory = 3.5
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    ## Runtime
    String docker_image = "ensemblorg/ensembl-vep:release_108.1"
    String VepDIR = "/opt/vep/src/ensembl-vep/vep"
    String FilterVEP = "/opt/vep/src/ensembl-vep/filter_vep"
    Int disk_size = ceil(5*size(cache, "GB")+size(input_vcf, "GB")+10)

    command {
        ## TODO GET CADD,FATHMM
        set -e
        mkdir -p ensembl_vep/homo_sapiens
        tar -xvzf ${cache}
        mv homo_sapiens/108_GRCh38 ensembl_vep/homo_sapiens/
    

        ## Getting CADD
        curl -O https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/108/CADD.pm 
        mv CADD.pm ensembl_vep/
        
        ls ensembl_vep

        ## Running VEP
        ${VepDIR} \
        --offline --cache \
        --species ${species} \
        --assembly ${assembly} \
        --dir_cache ensembl_vep \
        --dir_plugins ensembl_vep \
        --fork 16 \
        --no_escape \
        --fasta ${ref_fasta} \
        --format vcf \
        --vcf \
        --sift b \
        --polyphen b \
        --symbol \
        --af_gnomadg \
        --mane \
        --pick \
        --pick_order mane \
        --hgvs \
        --no_escape \
        --plugin CADD,${caddSnv},${caddIndel} \
        --custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNDN,CLNDISDB,CLNHGVS,CLNREVSTAT,CLNVC,ORIGIN \
        --custom ${arr_snp},ARR,vcf,exact,0,ARR_SNP \
        --custom ${cad_snp},CAD,vcf,exact,0,CAD_SNP \
        --custom ${twist_panel},TWIST,vcf,overlap,0,TWIST_GENE \
        --input_file ${input_vcf} \
        --output_file ${cohort_string}.annotate.vcf
    }

    output {
        File anno_vcf = "${cohort_string}.annotate.vcf"
    }

    runtime {
    	docker: docker_image
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
		cpu: cpu
        preemptible: 1
        maxRetries: 0
    }

}


task variant_to_table {
    input {
        File input_vcf
        String cohort_string
        Float memory = 3.5
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
   
    }
    Int disk_size = ceil(2*size(input_vcf, "GB")+ 10)

    command {
        vcf2table \
        --format GT \
        ~{input_vcf} \
        --filter \
        > ~{cohort_string}.annotate.txt
        
    }

    output {
        File output_txt = "${cohort_string}.annotate.txt"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: 1
        maxRetries: 1
    }
}

task convert_to_excel {
    input {
        File input_table
        File python_script
        String cohort_string
        Float memory = 3.5
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }

    Int disk_size = ceil(2*size(input_table, "GB") + 5)

    command {
        python3 ~{python_script} ~{input_table} ~{cohort_string}.annotate.xlsx
    }

    output {
        File output_xlsx = "${cohort_string}.annotate.xlsx"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: 1
        maxRetries: 1 
    }
}