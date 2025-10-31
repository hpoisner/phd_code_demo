version 1.0
# Example workflow

workflow printreads_immunelens {
    input {
        Array[File] cram_files
        Array[File] crai_files
        File fasta
        File fasta_index
        File dict 
        File immune_lens_rscript
    }
    
    call PrintReads {
        input:
        cram_files = cram_files,
        crai_files = crai_files,
        fasta = fasta,
        fasta_index = fasta_index,
        dict = dict 
    }

    scatter (i in range(length(PrintReads.tcra))) {
        call Immunelens {
        input:
            tcra = PrintReads.tcra[i],
            igh = PrintReads.igh[i],
            immune_lens_rscript = immune_lens_rscript
        }  
    }
    

    output {
        Array[File] cell_fractions = Immunelens.cell_fractions
    }

    meta {
        author: 'Hannah Poisner'
        email: 'hannah.m.poisner@vanderbilt.edu'
        description: "This workflow runs the ImmuneLens Program for TCRA and IGH"
    }
}


task PrintReads {
    input {
        Array[File] cram_files
        Array[File] crai_files
        File fasta
        File fasta_index
        File dict
        Int cpu = 4
        Int preemptible = 2
    }

    parameter_meta {
        cram_files: {
            description: "cram files",
            localization_optional: true
        }
        crai_files: {
            description: "crai files",
            localization_optional: true
        }
    }

    command <<<
    mkdir -p outputs
    cram_files_array=(~{sep=' ' cram_files})
    crai_files_array=(~{sep=' ' crai_files})

    for ((i=0; i<${#cram_files_array[@]}; i++)); do
        cram_file="${cram_files_array[i]}"
        crai_file="${crai_files_array[i]}"
        sample_id=$(basename "$cram_file" .cram)

        # Optional: Link index if needed
        ln -sf "$crai_file" "${cram_file}.crai"

        /gatk/gatk PrintReads \
            -I "$cram_file" \
            -L chr14:21621904-22752132 \
            -L chr14:105566277-106879844 \
            -R ~{fasta} \
            -O outputs/${sample_id}_cf.cram \
            --cloud-prefetch-buffer 0 \
            --cloud-index-prefetch-buffer 0

        samtools index outputs/${sample_id}_cf.cram outputs/${sample_id}_cf.cram.crai

        samtools depth -a --reference ~{fasta} -r chr14:21621904-22752132 outputs/${sample_id}_cf.cram > outputs/${sample_id}_tcra.txt
        samtools depth -a --reference ~{fasta} -r chr14:105566277-106879844 outputs/${sample_id}_cf.cram > outputs/${sample_id}_igh.txt
    done

    >>>

    output {
        Array[File] tcra = glob("outputs/*_tcra.txt")
        Array[File] igh = glob("outputs/*_igh.txt")
    }

    runtime {
    docker: 'gcr.io/bick-aps2/briansha/pileup_preprocessing:latest'
    memory: "20 GB"
    disks: "local-disk 500 HDD"
    preemptible: preemptible
    }    


}


task Immunelens {
    input {
        File tcra
        File igh
        File immune_lens_rscript
        Int cpu = 1
        Int preemptible = 2
    }
    String prefix = basename(tcra, "_tcra.txt")

    command <<<
        R < ~{immune_lens_rscript} --vanilla --args ~{tcra} ~{igh} ~{prefix}_il.txt   
        >>>
    output {
        File cell_fractions = "~{prefix}_il.txt"
    }
    runtime {
    docker: "briansha/immunelens:1.0.2"
    memory: "20 GB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
    }
}