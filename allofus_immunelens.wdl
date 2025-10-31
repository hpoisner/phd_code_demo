version 1.0
# Example workflow

workflow immunelens {
    input {
        Array[File] tcra_dp
        Array[File] igh_dp
        File immune_lens_rscript
    }
    

    scatter (tcra_igh in zip(tcra_dp, igh_dp)) {
        call Immunelens {
            input:
            tcra = tcra_igh.left,
            igh = tcra_igh.right,
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

task Immunelens{
    input {
        File tcra
        File igh
        File immune_lens_rscript
        Int cpu = 1
        Int preemptible = 2
    }
    String prefix = basename(tcra, ".txt")

    command {
        R < ${immune_lens_rscript} --vanilla --args ${tcra} ${igh} ${prefix}_il_v2.txt   
        }
    output {
        File cell_fractions = "${prefix}_il_v2.txt"
    }
    runtime {
    docker: "briansha/immunelens:1.0.2"
    memory: "8 GB"
    disks: "local-disk 50 HDD"
    bootDiskSizeGb: "50"
    preemptible: preemptible
    }
}
