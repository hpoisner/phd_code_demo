#A Run file for submitting DNANexus Regenie runs.

# User specified input files
pheno_file="project-J14qfbjJ21Qq0q2X50kBgXzb:file-J39G3x8J21Qk3pX8XV08BkxQ"
covar_file="${pheno_file}"
mask_file="project-J14qfbjJ21Qq0q2X50kBgXzb:file-J38jXxQJ21QV410jF8P3bfQQ"
pred_file="project-J14qfbjJ21Qq0q2X50kBgXzb:file-J3B3X80JYkJfJJy8x1gzFPyf"
loco_file="project-J14qfbjJ21Qq0q2X50kBgXzb:file-J3B3X80JYkJg1YPfqZ47JPf1"



# User specified input strings
phenocol="tp53"
covarcol="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
catcov="sex"
outfile_name="rvas_tp53_prev"

# User specified runtime
# 32 CPU job at 28 min
# 16 CPU job at 33 min
# 8 CPU at 47 min
# 4 CPU at 57 min
# 2 CPU at 1 hour 56 min
# Go for 16 CPU jobs at low priority. 
# Trying to minimize the chance we lose a spot and still maintain reasonable execution time efficiency.
docker_image="briansha/regenie:v4.0"
priority="normal"
instance_type="mem3_ssd1_v2_x16"
output_folder="pershy1/rare_variant_tp53_prev"
tag="wgs_tp53_prev"

# User specified TSV file
# Debugging: Make sure this file is not "ASCII text, with CRLF line terminators"
# - "file chip_to_ccus.txt"
# - If it is, change it. The presence of CRLF line terminators indicates that the file (chip_to_ccus.txt) uses Windows-style line endings, 
# which are incompatible with most Unix/Linux shell scripts that expect LF line endings. 
# These terminators can cause issues when processing the file, as they add unexpected \r characters that disrupt variable assignments and parsing.
# To fix: "sed -i '' 's/\r$//' chip_to_ccus.txt"
# - Check: "file chip_to_ccus.txt"
# - chip_to_ccus.txt: ASCII text
tsv_file="rare_files_other_chroms_rcv.txt"

# Dynamically set input file names
pheno_file_name=$(dx ls "${pheno_file}")
covar_file_name=$(dx ls "${covar_file}")
mask_file_name=$(dx ls "${mask_file}")
pred_file_name=$(dx ls "${pred_file}")
loco_file_name=$(dx ls "${loco_file}")

# Ensure the file ends with a newline
#echo >> "${tsv_file}"

# Skip the header line and iterate over the rows
tail -n +2 "${tsv_file}" | {
    while IFS=$'\t' read -r chr igenotype_bgens igenotype_bgis igenotype_samples istep2_anno_file istep2_gene_list; do
    
        # Dynamically set variables
        step2_bgen_file="${igenotype_bgens}"
        step2_bgen_file_index="${igenotype_bgis}"
        step2_sample_file="${igenotype_samples}"
        step2_anno_file="${istep2_anno_file}"
        step2_gene_list="${istep2_gene_list}"
        OUTPUT_FOLDER=${output_folder}

        echo ${step_2_gene_list}
        # Dynamically set variable names
        step2_bgen_file_name=$(dx ls "${step2_bgen_file}")
        step2_bgen_file_index_name=$(dx ls "${step2_bgen_file_index}")
        step2_sample_file_name=$(dx ls "${step2_sample_file}")
        step2_anno_file_name=$(dx ls "${step2_anno_file}")
        step2_gene_list_name=$(dx ls "${step2_gene_list}")

        # Dynamically set variable basenames
        step2_bgen_file_name_basename=$(basename "${step2_bgen_file_name}" | sed 's/\.[^.]*$//')

        # Dynamically set input strings
        step2_outfile="${outfile_name}_${step2_bgen_file_name}"

        # Echo all names for debugging
        echo "Processing Chromosome: ${chr}"
        echo "Input Files:"
        echo "  - Phenotype File: ${pheno_file_name}"
        echo "  - Covariate File: ${covar_file_name}"
        echo "  - Pred File: ${pred_file_name}"
        echo " - Loco File: ${loco_file_name}"
        echo "Chromosome-Specific Files:"
        echo "  - Sample File: ${step2_sample_file_name}"
        echo "  - BGEN File: ${step2_bgen_file_name}"
        echo "  - BGEN Index File: ${step2_bgen_file_index_name}"
        echo "  - Anno File (Step 2): ${step2_anno_file_name}"
        echo "  - Gene List (Step 2): ${step2_gene_list_name}"
        echo "Outputs:"
        echo "  - Step 2 Output File: ${step2_outfile}"
        echo "Destination: ${OUTPUT_FOLDER}"

        # Run dx command for each chromosome
        dx run swiss-army-knife \
            --priority "${priority}" \
            --instance-type "${instance_type}" \
            -iin="${pheno_file}" \
            -iin="${step2_anno_file}" \
            -iin="${step2_gene_list}" \
            -iin="${mask_file}" \
            -iin="${step2_sample_file}" \
            -iin="${step2_bgen_file}" \
            -iin="${pred_file}" \
            -iin="${loco_file}" \
            -icmd="
            # Step 2
            regenie \
                --step 2 \
                --bgen ${step2_bgen_file_name} \
                --phenoFile ${pheno_file_name} \
                --bsize 200 \
                --pred ${pred_file_name} \
                --sample ${step2_sample_file_name} \
                --anno-file ${step2_anno_file_name} \
                --set-list ${step2_gene_list_name} \
                --mask-def ${mask_file_name} \
                --covarFile ${covar_file_name} \
                --phenoColList ${phenocol} \
                --covarColList ${covarcol} \
                --catCovarList ${catcov} \
                --bt \
                --vc-tests skato \
                --out ${step2_outfile}
            " \
            -imount_inputs=true \
            -iimage="${docker_image}" \
            --destination "${OUTPUT_FOLDER}" \
            --tag "${tag}" \
            --tag "batch_${chr}" \
            -y \
            --brief
    done
}