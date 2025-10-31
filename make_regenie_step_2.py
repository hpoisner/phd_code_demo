import subprocess
import argparse

parser = argparse.ArgumentParser(description = 'Arguments for generating commands file')

parser.add_argument('--output_folder', '-of', type=str, required=True, help='Path to output folder, must use double quotes')
parser.add_argument('--output_file', '-ofi', type=str, required=True, help='Command File name, must use double quotes')

args = parser.parse_args()

OUTPUT_FOLDER = args.output_folder
output_file = args.output_file

with open("new_step_1_step2_batch_files.txt") as file, open(output_file, "w") as outfile:
    next(file)  # Skip header line
    for line in file:
        fields = line.strip().split('\t')
        chr, igenotype_pgens, igenotype_psams, igenotype_pvars = fields[:4]
        if chr!="X": 
            # Remove newline character
            pgen_proc=subprocess.Popen(f'dx ls "{igenotype_pgens}"', shell=True, stdout=subprocess.PIPE)
            psam_proc=subprocess.Popen(f'dx ls "{igenotype_psams}"', shell=True, stdout=subprocess.PIPE)
            pvar_proc=subprocess.Popen(f'dx ls "{igenotype_pvars}"', shell=True, stdout=subprocess.PIPE)
            pgen_name = pgen_proc.stdout.read().decode('ascii').strip('\n')
            psam_name = psam_proc.stdout.read().decode('ascii').strip('\n')
            pvar_name = pvar_proc.stdout.read().decode('ascii').strip('\n')
            # Construct dx run command
            cmd = f'''dx run app-swiss-army-knife \
            --priority normal \
            --instance-type "mem1_ssd1_v2_x36" \
            --tag="chr{chr}_make_gwas_step2" \
            -iin="{igenotype_pgens}" \
            -iin="{igenotype_psams}" \
            -iin="{igenotype_pvars}" \
            -icmd="plink2 --pgen {pgen_name} --psam {psam_name} --pvar {pvar_name} --no-pheno --chr {chr} --maf 0.01 --hwe 1e-15 --geno 0.1 --write-snplist --out chr{chr}_ukb_step2_maf_0_01_geno_0_1_hwe1e15" \
            --destination "{OUTPUT_FOLDER}" \
            -y --brief \
            --ignore-reuse\n'''
        
            # Write command to output file
            outfile.write(cmd)
        else:
            # Remove newline character                                                                                                                                                                                                    
            pgen_proc=subprocess.Popen(f'dx ls "{igenotype_pgens}"', shell=True, stdout=subprocess.PIPE)
            psam_proc=subprocess.Popen(f'dx ls "{igenotype_psams}"', shell=True, stdout=subprocess.PIPE)
            pvar_proc=subprocess.Popen(f'dx ls "{igenotype_pvars}"', shell=True, stdout=subprocess.PIPE)
            pgen_name = pgen_proc.stdout.read().decode('ascii').strip('\n')
            psam_name = psam_proc.stdout.read().decode('ascii').strip('\n')
            pvar_name = pvar_proc.stdout.read().decode('ascii').strip('\n')
            # Construct dx run command
            cmd = f'''dx run app-swiss-army-knife \
            --priority normal \
            --instance-type "mem1_ssd1_v2_x36" \
            --tag="chr{chr}_make_gwas_step1" \
            -iin="{igenotype_pgens}" \
            -iin="{igenotype_psams}" \
            -iin="{igenotype_pvars}" \
            -icmd="plink2 --pgen {pgen_name} --psam {psam_name} --pvar {pvar_name} --no-pheno --chr {chr} --maf 0.01 --hwe 1e-15 --geno 0.1 --split-par hg38 --write-snplist --out chr{chr}_ukb_step2_maf_0_01_geno_0_1_hwe1e15" \
            --destination "{OUTPUT_FOLDER}" \
            -y --brief \
            --ignore-reuse\n'''

            # Write command to output file
            outfile.write(cmd)

# Execute the commands from the file using subprocess if needed
# subprocess.run(f"bash {output_file}", shell=True)

# Execute outside of this script on the CLI
# sh commands.txt
