#!/bin/bash

#SBATCH --job-name=1900_busco_check_genome_completeness
#SBATCH --partition=cpu
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=wanghy@dicp.ac.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# module purge
# module load miniconda3
# source activate /lustre/home/acct-clslhz/clslhz/miniconda3/envs/busco
cd /lustre/home/acct-clslhz/clslhz/why_ssGEMs/result/check_assembly_completeness


(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAA_6.re.fna.out -m genome -e 1e-3
mv AAA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAB_6.re.fna.out -m genome -e 1e-3
mv AAB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAB_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAC_6.re.fna.out -m genome -e 1e-3
mv AAC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAD_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAD_6.re.fna.out -m genome -e 1e-3
mv AAD_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAD_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAD_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAE_6.re.fna.out -m genome -e 1e-3
mv AAE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAG_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAG_6.re.fna.out -m genome -e 1e-3
mv AAG_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAG_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAG_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAH_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAH_6.re.fna.out -m genome -e 1e-3
mv AAH_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAH_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAH_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAI_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAI_6.re.fna.out -m genome -e 1e-3
mv AAI_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAI_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAI_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAK_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAK_6.re.fna.out -m genome -e 1e-3
mv AAK_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAK_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAK_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAL_3.re.fna.out -m genome -e 1e-3
mv AAL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAL_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAM_3.re.fna.out -m genome -e 1e-3
mv AAM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAN_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAN_3.re.fna.out -m genome -e 1e-3
mv AAN_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAN_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAN_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAP_3.re.fna.out -m genome -e 1e-3
mv AAP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAQ_3.re.fna.out -m genome -e 1e-3
mv AAQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAQ_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAR_3.re.fna.out -m genome -e 1e-3
mv AAR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAS_3.re.fna.out -m genome -e 1e-3
mv AAS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAT_3.re.fna.out -m genome -e 1e-3
mv AAT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AAV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AAV_3.re.fna.out -m genome -e 1e-3
mv AAV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AAV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AAV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABA_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABA_3.re.fna.out -m genome -e 1e-3
mv ABA_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABA_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABA_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABB_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABB_5.re.fna.out -m genome -e 1e-3
mv ABB_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABB_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABB_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABC_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABC_5.re.fna.out -m genome -e 1e-3
mv ABC_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABC_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABC_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABD_5.re.fna.out -m genome -e 1e-3
mv ABD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABE_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABE_5.re.fna.out -m genome -e 1e-3
mv ABE_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABE_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABE_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABF_6.re.fna.out -m genome -e 1e-3
mv ABF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABG_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABG_5.re.fna.out -m genome -e 1e-3
mv ABG_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABG_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABG_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABH_5.re.fna.out -m genome -e 1e-3
mv ABH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABI_5.re.fna.out -m genome -e 1e-3
mv ABI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABI_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABK_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABK_5.re.fna.out -m genome -e 1e-3
mv ABK_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABK_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABK_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABL_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABL_5.re.fna.out -m genome -e 1e-3
mv ABL_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABL_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABL_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABM_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABM_5.re.fna.out -m genome -e 1e-3
mv ABM_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABM_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABM_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABP_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABP_6.re.fna.out -m genome -e 1e-3
mv ABP_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABP_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABP_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABQ_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABQ_6.re.fna.out -m genome -e 1e-3
mv ABQ_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABQ_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABQ_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABR_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABR_6.re.fna.out -m genome -e 1e-3
mv ABR_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABR_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABR_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABS_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABS_6.re.fna.out -m genome -e 1e-3
mv ABS_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABS_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABS_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABT_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABT_6.re.fna.out -m genome -e 1e-3
mv ABT_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABT_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABT_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ABV_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ABV_6.re.fna.out -m genome -e 1e-3
mv ABV_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ABV_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ABV_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACA_6.re.fna.out -m genome -e 1e-3
mv ACA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACB_6.re.fna.out -m genome -e 1e-3
mv ACB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACB_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACC_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACC_7.re.fna.out -m genome -e 1e-3
mv ACC_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACC_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACC_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACD_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACD_7.re.fna.out -m genome -e 1e-3
mv ACD_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACD_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACD_7.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACF_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACF_7.re.fna.out -m genome -e 1e-3
mv ACF_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACF_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACF_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACG_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACG_7.re.fna.out -m genome -e 1e-3
mv ACG_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACG_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACG_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACH_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACH_7.re.fna.out -m genome -e 1e-3
mv ACH_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACH_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACH_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACI_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACI_7.re.fna.out -m genome -e 1e-3
mv ACI_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACI_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACI_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACK_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACK_7.re.fna.out -m genome -e 1e-3
mv ACK_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACK_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACK_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACL_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACL_7.re.fna.out -m genome -e 1e-3
mv ACL_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACL_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACL_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACM_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACM_7.re.fna.out -m genome -e 1e-3
mv ACM_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACM_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACM_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACN_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACN_8.re.fna.out -m genome -e 1e-3
mv ACN_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACN_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACN_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACP_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACP_8.re.fna.out -m genome -e 1e-3
mv ACP_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACP_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACP_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACQ_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACQ_8.re.fna.out -m genome -e 1e-3
mv ACQ_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACQ_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACQ_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACR_2.re.fna.out -m genome -e 1e-3
mv ACR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACR_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACS_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACS_2.re.fna.out -m genome -e 1e-3
mv ACS_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACS_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACS_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACT_2.re.fna.out -m genome -e 1e-3
mv ACT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ACV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ACV_2.re.fna.out -m genome -e 1e-3
mv ACV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ACV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ACV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADA_2.re.fna.out -m genome -e 1e-3
mv ADA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADB_2.re.fna.out -m genome -e 1e-3
mv ADB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADC_2.re.fna.out -m genome -e 1e-3
mv ADC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADC_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADD_2.re.fna.out -m genome -e 1e-3
mv ADD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADE_2.re.fna.out -m genome -e 1e-3
mv ADE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADF_4.re.fna.out -m genome -e 1e-3
mv ADF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADF_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADG_4.re.fna.out -m genome -e 1e-3
mv ADG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADH_4.re.fna.out -m genome -e 1e-3
mv ADH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADI_4.re.fna.out -m genome -e 1e-3
mv ADI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADK_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADK_8.re.fna.out -m genome -e 1e-3
mv ADK_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADK_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADK_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADL_4.re.fna.out -m genome -e 1e-3
mv ADL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADM_2.re.fna.out -m genome -e 1e-3
mv ADM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADN_4.re.fna.out -m genome -e 1e-3
mv ADN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADP_3.re.fna.out -m genome -e 1e-3
mv ADP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADQ_4.re.fna.out -m genome -e 1e-3
mv ADQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADR_4.re.fna.out -m genome -e 1e-3
mv ADR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADR_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADS_4.re.fna.out -m genome -e 1e-3
mv ADS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADT_4.re.fna.out -m genome -e 1e-3
mv ADT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ADV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ADV_4.re.fna.out -m genome -e 1e-3
mv ADV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ADV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ADV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEA_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEA_8.re.fna.out -m genome -e 1e-3
mv AEA_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEA_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEA_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEB_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEB_8.re.fna.out -m genome -e 1e-3
mv AEB_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEB_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEB_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEC_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEC_8.re.fna.out -m genome -e 1e-3
mv AEC_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEC_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEC_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEE_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEE_3.re.fna.out -m genome -e 1e-3
mv AEE_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEE_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEE_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEF_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEF_8.re.fna.out -m genome -e 1e-3
mv AEF_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEF_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEF_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEG_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEG_8.re.fna.out -m genome -e 1e-3
mv AEG_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEG_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEG_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEH_3.re.fna.out -m genome -e 1e-3
mv AEH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEH_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEI_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEI_3.re.fna.out -m genome -e 1e-3
mv AEI_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEI_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEI_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEK_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEK_3.re.fna.out -m genome -e 1e-3
mv AEK_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEK_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEK_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEL_3.re.fna.out -m genome -e 1e-3
mv AEL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEM_3.re.fna.out -m genome -e 1e-3
mv AEM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEN_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEN_3.re.fna.out -m genome -e 1e-3
mv AEN_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEN_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEN_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEP_3.re.fna.out -m genome -e 1e-3
mv AEP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEQ_3.re.fna.out -m genome -e 1e-3
mv AEQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEQ_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AER_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AER_4.re.fna.out -m genome -e 1e-3
mv AER_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AER_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AER_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AES_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AES_4.re.fna.out -m genome -e 1e-3
mv AES_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AES_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AES_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AET_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AET_5.re.fna.out -m genome -e 1e-3
mv AET_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AET_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AET_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AEV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AEV_5.re.fna.out -m genome -e 1e-3
mv AEV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AEV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AEV_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFA_4.re.fna.out -m genome -e 1e-3
mv AFA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFB_4.re.fna.out -m genome -e 1e-3
mv AFB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFC_4.re.fna.out -m genome -e 1e-3
mv AFC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFD_4.re.fna.out -m genome -e 1e-3
mv AFD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFF_4.re.fna.out -m genome -e 1e-3
mv AFF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFG_4.re.fna.out -m genome -e 1e-3
mv AFG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFH_4.re.fna.out -m genome -e 1e-3
mv AFH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFI_5.re.fna.out -m genome -e 1e-3
mv AFI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFI_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFK_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFK_5.re.fna.out -m genome -e 1e-3
mv AFK_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFK_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFK_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFL_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFL_5.re.fna.out -m genome -e 1e-3
mv AFL_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFL_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFL_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFM_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFM_5.re.fna.out -m genome -e 1e-3
mv AFM_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFM_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFM_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFN_4.re.fna.out -m genome -e 1e-3
mv AFN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFP_1.re.fna.out -m genome -e 1e-3
mv AFP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFP_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFQ_1.re.fna.out -m genome -e 1e-3
mv AFQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFR_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFR_5.re.fna.out -m genome -e 1e-3
mv AFR_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFR_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFR_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFS_5.re.fna.out -m genome -e 1e-3
mv AFS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFT_5.re.fna.out -m genome -e 1e-3
mv AFT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AFV_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AFV_1.re.fna.out -m genome -e 1e-3
mv AFV_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AFV_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AFV_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGA_5.re.fna.out -m genome -e 1e-3
mv AGA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGA_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGB_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGB_1.re.fna.out -m genome -e 1e-3
mv AGB_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGB_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGB_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGC_1.re.fna.out -m genome -e 1e-3
mv AGC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGE_1.re.fna.out -m genome -e 1e-3
mv AGE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGE_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGF_1.re.fna.out -m genome -e 1e-3
mv AGF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGF_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGG_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGG_1.re.fna.out -m genome -e 1e-3
mv AGG_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGG_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGG_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGH_1.re.fna.out -m genome -e 1e-3
mv AGH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGI_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGI_1.re.fna.out -m genome -e 1e-3
mv AGI_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGI_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGI_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGK_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGK_1.re.fna.out -m genome -e 1e-3
mv AGK_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGK_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGK_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGL_2.re.fna.out -m genome -e 1e-3
mv AGL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGM_1.re.fna.out -m genome -e 1e-3
mv AGM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGM_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGN_2.re.fna.out -m genome -e 1e-3
mv AGN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGN_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGP_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGP_2.re.fna.out -m genome -e 1e-3
mv AGP_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGP_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGP_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGR_2.re.fna.out -m genome -e 1e-3
mv AGR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGR_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGS_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGS_2.re.fna.out -m genome -e 1e-3
mv AGS_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGS_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGS_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGT_2.re.fna.out -m genome -e 1e-3
mv AGT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AGV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AGV_2.re.fna.out -m genome -e 1e-3
mv AGV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AGV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AGV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHA_2.re.fna.out -m genome -e 1e-3
mv AHA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHB_2.re.fna.out -m genome -e 1e-3
mv AHB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHC_2.re.fna.out -m genome -e 1e-3
mv AHC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHC_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHD_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHD_1.re.fna.out -m genome -e 1e-3
mv AHD_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHD_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHD_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHE_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHE_3.re.fna.out -m genome -e 1e-3
mv AHE_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHE_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHE_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHF_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHF_3.re.fna.out -m genome -e 1e-3
mv AHF_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHF_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHF_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHG_3.re.fna.out -m genome -e 1e-3
mv AHG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHH_3.re.fna.out -m genome -e 1e-3
mv AHH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHI_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHI_3.re.fna.out -m genome -e 1e-3
mv AHI_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHI_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHI_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHK_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHK_3.re.fna.out -m genome -e 1e-3
mv AHK_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHK_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHK_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHL_3.re.fna.out -m genome -e 1e-3
mv AHL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHM_1.re.fna.out -m genome -e 1e-3
mv AHM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHM_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHN_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHN_3.re.fna.out -m genome -e 1e-3
mv AHN_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHN_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHN_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHP_1.re.fna.out -m genome -e 1e-3
mv AHP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHP_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHQ_1.re.fna.out -m genome -e 1e-3
mv AHQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHR_1.re.fna.out -m genome -e 1e-3
mv AHR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHS_3.re.fna.out -m genome -e 1e-3
mv AHS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHT_3.re.fna.out -m genome -e 1e-3
mv AHT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AHV_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AHV_1.re.fna.out -m genome -e 1e-3
mv AHV_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AHV_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AHV_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIA_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIA_1.re.fna.out -m genome -e 1e-3
mv AIA_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIA_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIA_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIB_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIB_1.re.fna.out -m genome -e 1e-3
mv AIB_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIB_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIB_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIC_1.re.fna.out -m genome -e 1e-3
mv AIC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AID_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AID_2.re.fna.out -m genome -e 1e-3
mv AID_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AID_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AID_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIE_1.re.fna.out -m genome -e 1e-3
mv AIE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIE_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIF_1.re.fna.out -m genome -e 1e-3
mv AIF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIF_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIG_2.re.fna.out -m genome -e 1e-3
mv AIG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AII_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AII_1.re.fna.out -m genome -e 1e-3
mv AII_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AII_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AII_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIK_2.re.fna.out -m genome -e 1e-3
mv AIK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIL_3.re.fna.out -m genome -e 1e-3
mv AIL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIM_2.re.fna.out -m genome -e 1e-3
mv AIM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIN_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIN_1.re.fna.out -m genome -e 1e-3
mv AIN_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIN_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIN_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIP_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIP_2.re.fna.out -m genome -e 1e-3
mv AIP_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIP_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIP_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIQ_1.re.fna.out -m genome -e 1e-3
mv AIQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIR_1.re.fna.out -m genome -e 1e-3
mv AIR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIR_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AIT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AIT_2.re.fna.out -m genome -e 1e-3
mv AIT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AIT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AIT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKB_2.re.fna.out -m genome -e 1e-3
mv AKB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKC_4.re.fna.out -m genome -e 1e-3
mv AKC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKD_2.re.fna.out -m genome -e 1e-3
mv AKD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKE_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKE_3.re.fna.out -m genome -e 1e-3
mv AKE_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKE_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKE_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKG_3.re.fna.out -m genome -e 1e-3
mv AKG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKH_3.re.fna.out -m genome -e 1e-3
mv AKH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKI_4.re.fna.out -m genome -e 1e-3
mv AKI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKL_3.re.fna.out -m genome -e 1e-3
mv AKL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKM_4.re.fna.out -m genome -e 1e-3
mv AKM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKM_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKN_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKN_3.re.fna.out -m genome -e 1e-3
mv AKN_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKN_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKN_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKP_3.re.fna.out -m genome -e 1e-3
mv AKP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKQ_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKQ_5.re.fna.out -m genome -e 1e-3
mv AKQ_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKQ_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKQ_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKR_3.re.fna.out -m genome -e 1e-3
mv AKR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKS_5.re.fna.out -m genome -e 1e-3
mv AKS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKT_4.re.fna.out -m genome -e 1e-3
mv AKT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AKV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AKV_3.re.fna.out -m genome -e 1e-3
mv AKV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AKV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AKV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALA_5.re.fna.out -m genome -e 1e-3
mv ALA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALA_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALB_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALB_3.re.fna.out -m genome -e 1e-3
mv ALB_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALB_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALB_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALC_4.re.fna.out -m genome -e 1e-3
mv ALC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALC_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALD_4.re.fna.out -m genome -e 1e-3
mv ALD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALE_4.re.fna.out -m genome -e 1e-3
mv ALE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALG_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALG_5.re.fna.out -m genome -e 1e-3
mv ALG_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALG_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALG_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALH_4.re.fna.out -m genome -e 1e-3
mv ALH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALI_5.re.fna.out -m genome -e 1e-3
mv ALI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALI_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALK_4.re.fna.out -m genome -e 1e-3
mv ALK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALL_4.re.fna.out -m genome -e 1e-3
mv ALL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALM_4.re.fna.out -m genome -e 1e-3
mv ALM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALN_4.re.fna.out -m genome -e 1e-3
mv ALN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALP_4.re.fna.out -m genome -e 1e-3
mv ALP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALP_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALQ_4.re.fna.out -m genome -e 1e-3
mv ALQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALR_4.re.fna.out -m genome -e 1e-3
mv ALR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALS_5.re.fna.out -m genome -e 1e-3
mv ALS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALT_5.re.fna.out -m genome -e 1e-3
mv ALT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ALV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ALV_2.re.fna.out -m genome -e 1e-3
mv ALV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ALV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ALV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMA_5.re.fna.out -m genome -e 1e-3
mv AMA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMA_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMB_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMB_5.re.fna.out -m genome -e 1e-3
mv AMB_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMB_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMB_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMC_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMC_5.re.fna.out -m genome -e 1e-3
mv AMC_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMC_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMC_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMD_5.re.fna.out -m genome -e 1e-3
mv AMD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMF_4.re.fna.out -m genome -e 1e-3
mv AMF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMF_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMG_4.re.fna.out -m genome -e 1e-3
mv AMG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMH_5.re.fna.out -m genome -e 1e-3
mv AMH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMI_4.re.fna.out -m genome -e 1e-3
mv AMI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMK_4.re.fna.out -m genome -e 1e-3
mv AMK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AML_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AML_4.re.fna.out -m genome -e 1e-3
mv AML_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AML_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AML_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMM_4.re.fna.out -m genome -e 1e-3
mv AMM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMN_5.re.fna.out -m genome -e 1e-3
mv AMN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMP_5.re.fna.out -m genome -e 1e-3
mv AMP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMP_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMQ_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMQ_5.re.fna.out -m genome -e 1e-3
mv AMQ_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMQ_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMQ_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMR_4.re.fna.out -m genome -e 1e-3
mv AMR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMR_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMS_5.re.fna.out -m genome -e 1e-3
mv AMS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AMT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AMT_5.re.fna.out -m genome -e 1e-3
mv AMT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AMT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AMT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANA_5.re.fna.out -m genome -e 1e-3
mv ANA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANA_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANC_4.re.fna.out -m genome -e 1e-3
mv ANC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AND_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AND_4.re.fna.out -m genome -e 1e-3
mv AND_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AND_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AND_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANE_4.re.fna.out -m genome -e 1e-3
mv ANE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANF_4.re.fna.out -m genome -e 1e-3
mv ANF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANG_4.re.fna.out -m genome -e 1e-3
mv ANG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANH_4.re.fna.out -m genome -e 1e-3
mv ANH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANI_4.re.fna.out -m genome -e 1e-3
mv ANI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANI_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANK_4.re.fna.out -m genome -e 1e-3
mv ANK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANL_2.re.fna.out -m genome -e 1e-3
mv ANL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANM_2.re.fna.out -m genome -e 1e-3
mv ANM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANN_4.re.fna.out -m genome -e 1e-3
mv ANN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANP_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANP_2.re.fna.out -m genome -e 1e-3
mv ANP_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANP_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANP_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANQ_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANQ_2.re.fna.out -m genome -e 1e-3
mv ANQ_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANQ_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANQ_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANR_2.re.fna.out -m genome -e 1e-3
mv ANR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANR_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANS_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANS_2.re.fna.out -m genome -e 1e-3
mv ANS_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANS_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANS_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANT_2.re.fna.out -m genome -e 1e-3
mv ANT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ANV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ANV_5.re.fna.out -m genome -e 1e-3
mv ANV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ANV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ANV_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APA_2.re.fna.out -m genome -e 1e-3
mv APA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.APA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APC_4.re.fna.out -m genome -e 1e-3
mv APC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APD_4.re.fna.out -m genome -e 1e-3
mv APD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APE_4.re.fna.out -m genome -e 1e-3
mv APE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APF_4.re.fna.out -m genome -e 1e-3
mv APF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APG_4.re.fna.out -m genome -e 1e-3
mv APG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APH_4.re.fna.out -m genome -e 1e-3
mv APH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/API_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o API_4.re.fna.out -m genome -e 1e-3
mv API_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.API_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r API_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APK_4.re.fna.out -m genome -e 1e-3
mv APK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APL_4.re.fna.out -m genome -e 1e-3
mv APL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APL_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APM_4.re.fna.out -m genome -e 1e-3
mv APM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APN_4.re.fna.out -m genome -e 1e-3
mv APN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APP_4.re.fna.out -m genome -e 1e-3
mv APP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APQ_4.re.fna.out -m genome -e 1e-3
mv APQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APR_4.re.fna.out -m genome -e 1e-3
mv APR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APS_4.re.fna.out -m genome -e 1e-3
mv APS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APT_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APT_7.re.fna.out -m genome -e 1e-3
mv APT_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.APT_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APT_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/APV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o APV_4.re.fna.out -m genome -e 1e-3
mv APV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.APV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r APV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQA_4.re.fna.out -m genome -e 1e-3
mv AQA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQB_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQB_7.re.fna.out -m genome -e 1e-3
mv AQB_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQB_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQB_7.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQC_4.re.fna.out -m genome -e 1e-3
mv AQC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQD_4.re.fna.out -m genome -e 1e-3
mv AQD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQE_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQE_5.re.fna.out -m genome -e 1e-3
mv AQE_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQE_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQE_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQG_4.re.fna.out -m genome -e 1e-3
mv AQG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQH_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQH_2.re.fna.out -m genome -e 1e-3
mv AQH_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQH_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQH_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQI_4.re.fna.out -m genome -e 1e-3
mv AQI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQM_4.re.fna.out -m genome -e 1e-3
mv AQM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQN_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQN_7.re.fna.out -m genome -e 1e-3
mv AQN_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQN_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQN_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQP_4.re.fna.out -m genome -e 1e-3
mv AQP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQQ_4.re.fna.out -m genome -e 1e-3
mv AQQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQQ_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQR_4.re.fna.out -m genome -e 1e-3
mv AQR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQS_4.re.fna.out -m genome -e 1e-3
mv AQS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AQT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AQT_4.re.fna.out -m genome -e 1e-3
mv AQT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.AQT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AQT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARA_4.re.fna.out -m genome -e 1e-3
mv ARA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARB_2.re.fna.out -m genome -e 1e-3
mv ARB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARC_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARC_5.re.fna.out -m genome -e 1e-3
mv ARC_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARC_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARC_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARE_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARE_5.re.fna.out -m genome -e 1e-3
mv ARE_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARE_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARE_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARG_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARG_5.re.fna.out -m genome -e 1e-3
mv ARG_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARG_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARG_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARH_5.re.fna.out -m genome -e 1e-3
mv ARH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARI_5.re.fna.out -m genome -e 1e-3
mv ARI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARI_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARK_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARK_5.re.fna.out -m genome -e 1e-3
mv ARK_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARK_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARK_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARL_2.re.fna.out -m genome -e 1e-3
mv ARL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARM_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARM_6.re.fna.out -m genome -e 1e-3
mv ARM_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARM_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARM_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARN_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARN_6.re.fna.out -m genome -e 1e-3
mv ARN_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARN_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARN_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARP_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARP_2.re.fna.out -m genome -e 1e-3
mv ARP_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARP_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARP_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARQ_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARQ_6.re.fna.out -m genome -e 1e-3
mv ARQ_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARQ_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARQ_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARR_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARR_6.re.fna.out -m genome -e 1e-3
mv ARR_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARR_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARR_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARS_1.re.fna.out -m genome -e 1e-3
mv ARS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ART_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ART_6.re.fna.out -m genome -e 1e-3
mv ART_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ART_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ART_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ARV_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ARV_6.re.fna.out -m genome -e 1e-3
mv ARV_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ARV_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ARV_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASA_6.re.fna.out -m genome -e 1e-3
mv ASA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASB_2.re.fna.out -m genome -e 1e-3
mv ASB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASC_6.re.fna.out -m genome -e 1e-3
mv ASC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASD_2.re.fna.out -m genome -e 1e-3
mv ASD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASE_2.re.fna.out -m genome -e 1e-3
mv ASE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASG_3.re.fna.out -m genome -e 1e-3
mv ASG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASH_3.re.fna.out -m genome -e 1e-3
mv ASH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASL_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASL_8.re.fna.out -m genome -e 1e-3
mv ASL_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASL_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASL_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASM_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASM_8.re.fna.out -m genome -e 1e-3
mv ASM_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASM_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASM_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASN_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASN_8.re.fna.out -m genome -e 1e-3
mv ASN_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASN_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASN_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASP_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASP_8.re.fna.out -m genome -e 1e-3
mv ASP_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASP_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASP_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASR_1.re.fna.out -m genome -e 1e-3
mv ASR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASS_1.re.fna.out -m genome -e 1e-3
mv ASS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AST_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AST_1.re.fna.out -m genome -e 1e-3
mv AST_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.AST_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AST_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ASV_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ASV_1.re.fna.out -m genome -e 1e-3
mv ASV_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ASV_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ASV_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATA_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATA_8.re.fna.out -m genome -e 1e-3
mv ATA_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATA_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATA_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATB_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATB_1.re.fna.out -m genome -e 1e-3
mv ATB_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATB_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATB_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATC_1.re.fna.out -m genome -e 1e-3
mv ATC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATD_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATD_1.re.fna.out -m genome -e 1e-3
mv ATD_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATD_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATD_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATE_1.re.fna.out -m genome -e 1e-3
mv ATE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATE_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATH_3.re.fna.out -m genome -e 1e-3
mv ATH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATM_4.re.fna.out -m genome -e 1e-3
mv ATM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATN_4.re.fna.out -m genome -e 1e-3
mv ATN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATP_4.re.fna.out -m genome -e 1e-3
mv ATP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATQ_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATQ_7.re.fna.out -m genome -e 1e-3
mv ATQ_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATQ_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATQ_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATR_4.re.fna.out -m genome -e 1e-3
mv ATR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATS_3.re.fna.out -m genome -e 1e-3
mv ATS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/ATV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o ATV_5.re.fna.out -m genome -e 1e-3
mv ATV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.ATV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r ATV_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVB_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVB_8.re.fna.out -m genome -e 1e-3
mv AVB_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVB_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVB_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVC_3.re.fna.out -m genome -e 1e-3
mv AVC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVC_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVD_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVD_3.re.fna.out -m genome -e 1e-3
mv AVD_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVD_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVD_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVH_5.re.fna.out -m genome -e 1e-3
mv AVH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVI_5.re.fna.out -m genome -e 1e-3
mv AVI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVI_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVK_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVK_6.re.fna.out -m genome -e 1e-3
mv AVK_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVK_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVK_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVL_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVL_6.re.fna.out -m genome -e 1e-3
mv AVL_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVL_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVL_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVM_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVM_5.re.fna.out -m genome -e 1e-3
mv AVM_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVM_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVM_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVN_5.re.fna.out -m genome -e 1e-3
mv AVN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVP_5.re.fna.out -m genome -e 1e-3
mv AVP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVP_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVQ_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVQ_5.re.fna.out -m genome -e 1e-3
mv AVQ_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVQ_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVQ_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVR_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVR_6.re.fna.out -m genome -e 1e-3
mv AVR_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVR_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVR_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVS_5.re.fna.out -m genome -e 1e-3
mv AVS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVT_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVT_6.re.fna.out -m genome -e 1e-3
mv AVT_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVT_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVT_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/AVV_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o AVV_8.re.fna.out -m genome -e 1e-3
mv AVV_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.AVV_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r AVV_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAA_6.re.fna.out -m genome -e 1e-3
mv BAA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAB_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAB_8.re.fna.out -m genome -e 1e-3
mv BAB_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAB_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAB_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAC_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAC_8.re.fna.out -m genome -e 1e-3
mv BAC_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAC_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAC_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAD_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAD_8.re.fna.out -m genome -e 1e-3
mv BAD_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAD_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAD_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAE_6.re.fna.out -m genome -e 1e-3
mv BAE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAF_6.re.fna.out -m genome -e 1e-3
mv BAF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAG_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAG_6.re.fna.out -m genome -e 1e-3
mv BAG_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAG_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAG_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAH_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAH_6.re.fna.out -m genome -e 1e-3
mv BAH_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAH_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAH_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAI_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAI_7.re.fna.out -m genome -e 1e-3
mv BAI_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAI_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAI_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAK_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAK_7.re.fna.out -m genome -e 1e-3
mv BAK_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAK_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAK_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAL_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAL_7.re.fna.out -m genome -e 1e-3
mv BAL_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAL_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAL_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAM_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAM_7.re.fna.out -m genome -e 1e-3
mv BAM_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAM_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAM_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAN_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAN_7.re.fna.out -m genome -e 1e-3
mv BAN_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAN_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAN_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAP_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAP_7.re.fna.out -m genome -e 1e-3
mv BAP_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAP_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAP_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAQ_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAQ_7.re.fna.out -m genome -e 1e-3
mv BAQ_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAQ_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAQ_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAR_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAR_7.re.fna.out -m genome -e 1e-3
mv BAR_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAR_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAR_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAS_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAS_8.re.fna.out -m genome -e 1e-3
mv BAS_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAS_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAS_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAT_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAT_8.re.fna.out -m genome -e 1e-3
mv BAT_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAT_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAT_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BAV_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BAV_7.re.fna.out -m genome -e 1e-3
mv BAV_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BAV_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BAV_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBA_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBA_8.re.fna.out -m genome -e 1e-3
mv BBA_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBA_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBA_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBB_6.re.fna.out -m genome -e 1e-3
mv BBB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBB_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBC_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBC_8.re.fna.out -m genome -e 1e-3
mv BBC_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBC_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBC_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBD_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBD_8.re.fna.out -m genome -e 1e-3
mv BBD_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBD_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBD_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBE_6.re.fna.out -m genome -e 1e-3
mv BBE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBF_6.re.fna.out -m genome -e 1e-3
mv BBF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBG_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBG_6.re.fna.out -m genome -e 1e-3
mv BBG_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBG_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBG_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBH_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBH_6.re.fna.out -m genome -e 1e-3
mv BBH_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBH_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBH_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBI_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBI_6.re.fna.out -m genome -e 1e-3
mv BBI_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBI_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBI_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBK_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBK_6.re.fna.out -m genome -e 1e-3
mv BBK_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBK_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBK_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBL_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBL_6.re.fna.out -m genome -e 1e-3
mv BBL_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBL_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBL_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBM_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBM_6.re.fna.out -m genome -e 1e-3
mv BBM_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBM_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBM_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBP_3.re.fna.out -m genome -e 1e-3
mv BBP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBQ_3.re.fna.out -m genome -e 1e-3
mv BBQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBQ_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBR_3.re.fna.out -m genome -e 1e-3
mv BBR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBS_3.re.fna.out -m genome -e 1e-3
mv BBS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBT_3.re.fna.out -m genome -e 1e-3
mv BBT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BBV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BBV_3.re.fna.out -m genome -e 1e-3
mv BBV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BBV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BBV_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCA_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCA_3.re.fna.out -m genome -e 1e-3
mv BCA_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCA_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCA_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCB_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCB_3.re.fna.out -m genome -e 1e-3
mv BCB_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCB_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCB_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCC_3.re.fna.out -m genome -e 1e-3
mv BCC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCC_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCD_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCD_8.re.fna.out -m genome -e 1e-3
mv BCD_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCD_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCD_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCE_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCE_8.re.fna.out -m genome -e 1e-3
mv BCE_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCE_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCE_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCF_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCF_8.re.fna.out -m genome -e 1e-3
mv BCF_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCF_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCF_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCG_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCG_8.re.fna.out -m genome -e 1e-3
mv BCG_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCG_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCG_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCH_5.re.fna.out -m genome -e 1e-3
mv BCH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCI_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCI_8.re.fna.out -m genome -e 1e-3
mv BCI_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCI_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCI_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCK_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCK_8.re.fna.out -m genome -e 1e-3
mv BCK_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCK_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCK_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCL_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCL_5.re.fna.out -m genome -e 1e-3
mv BCL_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCL_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCL_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCM_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCM_8.re.fna.out -m genome -e 1e-3
mv BCM_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCM_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCM_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCN_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCN_8.re.fna.out -m genome -e 1e-3
mv BCN_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCN_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCN_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCP_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCP_8.re.fna.out -m genome -e 1e-3
mv BCP_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCP_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCP_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCQ_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCQ_8.re.fna.out -m genome -e 1e-3
mv BCQ_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCQ_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCQ_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCR_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCR_8.re.fna.out -m genome -e 1e-3
mv BCR_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCR_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCR_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCS_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCS_8.re.fna.out -m genome -e 1e-3
mv BCS_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCS_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCS_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCT_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCT_8.re.fna.out -m genome -e 1e-3
mv BCT_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCT_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCT_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BCV_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BCV_8.re.fna.out -m genome -e 1e-3
mv BCV_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BCV_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BCV_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDA_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDA_8.re.fna.out -m genome -e 1e-3
mv BDA_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDA_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDA_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDB_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDB_8.re.fna.out -m genome -e 1e-3
mv BDB_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDB_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDB_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDC_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDC_8.re.fna.out -m genome -e 1e-3
mv BDC_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDC_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDC_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDD_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDD_8.re.fna.out -m genome -e 1e-3
mv BDD_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDD_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDD_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDE_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDE_8.re.fna.out -m genome -e 1e-3
mv BDE_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDE_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDE_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDF_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDF_8.re.fna.out -m genome -e 1e-3
mv BDF_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDF_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDF_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDG_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDG_5.re.fna.out -m genome -e 1e-3
mv BDG_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDG_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDG_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDH_5.re.fna.out -m genome -e 1e-3
mv BDH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDI_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDI_5.re.fna.out -m genome -e 1e-3
mv BDI_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDI_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDI_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDK_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDK_5.re.fna.out -m genome -e 1e-3
mv BDK_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDK_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDK_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDL_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDL_5.re.fna.out -m genome -e 1e-3
mv BDL_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDL_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDL_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDM_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDM_5.re.fna.out -m genome -e 1e-3
mv BDM_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDM_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDM_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDN_5.re.fna.out -m genome -e 1e-3
mv BDN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDP_5.re.fna.out -m genome -e 1e-3
mv BDP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDP_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDQ_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDQ_6.re.fna.out -m genome -e 1e-3
mv BDQ_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDQ_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDQ_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDR_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDR_6.re.fna.out -m genome -e 1e-3
mv BDR_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDR_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDR_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDS_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDS_6.re.fna.out -m genome -e 1e-3
mv BDS_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDS_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDS_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDT_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDT_6.re.fna.out -m genome -e 1e-3
mv BDT_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDT_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDT_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BDV_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BDV_6.re.fna.out -m genome -e 1e-3
mv BDV_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BDV_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BDV_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEA_6.re.fna.out -m genome -e 1e-3
mv BEA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEB_6.re.fna.out -m genome -e 1e-3
mv BEB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEB_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEC_6.re.fna.out -m genome -e 1e-3
mv BEC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BED_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BED_6.re.fna.out -m genome -e 1e-3
mv BED_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BED_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BED_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEE_6.re.fna.out -m genome -e 1e-3
mv BEE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEF_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEF_7.re.fna.out -m genome -e 1e-3
mv BEF_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEF_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEF_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEG_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEG_7.re.fna.out -m genome -e 1e-3
mv BEG_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEG_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEG_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEH_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEH_7.re.fna.out -m genome -e 1e-3
mv BEH_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEH_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEH_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEI_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEI_7.re.fna.out -m genome -e 1e-3
mv BEI_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEI_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEI_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEK_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEK_7.re.fna.out -m genome -e 1e-3
mv BEK_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEK_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEK_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEM_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEM_7.re.fna.out -m genome -e 1e-3
mv BEM_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEM_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEM_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEN_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEN_7.re.fna.out -m genome -e 1e-3
mv BEN_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEN_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEN_7.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEP_4.re.fna.out -m genome -e 1e-3
mv BEP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEQ_4.re.fna.out -m genome -e 1e-3
mv BEQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BER_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BER_2.re.fna.out -m genome -e 1e-3
mv BER_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BER_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BER_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BES_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BES_2.re.fna.out -m genome -e 1e-3
mv BES_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BES_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BES_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BET_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BET_1.re.fna.out -m genome -e 1e-3
mv BET_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BET_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BET_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BEV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BEV_3.re.fna.out -m genome -e 1e-3
mv BEV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BEV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BEV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFA_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFA_3.re.fna.out -m genome -e 1e-3
mv BFA_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFA_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFA_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFB_4.re.fna.out -m genome -e 1e-3
mv BFB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFC_3.re.fna.out -m genome -e 1e-3
mv BFC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFC_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFD_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFD_1.re.fna.out -m genome -e 1e-3
mv BFD_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFD_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFD_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFE_4.re.fna.out -m genome -e 1e-3
mv BFE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFG_4.re.fna.out -m genome -e 1e-3
mv BFG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFH_3.re.fna.out -m genome -e 1e-3
mv BFH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFK_2.re.fna.out -m genome -e 1e-3
mv BFK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFL_4.re.fna.out -m genome -e 1e-3
mv BFL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFM_3.re.fna.out -m genome -e 1e-3
mv BFM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFN_4.re.fna.out -m genome -e 1e-3
mv BFN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFP_3.re.fna.out -m genome -e 1e-3
mv BFP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFQ_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFQ_8.re.fna.out -m genome -e 1e-3
mv BFQ_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFQ_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFQ_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFR_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFR_8.re.fna.out -m genome -e 1e-3
mv BFR_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFR_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFR_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFS_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFS_8.re.fna.out -m genome -e 1e-3
mv BFS_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFS_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFS_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFT_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFT_8.re.fna.out -m genome -e 1e-3
mv BFT_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFT_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFT_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BFV_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BFV_8.re.fna.out -m genome -e 1e-3
mv BFV_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BFV_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BFV_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGA_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGA_7.re.fna.out -m genome -e 1e-3
mv BGA_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGA_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGA_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGB_7.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGB_7.re.fna.out -m genome -e 1e-3
mv BGB_7.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGB_7.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGB_7.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGC_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGC_8.re.fna.out -m genome -e 1e-3
mv BGC_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGC_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGC_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGD_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGD_8.re.fna.out -m genome -e 1e-3
mv BGD_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGD_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGD_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGE_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGE_8.re.fna.out -m genome -e 1e-3
mv BGE_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGE_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGE_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGF_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGF_8.re.fna.out -m genome -e 1e-3
mv BGF_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGF_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGF_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGG_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGG_8.re.fna.out -m genome -e 1e-3
mv BGG_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGG_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGG_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGH_1.re.fna.out -m genome -e 1e-3
mv BGH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGI_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGI_1.re.fna.out -m genome -e 1e-3
mv BGI_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGI_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGI_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGK_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGK_1.re.fna.out -m genome -e 1e-3
mv BGK_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGK_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGK_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGL_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGL_1.re.fna.out -m genome -e 1e-3
mv BGL_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGL_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGL_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGM_1.re.fna.out -m genome -e 1e-3
mv BGM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGM_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGN_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGN_1.re.fna.out -m genome -e 1e-3
mv BGN_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGN_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGN_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGP_1.re.fna.out -m genome -e 1e-3
mv BGP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGP_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGQ_1.re.fna.out -m genome -e 1e-3
mv BGQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGR_1.re.fna.out -m genome -e 1e-3
mv BGR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGS_1.re.fna.out -m genome -e 1e-3
mv BGS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGS_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGT_2.re.fna.out -m genome -e 1e-3
mv BGT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BGV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BGV_2.re.fna.out -m genome -e 1e-3
mv BGV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BGV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BGV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHA_2.re.fna.out -m genome -e 1e-3
mv BHA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHB_2.re.fna.out -m genome -e 1e-3
mv BHB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHC_2.re.fna.out -m genome -e 1e-3
mv BHC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHC_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHD_2.re.fna.out -m genome -e 1e-3
mv BHD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHE_2.re.fna.out -m genome -e 1e-3
mv BHE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHF_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHF_2.re.fna.out -m genome -e 1e-3
mv BHF_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHF_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHF_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHG_2.re.fna.out -m genome -e 1e-3
mv BHG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHH_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHH_2.re.fna.out -m genome -e 1e-3
mv BHH_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHH_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHH_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHI_4.re.fna.out -m genome -e 1e-3
mv BHI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHK_4.re.fna.out -m genome -e 1e-3
mv BHK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHL_4.re.fna.out -m genome -e 1e-3
mv BHL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHM_4.re.fna.out -m genome -e 1e-3
mv BHM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHN_4.re.fna.out -m genome -e 1e-3
mv BHN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHP_4.re.fna.out -m genome -e 1e-3
mv BHP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHQ_4.re.fna.out -m genome -e 1e-3
mv BHQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHR_4.re.fna.out -m genome -e 1e-3
mv BHR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHS_1.re.fna.out -m genome -e 1e-3
mv BHS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHT_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHT_1.re.fna.out -m genome -e 1e-3
mv BHT_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHT_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHT_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BHV_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BHV_1.re.fna.out -m genome -e 1e-3
mv BHV_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BHV_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BHV_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIA_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIA_1.re.fna.out -m genome -e 1e-3
mv BIA_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIA_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIA_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIB_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIB_1.re.fna.out -m genome -e 1e-3
mv BIB_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIB_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIB_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIC_1.re.fna.out -m genome -e 1e-3
mv BIC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BID_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BID_1.re.fna.out -m genome -e 1e-3
mv BID_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BID_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BID_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIE_1.re.fna.out -m genome -e 1e-3
mv BIE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIE_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIF_1.re.fna.out -m genome -e 1e-3
mv BIF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIF_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIG_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIG_1.re.fna.out -m genome -e 1e-3
mv BIG_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIG_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIG_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIH_4.re.fna.out -m genome -e 1e-3
mv BIH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BII_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BII_4.re.fna.out -m genome -e 1e-3
mv BII_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BII_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BII_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIK_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIK_3.re.fna.out -m genome -e 1e-3
mv BIK_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIK_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIK_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIL_3.re.fna.out -m genome -e 1e-3
mv BIL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIM_3.re.fna.out -m genome -e 1e-3
mv BIM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIN_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIN_3.re.fna.out -m genome -e 1e-3
mv BIN_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIN_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIN_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIP_3.re.fna.out -m genome -e 1e-3
mv BIP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIQ_3.re.fna.out -m genome -e 1e-3
mv BIQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIQ_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIR_2.re.fna.out -m genome -e 1e-3
mv BIR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIR_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIS_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIS_2.re.fna.out -m genome -e 1e-3
mv BIS_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIS_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIS_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BIT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BIT_2.re.fna.out -m genome -e 1e-3
mv BIT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BIT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BIT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKA_2.re.fna.out -m genome -e 1e-3
mv BKA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKA_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKC_1.re.fna.out -m genome -e 1e-3
mv BKC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKD_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKD_1.re.fna.out -m genome -e 1e-3
mv BKD_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKD_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKD_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKE_1.re.fna.out -m genome -e 1e-3
mv BKE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKE_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKF_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKF_3.re.fna.out -m genome -e 1e-3
mv BKF_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKF_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKF_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKG_3.re.fna.out -m genome -e 1e-3
mv BKG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKH_1.re.fna.out -m genome -e 1e-3
mv BKH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKI_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKI_2.re.fna.out -m genome -e 1e-3
mv BKI_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKI_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKI_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKK_2.re.fna.out -m genome -e 1e-3
mv BKK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKL_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKL_1.re.fna.out -m genome -e 1e-3
mv BKL_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKL_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKL_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKM_1.re.fna.out -m genome -e 1e-3
mv BKM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKM_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKN_2.re.fna.out -m genome -e 1e-3
mv BKN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKN_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKP_1.re.fna.out -m genome -e 1e-3
mv BKP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKP_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKQ_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKQ_2.re.fna.out -m genome -e 1e-3
mv BKQ_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKQ_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKQ_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKR_1.re.fna.out -m genome -e 1e-3
mv BKR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKS_1.re.fna.out -m genome -e 1e-3
mv BKS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BKT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BKT_2.re.fna.out -m genome -e 1e-3
mv BKT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BKT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BKT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLA_4.re.fna.out -m genome -e 1e-3
mv BLA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLB_4.re.fna.out -m genome -e 1e-3
mv BLB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLC_4.re.fna.out -m genome -e 1e-3
mv BLC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLD_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLD_3.re.fna.out -m genome -e 1e-3
mv BLD_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLD_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLD_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLE_2.re.fna.out -m genome -e 1e-3
mv BLE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLF_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLF_2.re.fna.out -m genome -e 1e-3
mv BLF_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLF_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLF_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLG_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLG_1.re.fna.out -m genome -e 1e-3
mv BLG_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLG_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLG_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLH_1.re.fna.out -m genome -e 1e-3
mv BLH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLI_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLI_3.re.fna.out -m genome -e 1e-3
mv BLI_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLI_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLI_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLK_2.re.fna.out -m genome -e 1e-3
mv BLK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLL_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLL_1.re.fna.out -m genome -e 1e-3
mv BLL_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLL_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLL_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLM_3.re.fna.out -m genome -e 1e-3
mv BLM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLN_2.re.fna.out -m genome -e 1e-3
mv BLN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLN_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLP_5.re.fna.out -m genome -e 1e-3
mv BLP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLP_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLQ_1.re.fna.out -m genome -e 1e-3
mv BLQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLR_3.re.fna.out -m genome -e 1e-3
mv BLR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLS_5.re.fna.out -m genome -e 1e-3
mv BLS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLT_2.re.fna.out -m genome -e 1e-3
mv BLT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BLV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BLV_3.re.fna.out -m genome -e 1e-3
mv BLV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BLV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BLV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMA_2.re.fna.out -m genome -e 1e-3
mv BMA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMB_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMB_1.re.fna.out -m genome -e 1e-3
mv BMB_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMB_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMB_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMC_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMC_1.re.fna.out -m genome -e 1e-3
mv BMC_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMC_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMC_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMD_5.re.fna.out -m genome -e 1e-3
mv BMD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BME_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BME_1.re.fna.out -m genome -e 1e-3
mv BME_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BME_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BME_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMF_1.re.fna.out -m genome -e 1e-3
mv BMF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMF_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMG_2.re.fna.out -m genome -e 1e-3
mv BMG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMH_5.re.fna.out -m genome -e 1e-3
mv BMH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMI_4.re.fna.out -m genome -e 1e-3
mv BMI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMK_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMK_5.re.fna.out -m genome -e 1e-3
mv BMK_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMK_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMK_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BML_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BML_5.re.fna.out -m genome -e 1e-3
mv BML_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BML_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BML_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMM_3.re.fna.out -m genome -e 1e-3
mv BMM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMN_4.re.fna.out -m genome -e 1e-3
mv BMN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMQ_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMQ_5.re.fna.out -m genome -e 1e-3
mv BMQ_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMQ_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMQ_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMR_2.re.fna.out -m genome -e 1e-3
mv BMR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMR_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BMV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BMV_2.re.fna.out -m genome -e 1e-3
mv BMV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BMV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BMV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNA_4.re.fna.out -m genome -e 1e-3
mv BNA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNB_4.re.fna.out -m genome -e 1e-3
mv BNB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNC_4.re.fna.out -m genome -e 1e-3
mv BNC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BND_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BND_3.re.fna.out -m genome -e 1e-3
mv BND_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BND_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BND_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNE_2.re.fna.out -m genome -e 1e-3
mv BNE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNF_4.re.fna.out -m genome -e 1e-3
mv BNF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNG_4.re.fna.out -m genome -e 1e-3
mv BNG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNH_3.re.fna.out -m genome -e 1e-3
mv BNH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNI_4.re.fna.out -m genome -e 1e-3
mv BNI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNI_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNK_4.re.fna.out -m genome -e 1e-3
mv BNK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNL_4.re.fna.out -m genome -e 1e-3
mv BNL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNM_4.re.fna.out -m genome -e 1e-3
mv BNM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNN_2.re.fna.out -m genome -e 1e-3
mv BNN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNN_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNP_4.re.fna.out -m genome -e 1e-3
mv BNP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNQ_4.re.fna.out -m genome -e 1e-3
mv BNQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNR_4.re.fna.out -m genome -e 1e-3
mv BNR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNT_4.re.fna.out -m genome -e 1e-3
mv BNT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BNV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BNV_2.re.fna.out -m genome -e 1e-3
mv BNV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BNV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BNV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPA_4.re.fna.out -m genome -e 1e-3
mv BPA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPA_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPB_4.re.fna.out -m genome -e 1e-3
mv BPB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPC_4.re.fna.out -m genome -e 1e-3
mv BPC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPF_4.re.fna.out -m genome -e 1e-3
mv BPF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPG_4.re.fna.out -m genome -e 1e-3
mv BPG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPH_4.re.fna.out -m genome -e 1e-3
mv BPH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPI_4.re.fna.out -m genome -e 1e-3
mv BPI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPK_4.re.fna.out -m genome -e 1e-3
mv BPK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPL_2.re.fna.out -m genome -e 1e-3
mv BPL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPM_4.re.fna.out -m genome -e 1e-3
mv BPM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPN_2.re.fna.out -m genome -e 1e-3
mv BPN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPN_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPP_4.re.fna.out -m genome -e 1e-3
mv BPP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPQ_4.re.fna.out -m genome -e 1e-3
mv BPQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPR_4.re.fna.out -m genome -e 1e-3
mv BPR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPT_5.re.fna.out -m genome -e 1e-3
mv BPT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BPV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BPV_4.re.fna.out -m genome -e 1e-3
mv BPV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BPV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BPV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQA_4.re.fna.out -m genome -e 1e-3
mv BQA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQB_4.re.fna.out -m genome -e 1e-3
mv BQB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQC_4.re.fna.out -m genome -e 1e-3
mv BQC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQE_4.re.fna.out -m genome -e 1e-3
mv BQE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQG_4.re.fna.out -m genome -e 1e-3
mv BQG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQG_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQH_4.re.fna.out -m genome -e 1e-3
mv BQH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQI_4.re.fna.out -m genome -e 1e-3
mv BQI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQM_4.re.fna.out -m genome -e 1e-3
mv BQM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQN_4.re.fna.out -m genome -e 1e-3
mv BQN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQP_4.re.fna.out -m genome -e 1e-3
mv BQP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQQ_4.re.fna.out -m genome -e 1e-3
mv BQQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQR_4.re.fna.out -m genome -e 1e-3
mv BQR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQS_4.re.fna.out -m genome -e 1e-3
mv BQS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BQT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BQT_4.re.fna.out -m genome -e 1e-3
mv BQT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BQT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BQT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRA_4.re.fna.out -m genome -e 1e-3
mv BRA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRA_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRB_4.re.fna.out -m genome -e 1e-3
mv BRB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRC_4.re.fna.out -m genome -e 1e-3
mv BRC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRD_4.re.fna.out -m genome -e 1e-3
mv BRD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRE_4.re.fna.out -m genome -e 1e-3
mv BRE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRF_4.re.fna.out -m genome -e 1e-3
mv BRF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRG_4.re.fna.out -m genome -e 1e-3
mv BRG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRH_4.re.fna.out -m genome -e 1e-3
mv BRH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRK_4.re.fna.out -m genome -e 1e-3
mv BRK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRL_4.re.fna.out -m genome -e 1e-3
mv BRL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRM_4.re.fna.out -m genome -e 1e-3
mv BRM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRM_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRP_4.re.fna.out -m genome -e 1e-3
mv BRP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRQ_4.re.fna.out -m genome -e 1e-3
mv BRQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BRV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BRV_4.re.fna.out -m genome -e 1e-3
mv BRV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BRV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BRV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSA_4.re.fna.out -m genome -e 1e-3
mv BSA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSB_4.re.fna.out -m genome -e 1e-3
mv BSB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSC_4.re.fna.out -m genome -e 1e-3
mv BSC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSD_4.re.fna.out -m genome -e 1e-3
mv BSD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSE_4.re.fna.out -m genome -e 1e-3
mv BSE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSG_4.re.fna.out -m genome -e 1e-3
mv BSG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSI_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSI_3.re.fna.out -m genome -e 1e-3
mv BSI_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSI_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSI_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSK_4.re.fna.out -m genome -e 1e-3
mv BSK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSL_4.re.fna.out -m genome -e 1e-3
mv BSL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSM_3.re.fna.out -m genome -e 1e-3
mv BSM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSN_5.re.fna.out -m genome -e 1e-3
mv BSN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSP_3.re.fna.out -m genome -e 1e-3
mv BSP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSR_4.re.fna.out -m genome -e 1e-3
mv BSR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BSS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BSS_4.re.fna.out -m genome -e 1e-3
mv BSS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BSS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BSS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BST_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BST_3.re.fna.out -m genome -e 1e-3
mv BST_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BST_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BST_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTA_5.re.fna.out -m genome -e 1e-3
mv BTA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTA_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTB_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTB_3.re.fna.out -m genome -e 1e-3
mv BTB_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTB_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTB_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTC_4.re.fna.out -m genome -e 1e-3
mv BTC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTD_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTD_1.re.fna.out -m genome -e 1e-3
mv BTD_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTD_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTD_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTE_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTE_3.re.fna.out -m genome -e 1e-3
mv BTE_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTE_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTE_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTF_4.re.fna.out -m genome -e 1e-3
mv BTF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTG_3.re.fna.out -m genome -e 1e-3
mv BTG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTH_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTH_3.re.fna.out -m genome -e 1e-3
mv BTH_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTH_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTH_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTI_4.re.fna.out -m genome -e 1e-3
mv BTI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTL_3.re.fna.out -m genome -e 1e-3
mv BTL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTN_4.re.fna.out -m genome -e 1e-3
mv BTN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTP_5.re.fna.out -m genome -e 1e-3
mv BTP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTP_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTR_4.re.fna.out -m genome -e 1e-3
mv BTR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTS_3.re.fna.out -m genome -e 1e-3
mv BTS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTT_3.re.fna.out -m genome -e 1e-3
mv BTT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BTV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BTV_5.re.fna.out -m genome -e 1e-3
mv BTV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BTV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BTV_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVA_4.re.fna.out -m genome -e 1e-3
mv BVA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVB_6.re.fna.out -m genome -e 1e-3
mv BVB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVB_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVC_6.re.fna.out -m genome -e 1e-3
mv BVC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVD_5.re.fna.out -m genome -e 1e-3
mv BVD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVE_4.re.fna.out -m genome -e 1e-3
mv BVE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVF_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVF_5.re.fna.out -m genome -e 1e-3
mv BVF_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVF_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVF_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVG_4.re.fna.out -m genome -e 1e-3
mv BVG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVH_1.re.fna.out -m genome -e 1e-3
mv BVH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVK_2.re.fna.out -m genome -e 1e-3
mv BVK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVM_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVM_3.re.fna.out -m genome -e 1e-3
mv BVM_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVM_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVM_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVT_3.re.fna.out -m genome -e 1e-3
mv BVT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/BVV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o BVV_4.re.fna.out -m genome -e 1e-3
mv BVV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.BVV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r BVV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAA_2.re.fna.out -m genome -e 1e-3
mv CAA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAB_4.re.fna.out -m genome -e 1e-3
mv CAB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAC_2.re.fna.out -m genome -e 1e-3
mv CAC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAC_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAD_4.re.fna.out -m genome -e 1e-3
mv CAD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAD_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAE_2.re.fna.out -m genome -e 1e-3
mv CAE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAF_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAF_2.re.fna.out -m genome -e 1e-3
mv CAF_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAF_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAF_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAH_4.re.fna.out -m genome -e 1e-3
mv CAH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAI_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAI_2.re.fna.out -m genome -e 1e-3
mv CAI_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAI_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAI_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAK_4.re.fna.out -m genome -e 1e-3
mv CAK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAL_3.re.fna.out -m genome -e 1e-3
mv CAL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAM_2.re.fna.out -m genome -e 1e-3
mv CAM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAN_4.re.fna.out -m genome -e 1e-3
mv CAN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAP_4.re.fna.out -m genome -e 1e-3
mv CAP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAQ_3.re.fna.out -m genome -e 1e-3
mv CAQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAQ_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAR_3.re.fna.out -m genome -e 1e-3
mv CAR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAS_1.re.fna.out -m genome -e 1e-3
mv CAS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAT_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAT_6.re.fna.out -m genome -e 1e-3
mv CAT_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAT_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAT_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CAV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CAV_4.re.fna.out -m genome -e 1e-3
mv CAV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CAV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CAV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBA_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBA_3.re.fna.out -m genome -e 1e-3
mv CBA_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBA_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBA_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBB_4.re.fna.out -m genome -e 1e-3
mv CBB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBC_6.re.fna.out -m genome -e 1e-3
mv CBC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBD_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBD_6.re.fna.out -m genome -e 1e-3
mv CBD_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBD_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBD_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBE_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBE_5.re.fna.out -m genome -e 1e-3
mv CBE_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBE_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBE_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBF_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBF_5.re.fna.out -m genome -e 1e-3
mv CBF_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBF_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBF_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBG_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBG_3.re.fna.out -m genome -e 1e-3
mv CBG_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBG_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBG_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBH_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBH_1.re.fna.out -m genome -e 1e-3
mv CBH_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBH_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBH_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBI_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBI_2.re.fna.out -m genome -e 1e-3
mv CBI_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBI_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBI_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBK_4.re.fna.out -m genome -e 1e-3
mv CBK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBL_2.re.fna.out -m genome -e 1e-3
mv CBL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBM_4.re.fna.out -m genome -e 1e-3
mv CBM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBN_5.re.fna.out -m genome -e 1e-3
mv CBN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CBR_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CBR_2.re.fna.out -m genome -e 1e-3
mv CBR_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CBR_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CBR_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCC_3.re.fna.out -m genome -e 1e-3
mv CCC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCC_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCD_4.re.fna.out -m genome -e 1e-3
mv CCD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCD_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCE_6.re.fna.out -m genome -e 1e-3
mv CCE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCF_6.re.fna.out -m genome -e 1e-3
mv CCF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCG_4.re.fna.out -m genome -e 1e-3
mv CCG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCH_4.re.fna.out -m genome -e 1e-3
mv CCH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCI_4.re.fna.out -m genome -e 1e-3
mv CCI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCL_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCL_1.re.fna.out -m genome -e 1e-3
mv CCL_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCL_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCL_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCM_1.re.fna.out -m genome -e 1e-3
mv CCM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCM_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCN_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCN_1.re.fna.out -m genome -e 1e-3
mv CCN_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCN_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCN_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCP_1.re.fna.out -m genome -e 1e-3
mv CCP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCP_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCQ_1.re.fna.out -m genome -e 1e-3
mv CCQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCQ_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCR_1.re.fna.out -m genome -e 1e-3
mv CCR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCS_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCS_1.re.fna.out -m genome -e 1e-3
mv CCS_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCS_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCS_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCT_2.re.fna.out -m genome -e 1e-3
mv CCT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CCV_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CCV_1.re.fna.out -m genome -e 1e-3
mv CCV_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CCV_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CCV_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDA_2.re.fna.out -m genome -e 1e-3
mv CDA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDB_2.re.fna.out -m genome -e 1e-3
mv CDB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDC_3.re.fna.out -m genome -e 1e-3
mv CDC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDC_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDD_2.re.fna.out -m genome -e 1e-3
mv CDD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDE_2.re.fna.out -m genome -e 1e-3
mv CDE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDF_1.re.fna.out -m genome -e 1e-3
mv CDF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDF_1.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDG_2.re.fna.out -m genome -e 1e-3
mv CDG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDH.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDH.re.fna.out -m genome -e 1e-3
mv CDH.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDH.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDH.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDI_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDI_1.re.fna.out -m genome -e 1e-3
mv CDI_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDI_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDI_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDK_4.re.fna.out -m genome -e 1e-3
mv CDK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDL_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDL_3.re.fna.out -m genome -e 1e-3
mv CDL_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDL_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDL_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDM_4.re.fna.out -m genome -e 1e-3
mv CDM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDN_4.re.fna.out -m genome -e 1e-3
mv CDN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDP_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDP_3.re.fna.out -m genome -e 1e-3
mv CDP_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDP_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDP_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDQ_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDQ_3.re.fna.out -m genome -e 1e-3
mv CDQ_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDQ_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDQ_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDR_3.re.fna.out -m genome -e 1e-3
mv CDR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDR_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDS_4.re.fna.out -m genome -e 1e-3
mv CDS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDT_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDT_2.re.fna.out -m genome -e 1e-3
mv CDT_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDT_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDT_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CDV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CDV_4.re.fna.out -m genome -e 1e-3
mv CDV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CDV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CDV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEA_4.re.fna.out -m genome -e 1e-3
mv CEA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEB_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEB_3.re.fna.out -m genome -e 1e-3
mv CEB_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEB_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEB_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEC_4.re.fna.out -m genome -e 1e-3
mv CEC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CED_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CED_4.re.fna.out -m genome -e 1e-3
mv CED_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CED_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CED_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEE_4.re.fna.out -m genome -e 1e-3
mv CEE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEF_4.re.fna.out -m genome -e 1e-3
mv CEF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEG_4.re.fna.out -m genome -e 1e-3
mv CEG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEG_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEH_4.re.fna.out -m genome -e 1e-3
mv CEH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEI_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEI_2.re.fna.out -m genome -e 1e-3
mv CEI_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEI_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEI_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEK_4.re.fna.out -m genome -e 1e-3
mv CEK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEL_4.re.fna.out -m genome -e 1e-3
mv CEL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEM_2.re.fna.out -m genome -e 1e-3
mv CEM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEN_4.re.fna.out -m genome -e 1e-3
mv CEN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEP_4.re.fna.out -m genome -e 1e-3
mv CEP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEQ_4.re.fna.out -m genome -e 1e-3
mv CEQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CER_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CER_4.re.fna.out -m genome -e 1e-3
mv CER_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CER_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CER_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CES_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CES_4.re.fna.out -m genome -e 1e-3
mv CES_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CES_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CES_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CET_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CET_4.re.fna.out -m genome -e 1e-3
mv CET_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CET_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CET_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CEV_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CEV_2.re.fna.out -m genome -e 1e-3
mv CEV_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CEV_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CEV_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFA_4.re.fna.out -m genome -e 1e-3
mv CFA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFB_4.re.fna.out -m genome -e 1e-3
mv CFB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFC_2.re.fna.out -m genome -e 1e-3
mv CFC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFC_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFD_4.re.fna.out -m genome -e 1e-3
mv CFD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFE_4.re.fna.out -m genome -e 1e-3
mv CFE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFF_4.re.fna.out -m genome -e 1e-3
mv CFF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFG_4.re.fna.out -m genome -e 1e-3
mv CFG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFH.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFH.re.fna.out -m genome -e 1e-3
mv CFH.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFH.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFH.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFI_4.re.fna.out -m genome -e 1e-3
mv CFI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFM_4.re.fna.out -m genome -e 1e-3
mv CFM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFN_4.re.fna.out -m genome -e 1e-3
mv CFN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFP_4.re.fna.out -m genome -e 1e-3
mv CFP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFQ_4.re.fna.out -m genome -e 1e-3
mv CFQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFR_3.re.fna.out -m genome -e 1e-3
mv CFR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFR_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFS_4.re.fna.out -m genome -e 1e-3
mv CFS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFT_5.re.fna.out -m genome -e 1e-3
mv CFT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CFV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CFV_4.re.fna.out -m genome -e 1e-3
mv CFV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CFV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CFV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGB_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGB_5.re.fna.out -m genome -e 1e-3
mv CGB_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGB_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGB_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGC_4.re.fna.out -m genome -e 1e-3
mv CGC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGD_5.re.fna.out -m genome -e 1e-3
mv CGD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGE_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGE_5.re.fna.out -m genome -e 1e-3
mv CGE_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGE_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGE_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGG_4.re.fna.out -m genome -e 1e-3
mv CGG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGH_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGH_5.re.fna.out -m genome -e 1e-3
mv CGH_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGH_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGH_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGI_4.re.fna.out -m genome -e 1e-3
mv CGI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGK_4.re.fna.out -m genome -e 1e-3
mv CGK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGL_4.re.fna.out -m genome -e 1e-3
mv CGL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGN_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGN_5.re.fna.out -m genome -e 1e-3
mv CGN_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGN_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGN_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGP_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGP_5.re.fna.out -m genome -e 1e-3
mv CGP_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGP_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGP_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGQ_4.re.fna.out -m genome -e 1e-3
mv CGQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGR_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGR_5.re.fna.out -m genome -e 1e-3
mv CGR_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGR_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGR_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGS_5.re.fna.out -m genome -e 1e-3
mv CGS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGS_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGT_5.re.fna.out -m genome -e 1e-3
mv CGT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CGV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CGV_3.re.fna.out -m genome -e 1e-3
mv CGV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CGV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CGV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHA_4.re.fna.out -m genome -e 1e-3
mv CHA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHB_4.re.fna.out -m genome -e 1e-3
mv CHB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHC_4.re.fna.out -m genome -e 1e-3
mv CHC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHD_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHD_3.re.fna.out -m genome -e 1e-3
mv CHD_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHD_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHD_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHE_4.re.fna.out -m genome -e 1e-3
mv CHE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHE_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHF_4.re.fna.out -m genome -e 1e-3
mv CHF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHG_4.re.fna.out -m genome -e 1e-3
mv CHG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHH_4.re.fna.out -m genome -e 1e-3
mv CHH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHI_4.re.fna.out -m genome -e 1e-3
mv CHI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHK_4.re.fna.out -m genome -e 1e-3
mv CHK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHL_4.re.fna.out -m genome -e 1e-3
mv CHL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHM_4.re.fna.out -m genome -e 1e-3
mv CHM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHN_4.re.fna.out -m genome -e 1e-3
mv CHN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHP_4.re.fna.out -m genome -e 1e-3
mv CHP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHQ_4.re.fna.out -m genome -e 1e-3
mv CHQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHQ_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHR_4.re.fna.out -m genome -e 1e-3
mv CHR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHS_4.re.fna.out -m genome -e 1e-3
mv CHS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHT_4.re.fna.out -m genome -e 1e-3
mv CHT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CHV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CHV_4.re.fna.out -m genome -e 1e-3
mv CHV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CHV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CHV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIA_4.re.fna.out -m genome -e 1e-3
mv CIA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIB_4.re.fna.out -m genome -e 1e-3
mv CIB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIC_4.re.fna.out -m genome -e 1e-3
mv CIC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CID_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CID_4.re.fna.out -m genome -e 1e-3
mv CID_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CID_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CID_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIE_4.re.fna.out -m genome -e 1e-3
mv CIE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIF_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIF_2.re.fna.out -m genome -e 1e-3
mv CIF_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIF_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIF_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIG_2.re.fna.out -m genome -e 1e-3
mv CIG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIH_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIH_2.re.fna.out -m genome -e 1e-3
mv CIH_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIH_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIH_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CII_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CII_2.re.fna.out -m genome -e 1e-3
mv CII_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CII_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CII_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIK_2.re.fna.out -m genome -e 1e-3
mv CIK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIL_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIL_2.re.fna.out -m genome -e 1e-3
mv CIL_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIL_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIL_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIM_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIM_2.re.fna.out -m genome -e 1e-3
mv CIM_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIM_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIM_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIN_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIN_2.re.fna.out -m genome -e 1e-3
mv CIN_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIN_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIN_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIP_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIP_2.re.fna.out -m genome -e 1e-3
mv CIP_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIP_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIP_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIQ_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIQ_2.re.fna.out -m genome -e 1e-3
mv CIQ_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIQ_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIQ_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIR_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIR_3.re.fna.out -m genome -e 1e-3
mv CIR_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIR_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIR_3.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIS_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIS_3.re.fna.out -m genome -e 1e-3
mv CIS_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIS_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIS_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIT_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIT_3.re.fna.out -m genome -e 1e-3
mv CIT_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIT_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIT_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CIV_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CIV_3.re.fna.out -m genome -e 1e-3
mv CIV_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CIV_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CIV_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKA_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKA_3.re.fna.out -m genome -e 1e-3
mv CKA_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKA_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKA_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKB_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKB_3.re.fna.out -m genome -e 1e-3
mv CKB_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKB_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKB_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKC_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKC_3.re.fna.out -m genome -e 1e-3
mv CKC_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKC_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKC_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKD_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKD_3.re.fna.out -m genome -e 1e-3
mv CKD_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKD_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKD_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKE_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKE_3.re.fna.out -m genome -e 1e-3
mv CKE_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKE_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKE_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKF_3.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKF_3.re.fna.out -m genome -e 1e-3
mv CKF_3.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKF_3.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKF_3.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKG_4.re.fna.out -m genome -e 1e-3
mv CKG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKG_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKH_4.re.fna.out -m genome -e 1e-3
mv CKH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKI_4.re.fna.out -m genome -e 1e-3
mv CKI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKK_4.re.fna.out -m genome -e 1e-3
mv CKK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKL_4.re.fna.out -m genome -e 1e-3
mv CKL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKM_4.re.fna.out -m genome -e 1e-3
mv CKM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKN_4.re.fna.out -m genome -e 1e-3
mv CKN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKP_4.re.fna.out -m genome -e 1e-3
mv CKP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKQ_4.re.fna.out -m genome -e 1e-3
mv CKQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKR_1.re.fna.out -m genome -e 1e-3
mv CKR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKS_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKS_5.re.fna.out -m genome -e 1e-3
mv CKS_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKS_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKS_5.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKT_5.re.fna.out -m genome -e 1e-3
mv CKT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CKV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CKV_5.re.fna.out -m genome -e 1e-3
mv CKV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CKV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CKV_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLA_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLA_5.re.fna.out -m genome -e 1e-3
mv CLA_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLA_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLA_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLB_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLB_5.re.fna.out -m genome -e 1e-3
mv CLB_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLB_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLB_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLC_4.re.fna.out -m genome -e 1e-3
mv CLC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLD_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLD_5.re.fna.out -m genome -e 1e-3
mv CLD_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLD_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLD_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLG_4.re.fna.out -m genome -e 1e-3
mv CLG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLG_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLH_4.re.fna.out -m genome -e 1e-3
mv CLH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLI_4.re.fna.out -m genome -e 1e-3
mv CLI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLK_4.re.fna.out -m genome -e 1e-3
mv CLK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLK_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLL_4.re.fna.out -m genome -e 1e-3
mv CLL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLM_4.re.fna.out -m genome -e 1e-3
mv CLM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLN_4.re.fna.out -m genome -e 1e-3
mv CLN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLP_4.re.fna.out -m genome -e 1e-3
mv CLP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLQ_4.re.fna.out -m genome -e 1e-3
mv CLQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLR_4.re.fna.out -m genome -e 1e-3
mv CLR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLS_4.re.fna.out -m genome -e 1e-3
mv CLS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLT_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLT_5.re.fna.out -m genome -e 1e-3
mv CLT_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLT_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLT_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CLV_5.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CLV_5.re.fna.out -m genome -e 1e-3
mv CLV_5.re.fna.out/short_summary.specific.saccharomycetes_odb10.CLV_5.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CLV_5.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMB_6.re.fna.out -m genome -e 1e-3
mv CMB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMB_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMC_6.re.fna.out -m genome -e 1e-3
mv CMC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMD_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMD_6.re.fna.out -m genome -e 1e-3
mv CMD_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMD_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMD_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CME_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CME_6.re.fna.out -m genome -e 1e-3
mv CME_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CME_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CME_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMF_6.re.fna.out -m genome -e 1e-3
mv CMF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMG_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMG_6.re.fna.out -m genome -e 1e-3
mv CMG_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMG_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMG_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMH_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMH_6.re.fna.out -m genome -e 1e-3
mv CMH_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMH_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMH_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMI_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMI_6.re.fna.out -m genome -e 1e-3
mv CMI_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMI_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMI_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMK_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMK_6.re.fna.out -m genome -e 1e-3
mv CMK_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMK_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMK_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CML_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CML_8.re.fna.out -m genome -e 1e-3
mv CML_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CML_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CML_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMM_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMM_8.re.fna.out -m genome -e 1e-3
mv CMM_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMM_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMM_8.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMN_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMN_8.re.fna.out -m genome -e 1e-3
mv CMN_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMN_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMN_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMP_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMP_8.re.fna.out -m genome -e 1e-3
mv CMP_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMP_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMP_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMQ_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMQ_8.re.fna.out -m genome -e 1e-3
mv CMQ_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMQ_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMQ_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMS_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMS_8.re.fna.out -m genome -e 1e-3
mv CMS_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMS_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMS_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CMT_8.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CMT_8.re.fna.out -m genome -e 1e-3
mv CMT_8.re.fna.out/short_summary.specific.saccharomycetes_odb10.CMT_8.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CMT_8.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNB_4.re.fna.out -m genome -e 1e-3
mv CNB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNC_4.re.fna.out -m genome -e 1e-3
mv CNC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNE_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNE_4.re.fna.out -m genome -e 1e-3
mv CNE_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNE_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNE_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNF_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNF_4.re.fna.out -m genome -e 1e-3
mv CNF_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNF_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNF_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNG_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNG_4.re.fna.out -m genome -e 1e-3
mv CNG_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNG_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNG_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNH_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNH_4.re.fna.out -m genome -e 1e-3
mv CNH_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNH_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNH_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNI_4.re.fna.out -m genome -e 1e-3
mv CNI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNI_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNK_4.re.fna.out -m genome -e 1e-3
mv CNK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNL_4.re.fna.out -m genome -e 1e-3
mv CNL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNM_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNM_1.re.fna.out -m genome -e 1e-3
mv CNM_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNM_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNM_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNN_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNN_1.re.fna.out -m genome -e 1e-3
mv CNN_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNN_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNN_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNP_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNP_1.re.fna.out -m genome -e 1e-3
mv CNP_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNP_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNP_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNQ_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNQ_1.re.fna.out -m genome -e 1e-3
mv CNQ_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNQ_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNQ_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNR_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNR_1.re.fna.out -m genome -e 1e-3
mv CNR_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNR_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNR_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNS_4.re.fna.out -m genome -e 1e-3
mv CNS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNS_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNT_4.re.fna.out -m genome -e 1e-3
mv CNT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CNV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CNV_4.re.fna.out -m genome -e 1e-3
mv CNV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CNV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CNV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPA_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPA_4.re.fna.out -m genome -e 1e-3
mv CPA_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPA_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPA_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPB_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPB_4.re.fna.out -m genome -e 1e-3
mv CPB_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPB_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPB_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPC_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPC_4.re.fna.out -m genome -e 1e-3
mv CPC_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPC_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPC_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPD_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPD_4.re.fna.out -m genome -e 1e-3
mv CPD_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPD_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPD_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPE_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPE_1.re.fna.out -m genome -e 1e-3
mv CPE_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPE_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPE_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPF_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPF_1.re.fna.out -m genome -e 1e-3
mv CPF_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPF_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPF_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPG_1.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPG_1.re.fna.out -m genome -e 1e-3
mv CPG_1.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPG_1.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPG_1.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPI_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPI_4.re.fna.out -m genome -e 1e-3
mv CPI_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPI_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPI_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPK_4.re.fna.out -m genome -e 1e-3
mv CPK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPL_4.re.fna.out -m genome -e 1e-3
mv CPL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPM_4.re.fna.out -m genome -e 1e-3
mv CPM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPM_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPN_4.re.fna.out -m genome -e 1e-3
mv CPN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPQ_4.re.fna.out -m genome -e 1e-3
mv CPQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPR_4.re.fna.out -m genome -e 1e-3
mv CPR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPS_4.re.fna.out -m genome -e 1e-3
mv CPS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CPT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CPT_4.re.fna.out -m genome -e 1e-3
mv CPT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CPT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CPT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQA_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQA_6.re.fna.out -m genome -e 1e-3
mv CQA_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQA_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQA_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQB_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQB_6.re.fna.out -m genome -e 1e-3
mv CQB_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQB_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQB_6.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQC_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQC_6.re.fna.out -m genome -e 1e-3
mv CQC_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQC_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQC_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQD_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQD_6.re.fna.out -m genome -e 1e-3
mv CQD_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQD_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQD_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQE_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQE_6.re.fna.out -m genome -e 1e-3
mv CQE_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQE_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQE_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQF_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQF_6.re.fna.out -m genome -e 1e-3
mv CQF_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQF_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQF_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQG_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQG_6.re.fna.out -m genome -e 1e-3
mv CQG_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQG_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQG_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQH_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQH_6.re.fna.out -m genome -e 1e-3
mv CQH_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQH_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQH_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQI_6.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQI_6.re.fna.out -m genome -e 1e-3
mv CQI_6.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQI_6.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQI_6.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQK_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQK_4.re.fna.out -m genome -e 1e-3
mv CQK_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQK_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQK_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQL_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQL_4.re.fna.out -m genome -e 1e-3
mv CQL_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQL_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQL_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQM_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQM_4.re.fna.out -m genome -e 1e-3
mv CQM_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQM_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQM_4.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQN_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQN_4.re.fna.out -m genome -e 1e-3
mv CQN_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQN_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQN_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQP_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQP_4.re.fna.out -m genome -e 1e-3
mv CQP_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQP_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQP_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQQ_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQQ_4.re.fna.out -m genome -e 1e-3
mv CQQ_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQQ_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQQ_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQR_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQR_4.re.fna.out -m genome -e 1e-3
mv CQR_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQR_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQR_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQS_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQS_4.re.fna.out -m genome -e 1e-3
mv CQS_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQS_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQS_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQT_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQT_4.re.fna.out -m genome -e 1e-3
mv CQT_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQT_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQT_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CQV_4.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CQV_4.re.fna.out -m genome -e 1e-3
mv CQV_4.re.fna.out/short_summary.specific.saccharomycetes_odb10.CQV_4.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CQV_4.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRA_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRA_2.re.fna.out -m genome -e 1e-3
mv CRA_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRA_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRA_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRB_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRB_2.re.fna.out -m genome -e 1e-3
mv CRB_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRB_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRB_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRC_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRC_2.re.fna.out -m genome -e 1e-3
mv CRC_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRC_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRC_2.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRD_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRD_2.re.fna.out -m genome -e 1e-3
mv CRD_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRD_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRD_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRE_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRE_2.re.fna.out -m genome -e 1e-3
mv CRE_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRE_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRE_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRG_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRG_2.re.fna.out -m genome -e 1e-3
mv CRG_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRG_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRG_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRH_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRH_2.re.fna.out -m genome -e 1e-3
mv CRH_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRH_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRH_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRI_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRI_2.re.fna.out -m genome -e 1e-3
mv CRI_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRI_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRI_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRK_2.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRK_2.re.fna.out -m genome -e 1e-3
mv CRK_2.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRK_2.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRK_2.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/CRL_.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o CRL_.re.fna.out -m genome -e 1e-3
mv CRL_.re.fna.out/short_summary.specific.saccharomycetes_odb10.CRL_.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r CRL_.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GAL.re.fna.out -m genome -e 1e-3
mv GAL.re.fna.out/short_summary.specific.saccharomycetes_odb10.GAL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GAL.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAP.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GAP.re.fna.out -m genome -e 1e-3
mv GAP.re.fna.out/short_summary.specific.saccharomycetes_odb10.GAP.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GAP.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAS.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GAS.re.fna.out -m genome -e 1e-3
mv GAS.re.fna.out/short_summary.specific.saccharomycetes_odb10.GAS.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GAS.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAT.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GAT.re.fna.out -m genome -e 1e-3
mv GAT.re.fna.out/short_summary.specific.saccharomycetes_odb10.GAT.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GAT.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GAV.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GAV.re.fna.out -m genome -e 1e-3
mv GAV.re.fna.out/short_summary.specific.saccharomycetes_odb10.GAV.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GAV.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000151485.1_ASM15148v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000151485.1_ASM15148v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000151485.1_ASM15148v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000151485.1_ASM15148v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000151485.1_ASM15148v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000181435.1_ASM18143v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000181435.1_ASM18143v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000181435.1_ASM18143v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000181435.1_ASM18143v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000181435.1_ASM18143v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182075.1_M22_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000182075.1_M22_genomic.fna.out -m genome -e 1e-3
mv GCA_000182075.1_M22_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000182075.1_M22_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000182075.1_M22_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182175.1_ASM18217v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000182175.1_ASM18217v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000182175.1_ASM18217v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000182175.1_ASM18217v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000182175.1_ASM18217v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000182315.2_ASM18231v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000182315.2_ASM18231v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000182315.2_ASM18231v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000182315.2_ASM18231v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000182315.2_ASM18231v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190195.1_ScAWRI796_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000190195.1_ScAWRI796_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000190195.1_ScAWRI796_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000190195.1_ScAWRI796_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000190195.1_ScAWRI796_v01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190215.1_ScVIN13_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000190215.1_ScVIN13_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000190215.1_ScVIN13_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000190215.1_ScVIN13_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000190215.1_ScVIN13_v01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190235.1_ScVL3_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000190235.1_ScVL3_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000190235.1_ScVL3_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000190235.1_ScVL3_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000190235.1_ScVL3_v01_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000190255.1_ScFostersB_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000190255.1_ScFostersB_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000190255.1_ScFostersB_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000190255.1_ScFostersB_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000190255.1_ScFostersB_v01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000192455.1_Saccharomyces_cerevisiae_CBS7960-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000192495.1_Saccharomyces_cerevisiae_CLIB324-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000192515.1_Saccharomyces_cerevisiae_YJM269-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000209265.1_Saccharomyces_cerevisiae_PW5-1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000209285.1_Saccharomyces_cerevisiae_UC5-1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000209305.1_Saccharomyces_cerevisiae_T73-1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000209345.1_Saccharomyces_cerevisiae_CLIB382-1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000209365.1_Saccharomyces_cerevisiae_T7-1.0_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000234495.1_EC9-8_1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000234495.1_EC9-8_1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000234495.1_EC9-8_1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000234495.1_EC9-8_1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000234495.1_EC9-8_1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000260735.1_K7_1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000260735.1_K7_1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000260735.1_K7_1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000260735.1_K7_1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000260735.1_K7_1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000275665.1_YJSH1_1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000275665.1_YJSH1_1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000275665.1_YJSH1_1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000275665.1_YJSH1_1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000275665.1_YJSH1_1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000308935.1_ZTW1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000308935.1_ZTW1_genomic.fna.out -m genome -e 1e-3
mv GCA_000308935.1_ZTW1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000308935.1_ZTW1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000308935.1_ZTW1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000325965.1_ScQA23_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000325965.1_ScQA23_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000325965.1_ScQA23_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000325965.1_ScQA23_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000325965.1_ScQA23_v01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000326005.1_ScFostersO_v01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000326005.1_ScFostersO_v01_genomic.fna.out -m genome -e 1e-3
mv GCA_000326005.1_ScFostersO_v01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000326005.1_ScFostersO_v01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000326005.1_ScFostersO_v01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000365045.1_Saccharomyces_cerevisiae_M3707_v1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000365065.1_Saccharomyces_cerevisiae_M3836_v1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000365085.1_Saccharomyces_cerevisiae_M3837_v1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000365105.1_Saccharomyces_cerevisiae_M3838_v1.0_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000365125.1_Saccharomyces_cerevisiae_M3839_v1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna.out -m genome -e 1e-3
mv GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000416405.1_S._cerevisiae_NY1308_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna.out -m genome -e 1e-3
mv GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000508805.2_ScNAM34-4C_assembly01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000508825.2_ScIR-2_assembly01_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000508825.2_ScIR-2_assembly01_genomic.fna.out -m genome -e 1e-3
mv GCA_000508825.2_ScIR-2_assembly01_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000508825.2_ScIR-2_assembly01_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000508825.2_ScIR-2_assembly01_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568005.1_R008_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000568005.1_R008_genomic.fna.out -m genome -e 1e-3
mv GCA_000568005.1_R008_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000568005.1_R008_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000568005.1_R008_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568055.1_p301_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000568055.1_p301_genomic.fna.out -m genome -e 1e-3
mv GCA_000568055.1_p301_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000568055.1_p301_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000568055.1_p301_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568295.1_P283_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000568295.1_P283_genomic.fna.out -m genome -e 1e-3
mv GCA_000568295.1_P283_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000568295.1_P283_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000568295.1_P283_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000568365.1_R103_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000568365.1_R103_genomic.fna.out -m genome -e 1e-3
mv GCA_000568365.1_R103_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000568365.1_R103_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000568365.1_R103_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000662435.2_Sc_YJM993_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000662435.2_Sc_YJM993_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000662435.2_Sc_YJM993_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000662435.2_Sc_YJM993_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000662435.2_Sc_YJM993_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000723645.1_ZhangWP_SECN85_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000723645.1_ZhangWP_SECN85_genomic.fna.out -m genome -e 1e-3
mv GCA_000723645.1_ZhangWP_SECN85_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000723645.1_ZhangWP_SECN85_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000723645.1_ZhangWP_SECN85_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000733235.3_ASM73323v3_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000733235.3_ASM73323v3_genomic.fna.out -m genome -e 1e-3
mv GCA_000733235.3_ASM73323v3_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000733235.3_ASM73323v3_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000733235.3_ASM73323v3_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000756235.1_ASM75623v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000756235.1_ASM75623v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000756235.1_ASM75623v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000756235.1_ASM75623v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000756235.1_ASM75623v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000756245.1_ASM75624v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000756245.1_ASM75624v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000756245.1_ASM75624v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000756245.1_ASM75624v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000756245.1_ASM75624v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766165.2_ASM76616v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766165.2_ASM76616v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766165.2_ASM76616v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766165.2_ASM76616v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766165.2_ASM76616v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766175.2_ASM76617v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766175.2_ASM76617v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766175.2_ASM76617v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766175.2_ASM76617v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766175.2_ASM76617v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766185.2_ASM76618v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766185.2_ASM76618v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766185.2_ASM76618v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766185.2_ASM76618v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766185.2_ASM76618v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766245.2_ASM76624v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766245.2_ASM76624v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766245.2_ASM76624v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766245.2_ASM76624v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766245.2_ASM76624v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766265.2_ASM76626v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766265.2_ASM76626v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766265.2_ASM76626v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766265.2_ASM76626v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766265.2_ASM76626v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766305.2_ASM76630v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766305.2_ASM76630v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766305.2_ASM76630v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766305.2_ASM76630v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766305.2_ASM76630v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766315.2_ASM76631v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766315.2_ASM76631v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766315.2_ASM76631v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766315.2_ASM76631v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766315.2_ASM76631v2_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766375.2_ASM76637v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766375.2_ASM76637v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766375.2_ASM76637v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766375.2_ASM76637v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766375.2_ASM76637v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766395.2_ASM76639v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766395.2_ASM76639v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766395.2_ASM76639v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766395.2_ASM76639v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766395.2_ASM76639v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766415.2_ASM76641v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766415.2_ASM76641v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766415.2_ASM76641v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766415.2_ASM76641v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766415.2_ASM76641v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766455.2_ASM76645v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766455.2_ASM76645v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766455.2_ASM76645v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766455.2_ASM76645v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766455.2_ASM76645v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766495.2_ASM76649v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766495.2_ASM76649v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766495.2_ASM76649v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766495.2_ASM76649v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766495.2_ASM76649v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766515.2_ASM76651v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766515.2_ASM76651v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766515.2_ASM76651v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766515.2_ASM76651v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766515.2_ASM76651v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766535.2_ASM76653v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766535.2_ASM76653v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766535.2_ASM76653v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766535.2_ASM76653v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766535.2_ASM76653v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766555.2_ASM76655v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766555.2_ASM76655v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766555.2_ASM76655v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766555.2_ASM76655v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766555.2_ASM76655v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766575.2_ASM76657v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766575.2_ASM76657v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766575.2_ASM76657v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766575.2_ASM76657v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766575.2_ASM76657v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766595.2_ASM76659v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766595.2_ASM76659v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766595.2_ASM76659v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766595.2_ASM76659v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766595.2_ASM76659v2_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766615.2_ASM76661v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766615.2_ASM76661v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766615.2_ASM76661v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766615.2_ASM76661v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766615.2_ASM76661v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000766635.2_ASM76663v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000766635.2_ASM76663v2_genomic.fna.out -m genome -e 1e-3
mv GCA_000766635.2_ASM76663v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000766635.2_ASM76663v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000766635.2_ASM76663v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000767965.1_ASM76796v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000767965.1_ASM76796v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000767965.1_ASM76796v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000767965.1_ASM76796v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000767965.1_ASM76796v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000768095.1_ASM76809v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000768095.1_ASM76809v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000768095.1_ASM76809v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000768095.1_ASM76809v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000768095.1_ASM76809v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000773925.1_ASM77392v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000773925.1_ASM77392v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000773925.1_ASM77392v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000773925.1_ASM77392v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000773925.1_ASM77392v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000775555.1_ASM77555v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000775555.1_ASM77555v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000775555.1_ASM77555v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000775555.1_ASM77555v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000775555.1_ASM77555v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975585.2_Sc_YJM195_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975585.2_Sc_YJM195_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975585.2_Sc_YJM195_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975585.2_Sc_YJM195_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975585.2_Sc_YJM195_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975615.2_Sc_YJM244_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975615.2_Sc_YJM244_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975615.2_Sc_YJM244_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975615.2_Sc_YJM244_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975615.2_Sc_YJM244_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975645.3_Sc_YJM1078_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975645.3_Sc_YJM1078_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975675.2_Sc_YJM1083_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975675.2_Sc_YJM1083_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975705.2_Sc_YJM1129_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975705.2_Sc_YJM1129_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975735.4_Sc_YJM189_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975735.4_Sc_YJM189_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975735.4_Sc_YJM189_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975735.4_Sc_YJM189_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975735.4_Sc_YJM189_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975765.2_Sc_YJM193_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975765.2_Sc_YJM193_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975765.2_Sc_YJM193_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975765.2_Sc_YJM193_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975765.2_Sc_YJM193_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975795.2_Sc_YJM248_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975795.2_Sc_YJM248_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975795.2_Sc_YJM248_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975795.2_Sc_YJM248_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975795.2_Sc_YJM248_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975825.2_Sc_YJM270_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975825.2_Sc_YJM270_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975825.2_Sc_YJM270_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975825.2_Sc_YJM270_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975825.2_Sc_YJM270_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975855.3_Sc_YJM271_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975855.3_Sc_YJM271_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975855.3_Sc_YJM271_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975855.3_Sc_YJM271_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975855.3_Sc_YJM271_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975885.2_Sc_YJM320_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975885.2_Sc_YJM320_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975885.2_Sc_YJM320_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975885.2_Sc_YJM320_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975885.2_Sc_YJM320_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975915.2_Sc_YJM326_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975915.2_Sc_YJM326_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975915.2_Sc_YJM326_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975915.2_Sc_YJM326_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975915.2_Sc_YJM326_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975945.2_Sc_YJM428_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975945.2_Sc_YJM428_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975945.2_Sc_YJM428_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975945.2_Sc_YJM428_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975945.2_Sc_YJM428_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000975975.3_Sc_YJM450_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000975975.3_Sc_YJM450_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000975975.3_Sc_YJM450_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000975975.3_Sc_YJM450_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000975975.3_Sc_YJM450_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976005.2_Sc_YJM451_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976005.2_Sc_YJM451_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976005.2_Sc_YJM451_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976005.2_Sc_YJM451_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976005.2_Sc_YJM451_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976035.2_Sc_YJM453_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976035.2_Sc_YJM453_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976035.2_Sc_YJM453_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976035.2_Sc_YJM453_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976035.2_Sc_YJM453_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976065.2_Sc_YJM456_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976065.2_Sc_YJM456_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976065.2_Sc_YJM456_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976065.2_Sc_YJM456_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976065.2_Sc_YJM456_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976095.2_Sc_YJM470_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976095.2_Sc_YJM470_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976095.2_Sc_YJM470_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976095.2_Sc_YJM470_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976095.2_Sc_YJM470_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976125.2_Sc_YJM541_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976125.2_Sc_YJM541_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976125.2_Sc_YJM541_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976125.2_Sc_YJM541_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976125.2_Sc_YJM541_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976155.4_Sc_YJM554_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976155.4_Sc_YJM554_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976155.4_Sc_YJM554_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976155.4_Sc_YJM554_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976155.4_Sc_YJM554_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976185.2_Sc_YJM555_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976185.2_Sc_YJM555_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976185.2_Sc_YJM555_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976185.2_Sc_YJM555_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976185.2_Sc_YJM555_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976215.2_Sc_YJM627_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976215.2_Sc_YJM627_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976215.2_Sc_YJM627_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976215.2_Sc_YJM627_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976215.2_Sc_YJM627_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976245.2_Sc_YJM681_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976245.2_Sc_YJM681_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976245.2_Sc_YJM681_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976245.2_Sc_YJM681_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976245.2_Sc_YJM681_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976275.4_Sc_YJM682_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976275.4_Sc_YJM682_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976275.4_Sc_YJM682_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976275.4_Sc_YJM682_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976275.4_Sc_YJM682_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976305.2_Sc_YJM683_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976305.2_Sc_YJM683_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976305.2_Sc_YJM683_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976305.2_Sc_YJM683_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976305.2_Sc_YJM683_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976335.2_Sc_YJM689_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976335.2_Sc_YJM689_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976335.2_Sc_YJM689_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976335.2_Sc_YJM689_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976335.2_Sc_YJM689_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976365.2_Sc_YJM693_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976365.2_Sc_YJM693_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976365.2_Sc_YJM693_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976365.2_Sc_YJM693_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976365.2_Sc_YJM693_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976395.2_Sc_YJM969_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976395.2_Sc_YJM969_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976395.2_Sc_YJM969_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976395.2_Sc_YJM969_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976395.2_Sc_YJM969_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976425.2_Sc_YJM972_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976425.2_Sc_YJM972_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976425.2_Sc_YJM972_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976425.2_Sc_YJM972_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976425.2_Sc_YJM972_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976455.3_Sc_YJM975_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976455.3_Sc_YJM975_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976455.3_Sc_YJM975_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976455.3_Sc_YJM975_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976455.3_Sc_YJM975_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976485.2_Sc_YJM978_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976485.2_Sc_YJM978_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976485.2_Sc_YJM978_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976485.2_Sc_YJM978_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976485.2_Sc_YJM978_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976515.2_Sc_YJM981_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976515.2_Sc_YJM981_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976515.2_Sc_YJM981_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976515.2_Sc_YJM981_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976515.2_Sc_YJM981_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976545.2_Sc_YJM984_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976545.2_Sc_YJM984_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976545.2_Sc_YJM984_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976545.2_Sc_YJM984_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976545.2_Sc_YJM984_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976575.4_Sc_YJM987_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976575.4_Sc_YJM987_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976575.4_Sc_YJM987_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976575.4_Sc_YJM987_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976575.4_Sc_YJM987_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976605.3_Sc_YJM990_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976605.3_Sc_YJM990_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976605.3_Sc_YJM990_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976605.3_Sc_YJM990_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976605.3_Sc_YJM990_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976665.2_Sc_YJM996_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976665.2_Sc_YJM996_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976665.2_Sc_YJM996_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976665.2_Sc_YJM996_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976665.2_Sc_YJM996_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976695.3_Sc_YJM1133_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976695.3_Sc_YJM1133_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976725.3_Sc_YJM1190_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976725.3_Sc_YJM1190_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976755.2_Sc_YJM1199_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976755.2_Sc_YJM1199_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976785.2_Sc_YJM1202_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976785.2_Sc_YJM1202_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976815.3_Sc_YJM1208_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976815.3_Sc_YJM1208_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976845.3_Sc_YJM1242_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976845.3_Sc_YJM1242_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976875.3_Sc_YJM1244_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976875.3_Sc_YJM1244_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976905.2_Sc_YJM1248_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976905.2_Sc_YJM1248_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976935.2_Sc_YJM1250_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976935.2_Sc_YJM1250_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976965.3_Sc_YJM1252_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976965.3_Sc_YJM1252_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000976995.2_Sc_YJM1273_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000976995.2_Sc_YJM1273_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977025.2_Sc_YJM1304_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977025.2_Sc_YJM1304_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977055.3_Sc_YJM1307_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977055.3_Sc_YJM1307_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977085.2_Sc_YJM1311_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977085.2_Sc_YJM1311_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977115.3_Sc_YJM1326_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977115.3_Sc_YJM1326_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977145.2_Sc_YJM1332_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977145.2_Sc_YJM1332_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977175.3_Sc_YJM1336_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977175.3_Sc_YJM1336_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977205.2_Sc_YJM1338_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977205.2_Sc_YJM1338_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977235.2_Sc_YJM1341_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977235.2_Sc_YJM1341_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977265.3_Sc_YJM1342_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977265.3_Sc_YJM1342_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977295.2_Sc_YJM1355_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977295.2_Sc_YJM1355_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977325.2_Sc_YJM1356_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977325.2_Sc_YJM1356_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977355.3_Sc_YJM1381_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977355.3_Sc_YJM1381_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977385.2_Sc_YJM1383_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977385.2_Sc_YJM1383_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977415.2_Sc_YJM1385_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977415.2_Sc_YJM1385_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977445.2_Sc_YJM1386_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977445.2_Sc_YJM1386_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977475.3_Sc_YJM1387_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977475.3_Sc_YJM1387_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977505.2_Sc_YJM1388_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977505.2_Sc_YJM1388_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977535.2_Sc_YJM1389_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977535.2_Sc_YJM1389_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977565.2_Sc_YJM1399_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977565.2_Sc_YJM1399_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977595.2_Sc_YJM1400_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977595.2_Sc_YJM1400_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977625.2_Sc_YJM1401_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977625.2_Sc_YJM1401_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977655.2_Sc_YJM1402_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977655.2_Sc_YJM1402_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977685.2_Sc_YJM1415_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977685.2_Sc_YJM1415_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977715.4_Sc_YJM1417_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977715.4_Sc_YJM1417_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977745.2_Sc_YJM1418_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977745.2_Sc_YJM1418_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977775.3_Sc_YJM1419_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977775.3_Sc_YJM1419_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977805.2_Sc_YJM1433_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977805.2_Sc_YJM1433_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977835.2_Sc_YJM1434_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977835.2_Sc_YJM1434_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977865.2_Sc_YJM1439_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977865.2_Sc_YJM1439_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977895.2_Sc_YJM1443_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977895.2_Sc_YJM1443_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977925.2_Sc_YJM1444_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977925.2_Sc_YJM1444_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977955.2_Sc_YJM1447_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977955.2_Sc_YJM1447_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000977985.2_Sc_YJM1450_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000977985.2_Sc_YJM1450_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978015.2_Sc_YJM1460_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978015.2_Sc_YJM1460_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978045.2_Sc_YJM1463_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978045.2_Sc_YJM1463_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978075.3_Sc_YJM1477_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978075.3_Sc_YJM1477_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978105.2_Sc_YJM1478_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978105.2_Sc_YJM1478_v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978135.2_Sc_YJM1479_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978135.2_Sc_YJM1479_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978165.2_Sc_YJM1526_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978165.2_Sc_YJM1526_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978195.2_Sc_YJM1527_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978195.2_Sc_YJM1527_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978225.2_Sc_YJM1549_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978225.2_Sc_YJM1549_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978255.2_Sc_YJM1573_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978255.2_Sc_YJM1573_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978285.2_Sc_YJM1574_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978285.2_Sc_YJM1574_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978315.2_Sc_YJM1592_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978315.2_Sc_YJM1592_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_000978345.2_Sc_YJM1615_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_000978345.2_Sc_YJM1615_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001006245.1_ScSp60_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001006245.1_ScSp60_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001006245.1_ScSp60_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001006245.1_ScSp60_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001006245.1_ScSp60_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001029075.1_ASM102907v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001029075.1_ASM102907v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001029075.1_ASM102907v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001029075.1_ASM102907v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001029075.1_ASM102907v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001029995.1_ASM102999v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001029995.1_ASM102999v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001029995.1_ASM102999v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001029995.1_ASM102999v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001029995.1_ASM102999v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001051215.1_ASM105121v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001051215.1_ASM105121v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001051215.1_ASM105121v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001051215.1_ASM105121v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001051215.1_ASM105121v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001078085.1_ASM107808v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001078085.1_ASM107808v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001078085.1_ASM107808v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001078085.1_ASM107808v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001078085.1_ASM107808v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001078105.1_ASM107810v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001078105.1_ASM107810v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001078105.1_ASM107810v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001078105.1_ASM107810v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001078105.1_ASM107810v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001282415.1_Y-12632v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001282415.1_Y-12632v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001282415.1_Y-12632v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001282415.1_Y-12632v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001282415.1_Y-12632v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001578265.1_S.cerevisiae-74-D694-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001578275.1_S.cerevisiae-15V-P4-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001578285.1_S.cerevisiae-1B-D1606-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001578295.1_25-25-2V-P3982-SPAdes-2.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna.out -m genome -e 1e-3
mv GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001578385.1_S.cerevisiae-6P-33G-D373-2.0_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001584535.1_ASM158453v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001584535.1_ASM158453v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001584535.1_ASM158453v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001584535.1_ASM158453v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001584535.1_ASM158453v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001592655.1_ASM159265v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001592655.1_ASM159265v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001592655.1_ASM159265v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001592655.1_ASM159265v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001592655.1_ASM159265v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001611835.1_ASM161183v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001611835.1_ASM161183v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001611835.1_ASM161183v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001611835.1_ASM161183v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001611835.1_ASM161183v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001634645.1_GLBRCY22-3_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001634645.1_GLBRCY22-3_genomic.fna.out -m genome -e 1e-3
mv GCA_001634645.1_GLBRCY22-3_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001634645.1_GLBRCY22-3_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001634645.1_GLBRCY22-3_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001669875.1_ASM166987v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001669875.1_ASM166987v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001669875.1_ASM166987v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001669875.1_ASM166987v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001669875.1_ASM166987v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001669905.1_ASM166990v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001669905.1_ASM166990v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001669905.1_ASM166990v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001669905.1_ASM166990v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001669905.1_ASM166990v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001669935.1_ASM166993v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001669935.1_ASM166993v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001669935.1_ASM166993v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001669935.1_ASM166993v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001669935.1_ASM166993v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001669965.1_ASM166996v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001669965.1_ASM166996v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001669965.1_ASM166996v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001669965.1_ASM166996v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001669965.1_ASM166996v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001669995.1_ASM166999v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001669995.1_ASM166999v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001669995.1_ASM166999v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001669995.1_ASM166999v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001669995.1_ASM166999v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670025.1_ASM167002v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670025.1_ASM167002v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670025.1_ASM167002v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670025.1_ASM167002v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670025.1_ASM167002v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670055.1_ASM167005v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670055.1_ASM167005v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670055.1_ASM167005v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670055.1_ASM167005v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670055.1_ASM167005v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670085.1_ASM167008v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670085.1_ASM167008v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670085.1_ASM167008v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670085.1_ASM167008v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670085.1_ASM167008v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670115.1_ASM167011v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670115.1_ASM167011v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670115.1_ASM167011v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670115.1_ASM167011v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670115.1_ASM167011v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670145.1_ASM167014v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670145.1_ASM167014v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670145.1_ASM167014v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670145.1_ASM167014v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670145.1_ASM167014v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670175.1_ASM167017v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670175.1_ASM167017v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670175.1_ASM167017v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670175.1_ASM167017v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670175.1_ASM167017v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670205.1_ASM167020v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670205.1_ASM167020v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670205.1_ASM167020v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670205.1_ASM167020v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670205.1_ASM167020v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670235.1_ASM167023v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670235.1_ASM167023v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670235.1_ASM167023v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670235.1_ASM167023v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670235.1_ASM167023v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670265.1_ASM167026v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670265.1_ASM167026v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670265.1_ASM167026v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670265.1_ASM167026v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670265.1_ASM167026v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670295.1_ASM167029v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670295.1_ASM167029v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670295.1_ASM167029v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670295.1_ASM167029v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670295.1_ASM167029v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670325.1_ASM167032v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670325.1_ASM167032v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670325.1_ASM167032v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670325.1_ASM167032v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670325.1_ASM167032v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670355.1_ASM167035v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670355.1_ASM167035v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670355.1_ASM167035v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670355.1_ASM167035v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670355.1_ASM167035v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670385.1_ASM167038v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670385.1_ASM167038v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670385.1_ASM167038v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670385.1_ASM167038v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670385.1_ASM167038v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670415.1_ASM167041v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670415.1_ASM167041v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670415.1_ASM167041v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670415.1_ASM167041v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670415.1_ASM167041v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670845.1_ASM167084v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670845.1_ASM167084v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670845.1_ASM167084v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670845.1_ASM167084v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670845.1_ASM167084v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670935.1_ASM167093v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670935.1_ASM167093v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670935.1_ASM167093v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670935.1_ASM167093v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670935.1_ASM167093v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670965.1_ASM167096v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670965.1_ASM167096v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670965.1_ASM167096v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670965.1_ASM167096v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670965.1_ASM167096v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001670995.1_ASM167099v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001670995.1_ASM167099v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001670995.1_ASM167099v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001670995.1_ASM167099v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001670995.1_ASM167099v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671025.1_ASM167102v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671025.1_ASM167102v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671025.1_ASM167102v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671025.1_ASM167102v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671025.1_ASM167102v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671055.1_ASM167105v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671055.1_ASM167105v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671055.1_ASM167105v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671055.1_ASM167105v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671055.1_ASM167105v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671085.1_ASM167108v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671085.1_ASM167108v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671085.1_ASM167108v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671085.1_ASM167108v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671085.1_ASM167108v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671115.1_ASM167111v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671115.1_ASM167111v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671115.1_ASM167111v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671115.1_ASM167111v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671115.1_ASM167111v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671165.1_ASM167116v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671165.1_ASM167116v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671165.1_ASM167116v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671165.1_ASM167116v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671165.1_ASM167116v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671195.1_ASM167119v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671195.1_ASM167119v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671195.1_ASM167119v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671195.1_ASM167119v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671195.1_ASM167119v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671225.1_ASM167122v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671225.1_ASM167122v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671225.1_ASM167122v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671225.1_ASM167122v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671225.1_ASM167122v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671255.1_ASM167125v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671255.1_ASM167125v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671255.1_ASM167125v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671255.1_ASM167125v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671255.1_ASM167125v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671285.1_ASM167128v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671285.1_ASM167128v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671285.1_ASM167128v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671285.1_ASM167128v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671285.1_ASM167128v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671315.1_ASM167131v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671315.1_ASM167131v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671315.1_ASM167131v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671315.1_ASM167131v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671315.1_ASM167131v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671345.1_ASM167134v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671345.1_ASM167134v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671345.1_ASM167134v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671345.1_ASM167134v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671345.1_ASM167134v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671375.1_ASM167137v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671375.1_ASM167137v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671375.1_ASM167137v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671375.1_ASM167137v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671375.1_ASM167137v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671405.1_ASM167140v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671405.1_ASM167140v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671405.1_ASM167140v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671405.1_ASM167140v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671405.1_ASM167140v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671435.1_ASM167143v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671435.1_ASM167143v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671435.1_ASM167143v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671435.1_ASM167143v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671435.1_ASM167143v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671565.1_ASM167156v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671565.1_ASM167156v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671565.1_ASM167156v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671565.1_ASM167156v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671565.1_ASM167156v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671595.1_ASM167159v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671595.1_ASM167159v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671595.1_ASM167159v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671595.1_ASM167159v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671595.1_ASM167159v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671625.1_ASM167162v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671625.1_ASM167162v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671625.1_ASM167162v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671625.1_ASM167162v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671625.1_ASM167162v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671655.1_ASM167165v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671655.1_ASM167165v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671655.1_ASM167165v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671655.1_ASM167165v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671655.1_ASM167165v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671685.1_ASM167168v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671685.1_ASM167168v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671685.1_ASM167168v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671685.1_ASM167168v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671685.1_ASM167168v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671715.1_ASM167171v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671715.1_ASM167171v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671715.1_ASM167171v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671715.1_ASM167171v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671715.1_ASM167171v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671745.1_ASM167174v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671745.1_ASM167174v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671745.1_ASM167174v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671745.1_ASM167174v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671745.1_ASM167174v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671775.1_ASM167177v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671775.1_ASM167177v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671775.1_ASM167177v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671775.1_ASM167177v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671775.1_ASM167177v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671805.1_ASM167180v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671805.1_ASM167180v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671805.1_ASM167180v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671805.1_ASM167180v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671805.1_ASM167180v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671835.1_ASM167183v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671835.1_ASM167183v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671835.1_ASM167183v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671835.1_ASM167183v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671835.1_ASM167183v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671865.1_ASM167186v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671865.1_ASM167186v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671865.1_ASM167186v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671865.1_ASM167186v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671865.1_ASM167186v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671895.1_ASM167189v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671895.1_ASM167189v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671895.1_ASM167189v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671895.1_ASM167189v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671895.1_ASM167189v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671925.1_ASM167192v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671925.1_ASM167192v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671925.1_ASM167192v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671925.1_ASM167192v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671925.1_ASM167192v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671955.1_ASM167195v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671955.1_ASM167195v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671955.1_ASM167195v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671955.1_ASM167195v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671955.1_ASM167195v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001671985.1_ASM167198v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001671985.1_ASM167198v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001671985.1_ASM167198v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001671985.1_ASM167198v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001671985.1_ASM167198v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001735825.1_ASM173582v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001735825.1_ASM173582v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001735825.1_ASM173582v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001735825.1_ASM173582v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001735825.1_ASM173582v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737145.2_ASM173714v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737145.2_ASM173714v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737145.2_ASM173714v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737145.2_ASM173714v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737145.2_ASM173714v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737155.1_ASM173715v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737155.1_ASM173715v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737155.1_ASM173715v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737155.1_ASM173715v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737155.1_ASM173715v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737165.1_ASM173716v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737165.1_ASM173716v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737165.1_ASM173716v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737165.1_ASM173716v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737165.1_ASM173716v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737175.1_ASM173717v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737175.1_ASM173717v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737175.1_ASM173717v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737175.1_ASM173717v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737175.1_ASM173717v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737265.1_ASM173726v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737265.1_ASM173726v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737265.1_ASM173726v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737265.1_ASM173726v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737265.1_ASM173726v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737275.1_ASM173727v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737275.1_ASM173727v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737275.1_ASM173727v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737275.1_ASM173727v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737275.1_ASM173727v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737285.1_ASM173728v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737285.1_ASM173728v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737285.1_ASM173728v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737285.1_ASM173728v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737285.1_ASM173728v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737295.2_ASM173729v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737295.2_ASM173729v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737295.2_ASM173729v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737295.2_ASM173729v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737295.2_ASM173729v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737385.2_ASM173738v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737385.2_ASM173738v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737385.2_ASM173738v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737385.2_ASM173738v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737385.2_ASM173738v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737395.2_ASM173739v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737395.2_ASM173739v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737395.2_ASM173739v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737395.2_ASM173739v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737395.2_ASM173739v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737415.2_ASM173741v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737415.2_ASM173741v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737415.2_ASM173741v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737415.2_ASM173741v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737415.2_ASM173741v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737445.2_ASM173744v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737445.2_ASM173744v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737445.2_ASM173744v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737445.2_ASM173744v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737445.2_ASM173744v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737505.2_ASM173750v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737505.2_ASM173750v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737505.2_ASM173750v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737505.2_ASM173750v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737505.2_ASM173750v2_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737515.2_ASM173751v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737515.2_ASM173751v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737515.2_ASM173751v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737515.2_ASM173751v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737515.2_ASM173751v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737535.2_ASM173753v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737535.2_ASM173753v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737535.2_ASM173753v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737535.2_ASM173753v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737535.2_ASM173753v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737545.2_ASM173754v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737545.2_ASM173754v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737545.2_ASM173754v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737545.2_ASM173754v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737545.2_ASM173754v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737625.2_ASM173762v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737625.2_ASM173762v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737625.2_ASM173762v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737625.2_ASM173762v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737625.2_ASM173762v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737635.2_ASM173763v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737635.2_ASM173763v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737635.2_ASM173763v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737635.2_ASM173763v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737635.2_ASM173763v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737645.2_ASM173764v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737645.2_ASM173764v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737645.2_ASM173764v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737645.2_ASM173764v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737645.2_ASM173764v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737655.2_ASM173765v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737655.2_ASM173765v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737655.2_ASM173765v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737655.2_ASM173765v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737655.2_ASM173765v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737745.2_ASM173774v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737745.2_ASM173774v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737745.2_ASM173774v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737745.2_ASM173774v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737745.2_ASM173774v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737755.2_ASM173775v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737755.2_ASM173775v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737755.2_ASM173775v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737755.2_ASM173775v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737755.2_ASM173775v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737765.2_ASM173776v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737765.2_ASM173776v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737765.2_ASM173776v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737765.2_ASM173776v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737765.2_ASM173776v2_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737785.2_ASM173778v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737785.2_ASM173778v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737785.2_ASM173778v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737785.2_ASM173778v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737785.2_ASM173778v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737865.2_ASM173786v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737865.2_ASM173786v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001737865.2_ASM173786v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737865.2_ASM173786v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737865.2_ASM173786v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737875.1_ASM173787v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737875.1_ASM173787v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737875.1_ASM173787v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737875.1_ASM173787v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737875.1_ASM173787v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737885.1_ASM173788v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737885.1_ASM173788v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737885.1_ASM173788v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737885.1_ASM173788v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737885.1_ASM173788v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737895.1_ASM173789v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737895.1_ASM173789v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737895.1_ASM173789v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737895.1_ASM173789v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737895.1_ASM173789v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737985.1_ASM173798v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737985.1_ASM173798v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737985.1_ASM173798v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737985.1_ASM173798v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737985.1_ASM173798v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001737995.1_ASM173799v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001737995.1_ASM173799v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001737995.1_ASM173799v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001737995.1_ASM173799v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001737995.1_ASM173799v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738005.1_ASM173800v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738005.1_ASM173800v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738005.1_ASM173800v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738005.1_ASM173800v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738005.1_ASM173800v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738015.1_ASM173801v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738015.1_ASM173801v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738015.1_ASM173801v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738015.1_ASM173801v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738015.1_ASM173801v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738105.1_ASM173810v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738105.1_ASM173810v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738105.1_ASM173810v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738105.1_ASM173810v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738105.1_ASM173810v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738115.1_ASM173811v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738115.1_ASM173811v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738115.1_ASM173811v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738115.1_ASM173811v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738115.1_ASM173811v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738125.1_ASM173812v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738125.1_ASM173812v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738125.1_ASM173812v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738125.1_ASM173812v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738125.1_ASM173812v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738135.2_ASM173813v2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738135.2_ASM173813v2_genomic.fna.out -m genome -e 1e-3
mv GCA_001738135.2_ASM173813v2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738135.2_ASM173813v2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738135.2_ASM173813v2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738225.1_ASM173822v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738225.1_ASM173822v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738225.1_ASM173822v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738225.1_ASM173822v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738225.1_ASM173822v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738235.1_ASM173823v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738235.1_ASM173823v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738235.1_ASM173823v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738235.1_ASM173823v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738235.1_ASM173823v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738255.1_ASM173825v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738255.1_ASM173825v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738255.1_ASM173825v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738255.1_ASM173825v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738255.1_ASM173825v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738265.1_ASM173826v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738265.1_ASM173826v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738265.1_ASM173826v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738265.1_ASM173826v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738265.1_ASM173826v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738345.1_ASM173834v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738345.1_ASM173834v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738345.1_ASM173834v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738345.1_ASM173834v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738345.1_ASM173834v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738355.1_ASM173835v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738355.1_ASM173835v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738355.1_ASM173835v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738355.1_ASM173835v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738355.1_ASM173835v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738375.1_ASM173837v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738375.1_ASM173837v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738375.1_ASM173837v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738375.1_ASM173837v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738375.1_ASM173837v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738395.1_ASM173839v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738395.1_ASM173839v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738395.1_ASM173839v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738395.1_ASM173839v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738395.1_ASM173839v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738465.1_ASM173846v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738465.1_ASM173846v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738465.1_ASM173846v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738465.1_ASM173846v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738465.1_ASM173846v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738485.1_ASM173848v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738485.1_ASM173848v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738485.1_ASM173848v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738485.1_ASM173848v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738485.1_ASM173848v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738495.1_ASM173849v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738495.1_ASM173849v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738495.1_ASM173849v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738495.1_ASM173849v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738495.1_ASM173849v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738515.1_ASM173851v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738515.1_ASM173851v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738515.1_ASM173851v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738515.1_ASM173851v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738515.1_ASM173851v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738585.1_ASM173858v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738585.1_ASM173858v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738585.1_ASM173858v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738585.1_ASM173858v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738585.1_ASM173858v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738595.1_ASM173859v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738595.1_ASM173859v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738595.1_ASM173859v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738595.1_ASM173859v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738595.1_ASM173859v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738605.1_ASM173860v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738605.1_ASM173860v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738605.1_ASM173860v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738605.1_ASM173860v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738605.1_ASM173860v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738615.1_ASM173861v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738615.1_ASM173861v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738615.1_ASM173861v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738615.1_ASM173861v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738615.1_ASM173861v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738705.1_ASM173870v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738705.1_ASM173870v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738705.1_ASM173870v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738705.1_ASM173870v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738705.1_ASM173870v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738715.1_ASM173871v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738715.1_ASM173871v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738715.1_ASM173871v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738715.1_ASM173871v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738715.1_ASM173871v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738725.1_ASM173872v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738725.1_ASM173872v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738725.1_ASM173872v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738725.1_ASM173872v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738725.1_ASM173872v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738755.1_ASM173875v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738755.1_ASM173875v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738755.1_ASM173875v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738755.1_ASM173875v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738755.1_ASM173875v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738825.1_ASM173882v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738825.1_ASM173882v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738825.1_ASM173882v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738825.1_ASM173882v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738825.1_ASM173882v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738835.1_ASM173883v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738835.1_ASM173883v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738835.1_ASM173883v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738835.1_ASM173883v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738835.1_ASM173883v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738845.1_ASM173884v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738845.1_ASM173884v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738845.1_ASM173884v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738845.1_ASM173884v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738845.1_ASM173884v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738855.1_ASM173885v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738855.1_ASM173885v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738855.1_ASM173885v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738855.1_ASM173885v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738855.1_ASM173885v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738945.1_ASM173894v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738945.1_ASM173894v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738945.1_ASM173894v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738945.1_ASM173894v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738945.1_ASM173894v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738955.1_ASM173895v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738955.1_ASM173895v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738955.1_ASM173895v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738955.1_ASM173895v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738955.1_ASM173895v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738965.1_ASM173896v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738965.1_ASM173896v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738965.1_ASM173896v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738965.1_ASM173896v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738965.1_ASM173896v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001738985.1_ASM173898v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001738985.1_ASM173898v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001738985.1_ASM173898v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001738985.1_ASM173898v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001738985.1_ASM173898v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739065.1_ASM173906v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739065.1_ASM173906v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739065.1_ASM173906v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739065.1_ASM173906v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739065.1_ASM173906v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739075.1_ASM173907v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739075.1_ASM173907v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739075.1_ASM173907v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739075.1_ASM173907v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739075.1_ASM173907v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739095.1_ASM173909v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739095.1_ASM173909v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739095.1_ASM173909v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739095.1_ASM173909v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739095.1_ASM173909v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739105.1_ASM173910v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739105.1_ASM173910v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739105.1_ASM173910v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739105.1_ASM173910v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739105.1_ASM173910v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739185.1_ASM173918v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739185.1_ASM173918v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739185.1_ASM173918v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739185.1_ASM173918v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739185.1_ASM173918v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739215.1_ASM173921v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739215.1_ASM173921v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739215.1_ASM173921v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739215.1_ASM173921v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739215.1_ASM173921v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739225.1_ASM173922v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739225.1_ASM173922v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739225.1_ASM173922v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739225.1_ASM173922v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739225.1_ASM173922v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739235.1_ASM173923v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739235.1_ASM173923v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739235.1_ASM173923v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739235.1_ASM173923v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739235.1_ASM173923v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739305.1_ASM173930v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739305.1_ASM173930v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739305.1_ASM173930v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739305.1_ASM173930v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739305.1_ASM173930v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739335.1_ASM173933v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739335.1_ASM173933v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739335.1_ASM173933v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739335.1_ASM173933v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739335.1_ASM173933v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739345.1_ASM173934v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739345.1_ASM173934v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739345.1_ASM173934v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739345.1_ASM173934v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739345.1_ASM173934v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739355.1_ASM173935v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739355.1_ASM173935v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739355.1_ASM173935v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739355.1_ASM173935v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739355.1_ASM173935v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739425.1_ASM173942v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739425.1_ASM173942v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739425.1_ASM173942v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739425.1_ASM173942v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739425.1_ASM173942v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739455.1_ASM173945v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739455.1_ASM173945v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739455.1_ASM173945v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739455.1_ASM173945v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739455.1_ASM173945v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739475.1_ASM173947v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739475.1_ASM173947v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739475.1_ASM173947v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739475.1_ASM173947v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739475.1_ASM173947v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739495.1_ASM173949v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739495.1_ASM173949v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739495.1_ASM173949v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739495.1_ASM173949v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739495.1_ASM173949v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739545.1_ASM173954v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739545.1_ASM173954v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739545.1_ASM173954v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739545.1_ASM173954v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739545.1_ASM173954v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739575.1_ASM173957v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739575.1_ASM173957v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739575.1_ASM173957v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739575.1_ASM173957v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739575.1_ASM173957v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739595.1_ASM173959v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739595.1_ASM173959v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739595.1_ASM173959v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739595.1_ASM173959v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739595.1_ASM173959v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739625.1_ASM173962v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739625.1_ASM173962v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739625.1_ASM173962v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739625.1_ASM173962v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739625.1_ASM173962v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739665.1_ASM173966v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739665.1_ASM173966v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739665.1_ASM173966v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739665.1_ASM173966v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739665.1_ASM173966v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739695.1_ASM173969v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739695.1_ASM173969v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739695.1_ASM173969v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739695.1_ASM173969v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739695.1_ASM173969v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739705.1_ASM173970v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739705.1_ASM173970v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739705.1_ASM173970v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739705.1_ASM173970v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739705.1_ASM173970v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739735.1_ASM173973v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739735.1_ASM173973v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739735.1_ASM173973v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739735.1_ASM173973v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739735.1_ASM173973v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739775.1_ASM173977v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739775.1_ASM173977v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739775.1_ASM173977v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739775.1_ASM173977v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739775.1_ASM173977v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739815.1_ASM173981v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739815.1_ASM173981v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739815.1_ASM173981v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739815.1_ASM173981v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739815.1_ASM173981v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739825.1_ASM173982v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739825.1_ASM173982v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739825.1_ASM173982v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739825.1_ASM173982v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739825.1_ASM173982v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739845.1_ASM173984v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739845.1_ASM173984v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739845.1_ASM173984v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739845.1_ASM173984v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739845.1_ASM173984v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739865.1_ASM173986v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739865.1_ASM173986v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739865.1_ASM173986v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739865.1_ASM173986v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739865.1_ASM173986v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739935.1_ASM173993v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739935.1_ASM173993v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739935.1_ASM173993v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739935.1_ASM173993v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739935.1_ASM173993v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739945.1_ASM173994v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739945.1_ASM173994v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739945.1_ASM173994v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739945.1_ASM173994v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739945.1_ASM173994v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739965.1_ASM173996v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739965.1_ASM173996v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739965.1_ASM173996v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739965.1_ASM173996v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739965.1_ASM173996v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001739995.1_ASM173999v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001739995.1_ASM173999v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001739995.1_ASM173999v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001739995.1_ASM173999v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001739995.1_ASM173999v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740055.1_ASM174005v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740055.1_ASM174005v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740055.1_ASM174005v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740055.1_ASM174005v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740055.1_ASM174005v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740075.1_ASM174007v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740075.1_ASM174007v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740075.1_ASM174007v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740075.1_ASM174007v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740075.1_ASM174007v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740095.1_ASM174009v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740095.1_ASM174009v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740095.1_ASM174009v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740095.1_ASM174009v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740095.1_ASM174009v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740115.1_ASM174011v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740115.1_ASM174011v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740115.1_ASM174011v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740115.1_ASM174011v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740115.1_ASM174011v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740175.1_ASM174017v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740175.1_ASM174017v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740175.1_ASM174017v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740175.1_ASM174017v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740175.1_ASM174017v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740185.1_ASM174018v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740185.1_ASM174018v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740185.1_ASM174018v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740185.1_ASM174018v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740185.1_ASM174018v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740205.1_ASM174020v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740205.1_ASM174020v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740205.1_ASM174020v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740205.1_ASM174020v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740205.1_ASM174020v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740225.1_ASM174022v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740225.1_ASM174022v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740225.1_ASM174022v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740225.1_ASM174022v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740225.1_ASM174022v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740295.1_ASM174029v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740295.1_ASM174029v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740295.1_ASM174029v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740295.1_ASM174029v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740295.1_ASM174029v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740305.1_ASM174030v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740305.1_ASM174030v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740305.1_ASM174030v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740305.1_ASM174030v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740305.1_ASM174030v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740325.1_ASM174032v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740325.1_ASM174032v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740325.1_ASM174032v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740325.1_ASM174032v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740325.1_ASM174032v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740345.1_ASM174034v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740345.1_ASM174034v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740345.1_ASM174034v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740345.1_ASM174034v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740345.1_ASM174034v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740415.1_ASM174041v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740415.1_ASM174041v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740415.1_ASM174041v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740415.1_ASM174041v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740415.1_ASM174041v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740425.1_ASM174042v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740425.1_ASM174042v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740425.1_ASM174042v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740425.1_ASM174042v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740425.1_ASM174042v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740445.1_ASM174044v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740445.1_ASM174044v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740445.1_ASM174044v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740445.1_ASM174044v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740445.1_ASM174044v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740455.1_ASM174045v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740455.1_ASM174045v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740455.1_ASM174045v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740455.1_ASM174045v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740455.1_ASM174045v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740535.1_ASM174053v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740535.1_ASM174053v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740535.1_ASM174053v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740535.1_ASM174053v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740535.1_ASM174053v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740545.1_ASM174054v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740545.1_ASM174054v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740545.1_ASM174054v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740545.1_ASM174054v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740545.1_ASM174054v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740565.1_ASM174056v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740565.1_ASM174056v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740565.1_ASM174056v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740565.1_ASM174056v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740565.1_ASM174056v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740575.1_ASM174057v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740575.1_ASM174057v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740575.1_ASM174057v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740575.1_ASM174057v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740575.1_ASM174057v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740655.1_ASM174065v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740655.1_ASM174065v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740655.1_ASM174065v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740655.1_ASM174065v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740655.1_ASM174065v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740665.1_ASM174066v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740665.1_ASM174066v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740665.1_ASM174066v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740665.1_ASM174066v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740665.1_ASM174066v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740675.1_ASM174067v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740675.1_ASM174067v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740675.1_ASM174067v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740675.1_ASM174067v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740675.1_ASM174067v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740705.1_ASM174070v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740705.1_ASM174070v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740705.1_ASM174070v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740705.1_ASM174070v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740705.1_ASM174070v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740775.1_ASM174077v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740775.1_ASM174077v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740775.1_ASM174077v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740775.1_ASM174077v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740775.1_ASM174077v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740785.1_ASM174078v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740785.1_ASM174078v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740785.1_ASM174078v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740785.1_ASM174078v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740785.1_ASM174078v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740795.1_ASM174079v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740795.1_ASM174079v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740795.1_ASM174079v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740795.1_ASM174079v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740795.1_ASM174079v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740805.1_ASM174080v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740805.1_ASM174080v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740805.1_ASM174080v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740805.1_ASM174080v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740805.1_ASM174080v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740895.1_ASM174089v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740895.1_ASM174089v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740895.1_ASM174089v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740895.1_ASM174089v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740895.1_ASM174089v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740905.1_ASM174090v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740905.1_ASM174090v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740905.1_ASM174090v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740905.1_ASM174090v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740905.1_ASM174090v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740915.1_ASM174091v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740915.1_ASM174091v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740915.1_ASM174091v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740915.1_ASM174091v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740915.1_ASM174091v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001740925.1_ASM174092v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001740925.1_ASM174092v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001740925.1_ASM174092v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001740925.1_ASM174092v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001740925.1_ASM174092v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741015.1_ASM174101v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741015.1_ASM174101v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741015.1_ASM174101v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741015.1_ASM174101v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741015.1_ASM174101v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741025.1_ASM174102v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741025.1_ASM174102v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741025.1_ASM174102v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741025.1_ASM174102v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741025.1_ASM174102v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741035.1_ASM174103v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741035.1_ASM174103v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741035.1_ASM174103v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741035.1_ASM174103v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741035.1_ASM174103v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741045.1_ASM174104v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741045.1_ASM174104v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741045.1_ASM174104v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741045.1_ASM174104v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741045.1_ASM174104v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741135.1_ASM174113v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741135.1_ASM174113v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741135.1_ASM174113v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741135.1_ASM174113v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741135.1_ASM174113v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741145.1_ASM174114v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741145.1_ASM174114v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741145.1_ASM174114v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741145.1_ASM174114v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741145.1_ASM174114v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741165.1_ASM174116v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741165.1_ASM174116v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741165.1_ASM174116v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741165.1_ASM174116v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741165.1_ASM174116v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741195.1_ASM174119v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741195.1_ASM174119v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741195.1_ASM174119v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741195.1_ASM174119v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741195.1_ASM174119v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741255.1_ASM174125v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741255.1_ASM174125v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741255.1_ASM174125v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741255.1_ASM174125v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741255.1_ASM174125v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741285.1_ASM174128v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741285.1_ASM174128v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741285.1_ASM174128v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741285.1_ASM174128v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741285.1_ASM174128v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741295.1_ASM174129v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741295.1_ASM174129v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741295.1_ASM174129v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741295.1_ASM174129v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741295.1_ASM174129v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741305.1_ASM174130v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741305.1_ASM174130v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741305.1_ASM174130v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741305.1_ASM174130v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741305.1_ASM174130v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741375.1_ASM174137v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741375.1_ASM174137v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741375.1_ASM174137v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741375.1_ASM174137v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741375.1_ASM174137v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741405.1_ASM174140v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741405.1_ASM174140v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741405.1_ASM174140v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741405.1_ASM174140v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741405.1_ASM174140v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741415.1_ASM174141v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741415.1_ASM174141v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741415.1_ASM174141v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741415.1_ASM174141v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741415.1_ASM174141v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741425.1_ASM174142v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741425.1_ASM174142v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741425.1_ASM174142v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741425.1_ASM174142v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741425.1_ASM174142v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741495.1_ASM174149v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741495.1_ASM174149v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741495.1_ASM174149v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741495.1_ASM174149v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741495.1_ASM174149v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741525.1_ASM174152v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741525.1_ASM174152v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741525.1_ASM174152v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741525.1_ASM174152v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741525.1_ASM174152v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741535.1_ASM174153v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741535.1_ASM174153v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741535.1_ASM174153v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741535.1_ASM174153v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741535.1_ASM174153v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741555.1_ASM174155v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741555.1_ASM174155v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741555.1_ASM174155v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741555.1_ASM174155v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741555.1_ASM174155v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741605.1_ASM174160v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741605.1_ASM174160v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741605.1_ASM174160v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741605.1_ASM174160v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741605.1_ASM174160v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741685.1_ASM174168v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741685.1_ASM174168v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741685.1_ASM174168v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741685.1_ASM174168v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741685.1_ASM174168v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741695.1_ASM174169v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741695.1_ASM174169v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741695.1_ASM174169v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741695.1_ASM174169v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741695.1_ASM174169v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741705.1_ASM174170v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741705.1_ASM174170v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741705.1_ASM174170v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741705.1_ASM174170v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741705.1_ASM174170v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741725.1_ASM174172v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741725.1_ASM174172v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741725.1_ASM174172v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741725.1_ASM174172v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741725.1_ASM174172v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741805.1_ASM174180v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741805.1_ASM174180v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741805.1_ASM174180v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741805.1_ASM174180v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741805.1_ASM174180v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001741815.1_ASM174181v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001741815.1_ASM174181v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001741815.1_ASM174181v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001741815.1_ASM174181v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001741815.1_ASM174181v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001750485.1_S_cerevisiae_GUJ105_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001932575.1_CTBEBG1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001932575.1_CTBEBG1_genomic.fna.out -m genome -e 1e-3
mv GCA_001932575.1_CTBEBG1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001932575.1_CTBEBG1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001932575.1_CTBEBG1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001937245.1_ASM193724v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001937245.1_ASM193724v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001937245.1_ASM193724v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001937245.1_ASM193724v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001937245.1_ASM193724v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_001983315.1_ASM198331v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_001983315.1_ASM198331v1_genomic.fna.out -m genome -e 1e-3
mv GCA_001983315.1_ASM198331v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_001983315.1_ASM198331v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_001983315.1_ASM198331v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057805.1_ASM205780v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002057805.1_ASM205780v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002057805.1_ASM205780v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002057805.1_ASM205780v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002057805.1_ASM205780v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057885.1_ASM205788v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002057885.1_ASM205788v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002057885.1_ASM205788v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002057885.1_ASM205788v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002057885.1_ASM205788v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002057995.1_ASM205799v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002057995.1_ASM205799v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002057995.1_ASM205799v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002057995.1_ASM205799v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002057995.1_ASM205799v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002058095.1_ASM205809v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002058095.1_ASM205809v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002058095.1_ASM205809v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002058095.1_ASM205809v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002058095.1_ASM205809v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002163515.1_ASM216351v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002163515.1_ASM216351v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002163515.1_ASM216351v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002163515.1_ASM216351v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002163515.1_ASM216351v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002214855.1_ASM221485v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002214855.1_ASM221485v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002214855.1_ASM221485v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002214855.1_ASM221485v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002214855.1_ASM221485v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002271905.1_S_cerevisiae_DBVPG6044_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002271945.1_S_cerevisiae_Y12_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna.out -m genome -e 1e-3
mv GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002571405.1_CENPK113_7D_MinION_1.0_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002804325.1_ASM280432v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002804325.1_ASM280432v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002804325.1_ASM280432v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002804325.1_ASM280432v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002804325.1_ASM280432v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002807125.1_S_cerevisiae_SPSC01_v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_002933895.1_ASM293389v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_002933895.1_ASM293389v1_genomic.fna.out -m genome -e 1e-3
mv GCA_002933895.1_ASM293389v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_002933895.1_ASM293389v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_002933895.1_ASM293389v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003020065.1_ASM302006v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003020065.1_ASM302006v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003020065.1_ASM302006v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003020065.1_ASM302006v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003020065.1_ASM302006v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046695.1_ASM304669v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003046695.1_ASM304669v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003046695.1_ASM304669v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003046695.1_ASM304669v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003046695.1_ASM304669v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046715.1_ASM304671v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003046715.1_ASM304671v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003046715.1_ASM304671v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003046715.1_ASM304671v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003046715.1_ASM304671v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003046745.1_ASM304674v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003046745.1_ASM304674v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003046745.1_ASM304674v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003046745.1_ASM304674v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003046745.1_ASM304674v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003086655.1_ASM308665v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003086655.1_ASM308665v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003086655.1_ASM308665v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003086655.1_ASM308665v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003086655.1_ASM308665v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270095.1_ASM327009v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270095.1_ASM327009v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270095.1_ASM327009v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270095.1_ASM327009v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270095.1_ASM327009v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270135.1_ASM327013v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270135.1_ASM327013v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270135.1_ASM327013v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270135.1_ASM327013v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270135.1_ASM327013v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270145.1_ASM327014v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270145.1_ASM327014v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270145.1_ASM327014v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270145.1_ASM327014v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270145.1_ASM327014v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270175.1_ASM327017v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270175.1_ASM327017v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270175.1_ASM327017v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270175.1_ASM327017v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270175.1_ASM327017v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270185.1_ASM327018v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270185.1_ASM327018v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270185.1_ASM327018v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270185.1_ASM327018v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270185.1_ASM327018v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270205.1_ASM327020v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270205.1_ASM327020v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270205.1_ASM327020v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270205.1_ASM327020v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270205.1_ASM327020v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270245.1_ASM327024v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270245.1_ASM327024v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270245.1_ASM327024v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270245.1_ASM327024v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270245.1_ASM327024v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270305.1_ASM327030v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270305.1_ASM327030v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270305.1_ASM327030v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270305.1_ASM327030v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270305.1_ASM327030v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270315.1_ASM327031v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270315.1_ASM327031v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270315.1_ASM327031v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270315.1_ASM327031v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270315.1_ASM327031v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270325.1_ASM327032v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270325.1_ASM327032v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270325.1_ASM327032v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270325.1_ASM327032v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270325.1_ASM327032v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270335.1_ASM327033v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270335.1_ASM327033v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270335.1_ASM327033v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270335.1_ASM327033v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270335.1_ASM327033v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270425.1_ASM327042v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270425.1_ASM327042v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270425.1_ASM327042v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270425.1_ASM327042v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270425.1_ASM327042v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270435.1_ASM327043v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270435.1_ASM327043v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270435.1_ASM327043v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270435.1_ASM327043v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270435.1_ASM327043v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270455.1_ASM327045v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270455.1_ASM327045v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270455.1_ASM327045v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270455.1_ASM327045v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270455.1_ASM327045v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270465.1_ASM327046v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270465.1_ASM327046v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270465.1_ASM327046v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270465.1_ASM327046v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270465.1_ASM327046v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270535.1_ASM327053v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270535.1_ASM327053v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270535.1_ASM327053v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270535.1_ASM327053v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270535.1_ASM327053v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270565.1_ASM327056v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270565.1_ASM327056v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270565.1_ASM327056v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270565.1_ASM327056v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270565.1_ASM327056v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270575.1_ASM327057v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270575.1_ASM327057v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270575.1_ASM327057v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270575.1_ASM327057v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270575.1_ASM327057v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270595.1_ASM327059v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270595.1_ASM327059v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270595.1_ASM327059v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270595.1_ASM327059v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270595.1_ASM327059v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270665.1_ASM327066v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270665.1_ASM327066v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270665.1_ASM327066v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270665.1_ASM327066v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270665.1_ASM327066v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270695.1_ASM327069v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270695.1_ASM327069v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270695.1_ASM327069v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270695.1_ASM327069v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270695.1_ASM327069v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270705.1_ASM327070v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270705.1_ASM327070v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270705.1_ASM327070v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270705.1_ASM327070v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270705.1_ASM327070v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270715.1_ASM327071v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270715.1_ASM327071v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270715.1_ASM327071v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270715.1_ASM327071v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270715.1_ASM327071v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270765.1_ASM327076v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270765.1_ASM327076v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270765.1_ASM327076v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270765.1_ASM327076v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270765.1_ASM327076v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270815.1_ASM327081v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270815.1_ASM327081v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270815.1_ASM327081v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270815.1_ASM327081v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270815.1_ASM327081v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270825.1_ASM327082v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270825.1_ASM327082v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270825.1_ASM327082v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270825.1_ASM327082v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270825.1_ASM327082v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270835.1_ASM327083v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270835.1_ASM327083v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270835.1_ASM327083v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270835.1_ASM327083v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270835.1_ASM327083v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270855.1_ASM327085v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270855.1_ASM327085v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270855.1_ASM327085v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270855.1_ASM327085v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270855.1_ASM327085v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270935.1_ASM327093v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270935.1_ASM327093v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270935.1_ASM327093v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270935.1_ASM327093v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270935.1_ASM327093v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270945.1_ASM327094v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270945.1_ASM327094v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270945.1_ASM327094v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270945.1_ASM327094v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270945.1_ASM327094v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270955.1_ASM327095v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270955.1_ASM327095v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270955.1_ASM327095v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270955.1_ASM327095v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270955.1_ASM327095v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003270985.1_ASM327098v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003270985.1_ASM327098v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003270985.1_ASM327098v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003270985.1_ASM327098v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003270985.1_ASM327098v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271045.1_ASM327104v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271045.1_ASM327104v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271045.1_ASM327104v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271045.1_ASM327104v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271045.1_ASM327104v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271085.1_ASM327108v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271085.1_ASM327108v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271085.1_ASM327108v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271085.1_ASM327108v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271085.1_ASM327108v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271115.1_ASM327111v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271115.1_ASM327111v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271115.1_ASM327111v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271115.1_ASM327111v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271115.1_ASM327111v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271125.1_ASM327112v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271125.1_ASM327112v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271125.1_ASM327112v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271125.1_ASM327112v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271125.1_ASM327112v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271135.1_ASM327113v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271135.1_ASM327113v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271135.1_ASM327113v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271135.1_ASM327113v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271135.1_ASM327113v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271185.1_ASM327118v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271185.1_ASM327118v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271185.1_ASM327118v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271185.1_ASM327118v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271185.1_ASM327118v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271225.1_ASM327122v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271225.1_ASM327122v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271225.1_ASM327122v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271225.1_ASM327122v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271225.1_ASM327122v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271245.1_ASM327124v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271245.1_ASM327124v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271245.1_ASM327124v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271245.1_ASM327124v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271245.1_ASM327124v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271255.1_ASM327125v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271255.1_ASM327125v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271255.1_ASM327125v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271255.1_ASM327125v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271255.1_ASM327125v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271315.1_ASM327131v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271315.1_ASM327131v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271315.1_ASM327131v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271315.1_ASM327131v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271315.1_ASM327131v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271355.1_ASM327135v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271355.1_ASM327135v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271355.1_ASM327135v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271355.1_ASM327135v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271355.1_ASM327135v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271365.1_ASM327136v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271365.1_ASM327136v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271365.1_ASM327136v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271365.1_ASM327136v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271365.1_ASM327136v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271395.1_ASM327139v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271395.1_ASM327139v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271395.1_ASM327139v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271395.1_ASM327139v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271395.1_ASM327139v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271425.1_ASM327142v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271425.1_ASM327142v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271425.1_ASM327142v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271425.1_ASM327142v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271425.1_ASM327142v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271475.1_ASM327147v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271475.1_ASM327147v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271475.1_ASM327147v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271475.1_ASM327147v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271475.1_ASM327147v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271505.1_ASM327150v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271505.1_ASM327150v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271505.1_ASM327150v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271505.1_ASM327150v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271505.1_ASM327150v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271535.1_ASM327153v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271535.1_ASM327153v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271535.1_ASM327153v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271535.1_ASM327153v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271535.1_ASM327153v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271545.1_ASM327154v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271545.1_ASM327154v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271545.1_ASM327154v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271545.1_ASM327154v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271545.1_ASM327154v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271555.1_ASM327155v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271555.1_ASM327155v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271555.1_ASM327155v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271555.1_ASM327155v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271555.1_ASM327155v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271565.1_ASM327156v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271565.1_ASM327156v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271565.1_ASM327156v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271565.1_ASM327156v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271565.1_ASM327156v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271645.1_ASM327164v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271645.1_ASM327164v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271645.1_ASM327164v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271645.1_ASM327164v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271645.1_ASM327164v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271675.1_ASM327167v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271675.1_ASM327167v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271675.1_ASM327167v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271675.1_ASM327167v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271675.1_ASM327167v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271695.1_ASM327169v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271695.1_ASM327169v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271695.1_ASM327169v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271695.1_ASM327169v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271695.1_ASM327169v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271705.1_ASM327170v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271705.1_ASM327170v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271705.1_ASM327170v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271705.1_ASM327170v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271705.1_ASM327170v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271775.1_ASM327177v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271775.1_ASM327177v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271775.1_ASM327177v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271775.1_ASM327177v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271775.1_ASM327177v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271805.1_ASM327180v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271805.1_ASM327180v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271805.1_ASM327180v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271805.1_ASM327180v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271805.1_ASM327180v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271815.1_ASM327181v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271815.1_ASM327181v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271815.1_ASM327181v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271815.1_ASM327181v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271815.1_ASM327181v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271825.1_ASM327182v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271825.1_ASM327182v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271825.1_ASM327182v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271825.1_ASM327182v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271825.1_ASM327182v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271835.1_ASM327183v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271835.1_ASM327183v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271835.1_ASM327183v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271835.1_ASM327183v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271835.1_ASM327183v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271925.1_ASM327192v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271925.1_ASM327192v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271925.1_ASM327192v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271925.1_ASM327192v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271925.1_ASM327192v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271935.1_ASM327193v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271935.1_ASM327193v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271935.1_ASM327193v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271935.1_ASM327193v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271935.1_ASM327193v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271945.1_ASM327194v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271945.1_ASM327194v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271945.1_ASM327194v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271945.1_ASM327194v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271945.1_ASM327194v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003271955.1_ASM327195v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003271955.1_ASM327195v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003271955.1_ASM327195v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003271955.1_ASM327195v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003271955.1_ASM327195v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272045.1_ASM327204v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272045.1_ASM327204v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272045.1_ASM327204v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272045.1_ASM327204v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272045.1_ASM327204v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272055.1_ASM327205v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272055.1_ASM327205v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272055.1_ASM327205v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272055.1_ASM327205v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272055.1_ASM327205v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272065.1_ASM327206v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272065.1_ASM327206v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272065.1_ASM327206v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272065.1_ASM327206v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272065.1_ASM327206v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272075.1_ASM327207v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272075.1_ASM327207v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272075.1_ASM327207v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272075.1_ASM327207v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272075.1_ASM327207v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272165.1_ASM327216v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272165.1_ASM327216v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272165.1_ASM327216v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272165.1_ASM327216v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272165.1_ASM327216v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272175.1_ASM327217v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272175.1_ASM327217v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272175.1_ASM327217v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272175.1_ASM327217v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272175.1_ASM327217v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272185.1_ASM327218v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272185.1_ASM327218v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272185.1_ASM327218v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272185.1_ASM327218v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272185.1_ASM327218v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272195.1_ASM327219v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272195.1_ASM327219v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272195.1_ASM327219v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272195.1_ASM327219v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272195.1_ASM327219v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272285.1_ASM327228v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272285.1_ASM327228v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272285.1_ASM327228v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272285.1_ASM327228v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272285.1_ASM327228v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272295.1_ASM327229v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272295.1_ASM327229v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272295.1_ASM327229v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272295.1_ASM327229v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272295.1_ASM327229v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272305.1_ASM327230v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272305.1_ASM327230v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272305.1_ASM327230v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272305.1_ASM327230v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272305.1_ASM327230v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272315.1_ASM327231v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272315.1_ASM327231v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272315.1_ASM327231v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272315.1_ASM327231v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272315.1_ASM327231v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272405.1_ASM327240v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272405.1_ASM327240v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272405.1_ASM327240v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272405.1_ASM327240v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272405.1_ASM327240v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272415.1_ASM327241v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272415.1_ASM327241v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272415.1_ASM327241v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272415.1_ASM327241v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272415.1_ASM327241v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272425.1_ASM327242v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272425.1_ASM327242v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272425.1_ASM327242v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272425.1_ASM327242v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272425.1_ASM327242v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272445.1_ASM327244v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272445.1_ASM327244v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272445.1_ASM327244v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272445.1_ASM327244v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272445.1_ASM327244v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272525.1_ASM327252v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272525.1_ASM327252v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272525.1_ASM327252v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272525.1_ASM327252v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272525.1_ASM327252v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272535.1_ASM327253v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272535.1_ASM327253v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272535.1_ASM327253v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272535.1_ASM327253v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272535.1_ASM327253v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272575.1_ASM327257v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272575.1_ASM327257v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272575.1_ASM327257v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272575.1_ASM327257v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272575.1_ASM327257v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272595.1_ASM327259v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272595.1_ASM327259v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272595.1_ASM327259v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272595.1_ASM327259v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272595.1_ASM327259v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272625.1_ASM327262v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272625.1_ASM327262v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272625.1_ASM327262v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272625.1_ASM327262v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272625.1_ASM327262v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272665.1_ASM327266v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272665.1_ASM327266v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272665.1_ASM327266v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272665.1_ASM327266v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272665.1_ASM327266v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272705.1_ASM327270v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272705.1_ASM327270v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272705.1_ASM327270v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272705.1_ASM327270v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272705.1_ASM327270v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272715.1_ASM327271v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272715.1_ASM327271v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272715.1_ASM327271v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272715.1_ASM327271v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272715.1_ASM327271v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272725.1_ASM327272v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272725.1_ASM327272v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272725.1_ASM327272v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272725.1_ASM327272v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272725.1_ASM327272v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272755.1_ASM327275v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272755.1_ASM327275v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272755.1_ASM327275v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272755.1_ASM327275v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272755.1_ASM327275v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272815.1_ASM327281v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272815.1_ASM327281v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272815.1_ASM327281v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272815.1_ASM327281v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272815.1_ASM327281v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272835.1_ASM327283v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272835.1_ASM327283v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272835.1_ASM327283v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272835.1_ASM327283v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272835.1_ASM327283v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272855.1_ASM327285v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272855.1_ASM327285v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272855.1_ASM327285v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272855.1_ASM327285v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272855.1_ASM327285v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272905.1_ASM327290v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272905.1_ASM327290v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272905.1_ASM327290v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272905.1_ASM327290v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272905.1_ASM327290v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272925.1_ASM327292v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272925.1_ASM327292v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272925.1_ASM327292v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272925.1_ASM327292v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272925.1_ASM327292v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272945.1_ASM327294v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272945.1_ASM327294v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272945.1_ASM327294v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272945.1_ASM327294v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272945.1_ASM327294v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003272975.1_ASM327297v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003272975.1_ASM327297v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003272975.1_ASM327297v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003272975.1_ASM327297v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003272975.1_ASM327297v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273025.1_ASM327302v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273025.1_ASM327302v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273025.1_ASM327302v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273025.1_ASM327302v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273025.1_ASM327302v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273035.1_ASM327303v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273035.1_ASM327303v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273035.1_ASM327303v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273035.1_ASM327303v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273035.1_ASM327303v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273045.1_ASM327304v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273045.1_ASM327304v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273045.1_ASM327304v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273045.1_ASM327304v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273045.1_ASM327304v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273075.1_ASM327307v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273075.1_ASM327307v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273075.1_ASM327307v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273075.1_ASM327307v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273075.1_ASM327307v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273145.1_ASM327314v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273145.1_ASM327314v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273145.1_ASM327314v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273145.1_ASM327314v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273145.1_ASM327314v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273155.1_ASM327315v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273155.1_ASM327315v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273155.1_ASM327315v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273155.1_ASM327315v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273155.1_ASM327315v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273165.1_ASM327316v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273165.1_ASM327316v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273165.1_ASM327316v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273165.1_ASM327316v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273165.1_ASM327316v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273175.1_ASM327317v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273175.1_ASM327317v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273175.1_ASM327317v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273175.1_ASM327317v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273175.1_ASM327317v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273265.1_ASM327326v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273265.1_ASM327326v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273265.1_ASM327326v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273265.1_ASM327326v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273265.1_ASM327326v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273275.1_ASM327327v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273275.1_ASM327327v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273275.1_ASM327327v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273275.1_ASM327327v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273275.1_ASM327327v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273285.1_ASM327328v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273285.1_ASM327328v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273285.1_ASM327328v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273285.1_ASM327328v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273285.1_ASM327328v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273295.1_ASM327329v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273295.1_ASM327329v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273295.1_ASM327329v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273295.1_ASM327329v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273295.1_ASM327329v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273385.1_ASM327338v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273385.1_ASM327338v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273385.1_ASM327338v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273385.1_ASM327338v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273385.1_ASM327338v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273395.1_ASM327339v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273395.1_ASM327339v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273395.1_ASM327339v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273395.1_ASM327339v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273395.1_ASM327339v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273405.1_ASM327340v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273405.1_ASM327340v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273405.1_ASM327340v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273405.1_ASM327340v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273405.1_ASM327340v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273415.1_ASM327341v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273415.1_ASM327341v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273415.1_ASM327341v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273415.1_ASM327341v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273415.1_ASM327341v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273495.1_ASM327349v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273495.1_ASM327349v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273495.1_ASM327349v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273495.1_ASM327349v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273495.1_ASM327349v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273535.1_ASM327353v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273535.1_ASM327353v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273535.1_ASM327353v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273535.1_ASM327353v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273535.1_ASM327353v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273555.1_ASM327355v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273555.1_ASM327355v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273555.1_ASM327355v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273555.1_ASM327355v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273555.1_ASM327355v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273575.1_ASM327357v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273575.1_ASM327357v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273575.1_ASM327357v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273575.1_ASM327357v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273575.1_ASM327357v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273585.1_ASM327358v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273585.1_ASM327358v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273585.1_ASM327358v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273585.1_ASM327358v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273585.1_ASM327358v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273625.1_ASM327362v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273625.1_ASM327362v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273625.1_ASM327362v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273625.1_ASM327362v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273625.1_ASM327362v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273685.1_ASM327368v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273685.1_ASM327368v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273685.1_ASM327368v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273685.1_ASM327368v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273685.1_ASM327368v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273695.1_ASM327369v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273695.1_ASM327369v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273695.1_ASM327369v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273695.1_ASM327369v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273695.1_ASM327369v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273705.1_ASM327370v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273705.1_ASM327370v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273705.1_ASM327370v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273705.1_ASM327370v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273705.1_ASM327370v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273725.1_ASM327372v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273725.1_ASM327372v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273725.1_ASM327372v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273725.1_ASM327372v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273725.1_ASM327372v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273805.1_ASM327380v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273805.1_ASM327380v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273805.1_ASM327380v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273805.1_ASM327380v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273805.1_ASM327380v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273815.1_ASM327381v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273815.1_ASM327381v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273815.1_ASM327381v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273815.1_ASM327381v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273815.1_ASM327381v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273825.1_ASM327382v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273825.1_ASM327382v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273825.1_ASM327382v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273825.1_ASM327382v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273825.1_ASM327382v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273835.1_ASM327383v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273835.1_ASM327383v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273835.1_ASM327383v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273835.1_ASM327383v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273835.1_ASM327383v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273925.1_ASM327392v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273925.1_ASM327392v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273925.1_ASM327392v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273925.1_ASM327392v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273925.1_ASM327392v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273935.1_ASM327393v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273935.1_ASM327393v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273935.1_ASM327393v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273935.1_ASM327393v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273935.1_ASM327393v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273955.1_ASM327395v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273955.1_ASM327395v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273955.1_ASM327395v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273955.1_ASM327395v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273955.1_ASM327395v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003273985.1_ASM327398v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003273985.1_ASM327398v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003273985.1_ASM327398v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003273985.1_ASM327398v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003273985.1_ASM327398v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274025.1_ASM327402v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274025.1_ASM327402v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274025.1_ASM327402v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274025.1_ASM327402v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274025.1_ASM327402v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274055.1_ASM327405v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274055.1_ASM327405v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274055.1_ASM327405v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274055.1_ASM327405v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274055.1_ASM327405v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274085.1_ASM327408v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274085.1_ASM327408v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274085.1_ASM327408v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274085.1_ASM327408v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274085.1_ASM327408v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274095.1_ASM327409v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274095.1_ASM327409v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274095.1_ASM327409v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274095.1_ASM327409v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274095.1_ASM327409v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274125.1_ASM327412v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274125.1_ASM327412v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274125.1_ASM327412v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274125.1_ASM327412v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274125.1_ASM327412v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274195.1_ASM327419v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274195.1_ASM327419v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274195.1_ASM327419v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274195.1_ASM327419v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274195.1_ASM327419v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274205.1_ASM327420v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274205.1_ASM327420v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274205.1_ASM327420v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274205.1_ASM327420v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274205.1_ASM327420v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274215.1_ASM327421v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274215.1_ASM327421v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274215.1_ASM327421v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274215.1_ASM327421v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274215.1_ASM327421v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274225.1_ASM327422v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274225.1_ASM327422v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274225.1_ASM327422v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274225.1_ASM327422v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274225.1_ASM327422v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274315.1_ASM327431v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274315.1_ASM327431v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274315.1_ASM327431v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274315.1_ASM327431v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274315.1_ASM327431v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274325.1_ASM327432v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274325.1_ASM327432v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274325.1_ASM327432v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274325.1_ASM327432v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274325.1_ASM327432v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274345.1_ASM327434v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274345.1_ASM327434v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274345.1_ASM327434v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274345.1_ASM327434v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274345.1_ASM327434v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274365.1_ASM327436v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274365.1_ASM327436v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274365.1_ASM327436v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274365.1_ASM327436v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274365.1_ASM327436v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274425.1_ASM327442v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274425.1_ASM327442v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274425.1_ASM327442v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274425.1_ASM327442v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274425.1_ASM327442v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274465.1_ASM327446v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274465.1_ASM327446v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274465.1_ASM327446v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274465.1_ASM327446v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274465.1_ASM327446v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274475.1_ASM327447v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274475.1_ASM327447v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274475.1_ASM327447v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274475.1_ASM327447v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274475.1_ASM327447v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274485.1_ASM327448v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274485.1_ASM327448v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274485.1_ASM327448v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274485.1_ASM327448v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274485.1_ASM327448v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274495.1_ASM327449v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274495.1_ASM327449v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274495.1_ASM327449v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274495.1_ASM327449v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274495.1_ASM327449v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274585.1_ASM327458v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274585.1_ASM327458v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274585.1_ASM327458v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274585.1_ASM327458v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274585.1_ASM327458v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274595.1_ASM327459v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274595.1_ASM327459v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274595.1_ASM327459v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274595.1_ASM327459v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274595.1_ASM327459v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274605.1_ASM327460v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274605.1_ASM327460v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274605.1_ASM327460v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274605.1_ASM327460v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274605.1_ASM327460v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274615.1_ASM327461v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274615.1_ASM327461v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274615.1_ASM327461v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274615.1_ASM327461v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274615.1_ASM327461v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274705.1_ASM327470v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274705.1_ASM327470v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274705.1_ASM327470v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274705.1_ASM327470v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274705.1_ASM327470v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274715.1_ASM327471v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274715.1_ASM327471v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274715.1_ASM327471v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274715.1_ASM327471v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274715.1_ASM327471v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274725.1_ASM327472v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274725.1_ASM327472v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274725.1_ASM327472v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274725.1_ASM327472v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274725.1_ASM327472v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274735.1_ASM327473v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274735.1_ASM327473v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274735.1_ASM327473v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274735.1_ASM327473v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274735.1_ASM327473v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274825.1_ASM327482v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274825.1_ASM327482v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274825.1_ASM327482v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274825.1_ASM327482v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274825.1_ASM327482v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274845.1_ASM327484v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274845.1_ASM327484v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274845.1_ASM327484v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274845.1_ASM327484v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274845.1_ASM327484v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274855.1_ASM327485v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274855.1_ASM327485v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274855.1_ASM327485v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274855.1_ASM327485v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274855.1_ASM327485v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274865.1_ASM327486v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274865.1_ASM327486v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274865.1_ASM327486v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274865.1_ASM327486v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274865.1_ASM327486v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274905.1_ASM327490v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274905.1_ASM327490v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274905.1_ASM327490v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274905.1_ASM327490v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274905.1_ASM327490v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274975.1_ASM327497v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274975.1_ASM327497v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274975.1_ASM327497v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274975.1_ASM327497v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274975.1_ASM327497v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274985.1_ASM327498v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274985.1_ASM327498v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274985.1_ASM327498v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274985.1_ASM327498v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274985.1_ASM327498v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003274995.1_ASM327499v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003274995.1_ASM327499v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003274995.1_ASM327499v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003274995.1_ASM327499v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003274995.1_ASM327499v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275005.1_ASM327500v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275005.1_ASM327500v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275005.1_ASM327500v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275005.1_ASM327500v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275005.1_ASM327500v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275095.1_ASM327509v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275095.1_ASM327509v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275095.1_ASM327509v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275095.1_ASM327509v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275095.1_ASM327509v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275105.1_ASM327510v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275105.1_ASM327510v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275105.1_ASM327510v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275105.1_ASM327510v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275105.1_ASM327510v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275115.1_ASM327511v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275115.1_ASM327511v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275115.1_ASM327511v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275115.1_ASM327511v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275115.1_ASM327511v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275125.1_ASM327512v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275125.1_ASM327512v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275125.1_ASM327512v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275125.1_ASM327512v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275125.1_ASM327512v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275215.1_ASM327521v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275215.1_ASM327521v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275215.1_ASM327521v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275215.1_ASM327521v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275215.1_ASM327521v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275225.1_ASM327522v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275225.1_ASM327522v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275225.1_ASM327522v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275225.1_ASM327522v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275225.1_ASM327522v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275235.1_ASM327523v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275235.1_ASM327523v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275235.1_ASM327523v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275235.1_ASM327523v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275235.1_ASM327523v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275245.1_ASM327524v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275245.1_ASM327524v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275245.1_ASM327524v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275245.1_ASM327524v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275245.1_ASM327524v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275335.1_ASM327533v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275335.1_ASM327533v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275335.1_ASM327533v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275335.1_ASM327533v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275335.1_ASM327533v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275345.1_ASM327534v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275345.1_ASM327534v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275345.1_ASM327534v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275345.1_ASM327534v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275345.1_ASM327534v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275365.1_ASM327536v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275365.1_ASM327536v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275365.1_ASM327536v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275365.1_ASM327536v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275365.1_ASM327536v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275385.1_ASM327538v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275385.1_ASM327538v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275385.1_ASM327538v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275385.1_ASM327538v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275385.1_ASM327538v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275455.1_ASM327545v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275455.1_ASM327545v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275455.1_ASM327545v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275455.1_ASM327545v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275455.1_ASM327545v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275465.1_ASM327546v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275465.1_ASM327546v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275465.1_ASM327546v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275465.1_ASM327546v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275465.1_ASM327546v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275475.1_ASM327547v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275475.1_ASM327547v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275475.1_ASM327547v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275475.1_ASM327547v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275475.1_ASM327547v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275505.1_ASM327550v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275505.1_ASM327550v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275505.1_ASM327550v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275505.1_ASM327550v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275505.1_ASM327550v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275575.1_ASM327557v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275575.1_ASM327557v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275575.1_ASM327557v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275575.1_ASM327557v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275575.1_ASM327557v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275585.1_ASM327558v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275585.1_ASM327558v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275585.1_ASM327558v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275585.1_ASM327558v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275585.1_ASM327558v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275605.1_ASM327560v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275605.1_ASM327560v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275605.1_ASM327560v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275605.1_ASM327560v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275605.1_ASM327560v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275615.1_ASM327561v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275615.1_ASM327561v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275615.1_ASM327561v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275615.1_ASM327561v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275615.1_ASM327561v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275685.1_ASM327568v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275685.1_ASM327568v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275685.1_ASM327568v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275685.1_ASM327568v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275685.1_ASM327568v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275715.1_ASM327571v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275715.1_ASM327571v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275715.1_ASM327571v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275715.1_ASM327571v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275715.1_ASM327571v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275745.1_ASM327574v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275745.1_ASM327574v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275745.1_ASM327574v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275745.1_ASM327574v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275745.1_ASM327574v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275755.1_ASM327575v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275755.1_ASM327575v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275755.1_ASM327575v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275755.1_ASM327575v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275755.1_ASM327575v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275785.1_ASM327578v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275785.1_ASM327578v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275785.1_ASM327578v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275785.1_ASM327578v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275785.1_ASM327578v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275825.1_ASM327582v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275825.1_ASM327582v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275825.1_ASM327582v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275825.1_ASM327582v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275825.1_ASM327582v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275835.1_ASM327583v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275835.1_ASM327583v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275835.1_ASM327583v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275835.1_ASM327583v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275835.1_ASM327583v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275875.1_ASM327587v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275875.1_ASM327587v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275875.1_ASM327587v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275875.1_ASM327587v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275875.1_ASM327587v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275885.1_ASM327588v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275885.1_ASM327588v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275885.1_ASM327588v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275885.1_ASM327588v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275885.1_ASM327588v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275955.1_ASM327595v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275955.1_ASM327595v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275955.1_ASM327595v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275955.1_ASM327595v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275955.1_ASM327595v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003275965.1_ASM327596v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003275965.1_ASM327596v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003275965.1_ASM327596v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003275965.1_ASM327596v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003275965.1_ASM327596v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276005.1_ASM327600v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276005.1_ASM327600v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276005.1_ASM327600v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276005.1_ASM327600v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276005.1_ASM327600v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276015.1_ASM327601v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276015.1_ASM327601v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276015.1_ASM327601v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276015.1_ASM327601v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276015.1_ASM327601v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276075.1_ASM327607v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276075.1_ASM327607v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276075.1_ASM327607v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276075.1_ASM327607v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276075.1_ASM327607v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276085.1_ASM327608v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276085.1_ASM327608v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276085.1_ASM327608v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276085.1_ASM327608v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276085.1_ASM327608v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276135.1_ASM327613v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276135.1_ASM327613v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276135.1_ASM327613v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276135.1_ASM327613v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276135.1_ASM327613v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276145.1_ASM327614v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276145.1_ASM327614v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276145.1_ASM327614v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276145.1_ASM327614v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276145.1_ASM327614v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276185.1_ASM327618v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276185.1_ASM327618v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276185.1_ASM327618v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276185.1_ASM327618v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276185.1_ASM327618v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276205.1_ASM327620v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276205.1_ASM327620v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276205.1_ASM327620v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276205.1_ASM327620v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276205.1_ASM327620v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276235.1_ASM327623v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276235.1_ASM327623v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276235.1_ASM327623v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276235.1_ASM327623v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276235.1_ASM327623v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276245.1_ASM327624v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276245.1_ASM327624v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276245.1_ASM327624v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276245.1_ASM327624v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276245.1_ASM327624v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276295.1_ASM327629v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276295.1_ASM327629v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276295.1_ASM327629v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276295.1_ASM327629v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276295.1_ASM327629v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276325.1_ASM327632v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276325.1_ASM327632v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276325.1_ASM327632v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276325.1_ASM327632v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276325.1_ASM327632v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276365.1_ASM327636v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276365.1_ASM327636v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276365.1_ASM327636v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276365.1_ASM327636v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276365.1_ASM327636v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276375.1_ASM327637v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276375.1_ASM327637v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276375.1_ASM327637v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276375.1_ASM327637v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276375.1_ASM327637v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276395.1_ASM327639v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276395.1_ASM327639v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276395.1_ASM327639v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276395.1_ASM327639v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276395.1_ASM327639v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276405.1_ASM327640v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276405.1_ASM327640v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276405.1_ASM327640v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276405.1_ASM327640v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276405.1_ASM327640v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276485.1_ASM327648v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276485.1_ASM327648v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276485.1_ASM327648v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276485.1_ASM327648v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276485.1_ASM327648v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276495.1_ASM327649v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276495.1_ASM327649v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276495.1_ASM327649v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276495.1_ASM327649v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276495.1_ASM327649v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276505.1_ASM327650v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276505.1_ASM327650v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276505.1_ASM327650v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276505.1_ASM327650v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276505.1_ASM327650v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276515.1_ASM327651v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276515.1_ASM327651v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276515.1_ASM327651v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276515.1_ASM327651v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276515.1_ASM327651v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276605.1_ASM327660v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276605.1_ASM327660v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276605.1_ASM327660v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276605.1_ASM327660v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276605.1_ASM327660v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276615.1_ASM327661v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276615.1_ASM327661v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276615.1_ASM327661v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276615.1_ASM327661v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276615.1_ASM327661v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276625.1_ASM327662v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276625.1_ASM327662v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276625.1_ASM327662v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276625.1_ASM327662v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276625.1_ASM327662v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276635.1_ASM327663v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276635.1_ASM327663v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276635.1_ASM327663v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276635.1_ASM327663v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276635.1_ASM327663v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276725.1_ASM327672v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276725.1_ASM327672v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276725.1_ASM327672v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276725.1_ASM327672v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276725.1_ASM327672v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276735.1_ASM327673v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276735.1_ASM327673v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276735.1_ASM327673v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276735.1_ASM327673v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276735.1_ASM327673v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276745.1_ASM327674v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276745.1_ASM327674v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276745.1_ASM327674v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276745.1_ASM327674v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276745.1_ASM327674v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276755.1_ASM327675v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276755.1_ASM327675v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276755.1_ASM327675v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276755.1_ASM327675v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276755.1_ASM327675v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276845.1_ASM327684v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276845.1_ASM327684v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276845.1_ASM327684v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276845.1_ASM327684v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276845.1_ASM327684v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276855.1_ASM327685v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276855.1_ASM327685v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276855.1_ASM327685v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276855.1_ASM327685v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276855.1_ASM327685v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276865.1_ASM327686v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276865.1_ASM327686v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276865.1_ASM327686v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276865.1_ASM327686v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276865.1_ASM327686v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276875.1_ASM327687v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276875.1_ASM327687v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276875.1_ASM327687v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276875.1_ASM327687v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276875.1_ASM327687v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276965.1_ASM327696v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276965.1_ASM327696v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276965.1_ASM327696v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276965.1_ASM327696v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276965.1_ASM327696v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276975.1_ASM327697v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276975.1_ASM327697v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276975.1_ASM327697v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276975.1_ASM327697v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276975.1_ASM327697v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276985.1_ASM327698v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276985.1_ASM327698v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276985.1_ASM327698v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276985.1_ASM327698v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276985.1_ASM327698v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003276995.1_ASM327699v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003276995.1_ASM327699v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003276995.1_ASM327699v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003276995.1_ASM327699v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003276995.1_ASM327699v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277085.1_ASM327708v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277085.1_ASM327708v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277085.1_ASM327708v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277085.1_ASM327708v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277085.1_ASM327708v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277095.1_ASM327709v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277095.1_ASM327709v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277095.1_ASM327709v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277095.1_ASM327709v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277095.1_ASM327709v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277105.1_ASM327710v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277105.1_ASM327710v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277105.1_ASM327710v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277105.1_ASM327710v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277105.1_ASM327710v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277115.1_ASM327711v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277115.1_ASM327711v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277115.1_ASM327711v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277115.1_ASM327711v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277115.1_ASM327711v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277205.1_ASM327720v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277205.1_ASM327720v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277205.1_ASM327720v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277205.1_ASM327720v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277205.1_ASM327720v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277235.1_ASM327723v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277235.1_ASM327723v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277235.1_ASM327723v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277235.1_ASM327723v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277235.1_ASM327723v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277245.1_ASM327724v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277245.1_ASM327724v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277245.1_ASM327724v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277245.1_ASM327724v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277245.1_ASM327724v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277255.1_ASM327725v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277255.1_ASM327725v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277255.1_ASM327725v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277255.1_ASM327725v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277255.1_ASM327725v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277295.1_ASM327729v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277295.1_ASM327729v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277295.1_ASM327729v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277295.1_ASM327729v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277295.1_ASM327729v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277355.1_ASM327735v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277355.1_ASM327735v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277355.1_ASM327735v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277355.1_ASM327735v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277355.1_ASM327735v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277365.1_ASM327736v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277365.1_ASM327736v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277365.1_ASM327736v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277365.1_ASM327736v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277365.1_ASM327736v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277375.1_ASM327737v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277375.1_ASM327737v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277375.1_ASM327737v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277375.1_ASM327737v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277375.1_ASM327737v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277385.1_ASM327738v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277385.1_ASM327738v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277385.1_ASM327738v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277385.1_ASM327738v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277385.1_ASM327738v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277475.1_ASM327747v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277475.1_ASM327747v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277475.1_ASM327747v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277475.1_ASM327747v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277475.1_ASM327747v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277485.1_ASM327748v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277485.1_ASM327748v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277485.1_ASM327748v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277485.1_ASM327748v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277485.1_ASM327748v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277495.1_ASM327749v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277495.1_ASM327749v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277495.1_ASM327749v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277495.1_ASM327749v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277495.1_ASM327749v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277505.1_ASM327750v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277505.1_ASM327750v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277505.1_ASM327750v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277505.1_ASM327750v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277505.1_ASM327750v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277595.1_ASM327759v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277595.1_ASM327759v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277595.1_ASM327759v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277595.1_ASM327759v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277595.1_ASM327759v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277615.1_ASM327761v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277615.1_ASM327761v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277615.1_ASM327761v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277615.1_ASM327761v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277615.1_ASM327761v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277635.1_ASM327763v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277635.1_ASM327763v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277635.1_ASM327763v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277635.1_ASM327763v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277635.1_ASM327763v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277645.1_ASM327764v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277645.1_ASM327764v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277645.1_ASM327764v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277645.1_ASM327764v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277645.1_ASM327764v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277715.1_ASM327771v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277715.1_ASM327771v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277715.1_ASM327771v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277715.1_ASM327771v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277715.1_ASM327771v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277725.1_ASM327772v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277725.1_ASM327772v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277725.1_ASM327772v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277725.1_ASM327772v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277725.1_ASM327772v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277735.1_ASM327773v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277735.1_ASM327773v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277735.1_ASM327773v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277735.1_ASM327773v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277735.1_ASM327773v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277755.1_ASM327775v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277755.1_ASM327775v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277755.1_ASM327775v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277755.1_ASM327775v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277755.1_ASM327775v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277815.1_ASM327781v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277815.1_ASM327781v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277815.1_ASM327781v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277815.1_ASM327781v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277815.1_ASM327781v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277865.1_ASM327786v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277865.1_ASM327786v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277865.1_ASM327786v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277865.1_ASM327786v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277865.1_ASM327786v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003277875.1_ASM327787v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003277875.1_ASM327787v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003277875.1_ASM327787v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003277875.1_ASM327787v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003277875.1_ASM327787v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003326515.1_ASM332651v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003326515.1_ASM332651v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003326515.1_ASM332651v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003326515.1_ASM332651v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003326515.1_ASM332651v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003326525.1_ASM332652v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003326525.1_ASM332652v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003326525.1_ASM332652v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003326525.1_ASM332652v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003326525.1_ASM332652v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003326555.1_ASM332655v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003326555.1_ASM332655v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003326555.1_ASM332655v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003326555.1_ASM332655v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003326555.1_ASM332655v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_003326565.1_ASM332656v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_003326565.1_ASM332656v1_genomic.fna.out -m genome -e 1e-3
mv GCA_003326565.1_ASM332656v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_003326565.1_ASM332656v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_003326565.1_ASM332656v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864775.1_ASM1786477v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864775.1_ASM1786477v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864775.1_ASM1786477v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864775.1_ASM1786477v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864775.1_ASM1786477v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864785.1_ASM1786478v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864785.1_ASM1786478v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864785.1_ASM1786478v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864785.1_ASM1786478v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864785.1_ASM1786478v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864795.1_ASM1786479v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864795.1_ASM1786479v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864795.1_ASM1786479v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864795.1_ASM1786479v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864795.1_ASM1786479v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864805.1_ASM1786480v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864805.1_ASM1786480v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864805.1_ASM1786480v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864805.1_ASM1786480v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864805.1_ASM1786480v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864825.1_ASM1786482v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864825.1_ASM1786482v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864825.1_ASM1786482v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864825.1_ASM1786482v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864825.1_ASM1786482v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864915.1_ASM1786491v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864915.1_ASM1786491v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864915.1_ASM1786491v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864915.1_ASM1786491v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864915.1_ASM1786491v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864925.1_ASM1786492v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864925.1_ASM1786492v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864925.1_ASM1786492v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864925.1_ASM1786492v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864925.1_ASM1786492v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864945.1_ASM1786494v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864945.1_ASM1786494v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864945.1_ASM1786494v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864945.1_ASM1786494v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864945.1_ASM1786494v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017864975.1_ASM1786497v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017864975.1_ASM1786497v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017864975.1_ASM1786497v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017864975.1_ASM1786497v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017864975.1_ASM1786497v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865125.1_ASM1786512v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865125.1_ASM1786512v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865125.1_ASM1786512v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865125.1_ASM1786512v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865125.1_ASM1786512v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865225.1_ASM1786522v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865225.1_ASM1786522v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865225.1_ASM1786522v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865225.1_ASM1786522v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865225.1_ASM1786522v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865245.1_ASM1786524v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865245.1_ASM1786524v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865245.1_ASM1786524v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865245.1_ASM1786524v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865245.1_ASM1786524v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865275.1_ASM1786527v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865275.1_ASM1786527v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865275.1_ASM1786527v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865275.1_ASM1786527v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865275.1_ASM1786527v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865285.1_ASM1786528v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865285.1_ASM1786528v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865285.1_ASM1786528v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865285.1_ASM1786528v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865285.1_ASM1786528v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865375.1_ASM1786537v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865375.1_ASM1786537v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865375.1_ASM1786537v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865375.1_ASM1786537v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865375.1_ASM1786537v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865415.1_ASM1786541v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865415.1_ASM1786541v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865415.1_ASM1786541v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865415.1_ASM1786541v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865415.1_ASM1786541v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865475.1_ASM1786547v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865475.1_ASM1786547v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865475.1_ASM1786547v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865475.1_ASM1786547v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865475.1_ASM1786547v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865645.1_ASM1786564v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865645.1_ASM1786564v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865645.1_ASM1786564v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865645.1_ASM1786564v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865645.1_ASM1786564v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865655.1_ASM1786565v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865655.1_ASM1786565v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865655.1_ASM1786565v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865655.1_ASM1786565v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865655.1_ASM1786565v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865695.1_ASM1786569v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865695.1_ASM1786569v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865695.1_ASM1786569v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865695.1_ASM1786569v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865695.1_ASM1786569v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865715.1_ASM1786571v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865715.1_ASM1786571v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865715.1_ASM1786571v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865715.1_ASM1786571v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865715.1_ASM1786571v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865745.1_ASM1786574v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865745.1_ASM1786574v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865745.1_ASM1786574v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865745.1_ASM1786574v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865745.1_ASM1786574v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865875.1_ASM1786587v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865875.1_ASM1786587v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865875.1_ASM1786587v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865875.1_ASM1786587v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865875.1_ASM1786587v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017865995.1_ASM1786599v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017865995.1_ASM1786599v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017865995.1_ASM1786599v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017865995.1_ASM1786599v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017865995.1_ASM1786599v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866045.1_ASM1786604v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866045.1_ASM1786604v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866045.1_ASM1786604v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866045.1_ASM1786604v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866045.1_ASM1786604v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866055.1_ASM1786605v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866055.1_ASM1786605v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866055.1_ASM1786605v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866055.1_ASM1786605v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866055.1_ASM1786605v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866065.1_ASM1786606v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866065.1_ASM1786606v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866065.1_ASM1786606v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866065.1_ASM1786606v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866065.1_ASM1786606v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866075.1_ASM1786607v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866075.1_ASM1786607v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866075.1_ASM1786607v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866075.1_ASM1786607v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866075.1_ASM1786607v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866085.1_ASM1786608v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866085.1_ASM1786608v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866085.1_ASM1786608v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866085.1_ASM1786608v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866085.1_ASM1786608v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866135.1_ASM1786613v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866135.1_ASM1786613v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866135.1_ASM1786613v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866135.1_ASM1786613v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866135.1_ASM1786613v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866155.1_ASM1786615v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866155.1_ASM1786615v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866155.1_ASM1786615v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866155.1_ASM1786615v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866155.1_ASM1786615v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866175.1_ASM1786617v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866175.1_ASM1786617v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866175.1_ASM1786617v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866175.1_ASM1786617v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866175.1_ASM1786617v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866185.1_ASM1786618v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866185.1_ASM1786618v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866185.1_ASM1786618v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866185.1_ASM1786618v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866185.1_ASM1786618v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866255.1_ASM1786625v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866255.1_ASM1786625v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866255.1_ASM1786625v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866255.1_ASM1786625v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866255.1_ASM1786625v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866275.1_ASM1786627v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866275.1_ASM1786627v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866275.1_ASM1786627v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866275.1_ASM1786627v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866275.1_ASM1786627v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866305.1_ASM1786630v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866305.1_ASM1786630v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866305.1_ASM1786630v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866305.1_ASM1786630v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866305.1_ASM1786630v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866315.1_ASM1786631v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866315.1_ASM1786631v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866315.1_ASM1786631v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866315.1_ASM1786631v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866315.1_ASM1786631v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866385.1_ASM1786638v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866385.1_ASM1786638v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866385.1_ASM1786638v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866385.1_ASM1786638v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866385.1_ASM1786638v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866405.1_ASM1786640v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866405.1_ASM1786640v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866405.1_ASM1786640v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866405.1_ASM1786640v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866405.1_ASM1786640v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866455.1_ASM1786645v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866455.1_ASM1786645v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866455.1_ASM1786645v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866455.1_ASM1786645v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866455.1_ASM1786645v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866465.1_ASM1786646v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866465.1_ASM1786646v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866465.1_ASM1786646v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866465.1_ASM1786646v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866465.1_ASM1786646v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866475.1_ASM1786647v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866475.1_ASM1786647v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866475.1_ASM1786647v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866475.1_ASM1786647v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866475.1_ASM1786647v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866485.1_ASM1786648v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866485.1_ASM1786648v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866485.1_ASM1786648v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866485.1_ASM1786648v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866485.1_ASM1786648v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866495.1_ASM1786649v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866495.1_ASM1786649v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866495.1_ASM1786649v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866495.1_ASM1786649v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866495.1_ASM1786649v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866565.1_ASM1786656v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866565.1_ASM1786656v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866565.1_ASM1786656v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866565.1_ASM1786656v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866565.1_ASM1786656v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866575.1_ASM1786657v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866575.1_ASM1786657v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866575.1_ASM1786657v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866575.1_ASM1786657v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866575.1_ASM1786657v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866745.1_ASM1786674v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866745.1_ASM1786674v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866745.1_ASM1786674v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866745.1_ASM1786674v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866745.1_ASM1786674v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866755.1_ASM1786675v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866755.1_ASM1786675v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866755.1_ASM1786675v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866755.1_ASM1786675v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866755.1_ASM1786675v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866845.1_ASM1786684v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866845.1_ASM1786684v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866845.1_ASM1786684v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866845.1_ASM1786684v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866845.1_ASM1786684v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866945.1_ASM1786694v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866945.1_ASM1786694v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866945.1_ASM1786694v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866945.1_ASM1786694v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866945.1_ASM1786694v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866955.1_ASM1786695v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866955.1_ASM1786695v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866955.1_ASM1786695v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866955.1_ASM1786695v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866955.1_ASM1786695v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017866995.1_ASM1786699v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017866995.1_ASM1786699v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017866995.1_ASM1786699v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017866995.1_ASM1786699v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017866995.1_ASM1786699v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867005.1_ASM1786700v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867005.1_ASM1786700v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867005.1_ASM1786700v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867005.1_ASM1786700v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867005.1_ASM1786700v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867015.1_ASM1786701v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867015.1_ASM1786701v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867015.1_ASM1786701v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867015.1_ASM1786701v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867015.1_ASM1786701v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867205.1_ASM1786720v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867205.1_ASM1786720v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867205.1_ASM1786720v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867205.1_ASM1786720v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867205.1_ASM1786720v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867255.1_ASM1786725v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867255.1_ASM1786725v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867255.1_ASM1786725v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867255.1_ASM1786725v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867255.1_ASM1786725v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867285.1_ASM1786728v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867285.1_ASM1786728v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867285.1_ASM1786728v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867285.1_ASM1786728v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867285.1_ASM1786728v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867295.1_ASM1786729v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867295.1_ASM1786729v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867295.1_ASM1786729v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867295.1_ASM1786729v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867295.1_ASM1786729v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867305.1_ASM1786730v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867305.1_ASM1786730v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867305.1_ASM1786730v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867305.1_ASM1786730v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867305.1_ASM1786730v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867515.1_ASM1786751v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867515.1_ASM1786751v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867515.1_ASM1786751v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867515.1_ASM1786751v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867515.1_ASM1786751v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867525.1_ASM1786752v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867525.1_ASM1786752v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867525.1_ASM1786752v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867525.1_ASM1786752v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867525.1_ASM1786752v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867555.1_ASM1786755v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867555.1_ASM1786755v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867555.1_ASM1786755v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867555.1_ASM1786755v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867555.1_ASM1786755v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867585.1_ASM1786758v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867585.1_ASM1786758v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867585.1_ASM1786758v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867585.1_ASM1786758v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867585.1_ASM1786758v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867595.1_ASM1786759v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867595.1_ASM1786759v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867595.1_ASM1786759v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867595.1_ASM1786759v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867595.1_ASM1786759v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867735.1_ASM1786773v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867735.1_ASM1786773v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867735.1_ASM1786773v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867735.1_ASM1786773v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867735.1_ASM1786773v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867755.1_ASM1786775v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867755.1_ASM1786775v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867755.1_ASM1786775v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867755.1_ASM1786775v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867755.1_ASM1786775v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867835.1_ASM1786783v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867835.1_ASM1786783v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867835.1_ASM1786783v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867835.1_ASM1786783v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867835.1_ASM1786783v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867885.1_ASM1786788v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867885.1_ASM1786788v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867885.1_ASM1786788v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867885.1_ASM1786788v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867885.1_ASM1786788v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867915.1_ASM1786791v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867915.1_ASM1786791v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867915.1_ASM1786791v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867915.1_ASM1786791v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867915.1_ASM1786791v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867955.1_ASM1786795v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867955.1_ASM1786795v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867955.1_ASM1786795v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867955.1_ASM1786795v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867955.1_ASM1786795v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017867995.1_ASM1786799v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017867995.1_ASM1786799v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017867995.1_ASM1786799v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017867995.1_ASM1786799v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017867995.1_ASM1786799v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868025.1_ASM1786802v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868025.1_ASM1786802v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868025.1_ASM1786802v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868025.1_ASM1786802v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868025.1_ASM1786802v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868065.1_ASM1786806v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868065.1_ASM1786806v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868065.1_ASM1786806v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868065.1_ASM1786806v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868065.1_ASM1786806v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868075.1_ASM1786807v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868075.1_ASM1786807v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868075.1_ASM1786807v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868075.1_ASM1786807v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868075.1_ASM1786807v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868145.1_ASM1786814v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868145.1_ASM1786814v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868145.1_ASM1786814v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868145.1_ASM1786814v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868145.1_ASM1786814v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868155.1_ASM1786815v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868155.1_ASM1786815v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868155.1_ASM1786815v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868155.1_ASM1786815v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868155.1_ASM1786815v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868165.1_ASM1786816v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868165.1_ASM1786816v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868165.1_ASM1786816v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868165.1_ASM1786816v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868165.1_ASM1786816v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868175.1_ASM1786817v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868175.1_ASM1786817v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868175.1_ASM1786817v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868175.1_ASM1786817v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868175.1_ASM1786817v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868185.1_ASM1786818v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868185.1_ASM1786818v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868185.1_ASM1786818v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868185.1_ASM1786818v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868185.1_ASM1786818v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868345.1_ASM1786834v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868345.1_ASM1786834v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868345.1_ASM1786834v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868345.1_ASM1786834v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868345.1_ASM1786834v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868355.1_ASM1786835v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868355.1_ASM1786835v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868355.1_ASM1786835v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868355.1_ASM1786835v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868355.1_ASM1786835v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868485.1_ASM1786848v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868485.1_ASM1786848v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868485.1_ASM1786848v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868485.1_ASM1786848v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868485.1_ASM1786848v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868655.1_ASM1786865v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868655.1_ASM1786865v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868655.1_ASM1786865v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868655.1_ASM1786865v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868655.1_ASM1786865v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868835.1_ASM1786883v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868835.1_ASM1786883v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868835.1_ASM1786883v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868835.1_ASM1786883v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868835.1_ASM1786883v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868915.1_ASM1786891v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868915.1_ASM1786891v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868915.1_ASM1786891v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868915.1_ASM1786891v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868915.1_ASM1786891v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868925.1_ASM1786892v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868925.1_ASM1786892v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868925.1_ASM1786892v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868925.1_ASM1786892v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868925.1_ASM1786892v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017868975.1_ASM1786897v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017868975.1_ASM1786897v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017868975.1_ASM1786897v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017868975.1_ASM1786897v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017868975.1_ASM1786897v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869095.1_ASM1786909v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869095.1_ASM1786909v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869095.1_ASM1786909v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869095.1_ASM1786909v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869095.1_ASM1786909v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869105.1_ASM1786910v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869105.1_ASM1786910v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869105.1_ASM1786910v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869105.1_ASM1786910v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869105.1_ASM1786910v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869165.1_ASM1786916v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869165.1_ASM1786916v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869165.1_ASM1786916v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869165.1_ASM1786916v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869165.1_ASM1786916v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869305.1_ASM1786930v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869305.1_ASM1786930v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869305.1_ASM1786930v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869305.1_ASM1786930v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869305.1_ASM1786930v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869315.1_ASM1786931v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869315.1_ASM1786931v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869315.1_ASM1786931v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869315.1_ASM1786931v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869315.1_ASM1786931v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869325.1_ASM1786932v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869325.1_ASM1786932v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869325.1_ASM1786932v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869325.1_ASM1786932v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869325.1_ASM1786932v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869335.1_ASM1786933v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869335.1_ASM1786933v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869335.1_ASM1786933v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869335.1_ASM1786933v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869335.1_ASM1786933v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869405.1_ASM1786940v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869405.1_ASM1786940v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869405.1_ASM1786940v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869405.1_ASM1786940v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869405.1_ASM1786940v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869445.1_ASM1786944v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869445.1_ASM1786944v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869445.1_ASM1786944v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869445.1_ASM1786944v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869445.1_ASM1786944v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869495.1_ASM1786949v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869495.1_ASM1786949v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869495.1_ASM1786949v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869495.1_ASM1786949v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869495.1_ASM1786949v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869545.1_ASM1786954v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869545.1_ASM1786954v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869545.1_ASM1786954v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869545.1_ASM1786954v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869545.1_ASM1786954v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869555.1_ASM1786955v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869555.1_ASM1786955v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869555.1_ASM1786955v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869555.1_ASM1786955v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869555.1_ASM1786955v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869615.1_ASM1786961v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869615.1_ASM1786961v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869615.1_ASM1786961v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869615.1_ASM1786961v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869615.1_ASM1786961v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869625.1_ASM1786962v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869625.1_ASM1786962v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869625.1_ASM1786962v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869625.1_ASM1786962v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869625.1_ASM1786962v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869635.1_ASM1786963v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869635.1_ASM1786963v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869635.1_ASM1786963v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869635.1_ASM1786963v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869635.1_ASM1786963v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869715.1_ASM1786971v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869715.1_ASM1786971v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869715.1_ASM1786971v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869715.1_ASM1786971v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869715.1_ASM1786971v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869735.1_ASM1786973v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869735.1_ASM1786973v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869735.1_ASM1786973v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869735.1_ASM1786973v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869735.1_ASM1786973v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869745.1_ASM1786974v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869745.1_ASM1786974v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869745.1_ASM1786974v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869745.1_ASM1786974v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869745.1_ASM1786974v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869765.1_ASM1786976v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869765.1_ASM1786976v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869765.1_ASM1786976v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869765.1_ASM1786976v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869765.1_ASM1786976v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869795.1_ASM1786979v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869795.1_ASM1786979v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869795.1_ASM1786979v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869795.1_ASM1786979v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869795.1_ASM1786979v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869965.1_ASM1786996v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869965.1_ASM1786996v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869965.1_ASM1786996v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869965.1_ASM1786996v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869965.1_ASM1786996v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017869985.1_ASM1786998v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017869985.1_ASM1786998v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017869985.1_ASM1786998v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017869985.1_ASM1786998v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017869985.1_ASM1786998v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870085.1_ASM1787008v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870085.1_ASM1787008v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870085.1_ASM1787008v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870085.1_ASM1787008v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870085.1_ASM1787008v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870115.1_ASM1787011v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870115.1_ASM1787011v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870115.1_ASM1787011v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870115.1_ASM1787011v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870115.1_ASM1787011v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870125.1_ASM1787012v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870125.1_ASM1787012v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870125.1_ASM1787012v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870125.1_ASM1787012v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870125.1_ASM1787012v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870425.1_ASM1787042v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870425.1_ASM1787042v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870425.1_ASM1787042v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870425.1_ASM1787042v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870425.1_ASM1787042v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870435.1_ASM1787043v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870435.1_ASM1787043v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870435.1_ASM1787043v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870435.1_ASM1787043v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870435.1_ASM1787043v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870445.1_ASM1787044v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870445.1_ASM1787044v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870445.1_ASM1787044v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870445.1_ASM1787044v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870445.1_ASM1787044v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870455.1_ASM1787045v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870455.1_ASM1787045v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870455.1_ASM1787045v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870455.1_ASM1787045v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870455.1_ASM1787045v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870595.1_ASM1787059v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870595.1_ASM1787059v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870595.1_ASM1787059v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870595.1_ASM1787059v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870595.1_ASM1787059v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870715.1_ASM1787071v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870715.1_ASM1787071v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870715.1_ASM1787071v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870715.1_ASM1787071v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870715.1_ASM1787071v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870815.1_ASM1787081v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870815.1_ASM1787081v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870815.1_ASM1787081v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870815.1_ASM1787081v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870815.1_ASM1787081v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870925.1_ASM1787092v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870925.1_ASM1787092v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870925.1_ASM1787092v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870925.1_ASM1787092v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870925.1_ASM1787092v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870945.1_ASM1787094v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870945.1_ASM1787094v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870945.1_ASM1787094v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870945.1_ASM1787094v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870945.1_ASM1787094v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017870955.1_ASM1787095v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017870955.1_ASM1787095v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017870955.1_ASM1787095v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017870955.1_ASM1787095v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017870955.1_ASM1787095v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017871015.1_ASM1787101v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017871015.1_ASM1787101v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017871015.1_ASM1787101v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017871015.1_ASM1787101v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017871015.1_ASM1787101v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017871025.1_ASM1787102v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017871025.1_ASM1787102v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017871025.1_ASM1787102v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017871025.1_ASM1787102v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017871025.1_ASM1787102v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017871065.1_ASM1787106v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017871065.1_ASM1787106v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017871065.1_ASM1787106v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017871065.1_ASM1787106v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017871065.1_ASM1787106v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_017871085.1_ASM1787108v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_017871085.1_ASM1787108v1_genomic.fna.out -m genome -e 1e-3
mv GCA_017871085.1_ASM1787108v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_017871085.1_ASM1787108v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_017871085.1_ASM1787108v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393555.1_ASM1939355v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393555.1_ASM1939355v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393555.1_ASM1939355v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393555.1_ASM1939355v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393555.1_ASM1939355v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393615.1_ASM1939361v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393615.1_ASM1939361v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393615.1_ASM1939361v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393615.1_ASM1939361v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393615.1_ASM1939361v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393625.1_ASM1939362v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393625.1_ASM1939362v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393625.1_ASM1939362v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393625.1_ASM1939362v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393625.1_ASM1939362v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393785.1_ASM1939378v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393785.1_ASM1939378v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393785.1_ASM1939378v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393785.1_ASM1939378v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393785.1_ASM1939378v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393865.1_ASM1939386v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393865.1_ASM1939386v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393865.1_ASM1939386v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393865.1_ASM1939386v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393865.1_ASM1939386v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393935.1_ASM1939393v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393935.1_ASM1939393v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393935.1_ASM1939393v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393935.1_ASM1939393v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393935.1_ASM1939393v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393945.1_ASM1939394v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393945.1_ASM1939394v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393945.1_ASM1939394v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393945.1_ASM1939394v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393945.1_ASM1939394v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393955.1_ASM1939395v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393955.1_ASM1939395v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393955.1_ASM1939395v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393955.1_ASM1939395v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393955.1_ASM1939395v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019393975.1_ASM1939397v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019393975.1_ASM1939397v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019393975.1_ASM1939397v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019393975.1_ASM1939397v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019393975.1_ASM1939397v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394065.1_ASM1939406v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394065.1_ASM1939406v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394065.1_ASM1939406v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394065.1_ASM1939406v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394065.1_ASM1939406v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394075.1_ASM1939407v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394075.1_ASM1939407v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394075.1_ASM1939407v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394075.1_ASM1939407v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394075.1_ASM1939407v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394085.1_ASM1939408v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394085.1_ASM1939408v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394085.1_ASM1939408v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394085.1_ASM1939408v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394085.1_ASM1939408v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394095.1_ASM1939409v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394095.1_ASM1939409v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394095.1_ASM1939409v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394095.1_ASM1939409v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394095.1_ASM1939409v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394105.1_ASM1939410v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394105.1_ASM1939410v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394105.1_ASM1939410v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394105.1_ASM1939410v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394105.1_ASM1939410v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394265.1_ASM1939426v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394265.1_ASM1939426v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394265.1_ASM1939426v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394265.1_ASM1939426v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394265.1_ASM1939426v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394275.1_ASM1939427v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394275.1_ASM1939427v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394275.1_ASM1939427v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394275.1_ASM1939427v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394275.1_ASM1939427v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394305.1_ASM1939430v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394305.1_ASM1939430v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394305.1_ASM1939430v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394305.1_ASM1939430v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394305.1_ASM1939430v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394315.1_ASM1939431v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394315.1_ASM1939431v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394315.1_ASM1939431v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394315.1_ASM1939431v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394315.1_ASM1939431v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394325.1_ASM1939432v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394325.1_ASM1939432v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394325.1_ASM1939432v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394325.1_ASM1939432v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394325.1_ASM1939432v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394485.1_ASM1939448v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394485.1_ASM1939448v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394485.1_ASM1939448v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394485.1_ASM1939448v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394485.1_ASM1939448v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394495.1_ASM1939449v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394495.1_ASM1939449v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394495.1_ASM1939449v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394495.1_ASM1939449v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394495.1_ASM1939449v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394525.1_ASM1939452v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394525.1_ASM1939452v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394525.1_ASM1939452v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394525.1_ASM1939452v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394525.1_ASM1939452v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394535.1_ASM1939453v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394535.1_ASM1939453v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394535.1_ASM1939453v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394535.1_ASM1939453v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394535.1_ASM1939453v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394585.1_ASM1939458v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394585.1_ASM1939458v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394585.1_ASM1939458v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394585.1_ASM1939458v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394585.1_ASM1939458v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394645.1_ASM1939464v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394645.1_ASM1939464v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394645.1_ASM1939464v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394645.1_ASM1939464v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394645.1_ASM1939464v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394655.1_ASM1939465v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394655.1_ASM1939465v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394655.1_ASM1939465v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394655.1_ASM1939465v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394655.1_ASM1939465v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394665.1_ASM1939466v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394665.1_ASM1939466v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394665.1_ASM1939466v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394665.1_ASM1939466v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394665.1_ASM1939466v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394675.1_ASM1939467v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394675.1_ASM1939467v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394675.1_ASM1939467v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394675.1_ASM1939467v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394675.1_ASM1939467v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394725.1_ASM1939472v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394725.1_ASM1939472v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394725.1_ASM1939472v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394725.1_ASM1939472v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394725.1_ASM1939472v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394775.1_ASM1939477v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394775.1_ASM1939477v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394775.1_ASM1939477v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394775.1_ASM1939477v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394775.1_ASM1939477v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394785.1_ASM1939478v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394785.1_ASM1939478v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394785.1_ASM1939478v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394785.1_ASM1939478v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394785.1_ASM1939478v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394805.1_ASM1939480v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394805.1_ASM1939480v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394805.1_ASM1939480v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394805.1_ASM1939480v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394805.1_ASM1939480v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394815.1_ASM1939481v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394815.1_ASM1939481v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394815.1_ASM1939481v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394815.1_ASM1939481v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394815.1_ASM1939481v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394885.1_ASM1939488v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394885.1_ASM1939488v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394885.1_ASM1939488v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394885.1_ASM1939488v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394885.1_ASM1939488v1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394895.1_ASM1939489v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394895.1_ASM1939489v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394895.1_ASM1939489v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394895.1_ASM1939489v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394895.1_ASM1939489v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_019394915.1_ASM1939491v1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_019394915.1_ASM1939491v1_genomic.fna.out -m genome -e 1e-3
mv GCA_019394915.1_ASM1939491v1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_019394915.1_ASM1939491v1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_019394915.1_ASM1939491v1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099165.1_FS2D_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099165.1_FS2D_genomic.fna.out -m genome -e 1e-3
mv GCA_900099165.1_FS2D_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099165.1_FS2D_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099165.1_FS2D_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099195.1_20B2_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099195.1_20B2_genomic.fna.out -m genome -e 1e-3
mv GCA_900099195.1_20B2_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099195.1_20B2_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099195.1_20B2_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099255.1_TS12-A7_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099255.1_TS12-A7_genomic.fna.out -m genome -e 1e-3
mv GCA_900099255.1_TS12-A7_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099255.1_TS12-A7_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099255.1_TS12-A7_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099265.1_F25_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099265.1_F25_genomic.fna.out -m genome -e 1e-3
mv GCA_900099265.1_F25_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099265.1_F25_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099265.1_F25_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099275.1_22A4_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099275.1_22A4_genomic.fna.out -m genome -e 1e-3
mv GCA_900099275.1_22A4_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099275.1_22A4_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099275.1_22A4_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099285.1_K1-280-2B_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099285.1_K1-280-2B_genomic.fna.out -m genome -e 1e-3
mv GCA_900099285.1_K1-280-2B_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099285.1_K1-280-2B_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099285.1_K1-280-2B_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099305.1_P3-D5_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099305.1_P3-D5_genomic.fna.out -m genome -e 1e-3
mv GCA_900099305.1_P3-D5_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099305.1_P3-D5_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099305.1_P3-D5_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099345.1_D47-60-3B_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099345.1_D47-60-3B_genomic.fna.out -m genome -e 1e-3
mv GCA_900099345.1_D47-60-3B_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099345.1_D47-60-3B_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099345.1_D47-60-3B_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099355.1_1014-F5_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099355.1_1014-F5_genomic.fna.out -m genome -e 1e-3
mv GCA_900099355.1_1014-F5_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099355.1_1014-F5_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099355.1_1014-F5_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099375.1_GE7-4A_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099375.1_GE7-4A_genomic.fna.out -m genome -e 1e-3
mv GCA_900099375.1_GE7-4A_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099375.1_GE7-4A_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099375.1_GE7-4A_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099395.1_F12-3B_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099395.1_F12-3B_genomic.fna.out -m genome -e 1e-3
mv GCA_900099395.1_F12-3B_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099395.1_F12-3B_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099395.1_F12-3B_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099425.1_GUF54-A1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099425.1_GUF54-A1_genomic.fna.out -m genome -e 1e-3
mv GCA_900099425.1_GUF54-A1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099425.1_GUF54-A1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099425.1_GUF54-A1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900099445.1_6320-A7_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900099445.1_6320-A7_genomic.fna.out -m genome -e 1e-3
mv GCA_900099445.1_6320-A7_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900099445.1_6320-A7_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900099445.1_6320-A7_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177835.1_L560_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177835.1_L560_genomic.fna.out -m genome -e 1e-3
mv GCA_900177835.1_L560_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177835.1_L560_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177835.1_L560_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177845.1_L175_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177845.1_L175_genomic.fna.out -m genome -e 1e-3
mv GCA_900177845.1_L175_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177845.1_L175_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177845.1_L175_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177855.1_GMY321_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177855.1_GMY321_genomic.fna.out -m genome -e 1e-3
mv GCA_900177855.1_GMY321_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177855.1_GMY321_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177855.1_GMY321_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177865.1_L761_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177865.1_L761_genomic.fna.out -m genome -e 1e-3
mv GCA_900177865.1_L761_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177865.1_L761_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177865.1_L761_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177875.1_L328_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177875.1_L328_genomic.fna.out -m genome -e 1e-3
mv GCA_900177875.1_L328_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177875.1_L328_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177875.1_L328_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177885.1_L709_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177885.1_L709_genomic.fna.out -m genome -e 1e-3
mv GCA_900177885.1_L709_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177885.1_L709_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177885.1_L709_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177905.1_L564_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177905.1_L564_genomic.fna.out -m genome -e 1e-3
mv GCA_900177905.1_L564_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177905.1_L564_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177905.1_L564_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177915.1_L461_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177915.1_L461_genomic.fna.out -m genome -e 1e-3
mv GCA_900177915.1_L461_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177915.1_L461_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177915.1_L461_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177925.1_L716_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177925.1_L716_genomic.fna.out -m genome -e 1e-3
mv GCA_900177925.1_L716_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177925.1_L716_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177925.1_L716_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177935.1_L710_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177935.1_L710_genomic.fna.out -m genome -e 1e-3
mv GCA_900177935.1_L710_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177935.1_L710_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177935.1_L710_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177945.1_L174_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177945.1_L174_genomic.fna.out -m genome -e 1e-3
mv GCA_900177945.1_L174_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177945.1_L174_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177945.1_L174_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177975.1_L292_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177975.1_L292_genomic.fna.out -m genome -e 1e-3
mv GCA_900177975.1_L292_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177975.1_L292_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177975.1_L292_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177985.1_L712_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177985.1_L712_genomic.fna.out -m genome -e 1e-3
mv GCA_900177985.1_L712_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177985.1_L712_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177985.1_L712_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900177995.1_L723_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900177995.1_L723_genomic.fna.out -m genome -e 1e-3
mv GCA_900177995.1_L723_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900177995.1_L723_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900177995.1_L723_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178015.1_L725_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178015.1_L725_genomic.fna.out -m genome -e 1e-3
mv GCA_900178015.1_L725_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178015.1_L725_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178015.1_L725_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178025.1_L561_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178025.1_L561_genomic.fna.out -m genome -e 1e-3
mv GCA_900178025.1_L561_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178025.1_L561_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178025.1_L561_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178035.1_L528_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178035.1_L528_genomic.fna.out -m genome -e 1e-3
mv GCA_900178035.1_L528_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178035.1_L528_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178035.1_L528_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178045.1_SD085_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178045.1_SD085_genomic.fna.out -m genome -e 1e-3
mv GCA_900178045.1_SD085_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178045.1_SD085_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178045.1_SD085_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178065.1_L711_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178065.1_L711_genomic.fna.out -m genome -e 1e-3
mv GCA_900178065.1_L711_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178065.1_L711_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178065.1_L711_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178075.1_L181_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178075.1_L181_genomic.fna.out -m genome -e 1e-3
mv GCA_900178075.1_L181_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178075.1_L181_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178075.1_L181_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178105.1_BA_Souche789_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178105.1_BA_Souche789_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178105.1_BA_Souche789_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178105.1_BA_Souche789_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178105.1_BA_Souche789_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178125.1_BH_Souche796_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178125.1_BH_Souche796_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178125.1_BH_Souche796_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178125.1_BH_Souche796_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178125.1_BH_Souche796_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178135.1_CD_SoucheMTF2377_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178145.1_CA_SoucheMTF2344_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178155.1_BD_Souche2333_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178155.1_BD_Souche2333_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178155.1_BD_Souche2333_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178155.1_BD_Souche2333_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178155.1_BD_Souche2333_V1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178165.1_EH_SoucheMTF2432_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178175.1_CC_SoucheMTF2359_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178185.1_DC_SoucheMTF2247_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178195.1_EF_SoucheMTF2416_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178205.1_BF_Souche2251_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178205.1_BF_Souche2251_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178205.1_BF_Souche2251_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178205.1_BF_Souche2251_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178205.1_BF_Souche2251_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178215.1_FA_SoucheMTF2421_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178225.1_BG_Souche2331_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178225.1_BG_Souche2331_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178225.1_BG_Souche2331_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178225.1_BG_Souche2331_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178225.1_BG_Souche2331_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178235.1_CE_SoucheMTF2380_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178245.1_DD_SoucheMTF2388_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178255.1_EB_SoucheMTF2370_V1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178265.1_ED_SoucheMTF2387_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178275.1_BE_Souche2334_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178275.1_BE_Souche2334_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178275.1_BE_Souche2334_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178275.1_BE_Souche2334_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178275.1_BE_Souche2334_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178285.1_DG_SoucheMTF889_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178295.1_EG_SoucheMTF1798_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178305.1_DE_SoucheMTF2350_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178315.1_DB_SoucheMTF2249_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178325.1_FF_SoucheMTF2404_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178335.1_FD_SoucheMTF2419_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178345.1_FE_SoucheMTF2414_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178355.1_DF_SoucheMTF2351_V1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178365.1_EE_SoucheMTF2397_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178375.1_FC_SoucheMTF2418_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178385.1_FG_SoucheMTF2376_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178395.1_CB_SoucheMTF2357_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178405.1_FH_SoucheMTF2453_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178415.1_FI_SoucheMTF2264_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178425.1_EC_SoucheMTF2383_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178435.1_CG_SoucheMTF2412_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178445.1_EA_SoucheMTF2368_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178455.1_DH_SoucheMTF2352_V1_genomic.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178465.1_CH_SoucheMTF2413_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178475.1_FB_SoucheMTF2417_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178485.1_CF_SoucheMTF2386_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178495.1_DA_SoucheMTF2248_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/GCA_900178505.1_BC_Souche1152_V1_genomic.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o GCA_900178505.1_BC_Souche1152_V1_genomic.fna.out -m genome -e 1e-3
mv GCA_900178505.1_BC_Souche1152_V1_genomic.fna.out/short_summary.specific.saccharomycetes_odb10.GCA_900178505.1_BC_Souche1152_V1_genomic.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r GCA_900178505.1_BC_Souche1152_V1_genomic.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAA.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o MAA.re.fna.out -m genome -e 1e-3
mv MAA.re.fna.out/short_summary.specific.saccharomycetes_odb10.MAA.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r MAA.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAB.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o MAB.re.fna.out -m genome -e 1e-3
mv MAB.re.fna.out/short_summary.specific.saccharomycetes_odb10.MAB.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r MAB.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/MAL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o MAL.re.fna.out -m genome -e 1e-3
mv MAL.re.fna.out/short_summary.specific.saccharomycetes_odb10.MAL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r MAL.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAA.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAA.re.fna.out -m genome -e 1e-3
mv YAA.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAA.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAA.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAB.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAB.re.fna.out -m genome -e 1e-3
mv YAB.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAB.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAB.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAC.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAC.re.fna.out -m genome -e 1e-3
mv YAC.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAC.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAC.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAD.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAD.re.fna.out -m genome -e 1e-3
mv YAD.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAD.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAD.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAE.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAE.re.fna.out -m genome -e 1e-3
mv YAE.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAE.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAE.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAF.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAF.re.fna.out -m genome -e 1e-3
mv YAF.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAF.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAF.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAG.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAG.re.fna.out -m genome -e 1e-3
mv YAG.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAG.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAG.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAK.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAK.re.fna.out -m genome -e 1e-3
mv YAK.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAK.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAK.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAL.re.fna.out -m genome -e 1e-3
mv YAL.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAL.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAM.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAM.re.fna.out -m genome -e 1e-3
mv YAM.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAM.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAM.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAN.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAN.re.fna.out -m genome -e 1e-3
mv YAN.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAN.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAN.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAP.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAP.re.fna.out -m genome -e 1e-3
mv YAP.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAP.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAP.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAQ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAQ.re.fna.out -m genome -e 1e-3
mv YAQ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAQ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAQ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAR.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAR.re.fna.out -m genome -e 1e-3
mv YAR.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAR.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAR.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAT.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAT.re.fna.out -m genome -e 1e-3
mv YAT.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAT.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAT.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAU.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAU.re.fna.out -m genome -e 1e-3
mv YAU.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAU.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAU.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAV.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAV.re.fna.out -m genome -e 1e-3
mv YAV.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAV.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAV.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAW.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAW.re.fna.out -m genome -e 1e-3
mv YAW.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAW.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAW.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAX.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAX.re.fna.out -m genome -e 1e-3
mv YAX.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAX.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAX.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAY.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAY.re.fna.out -m genome -e 1e-3
mv YAY.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAY.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAY.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YAZ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YAZ.re.fna.out -m genome -e 1e-3
mv YAZ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YAZ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YAZ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBA.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBA.re.fna.out -m genome -e 1e-3
mv YBA.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBA.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBA.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBB.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBB.re.fna.out -m genome -e 1e-3
mv YBB.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBB.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBB.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBE.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBE.re.fna.out -m genome -e 1e-3
mv YBE.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBE.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBE.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBF.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBF.re.fna.out -m genome -e 1e-3
mv YBF.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBF.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBF.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBG.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBG.re.fna.out -m genome -e 1e-3
mv YBG.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBG.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBG.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBH.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBH.re.fna.out -m genome -e 1e-3
mv YBH.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBH.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBH.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBJ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBJ.re.fna.out -m genome -e 1e-3
mv YBJ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBJ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBJ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBK.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBK.re.fna.out -m genome -e 1e-3
mv YBK.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBK.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBK.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBL.re.fna.out -m genome -e 1e-3
mv YBL.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBL.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBM.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBM.re.fna.out -m genome -e 1e-3
mv YBM.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBM.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBM.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBN.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBN.re.fna.out -m genome -e 1e-3
mv YBN.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBN.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBN.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBO.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBO.re.fna.out -m genome -e 1e-3
mv YBO.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBO.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBO.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBP.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBP.re.fna.out -m genome -e 1e-3
mv YBP.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBP.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBP.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBQ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBQ.re.fna.out -m genome -e 1e-3
mv YBQ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBQ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBQ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBR.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBR.re.fna.out -m genome -e 1e-3
mv YBR.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBR.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBR.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBS.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBS.re.fna.out -m genome -e 1e-3
mv YBS.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBS.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBS.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBT.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBT.re.fna.out -m genome -e 1e-3
mv YBT.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBT.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBT.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBU.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBU.re.fna.out -m genome -e 1e-3
mv YBU.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBU.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBU.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBV.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBV.re.fna.out -m genome -e 1e-3
mv YBV.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBV.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBV.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBW.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBW.re.fna.out -m genome -e 1e-3
mv YBW.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBW.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBW.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBX.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBX.re.fna.out -m genome -e 1e-3
mv YBX.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBX.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBX.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBY.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBY.re.fna.out -m genome -e 1e-3
mv YBY.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBY.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBY.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YBZ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YBZ.re.fna.out -m genome -e 1e-3
mv YBZ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YBZ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YBZ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCA.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCA.re.fna.out -m genome -e 1e-3
mv YCA.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCA.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCA.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCB.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCB.re.fna.out -m genome -e 1e-3
mv YCB.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCB.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCB.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCC.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCC.re.fna.out -m genome -e 1e-3
mv YCC.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCC.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCC.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCD.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCD.re.fna.out -m genome -e 1e-3
mv YCD.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCD.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCD.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCE.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCE.re.fna.out -m genome -e 1e-3
mv YCE.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCE.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCE.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCF.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCF.re.fna.out -m genome -e 1e-3
mv YCF.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCF.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCF.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCG.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCG.re.fna.out -m genome -e 1e-3
mv YCG.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCG.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCG.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCH.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCH.re.fna.out -m genome -e 1e-3
mv YCH.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCH.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCH.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCI.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCI.re.fna.out -m genome -e 1e-3
mv YCI.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCI.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCI.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCJ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCJ.re.fna.out -m genome -e 1e-3
mv YCJ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCJ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCJ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCK.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCK.re.fna.out -m genome -e 1e-3
mv YCK.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCK.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCK.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCL.re.fna.out -m genome -e 1e-3
mv YCL.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCL.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCM.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCM.re.fna.out -m genome -e 1e-3
mv YCM.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCM.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCM.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCN.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCN.re.fna.out -m genome -e 1e-3
mv YCN.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCN.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCN.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCO.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCO.re.fna.out -m genome -e 1e-3
mv YCO.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCO.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCO.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCP.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCP.re.fna.out -m genome -e 1e-3
mv YCP.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCP.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCP.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCQ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCQ.re.fna.out -m genome -e 1e-3
mv YCQ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCQ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCQ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCR.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCR.re.fna.out -m genome -e 1e-3
mv YCR.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCR.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCR.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCS.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCS.re.fna.out -m genome -e 1e-3
mv YCS.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCS.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCS.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCT.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCT.re.fna.out -m genome -e 1e-3
mv YCT.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCT.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCT.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCU.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCU.re.fna.out -m genome -e 1e-3
mv YCU.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCU.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCU.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCV.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCV.re.fna.out -m genome -e 1e-3
mv YCV.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCV.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCV.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCW.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCW.re.fna.out -m genome -e 1e-3
mv YCW.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCW.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCW.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCX.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCX.re.fna.out -m genome -e 1e-3
mv YCX.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCX.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCX.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCY.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCY.re.fna.out -m genome -e 1e-3
mv YCY.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCY.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCY.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YCZ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YCZ.re.fna.out -m genome -e 1e-3
mv YCZ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YCZ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YCZ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDA.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDA.re.fna.out -m genome -e 1e-3
mv YDA.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDA.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDA.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDB.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDB.re.fna.out -m genome -e 1e-3
mv YDB.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDB.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDB.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDC.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDC.re.fna.out -m genome -e 1e-3
mv YDC.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDC.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDC.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDD.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDD.re.fna.out -m genome -e 1e-3
mv YDD.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDD.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDD.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDE.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDE.re.fna.out -m genome -e 1e-3
mv YDE.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDE.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDE.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDF.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDF.re.fna.out -m genome -e 1e-3
mv YDF.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDF.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDF.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDG.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDG.re.fna.out -m genome -e 1e-3
mv YDG.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDG.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDG.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDH.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDH.re.fna.out -m genome -e 1e-3
mv YDH.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDH.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDH.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDI.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDI.re.fna.out -m genome -e 1e-3
mv YDI.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDI.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDI.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDJ.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDJ.re.fna.out -m genome -e 1e-3
mv YDJ.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDJ.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDJ.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDK.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDK.re.fna.out -m genome -e 1e-3
mv YDK.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDK.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDK.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDL.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDL.re.fna.out -m genome -e 1e-3
mv YDL.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDL.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDL.re.fna.out)&
wait

(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDM.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDM.re.fna.out -m genome -e 1e-3
mv YDM.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDM.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDM.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDN.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDN.re.fna.out -m genome -e 1e-3
mv YDN.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDN.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDN.re.fna.out)&
(busco -i /lustre/home/acct-clslhz/clslhz/why_ssGEMs/data/assembled_genome/YDO.re.fna -l /lustre/home/acct-clslhz/clslhz/why_ssGEMs/test/busco/busco_downloads/lineages/saccharomycetes_odb10 -o YDO.re.fna.out -m genome -e 1e-3
mv YDO.re.fna.out/short_summary.specific.saccharomycetes_odb10.YDO.re.fna.out.txt 1900_assemblies_BUSCO_results/
rm -r YDO.re.fna.out)&
