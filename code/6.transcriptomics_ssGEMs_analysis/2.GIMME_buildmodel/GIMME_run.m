% use GIMME method to integrate transcriptomics data into ssGEMs

model_dir='../output/gapfilled_ssGEMs';

% use absolute transcriptomics data
transcriptomics_path='../output/sce969_transcriptome_tpmMatrix.xlsx';
output_dir='../output/GIMME_0.75_ssGEMs';
threshold_fraction=0.75;   % Values corresponding to the 75% percentile in the reaction expression scores were applied as thresholds
GIMME_build_model(model_dir,transcriptomics_path,output_dir,threshold_fraction);


% use relative transcriptomics data(foldchange(log2))
% transcriptomics_path='output/sce969_transcriptome_foldchange.xlsx'
% output_dir='output/GIMME_fc_ssGEMs'
% threshold=0.5

% GIMME_build_model(model_dir,transcriptomics_path,output_dir,threshold)



