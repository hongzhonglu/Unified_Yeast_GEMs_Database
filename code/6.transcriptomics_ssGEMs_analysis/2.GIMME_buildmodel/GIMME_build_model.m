function  GIMME_build_model(model_dir,transcriptomics_path,output_dir,threshold_fraction)
    
    % integrate transcriptomics data into ssGEMs
    % Data: Caudal, E. et al. (2023) doi:10.1101/2023.05.17.541122 
    % method: GIMME

    % add cobratoolbox to path
    % addpath(genpath('/dssg/home/acct-clslhz/clslhz/why_ssGEM/biosoft/cobratoolbox'));

    % initCobraToolbox(false)
    % changeCobraSolver ('glpk', 'all');
    changeCobraSolver('gurobi','all');


    % set working directory
    % cd 'code/6.transcriptomics_ssGEMs_analysis'


    % load transcriptomics data
    [num,txt] = xlsread(transcriptomics_path);
    % [num,txt] = xlsread('../output/sce969_transcriptome_tpmMatrix.csv');

    strainList = txt(1, (2:end));
    exp.gene = txt((2:end), 1);


    panModel= readCbModel('../../../model/panYeast.xml');
    
    strainList2=sort(strainList);
    % run GIMME for each strain
    for i = 1 : length(strainList2)
        strainName= strainList2{i};
        % search the index of the strainName in the strainList
        index = find(ismember(strainList,strainName));
        % check if the model exists
        if ~isfile(strcat(model_dir, '/', strainName, '.xml'))
            continue
        end
        disp(strainName);
        % check if the model is already integrated
        if isfile(strcat(output_dir, '/', strainName, '.xml'))
            continue
        end
        model = readCbModel(strcat(model_dir, '/', strainName, '.xml')); % stratgy1: use the ssGEMs as the input model
        % model = panModel; % stratgy2: use the panYeast as the input model
        model= SCmedium(model);
        optimizeCbModel(model,'max').f
        exp.rawValue = num(:,index);
        exp.value = num(:,index);
        [expressionRxns, parsedGPR] = mapExpressionToReactions(model, exp,true);
        % ignore reactions without genes by fill NaN with 1000
        expressionRxns(isnan(expressionRxns)) = 1000;
        threshold = quantile(expressionRxns(expressionRxns~=1000), 1-threshold_fraction);
        fprintf('threshold: %f\n', threshold);
        options.expressionRxns = expressionRxns;
        options.threshold = threshold;
        options.solver = 'GIMME';
        options.obj_frac = 0.8;
        integrated_model = createTissueSpecificModel(model, options);
        modelFileName = convertStringsToChars(strcat(output_dir,'/',strainName , '.xml'));
        gr=optimizeCbModel(integrated_model,'max').f;
        fprintf('%s\t%f\n', strainName, gr);
        % print the gene number and reaction number
        fprintf('rxn number: %d\n', length(integrated_model.rxns));
        fprintf('gene number: %d\n', length(integrated_model.genes));
        subsystem=integrated_model.subSystems;
        for i = 1:numel(subsystem)
            if iscell(subsystem{i}) && numel(subsystem{i}) > 1
                subsystem{i} = subsystem{i}{1};
            end
            if iscell(subsystem{i}) && numel(subsystem{i}) == 1
                subsystem{i} = char(subsystem{i});
            end
        end

        integrated_model.subSystems=subsystem;
        % disp(integrated_model.subSystems);

        % write the model
        %writeCbModel(integrated_model, 'mat', modelFileName);
        % writeCbModel(integrated_model, 'sbml', modelFileName);
        % check if gr >0
    %     if gr > 0
    %          writeCbModel(integrated_model, 'sbml', modelFileName);
    %          % print the strainName and gr value
    %          fprintf('%s\t%f\n', strainName, gr);
    %     end
    end
end

