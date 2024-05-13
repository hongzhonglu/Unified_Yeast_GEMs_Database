function [model] = SCmedium(model)

    % start with a clean slate: set all exchange reactions to upper bound = 1000
    % and lower bound = 0 (ie, unconstrained excretion, no uptake)
    [~,exchangeRxns] = getExchangeRxns(model,'out');
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;

    % set minimal medium
    desiredExchanges = {'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };   % potassium exchange

    blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange

    % set minimal nutrient uptake
    for i = 1:length(desiredExchanges)
        model.lb(find(strcmp(desiredExchanges{i}, model.rxns))) = -1000;
    end

    % block some rxns
    for i = 1:length(blockedExchanges)
        model.lb(find(strcmp(blockedExchanges{i}, model.rxnNames))) = 0;
        model.ub(find(strcmp(blockedExchanges{i}, model.rxnNames))) = 0;
    end

    % switch medium from minimal medium to SC medium which add 17 amino acids into medium
    model.lb(find(strcmp('L-alanine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-arginine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-asparagine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-aspartate exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-cysteine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-glutamate exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-glutamine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-glycine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-histidine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-isoleucine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-leucine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-lysine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-methionine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-phenylalanine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-proline exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-serine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-threonine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-tryptophan exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-tyrosine exchange', model.rxnNames))) = -0.01;
    model.lb(find(strcmp('L-valine exchange', model.rxnNames))) = -0.01;

    % set glucose uptake rate to 1
    model.lb(find(strcmp('D-glucose exchange', model.rxnNames))) = -1;

    end