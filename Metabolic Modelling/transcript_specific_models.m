%%% Code for setting the model's bounds according to the transcriptomic data
% by Giuseppe

% gamma = [2, 1.5, 1, 0.5];                            
% threshold = [0,10,25,50];                        
% lb = [-50,-30,-10,-5];   
% type = ["to_one","to_zero","binary"];
%clear, clc
%feasTol =1e-02;

load('human1.mat');
model = human1;
model = changeRxnBounds(model,'HMR_9034',-0.63,'l'); 
model = changeRxnBounds(model,'HMR_9034',-0.44,'u'); 

model = changeRxnBounds(model,'HMR_9135',0.54,'1'); 
model = changeRxnBounds(model,'HMR_9135',0.90,'u');

model = changeRxnBounds(model,'HMR_9063',-0.16,'1'); 
model = changeRxnBounds(model,'HMR_9063',-0.09,'u');

model = changeRxnBounds(model,'HMR_9048',-2,'1'); 
model = changeRxnBounds(model,'HMR_9048',0,'u');

% model.c(13350) = 0; %  remove biomass from being the objective if you want to optimise another function

% define parameters over which to iterate

gamma = [2];
threshold = [25];

% addpath(genpath('PATH_TO_COBRATOOLBOX'));     % add here the path to your cobratoolbox folder
% initCobraToolbox

[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);

% load gene expression and gene ids
load('gene_exp.mat');             % patients on the rows, gene expression along the columns (it is a matrix of only numbers, no ids)
load('gene_ids.mat');                   % gene ids of the genes present in the gene expression dataset, sorted by the same order in the gene_exp matrix, it is a string array
load('patient_ids.mat');              % patient id, sorted by the same order in the gene_exp matrix
%gene_exp=load('gene-exp.mat');
genes = model.genes;
genes_in_dataset = gene_ids; 
C{size(gene_exp, 1), length(model.rxns)} = [];

GeneExpressionArray = ones(numel(genes),1); 

k = 1;
tic
for g = 1:numel(gamma)                
    changeCobraSolver('gurobi', 'all');
	gam = gamma(g);
    new_k = main(threshold, gam, gene_exp,genes_in_dataset,patient_ids,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k);
    k = new_k;
end
toc

function new_k = main(threshold, gam, gene_exp,genes_in_dataset,patient_ids,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k)
    gamma = gam;
    for tr = 1:numel(threshold)
        prc_cut = threshold(tr);                                                            
        % cut data with respect to a threshold  
         data = gene_exp;
        %data = gene_exp./mean(gene_exp);                  % get the fold change (If A is a matrix, then mean(A) returns a row vector containing the mean of each column.)
        % data = data.*(data>prctile(data,prc_cut,1));    % if you want to binaruse data according to a percentile                                        
                
        % applying the bounds
        fprintf("Iteration (k): %d, Gamma: %d, Threshold: %d\n",k,gamma(g),threshold(tr));
        for t=1:size(data,1)          % in here we select a unique profile
           	expr_profile = data(t,:);
            pos_genes_in_datas
            et = zeros(numel(genes),1);% gene in the model human 1
            for i=1:length(pos_genes_in_dataset)
                position = find(strcmp(genes{i},genes_in_dataset),1); 
                if ~isempty(position)                                   
                    pos_genes_in_dataset(i) = position(1);              
                    GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));         
                end
            end
            if or(sum(isinf(GeneExpressionArray)) >= 1, sum(isnan(GeneExpressionArray)) >= 1)
                fprintf("\nError in the gene expression data!");dataframe
            end
            [fluxes] = transcriptomic_bounds(gamma(g), GeneExpressionArray, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length);
            %fprintf("\nNon-zero reactions: %d, sample number: %d, Solution status: %d\n", length(find(fluxes.v)), t, fluxes.stat);
            %fprintf("\nNon-zero reactions: %d, sample number: %d, Solution status: %d\n", length(find(fluxes.v)), t, fluxes.stat);
            %C(k,:) = num2cell([patient_ids(k), transpose(fluxes.v)]);
            C(k,:) = num2cell([patient_ids(k), transpose(fluxes)]);
            k = k +1;
        end
    end
    new_k = k; 
    T = cell2table(C);
    T.Properties.VariableNames{1} = 'patient_id';  
    writetable(T,"D:\Noushin.e\FBA\flux1.csv");     % set the folder where you want the data to be saved
end


%Find the index of a reaction, e.g., ‘ATPM', in the model
%rxnList = 'ATPM';
%rxnID = findRxnIDs(modelEcore, rxnList)
%Set ATP demand at the real rate now
%valATPM = 1.07;
%modelCore = changeRxnBounds(modelCore,'DM_atp_c_',valATPM,'l');
%modelCore = changeRxnBounds(modelCore,'DM_atp_c_',1000,'u');
%
%