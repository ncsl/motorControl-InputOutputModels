% function [X, Y, beta, stats] = ppmForClusteredData(clusterData, outname, varargin)
% if isempty(outname)
%     fprintf('Please make sure to have a .mat filename to save the function output. Check function usage!');
% end
ugh=5
clearvars  -except ugh
close all

clusterData = load('trainData.mat');
kineData = load('trainKine.mat'); 

clusterData = clusterData.trainData;
kineData = kineData.trainKine;

winAfterMove = 100;

%
conditions = {'mallet','push','pull','sphere'};
areas = {'S1','PMv','PMd','M1'};
clusterAreaIndices = struct('S1', [], 'PMv', [], 'PMd', [], 'M1', []);
otherClusterArea = struct('S1', [], 'PMv', [], 'PMd', [], 'M1', []);

%% extract the parameters 
% EXAMPLE PARAMS
% params = struct('intrinsic', [], 'extrinsic', [], 'ensemble', []);
% params.intrinsic.self_q = 10;
% params.intrinsic.sum_window = 10;
% params.L = 100;
% params.neuronIndex = 13;
% params.clusterNum = 'Apc1'; % make sure this is correct for your cluster struct
% 
% params.ensemble.ensemble_q = 10;
% params.ensemble.ensemble_aread_q = 10;
% 
% params.extrinsic.extrinsic_q = 5;

%params = load('test_params_adam'); 
params = load('test_params_sep'); 
numParam = fieldnames(params.test_params); 

global ugh
params = params.test_params.(numParam{ugh});

clusterToAnalyze = params.clusterNum;
dataLength = params.L;
neuron_num = params.neuronIndex; % neuron passed in set{1,...,16} for Apc1

% get # of trials, and # of neighbor neurons in cluster
trials = fieldnames(clusterData.mallet);
N = size(clusterData.mallet.(trials{1}).(clusterToAnalyze), 1) - 1;

nghbr_indices = 1:N+1;
nghbr_indices(neuron_num) = [];

%- first intrinsic
ar_order = params.intrinsic.self_q;
ar_sum_size = params.intrinsic.sum_window;

%- second ensemble
ens_area_ar = params.ensemble.ensemble_aread_q;

%- third extrinsic
ext_ar_order = params.extrinsic.extrinsic_q;

%- fourth other cluster ensemble
try
ens_cluster_index = params.ensemble.cluster_index;
ens_cluster_ar = params.ensemble.ensemble_cluster_q;
catch e
    ens_cluster_index = 2;
    ens_cluster_ar = 10;
end


%% Set # of covariates per type of covariate
if ar_sum_size ~= 0
    self_avg_hist = round((100-ar_order) / ar_sum_size);
    num_int = self_avg_hist + ar_order;
else
    self_avg_hist = 0; 
    num_int = 0;
end

if ens_area_ar ~= 0 
    num_ens = (dataLength/(ens_area_ar)) * length(areas); 
else
    num_ens = 0; 
end

num_ens_cluster = ens_cluster_ar*length(areas);

% set number of external covariates
num_ext = 4 * (ext_ar_order + 1);

%- total number of covariates
num_covar = 1 + num_int + num_ens + num_ext + num_ens_cluster;

%% Initialize Matrices X/Y for speed
count = 1;
for iCond=1:length(conditions)
    trials = fieldnames(clusterData.(conditions{iCond}));
    conditionData = clusterData.(conditions{iCond});
    for iTrial=1:length(trials)
        trialData = conditionData.(trials{iTrial});
    
        % set start and end points based on led cue time and movement
        t_start = 501 - round(1000*trialData.info.led);
        t_end = 501 + round(1000*trialData.info.close);
        % takes time points only between led cue and stop time
        if (t_end - 501 > 100)
            t_end = 501 + 100;
        end
        if t_start < 101
            t_start = 101; 
        end
        
        for iTime=t_start:t_end % loop through time windows
            count = count+1;
        end
    end
end
X = nan(count-1, num_covar);
Y = nan(count-1,1);

%% Build Covariate Matrix
count = 1; % for indexing through rows of matrix
for iCond=1:length(conditions) % loop through task-conditions
    trials = fieldnames(clusterData.(conditions{iCond}));
    conditionData = clusterData.(conditions{iCond});
    
    for iTrial=1:length(trials) % loop through trials        
        trialData = conditionData.(trials{iTrial});
        kinetrialData = kineData.(conditions{iCond}).(trials{iTrial});
        cluster = fieldnames(trialData);
        
        %% SET NEURONS IN CLUSTER AND ENSEMBLE CLUSTER
        % find neurons in each area
        clusterIndex = ~cellfun('isempty', (strfind(cluster, clusterToAnalyze)));
        clusterNeurons = trialData.info.areaM(:, clusterIndex);
        clusterNeurons = clusterNeurons(~cellfun('isempty', clusterNeurons));
        for iArea=1:length(areas)
            % find indices of neuron in each area
            area = areas{iArea};
            areaIndices = find(~cellfun('isempty', strfind(clusterNeurons, area)));
            if ismember(neuron_num, areaIndices)
                areaIndices(areaIndices==neuron_num) = [];
            end
            clusterAreaIndices.(area) = areaIndices;
        end
        
        % get neuron area data for other cluster
        otherClusters = ~cellfun('isempty', (strfind(cluster, cluster{ens_cluster_index})));
        otherNeurons = trialData.info.areaM(:, otherClusters);
        otherNeurons = otherNeurons(~cellfun('isempty', otherNeurons));
        for iArea=1:length(areas)
            % find indices of neuron in each area
            area = areas{iArea};
            areaIndices = find(~cellfun('isempty', strfind(otherNeurons, area)));
            otherClusterArea.(area) = areaIndices;
        end
        
        %% SET TIME 
        % set start and end points based on led cue time and movement
        t_start = 501 - round(1000*trialData.info.led);
        t_end = 501 + round(1000*trialData.info.close);
        % takes time points only between led cue and stop time
        if (t_end - 501 > 100)
            t_end = 501 + 100;
        end
        if t_start < 101
            t_start = 101; 
        end
        
        %establish kinematic matrix for trial
        [xPos, yPos, zPos] = kinetrialData.metric1{:}; 
        vel = kinetrialData.metric2{1}; 
        trialKine = [xPos, yPos, zPos, vel]; 
        
        for iTime=t_start:t_end % loop through time windows
            covMatIdx = 0; % index through covariate matrix columnns
            
            intMat = [];
            ensMat = [];
            extMat = [];
            ensClustMat = [];
            if num_int > 0
                %% Build Intrinsic Covariate Matrix
                intMat = nan(1, num_int);
                intMat(1) = iCond; % set the first variable to be condition
                covMatIdx = 1;

                %- construct self-history based on AR order
                self_hist = iTime - ar_order: iTime - 1;
                intMat(covMatIdx+1:covMatIdx+length(self_hist)) = trialData.(clusterToAnalyze)(neuron_num, self_hist);
                
                covMatIdx = covMatIdx + ar_order; % account for begginin 
                
                %- construct self-sum-history based on AR order
                self_sum_hist = iTime - 100 : iTime - ar_order -1;
                sum_hist = trialData.(clusterToAnalyze)(neuron_num, self_sum_hist);
                
                try
                    sum_hist = reshape(sum_hist, length(sum_hist) / self_avg_hist, self_avg_hist);
                catch e
                    disp(e);
                    disp(['Possibly window size is not able to be done...']);
                end
                sum_hist = sum(sum_hist, 1);
                intMat(covMatIdx+1:covMatIdx+self_avg_hist) = sum_hist;
                
                covMatIdx = covMatIdx + length(sum_hist);
            end
            
%             clear intMat
            
            if num_ens > 0
               %% Build Ensemble Covariate Matrix
               ensMat = nan(1, num_ens);
               covMatIdx = 0; % index through covariate matrix columnns
               
               if num_int == 0
                  ensMat(1) = iCond;
                  covMatIdx = 1; % index through covariate matrix columnns
               end

               %- construct spiking history of ensembles
               for iArea=1:length(fieldnames(clusterAreaIndices)) % loop through each neighbor
                   % get the area's neurons in this cluster
                   ens_hist = iTime - dataLength:iTime - 1;
                   currIndices = clusterAreaIndices.(areas{iArea});
                   
                   neuronData = trialData.(clusterToAnalyze)(currIndices, ens_hist); 
                   neuronData = sum(neuronData, 1);
                   
                   try
                       sum_hist = reshape(neuronData, ens_area_ar, length(neuronData) / ens_area_ar);
                   catch e
                       disp(e);
                       disp(['Possibly window size is not able to be done...']);
                   end
                   
                   sum_hist = sum(sum_hist, 1);
                   ensMat(covMatIdx+1 : covMatIdx+(dataLength/ens_area_ar)) = sum_hist;
                
                   covMatIdx = covMatIdx + length(sum_hist);
               end
            end
            
            %% Build Ensemble Covariate Matrix from Other Clusters
            if num_ens_cluster > 0
                ensClustMat = nan(1, num_ens_cluster);
                covMatIdx = 0;
                
                if num_int ==0 && num_ens ==0
                   ensClustMat(1) = iCond;
                   covMatIdx = 1;
                end
                
                %- construct spiking history of different cluster ensemble
                %- otherClusterArea -> contains information about cluster
                for iArea=1:length(fieldnames(otherClusterArea)) % loop through each neighbor
                   % get the area's neurons in this cluster
                   ens_hist = iTime - dataLength:iTime - 1;
                   currIndices = otherClusterArea.(areas{iArea});
                   
                   neuronData = trialData.(cluster{ens_cluster_index})(currIndices, ens_hist); 
                   neuronData = sum(neuronData, 1);
                   
                   try
                       sum_hist = reshape(neuronData, ens_cluster_ar, length(neuronData) / ens_cluster_ar);
                   catch e
                       disp(e);
                       disp(['Possibly window size is not able to be done...']);
                   end
                   
                   sum_hist = sum(sum_hist, 1);
                   ensClustMat(covMatIdx+1 : covMatIdx+(dataLength/ens_cluster_ar)) = sum_hist;
                
                   covMatIdx = covMatIdx + length(sum_hist);
               end
                
            end
            
            if num_ext > 0
                %% Build Extrinsic Covariate Matrix
                extMat = nan(1, num_ext);
                covMatIdx = 0;
                
                if num_int == 0 && num_ens == 0
                    extMat(1) = iCond;
                    covMatIdx = 1;
                end
                
                ext_hist = iTime - ext_ar_order : iTime;
                extMat(covMatIdx+1:covMatIdx+end) = [trialKine(ext_hist,1)', ...                    % [xPos(-5),xPos(-4),...xPos(0), yPos(), zPos(), vel()]
                    trialKine(ext_hist,2)',trialKine(ext_hist,3)',trialKine(ext_hist,4)'];

            end
            
            % build entire covMatrix for this time
            X(count,:) = [intMat, ensMat, ensClustMat, extMat];
            
            % build observation Matrix
            Y(count) = trialData.(clusterToAnalyze)(neuron_num, iTime); 
            count = count+1;
        end % time loop
    end % trial loop
end % condition loop

[beta, dev, stats] = glmfit(X, Y, 'poisson');
plotglm(dev, stats, [ar_order, self_avg_hist, num_ens, num_ens_cluster,num_ext]); 
plotks(X,Y,stats);

% figure;
% subplot(211);
% plot(stats.beta);
% subplot(212);
% plot(stats.p);
% [beta, dev, stats] = glmfit(X, Y, 'poisson');
% figure;
% subplot(211);
% plot(stats.beta);
% subplot(212);
% plot(stats.p);
% end