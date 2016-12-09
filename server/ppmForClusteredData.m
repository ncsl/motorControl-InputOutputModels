% function [X, Y, beta, stats] = ppmForClusteredData(clusterData, outname, varargin)
% if isempty(outname)
%     fprintf('Please make sure to have a .mat filename to save the function output. Check function usage!');
% end

clusterData = load('trainData.mat');
clusterData = clusterData.trainData;
params = struct('intrinsic', [], 'extrinsic', [], 'ensemble', []);
params.intrinsic.self_q = 10;
params.intrinsic.sum_window = 10;
params.L = 100;
params.neuronIndex = 13;
params.clusterNum = 'Apc1'; % make sure this is correct for your cluster struct

params.ensemble.ensemble_q = 10;
params.ensemble.ensemble_aread_q = 0;

winAfterMove = 100;

%
conditions = {'mallet','push','pull','sphere'};
areas = {'S1','PMv','PMd','M1'};
clusterAreaIndices = struct('S1', [], 'PMv', [], 'PMd', [], 'M1', []);

%% extract the parameters 
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
ens_ar_order = params.ensemble.ensemble_q;
ens_area_ar = params.ensemble.ensemble_aread_q;

%- third extrinsic


%%- Set # of covariates per type of covariate
if ar_sum_size ~= 0
    self_avg_hist = round((100-ar_order) / ar_sum_size);
    num_int = self_avg_hist + ar_order;
else
    num_int = 0;
end
num_ens = ens_area_ar*length(areas);
num_ext = 0;
num_covar = 1 + num_int + num_ens + num_ext;

%% Initialize Matrices X/Y for speed
count = 1;
for iCond=1:length(conditions)
    trials = fieldnames(clusterData.(conditions{iCond}));
    conditionData = clusterData.(conditions{iCond});
    for iTrial=1:length(trials)
        trialData = conditionData.(trials{iTrial});
    
        % set start and end points based on led cue time and movement
        t_start = 501 - round(100*trialData.info.led);
        t_end = 501 + round(100*trialData.info.close);
        % takes time points only between led cue and stop time
        if (t_end - 501 > 100)
            t_end = 501 + 100;
        end
        if t_start < 100
            fprintf('Invalid start time. \n');
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
        cluster = fieldnames(trialData);
        
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
        
        % set start and end points based on led cue time and movement
        t_start = 501 - round(100*trialData.info.led);
        t_end = 501 + round(100*trialData.info.close);
        % takes time points only between led cue and stop time
        if (t_end - 501 > 100)
            t_end = 501 + 100;
        end
        if t_start < 100
            fprintf('Invalid start time. \n');
        end
        
        for iTime=t_start:t_end % loop through time windows
            covMatIdx = 0; % index through covariate matrix columnns
            
            intMat = [];
            ensMat = [];
            extMat = [];
            if num_int > 0
                %% Build Intrinsic Covariate Matrix
                intMat = nan(1, num_int);
                intMat(1) = iCond; % set the first variable to be condition
                covMatIdx = 1;

                %- construct self-history based on AR order
                self_hist = iTime - ar_order: iTime-1;
                intMat(covMatIdx+1:covMatIdx+length(self_hist)) = trialData.(clusterToAnalyze)(neuron_num, self_hist);
                
                covMatIdx = covMatIdx + ar_order; % account for begginin 
                
                %- construct self-sum-history based on AR order
                self_sum_hist = iTime - (ar_order*self_avg_hist)-ar_order : iTime - ar_order-1;
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

               ens_hist = iTime - ens_area_ar : iTime-1;
               %- construct spiking history of ensembles
               for iArea=1:length(fieldnames(clusterAreaIndices)) % loop through each neighbor
                   % get the area's neurons in this cluster
                   currIndices = clusterAreaIndices.(areas{iArea});
                   
                   % sum for all areas
                   neuronData = trialData.(clusterToAnalyze)(currIndices, ens_hist); 
                   neuronData = sum(neuronData, 1);
                   
                   ensMat(covMatIdx+1:covMatIdx+length(neuronData)) = neuronData;
                   
                   covMatIdx = covMatIdx + length(neuronData);
               end
            end
            
            if num_ext > 0
                %% Build Extrinsic Covariate Matrix
                extMat = nan(1, num_ext);
                covMatIdx = 0; % to keep track of where you are in ext covariate matrix
                if num_int == 0 && num_ens == 0
                    extMat(1) = iCond;
                    covMatIdx = 1;
                end
                
                
            end
            
            % build entire covMatrix for this time
            X(count,:) = [intMat, ensMat, extMat];
            
            % build observation Matrix
            Y(count) = trialData.(clusterToAnalyze)(neuron_num, iTime); 
            count = count+1;
        end % time loop
    end % trial loop
end % condition loop

[beta, dev, stats] = glmfit(X, Y, 'binomial');
figure;
subplot(211);
plot(stats.beta);
subplot(212);
plot(stats.p);
% [beta, dev, stats] = glmfit(X, Y, 'poisson');
% figure;
% subplot(211);
% plot(stats.beta);
% subplot(212);
% plot(stats.p);
% end