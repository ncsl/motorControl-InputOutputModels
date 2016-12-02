clear all
close all
clc

graphDir = './Monkey X/graphs/';
figDir = './figures/pearsonR/';
clusterDir = './Monkey X/clusters/';

% create the cluster data structure
binnedDataDir = './Monkey X/binned neurons 1ms/';
cluster_struct = struct('push', [], 'pull', [], 'mallet', [], 'sphere', []);
conditions = fieldnames(cluster_struct);
neurons = dir(strcat(binnedDataDir, '*.mat'));
neurons = {neurons.name};

load(fullfile(graphDir, 'pearsonRGraph_allconditions.mat'));
cluster_struct = struct();
[C, L, U] = SpectralClustering(graph, 4, 1);
for iClust=1:size(C,2) % loop through number of clusters
    clusterName = strcat('cluster', num2str(iClust));
    currentCluster = C(:, iClust); % get the cluster indices for current cluster
    neuronIndices = full(currentCluster == 1);
    neuronsInCluster = neurons(neuronIndices);

    for iNeuron=1:length(neuronsInCluster) % loop through each neuron in cluster
        neuronName = neuronsInCluster{iNeuron};
        neuronName = neuronName(1:end-4);
        currentNeuron = fullfile(binnedDataDir, neuronsInCluster{iNeuron});
        
        % load in binned neuron and assign it to the cluster structure
        load(currentNeuron);
        for iCond=1:length(conditions)
            %%- for a condition, load all neuron's activity across all trials
            condition = binned_neuron.(conditions{iCond}); % get the data for this current condition
            trials = fieldnames(condition);
        
            for iTrial=1:length(trials)
                % get the trial data struct
                trialData = condition.(trials{iTrial});

                if isfield(trialData.eventIndices, 'move_time')
                    instruction_index = trialData.eventIndices.instruction_time;
                    move_index = trialData.eventIndices.move_time;

                    flag = 1;
                    if (move_index - instruction_index)/1000 > 0.5 % move time was too slow
                        flag = 0;
                    end

                    if flag % only append spiking patterns for trials to include
                        cluster_struct.(clusterName).(neuronName).(conditions{iCond}).(trials{iTrial}) = ...
                            binned_neuron.(conditions{iCond}).(trials{iTrial}).spikeHist(move_index-500:move_index+499);
                    end
                end
            end % end of loop thru trials 
        end % end of loop throug conditions
    end % loop through each neuron in cluster
end  

% figure;
% for iCond=1:length(conditions)
%     currentGraph = graph.(conditions{iCond});
%     [C, L, U] = SpectralClustering(currentGraph, 4, 1);
%     
%     for iClust=1:size(C,2) % loop through number of clusters
%         clusterName = strcat('cluster', num2str(iClust));
%         currentCluster = C(:, iClust); % get the cluster indices for current cluster
%         neuronIndices = full(currentCluster == 1);
%         neuronsInCluster = neurons(neuronIndices);
%         
%         for iNeuron=1:length(neuronsInCluster) % loop through each neuron in cluster
%             % load in binned neuron and assign it to the cluster structure
%             load(fullfile(binnedDataDir, neuronsInCluster{iNeuron}));
%             cluster_struct.(conditions{iCond}).(clusterName) = binned_neuron.(conditions{iCond})
%         end
%     end
%     
%     subplot(2,2,iCond);
%     imagesc(full(C))
%     title(['Looking at ', conditions{iCond}]);
%     xlabel('Cluster Index');
%     ylabel('Counts');
% end

if ~exist(figDir, 'dir')
    mkdir(figDir);
end
if ~exist(clusterDir, 'dir')
    mkdir(clusterDir);
end

% save the cluster data struct
save(fullfile(clusterDir, 'pearsonRSpectralClustered_allconditions'), 'cluster_struct', '-v7.3');

% print(fullfile(figDir, 'unnormalizedwith4_allconditions'), '-dpng', '-r0')


