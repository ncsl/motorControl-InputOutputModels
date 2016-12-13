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
% cluster_struct = struct();
[C, L, U] = SpectralClustering(graph, 4, 1);

tempclusters = C;
clusters = zeros(length(tempclusters),1);
for i=1:length(tempclusters)
    clusterIndex = find(tempclusters(i,:)==1);
    clusters(i) = clusterIndex;
end

save('cindices.mat', 'clusters')
% initialize struct
for iCond=1:length(conditions)
    cluster_struct.(conditions{iCond}) = struct();
    load(fullfile(binnedDataDir, neurons{1}));
    trials = fieldnames(binned_neuron.(conditions{iCond}));
    for iNeuron=1:length(neurons)
        for iTrial=1:length(trials)
            if ~isfield( cluster_struct.(conditions{iCond}), trials{iTrial})
                cluster_struct.(conditions{iCond}).(trials{iTrial}) = struct();
            end
            
            for iClust=1:size(C,2)
                clustername = strcat('cluster', num2str(iClust));
                if ~isfield( cluster_struct.(conditions{iCond}).(trials{iTrial}), clustername)
                    cluster_struct.(conditions{iCond}).(trials{iTrial}).(clustername) = [];
                end
            end
        end
    end
end

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
                info = trialData.info;
                info.trialnum = trialData.trialnum;
                info.neuronsInCluster = neuronsInCluster;
                
                % add meta data to cluster condition trial
                cluster_struct.(conditions{iCond}).(trials{iTrial}).info = info;
                
                currentClusterSpikes = cluster_struct.(conditions{iCond}).(trials{iTrial}).(clusterName);
                if isempty(cluster_struct.(conditions{iCond}).(trials{iTrial}).(clusterName))
                    cluster_struct.(conditions{iCond}).(trials{iTrial}).(clusterName) = trialData.spikeHist;
                else
                    cluster_struct.(conditions{iCond}).(trials{iTrial}).(clusterName) = ...
                        cat(2, currentClusterSpikes, trialData.spikeHist);
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


