clear all
close all
clc

graphDir = './Monkey X/graphs/';
figDir = './figures/pearsonR/';
clusterDir = './Monkey X/clusters/';

% create the cluster data structure
binnedDataDir = './Monkey X/binned neurons 1ms/';
cluster_struct = struct('push', [], 'pull', [], 'mallet', [], 'sphere', []);

neurons = dir(strcat(binnedDataDir, '*.mat'));
neurons = {neurons.name};

load(fullfile(graphDir, 'pearsonRGraph.mat'));
conditions = fieldnames(graph);

figure;
for iCond=1:length(conditions)
    currentGraph = graph.(conditions{iCond});
    [C, L, U] = SpectralClustering(currentGraph, 4, 1);
    
    for iClust=1:size(C,2) % loop through number of clusters
        clusterName = strcat('cluster', num2str(iClust));
        currentCluster = C(:, iClust); % get the cluster indices for current cluster
        neuronIndices = full(currentCluster == 1);
        neuronsInCluster = neurons(neuronIndices);
        
        for iNeuron=1:length(neuronsInCluster) % loop through each neuron in cluster
            % load in binned neuron and assign it to the cluster structure
            load(fullfile(binnedDataDir, neuronsInCluster{iNeuron}));
            cluster_struct.(conditions{iCond}).(clusterName) = binned_neuron.(conditions{iCond})
        end
    end
    
    subplot(2,2,iCond);
    imagesc(full(C))
    title(['Looking at ', conditions{iCond}]);
    xlabel('Cluster Index');
    ylabel('Counts');
end

if ~exist(figDir, 'dir')
    mkdir(figDir);
end
if ~exist(clusterDir, 'dir')
    mkdir(clusterDir);
end

% save the cluster data struct
save(fullfile(clusterDir, 'pearsonRSpectralClustered'), 'cluster_struct');

print(fullfile(figDir, 'unnormalizedwith4'), '-dpng', '-r0')


