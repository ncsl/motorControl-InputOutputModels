clear all
close all
clc

graphDir = './Monkey X/graphs/';
figDir = './figures/pearsonR/';

load(fullfile(graphDir, 'pearsonRGraph.mat'));
conditions = fieldnames(graph);

figure;
for iCond=1:length(conditions)
    currentGraph = graph.(conditions{iCond});
    [C, L, U] = SpectralClustering(currentGraph, 4, 1);
    
    subplot(2,2,iCond);
    imagesc(full(C))
    title(['Looking at ', conditions{iCond}]);
    xlabel('Cluster Index');
    ylabel('Counts');
end

if ~exist(figDir, 'dir')
    mkdir(figDir);
end

print(fullfile(figDir, 'unnormalizedwith4'), '-dpng', '-r0')