% script to create connectivity across all neurons regardless of brain
% region

clear all
clc
close all

addpath('./library/');

% data files we want to examine
dataDir = 'Monkey X/binned neurons 1ms/';
graphDir = fullfile('Monkey X/graphs/');
if ~exist(graphDir, 'dir')
    mkdir(graphDir);
end

matFiles = dir(strcat(dataDir, '*.mat'));
matFiles = {matFiles.name};

areas = {'S1','PMv','PMd','M1'};

organizedData = struct('mallet',[],'pull',[],'push',[],'sphere',[]); 
letters = 'a':'z';  
conditions = fieldnames(organizedData); 
regions = struct('PMv',[1,16],'S1',[17,32],'PMd',[33,48],'M1',[49,112],'M1medial',[113,128]);
regionNeurons = struct('PMv', [], 'S1', [], 'PMd', [], 'M1', []);%, 'M1medial', []);

countNeuron = zeros([length(areas) 1]);
trialsToInclude = cell(length(matFiles), 120);
%%- 1. Connect all neurons
for iNeuron=1:length(matFiles)
    currentNeuron = fullfile(dataDir, matFiles{iNeuron});
    neuronName = matFiles{iNeuron};
    neuronName = neuronName(1:end-4);
    load(currentNeuron);
    
    for iCond=1:length(conditions) % loop through all conditions
        %%- for a condition, load all neuron's activity across all trials
        condition = binned_neuron.(conditions{iCond}); % get the data for this current condition
        trials = fieldnames(condition);
        spiking_vector = [];
        for iTrial=1:length(trials)
            % get the trial data struct
            trialData = condition.(trials{iTrial});
            spikeHist = trialData.spikeHist;

            % get the entire spiking parttern for entire trial
            if isempty(spiking_vector)
                spiking_vector = spikeHist;
            else
                spiking_vector = cat(1, spiking_vector, spikeHist);
            end
        end

        organizedData.(conditions{iCond}).(neuronName) = spiking_vector;
    end    
end

%%- create graph irrespective of condition
data = [];
neurons = fieldnames(organizedData.mallet);
for iNeuron=1:length(neurons)
    neuronData = [];
    for iCond=1:length(conditions)
        if isempty(neuronData)
            neuronData = organizedData.(conditions{iCond}).(neurons{iNeuron});
        else
            neuronData = cat(1, neuronData, organizedData.(conditions{iCond}).(neurons{iNeuron}));
        end
    end
    if isempty(data)
        data = neuronData;
    else
        data = cat(2, data, neuronData);
    end
end
graph = corr(data);
save(fullfile(graphDir, 'pearsonRGraph_allconditions'), 'graph');

%%- create graph based on pearson R correlation
% graph = struct('mallet',[],'pull',[],'push',[],'sphere',[]);
% for iCond=1:length(conditions)
%     data = [];
%     
%     conditionData = organizedData.(conditions{iCond});
%     neurons = fieldnames(conditionData);
%     for iNeuron=1:length(neurons)
%         neuronData = conditionData.(neurons{iNeuron});
%         if isempty(data)
%             data = neuronData;
%         else
%             data = cat(2, data, neuronData);
%         end
%     end
%     
%     graph.(conditions{iCond}) = corr(data);
% end
% 
% save(fullfile(graphDir, 'pearsonRGraph'), 'graph');


%%- 2. Connect by areas
% loop through each file
% channel corresponds with brain region recording
% for num = 1:length(areas)
%     channel_range = regions.(areas{num});
% 
%     regionNeurons.(areas{num}) = {};
%     % iterates through all desired channels 
%     for i = channel_range(1):channel_range(2)
%         flag = 1;
%         while flag > 0 
%             % itereation through all letters for each channel       
%             for letter = 1:length(letters)                                       
%                 n = letters(letter);
% 
%                 % establishes neuron number
%                 neuro = sig(i,n);       % helper function 
%                 filetoload = fullfile(dataDir, strcat(neuro, '.mat'));
%                 
%                 if exist(filetoload, 'file') % only iterates through conditions if neuron actually exists
% %                     load(filetoload);                        % load neuron struct
%                     countNeuron(num) = countNeuron(num) + 1; % keep track of number of neurons in each region
%                     
%                     regionNeurons.(areas{num}){end+1} = neuro;
%                 else 
%                     flag = 0; %will terminate while loop once for loop is broken
%                     break; %breaks for loop if neuron does not exist
%                 end
%                 
%             end
%         end
%     end
% end   
% 
% %%- 2a. now loop through each region and construct connectivity graph
% for num=1:length(areas)
%     neurons = regionNeurons.(areas{num});
%     
%     for iNeuron=1:length(neurons)
%         filetoload = fullfile(dataDir, strcat(neurons{iNeuron}, '.mat'));
%         load(filetoload);        
%         
%         for iCond=1:length(conditions) % loop through all conditions
%             %%- for a condition, load all neuron's activity across all trials
%             condition = binned_neuron.(conditions{iCond}); % get the data for this current condition
%             trials = fieldnames(condition);
%             spiking_vector = [];
%             for iTrial=1:length(trials)
%                 % get the entire spiking parttern for entire trial
%                 trialData = condition.(trials{iTrial}).spikeHist; 
%                 if isempty(spiking_vector)
%                     spiking_vector = trialData;
%                 else
%                     spiking_vector = cat(1, spiking_vector, trialData);
%                 end
%             end
% 
%             organizedData.(conditions{iCond}).(areas{num}).(neurons{iNeuron}) = spiking_vector;
%         end    
%     end
% end