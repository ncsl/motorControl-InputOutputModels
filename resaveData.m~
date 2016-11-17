% script to resave data locally to make for easier data analysis
% By: Adam Li
% 11/11/16

% load in the data and set to certain monkey
neuron_names = load('neuron_names.mat');
neuron_names = neuron_names.X0918;
dataDir = './Monkey X/';

% directories for 4 different experimental conditions
pullDir = fullfile(dataDir, 'pull');    
pushDir = fullfile(dataDir, 'push');
malletDir = fullfile(dataDir, 'mallet');
sphereDir = fullfile(dataDir, 'sphere');

if ~exist(pullDir, 'dir')
    mkdir(pullDir);
    mkdir(pushDir);
    mkdir(malletDir);
    mkdir(sphereDir);
end

load('data_X0918.mat');
neuron_list = fieldnames(neurons);
for i=1:length(neuron_list)
    neuron = neurons.(neuron_list{i});
    
    push = neuron.push;
    mallet = neuron.mallet;
    sphere = neuron.sphere;
    pull = neuron.pull;
    
    
end

% loop through all neurons and save the corresponding extracted data
for i=1:length(neuron_names)
    [trains, kipos] = extractdata('X0918', neuron_names{i}, 'push');
    file = fullfile(pushDir, neuron_names{i});
    save(file, 'trains', 'kipos');
    
    [trains, kipos] = extractdata('X0918', neuron_names{i}, 'pull');
    file = fullfile(pullDir, neuron_names{i});
    save(file, 'trains', 'kipos');
    
    [trains, kipos] = extractdata('X0918', neuron_names{i}, 'mallet');
    file = fullfile(malletDir, neuron_names{i});
    save(file, 'trains', 'kipos');
    
    [trains, kipos] = extractdata('X0918', neuron_names{i}, 'sphere');
    file = fullfile(sphereDir, neuron_names{i});
    save(file, 'trains', 'kipos');
end

