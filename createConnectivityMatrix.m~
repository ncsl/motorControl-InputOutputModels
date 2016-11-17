% data = load('organizedData.mat');
% data = data.organizedData;
% 
% mallet = data.mallet;
% push = data.push;
% sphere = data.sphere;
% pull = data.pull;
% 
% experiment = mallet;
% trials = fieldnames(experiment);
% % loop through all trials
% for i=1:length(trials)
%     trialname = trials{i};
%     trial = experiment.(trialname);
%     
%     S1 = trial.S1;
%     Pmv = trial.PMv;
%     Pmd = trial.PMd;
%     M1 = trial.M1;
%     
% end

load('data_X0918.mat');
neuron_list = fieldnames(neurons);

organizedData = struct('mallet',[],'pull',[],'push',[],'sphere',[]); 
letters = 'a':'c';  
condition = fieldnames(organizedData); 
region = struct('PMv',[1,16],'S1',[17,32],'PMd',[33,48],'M1',[49,112],'M1medial',[113,128]);

xmax = 0.5;                             %500ms after epoch
xmin = -0.5;                            %500ms before epoch

areas = {'S1','PMv','PMd','M1'};

for iA=1:length(areas) % looping through each brain recording region
    channel = region.(areas{iA});
    
    for iChan=channel(1):channel(2) % loop thruogh range of desired channels
         count=1;
         while count > 0
             for letter=1:length(letters)
                 n = letters(letter);
                 
                 neuro = sig(iChan, n); % get neuron in sig<channel><letter> form
                 neuronID = sprintf('neurons.%s',neuro);
                 if isfield(neurons,neuro)
                    countNeuron(num) = countNeuron(num) + 1;
                    
                    for k = 1:length(condition) 
                        trials = fieldnames(eval(sprintf('%s.%s',neuronID,condition{k})));
                        ntrials = length(trials); 

                    end
                 end
             end
         end
    end
end

for i=1:length(neuron_list)
    neuron = neurons.(neuron_list{i});
    
    push = neuron.push;
    mallet = neuron.mallet;
    sphere = neuron.sphere;
    pull = neuron.pull;
    
    
end