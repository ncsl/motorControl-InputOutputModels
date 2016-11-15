function [organizedData] = mon_population()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ex: [organizedData] = mon_population()

% function: [organizedData] = mon_population()
%
%-----------------------------------------------------------------------------------------
%
% Description:  Reorganizes the neurons structure from the master data_*****.mat
%               file to a new structure that organizes together neurons of the same trial under
%               an all conditions structure.  
%
%               Does not analyze trials where the difference between led
%               and movement onset and movement and closure are greater than 500 msec. 
%               NOTE: need MATLAB to be in ..\Data folder to run 'pwd' properly.
%
%-----------------------------------------------------------------------------------------
%   
%   Input:   area        -   A cell of strings that establish which region we are
%                            attempting to iterate through. It must be amongst "PMv", "PMd", 
%                            "S1", "M1", and "M1medial". 
%                            If using more than one electrode region than input needs to be cell array. 
%
% 
%   Output:  rganizedData  -  
%            

%                          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 11/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%----------------------------------------------------------------------------------------%
% 1. Check input variables
%----------------------------------------------------------------------------------------%
%
% Initialization
global area
 
% checks if area is a valid string of characters    
if ~iscell(area)  
    area = {area};
end

for i = 1:length(area)
    if (~ischar(area{i}) ||...
        (~strcmp(area{i},'PMv') && ~strcmp(area{i},'S1') && ~strcmp(area{i},'PMd') &&~strcmp(area{i},'M1') && ~strcmp(area{i},'M1medial')))
        error('Error: brain region did not have data collected');
    end 
end
    
% load the data
load(sprintf('%s/data_%s.mat',pwd,'X0918'));


%%
%----------------------------------------------------------------------------------------%
% 2. Extract the neural data
%----------------------------------------------------------------------------------------%
%
% initialize the output struct array and other required variables

organizedData = struct('mallet',[],'pull',[],'push',[],'sphere',[]); 


letters = 'a':'z';  
condition = fieldnames(organizedData); 
region = struct('PMv',[1,16],'S1',[17,32],'PMd',[33,48],'M1',[49,112],'M1medial',[113,128]);

countNeuron = zeros([length(area) 1]);

xmax = 0.5;                             %500ms after epoch
xmin = -0.5;                            %500ms before epoch


%channel corresponds with brain region recording
for num = 1:length(area)
    channel = region.(area{num});

    % iterates through all desired channels 
    for i = channel(1):channel(2)
        count = 1;
        while count > 0 
    % itereation through all letters for each channel       
            for letter = 1:length(letters)                                       
                n = letters(letter);

    % establishes neuron number
                neuro = sig(i,n);       % helper function                                        
                neuronID = sprintf('neurons.%s',neuro);

    % only iterates through conditions if neuron actually exists
                if isfield(neurons,neuro)
                    countNeuron(num) = countNeuron(num) + 1;
                    
                    for k = 1:length(condition) 
                        trials = fieldnames(eval(sprintf('%s.%s',neuronID,condition{k})));
                        ntrials = length(trials); 


                        for j = 1:ntrials 
                            currenttrial = char(trials(j));

                            neuron = eval(sprintf('%s.%s.%s', neuronID, condition{k}, currenttrial));
                            spikes = num2cell(neuron.alldata, 1);

                            [aligned, differences] = mon_alignData(spikes, neuron.eventts, neuron.eventdata);              % normalizes spike data to move epoch

                            if k == 1
                                if j == 81 || j == 82
                                    differences{2} = 0;     % error trials 
                                end
                            elseif k == 4
                                if j == 34 
                                    differences{2} = 0; 
                                end
                            end
                      
                            if differences{1} < 0.5 && differences{2} < 0.5 && differences{2} > 0.01        % only add trials where monkey responds faster than 0.5 sec, and closes between 0.01 and 0.5 seconds 
                                alignMatrix = cell2mat(aligned); 
                                alignMatrix(alignMatrix < xmin | alignMatrix > xmax)=  NaN;    % Move between -0.5s and 0.5s
                                align = alignMatrix(~isnan(alignMatrix(:)));
                              
                                if ~isempty(align);
                                    organizedData.(condition{k}).(currenttrial).(area{num}){countNeuron(num), 1} = align;    % ntrials x 1 cell of a neuron's firing times normalized to move 
                                    if ~isfield(organizedData.(condition{k}).(currenttrial), 'info');
                                        organizedData.(condition{k}).(currenttrial).info.led = differences{1};
                                        organizedData.(condition{k}).(currenttrial).info.close = differences{2};
                                    end
                                end
                                

                            end
                        end
                    end

                else 
                    count = 0; %will terminate while loop once for loop is broken
                    break; %breaks for loop if neuron does not exist
                end
            end
        end
    end
end   

end