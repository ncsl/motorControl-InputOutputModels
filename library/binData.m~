addpath('../');

load('../Monkey X/data_X0918.mat'); % load in the raw data file
clear lfp kinematics
binsize = 1; % in milliseconds
stepsize = binsize/1000; % in seconds

xmax = 0.5;                                         %1500ms after epoch
xmin = -0.5;                                        %500ms before epoch
nbins = round((xmax - xmin)*(1000 / binsize));      %1ms non-overlapping bins

% variables about the neurons
num_neurons = length(fieldnames(neurons));
neuron_names = fieldnames(neurons);
conditions = {'push', 'pull', 'mallet', 'sphere'};

binned_neurons = struct(); % data to be saved
%% For Each Neuron Process and Bin Data
for iNeuron=1:num_neurons
    current_neuron = neurons.(neuron_names{iNeuron})
    binned_neuron = struct('push', [], 'pull', [], 'mallet', [], 'sphere', []);
    
    
    %- loop through each condition
    for iCond=1:length(conditions)
        condition = current_neuron.(conditions{iCond});
        new_trials = struct();
        
        %- loop through each trial for this condition
        trials = fieldnames(condition);

        tic;
        for iTrial=1:length(trials)
            current_trial = condition.(trials{iTrial});
            
            % log event marker data to align data
            start_trial = current_trial.eventts(1);
            end_trial = current_trial.eventts(end);
            trialnum = current_trial.trialnum;
            
            %%- STORE ADDITIONAL MARKERS HERE IF YOU WANT TO LOG THEM
            instruction_index = find(current_trial.eventdata == 245, 1, 'last'); %,1,'last'
            move_index = find(current_trial.eventdata == 64);    % 
            close_index = find(current_trial.eventdata == 248, 1, 'last'); %248,1,'last'
            move_time = current_trial.eventts(move_index);
            instruction_time = current_trial.eventts(instruction_index);
            close_time = current_trial.eventts(close_index);
            
            % perform binning throughout entire trial
            current_Time = start_trial; % indexing through time
            index = 1;
            spike_times = current_trial.alldata;
            timeIndices = zeros(round((end_trial-start_trial)/(binsize/1000)),2);
            eventIndices = struct();
            spikes_binned = zeros(round((end_trial-start_trial)/(binsize/1000)),1);
%             tic;
            while current_Time < end_trial
                % increment in bin size
                num_spikes = length(find(spike_times > current_Time & spike_times <= (current_Time + binsize/1000)));
                
                % store data
                timeIndices(index,:) = [current_Time, current_Time+stepsize];
                spikes_binned(index) = num_spikes;
                try
                    if (current_Time-instruction_time <= stepsize && current_Time-instruction_time > 0)
                        eventIndices.instruction_time = index;
                    end
                catch e
                    eventIndices.instruction_time = e;
                end
                try
                    if (current_Time-move_time <= stepsize && current_Time-move_time > 0)
                        eventIndices.move_time = index;
                    end
                catch e
                    eventIndices.endhold_time = e;
                end
                try
                    if (current_Time-close_time <= stepsize && current_Time-close_time > 0)
                        eventIndices.endhold_time = index;
                    end
                catch e
                    eventIndices.endhold_time = e;
                end
                % increment time and indices
                current_Time = current_Time + stepsize;
                index = index+1;
            end
%             toc;
            
            % store data in trial struct
            new_trials.(trials{iTrial}) = struct();
            new_trials.(trials{iTrial}).spikeHist = spikes_binned;
            new_trials.(trials{iTrial}).eventIndices = eventIndices;
            new_trials.(trials{iTrial}).timeIndices = timeIndices;
            new_trials.(trials{iTrial}).trialnum = trialnum; 
%             new_trials(iTrial).spikeHist = spikes_binned;
%             new_trials(iTrial).eventIndices = eventIndices;
%             new_trials(iTrial).timeIndices = timeIndices;
%             new_trials(iTrial).trialnum = trialnum; 
            
            clear spikes_binned eventIndices timeIndices
        end
        toc;
        % store data in condition struct
        binned_neuron.(conditions{iCond}) = new_trials;
        clear new_trials
    end
%     binned_neurons.(neuron_names{iNeuron}) = binned_neuron;
    save(neuron_names{iNeuron}, 'binned_neuron');
    clear binned_neuron
end