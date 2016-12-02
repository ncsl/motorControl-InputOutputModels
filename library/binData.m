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
            spike_times = current_trial.alldata;
            
            % log event marker data to align data
            start_trial = current_trial.eventts(1);
            end_trial = current_trial.eventts(end);
            trialnum = current_trial.trialnum;
            
            %%- STORE ADDITIONAL MARKERS HERE IF YOU WANT TO LOG THEM
<<<<<<< HEAD
            % indices of the markers
            instruction_index = find(current_trial.eventdata == 245,1,'last');
            move_index = find(current_trial.eventdata == 64);    % find move_indices
            if length(move_index) > 1
                move_index = move_index(2);
            end
            close_index = find(current_trial.eventdata == 248,1,'last'); %248,1,'last'
            
            % actual times of the markers
=======
            instruction_index = find(current_trial.eventdata == 245, 1, 'last'); %,1,'last'
            move_index = find(current_trial.eventdata == 64);    % 
            close_index = find(current_trial.eventdata == 248, 1, 'last'); %248,1,'last'
>>>>>>> 6920b92a9594267d2c614b4ae2c1cded666a2920
            move_time = current_trial.eventts(move_index);
            instruction_time = current_trial.eventts(instruction_index);
            close_time = current_trial.eventts(close_index);
            
            %%- time lock to MOVEMENT ONSET
            instruction_time = instruction_time - move_time;
            close_time = close_time - move_time;
            spike_times = spike_times - move_time;
            move_time = 0;
            
            flag = 0;
            % run a check on range of close and opening
            if abs(instruction_time) > 0.5% movement happens too slowly
                flag = 1; % ignore the trial
            end

            % only bin this trial if 
            if ~flag
                % keep spike_times within xmin and xmax
                spike_times = spike_times(spike_times>xmin & spike_times<xmax);
                
<<<<<<< HEAD
                 % initialize: perform binning throughout entire trial
                spikes_binned = zeros((xmax - xmin)/(binsize/1000), 1);
                
                % only keep spike_indices within xmax after move_time
                spike_indices = ceil(spike_times*(1000/binsize) + abs(xmin)*1000);
                spikes_binned(spike_indices) = 1;
                
                info.close_time = close_time;
                info.instruction_time = instruction_time;
                
                % store data in trial struct
                new_trials.(trials{iTrial}) = struct();
                new_trials.(trials{iTrial}).spikeHist = spikes_binned;
                new_trials.(trials{iTrial}).info = info;
                new_trials.(trials{iTrial}).trialnum = trialnum; 
=======
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
>>>>>>> 6920b92a9594267d2c614b4ae2c1cded666a2920
            end
            clear spikes_binned eventIndices timeIndices
        end
        toc;
        % store data in condition struct
        binned_neuron.(conditions{iCond}) = new_trials;
        clear new_trials
    end
%     binned_neurons.(neuron_names{iNeuron}) = binned_neuron;
%     binned_neuron.mallet
%     binned_neuron.pull
%     binned_neuron.push
%     binned_neuron.sphere

    save(neuron_names{iNeuron}, 'binned_neuron');
    clear binned_neuron
end