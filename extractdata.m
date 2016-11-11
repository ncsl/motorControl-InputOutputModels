function [trains,kinpos] = extractdata(sessionID,neuronID,objID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function: [trains,kinpos] = extractdata(sessionID,neuronID,objID)
%
%-----------------------------------------------------------------------------------------
%
% Description: It extracts the spike trains collected from a single neuron across several
%              trials of a reach-out task and save them in a struct array. It extracts the
%              kinematic variables (i.e., [x,y,z] positions of the optical markers in the
%              Cartesian coordinates) associated with the spiking activity in each trial
%              and save them in a struct array.
%
%-----------------------------------------------------------------------------------------
%
% Input:    sessionID - A string that contains the name of the experimental session to be
%                       accessed. It must be one among "X0918", "X1002", and "Y0210".
%
%
%           neuronID  - A string that contains the name of the neuron whose spike trains
%                       must be extracted. A list of valid names is reported in the file
%                       "Description_Dataset.pdf".
%
%
%           objID     - A string that contains the name of the target object used in the
%                       reach-out task of interest. It must be one among "sphere", "push"
%                       (for the tasks involving the pushbutton), "pull" (for the tasks
%                       involving the cylinder), and "mallet".
%
%-----------------------------------------------------------------------------------------
%
% Output:   trains    - A struct array with N entries (one per trial). Each entry has the
%                       following fields:
%
%                          ** spikes  = 1xM array of spike times (in s) collected during
%                                       the current trial (between behavioral marker codes
%                                       251 and 253);
%
%                          ** avgrate = average spiking rate (in spikes/s) of the neuron
%                                       in the current trial (between behavioral marker
%                                       codes 251 and 253);
%
%                          ** markers = 1x6 array of time stamps (in s) that correspond to
%                                       the event markers: 
%                                         - Trial Onset (code: 251);
%                                         - Cue ON (code: 245);
%                                         - Movement Onset (code: 64);
%                                         - Switch Closed (code: 248);
%                                         - End of FH (code: 70);
%                                         - End of Trial (code: 253).
%
%
%           kinpos    - A struct array with N entries (one per trial). Each entry has the
%                       following fields:
%
%                          ** x       = NxT matrix of int32 numbers denoting the position
%                                       (in micron) along the x-axis in the Cartesian
%                                       coordinates of each of the N optical markers at
%                                       each of the T frames collected during the current
%                                       trial (between behavioral marker codes 251 and
%                                       253);
%
%                          ** y       = NxT matrix of int32 numbers denoting the position
%                                       (in micron) along the y-axis in the Cartesian
%                                       coordinates of each of the N optical markers at
%                                       each of the T frames collected during the current
%                                       trial;
%
%                          ** z       = NxT matrix of int32 numbers denoting the position
%                                       (in micron) along the z-axis in the Cartesian
%                                       coordinates of each of the N optical markers at
%                                       each of the T frames collected during the current
%                                       trial;
%
%                          ** time    = 1xT array of time stamps (in s) that correspond to
%                                       the T frames collected during the current trial;
%
%                          ** markers = 1x6 array of time stamps (in s) that correspond to
%                                       the event markers: 
%                                         - Trial Onset (code: 251);
%                                         - Cue ON (code: 245);
%                                         - Movement Onset (code: 64);
%                                         - Switch Closed (code: 248);
%                                         - End of FH (code: 70);
%                                         - End of Trial (code: 253).
%
%-----------------------------------------------------------------------------------------
%
% List of invoked MATLAB functions: NONE
%
%-----------------------------------------------------------------------------------------
%
% Author: S. Santaniello
%
% Ver.: 1.0 - Date: 03/24/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%----------------------------------------------------------------------------------------%
% 1. Check input variables
%----------------------------------------------------------------------------------------%
%
% initialization
if (nargin<3), error('Error: Not enough input'); end

% check if sessionID is a valid string of characters
if (~ischar(sessionID) ||...
        (~strcmp(sessionID,'X0918') && ~strcmp(sessionID,'X1002') && ~strcmp(sessionID,'Y0210')))
    error('Error: session does not exist');
else
    % load the data
    if strcmp(sessionID,'Y0210')
        load(sprintf('%s/Monkey Y/data_%s.mat',pwd,sessionID));
    else
        load(sprintf('%s/Monkey X/data_%s.mat',pwd,sessionID));
    end
end

% check if neuronID is a valid string of characters
if (~ischar(neuronID) || ~isfield(neurons,neuronID))
    error('Error: neuron does not exist');
end

% check if objID is a valid string of characters
if (~ischar(objID) ||...
        (~strcmp(objID,'sphere') && ~strcmp(objID,'mallet') && ~strcmp(objID,'pull') && ~strcmp(objID,'push')))
    error('Error: target object not encoded');
else
    % extract the number of trials
    ntrials = length(fieldnames(eval(sprintf('neurons.%s.%s',neuronID,objID))));
end


%----------------------------------------------------------------------------------------%
% 2. Extract the data and format in the new struct arrays
%----------------------------------------------------------------------------------------%
%
% initialize the output struct arrays
trains = struct('spikes',[],'markers',[],'avgrate',[]);
kinpos = struct('x',[],'y',[],'z',[],'time',[],'markers',[]);

% for each trial...
for i=1:ntrials
    
    % load the spike train, behavioral markers, and kinematic variables for the current
    % trial
    tmp1 = eval(sprintf('kinematics.%s.trial%d',objID,i));
    tmp2 = eval(sprintf('neurons.%s.%s.trial%d',neuronID,objID,i));
    if (tmp1.trialnum~=tmp2.trialnum)
        error('Error: trial numbers do not match');
    else
        % kinematic variables: x, y, z value for each optimal marker
        kinpos(i).x = reshape(tmp1.orig_data(1,:,:),size(tmp1.orig_data,2),size(tmp1.orig_data,3));
        kinpos(i).y = reshape(tmp1.orig_data(2,:,:),size(tmp1.orig_data,2),size(tmp1.orig_data,3));
        kinpos(i).z = reshape(tmp1.orig_data(3,:,:),size(tmp1.orig_data,2),size(tmp1.orig_data,3));
        
        % kinematic variables: time stamp of the behavioral markers and position values
        m0 = find(tmp1.eventdata==253,1,'last');
        m1 = find(tmp1.eventdata==70,1,'last');
        m2 = find(tmp1.eventdata==248,1,'last');
        m3 = find(tmp1.eventdata==245,1,'last');
        m4 = find(tmp1.eventdata(m3:m2)==64); 
        if (length(m4)<3), m4 = m4(end)+m3-1; else m4 = m4(2)+m3-1; end
        m5 = find(tmp1.eventdata==251,1,'first');
        
        kinpos(i).markers = tmp1.eventts([m5 m4 m3 m2 m1 m0])-tmp1.eventts(m5);
        kinpos(i).time = (1:size(tmp1.orig_data,3))./200;
        clear m0 m1 m2 m3 m4 m5
        
        % spike train: time stamp of the behavioral markers and spikes
        m0 = find(tmp2.eventdata==253,1,'last');
        m1 = find(tmp2.eventdata==70,1,'last');
        m2 = find(tmp2.eventdata==248,1,'last');
        m3 = find(tmp2.eventdata==245,1,'last');
        m4 = find(tmp2.eventdata(m3:m2)==64); 
        if (length(m4)<3), m4 = m4(end)+m3-1; else m4 = m4(2)+m3-1; end
        m5 = find(tmp2.eventdata==251,1,'first');
        
        trains(i).markers = tmp2.eventts([m5 m4 m3 m2 m1 m0])-tmp2.eventts(m5);
        if (~isempty(tmp2.alldata))
            trains(i).spikes = tmp2.alldata(:)-tmp2.eventts(m5);
            trains(i).avgrate = length(trains(i).spikes)/(tmp2.eventts(m0)-tmp2.eventts(m5));
        else
            trains(i).spikes = [];
            trains(i).avgrate = 0;
        end
        clear m0 m1 m2 m3 m4 m5
    end
    clear tmp1 tmp2
end
