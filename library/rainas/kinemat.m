function [kine] = kinemat(sessionID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [kine] = kinematics(sessionID)
%-----------------------------------------------------------------------------------------
%
% Description:  Creates structure with only trials where monkey responds faster than 0.5 sec, 
%               and closes between 0.01 and 0.5 seconds
%
%-----------------------------------------------------------------------------------------
%   
%   Input:    sessionID  -   A string that contains the name of the experimental session to be
%                            accessed. It must be one among "X0918", "X1002", and "Y0210".
% 
%   Output:   kine       -   A structure  containing kinematic information from all trials within each
%                            condition, aligned -750ms to 750ms around movement epoch. Each trial 
%                            contains 3x1 cell of 30 x 300 matrix for x, y, z sensor coordinates
%                          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 07/11/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%----------------------------------------------------------------------------------------%
% 1. Check input variables
%----------------------------------------------------------------------------------------%
%
% Initialization

if (nargin<1), error('Error: Not enough input'); end

% check if sessionID is a valid string of characters
if (~ischar(sessionID) ||...
        (~strcmp(sessionID,'X0918') && ~strcmp(sessionID,'X1002') && ~strcmp(sessionID,'Y0210')))
    error('Error: session does not exist');   
end
    
% load the data
if strcmp(sessionID,'Y0210')
    load(sprintf('%s/Monkey Y/data_%s.mat',pwd,sessionID));
else
    load(sprintf('%s/Monkey X/data_%s.mat',pwd,sessionID));
end


%%
%----------------------------------------------------------------------------------------%
%2. Extract the kinematic data as cell array and format to the new struct array
%----------------------------------------------------------------------------------------%

kine = struct('mallet',[],'pull',[],'push',[],'sphere',[]); 
condition = fieldnames(kine); 
kinpos = struct('x',[],'y',[],'z',[]);

%adds cell array of kinematic data to trial structure
for i = 1:length(condition)
    ntrials = length(fieldnames(eval(sprintf('kinematics.%s',condition{i})))) - 3;

    for j = 1:ntrials           
                
        currenttrial = sprintf('trial%d',j);
        tmp1 = eval(sprintf('kinematics.%s.%s',condition{i},currenttrial));
        
        eventdata = tmp1.eventdata; 
        eventts = tmp1.eventts;
        
        led = find(eventdata == 245,1,'last');
        move = find(eventdata == 64);
        close = find(eventdata == 248,1,'last');
        
        if length(move) > 1
            move = move(2); % mallet also uses 64
        end
               
        difference{1} = eventts(move) - eventts(led);
        difference{2} = eventts(close) - eventts(move); 
        difference{3} = eventts(move) - eventts(1); 
        
        start = round((difference{3} - 0.75) * 200);
        stop = 300 + start;     % -750ms to 750ms around move epoch
        
        if start > 0 && difference{1} < 0.5 && difference{2} < 0.5 && difference{2} > 0.01  % only add trials where monkey responds faster than 0.5 sec, and closes between 0.01 and 0.5 seconds              
        
            % kinematic variables: x, y, z value for each optimal marker
            kinpos.x = reshape(tmp1.orig_data(1, :, start:stop - 1), size(tmp1.orig_data,2), []);
            kinpos.y = reshape(tmp1.orig_data(2, :, start:stop - 1), size(tmp1.orig_data,2), []);
            kinpos.z = reshape(tmp1.orig_data(3, :, start:stop - 1), size(tmp1.orig_data,2), []);

            kine.(condition{i}).(currenttrial) = struct2cell(kinpos);    
        end
    end
end


