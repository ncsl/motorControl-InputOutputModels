function [ pcakine ] = pcaKine( kinematics )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [ pcakine ] = pcaKine( kinematics )
%-----------------------------------------------------------------------------------------
%
% Description:  Establishes kinematic space 
%
%-----------------------------------------------------------------------------------------
%   
%   Input:    kine       -   A structure  containing kinematic information from all trials within each
%                            condition, aligned -750ms to 750ms around movement epoch. Each trial 
%                            contains 3x1 cell of 30 x 300 matrix for x, y, z sensor coordinates
% 
%   Output:   pcakine    -   A structure  containing reduced dimension kinematic information from all trials within each
%                            condition, aligned -750ms to 750ms around movement epoch. Each condition contains a ntrial x 4
%                            cell, where all 30 sensor coordinates have been reduced to a 1x300 matrix for x, y, z.
%                            The fourth column contains the trial name string.
%                          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 08/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%----------------------------------------------------------------------------------------%
% 1. Reduce sensor dimensionality in consistent three dimensional space.
%----------------------------------------------------------------------------------------%
%

pcakine = struct('mallet',[],'pull',[],'push',[],'sphere',[]);
condition = fieldnames(kinematics); 
dim = cell([3 1]); 

tlength = 300;

% all x,y,z sensor information for all trials and all conditions concatinated. 
for i = 1 : length(condition) 
    trials = fieldnames(kinematics.(condition{i}));
    ntrials = length(trials); 

    allScore = cell([3 1]);

    
    for j = 1 : ntrials 
        tmp = kinematics.(condition{i}).(trials{j}); 
        for x = 1:3
            dim{x} = horzcat(dim{x}, tmp{x});
        end
    end
end    
  

% PCA for each three dimension done seperatly. 
% Only want to reduce sensor dimensionality
for x = 1:3
    [coeff, score,~,~,ex] = pca(double(transpose(dim{x})));   
    allScore{x} = transpose(score);
    
end

% reorganizes structure to provide condition and trial seperation for kinematic data
for i = 1 : length(condition) 
    trials = fieldnames(kinematics.(condition{i}));

    if i == 1
        p = 0; 
    else 
        p = p + (ntrials * tlength); 
    end
    
    ntrials = length(trials);

    for j = 1 : ntrials;
        check = [(j - 1) * tlength + 1 + p , j * tlength + p];
        for x = 1:3      
            pcakine.(condition{i}){j, x} = allScore{x}(1, check(1):check(2));
        end     
        pcakine.(condition{i}){j, x + 1} = trials{j};
    end   
end  

end

