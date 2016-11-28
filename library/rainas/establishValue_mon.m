function [ valuemove ] = establishValue_mon(allmove)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[value] = establishValue(allmove)
% function: [ value ] = establishValue(allmove)
%
%-----------------------------------------------------------------------------------------
%
% Description:  Creates matrix of hist for all neurons trials. 
%
%-----------------------------------------------------------------------------------------
%   
%   Input:      allmove    -   ntrials x 1 cell containing all of a single neuron's
%                              firing times for a specific conditions
%                              aligned to movement onset of each trial.
%                              
%
% 
%   Output:     valuemove      -   1 x N matrix containing spike counts in 1ms time bins.
%                              aligned to movement onset
%                         
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 06/20/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

xmax = 0.5;                             %1500ms after epoch
xmin = -0.5;                            %500ms before epoch
nbins = round((xmax - xmin)*1000);      %1ms non-overlapping bins
%nbins = round((xmax - xmin)*100);       %10ms non-overlapping bins

allmove{length(allmove)+1} = [xmin, xmax]; % Move between -0.75s and 1.25s


nx = max(cellfun(@numel,allmove));                                                     % determines max number of neural activity of the cell
spikeMatrixmove = cell2mat(cellfun(@(x)[x,nan(1,nx-numel(x))],allmove,'uni',false)); 

valuemove = histcounts(spikeMatrixmove, nbins);         % Creates matrix of neurons 

valuemove(1) = valuemove(1) - 1; 
valuemove(length(valuemove)) = valuemove(length(valuemove)) - 1;

end

