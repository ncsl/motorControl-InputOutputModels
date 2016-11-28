function [ alignmovespikes, difference ] = mon_alignData(spikes, eventts, eventdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function: [ alignmovespikes, difference ] = mon_alignData(spikes, eventts, eventdata)
%
%-----------------------------------------------------------------------------------------
%
% Description:  aligns Mx1 cell array data to epoch 'move'
%
%-----------------------------------------------------------------------------------------
%   
%   Input:    spikes            -   A Mx1 cell array where each cell contains 1xN
%                                   vector representing a single neuron's spike train. 
%                                   Represents population activity within a single trial.
%
%             eventts           -   A vector that indicates times of specific events within the trial.
%
%             eventdata         -   A vector that relates experimental events to trial times.
%
% 
%   Output:   alignmovespikes   -   A cell array that represents the aligned spike train data, 
%                                   aligned{1} to move.                        
%                                                          
%             difference        -   A cell array that represents time
%                                   between led cue and movement onset {1} and time between
%                                   movement onset and switch close {2}
%          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 11/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

led = find(eventdata == 245,1,'last');
move = find(eventdata == 64);
close = find(eventdata == 248,1,'last');

if length(move) > 1
    move = move(2); % mallet also uses 64
end

difference{1} = eventts(move) - eventts(led);
difference{2} = eventts(close) - eventts(move); 

alignmove = eventts(move);

alignmovespikes = cellfun(@(x) x - alignmove,spikes,'UniformOutput',false);

end

