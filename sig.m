function [neuronID] = sig(channel, letter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function: [neuronID] = sig(channel, letter)
%
%-----------------------------------------------------------------------------------------
%
% Description: establishes neuron from given channel and letter
%
%----------------------------------------------------------------------------------------- 
%   Input:    channel   -   An integer that represents channel associated with
%                           specific brain region of interest.
%
%             letter    -   a:z specifies new neuron on channel. Multiple neurons per channel, 
%                           but letters don't translate across channels.
% 
%   Output:   neuronID  -   A string that represents neuronID 
%                           ex: sig001a
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 04/11/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if channel < 10
    neuronID = sprintf('sig00%d%s',channel,letter);
elseif channel < 100
    neuronID = sprintf('sig0%d%s',channel,letter);
else
    neuronID = sprintf('sig%d%s',channel,letter);
end


