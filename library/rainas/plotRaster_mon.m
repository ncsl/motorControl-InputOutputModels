function [  ] = plotRaster_mon(spikes, info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ex: [  ] = plotRaster_mon(organizedData)
% function: [  ] = plotRaster_mon(organizedData)
%
%-----------------------------------------------------------------------------------------
%
% Description:  Plots raster plot from 'plotSpikeRaster.m' function
%               downloaded to Dropbox\Data\Raster.   
%
%-----------------------------------------------------------------------------------------
%   
%   Input:    organizedData
%
% 
%   Output:   Figure is the output.
%                                                          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 11/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%----------------------------------------------------------------------------------------%
% 1. Raster Plot 
%----------------------------------------------------------------------------------------%
%
global area

fig = figure(1);
set(fig, 'Position',[400 100 1300 900]);
clf('reset')
num = length(area); % 4 

for i = 1 : num
    subplot(2, 2, i);
    
    plotSpikeRaster(spikes{i},'PlotType','vertline');    %found in ..\Dropbox\Data\Raster
    title(area{i});
    if i == 1
        ylabel('Neuron #'); 
        xlabel('Time (s)');
    end

    hold on
    line([-info.led -info.led], ylim,'Color', 'b');       % led cue 
    line([0 0], ylim,'Color', 'g');                     % movement onset
    line([info.close info.close], ylim,'Color', 'm');   % switch close

end

