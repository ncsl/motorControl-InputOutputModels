
global area
area = {'S1','PMv','PMd','M1'};

organizedData = mon_population();   % align spike data to movement onset 



%% Raw Data Raster

trial = 'trial50';
condition = 'mallet';

for i = 1:length(area)
    spike = organizedData.(condition).(trial).(area{i});
    spike = spike(~cellfun('isempty',spike));
    x = cellfun('size',spike, 2); 
    [~, Index] = sort(x, 'ascend');
    
    spikes{i,1} = spike(Index);
    
    
    % value = establishValue_mon(spikes{i,1});   %% dummy example for function
end
info = organizedData.(condition).(trial).info;


plotRaster_mon(spikes, info); 
suptitle(sprintf('%s %s %s', [area{:}], trial, condition))



%% Kinematics 

fig = figure(1);
set(fig, 'Position',[400 100 1300 900]); 


[kine] = kinemat(monkey);
[pcakine] = pcaKine(kine);

clf('reset');
hold on 
plotKine(pcakine);
suptitle(sprintf('%s', monkey));  
hold off 
