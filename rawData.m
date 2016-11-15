
global area
area = {'S1','PMv','PMd','M1'};

organizedData = mon_population();   % align spike data to movement onset 



%% 

trial = 'trial50';
condition = 'mallet';

for i = 1:length(area)
    spike = organizedData.(condition).(trial).(area{i});
    spike = spike(~cellfun('isempty',spike));
    x = cellfun('size',spike, 2); 
    [~, Index] = sort(x, 'ascend');
    
    spikes{i,1} = spike(Index);
end
info = organizedData.(condition).(trial).info;


plotRaster_mon(spikes, info); 
suptitle(sprintf('%s %s %s', [area{:}], trial, condition))