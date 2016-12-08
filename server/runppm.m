function runppm(index)
    for index=1:40
    params = load('test_params.mat');
    load('trainData.mat');
    
    try
        index = str2num(index);
    catch e
        disp(e);
    end
    disp(['index is ', num2str(index)]);
    
    clusterData = trainData;
    allparams = params.test_params;
    outname = fullfile('glm_fits', strcat('m1neuron_', num2str(index)));
    paramlist = fieldnames(allparams);
    params = allparams.(paramlist{index});

    % run actual ppm
    [X,Y,beta,stats] = pointProcessModel(clusterData, outname, params);
    
    disp(['Finished ', num2str(index)]);
    end
end
