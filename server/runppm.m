function runppm(index)
    params = load('test_params.mat');
    load('trainData.mat');

    clusterData = trainData;
    allparams = params.test_params;
    outname = fullfile('glm_fits', strcat('m1neuron_', num2str(index)));
    paramlist = fieldnames(allparams);
    params = allparams.(paramlist{index});

    % run actual ppm
    [X,Y,beta,stats] = pointProcessModel(clusterData, outname, params);
end
