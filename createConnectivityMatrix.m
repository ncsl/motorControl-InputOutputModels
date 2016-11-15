data = load('organizedData.mat');
data = data.organizedData;

mallet = data.mallet;
push = data.push;
sphere = data.sphere;
pull = data.pull;

experiment = mallet;
trials = fieldnames(experiment);
% loop through all trials
for i=1:length(trials)
    trialname = trials{i};
    trial = experiment.(trialname);
    
    S1 = trial.S1;
    Pmv = trial.PMv;
    Pmd = trial.PMd;
    M1 = trial.M1;
    
end