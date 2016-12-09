function [X,Y,beta,stats] = pointProcessModel(clusterData, outname, varargin)
 
% (Please read the input format below on how to use this function.)
% -------------------------------------------------------------------------
% Usage:
%   [X,Y,beta,stats] = pointProcessModel_cluster(clusterData, kinematics, varargin)
%
% Inputs :
% ======
% clusterData : type 'struct' with format:
% clusterData.(task condition).(trial#).(cluster#) is a matrix of dimension
% [#neurons, #observations]
% clusterData.(tast condition).(trial#).info (contains led cue time and
% switch close time)
%
% varargin:
% params : type 'struct'    
%     params.L -> total length of time for which combined cluster actvity is
%                 considered
%     params.own_hist -> #contiguous self history time bins
%     params.own_avg_hist -> window size for self history in t-100 to
%                            t-own_hist bins
%     params.cluster1_hist -> window size for cluster-1 history in t-100 to
%                            t-own_hist bins
%     params.cluster2_hist -> window size for cluster-1 history in t-100 to
%                            t-own_hist bins
%     params.cluster3_hist -> window size for cluster-1 history in t-100 to
%                            t-own_hist bins
%     params.cluster4_hist -> window size for cluster-1 history in t-100 to
%                            t-own_hist bins
%
% Outputs:
% =======
% X : covariate matrix (sparse)
% Y : observation vector (0 or 1)
% beta : glm parameters
% stats : statistics obtained from glmfit
%
% Author: Sanjukta Nandi Bose
% Last modified : December 06, 2016
% -------------------------------------------------------------------------

global area 

condition = {'mallet','push','pull','sphere'};
area = {'S1','PMv','PMd','M1'};

if isempty(outname)
    fprintf('Please make sure to have a .mat filename to save the function output. Check function usage!');
end

if nargin>2
    params = varargin{1};
    % options is a structure that specifies parameters associated with the
    % history of each neuron/cluster of neurons.
    own_hist = params.own_hist ; % maximum own spiking history for i-th time point is (i-own_hist_max)th bin
    own_avg_hist = params.own_avg_hist ;
    L = params.L;
    cluster1_hist = params.cluster1_hist; 
    cluster2_hist = params.cluster2_hist;
    cluster3_hist = params.cluster3_hist;
    cluster4_hist = params.cluster4_hist;
else
    % time bins are of length 1 ms
    % Default params:
    L = 100; % total time length for which combined cluster spiking activity is considered
    own_hist = 15 ; % maximum own spiking history for i-th time point is (i-own_hist_max)th bin
    own_avg_hist = 10;
    cluster1_hist = 10; 
    cluster2_hist = 10;
    cluster3_hist = 10;
    cluster4_hist = 10;
end

% AR order for each cluster?
n1 = cluster1_hist; 
n2 = cluster2_hist; 
n3 = cluster3_hist; 
n4 = cluster4_hist;

trials = fieldnames(clusterData.mallet);
N = size(clusterData.mallet.(trials{1}).Apc1,1)-1; % Number of neighboring neurons in the same cluster

p0 = round((100-own_hist)/own_avg_hist);
p1 = round(L/cluster1_hist);
p2 = round(L/cluster2_hist);
p3 = round(L/cluster3_hist);
p4 = round(L/cluster4_hist);
num_covar = 1 + own_hist + p0 + p1 +  p2  +  p3 + p4;
X = sparse(1,num_covar); % initialize column size for input covariate matrix for glm
count =1;

% Pick a neuron number.
% default: an M1 neuron to be chosen.
% numNeuronsInCluster = size(clusterData.mallet.(trials{1}).Apc1,1);
% clusterVec = clusterData.mallet.(trials{1}).info.areaM(:,1);
% clusterVec = clusterVec(~cellfun('isempty', clusterVec));
% neuron_num = cellfun(@(x) strcmp(x,'M1'), clusterVec, 'UniformOutput', 1);
% neuron_num = find(neuron_num==1, neuronindex);
% neuron_num = find(ismember(clusterData.mallet.(trials{1}).info.areaM(numNeuronsInCluster,1),'M1')==1,2); % finds the 1st M1 neuron in cluster Apc1 (cluster1) : Target neuron
% neuron_num = find(ismember(clusterData.mallet.(trials{1}).info.areaM(size(clusterData.mallet.(trials{1}).Apc1,1),1),'M1')==1,1); % finds the 1st M1 neuron in cluster Apc1 (cluster1) : Target neuron

% neuron index passed in set{1,...16} for APC1
neuron_num = params.neuronIndex;

cluster_ngbrs = 1:N+1; % cluster_ngbrs,indices of neurons other than target neuron

cluster_ngbrs(neuron_num) = [];

for i = 1:length(condition) % loops through task-conditions
    trials = fieldnames(clusterData.(condition{i}));
    num_trials = length(trials);
    
    for j = 1:num_trials %1:50 % % loops through trials
        try
            cluster = fieldnames(clusterData.(condition{i}).(trials{j}));
        catch
            fprintf('%s.%s.%s does not exist. \n',clusterData, condition{i}, trials{j});
            continue;
        end
        
        num_cluster = length(cluster);
        t_start = 501 - round(100*clusterData.(condition{i}).(trials{j}).info.led);
        t_end = 501 + round(100*clusterData.(condition{i}).(trials{j}).info.close);
        % takes time points only between led cue and stop time
        if (t_end - 501 > 100)
            t_end = 501 + 100;
        end
        if t_start < 100
            fprintf('Invalid start time. \n');
        end
        
        for k = 1  % current cluster number
            
            for t = t_start : t_end
                % 'neuron_num's own history from 0 -> param
                test = sparse(1,num_covar);
                test(1,1) = i;
                test(1,2:own_hist+1) = clusterData.(condition{i}).(trials{j}).(cluster{k})(neuron_num,t-own_hist:t-1);
                idx = own_hist + 1;
               
                % 'neuron_num's own history summed history
                temp = reshape(clusterData.(condition{i}).(trials{j}).(cluster{k})(neuron_num,t-(own_avg_hist*p0)-own_hist:t-own_hist-1),[own_avg_hist,p0]);
                temp = sum(temp,1);
                test(idx+1:idx+p0) = temp;
                clear temp;
                
                % cluster1_avg_hist
                idx = idx+ p0;
                temp1 = [];
                for p = 1:size(clusterData.(condition{i}).(trials{j}).(cluster{1}),1)
                    if p == neuron_num
                        continue;
                    else
                        temp = reshape(clusterData.(condition{i}).(trials{j}).(cluster{1})(p,t-(n1*p1):t-1),[n1,p1]);
                        temp = sum(temp,1);
                        temp1 = [temp1; temp];
                        clear temp;
                    end
                end
                test(idx+1:idx+p1) = sum(temp1,1);
                clear temp1;
                
                % cluster2_avg_hist
                idx = idx+ p1;
                temp1 = [];
                numNeuronsClust2 = size(clusterData.(condition{i}).(trials{j}).(cluster{2}),1);
                for p = 1:numNeuronsClust2 % loop through number of clust 2 neurons
                    temp = reshape(clusterData.(condition{i}).(trials{j}).(cluster{2})(p,t-(n2*p2):t-1),[n2,p2]);
                    temp = sum(temp,1);
                    temp1 = [temp1; temp];
                    clear temp;
                end
                test(idx+1:idx+p2) = sum(temp1,1);
                clear temp1;
                
                % cluster3_avg_hist
                idx = idx+ p2;
                temp1 = [];
                for p = 1:size(clusterData.(condition{i}).(trials{j}).(cluster{3}),1)
                    
                    temp = reshape(clusterData.(condition{i}).(trials{j}).(cluster{3})(p,t-(n3*p3):t-1),[n3,p3]);
                    temp = sum(temp,1);
                    temp1 = [temp1; temp];
                    clear temp;
                    
                end
                test(idx+1:idx+p3) = sum(temp1,1);
                clear temp1;
                
                % cluster4_avg_hist
                idx = idx+ p3;
                temp1 = [];
                for p = 1:size(clusterData.(condition{i}).(trials{j}).(cluster{4}),1)
                    
                    temp = reshape(clusterData.(condition{i}).(trials{j}).(cluster{4})(p,t-(n4*p4):t-1),[n4,p4]);
                    temp = sum(temp,1);
                    temp1 = [temp1; temp];
                    clear temp;
                    
                end
                test(idx+1:idx+p4) = sum(temp1,1);
                clear temp1;
                
                if count ==1
                    X = test;
                else
                    X = [X; test];
                end
                Y(count) = clusterData.(condition{i}).(trials{j}).(cluster{k})(neuron_num,t);
                clear test;
                count = count + 1;
                
                % kinematic data for you to add...
                
            end % loop through time index
        end % loop thru cluster
    end % loop thru trial
end % loop thru condition

% Fit glm model to obtain the point process model with link function
% log(lambda) = X * beta; (if Poisson dist.), beta: parameters of the glm model
% log(p_hat/(1-p_hat)) = X * beta (if binomial dist.)


X = full(X);

[beta,dev,stats] = glmfit(X,Y(:),'Binomial');

% X= sparse(X);
% Y = sparse(Y(:));

% save(outname,'X','Y','beta','stats','-v7.3');
save(outname,'params','stats','-v7.3');
end
