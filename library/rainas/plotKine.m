function [  ] = plotKine( kine )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [  ] = plotKine( kine )
%-----------------------------------------------------------------------------------------
%
% Description:  Creates figure represeting kinematic trial information in
%               reduced dimension principal component space. Plots each
%               conditions mean as well as error trials that exceed three standard
%               deviations.
%
%-----------------------------------------------------------------------------------------
%   
%   Input:    kine       -   A structure  containing reduced dimension kinematic information from all trials within each
%                            condition, aligned -750ms to 750ms around movement epoch. Each condition contains a ntrial x 4
%                            cell, where all 30 sensor coordinates have been reduced to a 1x300 matrix for x, y, z.
%                            The fourth column contains the trial name string.
% 
%   Output:   Figure is the output
%                          
%-----------------------------------------------------------------------------------------
% Author: R D'Aleo
%
% Ver.: 1.0 - Date: 08/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%----------------------------------------------------------------------------------------%
% 1. Create Figure Guidelines 
%----------------------------------------------------------------------------------------%
%

condition = fieldnames(kine);

clf('reset');
fig = figure(1);
set(fig, 'Position',[400 100 1300 900]); 
set(fig, 'color', [1 1 1]);
hold on;

color = colormap(cool(4)); 
nms = {};

%%
subplot(2,3,[1,3]);
hold on
ax = gca; 
ax.XTickLabel = {}; 
ax.YTickLabel = {};
ax.ZTickLabel = {}; 
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

% subplot(2,3,3);
% hold on
% ax = gca; 
% ax.XTickLabel = {}; 
% ax.YTickLabel = {}; 
% ax.ZTickLabel = {}; 

xlabel('x');
ylabel('y');
zlabel('z');
view(3);

subplot(2,3,4);
hold on
ax = gca; 
ax.XTickLabel = {}; 
ax.YTickLabel = {}; 

xlabel('x');
ylabel('y');

subplot(2,3,5);
hold on
ax = gca; 
ax.XTickLabel = {}; 
ax.YTickLabel = {}; 

xlabel('x');
ylabel('z');

subplot(2,3,6);
hold on
ax = gca; 
ax.XTickLabel = {}; 
ax.YTickLabel = {}; 

xlabel('y');
ylabel('z');

%%
%----------------------------------------------------------------------------------------%
% 2. Establish condition mean and error trials 
%----------------------------------------------------------------------------------------%
%

for i =  1 : length(condition); 
    colr = color(i,:);
    corr = colr+(1-colr)*(0.5);
    count = 1;

    
    tmpx = cell2mat(kine.(condition{i})(:,1));
    tmpy = cell2mat(kine.(condition{i})(:,2));
    tmpz = cell2mat(kine.(condition{i})(:,3));
    
    mux = mean(tmpx);
    muy = mean(tmpy);
    muz = mean(tmpz);

    stdx = std(tmpx);
    stdy = std(tmpy);
    stdz = std(tmpz);
    
    trialnames = kine.(condition{i})(:,4);
    ntrials = size(trialnames, 1);
% 
%     for j = 1 : ntrials;
%         
%         outlier(1,:) = ((mux(150:end) + 3*stdx(150:end)) - tmpx(j,150:end));   %Error occurs after movement onset
%         outlier(2,:) = (tmpx(j,150:end) - (mux(150:end) - 3*stdx(150:end)));
%         outlier(3,:) = ((muy(150:end) + 3*stdy(150:end)) - tmpy(j,150:end));
%         outlier(4,:) = (tmpy(j,150:end) - (muy(150:end) - 3*stdy(150:end)));
%         outlier(5,:) = ((muz(150:end) + 3*stdz(150:end)) - tmpz(j,150:end));
%         outlier(6,:) = (tmpz(j,150:end) - (muz(150:end) - 3*stdz(150:end)));
% % 
% %         outlier(1,:) = ((mux + 3*stdx) - tmpx(j,:));                         %Error occurs before or after movement onset
% %         outlier(2,:) = (tmpx(j,:) - (mux - 3*stdx));
% %         outlier(3,:) = ((muy + 3*stdy) - tmpy(j,:));
% %         outlier(4,:) = (tmpy(j,:) - (muy- 3*stdy));
% %         outlier(5,:) = ((muz + 3*stdz) - tmpz(j,:));
% %         outlier(6,:) = (tmpz(j,:) - (muz - 3*stdz));
% 
% %         outlier(1,:) = ((mux(1:150) + 3*stdx(1:150)) - tmpx(j,1:150));       %Error occurs before movement onset
% %         outlier(2,:) = (tmpx(j,1:150) - (mux(1:150) - 3*stdx(1:150)));
% %         outlier(3,:) = ((muy(1:150) + 3*stdy(1:150)) - tmpy(j,1:150));
% %         outlier(4,:) = (tmpy(j,1:150) - (muy(1:150) - 3*stdy(1:150)));
% %         outlier(5,:) = ((muz(1:150) + 3*stdz(1:150)) - tmpz(j,1:150));
% %         outlier(6,:) = (tmpz(j,1:150) - (muz(1:150) - 3*stdz(1:150)));
%         
%         if ~isempty(find(outlier < 0,1))  % Plot error trials 
%             nms{i,count} = trialnames{j}(6:end);
%             count = count + 1;
%             
%             subplot(2,3,4);
%             plot(tmpx(j,:), tmpy(j,:), 'Color', corr,'LineWidth', 1);
%             
%             subplot(2,3,5);
%             plot(tmpx(j,:), tmpz(j,:), 'Color', corr,'LineWidth', 1);
%             
%             subplot(2,3,6);
%             plot(tmpy(j,:), tmpz(j,:), 'Color', corr,'LineWidth', 1);
%         end
%      
%     end

    % Plot condition mean
    
    subplot(2,3,[1,3]); 
    line(i) = plot3(mux, muy, muz, 'Color', color(i,:),'LineWidth', 3);
    plot3(mux(150),muy(150),muz(150), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'g');

%     subplot(2,3,3);
%     plot3(mux(1:152), muy(1:152), muz(1:152), 'Color', colr,'LineWidth', 3);
%     plot3(mux(150),muy(150),muz(150), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'g');
%     plot3(mux(150:151),muy(150:151),muz(150:151), 'Color', 'g','LineWidth', 1); 
%     plot3(mux(100),muy(100),muz(100), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'b');
 

    subplot(2,3,4);
    ax.XTickLabel = {}; 
    plot(mux, muy, 'Color', colr,'LineWidth', 2);
    %shadedErrorBar3D(1:size(mux,2),[mux;muy], [stdx;stdy], {'Color', color(i,:)}); 
    plot(mux(150),muy(150), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'g');

    subplot(2,3,5);
    plot(mux, muz, 'Color', colr,'LineWidth', 2);
    %shadedErrorBar3D(1:size(mux,2),[mux;muz], [stdx;stdz], {'Color', color(i,:)}); 
    plot(mux(150),muz(150), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'g');

    subplot(2,3,6);
    plot(muy, muz, 'Color', colr,'LineWidth', 2);
    %shadedErrorBar3D(1:size(mux,2),[muy;muz], [stdy;stdz], {'Color', color(i,:)}); 
    plot(muy(150),muz(150), 'o', 'MarkerEdgeColor', colr, 'MarkerFaceColor', 'g');
   
end

% subplot(2,3,3);
% axis tight; 
subplot(2,3,[1,3]);
axis tight
    

legend(line, condition{:},'Position',[-110 490 550 300]); % 'Location','bestoutside');

subplot(2,3,6);
x = xlim;
y = ylim;
%text(x(1) + x(1)/-3 , y(2) + y(2)/40 , 'Error Trials')
for i = 1 : size(nms,1)
    text(x(1), y(2) + y(2)/40 - (y(2)*3/40)*i,{sprintf('%s: %s',condition{i}, sprintf('%s ', nms{i,:}))});
end

