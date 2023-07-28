%% summary all data 

% acute recordings trained 


% chronic recordings


% acute recordings multi-image world (but this was two-sided, a whole
% thing)


% implant stability
%JF067, JF082 (vs JF078, JF084)

% learning phase 2 : JF082/84 against the world 
animalsAP = { ...
    'AP100','AP101','AP103','AP104','AP105', ...
    'AP106','AP107','AP108','AP109','AP111', ...
    'AP113','AP114','AP115'};

animalsJF_acute = { ...
    'JF070', 'JF072', 'JF073'};

animalsJF_chronic = { ...
    'JF067', 'JF078', 'JF082', 'JF084'};
%plot = False; 
bhvData = noGoWorld_behavior_debug(animalsJF_chronic);
bhvData = noGoWorld_behavior_debug(animalsJF_acute);
bhvData = noGoWorld_behavior_debug(animalsAP);

timeToLearn_chronic = [5, 4];
timeToLearn_acuteRec = [6,6,7];
timeToLearn_wierdos = [13, 16];

figure();
swarmchart(ones(13,1), [3,5,6,5,5,4,4,5,6,3,5,4,4], 'filled'); hold on;
swarmchart([3,3], timeToLearn_chronic, 'filled'); hold on;
swarmchart([2,2,2], timeToLearn_acuteRec, 'filled'); hold on;
swarmchart([4,4], timeToLearn_wierdos, 'filled')
ylabel('days to learn go 1')
xticks([1, 2, 3, 4])
xticklabels({'AP mice', 'acute', 'chronic other', 'wierdos'})
