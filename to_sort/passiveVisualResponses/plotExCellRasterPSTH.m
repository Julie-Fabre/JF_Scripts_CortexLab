function [pltTime1, pltY1, pltTime2, pltY2,pltTime3, pltY3, pltCol3, pltTime4, pltY4, pltCol4]= plotExCellRasterPSTH(stimIDs,  conditionsNew, colorsO, stimOn_times,spike_templates,...
    spike_times_timeline,iUnit,alphaV,szV, figH, sortThis, plotThis)
if plotThis
figure(figH)
clf;
end
if size(conditionsNew,2)==3
    allOri = unique(conditionsNew(:, 3));
else
    allOri = conditionsNew;
end

raster_window = [-0.25, 0.4];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
timeT = t_bins(1:end-1) + diff(t_bins) ./ 2;

use_align = reshape(stimOn_times(:), [], 1);
t_peri_event = use_align + t_bins;

% A. raster 
curr_spikes_idx = ismember(spike_templates,iUnit);
curr_raster_spike_times =spike_times_timeline(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
    curr_raster_spike_times > max(t_peri_event(:))) = [];

if ~any(diff(reshape(t_peri_event',[],1)) < 0)
    % (if no backward time jumps, can do long bin and cut out in-between, faster)
    curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
        reshape(t_peri_event',[],1)),NaN],size(t_peri_event'))';
    curr_raster = curr_raster_continuous(:,1:end-1);   
else
    % (otherwise, bin trial-by-trial)
    curr_raster = cell2mat(arrayfun(@(x) ...
        histcounts(curr_raster_spike_times,t_peri_event(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false));
end
[raster_y,raster_x] = find(curr_raster);
pltTime1= timeT(raster_x);
pltY1 = raster_y;
if plotThis
clf;

if max(raster_y)<300
    subplot(4,2,[5:6])
    cla;
else
    subplot(4,2,[1:6])
    cla;
end

sc = scatter(timeT(raster_x), raster_y, szV, 'k', 'filled','MarkerFaceAlpha',alphaV,'MarkerEdgeAlpha',alphaV);
xlim([timeT(1), timeT(end)])
ylim([0, size(t_peri_event, 1)])
line([0, 0], [0, size(t_peri_event, 1)],'Color','r')
ylabel('unsorted trial #')
set(gca,'xtick',[])
makepretty;
end
% B. PSTH
smooth_size = 51;
gw = gausswin(smooth_size, 3)';
smWin = gw ./ sum(gw);
bin_t = mean(diff(timeT));

curr_psth = nanmean(curr_raster);
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
    smWin, 'valid') ./ bin_t;
pltTime2 = timeT;
pltY2 = curr_smoothed_psth;
if plotThis
subplot(4,2,7:8)
cla;


plot(timeT,curr_smoothed_psth,'Color', 'k'); hold on;
%plotshaded(t,[curr_smoothed_psth - nanstd(],
xlabel('time from stim onset (s)')
ylim([min(curr_smoothed_psth), max(curr_smoothed_psth)]);
line([0, 0], [min(curr_smoothed_psth), max(curr_smoothed_psth)],'Color','r')
%line([0.3, 0.3], [min(curr_smoothed_psth), max(curr_smoothed_psth)],'Color','r')
xlim([timeT(1), timeT(end)])
ylabel('sp/s')
makepretty;
set(gcf, 'color','white')
end

if plotThis
figure(figH+1)
cla
end
sID = ones(size(stimIDs, 1), 1);
for iOri = 1:size(allOri,1)
%     if iOri==1
%         raster_x
%     end
% if size(conditionsNew,2)==3
%     allOriRows = find(conditionsNew(:, 3) == allOri(iOri));
% else
%     allOriRows = find(conditionsNew == allOri(iOri));
% end
%     
    sID(ismember(stimIDs, allOri(iOri))) = iOri;
end
% B. PSTH
smooth_size = 51;


gw = gausswin(smooth_size, 3)';
smWin = gw ./ sum(gw);
bin_t = mean(diff(timeT));

curr_psth = grpstats(curr_raster,sID,@(x) mean(x,1));
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0, floor(length(smWin)/2)], 'replicate', 'both'), ...
    smWin, 'valid') ./ bin_t;
if sortThis
    [~,sid_sort] = sort(nanmean(curr_smoothed_psth(:,find(timeT>0.1, 1):find(timeT<0.3, 1,'last')),2));
    for iSID=1:size(sid_sort,1)
        sID_new(stimIDs==iSID)=sid_sort(iSID);
    end
    [~,trial_sort] = sort(sID_new);
    
elseif size(conditionsNew,2)==3
    [~,trial_sort] = sort(sID);
else
    
    [~,trial_sort] = sort(stimIDs);
end

curr_raster_sorted = curr_raster(trial_sort,:);
[raster_y,raster_x] = find(curr_raster_sorted);
[~,~,row_group] = unique(sID(trial_sort),'sorted');
raster_dot_color = colorsO(row_group(raster_y),:);
pltTime3 = timeT(raster_x);
pltY3 = raster_y;
pltCol3 = raster_dot_color;
pltTime4 = timeT;
pltY4 = curr_smoothed_psth;
pltCol4 =  colorsO; 
if plotThis
if max(raster_y)<300
    subplot(4,2,[5:6])
    cla;
else
    subplot(4,2,[1:6])
    cla;
end
% colorsO = [rgb('FireBrick'); rgb('Red'); rgb('DarkOrange');...
%      rgb('YellowGreen'); rgb('Green'); rgb('DarkTurquoise'); rgb('Blue'); rgb('Purple')];
%    


raster_dots = scatter(timeT(raster_x), raster_y, szV, 'k', 'filled','MarkerFaceAlpha',alphaV,'MarkerEdgeAlpha',alphaV);
xlim([timeT(1), timeT(end)])
ylim([0, size(t_peri_event, 1)])
line([0, 0], [0, size(t_peri_event, 1)],'Color','k')

set(raster_dots,'CData',raster_dot_color);
% for iO=2:size(allOri,1)
% line([-0.2, 0.5], [size(curr_raster_sorted,1)/size(allOri,1)*(iO-1)+1, size(curr_raster_sorted,1)/size(allOri,1)*(iO-1)+1],'Color','k')
% end

ylabel('sorted trial #')
set(gca,'xtick',[])
makepretty;

subplot(4,2,7:8)
cla;

psth_lines = plot(timeT,curr_smoothed_psth); hold on;
arrayfun(@(align_group) set(psth_lines(align_group), ...
    'XData',timeT,'YData',curr_smoothed_psth(align_group,:), ...
    'Color',colorsO(align_group,:)),1:size(curr_psth,1));
xlabel('time from stim onset (s)')
ylabel('sp/s')
xlim([timeT(1), timeT(end)])
ylim([min(min(curr_smoothed_psth)), max(max(curr_smoothed_psth))]);
line([0, 0], [min(min(curr_smoothed_psth)), max(max(curr_smoothed_psth))],'Color','k')

makepretty;
set(gcf, 'color','white')
end
% get reposonive boolean
%w = gausswin(10);


%                     forSTb = nanmean(filter(w, 1, (curr_raster(1:1:end, 1:150) ./ 0.01)), 2);
%                     forST1 = nanmean(filter(w, 1, (curr_raster(1:1:end, 201:350) ./ 0.01)), 2);
%                    %remove zero rows, skews everyhtng (?)
% 
%                     rm =  ~isnan(forSTb) & ~isnan(forST1) & ~isinf(forSTb) & ~isinf(forST1);
%                     
% 
%                     if sum(rm) > 1 
%                          resp = [(signrank(squeeze(forSTb(rm)), squeeze(forST1(rm))))]
%                          oo = nanmean(curr_smoothed_psth(:,150:400),2)
%                         selIdx = (nanmax(oo)-nanmean(oo))./ nanmax(oo);
%                     end

% % get selectivity index 
% figure(figH)
% if max(raster_y)<300
%     subplot(4,2,[5:6])
%     
% else
%     subplot(4,2,[1:6])
%     
% end
% title(['resp = ' num2str(resp(1))], 'Fontsize',12)
% 
% figure(figH+1)
% if max(raster_y)<300
%     subplot(4,2,[5:6])
%     
% else
%     subplot(4,2,[1:6])
%     
% end
% title(['sel idx = ' num2str(selIdx)], 'Fontsize',12)
end