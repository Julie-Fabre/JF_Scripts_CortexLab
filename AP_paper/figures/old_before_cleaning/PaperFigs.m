
%% ANDY PAPER FIGURES - JF

%% ~~1. CLASSIFICATION~~
clear all;
naiveComp = 0;
data_fn = {'E:/trial_activity_muscimol_passiveCorr','E:/trial_activity_ctx_passiveCorr','E:\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld'};


loadAllTask; 
SUAbugsBgone; 
trained = 0;
naiveComp = 0;
loadingTrialDataUnpack;
%% ~METRICS VENN DIAGRAM~
% % missing, number spikes, duration, peaksTroughd,rpv,
% somatic, amplitude.
param = struct;
param.minNumSpikes = 300;
param.maxNumPeak = 3;
param.minAmpli = 15;
param.maxPercMissing = 30;
param.maxRPV = 2;
param.somaCluster = 1;

numPeaksTroughs_all = vertcat(numPeaksTroughs{:});
numPeaksTroughs_all = cellfun(@transpose, numPeaksTroughs_all, 'UniformOutput', false);
numPeaksTroughs_all = vertcat(numPeaksTroughs_all{:});
numPeaksTroughs_all = numPeaksTroughs_all > param.maxNumPeak;

numSpikes_all = vertcat(numSpikes{:});
numSpikes_all = cellfun(@transpose, numSpikes_all, 'UniformOutput', false);
numSpikes_all = vertcat(numSpikes_all{:});
numSpikes_all = numSpikes_all < param.minNumSpikes;

percent_missing_ndtr_all = vertcat(percent_missing_ndtr{:});
thesemins = [];
for iRec = 1:size(percent_missing_ndtr_all, 1)
    pp = percent_missing_ndtr_all{iRec};
    pp = pp(2:end, :);
    thesemins = [thesemins, nanmin(pp)];
end
percent_missing_ndtr_all = thesemins > 30;

waveformRawAmpli_all = vertcat(waveformRawAmpli{:});
waveformRawAmpli_all = cellfun(@transpose, waveformRawAmpli_all, 'UniformOutput', false);
waveformRawAmpli_all = vertcat(waveformRawAmpli_all{:});
waveformRawAmpli_all = waveformRawAmpli_all .* 0.195 < param.minAmpli;

fractionRPVchunk_all = vertcat(fractionRPVchunk{:});
fractionRPVchunk_all = cellfun(@transpose, fractionRPVchunk_all, 'UniformOutput', false);
fractionRPVchunk_all = vertcat(fractionRPVchunk_all{:});
fractionRPVchunk_all = fractionRPVchunk_all > param.maxRPV;

templateDuration_all = vertcat(templateDuration{:});
templateDuration_all = cellfun(@transpose, templateDuration_all, 'UniformOutput', false);
templateDuration_all = vertcat(templateDuration_all{:});
templateDuration_all = templateDuration_all > 1000;

somatic_all = vertcat(somatic{:});
somatic_all = cellfun(@transpose, somatic_all, 'UniformOutput', false);
somatic_all = vertcat(somatic_all{:});
somatic_all = somatic_all == 0;

%A num peaks
%B num spikes
%C ampli
%D rpv
%E duration
%F somatic
%G perc. missing

T = table(numPeaksTroughs_all, numSpikes_all, waveformRawAmpli_all, fractionRPVchunk_all, templateDuration_all, somatic_all, percent_missing_ndtr_all');
writetable(T, 'C:/Users/Julie/Dropbox/Paper/metricsTable.csv');

%-> loaded into R with script qualityMetricVenn.R

%% ~METRICS: KEPT UNIT NUMBERS + GRAPH~
%% mean +/- sem
prop_allcat = vertcat(prop_long_isi{:});
dd = cellfun(@transpose, depth_aligned, 'UniformOutput', false);
depthThis_allcat = vertcat(dd{:});
good_allcat = vertcat(goodUnits{:});

for iRec = 1:45

    thisval(iRec) = sum(good_allcat{iRec}) / numel(prop_allcat{iRec});


end
disp(strcat(num2str(mean(thisval))), '+/-', num2str(AP_sem(thisval, 2)))))
disp(strcat(num2str(mean(thisval)), '+/-', num2str(std(thisval, 2)))
%% number of kept/number total
propall_allcat = vertcat(prop_long_isi{:});
propall_allcat = cellfun(@transpose, propall_allcat, 'UniformOutput', false);
propall_allcat = vertcat(propall_allcat{:});


goodall_allcat = vertcat(goodUnits{:});
goodall_allcat = cellfun(@transpose, goodall_allcat, 'UniformOutput', false);
goodall_allcat = vertcat(goodall_allcat{:});

disp(strcat(num2str(sum(goodall_allcat)), '/', num2str(numel(propall_allcat))))
%% graph
cell_clean = celltype_allcat(goodUnits_allcat);
domain_clean = domain_aligned_allcat(goodUnits_allcat);
domain_clean(cell_clean==5)=[];
cell_clean(cell_clean==5)=[];
cell_clean(cell_clean==6)=3;

domain_type_count = accumarray([domain_clean, cell_clean], 1);

figure;
H=bar(domain_type_count, 'stacked');
H(1).FaceColor = 'flat';
H(1).CData = rgb('DarkRed');
H(2).FaceColor = 'flat';
H(2).CData = rgb('RoyalBlue');
H(3).FaceColor = 'flat';
H(3).CData = rgb('Lime');
H(4).FaceColor = 'flat';
H(4).CData = rgb('Pink');

xlabel('Domain');
ylabel('Count');
legend({'MSN', 'FSI', 'TAN', 'UIN'});
makepretty;

cell_clean = celltype_allcat;
domain_clean = domain_aligned_allcat;
domain_clean(cell_clean==5)=[];
cell_clean(cell_clean==5)=[];
cell_clean(cell_clean==6)=3;

domain_type_count = accumarray([domain_clean, cell_clean], 1);

figure;
H=bar(domain_type_count, 'stacked');
H(1).FaceColor = 'flat';
H(1).CData = rgb('DarkRed');
H(2).FaceColor = 'flat';
H(2).CData = rgb('RoyalBlue');
H(3).FaceColor = 'flat';
H(3).CData = rgb('Lime');
H(4).FaceColor = 'flat';
H(4).CData = rgb('Pink');

xlabel('Domain');
ylabel('Count');
legend({'MSN', 'FSI', 'TAN', 'UIN'});
makepretty;
%% ~EXAMPLES OF GOOD AND BAD UNITS~

%% rpv

%% spikes missing

%% # spikes

%% wav

%% amplitude



%% ~EX PLOT~

prop_allcat = vertcat(prop_long_isi{:});
prop_allcat = cellfun(@transpose, prop_allcat, 'UniformOutput', false);
prop_allcat = vertcat(prop_allcat{:});

for iCell = 1:size(acg_allcat, 1)
    pss_allcat2temp = find(acg_allcat(iCell, 500:1000) >= nanmean(acg_allcat(iCell, 600:900)));
    pss_allcat2(iCell) = pss_allcat2temp(1);
end
tanMethod1 = goodUnits_allcat & pss_allcat > 20;
tanMethod2 = goodUnits_allcat & pss_allcat2' > 40 & templateDuration_us<850 & templateDuration_us>400;

[~, waveform_trough] = min(wv_allcat, [], 2);
[~, waveform_peak_rel] = arrayfun(@(x) ...
    max(wv_allcat(x, waveform_trough(x):end), [], 2), ...
    transpose(1:size(wv_allcat, 1)));
waveform_peak = waveform_peak_rel + waveform_trough;

templateDuration = waveform_peak - waveform_trough;
templateDuration_us = (templateDuration / 30000) * 1e6;

recording_neurons = find(celltype_allcat == 1); %recordings ordered by cell type
recordingsCidx = [diff(recording_neurons')~=1, true]; %[diff(recording_neurons')~=1, true]
recording_neurons = recording_neurons(recordingsCidx);

theseNeu = zeros(size(depth_aligned_allcat, 1), 1);
for iRec = 1:45
    if iRec == 1
        theseNeu(1:recording_neurons(iRec)) = iRec;
    else
        theseNeu(recording_neurons(iRec-1)+1:recording_neurons(iRec)) = iRec;
    end
end

iRec = 12;
figure();
clf;
scatter(pss_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 1), ...
    templateDuration_us(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 1), [], rgb('DarkRed'), 'filled')
hold on;
scatter(pss_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 2), ...
    templateDuration_us(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 2), [], rgb('RoyalBlue'), 'filled')
hold on;
scatter(pss_allcat(theseNeu == iRec & goodUnits_allcat & tanMethod1), ...
    templateDuration_us(theseNeu == iRec & goodUnits_allcat & tanMethod1), [], rgb('Lime'), 'filled')
hold on;
scatter(pss_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 4), ...
    templateDuration_us(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 4), [], rgb('HotPink'), 'filled')
xlabel('post spike suppression (ms)')
ylabel('peak to trough (\mus)')
ylim([150, 800])
makepretty;

figure();


scatter(prop_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 1), ...
    nanmean(permute(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 1:18, theseNeu == iRec & ...
    goodUnits_allcat & celltype_allcat == 1))), [2, 1]), 2), [], rgb('DarkRed'), 'filled');
hold on;
scatter(prop_allcat(theseNeu == iRec & goodUnits_allcat & tanMethod1), ...
    nanmean(permute(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 1:18, theseNeu == iRec & ...
    goodUnits_allcat & tanMethod1))), [2, 1]), 2), [], rgb('Lime'), 'filled');
scatter(prop_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 2), ...
    nanmean(permute(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 1:18, theseNeu == iRec & ...
    goodUnits_allcat & celltype_allcat == 2))), [2, 1]), 2), [], rgb('RoyalBlue'), 'filled');

scatter(prop_allcat(theseNeu == iRec & goodUnits_allcat & celltype_allcat == 4), ...
    nanmean(permute(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 1:18, theseNeu == iRec & ...
    goodUnits_allcat & celltype_allcat == 4))), [2, 1]), 2), [], rgb('HotPink'), 'filled');

xlabel('Prop. ISI>2s')
ylabel('mean ITI f.r. (Hz)')
makepretty;

%% ~FULL PLOT~
%% HM template dur 
% pssBinned(:,1,1) = histcounts(pss_allcat2(goodUnits_allcat & celltype_allcat == 1), [0:10:185]);
% pssBinned(:,1,1) = pssBinned(:,1,1)./sum(pssBinned(:,1,1));
% pssBinned(:,1,2) = histcounts(pss_allcat2(goodUnits_allcat& celltype_allcat == 2), [0:10:185]);
% pssBinned(:,1,2) = pssBinned(:,1,2)./sum(pssBinned(:,1,2));
% pssBinned(:,1,3) = histcounts(pss_allcat2(goodUnits_allcat& (celltype_allcat == 3 |  celltype_allcat == 6)), [0:10:185]);
% pssBinned(:,1,3) = pssBinned(:,1,3)./sum(pssBinned(:,1,3));
% pssBinned(:,1,4) = histcounts(pss_allcat2(goodUnits_allcat& celltype_allcat == 4), [0:10:185]);
% pssBinned(:,1,4) = pssBinned(:,1,4)./sum(pssBinned(:,1,4));
% 
% pssBinned(:,2,1) = histcounts(templateDuration_us(goodUnits_allcat & celltype_allcat == 1), [0:54:1000]);
% pssBinned(:,2,1) = pssBinned(:,1,1)./sum(pssBinned(:,1,1));
% pssBinned(:,2,2) = histcounts(templateDuration_us(goodUnits_allcat& celltype_allcat == 2), [0:54:1000]);
% pssBinned(:,2,2) = pssBinned(:,1,2)./sum(pssBinned(:,1,2));
% pssBinned(:,2,3) = histcounts(templateDuration_us(goodUnits_allcat& (celltype_allcat == 3 |  celltype_allcat == 6)), [0:54:1000]);
% pssBinned(:,2,3) = pssBinned(:,1,3)./sum(pssBinned(:,1,3));
% pssBinned(:,2,4) = histcounts(templateDuration_us(goodUnits_allcat& celltype_allcat == 4), [0:54:1000]);
% pssBinned(:,2,4) = pssBinned(:,1,4)./sum(pssBinned(:,1,4));

%catPssDur = cat(2, pssBinned, durBinned));
clearvars pDcounts
[pDcounts(:,:,1), xedges, yedges] = histcounts2(pss_allcat2(goodUnits_allcat & celltype_allcat == 1), templateDuration_us(goodUnits_allcat & celltype_allcat == 1)',[0:7:185], [0:34:1000]);
pDcounts(:,:,1) = pDcounts(:,:,1)./sum(sum(pDcounts(:,:,1),1),2);
pDcounts(:,:,2) = histcounts2(pss_allcat2(goodUnits_allcat & celltype_allcat == 2), templateDuration_us(goodUnits_allcat & celltype_allcat == 2)',[0:7:185], [0:34:1000]);
pDcounts(:,:,2) = pDcounts(:,:,2)./sum(sum(pDcounts(:,:,2),1),2);
pDcounts(:,:,3) = histcounts2(pss_allcat2(tanMethod2), ...
    templateDuration_us(tanMethod2)',[0:7:185], [0:34:1000]);
pDcounts(:,:,3) = pDcounts(:,:,3)./sum(sum(pDcounts(:,:,3),1),2);
pDcounts(:,:,4) = histcounts2(pss_allcat2(goodUnits_allcat & celltype_allcat == 4), templateDuration_us(goodUnits_allcat & celltype_allcat == 4)',[0:7:185], [0:34:1000]);
pDcounts(:,:,4) = pDcounts(:,:,4)./sum(sum(pDcounts(:,:,4),1),2);
figure(); 
imagesc(yedges,xedges,pDcounts(:,:,1)+ pDcounts(:,:,2)+ pDcounts(:,:,3)+ pDcounts(:,:,4))
hold on; 
scatter( templateDuration_us(tanMethod2)',...
    pss_allcat2(tanMethod2),'k')
xlabel('template duration (\mus)')

ylabel('post spike suppression (ms)')
makepretty;
% n bins x n bins x 3 cell type matrix 'counts', you can do something like image(counts./sum(sum(counts,1),2)) and that would show up as a normalized and color-coded heatmap
%% HM 
%% template duration
figure();
hist(templateDuration_us(goodUnits_allcat & templateDuration_us<1000), 23)
xlabel('peak to trough (\mus)')
ylabel('# cells')
h = findobj(gca, 'Type', 'patch');
h.FaceColor = rgb('DarkOrange');
makepretty;

figure();
hist(templateDuration_us(:), 35)
xlabel('peak to trough (\mus)')
ylabel('# cells')
h = findobj(gca, 'Type', 'patch');
h.FaceColor = rgb('DarkOrange');
makepretty;

%% post spike supression
figure();
histogram(pss_allcat2(goodUnits_allcat), 40, 'FaceColor', rgb('DarkOrange'))
xlabel('post spike suppression (ms)')
ylabel('# cells')
% h = findobj(gca,'Type','patch');
% h.FaceColor =  rgb('DarkOrange');
ylim([0, 3400])
breakyaxis([200, 3300])
%set(gca, 'YScale', 'log')
makepretty;

figure();
histogram(pss_allcat2(:), 40, 'FaceColor', rgb('DarkOrange'))
xlabel('post spike suppression (ms)')
ylabel('# cells')
% h = findobj(gca,'Type','patch');
% h.FaceColor =  rgb('DarkOrange');
ylim([0, 12400])
xlim([0, 250])
breakyaxis([1000, 12200])

%set(gca, 'YScale', 'log')
makepretty;

%% prop ISI
figure();
histogram(prop_allcat(goodUnits_allcat & (celltype_allcat == 2 | celltype_allcat == 4)), 50,'FaceColor', rgb('DarkOrange'))
xlabel('prop ISI > 2s')
ylabel('# of cells')
makeprettyLite;

figure();
histogram(prop_allcat( (celltype_allcat == 2 | celltype_allcat == 4)), 50,'FaceColor', rgb('DarkOrange'))
xlabel('prop ISI > 2s')
ylabel('# of cells')
makeprettyLite;


figure();
scatter(nanmean(permute(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 1:18, goodUnits_allcat & (celltype_allcat == 2 | celltype_allcat == 4)))), [2, 1]), 2), ...
    prop_allcat(goodUnits_allcat & (celltype_allcat == 2 | celltype_allcat == 4)))
ylabel('prop ISI > 2s')
xlabel('baseline f.r. (Hz)')
makeprettyLite;

%% ~WV, ACG, F.R~
%wv
waveform_t = 1e3*((0:size(wv_allcat,2)-1)/30000);
figure();
title('Mean waveforms');hold on;
P(1)=plot(waveform_t, mean(wv_allcat(goodUnits_allcat & celltype_allcat==2,:)),':','Color',rgb('RoyalBlue'),'LineWidth',2);
plotshaded(waveform_t,[-std(wv_allcat(goodUnits_allcat & celltype_allcat==2,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==2,:)); ...
    std(wv_allcat(goodUnits_allcat & celltype_allcat==2,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==2,:))],rgb('RoyalBlue'));
       
hold on;
P(2)=plot(waveform_t,mean(wv_allcat(goodUnits_allcat & celltype_allcat==1,:)), '-.','Color', rgb('DarkRed'),'LineWidth',2);
plotshaded(waveform_t,[-std(wv_allcat(goodUnits_allcat & celltype_allcat==1,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==1,:)); ...
    std(wv_allcat(goodUnits_allcat & celltype_allcat==1,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==1,:))],rgb('DarkRed'));

hold on;
P(3)=plot(waveform_t,mean(wv_allcat(tanMethod1,:)),'Color',rgb('Lime'),'LineWidth',2);
plotshaded(waveform_t,[-std(wv_allcat(tanMethod1,:)) + mean(wv_allcat(tanMethod1,:)); std(wv_allcat(tanMethod1,:)) + mean(wv_allcat(tanMethod1,:))],rgb('Lime'));
xlim([0 waveform_t(end)])
makepretty;

hold on;
P(4)=plot(waveform_t,mean(wv_allcat(goodUnits_allcat & celltype_allcat==4,:)), '--','Color', rgb('HotPink'),'LineWidth',2);
plotshaded(waveform_t,[-std(wv_allcat(goodUnits_allcat & celltype_allcat==4,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==4,:));...
    std(wv_allcat(goodUnits_allcat & celltype_allcat==4,:)) + mean(wv_allcat(goodUnits_allcat & celltype_allcat==4,:))],rgb('HotPink'));
xlim([0 waveform_t(end)])
makepretty;

ylabel('Amplitude (\muV)')
xlabel('Time (ms)')
legend(P,{'FSI mean +/- std','MSN mean +/- std','TAN mean +/- std', 'UIN mean +/- std'});
makepretty;

%acg 
% acgmean=smoothdata(mean(acg(fsi,500:end)./50),'gaussian', [0 5]);
% acgstd=smoothdata(std(acg(fsi,500:end)./50),'gaussian', [200 200]);
figure();
subplot(141)
cla();
title('FSIs');hold on;
plot(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat & celltype_allcat==2,501:end)), 'Color', rgb('RoyalBlue'))
%plot(0:0.001:0.5,mean(acg(fsi,500:end)),'b');%,'FaceColor','b');
hold on;
plotshaded(0:0.001:0.5,[-nanstd(acg_allcat(goodUnits_allcat & celltype_allcat==2,501:end))+ nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==2,501:end));...
    nanstd(acg_allcat(goodUnits_allcat &celltype_allcat==2,501:end)) + nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==2,501:end))],rgb('RoyalBlue'));
area(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==2,501:end)),'FaceColor',rgb('RoyalBlue'));
xlim([0 0.5])
ylim([0 50])
ylabel('sp/s')
xlabel('Time (s)')
makepretty;

subplot(142)
cla();
title('MSNs');hold on;
plot(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end)), 'Color', rgb('DarkRed'))
%plot(0:0.001:0.5,mean(acg(fsi,500:end)),'b');%,'FaceColor','b');
hold on;
plotshaded(0:0.001:0.5,[- smoothdata(nanstd(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end)),'gaussian', [0 100])+ nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end));...
    smoothdata(nanstd(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end)),'gaussian', [0 100]) + nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end))],rgb('DarkRed'));
area(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==1,501:end)),'FaceColor',rgb('DarkRed'));
xlim([0 0.5])
ylim([0 50])
%ylim([0 150])
xlabel('Time (s)')
set(gca,'YTickLabel',[]);
makepretty;

subplot(143)
cla();
title('TANs');hold on;
plot(0:0.001:0.5,nanmean(acg_allcat(tanMethod1,501:end)), 'Color', rgb('Lime'))
%plot(0:0.001:0.5,mean(acg(fsi,500:end)),'b');%,'FaceColor','b');
hold on;
plotshaded(0:0.001:0.5,[-nanstd(acg_allcat(tanMethod1,501:end))+ nanmean(acg_allcat(tanMethod1,501:end));...
    nanstd(acg_allcat(tanMethod1,501:end)) + nanmean(acg_allcat(tanMethod1,501:end))],rgb('Lime'));
area(0:0.001:0.5,nanmean(acg_allcat(tanMethod1,501:end)),'FaceColor',rgb('Lime'));
xlim([0 0.5])
ylim([0 50])
%ylim([0 150])
xlabel('Time (s)')
set(gca,'YTickLabel',[]);
makepretty;

subplot(144)
cla();
title('UINs');hold on;
plot(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end)), 'Color', rgb('HotPink'))
%plot(0:0.001:0.5,mean(acg(fsi,500:end)),'b');%,'FaceColor','b');
hold on;
plotshaded(0:0.001:0.5,[-nanstd(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end))+ nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end));...
    nanstd(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end)) + nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end))],rgb('HotPink'));
area(0:0.001:0.5,nanmean(acg_allcat(goodUnits_allcat &celltype_allcat==4,501:end)),'FaceColor',rgb('HotPink'));
xlim([0 0.5])
ylim([0 50])
%ylim([0 150])
xlabel('Time (s)')
set(gca,'YTickLabel',[]);
makepretty;

%cv2, f.r:distribution plots 

msnB = nanmean(permute(squeeze(nanmax(mua_allcat_stimalign_exp_pad2(:, 1:18,goodUnits_allcat & celltype_allcat==1))), [2, 1]),2); 
fsiB = nanmean(permute(squeeze(nanmax(mua_allcat_stimalign_exp_pad2(:, 1:18,goodUnits_allcat & celltype_allcat==2))), [2, 1]),2); 
tanB = nanmean(permute(squeeze(nanmax(mua_allcat_stimalign_exp_pad2(:, 1:18,goodUnits_allcat & tanMethod1))), [2, 1]),2); 
uinB = nanmean(permute(squeeze(nanmax(mua_allcat_stimalign_exp_pad2(:, 1:18,goodUnits_allcat & celltype_allcat==4))), [2, 1]),2); 
figure(); 
distributionPlot([msnB;fsiB;tanB;uinB], 'addSpread', true, 'showMM', false,'groups', [ones(size(msnB,1),1); ones(size(fsiB,1),1).*2;...
    ones(size(tanB,1),1).*3;ones(size(uinB,1),1).*4])

figure();
boxplot([fsiB; msnB; tanB; uinB], [ones(size(fsiB,1),1); ones(size(msnB,1),1).*2; ...
    ones(size(tanB,1),1).*3; ones(size(uinB,1),1).*4],'Notch', 'on', 'Symbol', 'wx', 'labels', { 'FSI','MSN', 'TAN', 'UIN'},'Colors', ...
    [ rgb('RoyalBlue');rgb('DarkRed'); rgb('Lime'); rgb('HotPink')])
[p12, h12]=ranksum(msnB, fsiB);
[p34, h34]=ranksum(tanB, uinB);
[p24, h24]=ranksum(fsiB, uinB);
[p14, h14]=ranksum(msnB, uinB);
ylim([0 100])
siglineJF([ 1 2], strcat('p < 0.0001'))
siglineJF([3 4], strcat('p < 0.0001'))
ylim([0 120])
siglineJF([2 4], strcat('p < 0.0001'))
ylim([0 140])
siglineJF([1 4], strcat('p < 0.0001'))
ylim([0 170])
ylabel('baseline f.r. (Hz)')
makepretty;

%% ~EXAMPLE CELLS: WAVEFORM AND ACG~ 

%% ~~MUA~~

%% ~PSTH~

plot_trials = cellfun(@(stim) ...
    stim >= 0 | stim < 0, trial_stim_allcat_exp, 'uni', false);
maxtrials=max( cellfun(@(x) sum(x), plot_trials));
mmm = mua_allcat_stimalign_exp;
mua_allcat_stimalign_exp_padfull = cell2mat(permute(cellfun(@(act, trials_1) ...
padarray(act([find(trials_1)], :, :), ...
[maxtrials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
mmm, plot_trials, 'uni', false), [2, 3, 1]));

goodUnitsRec = vertcat(goodUnits{:}); 
allGroupsRec = vertcat(allGroups{:});
depthsRec = vertcat(allDepths{:});

celltypeL = {'MSN', 'FSI', 'TAN', 'UIN'}; 

figure();
for iDepth = 1:3
    for iCellType = 1:4 
        
    fullIM=[];
    fullMove=[];
    fullContrast=[];
    
        for iRec =1:45
            str = ~isnan(depthsRec{iRec});
            thisGroup = allGroupsRec{iRec}; 
            if iCellType==3
                 theseNeurons =(( thisGroup(str)'==iDepth+3*(iCellType-1) )| (thisGroup(str)'==iDepth+15))& goodUnitsRec{iRec}==1;
            else
                 theseNeurons = thisGroup(str)'==iDepth+3*(iCellType-1)& goodUnitsRec{iRec}==1;
            end
           
            theseTrials = trial_stim_allcat_exp{iRec}>0 & trial_outcome_allcat_exp{iRec} ==1;
            thisMoveT = move_t_exp{iRec};
            thismua = mua_allcat_stimalign_exp{iRec};
            thisIM = (permute(squeeze(nanmean(thismua(theseTrials, :, theseNeurons),3)), [2, 1])) ;    
            fullIM=[fullIM, thisIM];
            fullMove=[fullMove; thisMoveT(theseTrials)];
            thisContrast = trial_stim_allcat_exp{iRec}; 
            fullContrast = [fullContrast; thisContrast(theseTrials)]; 
        end

        [sortedMove, sortedMoveIdx] = sort(fullMove); %sort by move time 
        [sortedCon, sortedConIdx] = sort(fullContrast(sortedMoveIdx)); %sort by contrast 
        
        IMfinal = permute(fullIM(:,sortedMoveIdx(sortedConIdx)), [2, 1]);
       % smooth_filt = [round(size(IMfinal,1)*0.01),1]; 
       smooth_filt = [50,1]; % (trials x frames)
        curr_plot_smooth = nanconv(IMfinal,ones(smooth_filt),'same');
        
        subplot(3, 5, iCellType + (iDepth-1)*5)
        imagesc(t, [],curr_plot_smooth);
        xlim([-0.5, 1])
        title(strcat(celltypeL{iCellType},', str ', num2str(iDepth)))
        if iDepth==1
            xlabel('time from stim (ms)')
            ylabel('trial #')
      
        else
            set(gca,'YTickLabel',[]);
        end
        
        hold on;
        line([0, 0], ylim, 'linestyle', '--', 'color', 'r');
        line([0.5, 0.5], ylim, 'linestyle', '--', 'color', rgb('Gold'));
        colormap(brewermap([],'Greys'));
        yl=ylim;
        line(sortedMove(sortedConIdx), yl(1)+0.5:yl(2)-0.5, 'linestyle', '--', 'color', rgb('Purple'));
        makepretty;
    end
end


figure();
for iDepth = 1:3
    for iCellType = 1:4 
    fullIM=[];
    fullMove=[];
    
    fullOutcome=[];
    fullDirection =[];
    
        for iRec =1:45
            str = ~isnan(depthsRec{iRec});
            thisGroup = allGroupsRec{iRec}; 
            if iCellType==3
                 theseNeurons =(( thisGroup(str)'==iDepth+3*(iCellType-1) )| (thisGroup(str)'==iDepth+15))& goodUnitsRec{iRec}==1;
            else
                 theseNeurons = thisGroup(str)'==iDepth+3*(iCellType-1)& goodUnitsRec{iRec}==1;
            end
            theseTrials = trial_outcome_allcat_exp{iRec} ==1 | trial_outcome_allcat_exp{iRec} ==-1;
            thisMoveT = move_t_exp{iRec}-t(outcome_idx_exp{iRec})';
            if ~isempty(find(thisMoveT>0))
                theseTrials = theseTrials(thisMoveT<=0);
                thisMoveT = thisMoveT(thisMoveT<=0);
                
            end
            thismua = mua_allcat_outcomealign_exp{iRec};
            thisIM = (permute(squeeze(nanmean(thismua(theseTrials, :, theseNeurons),3)), [2, 1])) ;    
            fullIM=[fullIM, thisIM];
            fullMove=[fullMove; thisMoveT(theseTrials)];
            %t(outcome_idx_exp)
            thisOutcome = trial_outcome_allcat_exp{iRec} ; 
            fullOutcome = [fullOutcome; thisOutcome(theseTrials)]; 
            
            thisDirection =  (trial_outcome_allcat_exp{iRec} == 1 &  trial_stim_allcat_exp{iRec}>0)...
                | (trial_outcome_allcat_exp{iRec} == -1 &  trial_stim_allcat_exp{iRec}<0);
            thisDirection(find(thisDirection))=2;
            fullDirection = [fullDirection; thisDirection(theseTrials)]; 
        end

         [sortedMove, sortedMoveIdx] = sort(fullMove); %sort by move time 
        [sortedOut, sortedOutIdx] = sort(fullOutcome(sortedMoveIdx)); %sort by outcometype
        [sortedDir, sortedDirIdx] = sort(fullDirection(sortedOutIdx(sortedMoveIdx))); %sort by movedir
        
        IMfinal = permute(fullIM(:,sortedMoveIdx(sortedOutIdx(sortedDirIdx))), [2, 1]);
        % smooth_filt = [round(size(IMfinal,1)*0.01),1]; 
        smooth_filt = [50,1]; % (trials x frames)
        curr_plot_smooth = nanconv(IMfinal,ones(smooth_filt),'same');
        
        subplot(3, 5, iCellType + (iDepth-1)*5)
        imagesc(t, [],curr_plot_smooth);
        xlim([-0.5, 1])
        title(strcat(celltypeL{iCellType},', str ', num2str(iDepth)))
        if iDepth==1
            xlabel('time from rew (ms)')
            ylabel('trial #')
        else
            set(gca,'YTickLabel',[]);
        end
        
        hold on;
        subplot(3, 5, iCellType + (iDepth-1)*5)
        set(gca, 'YLim', [0, size(curr_plot_smooth,1)])
        line([0, 0], [0, size(curr_plot_smooth,1)], 'linestyle', '--', 'color', rgb('LightBlue'));
        colormap(brewermap([],'Greys'));
        yl=[0, size(curr_plot_smooth,1)];
        line(sortedMove(sortedOutIdx(sortedDirIdx)), yl(1)+0.5:yl(2)-0.5, 'linestyle', '--', 'color', rgb('Purple'));
        makepretty;
    end
end

%% ~~SUA~~

%% ~GOOD LOOKING EXAMPLE CELLS~
clear all;
naiveComp = 0;
data_fn = 'E:\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld.mat';
loadingTRialData;

plot_trials = cellfun(@(stim, move) ...
        stim > 0 & move > 0.2 , trial_stim_allcat_exp, move_t_exp, ...
        'uni', false); %correct trials + late movement 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_stimalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_stimalign_exp, plot_trials, 'uni', false), [2, 3, 1]));


outcome_idx_exp = mat2cell(outcome_idx, use_split, 1);
mua_allcat_outcomealign_exp2 = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp2)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp2{curr_exp}, 1)
        mua_allcat_outcomealign_exp2{curr_exp}(curr_trial, :, :) = ...
            circshift(mua_allcat_outcomealign_exp2{curr_exp}(curr_trial, :, :), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples, 2);
    end
end
plot_trials = cellfun(@(out, out_idx) ...
        out == 1 & t(out_idx)' > 0.55 , trial_outcome_allcat_exp, outcome_idx_exp,...
        'uni', false); %correct trials, and go cue before reward
clearvars mua_allcat_outalign_exp_pad2
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_outalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_outcomealign_exp, plot_trials, 'uni', false), [2, 3, 1]));
    

move_idx_exp = mat2cell(move_idx, use_split, 1);
mua_allcat_movealign_exp2 = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp2)
    for curr_trial = 1:size(mua_allcat_movealign_exp2{curr_exp}, 1)
        mua_allcat_movealign_exp2{curr_exp}(curr_trial, :, :) = ...
            circshift(mua_allcat_movealign_exp2{curr_exp}(curr_trial, :, :), ...
            -move_idx_exp{curr_exp}(curr_trial)+leeway_samples, 2);
    end
end
%outcome_t_exp = t(outcome_idx_exp); 
plot_trials = cellfun(@(out, stim, out_idx, move_t ) ...
        ((out == 1 & stim >0) | (out ==-1&stim<0) )& t(out_idx)' - move_t > 0.2 & (move_t > 0.55 | move_t < 0.3)...
        , trial_outcome_allcat_exp,...
        trial_stim_allcat_exp,outcome_idx_exp,move_t_exp,...
        'uni', false);%contralateral movement , reward after 200ms of movement start, gocue before or 200ms after 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_movealign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_movealign_exp2, plot_trials, 'uni', false), [2, 3, 1]));
move_idx_exp = mat2cell(move_idx, use_split, 1);
mua_allcat_movealign_exp2 = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp2)
    for curr_trial = 1:size(mua_allcat_movealign_exp2{curr_exp}, 1)
        move_t_stimal{curr_exp}(curr_trial) = move_idx_exp{curr_exp}(curr_trial);
    end
end
move_t_stimal = cellfun(@transpose, move_t_stimal, 'UniformOutput', false);
move_t_stimal = move_t_stimal'; 
plot_trials = cellfun(@(stim, move) ...
        stim > 0 & move > 0.2 , trial_stim_allcat_exp, move_t_exp, ...
        'uni', false); %correct trials + late movement 
max_trials = max(cellfun(@(x) sum(x), plot_trials));

move_t_stimal_pad = cell2mat(permute(cellfun(@(move, trials_1) ...
    padarray(move([find(trials_1)]), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
   move_t_stimal, plot_trials, 'uni', false), [2, 3, 1]));

%move to rew, use move plot trials 
movetoout_idx_exp = mat2cell(outcome_idx, use_split, 1);
mua_allcat_outalign_exp2 = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outalign_exp2)
    for curr_trial = 1:size(mua_allcat_outalign_exp2{curr_exp}, 1)
        movetoout_t_stimal{curr_exp}(curr_trial) = movetoout_idx_exp{curr_exp}(curr_trial)- move_idx_exp{curr_exp}(curr_trial)+leeway_samples;
    end
end
movetoout_t_stimal = cellfun(@transpose, movetoout_t_stimal, 'UniformOutput', false);
movetoout_t_stimal = movetoout_t_stimal'; 
plot_trials = cellfun(@(out, stim, out_idx, move_t ) ...
        ((out == 1 & stim >0) | (out ==-1&stim<0) )& t(out_idx)' - move_t > 0.2 & (move_t > 0.55 | move_t < 0.3)...
        , trial_outcome_allcat_exp,...
        trial_stim_allcat_exp,outcome_idx_exp,move_t_exp,...
        'uni', false);%contralateral movement , reward after 200ms of movement start, gocue before or 200ms after 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
movetoout_t_stimal_pad = cell2mat(permute(cellfun(@(out, trials_1) ...
    padarray(out([find(trials_1)]), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
   movetoout_t_stimal, plot_trials, 'uni', false), [2, 3, 1]));

%move to go cue, use move plot trials 
movetogocue_idx_exp = mat2cell(outcome_idx, use_split, 1);
mua_allcat_outalign_exp2 = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outalign_exp2)
    for curr_trial = 1:size(mua_allcat_outalign_exp2{curr_exp}, 1)
        movetogocue_t_stimal{curr_exp}(curr_trial) = 0.5 - move_idx_exp{curr_exp}(curr_trial)+leeway_samples;
    end
end
movetogocue_t_stimal = cellfun(@transpose, movetogocue_t_stimal, 'UniformOutput', false);
movetogocue_t_stimal =movetogocue_t_stimal'; 
plot_trials = cellfun(@(out, stim, out_idx, move_t ) ...
        ((out == 1 & stim >0) | (out ==-1&stim<0) )& t(out_idx)' - move_t > 0.2 & (move_t > 0.55 | move_t < 0.3)...
        , trial_outcome_allcat_exp,...
        trial_stim_allcat_exp,outcome_idx_exp,move_t_exp,...
        'uni', false);%contralateral movement , reward after 200ms of movement start, gocue before or 200ms after  
max_trials = max(cellfun(@(x) sum(x), plot_trials));
movetogocue_t_stimal_pad = cell2mat(permute(cellfun(@(out, trials_1) ...
    padarray(out([find(trials_1)]), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
   movetoout_t_stimal, plot_trials, 'uni', false), [2, 3, 1]));

%rew to last move, use rew plot trials 
active_trials = any(abs(wheel_allrec{iRec}(:, t >= 0 & t <= 0.3)) > wheel_thresh, 2);
for iTrial = 1:size(wheel_allcat,1)
    if max(wheel_allcat(iTrial,:) > wheel_thresh)
        mE = find(wheel_allcat(iTrial,:) > wheel_thresh, 1, 'last');
        moveEnd(iTrial) = mE(end); 
    else
        moveEnd(iTrial) = NaN; 
    end
end
moveEnd = moveEnd(trial_outcome_allcat==1); %reward only trials 

%which cell in which rec
recording_neurons = find(celltype_allcat == 1); %recordings ordered by cell type
recordingsCidx = [diff(recording_neurons')~=1, true];%[diff(recording_neurons')~=1, true]
recording_neurons = recording_neurons(recordingsCidx);

theseNeu = zeros(size(depth_aligned_allcat, 1), 1);
for iRec = 1:45
    if iRec == 1
        theseNeu(1:recording_neurons(iRec)) = iRec;
    else
        theseNeu(recording_neurons(iRec-1)+1:recording_neurons(iRec)) = iRec;
    end
end

%which trials in which rec 
trials_rec = [1;cumsum(trials_recording)]; 

%units to plot
stimSpikes = mean(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, :)))); 
%tt =  goodUnits_allcat & stimSpikes'>10 & celltype_allcat == 1;
tt =  celltype_allcat==2 &domain_allcat==1 ;
cc=celltype_allcat(find(tt));
ccc=domain_aligned_allcat(find(tt));

celltype_labels={'MSN','FSI','TAN','UIN','CRAP','TAN'};
%cell string 
cell_string = arrayfun(@(x) ['Cell: ', num2str(x), ...
    ', Domain: ', num2str(ccc(x)), ...
    ', ', celltype_labels{cc(x)}], 1:size(ccc,1), 'uni', false);

cell_string = arrayfun(@(x) ['Recording', , 'Cell', ,    ' Domain: ', num2str(ccc(x)), ...
    ', ', celltype_labels{cc(x)}],  1:size(ccc,1), 'uni', false);

%dynamic plot 
dynamicCellPlot_rastermean_overlaid(squeeze(nanmean((mua_allcat_stimalign_exp_pad2(1:100, :, find(tt))))),...
    squeeze(nanmean((mua_allcat_movealign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanmean((mua_allcat_outalign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_stimalign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_movealign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_outalign_exp_pad2(:, :, find(tt))))),...
    mua_allcat_stimalign_exp_pad2(:, :, find(tt)), ...
    mua_allcat_movealign_exp_pad2(:, :, find(tt)), ...
    mua_allcat_outalign_exp_pad2(:, :, find(tt)), ...
    move_t_stimal_pad, ...
    movetoout_t_stimal_pad,...
    moveEnd,...
    trials_rec,...
    theseNeu(find(tt)), {cell_string{:}}, t )

%% ~~TANS: stim vs reward~~


%% above/below unity line contrasts ==1 
plot_trials = cellfun(@(stim, move) ...
        stim ==1 , trial_stim_allcat_exp, move_t_exp, ...
        'uni', false); %correct trials + late movement 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_stimalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_stimalign_exp, plot_trials, 'uni', false), [2, 3, 1]));

tt = (celltype_allcat == 3 | celltype_allcat == 6) & goodUnits_allcat;

tanStimR1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 1)));
tanOutR1 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 1)));
tanStimB1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 1)));
tanOutB1 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 1)));

tanStimR2 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 2)));
tanOutR2 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 2)));
tanStimB2 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 2)));
tanOutB2 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 2)));

tanStimR3 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 3)));
tanOutR3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 3)));
tanStimB3 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 3)));
tanOutB3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 3)));

figure();
scatter(nanmax((tanStimR1 - nanmean(tanStimB1))./nanstd(tanStimB1)), ...
    nanmax((tanOutR1 - nanmean(tanOutB1))./nanstd(tanOutB1)), [], rgb('Orange'));
hold on;

scatter(nanmax((tanStimR2 - nanmean(tanStimB2))./nanstd(tanStimB2)), ...
    nanmax((tanOutR2 - nanmean(tanOutB2))./nanstd(tanOutB2)), [], rgb('Blue'));
hold on;
scatter(nanmax((tanStimR3 - nanmean(tanStimB3))./nanstd(tanStimB3)), ...
    nanmax((tanOutR3 - nanmean(tanOutB3))./nanstd(tanOutB3)), [], rgb('Green'));
hold on;
line([-2:38], [-2:38], 'Color', 'r')
legend({'1', '2', '3'})
xlim([-2 38])
ylim([-2 38])
xlabel('std inc after stim')
ylabel('std inc after rew')
makeprettyLite;

tanStimR = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt )));
tanOutR = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt)));
tanStimB = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt )));
tanOutB = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt )));

figure();
scatter(nanmax((tanStimR - nanmean(tanStimB))./nanstd(tanStimB)), ...
    nanmax((tanOutR - nanmean(tanOutB))./nanstd(tanOutB)), [], rgb('Orange'));
hold on;


hold on;
line([-2:38], [-2:38], 'Color', 'r')

xlim([-2 38])
ylim([-2 38])
xlabel('std inc after stim')
ylabel('std inc after rew')
makeprettyLite;

%% get average activity(stim or reward) 
stimR = nanmax((tanStimR - nanmean(tanStimB))./nanstd(tanStimB));
outR=    nanmax((tanOutR - nanmean(tanOutB))./nanstd(tanOutB));

stimActTANs = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:36, tt )));
stimActTANs = stimActTANs(:,stimR>outR);

outActTANs = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:36, tt )));
outActTANs = outActTANs(:,stimR<outR);

%% correlate
stimCorrs = corr(stimActTANs);
stimCorrs = triu(stimCorrs,1);
stimCorrs(stimCorrs==0)=[];

outCorrs = corr(outActTANs);
outCorrs = triu(outCorrs,1);
outCorrs(outCorrs==0)=[];

figure();
distributionPlot([stimCorrs'; outCorrs'],  'addSpread', true,'showMM', false, 'groups',...
    [ones(size(stimCorrs,2),1); ones(size(outCorrs,2),1).*2])
xticklabels({'stim', 'rew'})
ylabel('correlation rho')
makeprettyLite;

%% %% above/below unity line contrasts > 0 
plot_trials = cellfun(@(stim, move) ...
        stim >0 , trial_stim_allcat_exp, move_t_exp, ...
        'uni', false); %correct trials + late movement 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_stimalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_stimalign_exp, plot_trials, 'uni', false), [2, 3, 1]));

tt = (celltype_allcat == 3 | celltype_allcat == 6) & goodUnits_allcat;

tanStimR1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 1)));
tanOutR1 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 1)));
tanStimB1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 1)));
tanOutB1 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 1)));

tanStimR2 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 2)));
tanOutR2 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 2)));
tanStimB2 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 2)));
tanOutB2 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 2)));

tanStimR3 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 3)));
tanOutR3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt & domain_aligned_allcat == 3)));
tanStimB3 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 3)));
tanOutB3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt & domain_aligned_allcat == 3)));

figure();
scatter(nanmax((tanStimR1 - nanmean(tanStimB1))./nanstd(tanStimB1)), ...
    nanmax((tanOutR1 - nanmean(tanOutB1))./nanstd(tanOutB1)), [], rgb('Orange'));
hold on;

scatter(nanmax((tanStimR2 - nanmean(tanStimB2))./nanstd(tanStimB2)), ...
    nanmax((tanOutR2 - nanmean(tanOutB2))./nanstd(tanOutB2)), [], rgb('Blue'));
hold on;
scatter(nanmax((tanStimR3 - nanmean(tanStimB3))./nanstd(tanStimB3)), ...
    nanmax((tanOutR3 - nanmean(tanOutB3))./nanstd(tanOutB3)), [], rgb('Green'));
hold on;
line([-2:38], [-2:38], 'Color', 'r')
legend({'1', '2', '3'})
xlim([-2 38])
ylim([-2 38])
xlabel('std inc after stim')
ylabel('std inc after rew')
makeprettyLite;

tanStimR = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, tt )));
tanOutR = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:22, tt)));
tanStimB = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:18, tt )));
tanOutB = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:18, tt )));

figure();
scatter(nanmax((tanStimR - nanmean(tanStimB))./nanstd(tanStimB)), ...
    nanmax((tanOutR - nanmean(tanOutB))./nanstd(tanOutB)), [], rgb('Orange'));
hold on;


hold on;
line([-2:38], [-2:38], 'Color', 'r')

xlim([-2 38])
ylim([-2 38])
xlabel('std inc after stim')
ylabel('std inc after rew')
makeprettyLite;

%% get average activity(stim or reward) 
stimR = nanmax((tanStimR - nanmean(tanStimB))./nanstd(tanStimB));
outR=    nanmax((tanOutR - nanmean(tanOutB))./nanstd(tanOutB));

stimActTANs = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:32, tt )));
stimActTANs = stimActTANs(:,stimR>outR);

outActTANs = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:32, tt )));
outActTANs = outActTANs(:,stimR<outR);

%% correlate
stimCorrs = corr(stimActTANs);
stimCorrs = triu(stimCorrs,1);
stimCorrs(stimCorrs==0)=[];

outCorrs = corr(outActTANs);
outCorrs = triu(outCorrs,1);
outCorrs(outCorrs==0)=[];

figure();
distributionPlot([stimCorrs'; outCorrs'],  'addSpread', true,'showMM', false, 'groups',...
    [ones(size(stimCorrs,2),1); ones(size(outCorrs,2),1).*2])
xticklabels({'stim', 'rew'})
ylabel('correlation rho')
makeprettyLite;
%% correlate by rec/animal -> too few TANs for this..

%% plot average traces TANs 
stimActTANs = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 11:32, tt )));
stimActTANs = stimActTANs(:,stimR>outR);

outActTANs = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 11:32, tt )));
outActTANs = outActTANs(:,stimR<outR);

figure(); 
plot(t(11:32),stimActTANs)

figure(); 
plot(t(11:32),outActTANs)

%% correlate MSNs and FSIs by rec. cell-cell. -0.1-0.3 -allcells
plot_trials = cellfun(@(stim, move) ...
        stim >0 &move>0.3, trial_stim_allcat_exp, move_t_exp, ...
        'uni', false); %correct trials + late movement 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_stimalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_stimalign_exp, plot_trials, 'uni', false), [2, 3, 1]));

plot_trials = cellfun(@(out, stim, out_idx, move_t ) ...
        ((out == 1 & stim >0) | (out ==-1&stim<0) )& t(out_idx)' - move_t > 0.2 & (move_t > 0.55 | move_t < 0.3)...
        , trial_outcome_allcat_exp,...
        trial_stim_allcat_exp,outcome_idx_exp,move_t_exp,...
        'uni', false);%contralateral movement , reward after 200ms of movement start, gocue before or 200ms after 
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_movealign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_movealign_exp2, plot_trials, 'uni', false), [2, 3, 1]));

plot_trials = cellfun(@(out, out_idx) ...
        out == 1 & t(out_idx)' > 0.55 , trial_outcome_allcat_exp, outcome_idx_exp,...
        'uni', false); %correct trials, and go cue before reward
clearvars mua_allcat_outalign_exp_pad2
max_trials = max(cellfun(@(x) sum(x), plot_trials));
mua_allcat_outalign_exp_pad2 = cell2mat(permute(cellfun(@(act, trials_1) ...
    padarray(act([find(trials_1)], :, :), ...
    [max_trials - (sum(trials_1)), 0, 0], NaN, 'post'), ...
    mua_allcat_outcomealign_exp, plot_trials, 'uni', false), [2, 3, 1]));


tt = (celltype_allcat == 1) & goodUnits_allcat;
msnStr1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 1)));
msnStr2 = squeeze(nanmean(mua_allcat_movealign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 2)));
msnStr3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 3)));

msn1Corrs = corr(msnStr1);
msn1Corrs = triu(msn1Corrs,1);
msn1Corrs(msn1Corrs==0)=[];

msn2Corrs = corr(msnStr2);
msn2Corrs = triu(msn2Corrs,1);
msn2Corrs(msn2Corrs==0)=[];

msn3Corrs = corr(msnStr3);
msn3Corrs = triu(msn3Corrs,1);
msn3Corrs(msn3Corrs==0)=[];



tt = (celltype_allcat == 2) & goodUnits_allcat;
fsiStr1 = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 1)));
fsiStr2 = squeeze(nanmean(mua_allcat_movealign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 2)));
fsiStr3 = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, 19:29, tt & domain_aligned_allcat == 3)));

fsi1Corrs = corr(fsiStr1);
fsi1Corrs = triu(fsi1Corrs,1);
fsi1Corrs(fsi1Corrs==0)=[];

fsi2Corrs = corr(fsiStr2);
fsi2Corrs = triu(fsi2Corrs,1);
fsi2Corrs(fsi2Corrs==0)=[];

fsi3Corrs = corr(fsiStr3);
fsi3Corrs = triu(fsi3Corrs,1);
fsi3Corrs(fsi3Corrs==0)=[];

figure();
distributionPlot([msn1Corrs'; msn2Corrs'; msn3Corrs'; fsi1Corrs'; fsi2Corrs'; fsi3Corrs'],  'showMM', false, 'groups',...
    [ones(size(msn1Corrs,2),1); ones(size(msn2Corrs,2),1).*2; ones(size(msn3Corrs,2),1).*3; ones(size(fsi1Corrs,2),1).*4; ...
    ones(size(fsi2Corrs,2),1).*5; ones(size(fsi3Corrs,2),1).*6 ])
xticklabels({'str1 MSN', 'str2MSN', 'str3MSN', 'str1FSI', 'str2FSI', 'str3FSI'})
ylabel('correlation rho')
makeprettyLite;

%% by rec MSNs and FSIs 
goodUnitsRec = vertcat(goodUnits{:}); 
allGroupsRec = vertcat(allGroups{:});
depthsRec = vertcat(allDepths{:});
outcomeRec = vertcat(outcome_all{:});

tStart = 19;
tStop = 30;
for iRec = 1:45
    %get stim r msn
            
    str = ~isnan(depthsRec{iRec});
    thisGroup = allGroupsRec{iRec}; 
    theseNeurons = thisGroup(str)'==1& goodUnitsRec{iRec}==1;

    theseTrials = trial_stim_allcat_exp{iRec}>0 & move_t_exp{iRec} >0.3; 
    thismua = mua_allcat_stimalign_exp{iRec};
    IMstr1MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if ~isempty(IMstr1MSN)
            msn1RecCorrs = corr(IMstr1MSN );
    msn1RecCorrs = triu(msn1RecCorrs,1);
    msn1RecCorrs(msn1RecCorrs==0)=[];
    
    allmsn1RecCorrs(iRec) = nanmean(msn1RecCorrs);
    else
        allmsn1RecCorrs(iRec) =NaN; 
    end

           
    %str1 fsi 
    theseNeurons = thisGroup(str)'==4& goodUnitsRec{iRec}==1;
    
    IMstr1FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    
    if ~isempty(IMstr1FSI)
            fsi1RecCorrs = corr(IMstr1FSI);
    fsi1RecCorrs = triu(fsi1RecCorrs,1);
    fsi1RecCorrs(fsi1RecCorrs==0)=[];
    
    allfsi1RecCorrs(iRec) = nanmean(fsi1RecCorrs);
    else
        allfsi1RecCorrs(iRec) = NaN; 
    end
    
    %get move r
    theseNeurons = thisGroup(str)'==2& goodUnitsRec{iRec}==1;

    theseTrials = (trial_stim_allcat_exp{iRec}>0 & trial_outcome_allcat_exp{iRec} ==1) | (trial_stim_allcat_exp{iRec}<0 & trial_outcome_allcat_exp{iRec} ==-1); 
    thismua = mua_allcat_movealign_exp2{iRec};
    IMstr2MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
   
     if ~isempty(IMstr2MSN)
            msn2RecCorrs = corr(IMstr2MSN );
    msn2RecCorrs = triu(msn2RecCorrs,1);
    msn2RecCorrs(msn2RecCorrs==0)=[];
    
    allmsn2RecCorrs(iRec) = nanmean(msn2RecCorrs);
    else
        allmsn2RecCorrs(iRec) = NaN; 
     end
    
           
    %str2 fsi 
    theseNeurons = thisGroup(str)'==5& goodUnitsRec{iRec}==1;
    
    IMstr2FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
  
    
     if ~isempty(IMstr2FSI)
         fsi2RecCorrs = corr(IMstr2FSI );
    fsi2RecCorrs = triu(fsi2RecCorrs,1);
    fsi2RecCorrs(fsi2RecCorrs==0)=[];
    
    allfsi2RecCorrs(iRec) = nanmean(fsi2RecCorrs);
    else
       allfsi2RecCorrs(iRec) = NaN; 
     end
    
    %get rew r
    theseNeurons = thisGroup(str)'==3& goodUnitsRec{iRec}==1;

    theseTrials = trial_outcome_allcat_exp{iRec} ==1; 
    thismua = mua_allcat_outcomealign_exp{iRec};
    IMstr3MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  

    
     if ~isempty(IMstr3MSN)
          msn3RecCorrs = corr(IMstr3MSN );
    msn3RecCorrs = triu(msn3RecCorrs,1);
    msn3RecCorrs(msn3RecCorrs==0)=[];
    
    allmsn3RecCorrs(iRec) = nanmean(msn3RecCorrs);
    else
       allmsn3RecCorrs(iRec) = NaN; 
     end
           
    %str3 fsi 
    theseNeurons = thisGroup(str)'==6& goodUnitsRec{iRec}==1;
    
    IMstr3FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    
    if ~isempty(IMstr3FSI)
         fsi3RecCorrs = corr(IMstr3FSI);
    fsi3RecCorrs = triu(fsi3RecCorrs,1);
    fsi3RecCorrs(fsi3RecCorrs==0)=[];
    
    allfsi3RecCorrs(iRec) = nanmean(fsi3RecCorrs);
    else
       allfsi3RecCorrs(iRec) =  NaN; 
    end
     
end

%average corr for each and std. 
figure();
e=errorbar([1:2],[nanmean(allmsn1RecCorrs), nanmean(allfsi1RecCorrs)],[nanstd(allmsn1RecCorrs), nanstd(allfsi1RecCorrs)],'LineWidth', 2);
e.Color = 'red';
xlim([0.5 2.5])
hold on;
makepretty;
e2=errorbar([1:2],[nanmean(allmsn2RecCorrs), nanmean(allfsi2RecCorrs)],[nanstd(allmsn2RecCorrs), nanstd(allfsi2RecCorrs)],'LineWidth', 2);
e2.Color = rgb('Purple');
makepretty;
e3=errorbar([1:2],[nanmean(allmsn3RecCorrs), nanmean(allfsi3RecCorrs)],[nanstd(allmsn3RecCorrs), nanstd(allfsi3RecCorrs)],'LineWidth', 2);
e3.Color = rgb('Blue');
makepretty;
ylabel('correlation rho')
xticklabels({'', 'MSN', '', 'FSI', ''})
makepretty;

%% TANs, FSIs, MSNs, inter-trial cell corrs-messy and not well written code
%get all traces
%and rec number


tStart = 19;
tStop = 36;
% tStart = 19; 
% tStop=26;
allIMfsi1 = [];
allRecfsi1 = [];
allIMfsi2 = [];
allRecfsi2 = [];
allIMfsi3 = [];
allRecfsi3 = [];

allIMmsn1 = [];
allRecmsn1 = [];
allIMmsn2 = [];
allRecmsn2 = [];
allIMmsn3 = [];
allRecmsn3 = [];

for iRec = 1:45
    %get stim r msn
            
    str = ~isnan(depthsRec{iRec});
    thisGroup = allGroupsRec{iRec}; 
    theseNeurons = thisGroup(str)'==1& goodUnitsRec{iRec}==1;

    theseTrials = trial_stim_allcat_exp{iRec}>0 & move_t_exp{iRec} >0.3; 
    thismua = mua_allcat_stimalign_exp{iRec};
    IMstr1MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if numel(find(theseNeurons))==1
        allIMmsn1 = [allIMmsn1, IMstr1MSN'];
        allRecmsn1 = [allRecmsn1; ones(size(IMstr1MSN,1),1).*iRec];
    else
        allIMmsn1 = [allIMmsn1, IMstr1MSN];
        allRecmsn1 = [allRecmsn1; ones(size(IMstr1MSN,2),1).*iRec];
    end
    %theseFRmsn1 = 
           
    %str1 fsi 
    theseNeurons = thisGroup(str)'==4& goodUnitsRec{iRec}==1;
    
    IMstr1FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if numel(find(theseNeurons))==1
        allIMfsi1 = [allIMfsi1, IMstr1FSI'];
        allRecfsi1 = [allRecfsi1; ones(size(IMstr1FSI,1),1).*iRec];
    else
        allIMfsi1 = [allIMfsi1, IMstr1FSI];
        allRecfsi1 = [allRecfsi1; ones(size(IMstr1FSI,2),1).*iRec];
    end
    
    
    %get move r
    theseNeurons = thisGroup(str)'==2& goodUnitsRec{iRec}==1;

    theseTrials = (trial_stim_allcat_exp{iRec}>0 & trial_outcome_allcat_exp{iRec} ==1) | (trial_stim_allcat_exp{iRec}<0 & trial_outcome_allcat_exp{iRec} ==-1); 
    thismua = mua_allcat_movealign_exp2{iRec};
    IMstr2MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if numel(find(theseNeurons))==1
        allIMmsn2 = [allIMmsn2, IMstr2MSN'];
        allRecmsn2 = [allRecmsn2; ones(size(IMstr2MSN,1),1).*iRec];
    else
        allIMmsn2 = [allIMmsn2, IMstr2MSN];
        allRecmsn2 = [allRecmsn2; ones(size(IMstr2MSN,2),1).*iRec];
    end
    
    
           
    %str2 fsi 
    theseNeurons = thisGroup(str)'==5& goodUnitsRec{iRec}==1;
    
    IMstr2FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
  
    
    if numel(find(theseNeurons))==1
        allIMfsi2 = [allIMfsi2, IMstr2FSI'];
        allRecfsi2 = [allRecfsi2; ones(size(IMstr2FSI,1),1).*iRec];
    else
        allIMfsi2 = [allIMfsi2, IMstr2FSI];
        allRecfsi2 = [allRecfsi2; ones(size(IMstr2FSI,2),1).*iRec];
    end
    
    
    %get rew r
    theseNeurons = thisGroup(str)'==3& goodUnitsRec{iRec}==1;

    theseTrials = trial_outcome_allcat_exp{iRec} ==1; 
    thismua = mua_allcat_outcomealign_exp{iRec};
    IMstr3MSN = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if numel(find(theseNeurons))==1
        allIMmsn3 = [allIMmsn3, IMstr3MSN'];
        allRecmsn3 = [allRecmsn3; ones(size(IMstr3MSN,1),1).*iRec];
    else
        allIMmsn3 = [allIMmsn3, IMstr3MSN];
        allRecmsn3 = [allRecmsn3; ones(size(IMstr3MSN,2),1).*iRec];
    end
    
           
    %str3 fsi 
    theseNeurons = thisGroup(str)'==6& goodUnitsRec{iRec}==1;
    
    IMstr3FSI = (squeeze(nanmean(thismua(theseTrials, tStart:tStop, theseNeurons)))) ;  
    if numel(find(theseNeurons))==1
        allIMfsi3 = [allIMfsi3, IMstr3FSI'];
         allRecfsi3 = [allRecfsi3; ones(size(IMstr3FSI,1),1).*iRec];
    else
        allIMfsi3 = [allIMfsi3, IMstr3FSI];
         allRecfsi3 = [allRecfsi3; ones(size(IMstr3FSI,2),1).*iRec];
    end
    
   
    
    clearvars IMstr3FSI IMstr2FSI IMstr1FSI IMstr2MSN IMstr1MSN IMstr3MSN
end

%cycle trhough again, get cell corrs for cells not in rec 
size(allRecmsn1)
size(allIMmsn1)
% for iRec = 1:45 
%     thesecorrs = corr(allIMmsn1(:, allRecmsn1==iRec), allIMmsn1(:, allRecmsn1~=iRec));
%     stimCorrs = triu(stimCorrs,1);
%     stimCorrs(stimCorrs==0)=[];
% end
clearvars FRminMSN1 FRminMSN2 FRminMSN3 FRminFSI1  FRminFSI2  FRminFSI3 FRoutCorrs FRstimCorrs...
    thiscorrMSN1 
thiscorrMSN1 = corr(allIMmsn1, allIMmsn1);
thiscorrMSN1 = thiscorrMSN1.';
m  = tril(true(size(thiscorrMSN1)),-1);
thiscorrMSN1 = thiscorrMSN1(m).';

for iCell1=1:size(allIMmsn1,2)
    for iCell2 = 1:size(allIMmsn1,2)
       FRminMSN1(iCell1,iCell2) = min(nanmean(allIMmsn1(:,iCell1)),nanmean(allIMmsn1(:,iCell2)));
    end
end
FRminMSN1 = FRminMSN1.';
m  = tril(true(size(FRminMSN1)),-1);
FRminMSN1 = FRminMSN1(m).';


thiscorrMSN2 = corr(allIMmsn2, allIMmsn2);
thiscorrMSN2 = thiscorrMSN2.';
m  = tril(true(size(thiscorrMSN2)),-1);
thiscorrMSN2 = thiscorrMSN2(m).';
for iCell1=1:size(allIMmsn2,2)
    for iCell2 = 1:size(allIMmsn2,2)
       FRminMSN2(iCell1,iCell2) = min(nanmean(allIMmsn2(:,iCell1)),nanmean(allIMmsn2(:,iCell2)));
    end
end
FRminMSN2 = FRminMSN2.';
m  = tril(true(size(FRminMSN2)),-1);
FRminMSN2 = FRminMSN2(m).';

thiscorrMSN3 = corr(allIMmsn3, allIMmsn3);
thiscorrMSN3 = thiscorrMSN3.';
m  = tril(true(size(thiscorrMSN3)),-1);
thiscorrMSN3 = thiscorrMSN3(m).';
for iCell1=1:size(allIMmsn3,2)
    for iCell2 = 1:size(allIMmsn3,2)
       FRminMSN3(iCell1,iCell2) = min(nanmean(allIMmsn3(:,iCell1)),nanmean(allIMmsn3(:,iCell2)));
    end
end
FRminMSN3 = FRminMSN3.';
m  = tril(true(size(FRminMSN3)),-1);
FRminMSN3 = FRminMSN3(m).';



thiscorrFSI1 = corr(allIMfsi1, allIMfsi1);
thiscorrFSI1 = thiscorrFSI1.';
m  = tril(true(size(thiscorrFSI1)),-1);
thiscorrFSI1 = thiscorrFSI1(m).';

for iCell1=1:size(allIMfsi1,2)
    for iCell2 = 1:size(allIMfsi1,2)
       FRminFSI1(iCell1,iCell2) = min(nanmean(allIMfsi1(:,iCell1)),nanmean(allIMfsi1(:,iCell2)));
    end
end
FRminFSI1 = FRminFSI1.';
m  = tril(true(size(FRminFSI1)),-1);
FRminFSI1 = FRminFSI1(m).';

thiscorrFSI2 = corr(allIMfsi2, allIMfsi2);
thiscorrFSI2 = thiscorrFSI2.';
m  = tril(true(size(thiscorrFSI2)));
thiscorrFSI2 = thiscorrFSI2(m).';

for iCell1=1:size(allIMfsi2,2)
    for iCell2 = 1:size(allIMfsi2,2)
       FRminFSI2(iCell1,iCell2) = min(nanmean(allIMfsi2(:,iCell1)),nanmean(allIMfsi2(:,iCell2)));
    end
end
FRminFSI2 = FRminFSI2.';
m  = tril(true(size(FRminFSI2)),-1);
FRminFSI2 = FRminFSI2(m).';

thiscorrFSI3 = corr(allIMfsi3, allIMfsi3);
thiscorrFSI3 = thiscorrFSI3.';
m  = tril(true(size(thiscorrFSI3)),-1);
thiscorrFSI3 = thiscorrFSI3(m).';

for iCell1=1:size(allIMfsi3,2)
    for iCell2 = 1:size(allIMfsi3,2)
       FRminFSI3(iCell1,iCell2) = min(nanmean(allIMfsi3(:,iCell1)),nanmean(allIMfsi3(:,iCell2)));
    end
end
FRminFSI3 = FRminFSI3.';
m  = tril(true(size(FRminFSI3)),-1);
FRminFSI3 = FRminFSI3(m).';


% for iRec=1:45
%     thiscorrMSN1(allRecmsn1==iRec,allRecmsn1==iRec)=0;
%     thiscorrMSN2(allRecmsn2==iRec,allRecmsn2==iRec)=0;
%     thiscorrMSN3(allRecmsn3==iRec,allRecmsn3==iRec)=0;
%     
%     thiscorrFSI1(allRecfsi1==iRec,allRecfsi1==iRec)=0;
%     thiscorrFSI2(allRecfsi2==iRec,allRecfsi2==iRec)=0;
%     thiscorrFSI3(allRecfsi3==iRec,allRecfsi3==iRec)=0;
%     
% end

%tans
tt = (celltype_allcat == 3 | celltype_allcat== 6 ) & goodUnits_allcat;

stimActTANs = squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, tStart:tStop, tt )));
stimActTANs = stimActTANs(:,stimR>outR);

outActTANs = squeeze(nanmean(mua_allcat_outalign_exp_pad2(:, tStart:tStop, tt )));
outActTANs = outActTANs(:,stimR<outR);

stimCorrs = corr(stimActTANs);
stimCorrs = stimCorrs.';
m  = tril(true(size(stimCorrs)),-1);
stimCorrs = stimCorrs(m).';

for iCell1=1:size(stimActTANs,2)
    for iCell2 = 1:size(stimActTANs,2)
       FRstimCorrs(iCell1,iCell2) = min(nanmean(stimActTANs(:,iCell1)),nanmean(stimActTANs(:,iCell2)));
    end
end
FRstimCorrs = FRstimCorrs.';
m  = tril(true(size(FRstimCorrs)),-1);
FRstimCorrs = FRstimCorrs(m).';

outCorrs = corr(outActTANs);
outCorrs = outCorrs.';
m  = tril(true(size(outCorrs)),-1);
outCorrs = outCorrs(m).';

for iCell1=1:size(outActTANs,2)
    for iCell2 = 1:size(outActTANs,2)
       FRoutCorrs(iCell1,iCell2) = min(nanmean(outActTANs(:,iCell1)),nanmean(outActTANs(:,iCell2)));
    end
end
FRoutCorrs = FRoutCorrs.';
m  = tril(true(size(FRoutCorrs)),-1);
FRoutCorrs = FRoutCorrs(m).';

frSeq = 0:0.5:18;
clearvars msn1BIN fsi1BIN msn2BIN fsi2BIN msn3BIN FRoutCorrsBIN FRstimCorrsBIN
for iFR = 2:numel(frSeq)
    msn1BIN(iFR) = nanmean(thiscorrMSN1(FRminMSN1>frSeq(iFR-1) & FRminMSN1<frSeq(iFR))); 
    fsi1BIN(iFR) = nanmean(thiscorrFSI1(FRminFSI1>frSeq(iFR-1) & FRminFSI1<frSeq(iFR))); 
    msn2BIN(iFR) = nanmean(thiscorrMSN2(FRminMSN2>frSeq(iFR-1) & FRminMSN2<frSeq(iFR))); 
    fsi2BIN(iFR) = nanmean(thiscorrFSI2(FRminFSI2>frSeq(iFR-1) & FRminFSI2<frSeq(iFR))); 
    msn3BIN(iFR) = nanmean(thiscorrMSN3(FRminMSN3>frSeq(iFR-1) & FRminMSN3<frSeq(iFR))); 
    fsi3BIN(iFR) = nanmean(thiscorrFSI3(FRminFSI3>frSeq(iFR-1) & FRminFSI3<frSeq(iFR))); 
    FRoutCorrsBIN(iFR) = nanmean(outCorrs(FRoutCorrs>frSeq(iFR-1) & FRoutCorrs<frSeq(iFR))); 
    FRstimCorrsBIN(iFR) = nanmean(stimCorrs(FRstimCorrs>frSeq(iFR-1) & FRstimCorrs<frSeq(iFR))); 
end

figure();
plot(frSeq,msn1BIN,'Color',rgb('DarkRed'))
hold on; 
plot(frSeq,fsi1BIN,'Color',rgb('MidnightBlue'))
plot(frSeq,FRstimCorrsBIN,'Color', rgb('DarkGreen'))
plot(frSeq,msn2BIN,'Color',rgb('Red'))
plot(frSeq,fsi2BIN,'Color',rgb('Blue'))
plot(frSeq,msn3BIN,'Color',rgb('Salmon'))
plot(frSeq,fsi3BIN,'Color',rgb('SkyBlue'))
plot(frSeq,FRoutCorrsBIN,'Color', rgb('LimeGreen'))
xlabel('min pairwise f.r. (Hz)')
ylabel('mean correlation rho')
legend({'str 1MSN', 'str1 FSI', 'stim TAN', 'str2 MSN', 'str2 FSI', 'str3 MSN', 'str3 FSI' ,'rew TAN'})
makepretty;


%final plot 
figure();
distributionPlot([thiscorrMSN1'; thiscorrMSN2'; thiscorrMSN3';thiscorrFSI1'; thiscorrFSI2'; thiscorrFSI3';stimCorrs';outCorrs'],  'showMM',1, 'groups',...
    [ones(size(thiscorrMSN1,2),1); ones(size(thiscorrMSN2,2),1).*2; ones(size(thiscorrMSN3,2),1).*3; ones(size(thiscorrFSI1,2),1).*4; ...
    ones(size(thiscorrFSI2,2),1).*5; ones(size(thiscorrFSI3,2),1).*6 ; ones(size(stimCorrs,2),1).*7 ;ones(size(outCorrs,2),1).*8 ])
xticklabels({'str1 MSN', 'str2MSN', 'str3MSN', 'str1FSI', 'str2FSI', 'str3FSI', 'stim TAN', 'rew TAN'})
ylabel('correlation rho')
makeprettyLite;

%% spike time tiling coefficient - give up
for iCell1=1:size(stimActTANs,2)
    for iCell2 = 1:size(stimActTANs,2)
        sstc(iCell1,iCell2) = STTC(t(tStop) -t(tStart), t(2)-t(1), stimActTANs(:,iCell1), stimActTANs(:,iCell2),  0.01);
    end
end
sstc = triu(sstc,1);
sstc(sstc==0)=[];

for iCell1=1:size(outActTANs,2)
    for iCell2 = 1:size(outActTANs,2)
        sstcO(iCell1,iCell2) = STTC(t(tStop) -t(tStart), t(2)-t(1), outActTANs(:,iCell1), outActTANs(:,iCell2),  0.01);
    end
end
sstcO = triu(sstcO,1);
sstcO(sstcO==0)=[];

figure();
distributionPlot([sstc'; sstcO'],  'showMM',1, 'groups',...
    [ones(size(sstc,2),1); ones(size(sstcO,2),1).*2])
xticklabels({'stim TAN', 'rew TAN'})
ylabel('correlation rho')
makeprettyLite;

%% firing rate 
figure(); 
ndhist(FRminMSN1, thiscorrMSN1)
xlabel('min pairwise f.r.')
ylabel('correlation rho')
xlim([0 1])

figure(); 
ndhist(FRminFSI1, thiscorrFSI1,'bins',3)
xlabel('min pairwise f.r.')
ylabel('correlation rho')

%binned line plots : for each firing rate, get median correlation 

frSeq = 0:0.5:18;
clearvars msn1BIN fsi1BIN msn2BIN fsi2BIN msn3BIN FRoutCorrsBIN FRstimCorrsBIN
for iFR = 2:numel(frSeq)
    msn1BIN(iFR) = nanmean(thiscorrMSN1(FRminMSN1>frSeq(iFR-1) & FRminMSN1<frSeq(iFR))); 
    fsi1BIN(iFR) = nanmean(thiscorrFSI1(FRminFSI1>frSeq(iFR-1) & FRminFSI1<frSeq(iFR))); 
    msn2BIN(iFR) = nanmean(thiscorrMSN2(FRminMSN2>frSeq(iFR-1) & FRminMSN2<frSeq(iFR))); 
    fsi2BIN(iFR) = nanmean(thiscorrFSI2(FRminFSI2>frSeq(iFR-1) & FRminFSI2<frSeq(iFR))); 
    msn3BIN(iFR) = nanmean(thiscorrMSN3(FRminMSN3>frSeq(iFR-1) & FRminMSN3<frSeq(iFR))); 
    fsi3BIN(iFR) = nanmean(thiscorrFSI3(FRminFSI3>frSeq(iFR-1) & FRminFSI3<frSeq(iFR))); 
    FRoutCorrsBIN(iFR) = nanmean(outCorrs(FRoutCorrs>frSeq(iFR-1) & FRoutCorrs<frSeq(iFR))); 
    FRstimCorrsBIN(iFR) = nanmean(stimCorrs(FRstimCorrs>frSeq(iFR-1) & FRstimCorrs<frSeq(iFR))); 
end

figure();
plot(frSeq,msn1BIN,'Color',rgb('DarkRed'))
hold on; 
plot(frSeq,fsi1BIN,'Color',rgb('MidnightBlue'))
plot(frSeq,FRstimCorrsBIN,'Color', rgb('DarkGreen'))
plot(frSeq,msn2BIN,'Color',rgb('Red'))
plot(frSeq,fsi2BIN,'Color',rgb('Blue'))
plot(frSeq,msn3BIN,'Color',rgb('Salmon'))
plot(frSeq,fsi3BIN,'Color',rgb('SkyBlue'))
plot(frSeq,FRoutCorrsBIN,'Color', rgb('LimeGreen'))
xlabel('min pairwise f.r. (Hz)')
ylabel('mean correlation rho')
legend({'str 1MSN', 'str1 FSI', 'stim TAN', 'str2 MSN', 'str2 FSI', 'str3 MSN', 'str3 FSI' ,'rew TAN'})
makepretty;

%% |RMSQ 


%% example cell;s
tt = domain_allcat==1 |  domain_allcat==2 | domain_allcat==3;
cc=celltype_allcat(find(tt));
ccc=domain_aligned_allcat(find(tt));

neuAnimal = [];
cellNum=[];
fullRec=[];
neuRec=[];
recC=1;
cellNum2 = {};
for iAnimal =1:6
    thisN = cellfun(@transpose, allNeurons{iAnimal}, 'UniformOutput', false);
    thisM = vertcat(thisN{:});
    neuAnimal = [neuAnimal; ones(size(thisM,1),1).*iAnimal];
    for iRec = 1:size(thisN,1)
        neuRec = [neuRec; ones(size(thisN{iRec},1),1).*iRec];
        cellNum= [cellNum, 1:size(thisN{iRec},1)];
        cellNum2{recC}=  1:size(thisN{iRec},1);
        if iAnimal == 1 && iRec == 1
            oo=ones(size(thisN{iRec},1),1);
        else
            oo=ones(size(thisN{iRec},1),1)+max(fullRec);
        end
        fullRec = [fullRec; oo];
        recC = recC+1;
    end
end
%animals
tt =goodUnits_allcat & ( celltype_allcat==2 | celltype_allcat==2)&domain_allcat==1;
cc=celltype_allcat(find(tt));
ccc=domain_aligned_allcat(find(tt));
neuAcc = neuAnimal(find(tt)); 
neuRcc = neuRec(find(tt));
neuCcc = cellNum(find(tt)); 
cell_string = arrayfun(@(x) [  'Cell: ', num2str(neuCcc(x)),  ' Domain: ', num2str(ccc(x)), ...
    ', ', celltype_labels{cc(x)}, 'animal: ', animals{neuAcc(x)}, 'rec: ', num2str(neuRcc(x))],  1:size(ccc,1), 'uni', false);

%dynamic plot 
dynamicCellPlot_rastermean_overlaid(squeeze(nanmean((mua_allcat_stimalign_exp_pad2(1:100, :, find(tt))))),...
    squeeze(nanmean((mua_allcat_movealign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanmean((mua_allcat_outalign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_stimalign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_movealign_exp_pad2(:, :, find(tt))))),...
    squeeze(nanstd((mua_allcat_outalign_exp_pad2(:, :, find(tt))))),...
    mua_allcat_stimalign_exp_pad2(:, :, find(tt)), ...
    mua_allcat_movealign_exp_pad2(:, :, find(tt)), ...
    mua_allcat_outalign_exp_pad2(:, :, find(tt)), ...
    move_t_stimal_pad, ...
    movetoout_t_stimal_pad,...
    moveEnd,...
    trials_rec,...
    theseNeu(find(tt)), cell_string, t , acg_allcat(find(tt),:), wv_allcat(find(tt),:))

stimSpikes = nanmax(squeeze(nanmean(mua_allcat_stimalign_exp_pad2(:, 19:22, :)))); 
tt = domain_allcat==1 |  domain_allcat==2 |  domain_allcat==3 ;
cc=celltype_allcat(find(tt));
ccc=domain_aligned_allcat(find(tt));
neuAcc = neuAnimal(find(tt)); 
neuRcc = neuRec(find(tt)); 
neuCcc = cellNum(find(tt)); 
fullrecCcc = fullRec(find(tt)); 
cell_string = arrayfun(@(x) [  'Cell: ', num2str(neuCcc(x)),  ' Domain: ', num2str(ccc(x)), ...
    ', ', celltype_labels{cc(x)}, 'animal: ', animals{neuAcc(x)}, 'rec: ', num2str(neuRcc(x))],  1:size(ccc,1), 'uni', false);

dynamicCellPlot_rastermean_overlaidRec(mua_allcat_stimalign_exp, mua_allcat_movealign_exp2, mua_allcat_outcomealign_exp, ...
    fullrecCcc,neuCcc, move_t_exp, {cell_string{:}}, t )

%% checking 
waveform_t = 1e3*((0:size(wv_allcat,2)-1)/30000);

[187, 1, 8;
91, 1, 2;
142, 1, 1;
126, 1, 1]; %MSNs

[412,6,1;
107,1,2;
128,1,3];%FSIs

[325,4,5;
482,3,2;
44,1,2];%UINs

467
463
366
495
491
297

[364,5,1;%wonky
310,2,2;%a bit wonky
517,5,4;%=/=
290,5,1;%wonky
480,3,2;% abit wonky
184,2,1];%a bit wonly 

figure(2);
thiscell =297;
thisanimal = 2;
thisrec = 1; 
find(allNeurons{thisanimal}{thisrec}==thiscell)
allNeurons{thisanimal}{thisrec}(thiscell)



thism = mua_all{thisanimal}{thisrec};
clf;
subplot(131)
plot(t,squeeze(mean(thism(:,:, thiscell))))
xlabel('time from stim (s)')
subplot(132)
title(['dur = ', num2str(templateDuration{thisanimal}{thisrec}(thiscell)), 'ms'])
hold on;
plot(waveform_t,wv{thisanimal}{thisrec}(thiscell,:))
subplot(133)
title(['dur = ', num2str(templateDuration{thisanimal}{thisrec}(thiscell)), 'ms'])
hold on;
plot(acg{thisanimal}{thisrec}(thiscell,:))
xlim([0 1000])

