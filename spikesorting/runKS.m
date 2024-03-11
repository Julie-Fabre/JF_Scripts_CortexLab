% JF_runKSandBombcell(animal, date, site, recording, rerunKS, rerunQM)

rerunQM = 1;
rerunKS = 1;

%JF_runKSandBombcell('JF115', '2024-02-15', 3, [],rerunKS,0)%error
%JF_runKSandBombcell('JF123', '2024-02-23', 2, [],rerunKS,rerunQM )%error

thisDay = '2024-02-28';
JF_runKSandBombcell('JF123', thisDay, 1, [], rerunKS, rerunQM)
JF_extract_sync({'JF123'});

%% chronic
rerunQM = 1;
rerunKS = 0;
% JF082
dates = {'2022-08-05', '2022-08-06', '2022-08-07', '2022-08-08', '2022-08-09', ...
    '2022-08-10', '2022-08-11', '2022-08-12', '2022-08-15', '2022-08-16', '2022-08-18', '2022-08-22'};
for iDate = 6:size(dates, 2)
    for iSite = 1:2
        JF_runKSandBombcell('JF082', dates{iDate}, iSite, [], rerunKS, rerunQM)
    end
end
%JF_extract_sync({'JF082'});


% JF067

dates = {'2022-01-25', '2022-01-26', '2022-01-27', '2022-01-28', '2022-01-30', ...
    '2022-01-31', '2022-02-01', '2022-02-02', '2022-02-03', '2022-02-04', '2022-02-05', ...
    '2022-02-06', '2022-02-07', '2022-02-08', '2022-02-09', '2022-02-10'};
for iDate = 1:size(dates, 2)
    for iSite = 1
        JF_runKSandBombcell('JF067', dates{iDate}, iSite, [], rerunKS, rerunQM)
    end
end
%JF_extract_sync({'JF067'});


% JF084
dates = {'2022-08-05', '2022-08-06', '2022-08-07', '2022-08-08', '2022-08-09', ...
    '2022-08-10', '2022-08-11', '2022-08-12', '2022-08-15', '2022-08-16', '2022-08-18', '2022-08-22'};
for iDate = 1:size(dates, 2)
    for iSite = 1:4
        JF_runKSandBombcell('JF084', dates{iDate}, iSite, [], rerunKS, rerunQM)
    end
end
%JF_extract_sync({'JF084'});

% JF078
dates = {'2022-05-21', '2022-05-22', '2022-05-23', '2022-05-25', '2022-05-26', ...
    '2022-05-28', '2022-05-29', '2022-05-30', '2022-05-31', '2022-06-01'};
for iDate = 1:size(dates, 2)
    for iSite = 1:3
        JF_runKSandBombcell('JF078', dates{iDate}, iSite, [], rerunKS, rerunQM)
    end
end
%JF_extract_sync({'JF078'});


um_pipelineJF_v2;
% re extract raw aveformns inclusing full ones? 
%% ongoing
rerunKS = 1;
%JF052, JF045, JF053, JF054

JF_runKSandBombcell('JF052', '2021-08-31', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF052', '2021-08-31', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF052', '2021-08-31', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF052', '2021-09-01', 1, [], rerunKS, 0)

JF_runKSandBombcell('JF053', '2021-09-27', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF053', '2021-09-27', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF053', '2021-09-27', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF053', '2021-09-27', 4, [], rerunKS, 0)

JF_runKSandBombcell('JF054', '2021-09-25', 1, [], rerunKS, 0) %error
JF_runKSandBombcell('JF054', '2021-09-25', 2, [], rerunKS, 0) %error
JF_runKSandBombcell('JF054', '2021-09-26', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF054', '2021-09-26', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF054', '2021-09-26', 3, [], rerunKS, 0) %error
JF_runKSandBombcell('JF054', '2021-09-26', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF054', '2021-09-26', 5, [], rerunKS, 0)
JF_runKSandBombcell('JF054', '2021-09-26', 6, [], rerunKS, 0)

%JF_runKSandBombcell('JF117', '2024-01-21', 2, [],rerunKS,0)

JF_extract_sync;

% JF089
%JF_runKSandBombcell('JF089', '2022-11-15', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-15', 2, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-16', 1, [],rerunKS,0)
JF_runKSandBombcell('JF089', '2022-11-16', 2, [], rerunKS, 0) %error
%JF_runKSandBombcell('JF089', '2022-11-17', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-17', 2, [],rerunKS,0)


% JF090
%JF_runKSandBombcell('JF090', '2022-11-15', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF090', '2022-11-16', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF090', '2022-11-17', 1, [],rerunKS,0)

% JF059
%JF_runKSandBombcell('JF059', '2022-01-13', 1, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-13', 2, [], rerunKS, 0) % error
JF_runKSandBombcell('JF059', '2022-01-14', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-14', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-14', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-14', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-14', 5, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-14', 6, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-15', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-15', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF059', '2022-01-15', 3, [], rerunKS, 0)

% JF062
JF_runKSandBombcell('JF062', '2021-11-05', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-05', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-06', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-07', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-07', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-07', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF062', '2021-11-07', 4, [], rerunKS, 0)

% JF063 - no histo
JF_runKSandBombcell('JF063', '2021-11-09', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-09', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-09', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-09', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-09', 5, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-09', 6, [], rerunKS, 0)


% JF065 - no histo
JF_runKSandBombcell('JF065', '2022-01-18', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-18', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-18', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-18', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-19', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-19', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-20', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-20', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-20', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-20', 4, [], rerunKS, 0)

% JF066 - no histo
JF_runKSandBombcell('JF066', '2022-02-03', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF066', '2022-02-03', 2, [], rerunKS, 0)

% JF068 - no histo
JF_runKSandBombcell('JF068', '2022-01-27', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF068', '2022-01-27', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF068', '2022-01-27', 3, [], rerunKS, 0)
JF_runKSandBombcell('JF068', '2022-01-27', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF068', '2022-01-28', 1, [], rerunKS, 0) % not sure abut the rest that day

%%
% JF086 - no histo
JF_runKSandBombcell('JF086', '2022-11-01', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF086', '2022-11-02', 1, [], rerunKS, 0)

% JF088
JF_runKSandBombcell('JF088', '2022-11-29', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF088', '2022-11-30', 1, [], rerunKS, 0)

% JF092
JF_runKSandBombcell('JF092', '2022-12-12', 1, [], rerunKS, 0)
JF_runKSandBombcell('JF092', '2022-12-12', 2, [], rerunKS, 0)


set(0, 'DefaultFigureVisible', 'on')

%%
JF_runKSandBombcell('JF063', '2021-11-10', 2, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-10', 4, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-10', 5, [], rerunKS, 0)
JF_runKSandBombcell('JF063', '2021-11-10', 6, [], rerunKS, 0)
JF_runKSandBombcell('JF065', '2022-01-18', 3, [], rerunKS, 0)
