% JF_runKSandBombcell(animal, date, site, recording, rerunKS, rerunQM)

rerunQM = 1;
rerunKS = 1;

JF_runKSandBombcell('JF115', '2024-02-15', 3, [],rerunKS,0)%error


JF_runKSandBombcell('JF115', '2024-02-16', 1, [],rerunKS,0)
JF_runKSandBombcell('JF115', '2024-02-16', 2, [],rerunKS,0)
JF_runKSandBombcell('JF115', '2024-02-16', 3, [],rerunKS,0)


JF_extract_sync;

set(0,'DefaultFigureVisible','off')
%% ongoing
rerunKS = 1;
%JF052, JF045, JF053, JF054

JF_runKSandBombcell('JF052', '2021-08-31', 1, [],rerunKS,0)
JF_runKSandBombcell('JF052', '2021-08-31', 2, [],rerunKS,0)
JF_runKSandBombcell('JF052', '2021-08-31', 3, [],rerunKS,0)
JF_runKSandBombcell('JF052', '2021-09-01', 1, [],rerunKS,0)

JF_runKSandBombcell('JF053', '2021-09-27', 1, [],rerunKS,0)
JF_runKSandBombcell('JF053', '2021-09-27', 2, [],rerunKS,0)
JF_runKSandBombcell('JF053', '2021-09-27', 3, [],rerunKS,0)
JF_runKSandBombcell('JF053', '2021-09-27', 4, [],rerunKS,0)

JF_runKSandBombcell('JF054', '2021-09-25', 1, [],rerunKS,0)%error
JF_runKSandBombcell('JF054', '2021-09-25', 2, [],rerunKS,0)%error
JF_runKSandBombcell('JF054', '2021-09-26', 1, [],rerunKS,0)
JF_runKSandBombcell('JF054', '2021-09-26', 2, [],rerunKS,0)
JF_runKSandBombcell('JF054', '2021-09-26', 3, [],rerunKS,0)%error
JF_runKSandBombcell('JF054', '2021-09-26', 4, [],rerunKS,0)
JF_runKSandBombcell('JF054', '2021-09-26', 5, [],rerunKS,0)
JF_runKSandBombcell('JF054', '2021-09-26', 6, [],rerunKS,0)

%JF_runKSandBombcell('JF117', '2024-01-21', 2, [],rerunKS,0)

JF_extract_sync;

% JF089
%JF_runKSandBombcell('JF089', '2022-11-15', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-15', 2, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-16', 1, [],rerunKS,0)
JF_runKSandBombcell('JF089', '2022-11-16', 2, [],rerunKS,0)%error
%JF_runKSandBombcell('JF089', '2022-11-17', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF089', '2022-11-17', 2, [],rerunKS,0)


% JF090
%JF_runKSandBombcell('JF090', '2022-11-15', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF090', '2022-11-16', 1, [],rerunKS,0)
%JF_runKSandBombcell('JF090', '2022-11-17', 1, [],rerunKS,0)

% JF059
%JF_runKSandBombcell('JF059', '2022-01-13', 1, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-13', 2, [],rerunKS,0)% error
JF_runKSandBombcell('JF059', '2022-01-14', 1, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-14', 2, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-14', 3, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-14', 4, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-14', 5, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-14', 6, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-15', 1, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-15', 2, [],rerunKS,0)
JF_runKSandBombcell('JF059', '2022-01-15', 3, [],rerunKS,0)

% JF062
JF_runKSandBombcell('JF062', '2021-11-05', 1, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-05', 2, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-06', 1, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-07', 1, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-07', 2, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-07', 3, [],rerunKS,0)
JF_runKSandBombcell('JF062', '2021-11-07', 4, [],rerunKS,0)

% JF063 - no histo 
JF_runKSandBombcell('JF063', '2021-11-09', 1, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-09', 2, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-09', 3, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-09', 4, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-09', 5, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-09', 6, [],rerunKS,0)




% JF065 - no histo
JF_runKSandBombcell('JF065', '2022-01-18', 1, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-18', 2, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-18', 3, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-18', 4, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-19', 1, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-19', 2, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-20', 1, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-20', 2, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-20', 3, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-20', 4, [],rerunKS,0)

% JF066 - no histo 
JF_runKSandBombcell('JF066', '2022-02-03', 1, [],rerunKS,0)
JF_runKSandBombcell('JF066', '2022-02-03', 2, [],rerunKS,0)

% JF068 - no histo 
JF_runKSandBombcell('JF068', '2022-01-27', 1, [],rerunKS,0)
JF_runKSandBombcell('JF068', '2022-01-27', 2, [],rerunKS,0)
JF_runKSandBombcell('JF068', '2022-01-27', 3, [],rerunKS,0)
JF_runKSandBombcell('JF068', '2022-01-27', 4, [],rerunKS,0)
JF_runKSandBombcell('JF068', '2022-01-28', 1, [],rerunKS,0)% not sure abut the rest that day 

%%
% JF086 - no histo
JF_runKSandBombcell('JF086', '2022-11-01', 1, [],rerunKS,0)
JF_runKSandBombcell('JF086', '2022-11-02', 1, [],rerunKS,0)

% JF088
JF_runKSandBombcell('JF088', '2022-11-29', 1, [],rerunKS,0)
JF_runKSandBombcell('JF088', '2022-11-30', 1, [],rerunKS,0)

% JF092 
JF_runKSandBombcell('JF092', '2022-12-12', 1, [],rerunKS,0)
JF_runKSandBombcell('JF092', '2022-12-12', 2, [],rerunKS,0)



set(0,'DefaultFigureVisible','on')

%%
JF_runKSandBombcell('JF063', '2021-11-10', 2, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-10', 4, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-10', 5, [],rerunKS,0)
JF_runKSandBombcell('JF063', '2021-11-10', 6, [],rerunKS,0)
JF_runKSandBombcell('JF065', '2022-01-18', 3, [],rerunKS,0)