
clear all;

animal = 'JF020';
day = '2020-11-14';
experiment = 3;
lilrig_load_experimentJF;
lilrig_retinotopyJF;
set(gcf,'color','w');
saveas(gcf, ['//znas.cortexlab.net/Subjects' filesep animal filesep day filesep 'retinotopy_' animal '_' day])
%AP_reference_outline('grid')