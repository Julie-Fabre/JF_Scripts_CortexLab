%get all V1 experiments
expIDs = findAllenExperiments('injection', 'VISp', 'line', '0', 'primary', true);
proj = getProjectionDataFromExperiment(expIDs(1));
%plot proj to caudoputamen, projection_intensity, projection volume,
%structure_id 

[~, brain_outline] = plotBrainGrid([],[]);
hold on;
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])