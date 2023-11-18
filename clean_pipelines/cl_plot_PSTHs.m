contra = 1;
center = 0;
keepUnits = [1,2];
keepVis = 1;
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full2.mat');
        idx = 5;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo2.mat');
        idx = 2;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo2.mat');
        idx = 2;
    end
cl_plot_average_task(task_data, idx, contra, center,  keepUnits, keepVis)
end