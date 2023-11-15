contra = 1;
center = 0;
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat');
        idx = 5;
        passive = 1;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
        idx = 2;
        passive = 0;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo3.mat');
        idx = 2;
        passive = 0;

    end
cl_plot_average_task(task_data, idx, contra, center, passive)
end