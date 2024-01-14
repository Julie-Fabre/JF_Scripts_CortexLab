
function cl_plot_PSTHs(contra, center)
keepUnits = [1,2];
keepVis = 1;
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive_full.mat');
        idx = 5;
        %task_data = load('/home/julie/Dropbox/MATLAB/task_data_passive2.mat');
        %idx = 4;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_gogogo.mat');
        idx = 2;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/task_data_goNogo.mat');
        idx = 2;
    end
cl_plot_average_task(task_data, idx, contra, center,  keepUnits, keepVis)
end
end