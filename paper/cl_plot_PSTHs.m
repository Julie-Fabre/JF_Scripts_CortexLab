
function cl_plot_PSTHs(contra, center)
keepUnits = [1,2];
keepVis = 1;
for iTask = 1:3
    if iTask == 1 %passive
        task_data = load('/home/julie/Dropbox/MATLAB/naive_data5.mat');
        idx = 5;
    elseif iTask == 2
        task_data = load('/home/julie/Dropbox/MATLAB/gogogo_data2.mat');
        idx = 2;
    elseif iTask == 3
        task_data = load('/home/julie/Dropbox/MATLAB/goNogo_data2.mat');
        idx = 2;
    end
cl_plot_average_task(task_data, idx, contra, center,  keepUnits, keepVis)

cl_plot_average_task_diff(task_data, idx, contra, center,  keepUnits, keepVis)

end
end