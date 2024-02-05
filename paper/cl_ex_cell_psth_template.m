



% 
% cl_plotExCell_psth(animal, experiment, iProbe, unit, align_type, group_type, ...
%     raster_window, psth_bin_size, plot_raster, plot_psth)

cl_plotExCell_psth('JF090', 1, 8, 3, 'stimOn_noMove', 2,...
    [-0.2, 0.6], 0.001, 1, 1, 'spatialFreq')

cl_plotExCell_psth('JF090', 1, 8, 3, 'stimOn_noMove', 3,...
    [-0.2, 0.6], 0.001, 1, 1, 'ori')

cl_plotExCell_psth('JF090', 2, 8, 3, 'stimOn_noMove', 1,...
    [-0.2, 0.6], 0.001, 1, 1, 'locations')

cl_plotExCell_psth('JF090', 3, 8, 3, 'stimOn_noMove', 1,...
    [-0.2, 0.6], 0.001, 1, 1, 'natImg')