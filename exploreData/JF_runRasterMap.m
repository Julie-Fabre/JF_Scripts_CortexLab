function [cell_traces, isort1, isort2] = JF_runRasterMap(spike_times, spike_templates, ops)

% sorts the matrix S (neurons by time) along the first axis
% ops.nC = 30, number of clusters to use 
% ops.iPC = 1:100, number of PCs to use 
% ops.isort = [], initial sorting, otherwise will be the top PC sort
% ops.useGPU = 0, whether to use the GPU
% ops.upsamp = 100, upsampling factor for the embedding position
% ops.sigUp = 1, % standard deviation for upsampling


% Parameters
% 
% Rastermap first takes the specified PCs of the data, and then embeds them into n_X clusters. It returns upsampled cluster identities (n_X x upsamp). Clusters are also computed across Y (n_Y) and smoothed, to help with fitting.
% 
%     n_components : int, optional (default: 2) dimension of the embedding space
%     n_X : int, optional (default: 40) size of the grid on which the Fourier modes are rasterized
%     nPC : nparray, int, optional (default: 400) how many of the top PCs to use during optimization
%     alpha : float, optional (default: 1.0) exponent of the power law enforced on component n as: 1/(K+n)^alpha
%     K : float, optional (default: 1.0) additive offset of the power law enforced on component n as: 1/(K+n)^alpha
%     init : initialization of algorithm (default: 'pca') can use 'pca', 'random', or a matrix n_samples x n_components
% 
% Outputs
% 
% Rastermap model has the following attributes after running 'fit':
% 
%     embedding : array-like, shape (n_samples, n_components) Stores the embedding vectors.
%     u,sv,v : singular value decomposition of data S, potentially with smoothing
%     isort1 : sorting along first dimension (n_samples) of matrix
%     cmap : correlation of each item with all locations in the embedding map (before upsampling)
%     A : PC coefficients of each Fourier mode


% remove zero rows 

% get unit matrix 
all_templates = unique(spike_templates);
thisTime = [min(spike_times), max(spike_times)];
raster_window = [thisTime(1),thisTime(1)+2000] ;
psth_bin_size = 0.1;
timeVector = thisTime(1):psth_bin_size:thisTime(1)+2000; % 2000 seconds 
cell_traces = zeros(length(all_templates), length(timeVector)-1);


for iUnit = 1:length(all_templates)

curr_raster_spike_times = spike_times(spike_templates == all_templates(iUnit));
curr_raster_spike_times(curr_raster_spike_times > timeVector(end))= [];
cell_traces(iUnit,:) = histcounts(curr_raster_spike_times, timeVector);

end

% remove no firing cells
%firing_rate = sum(cell_traces,2) ./1000;
%removeMe = firing_rate
% run rastermap
[isort1, isort2, Sm] = mapTmap(cell_traces, ops);
imagesc(timeVector, [], cell_traces(isort1,:))
xlabel('time (s)')
ylabel('neuron #')
makepretty
colorbar;
% dot for each spike ? only good cells? 
end
