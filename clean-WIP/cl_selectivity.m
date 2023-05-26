%% cl_selectivity 

% from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8114876/:
% Might a calculation that is more robust to trial-to-trial variability reduce 
% the sensitivity of measurements to inclusion criteria or CV? We recalculated 
% OSI, DSI and TF with cross-validation, using half of the trials to identify 
% the stimulus condition that evoked the largest mean responses (grating 
% direction and TF) and then calculated OSI, DSI and TF for these preferred 
% conditions from the other half of the trials. The overall effect of 
% including more neurons based on their CV on the cross-validated metrics 
% across different areas was similar to that on the non-cross-validated 
% metrics (Fig. 5). The notable difference is that the noisy neurons in 
% the lowest decile of robustness no longer have high DSI or OSI values, 
% but are shifted to much lower values (Fig. 5F,J). This difference is 
% also reflected in the fact that the overall curves are shifted to lower 
% values (compare Figs. 5E,I and Fig. 4E,I). Thus, while more statistically 
% robust metrics calculated through cross-validation likely better reflect 
% the true values of the population, they do not reduce the impact of 
% selection on those metrics.

% for each unit, get maximum response on half trials, and then calculate
% selectivity index based on this 
% matrix of neurons mean response to each condition. use absolute value
% (mean subtracted) 

% for each condition, measure 