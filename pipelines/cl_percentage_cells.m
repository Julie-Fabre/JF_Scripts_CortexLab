% from https://www.nature.com/articles/s41467-021-25436-3
% A neuronfs responsiveness to a stimulus was determined based on a collective 
% measure of the reliability of the neurons in a given field using time-shuffled 
% data. First, a neuron’s activity on each trial was circularly shuffled by 
% a random amount. Next, a reliability value was calculated using this shuffled 
% data. This was repeated 1000 times to yield a distribution of reliability 
% values, and the 99th percentile of this distribution was stored. This 99th 
% percentile threshold was found for every neuron. If a neuron’s actual 
% average reliability across sessions was statistically greater than the 
% average of these 99th percentile values (two-tailed one-sample t-test), 
% it was classified as responsive to the stimulus. 