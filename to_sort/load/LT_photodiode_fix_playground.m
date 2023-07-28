%% leo 's photodiode fix playground 
%% Fix the weird photodiode trace %%
clear; close all; clc
load('\\zaru.cortexlab.net\Subjects\JF086\2022-11-02\2\2022-11-02_2_JF086_Timeline.mat')
load('\\zaru.cortexlab.net\Subjects\JF086\2022-11-02\2\2022-11-02_2_JF086_Block.mat')
​
%% Load the raw photodiode trace alongside its timestamps
t_points = Timeline.rawDAQTimestamps;
pd_trace = Timeline.rawDAQData(:,14);
​
%% Find the obvious ON/OFF flip indices
% Find the photodiode trace differential and gather the indices when 
% stimulus was flipped on and off by finding the local maxima and minima
pd_diff = diff(pd_trace);
[~,on_ind] = findpeaks(pd_diff,'MinPeakHeight',0.2); 
[~,off_ind] = findpeaks(-pd_diff,'MinPeakHeight',0.2); 
​
% Remove the first and last ON blocks as they are unusual
on_ind2 = on_ind(2:end-1);
​
% The photodiode trace is lowpass filtered to later detect smaller ON
% blocks that can be confused with the similar amplitude noise
if max(diff(t_points) == min(diff(t_points)))
    disp('Consistent sampling interval at 0.001 s')
end
% Using max(diff(t_points) and min(diff(t_points)), it was found that
% the sampling rate was consistently 1000 Hz (0.001s sampling interval).
pd_filt = lowpass(pd_trace,1.5,1000,'Steepness',0.9);
% The passband frequency was set to 2 Hz as that was slightly greater than
% the minimum stimulus presentation interval
pd_filtdiff = diff(pd_filt);
​
%% Fix the unusual ON-OFF-OFF patterns
on_ind = zeros(size(off_ind)); i1 = 1;
for i2 = 1:numel(off_ind)-1
    if on_ind2(i1) > off_ind(i2) && on_ind2(i1) < off_ind(i2+1)
        on_ind(i2) = on_ind2(i1);
        i1 = i1 + 1;
    else
        t1 = round(off_ind(i2) + (off_ind(i2+1)-off_ind(i2))/4);
        t2 = round(off_ind(i2) + (off_ind(i2+1)-off_ind(i2))*3/4);
        % Unnoticed deflection should be ~1/2 between the two OFF
        [~,ind] = min(pd_filtdiff(t1:t2));
        on_ind(i2) = ind + t1 - 1;
    end
end
​
on_ind = on_ind(1:end-1);
​
clear i1 i2 ind on_ind2 t1 t2
​
%% Fix the ON periods that are too long
stim_max = max(block.events.stimITIsValues); % maximum stimulation interval
% Convert seconds unit to index units (sampling interval 0.001 s)
stim_max = stim_max/0.001;
​
% Flag the linear indices of on_ind that violates the stimulus period
flag = find(off_ind(2:end)-on_ind > stim_max + 200); % 200 units of freedom
​
% Find the missing ON/OFF timestamps within the long ON blocks
off_miss = zeros(size(flag));
on_miss = zeros(size(flag));
for i = 1:numel(flag)
    % Find the range in which the ON/OFF timestamps would lie
    t1 = on_ind(flag(i)); t2 = off_ind(flag(i)+1); % long ON block [t1,t2]
    off_t = round(t1 + (t2-t1)/6):round(t1 + (t2-t1)/2); % ~1/3 of long ON block
    on_t = round(t1 + (t2-t1)/2):round(t1 + (t2-t1)*5/6); % ~2/3 of long ON block
    % Find the local min/max in the lowpass filtered differential
    [~,localmin] = min(pd_filtdiff(off_t));
    [~,localmax] = max(pd_filtdiff(on_t));
    % Store the ON/OFF timestamps
    off_miss(i) = localmin + off_t(1);
    on_miss(i) = localmax + on_t(1);
end
​
% Create ON/OFF index arrays that includes the missing ON blocks
on_indfin = zeros(numel(on_ind)+numel(flag),1);
off_indfin = zeros(numel(off_ind)+numel(flag),1);
​
% Insert the missing ON/OFF timestamps after the flagged ON/OFF timestamps
for i = 1:numel(flag)
    if i == 1
        on_indfin(1:flag(i)) = on_ind(1:flag(i));
        off_indfin(1:flag(i)) = off_ind(1:flag(i));
    else
        on_indfin(flag(i-1)+i:flag(i)+i-1) = on_ind(flag(i-1)+1:flag(i));
        off_indfin(flag(i-1)+i:flag(i)+i-1) = off_ind(flag(i-1)+1:flag(i));
    end
    on_indfin(flag(i)+i) = on_miss(i);
    off_indfin(flag(i)+i) = off_miss(i);
    if i == numel(flag)
        on_indfin(flag(i)+i+1:end) = on_ind(flag(i)+1:end);
        off_indfin(flag(i)+i+1:end) = off_ind(flag(i)+1:end);
    end
end
​
clear stim_max i on_t off_t t1 t2 block flag localmax localmin
clear off_ind off_miss on_ind on_miss
​
%% Create a boolean vector representing the photodiode ON/OFF periods
pd_btrace = logical(zeros(size(pd_trace)));
for i = 1:numel(on_indfin)
    pd_btrace(on_indfin(i):off_indfin(i+1)) = 1;
end
​
figure();
title('Abnormal Photodiode v Rescued Stim Period'); hold on
plot(t_points,pd_trace); hold on
plot(t_points,pd_btrace*5.5); hold on
xlabel('Time (s)'); ylabel('Photodiode Output');
​
​
clear i on_indfin off_indfin
clear pd_diff pd_filt pd_filtdiff Timeline