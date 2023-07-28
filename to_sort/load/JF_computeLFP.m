
function lfp = JF_computeLFP(dat, Fs, fact)

filtDat = dat;


lowPassCutoff = 300; % Hz
[b1, a1] = butter(5, lowPassCutoff/Fs, 'low');
filtDat = filtfilt(b1,a1, filtDat);



NT = size(filtDat,1);
lfp = permute(mean(reshape(filtDat, fact, NT/fact, []), 1), [2 3 1]);

lfp = -lfp;
