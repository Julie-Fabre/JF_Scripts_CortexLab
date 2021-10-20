
%% get LFP from 2.0 probes
function lfp = JF_get_NPX2_LFP(ephys_path, start_time, stop_time)
    myPaths;
    binFile = dir([ephys_path(1:end-15) ephys_path(end-4:end)  filesep '*.meta']);
    if isempty(binFile)
                binFile = dir([ephys_path(1:end-15) ephys_path(end-4:end)  filesep 'experiment1/*.meta']);
            end
    meta = ReadMeta_GLX(binFile.name, binFile.folder);
    ephysAPfile = [ephys_path(1:end-15) ephys_path(end-4:end)];
    if contains(meta.imRoFile, 'NPtype24_hStripe_botRow0_ref0.imro') %2.0 4 shank, bottom stripe
        chanMapData = load([dropboxPath, filesep, 'MATLAB/JF_scripts_CortexLab/kilosort/chanMapNP2_4Shank_bottRow_flipper.mat']);
        chanMap = chanMapData.chanMap;
        connected = chanMapData.connected;
    else %1.0 bottom row
        chanMapData = load([dropboxPath, filesep, 'MATLAB/JF_scripts_CortexLab/kilosort/chanMapNP2_1Shank_flipper.mat']);
        chanMap = chanMapData.chanMap;
        connected = chanMapData.connected;
    end
    ichn = chanMap(connected > 0);
    ichn = sort(ichn); %
    ops.numChannels = 385; % to edit for sure;
    ops.apSampleRate = 30000;
    ops.fsubsamp = 100;
    fsubsamp = ops.fsubsamp;
    ops.ichn = ichn;
    chunkSize = ops.apSampleRate * fsubsamp; % this number MUST be a multiple of fsubsamp
    chunkSize_t = 1 / chunkSize;
    AP_filename = dir([ephysAPfile filesep '*.bin']);
    nSamps = AP_filename.bytes / 2 / ops.numChannels;
    nChunksTotal = ceil(nSamps/chunkSize);
    lfp = zeros(ceil(1/fsubsamp), size(ichn, 2));

    mmf = memmapfile([AP_filename.folder, filesep, AP_filename.name], 'Format', {'int16', [ops.numChannels, nSamps], 'x'});

    for chunkID = 1%:nChunksTotal


        start = chunkSize * (chunkID - 1) + 1;
        startdownsized = (start - 1) / 100 + 1;

        if chunkID == nChunksTotal
            stop = nSamps;
            stopdownsized = size(lfp, 1);
        else
            stop = chunkSize * chunkID;
            stopdownsized = stop / 100;
        end


        dat = mmf.Data.x(1:384, start:stop);
        dat = double(dat);

        dat = double(permute(mean(reshape(dat(ichn, :), [size(ichn), size(dat, 2)]), 1), [3, 2, 1]));
        dat(fsubsamp*ceil(size(dat, 1)/fsubsamp), :) = 0;
        mua0 = JF_computeLFP(dat, ops.apSampleRate, fsubsamp);

        lfp(startdownsized:stopdownsized, :) = mua0;

    end
end