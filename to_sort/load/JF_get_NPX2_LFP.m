
%% get LFP from 2.0 probes
function lfp = JF_get_NPX2_LFP(ephys_path, start_time, stop_time)
    cl_myPaths;
    binFile = dir([ephys_path(1:end-15) ephys_path(end-4:end)  filesep '*.meta']);
    ephysAPfile = [ephys_path(1:end-15) ephys_path(end-4:end)];
    if isempty(binFile)
       binFile = dir([ephys_path(1:end-15) ephys_path(end-4:end)  filesep 'experiment1/*.meta']);
       ephysAPfile = [ephys_path(1:end-15) ephys_path(end-4:end)  filesep 'experiment1/'];
    end
    if isempty(binFile)
       binFile = dir([ephys_path(1:end-23) ephys_path(end-12:end-7)  filesep '*.meta']);
       ephysAPfile = [ephys_path(1:end-23) ephys_path(end-12:end-7)  ];
    end
    if size(binFile,1) > 1
        binFile = binFile(1); %lfp one
    end
    meta = ReadMeta_GLX(binFile.name, binFile.folder);
    
    [chanMapFile, xcoords, ycoords] = SGLXMetaToCoords_JF(binFile.name, binFile.folder, '');
   
    ops.numChannels = length(xcoords)+1; % to edit for sure;
    ops.apSampleRate = 30000;
    fsubsamp = 100;
    chunkSize = ops.apSampleRate * fsubsamp; % this number MUST be a multiple of fsubsamp
    AP_filename = dir([ephysAPfile filesep '*.bin']);
    nSamps = AP_filename.bytes / 2 / ops.numChannels;
    nChunksTotal = ceil(nSamps/chunkSize);
    lfp = zeros(ceil(1/fsubsamp), ops.numChannels);

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

        %dat = double(permute(mean(reshape(dat(ops.numChannels-1, :), [ops.numChannels-1, size(dat, 2)]), 1), [3, 2, 1]));
        dat(fsubsamp*ceil(size(dat, 1)/fsubsamp), :) = 0;
        mua0 = JF_computeLFP(dat, ops.apSampleRate, fsubsamp);

        lfp(startdownsized:stopdownsized, :) = mua0;

    end
end