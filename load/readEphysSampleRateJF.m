function ephys_sample_rate = readEphysSampleRateJF(ephys_path)
% for open ephys
header_path = [ephys_path, 'dat_params.txt'];
header_fid = fopen(header_path);
header_info = textscan(header_fid, '%s %s', 'delimiter', {' = '});
fclose(header_fid);

header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

% Load spike data
if isfield(header, 'sample_rate')
    ephys_sample_rate = str2num(header.sample_rate);
elseif isfield(header, 'ap_sample_rate')
    ephys_sample_rate = str2num(header.ap_sample_rate);
end

end