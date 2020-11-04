function cluster_groups = readClusterGroupsJF(ephys_path)

cluster_filepattern = [ephys_path, 'cluster_group*'];
cluster_filedir = dir(cluster_filepattern);
if ~isempty(cluster_filedir)
    cluster_filename = [ephys_path, filesep, cluster_filedir.name];
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid, '%d%s', 'HeaderLines', 1);
    fclose(fid);
end

end