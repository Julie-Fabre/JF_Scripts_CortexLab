
%% lab servers

if ispc % eww. 
    zaruPath = '\\zaru.cortexlab.net\Subjects'; %\\zaru.cortexlab.net\
    znasPath = '\\znas.cortexlab.net\Subjects';
    zserverPath = '\\zserver.cortexlab.net\Subjects'; %
    zserverLabPath = '\\znas.cortexlab.net\Lab'; %
    zinuPath = '\\zinu.cortexlab.net\Subjects';
    zortexPath = '\\zortex.cortexlab.net\Subjects';

    %% extra/backup drives
    dropboxPath = 'C:\Users\julie\Dropbox';
    %fabreServerPath = '/home/netshare/fabreServer/';

    extraHDPath = 'D:\';
    % cortexLabBackup1 = '/media/julie/CortexLabBU1/';
    % cortexLabBackup2 = '/media/julie/CortexLabBU2/';
    % cortexLabBackup3 = '/media/julie/CortexLabBU3/';
    localExtHdPath = 'F:\';

    %% specific paths
    csvPath = 'C:\Users\julie\Dropbox\PhD_summary\';
    % atlasBrainRegLocation = '/home/julie/.brainglobe/allen_mouse_25um_v1.2';
    %
    brainsawPath = 'U:\';
    allenAtlasPath = 'C:\Users\julie\Dropbox\Atlas\allenCCF';
    % regParamsPath = '/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/histology/';
    tempPath = 'C:\temp\';
    JF_Scripts_CortexLabPath = 'C:\Users\julie\Documents\MATLAB\JF_Scripts_Cortexlab\';
    brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives

else

    %% lab servers
    zaruPath = '/home/netshare/zaru/';
    znasPath = '/home/netshare/znas/';
    zserverPath = '/home/netshare/zserver-data/'; %
    zserverLabPath = '/home/netshare/zserver-lab/'; %
    zinuPath = '/home/netshare/zinu/';
    zortexPath = '/home/netshare/zortex';

    %% extra/backup drives
    dropboxPath = '/home/julie/Dropbox/';
    fabreServerPath = '/home/netshare/fabreServer/';

    extraHDPath = '/media/julie/Expansion';%'/media/julie/ExtraHD/';
    cortexLabBackup1 = '/media/julie/CortexLabBU1/';
    cortexLabBackup2 = '/media/julie/CortexLabBU2/';
    cortexLabBackup3 = '/media/julie/CortexLabBU3/';
    localExtHdPath = '/media/julie/Expansion';

    %% specific paths
    atlasBrainRegLocation = '/home/julie/.brainglobe/allen_mouse_25um_v1.2';
    csvPath = '/home/julie/Dropbox/PhD_summary/';
    brainsawPath = '/home/netshare/znas-brainsaw/';
    allenAtlasPath = '/home/julie/Dropbox/Atlas/allenCCF';
    regParamsPath = '/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/histology/';
    tempPath = '/home/julie/temp/';
    JF_Scripts_CortexLabPath = '/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/';
    brainglobeLocation = '/home/julie/.brainglobe/'; % where your brainglobe data lives
end