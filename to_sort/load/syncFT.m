function sync = syncFT(AP_filename, nChansTotal, outputDir) 

if size(AP_filename,2) > 1 && iscell(AP_filename)
    AP_filename = AP_filename{1};
end
      
d=dir(AP_filename);

nSamps = d.bytes/2/nChansTotal;
mmf = memmapfile(AP_filename,'Format',{'int16', [nChansTotal nSamps],'x'});

disp('extracting sync...'); 
sync = mmf.Data.x(385,:);    % save sync data    
try
    save(sprintf('%s//sync.mat', outputDir),'sync');
catch 
    try
        warning('saving to local disk')
        if ~isempty(outputDir)
            folderHierarchy = strsplit(outputDir, filesep);
            mkdir(fullfile('/media/julie/ExtraHD/', folderHierarchy{5:end}))
            save(sprintf('%s//sync.mat',fullfile('/media/julie/ExtraHD/', folderHierarchy{5:end})),'sync');
        else
            folderHierarchy = strsplit(AP_filename, filesep);
            mkdir(fullfile('/media/julie/ExtraHD/', folderHierarchy{5:end-1}))
            save(sprintf('%s//sync.mat',fullfile('/media/julie/ExtraHD/', folderHierarchy{5:end-1})),'sync');
        end
    catch
        warning('saving error')
    end
end
end 
