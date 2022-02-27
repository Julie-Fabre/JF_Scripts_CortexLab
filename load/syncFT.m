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
        mkdir(['/media/experiment/HDD/' AP_filename '/'])
        save(sprintf('%s//sync.mat',['/media/experiment/HDD/' AP_filename '/']),'sync');
    catch
        warning('saving error')
    end
end
end 
