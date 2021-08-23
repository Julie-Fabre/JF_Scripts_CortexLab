function sync = syncFT(AP_filename, nChansTotal, outputDir) 
d=dir(AP_filename);
nSamps = d.bytes/2/nChansTotal;
mmf = memmapfile(AP_filename,'Format',{'int16', [nChansTotal nSamps],'x'});

disp('extracting sync...'); 
sync=mmf.Data.x(385,:);    % save sync data---    
try
    save(sprintf('%s//sync.mat', outputDir),'sync');
catch 
    warning('saving error')
end
end 
