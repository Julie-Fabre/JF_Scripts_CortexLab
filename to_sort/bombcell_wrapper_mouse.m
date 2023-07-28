function bombcell_wrapper_mouse(animal, rerunKS,rerunBC)

protocol = ''; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);

for iExp = 1:length(experiments)
    site_dirs = dir(['/home/netshare/zaru', filesep, animal, filesep,  experiments(iExp).day, filesep, 'ephys', filesep, 'site*']);
    
    for iSite =1:size(site_dirs,1)
        this_site = str2num(strrep(site_dirs(iSite).name, 'site', ''));
        try
            JF_runKSandBombcell(animal,experiments(iExp).day,  this_site, [],rerunKS,rerunBC)%
        catch
            disp(['error ', animal, ', ', experiments(iExp).day, ', ' this_site])
            %keyboard;
        end
    end
end