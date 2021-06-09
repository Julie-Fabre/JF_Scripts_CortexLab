%all Units
iRegion = 4;
regions = {'CP', 'STN', 'GPe', 'SNr', 'GPi'};
locCat = cat(1, ephysData.location);
siteCat = cat(1, ephysData.site);
dayCat = cat(1, ephysData.date);
currDayCat = cat(1, ephysData.curr_day);
animalsCat = cat(1, ephysData.animal);


ind = find(ismember(locCat, regions{iRegion}));
%depthsCat = cat(1,ephysData(ind).templateDepths); 
rec = table;
rec.animals = animalsCat(ind,:);
rec.currDays = currDayCat(ind,:);
rec.sites = siteCat(ind);
rec.days = dayCat(ind,:);
clearvars depthsCat
for iind = 1:size(ind,1)
    if ~isempty(ephysData(ind(iind)).template_depths)
    depthsCat(iind,:) =3840-[min(ephysData(ind(iind)).template_depths), max(ephysData(ind(iind)).template_depths)];
    end
end
rec.depthsCat = round(depthsCat);
uniqueRec=unique(rec, 'rows');
