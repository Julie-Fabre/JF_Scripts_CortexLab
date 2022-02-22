clearvars single_best_match

single_best_match = zeros(size(all_corrs,1),1);


excludeThis = [];
all_corrs_excluded = sum(all_corrs,3);

 for iMax = 1:size(all_corrs,1)
     %iMax=iMax+1
     lin_corrs = reshape(all_corrs_excluded,1,[]);
     lin_corrs = lin_corrs(~isnan(lin_corrs));
     [~,sorted_corrs] = sort(lin_corrs);
     theseMaxCorrs = lin_corrs(sorted_corrs);
     
     thisMaxCorr = theseMaxCorrs(end);
    
     [thisMaxX, thisMaxY] = find(all_corrs_excluded == thisMaxCorr);
     n = 0; 
     try 
    while abs(sum(template_positions(thisMaxX, :, 1)-template_positions(thisMaxY, :, 2), 2)) > 30
        thisMaxCorr = theseMaxCorrs(end -n);
        [thisMaxX, thisMaxY] = find(all_corrs_excluded == thisMaxCorr);
        n = n+1;
    end
     catch 
         disp ('no match')
         
     end 
     if length(thisMaxX)>1
         disp('warning: 2 identical values')
     end
     single_best_match(thisMaxX) = thisMaxY;
     disp(single_best_match(thisMaxX))
     all_corrs_excluded(:,thisMaxY) = NaN;% ensure you more than one n day unit can't be matched with the same n+1 day unit
     all_corrs_excluded(thisMaxX,:) = NaN;% ensure you more than one n+1 day unit can't be matched with the same n day unit
 end