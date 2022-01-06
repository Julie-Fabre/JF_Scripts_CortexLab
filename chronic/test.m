    for curr_template = 1:n_clu
        unique_clu2 = dat.d2.clu.ID;
        n_clu2 = size(unique_clu, 1);
       for curr_template2 = 1:n_clu2
           ss = smoothdata(all_ccg{1}(curr_template, :), 'gaussian', 40);
           ss2=smoothdata(all_ccg{2}(curr_template2, :), 'gaussian', 40);
           ccc = dtw(ss', ...
              ss2');
           all_corrs(curr_template, curr_template2, 1) = ccc;
           cc = dtw(waveforms{:, :, 1}(curr_template,:), waveforms{:, :, 2}(curr_template2,:));
           all_corrs(curr_template, curr_template2, 2) = cc;
           wv_duration1(curr_template) = find(waveforms{:, :, 1}(curr_template,:) == max(waveforms{:, :, 1}(curr_template,:) )) - ...
            find(waveforms{:, :, 1}(curr_template,:) == min(waveforms{:, :, 1}(curr_template,:) ));
        wv_duration2(curr_template2) = find(waveforms{:, :, 2}(curr_template2,:) == max(waveforms{:, :, 2}(curr_template2,:) )) - ...
            find(waveforms{:, :, 2}(curr_template2,:) == min(waveforms{:, :, 2}(curr_template2,:) ));
        acg_peak_end_ratio1(curr_template) = max(smoothdata(all_ccg{1}(curr_template, :), 'gaussian', 40))/...
            nanmean(all_ccg{1}(curr_template, :));
       acg_peak_end_ratio2(curr_template2) = max(smoothdata(all_ccg{2}(curr_template2, :), 'gaussian', 40))/...
            nanmean(all_ccg{2}(curr_template2, :));
        acg_pss1(curr_template) = find(smoothdata(all_ccg{1}(curr_template, :), 'gaussian', 40) ...
            >= nanmean(all_ccg{1}(curr_template, 450:500)),1,'first');
        acg_pss2(curr_template2) = find(smoothdata(all_ccg{2}(curr_template2, :), 'gaussian', 40) ...
            >= nanmean(all_ccg{2}(curr_template2, 450:500)),1,'first');
           
       end
       best_match(curr_template) = find(sum([squeeze(all_corrs(curr_template, :, :))],2) == ...
           max(sum([squeeze(all_corrs(curr_template, :, :))],2)));
       best_match1(curr_template) = find([squeeze(all_corrs(curr_template, :, 1))] == ...
           max([squeeze(all_corrs(curr_template, :, 1))]));
       best_match2(curr_template) = find([squeeze(all_corrs(curr_template, :, 2))] == ...
           max([squeeze(all_corrs(curr_template, :, 2))]));
       best_match3(curr_template) = find(sum([squeeze(all_corrs(curr_template, euclidian_dist<=50, :))],2) == ...
           max(sum([squeeze(all_corrs(curr_template, euclidian_dist<=50, :))],2))); 
       
        

       best_match4(curr_template) = find(sum([abs(wv_duration1 - wv_duration2);abs(acg_peak_end_ratio1 - acg_peak_end_ratio2); ...
           abs(acg_pss1 - acg_pss2)])==min(sum([abs(wv_duration1 - wv_duration2);abs(acg_peak_end_ratio1 - acg_peak_end_ratio2); ...
           abs(acg_pss1 - acg_pss2)])),1,'first');
    end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    