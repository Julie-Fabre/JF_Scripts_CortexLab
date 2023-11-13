close all;
clear all;
keepVis = 0;
keepVis2 = 0;
for ppp = 1:3
    %close all;
    for rrr = 1:3
        %close all;
        passiveYes = 0;
        if keepVis == 0 && passiveYes

            %% passive
            passive = 1;
            goNogo = 0;
            cl_plot_dprime;

            % plot d prime
            figure(100);
            for iRegion = 1:size(plot_regions, 2)
                for iPair = 1:3
                    subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)));
                    hold on;
                    % if 0, remove
                    %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
                    kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))));

                    histogram(d_prime{iRegion}(kp, iPair), [0:0.05:2], 'Normalization', 'probability', 'FaceAlpha', 0.2) % 1.4142   -1.4142 = one of stims has 0 spikes
                    %nanmedian(d_prime(kp,iPair,iRegion))
                    %nanmean(d_prime(kp,iPair,iRegion))
                    %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime{iRegion}(kp,iPair))))])
                    xlabel('d-prime')
                    ylabel('# of neurons')
                    xlim([0, 2])
                    %makepretty_lite;
                end
            end

            figure('name', 'passive');
            hold on;

            for iRegion = 1:3
                num_recs = size(unique(d_prime_session_num{iRegion}(:, 1)), 1);
                theseColors = [lines(num_recs)];
                theseColors = repelem(theseColors, 3, 1);

                kp1 = find(abs(d_prime{iRegion}(:, 1)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 1))) & ~isnan(abs(d_prime{iRegion}(:, 1))));
                kp2 = find(abs(d_prime{iRegion}(:, 2)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 2))) & ~isnan(abs(d_prime{iRegion}(:, 2))));
                kp3 = find(abs(d_prime{iRegion}(:, 3)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 3))) & ~isnan(abs(d_prime{iRegion}(:, 3))));

                subplot(3, 1, iRegion)
                violinplot([d_prime{iRegion}(kp1, 1); d_prime{iRegion}(kp2, 2); d_prime{iRegion}(kp3, 3)], ...
                    [d_prime_session_num{iRegion}(kp1, 1); d_prime_session_num{iRegion}(kp2, 2) + 0.2; d_prime_session_num{iRegion}(kp3, 3) + 0.4], ...
                    'ViolinColor', theseColors);
                ylim([0, 2.6])
                uu = unique(d_prime_session_num{iRegion});
                mouseys = [d_prime_animal_num{iRegion, ...
                    uu(~isnan(unique(d_prime_session_num{iRegion}(:, 1))), 1) ...
                    }];
                title(num2str(mouseys))

            end

            save('d_prime_session_fraction_passive.mat', 'd_prime_session_fraction')
            save('d_prime_session_fraction_passive_median.mat', 'd_prime_session_fraction_median')
        end

        %% gogogo
        passive = 0;
        goNogo = 1;
        cl_plot_dprime;
        %lmeTable.stimulus
        iRegion = rrr;
        iPair = ppp;

        lmeTable = table;
        %for iRegion = 1%:3
        kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))) & ~isnan(abs(d_prime{iRegion}(:, iPair))));
        uu = unique(d_prime_session_num{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{iRegion}(:, iPair))));
        mouseys = [d_prime_animal_num{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{iRegion}(:, iPair))), 1) ...
            }];

        lmeTable.dprime = d_prime{iRegion}(kp, iPair);
        lmeTable.sessionType = ones(size(d_prime{iRegion}(kp, iPair), 1), 1);
        lmeTable.session = d_prime_session_num{iRegion}(kp, iPair);
        mousee = d_prime_session_num{iRegion}(kp, iPair);
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        lmeTable.mouse = mousee;
        lmeTable.pss = pss{iRegion}(kp, iPair);
        lmeTable.tempDur = tempDur{iRegion}(kp, iPair);
        lmeTable.propISI = propISI{iRegion}(kp, iPair);
        lmeTable.fr = fr{iRegion}(kp, iPair);
        %end

        %iRegion = rrr;
        % if ppp ==1
        diff_go_nogo = table;
        diff_go_nogo.value = [d_prime{1}(:, 3) - d_prime{1}(:, 1); ...
            d_prime{2}(:, 3) - d_prime{2}(:, 1); ...
            d_prime{3}(:, 3) - d_prime{3}(:, 1)];
        diff_go_nogo.session = [d_prime_session_num{1}(:, 1); d_prime_session_num{2}(:, 1); d_prime_session_num{3}(:, 1)];
        diff_go_nogo.region_num = [ones(size(d_prime{1}(:, 3), 1), 1); ones(size(d_prime{2}(:, 3), 1), 1) .* 2; ones(size(d_prime{3}(:, 3), 1), 1) .* 3];
        %diff_go_nogo
        figure();
        violinplot(diff_go_nogo.value, diff_go_nogo.region_num)


        % plot d prime
        figure(100);
        for iRegion = 1:size(plot_regions, 2)
            for iPair = 1:3
                subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)));
                hold on;
                % if 0, remove
                %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
                kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))));

                histogram(d_prime{iRegion}(kp, iPair), [0:0.05:2], 'Normalization', 'probability', 'FaceAlpha', 0.2) % 1.4142   -1.4142 = one of stims has 0 spikes
                %nanmedian(d_prime(kp,iPair,iRegion))
                %nanmean(d_prime(kp,iPair,iRegion))
                %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime{iRegion}(kp,iPair))))])
                xlabel('d-prime')
                ylabel('# of neurons')
                xlim([0, 2])
                %makepretty_lite;
            end
        end

        figure('name', 'gonogo');
        hold on;
        for iRegion = 1:3
            num_recs = size(unique(d_prime_session_num{iRegion}(:, 1)), 1);
            theseColors = [lines(num_recs)];
            theseColors = repelem(theseColors, 3, 1);

            kp1 = find(abs(d_prime{iRegion}(:, 1)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 1))) & ~isnan(abs(d_prime{iRegion}(:, 1))));
            kp2 = find(abs(d_prime{iRegion}(:, 2)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 2))) & ~isnan(abs(d_prime{iRegion}(:, 2))));
            kp3 = find(abs(d_prime{iRegion}(:, 3)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 3))) & ~isnan(abs(d_prime{iRegion}(:, 3))));

            subplot(3, 1, iRegion)
            violinplot([d_prime{iRegion}(kp1, 1); d_prime{iRegion}(kp2, 2); d_prime{iRegion}(kp3, 3)], ...
                [d_prime_session_num{iRegion}(kp1, 1); d_prime_session_num{iRegion}(kp2, 2) + 0.2; d_prime_session_num{iRegion}(kp3, 3) + 0.4], ...
                'ViolinColor', theseColors);
            ylim([0, 2.6])
            uu = unique(d_prime_session_num{iRegion});
            mouseys = [d_prime_animal_num{iRegion, ...
                uu(~isnan(unique(d_prime_session_num{iRegion}(:, 1))), 1) ...
                }];
            title(num2str(mouseys))
            %for iMouse = 1:size(unique(mouseys),2)
            %    line([], [])
            %end

        end

        save('d_prime_session_fraction_gonogo.mat', 'd_prime_session_fraction')
        save('d_prime_session_fraction_gonogo_median.mat', 'd_prime_session_fraction_median')

        %% go no go
        passive = 0;
        goNogo = 0;
        cl_plot_dprime;

        % plot d prime
        figure(100);
        for iRegion = 1:size(plot_regions, 2)
            for iPair = 1:3
                subplot(3, size(plot_regions, 2), iPair+(iRegion - 1)*(size(plot_regions, 2)));
                hold on;
                % if 0, remove
                %kp = find(abs(d_prime{iRegion}(:,iPair)) < 1.4141 & abs(d_prime{iRegion}(:,iPair)) ~= 0);
                kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))));

                histogram(d_prime{iRegion}(kp, iPair), [0:0.05:2], 'Normalization', 'probability', 'FaceAlpha', 0.2) % 1.4142   -1.4142 = one of stims has 0 spikes
                %nanmedian(d_prime(kp,iPair,iRegion))
                %nanmean(d_prime(kp,iPair,iRegion))
                %title([num2str(thesePairs(iPair,1)) ' vs ' num2str(thesePairs(iPair,2)) ', mean= ' num2str(nanmean(abs(d_prime{iRegion}(kp,iPair))))])
                xlabel('d-prime')
                ylabel('fraction of neurons')
                xlim([0, 2])
                %makepretty_lite;
            end
        end

        figure('name', 'gogogo');
        hold on;
        for iRegion = 1:3
            num_recs = size(unique(d_prime_session_num{iRegion}(:, 1)), 1);
            theseColors = [lines(num_recs)];
            theseColors = repelem(theseColors, 3, 1);

            kp1 = find(abs(d_prime{iRegion}(:, 1)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 1))) & ~isnan(abs(d_prime{iRegion}(:, 1))));
            kp2 = find(abs(d_prime{iRegion}(:, 2)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 2))) & ~isnan(abs(d_prime{iRegion}(:, 2))));
            kp3 = find(abs(d_prime{iRegion}(:, 3)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, 3))) & ~isnan(abs(d_prime{iRegion}(:, 3))));

            subplot(3, 1, iRegion)
            violinplot([d_prime{iRegion}(kp1, 1); d_prime{iRegion}(kp2, 2); d_prime{iRegion}(kp3, 3)], ...
                [d_prime_session_num{iRegion}(kp1, 1); d_prime_session_num{iRegion}(kp2, 2) + 0.2; d_prime_session_num{iRegion}(kp3, 3) + 0.4], ...
                'ViolinColor', theseColors);
            ylim([0, 2.6])
            uu = unique(d_prime_session_num{iRegion});
            mouseys = [d_prime_animal_num{iRegion, ...
                uu(~isnan(unique(d_prime_session_num{iRegion}(:, 1))), 1) ...
                }];
            title(num2str(mouseys))

        end
        save('d_prime_session_fraction_gogogo.mat', 'd_prime_session_fraction')
        save('d_prime_session_fraction_gogogo_median.mat', 'd_prime_session_fraction_median')

        iPair = ppp;
        iRegion = rrr;
        kp = find(abs(d_prime{iRegion}(:, iPair)) ~= 0 & ~isinf(abs(d_prime{iRegion}(:, iPair))) & ~isnan(abs(d_prime{iRegion}(:, iPair))));
        uu = unique(d_prime_session_num{iRegion});
        uu = uu(~isnan(unique(d_prime_session_num{iRegion}(:, 1))));
        mouseys = [d_prime_animal_num{iRegion, ...
            uu(~isnan(unique(d_prime_session_num{iRegion}(:, 1))), 1) ...
            }];
        t2 = table;
        t2.dprime = d_prime{iRegion}(kp, iPair);
        t2.sessionType = ones(size(d_prime{iRegion}(kp, iPair), 1), 1) * 2;
        t2.session = d_prime_session_num{iRegion}(kp, iPair, 1);
        mousee = d_prime_session_num{iRegion}(kp, iPair);
        for iS = 1:size(uu, 1)
            mousee(mousee == uu(iS)) = mouseys(iS);
        end
        t2.mouse = mousee;
        t2.pss = pss{iRegion}(kp, iPair);
        t2.tempDur = tempDur{iRegion}(kp, iPair);
        t2.propISI = propISI{iRegion}(kp, iPair);
        t2.fr = fr{iRegion}(kp, iPair);

        lmeTable = [lmeTable; t2];

        %% linear mixed -effects model
        %lme = fitlme(tbl,'CityMPG~Horsepower+(1|EngineType)+(Horsepower-1|EngineType)')
        lme = fitglme(lmeTable, 'dprime ~ sessionType  + (1|mouse) + (1|mouse:session)') %what if I add 1|neuron?
        %         lme = fitglme(lmeTable, 'dprime ~ sessionType  + (1|neuron) cd(1|mouse) + (1|mouse:session)')
        %lme = fitglme(lmeTable, 'dprime ~ sessionType  +  (1|mouse)')
        %lme = fitglme(lmeTable, 'dprime ~ sessionType  +  (1|session)')

        %% violin plot
        figure(300)
        subplot(3, 3, ppp+3*(rrr - 1))
        violinplot(lmeTable.dprime, lmeTable.sessionType);
        ylim([0, 4])
        ylabel('dprime(nogo - go2)')
        xticks([1, 2])
        xticklabels({'go/no go', 'go go go'})
        xlabel('task')

        figure('name', ['region', num2str(rrr)]);
        subplot(1, 4, 1)
        scatter(lmeTable.tempDur, lmeTable.pss, 4, 'filled')
        xlabel('peak to trough duration (us)')
        ylabel('post spike suppression (ms)')

        subplot(1, 4, 4)
        stem3(lmeTable.tempDur, lmeTable.pss, lmeTable.fr)
        xlabel('peak to trough duration (us)')
        ylabel('post spike suppression (ms)')
        zlabel('firing rate')

        subplot(1, 4, 2)
        scatter(lmeTable.tempDur, lmeTable.propISI, 4, 'filled')
        xlabel('peak to trough duration (us)')
        ylabel('prop. ISI > 2s')

        subplot(1, 4, 3)
        scatter(lmeTable.propISI, lmeTable.pss, 4, 'filled')
        xlabel('prop. ISI > 2s')
        ylabel('post spike suppression (ms)')

        prettify_plot;

        figure(301)
        subplot(3, 3, ppp+3*(rrr - 1))
        if rrr == 1
            % msns
            violinplot(lmeTable.dprime(lmeTable.tempDur >= 500 & lmeTable.pss <= 20), ...
                lmeTable.sessionType((lmeTable.tempDur >= 500 & lmeTable.pss <= 20)));

            %fsi

            %tans

        else
            violinplot(lmeTable.dprime, lmeTable.sessionType);

        end
        ylim([0, 4])
        ylabel('dprime(nogo - go2)')
        xticks([1, 2])
        xticklabels({'go/no go', 'go go go'})
        xlabel('task')


        if rrr == 1
            figure(302)
            subplot(3, 3, (ppp - 1)*3+1)
            % msns
            violinplot(lmeTable.dprime(lmeTable.tempDur >= 500 & lmeTable.pss <= 20), ...
                lmeTable.sessionType((lmeTable.tempDur >= 500 & lmeTable.pss <= 20)));
            title('MSNs')
            ylim([0, 4.1])
            %fsi
            subplot(3, 3, (ppp - 1)*3+2)
            % msns
            violinplot(lmeTable.dprime(lmeTable.tempDur < 500), ...
                lmeTable.sessionType((lmeTable.tempDur < 500)));
            title('FSIs')
            ylim([0, 4.1])
            %tans
            subplot(3, 3, (ppp - 1)*3+3)
            % msns
            violinplot(lmeTable.dprime(lmeTable.tempDur >= 500 & lmeTable.pss > 20), ...
                lmeTable.sessionType((lmeTable.tempDur >= 500 & lmeTable.pss > 20)));
            title('TANs')
            ylim([0, 4.1])
        end

        %%

    end
end
