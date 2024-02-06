%% [FIG 1B]: example performance

animal = 'JF117';
protocol = '';
experiments = AP_find_experimentsJF(animal,protocol);
plot_days = 27;

figure;
h = tiledlayout(length(plot_days),1);
for curr_day = plot_days
 
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(2);

    load_parts.imaging = false;
    JF_load_experiment;
    
    t = Timeline.rawDAQTimestamps;
    
    % Convert wheel velocity from clicks/s to mm/s
    % (mm in clicks from +hw.DaqRotaryEncoder, Lilrig encoder = 100)
    wheel_click2mm = 0.4869;
    wheel_velocity_mm = wheel_velocity*wheel_click2mm;

    % (minimum ITI: new trial + trial quiescence)
    min_iti_t = signals_events.newTrialTimes;% + ...
        %signals_events.trialQuiescenceValues;
    
    plot_t = [1980,2040];%[1800,1860];
    plot_t_idx = t > plot_t(1) & t < plot_t(2);
    plot_stim_idx = find(stimOn_times > plot_t(1) & stimOn_times < plot_t(2))';
    plot_min_iti_t_idx = find(min_iti_t > plot_t(1) & min_iti_t < plot_t(2));
    plot_reward_idx = find(reward_t_timeline > plot_t(1) & reward_t_timeline < plot_t(2));
    
    nexttile; hold on;
    
    % Plot stim and rewards
    yyaxis left; hold on

    for i = plot_stim_idx
        line(repmat(stimOn_times(i),2,1), ...
            [0,1],'color',[0.9,0.2,0.2],'linewidth',2);
    end

    for i = plot_reward_idx
        line(repmat(reward_t_timeline(i),2,1), ...
            [0,1],'color','b','linewidth',2);
    end

    % Plot wheel velocity
    yyaxis right;
    plot(t(plot_t_idx),abs(wheel_velocity_mm(plot_t_idx)),'k');
    
    axis off;
    title(sprintf('Day %d',curr_day));
    
end

h_ax = flipud(allchild(h));
linkaxes(h_ax,'y');

% Draw scalebars
t_scale = 5;
vel_scale = 100;
AP_scalebar(t_scale,vel_scale);