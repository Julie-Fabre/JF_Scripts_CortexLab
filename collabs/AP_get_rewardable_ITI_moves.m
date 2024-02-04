% Grab all movements separately (in degrees from move start)
wheel_moves_deg = arrayfun(@(x) wheel_position_deg( ...
    Timeline.rawDAQTimestamps >= wheel_starts(x) & ...
    Timeline.rawDAQTimestamps <= wheel_stops(x)) - ...
    wheel_position_deg(find(Timeline.rawDAQTimestamps >= wheel_starts(x),1)), ...
    1:length(wheel_starts),'uni',false);

% Find movements that hit reward limit (and don't hit punish limit)
deg_reward = -90;
deg_punish = 90;

wheel_moves_deg_rewardlimit = find(cellfun(@(x) ...
    any(x <= deg_reward) && ...
    ~any(x >= deg_punish),wheel_moves_deg));

wheel_move_nostim_rewardable_idx = ...
    intersect(wheel_move_nostim_idx,wheel_moves_deg_rewardlimit);

% Get move no-stim align times
move_nostim_rewardable_align = wheel_starts(wheel_move_nostim_rewardable_idx);