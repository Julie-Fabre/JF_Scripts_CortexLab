
%% Choiceworld trial activity (striatum domain)
%clear all
disp('Choiceworld trial activity (striatum domain)')

trained_passive_animals = { 'AP028', 'AP029'};
trained_animals = {'AP024', 'AP025', 'AP026', 'AP027', 'AP028', 'AP029'};
naive_animals = {'AP032', 'AP033', 'AP034', 'AP035', 'AP036'}; %'AP034',
ctx_passive_animals =  {'AP043', 'AP060', 'AP061'};
trained = 1; %use trained or naive animals
passive = 0;
kal = 0;
ctxpassive=0;
if trained == 1 && passive == 1
    animals  = trained_passive_animals;
    protocol  = 'AP_choiceWorldStimPassive';
    protocol2 = 'AP_choiceWorldStimPassive';
elseif trained == 1
    animals = trained_animals;
    protocol  = 'vanillaChoiceworld';
    protocol2 = 'vanillaChoiceworld';
elseif kal == 1
    animals = naive_animals;
    protocol  = 'stimKalatsky';
    protocol2 = 'stimKalatsky';
elseif ctxpassive==1
    animals = ctx_passive_animals;
    protocol  = 'AP_lcrGratingPassive';
    protocol2 = 'AP_lcrGratingPassive';
    
else
    animals = naive_animals;
    protocol  = 'AP_choiceWorldStimPassive';
    protocol2 = 'AP_choiceWorldStimPassive';
end
% Initialize save variable
trial_data_all = struct;
corona = 1;
previous = 0; %use previously claulated metrics
if previous
    load('C:\Users\Julie\Dropbox\Analysis\andyData\ephysData')
    load('C:\Users\Julie\Dropbox\Analysis\andyData\qMetrics')
    load('C:\Users\Julie\Dropbox\Analysis\andyData\qMetricsAdd')
    load('C:\Users\Julie\Dropbox\Analysis\andyData\qMetricsAddI')
    load('C:\Users\Julie\Dropbox\Analysis\andyData\ephysParams')
    load('C:\Users\Julie\Dropbox\Analysis\andyData\cellTypes')
end

thisCount = 0;
redo = 0; %re-calculate quality metrics and ephys params for cell type classification

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    
    if ctxpassive ==1
        experiments = cl_find_experiments(animal, protocol, protocol2);
    else
        experiments = cl_find_experiments(animal, protocol, protocol2);
    end
    experiments = experiments([experiments.imaging] & [experiments.ephys]);

    disp(['Loading ', animal]);
    allExps{curr_animal, :} = experiments;
    
    for curr_day = 1:length(experiments)
        thisCount = thisCount + 1;
        thisDate = experiments(curr_day).thisDate;
        experiment = experiments(curr_day).experiment;
        % Load experiment
        n_aligned_depths = 3;
        str_align = 'kernel';
        load_parts.ephys = true;
        load_parts.cam = false;
        load_parts.imaging = true;
        verbose = false;
        [timeline_filename, timeline_exists] = AP_cortexlab_filename(animal, thisDate, experiment, 'timeline');
        [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,thisDate,experiment,'protocol');

        if timeline_exists
           AP_load_experiment;
           
            if ephys_exists
                
                    AP_ctx_str_grab_trial_dataJF_SUA;
                    % Store trial data into master structure
                    trial_data_fieldnames = fieldnames(trial_data);

                    for curr_trial_data_field = trial_data_fieldnames'

                        trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal, 1}{curr_day, 1} = ...
                            trial_data.(cell2mat(curr_trial_data_field));

                    end
                    % Store general info
                    trial_data_all.animals = animals;
                    trial_data_all.t = t;
                    if trained == 1 && passive == 0
                        trial_data_all.task_regressor_labels = task_regressor_labels;
                        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
                    end

                    AP_print_progress_fraction(curr_day, length(experiments));
                    % Clear for next loop
                
                clearvars -except allExps corona animals curr_animal animal protocol experiments curr_day ...
                    trial_data_all ephysData qMetrics qMetricsAddI qMetricsAdd cellTypes ephysParams thisCount redo ...
                    trained protocol passive protocol2 kal ctxpassive
            end
        end
    end
end

clearvars -except trial_data_all redo trained protocol passive protocol2 kal ctxpassive
disp('Finished loading all')
% Save
save_path = '/home/julie/Dropbox/MATLAB/onPaths/PhD_thesis/andyData';
if trained == 1 && passive == 1
    save_filename = [save_path, filesep, 'trial_activity_trainedPassive'];
elseif trained == 1
    save_filename = [save_path, filesep, 'trial_activity_trained'];
elseif kal == 1
    save_filename = [save_path, filesep, 'trial_activity_naiveKalatsky'];
elseif ctxpassive == 1
    save_filename = [save_path, filesep, 'ctx_passive'];
else
    
    save_filename = [save_path, filesep, 'trial_activity_naive'];
end

save(save_filename, '-v7.3');
disp(['Saved ', save_filename]);
