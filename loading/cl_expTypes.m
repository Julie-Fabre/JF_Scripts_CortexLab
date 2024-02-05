function expTypes = cl_expTypes(load_type, exp_type)

    if contains(load_type, 'naive')
        expTypes = {'JF_GratingPassiveVarITI', 'JF_locations', 'JF_natural_imagesFit', 'JF_choiceworldStimuli_onlyTask', ...
        'JF_choiceworldStimuli', 'JF_natural_images'};
        if contains(load_type, 'gratings')

        elseif contains(load_type, 'locations')
        elseif contains(load_type, 'cwStims')
        end
    else
        expTypes = {'noGo_', 'JF_choiceworldStimuli'};
    end

end