function expTypes = cl_expTypes(load_type)

    if contains(load_type, 'passive')
    expTypes = {'JF_GratingPassiveVarITI', 'JF_locations', 'JF_natural_imagesFit', 'JF_choiceworldStimuli_onlyTask', ...
        'JF_choiceworldStimuli', 'JF_natural_images'};
    else
    expTypes = {'noGo_', 'JF_choiceworldStimuli'};
    end

end