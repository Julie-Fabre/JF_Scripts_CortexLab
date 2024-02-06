function dataSet_combined = cl_combineDataSets(dataSet1, dataSet2, index1, index2)
    
try
    dataSet_combined.psth{1} = [dataSet1.psth{index1}; dataSet2.psth{index2}]; 
catch
    warning('psth size incompatible')
end
    dataSet_combined.unitType = [dataSet1.unitType, dataSet2.unitType]; 
    dataSet_combined.unit_coords = [dataSet1.unit_coords; dataSet2.unit_coords]; 
    dataSet_combined.unit_area = [dataSet1.unit_area; dataSet2.unit_area]; 
    dataSet_combined.pss = [dataSet1.pss, dataSet2.pss]; 
    dataSet_combined.propLongISI = [dataSet1.propLongISI, dataSet2.propLongISI]; 
    dataSet_combined.templateDuration = [dataSet1.templateDuration, dataSet2.templateDuration]; 
    dataSet_combined.fr = [dataSet1.fr, dataSet2.fr]; 
end