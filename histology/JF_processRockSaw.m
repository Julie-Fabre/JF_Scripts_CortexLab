%% rocksaw histology process 

%% load images 


%% automatic registering using elastix 
 reg=elastix(lenaTrans,lena,[],'elastix_default.yml','paramstruct',p);
 
 %% user xxx 