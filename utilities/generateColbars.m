figure()
subplot(1,3,1)
originalColormap = brewermap([], '*RdBu');
        middleIndex = ceil(size(originalColormap, 1)/2);
        whiteColor = [1, 1, 1];
        grayColor = [0.5, 0.5, 0.5];
        modifiedColormap = [grayColor; grayColor; grayColor; grayColor; grayColor; grayColor; grayColor; grayColor;...
            grayColor; grayColor; grayColor; grayColor;originalColormap(1:middleIndex-1, :); whiteColor; whiteColor; originalColormap(middleIndex+1:end, :)];

        % Apply the modified colormap
        colormap(modifiedColormap);
 

caxis([-0.2, 0.2]);
cb =colorbar; 
cb.Label.String = '|ΔFR|/FR'


subplot(1,3,2)
originalColormap = brewermap([], '*RdBu');
        middleIndex = ceil(size(originalColormap, 1)/2);
        whiteColor = [1, 1, 1];
        grayColor = [0.5, 0.5, 0.5];
        modifiedColormap = [grayColor; grayColor; grayColor; grayColor; grayColor; grayColor; grayColor; grayColor;...
            grayColor; grayColor; grayColor; grayColor;originalColormap(1:middleIndex-1, :); whiteColor; whiteColor; originalColormap(middleIndex+1:end, :)];

        % Apply the modified colormap
        colormap(modifiedColormap);
 

caxis([-0.15, 0.15]);
cb2 =colorbar; 
cb2.Label.String = '|ΔFR|/FR'




figure()
subplot(1,3,3)
originalColormap = brewermap([], 'Reds');
 
% Apply the modified colormap
colormap(originalColormap);
 

caxis([0, 1]);
cb =colorbar; 
cb.Label.String = ['Normalized' newline 'projection' newline 'strength'];
