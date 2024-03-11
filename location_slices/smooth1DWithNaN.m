function smoothedSeries = smooth1DWithNaN(pixelSeries, kernel)
    % Check if the series is entirely NaN or has sufficient non-NaN points
    if all(isnan(pixelSeries)) || numel(find(~isnan(pixelSeries))) < 2
        % If all values are NaN or less than 2 non-NaN values exist,
        % return the original series without modification
        smoothedSeries = pixelSeries;
        return;
    end
    
    % Find NaN indices
    nanIdx = isnan(pixelSeries);

    % Interpolate NaN values linearly for smoothing purposes
    pixelSeriesInterpolated = pixelSeries;
    
    % Protect against cases with insufficient non-NaN data
    try
        pixelSeriesInterpolated(nanIdx) = interp1(find(~nanIdx), pixelSeries(~nanIdx), find(nanIdx), 'linear', 'extrap');
    catch
        % In case interpolation still fails, default to not altering the series
        pixelSeriesInterpolated(nanIdx) = NaN;
    end

    % Apply 1D smoothing
    smoothedSeriesInterpolated = conv(pixelSeriesInterpolated, kernel, 'same');

    % Reinstate NaN values to their original positions
    smoothedSeries = smoothedSeriesInterpolated;
    smoothedSeries(nanIdx) = NaN;
end
