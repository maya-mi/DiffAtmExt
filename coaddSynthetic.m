nPoints = 100;
wvGrid = linspace(-w*1.75, w*1.75, nPoints);

interpedData = zeros(nDays, nPoints);
for i = 1:nDays
    dayYAll = zeros(nObs, nPoints);
    for j = 1:nObs
        ind = j + nObs*(0:nPix -1);
        dayX = xData(i, ind);
        dayY = yData(i, ind);
        dayYAll(j, :) = interp1(dayX, dayY, wvGrid, 'spline', 1);
    end
    interpedData(i, :) = mean(dayYAll);
end