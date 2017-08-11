function [wvGrid, interpedData, f] = coadd(xData, yData, nPix, nPoints, nObs,nDays, w)
wvGrid = linspace(-w*1.5, w*1.5, nPoints);
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/b(4)).^2 * .5));
nParam = 4; 

interpedData = zeros(nDays, nPoints);
f = zeros(nDays, nParam);
for i = 1:nDays
    dayYAll = zeros(nObs, nPoints);
    for j = 1:nObs
        ind = j + nObs*(0:nPix -1);
        dayX = xData(i, ind);
        dayY = yData(i, ind);
        dayYAll(j, :) = interp1(dayX, dayY, wvGrid, 'spline', 1);
    end
    fitY = mean(dayYAll);
    interpedData(i, :) = fitY;
    b0 = [1, 1 - min(fitY), 0, 0.075];
    f(i, :) = nlinfit(wvGrid, fitY, mfit, b0);
end


end

