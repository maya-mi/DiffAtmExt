load results100.mat
nPoints = 50;
[~, ~, f] = coadd(xData, ys, nPix, nPoints, nObs, nDays, w);
[~, ~, fS] = coadd(xData, ysSkew, nPix, nPoints, nObs, nDays, w);
[~, ~, fF] = coadd(xData, ysF, nPix, nPoints, nObs, nDays, w);
[~, ~, fFS] = coadd(xData, ysFSkew, nPix, nPoints, nObs, nDays, w);