load lineParams.mat
e2dsH = 12.5e+4;
e2dsHstd = stdH*e2dsH;

e2dsD = e2dsH * d;
e2dsDstd = stdD*e2dsH;

nPix = 18;
x = linspace(-w*3.5, w*3.5, nPix);
xOver = [];
nObs = 50;
overCount = 3;
clusterWidth = .005;
for j = 1:nPix
    xOver = [xOver, linspace(x(j) - clusterWidth, x(j) + clusterWidth, nObs * overCount)];
end

baseline = sort(repmat(nObs*overCount*((1:nPix) - 1), 1, nObs));

nPoints = 100;
xSimple = linspace(-w*4, w*4, nPoints);
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/b(4)).^2 * .5));

hs = normrnd(e2dsH, e2dsHstd, nDays, 1);
cs = normrnd(c, stdC, nDays, 1);
ds = abs(normrnd(e2dsD, e2dsDstd, nDays, 1));
ws = normrnd(w, stdW, nDays, 1);


xData = zeros(nDays, nPix*nObs);
yData = zeros(nDays, nPoints);
yDataSkew = zeros(nDays, nPoints);

avgRV = avgJPL(:, 2);
c = 299792.458; %speed of light in km/sec
offsets = avgRV/ c *6173.3;

for i = 1:nDays
    idx = baseline + ceil(rand(1, nPix*nObs)* nObs*overCount);
    dayX = xOver(idx) + offsets(i);
    dayY = mfit([hs(i) ds(i) cs(i) ws(i)], xSimple);
    yDataSkew(i, :) = dayY - 10000*xSimple;
    xData(i, :) = dayX;
    yData(i, :) = dayY;
%     yData(i, :) = poissrnd(dayY);
%     yDataSkew(i, :) = poissrnd(dayYSkew);
%     normConst = sort(dayY);
%     yData(i, :) = dayY ./ normConst(round(nObs * nPix * .95));
%     normConstSkew = sort(dayYSkew);
%     yDataSkew(i, :) = dayYSkew ./ normConst(round(nObs*nPix*.95));
        
end
    

 
