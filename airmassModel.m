r = 100;
x = -r + 1:r;
y = x;

c = 299792.458; %speed of light in km/sec
lambdaI = 6173.3;

mask = x.^2 + (y.^2)' < r^2;
%angle from center of sun/line of sight path to point on sun surface
rOut = sqrt(x.^2 + (y.^2)')/r; 
mu = .6;
limbDark = 1 - mu + mu*sqrt(1 - rOut.^2);
%weight lines based on limb darkening
limbDark(~mask) = 0;

%Accounting for differential solar rotation rate
%NOTE: got constants from wiki in deg/day-- converted to rad/sec
A = 0.244346 / 86400;
B = -0.041818089 / 86400;
C = -0.031189034 / 86400;
R = 695700; %radius of sun in km

latAng = abs(asin(repmat(y', 1, 2*r) / r));
longAng = acos(repmat(x, 2*r, 1)/ r) - pi/2;
omega = A + B*sin(latAng).^2 + C*sin(latAng).^4; %in rad/sec
%final RV mask, assuming a vertical axis of rotation
%R in km, so units of km/sec
vProj = R*cos(latAng).*omega.*sin(longAng);
vProj(~mask) = 0;


load avgJPL
nDays = length(edges);
lat = 28.7541; %latitude of la palma, in degrees
lng = 17.8892; %longitude, degrees
dec = avgJPL(:, 3);
ra = avgJPL(:, 2);
%z = deg2rad(90 + lat - avgJPL(:, 3)) - pi/2; %angle from zenith to mid of sun
a = asin(sind(dec) .* sind(lat) + cosd(dec).*cosd(lat).*cosd(avgHrAng));
z = pi/2 - a;
d = avgJPL(:, 1);

eta = asin(sind(avgHrAng).*cosd(lat)./cos(a));
axRot = eta - deg2rad(avgNPAng);
dz = tan(R./d);

%using linear approx
ext = avgJPL(:, 4); %.13? Nervous about this: not sure this is the ext coeff I want?
amass = avgJPL(:, 5);

%Generate synthetic data
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
cs = c * ones(nDays, 1); %normrnd(c, stdC, nDays, 1); HOLDING CONSTANT
ds = abs(normrnd(e2dsD, e2dsDstd, nDays, 1));
ws = normrnd(w, stdW, nDays, 1);


xData = zeros(nDays, nPix*nObs);
yData = zeros(nDays, nPoints);
yDataSkew = zeros(nDays, nPoints);
yFinal = zeros(nDays, nPix*nObs);
ySkewFinal = zeros(nDays, nPix*nObs);


c = 299792.458; %speed of light in km/sec
offsets = avgRV/ c *6173.3;

for i = 1:nDays
    idx = baseline + ceil(rand(1, nPix*nObs)* nObs*overCount);
    dayX = xOver(idx) + offsets(i);
    dayY = mfit([hs(i) ds(i) cs(i) ws(i)], xSimple);
    yDataSkew(i, :) = dayY + 10000*xSimple;
    xData(i, :) = dayX;
    yData(i, :) = dayY;
    final = mfit([hs(i) ds(i) cs(i) ws(i)], dayX);
    yFinal(i, :) = poissrnd(final);
    ySkewFinal(i, :) = poissrnd(final - 10000*dayX);
        
end

yFull = zeros(nDays, nPix*nObs);
ySkewFull = zeros(nDays, nPix*nObs);
test = zeros(nDays, 1);
for i = 1:nDays
    AM_grad = repmat(linspace(sec(z(i) - dz(i)), sec(z(i)+dz(i)), 2*r)' ./ sec(z(i)), 1, 2*r);
    AM_eff = amass(i) * ext(i) * AM_grad;
    weight = (1 - AM_eff) .* limbDark;
    weight = weight ./ sum(sum(weight));
    %Note: might need to flip the sign here?
    rotatedVProj = zeros(2*r);
    for j = 1:2*r
        for k = 1:2*r
            if (mask(j, k))
                try
                    xNew = round(cos(axRot(i))* (j - r) - sin(axRot(i))*(k - r)) + r;
                    yNew = round(sin(axRot(i))*(j - r) + cos(axRot(i)) * (k - r)) + r;
                    rotatedVProj(j, k) = vProj(xNew, yNew);
                    
                catch
                end
            end
        end
    end
    lambdaShift = lambdaI*(1 - rotatedVProj ./ c) - lambdaI;
    test(i) = sum(sum(lambdaShift .* weight));
    for j = 1:2*r
        for k = 1:2*r
            if weight(j, k) > 0
                ys = interp1(xSimple + lambdaShift(j, k), yData(i, :), xData(i, :), 'linear', e2dsH);
                ySks = interp1(xSimple + lambdaShift(j, k), yDataSkew(i, :), xData(i, :), 'linear', e2dsH);
                yFull(i, :) = yFull(i, :) + ys * weight(j, k);
                ySkewFull(i, :) = ySkewFull(i, :) + ySks * weight(j, k);
            end
        end
    end
end
    
yFull = poissrnd(yFull);
ySkewFull = poissrnd(ySkewFull);

rats = sort(yFull, 2);
ratsSkew = sort(ySkewFull, 2);
scalars = mean(rats(:, round(nObs * nPix * 2 / 3):end), 2);
scalarsSkew = mean(ratsSkew(:, round(nObs * nPix * 2 / 3):end), 2);
ys = yFull ./ scalars;
ysSkew = ySkewFull ./ scalarsSkew;


ratsF = sort(yFinal, 2);
ratsFSkew = sort(ySkewFinal, 2);
scalarsF = mean(ratsF(:, round(nObs * nPix * 2 / 3):end), 2);
scalarsFSkew = mean(ratsFSkew(:, round(nObs * nPix * 2 / 3):end), 2);
ysF = yFinal ./ scalarsF;
ysFSkew = ySkewFinal ./ scalarsFSkew;

fCen = zeros(nDays, 4);
fSkewCen = zeros(nDays, 4);
fFCen = zeros(nDays, 4);
fFSkewCen = zeros(nDays, 4);

ind = zeros(nDays, 2);
indSkew = zeros(nDays, 2);
for i = 1:nDays
    [~, lI] = min(abs(xData(i, :) + 1.75 * w));
    [~, rI] = min(abs(xData(i, :) - 1.75 * w));
    fitX = xData(i, lI:rI);
    fitY = ys(i, lI:rI);
    b0 = [1, 1 - min(fitY), 0, 0.075];
    fCen(i, :) = nlinfit(fitX, fitY, mfit, b0);
    ind(i, :) = [lI rI];

    fitYSkew = ysSkew(i, lI:rI);
    b0Skew = [1, 1 - min(fitYSkew), 0, 0.075];
    fSkewCen(i, :) = nlinfit(fitX, fitYSkew, mfit, b0Skew);
    indSkew(i, :) = [lI rI];
    
    
    fitYF = ysF(i, lI:rI);
    b0F = [1, 1 - min(fitYF), 0, 0.075];
    fFCen(i, :) = nlinfit(fitX, fitYF, mfit, b0F);
    
    
    fitYFSkew = ysFSkew(i, lI:rI);
    b0FSkew = [1, 1 - min(fitYFSkew), 0, 0.075];
    fFSkewCen(i, :) = nlinfit(fitX, fitYFSkew, mfit, b0FSkew);
   
end

