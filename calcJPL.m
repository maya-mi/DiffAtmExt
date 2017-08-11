% Loading in jpl info

load jpl.mat
load cloudCuts

goodObs = datefind(cloudOK, jplTimes, minutes(1));
jplTimes = jplTimes(goodObs);
jplData = jplData(goodObs, :);
jplFields = {'delta', 'raObs', 'decObs', 'ext', 'amass'};
nJPL = 5;

[D, full_edges] = discretize(jplTimes, days);
nDays = length(full_edges);
avgJPL = zeros(nDays, nJPL);
avgT = zeros(nDays, 1);
for i = 1:nDays
    avgJPL(i, :) = mean(jplData(D == i, :), 1);
    avgT(i, :) = datenum(mean(jplTimes(D == i)));
end
 
load jplSupp.mat
goodObs = datefind(cloudOK, hrTimes, hours(1));
hrTimes = hrTimes(goodObs);
hrAng = hrAng(goodObs);
npAng = npAng(goodObs);

avgHrAng = zeros(nDays, 1);
avgNPAng = zeros(nDays, 1);

[D1, full_edges1] = discretize(hrTimes, days);

for i = 1:length(full_edges1)
    avgHrAng(i) = mean(hrAng(i == D1));
    avgNPAng(i) = mean(npAng(i == D1));
end

load edges.mat
avgHrAng = avgHrAng(ismember(full_edges1, edges));
avgNPAng = avgNPAng(ismember(full_edges1, edges));
avgJPL = avgJPL(ismember(full_edges, edges), :);
avgT = datetime(avgT(ismember(full_edges, edges), :), 'ConvertFrom', 'datenum');
T = (juliandate(avgT) - 2451545)/36525;
gmst = 24110.54841 + 8640184.812866*T + 0.093104*T.^2 - .0000062*T.^3;

load avgRV.mat
save('avgJPL.mat', 'avgJPL', 'avgRV', 'edges', 'avgT', 'T', 'gmst', 'avgHrAng', 'avgNPAng')