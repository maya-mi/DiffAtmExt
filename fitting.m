load('processed.mat')
load('propError.mat')
[nObs, nLines, nPix] = size(normOrders);
finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/b(4)).^2 * .5));
notEmpty = nums > 5;
%Pre-designating fit vars
f = zeros(sum(notEmpty), nLines, 4);
errFit = zeros(sum(notEmpty), nLines, 4);
chisq = zeros(sum(notEmpty), nLines);
reduced = zeros(sum(notEmpty), nLines);


index = 0;
for i = 1:max(D)
    if notEmpty(i) %Enforcing min number of daily exposures
        index = index + 1;
        for j = 1:nLines
            fitX = squeeze(wavelengths(i == D, j, 11:23)) - ironA(j); %hard-coding in appropriate range
            fitX = reshape(fitX, 1, numel(fitX));
            [fitX, I] = sort(fitX); %sorting to make bisectors easier
            fitY = squeeze(normOrders(i == D, j, 11:23)); 
            fitY = reshape(fitY, 1, numel(fitY));
            fitY = fitY(I);
            plot(fitX, fitY)
            pause(2)
            err = squeeze(propError(i == D, j, 11:23));
            err = reshape(err, 1, numel(err));
            errSq = err(I) .^2;
            b0 = [1, 1 - min(fitY), 0, 0.075];
            [f(index, j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, mfit, b0, 'Weights', (1 ./ errSq));
            errFit(index, j, :) = sqrt(diag(cov));
            model = mfit(f(index, j, :), fitX);
            chisq(index, j) = sum((fitY - model) .^2 ./ errSq);
            reduced(index, j) = chisq(index, j) ./ (length(fitX) - 4);
        end
               
    end
end

f = f(hasSDO, :, :);
errFit = errFit(hasSDO, :, :);
chisq = chisq(hasSDO, :, :);
reduced = reduced(hasSDO, :, :);

save('fitResults.mat', 'f', 'errFit', 'chisq', 'reduced')