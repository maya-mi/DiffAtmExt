load('processed.mat')
load('propError.mat')
load('fitResults.mat')
pastFits = f;
[nObs, nLines, nPix] = size(normOrders);
finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/b(4)).^2 * .5));
notEmpty = nums > 5;
%Pre-designating fit vars
f = zeros(sum(hasSDO), nLines, 4);
errFit = zeros(sum(hasSDO), nLines, 4);
chisq = zeros(sum(hasSDO), nLines);
reduced = zeros(sum(hasSDO), nLines);
caughtCases = zeros(nLines, 1);
sdoI = 0;
index = 0;
for i = 1:max(D)
    if notEmpty(i) %Enforcing min number of daily exposures
        sdoI = sdoI + 1;
        if hasSDO(sdoI)
            index = index + 1;
            for j = 1:nLines
                fitX = squeeze(wavelengths(i == D, j, :)) - ironA(j) - mean(squeeze(pastFits(:, j, 3))); %hard-coding in appropriate range
                [~, lI] = max(find(mean(fitX) < - 1.7 * squeeze(mean(pastFits(:, j, 4)))));
                [~, rI] = max(find(mean(fitX) < 1.7 * squeeze(mean(pastFits(:, j, 4)))));
                try
                    fitX = fitX(:, lI:rI);
                    lB = mean(fitX(:, 1));
                    rB = mean(fitX(:, end));
                    fitX = reshape(fitX, 1, numel(fitX));
                    [fitX, I] = sort(fitX); %sorting to make bisectors easier
                    fitY = squeeze(normOrders(i == D, j, lI:rI)); 
                    fitY = reshape(fitY, 1, numel(fitY));
                    fitY = fitY(I);
                    err = squeeze(propError(i == D, j, lI:rI));
                    err = reshape(err, 1, numel(err));
                    errSq = err(I) .^2;
                    b0 = [1, 1 - min(fitY), 0, 0.075];
                    [f(index, j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, mfit, b0, 'Weights', (1 ./ errSq));
                    errFit(index, j, :) = sqrt(diag(cov));
                    model = mfit(f(index, j, :), fitX);
                    chisq(index, j) = sum((fitY - model) .^2 ./ errSq);
                    reduced(index, j) = chisq(index, j) ./ (length(fitX) - 4);
                catch
                    caughtCases(j) = caughtCases(j) + 1;
                    fitX = squeeze(wavelengths(i == D, j, 11:23)) - ironA(j); %hard-coding in appropriate range
                    fitX = reshape(fitX, 1, numel(fitX));
                    [fitX, I] = sort(fitX); %sorting to make bisectors easier
                    fitY = squeeze(normOrders(i == D, j, 11:23)); 
                    fitY = reshape(fitY, 1, numel(fitY));
                    fitY = fitY(I);
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
               
    end
end

save('centeredFitResults.mat', 'f', 'errFit', 'chisq', 'reduced', 'caughtCases', 'ironA')

