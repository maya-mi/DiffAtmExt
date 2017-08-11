load edges
finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
notEmpty = nums > 5;
aV = zeros(sum(notEmpty), 1);

aM = fitsread('order_airmass.fits');
index = 0;
for i = 1:max(D)
    if notEmpty(i) %Enforcing min number of daily exposures
        index = index + 1;
        a = aM(i == D);
        aV(index) = max(a) - min(a);         
    end
end

aV = aV(hasSDO);

%use 24, 235

