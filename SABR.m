Settle = '11-Mar-2019';
ExerciseDate = '18-SEP-2020';
MarketStrikes = [91 95 98 107 107.5 108 108.5 109 111 114 122 125.5 137.5 160]'/100;
MarketVolatilities = [21.868 20.866 20.151 18.244 18.157 18.065 17.966 17.862 17.495 16.971 15.713 15.231 14.080 13.444]'/100;
CurrentForwardValue = MarketStrikes(4);
ATMVolatility = MarketVolatilities(4);
Beta1 = 0.5;
% Calibrate Alpha, Rho, and Nu
objFun = @(X) MarketVolatilities - ...
    blackvolbysabr(X(1), Beta1, X(2), X(3), Settle, ...
    ExerciseDate, CurrentForwardValue, MarketStrikes);

X = lsqnonlin(objFun, [0.5 0 0.5], [0 -1 0], [Inf 1 Inf]);

Alpha1 = X(1);
Rho1 = X(2);
Nu1 = X(3);
%%METHOD 2 Calibrate Rho and Nu by Implying Alpha from At-The-Money Volatility
Beta2 = 0.5;
% Year fraction from Settle to option maturity
T = yearfrac(Settle, ExerciseDate, 1);

% This function solves the SABR at-the-money volatility equation as a
% polynomial of Alpha
alpharoots = @(Rho,Nu) roots([...
    (1 - Beta2)^2*T/24/CurrentForwardValue^(2 - 2*Beta2) ...
    Rho*Beta2*Nu*T/4/CurrentForwardValue^(1 - Beta2) ...
    (1 + (2 - 3*Rho^2)*Nu^2*T/24) ...
    -ATMVolatility*CurrentForwardValue^(1 - Beta2)]);

% This function converts at-the-money volatility into Alpha by picking the
% smallest positive real root 
atmVol2SabrAlpha = @(Rho,Nu) min(real(arrayfun(@(x) ...
    x*(x>0) + realmax*(x<0 || abs(imag(x))>1e-6), alpharoots(Rho,Nu))));
% Calibrate Rho and Nu (while converting at-the-money volatility into Alpha
% using atmVol2SabrAlpha)
objFun = @(X) MarketVolatilities - ...
    blackvolbysabr(atmVol2SabrAlpha(X(1), X(2)), ...
    Beta2, X(1), X(2), Settle, ExerciseDate, CurrentForwardValue, ...
    MarketStrikes);

X = lsqnonlin(objFun, [0 0.5], [-1 0], [1 Inf]);

Rho2 = X(1);
Nu2 = X(2);
% Obtain final Alpha from at-the-money volatility using calibrated parameters
Alpha2 = atmVol2SabrAlpha(Rho2, Nu2);

% Display calibrated parameters
C = {Alpha1 Beta1 Rho1 Nu1;Alpha2 Beta2 Rho2 Nu2};
CalibratedPrameters = cell2table(C,...
    'VariableNames',{'Alpha' 'Beta' 'Rho' 'Nu'},...
    'RowNames',{'Method 1';'Method 2'})


PlottingStrikes = (80:0.1:170)'/100;

% Compute volatilities for model calibrated by Method 1
ComputedVols1 = blackvolbysabr(Alpha1, Beta1, Rho1, Nu1, Settle, ...
    ExerciseDate, CurrentForwardValue, PlottingStrikes);

% Compute volatilities for model calibrated by Method 2
ComputedVols2 = blackvolbysabr(Alpha2, Beta2, Rho2, Nu2, Settle, ...
    ExerciseDate, CurrentForwardValue, PlottingStrikes);

figure;
plot(MarketStrikes,MarketVolatilities,'xk',...
    PlottingStrikes,ComputedVols1,'b', ...
    PlottingStrikes,ComputedVols2,'r', ...
    CurrentForwardValue,ATMVolatility,'ok',...
    'MarkerSize',10);
xlim([0.5 2]);
ylim([0.1 0.25]);
xlabel('Strike', 'FontWeight', 'bold');
ylabel('Implied Black Volatility', 'FontWeight', 'bold');
legend('Market Volatilities', 'SABR Model (Method 1)',...
    'SABR Model (Method 2)', 'At-the-money volatility');