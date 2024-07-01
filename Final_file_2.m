clc;
clear;

%% Generating returns from the closing prices

x =  readmatrix('ClosingPrices.csv');
returnT = randi([0 0], 1255, 55);
for i = 2 : 1256
  for j = 2 : 56
      a = x(i, j);
      b = x(i - 1, j);
      returnT(i - 1, j) = (a - b) / b;
  end
end
writematrix(returnT, "Returns.csv");

%% Minimizing MAD and finding out the weights

% This program has the returns of first 926 days of the closing prices
% table

T = 925;
n = 55;

Returns1 = readmatrix('Returns.csv');
Returns = Returns1(1 : T, 2 : 56);

% Trying to get the coefficients of the objective function first
% Number of 'd' variables = 'T'(Total number of days)
% Number of 'x' variables = 'n'(Total number of companies)

% ObjFunc has to be in a row vector form, hence transpose is taken in line
% 17

ObjFunc = cat(1, (1 / T) * (ones(T, 1)), zeros(55, 1))';

AvgOfComp = (1 / T) * sum(Returns, 1);

% Coming to constraints, there will be 750 + 750 + 1 constraints(and
% 1 equality constraint)

% First 750 constraints will have + symbol
% Then the next 750 constraints will have - symbol

PositiveConst = zeros(T, T + n);
NegativeConst = PositiveConst;

PositiveConst(1 : T, 1 : T) = (-1) * eye(T);
NegativeConst(1 : T, 1 : T) = (-1) * eye(T);

% column number 751 to 805
% row number 1 to 750

for i = 1 : T
    for j = T + 1 : T + n
       PositiveConst(i, j) = Returns(i, j - T) - AvgOfComp(j - T);
       NegativeConst(i, j) = (-1) * Returns(i, j - T) + AvgOfComp(j - T);
    end
end

MeanConst = cat(2, zeros(1, T), (-1) * AvgOfComp);

A = cat(1, PositiveConst, NegativeConst, MeanConst);
b = zeros(1, 2 * T + 1)';
b(height(b), 1) = -0.0010;

Aeq = cat(2, zeros(1, T), ones(1, n));
Beq = ones(1, 1);

lb = zeros(1, T + n);
ub = [];

[X, Z] = linprog(ObjFunc, A, b, Aeq, Beq, lb, ub);
Weights = X(T + 1 : T + n, 1);

%% Finding out the Mean, Standard Deviation, Sharpe Ratio, VaR, CVaR and Starr Ratio of the portfolio returns

LeftoverData = Returns1(T + 1 : 1255, 2 : 56);
PortfolioReturns = LeftoverData * Weights;

MeanPortfolioReturns = mean(PortfolioReturns);
StdDev = std(PortfolioReturns, 0, 1);
SharpeRatio = MeanPortfolioReturns / StdDev;

SortedPortfolioReturns = sort(PortfolioReturns);

[VAR, Index] = computeHistoricalVaR(PortfolioReturns, 0.95);
CVAR = mean(SortedPortfolioReturns(1 : Index));
StarrRatio = mean(PortfolioReturns) / (-1) * CVAR;


% NAIVE
EqualWts = (1 / 55) * ones(55, 1);
NaivePortfolioReturns = LeftoverData * EqualWts;

MeanNaive = mean(NaivePortfolioReturns);
StdDevNaive = std(NaivePortfolioReturns, 0, 1);
SharpeRatioNaive = MeanNaive / StdDevNaive;

NaiveSorted = sort(NaivePortfolioReturns);

[VARNaive, IndexNaive] = computeHistoricalVaR(NaivePortfolioReturns, 0.95);
CVARNaive = mean(NaiveSorted(1 : IndexNaive));
StarrRatioNaive = mean(NaivePortfolioReturns) / (-1) * CVARNaive;

Portfolios = ["Using LINPROG"; "Using NAIVE"];
Means = [MeanPortfolioReturns; MeanNaive];
StandardDev = [StdDev; StdDevNaive];
SharpeRat = [SharpeRatio; SharpeRatioNaive];
VaR = [VAR; VARNaive];
CVaR = [CVAR; CVARNaive];
StarrRat = [StarrRatio; StarrRatioNaive];

% The comparison of the values is stored in the "Comparison" table
Comparison = table(Portfolios, Means, StandardDev, SharpeRat, VaR, CVaR, StarrRat);

function [VaR, VaR_index] = computeHistoricalVaR(returns,confidence_level)
    sorted_returns = sort(returns);
    num_returns = numel(returns);
    VaR_index = ceil((1-confidence_level)*num_returns);
    VaR = sorted_returns(VaR_index);
end