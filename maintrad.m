% This m file replicates the traditional momentum strategy 
% discussed in the emphirical study on skewness and momentum 
% by Jacobs, Regel and Weber (2015)
% ==================================================================
% Name:         Victoria Xie Qingtong
% Organization: Imperial College London
% ==================================================================

%% load data
cd('/Users/victoriaxie/Documents/Study/MSc Finance/Autumn Term/Econometrics/CW/code');
clear;
clc;
load('assignment2016.mat');

%% calculate monthly cumulative return
num_stocks = 16150;
num_tseries = 6306;
[year, month] = datevec(data(1:num_tseries, 2));
minyear = min(year);
maxyear = max(year);
mthret = NaN((maxyear-minyear+1)*12, num_stocks);
for k = 1 : num_stocks
    dayret = data((k-1)*num_tseries+1 : k*num_tseries, 4);
    for i = minyear : maxyear
        for j = 1 : 12
            mthret((i-minyear)*12+j, k) = sum(dayret(year==i & month==j));
        end
    end
end

%% calculate past 12-month cumulative return for 1991-Feb to 2014-Dec
cumret = NaN((maxyear-minyear+1)*12-13, num_stocks);
for t = 1 : size(cumret, 1)
    cumret(t, :) = sum(mthret(t:t+11, :));
end

%% form the traditional momentum portfolio
weight = zeros((maxyear-minyear+1)*12-13, num_stocks);
for t = 1 : size(cumret, 1)
    location = find(~isnan(cumret(t, :)));
    [~, idx] = sort(cumret(t, :), 'descend');
    sort_idx = intersect(idx, location, 'stable');
    decile = round(length(sort_idx)/10);
    long_idx = sort_idx(1:decile);
    short_idx = sort_idx(length(sort_idx)-decile+1:length(sort_idx));
    weight(t, long_idx) = 1/(2*decile);
    weight(t, short_idx) = -1/(2*decile);
end

%% calculate the return series of the traditional momentum strategy
furret = mthret(14:end, :);
furret(isnan(furret)) = 0;
tradret = diag(furret * weight');

%% report summary statistics of the traditional momentum strategy
avgret = mean(tradret);
[h,p,ci,stats] = ttest(tradret);
sharpe = sqrt(252)*avgret/stats.sd;