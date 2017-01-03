% This m file replicates the skewness-enhanced momentum strategy 
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
breakpoints = csvread('mebreakpoints.csv', 0, 1)*10^6;

%% filter out stocks below 10% breakpoint
num_stocks = 16150;
num_tseries = 6306;
[year, month] = datevec(data(1:num_tseries, 2));
minyear = min(year);
maxyear = max(year);
daymktcap = NaN(num_tseries, num_stocks);
dayret = NaN(num_tseries, num_stocks);
for k = 1 : num_stocks
    idx_start = (k-1)*num_tseries+1;
    idx_end = k*num_tseries;
    daymktcap(:, k) = data(idx_start : idx_end, 3);
    dayret(:, k) = data(idx_start : idx_end, 4);
end

for i = 1 : num_tseries
    idx = (year(i)-minyear)*12+month(i);
    filter = daymktcap(i, :)>breakpoints(idx, 2);
    %filter = ones(1, num_stocks); % fake filter
    dayret(i, :) = dayret(i, :)./filter;
end

%% caculate monthly cumulative return and maximum return
mthret = NaN((maxyear-minyear+1)*12, num_stocks);
maxret = NaN((maxyear-minyear+1)*12, num_stocks);
for i = minyear : maxyear
    for j = 1 : 12
        mthret((i-minyear)*12+j, :) = sum(dayret(year==i & month==j, :));
        maxret((i-minyear)*12+j, :) = max(dayret(year==i & month==j, :));
    end
end

%% calculate proxy for sknewness measure and past 12mth cumulative return
num_tseries = (maxyear-minyear+1)*12-13;
proxy = NaN(num_tseries, num_stocks);
cumret = NaN(num_tseries, num_stocks);
for t = 1 : num_tseries
    proxy(t, :) = maxret(t+12, :);
    cumret(t, :) = sum(mthret(t:t+11, :));
end

%% form skewness-enhanced momemtum portfolio
weight = zeros((maxyear-minyear+1)*12-13, num_stocks);
for t = 1 : size(cumret, 1)
    location = find(~isnan(cumret(t, :)) & ~isinf(cumret(t, :)));
    [~, idx] = sort(cumret(t, :), 'descend');
    momt_idx = intersect(idx, location, 'stable');
    location = find(~isnan(proxy(t, :)) & ~isinf(proxy(t, :)));
    [~, idx] = sort(proxy(t, :), 'ascend');
    skew_idx = intersect(idx, location, 'stable');
    momt_decile = round(length(momt_idx)/10);
    skew_decile = round(length(skew_idx)/10);
    momt_firstgrp = momt_idx(1 : momt_decile);
    skew_firstgrp = skew_idx(1 : skew_decile);
    momt_lastgrp = momt_idx(length(momt_idx)-momt_decile+1 : length(momt_idx));
    skew_lastgrp = skew_idx(length(skew_idx)-skew_decile+1 : length(skew_idx));
    long_idx = intersect(momt_firstgrp, skew_firstgrp, 'stable');
    short_idx = intersect(momt_lastgrp, skew_lastgrp, 'stable');
    weight(t, long_idx) = 1/length(long_idx);
    weight(t, short_idx) = -1/length(short_idx);
end

%% calculate the return series of the skewness-enhanced momentum strategy
furret = mthret(14:end, :);
furret(isnan(furret)) = 0;
furret(isinf(furret)) = 0;
tradret = diag(furret * weight');

%% report summary statistics of the skewness-enhanced momentum strategy
avgret = mean(tradret);
[h,p,ci,stats] = ttest(tradret);
sharpe = sqrt(252)*avgret/stats.sd;