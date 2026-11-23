%   price disovery listing exchange vs non-listing exchanges in EVENT time
%   structured to allow indexing (over time and symbols) using environment variables
%   on HPC use 'submit -a 1:10 LexExl.s'
clc; clf
clear all; close all
format; format compact
addpath('./mClasses','./mFiles')
diary off;

inDate=20161003;    %   default date
symbol='IBM';      %   default symbol
saveDir = './EventTime/save/';
outDir = './EventTime/listing/';

%%  System-specific setup (are we running locally or in a batch environment)
switch computer
    case 'PCWIN64'
        hpc = false;
        pIndex = [];
    otherwise   %   running on NYU's high performance cluster
        hpc = true;
        pIndex = getenv('SLURM_ARRAY_TASK_ID');   %   Used in array jobs. Runs from 1 to 60 (1-30 are IBM; 31-60 are NVDA)
        if ~isempty(pIndex)
            k = str2double(pIndex);
            if k>0 && k<=60
                kSymbol = fix((k-1)/30) + 1;
                kDate = mod(k-1,30) + 1;
            end
        end
end

%%  Misc setup
inCsvDir = '../wrdsSasDatasets/';
progName = 'MVARiLexExl01ET';
isGroups = {{'nbbLex','nboLex'},{'nbbExl','nboExl'}};
%   Find all symbols and inDates
sdFile = 'symDates.mat';
%   Build symbol, date arrays.
% inFile = [inCsvDir 'part_nbb.csv'];
% x = readtable(inFile);
% allSymbols=unique(x.SYM_ROOT)
% allInDates = int64(unique(x.DATE));
% stopDate = allInDates(31);
% save(sdFile,'allSymbols','allInDates')
load(sdFile)
if exist('kSymbol','var'); symbol = allSymbols{kSymbol}; end
if exist('kDate','var'); inDate = allInDates(kDate); end

sdTimer = tic;
inDateVec = datetime(fix(inDate/10000), fix(mod(inDate/100,100)), mod(inDate,100));
fileRoot = sprintf('%s_%s_%s',progName,symbol,int2str(inDate));
diaryFile = [outDir fileRoot '.txt'];
% delete(fullfile(diaryFile));
saveFile = [saveDir fileRoot '.mat'];
diary(diaryFile);
fprintf('%s starting at %s\n',fileRoot,datetime);

mv = MVARi;

mv.symbol = symbol;
mv.inDate = inDate;

%%  read data

fprintf('Reading data for event-time analysis (%s)\n',datetime);
iTimer=tic;
% inFile = [inCsvDir 'part_nbb.csv'];
% opts = detectImportOptions(inFile)
t1 = getPriceET(readtable([inCsvDir 'lex_nbb.csv']),'nbb','timeMid',symbol,inDate);
t2 = getPriceET(readtable([inCsvDir 'lex_nbo.csv']),'nbo','timeMid',symbol,inDate);
t3 = getPriceET(readtable([inCsvDir 'exl_nbb.csv']),'nbb','timeMid',symbol,inDate);
t4 = getPriceET(readtable([inCsvDir 'exl_nbo.csv']),'nbo','timeMid',symbol,inDate);
cnan = nan(height(t1),1);
x = [t1.time t1.price cnan cnan cnan];
cnan = nan(height(t2),1);
x = vertcat(x, [t2.time cnan t2.price cnan cnan]);
cnan = nan(height(t3),1);
x = vertcat(x, [t3.time cnan cnan t3.price cnan]);
cnan = nan(height(t4),1);
x = vertcat(x, [t4.time cnan cnan cnan t4.price]);
seriesNames = {'nbbLex','nboLex','nbbExl', 'nboExl'};
[b,k] = sort(x(:,1));
xs = x(k,:);
n = size(xs,1);
for i=2:n   %   propagate missing prices
    for j=2:5
        if isnan(xs(i,j)); xs(i,j) = xs(i-1,j); end
    end
end
k = ~any(isnan(xs),2);
xs = xs(k,:);
k = xs(:,1)>=seconds(duration(9,45,0)) & xs(:,1)<=seconds(duration(16,0,0));
xs = xs(k,:);
disptable(xs(1:5,:),[],[],'%15.6f')
disptable(xs(length(xs)-5:end,:),[],[],'%15.6f')

%%  VECM in event time
mv.vecmNobs = length(x);
P1 = 10;
n = size(xs,2)-1;
%   set up model
mdl = vecm(n,3,P1);
B = [ones(n-1,1) -eye(n-1)];
mdl.Cointegration = B';
mdl.Constant = zeros(n,1);
mdl.Trend = zeros(n,1);
mdl.SeriesNames = seriesNames;
xs = xs - repmat(mean(xs),size(xs,1),1);    % demean the data.
[est,se] = estimate(mdl,xs(:,2:end),'Display','full','Model','H2','MaxIterations',1);
mv.vecmEst = est;
mv.vecmSE = se;

mv.rwd = randomWalkDecomp;
r = mv.rwd;
r.initm(est);
r.isBoundsGrouped(isGroups);
r.rDisplay;
r.isBounds(true);

save(saveFile,'mv')
et = etDuration(toc(sdTimer));
fprintf('\nEnding %s %d analysis at %s. Elapsed time is %s.\n\n',symbol, inDate, datetime, et)
diary off
