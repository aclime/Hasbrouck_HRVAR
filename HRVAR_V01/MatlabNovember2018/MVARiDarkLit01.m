%   price disovery listing exchange vs non-listing exchange
%   structured to allow indexing (over time and symbols) using environment variables
%   on HPC use 'submit -a 1:10 DarkLit.s'
clc; clear all; close all
format; format compact
addpath('./mClasses','./mFiles')
diary off;

%%  General setup
secondFactors = 10.^(0:5);
nSecondFactors = length(secondFactors);
mva(nSecondFactors) = MVARi;
inDate=20161003;    %   default date
symbol='IBM';      %   default symbol
saveDir = './DarkLit01/save/';
outDir = './DarkLit01/listing/';

%%  System-specific setup (are we running locally or in a batch environment)
switch computer
    case 'PCWIN64'
        hpc = false;
        nCPU = 0;
        pIndex = zeros(0,1);
    otherwise   %   running on NYU's high performance cluster
        hpc = true;
        nCPU = 8;
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
progName = 'MVARiDarkLit01';
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

%%  secondFactor loop (analyze at multiple time resolutions)
for iSec=1:nSecondFactors
    %%   secondFactor loop initializations
    secondFactor = secondFactors(iSec);
    fprintf('\n%s secondFactor=%d\n',fileRoot,secondFactor)
    mv = mva(iSec);
    fprintf('Setting polynomials\n')
    lagPDLmax = 10^6;
    mv.polys = polynom(0,1,'p0',0);
    for i=0:min([log10(secondFactor) log10(lagPDLmax)-1])
        newPoly = polynom(0,9*10^i,['p' int2str(i+1)],10^i);
        mv.polys = horzcat(mv.polys, newPoly);
    end
    for i=1:length(mv.polys); fprintf('%s\n', mv.polys(i).toString()); end
    mv.symbol = symbol;
    mv.inDate = inDate;
    mv.secondFactor = secondFactor;
    T=18*3600*secondFactor;
    spd.setgetMax(T);
    
    %%  read data
    fprintf('Reading data for secondFactor=%d (%s)\n',secondFactor,datetime);
    iTimer=tic;
    % inFile = [inCsvDir 'part_nbb.csv'];
    % opts = detectImportOptions(inFile);
    % x = readtable(inFile);
    mDisplay = true;
    t = getPrice(readtable([inCsvDir 'part_nbb.csv']),'nbb','timeMid',secondFactor,mDisplay,symbol,inDate);
    prices(1)=spl(t.time2, t.price,[],[],'nbbPart');
    t = getPrice(readtable([inCsvDir 'part_nbo.csv']),'nbo','timeMid',secondFactor,mDisplay,symbol,inDate);
    prices(2)=spl(t.time2, t.price,[],[],'nboPart');
    t = getPrice(readtable([inCsvDir 'ctLit.csv']),'PRICE','PART_TIME',secondFactor,mDisplay,symbol,inDate);
    prices(3) = spl(t.time2, t.price, [], [], 'trLit');
    t = getPrice(readtable([inCsvDir 'ctDark.csv']),'PRICE','PART_TIME',secondFactor,mDisplay,symbol,inDate);
    prices(4) = spl(t.time2, t.price, [], [], 'trDark');
    isGroups = {{'nbbPart','nboPart'},{'trLit'},{'trDark'}};    %   groupings used to compute information shares
    mv.prices = copy(prices);
    mv.priceNames = cellstr( arrayfun(@(x) x.name,mv.prices,'UniformOutput',false) );
    save(saveFile,'mva');
    et = etDuration(toc(iTimer));
    fprintf('Data read finished for %s %d. secondFactor=%d. Elapsed time is %s.\n',symbol,inDate,secondFactor,et);
    % mv.pricePlot;
    
    %%  Estimate
    iTimer = tic;
    load(saveFile);
    mv = mva(iSec);
    fprintf('starting estimates secondFactor=%d. %s\n(%s)\n',secondFactor,datetime,strjoin(mv.priceNames));
    mv.intercept = false;
    mv.ecm = true;
    spd.setgetMax(mv.prices(1).maxSize);
    mv.setup;
    mv.setNamesPDL();
    mv.zpzLayoutPDL();
    mv.setValid(mv.polys.maxLength,duration(9,45,0),duration(16,0,0),true);
    mv.eVecDemean(true);
    zpzDisplay = 1;
    mv.buildzpzSparsePDLpar(zpzDisplay,nCPU);
    mv.estimateSparse();
    mv.dispEstimates
    save(saveFile,'mva');
    et = etDuration(toc(iTimer));
    fprintf('Estimation finished for %s %d. secondFactor=%d. Elapsed time is %s.\n',symbol,inDate,secondFactor,et);
    
    %% Impulse response functions
    iTimer = tic;
    load(saveFile);
    mv = mva(iSec);
    fprintf('starting irf computation for secondFactor=%d. %s\n',mv.secondFactor,datetime);
    nAhead = min(500*secondFactor,500*10^4);
    irfDisplayLevel = 2;
    mv.irfPDLpacked(nAhead, irfDisplayLevel, nCPU);
    save(saveFile,'mva')
    et = toc(iTimer);
    fprintf('%s %d irf finished. elapsed time is %3.2f seconds.\n',symbol,inDate,et);
    et = etDuration(toc(iTimer));
    fprintf('IRF calculation finished for %s %d. secondFactor=%d. Elapsed time is %s.\n',symbol,inDate,secondFactor,et);
    
    %% Random-walk decompositions
    load(saveFile);
    mv = mva(iSec);
    fprintf('Random walk decomposition for secondFactor=%d. %s\n',mv.secondFactor,datetime);
    mv.rwd = randomWalkDecomp;
    r = mv.rwd;
    r.init(mv);
    r.isBoundsGrouped(isGroups);
    r.rDisplay;
    r.isBounds(true);
    save(saveFile,'mva');
end
et = etDuration(toc(sdTimer));
fprintf('\nEnding %s %d analysis at %s. Elapsed time is %s.\n\n',symbol, inDate, datetime, et)
diary off
return
