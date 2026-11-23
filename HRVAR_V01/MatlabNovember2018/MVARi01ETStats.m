%   Compute statistics across a series of MVARi save (mat) files
clc
clear all
close all
format
format compact
addpath('./mClasses','./mFiles')
diary off;
delete('allStatsET.txt')
diary('allStatsET.txt')

%%  setups
topDir = './HPC results/';
Dir2 = 'EventTime/';
nameRoot={'MVARiPartSIP01ET','MVARiLexExl01ET','MVARiDarkLit01ET'};

sdFile = 'symDates.mat';
isGroups = {{'nbbPart','nboPart'},{'nbbSIP','nboSIP'}};
%   Build symbol, date arrays.
% inFile = [inCsvDir 'part_nbb.csv'];
% x = readtable(inFile);
% allSymbols=unique(x.SYM_ROOT)
% allInDates = int64(unique(x.DATE));
% stopDate = allInDates(31);
% save(sdFile,'allSymbols','allInDates')
load(sdFile)
allInDates = repmat(allInDates(1:30),2,1);
allSymbols = reshape(repmat(allSymbols',30,1),[],1);

%%   check to see if anything is missing
for ir=1:length(nameRoot)
    %     nameRoot{ir}
    matDir = [topDir Dir2 'save/'];
    inFiles = dir(fullfile(matDir,[nameRoot{ir} '*.mat']));
    nFiles = length(inFiles);
    fprintf('\n%s\n%s %s: %d files found in %s\n',repmat('=',1,100),datetime,nameRoot{ir},nFiles,matDir)
    inSymbols = cell(0);
    inDates = [];
    for i=1:nFiles
        fn = inFiles(i).name;
        k = strfind(fn,"_");
        inSymbols(i) = extractBetween(fn,k(1)+1,k(2)-1);
        kp = strfind(fn,".");
        inDates(i) = str2double(string(extractBetween(fn,k(2)+1,kp-1)));
    end
    for i=1:60
        k = strcmp(allSymbols(i),inSymbols) & allInDates(i)==inDates;
        if ~any(k); fprintf('%2d (%s, %d) not found in mat files\n',i,allSymbols{i},allInDates(i)); end
    end
end

%%  Compute statistics over all days
for ir=1:length(nameRoot)
    %     nameRoot{ir}
    matDir = [topDir Dir2 'save/'];
    inFiles = dir(fullfile(matDir,[nameRoot{ir} '*.mat']));
    nFiles = length(inFiles);
    fprintf('\n\n%s\n%s %s: %d files in %s\n',repmat('=',1,100),datetime,nameRoot{ir},nFiles,matDir)
    for i=1:nFiles
        fn = inFiles(i).name;
        k = strfind(fn,"_");
        symbols(i) = extractBetween(fn,k(1)+1,k(2)-1);
        kp = strfind(fn,".");
        dates(i) = extractBetween(fn,k(2)+1,kp-1);
    end
    allSymbols = sort(unique(symbols));
    allDates = sort(unique(dates));
    for i=1:length(allSymbols)
        k = strcmp(allSymbols{i},symbols);
        sDates{i} = dates(k);
    end
    % for i=1:length(allDates)
    %     fprintf('%4d %s\n', i,allDates{i});
    % end
    for iSymbol=1:length(allSymbols)
        fprintf('%s\nStarting %s\n',repmat('-',1,100),allSymbols{iSymbol})
        sd = sDates{iSymbol};
        nr = length(sd);
        for iDate=1:nr
            k =find(strcmp(allSymbols{iSymbol},symbols) & strcmp(dates,sd{iDate}));
            load(fullfile(matDir,inFiles(k).name),'mv')
            fn2(iDate) = string(fullfile(matDir,inFiles(k).name));
            mvs(iDate) = mv;
        end
        nr = length(mvs);
        fprintf('Results for %d dates read.\n',nr)
        
        disp('Results for October 3, 2016')
        r = mvs(1).rwd;
        if isempty(r); 
            fprintf('Empty rwd for %s.',fn2(1));
            continue;
        end
            
        r.rDisplay
        columnLabels = {'Min','Max'};
   
        s = size(r.isGroupBounds);
        
        isgs=zeros(s(1),s(2),0);
        for i=1:nr
            if isempty(mvs(i).rwd) | isempty(mvs(i).rwd.isGroupBounds)
                fprintf('Note: For %s %d, rwd or rwd.isGroupBounds is empty.\n',mvs(i).symbol,mvs(i).inDate)
            else
                isgs = cat(3,isgs,mvs(i).rwd.isGroupBounds);
            end
        end
        nObs = size(isgs,3);
        fprintf('Means and SEMs over first %d days in sample (%d valid observations)\n',nr,nObs);
        m = mean(isgs,3);
        disptable(m,columnLabels,r.isGroupLabels,'%8.3f');
        sem = std(isgs,1,3)/sqrt(nr);
        disptable(sem,columnLabels,r.isGroupLabels,'%8.3f');
        continue
    end
end
fprintf('\n\n\n')
sdFile = 'symDates.mat';
%   Build symbol, date arrays.
% inFile = [inCsvDir 'part_nbb.csv'];
% x = readtable(inFile);
% allSymbols=unique(x.SYM_ROOT)
% allInDates = int64(unique(x.DATE));
% stopDate = allInDates(31);
% save(sdFile,'allSymbols','allInDates')
load(sdFile)
for i=1:length(allInDates)
    fprintf('%2d %10d\n',i,allInDates(i));
end
diary off;
