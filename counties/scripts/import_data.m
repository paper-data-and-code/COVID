% Import county level data (movement, population, demographic, case, area)

clear all
close all
clc


cityNames = {'CHI' 'COLO' 'DAL' 'DET' 'FWTX' 'HOU' 'IND' 'LA' 'LV' 'MIA' 'NY' 'PHI' 'PHX' 'SD' 'SEA'};

for j = 1:length(cityNames)
    cityName = cityNames{j}

    [parentdir,~,~]=fileparts(pwd);
    params_load.fileNameBase = [parentdir '/data_raw/' cityName '/'];
    [W N Z2G G groups] = LoadNodeData(params_load);
    C = LoadCityData(params_load);
    [C_rep C_est] = LoadCase_RepEst(params_load);
    [A ZIP_name] = LoadArea(params_load);

    save([parentdir '/data_imported/' cityName '/data.mat']);  
end



%% Load ZIP level data
function [W N Z2G G groups] = LoadNodeData(params_load)

    fileNameBase = params_load.fileNameBase;
    
    % Movement data
    fileName = [fileNameBase 'movement.xlsx'];
    opts = detectImportOptions(fileName);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        W(:,:,index) = table2array(M);
    end

    % Population size data
    fileName = [fileNameBase 'popsize.xlsx'];
    opts = detectImportOptions(fileName);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        N(:,index) = table2array(M(:,2));
    end

    % ZIP to group data
    fileName = [fileNameBase 'zips_to_groups.xlsx'];
    opts = detectImportOptions(fileName);
    opts.DataRange = 'A2';
    M = readtable(fileName,opts);
    groups = table2cell(unique(M(:,2)));
    for i = 1:length(groups)
        index = find(strcmp(string(groups(i,:)),table2cell(M(:,2))));
        Z2G(index,1) = i;
        G{i} = index;
    end
end


%% Load city level data
function C = LoadCityData(params_load)
    
    fileNameBase = params_load.fileNameBase;

    % City level data
    fileName = [fileNameBase 'count.xlsx'];
    opts = detectImportOptions(fileName);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        C(index) = table2array(M(1,4));
    end
end

%% Load city estimate data
function [C_rep C_est] = LoadCase_RepEst(params_load)
    fileNameBase = params_load.fileNameBase;
    fileName = [fileNameBase 'est_count.xlsx'];
    
    opts = detectImportOptions(fileName);
    M = readtable(fileName,opts);
    
    C_rep = table2array(M(:,2));
    C_est = table2array(M(:,3));
end

%% Load area data
function [A ZIP_name] = LoadArea(params_load)
    fileNameBase = params_load.fileNameBase;
    fileName = [fileNameBase 'area.xlsx'];

    opts = detectImportOptions(fileName);
    M = readtable(fileName,opts);
    
    ZIP_name = table2array(M(:,1));
    A = str2double(table2array(M(:,2)));
    
end