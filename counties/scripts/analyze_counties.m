% analyze movement patterns

clear all
close all
clc

colors.group = { '#073b4c' '#ef476f' '#ffd166' '#06d6a0' '#118ab2' '#adb5bd'};
colors.anal  = {'#2d00f7' '#ff5400' '#f20089'};

groupNames = {'asian' 'black' 'latinx' 'white' 'mixed'};
cityNames = {'CHI' 'COLO' 'DAL' 'DET' 'FWTX' 'HOU' 'IND' 'LA' 'LV' 'MIA' 'NY' 'PHI' 'PHX' 'SD' 'SEA'};

% create figures for plotting analysis results
fig_edges = figure('Position',[0 0 700 1000]);
tlt = tiledlayout(5, 3);
tlt.TileSpacing = 'compact';
tlt.Padding = 'none';

fig_homophily = figure('Position',[0 0 700 1000]);
tlt = tiledlayout(5, 3);
tlt.TileSpacing = 'compact';
tlt.Padding = 'none';

fig_groups = figure('Position',[0 0 700 1000]);
tlt = tiledlayout(5, 3);
tlt.TileSpacing = 'compact';
tlt.Padding = 'none';

for i = 1:length(cityNames)
    
    % load city specific data
    cityName = cityNames{i};
    [parentdir,~,~]=fileparts(pwd);
    RAW   = load([parentdir '/data_imported/' cityName '/data.mat']);
    W_mean = mean(RAW.W,3);
    ZIP.pop    = RAW.N;
    ZIP.W      = RAW.W;
    ZIP.hetero = 100 - 100*diag(W_mean)./sum(W_mean,2);
    
    W = RAW.W;
    
    % trip distribution (self vs average non-self)
    w_frac = [];
    for t = 1:size(W,3)
       dW = triu(W(:,:,t),1)+tril(W(:,:,t),-1);
       row_mean = sum(dW,2)/(size(dW,2)-1);
       w_frac(:,t) = diag(W(:,:,t))./row_mean;
    end

    figure(fig_edges)
    nexttile(i)
    [~,edges] = histcounts(log10(w_frac(:)),250);
    histogram(w_frac(:),10.^edges,'EdgeColor','none','FaceColor',hex2rgb(colors.group{6}))
    set(gca,'XScale','log')
    axis square
    grid on
    xlabel('')
    ylabel('#')
    xlim([1e0 1e4])
    xticks([1e0 1e1 1e2 1e3 1e4])
    title(cityNames{i})
    
    % compute group and city level heterophily
    city.hetero(i) = 0;
    for g = 1:5
        index = find(contains(RAW.groups,groupNames{g}));
        if isempty(index)
            group.ZIPs{g} = [];
        else
            group.ZIPs{g} = RAW.G{index};
        end
        group.pop(g,:) = sum(RAW.N(group.ZIPs{g},:),1);
        group.names{g} = [upper(groupNames{g}(1)) groupNames{g}(2:end)];
        
        group.hetero(g) = 100 - 100*sum(W_mean(group.ZIPs{g},group.ZIPs{g}),[1 2])/sum(W_mean(group.ZIPs{g},:),[1 2]);
        city.hetero(i) = city.hetero(i) + sum(W_mean(group.ZIPs{g},group.ZIPs{g}),[1 2]);
    end
    city.hetero(i) = 100-100*city.hetero(i)/sum(W_mean,[1 2]);
       
    % scatter plot for ZIP level pop size and homophily
    figure(fig_homophily)
    nexttile(i)
    semilogx(0,0)
    hold on
    for g = 1:5
        index = group.ZIPs{g};
        s = scatter(mean(ZIP.pop(index,:),2),100-ZIP.hetero(index),50,'filled');
        s.MarkerFaceColor = colors.group{g};
    end
    axis square
    grid on
    xlim([1e2 1e6])
    xticks([1e2 1e3 1e4 1e5 1e6])
    ylim([0 100])
    xlabel('population')
    ylabel('homophily (%)')
    title(cityNames{i})
        
    % scatter plot connecting pop frac and homophily for each group
    figure(fig_groups)
    nexttile(i)
    hold on
    for g = 1:5
        s = scatter(100*sum(group.pop(g,:),2)/sum(group.pop,[1,2]),100-group.hetero(g),150,'filled');
        s.MarkerFaceColor = colors.group{g};
    end
    axis square
    grid on
    xlim([0 100])
    ylim([0 100])
    xlabel('population (%)')
    ylabel('homophily (%)')
    title(cityNames{i})
    
    % save heterophily and pop frac data for city level comparison
    overall.pop_frac(i,:) = 100*sum(group.pop,2)/sum(group.pop,[1,2]);
    overall.hetero(i,:)   = group.hetero;
    overall.city_het(i)   = city.hetero(i);
end

figure(fig_edges)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
print(gcf,[parentdir '/figures/cities_edges.eps'],'-depsc')

figure(fig_homophily)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
print(gcf,[parentdir '/figures/cities_homophily.eps'],'-depsc')

figure(fig_groups)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
print(gcf,[parentdir '/figures/cities_groups.eps'],'-depsc')

%% Radar plot for homophily, dot size is pop frac, homophily
figure('Position',[0 0 700 1000])
tlt = tiledlayout(5, 3);
tlt.TileSpacing = 'compact';
tlt.Padding = 'none';
for i = 1:length(cityNames)
    nexttile
    polarplot(0,0)
    hold on
    for g = 1:5
        if overall.pop_frac(i,g)>0
            polarplot(-pi/3 + g*pi/3,100-overall.hetero(i,g),'o','MarkerFaceColor',colors.group{g},'MarkerEdgeColor','none','MarkerSize',3*sqrt(overall.pop_frac(i,g)))   
        end
    end
    polarplot(-pi/3 + 2*pi,100-overall.city_het(i),'p','MarkerFaceColor',colors.group{6},'MarkerEdgeColor','none','MarkerSize',20)       
    thetaticks([0:60:360])
    thetaticklabels([group.names 'City'])
    rlim([0 100])
    rticks([0:25:100])
    rticklabels({'','','','',''})
    text(-0.5*pi, 50, '50%', 'horiz', 'center')
    text(-0.5*pi, 100, '100%', 'horiz', 'center')
    title(cityNames{i})
end
set(findall(gcf,'-property','FontSize'),'FontSize',10)
print(gcf,[parentdir '/figures/cityComparison_homo.eps'],'-depsc')







%% hex2rgb
function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary: 
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1; 
end
%% Convert from hex to rgb: 
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
end