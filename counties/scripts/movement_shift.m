% calculate homophily for during the two phases each 
clear all
close all
clc

groupNames = {'asian' 'black' 'latinx' 'white' 'mixed'};
cityNames = {'CHI' 'COLO' 'DAL' 'DET' 'FWTX' 'HOU' 'IND' 'LA' 'LV' 'MIA' 'NY' 'PHI' 'PHX' 'SD' 'SEA'};

% two phases
phase1 = [1:60]; 
phase2 = [61:360];

for i = 1:length(cityNames)
    i
    % load city specific data
    cityName = cityNames{i};
    [parentdir,~,~]=fileparts(pwd);
    RAW   = load([parentdir '/data_imported/' cityName '/data.mat']);
    
    for g = 1:5
        index = find(contains(RAW.groups,groupNames{g}));
        if isempty(index)
            group.ZIPs{g} = [];
        else
            group.ZIPs{g} = RAW.G{index};
        end
    end
    
    W1 = squeeze(sum(RAW.W(:,:,phase1),3));
    W2 = squeeze(sum(RAW.W(:,:,phase2),3));
    
    homo1{i} = Homophily(W1,group.ZIPs);
    homo2{i} = Homophily(W2,group.ZIPs);
    
    
    HOMO_SELF{1}(i,:) = homo1{i}.self; % ZIP level homophily during phase 1
    HOMO_SELF{2}(i,:) = homo2{i}.self; % ZIP level homophily during phase 2
    
    HOMO_GROUP{1}(i,:) = homo1{i}.group; % group level homophily during phase 1
    HOMO_GROUP{2}(i,:) = homo2{i}.group; % group level homophily during phase 22
end

% calculate homophily
function homo = Homophily(W,group_ZIPs)   
    for g = 1:5
        group_index = group_ZIPs{g};
        if isempty(group_index)
            homo.self(g) = 0;
            homo.group(g) = 0;
        else
            D = diag(W);
            % ZIP level homophily
            homo.self(g) = 100*sum(D(group_index))/sum(W(group_index,:),[1 2]);

            % group level homophily
            homo.group(g) = 100*sum(W(group_index,group_index),[1 2])/sum(W(group_index,:),[1 2]);
        end
    end
end