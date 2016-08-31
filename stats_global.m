%% sourisfolle() - experimental economics of the feeling of control

%% Conditions :
  % 1 = Slow target & High control
  % 2 = Fast target & High control
  % 3 = Slow target & Low control
  % 4 = Fast target & Low control
  
  
%% Parametres :
    inter = 1;     %inter=1 : interversion des croyances telles que les performances aux conditions 1 et 3 soient superieures respectivement aux conditions 2 et 4 (correction des donnes en cas d'erreur du sujet) // inter = 0 : pas d'interversion
    compareToSubj = 0;  % Estimation de la preference pour le controle en fonction des croyances (performances subjectives). Mettre 0 pour estimer en fonction des performances objectives. 

        
    pref_control = 0;    % [Variable locale] Comptage des sujets qui expriment une preference pour le controle 
	pref_control_global = 0;    % [Variable globale] Comptage des sujets qui expriment une preference pour le controle 
	excluded_values0 = [4, 7, 13, 20];   % Pilote : Sujets 4 et 13 : absents // sujets 7 et 20 : outliers (switching points aberranst)
	excluded_values1 = [6, 10];   % Expe1 : Sujet6 : absent // sujet 10 : outlier (switching points aberrants)
	excluded_values2 = [7, 12, 15, 18, 19];   % Expe2 : Sujet 7, 12, 18 : absents // sujets 15, 19 : outliers (switching points aberrants)
    excluded_values3 = [1, 4, 6];   % Expe3 : Sujet 4 : absents // sujets 1 et 6: outliers (switching points aberrants)    

    
%% Extraction et traitement des donnees    
% Pilote
clear all
close all
clc
kk=0;
% risk_aversion = [];
basedir = 'G:\client_sourisfolle\stats_globales';
dirlist = {'180716','250716_13h','250716_14h30','280716'};
exclist = {[4, 13],[6],[7, 12, 18],[4]}; % absents
for d=1:4    
    for k=1:20
        % take a new participant
        if sum(k==exclist{d})==0
            filename = fullfile(basedir,['DataAt' dirlist{d}],['sourisfolle_' num2str(k) '.mat']);
            load(filename)
        
            kk=kk+1;
            datas{kk}=data;
            
            % aggregate replay data
            replay_suj = rmfield(data.replay.trials,'xy');
            replay_suj = catstruct(replay_suj);
            replay_suj.suj = ones(size(replay_suj.condition))*kk;
            if kk==1
                replay_all = replay_suj;
            else
                replay_all = catstruct(replay_all,replay_suj);
            end
            
            % aggregate main task data
            
            maintask_suj = catstruct(data.maintask.trials);
            maintask_suj.suj = ones(size(maintask_suj.condition))*kk;
            if kk==1
                maintask_all = maintask_suj;
            else
                maintask_all = catstruct(maintask_all,maintask_suj);
            end
            
            % aggregate subjective performance data
            
            subjective_suj = catstruct(data.successrate);
            subjective_suj.suj = ones(size(subjective_suj))*kk;
            subjective_suj.condition = [1 2 3 4];
            if kk==1
                subjective_all = subjective_suj;
            else
                subjective_all = catstruct(subjective_all,subjective_suj);
            end
            
            
        end % end of participant
        
    end  % end of group

    % now load questionnaire data for the group
%     tabAdress = fullfile(basedir,['Quest' dirlist{d} '.xlsx']);
%     [num, txt, tab] = xlsread(tabAdress);
%     for i=1:length(txt)
%         if strcmp(txt(i,14),'A')
%             risk_aversion = [risk_aversion, txt(i,14)];
%         end
%         if strcmp(txt(i,14),'B')
%             risk_aversion = [risk_aversion, txt(i,14)];
%         end
%     end 




end

% condition_subj = [];
% for j=0:3
%     for i=1:73
%         condition_subj = [condition_subj,1+j];
%     end
% end
% subjective_all.condition = condition_subj;

%% ANALYSIS
% identify outliers


% remove outliers


% aggregate variables of interest


bySub.caught_byCond = aggreg1fac([replay_all.caught],'mean',replay_all.suj,replay_all.condition);
bySub.rtReplay_byCond = aggreg1fac([replay_all.rt],'nanmedian',replay_all.suj,replay_all.condition);
bySub.speedmouse_byCond = aggreg1fac([replay_all.speedmouse],'mean',replay_all.suj,replay_all.condition);
bySub.speedtarget_byCond = aggreg1fac([replay_all.speedtarget],'mean',replay_all.suj,replay_all.condition);
bySub.subjective_byCond = aggreg1fac([subjective_all.r],'mean',subjective_all.suj,subjective_all.condition);



% regressions



 
% display(['Sujet numero ' num2str(k) ' :']);

%Performances objectives :
    % Calcul des performances objectives
    perf = objective_performance(data);
    % Condition 1 
    perf1 = perf(1);
    % Condition 2 
    perf2 = perf(2);
    % Condition 3
    perf3 = perf(3);
    % Condition 4;
    perf4 = perf(4);
       
    % Ecart condition 2 et 3
    ecart_perf = perf(2) - perf(3);

    Mdiff_perf = [Mdiff_perf, ecart_perf];
    
 % Maintask : calcul switching point (choice=0 : EC / choice=1 : jeu)
    % Condition 1
    sum1=0;
    for i=22:31
        sum1 = sum1 + data.maintask.trials(i).choice;
    end
    ec1 = data.maintask.trials(23).price + sum1 - 0.5;
    
	% Condition 2
    sum2=0;
    for i=53:62
        sum2 = sum2 + data.maintask.trials(i).choice;
    end
    ec2 = data.maintask.trials(54).price + sum2 - 0.5;
    
    % Condition 3
    sum3=0;
    for i=84:93
        sum3 = sum3 + data.maintask.trials(i).choice;
    end
    ec3 = data.maintask.trials(85).price + sum3 - 0.5;
    
    % Condition 4
    sum4=0;
    for i=115:124
        sum4 = sum4 + data.maintask.trials(i).choice;
    end
    ec4 = data.maintask.trials(116).price + sum4 - 0.5;
    
    Mec0_temp{k} = [ec1, ec2, ec3, ec4];
    
    % Ecart condition 2 et 3
    ecart_ec = ec2 - ec3;
    
    Mdiff_ec = [Mdiff_ec, ecart_ec];
    
    % preference pour le controle (si convergence des performances) :
    ec2 > ec3;

% Performances subjectives 
    % Condition 1
    subj1 = data.successrate(1).r;
    % Condition 2
    subj2 = data.successrate(2).r;
    % Condition 3
    subj3 = data.successrate(3).r;
    % Condition 4 
    subj4 = data.successrate(4).r;
    
% Verification : inversion des conditions 1&3 et 2&4
    subj1 > subj3;
    subj2 > subj4;
    % Interversion en cas d'erreur
    if inter == 1
        if subj1 < subj3
            temp = subj1;
            subj1 = subj3;
            subj3 = temp;
        end
        if subj2 < subj4
            temp = subj2;
            subj2 = subj4;
            subj4 = temp;
        end
    end
    
    % Ecart condition 2 et 3
    ecart_subj = subj2 - subj3;
    
    Mdiff_subj = [Mdiff_subj, ecart_subj]; 
        
% Ecarts performances objectives & performances subjectives (signe attendu
% aux condition 2 et 3 : -, car sous-estimation en general)
    % Condition 1
    ecart1 = perf(1) - data.successrate(1).r;
    % Condition 2
    ecart2 = perf(2) - data.successrate(2).r;
    % Condition 3
    ecart3 = perf(3) - data.successrate(3).r;
    % Condition 4 
    ecart4 = perf(4) - data.successrate(4).r;
    
 % preference pour le controle avec non convergence des performances
 %subjectives
    if compareToSubj == 1
        if subj2*100 - subj3*100 > ec2 - ec3
            pref_control = pref_control + 1;
            Mpref_control = [Mpref_control, 1];
        else
            Mpref_control = [Mpref_control, 0];
        end
    
 % preference pour le controle avec non convergence des performances
 %objectives
    else
        if perf2*100 - perf3*100 > ec2 - ec3
            pref_control = pref_control + 1;
            Mpref_control = [Mpref_control, 1];
        else
            Mpref_control = [Mpref_control, 0];
        end
    end

% 
% % Suppression des éléments vides de Mec0_temp
% for i = 1:length(Mec0_temp)
%     if length(Mec0_temp{i}) ~= 0
%         Mec0 = [Mec0, {Mec0_temp{i}}];
%     end
% end
% 
% pref_control_global = pref_control;
% Mpref_control_global = [Mpref_control];
% 
% Mpref_control = [];
% pref_control = 0;

pref_control_global = pref_control_global + pref_control;
Mpref_control_global = [Mpref_control_global, Mpref_control];
Mec_global = [Mec0, Mec1, Mec2, Mec3];
    
display(pref_control_global);  % Nombre de sujet ayant une preference pour le controle
pref_control_relative = pref_control_global/(80 - length(excluded_values0) - length(excluded_values1) - length(excluded_values2) - length(excluded_values3));    %Proportion de sujets ayant une preference pour le controle
display(pref_control_relative);  %Proportion de sujets sur l'echantillon exploitable ayant une preference pour le controle
display(Mpref_control_global);     % Vecteur dont le jeme element vaut 1 si l'individu j exprime une preference pour le controle
display(Mdiff_perf);     % Vecteur dont le j-ieme element est la différence des performances objectives de l'individu j : perf2 - perf3
display(Mdiff_subj);     % Vecteur dont le j-ieme element est la différence des performances subjectives de l'individu j : subj2 - subj3
display(Mdiff_ec);      % [Vecteur global] Liste dont le j-ieme element est la différence des equivalents certains de l'individu j : ec2 - ec3

% Vecteur aversion pour le risque : construction d'un vecteur dont le j-ieme
% element est le score d'aversion pour le risque. Le score va de 0 a  10, 10
% etant une grande aversion, 0 un guot pour le risque, 5 une indifference au
% risque. 

%Pilote
[num, txt, tab] = xlsread('Quest18072016.xlsx');

for i=1:length(tab)
    if strcmp(tab{i,11},'A')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if strcmp(tab{i,11},'B')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if length(temp_risk_aversion) == 11
        temp_risk_aversion(length(temp_risk_aversion)) = []; % Enleve dernier element de la liste (choix effectue par le sujet sur la ligne tiree aleatoirement par le programme en vue de la remuneration)
        for k=1:length(temp_risk_aversion)
            if strcmp(temp_risk_aversion(k),'A')
               risk_aversion = risk_aversion + 1;
            end
        end
        Mrisk_aversion0 = [Mrisk_aversion0, risk_aversion];
        temp_risk_aversion = [];
        risk_aversion = 0;
   end
end 
Mrisk_aversion0(18) = []; Mrisk_aversion0(6) = [];    % Suppression des outliers (sujet 7 : 6eme element du vecteur car sujet 4 absent ; sujet 20 : 18eme element du vecteur car sujet 13 absent)

% Vecteur homme=1/femme=0
for i=1:length(tab)
    if strcmp(tab{i,11},'Homme')
        Msexe0 = [Msexe0, 1];
    end
    if strcmp(tab{i,11},'Femme')
        Msexe0 = [Msexe0, 0];
    end
end
Msexe0(18) = []; Msexe0(6) = [];    % Suppression des outliers (sujet 7 : 6eme element du vecteur car sujet 4 absent ; sujet 20 : 18eme element du vecteur car sujet 13 absent)

% Expe1
[num, txt, tab] = xlsread('QuestGains25072016_13h.xlsx');
temp_risk_aversion = [];

for i=1:length(tab)
    if strcmp(tab{i,11},'A')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if strcmp(tab{i,11},'B')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if length(temp_risk_aversion) == 11
        temp_risk_aversion(length(temp_risk_aversion)) = []; % Enlève dernier élément de la liste (choix effectué par le sujet sur la ligne tirée aléatoirement par le programme en vue de la rémunération)
        for k=1:length(temp_risk_aversion)
            if strcmp(temp_risk_aversion(k),'A')
               risk_aversion = risk_aversion + 1;
            end
        end
        Mrisk_aversion1 = [Mrisk_aversion1, risk_aversion];
        temp_risk_aversion = [];
        risk_aversion = 0;
   end
end 
Mrisk_aversion1(9) = [];   % Suppression des outliers (sujet 10 : 9 ème élément du vecteur car sujet 6 absent)

% Vecteur homme=1/femme=0
for i=1:length(tab)
    if strcmp(tab{i,11},'Homme')
        Msexe1 = [Msexe1, 1];
    end
    if strcmp(tab{i,11},'Femme')
        Msexe1 = [Msexe1, 0];
    end
end
Msexe1(9) = [];  % Suppression des outliers (sujet 10 : 9 ème élément du vecteur car sujet 6 absent)

% Expe2
[num, txt, tab] = xlsread('QuestGains25072016_14h30.xlsx');
temp_risk_aversion = [];

for i=1:length(tab)
    if strcmp(tab{i,11},'A')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if strcmp(tab{i,11},'B')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if length(temp_risk_aversion) == 11
        temp_risk_aversion(length(temp_risk_aversion)) = []; % Enlève dernier élément de la liste (choix effectué par le sujet sur la ligne tirée aléatoirement par le programme en vue de la rémunération)
        for k=1:length(temp_risk_aversion)
            if strcmp(temp_risk_aversion(k),'A')
               risk_aversion = risk_aversion + 1;
            end
        end
        Mrisk_aversion2 = [Mrisk_aversion2, risk_aversion];
        temp_risk_aversion = [];
        risk_aversion = 0;
   end
end 
Mrisk_aversion2(16) = []; Mrisk_aversion2(13) = [];    % Suppression des outliers (sujet 15 : 13ème élément du vecteur car sujets 7 et 12 absents ; sujet 19 : 16ème élément du vecteur car sujets 7, 12 et 18 absents)

% Vecteur homme=1/femme=0
for i=1:length(tab)
    if strcmp(tab{i,11},'Homme')
        Msexe2 = [Msexe2, 1];
    end
    if strcmp(tab{i,11},'Femme')
        Msexe2 = [Msexe2, 0];
    end
end
Msexe2(16) = []; Msexe2(13) = [];  % Suppression des outliers (sujet 15 : 13ème élément du vecteur car sujets 7 et 12 absents ; sujet 19 : 16ème élément du vecteur car sujets 7, 12 et 18 absents)

% Expe3
[num, txt, tab] = xlsread('QuestGains28072016.xlsx');
temp_risk_aversion = [];

for i=1:length(tab)
    if strcmp(tab{i,11},'A')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if strcmp(tab{i,11},'B')
        temp_risk_aversion = [temp_risk_aversion, tab{i,11}];
    end
    if length(temp_risk_aversion) == 11
        temp_risk_aversion(length(temp_risk_aversion)) = []; % Enlève dernier élément de la liste (choix effectué par le sujet sur la ligne tirée aléatoirement par le programme en vue de la rémunération)
        for k=1:length(temp_risk_aversion)
            if strcmp(temp_risk_aversion(k),'A')
               risk_aversion = risk_aversion + 1;
            end
        end
        Mrisk_aversion3 = [Mrisk_aversion3, risk_aversion];
        temp_risk_aversion = [];
        risk_aversion = 0;
   end
end 
Mrisk_aversion3(5) = []; Mrisk_aversion3(1) = [];    % Suppression des outliers (sujet 6 : 5ème élément du vecteur car sujet 4 absent ; sujet 1 : 1er élément du vecteur)

% Vecteur homme=1/femme=0
for i=1:length(tab)
    if strcmp(tab{i,11},'Homme')
        Msexe3 = [Msexe3, 1];
    end
    if strcmp(tab{i,11},'Femme')
        Msexe3 = [Msexe3, 0];
    end
end
Msexe3(5) = []; Msexe3(1) = [];  % Suppression des outliers (sujet 6 : 5ème élément du vecteur car sujet 4 absent ; sujet 1 : 1er élément du vecteur)


Mrisk_aversion_global = [Mrisk_aversion0, Mrisk_aversion1, Mrisk_aversion2, Mrisk_aversion3];
display(Mrisk_aversion_global);
Msexe_global = [Msexe0, Msexe1, Msexe2, Msexe3];
display(Msexe_global);


%%
save('sourisfolle_globaldata.mat','Mdiff_ec','Mdiff_perf')

%% ANALYSES
clear all
close all
clc
cd('G:\client_sourisfolle\stats_globales')
load('sourisfolle_globaldata.mat')
addpath('tools')


%%

% Probit : variable dependante : Mpref_control / variables explicatives :
% Mrisk_aversion, Msexe
[b,dev,stats] = glmfit([Mrisk_aversion_global', Msexe_global'],Mpref_control_global','binomial','link','probit');


%% equivalent certain
Mec_global_mat = reshape([Mec_global{:}],[4 66])';
hist( Mec_global_mat )
hist(Mec_global_mat(:,2)- Mec_global_mat(:,3))

%% 
figure
plot(Mdiff_perf, Mdiff_subj, '.')
corr(Mdiff_perf', Mdiff_subj')

% refline % pour ajouter une ligne sur un plot

%%
for k=1:66
    if mod(k,20)==1
        figure
    end
    subplot(4,5,mod(k-1,20)+1)
    plot(x,y)
end





















