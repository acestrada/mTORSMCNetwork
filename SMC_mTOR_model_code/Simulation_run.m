
TSCKO = [1 0];
mTORC1block = [1 0];
tspan = [0 2000];
options = [];

n_val = 1.4; %Hill exponent
EC50val = 0.52; %Half-maximal activation
input_weight_stress = 0.24054; %Input stress

contractile_phenotype_species = [67 65 64]; %SMMHC, SMA, SM22
synthetic_phenotype_species = [10 16 73]; %Col3a1, Eln, TIMP
degradative_phenotype_species = [34 42 61 40 8]; %LAMP2, MMP2, S6, MITF, betaCatenin

%Qualitative comparison:
%Col3a1, Eln, S6K, S6, p4EBP1, AKT, betaCatenin, MITF, LAMP2, MMP2, SMMHC,
%SMA, SM22, Elastin
qualitative_interest_species = [10 16 62 61 49 1 8 40 34 42 67 65 64 14]; 

% Baseline condition
[params_base,y0_base] = SMC_mTOR_ODE_loadParams(n_val,EC50val,input_weight_stress,TSCKO(1),mTORC1block(1));
[t_base,y_base] = ode45(@SMC_mTOR_ODE,tspan,y0_base,options,params_base);
yf_base = y_base(end,:);

for nn = 1:length(t_base)
    Proliferation_base(nn) = mean([real(y_base(nn,8)) (1-real(y_base(nn,48))) real(y_base(nn,61))]);
    Phagocytosis_base(nn) = mean(real([y_base(nn,12) y_base(nn,34)]));
end

% Tsc1 KO condition
[params_KO,y0_KO] = SMC_mTOR_ODE_loadParams(n_val,EC50val,input_weight_stress,TSCKO(2),mTORC1block(1));
[t_KO,y_KO] = ode45(@SMC_mTOR_ODE,tspan,yf_base,options,params_KO);
yf_KO = y_KO(end,:);

for mm = 1:length(t_KO)
    Proliferation_KO(mm) = mean([real(y_KO(mm,8)) (1-real(y_KO(mm,48))) real(y_KO(mm,61))]);
    Phagocytosis_KO(mm) = mean(real([y_KO(mm,12) y_KO(mm,34)]));
end

% Rapamycin condition
[params_rapa,y0_rapa] = SMC_mTOR_ODE_loadParams(n_val,EC50val,input_weight_stress,TSCKO(2),mTORC1block(2));
[t_rapa,y_rapa] = ode45(@SMC_mTOR_ODE,tspan,yf_base,options,params_rapa);
yf_rapa = y_rapa(end,:);

for qq = 1:length(t_rapa)
    Proliferation_rapa(qq) = mean([real(y_rapa(qq,8)) (1-real(y_rapa(qq,48))) real(y_rapa(qq,61))]);
    Phagocytosis_rapa(qq) = mean(real([y_rapa(qq,12) y_rapa(qq,34)]));
end

Base_phenotypes = real([mean(yf_base(contractile_phenotype_species)); mean(yf_base(synthetic_phenotype_species));mean(yf_base(degradative_phenotype_species))]);
KO_phenotypes = real([mean(yf_KO(contractile_phenotype_species)); mean(yf_KO(synthetic_phenotype_species));mean(yf_KO(degradative_phenotype_species))]);
Rapamycin_phenotypes = real([mean(yf_rapa(contractile_phenotype_species)); mean(yf_rapa(synthetic_phenotype_species));mean(yf_rapa(degradative_phenotype_species))]);

% Qualitative comparison plotting
% Colormap
mymap = [5,48,97
    33,102,172
    67,147,195
    146,197,222
    209,229,240
    247,247,247
    253,219,199
    244,165,130
    214,96,77
    178,24,43
    103,0,31]./255;

qualitative_comparison_base = [real(yf_base(qualitative_interest_species)) Proliferation_base(end) real(yf_base(5)) Phagocytosis_base(end)];
qualitative_comparison_KO = [real(yf_KO(qualitative_interest_species)) Proliferation_KO(end) real(yf_KO(5)) Phagocytosis_KO(end)];
qualitative_comparison_rapa = [real(yf_rapa(qualitative_interest_species)) Proliferation_rapa(end) real(yf_rapa(5)) Phagocytosis_rapa(end)];

qualitative_comparison_labels = {'Col3a1','Eln','p-S6K','p-S6','p-4EBP1','p-AKT', '\beta-Catenin','MITF','LAMP2','MMP2','SMMHC','SMA','SM22','Elastin','Proliferation','Apoptosis','Phagocytosis'};

figure
imagesc([qualitative_comparison_KO-qualitative_comparison_base])
colormap(mymap)
caxis([-1 1])
xticks(1:length(qualitative_comparison_base))
xticklabels(qualitative_comparison_labels)
ylabel('\Delta(KO - Baseline)')

figure
imagesc([qualitative_comparison_rapa-qualitative_comparison_KO])
colormap(mymap)
caxis([-1 1])
xticks(1:length(qualitative_comparison_KO))
xticklabels(qualitative_comparison_labels)
ylabel('\Delta(Rapamycin - KO)')

figure
imagesc([qualitative_comparison_rapa-qualitative_comparison_base])
colormap(mymap)
caxis([-1 1])
xticks(1:length(qualitative_comparison_base))
xticklabels(qualitative_comparison_labels)
ylabel('\Delta(Rapamycin - Baseline)')


% Phenotypes plotting
X_base = [0 sqrt((Base_phenotypes(2)^2)/2) -sqrt((Base_phenotypes(3)^2)/2)];
Y_base = [Base_phenotypes(1) -sqrt((Base_phenotypes(2)^2)/2) -sqrt((Base_phenotypes(3)^2)/2)];
V_base = [X_base(1) Y_base(1);X_base(2) Y_base(2);X_base(3) Y_base(3)];
F_base = [1 2 3];
X_KO = [0 sqrt((KO_phenotypes(2)^2)/2) -sqrt((KO_phenotypes(3)^2)/2)];
Y_KO = [KO_phenotypes(1) -sqrt((KO_phenotypes(2)^2)/2) -sqrt((KO_phenotypes(3)^2)/2)];
V_KO = [X_KO(1) Y_KO(1);X_KO(2) Y_KO(2);X_KO(3) Y_KO(3)];

figure
hold on
patch('Faces',F_base,'Vertices',V_base,'FaceColor','blue','FaceAlpha',.3)
patch('Faces',F_base,'Vertices',V_KO,'FaceColor','red','FaceAlpha',.3)
legend('Baseline','KO')





