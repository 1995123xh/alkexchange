%% natural_exchange
% Simulating H isotope exchange in natural samples
% initialization
close all
clear all

load("all_betas.mat")
molecules = readtable('MolDescriptors_metal.xlsx', 'Range','E2:K28');
radicals = readtable('MolDescriptors_metal.xlsx','Range','M2:Q18');
%kmatrix=xlsread('MolDescriptors_metal.xlsx',2,'C4:M7');
HabstractFrac = xlsread('MolDescriptors_metal.xlsx',1,'W3:W17');
DabstractFrac = xlsread('MolDescriptors_metal.xlsx',1,'AB3:AB17');
Habstract = xlsread('MolDescriptors_metal.xlsx',1,'T3:U17');
Habstract_deg = xlsread('MolDescriptors_metal.xlsx',1,'V3:V17');
Dabstract = xlsread('MolDescriptors_metal.xlsx',1,'Y3:Z17');
Dabstract_deg = xlsread('MolDescriptors_metal.xlsx',1,'AA3:AA17');

% EIE at 100c
EIE=molecules.beta_100; 
allSpecies=height(molecules)+height(radicals)+2;  %add last two: metal H and metal D

% diagnal n*n matrix for making reaction balance
diagM=zeros(allSpecies,allSpecies);  
for i=1:allSpecies
    diagM(i,i)=1;
end

Habstract(:,2)=Habstract(:,2)+height(molecules);
Dabstract(:,2)=Dabstract(:,2)+height(molecules);

r=[];
n=1;

% row: species column:reaction
Mabsbalance=zeros(allSpecies, length(Habstract)*2+length(Dabstract)*2); 
reactantsList=zeros(4,length(Habstract)*2+length(Dabstract)*2);
for abst=1:length(Habstract)     %H-abstraction reaction
            Mabsbalance(:, n)=-diagM(:,Habstract(abst,1))+diagM(:,Habstract(abst,2))+diagM(:,allSpecies-1);
            reactantsList(1:2, n) = [Habstract(abst,1); Habstract(abst,2)];
            reactantsList(3, n) = Habstract_deg(abst);   %reactant degrees
            reactantsList(4, n) = HabstractFrac(abst);   %Fraction of reaction going down this path
            n=n+1;
end
for abst=1:length(Dabstract)     %D-abstraction reaction
            Mabsbalance(:, n)=-diagM(:,Dabstract(abst,1))+diagM(:,Dabstract(abst,2))+diagM(:,allSpecies);
            reactantsList(1:2, n) = [Dabstract(abst,1); Dabstract(abst,2)];
            reactantsList(3, n) = Dabstract_deg(abst);   %reactant degrees
            reactantsList(4, n) = DabstractFrac(abst);   %Fraction of reaction going down this path
            n=n+1;
end

for abst=1:length(Habstract)     %H-addition reaction
            Mabsbalance(:, n)=+diagM(:,Habstract(abst,1))-diagM(:,Habstract(abst,2))-diagM(:,allSpecies-1);
            reactantsList(1:2, n) = [Habstract(abst,1); Habstract(abst,2)];
            reactantsList(3, n) = Habstract_deg(abst);   %reactant degrees
            reactantsList(4, n) = HabstractFrac(abst);   %Fraction of reaction going down this path
            n=n+1;
end

for abst=1:length(Dabstract)     %H-addition reaction
            Mabsbalance(:, n)=+diagM(:,Dabstract(abst,1))-diagM(:,Dabstract(abst,2))-diagM(:,allSpecies);
            reactantsList(1:2, n) = [Dabstract(abst,1); Dabstract(abst,2)];
            reactantsList(3, n) = Dabstract_deg(abst);   %reactant degrees
            reactantsList(4, n) = DabstractFrac(abst);   %Fraction of reaction going down this path
            n=n+1;
end




genesis_conc = [84.37	7.52	4.33	1.39  0.75		0.354	0.389]; % 1,2,3,n4,i4,n5,i5; Genesis A15 ST1 from Thiagarajan 2020 paper
genesis_iso = [-209.1	-150	-122.7	-116  -120		-107.3  -109.3	]; % 1,2,3,n4,i4,n5,i5; Genesis A15 ST1 from Thiagarajan 2020 paper


Temp = 373;

jw_conc = genesis_conc;
jw_iso = genesis_iso;

c0=zeros(allSpecies,1); %initial condition
compo = [jw_conc(1), 0, jw_conc(2:end)];
compo = compo/sum(compo);
Dratio_ini=0.0001567; %D/H ratio
Dmodifier = [jw_iso(1), 0, jw_iso(2:end)];
Dmodifier = 1 + Dmodifier/1000;
Cratio=0.01;
atm=10*2.5e19;  %molecuels/cm3
% compo=[0.00,0.000, 0.3, 0.4,0.5,0.5,0.6,0.6];

c0(1)=atm*compo(1);
c0(2)=c0(1)*Dratio_ini*Dmodifier(1)*4;
c0(3)=c0(1)*(Dratio_ini*Dmodifier(1))^2*6;
c0(4)=c0(1)*Cratio;
c0(5)=c0(1)*Cratio*Dratio_ini*Dmodifier(1)*4;
c0(6)=c0(1)*Cratio*(Dratio_ini*Dmodifier(1))^2*6;
CSnum=2;
for i=7:height(molecules)
    if molecules.CStructure(i)~=CSnum  %change to the next molecule
        CSnum=CSnum+1;
        c0(i)=compo(CSnum)*atm;
        basesym=molecules.x__sym_(i);
    else
        c0(i)=compo(CSnum)*atm*basesym/molecules.x__sym_(i)*Dratio_ini*Dmodifier(CSnum);
    end
end
c0(3)=c0(3)*1.05;
c0(5)=c0(5)*1.05;

sf=0.23;
k=1e-5;

c0(height(molecules)+1:height(molecules)+height(radicals))=1e10;
c0(allSpecies-1)=1e15;
c0(allSpecies)=c0(allSpecies-1)*Dratio_ini;


%% run the model

tspan=[0 5e8];

options1 = odeset('Refine',1);
options2 = odeset(options1,'NonNegative',1); 
options3 = odeset(options2, 'RelTol', 1e-10);
options4 = odeset(options3, 'AbsTol', 1e-15);
options5 = odeset(options4);
options6 = odeset(options5,'MaxStep',5e8);


[time, c_out]=ode15s(@(t,c) dcdt_metal(t,c, Mabsbalance, reactantsList, EIE, sf,k, Temp), tspan,c0, options6);

timec=0.00001;%time compressor
time=time*timec;

%% plot

methanedD=1/4*c_out(2,:)/c_out(1,:)/Dratio_ini*1000-1000;

CSnum=0;
Hnumbers=[4,4,6,8,10,10,12,12];
Dratio=[];
subscount=0;
for i=1:height(molecules)
    if molecules.CStructure(i)~=CSnum  %change to the next molecule
        if CSnum>0
        newDratio=subscount./(Basemolecules*Hnumbers(CSnum));
        Dratio=[Dratio, newDratio];
        end
        CSnum=CSnum+1;
        subscount=zeros(size(c_out,1),1);
        Basemolecules=c_out(:,i);
    elseif i==3 || i==6
        subscount=subscount+c_out(:,i)*2;  %%doubly deuterated;
        Basemolecules=Basemolecules-1;
    else
        subscount=subscount+c_out(:,i);
    end
end
newDratio=subscount./(Basemolecules*Hnumbers(CSnum));
Dratio=[Dratio, newDratio];

figure()
hold on
for i=1:size(Dratio,2)
    plot(time,Dratio(:,i)/Dratio_ini*1000-1000)
    hold on
end
legend('C1','C1(13C)','C2','C3','nC4','iC4','nC5','iC5')
ax=gca;
% ax.XScale='log';

deltaD = Dratio/Dratio_ini*1000-1000;


% calculate the equi. curves
equil_temp = 101;
C2_C1_alphaD_200 = CSeps.C2_C1_D(equil_temp)/1000 + 1;
C3_C2_alphaD_200 = (CSeps.C3_C2_D(equil_temp)/1000 + 1);
C4n_C3_alphaD_200 = (CSeps.C4n_all_C1_D(equil_temp)/1000 + 1)./(CSeps.C3_C1_D(equil_temp)/1000 + 1);
C4i_C4n_alphaD_200 = (CSeps.C4i_C1_D(equil_temp)/1000 + 1)./(CSeps.C4n_all_C1_D(equil_temp)/1000 + 1);
C5i_C5n_alphaD_200 = (CSeps.C5i_all_C1_D(equil_temp)/1000 + 1)./(CSeps.C5n_all_C1_D(equil_temp)/1000 + 1);

alphas = [C2_C1_alphaD_200, C3_C2_alphaD_200, C4i_C4n_alphaD_200, C5i_C5n_alphaD_200];



figure();

% Define the column pairs for each subplot
columnPairs = [1, 3; 3, 4; 5, 6; 7, 8];

% Define the labels for each axis
xLabels = {'C1', 'C2', '{\it n}C4', '{\it n}C5'};
yLabels = {'C2', 'C3', '{\it i}C4', '{\it i}C5'};

% indices for ploting gom data
columnPairs_gom = [11, 12; 12, 13; 15, 14; 17, 16];

% Loop through each subplot
for i = 1:4
    subplot(2, 2, i);
    s = scatter(deltaD(1, columnPairs(i, 1)), deltaD(1, columnPairs(i, 2)), 'filled');
    s.MarkerEdgeColor = 'k';
    s.SizeData = 50;
    hold on;
    plot(deltaD(:, columnPairs(i, 1)), deltaD(:, columnPairs(i, 2)), 'LineWidth', 1.7);
    xlabel(xLabels{i});
    ylabel(yLabels{i});
    ax = gca;
    xlims_backup = [ax.XLim(1)-20, ax.XLim(2)+20];
    ylims_backup = [ax.YLim(1)-20, ax.YLim(2)+20];
    equicurve_x = linspace(min(deltaD(:, columnPairs(i, 1)))-30, max(deltaD(:, columnPairs(i, 1)))+30, 100);
    equicurve_y = (equicurve_x/1000 + 1) * alphas(i)*1000 - 1000;
    plot(equicurve_x, equicurve_y, 'k-.', LineWidth = 0.8)
    box on
    ax.LineWidth = 1;
    ax.XLim = xlims_backup;
    ax.YLim = ylims_backup;
end
