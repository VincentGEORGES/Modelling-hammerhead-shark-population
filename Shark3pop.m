%%Projet MEM 
%%
%Modélisation des Populations de Requins marteau halicornes ; approche
%matricielle
%
%Anais Martin - Vincent GEORGES
clear all
close all
% Vecteur temps
ti = 1;
dt = 5;
tf = 300;
time = ti:dt:tf;
tt = size(time, 2);

%Matrices de stockage
N = zeros(21,tt); %Abondances
test = zeros(10,tt);%Debugging

%Conditions initiales
N(:,1) = [5000 500 50 5 0 0 0 5000 500 50 5 0 0 0 5000 500 50 5 0 0 0];

%Paramètres basiques

%Facteur K logistique : par lagon 
Ka= 50000;
Kb = 50000;
Kc = 50000;
%Pop a

Pa = 1; % Impact de pèche sur a(0:1), multiplie le taux de survie : 1: pas de pèche, 0 Eradication par la pèche.
sa = 0.7; % Taux de survie basique annuel
sdS1a = 0.02; % ecart type de variation stochastique des Juvéniles
sdSxa = 0.05; %des autres

Pb = 1; % Pareil pour b
sb = 0.6;
sdS1b = 0.05;
sdSxb = 0.01;

Pc = 1;
sc = 0.5;
sdS1c = 0.05;
sdSxc = 0.01;

%taux de survie théorique : s'implante dans la boucle, puissance 5 car on
%intégre sur 5 ans 
s1a_t =(sa^5)*0.5 ;  %mortalité juvénile
s2a_t = sa^5 *Pa;
s3a_t = sa^5 *Pa;
s4a_t = sa^5*1.1*Pa; %1.1 mortalité diminue chez les gros
s5a_t = sa^5*1.1*Pa; 
s6a_t = sa^5*1.1*Pa;
s7a_t = sa^5*1.1*Pa;

s1b_t=(sb^5)*0.5;
s2b_t= sb^5*Pb;
s3b_t= sb^5*Pb;
s4b_t= sb^5*1.1*Pb;
s5b_t= sb^5*1.1*Pb;
s6b_t= sb^5*1.1*Pb;
s7b_t= sb^5*1.1*Pb;

s1c_t=(sc^5)*0.5;
s2c_t= sc^5*Pc;
s3c_t= sc^5*Pc;
s4c_t= sc^5*1.1*Pc;
s5c_t= sc^5*1.1*Pc;
s6c_t= sc^5*1.1*Pc;
s7c_t= sc^5*1.1*Pc;

% fecundity a
fa = 1.8; %fécondité basique

f2a_t= fa*5;
f3a_t= fa*5; % *5 car intégre sur 5 ans
f4a_t= fa*5;
f5a_t= fa*5*1.1;
f6a_t= fa*5*1.1;
f7a_t= fa*5*1.1;
sdfa = 0.005;

%fécondité  b
fb = 1.8;

f2b_t= fb*5;
f3b_t= fb*5;
f4b_t= fb*5;
f5b_t= fb*5*1.1;
f6b_t= fb*5*1.1;
f7b_t= fb*5*1.2;

fc = 1.8;

f2c_t= fc*5;
f3c_t= fc*5;
f4c_t= fc*5;
f5c_t= fc*5*1.1;
f6c_t= fc*5*1.1;
f7c_t= fc*5*1.2;

%Taux de migration dans les deux directions
mab = 0;
mba = 0;
mac = 0;
mca = 0;
mbc = 0;
mcb = 0;

mabJ= mab*0.8;   %Plus grande migration chez les individus agés
mabV= mab*1.2;

mbaJ= mba*0.8;
mbaV= mab*1.2;

macJ = mac*0.8;
macV = mac*1.2;

mcaJ = mca*0.8;
mcaV = mca*1.2;

mbcJ = mbc*0.8;
mbcV = mbc*1.2;

mcbJ = mcb*0.8;
mcbV = mcb*1.2;


%%
for i=1:(tt-1)

%Facteur logistique : limitation par les juvéniles dans les lagons
%Facteur limitant dans ce modèle

popJuva = N(1,i);
popJuvb = N(8,i);
popJuvc = N(15,i);

%Calcul des taux de survies réels; facteur stochastique groupée par
%classes

s1a = max(estoch(s1a_t,sdS1a,0.00001,1)* (1-(popJuva/Ka)),0); %application de la logistique à cette classe
s2a = max(estoch(s2a_t,sdSxa,0.00001,1) ,0);
s3a = s2a;
s4a = s2a;
s5a = max(estoch(s5a_t,sdSxa,0.01,1),0);
s6a = s5a;
s7a = s5a;

%Survie B
s1b = max(estoch(s1b_t,sdS1b,0.0001,1)*(1-(popJuvb/Kb)),0);
s2b = max(estoch(s2b_t,sdSxb,0.0001,1),0);
s3b = s2b;
s4b = s2b;
s5b = max(estoch(s5b_t,sdSxb,0.0001,1),0);
s6b = s5b;
s7b = s5b;

%Survie C

s1c = max(estoch(s1c_t,sdS1c,0.0001,1)*(1-(popJuvc/Kc)),0);
s2c = max(estoch(s2c_t,sdSxc,0.0001,1),0);
s3c = s2c;
s4c = s2c;
s5c = max(estoch(s5c_t,sdSxc,0.0001,1),0);
s6c = s5c;
s7c = s5c;

% fecondité a

f2a= f2a_t;
f3a= f3a_t;
f4a= f4a_t;
f5a= f5a_t;
f6a= f6a_t;
f7a= f7a_t;

%fécondité  b

f2b= f2b_t;
f3b= f3b_t;
f4b= f4b_t;
f5b= f5b_t;
f6b= f6b_t;
f7b= f7b_t;

%fécondité  C

f2c= f2c_t;
f3c= f3c_t;
f4c= f4c_t;
f5c= f5c_t;
f6c= f6c_t;
f7c= f7c_t;


%Debug matrix
%Garde en mémoires les taux tirés au hasards ; permet de surveiller en cas
%d'emballement soudain de la machine
test(1,i) = s1a;
test(2,i) = s2a;
test(3,i) = s5a;
test(4,i) = s1b;
test(5,i) = s2b;
test(6,i) = s5b;
test(7,i) = s2b_t;
test(8,i) = s1c;
test(9,i) = s2c;
test(10,i)= s5c; 


%Création des matrices de transition

vFa =[0 f2a f3a f4a f5a f6a f7a]; % Vecteur fecondité
vSa =[s1a (1-macJ)*(1-mabJ)*s2a (1-macJ)*(1-mabJ)*s3a (1-macJ)*(1-mabJ)*s4a (1-macV)*(1-mabV)*s5a (1-macV)*(1-mabV)*s6a ]; %Vecteur survie
Ma = diag(vSa,-1);
Ma(1,:) = vFa;

vFb =[0 f2b f3b f4b f5b f6b f7b];
vSb =[s1b (1-mbcJ)*(1-mbaJ)*s2b (1-mbcJ)*(1-mbaJ)*s3b (1-mbcJ)*(1-mbaJ)*s4b (1-mbcV)*(1-mbaV)*s5b (1-mbcV)*(1-mbaV)*s6b ];
Mb = diag(vSb,-1);
Mb(1,:) = vFb;

vFc =[0 f2c f3c f4c f5c f6c f7c];
vSc =[s1c (1-mcbJ)*(1-mcaJ)*s2c (1-mcbJ)*(1-mcaJ)*s3c (1-mcbJ)*(1-mcaJ)*s4c (1-mcbV)*(1-mcaV)*s5c (1-mcbV)*(1-mcaV)*s6c];
Mc = diag(vSc,-1);
Mc(1,:) = vFc;

vFba = [0 0 0 0 0 0 0];
vSba = [0 mbaJ*s2b mbaJ*s3b mbaJ*s4b mbaV*s5b mbaV*s6b mbaV*s7b];
Mba = diag(vSba,-1);
Mba(:,8) = []; 
Mba(1,:) = [];
Mba(1,:) = vFba;

vFbc = [0 0 0 0 0 0 0];
vSbc = [0 mbcJ*s2b mbcJ*s3b mbcJ*s4b mbcV*s5b mbcV*s6b mbcV*s7b];
Mbc = diag(vSbc,-1);
Mbc(:,8) = []; 
Mbc(1,:) = [];
Mbc(1,:) = vFbc;


vFab = [0 0 0 0 0 0 0];
vSab = [0 mabJ*s2a mabJ*s3a mabJ*s4a mabV*s5a mabV*s6a mabV*s7a];
Mab = diag(vSba,-1);
Mab(:,8) = []; 
Mab(1,:) = [];
Mab(1,:) = vFab;

vFac = [0 0 0 0 0 0 0];
vSac = [0 macJ*s2b macJ*s3b macJ*s4b macV*s5b macV*s6b macV*s7b];
Mac = diag(vSac,-1);
Mac(:,8) = []; 
Mac(1,:) = [];
Mac(1,:) = vFac;

vFca = [0 0 0 0 0 0 0];
vSca = [0 mcaJ*s2c mcaJ*s3c mcaJ*s4c mcaV*s5c mcaV*s6c mcaV*s7c];
Mca = diag(vSca,-1);
Mca(:,8) = []; 
Mca(1,:) = [];
Mca(1,:) = vFca;

vFcb = [0 0 0 0 0 0 0];
vScb = [0 mcbJ*s2c mcbJ*s3c mcbJ*s4c mcbV*s5c mcbV*s6c mcbV*s7c];
Mcb = diag(vScb,-1);
Mcb(:,8) = []; 
Mcb(1,:) = [];
Mcb(1,:) = vFcb;

M = [Ma Mba Mca; %Assemblage
     Mab Mb Mcb;
     Mac Mbc Mc];
%Integration
%sur le temps
N(:,i+1) = M * N(:,i);

%Petit message d'avancement
X = [num2str((tt-1)-i),' left to compute'];
    disp(X)


end

%Plots

figure(1)
subplot(3,1,1)
plot(time, N(1,:),'r-+')
hold on
plot(time, sum(N(2:4,:)),'r-o')
hold on
plot(time, sum(N(5:7,:)),'r-s')
hold on
plot(time, N(8,:),'g-+')
hold on
plot(time, sum(N(9:11,:)),'g-o')
hold on
plot(time, sum(N(12:14,:)),'g-s')
hold on
plot(time, N(15,:),'b-+')
hold on
plot(time, sum(N(16:18,:)),'b-o')
hold on
plot(time, sum(N(19:21,:)),'b-s')

legend('Bébés a','Jeunes a','Vieux a','Bébés b','Jeunes b','Vieux b','Bébés c','Jeunes c','Vieux c')
xlabel('Years');
ylabel('Abundance');
x = [' mab = ',num2str(mab),' mba =',num2str(mba),' Pa = ',num2str(Pa),' Pb = ',num2str(Pb)];
title(x)

subplot(3,1,2)
semilogy(time, N(1,:),'r-+')
hold on
semilogy(time, sum(N(2:4,:)),'r-o')
hold on
semilogy(time, sum(N(5:7,:)),'r-*')
hold on
semilogy(time, N(8,:),'g-+')
hold on
semilogy(time, sum(N(9:11,:)),'g-o')
hold on
semilogy(time, sum(N(12:14,:)),'g-*')
hold on
semilogy(time, N(15,:),'b-+')
hold on
semilogy(time, sum(N(16:18,:)),'b-o')
hold on
semilogy(time, sum(N(19:21,:)),'b-*')
legend('Bébés a','Jeunes a','Vieux a','Bébés b','Jeunes b','Vieux b','Bébés c','Jeunes c','Vieux c')
xlabel('Years');
ylabel('Abundance');
x = [' mab = ',num2str(mab),' mba =',num2str(mba),' Pa = ',num2str(Pa),' Pb = ',num2str(Pb)];
title(x)

subplot(3,1,3)
plot(time,sum(N))
hold on
plot(time,sum(N(1:7,:)),'r-+')
hold on
plot(time,sum(N(8:14,:)),'g-o')
hold on 
plot(time,sum(N(15:21,:)),'b-*')

legend('total abundance','pop a ', 'popb', 'pop c')

figure(2)
plot(time,test)
legend('s1a','s2a','s5a','s1b','s2b','s5b','s1c','s2c','s5c')



%Fonction de stochasticité
function y = estoch(mu, sigma, lower, upper)
% ESTOCH simulate environmental stochasticity for a survival/fecundity rate
%   y = estoch(mu, sigma, lower, upper) generates the new value which is
%   centered around the old value, mu, with a variability sigma and bounded
%   in [lower, upper]. For survival lower=0, upper=1, for fecundity
%   lower=0, upper=Inf
    normd = makedist('Normal', mu, sigma);
    tnormd = truncate(normd, lower, upper);
    x = random(tnormd,1, 1);
    if x > 1
        y = 1;
    else
          y = x;
    end
end