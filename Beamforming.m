%% ----------------- Projet de traitement du son : Groupe 11 ---------------------- %%
% Sujet : Localisation de sources sonores
% Réalisé par : 
%   - ALOULEN Bilal
%   - BANI Mohamed-Iadh
%   - BEN FTIMA Mohamed
%   - DA ROCHA Martins Mickaël


%% Fermeture des figure et suppression des variables
clc
clear all
close all

%% Initialisation des constantes
c = 334;         % Vitesse du son
fs = 1000;       % Fréquence d'échantillonnage
T = 0.1;         % Période
t = 0:1/fs:T;    % Intervalle de temps [0, 0.1]
L = length(t);   % Longueur du vecteur temps = 101
f = 500;         % Fréquence des signaux
w = 2*pi*f;      % Pulsation
k = w/c;         % Nombre d'onde k

%% Coordonnées des microphones
M = 13;          % Nombre des microphones
yi = zeros(M,1); % Définition d'une matrice de zéros dimension M*1
% Coordonnées en z des micros
zi = [9;15;12;12;12;12;12;11;13;10;14;12;12];
% Coordonnées en x des micros
xi = [12;12;9;15;12;11;13;12;12;12;12;10;14];

% Figure montrant l'emplacement des microphones
figure(1)
plot(xi,zi,'b*');
xlabel('x'),ylabel('z')
title('Réseau de microphones')


%% Coordonnées des sources sonores
% Coordonnées de la source sonore
x1 = 13;
y1 = 20;
z1 = 13;
% Coordonnées de l'antenne constituée de l'ensemble des microphones
x2 = 12;
y2 = 0;
z2 = 12;
 
Ric1 = sqrt((x1-xi).^2+(y1-yi).^2+(z1-zi).^2);     % Distance entre la source sonore et chaque microphone
Ric2 = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);     % Distance entre la source et le centre de l'antenne
Rn1 = Ric1 - Ric2;                                 

s1 = cos(2*w*t);                                   % Le signal reçu par le microphone de référence
Am = 10^(-1);                                      % Amplitude du signal
n1 = Am * (randn(M, L) + j*randn(M, L));           % Bruit blanc gaussien de chaque micro
p1 = zeros(M,L);


%% Calcul du retard de chaque microphone
% Calcul de la matrice de pression acoustique reçue par chaque microphone
for k1 = 1:M
    p1(k1,:) = Ric2/Ric1(k1) * s1.*exp(-j*w*Rn1(k1)/c);
end

p = p1+n1;  % Matrice des signaux de pression acoustique reçus par chaque microphone
R = p*p'/L; % Matrice d'autocovariance des données reçues  


%% Plage de balayage
step_x = 0.1;       % Taille de pas
step_z = 0.1;
y = y1;
x = (9:step_x:15); 
z = (9:step_z:15); 

for k1=1:length(z)
    for k2=1:length(x)
        Ri = sqrt((x(k2)-xi).^2+(y-yi).^2+(z(k1)-zi).^2);  % Vecteur de distance
        Ri2 = sqrt((x(k2)-x2).^2+(y-y2).^2+(z(k1)-z2).^2); 
        Rn = Ri-Ri2;   % Vecteur de différence du point de balayage à chaque élément et l'élément de référence
        b = exp(-j*w*Rn/c); % Vecteur de direction de mise au point de pression acoustique
        Pcbf(k1,k2) = abs(b'*R*b);
    end
end


%% Normalisation
for k1 = 1:length(z)
    pp(k1) = max(Pcbf(k1,:)); % Pcbf La valeur du plus grand élément de la ligne k1
end

Pcbf = Pcbf/max(pp);  % Division de tous les éléments par leur valeurmaximale


%% Affichages
figure(2)
surf(x,z,Pcbf);
xlabel('x(m)'),ylabel('z(m)')
title('Diagramme de source sonore 3D')
colorbar
 
figure(3)
pcolor(x,z,Pcbf);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
title('Diagramme de source sonore 2D')
colorbar