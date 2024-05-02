%clear all;
format long;
clc
%%%%%%%%%%%%
% PARAMÈTRES
%%%%%%%%%%%%

% taille de la matrice symétrique
for n = 200:200:600

% type de la matrice (voir matgen_csad)
% imat == 1 valeurs propres D(i) = i
% imat == 2 valeurs propres D(i) = random(1/cond, 1) avec leur logarithmes
%                                  uniformément répartie, cond = 1e10
% imat == 3 valeurs propres D(i) = cond**(-(i-1)/(n-1)) avec cond = 1e5
% imat == 4 valeurs propres D(i) = 1 - ((i-1)/(n-1))*(1 - 1/cond) avec cond = 1e2
for imat = 1:4

% tolérance
eps = 1e-8;
% nombre d'itérations max pour atteindre la convergence
maxit = 1000;

% on génère la matrice (1) ou on lit dans un fichier (0)
genere = 0;

% méthode de calcul
v = 0; % subspace iteration v0

% nombre de valeurs propres cherchées (v0)
m = 20;

[W, V, flag, q, qv] = eigen_2024(imat, n, v, m, eps, maxit, [], [], genere);
fprintf('Qualité des couples propres (par rapport au critère d''arrêt) = [%0.3e , %0.3e]\n', min(qv), max(qv));
fprintf('Qualité des valeurs propres (par rapport au spectre de la matrice) = [%0.3e , %0.3e] \n', min(q), max(q));

% méthode de calcul
v = 1;

% pourcentage de la trace que l'on veut atteindre (v1, v2)
percentage = .1;

[W, V, flag, q, qv] = eigen_2024(imat, n, v, m, eps, maxit, percentage, [], genere);
fprintf('Qualité des couples propres (par rapport au critère d''arrêt) = [%0.3e , %0.3e]\n', min(qv), max(qv));
fprintf('Qualité des valeurs propres (par rapport au spectre de la matrice) = [%0.3e , %0.3e] \n', min(q), max(q));

% méthode de calcul
v = 2;

p = 8; %8 etant la meilleure valeur obtenue experimentalement
[W, V, flag, q, qv] = eigen_2024(imat, n, v, m, eps, maxit, percentage, p, genere);
fprintf('Qualité des couples propres (par rapport au critère d''arrêt) = [%0.3e , %0.3e]\n', min(qv), max(qv));
fprintf('Qualité des valeurs propres (par rapport au spectre de la matrice) = [%0.3e , %0.3e] \n', min(q), max(q));


% méthode de calcul
v = 3;

[W, V, flag, q, qv] = eigen_2024(imat, n, v, m, eps, maxit, percentage,p , genere);
fprintf('Qualité des couples propres (par rapport au critère d''arrêt) = [%0.3e , %0.3e]\n', min(qv), max(qv));
fprintf('Qualité des valeurs propres (par rapport au spectre de la matrice) = [%0.3e , %0.3e] \n', min(q), max(q));
end
end
