%%  Application de la SVD : compression d'images

clear all
close all

% Lecture de l'image
I = imread('BD_Asterix_1.png');
I = rgb2gray(I);
I = double(I);

[q, p] = size(I)

% Décomposition par SVD
fprintf('Décomposition en valeurs singulières\n')
tic
[U, S, V] = svd(I);
toc

l = min(p, q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On choisit de ne considérer que 200 vecteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vecteur pour stocker la différence entre l'image et l'image reconstuite
inter = 1:40:(200 + 40);
inter(end) = 200;
differenceSVD = zeros(size(inter, 2), 1);

% images reconstruites en utilisant de 1 à 200 vecteurs (avec un pas de 40)
ti = 0;
td = 0;
for k = inter

    % Calcul de l'image de rang k
    Im_k = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';

    % Affichage de l'image reconstruite
    ti = ti + 1;
    figure(ti)
    colormap('gray')
    imagesc(Im_k), axis equal

    % Calcul de la différence entre les 2 images
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I - Im_k).^2)));
    % pause
end

% Figure des différences entre image réelle et image reconstruite
ti = ti + 1;
figure(ti)
hold on
plot(inter, differenceSVD, 'rx')
ylabel('RMSE')
xlabel('rank k')
pause


% Plugger les différentes méthodes : eig, puissance itérée et les 4 versions de la "subspace iteration method"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUELQUES VALEURS PAR DÉFAUT DE PARAMÈTRES,
% VALEURS QUE VOUS POUVEZ/DEVEZ FAIRE ÉVOLUER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolérance
eps = 1e-3;
% nombre d'itérations max pour atteindre la convergence
maxit = 10000;

% taille de l'espace de recherche (m)
search_space = 400;

% pourcentage que l'on se fixe
percentage = 0.995;

% p pour les versions 2 et 3 (attention p déjà utilisé comme taille)
puiss = 1;

% version a utiliser
v = 1;

%%%%%%%%%%%%%
% À COMPLÉTER
%%%%%%%%%%%%%

pause
close all

figure
tic

M = I * I';

t = tiledlayout(3, 3);

%%
% calcul des couples propres
switch v
    case 0

        %% eig method
        [U, D] = eig(M);
        [D, idx] = sort(diag(D), 'descend');

        U = U(:, idx);
        U = U(:, 1:search_space);

        fprintf('\nTemps EIG\n');
        title(t, 'Results using EIG');
    case 1
        %% power method

        [ U, D, n_ev, itv, flag, ] = power_v12(M, search_space, percentage, eps, maxit);
        U = U(:,1:n_ev);

        fprintf('\nTemps Power Method\n');
        title(t, 'Results using the Power Method');
end

%% calcul des valeurs singulières

if isdiag(D)
    Sigma = sqrt(D);
else
    Sigma = diag(sqrt(D(1:search_space)));
end

%% calcul de l'autre ensemble de vecteurs

V = I' * U .* repmat(1./diag(Sigma)', p, 1);

%% calcul des meilleures approximations de rang faible
difference = zeros(size(inter, 2), 1);

td = 0;

toc 

for k = inter

    % Calcul de l'image de rang k
    Im_k = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';

    % Affichage de l'image reconstruite
    ti = ti + 1;
    nexttile;
    colormap('gray')
    imagesc(Im_k), axis equal
    title(['k = ' num2str(k)]);

    % Calcul de la différence entre les 2 images
    td = td + 1;
    difference(td) = sqrt(sum(sum((I - Im_k).^2)));
    %pause
end

% Figure des différences entre image réelle et image reconstruite
figure;
hold on
plot(inter, difference, 'rx')
ylabel('RMSE')
xlabel('rank k')
pause
