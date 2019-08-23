%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%SEGMENTAÇÃO POR CLASSIFICAÇÃO (K-MEDIAS)%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
erro=0;

%% Carregar a imagem
img_1 = imread('palito.tif');
%img_1 = rgb2gray(img_1);
img_2 = double(img_1);
[lin1, col1] = size(img_2);

%% Adicionar Ruídos na Imagem
ruido = random('norm',2,1,[lin1, col1]);
img_ruido = ruido + img_2;

%% Inicialização dos K-centroides
k = 3; % Número de Clusters
for i=1:k
c(i) = img_ruido(randi(lin1,1,1),randi(lin1,1,1));
end

%% Cálculo das Distâncias pro Centroides

dc1 = arrayfun(@(x) abs(x - c(1)), img_ruido);
dc2 = arrayfun(@(x) abs(x - c(2)), img_ruido);
dc3 = arrayfun(@(x) abs(x - c(3)), img_ruido);


% Matriz Distâncias
dist = [dc1(:), dc2(:), dc3(:)];
[~, ind] = min(dist, [], 2);

for l=1:k
[L(l), C(l)] = ind2sub(size(img_ruido), find(ind == l));
end


%% Inicialização
g1 = zeros(lin1,col1);



figure(1),imshow(img_1)
title('\fontsize{12}Imagem  Original')

figure(2),imshow(ruido)
title('\fontsize{12}Ruido Adicional')

figure(3),imshow(g1)
title('\fontsize{12}Imagem  Segmentada')
