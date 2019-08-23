%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  CLASSIFICAÇÃO (K-MEDIAS)  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all

%% Carrega Dado
X = imread ('palito.tif');
[Lin Col] = size(X);
k = 3; %Número de Clusters

M = zeros(k,Col);
for i=1:k 
  M(i,:)= X(i,:) +i ;    
end  

%%
Mnew1 = M;
E_distance = zeros(Lin,k);


for i = 1:k
Mnew(:,:,i) = repmat(Mnew1(i,:),Lin,1);
E_distance(:,i) = sqrt(sum((double(X) - Mnew(:,:,i)).^2,2));
end 

S = zeros(Lin,k);
[M , I] = (min(E_distance(:,:),[],2));  %Retorna os min val com os ind. da col

%% Número de amostras por clusters
NumberOfsamples = zeros(k,1);
for i = 1: Lin
    for j = 1: k 
      if (I(i,1) == j)
        NumberOfsamples(j) = NumberOfsamples(j)+1;    
        end
    end
end

for i=1:k
Indx{i} = zeros(NumberOfsamples(i),1);
end

for j = 1:k  
Indx{j} = find(I(:,1)==j);      
end

%% calculation for the sum of samples from same cluster
Mean = zeros(k,Col); 
for j= 1:k
temp = Indx{j};

    for i=1:NumberOfsamples(j)
    temp1(i,:) = X(temp(i),:);
    Mean(j,:) = double(temp1(i,:)) + Mean(j,:);    
    end
end

for j = 1:k
    Mnew2(j,:)= (1/NumberOfsamples(j)).*(Mean(j,:));
end

Mnew1 = Mnew2 ;


Means = Mnew1 ;
Samplesfromeachclass = NumberOfsamples;
