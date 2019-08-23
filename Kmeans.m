clc;
clear;
close all;
%% Carregar Dados
tic;
dado = load('grupo5.txt');


%% Separar 2 dados de entrada para plot
Xtrain = dado(:,1:8);
figure(1);
hold off
plot(Xtrain(:,4),Xtrain(:,5),'ro'); 
title 'ICC x Isquêmica x Chagas x Saudáveis';
xlabel 'SDNN';
ylabel 'SDANN';
Lambda = 1;

%% Número de Clusters
K = 2; %input('enter the value of K =');

N = size(Xtrain,1);
kernel = zeros(N); 

%% Cálculo da função Kernel
for i=1:N
    for j=1:N
        kernel(i,j) = exp(-Lambda*sum((Xtrain(i,:)-Xtrain(j,:)).^2));
    end
end

%% Separar a saída dado em binário 
converged = 0;
cont = 0;
Z_or = zeros(N,K);
Z_or(1:32,1) = 1;
Z_or(33:42,2) = 1;
%Z_or(30:39,3) = 1;
%Z_or(40:42,4) = 1;
Z = Z_or;

distance = zeros(N,K);


while ~converged 
    
    for k = 1:K
    Qtd_Saida(k)= sum(Z(:,k)); % Qtd de dados por padrão(saída)
    % Compute kernelised distance
    distance(:,k) = diag(kernel) - (2/(Qtd_Saida(k)))*sum(repmat(Z(:,k)',N,1).*kernel,2) + ...
     Qtd_Saida(k)^(-2)*sum(sum((Z(:,k)*Z(:,k)').*kernel));
    end
    oldZ = Z;
    Z = (distance == repmat(min(distance,[],2),1,K)); % Transf. saída em binário
    Z = 1.0 * Z ;   
    cont = cont + 1;

   if sum(sum(oldZ~=Z))==0
       converged = 1;
   end
end   

TrainClusterAssignments = Z ;

 for j = 1:K  
   Indx{j} = find(Z(:,j)==1);      
end

figure(2);
if (K==2)
%%Plot the data with the means
plot(Xtrain(Indx{1,1},4),Xtrain(Indx{1,1},5),'b.','MarkerSize',14 );
hold on
plot(Xtrain(Indx{1,2},4),Xtrain(Indx{1,2},5),'g.','MarkerSize',14 );
title 'ICC x Saudáveis';
xlabel 'SDNN';
ylabel 'SDANN';
legend('Cardiopatas','Saudável','Location','NW');
title ('K-means');
end
% 
% if (K==3)
% plot(Xtrain(Indx{1,1},4),Xtrain(Indx{1,1},5),'r.','MarkerSize',12 );
% hold on
% plot(Xtrain(Indx{1,2},4),Xtrain(Indx{1,2},5),'b.','MarkerSize',12 );
% hold on
% plot(Xtrain(Indx{1,3},4),Xtrain(Indx{1,3},5),'g.','MarkerSize' ,12);
% hold on
% xlabel 'SDNN';
% ylabel 'SDANN';
% legend('Cluster 1','Cluster 2','Cluster 3','Location','NW');
% title ('K-means');
% end
% 
if (K==4)
plot(Xtrain(Indx{1,1},4),Xtrain(Indx{1,1},5),'r.','MarkerSize',14);
hold on
plot(Xtrain(Indx{1,2},4),Xtrain(Indx{1,2},5),'b.','MarkerSize',14);
hold on
plot(Xtrain(Indx{1,3},4),Xtrain(Indx{1,3},5),'g.','MarkerSize' ,14);
hold on
plot(Xtrain(Indx{1,4},4),Xtrain(Indx{1,4},5),'m.','MarkerSize' ,14);
hold on
xlabel 'SDNN';
ylabel 'SDANN';
legend('Idiopática','Idiopática/Isquêmica','Saudáveis','Idiopática/Chagásica','Location','NW');
title ('K-means');
end

 convergenceStatus = converged;
 
 %% Resultados
cont
erro_saida = Z_or(:,1) - Z(:,1);
erroc = sum(erro_saida);
erro_quad = erroc/N;
tx_acerto = 100 - ((erroc/N)*100)
desviop = std(erro_saida)
tempo = toc