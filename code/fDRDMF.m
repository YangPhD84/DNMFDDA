function [M, Z, W] = fDRDMF(A,data,data_Similarity,layers,alpha,beta,lammad,gamma,MaxIter,tol1,tol2,T,Dataset)

%% initialization
%initialize—Z、W、Lr
[Z,W] = fdrinit(data, layers);

num_of_layers = numel(layers);
numOfView = numel(data);

%initialize—基于原始图的拉普拉斯矩阵L_W
K_W = cell(1,numOfView);
O_W = cell(1,numOfView);
L_W = cell(1,numOfView);%构建关于基矩阵W的拉普拉斯矩阵L_W

for v_ind = 1:numOfView
    K_W{1,v_ind} = diag(sum(data_Similarity{v_ind},2));
    O_W{1,v_ind} = data_Similarity{v_ind};
    L_W{1,v_ind} = K_W{v_ind} - O_W{v_ind};
end


WH = A;
stop1 = 1;
stop2 = 1;

%%%%
for iter = 1:MaxIter
%     iter
     for v = 1:numOfView 
         
        for i = 1:num_of_layers
            %Zupdate
            Z{v,i}=Zupdate(data{1,v},i,v,num_of_layers,Z,W,beta,gamma);
        
            %Wupdate
            W{1,i}=Wupdate(data,Z,W,i,v,num_of_layers,alpha,lammad,L_W,gamma);
        end
     end
    
     %得出最终关联矩阵
    [row,col] = size(data{1,1});
    M = zeros(row,col);
    for i = 1:numOfView
        Z_all = Z{i,num_of_layers};
        for j = num_of_layers-1:-1:1
            Z_all = Z_all*Z{i,j};
        end
        M_result{i,1} = W{1,num_of_layers}*Z_all;
        M = M + M_result{i,1};
    end
    M = M./numOfView;
     %% 和multiGMF同样的判定迭代终止条件
     if Dataset == 1
         if T == 1
             stop1_0 = stop1;
             stop1 = norm(M(:,594:end)'-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,594:end)';
             if stop1<tol1 & stop2<tol2
                 break
             end
         elseif T == 0
             stop1_0 = stop1;
             stop1 = norm(M(:,314:end)-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,314:end);
             if stop1<tol1 & stop2<tol2
                 break
             end
         end
     elseif Dataset == 2
         if T == 1
             stop1_0 = stop1;
             stop1 = norm(M(:,664:end)'-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,664:end)';
             if stop1<tol1 & stop2<tol2
                 break
             end
         elseif T == 0
             stop1_0 = stop1;
             stop1 = norm(M(:,410:end)-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,410:end);
             if stop1<tol1 & stop2<tol2
                 break
             end
         end
     elseif Dataset == 3
         if T == 1
             stop1_0 = stop1;
             stop1 = norm(M(:,1238:end)'-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,1238:end)';
             if stop1<tol1 & stop2<tol2
                 break
             end
         elseif T == 0
             stop1_0 = stop1;
             stop1 = norm(M(:,279:end)-WH,'fro')/norm(WH,'fro');
             stop2 = abs(stop1 - stop1_0)/max(1, abs(stop1_0));
             WH = M(:,279:end);
             if stop1<tol1 & stop2<tol2
                 break
             end
         end
     end
end

end