function [Z,W] = fdrinit(X, layers)

num_of_layers = numel(layers);
numOfView = numel(X);
numOfSample = size(X{1,1},2);
Z = cell(numOfView, num_of_layers);
W = cell(1, num_of_layers);

for v_ind = 1:numOfView
    
    fac=X{1,v_ind};
    for j = 1:length(layers)
        if v_ind==1
%             [Z{v_ind,j},W{1,j}, ~] = seminmf(fac,layers(j));
            
            [W{1,j},Z{v_ind,j}, ~] = seminmf(fac,layers(j));
            fac=W{1,j};
%             W_i{v_ind,j} = W{1,j};
        else
%             [Z{v_ind,j},Wmint1, ~] = seminmf(fac,layers(j));
            
            [Wmint1,Z{v_ind,j}, ~] = seminmf(fac,layers(j));
%             [Z{v_ind,j},Wmint1]=SVDinitSemiNMF(fac,layers(j),1);
            W{1,j}=W{1,j}+Wmint1;
            fac=Wmint1;
%             W_i{v_ind,j} = Wmint1;
        end
    end
end

for X_v = 1:num_of_layers
    W{1,X_v}=W{1,X_v}/numOfView;
end


end