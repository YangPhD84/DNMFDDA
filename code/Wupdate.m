function Wi=Wupdate(data,Z,W,i,v,m,alpha,lammad,L_W,gamma)
%data表示所有视图数据，Z表示所有深度矩阵分解项，Wi表示第i层矩阵，i表示第i层，v表示第v个视图，m表示总层数

Wi = W{1,i};
Isize=size(Z{v,i},1);
phi = eye(Isize);
X_phi = 0;
W_phi_phi = 0;
W_Z = 0;
W_Z_Z = 0;
W_Z_t = 0;
for v_ind = 1:v
    for j_pi = 1:i
        phi = phi*Z{v_ind,i-j_pi+1};
    end
    phi_tran = phi';
    
    p1 = data{1,v_ind}*phi_tran;
    X_phi = X_phi + p1;
    p2 = Wi*phi*phi_tran;
    W_phi_phi = W_phi_phi + p2;
    phi = eye(Isize);
end
X_phi_p = (abs(X_phi)+X_phi)/2;
X_phi_n = (abs(X_phi)-X_phi)/2;
W_phi_phi_p = (abs(W_phi_phi)+W_phi_phi)/2;
W_phi_phi_n = (abs(W_phi_phi)-W_phi_phi)/2;

% %2,1
% D = eye(size(Wi,1));
% for D_i = 1: size(Wi,1)
%     Wi_norm = norm(Wi(D_i,:));
%     D(D_i,D_i) = 1/Wi_norm;
% end
% D_Wi = D*Wi;
% D_Wi_p = (abs(D_Wi)+D_Wi)/2;
% D_Wi_n = (abs(D_Wi)-D_Wi)/2;

%F
W_p = (abs(Wi)+Wi)/2;
W_n = (abs(Wi)-Wi)/2;

%GrLap
GrL_W = L_W{1,v} * Wi;
L_W_p = (abs(GrL_W)+GrL_W)/2;
L_W_n = (abs(GrL_W)-GrL_W)/2;

%%
if i == 1 & m == 1
    
    Wi = Wi.*sqrt((X_phi_p+W_phi_phi_n+alpha*L_W_n+lammad*W_n)./max((X_phi_n+W_phi_phi_p+alpha*L_W_p+lammad*W_p),1e-10));
    
elseif i == 1 & m ~= 1
    %
    for v_ind = 1:v
        p3 = W{1,i+1}*Z{v_ind,i+1};
        W_Z = W_Z + p3;
    end
    W_Z_p = (abs(W_Z)+W_Z)/2;
    W_Z_n = (abs(W_Z)-W_Z)/2;
    
    Wi = Wi.*sqrt((X_phi_p+W_phi_phi_n+alpha*L_W_n+(lammad+gamma)*W_n+gamma*W_Z_p)./max((X_phi_n+W_phi_phi_p+alpha*L_W_p+(lammad+gamma)*W_p+gamma*W_Z_n),1e-10));

elseif i == m
    %
    for v_ind = 1:v
        p4 = W{1,i}*Z{v_ind,i}*Z{v_ind,i}';
        W_Z_Z = W_Z_Z + p4;
        
        p5 = W{1,i-1}*Z{v_ind,i}';
        W_Z_t = W_Z_t + p5;
    end
    W_Z_Z_p = (abs(W_Z_Z)+W_Z_Z)/2;
    W_Z_Z_n = (abs(W_Z_Z)-W_Z_Z)/2;
    
    W_Z_t_p = (abs(W_Z_t)+W_Z_t)/2;
    W_Z_t_n = (abs(W_Z_t)-W_Z_t)/2;
    
     Wi = Wi.*sqrt((X_phi_p+W_phi_phi_n+alpha*L_W_n+lammad*W_n+gamma*W_Z_t_p+gamma*W_Z_Z_n)./max((X_phi_n+W_phi_phi_p+alpha*L_W_p+lammad*W_p+gamma*W_Z_t_n+gamma*W_Z_Z_p),1e-10));

elseif i ~= m
    %
     for v_ind = 1:v
        p3 = W{1,i+1}*Z{v_ind,i+1};
        W_Z = W_Z + p3;
        
        p4 = W{1,i}*Z{v_ind,i}*Z{v_ind,i}';
        W_Z_Z = W_Z_Z + p4;
        
        p5 = W{1,i-1}*Z{v_ind,i}';
        W_Z_t = W_Z_t + p5;
    end
    W_Z_p = (abs(W_Z)+W_Z)/2;
    W_Z_n = (abs(W_Z)-W_Z)/2;
    
    W_Z_Z_p = (abs(W_Z_Z)+W_Z_Z)/2;
    W_Z_Z_n = (abs(W_Z_Z)-W_Z_Z)/2;
    
    W_Z_t_p = (abs(W_Z_t)+W_Z_t)/2;
    W_Z_t_n = (abs(W_Z_t)-W_Z_t)/2;
    
    Wi = Wi.*sqrt((X_phi_p+W_phi_phi_n+alpha*L_W_n+(lammad+gamma)*W_n+gamma*W_Z_t_p+gamma*W_Z_Z_n+gamma*W_Z_p)./max((X_phi_n+W_phi_phi_p+alpha*L_W_p+(lammad+gamma)*W_p+gamma*W_Z_t_n+gamma*W_Z_Z_p+gamma*W_Z_n),1e-10));

end

end