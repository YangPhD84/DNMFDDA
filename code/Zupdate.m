function Zi=Zupdate(datav,i,v,m,Z,W,beta,gamma)

Zi = Z{v,i};
if i>1
    Rsize=size(Z{v,i-1},1);
    phi=eye(Rsize);
    for j_pi = 1:i-1
        phi = phi*Z{v,i-j_pi};
    end

    if any(isnan(phi), 'all') || any(isinf(phi), 'all')
        return;
    end
end
%
Z_p = (abs(Z{v,i})+Z{v,i})/2;
Z_n = (abs(Z{v,i})-Z{v,i})/2;

%
W_W_Z = W{1,i}'*W{1,i}*Z{v,i};
W_W_Z_p = (abs(W_W_Z)+W_W_Z)/2;
W_W_Z_n = (abs(W_W_Z)-W_W_Z)/2;

%%
if i == 1
    
    W_X = W{1,i}'*datav;
    W_X_p = (abs(W_X)+W_X)/2;
    W_X_n = (abs(W_X)-W_X)/2;
    
    Zi = Zi.*sqrt((W_W_Z_n+W_X_p+beta*Z_n)./max((W_W_Z_p+W_X_n+beta*Z_p),1e-10));
else
    %      
    W_Z_phi = W{1,i}'*W{1,i}*Z{v,i}*phi*phi';
    W_Z_phi_p = (abs(W_Z_phi)+W_Z_phi)/2;
    W_Z_phi_n = (abs(W_Z_phi)-W_Z_phi)/2;
    
    W_X_phi = W{1,i}'*datav*phi';
    W_X_phi_p = (abs(W_X_phi)+W_X_phi)/2;
    W_X_phi_n = (abs(W_X_phi)-W_X_phi)/2;
    
    W_W = W{1,i}'*W{1,i-1};
    W_W_p = (abs(W_W)+W_W)/2;
    W_W_n = (abs(W_W)-W_W)/2;
    
    Zi = Zi.*sqrt((W_Z_phi_n+W_X_phi_p+beta*Z_n+gamma*W_W_p+gamma*W_W_Z_n)./max((W_Z_phi_p+W_X_phi_n+beta*Z_p+gamma*W_W_n+gamma*W_W_Z_p),1e-10));
    
end
end
