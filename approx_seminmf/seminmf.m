function [Z, H, dnorm] = seminmf( X, k)
% Matrix sizes
% X: m x n
% Z: m x num_of_components
% H: num_of_components x num_of_components

% Process optional arguments
pnames = {'Z0' 'W0' 'bUpdateW' 'maxiter' 'TolFun' 'bUpdateZ' 'save' 'fast'};

W0 = rand(k, size(X, 2)); 
Z0 = rand(size(X, 1),k ); 

dflts  = {Z0, W0, 1, 100,  1e-5, 1, 0, 0};

[Z, H, bUpdateW, max_iter, tolfun, bUpdateZ, doSave, fastapprox] = ...
        internal.stats.parseArgs(pnames,dflts);

if fastapprox
    if k < size(X,1)&& k < size(X,2)
        H = LPinitSemiNMF(X, k);
        %
        Z = LPinitSemiNMF(X, size(X, 1));
    else
        H = rand(k, size(X, 2)); 
        %
        Z = rand(k, size(X, 1)); 
    end
end

key = generate_checksum(X, k);

if ispc
    path = ['\\fs-vol-hci2.doc.ic.ac.uk\hci2\projects\trigeorgis\nmf\seminmf_cache\' key '.mat'];
else
    path = ['/home/zhao.han/HierarchicalMVC/Deep-Semi-NMF-master/seminmf_cache/' key '.mat'];
end

if  doSave && exist(path, 'file') ~= 0
    load(path);
    return;
end

X_H = X*H';
X_H_p = (abs(X_H)+X_H)./2;
X_H_n = (abs(X_H)-X_H)./2;
Z_HH = H*H';
Z_HH_p = (abs(Z_HH)+Z_HH)./2;
Z_HH_n = (abs(Z_HH)-Z_HH)./2;

if length(Z) == 1
%     Z = X * pinv(H);
    %
    Z = Z .* sqrt((X_H_p + Z*Z_HH_n) ./ max(X_H_n + Z*Z_HH_p, eps));
end

dnorm = norm(X - Z * H, 'fro');

for i = 1:max_iter
    
    if bUpdateZ
        try
%             Z = X * pinv(H);
            
            %
            Z = Z .* sqrt((X_H_p + Z*Z_HH_n) ./ max(X_H_n + Z*Z_HH_p, eps));
        catch
            display('Error inverting');
        end
    end
    
    A = Z' * X;
    Ap = (abs(A)+A)./2;
    An = (abs(A)-A)./2;
    
    B = Z' * Z;
    Bp = (abs(B)+B)./2;
    Bn = (abs(B)-B)./2;
    
    if bUpdateW
        H = H .* sqrt((Ap + Bn * H) ./ max(An + Bp * H, eps));
    end
      
    if mod(i, 10) == 0 || mod(i+1, 10) == 0 
        
        s = X - Z * H;
        dnorm = sqrt(sum(s(:).^2));
        % dnorm = norm(gX - Z * H, 'fro');
        
        if mod(i+1, 10) == 0
            dnorm0 = dnorm;
            continue
        end

        if 0 && exist('dnorm0')
            assert(dnorm <= dnorm0, sprintf('Rec. error increasing! From %f to %f. (%d)', dnorm0, dnorm, k));
        end

        % Check for convergence
        if exist('dnorm0') && dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        end
     
    end
end

if doSave
dnorm = norm(X - Z * H, 'fro');
save(path, 'Z', 'H', 'dnorm');
end
