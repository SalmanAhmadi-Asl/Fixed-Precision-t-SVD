function [U, S, V, ranks] = t_SVDbasic_auto(X, eps, b)
    
    if ndims(X) ~= 3
        error("X should be 3rd-order tensor")
    end
    
    [n1, n2, n3] = size(X);
    xf = fft(X, [], 3);
    
    ranks = zeros(1,n3);
    neps = eps / n3^.5;
    uf = zeros(0,0,n3);
    vf = zeros(0,0,n3);
    sf = zeros(0,0,n3);
    
    for i=1:n3
        
        b = i;
        [u1, s1, v1] = svd2d(xf(:,:,i), neps, b);
        k = size(s1, 1);
        ranks(i) = k;
        
        if i == 1
            maxk = k;
        end
        
        if k < maxk
            u1(n1, maxk) = 0;
            v1(n2, maxk) = 0;
            s1(maxk, maxk) = 0;
        end
        
        if k > maxk
            maxk = k;
            uf(n1, maxk, i) = 0;
            vf(n2, maxk, i) = 0;
            sf(maxk, maxk, i) = 0;
        end
        
        uf(:,:,i) = u1;
        vf(:,:,i) = v1;
        sf(:,:,i) = s1;
    end
    
    U = ifft(uf, [], 3);
    V = ifft(vf, [], 3);
    S = ifft(sf, [], 3);
    
end

function [u, s, v] = svd2d(x, eps, b)
    [u, s, v] = svd(x);
end