function [U,S,V] = t_rSVD_auto(A, relerr, b, P)
    [n1, n2, n3] = size(A);
    
    [Q, B] = t_rQB_auto_(A, relerr, b, P);
    %Q = conj(permute(fft(Q, [], 3), [2,1,3]));
    B = fft(B, [], 3);
    
    kmax = size(B, 1);
    
    U = zeros(kmax, kmax, n3);
    S = zeros(kmax, kmax, n3);
    V = zeros(n2, kmax, n3);
    
    for i=1:n3
       [u1,s1,v1] = svd(B(:,:,i), 'econ');
       k = size(u1, 2);
        if k < kmax
            u1(kmax,kmax) = 0;
            s1(kmax,kmax) = 0;
            v1(n2,kmax) = 0;
        end
        U(:,:,i) = u1;
        S(:,:,i) = s1;
        V(:,:,i) = v1;
    end
    
    U = ifft(U, [], 3);
    U = t_prod(Q, U);
    S = ifft(S, [], 3);
    V = ifft(V, [], 3);
end