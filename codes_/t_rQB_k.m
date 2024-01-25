function [Q, B] = t_rQB_k(A, k, b, P, opt, flag)
% in every slices it used randQB_EI_k from fSVT/randQB
% k is fixed rank parameter 
% P is the power parameter, b is rank-increase step (usually a factor of k).
% opt is a switch for output errors (can be any value).
% errs is the error ||A-QB||_F, errQ is ||I-Q'*Q||_inf.
% err_id is the error indicator, ||A||^2-||B||^2.

[n1, n2, n3] = size(A);

Q = zeros(n1, k, n3);
B = zeros(k, n2, n3);


if flag 
    A = fft(A, [], 3);
    
    for i=1:n3
        [Q1, B1, ~, ~, ~] = randQB_EI_k(A(:,:,i), k+1, b, P, opt);
        Q(:,:,i) = Q1;
        B(:,:,i) = B1;
        
    end
    
    Q = ifft(Q, [], 3);
    B = ifft(B, [], 3);
else
    for i=1:n3
        A1 = fft(A(:,:,i));
        [Q1, B1, ~, ~, ~] = randQB_EI_k(A1, k+1, b, P, opt);
        Q(:,:,i) = ifft(Q1);
        B(:,:,i) = ifft(B1);
    end
end

end 