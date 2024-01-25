function Q = t_fft_orth(X)
    [n1, ~, n3] = size(X);
    Q = zeros(n1, 0, n3);
    kmax = 0;
    for i=1:n3
        [Qi, ~] = qr(X(:,:,i), 0);
        k = size(Qi, 2);
        if k < kmax
            Qi(n1,kmax) = 0;
        end
        if k > kmax
            kmax = k;
            Q(n1,kmax,n3) = 0;
        end
        Q(:,:,i) = Qi;
    end
end