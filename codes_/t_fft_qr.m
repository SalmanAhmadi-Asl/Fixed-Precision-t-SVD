function [Q, R] = t_fft_qr(X)
    [n1, n2, n3] = size(X);
    Q = zeros(n1, 0, n3);
    R = zeros(0, n2, n3);
    kmax = 0;
    for i=1:n3
        [Qi, Ri] = qr(X(:,:,i), 0);
        k = size(Qi, 2);
        if k < kmax
            Qi(n1,kmax) = 0;
            Ri(kmax, n2) = 0;
        end
        if k > kmax
            kmax = k;
            Q(n1,kmax,n3) = 0;
            R(kmax,n2,n3) = 0;
        end
        Q(:,:,i) = Qi;
        R(:,:,i) = Ri;
    end
end