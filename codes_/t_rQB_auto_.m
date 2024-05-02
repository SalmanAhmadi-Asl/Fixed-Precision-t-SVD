function [Q, B, k] = t_rQB_auto_(A, relerr, b, P)
% rank-revealing tubal QB
    [n1, n2, n3] = size(A);
    
    E = norm(A(:), 'fro')^2 ; E0 = E; % !! All norms are in Signal space.
    threshold = relerr^2 * E;
    maxiter = 1000;
    maxiter = min(maxiter, ceil(n2/b));
    flag = false;
    
    %A = fft(A, [], 3);
    l = 0;
    Q = zeros(n1, l, n3);
    B = zeros(l, n2, n3);
    
    for i=1:maxiter
        Omg = randn(n2, b, n3);
        Y = t_prod(A, Omg) - t_prod(Q, t_prod(B, Omg)); 
        
         Qi = QR_tubal(Y); 
        for j = 1:P
            Qi = QR_tubal(t_prod(t_trans(A), Qi)-t_prod(t_prod(t_trans(B),t_trans(Q)),Qi) );
            Qi = QR_tubal(t_prod(A, Qi)-t_prod(t_prod(Q,B),Qi) );
        end
        Qi = QR_tubal(Qi - t_prod(Q, (t_prod(t_trans(Q), Qi)))); 

        Bi=t_prod(t_trans(Qi),A);

        Q = cat(2, Q, Qi); 
        B = cat(1, B, Bi); 
        
        temp = E - norm(Bi(:), 'fro')^2; 

        if temp < threshold   % precise rank determination 
            for j=1:b
                btmp = Bi(j,:,:);
                E = E - norm(btmp(:), 'fro')^2;
                if E < threshold
                    flag = true;
                    break;
                end
            end
        else
            E = temp;
        end
        if flag
            k = (i - 1) * b + j;
            Q = Q(:,1:k,:);
            B = B(1:k,:,:);
            break;
        end
        
    end

    if ~flag,
        k = i*b;
        fprintf('E = %f. Fail to converge within maxiter!\n\n', sqrt(E/E0));
    end

end
