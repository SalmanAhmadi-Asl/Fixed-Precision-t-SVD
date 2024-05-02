function at = t_trans(a)
    % tubal transpose  of 3d-tensor
    
    if ndims(a) ~= 3
       error("dimension of tensor a should be equal 3")
    end
    
    [n1, n2, n3] = size(a);
    at = zeros([n2, n1, n3]);
    
    at(:,:,1) = a(:,:,1)';
    for i=2:n3
        at(:,:,n3-i+2) = a(:,:,i)';
    end   

end
