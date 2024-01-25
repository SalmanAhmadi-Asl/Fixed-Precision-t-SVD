function [Q_1,R_1]= fastQR_tubal(X)
[~,~,n3]=size(X); 
X=fft(X,[],3);
[n1,n2,n3] = size(X);
halfn3 = ceil((n3+1)/2);
for i=1:halfn3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
end
for i=halfn3+1:n3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,n3+2-i),0);
end
Q_1=ifft(Q,[],3);
R_1=ifft(R,[],3);