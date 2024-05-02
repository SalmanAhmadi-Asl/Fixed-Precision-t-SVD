function [Q_1,R_1]= QR_tubal(X)
[~,~,n3]=size(X); 
X=fft(X,[],3);
for i=1:n3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
end
Q_1=ifft(Q,[],3);
R_1=ifft(R,[],3);