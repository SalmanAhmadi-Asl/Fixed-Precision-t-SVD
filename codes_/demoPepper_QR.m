clc;clear all
randn('state',1); rand('state',1); 


A1=imread('peppers.png');

% A1(:,:,2) = zeros(256)*255;
% A1(:,:,1) = zeros(256)*255;

%A1=reshape(A1,[64,48,64]);
ad = double(A1) / 255;

relerr = 2e-1;
b = 4;
P = 2;
[Q, B, k] = t_rQB_auto_(ad, relerr, b, P);
disp(k);

fprintf("Compression rate: %f\n", (numel(A1)/(numel(Q)+numel(B)) ))
Anew  = t_prod(Q, B) * 255;
Anew=reshape(Anew,[256,256,3]);
A1=reshape(A1,[256,256,3]);

subplot(1,2,1)
imshow(A1)
subplot(1,2,2)
imshow(uint8(Anew))
