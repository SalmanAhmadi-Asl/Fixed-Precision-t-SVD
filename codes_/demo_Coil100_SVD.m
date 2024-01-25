clc;clear all
randn('state',1); rand('state',1); 

load('coil100.mat');
A1=double(X);
A1=A1(:,:,:,1:30,1:30);
%ad = reshape(A1,[128,128,2700]) / 255;
ad=reshape(A1,[384,384,300]) / 255;


relerr = 5e-1;
b = 50;
P =1;
tic
[U, S, V] = t_rSVD_auto(ad, relerr, b, P);
toc

t_prod(t_prod(U, S), t_trans(V)) * 255;

tic
[U,S,V] = tensor_t_svd(ad,201);
toc

Anew = t_prod(t_prod(U, S), t_trans(V)) * 255;
Anew = reshape(Anew,[128,128,3,30,30]);

n_obj = 11; n_angle = 10;
subplot(1,2,1)
imshow(uint8(A1(:,:,:,n_obj,n_angle)))
subplot(1,2,2)
imshow(uint8(Anew(:,:,:,n_obj,n_angle)))

fprintf("Compression rate: %f\n", (numel(A1)/(numel(U)+numel(V) + numel(S)) ))