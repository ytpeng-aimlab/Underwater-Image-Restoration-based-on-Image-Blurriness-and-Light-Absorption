%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    File: UIR_IBLA.m                                        %
%    Author: Jerry Peng                                      %
%    Date: Aug/2015                                          %
%                                                            %
%    Underwater Image Restoration based on Image             %
%    Blurriness and Light Absorption                         %
%------------------------------------------------------------%
%    MODIFICATIONS                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

beta = [1/7; 1/25 ; 1/28];
t1 = 1;
t0 = 0.1;
D = 8;
B_thr_high= 0.5;
r_thr_high= 0.1;
as = 32;
win = 7;
use_lap = 1;
wl = [620 540 450];
Strch = @(x) (x-min(x(:))).*(1/(max(x(:))-min(x(:))));
%%
[fn,pn,fi]=uigetfile('*.bmp;*.jpg;*.png','open files');
I = im2double(imread([pn fn]));
[height, width, ~] = size(I);

if width*height < 180000
    lambda = 10e-6;
else
    lambda = 10e-3;
end

[trans_mip, trans_red] = estTransRed(I, win);

%% Blurriness
t_b_est = estBlur(I, win);

% if use_lap == 0
%     t_blur_est = imguidedfilter(t_b_est, I, 'NeighborhoodSize',[7 7]);
% else
    t_blur_est = laplacian_matting( t_b_est, I, lambda);
% end

B = estBacklight(I, t_blur_est);
mean_b = mean(B);
mean_r = mean(reshape(I(:,:,1), height*width, 1));

alpha = sigmf(mean_b,[as B_thr_high]);
alpha2= sigmf(mean_r,[as r_thr_high]);

t_pro_a_est = trans_mip*alpha+trans_red*(1-alpha);
t_strch_blur_est = Strch(t_b_est);

if use_lap == 0
    t_pro_est = imguidedfilter(t_pro_a_est*alpha2+t_strch_blur_est*(1-alpha2), I, 'NeighborhoodSize',[win win]);
else
    t_pro_est = laplacian_matting(t_pro_a_est*alpha2+t_strch_blur_est*(1-alpha2), I, lambda);
end

%% estimate Distance
BLmap = zeros(size(I));
for ind = 1:3
    BLmap(:,:,ind) = B(ind) * ones(size(height, width));
end
diff_BL_I= abs(BLmap-I);
[max_diff_bl_i, pos] = max(diff_BL_I(:));
R = max_diff_bl_i/max(BLmap(pos), 1-BLmap(pos));

%% Trans

BL = max(B, 0.1);
b_sf  = -0.00113*wl+1.62517;
cg2cr = (b_sf(2)*BL(1))/(b_sf(1)*BL(2));
cb2cr = (b_sf(3)*BL(1))/(b_sf(1)*BL(3));

trans = zeros(size(I));

beta(2) = beta(1)*cg2cr;
beta(3) = beta(1)*cb2cr;

dist = D*(1-t_pro_est+(1-R));
for ind = 1:3
    trans(:,:,ind) = exp(-beta(ind).*dist);
end

J_pro = zeros(size(I));
for ind = 1:3
    J_pro(:,:,ind) =B(ind)+(I(:,:,ind)-B(ind))./max(trans(:,:,ind), t0);
end
J_pro(J_pro < 0) = 0;
J_pro(J_pro > 1) = 1;

figure, imshow([I J_pro], 'Border','tight');
