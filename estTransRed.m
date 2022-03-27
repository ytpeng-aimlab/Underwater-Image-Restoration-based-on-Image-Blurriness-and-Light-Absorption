%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    File: estTransRed.m                                     %
%    Author: Jerry Peng                                      %
%    Date: Aug/2015                                          %
%                                                            %
%    Transmission estimation based on                        %
%    attenuation of red channel                              %
%------------------------------------------------------------%
%    MODIFICATIONS                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trans_mip, trans_red] = estTransRed(I, win)
Strch = @(x) (x-min(x(:))).*(1/(max(x(:))-min(x(:))));

max_r_im = ordfilt2(I(:,:,1),win^2,ones(win,win),zeros(win,win),'symmetric');
max_g_im = ordfilt2(I(:,:,2),win^2,ones(win,win),zeros(win,win),'symmetric');
max_b_im = ordfilt2(I(:,:,3),win^2,ones(win,win),zeros(win,win),'symmetric');

gb_im = cat(3,max_g_im,max_b_im);
max_gb_im = max(gb_im,[],3);
diff = max_r_im - max_gb_im;
trans_mip= Strch(diff);
trans_red= Strch(max_r_im);