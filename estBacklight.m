%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    File: estBacklight.m                                    %
%    Author: Jerry Peng                                      %
%    Date: Aug/2015                                          %
%                                                            %
%    background light estimation                             %
% -----------------------------------------------------------%
%   MODIFICATIONS                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B_pro = estBacklight(im, blur)
I = im;
I_s = I;
as = 32;
thr_sig = 0.2;
[h, w, ~] = size(I);
s = h*w;
thr_s = floor(s/1024);

grI = rgb2gray(I);
fin_lef = zeros(2, 1);
fin_rig = zeros(2, 1);
fin_top = zeros(2, 1);
fin_bot = zeros(2, 1);
var_blur= zeros(4, 1);
vb = zeros(2, 1);
bv = zeros(2, 1);

for f = 1:2
lef =1;
rig =w;
top =1;
bot =h;
min_res = 1e8;
ss = 1e8;
sub_lef = zeros(4, 1);
sub_rig = zeros(4, 1);
sub_top = zeros(4, 1);
sub_bot = zeros(4, 1);
while ss > thr_s
%tl part
sub_lef(1) = lef;
sub_rig(1) = lef+floor((rig-lef+1)/2);
sub_top(1) = top;
sub_bot(1) = top+floor((bot-top+1)/2);
%tr part
sub_lef(2) = sub_rig(1)+1;
sub_rig(2) = rig;
sub_top(2) = top;
sub_bot(2) = top+floor((bot-top+1)/2);
%bl part
sub_lef(3) = lef;
sub_rig(3) = lef+floor((rig-lef+1)/2);
sub_top(3) = sub_bot(1)+1;
sub_bot(3) = bot;
%br part
sub_lef(4) = sub_rig(3)+1;
sub_rig(4) = rig;
sub_top(4) = sub_bot(2)+1;
sub_bot(4) = bot;

ls = 3;
if f == 1
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 1) = 0;
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 2) = 0;
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 3) = 255;

I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 1) = 0;
I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 2) = 0;
I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 3) = 255;
else
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 1) = 255;
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 2) = 0;
I_s(sub_bot(1):sub_bot(1)+ls, sub_lef(1):sub_rig(2), 3) = 0;

I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 1) = 255;
I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 2) = 0;
I_s(sub_top(1):sub_bot(3), sub_rig(1):sub_rig(1)+ls, 3) = 0;
end    
res = zeros(4, 1);
for ind = 1:4
    if f == 1 %var
        res(ind) = var(double(reshape(grI(sub_top(ind):sub_bot(ind), sub_lef(ind):sub_rig(ind)), [], 1)), 1);
        var_blur(ind) = mean(double(reshape(blur(sub_top(ind):sub_bot(ind), sub_lef(ind):sub_rig(ind)), [], 1)), 1);
    else %blur
        res(ind) = mean(double(reshape(blur(sub_top(ind):sub_bot(ind), sub_lef(ind):sub_rig(ind)), [], 1)), 1);
        var_blur(ind) = var(double(reshape(grI(sub_top(ind):sub_bot(ind), sub_lef(ind):sub_rig(ind)), [], 1)), 1);
    end
end

min_res = min(res(:));
for ind = 1:4
    if res(ind) == min_res
        lef = sub_lef(ind);
        rig = sub_rig(ind);
        top = sub_top(ind);
        bot = sub_bot(ind);
        if f == 1 % var
            vb(1) = min_res;
            vb(2) = var_blur(ind);
        else
            bv(1) = min_res;
            bv(2) = var_blur(ind); 
        end
        break;
    end
end
ww = abs(rig-lef);
hh = abs(bot-top);
ss = ww*hh;
end

fin_lef(f) = lef;
fin_rig(f) = rig;
fin_top(f) = top;
fin_bot(f) = bot;
for i=1:3
    B_cand(f, i) = mean(reshape(I(top:bot, lef:rig, i), [], 1));
end

end % f = 1 (var) / 2 (blur)

for idx = 1:2
    I_s(fin_top(idx):fin_bot(idx), fin_lef(idx):fin_rig(idx), :) = 255;
end

imsize = w * h;
JDark = 1-blur;
numpx = floor(imsize/1000); % accomodate for small images
JDarkVec = reshape(JDark,imsize,1); % a vector of pixels in JDark
ImVec = reshape(I,imsize,3);  % a vector of pixels in my image

[~, indices] = sort(JDarkVec); %sort
indices = indices(imsize-numpx+1:end); % need the last few pixels because those are closest to 1

atmSum = zeros(1,3);
for ind = 1:numpx
    atmSum = atmSum + ImVec(indices(ind),:);
end
A = atmSum / numpx;
B_pro = zeros(1, 3);

for ind = 1:3
    theta = sigmf(numel(find(I(:,:,ind) > 0.5))/s,[as thr_sig]);
    B_min = min([B_cand(1, ind), B_cand(2, ind), A(ind)]);
    B_max = max([B_cand(1, ind), B_cand(2, ind), A(ind)]);
    B_pro(ind) = B_min*(1-theta)+B_max*theta;
end