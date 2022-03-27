%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    File: estBlur.m                                         %
%    Author: Jerry Peng                                      %
%    Date: Nov/2014                                          %
%                                                            %
%    Blurriness map estimation                               %
% -----------------------------------------------------------%
%   MODIFICATIONS                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Blur = estBlur(I, win)
I = im2double(I);
[height, width, ~] = size(I);
imYUV = rgb2ycbcr(I);
imY = imYUV(:,:,1);
radius = [9, 17, 33, 65];
DiffImage = zeros(height, width, 4);

for idx = 1:length(radius)
    r = radius(idx);
    sigma = r;
    GFImage = imfilter(imY, fspecial('gaussian', r, sigma), 'replicate');
    DiffImage(:,:,idx) = abs(imY-GFImage);
end
roughBlurMap = mean(DiffImage, 3);

% dilation and hole-filling
Blur = imfill(imdilate(roughBlurMap, strel('square', win)), 8, 'holes');