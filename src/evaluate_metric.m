function d = evaluate_metric(It, Io)
%EVALUATE_METRIC Compute the gamut distance between the target image
%   and the output image.
%
%   Input
%   -------
%   It: the target image
%   Io: the output image
%
%   Output 
%   -------
%   d: the gamut distance between the target and output images
%
%   Reference 
%   -------
%   Rang M.H. Nguyen, Seon Joo Kim, Michael S. Brown
%   Illuminant Aware Gamut-Based Color Transfer
%   Computer Graphics Forum 2014
%
%   Date
%   -------
%   Nov. 24, 2014

Io = reshape(Io, [], 3);
It = reshape(It, [], 3);

[cho, Vo] = convhull(Io, 'simplify', true);
[cht, Vt] = convhull(It, 'simplify', true);

ido = unique(cho(:));
idt = unique(cht(:));

[~, Vtotal] = convhull([Io(ido,:); It(idt,:)]);

d = (Vtotal - Vo) + (Vtotal - Vt);

