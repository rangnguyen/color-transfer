function Io = color_transfer(Is, It)
%COLOR_TRANSFER Transfer the color style of the target image into
%   the source image. As a result, these both images share the same
%   'look and feel'.
%
%   Input
%   -------
%   Is: the source image
%   It: the target image
%
%   Output 
%   -------
%   Io: the output image
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

addpath('weighted_GE');
global pi pt Vt
[H, W, ~] = size(Is);
[Ht,Wt,~] = size(It);

%% STEP 1: White-balancing and rotating 
% Grey-Egde algorithm to estimate illuminations of the source and target
mink_norm = 5;
sigma = 2;
kappa = 10; 
[wRs, wGs, wBs, ~] = weightedGE(Is, kappa, mink_norm, sigma);
WBs = [wRs wGs wBs];
[wRt, wGt, wBt, ~] = weightedGE(It, kappa, mink_norm, sigma);
WBt = [wRt wGt wBt];
WBs = sqrt(3)*WBs/sqrt(sum(WBs.^2));
WBt = sqrt(3)*WBt/sqrt(sum(WBt.^2));
Is = reshape(Is, [], 3)';
It = reshape(It, [], 3)';
Is = diag(1./WBs) * Is;
It = diag(1./WBt) * It;

% Rotate 2 images to B axis
Is = rotate_to_Zaxis(Is, [1 1 1]);
It = rotate_to_Zaxis(It, [1 1 1]);

%% STEP 2: Matching luminance
Is = Is';
It = It';
Is(:,3) = normalizeIntensity(Is(:,3), It(:,3), H, W);

%% STEP 3: Aligning gamuts
Ms = mean(Is, 1);
Mt = mean(It, 1);
% Shift the means to the origin
Is = Is - repmat(Ms, H*W, 1);
It = It - repmat(Mt, Ht*Wt, 1);
% Compute the convex-hull
[CHi, ~] = convhull(Is, 'simplify', true);
[CHt, Vt] = convhull(It, 'simplify', true);
idi = unique(CHi(:));
idt = unique(CHt(:));
pi = Is(idi,:);
pt = It(idt,:);
% Compute the optimal matrix to align two gamuts
x0 = [0 1 1];
options = optimset('MaxIter' ,50, 'Display', 'iter', 'DiffMinChange', 0.01);
x = fminunc(@myfunoptimal,x0, options);
T = [x(2)*cos(x(1)) -x(2)*sin(x(1)) 0;
     x(3)*sin(x(1))  x(3)*cos(x(1)) 0;
     0               0              1];
disp('The optimal matrix to align two gamuts:');
disp(T);
% Align two gamuts
Io = T*Is';
Mt(3) = Ms(3);
Io = Io + repmat(Mt',1, H*W);

%% STEP 4: Rotate back and undo white-balancing
Io = rotate_back_from_Zaxis(Io, [1 1 1]);

Io = diag(WBt) * Io;

Io(Io < 0) = 0;
Io(Io > 1) = 1;
Io = Io';
Io = reshape(Io, H, W, 3);


function Io = normalizeIntensity(Is, It, ws, hs)

% create matrix
lambda = 1;
DX = createDXForward(ws,hs);
DX = DX'*DX;

DY = createDYForward(ws,hs);
DY = DY'*DY;

D = lambda*(DX + DY);
A = speye(ws*hs) + D;
clear DX DY

t = sqrt(3);
Is = Is / t;
It = It / t;
If = histeq(Is, imhist(It));


% solve the sparse matrix
b = If + D*Is;
Io = A \ b;
Io = t * Io;


function DX = createDXForward(N, M)
K = N*M;
DX = spdiags([-ones(K,1) -2*ones(K,1) -ones(K,1) ones(K,1) 2*ones(K,1) ones(K,1)],...
              [-N -N+1 -N+2 N N+1 N+2],K,K);
DX(1:N,:) = 0;
DX(K-N+1:end,:) = 0;



function DY = createDYForward(N,M)
K = N*M;
DY = spdiags([-ones(K,1) ones(K,1) -2*ones(K,1) 2*ones(K,1) -ones(K,1) ones(K,1)],...
              [-N -N+2 -1 1 N N+2],K,K);
DY(1:N,:) = 0;
DY(K-N+1:end,:) = 0;

function f = myfunoptimal(x)
global pi pt Vt

T = [x(2)*cos(x(1)) -x(2)*sin(x(1)) 0;
     x(3)*sin(x(1))  x(3)*cos(x(1)) 0;
     0              0               1];
    
Io = pi * T';

[~, Vo] = convhull(Io);
[~, Vtotal] = convhull([Io; pt]);
f = (Vtotal - Vt) + (Vtotal - Vo);


function I = rotate_to_Zaxis(I, a)
b = sqrt(a(1)^2 + a(2)^2);

Txz = [a(1)/b    a(2)/b   0;
      -a(2)/b    a(1)/b   0;
       0         0        1];

c = sqrt(a(1)^2+a(2)^2+a(3)^2);
Tz =  [a(3)/c    0        -b/c;
      0          1         0;
      b/c        0        a(3)/c];

T = Tz*Txz;
I = T*I;

function I = rotate_back_from_Zaxis(I, a)
b = sqrt(a(1)^2 + a(2)^2);

iTxz = [a(1)/b    -a(2)/b   0;
      a(2)/b    a(1)/b   0;
       0         0        1];

c = sqrt(a(1)^2+a(2)^2+a(3)^2);
iTz =  [a(3)/c    0        b/c;
      0          1         0;
      -b/c        0        a(3)/c];

iT = iTxz*iTz;
I = iT * I;

