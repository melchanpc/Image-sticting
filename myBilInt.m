function [zp] = myBilInt(zin,Xp)
% myBilInt Bilinear interpolation
%
% Input arguments:
%   zin : the input image - a matrix
%   Xp  : a n*2 matrix of points to be interpolated
% Output argument:
%   zp : the interpolated values

% Get some info from the inputs
N = length(Xp);
sz = size(zin);

% Remove input points that exceed the boundaries of the input image
badX_1 = find(Xp(:,1) > sz(2)-1);       % get the indices
badX_2 = find(Xp(:,1) < 1);
badY_1 = find(Xp(:,2) > sz(1)-1);
badY_2 = find(Xp(:,2) < 1);
badRowIdx = [badX_1; badX_2; badY_1; badY_2];
Xp(badRowIdx,:) =[];               % remove them from Xp

% x1, y1 are the lower bound of the unit square
% x2, y2 are the upper bound of the unit square
x1 = floor(Xp(:,1)); x2 = x1+1;
y1 = floor(Xp(:,2)); y2 = y1+1;

% Get the linear indices of all four points on the unit square for all the
% N points to be interpolated
idx = sub2ind(sz,[y1;y1;y2;y2],[x1;x2;x1;x2]);

% Get the pixel values at all these points
F = double(zin(idx));     

% Convert into N*4 array so that the values are arranged in the order
% f(x_i,y_i), f(x_i+1,y_i), f(x_i,y_i+1) f(x_i+1,y_i+1)
F = reshape(F,[],4);

% Here x and y are the distances from the lower corner of the unit square
x = Xp(:,1) - x1;  y = Xp(:,2) - y1;

% Populate A  matrix A of weights for the points
A = [(1-x).*(1-y) x.*(1-y) (1-x).*y x.*y];

% Calculate points - here temp because have not yet include the points that
% are out of bounds
tempZp = sum(A.*F,2);

% Assign 0 value to the out of bound points 
totIdx = 1:N;
zp = zeros(N,1);
totIdx(badRowIdx) = [];
zp(totIdx) = tempZp;        % finally the output

end

