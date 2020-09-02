%% ECE4076 - Lab 3: Image stitching by homography
clear; clc; close all;

USE_ALPHA = 1;
USE_VAR_SEAM = 1;

%% Task 1: Draw test points

% Read the images
imL = imread('left.jpg');
imR = imread('right.jpg');

% Get the size of the image
szIm = size(imL);

% Marked points for left image -- given
ptsL = [338,197,1; 468,290,1; 253,170,1; 263,256,1; 242,136,1];

% Mark the points on left image and show it
imL_marked = insertMarker(imL,ptsL(:,1:2),'x','color','r');
figure(1); subplot(1,2,1); imshow(imL_marked); title('Left image');

%% Task 2: Use Homography to find right image points
% The given homography matrix
H = [1.6010 -0.0300 -317.9341; 0.1279 1.5325 -22.5847; 0.0007 0 1.2865];

% Calculates the corresponding points on the right image and rescale
ptsR = (H*ptsL')';
ptsR = ptsR./ptsR(:,3);

% Mark the points on the right image and show it
imR_marked = insertMarker(imR,ptsR(:,1:2),'x','color','r');
figure(1); subplot(1,2,2);imshow(imR_marked); title('Right image');
set(gcf, 'Position', get(0, 'Screensize'));

%% Normalize both images
% Normalize the images by obtaining the average value of windows of a fixed
% size centred around the given test points.
% This average is compared between the images and a scaling factor for each
% image is derived

szWin = 30;          % window size
L = length(ptsL);    % number of test points
avgWin = zeros(L,2); % mean value of each window

% Repeat for all points and for both image
for i = 1:L
    for j = 1:2
        if j==1
            pts = ptsL; im = imL;
        else
            pts = round(ptsR); im = imR;
        end
        % some pre-processing to get the indices of the image
        x = pts(i,1)-szWin:pts(i,1)+szWin;
        y = pts(i,2)-szWin:pts(i,2)+szWin;
        [X, Y] = meshgrid(x,y);
        idx = sub2ind(szIm,Y,X);
                
        win = im(idx);                  % the window
        avgWin(i,j) = mean(mean(win));  % calculate the mean
    end
end

% For each point, calculate the ratio of each image's mean to the average
% of both images' mean
k = avgWin./mean(avgWin,2);
k = mean(k);        % then take the mean of this for all the points


% Apply the scale
imL = k(2)*imL;
imR = k(1)*imR;


% Display the scaled images
figure(2); subplot(1,2,1); imshow(imL); title('Left image (scaled)');
subplot(1,2,2); imshow(imR); title('Right image (scaled)');
set(gcf, 'Position', get(0, 'Screensize'));

%% Task 3: Bilinear interpolation of the right image
% Calculate the interpolated values for the given points
int = myBilInt(imR,ptsR(:,1:2));

% Displays the result to the command window
disp('The interpolated intensity value for the transformed points are:');
for i = 1:length(int)
    fprintf('  %.2f  ', int(i));
end
fprintf('\n');


%% Rough estimate of seam column
% Aim is to get the first column from the left image whose corresponding
% point exists on the right image (computed using the homography matrix)

% Get the set of points for each column of the first row from the left
% image
setPts = [(1:szIm(2))' ones(szIm(2),2)];

% Calculates the corresponding points on the right image
ptsR = (H*setPts')'; 
ptsR = ptsR./ptsR(:,3);

% The estimated seam point is the midpoint between the first overlap column
% and the end (width) of the image
firstOverlap = find(ptsR(:,1)>1,1);  
seamPtEst = round((firstOverlap + 1.15*szIm(2))/2);

%% Task 4: Image stitching

% some processing to get all the points to be stitched
% aim is to list out the points in an array of size N*3,
% where N = szIm(1)*szIm(2), and the 3rd column is 1 for use of homography
% matrix to get the corresponding right image points
ptsSt_x = seamPtEst+1:2*szIm(2);      % x-coordinates of the stitched image
ptsSt_y = 1:szIm(1);                % y-coordinates are same as original
[ptsSt_X, ptsSt_Y] = meshgrid(ptsSt_x, ptsSt_y);
ptsSt_X = reshape(ptsSt_X,[],1);
ptsSt_Y = reshape(ptsSt_Y,[],1);
ptsSt = [ptsSt_X ptsSt_Y ones(length(ptsSt_X),1)];
% clears the temporary variables to prevent memory wastage 
clear ptsSt_x ptsSt_y ptsSt_X ptsSt_Y;  
 
% Evaluate the equivalent right image points using the homography matrix
ptsR = (H*ptsSt')'; 
ptsR = ptsR./ptsR(:,3);

% Calls the billinear function to get the RHS of the stitched image
imStR = myBilInt(imR,ptsR(:,1:2));  % result is in vector form
imStR = reshape(imStR,szIm(1),[]);  % reshape it to get back to 2D
imStR = uint8(imStR); 

%% Continue with seam point evaluation
% For each row of the image, get the column which gives the minimum
% difference between the corresponding points on the left and right images
% This column would be the seam point for that row 
% (each row has its own seam column)

if (USE_VAR_SEAM)
    % Get the two matrices for comparison
    compR = imStR(:,1:szIm(2)-seamPtEst);
    compL = imL(:,seamPtEst+1:end);

    % Calculates absolute difference and get the corresponding smallest indices
    err = double(abs(compL-compR));
    [minErr, seamCol] = min(err,[],2);
else
    seamCol = ones(szIm(1),1);
end

% Now stitch the image

if (USE_ALPHA)
    winWidth = 15;          % window width of alpha blending
    imStitch = uint8(zeros(szIm(1),2*szIm(2)-winWidth));  % initialize the stitched image

    winHeight = 11;          % window height of alpha blending

    weightsRight = linspace(0,1,winWidth);  % weights for blending
    weightsRight = sin(weightsRight*pi/2);
    weightsLeft= 1-weightsRight;

    % Repeat for all windows
    for i = ceil(winHeight/2):winHeight:szIm(1)
        idxL = i - floor(winHeight/2);      % lower row index
        idxU = i + floor(winHeight/2);      % upper row index
        if ((i+winHeight)>szIm(1))
            idxU = szIm(1);
        end
        winL = seamPtEst+seamCol(i) - 1 - winWidth;  % pixel col to the left of window
        winR = seamCol(i)+winWidth;            % pixel col to the right of window
        midPatch = bsxfun(@times,double(imL(idxL:idxU,winL+1:winL+winWidth)),weightsLeft) + ...
                     bsxfun(@times,double(imStR(idxL:idxU,seamCol(i):winR-1)),weightsRight) ;
        % stitched image consists of left, overlapped portion and right image
        imStitch(idxL:idxU,:) = [imL(idxL:idxU,1:winL) midPatch imStR(idxL:idxU,winR:end)];
    end
else
    imStitch = uint8(zeros(szIm(1),2*szIm(2)));
    for i = 1:szIm(1)
        imStitch(i,:) = [imL(i,1:seamPtEst+seamCol(i)-1) imStR(i,seamCol(i):end)];
    end
end

%% Remove black pixels at the RHS
%  Find the leftmost column where there are more than half the pixels that
%  are black
%  Then remove all of the columns starting from this one.
for i = szIm(2)+1:2*szIm(2)
    blackPx = find(imStitch(:,i) == 0);
    if (length(blackPx) > szIm(1)/2)
        break;
    end
end

% Here i corresponds to the first column to be removed
imStitch = imStitch(:,1:i-1);

% FINALLY ... displays the stitched image
figure(3); imshow(imStitch);  title('Stitched image');
set(gcf, 'Position', get(0, 'Screensize'));

% END OF SCRIPT :)