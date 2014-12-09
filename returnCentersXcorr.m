function [X,Y,R] = returnCentersXcorr(image,coneRadius,cannyThreshold,radiusThreshold,imageThreshold,plotOn,dividingLine,iteration,totalIterations,adapt,points)
% Written by Jennifer Rieser, University of Pennsylvania 2014


%%% Output:
% X = particle X positions
% Y = particle Y positions
% R = heights of distance map at (X,Y)


%%% Input:
% image = grayscale uint8 image

% coneRadius =   radius of right-circular cone (also height of cone) should be a
%           slightly less or equal to particle size. This cone will be
%           convolved with distance map to improve peaks and give them a paraboloid height profile.

% cannyThreshold = threshold used for canny edge detection.
%
% radiusThreshold = minimum feature size.  Any peaks in distance map with
%                    height less than this value will be excluded from final result.

% imageThreshold = Pixel intensity value corresponding to background.
%             Any peaks in the distance map will be excluded if the
%             original image intensity exceeds this threshold. Note: this
%             is assuming particles are dark and background is light.  If
%             this is not the case, then input imcomplement(image) instead
%             of image.

% plotOn = true will produce a plot of identified particles overlaid on top of original image.

% dividingLine = only matters for plotting -- plotted particles will be color coded based on
%                big = larger than dividingLine, small = smaller than dividingLine

% iteration = 1 unless this is enclosed in a for-loop.  Only used for plotting.

% totalIterations = 1 unless this is enclosed in a for-loop.

% adapt = 0 for no image histogram equalization, 1 otherwise.  Sometimes this improves feature finding.

% points = number of points to keep on each side of maximum value for sub-pixel center identification

%%% Default values... %%%
if nargin < 2 || isempty(coneRadius)
    coneRadius = 10;   
end

[xx,yy] = meshgrid(size(image));
coneOut = generateRightCircularCone(xx,yy,round(mean(size(xx,1))),round(mean(size(xx,2))),coneRadius,coneRadius);
coneOut = coneOut';

% threshold for canny edge detection
if nargin < 3 || isempty(cannyThreshold)
    cannyThreshold = 0.085;
end


% minimum radius value required to be considered a particle
if nargin < 4 || isempty(radiusThreshold)
    radiusThreshold = 7;
end

maxRadiusThreshold = 3*radiusThreshold;

% pixel intensity value above which a particle is ignored
if nargin < 5 || isempty(imageThreshold)
    imageThreshold = 100;
end

% true if plot is desired to evaluate accuracy... does make code run
% slower, especially for large runs
if nargin < 6 || isempty(plotOn)
    plotOn = true;
end


% for bidisperse, mean of radii
if nargin < 7 || isempty(dividingLine)
    dividingLine = 14;
end

% current iteration
if nargin < 8 || isempty(iteration)
    iteration = 1;
end

if nargin < 9 || isempty(totalIterations)
    totalIterations = 1;
    plotOn = true;
end

if nargin < 10 || isempty(adapt)
    adapt = 0;
end

if nargin < 11 || isempty(points)
    points = 8;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using a 3x3 grid near pixelated center to fit a paraboloid ... this is
% the coefficient matrix
aMat = [1 1 -1 -1 1 1; 0 1 0 0 1 1;1 1 1 1 1 1; 1 0 0 -1 0 1; 0 0 0 0 0 1; 1 0 0 1 0 1; 1 1 1 -1 -1 1; 0 1 0 0 -1 1; 1 1 -1 1 -1 1];
combinations = [2 3 4 5 6 8; 1 2 4 5 6 8; 2 4 5 6 7 8; 2 4 5 6 8 9]; %points to include in the calculation: grid numbering: [1 2 3; 4 5 6; 7 8 9]
aMatInv = cell(length(combinations(:,1)),1);

for j = 1:length(combinations(:,1));
    aMatInv{j} = inv(aMat(combinations(j,:),:));
end

s = size(image);

%histogram equalization
if adapt
    adapthisteq(image);
end

% image gradient
[gradmag,~] = imgradient(image,'prewitt');
gm8 = uint8(gradmag);

% edge detection
E = edge(gm8,'canny',cannyThreshold);

% compute distance map
distmap = bwdist(E);

% throw away features that are too bright
distmap(image > imageThreshold) = 0;
test = fftshift(real(ifft2(fft2(double(coneOut)).*fft2(double(distmap)))))/(s(1)*s(2));


% find maxima of distance map
[x,y] = find(imregionalmax(test) & distmap > radiusThreshold & distmap < maxRadiusThreshold);
vals = distmap(imregionalmax(test) & distmap > radiusThreshold & distmap < maxRadiusThreshold);


% only keep features at least radiusThreshold away from boundaries 
idx = find(x >= radiusThreshold & x <= s(1)-radiusThreshold & y >= radiusThreshold & y <= s(2)-radiusThreshold);
x = x(idx);
y = y(idx);
vals = vals(idx);


L = length(x);
D = zeros(L,L);
% find distances between all maxima locations
for i=1:L
    for j=(i+1):L
        D(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        D(j,i) = D(i,j);
    end
end

% eliminate multiple points from the same distance map peak  
for i=1:length(x)
    idx = find(D(i,:)<2.5*radiusThreshold);
    idxVals = vals(idx);
    [~,maxIdx] = max(idxVals);
    if length(maxIdx) > 1
        maxIdx = maxIdx(1);
    end
    getRidOf = setdiff(idx,idx(maxIdx));
    vals(getRidOf) = -1;
end

x = x(vals>0);
y = y(vals>0);
r = zeros(size(x));

vals = vals(vals>0);
isLarge = vals > dividingLine;

for i=1:length(x)
    
    % paraboloid determination and averaging
    x1 = max(x(i)-points,1);
    x2 = min(x(i)+points,s(1));
    y1 = max(y(i)-points,1);
    y2 = min(y(i)+points,s(2));
    xVec = double(y1:y2);
    yVec = double(x1:x2);
    zMat = double(test(yVec,xVec));
    [xx,yy] = meshgrid(xVec,yVec);
    
    zzz = zMat(points:points+2,points:points+2);
    zVec = [zzz(3,:) zzz(2,:) zzz(1,:)]';
    
    xmins = zeros(length(combinations(:,1)),1);
    ymins = xmins;
    for j = 1:length(combinations(:,1));
        ajvec = aMatInv{j}*zVec(combinations(j,:));
        [xmins(j),ymins(j)] = paraboloidMin(ajvec);
    end
    
    x(i) = xx(points+1,points+1)+mean(xmins);
    y(i) = yy(points+1,points+1)+mean(ymins);
    r(i) = vals(i);
    
end

% there is a systematic shift in each of these that comes from the edge
% detection on the image gradient.

X = x+1;
Y = y+1;
R = r+1;

% plot the identified particles, color coded based on size
if plotOn
    thetas = linspace(0,2*pi,200);
    figure(6),clf
    imshow(image)
    hold on,title(num2str(iteration))
    plot(X(isLarge),Y(isLarge),'m.')
    plot(X(~isLarge),Y(~isLarge),'g.')
    for i=1:length(X)
        if isLarge(i)
            plot(X(i)+r(i)*cos(thetas),Y(i)+r(i)*sin(thetas),'-m','linewidth',2)
        else
            plot(X(i)+r(i)*cos(thetas),Y(i)+r(i)*sin(thetas),'g-','linewidth',2)
        end
    end
    drawnow;
    hold off
end
