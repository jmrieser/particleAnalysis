function coneOut = generateRightCircularCone(x,y,x0,y0,r,h)

% x and y should be outputs from mesgrid
% (x0,y0) is the location of the center of the cone
% r = radius of the cone
% h = height of the cone

coneOut = h - h/r*sqrt((x-x0).^2 + (y-y0).^2);

coneOut(coneOut < 0) = 0;