function deformation = triangleStrains(triangleXs,triangleYs)
% Written by Jennifer Rieser, University of Pennsylvania 2014


%%% inputs:
% triangleXs and triangleYs are Nx1 cells, where N = total number of triangles
% triangleXs{n} is a T x 3 matrix of the x-coordinates of each triangle vertex
% triangleYs{n} is a T x 3 matrix of the y-coordinates of each triangle vertex

%%% output
% deformation = structure containing NxT matrices of the following quantities:
% deformation.eps_xx = xx-component of strain tensor
% deformation.eps_xy =  xy-component of strain tensor
% deformation.eps_yy =  yy-component of strain tensor
% deformation.rotation = rotation of triangle
% deformation.pureShear = pure shear of triangle
% deformation.simpleShear = simple shear of triangle
% deformation.dilation = dilation of triangle -- 1st invariant of strain tensor
% deformation.deviatoricStrain = deviatoricStrain -- invariant of strain tensor
% deformation.maxEig = maximum eigenvalue
% deformation.minEig = minimum eigenvalue
% deformation.maxEigVecX = x-component of eigenvector associated with maximum eigenvalue
% deformation.maxEigVecY = y-component of eigenvector associated with maximum eigenvalue
% deformation.minEigVecX = x-component of eigenvector associated with minimum eigenvalue 
% deformation.minEigVecY = y-component of eigenvector associated with minimum eigenvalue
% deformation.thetaMax = angle of max eigenvector
% deformation.maxShear = maximum shear strain
% deformation.thetaMaxShear = maximum shear direction


%%% Description of calculation
% This function calculates the strains for each triangle as s function of time,
% assuming that the displacements are linear in x and y.  Let u =
% x-displacement and v = y-displacement, then:

% u(x,y) = constant1 + a*x + b*y (1)
% v(x,y) = constant2 + c*x + d*y (2)

% constant1 and constant2 describe triangle centroid displacement
% [a b; c d] = displacement gradient tensor, where
% a = du/dx; b = du/dy; c = dv/dx; d = dv/dy;

% putting the x and y displacements of the three vertices into (1) and (2),
% there are six equations and six unknowns, which yields the following
% matrix equations:

% shift vertices so that original triangle centroid = 0. (If the triangle
% is not shifted, a,b,c,d are unchanged, but constant1 and constant2 are
% changed and no longer equal to centroid shift)

% Let xis = xi - xcm; yis = yi-ycm; xcm = mean(xi) and ycm = mean(yi) for
% triangles.

% [u1]   [ 1 x1s y1s ] [constant1]
% |u2| = | 1 x2s y2s | |    a    |
% [u3]   [ 1 x3s y3s ] [    b    ]

% [v1]   [ 1 x1s y1s ] [constant2]
% |v2| = | 1 x2s y2s | |    c    |
% [v3]   [ 1 x3s y3s ] [    d    ]

% to solve for a,b,c,d,constant1,constant2, need to invert the following
% matrix:

%         [ 1 x1s y1s ]
%  X =    | 1 x2s y2s |
%         [ 1 x3s y3s ]

% From Mathematica, 
%                                                                                [    2*Area/3       2*Area/3       2*Area/3     ]
%  X^{-1} = (-x2s y1s + x3s y1s + x1s y2s - x3s y2s - x1s y3s + x2s y3s)^(-1) *  |   (y2s - y3s)    (y3s - y1s)    (y1s - y2s)   |
%                                                                                [   (x3s - x2s)    (x1s - x3s)    (x2s - x1s)   ]

% where prefactor = 1/(2*Area) [see http://mathworld.wolfram.com/PolygonArea.html;http://mathworld.wolfram.com/TriangleArea.html eq(17)].

% Define alpha, beta, and gamma in the following way:

%                         [alpha1  alpha2  alpha3]
%  X^{-1} = 1/(2*Area) *  |beta1   beta2   beta3 |
%                         [gamma1  gamma2  gamma3]

% where prefactor = 1/(2*Area).  Define alpha, beta, and gamma in the following way:


N = length(triangleXs);
[T,~] = size(triangleXs{1});

eps_xx = zeros(N,T);
gamma_xy = zeros(N,T);
eps_yy = zeros(N,T);
omega = zeros(N,T);
deviatoricStrain = zeros(N,T);
maxEig = zeros(N,T);
minEig = zeros(N,T);
maxShear = zeros(N,T);
maxShearTheta = zeros(N,T);

thetaMax = zeros(N,T);
maxevx = zeros(N,T);
maxevy = zeros(N,T);

minevx = zeros(N,T);
minevy = zeros(N,T);


for n = 1:N
    for t = 1:T-1
        
        
        
        x = triangleXs{n}(t,1:3);
        y = triangleYs{n}(t,1:3);
        
        xnew = triangleXs{n}(t+1,1:3);
        ynew = triangleYs{n}(t+1,1:3);
        
        
        % original triangle needs to have centroid located at zero, so
        % shift both triangles by initial centroid
        
        xCentroid = mean(x);
        yCentroid = mean(y);
        
        x = x-xCentroid;
        y = y-yCentroid;
        
        xnew = xnew - xCentroid;
        ynew = ynew - yCentroid;
        
        %displacements in x and y
        u = xnew-x;
        v = ynew-y;
        
        % define betas, and gammas
        % i,j,k refer to vertices
        beta(1) =  y(2)-y(3);
        beta(2) =  y(3)-y(1);
        beta(3) =  y(1)-y(2);
        
        gamma(1) =  x(3)-x(2);
        gamma(2) =  x(1)-x(3);
        gamma(3) =  x(2)-x(1);
        
        % note: this definition of area is signed, so ordering of vertices
        % does not have to all be the same.  twoA > 0 if indices are
        % ccw-ordered and twoA < 0 if indices are cw-ordered.
        
        twoA = sum(x.*beta);
        
        %dudx
        eps_xx(n,t) = 1/twoA*(sum(beta.*u));
        
        %dvdy
        eps_yy(n,t) = 1/twoA*(sum(gamma.*v));
        
        %dudy
        dudy = 1/twoA*(sum(gamma.*u));
        
        %dvdx
        dvdx = 1/twoA*(sum(beta.*v));
        
        % for small deformations, strain tensor is approximately symmetric part of
        % dispacement gradient tensor and rotation is anti-symmetric part:
        gamma_xy(n,t) = 1/2*(dudy + dvdx);
        omega(n,t) = 1/2*(dvdx - dudy);
        em = 0.5*(eps_xx(n,t)+eps_yy(n,t));
        strain = [eps_xx(n,t) gamma_xy(n,t);gamma_xy(n,t) eps_yy(n,t)];
        deviatoricStrain(n,t) = sqrt(0.5*abs(sum(sum((strain-em*eye(2,2)).*(strain-em*eye(2,2))))));
        maxShear(n,t) = sqrt((eps_xx(n,t)-eps_yy(n,t)).^2 +gamma_xy(n,t).^2);
        maxShearTheta(n,t) = mod(atan(-(eps_xx(n,t)-eps_yy(n,t))/gamma_xy(n,t)),pi);
        
        [vecs,vals] = eig(strain);
        
       [maxeig,ix] = max(real(diag(vals)));
        maxEig(n,t) = maxeig;
        thetaMax(n,t) = mod(atan2(vecs(2,ix),vecs(1,ix)),pi);
        maxevx(n,t) = real(vecs(1,ix));
        maxevy(n,t) = real(vecs(2,ix));
        
        [mineig,ix] = min(real(diag(vals)));
        minEig(n,t) = mineig;
        minevx(n,t) = real(vecs(1,ix));
        minevy(n,t) = real(vecs(2,ix));
    end
end

deformation.eps_xx = eps_xx;
deformation.eps_xy = gamma_xy;
deformation.eps_yy = eps_yy;
deformation.rotation = omega;
deformation.pureShear = eps_xx - eps_yy;
deformation.simpleShear = 2*gamma_xy;
deformation.dilation = eps_xx+eps_yy;
deformation.deviatoricStrain = deviatoricStrain;
deformation.maxEig = maxEig;
deformation.minEig = minEig;
deformation.maxEigVecX = maxevx; 
deformation.maxEigVecY = maxevy;
deformation.minEigVecX = minevx; 
deformation.minEigVecY = minevy;
deformation.thetaMax = thetaMax;
deformation.maxShear = maxShear;
deformation.thetaMaxShear = thetaMaxShear;