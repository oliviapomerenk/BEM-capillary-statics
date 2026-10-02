function [x,y,N,perimeter,Rmin,ds] = processDrawnParticleContour(pos,L,min_edge_distance, C_BEM, C_curv, Nmin, Nmax, devFactor)

% Also used ChatGPT


% ============================================================
% Process a user-drawn particle contour
%
% The raw contour is intentionally reduced to a coarse,
% rounded representation before determining the BEM resolution.
%
% Inputs:
%   pos                - Nx2 array of points from drawfreehand
%   L                  - characteristic particle radius
%   min_edge_distance  - minimum distance between particles
%
% Outputs:
%   x,y                - equally spaced boundary points for BEM
%   N                  - number of BEM boundary points
%   perimeter          - physical particle perimeter
%   Rmin               - minimum radius of curvature
%   ds                 - BEM boundary point spacing
%
% ============================================================


% This controls how much of the small-scale detail in the
% hand-drawn contour is retained (smaller = coarser)
nShape = 20;


%% ------------------------------------------------------------
% 1. Clean up raw points
% ------------------------------------------------------------

xRaw = pos(:,1);
yRaw = pos(:,2);

xRaw = xRaw(:);
yRaw = yRaw(:);

% Remove consecutive duplicate points
keep = [true; diff(xRaw) ~= 0 | diff(yRaw) ~= 0];

xRaw = xRaw(keep);
yRaw = yRaw(keep);

if length(xRaw) < 3
    error('The particle contour must contain at least three distinct points.');
end

% Close contour
if xRaw(1) ~= xRaw(end) || yRaw(1) ~= yRaw(end)
    xRaw(end+1) = xRaw(1);
    yRaw(end+1) = yRaw(1);
end


%% ------------------------------------------------------------
% 2. Fit a smooth closed spline
% ------------------------------------------------------------

pp = cscvn([xRaw'; yRaw']);


%% ------------------------------------------------------------
% 3. Sample the spline COARSELY
% ------------------------------------------------------------


tShape = linspace(pp.breaks(1),pp.breaks(end),nShape+1);

xyShape = fnval(pp,tShape);

xShape = xyShape(1,:)';
yShape = xyShape(2,:)';

% Remove duplicate endpoint
xShape(end) = [];
yShape(end) = [];


%% ------------------------------------------------------------
% 4. Scale the coarse shape so its maximum diameter is 2L
% ------------------------------------------------------------

% Maximum diameter of the coarse shape

maxDistance = 0;

for i = 1:nShape

    distances = sqrt( ...
        (xShape(i) - xShape).^2 + ...
        (yShape(i) - yShape).^2 );

    maxDistance = max(maxDistance,max(distances));

end

if maxDistance <= 0
    error('The particle contour has zero diameter.');
end

scaleFactor = 2*L/maxDistance;

xShape = scaleFactor*xShape;
yShape = scaleFactor*yShape;


%% ------------------------------------------------------------
% 5. Fit a spline to the COARSE shape
%
% This is important: from this point onward, the geometry is
% the rounded/coarse version, not the original hand-drawn curve.
% ------------------------------------------------------------

% Close the coarse contour for cscvn
xClosed = [xShape; xShape(1)];
yClosed = [yShape; yShape(1)];

ppShape = cscvn([xClosed'; yClosed']);


%% ------------------------------------------------------------
% 6. Sample the coarse shape for geometry calculations
% ------------------------------------------------------------

nCalc = 500;

tCalc = linspace(ppShape.breaks(1), ...
                 ppShape.breaks(end), ...
                 nCalc);

xyCalc = fnval(ppShape,tCalc);

xCalc = xyCalc(1,:)';
yCalc = xyCalc(2,:)';


%% ------------------------------------------------------------
% 7. Calculate perimeter
% ------------------------------------------------------------

dsCalc = sqrt(diff(xCalc).^2 + diff(yCalc).^2);

sCalc = [0; cumsum(dsCalc)];

perimeter = sCalc(end);


%% ------------------------------------------------------------
% 8. Calculate minimum radius of curvature
% ------------------------------------------------------------

pp1 = fnder(ppShape,1);
pp2 = fnder(ppShape,2);

d1 = fnval(pp1,tCalc);
d2 = fnval(pp2,tCalc);

dxdt = d1(1,:)';
dydt = d1(2,:)';

d2xdt2 = d2(1,:)';
d2ydt2 = d2(2,:)';

numerator = abs(dxdt.*d2ydt2 - dydt.*d2xdt2);

denominator = (dxdt.^2 + dydt.^2).^(3/2);

kappa = numerator ./ denominator;

kappa = kappa(isfinite(kappa) & kappa > 0);

if isempty(kappa)
    error('Could not calculate the curvature of the particle contour.');
end

kappaMax = max(kappa);

Rmin = 1/kappaMax;


%% ------------------------------------------------------------
% 9. Determine target BEM spacing
% ------------------------------------------------------------

% if isfinite(min_edge_distance)
% 
%     smallestLengthScale = min(Rmin,min_edge_distance);
% 
% else
% 
%     smallestLengthScale = Rmin;
% 
% end

if isfinite(min_edge_distance)

    smallestLengthScale = min(L,min_edge_distance);

else

    smallestLengthScale = L;

end

ds = min(C_BEM*smallestLengthScale, C_curv*Rmin);


%% ------------------------------------------------------------
% 10. Determine number of BEM points
% ------------------------------------------------------------


N = ceil(devFactor * perimeter/ds);

N = min(Nmax,max(N,Nmin));


%% ------------------------------------------------------------
% 11. Resample rounded shape uniformly in arclength
% ------------------------------------------------------------

sNew = (0:N-1)' * perimeter/N;

x = interp1(sCalc,xCalc,sNew,'pchip');
y = interp1(sCalc,yCalc,sNew,'pchip');


%% ------------------------------------------------------------
% 12. Center particle at (0,0)
% ------------------------------------------------------------

lastwarn('', ''); 

P = polyshape(x,y);

% check for self-intersections (seems to be the only reason a warning is
% thrown)
[~, warningID] = lastwarn;
if ~isempty(warningID)
    error('Custom contour is not a smooth, closed contour (likely self-intersection).');
end


[xcentroid,ycentroid] = centroid(P);

x = x - xcentroid;
y = y - ycentroid;


end