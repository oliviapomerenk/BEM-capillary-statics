function [Tx,Ty,Tz] = CalcTorque(mu,w,sigma,r_x,r_y,xcm,ycm,nx,ny,H_coll,zcm_height,hs,Tz_baseline)

mu = mu(:);
hs = hs(:);
w = w(:);
r_x = r_x(:);
r_y = r_y(:);
nx = nx(:);
ny = ny(:);
H_coll = H_coll(:);

% Normals and Tangents
nfx = -nx;
nfy = -ny;

tfx = -nfy;
tfy = nfx;

hn = -mu;

% O(h^2) expansion again
dFx = sigma * ( (1 + 0.5*(hs.^2 - hn.^2)).*nfx - (hn.*hs).*tfx );
dFy = sigma * ( (1 + 0.5*(hs.^2 - hn.^2)).*nfy - (hn.*hs).*tfy );
dFz = sigma * hn;

% lever arms 
Rx = r_x - xcm;
Ry = r_y - ycm;
Rz = H_coll - zcm_height;

% torques: T = sum( w .* (R x dF) )
Tx = sum( w .* ( Ry .* dFz - Rz .* dFy ) );
Ty = sum( w .* ( Rz .* dFx - Rx .* dFz ) );
Tz = sum( w .* ( Rx .* dFy - Ry .* dFx ) );

% baseline subtraction
Tz = Tz - Tz_baseline;



end

