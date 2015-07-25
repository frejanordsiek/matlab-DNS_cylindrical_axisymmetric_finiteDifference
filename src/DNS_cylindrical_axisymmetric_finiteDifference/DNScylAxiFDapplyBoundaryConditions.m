function [Gamma eta psi] = DNScylAxiFDapplyBoundaryConditions(Gamma, eta, psi, t, boundaryInformation)
%NSCYLAXIFDAPPLYBOUNDARYCONDITIONS applies the boundary conditions.
% [GAMMA ETA PSI] = NSCYLAXIFDAPPLYBOUNDARYCONDITIONS(GAMMA,ETA,PSI,T,BOUNDARYINFORMATION)
% Applies the boundary conditions at time T using the boundary information
% in the boundaryInformation structure BOUNDARYINFORMATION. GAMMA, ETA, and
% PSI are the angular velocity, theta component of vorticity, and the 2D
% streamfunction respectively. A boundaryInformation structure has five
% fields:
%
% periodicZ : Whether the system is axially periodic (true) or not (false)
%
% rMinBoundary, rMaxBoundary, zMinBoundary, zMaxBoundary : structures that
%   describe the boundaries at the minimum r-value, maximum r-value,
%   minimum z-value, and maximum z-value respectively. Each has a mandatory
%   field 'Type' which tells what type of boundary it it. It must be 'wall'
%   for all boundaries except the minimum r boundary where it can be 'axis'
%   if it is the axis of rotation (its r-value had better be zero or funny
%   results will occur). If it is a 'wall', then the field 'GammaValues'
%   must be a function handle that takes a single argument, the time T, and
%   returns the Gamma values that are the wall is supposed to have at each
%   grid point on the wall at that time.
%
% References:
%   1)  JM Lopez and J Shen. An Efficient Spectral-Projection Method for
%       the Navier-Stokes Equations in Cylindrical Geometries - Part I -
%       Axisymmetric Cases. Journal of Computational Physics. 139, 308-326
%       (1998)
%       DOI: 10.1006/jcph.1997.5872
%
%
%
% Copyright 2012 Freja Nordsiek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% psi = 0 on the walls and axis, normal and tangential derivatives of psi
% are zero on walls, eta = 0 on the axis, Gamma = 0 the axis, and
% Gamma = r*v on walls where v is the azimuthal velocity of the wall at the
% given point (Lopez and Shen 1998).

% Do axial boundary conditions first, then radial.

% If the system is not axially periodic, then the boundary conditions need
% to be enforced on the ends. Otherwise, there is nothing that needs to be
% done.

if ~boundaryInformation.periodicZ
    
    % The boundary at the minimum z must be a wall.
    
    if ~strcmpi(boundaryInformation.zMinBoundary.Type,'wall')
        error('z-min boundary type not supported.');
    end
    
    % Since the boundary is a wall, set psi equal to zero at it and the
    % next value of z in and set Gamma to what the GammaValues function
    % handle returns.

    Gamma(end,:) = boundaryInformation.zMinBoundary.GammaValues(t);
    psi((end-1):end,:) = 0;
%     psi(end,:) = 0;
    
    % The boundary at the maximum z must be a wall.
    
    if ~strcmpi(boundaryInformation.zMaxBoundary.Type,'wall')
        error('z-max boundary type not supported.');
    end
    
    % Since the boundary is a wall, set psi equal to zero at it and the
    % next value of z in and set Gamma to what the GammaValues function
    % handle returns.

    Gamma(1,:) = boundaryInformation.zMaxBoundary.GammaValues(t);
    psi(1:2,:) = 0;
%     psi(1,:) = 0;
    
end

% The boundary at the minimum r can either be the axis or a wall. If it is
% an axis; Gamma, eta, and psi are set to zero on it. If it is a wall, psi
% is set to zero at it and the next r (makes the normal and tangential
% derivatives of it zero) and Gamma is set to the values gotten from the
% GammaValues function handle at time t.

switch lower(boundaryInformation.rMinBoundary.Type)
    case 'axis'
        Gamma(:,1) = 0;
        eta(:,1) = 0;
        psi(:,1) = 0;
    case 'wall'
        Gamma(:,1) = boundaryInformation.rMinBoundary.GammaValues(t);
        psi(:,1:2) = 0;
%         psi(:,1) = 0;
    otherwise
        error('r-min boundary type not supported.');
end

% The boundary at the maximum r must be a wall.

if ~strcmpi(boundaryInformation.rMaxBoundary.Type,'wall')
    error('r-max boundary type not supported.');
end

% Since the boundary is a wall, set psi equal to zero at it and the next
% value of r in and set Gamma to what the GammaValues function handle
% returns.

Gamma(:,end) = boundaryInformation.rMaxBoundary.GammaValues(t);
psi(:,(end-1):end) = 0;
% psi(:,end) = 0;
