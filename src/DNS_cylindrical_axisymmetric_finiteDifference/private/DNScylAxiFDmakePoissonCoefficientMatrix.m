function PoissonCoefficientMatrix = DNScylAxiFDmakePoissonCoefficientMatrix(R,Z,periodicZ)
%NSCYLAXIFDMAKEPOISSONCOEFFICIENTMATRIX makes the coefficient matrix to do the poisson equation in equation A.6.
% POISSONCOEFFICIENTMATRIX = NSCYLAXIFDMAKEPOISSONCOEFFICIENTMATRIX(R,Z,PERIODICZ)
% Makes the coefficient matrix, POISSONCOEFFICIENTMATRIX, to do the poisson
% equation for Gamma in equation A.6 in Appendix Z (Lopez and Shen 1998)
% that relates it to eta and R. R and Z are the coordinate matrices (made
% by meshgrid) where coordinates are in (z,r) order. PERIODICZ is a flag
% indicating whether the axial boundary condition is periodic (true) or not
% (false). POISSONCOEFFICIENTMATRIX is constructed such that
% R(:).*eta(:) = POISSONCOEFFICIENTMATRIX * psi(:).
% POISSONCOEFFICIENTMATRIX is sparse and only deals with the interior
% points.
%
% References:
%   1)  JM Lopez and J Shen. An Efficient Spectral-Projection Method for
%       the Navier-Stokes Equations in Cylindrical Geometries - Part I -
%       Axisymmetric Cases. Journal of Computational Physics. 139, 308-326
%       (1998)
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

% Equation A.6 works out to be
%
% -psi_rr + (1/r)*psi_r - psi_zz = r*eta. PoissonCoefficientMatrix will
% describe the right hand side.

% The coefficient matrices for each term will be made separately because
% the second term requires special handling due to the division by r (don't
% want any 0/0 -> NaN's filling the matrix) and the third term requires
% special handling depending on whether the axial boundary condition is
% periodic or not.

% Make the first term.

negPsi_rr_coeffs = -d2fdx2FiniteDifferenceNonUniformMatrix(R,R(1,:),2);

% Make the second term without dividing by R first. Then, only divide the
% non-zero entries by R at tthose points (doing it this way avoids dividing
% any zero entries by R that would then either stay zero or become NaN).
% This has to be done using a for loop.

psi_r_overR_coeffs = dfdxFiniteDifferenceNonUniformMatrix(R,R(1,:),2);
for ii=find(psi_r_overR_coeffs)'
    [row col] = ind2sub(size(psi_r_overR_coeffs),ii);
    psi_r_overR_coeffs(row,col) = psi_r_overR_coeffs(row,col) / R(col);
end

% Make the third term. The way that this is done depends on whether the
% axial boundary conditions are periodic or not.

if periodicZ
    dz = diff(Z(1:2,1));
    negPsi_zz_coeffs = -d2fdx2FiniteDifferencePeriodicMatrix(R,dz,2);
else
    negPsi_zz_coeffs = -d2fdx2FiniteDifferenceNonUniformMatrix(R,Z(:,1),1);
end

% Make PoissonCoefficientMatrix by adding the different terms together.

PoissonCoefficientMatrix =  negPsi_rr_coeffs + psi_r_overR_coeffs + negPsi_zz_coeffs;

% Now, some points not in the interior were included and even for those in
% the interior, some may be using points out of the interior. These points
% need to be eliminated. So, we take all the nonzero entries and get the
% row and column indices. If the r-values (and z-values if the system is
% not axially periodic) of either point (row or column) are on the borders,
% they are eliminated by setting their coefficient to zero.

indices = find(PoissonCoefficientMatrix);
[rows cols] = ind2sub(size(PoissonCoefficientMatrix),indices);

if ~periodicZ
%     indices2 = find(R(rows) == R(1,1) | R(rows) == R(1,end) | Z(rows) == Z(1,1) | Z(rows) == Z(end,1));
    indices2 = find(R(rows) == R(1,1) | R(rows) == R(1,end) | Z(rows) == Z(1,1) | Z(rows) == Z(end,1) ...
        | R(cols) == R(1,1) | R(cols) == R(1,end) | Z(cols) == Z(1,1) | Z(cols) == Z(end,1));
else
%     indices2 = find(R(rows) == R(1,1) | R(rows) == R(1,end));
    indices2 = find(R(rows) == R(1,1) | R(rows) == R(1,end) ...
        | R(cols) == R(1,1) | R(cols) == R(1,end));
end
PoissonCoefficientMatrix(indices(indices2)) = 0;
