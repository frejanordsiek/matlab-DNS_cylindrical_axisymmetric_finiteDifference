function D2fDx2 = d2fdx2FiniteDifferencePeriodic(f,h,dim)
%D2FDX2FINITEDIFFERENCEPERIODIC does second finite difference on a periodic uniform grid.
% D2FDX2 = D2FDX2FINITEDIFFERENCEPERIODIC(F, H, DIM)
% Does the second centered finite difference on a periodic uniform grid
% of F along direction DIM (1 is row, 2 is column, and 3 is page) where H
% is the grid spacing. The formula is, for a one dimensional F (extends to
% two and three dimensions by just re-application of it to each
% row/col/page/whatever), given by
%
% D2f_i/Dx2 = (f_(i+1) - 2*f_i + f_(i-1)) / h^2
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

% All we have to do is simply add a grid points on each that are the grid
% points on the other side and just pass it onto the function that handles
% the non-periodic case.

switch dim
    case 1
        D2fDx2 = d2fdx2FiniteDifference(cat(1,f(end,:,:),f,f(1,:,:)),h,dim);
    case 2
        D2fDx2 = d2fdx2FiniteDifference(cat(2,f(:,end,:),f,f(:,1,:)),h,dim);
    case 3
        D2fDx2 = d2fdx2FiniteDifference(cat(3,f(:,:,end),f,f(:,:,1)),h,dim);
    otherwise
        error('DIM can only be 1, 2, or 3.');
end
