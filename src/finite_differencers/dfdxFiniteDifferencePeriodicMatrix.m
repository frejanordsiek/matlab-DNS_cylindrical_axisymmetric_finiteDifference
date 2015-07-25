function coeffMatrix = dfdxFiniteDifferencePeriodicMatrix(f,h,dim)
%DFDXFINITEDIFFERENCEPERIODICMATRIX finds the coeff matrix for the first finite difference on a periodic uniform grid.
% COEFFMATRIX = DFDXFINITEDIFFERENCEPERIODICMATRIX(F, H, DIM)
% Creates the coefficient matrix, COEFFMATRIX which is sparse, to take the
% first centered finite difference of F on a periodic uniform grid along
% direction DIM (1 is row, 2 is column, and 3 is page) where H is the grid
% spacing. Then, the finite difference is calculated by
% dF(:)=COEFFMATRIX*F(:). The difference formula is, for a one dimensional
% F (extends to two and three dimensions by just re-application of it to
% each row/col/page/whatever), given by
%
% Df_i/Dx = (f_(i+1) - f_(i-1)) / 2*h
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

% Get the size vector of f and pad any missing column or page sizes with
% one.

siz = size(f);
siz((numel(siz)+1):3) = 1;

% Three vectors will be used in constructing the coefficient matrix. rows
% will hold the indices of the grid points that the difference is being
% evaluated at, cols will hold the indices of the grid points in f that
% will be used to evaluate the difference, and coeffs will hold the
% coefficients applied tot eh grid points specified in cols to find the
% differences at the grid points specified in rows.

% The total number of entries will be the product fo the dimensions of f.

rows = zeros(prod(siz),1);
cols = rows;
coeffs = rows;

% There are some slight differences to how it needs to be handled for each
% dimension type, so they are handled separately; but the basic algorithm
% is the same.
%
% Each row, column, and page are iterated over except the end ones along
% dimension dim. Then, as two entries are to be added in all three vectors,
% we simplify things can calculate the index to which (1:2) must be added
% to reach the proper entries in the vectors. This is simply the page index
% minus the starting plus the column index minus the starting times the
% number of pages we are looping over plus the row index minus the starting
% times the number of pages we are looping over times the number of columns
% we are looping over. Then, the two entries into rows are just the index
% associated with (ii,jj,kk). The two entries in cols are just the indices
% associated with (ii,jj,kk) but +/- 1 in direction dim. The coefficients
% are simply +/- 1/(2*h). The two boundaries and the interior are then
% handled separately.

switch dim
    case 1
        for ii=1:siz(1)
            for jj=1:siz(2)
                for kk=1:siz(3)
                    baseIndex = 2*((kk-1) + (jj-1)*siz(3) + (ii-1)*siz(3)*siz(2));
                    rows(baseIndex + (1:2)) = sub2ind(siz,ii,jj,kk)+[0;0];
                    switch ii
                        case 1
                            cols(baseIndex + (1:2)) = sub2ind(siz,[2;siz(1)],jj+[0;0],kk+[0;0]);
                        case siz(1)
                            cols(baseIndex + (1:2)) = sub2ind(siz,[1;siz(1)-1],jj+[0;0],kk+[0;0]);
                        otherwise
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[1;-1],jj+[0;0],kk+[0;0]);
                    end
                    coeffs(baseIndex + (1:2)) = [1; -1]/(2*h);
                end
            end
        end
    case 2
        for ii=1:siz(1)
            for jj=1:siz(2)
                for kk=1:siz(3)
                    baseIndex = 2*((kk-1) + (jj-1)*siz(3) + (ii-1)*siz(3)*siz(2));
                    rows(baseIndex + (1:2)) = sub2ind(siz,ii,jj,kk)+[0;0];
                    switch jj
                        case 1
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],[2;siz(2)],kk+[0;0]);
                        case siz(2)
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],[1;siz(2)-1],kk+[0;0]);
                        otherwise
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],jj+[1;-1],kk+[0;0]);
                    end
                    coeffs(baseIndex + (1:2)) = [1; -1]/(2*h);
                end
            end
        end
    case 3
        for ii=1:siz(1)
            for jj=1:siz(2)
                for kk=1:siz(3)
                    baseIndex = 2*((kk-1) + (jj-1)*siz(3) + (ii-1)*siz(3)*siz(2));
                    rows(baseIndex + (1:2)) = sub2ind(siz,ii,jj,kk)+[0;0];
                    switch kk
                        case 1
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],jj+[0;0],[2;siz(3)]);
                        case siz(3)
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],jj+[0;0],[1;siz(3)-1]);
                        otherwise
                            cols(baseIndex + (1:2)) = sub2ind(siz,ii+[0;0],jj+[0;0],kk+[1;-1]);
                    end
                    coeffs(baseIndex + (1:2)) = [1; -1]/(2*h);
                end
            end
        end
    otherwise
        error('DIM can only be 1, 2, or 3.');
end

% To make the sparse matrix a little faster to access, the entries will be
% sorted by rows and then cols. The resulting sorted vectors are then fed
% into the sparse command to make a sparse matrix with a number of rows and
% columns equal to the total number of elements of f with values coeffs at
% (rows,cols).

junk = sortrows([rows,cols,coeffs]);

coeffMatrix = sparse(junk(:,1),junk(:,2),junk(:,3),numel(f),numel(f));
