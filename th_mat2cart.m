function [xy] = th_mat2cart(ij,n)

% This function transforms matrix coordinates (i,j) to Cartesian coordinates
% (x,y) to allow for straightforward plotting of matrix space properties.
% It implements the mapping
%
%  f:N_n^2 -> R^2, (i,j) |-> f(i,j) := (x,y) := (j - (n+1)/2,-i +(n+1)/2) 
%
% where (i,j) refer to matrix row and column indices and (x,y) to Cartesian
% coordinates for a given square matrix size of size n x n. The Cartesian 
% coordinates are zero centered at the center of the matrix
%
% Inputs          
%       ij: (2 x 1) array of matrix row index i and column index j
%       n:  matrix size 
%
% Outputs
%      xy: (2 x 1) array of Cartesian x-ordinate x and y-ordinate y 
%
% Copyright (C) Dirk Ostwald, Lilla Horvath
% -------------------------------------------------------------------------

xy(1) =  ij(2) - ((n + 1)/2);
xy(2) = -ij(1) + ((n + 1)/2);

end
