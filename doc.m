%% meshed_sinusoidal_icosahedron
%
% Function to compute and display
% a meshed sinusoidal icosahedron.
%
% Author : nicolas.douillet (at) free.fr 2016-2024.
%
%% Syntax
%
% meshed_sinusoidal_icosahedron;
%
% meshed_sinusoidal_icosahedron(sampling);
%
% meshed_sinusoidal_icosahedron(sampling, w);
%
% meshed_sinusoidal_icosahedron(sampling, w, option_display);
%
% [M, T] = meshed_sinusoidal_icosahedron(sampling, w, option_display);
%
%% Description
%
% meshed_sinusoidal_icosahedron computes and displays a meshed sinusoidal
% icosahedron with parameters sampling = 60, w = 1, and option_display =
% true, by default.
%
% meshed_sinusoidal_icosahedron(sampling) samples at the value sampling
% each icosahedron basis triangle (20).
%
% meshed_sinusoidal_icosahedron(sampling, w) uses the given shape parameter
% w.
%
% meshed_sinusoidal_icosahedron(sampling, w, option_display) displays the
% result when option_display = *true/*1 and doesn't when option_display =
% false/0.
%
% [M, T] = meshed_sinusoidal_icosahedron(sampling, w, option_display)
% stores the point set coordinates and its corresponding triangulation
% respectively in M and T.
%
%% See also
%
% <https://fr.mathworks.com/help/matlab/ref/sphere.html sphere> |
% <https://fr.mathworks.com/help/matlab/ref/mesh.html mesh> |
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html trimesh> |
% <https://fr.mathworks.com/help/symbolic/mupad_ref/plot-icosahedron.html plot::icosahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73159-meshed-reuleaux-tetrahedron meshed_reuleaux_tetrahedron> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/69212-n-level-geoid n_level_geoid>
%
%% Input arguments
%
% - sampling : positive integer scalar double, sampling > 2. Remarkable value : sampling = 3 gives an icosahedron
%   
% - w : real scalar double, the shape parameter. Remarkable value : w = 0 gives a geoid.
%
% - option_display : logical *true (1) / false (0).
%
%% Output arguments
%
%        [|  |  | ]
% - M = [Mx My Mz], real matrix double, the coordinates of the generated point set. Size(M) = [nb_vertices,3].
%        [|  |  | ]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Default parameters values
meshed_sinusoidal_icosahedron;

%% Example #2
% Minimum sampling step
meshed_sinusoidal_icosahedron(6);

%% Example #3
% Negative shape parameter value
meshed_sinusoidal_icosahedron(60,-1);

%% Example #4
% Large shape parameter value
meshed_sinusoidal_icosahedron(60,3);