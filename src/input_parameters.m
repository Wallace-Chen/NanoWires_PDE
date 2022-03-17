% this file sets the various input parameters to be used.
global degree_of_polygon;
global vector_of_side_lengths;
global vector_of_V0;
global vector_of_masses;
global vector_of_eps;
global number_of_triangles;
global number_of_interfaces;
global outer_shell_side_length;
% a number indicating the degree of the polygon, e.g. 3 means triangle, 6
% means hexagon, 0 means circle
degree_of_polygon = 6;

% a vector containing the side lengths of each of the interfaces in units
% of l0 nm [see constat.m for l0]. This vector is ordered as follows: the
% first element contains the side length of the outermost shell, and the
% subsequent elements contains side lengths of the remaining interfaces as
% one "works inwards" towards the innermost core region
vector_of_side_lengths = [4.5 3];

% a vector containing the bandedge energies in units of eV [see constat.m
% for eV]. This vector is ordered as follows: the first element contains
% the bandedge energy of the outermost region, and the subsequent elements
% contain bandedge energies of the remaining regions as one "works inwards"
% towards the innermost core region
vector_of_V0 = [0.5 0.0];

% a vector containing the effective mass of each of the nanowire regions in
% units of m_0 [see constat.m for m0]. This vector is ordered as follows:
% the first element contains the effective mass of the outermost region,
% and the subsequent elements contain effective masses of the remaining
% regions as one "works inwards" towards the innermost core region
vector_of_masses = [0.2-0.12*0.3 0.2];

% a vector containing the static dielectric of each of the nanowire
% regions. This vector is ordered as follows: the first element contains
% the static dielectric of the outermost region, and the subsequent
% elements contain static dielectrics of the remaining regions as one
% "works inwards" towards the innermost core region
vector_of_eps = [9.28-0.61*0.3 9.28];

% number of triangles to be used in the finite element mesh
number_of_triangles = 500;

number_of_interfaces = length(vector_of_side_lengths);
outer_shell_side_length = vector_of_side_lengths(1);