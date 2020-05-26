% bvp
% Copy the description from 1.1 here
%
% Parameters
% ==========
%    c A vector made up of the coefficients to the derivatives in the 2nd order 
%      LODE
%    x_int Input values to the boundary values
%    u_int Output values to the boundary values 
%
%    g The function isolated on the right side without derivative terms related %      to the function u(x).
%
%    n The number of that points that will evenly divide the interval.
%
% Return Values
% =============
%    x The vector xout is made up of n elements that correspond to the x value 
%      where the approximation for the solution u(x) is evaluated.
%    u The vector uout is made up of n elements that correspond to the 
%      evaluation of u(x) at xout

function [x_values, u_values] = bvp(c, x_int, u_int, g, n)
% Argument Checking

if ~all( size( c ) == [1, 3] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a 3-dimensional row vector' ) );
end

if ~all( size( x_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_int is not a 2-dimensional row vector' ) );
end

if ~all( size( u_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_int is not a 2-dimensional row vector' ) );
end

if ~isa( g, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument g is not a function handle' ) );
end

if ~isscalar( n ) || ( n ~= round( n ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n is not an integer' ) );
end

% define the step size h. h is needed in calculating the finite 
% difference formula

h = (x_int(2) - x_int(1)) / (n-1);

% create variables d1, d2, d3 to simplify the coefficients

d1 = (2 * c(1)) - (h * c(2));
d2 = (2 * (h ^ 2) * c(3)) - (4 * c(1));
d3 = ((2 * c(1)) + (h * c(2)));

% a vector of x values needs to be created which will be evaluated in the
% approximation

x_values = (linspace(x_int(1), x_int(2), n)');

% create the tri diagonal matrix of d coefficients

d_matrix = (diag( d2 * ones((n-2), 1)) + diag( d3 * ones((n-3), 1), 1) + diag( d1 * ones((n-3), 1), -1));

% create a vector of x_values that will be evaluated at solution function

mysterious_x_values = x_values(2: end -1);

% create target vector with forcing function involved

target_vector = (2 .* (h ^ 2) .* (g(mysterious_x_values)));

intermediate_zeroes_vector = zeros((n-2), 1);
intermediate_zeroes_vector(1) = (-d1 * u_int(1));
intermediate_zeroes_vector(end) = (-d3 * u_int(2));

final_target_vector = target_vector + intermediate_zeroes_vector;

% solve for the final vector u

u = d_matrix \ final_target_vector;

% add the known u terms to the final u vector

u_values = zeros(n, 1);
u_values(1) = u_int(1);
u_values(end) = u_int(2);
u_values(2:end - 1) = u;

% plot the solution

plot(x_values, u_values, 'r.');
title( 'ey3chen and mbmoin' );


