function [ x, y, seed ] = arfima_model( d1, d2, W, len, seed )
	% generates two-component ARFIMA process. 
    % Written by H. Hennig (2013)
	% Based on article by Podobnik et al., Physica A, 387, 3954 (2008)
    % Input example:
    %   [x,y] = arfima(d1, d2, W);
    % Parameters:
    %   d1: scaling parameter for autocorrelations of X. fractal (Hurst) exponent H1 = 0.5+d1
	%   d2: scaling parameter for autocorrelations of Y. fractal (Hurst) exponent H2 = 0.5+d2
	%   W:  coupling between time series X and Y. 0.5<W<1.
	%       W=0.5: max coupling. W=1: no coupling
	% For 1/f^alpha noise we find alpha = 2*H-1 where H is the Hurst exponent
	
    % check input
	if ~exist('d1','var') || isempty(d1), d1 = 0.4; end
	if ~exist('d2','var') || isempty(d2), d2 = 0.4; end
	if ~exist('W','var')  || isempty(W), W = 0.5; end
	if ~exist('len','var') || isempty(len), len = 1e3; end
	if ~exist('seed','var') || isempty(seed), seed = floor(1e6*rand); end
	
	% initialize
	weights_x = zeros(len,1);
	weights_y = zeros(len,1);
	zx = zeros(len,1);
	zy = zeros(len,1);
	
	% independent and identically distributed (i.i.d.) Gaussian noise
	rng(seed) % seeds the random number generator so that rand and randn produce a predictable sequence of numbers
	r = randn(2*len,1);
	e_x = r(1:len);
	e_y = r(len+1:end);
	
	% calc weights
	for n = 1:len
		weights_x(n) = weights(n,d1);
		weights_y(n) = weights(n,d2);
	end
	
	% first point in time series: no history yet to consider
	x(1) = e_x(1);
	y(1) = e_y(1);
	
	zx(1) = e_x(1);
	fprintf('Generated two component ARFIMA process with Hurst exponent H1=%1.2f, H2=%1.2f\n',d1+0.5,d2+0.5);
	
    % Two/component ARFIMA, see definition in [1]	
	for t = 2:len
       		
        % calc stochastic variables zx and zy
        zx(t-1) = W*x(t-1) + (1-W)*y(t-1) + e_x(t);
		
        % calc x and y, i.e., weighted sum over all past values (to generate memory)
        zy(t-1) = (1-W)*x(t-1) + W*y(t-1) + e_y(t);

        % DCCA exponent delta = (H1+H2)/2
        x(t) = weights_x(1:t-1)' * zx(t-1:-1:1);
		y(t) = weights_y(1:t-1)' * zy(t-1:-1:1);        		
	end
	
	
	
	
	%--------------------------------------------------------------
	% weights a_n(d)
	
function [a]= weights(n,d)
	
	%a = d * gamma(n-d) / ( gamma(1-d) * gamma(n+1) );
	loga = log(d) + gammaln(n-d) - gammaln(1-d) - gammaln(n+1);
	a = exp(loga);
