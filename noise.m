function [x, seed] = noise(alpha, N, seed, verbose)
	% generate 1/f^alpha noise with N data points. 
    % Written by Holger Hennig (2008)
	% For 0<alpha<2 the time series exhibits long-range correlations.
	% Input: x=noise(alpha,N);
	% Default: alpha=1, N=1000, mean 0 and standard deviation 1.
	
    % Creates 1/f^alpha noise from inverse Fast Fourier Transformation (FFT) 
    % of 1/f power spectral density (PSD).
	% Steps:
	%	1. generate PSD with s(f)=1/f^alpha for f=0,...,1.
	%	2. Randomize phases
	%       a) Generate uniformly distributed random phases phi=0,...,2*pi
	%       b) y(f) = sqrt(s)*exp(i*phi), y(1) = 0; (zero mean)
	%	3. Symmetrize the PSD. y(N-k+2))=conj(y(k)) for k=2...N
	%	4. inverse FFT of y gives time series x. Note that x is real due to
	%	   symmetrizing the PSD (step 3).
	
	if nargin < 1 || isempty(alpha), alpha = 1; end
	if nargin < 2 || isempty(N), N = 1000; end
	if ~exist('seed','var') || isempty(seed), seed = floor(1e6*rand); end
	if ~exist('verbose','var') || isempty(verbose), verbose = 1; end
	
	% n phases and n amplitudes result in 2*n real data points x(t) after FFT.
	n=ceil(N/2);
	
	% Frequencies in Interval [0, 1/2]
	f=linspace(0,1/2,n);
	
	% 1/f^alpha power spectral density
	s=1./f.^alpha;
	
	% zero mean
	s(1)=0;
	
	% random phases uniformly distributed in [0, 2pi]
	rng(seed) % seeds the random number generator so that rand and randn produce a predictable sequence of numbers
	phi = rand(1,n)*2*pi;
	amp = sqrt(s);
	y = amp.*exp(1i.*phi);
	
	% symmetrize Y(N-k+2)) = conj(Y(k)) for k=2...N. "conj" is the complex conjugate
	Y = [y,0,fliplr(conj(y(2:end)))]; % fliplr returns a vector of the same length with the order of its elements reversed
	x = ifft(Y); % inverse FFT
	x = x(1:N); %for odd N we generated N+1 data points.
	
	% set standard deviation to 1
	x = x/std(x);
	
	if verbose >0.5, fprintf('Generated 1/f^alpha noise with alpha=%2.1f\n',alpha); end
	
end

