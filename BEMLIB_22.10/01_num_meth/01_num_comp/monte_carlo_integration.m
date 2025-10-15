clear all

m = 4;
alpha = 0.3;
N = 2^20;
sqrt(pi/2) / (2*alpha)^(m+0.5) * mean(randn(1,N).^(2*m))
exact =  1*3*5*7/2/(2*alpha)^m * sqrt(pi/alpha)

