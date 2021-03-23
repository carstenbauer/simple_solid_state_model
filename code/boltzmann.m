 function out=boltzmann(coeff,X,Y)
 lambda = coeff(1);
 N = 1-exp(-lambda*X);
 diff = N - Y; 
 square_diff = diff.^2;
 out = sum(square_diff);