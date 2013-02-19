function [cmap] = blankstep()

N = 200; N2 = (N-102)/2;
cmap = [0.9*ones(N2,1); ones(N-2*N2,1); 0.7*ones(N2,1)]*[1 1 1];

