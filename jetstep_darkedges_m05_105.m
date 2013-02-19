function [cmap] = jetstep_darkedges_m05_105()

M = 10;

cc = colormap(jet); % get jet colormap
N = length(cc);
ind = round([1:2:2*M-1]/(2*M)*N); 
cclist = cc(ind,:);   % choose six colors from colormap, leaving off ends

cmap = [ones(8,1)*[0.3 0.3 0.3];...   %-0.05:-0.01
        ones(3,1)*[1 1 1]];   %-0.00:-0.00

M = size(cclist,1);
cmap = [cmap; ones(19,1)*cclist(1,:)];
for n = 2:M-1
   cmap = [cmap; ones(20,1)*cclist(n,:)];
end

cmap = [cmap; 
        ones(24,1)*cclist(M,:);     % 0.90:1.00
        ones(8,1)*[0.6 0.6 0.6]];   % 1.01:1.05


