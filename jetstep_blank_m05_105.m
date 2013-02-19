function [cmap] = jetstep_blank_m05_105()

M = 10;

cc = colormap(jet); % get jet colormap
N = length(cc);
ind = round([1:2:2*M-1]/(2*M)*N); 
cclist = cc(ind,:);   % choose six colors from colormap, leaving off ends

cmap = ones(6,1)*[0.9 0.9 0.9];   %-0.05:-0.00

M = size(cclist,1);
cmap = [cmap; ones(9,1)*cclist(1,:)];
for n = 2:M-1
   cmap = [cmap; ones(10,1)*cclist(n,:)];
end

cmap = [cmap; 
        ones(11,1)*cclist(M,:);     % 0.90:1.00
        ones(5,1)*[0.9 0.9 0.9]];   % 1.01:1.05


