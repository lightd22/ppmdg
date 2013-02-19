% Plot Advection Tests using plot_2dadv.m
% By: Devin Light
% ------

clear all;
close all;
clc;

tests = {'adv_sine', ... % Uniform adv of sine^4
         'def_cosinebell', ... % LeVeque deformation test cosinebell
         'def_smooth_cosinebell', ... % Smoother version of LeVeque test
         };
res = {'1','2','3'};

which_test = tests(1);
which_res = res(1);


ncfilename = strcat('weno2d_' ,which_test{1}, '.nc');
for n=1:2
    if n==1
        methname = 'PPM, No Limiting';
        nc = ['pfctnon/' ncfilename];
        file = ['figures/ppm' which_test{1}];
        [Q1,x1,y1,t1] = plot_2dadv(methname,nc,which_res,file);
    elseif n==2
        methname = 'PPMDG, No Limiting';
        nc = ['ppmdghy/', ncfilename];
        file = ['figures/ppmdg' which_test{1}];
        [Q2,x2,y2,t2] = plot_2dadv(methname,nc,which_res,file);
    end
    
end

% Make combined animation
% ---


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
subplot(1,2,1)
%caxis manual; % allow subsequent plots to use the same color limits
%caxis([-0.1 1]); % set the color axis scaling to your min and max color limits

subplot(1,2,2)
set(gca, 'nextplot', 'replacechildren'); 
%caxis manual; % allow subsequent plots to use the same color limits
%caxis([-0.1 1]); % set the color axis scaling to your min and max color limits
pos = get(gcf, 'Position');
width = pos(3); height = pos(4);
mov = zeros(height, width, 1, length(t1), 'uint8');

for id = 1:length(t1)
        subplot(1,2,1)
        tmp1 = squeeze(Q1(id,:,:));
        contourf(x1,y1,tmp1)
        colorbar('location','EastOutside')
        ftitle = sprintf('PPM, No Limit - Time: %0.2f sec',t1(id));
        title(ftitle);
        
        subplot(1,2,2)
        tmp2 = squeeze(Q2(id,:,:));
        contourf(x2,y2,tmp2)
        colorbar('location','EastOutside')
        ftitle = sprintf('PPMDG - Time: %0.2f sec',t2(id));
        title(ftitle);

        % Get frame as an image
        f = getframe(gcf);

        % Create a colormap for the first frame. For the rest of the frames,
        % use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
        else
            mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
        end
end

% Create animated GIF
imwrite(mov, map, 'combo.gif', 'DelayTime', 0, 'LoopCount', inf);

