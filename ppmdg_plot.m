% Advection Test Plotting
%
% ---

clear all;
clf;
clc;
close all;

scrsz = get(0,'ScreenSize');

ncfilename = 'ppmdghy/weno2d_def_cosinebell.nc';

Q1 = nc_varget(ncfilename, 'Q1');
Q2 = nc_varget(ncfilename, 'Q2');
Q3 = nc_varget(ncfilename, 'Q3');
x1 = nc_varget(ncfilename, 'x1');
x2 = nc_varget(ncfilename, 'x2');
x3 = nc_varget(ncfilename, 'x3');
y1 = nc_varget(ncfilename, 'y1');
y2 = nc_varget(ncfilename, 'y2');
y3 = nc_varget(ncfilename, 'y3');
t = nc_varget(ncfilename, 'time');

nt = size(Q1,1);
nx = size(Q1,2);
nelem = nx/5;
elemdx = 1/nelem;

for n=1:1
    figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])%'Position', [100, 100, 2*256, 2*256 ])
    tmp = squeeze(Q3(n,:,:));
    %hold on
    hh1 = contourf(x3,y3,tmp);
    %for j=1:nelem
    %    line([j*elemdx j*elemdx],[0 1],'Marker','.','LineStyle','-')
    %end
    %hold off
    %axis([0.4 0.5 0 0.5])
    ht=title(sprintf('Time: %0.2f sec', t(n)));
end

% Get figure size
pos = get(gcf, 'Position');
width = pos(3); height = pos(4);

% Preallocate data (for storing frame data)
mov = zeros(height, width, 1, length(t), 'uint8');

% Loop through by changing XData and YData
for id = 1:length(t)
    % Update graphics data. This is more efficient than recreating plots.
    tmp = squeeze(Q2(id,:,:));
    contourf(x2,y2,tmp)
    title(sprintf('Time: %0.2f sec',t(id)));

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
imwrite(mov, map, 'dganimation.gif', 'DelayTime', 0, 'LoopCount', inf);

% ---
% Loop through ppm unlimited data
% ---

ncfilename = 'pfctnon/weno2d_def_cosinebell.nc';

Q1 = nc_varget(ncfilename, 'Q1');
Q2 = nc_varget(ncfilename, 'Q2');
Q3 = nc_varget(ncfilename, 'Q3');
x1 = nc_varget(ncfilename, 'x1');
x2 = nc_varget(ncfilename, 'x2');
x3 = nc_varget(ncfilename, 'x3');
y1 = nc_varget(ncfilename, 'y1');
y2 = nc_varget(ncfilename, 'y2');
y3 = nc_varget(ncfilename, 'y3');
t = nc_varget(ncfilename, 'time');

for n=1:1
    figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])%'Position', [100, 100, 2*256, 2*256 ])
    tmp = squeeze(Q3(n,:,:));
    %hold on
    hh1 = contourf(x3,y3,tmp);
    %for j=1:nelem
    %    line([j*elemdx j*elemdx],[0 1],'Marker','.','LineStyle','-')
    %end
    %hold off
    %axis([0.4 0.5 0 0.5])
    ht=title(sprintf('Time: %0.2f sec', t(n)));
end

% Get figure size
pos = get(gcf, 'Position');
width = pos(3); height = pos(4);

% Preallocate data (for storing frame data)
mov = zeros(height, width, 1, length(t), 'uint8');

% Loop through by changing XData and YData
for id = 1:length(t)
    % Update graphics data. This is more efficient than recreating plots.
    tmp = squeeze(Q2(id,:,:));
    contourf(x2,y2,tmp)
    title(sprintf('Time: %0.2f sec',t(id)));

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
imwrite(mov, map, 'ppmanimation.gif', 'DelayTime', 0, 'LoopCount', inf);


ncfilename = 'pfctnon/weno2d_adv_sine.nc';

close all
