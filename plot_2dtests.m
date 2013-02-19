clear all

which_test = {'adv_sine', ... % uniform advection of  sin^4
              'def_smooth_cosinebell', ... % Leveque deformation test w/smoother init cond
              'fdefvel_smooth_cosbell', ... % Blossey-Durran deformation w/smoother init cond
              'fdef_sqwave'}; % Blossey-Durran deformation on a discontinuous initial cond.
Ntest = length(which_test);
for mm = 1:1
  ncfilename = ['weno2d_' which_test{mm} '.nc'];
  
  which_res = {'1','2','3'};
  Nres = length(which_res);
  for nn = 1:Nres
    resname = which_res{nn};

    figure(10+mm*Nres+nn); clf

    for n = 1:4
      if n==1
        name = 'PPM, No limiting';
        nc = ['pfctnon/' ncfilename];
        xm = 0.285; ym = 0.71; subfiglabel = 'a';
      elseif n==2
        name = 'PPM, Limiting for non-negativity';
        nc = ['pfctpos/' ncfilename];
        xm = 0.51; ym = 0.71; subfiglabel = 'b';
      elseif n==3
        name = 'PPM, Selective (FCT)';
        nc = ['pfctsel/' ncfilename];
        xm = 0.285; ym = 0.41; subfiglabel = 'c';
      elseif n==4
        name = 'PPM, Selective, Positive (PMod)';
        nc = ['ppmpscl/' ncfilename];
        xm = 0.51; ym = 0.41; subfiglabel = 'd';
%      elseif n==5
%        name = 'PPMDG Hybrid';
%        nc = ['ppmdg/' ncfilename]; %['ppmdg/' ncfilename];
%        xm = 0.285; ym = 0.71; subfiglabel = 'e';
      end
      plotfinal2d(name,nc,resname,xm,ym,subfiglabel);
    end
    
    % $$$     print -depsc2 figure_2d_fdef_sqwave.eps
    eval(['print -dpng -r200 figure_2d_' which_test{mm} ...
          '_res' which_res{nn} '.png']);
      
    % Plot second set of figures  
    figure(20+mm*Nres+nn); clf

    for n = 5:5
      if n==5
        name = 'PPMDG Hybrid';
        nc = ['ppmdg/' ncfilename]; %['ppmdg/' ncfilename];
        xm = 0.285; ym = 0.71; subfiglabel = 'e';
      end
      plotfinal2d(name,nc,resname,xm,ym,subfiglabel);
    end
    
    

% $$$     print -depsc2 figure_2d_fdef_sqwave.eps
    eval(['print -dpng -r200 sec_figure_2d_' which_test{mm} ...
          '_res' which_res{nn} '.png']);
  end
  
end

close all