% ####################################################################
% # Author: Shan Qu, Delft University of Technology                ###
% # First created: May 2019                                        ###
% # This code is propriatary under the Delphi Research Consortium  ###
% #                                                                ###
% # product: ML workflow for event detection                       ###
% # generate training data and label
% ####################################################################

clear;clc;close;
%%
% load mic_syn.mat
load Data/data0_syn1.mat

dt = 0.001;
T = 3.1;
nt = size(d, 1);
nx = size(d, 2);
dx = 7.5;
x = 0 : dx : (nx-1)*dx;

%% interpolation
dt_fine = 0.0005;
t = 0 : dt : T;
t_fine = 0 : dt_fine : T;
nt_fine = length(t_fine);
d_fine = zeros(nt_fine, nx);
for ix = 1 : nx
    d_fine(:, ix) = interp1(t, d(:, ix), t_fine, 'spline');
end

dt = dt_fine;
d = d_fine;
t = t_fine;
nt = nt_fine;
    
%% add random noise
d_noise = awgn(d, 0.05,'measured','linear');

figure;
imagesc(d_noise);
% figure;
% imagesc(d, [-0.2 0.2]);

%% labeling data via mask
% d_mask = d;
% figure;
% imagesc(d_mask, [-0.1 0.1]);
% alpha(100);
% n=10;
% while strcmp(input('More?','s'),'y')
%     sig = input('D=');
%     d_mask(roipoly) = sig;
% %     imcontour(dModel);
% %     n=n+1;
% %     filenm=['den_Merger_' num2str(n) '.dat'];
% %     save(filenm,'dModel','-ascii','-double');
% end
% 
% figure;
% imagesc(d_mask);

%% labeling data via threshold for synthetic data
d_mask = zeros(size(d));
for ix = 1 : nx
    for it = 1 : nt
        if abs(d(it, ix)) < 0.05
            d_mask(it, ix) = false;
        else
            d_mask(it, ix) = true;
        end
    end
end

figure;
imagesc(d_mask);

%% create segment

d_seg = 58*2;
N_seg_pertr = length(1 : d_seg : nt);
N_seg = N_seg_pertr * nx;

data = zeros(d_seg, N_seg);
label = zeros(N_seg, 1);
label_toshow = zeros(N_seg_pertr, nx);
tmp = 1 : d_seg : nt;
sp = length(tmp(end) : nt);
for i_seg = 1 : N_seg
    ix = floor((i_seg-1) / N_seg_pertr) + 1;
    i_seg_pertr = mod(i_seg-1, N_seg_pertr)+1;
    it = (i_seg_pertr-1) * d_seg + 1;
    
    if i_seg_pertr == N_seg_pertr
        data(1:sp, i_seg) = d_noise(it:end, ix);

        if ((sum(d_mask(it:end, ix))/length(d_mask(it:end, ix)) * 100) < 20)
            label_toshow(i_seg_pertr, ix) = 0;
            label(i_seg) = 0;
        else
            label_toshow(i_seg_pertr, ix) = 1;
            label(i_seg) = 1;
        end
        
    else
        data(:, i_seg) = d_noise(it:it+d_seg-1, ix);
        if ((sum(d_mask(it:it+d_seg-1, ix))/length(d_mask(it:it+d_seg-1, ix)) * 100) < 20)
            label_toshow(i_seg_pertr, ix) = 0;
            label(i_seg) = 0;
        else
            label_toshow(i_seg_pertr, ix) = 1;
            label(i_seg) = 1;
        end
    end
end


figure;imagesc(label_toshow)


label_pred = label;
N_seg_pertr = length(1 : d_seg : nt);
N_seg = N_seg_pertr * nx;

% label = zeros(N_seg, 1);
label_toshow = zeros(size(d_noise));
true_label_toshow = zeros(size(d_noise));
tmp = 1 : d_seg : nt;
sp = length(tmp(end) : nt);
for i_seg = 1 : N_seg
    ix = floor((i_seg-1) / N_seg_pertr) + 1;
    i_seg_pertr = mod(i_seg-1, N_seg_pertr)+1;
    it = (i_seg_pertr-1) * d_seg + 1;
    
    if i_seg_pertr == N_seg_pertr
        if label_pred(i_seg)
            label_toshow(it:end, ix) = 1;
        else
            label_toshow(it:end, ix) = 0;
        end
%         if label(i_seg)
%             true_label_toshow(it:end, ix) = 1;
%         else
%             true_label_toshow(it:end, ix) = 0;
%         end
        
    else
        if label_pred(i_seg)
            label_toshow(it:it+d_seg-1, ix) = 1;
        else
            label_toshow(it:it+d_seg-1, ix) = 0;
        end
%         if label(i_seg)
%             true_label_toshow(it:it+d_seg-1, ix) = 1;
%         else
%             true_label_toshow(it:it+d_seg-1, ix) = 0;
%         end
    end
end

load cgray.mat
load cmap.mat


fig=figure(1)
clf
%ha = tight_subplot(Nh, Nw, gap, marg_ver, marg_hor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha = tight_subplot(1,1,[0.01 0.01],[0.15 0.05],[0.2 0.1]);
axes(ha(1)); imagesc(x(6:end-5),t,d_noise(:, 6:end-5),  [-0.15 0.15]); colorbar('eastoutside');colormap(cgray);


axes(ha(1));
ylabel('Time [s]')
xlabel('x [m]')
set(gca,'fontsize',30); 

set(fig,'paperposition',[0 0.1 10 10])


print -depsc Fig/d1_noise_n005.eps
print -djpeg Fig/d1_noise_n005.jpg

fig=figure(1)
clf
%ha = tight_subplot(Nh, Nw, gap, marg_ver, marg_hor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha = tight_subplot(1,1,[0.01 0.01],[0.15 0.05],[0.2 0.1]);
axes(ha(1)); imagesc(x(6:end-5),t,d_noise(:, 6:end-5).*label_toshow(:, 6:end-5), [-0.15 0.15]);colorbar('eastoutside');colormap(cgray)

axes(ha(1));
ylabel('Time [s]')
xlabel('x [m]')
set(gca,'fontsize',30); 

set(fig,'paperposition',[0 0.1 10 10])


print -depsc Fig/d1_label_noise_n005.eps
print -djpeg Fig/d1_label_noise_n005.jpg


%%
save Data/d1_noise_n005.mat d_noise
save Data/data1_n005.mat data
save Data/label1_n005.mat label


        
