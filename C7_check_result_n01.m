% ####################################################################
% # Author: Shan Qu, Delft University of Technology                ###
% # First created: May 2019                                        ###
% # This code is propriatary under the Delphi Research Consortium  ###
% #                                                                ###
% # product: ML workflow for event detection                       ###
% # check predicted results
% ####################################################################


clear;clc;close
%%
load Data/d2_noise_n01.mat

dt = 0.0005;
T = 3.1;
nt = size(d_noise, 1);
nx = size(d_noise, 2);
dx = 7.5;
x = 0 : dx : (nx-1)*dx;
t = 0 : dt : T;

figure;
imagesc(d_noise, [-0.2 0.2]);

%%
load Data/label_pred2_n01.mat
% load label2_20.mat

d_seg = 58 * 2.0
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


print -depsc Fig/d2_noise_n01.eps
print -djpeg Fig/d2_noise_n01.jpg


fig=figure
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


print -depsc Fig/d2_pred_noise_n01.eps
print -djpeg Fig/d2_pred_noise_n01.jpg

% figure;
% plot(d_noise(:, 150));
% view([90 90])
% figure;
% plot(d(:, 150))
% view([90 90])


%%

% d_new = d_noise .* label_toshow;
% load Data/d.mat
% 
% fig=figure(1)
% ha = tight_subplot(1,1,[0.01 0.01],[0.05 0.05],[0.01 0.05]);
% 
% axes(ha(1));
% subplot(1,5,1);
% a1=plot(t,d_noise(:, 21),'LineWidth',3); 
% hold on
% a2=plot(t,d(:, 21),'LineWidth',3); 
% hold on
% a3=plot(t,label_toshow(:, 21),'LineWidth',3); 
% view([90 90])
% set(gca,'YLim',[-1.2 1.2]);
% set(gca,'XLim',[0 3])
% xlabel('Time [s]')
% ylabel('x [m]')
% set(gca,'fontsize',30); 
% 
% subplot(1,5,2);
% a1=plot(t,d_noise(:, 71),'LineWidth',3); 
% hold on
% a2=plot(t,d(:, 71),'LineWidth',3); 
% hold on
% a3=plot(t,label_toshow(:, 71),'LineWidth',3); 
% view([90 90])
% set(gca,'YLim',[-1.2 1.2]);
% set(gca,'XLim',[0 3])
% set(gca,'xtick',[])
% ylabel('x [m]')
% set(gca,'fontsize',30); 
% 
% subplot(1,5,3);
% a1=plot(t,d_noise(:, 121),'LineWidth',3); 
% hold on
% a2=plot(t,d(:, 121),'LineWidth',3); 
% hold on
% a3=plot(t,label_toshow(:, 121),'LineWidth', 3); 
% view([90 90])
% set(gca,'YLim',[-1.2 1.2]);
% set(gca,'XLim',[0 3])
% set(gca,'xtick',[])
% ylabel('x [m]')
% set(gca,'fontsize',30); 
% 
% subplot(1,5,4);
% a1=plot(t,d_noise(:, 171),'LineWidth', 3); 
% hold on
% a2=plot(t,d(:, 171),'LineWidth', 3); 
% hold on
% a3=plot(t,label_toshow(:, 171),'LineWidth', 3); 
% view([90 90])
% set(gca,'YLim',[-1.2 1.2]);
% set(gca,'XLim',[0 3])
% set(gca,'xtick',[])
% ylabel('x [m]')
% set(gca,'fontsize',30); 
% 
% subplot(1,5,5);
% a1=plot(t,d_noise(:, 221),'LineWidth', 3); 
% hold on
% a2=plot(t,d(:, 221),'LineWidth', 3); 
% hold on
% a3=plot(t,label_toshow(:, 221),'LineWidth', 3); 
% view([90 90])
% set(gca,'YLim',[-1.2 1.2]);
% set(gca,'XLim',[0 3])
% set(gca,'xtick',[])
% ylabel('x [m]')
% set(gca,'fontsize',30); 
% 
% set(fig,'paperposition',[0 0.1 20 20])
% 
% print -depsc Fig/profile.eps
% print -djpeg Fig/profile.jpg


%%

load Data/label_pred2_n01_no2d.mat
% load label2_20.mat

d_seg = 58 * 2.0
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



% fig=figure(1)
% clf
% %ha = tight_subplot(Nh, Nw, gap, marg_ver, marg_hor)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ha = tight_subplot(1,1,[0.01 0.01],[0.15 0.05],[0.2 0.1]);
% axes(ha(1)); imagesc(x(6:end-5),t,d_noise(:, 6:end-5),  [-0.15 0.15]); colorbar('eastoutside');colormap(cgray);
% 
% 
% axes(ha(1));
% ylabel('Time [s]')
% xlabel('x [m]')
% set(gca,'fontsize',30); 
% 
% set(fig,'paperposition',[0 0.1 10 10])
% 
% 
% print -depsc Fig/d2_noise_n12_no2d.eps
% print -djpeg Fig/d2_noise_n12_no2d.jpg


fig=figure
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


print -depsc Fig/d2_pred_noise_n01_no2d.eps
print -djpeg Fig/d2_pred_noise_n01_no2d.jpg

%%
% load Data/label_pred2_rf_n01.mat
% % load label2_20.mat
% 
% d_seg = 58 * 2.0
% N_seg_pertr = length(1 : d_seg : nt);
% N_seg = N_seg_pertr * nx;
% 
% % label = zeros(N_seg, 1);
% label_toshow = zeros(size(d_noise));
% true_label_toshow = zeros(size(d_noise));
% tmp = 1 : d_seg : nt;
% sp = length(tmp(end) : nt);
% for i_seg = 1 : N_seg
%     ix = floor((i_seg-1) / N_seg_pertr) + 1;
%     i_seg_pertr = mod(i_seg-1, N_seg_pertr)+1;
%     it = (i_seg_pertr-1) * d_seg + 1;
%     
%     if i_seg_pertr == N_seg_pertr
%         if label_pred(i_seg)
%             label_toshow(it:end, ix) = 1;
%         else
%             label_toshow(it:end, ix) = 0;
%         end
% %         if label(i_seg)
% %             true_label_toshow(it:end, ix) = 1;
% %         else
% %             true_label_toshow(it:end, ix) = 0;
% %         end
%         
%     else
%         if label_pred(i_seg)
%             label_toshow(it:it+d_seg-1, ix) = 1;
%         else
%             label_toshow(it:it+d_seg-1, ix) = 0;
%         end
% %         if label(i_seg)
% %             true_label_toshow(it:it+d_seg-1, ix) = 1;
% %         else
% %             true_label_toshow(it:it+d_seg-1, ix) = 0;
% %         end
%     end
% end
% 
% load cgray.mat
% load cmap.mat
% 
% 
% fig=figure
% clf
% %ha = tight_subplot(Nh, Nw, gap, marg_ver, marg_hor)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ha = tight_subplot(1,1,[0.01 0.01],[0.15 0.05],[0.2 0.1]);
% axes(ha(1)); imagesc(x(6:end-5),t,d_noise(:, 6:end-5).*label_toshow(:, 6:end-5), [-0.15 0.15]);colorbar('eastoutside');colormap(cgray)
% 
% 
% axes(ha(1));
% ylabel('Time [s]')
% xlabel('x [m]')
% set(gca,'fontsize',30); 
% 
% set(fig,'paperposition',[0 0.1 10 10])
% 
% 
% print -depsc Fig/d2_pred_rf_noise_n01.eps
% print -djpeg Fig/d2_pred_rf_noise_n01.jpg

% figure;
% plot(d_noise(:, 150));
% view([90 90])
% figure;
% plot(d(:, 150))
% view([90 90])

