% ####################################################################
% # Author: Shan Qu, Delft University of Technology                ###
% # First created: May 2019                                        ###
% # This code is propriatary under the Delphi Research Consortium  ###
% #                                                                ###
% # product: ML workflow for event detection                       ###
% # generate 2D features
% ####################################################################

clear;clc;close;
%%
load Data/d2_noise_n01.mat

dt = 0.0005;
T = 3.1;
nt = size(d_noise, 1);
nx = size(d_noise, 2);
dx = 7.5;
x = 0 : dx : (nx-1)*dx;
t = 0 : dt : T;
%%

% figure;
% imagesc(d_noise, [-0.7 0.7]);
% figure;
% imagesc(d, [-0.2 0.2]);


d_noise_grayscale = mat2gray(d_noise, 0.5*[-0.5 0.5]);
%%
d_seg = 58*2;
N_seg_pertr = length(1 : d_seg : nt);
N_seg = N_seg_pertr * nx;

feature_glcm = zeros(N_seg, 128);
feature_glcm_toshow = zeros(N_seg_pertr, nx, 128);
for i_seg = 1 : N_seg
    ix = floor((i_seg-1) / N_seg_pertr) + 1;
    i_seg_pertr = mod(i_seg-1, N_seg_pertr)+1;
    it = (i_seg_pertr-1) * d_seg + 1;
    
    if (ix >= 9) && (ix <= nx-8)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-8:ix+8);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 8; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; -7 7; -8 8; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; -7 -7; -8 -8; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; -7 0; -8 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-8:ix+8);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 8; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; -7 7; -8 8; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; -7 -7; -8 -8; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; -7 0; -8 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end
    
    if (ix == 8) || (ix == nx-7)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-7:ix+7);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; -7 7; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; -7 -7; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; -7 0; 0 0;],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 8 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-7:ix+7);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; -7 7; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; -7 -7; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; -7 0; 0 0;],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 8 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end

    if (ix == 7) || (ix == nx-6)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-6:ix+6);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 7 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-6:ix+6);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; -6 6; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; -6 -6; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; -6 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 7 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end

    if (ix == 6) || (ix == nx-5)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-5:ix+5);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 6 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-5:ix+5);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 5; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; -5 5; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; -5 -5; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; -5 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 6 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end
    
    if (ix == 5) || (ix == nx-4)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-4:ix+4);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 5 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-4:ix+4);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 4; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; -4 4; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; -4 -4; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; -4 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 5 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end
    
    if (ix == 4) || (ix == nx-3)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-3:ix+3);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; 0 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 4 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end     
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-3:ix+3);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; -3 3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; -3 -3; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; -3 0; 0 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 4 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end    
    if (ix == 3) || (ix == nx-2)
        if i_seg_pertr == N_seg_pertr
%             data(1:sp, i_seg) = d_noise(it:end, ix);
            tmp = d_noise_grayscale(it:end, ix-2:ix+2);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 3 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        else
            tmp = d_noise_grayscale(it:it+d_seg-1, ix-2:ix+2);
            glcm = graycomatrix(tmp,'Offset',[ 0 1; 0 2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 1; -2 2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 -1; -2 -2; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; ...
                -1 0; -2 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],'Symmetric', true);
            stats = graycoprops(glcm);
            for i = 1 : 4
                for j = 3 : 8
                    stats.Contrast((i-1)*8 + j) = 0.0;
                    stats.Correlation((i-1)*8 + j) = 0.0;
                    stats.Energy((i-1)*8 + j) = 0.0;
                    stats.Homogeneity((i-1)*8 + j) = 0.0;
                end
            end            
            for i = 1 : 32
                feature_glcm(i_seg, 1 + i-1) = stats.Contrast(i);
                feature_glcm(i_seg, 32+1 + i-1) = stats.Correlation(i);
                feature_glcm(i_seg, 32+32+1 + i-1) = stats.Energy(i);
                feature_glcm(i_seg, 32+32+32+1 + i-1) = stats.Homogeneity(i);
            end
            
            feature_glcm_toshow(i_seg_pertr, ix, :) = feature_glcm(i_seg, :);
        end
    end  
    
end

figure;
i = 2;
subplot(4,4,1);
imagesc(squeeze(feature_glcm_toshow(:,:,1 + i-1)));
subplot(4,4,2);
imagesc(squeeze(feature_glcm_toshow(:,:,32+1 + i-1)));
subplot(4,4,3);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));
subplot(4,4,4);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));

i = i + 8;
subplot(4,4,5);
imagesc(squeeze(feature_glcm_toshow(:,:,1 + i-1)));
subplot(4,4,6);
imagesc(squeeze(feature_glcm_toshow(:,:,32+1 + i-1)));
subplot(4,4,7);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));
subplot(4,4,8);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));

i = i + 8;
subplot(4,4,9);
imagesc(squeeze(feature_glcm_toshow(:,:,1 + i-1)));
subplot(4,4,10);
imagesc(squeeze(feature_glcm_toshow(:,:,32+1 + i-1)));
subplot(4,4,11);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));
subplot(4,4,12);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));

i = i + 8;
subplot(4,4,13);
imagesc(squeeze(feature_glcm_toshow(:,:,1 + i-1)));
subplot(4,4,14);
imagesc(squeeze(feature_glcm_toshow(:,:,32+1 + i-1)));
subplot(4,4,15);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));
subplot(4,4,16);
imagesc(squeeze(feature_glcm_toshow(:,:,32+32+1 + i-1)));




save Data/testdata_2dfeature_n01.mat feature_glcm

















        
