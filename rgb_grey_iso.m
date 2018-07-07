%% Reads grayscale jpg files finds their average luminance and adjusts luminace in the stimulus set
clc
clear all
close all
path = pwd;

path_rgb = [path,filesep,'RGB'];
path_gs = [path,filesep,'Grayscale'];
path_gsi = [path,filesep,'Grayscale_iso'];

cd(path_rgb)

% read jpg files
img_files = dir('*.jpg');
img_files = {img_files(:).name}';

%% Read images and get intensities of non-white areas
th = 245;

mean_l = zeros(length(img_files),1);
for i = 1:length(img_files)
    
    img_org = rgb2gray(imread(img_files{i})); % load image
    %[~,name,ext] = fileparts(img_files{i});
    imwrite(img_org,fullfile(path_gs,img_files{i}));
    
    pic = img_org(img_org <= th);
    mean_l(i) = mean(pic);
    
end

mean_l_tot = mean(mean_l);
mean_l_d = mean_l - mean_l_tot;

%% Apply difference to intensity mean to images

for i = 1:length(img_files)
    
    img_org = rgb2gray(imread(img_files{i})); % load image
    %[~,name,ext] = fileparts(img_files{i});
    
    [x,y] = find(img_org<= th);
    
    for j = 1:length(x)
        img_org(x(j),y(j)) = img_org(x(j),y(j)) - mean_l_d(i);
    end
    
    img_org(img_org<0) = 0;
    img_org(img_org>255) = 255;

    imwrite(img_org,fullfile(path_gsi,img_files{i}));
    
end

cd(path)
