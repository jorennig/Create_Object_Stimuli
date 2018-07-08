function std_img = normalize_img(img);
% normalize the img to have

img = double(img);% convert in case of uint8 values
std_img = (img-mean(img(:)))/std(img(:));