function img_norm = make_isoluminant(img, mu, sigma);
% takes the image img and normalizes it to have mean mu and std sigma

img_norm = normalize_img(img);

img_norm = img_norm*sigma+mu;