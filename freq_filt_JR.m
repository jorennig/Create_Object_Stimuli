clc
clear all
close all

nd = 'img_filt_neu';
mkdir(nd);

%%% parameters
deg_image = 1;% angle of view subtended by the image

cpdLP = 6;% low pass cut off in cycles per degree
cpdHP = 7;%  high pass cut off in cycles per degree

cut_offLP = cpdLP/deg_image;% cycles per image
cut_offHP = cpdHP/deg_image;% cycles per image

filt_type = 'gaussian';% one in {'gaussian','ideal','btw'}

btw_order = 4;% for the Butterworth filter specify also order

filt_flag = 0;% put to 1 to visualize the filter in space and frequency

[FileName,PathName] = uigetfile('*.*','Select pictures to filter','MultiSelect', 'on' );

for i=1:max(1,iscell(FileName)*length(FileName))
    
    % read image
    if ~iscell(FileName)
        img = imread([PathName FileName]);
    else
        img = imread([PathName FileName{i}]);
    end
    
    if numel(size(img))>2; % rgb image, convert to gray
        img = rgb2gray(img);
    end
    
%     % Bild verkleinern
%     fk = 20;
%     % left right
%     img(1:fk,:)=[];
%     img((size(img,1)+1)-fk:size(img,1),:)=[];
%     
%     % up down
%     img(:,1:fk)=[];
%     img(:,(size(img,2)+1)-fk:size(img,2))=[];
    
    % determine padding size and cutoff based on image size
    PQ = paddedsize(size(img));
    D0LP = cut_offLP/100*PQ(1);%size(img,1);% cutoff adapted to img size
    D0HP = cut_offHP/100*PQ(1);%size(img,1);% cutoff adapted to img size
    
    %%% low-pass/high-pass filters, designed in frequency! %%%%
    
    LP = lpfilter(filt_type, PQ(1), PQ(2), D0LP, btw_order); % Calculate the LPF
    HP = hpfilter(filt_type, PQ(1), PQ(2), D0HP, btw_order);
    
    
    if (filt_flag)% show the filter only the first time
        
        
        filt_spLP = real(fftshift(ifft2(LP)));% antitransform to get the filter in space
        filt_spHP = real(fftshift(ifft2(HP)));% antitransform to get the filter in space
        
        figure,
        subplot(1,2,1),
        imagesc(filt_spLP),title('low-pass filter in space')
        subplot(1,2,2),
        %freqz2(filt_spLP,size(filt_spLP)), title('filter response in frequency')
        imagesc(abs(fftshift(LP)))
        
        figure,
        subplot(1,2,1),
        imagesc(filt_spHP),title('high-pass filter in space')
        subplot(1,2,2),
        %freqz2(filt_spHP,size(filt_spHP)), title('filter response in frequency')
        imagesc(abs(fftshift(HP)))
        pause
        close
    end
    
    % Fourier transform the image
    F=fft2(double(img),size(LP,1),size(LP,2)); % Calculate the discrete Fourier
    
    %%%%  frequency filtering -> element-wise product of filter and img fourier
    % transformed
    
    img_LP=real(ifft2(LP.*F)); % multiply the Fourier spectrum by the LPF and apply the inverse, discrete Fourier transform
    
    img_HP=real(ifft2(HP.*F));
    
    img_LP=img_LP(1:size(img,1), 1:size(img,2)); % Resize the image to undo padding
    img_HP=img_HP(1:size(img,1), 1:size(img,2)); % Resize the image to undo padding
    
    % plot images and Fourier Spectra
    % transform image and filtered images to be iso-luminant
    mu = mean(img_LP(:));
    sigma = std(img_HP(:));
    
    img_iso = make_isoluminant(img,mu,sigma);
    img_LPiso = make_isoluminant(img_LP,mu,sigma);
    img_HPiso = make_isoluminant(img_HP,mu,sigma);
    
    %Rand entfernen
    cv = 20; % get color from the defined area
    val_iso = mean2(img_iso(cv:cv*2,cv:cv*2));
    val_LPiso = mean2(img_LPiso(cv:cv*2,cv:cv*2));
    val_HPiso = mean2(img_HPiso(cv:cv*2,cv:cv*2));

    pc = 15; % frame size px
    
%     % Rand vertikal entfernen
%     for p = 1:size(img_iso,1)
%         for q = 1:pc
%         
%             img_LPiso(p,q) = val_LPiso;
%             img_LPiso(p,(size(img_iso,2)+1)-q) = val_LPiso;
% 
%             img_HPiso(p,q) = val_HPiso;
%             img_HPiso(p,(size(img_iso,2)+1)-q) = val_HPiso;
%             
%         end
%     end
%     
%     % Rand horizontal entfernen
%     for p = 1:pc
%         for q = 1:size(img_iso,2)
%         
%             img_LPiso(p,q) = val_LPiso;
%             img_LPiso((size(img_iso,1)+1)-p,q) = val_LPiso;
%             
%             img_HPiso(p,q) = val_HPiso;
%             img_HPiso((size(img_iso,1)+1)-p,q) = val_HPiso;
% 
%         end
%     end    
%     
    % Images darker
    lf = 40; % luminance factor
    
    img_iso = img_iso-lf;
    img_LPiso = img_LPiso-lf;
    img_HPiso = img_HPiso-lf;
    
    % Plot & save
    name = FileName{i};
    name = name(1:end-4);
    
    lt = 0;
    ht = 255;
    
    res = 120; % resolution in dpi
    path = [PathName nd];
    
    %% White/black background
%     img_iso_b = img_iso;
%     img_iso_w = img_iso;
%     
%     ref_grey = img_iso(1,1);
%     
%     for p = 1:size(img_iso,1)
%         for q = 1:size(img_iso,2)
%             
%             if img_iso_b(p,q) == ref_grey;
%                 
%                 img_iso_b(p,q) = 0;
%                 
%             end
%         
%         end    
%     end
    
    % Save filtered images
    % IMG Original    
    figure(1)
    f = imshow(img_iso,[lt ht]);
    
    filename = ['BB_' name];    
    export_fig(gcf,fullfile(path,filename),'-jpg',['-r',num2str(res)]);
    
    close all    
    
    % LP
    figure(2)
    g = imshow(img_LPiso,[lt ht]);
    
    filename = ['LP_' name];    
    export_fig(gcf,fullfile(path,filename),'-jpg',['-r',num2str(res)]);
    
    close all    
    
    % HP
    figure(3)
    h = imshow(img_HPiso,[lt ht]);
    
    filename = ['HP_' name];    
    export_fig(gcf,fullfile(path,filename),'-jpg',['-r',num2str(res)]);
    
    close all
    
%     set(gca,'LooseInset',get(gca,'TightInset')) % ohne weiﬂen Rand
%     
%     f = gcf; % f is the handle of the figure you want to export
%     figpos = getpixelposition(f);
%     resolution = get(0,'ScreenPixelsPerInch');
%     set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
%     path = [PathName nd];
%     name = ['HP_' FileName{i}];
%     print(f,fullfile(path,name),'-dbmp',['-r',num2str(rez)],'-opengl') % save file 
        
%     figure(2),
%     subplot(1,3,1), imshow(imgiso,[0 255]), title('original image')
%     subplot(1,3,2), imshow(img_LPiso,[0 255]), title('low-pass filtered image')
%     subplot(1,3,3), imshow(img_HPiso,[0 255]), title('high-pass filtered image')
%     
%     
%     
%     % spectra
%     S1 = log(1+abs(fftshift(F)));
%     S2 = log(1+abs(fftshift(fft2(img_LP,size(LP,1),size(LP,2))))); % use abs to compute the magnitude (handling imaginary) and use log to brighten display
%     S3 = log(1+abs(fftshift(fft2(img_HP,size(LP,1),size(LP,2)))));
%     
%     figure(3),
%     subplot(1,3,1), imagesc(S1), title('Freq. spectrum original img')
%     axis('square')
%     subplot(1,3,2), imagesc(S2),  title('Freq. spectrum lp img')
%     axis('square')
%     subplot(1,3,3), imagesc(S3),  title('Freq. spectrum hp img')
%     axis('square')
%     
%     pause
%     close all
    
end