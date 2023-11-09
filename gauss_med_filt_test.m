function [denoised_img] = gauss_med_filt_test(img, n, B)
%%%%%%%FFT2 of Image%%%%%%%%%%%
    figure('name','Noisy Image')
    imshow(img)
    fft_img=fft2(double(img))/size(img,1)/size(img,2);
    shifted_img =fftshift(fft_img); %power spectrum
 %   figure('name','Power Spectrum')
 %   imshow(log((255 * real(shifted_img))), []); %visualization of shifted image
    C = shifted_img;
 %   figure('name','Magnitude of Noisy Image')
 %   imagesc(abs(C))
    %%%%%%%%%%% Noisy Image Surface Visualization%%%%%%%%%%%%
 %   C_viz = log(C); %log transform for added variance + visualization
 %   G_size = size(C_viz);
 %   G_freq = 75; % frequency-filter width threshold
 %   [X,Y] = meshgrid(1:G_size(1),1:G_size(2));
 %   i = min(X-1,G_size(1)-X+1);
 %   j = min(Y-1,G_size(2)-Y+1);
 %   H = exp(-.5*(i.^2+j.^2)/G_freq^2);
 %   G = real(ifft2(H.*fft2(C_viz)));
 %   figure('name','Noisy Image Fourier Surface Visualization')
 %   surf(X,Y,G,'edgecolor','none');
 %   colorbar;
 %   colormap default;
    %%%%%%%Surface Params%%%%%%%%%%%
    edge = floor(n/2); %this is how you get the floor div.
    e = 1; %scaling coefficient
    g = zeros(n);
    center = ceil(n/2);
    k1 = floor((n-1)/2);
    %%%%%%%create surface%%%%%%%%%%%
    for u = -edge:edge % this is floor division of -n/2
        for v = -edge:edge% this is floor division of n/2
            g(u+center,v+center) = 1-e * exp(-B*((u^2+k1^2) + (v^2+k1^2)));
        end
    end
    figure('name','Gaussian Surface')
    mesh(g)
    %%%%%%%%%%%%%%%%%%%%Pixel Iteration%%%%%%%%%%%%%%%%%%%
%    [row,col] = size(C);
%    C_padded = mirroring2(C,n);
%    thold = 9; %trying 3, 8. 6 causes black stripes. Which means the freq. peaks are replaced with zeros.
%    gauss_filt_img = zeros(row,col);
%    iter = n-1;
%    for i = 1:row
%        for j = 1:col
%            target = C_padded(i:i+iter, j:j+iter);
%            med = median(abs(target), "all");
%            %test_val1 = abs(C_padded(i,j));
%            test_val = abs(target(ceil(numel(target)/2)));
%            comp = test_val/med;
%            if comp <= thold
%                gauss_filt_img(i,j) = C_padded(i,j);
%            elseif comp >= thold
%                new_target = target.*g; 
%                gauss_filt_img(i:i+(n-1), j:j+(n-1)) = new_target; %issue might be here?
%            end
%        end
%    end
    
    [row,col] = size(C);
    C_padded = mirroring2(C,n);
    thold = 5; %trying 3, 8. 6 causes black stripes. Which means the freq. peaks are replaced with zeros.
    % The resulting Fourier transform set to be equal to "noisy" extended
    % transform
    gauss_filt_img = C_padded;
 
    
    m = floor(n/2);
    i1 = 1; % indexes for the extended (mirrored) Fourier transform
    j1 = 1; % indexes for the extended (mirrored) Fourier transform
    % main loops over actual elements of "noisy" Foirier Transform with
    % n x n area around zero frequency excluded from the processing
    for i = 1:row
        for j = 1:col
            
        if ((i < floor(row/2)-m) || (i > floor(row/2)+m)) && ((j < floor(col/2)-m) || (j > floor(col/2)+m)) 
           % targeting window 
           target = C_padded(i1 : i1+n-1, j1 : j1+n-1); 
           % median over magnitude of the targeting window
           med = median(abs(target), "all");
           % magnitude in the center of the window (targeting value)
           test_val = abs(C(i,j));
           % ratio used in detector
           comp = test_val/med;
           % is center of the window a peak
            if comp >= thold
                % if so, multiply a targeting window by the Gaussian
                % surface
                new_target = target.*g; 
                % updating resulting Fourier transofrm by putting a window
                % with a peak and surrounding area suppressed
                gauss_filt_img(i1 : i1+n-1, j1 : j1+n-1) = new_target; 
            end
            
        end
        j1 = j1+1;    
        end
        i1 = i1+1;
        j1 = 1;
    end
    
    % cutting off extended edges of the resulting Fourier transform
    gauss_filt_img = gauss_filt_img(m+1:m+row, m+1: m+col);
    
denoised_img = uint8(real(ifft2(ifftshift(gauss_filt_img*size(img,1)*size(img,2))))); %is this right????
%figure('name','Image without noise')
%imshow(denoised_img);
%subplot(2,1,1);
%imshow(img);
%title('Original Image');
%axis tight
%subplot(2,1,2);
imshow(uint8(denoised_img));
%title('Gaussian Notch Filtered Image');
%axis tight
%%%%%end of function%%%%%%%      
end