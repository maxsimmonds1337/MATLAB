%% read in the image
I = imread('C:\Users\Max Simmonds\Desktop\MI5.png');

% set(subplot(1,3,1),'fontname','times')  % Set it to times
% set(gca, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [])
% title('Image obtained using home-made V-dipole antenna', 'FontSize', 14)
% 
% I = imread('earth_colour_no_strips.jpg');
% subplot(1,3,2)
% imshow(I)
% 
% set(gca,'fontname','times')  % Set it to times
% set(gca, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [])
% title('Image with MCIR colour map', 'FontSize', 14)

[row, col] = size(I);
 % 
% Imask = I; %% copy the image, to use it as a mask
% 
% %% loop through every pixel in the image, and convert to hex for each layer (red/green/blue)
for i = 1:row
    for j = 1:col/3
%         
%         %% extract the RBG colours
         R = I(i,j,1);
         G = I(i,j,2);
         B = I(i,j,3);
%         
%         %% if any of them are greater than the threshold, assume it's a cloud - make the pixel red
%         if (R > cloud_thresh)
%             
%             %% count the pixel, to generate the CCI
%             cloud_count = cloud_count + 1;
%             
%             Imask(i,j,1) = 225;
%             Imask(i,j,2) = 0;
%             Imask(i,j,3) = 225;            
    end
 end
% 
% subplot(1,3,3)
% imshow(Imask)
% set(gca,'fontname','times')  % Set it to times
% set(gca, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [])
% title('Processed image with cloud detection algorithm', 'FontSize', 14)
% 
% CCI = ((cloud_count)/(row*col)) * 100
% suptitle(sprintf('Cloud Coverage Index (CCI) = %0.2f%%', CCI));
% set(gca, 'color', 'k')