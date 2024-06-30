%% main file
clc
clear
c              = 3e8;   
z_trim     = 500;  
z_offset = 0;
width = 0.4;
bin_resolution = 5e-12; 
isbackprop = 0;         
isdiffuse  = 0;                  
snr        = 8e-1; 
load('Z.mat')
tof_data(:,:,1:500)=0;

load('t_d.mat')

Q = 10:1:60;
[~,n] = size(Q);
NLOSvol = {}; 
[M,N,Z] = size(tof_data);
tof_G = zeros(M,N,Z);
tof_FP = zeros(M,N,Z);
tofssim1 = zeros(size(Q));


for q = 1:n
 
   for i = 1:M
        for j = 1:N
            a = tof_data(i,j,:);
            a = a(:);
            B = imgaussfilt(a,Q(q));
            tof_G(i,j,:) = B;
        end
   end
   
   tof = tof_data - tof_G;
   
   vol = back_projection(tof,width,bin_resolution,z_offset,z_trim,isdiffuse);
   volb = permute(vol, [2 3 1]);
   volb = volb./max(volb(:));
   vola = max(volb,[],3);
   NLOSvol{q}=vola; 
   
   tof_FP = forward_projection(vol,t_d);
   tof_FP = tof_FP./max(tof_FP(:));
   tof_data = tof_data./max(tof_data(:));
   tofssim1(q) = ssim(tof_FP,tof_data);
   
   display(q);


end

[tofssim1_max,tofssim1_index]=max(tofssim1);
tofssim1 = tofssim1 -min(tofssim1);
tofssim1 = tofssim1./max(tofssim1);
display(Q(tofssim1_index));

figure(3)
plot(Q.*5,tofssim1);
text(Q(tofssim1_index).*5,tofssim1_max+0.005,'o','color','red')


figure(4)
subplot(1,3,1)
imagesc(NLOSvol{tofssim1_index});
axis off
subplot(1,3,2)
imagesc(NLOSvol{1});
axis off
subplot(1,3,3)
imagesc(NLOSvol{q});
axis off
colormap('hot');




