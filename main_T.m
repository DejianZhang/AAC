%% main file
clc
clear

c              = 3e8;   
z_trim     = 500;  
z_offset = 0;
width = 0.4;
bin_resolution = 5e-12; 
s_lamda_limit = 100 * 0.8/63;
sampling_coeff = 3;
isbackprop = 0;         
isdiffuse  = 0;         
K          = 2;         
snr        = 8e-1; 
load('Groundtruth_T.mat')
volc = TOF./max(TOF(:));
load('Lambertian_T.mat')
tof_data = TOF;

Q = 1:1:100;
[~,n] = size(Q);
NLOSvol = {}; 
[M,N,Z] = size(tof_data);
tof_G = zeros(M,N,Z);
tof_FP = zeros(M,N,Z);
ssim1 = zeros(size(Q));
corr1 = zeros(size(Q));
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
   
   ssim1(q)=ssim(vola,volc);
   corr1(q)=corr2(vola,volc);
   
   tof_FP = forward_projection_sim(vol);
   tof_FP = tof_FP./max(tof_FP(:));
   tof_data(:,:,1:500)=0;
   tof_data = tof_data./max(tof_data(:));

   tofssim1(q) = ssim(tof_FP,tof_data);
   
   display(q);
end
 
ssim1 = ssim1 - min(ssim1);
ssim1 = ssim1/max(ssim1);
corr1 = corr1 -min(corr1);
corr1 = corr1./max(corr1);
tofssim1 = tofssim1 -min(tofssim1);
tofssim1 = tofssim1./max(tofssim1);

[ssim1_max,ssim1_index]=max(ssim1);
[corr1_max,corr1_index]=max(corr1);
[tofssim1_max,tofssim1_index]=max(tofssim1);


display(Q(ssim1_index));
display(Q(corr1_index));
display(Q(tofssim1_index));

figure
set(gcf,"Position",[150,150,600,600])
hold on
plot(Q.*5,ssim1);
plot(Q.*5,corr1);
plot(Q.*5,tofssim1);

text(Q(ssim1_index).*5,ssim1_max+0.005,'o','color','black')
text(Q(corr1_index).*5,corr1_max+0.005,'o','color','blue')
text(Q(tofssim1_index).*5,tofssim1_max+0.005,'o','color','red')
hold off

figure
set(gcf,'position',[150 150 800 500])
subplot(2,3,1)
imagesc(NLOSvol{ssim1_index});
axis off

subplot(2,3,2)
imagesc(NLOSvol{corr1_index});
axis off

subplot(2,3,3)
imagesc(NLOSvol{tofssim1_index});
axis off

subplot(2,3,4)
imagesc(volc);
axis off

subplot(2,3,5)
imagesc(NLOSvol{1});
axis off

subplot(2,3,6)
imagesc(NLOSvol{q});
axis off
colormap('hot');


