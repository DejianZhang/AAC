%% forward projection, fro experimetal data
function tof_FP = forward_projection(vol,t_d)
tof_FP = zeros(64,64,2048);
bin_resolution = 5e-12;
c              = 3e8; 
width=0.8;
X_voxel=width/63;         
Y_voxel=width/63;         
Z_voxel=1.5/1000./2; 


vol_BP = permute(vol, [2 3 1]);
for i = 1:64
    for j = 1:64
        vol_ABP = vol_BP(i,j,:);
        vol_ABP= vol_ABP(:);
        [maxVal,~] = max(vol_ABP);
        vol_ABP(vol_ABP<maxVal) = 0;
        vol_BP(i,j,:) = vol_ABP;
    end
end
[MAXC,INDEX]=max(vol_BP,[],3,"linear");
[X,Y,Z] = ind2sub(size(vol_BP), INDEX);
X = X(:);
Y = Y(:);
Z = Z(:);
P = MAXC(:); 
P =P./max(P);



for x=1:64
    for y=1:64 
        distance = ((X_voxel.*(X-x)).^2+(Y_voxel.*(Y-y)).^2+(Z_voxel.*Z).^2).^0.5;
        distance = repelem(distance,10) + t_d(1:40960); %add some timing jitter
        optical_path = distance(:).*2;
        cosa=Z_voxel.*repelem(Z,10)./distance;
        D = round(repelem(P,10).*10.*(cosa(:).^2)./(distance(:).^2)); %retroreflective object
        D(D>500)=500;
        optical_path = repelem(optical_path,D);
       % optical_path = optical_path + t_d(1:sum(D));
        t = optical_path/c;       
        delta = 0:bin_resolution:2048*bin_resolution;   
        t_hist = histcounts(t,delta);
        tof_FP(x,y,:)=t_hist;  
    end
end


