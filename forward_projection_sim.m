%% forward projection, for simulation data
function tof_FP = forward_projection_sim(vol)
tof_FP = zeros(64,64,2048);
bin_resolution = 5e-12;
c              = 3e8; 
width=0.8;
X_voxel=width/64;         
Y_voxel=width/64;         
Z_voxel=1.5/1000; 

vol_BP = vol(503:988,:,:);
vol_BP = permute(vol_BP, [2 3 1]);
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
Z = Z+502;
P = MAXC(:); 
P =P./max(P);

for x=1:64
    for y=1:64 
        distance = ((X_voxel.*(X-x)).^2+(Y_voxel.*(Y-y)).^2+(Z_voxel.*Z./2).^2).^0.5;
        optical_path = distance(:).*2;
        cosa=(Z_voxel.*Z)./distance;
        cosa=cosa(:);
        D = round((P.*(cosa.^4)./(distance.^4))); %lambertian object
        optical_path = repelem(optical_path,D);
        t = optical_path/c;       
        delta = 0:bin_resolution:2047*bin_resolution;
        t_hist = hist(t,delta);
        tof_FP(x,y,:)=t_hist;    
    end
end


