%% back projection, adopted from LCT NLOS.
function vol = back_projection(tof_data,width,bin_resolution,z_offset,z_trim,isdiffuse)
% Confocal Non-Light-of-Sight (C-NLOS) reconstruction procedure for paper
% titled "Confocal Non-Line-of-Sight Imaging Based on the Light Cone Transform"
% by Matthew O'Toole, David B. Lindell, and Gordon Wetzstein.
%
% 

% Constants
c              = 3e8;   % Speed of light (meters per second)
% Adjustable parameters
isbackprop = 1;         % Toggle backprojection
K          = 0;         % Downsample data to (4 ps) * 2^K = 16 ps for K = 2


% Load scene & set visualization parameter

        
% Because the scene is diffuse, toggle the diffuse flag and 
% adjust SNR value correspondingly.


N = size(tof_data,1);        % Spatial resolution of data
M = size(tof_data,3);        % Temporal resolution of data
range = M.*c.*bin_resolution; % Maximum range for histogram
    
% Downsample data to 16 picoseconds
for k = 1:K
    M = M./2;
    bin_resolution = 2*bin_resolution;
    tof_data = tof_data(:,:,1:2:end) + tof_data(:,:,2:2:end);
    z_trim = round(z_trim./2);
    z_offset = round(z_offset./2);
end
    
% Set first group of histogram bins to zero (to remove direct component)
tof_data(:,:,1:z_trim) = 0;
    
% Define NLOS blur kernel 
psf = definePsf(N,M,width./range);

% Compute inverse filter of NLOS blur kernel
fpsf = fftn(psf);
if (~isbackprop)
    invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./snr);
else
    invpsf = conj(fpsf);
end

% Define transform operators
[mtx,mtxi] = resamplingOperator(M);

% Permute data dimensions
data = permute(tof_data,[3 2 1]);

% Define volume representing voxel distance from wall
grid_z = repmat(linspace(0,1,M)',[1 N N]);


% Step 1: Scale radiometric component
if (isdiffuse)
    data = data.*(grid_z.^4);
else
    data = data.*(grid_z.^2);
end

% Step 2: Resample time axis and pad result
tdata = zeros(2.*M,2.*N,2.*N);
tdata(1:end./2,1:end./2,1:end./2)  = reshape(mtx*data(:,:),[M N N]);

% Step 3: Convolve with inverse filter and unpad result
tvol = ifftn(fftn(tdata).*invpsf);
tvol = tvol(1:end./2,1:end./2,1:end./2);

% Step 4: Resample depth axis and clamp results
vol  = reshape(mtxi*tvol(:,:),[M N N]);
vol  = max(real(vol),0);

% Crop and flip reconstructed volume for visualization
ind = round(M.*2.*width./(range./2));
vol = vol(:,:,end:-1:1);
vol = vol((1:ind)+z_offset,:,:);


function psf = definePsf(U,V,slope)
% Local function to computeD NLOS blur kernel

x = linspace(-1,1,2.*U);
y = linspace(-1,1,2.*U);
z = linspace(0,2,2.*V);
[grid_z,grid_y,grid_x] = ndgrid(z,y,x);

% Define PSF
psf = abs(((4.*slope).^2).*(grid_x.^2 + grid_y.^2) - grid_z);
psf = double(psf == repmat(min(psf,[],1),[2.*V 1 1]));
psf = psf./sum(psf(:,U,U));
psf = psf./norm(psf(:));
psf = circshift(psf,[0 U U]);


function [mtx,mtxi] = resamplingOperator(M)
% Local function that defines resampling operators

mtx = sparse([],[],[],M.^2,M,M.^2);

x = 1:M.^2;
mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
mtx  = spdiags(1./sqrt(x)',0,M.^2,M.^2)*mtx;
mtxi = mtx';

K = log(M)./log(2);
for k = 1:round(K)
    mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
    mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
end
