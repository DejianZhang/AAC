# Non-line-of-sight imaging with adaptive artifact cancellation Code & Datasets

This repository contains code for the paper _Non-line-of-sight imaging with adaptive artifact cancellation_ by Hongyuan Zhou, Ziyang Chen, Jumin Qiu, Sijia Zhong, Dejian Zhang, Tongbiao Wang, Qiegen Liu, and Tianbao Yu.

## Data

### "LT" and "Z" 

- Description: TOF data of retroreflective letters placed at a distance of approximately 0.8m from the wall.
- Resolution: 64 x 64
- Scanned Area: 0.8 m x 0.8 m planar wall
- Acquisition Method: Experiment 

### "bunny" and "T" 
- Description: TOF data of Lambertian objects placed at a distance of approximately 0.8m from the wall. 
- Resolution: 64 x 64
- Scanned Area: 0.8 m x 0.8 m planar wall
- Acquisition method: Simulation

## Code

### forward projection
The code implements the forward projection process: input the voxel information of object and obtains the TOF histogram. It is worth noting that for experimental data, forward projection needs to simulate the timing jitter of the experimental system. The variable t_d in the code is generated based on the timing jitter of the experimental system, and only one decay is performed when calculating the number of photons, as the experiment used a retroreflective object. For simulation data, the forward projection performs two rounds of attenuation when calculating the number of photons, as it uses Lambertian objects.

### main
These codes implement the algorithm of adaptive artifact cancellation:(1) subtract the TOF histogram convolved with Gaussian kernel from the original TOF histogram to generate a new sequence of flight time TOF histograms; (2) perform backprojection to obtain reconstructed voxels; (3) perform forward projection to obtain the TOF histogram of the reconstructed voxels; (4)calculate the SSIM between this histogram and the original histogram, change the sigma value of Gaussian convolution kernel, and repeat this process (here we plot the curve of SSIM against sigma, of course, other methods can be used to improve efficiency, such as genetic algorithm), (5) find the optimal sigma, and draw the results.