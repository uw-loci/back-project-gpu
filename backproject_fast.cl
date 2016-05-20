#define LASER_CHUNCK 32
#define CAM_CHUNCK 32
		
__kernel void backproject_fast(__global double *d1l, __global double *laserpos,
	__global double *lasernormal, __global double *voxels, __global double *cpos,
	__global double *cameranormal, __global double *d4l, __global double *sI, 
	double shift, double tpp, int nlasers, int nx, int nt, int nvoxels,
	double intensity_correction, __global double *output)
{
	__local double sum[LASER_CHUNCK * CAM_CHUNCK];
	int voxelIdx = get_global_id(0);
	int laserIdx = get_local_id(1);
	int camIdx = get_local_id(2);
	int voxelSubIdx = laserIdx * CAM_CHUNCK + camIdx;
	
	double3 l = (double3)(laserpos[laserIdx*3], laserpos[laserIdx*3+1], laserpos[laserIdx*3+2]);
	double3 ln = (double3)(lasernormal[laserIdx*3], lasernormal[laserIdx*3+1], lasernormal[laserIdx*3+2]);
	double3 v = (double3)(voxels[voxelIdx*3], voxels[voxelIdx*3+1], voxels[voxelIdx*3+2]);
	double3 c = (double3)(cpos[camIdx*3], cpos[camIdx*3+1], cpos[camIdx*3+2]);
	double3 cn = (double3)(cameranormal[camIdx*3], cameranormal[camIdx*3+1], cameranormal[camIdx*3+2]);
	
	double d2 = distance(l, v);
	double d3 = distance(v, c);
	double3 vlv = (v - l)/d2;
	double dotlv = dot(vlv, ln);
	double3 vcv = (v - c)/d3;
	double dotcv = dot(vcv, cn);
	double d = d1l[laserIdx] + d2 + d3 + d4l[camIdx];
	
	int tindex = (d-(shift))/(tpp) + 0.5;
	if ((tindex>=0) && (tindex<nt) && (dotlv > 0) && (dotcv > 0))
	{                          
		int index = laserIdx*nx*nt + tindex + camIdx*nt;
		sum[voxelSubIdx] = sI[index] * (intensity_correction);
	} else {
		sum[voxelSubIdx] = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (laserIdx != 0 || camIdx != 0)
		return;
	double result = 0;
	for (int i = 0; i < nlasers; i++) {
		for (int j = 0; j < nx; j++) {
			result += sum[i*CAM_CHUNCK+j];
		}
	}
	output[voxelIdx] = result;
}