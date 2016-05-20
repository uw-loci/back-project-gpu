/* backprojects */
#include <math.h>
#include <matrix.h>
#include <mex.h>

#define MAX_LASERS 100
#define MAX_NX 20000

#define dist(v1, v2) \
	sqrt(((v1)[0]-(v2)[0])*((v1)[0]-(v2)[0])\
		+((v1)[1]-(v2)[1])*((v1)[1]-(v2)[1])\
		+((v1)[2]-(v2)[2])*((v1)[2]-(v2)[2]))

#define VOXEL_CHUNCK 32768
#define LASER_CHUNCK 32
#define CAM_CHUNCK 32

void __global__ calChunck(double *d1l, double *laserpos, double *voxels,
	double *cpos, double *d4l, double shift, double tpp, int nlasers, int nx,
	int nt, int nvoxels, double *sI, double intensity_correction, double *output)
{
	__shared__ double sum[LASER_CHUNCK * CAM_CHUNCK];
	int voxelIdx = threadIdx.x * CAM_CHUNCK + threadIdx.y;
	if (blockIdx.x >= nvoxels || threadIdx.x >= nlasers || threadIdx.y >= nx)
		return;
	double d2 = dist(&laserpos[threadIdx.x*3], &voxels[blockIdx.x*3]);
	double d3 = dist(&voxels[blockIdx.x*3], &cpos[threadIdx.y*3]);
	double d = d1l[threadIdx.x] + d2 + d3 + d4l[threadIdx.y];
	int tindex = (d-(shift))/(tpp) + 0.5;
	if ((tindex>=0) && (tindex<nt))
	{                          
		int index = threadIdx.x*nx*nt + tindex +threadIdx.y*nt;
		sum[voxelIdx] = sI[index%50000] * (intensity_correction);
	} else {
		sum[voxelIdx] = 0;
	}
	__syncthreads();
	if (threadIdx.x != 0 || threadIdx.y != 0)
		return;
	double result = 0;
	for (int i = 0; i < nlasers; i++) {
		for (int j = 0; j < nx; j++) {
			result += sum[i*CAM_CHUNCK+j];
		}
	}
	output[blockIdx.x] += result;
}

//Load Variables
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* Macros for the ouput and input arguments */
#define xsI prhs[0]
#define xlaserpos prhs[1]
#define xvoxels prhs[2]
#define xcpos prhs[3]
#define xlcop prhs[4]
#define xccop prhs[5]
#define xtpp prhs[6]
#define xshift prhs[7]
#define xlasernorm prhs[8]
#define xcameranorm prhs[9]

#define xoutput plhs[0]

cudaError_t rc;

double *sI = mxGetPr(xsI);
double *laserpos = mxGetPr(xlaserpos);
double *voxels = mxGetPr(xvoxels);
double *cpos = mxGetPr(xcpos);
double *geo_laser_cop = mxGetPr(xlcop);
double *geo_camera_cop = mxGetPr(xccop);
double *cameranormal = mxGetPr(xcameranorm);
double *lasernormal = mxGetPr(xlasernorm);

double tpp = *((double *)mxGetData(xtpp)); 
double shift = *((double *)mxGetData(xshift));

int msI = mxGetM(xsI);
int nlasers = mxGetN(xsI);
int nx = mxGetN(xcpos);
int nt = mxGetM(xsI)/nx;
int nvoxels = mxGetN(xvoxels);

double *d_sI, *d_laserpos, *d_voxels, *d_cpos;

rc = cudaMalloc((void **)&d_sI, sizeof(double)*msI*nlasers);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
rc = cudaMemcpy(d_sI, sI, sizeof(double)*msI*nlasers, cudaMemcpyHostToDevice);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));

rc = cudaMalloc((void **)&d_laserpos, sizeof(double)*3*nlasers);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
rc = cudaMemcpy(d_laserpos, laserpos, sizeof(double)*3*nlasers, cudaMemcpyHostToDevice);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));

rc = cudaMalloc((void **)&d_cpos, sizeof(double)*3*nx);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
rc = cudaMemcpy(d_cpos, cpos, sizeof(double)*3*nx, cudaMemcpyHostToDevice);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));

int x = 0,p=0;
int tpos = 0;

mexPrintf("nlasers: %d\n", nlasers);
mexPrintf("nvoxels: %d\n", nvoxels);
mexPrintf("nx: %d\n", nx);
mexPrintf("nt: %d\n", nt);
mexPrintf("tpp: %f\n", tpp);
mexPrintf("shift: %f\n", shift);
mexPrintf("First voxel: %f %f %f\n", voxels[0], voxels[1], voxels[2]);
mexPrintf("Laser cop: %f %f %f\n", geo_laser_cop[0], geo_laser_cop[1], geo_laser_cop[2]);


//Start Backproject
mxArray *out_array = xoutput = mxCreateDoubleMatrix(nvoxels,1,mxREAL);

double *output = mxGetPr(out_array);
double* d1l=new double[nlasers];
double* d4l=new double[nx];


for( tpos = 0;tpos<nlasers;tpos++) {
  d1l[tpos] = dist(geo_laser_cop,&laserpos[tpos*3]);
//  mexPrintf("geo_laser_cop [%f %f %f]\n", geo_laser_cop[0], geo_laser_cop[1], geo_laser_cop[2]);
}

for (x=0;x<nx;x++) {
  d4l[x] = dist(&cpos[x*3],geo_camera_cop);
}

double *d_output;
double *d_d1l, *d_d4l;

rc = cudaMalloc((void **)&d_d1l, sizeof(double)*nlasers);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
rc = cudaMemcpy(d_d1l, d1l, sizeof(double)*nlasers, cudaMemcpyHostToDevice);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));

rc = cudaMalloc((void **)&d_d4l, sizeof(double)*nx);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
rc = cudaMemcpy(d_d4l, d4l, sizeof(double)*nx, cudaMemcpyHostToDevice);
if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));

for (p = 0; p < nvoxels; p += VOXEL_CHUNCK) {
	mexPrintf("%d percent done\n", (p / (nvoxels / 10)) * 10);
	// mexPrintf("%d %d\n",p, nvoxels);
	mexEvalString("pause(.001);"); // to dump string.
	
	int vchunck = min(VOXEL_CHUNCK, nvoxels - p);
	
	rc = cudaMalloc((void **)&d_voxels, sizeof(double)*3*vchunck);
	if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
	rc = cudaMemcpy(d_voxels, &voxels[3*p], sizeof(double)*3*vchunck, cudaMemcpyHostToDevice);
	if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
	
	rc = cudaMalloc((void **)&d_output, sizeof(double)*vchunck);
	if (rc != cudaSuccess) printf(	"ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
	rc = cudaMemset(d_output, 0, sizeof(double)*vchunck);
	if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
	
	for (tpos = 0; tpos < nlasers; tpos += LASER_CHUNCK) {
		int lchunck = min(LASER_CHUNCK, nlasers - tpos);
		
		for (x = 0; x < nx; x += CAM_CHUNCK) {
			int xchunck = min(CAM_CHUNCK, nx - x);
	
			calChunck<<<VOXEL_CHUNCK, dim3(LASER_CHUNCK, CAM_CHUNCK)>>>(&d_d1l[tpos], &d_laserpos[3*tpos],
				d_voxels, &d_cpos[3*x], &d_d4l[x], shift, tpp, lchunck, xchunck, nt, vchunck, d_sI, 1.0, d_output);
		}
	}
	
	rc = cudaMemcpy(&output[p], d_output, sizeof(double)*vchunck, cudaMemcpyDeviceToHost);
	if (rc != cudaSuccess) printf("ERROR ON CUDA: %s\n", cudaGetErrorString(rc));
}

cudaFree(d_d1l);
cudaFree(d_d4l);
cudaFree(d_sI);
cudaFree(d_voxels);
cudaFree(d_laserpos);
cudaFree(d_cpos);
cudaFree(d_output);

/*
for(p=0;p<nvoxels;p++)
{
	break;
  if (p % (nvoxels / 10) == 0)
  {
      mexPrintf("%d percent done\n", (p / (nvoxels / 10)) * 10);
      // mexPrintf("%d %d\n",p, nvoxels);
      mexEvalString("pause(.001);"); // to dump string.
  }
    
  double thesum=0;
  double * voxel1 = &voxels[p*3];
  
  for( tpos = 0;tpos<nlasers;tpos++) {
      double * laserpos1 = &laserpos[tpos*3];
      double * lasernorm1 = &lasernormal[tpos*3];
      double d1 = d1l[tpos];
      double d2 = distance2(laserpos1,voxel1);
      
      for (x=0;x<nx;x++) {
          double * cpos1 = &cpos[x*3];
          double * cameranorm1 = &cameranormal[x*3];
          double d3 = distance2(voxel1,cpos1);
          double d4 = d4l[x];
          double d=d1+d2+d3+d4;
          
          //mexPrintf("d1: %f, d2: %f d3: %f, d4: %f\n", d1, d2, d3, d4);
          
          double vlv[3];
          double vcv[3];
          vlv[0]= (voxel1[0]-laserpos1[0])/d2;
          vlv[1]= (voxel1[1]-laserpos1[1])/d2;
          vlv[2]= (voxel1[2]-laserpos1[2])/d2;
          double dotlv = vlv[0]*lasernorm1[0]+vlv[1]*lasernorm1[1]+vlv[2]*lasernorm1[2];
//           mexPrintf("vlv: %f\n" , vlv[2]);
//           mexPrintf("laser normal: %f %f %f\n", lasernormal[0], lasernormal[1], lasernormal[2]);
          
          vcv[0]= (voxel1[0]-cpos1[0])/d3;
          vcv[1]= (voxel1[1]-cpos1[1])/d3;
          vcv[2]= (voxel1[2]-cpos1[2])/d3;
          double dotcv = vcv[0]*cameranorm1[0]+vcv[1]*cameranorm1[1]+vcv[2]*cameranorm1[2];
//           mexPrintf("vcv: %f\n", vcv[2]);
//           mexPrintf("camera normal: %f %f %f\n", cameranormal[0], cameranormal[1], cameranormal[2]);
          
          double intensity_correction = 1.0;          
//           intensity_correction =  sqrt(d2*d3);
//           intensity_correction =  d2*d3;
          
          int tindex = round((d-(shift))/tpp);
                  
          int index = 0;
          double tol=0.3;
          if(voxel1[0]>-tol && voxel1[0]<tol && voxel1[1]>46-tol && voxel1[1]<46+tol && voxel1[2]>-40-tol && voxel1[2]<-40+tol) {
                mexPrintf("---------------------------------------------------------\n");
                mexPrintf("Voxel 4000: [%f %f %f]\n", voxel1[0], voxel1[1], voxel1[2]);
                mexPrintf("cpos: [%f %f %f]\n", cpos1[0], cpos1[1], cpos1[2]);
                mexPrintf("lpos: [%f %f %f]\n", laserpos1[0], laserpos1[1], laserpos1[2]);
                mexPrintf("d1 %f, d2 %f, d3 %f, d4 %f\n", d1, d2, d3, d4);
                mexPrintf("INDEX: %d\n", tindex);
                mexPrintf("d %f\n", d);
                }
                      
          if ((tindex>=0) && (tindex<nt) && (dotlv>0) && (dotcv>0))
          //if ((tindex>=0) && (tindex<nt))
          {                          
              index = tpos*nx*nt + tindex +x*nt;                            
              thesum = thesum + sI[index] *intensity_correction;
          }

          /*
          double tindexd = ((d-(shift))/tpp);
          int tindexl = floor(tindexd), tindexu=ceil(tindexd);  
         
          if ((tindexl>=0) && (tindexu<nt))
          {            
              double w = tindexd-tindexl;
              int indexl = tpos*nx*nt + tindexl +x*nt, indexu = tpos*nx*nt + tindexu +x*nt;                                       
              thesum = thesum + ( (1-w) * sI[indexl] + w * sI[indexu])*intensity_correction;
              //sIout[index] = -1;
          }*//*
         
           
          if (x == 250 && p == 400 && tpos==0)

          {
              mexPrintf("tindex: %d\n", tindex);
              mexPrintf("d1: %f\n", d1);
              mexPrintf("d2: %f\n", d2);
              mexPrintf("d3: %f\n", d3);
              mexPrintf("d4: %f\n", d4);
              mexPrintf("Laser pos: %f %f %f\n", laserpos1[0], laserpos1[1], laserpos1[2]);
              mexPrintf("Point pos: %f %f %f\n", voxel1[0], voxel1[1], voxel1[2]);
              mexPrintf("Cam pos: %f %f %f\n", cpos1[0], cpos1[1], cpos1[2]);
              mexPrintf("x: %d\n", x);
              mexPrintf("pixel index: %d\n",  tindex +x*nt);
              //mexPrintf("pixel value: %f\n", sI[tpos*nx*nt + tindex +x*nt]);
              mexPrintf("pixel value: %f\n", sI[index]);
          }
          
      }
  }
  output[p] = thesum; 
}*/
mexPrintf("100 percent done\n");
delete[] d1l;
delete[] d4l;
}
