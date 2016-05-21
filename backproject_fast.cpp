/* backprojects */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>

#ifdef __APPLE__
    #include "OpenCL/opencl.h"
#else
    #include "CL/cl.h"
#endif

#define MAX_LASERS 100
#define MAX_NX 20000

#define dist(v1, v2) \
	sqrt(((v1)[0]-(v2)[0])*((v1)[0]-(v2)[0])\
		+((v1)[1]-(v2)[1])*((v1)[1]-(v2)[1])\
		+((v1)[2]-(v2)[2])*((v1)[2]-(v2)[2]))
		
#define min(a, b) ((a) < (b) ? (a) : (b))

#define VOXEL_CHUNCK 65536
#define LASER_CHUNCK 32
#define CAM_CHUNCK 32

inline void checkRet(int ret) {
	if (ret != 0) {
		printf("ret = %d\n", ret);
	}
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

	/* read source from file */
#define MAX_SOURCE_SIZE (0x100000)
	FILE *source;
	char *source_str;
	size_t source_size;

	source = fopen("./backproject_fast.cl", "r");
	if (!source) {
		fprintf(stderr, "Failed to load kernel.\n");
		return;
	}
	source_str = (char*)malloc(MAX_SOURCE_SIZE);
	source_size = fread(source_str, 1, MAX_SOURCE_SIZE, source);
	fclose(source);

	cl_int ret;

	cl_uint platformIdx;
	cl_platform_id *platformIds;
	cl_uint platformIdCount;
	ret = clGetPlatformIDs(0, NULL, &platformIdCount); checkRet(ret);
	platformIds = new cl_platform_id[platformIdCount];
	ret = clGetPlatformIDs(platformIdCount, platformIds, NULL); checkRet(ret);

	for (platformIdx = 0; platformIdx < platformIdCount; platformIdx++) {
		char name[128];
		clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_VENDOR, 128, name, NULL);
		// only use Nvidia or AMD graphics cards
		if (strncmp("NVIDIA", name, 6) == 0 || strncmp("AMD", name, 3) == 0)
			break;
	}

	if (platformIdx == platformIdCount) {
		printf("No proper graphic card could be found.\n");
		return;
	}
	
	{
		char info[1024];
		printf("Platform info:\n");
		ret = clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_PROFILE, 1024, info, NULL);
		printf("    profile: %s\n", info);
		ret = clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_VERSION, 1024, info, NULL);
		printf("    version: %s\n", info);
		ret = clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_NAME, 1024, info, NULL);
		printf("    name: %s\n", info);
		ret = clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_VENDOR, 1024, info, NULL);
		printf("    vendor: %s\n", info);
		ret = clGetPlatformInfo(platformIds[platformIdx], CL_PLATFORM_EXTENSIONS, 1024, info, NULL);
		printf("    extensions: %s\n", info);
		printf("\n");
	}

	/* use the first device */
	cl_device_id deviceId;
	ret = clGetDeviceIDs(platformIds[platformIdx], CL_DEVICE_TYPE_DEFAULT, 1, &deviceId, NULL);
	{
		char info[1024];
		printf("Device info:\n");
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_ADDRESS_BITS, 1024, info, NULL);
		printf("    address bits: %d\n", *(cl_uint *)info);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_EXTENSIONS, 1024, info, NULL);
		printf("    extensions: %s\n", info);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_GLOBAL_MEM_SIZE, 1024, info, NULL);
		printf("    global mem size: %u\n", *(cl_ulong *)info);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_LOCAL_MEM_SIZE, 1024, info, NULL);
		printf("    local mem size: %u\n", *(cl_ulong *)info);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_MAX_COMPUTE_UNITS, 1024, info, NULL);
		printf("    max compute units: %d\n", *(cl_uint *)info);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_MAX_WORK_GROUP_SIZE, 1024, info, NULL);
		printf("    max work group size: %ld\n", *(size_t *)info);
		
		cl_uint tmp;
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(tmp), &tmp, NULL);
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_MAX_WORK_ITEM_SIZES, 1024, info, NULL);
		printf("    max work item sizes: ");
		for (int i = 0; i < tmp; i++) {
			printf("%ld  ", ((size_t *)info)[i]);
		}
		printf("\n");
		
		ret = clGetDeviceInfo(deviceId, CL_DEVICE_NAME, 1024, info, NULL);
		printf("    device name: %s\n", info);
		
		printf("\n");
	}
	delete[] platformIds;

	cl_context context = clCreateContext(NULL, 1, &deviceId, NULL, NULL, &ret); checkRet(ret);
	cl_command_queue commandQueue = clCreateCommandQueue(context, deviceId, 0, &ret); checkRet(ret);
	
	cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
		(const size_t *)&source_size, &ret); checkRet(ret);
	ret = clBuildProgram(program, 1, &deviceId, NULL, NULL, NULL); checkRet(ret);
	if (ret != 0) {
		size_t len;
		clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, NULL, NULL, &len);
		char *log = new char[len];
		clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, len, log, NULL);
		printf("%s\n", log);
		delete log;
		free(source_str);
		return;
	}
	cl_kernel kernel = clCreateKernel(program, "backproject_fast", &ret); checkRet(ret);
	free(source_str);

	cl_mem d_sI, d_laserpos, d_lasernormal, d_voxels, d_cpos, d_cameranormal;

	d_sI = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*msI*nlasers, sI, &ret); checkRet(ret);

	d_laserpos = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*3*nlasers, laserpos, &ret); checkRet(ret);

	d_lasernormal = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*3*nlasers, lasernormal, &ret); checkRet(ret);

	d_cpos = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*3*nx, cpos, &ret); checkRet(ret);

	d_cameranormal = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*3*nx, cameranormal, &ret); checkRet(ret);

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


	cl_mem d_output;
	cl_mem d_d1l, d_d4l;

	d_d1l = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		 sizeof(double)*nlasers, d1l, &ret); checkRet(ret);

	d_d4l = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		 sizeof(double)*nx, d4l, &ret); checkRet(ret);

	// argument indices
#define IDX_D1L                  0
#define IDX_LASERPOS             1
#define IDX_LASERNORMAL          2
#define IDX_VOXELS               3
#define IDX_CPOS                 4
#define IDX_CAMERANORMAL         5
#define IDX_D4L                  6
#define IDX_SI                   7
#define IDX_SHIFT                8
#define IDX_TPP                  9
#define IDX_NLASERS              10
#define IDX_NX                   11
#define IDX_NT                   12
#define IDX_NVOXELS              13
#define IDX_INTENSITY_CORRECTION 14
#define IDX_OUTPUT               15

	ret = clSetKernelArg(kernel, IDX_SHIFT, sizeof(shift), (void*)&shift); checkRet(ret);
	ret = clSetKernelArg(kernel, IDX_TPP, sizeof(tpp), (void*)&tpp); checkRet(ret);
	ret = clSetKernelArg(kernel, IDX_NT, sizeof(nt), (void*)&nt); checkRet(ret);
	ret = clSetKernelArg(kernel, IDX_SI, sizeof(d_sI), (void*)&d_sI); checkRet(ret);
	double intensity_correction = 1.0;
	ret = clSetKernelArg(kernel, IDX_INTENSITY_CORRECTION, sizeof(intensity_correction), (void*)&intensity_correction); checkRet(ret);

	size_t offset[] = {0, 0, 0};
	size_t globalWorkSize[3], localWorkSize[3];
	localWorkSize[0] = 1;

	for (p = 0; p < nvoxels; p += VOXEL_CHUNCK) {
		mexPrintf("%d percent done\n", p * 100 / nvoxels);
		// mexPrintf("%d %d\n",p, nvoxels);
		mexEvalString("pause(.001);"); // to dump string.
		
		int vchunck = min(VOXEL_CHUNCK, nvoxels - p);
		globalWorkSize[0] = vchunck;
		
		d_voxels = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(double)*3*vchunck, &voxels[3*p], &ret); checkRet(ret);
		
		d_output = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double)*vchunck,
			NULL, &ret); checkRet(ret);
		double zero = 0.0;
		ret = clEnqueueFillBuffer(commandQueue, d_output, &zero, sizeof(zero), 0,
			sizeof(double)*vchunck, 0, NULL, NULL); checkRet(ret);
			
		ret = clSetKernelArg(kernel, IDX_VOXELS, sizeof(d_voxels), (void*)&d_voxels); checkRet(ret);
		ret = clSetKernelArg(kernel, IDX_NVOXELS, sizeof(vchunck), (void*)&vchunck); checkRet(ret);
		ret = clSetKernelArg(kernel, IDX_OUTPUT, sizeof(d_output), (void*)&d_output); checkRet(ret);
		
		for (tpos = 0; tpos < nlasers; tpos += LASER_CHUNCK) {
			int lchunck = min(LASER_CHUNCK, nlasers - tpos);
			globalWorkSize[1] = localWorkSize[1] = lchunck;
			
			cl_buffer_region d_d1l_off_Info = {sizeof(double)*tpos, sizeof(double)*lchunck};
			cl_mem d_d1l_off = clCreateSubBuffer(d_d1l, CL_MEM_READ_ONLY,
				CL_BUFFER_CREATE_TYPE_REGION, &d_d1l_off_Info, &ret); checkRet(ret);
			ret = clSetKernelArg(kernel, IDX_D1L, sizeof(d_d1l_off), (void*)&d_d1l_off); checkRet(ret);
			
			cl_buffer_region d_laserpos_off_info = {sizeof(double)*3*tpos, sizeof(double)*3*lchunck};
			cl_mem d_laserpos_off = clCreateSubBuffer(d_laserpos, CL_MEM_READ_ONLY,
				CL_BUFFER_CREATE_TYPE_REGION, &d_laserpos_off_info, &ret); checkRet(ret);
			ret = clSetKernelArg(kernel, IDX_LASERPOS, sizeof(d_laserpos_off), (void*)&d_laserpos_off); checkRet(ret);
			
			cl_buffer_region d_lasernormal_off_info = {sizeof(double)*3*tpos, sizeof(double)*3*lchunck};
			cl_mem d_lasernormal_off = clCreateSubBuffer(d_lasernormal, CL_MEM_READ_ONLY,
				CL_BUFFER_CREATE_TYPE_REGION, &d_lasernormal_off_info, &ret); checkRet(ret);
			ret = clSetKernelArg(kernel, IDX_LASERNORMAL, sizeof(d_lasernormal_off), (void*)&d_lasernormal_off); checkRet(ret);
			
			ret = clSetKernelArg(kernel, IDX_NLASERS, sizeof(lchunck), (void*)&lchunck); checkRet(ret);
			
			for (x = 0; x < nx; x += CAM_CHUNCK) {
				int xchunck = min(CAM_CHUNCK, nx - x);
				globalWorkSize[2] = localWorkSize[2] = xchunck;
				
				cl_buffer_region d_cpos_off_info = {sizeof(double)*3*x, sizeof(double)*3*xchunck};
				cl_mem d_cpos_off = clCreateSubBuffer(d_cpos, CL_MEM_READ_ONLY,
					CL_BUFFER_CREATE_TYPE_REGION, &d_cpos_off_info, &ret); checkRet(ret);
				ret = clSetKernelArg(kernel, IDX_CPOS, sizeof(d_cpos_off), (void*)&d_cpos_off); checkRet(ret);
				
				cl_buffer_region d_cameranormal_off_info = {sizeof(double)*3*x, sizeof(double)*3*xchunck};
				cl_mem d_cameranormal_off = clCreateSubBuffer(d_cameranormal, CL_MEM_READ_ONLY,
					CL_BUFFER_CREATE_TYPE_REGION, &d_cameranormal_off_info, &ret); checkRet(ret);
				ret = clSetKernelArg(kernel, IDX_CAMERANORMAL, sizeof(d_cameranormal_off), (void*)&d_cameranormal_off); checkRet(ret);
				
				cl_buffer_region d_d4l_off_info = {sizeof(double)*x, sizeof(double)*xchunck};
				cl_mem d_d4l_off = clCreateSubBuffer(d_d4l, CL_MEM_READ_ONLY,
					CL_BUFFER_CREATE_TYPE_REGION, &d_d4l_off_info, &ret); checkRet(ret);
				ret = clSetKernelArg(kernel, IDX_D4L, sizeof(d_d4l_off), (void*)&d_d4l_off); checkRet(ret);
				
				ret = clSetKernelArg(kernel, IDX_NX, sizeof(xchunck), (void*)&xchunck); checkRet(ret);
		
				ret = clEnqueueNDRangeKernel(commandQueue, kernel, 3, offset, globalWorkSize,
					localWorkSize, 0, NULL, NULL);
				checkRet(ret);
			}
		}
		ret = clEnqueueReadBuffer(commandQueue, d_output, CL_TRUE, 0, sizeof(double)*vchunck,
			&output[p], 0, NULL, NULL); checkRet(ret);
	}

	mexPrintf("100 percent done\n");
	delete[] d1l;
	delete[] d4l;
}
