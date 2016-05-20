/* backprojects */
#include <math.h>
#include <matrix.h>
#include <mex.h>

#define MAX_LASERS 200
#define MAX_NX 20000

inline double distance2(double* v1,double *v2)
{
  double a = v1[0]-v2[0];
  double b = v1[1]-v2[1];
  double c = v1[2]-v2[2];
  
  return sqrt(a*a+b*b+c*c);
}

inline double round(double x) { return floor(x+0.5); }

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

#define xoutput plhs[0]
//#define xdebug plhs[1]

double *sI = mxGetPr(xsI);
double *laserpos = mxGetPr(xlaserpos);
double *voxels = mxGetPr(xvoxels);
double *cpos = mxGetPr(xcpos);
double *geo_laser_cop = mxGetPr(xlcop);
double *geo_camera_cop = mxGetPr(xccop);

double tpp = *((double *)mxGetData(xtpp)); 
double shift = *((double *)mxGetData(xshift));

int nlasers = mxGetN(xsI);
int nx = mxGetN(xcpos);
int nt = mxGetM(xsI)/nx;
int nvoxels = mxGetN(xvoxels);

mexPrintf("nlasers: %d\n", nlasers);
mexPrintf("nvoxels: %d\n", nvoxels);
mexPrintf("nx: %d\n", nx);
mexPrintf("nt: %d\n", nt);
mexPrintf("tpp: %f\n", tpp);
mexPrintf("shift: %f\n", shift);
mexPrintf("First voxel: %f %f %f\n", voxels[0], voxels[1], voxels[2]);
mexPrintf("Laser cop: %f %f %f\n", geo_laser_cop[0], geo_laser_cop[1], geo_laser_cop[2]);

//return;
//double * cpos1l;
/*for(int i=0; i<nx; i++) {
    cpos1l = &cpos[i*3];
    mexPrintf("cpos %d: %f %f %f\n", i, cpos1l[0], cpos1l[1], cpos1l[2]);
}*/

int x=0,p=0;

mxArray *out_array = xoutput = mxCreateDoubleMatrix(nvoxels,1,mxREAL);
double *output = mxGetPr(out_array);
//mxArray *sIdebug = xdebug = mxCreateDoubleMatrix(nt,nx,mxREAL);
//double *sIout = mxGetPr(sIdebug);

//for(x=0;x<nx*nt;x++)
//  sIout[x]=sI[x];

int tpos = 0;
double d1l[MAX_LASERS];
double d4l[MAX_NX];
int offindex=0;

for( tpos = 0;tpos<nlasers;tpos++) {
  d1l[tpos] = distance2(geo_laser_cop,&laserpos[tpos*3]);
//  mexPrintf("geo_laser_cop [%f %f %f]\n", geo_laser_cop[0], geo_laser_cop[1], geo_laser_cop[2]);
//  mexPrintf("d1,tpos %d: %f\n", tpos, d1l[tpos]);
}

for (x=0;x<nx;x++) {
  d4l[x] = distance2(&cpos[x*3],geo_camera_cop);
}

for(p=0;p<nvoxels;p++)
{
  if (p % (nvoxels / 10) == 0)
  {
      mexPrintf("%d percent done\n", (p / (nvoxels / 10)) * 10);
      mexEvalString("pause(.001);"); // to dump string.
  }
    
  double thesum=0;
  double * voxel1 = &voxels[p*3];
  
  for( tpos = 0;tpos<nlasers;tpos++) {
      double * laserpos1 = &laserpos[tpos*3];
      //double d1 = distance2(geo_laser_cop,laserpos1);
      double d1 = d1l[tpos];
      double d2 = distance2(laserpos1,voxel1);
      //mexPrintf("d1: %f, d2: %f \n", d1, d2);
      for (x=0;x<nx;x++) {
          double * cpos1 = &cpos[x*3];
          double d3 = distance2(voxel1,cpos1);
          //double d4 = distance2(cpos1,geo_camera_cop);
          double d4 = d4l[x];
          double d=d1+d2+d3+d4;
          double intensity_correction = 1.0;

          if (true) // do intensity correction?
              intensity_correction =  d2*d3;    
          
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
          
          if ((tindex>=0) && (tindex<nt))
          {                          
              index = tpos*nx*nt + tindex + x*nt;
//               if(x==1)
//                   mexPrintf("%d : %f %f %f %f\n", index, d1, d2, d3, d4);
              
              thesum = thesum + sI[index] *intensity_correction;
              //sIout[index] = -1;
          }
          else {
                    if(tindex<offindex || offindex==0) offindex=tindex;
 //                 thesum=thesum+0.5*intensity_correction;
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
          }*/
         
           
          if (x == 250 && p == 400 && tpos==0)
//          if (sI[index]!=0)
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
}

  mexPrintf("offindex: %i\n", offindex);
}
