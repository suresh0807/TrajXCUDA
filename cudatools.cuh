/*
################# TrajXCUDA #########################
######### Suresh Kondati Natarajan ##################
##### Lehrstuhl fuer theoretische chemie ############
######## Ruhr Universitaet Bochum ###################
*/

#ifndef CUDATOOLS_CUH
#define CUDATOOLS_CUH


#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<cuda.h>
#include <complex.h>
#include <fftw3.h>

#include<iostream>
#include<string>
#include<fstream>

using namespace std;

extern int natoms,nstruct;
extern float *A, *Aint, *B, *VEL,*METAL; //Host
extern float *Aavg;
extern char c;
extern int *exch, *dev_exch;
extern float *rdf, *SD, *SDsum;
extern int *Ofrag, *Hbfrag , *dev_Ofrag, *dev_Hbfrag;
extern int *Ofragsum, *Hbfragsum;
extern float *dev_rdf, *dev_SD, *dev_SDsum;
extern float *dev_A,*dev_B,*dev_VEL; //Device
extern string Dummy;
extern string type;
extern string order;
//extern char Dummy[100];
//extern char type[50];
//extern char order[20];
extern float xvec,yvec,zvec,yxvec;
extern float *lattice, *dev_lattice;
//extern char felt[2];
extern string felt;
extern string keyword, unwrap;
extern string vacuum,choose_atoms,choose_dirn;
extern int bin;
extern int dircover;
extern float bin_rad;
//extern char Elt1[2],Elt2[2];
extern string Elt1,Elt2;
extern int nElt1,nElt2;
extern float pi;
extern float *rdfavgf, *SDavgf;
extern float rho,minbondist,maxbondist;
extern float *dev_rdfavgf, *dev_SDavgf;
extern float *rdfavg, *SDavg;
//extern char average[10];
extern string average;
extern float *dev_rdfavg, *dev_SDavg;
extern float minx,maxx,minz,maxz,miny,maxy,set_maxz,set_minz,set_maxx,set_minx,set_maxy,set_miny;
extern string set_max_x,set_max_y,set_max_z;
extern int count_metal;
extern float rdf_max_rad;
extern string rdf_metal_exclude;
extern float metal_atom_volume;
extern float minix,maxix,miniz,maxiz,miniy,maxiy;
extern float xrange,yrange,zrange;
extern float xint,yint,zint;
extern float *xtick,*ytick,*ztick,*dev_xtick,*dev_ytick,*dev_ztick;
extern int xsplit,ysplit,zsplit;
extern float mintop,maxtop;
extern int *density, *dev_density;
extern float *velocity, *dev_velocity,*Vdensity, *dev_Vdensity;
extern int skip, skips, check_index;
extern int SIZER, *dev_MSDchaos;
extern float *SD1, *SDsum1,  *MSDtics;
extern float *dev_SD1 , *dev_SDsum1;
extern float *SDavg1;
extern int *MSDchaos, MSDsteps;
extern float *dev_SDavg1;
//extern char command[50];
extern string command,filter_density,choose_atoms_from,HBlifetime,unwrapout,HB_histograms,density_plot,strict_dirn,file_path;
extern int symbol,color;
extern float size,bondistOH,bondistOM;
extern float bondist, Mbondist,bondist_int_bulk,bondist_int_bulk2;
extern int switchint, maxneigh, minneigh;
extern float timestep, diffuse_time, total_time, Hbond_angle_dev, max_O_O;
extern int num_bins, SD_store, origins;
extern string veloc,diffusion_kernel,transient_HB_definition,HB_crossing;
extern string cell_type, metal_species,write_cube,scope,plot,HB_for,msd_for,diffco_for,rdf_for,vdos_for,include_transient_bonds,dipole_orient,atop_orient;
extern int nA_int, nB_int, nA_bulk, nB_bulk,max_frames, *Aintnum, *Bintnum, *Abulknum, *Bbulknum;
extern float *A_int, *B_int, *A_bulk, *B_bulk;
extern string diffuse_direction,HB_criterium_set,rdf_between,HB_lifestyle,lifetime_for;
extern int xsrt, xend, xski,start_frame,transtime;
extern int *MSDcounter, *dev_MSDcounter, *dev_MSDchaos;
extern float MSDint;
extern float *A_all, *B_all;
extern float Coulomb_cutoff;
extern string ext_charge;

void read_inputs();
void density_top_avg_pos();
void density_side_avg_pos();
void compute_densityvelocitygrid();
void compute_densitygrid();
void compute_densitygrid1();
void compute_velocitygrid();
void compute_exchange();
void compute_Vdist(void);
void compute_coverage(void);
void compute_resident_time(void);
void compute_angledist();
void compute_adatoms(void);
void compute_lifetime(int *A, int *B, float time);
void compute_rdf(void);
void compute_rdf_frames(void);
void compute_tetrawater(void);
void compute_Hbonds(void);
void compute_Hbonds1(void);
void compute_Hbonds2(void);
void compute_HbondsSSP(void);
void compute_diffusion(void);
void compute_diffusion1(void);
void compute_diffusion_lifetime();
void compute_diffusion_z();
void select_atoms(string whichmat);
void compute_Zdist(void);
void compute_VDOS(void);
void compute_Ofrag(void);
void compute_coulomb_charge(void);
float square(float x);
float angle(float x, float y, float z, float x1, float y1, float z1, float x2, float y2, float z2);
float angle(float x1, float y1, float z1, float x2, float y2, float z2,float x3, float y3, float z3, float x4, float y4, float z4);
float minimum(float *a,int natoms,int stride);
float maximum(float *a,int natoms,int stride);
__global__ void density_calc(float *A,int *density,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick);
__global__ void density_calc(float *A,int *density,int *exch,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick);
__global__ void reduce(float *a, float *b, int bin, int nstruct);
__global__ void reducef(float *a, float *b, int bin, int nstruct, int nElt1);
__global__ void distmat(float *x, float*y, float *a, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec,float yxvec, float rdf_max_rad);
__global__ void distmatmono(float *x, float*y, float *a, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec,float yxvec, float rdf_max_rad);
__global__ void velocity_calc(float *A,float *B,float *velocity,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick);
__global__ void Ofragmat(float *x, float*y, int *a,float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec);
__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec);
__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float *latti);
__global__ void SD_calc(float *A,float *B,float *C,int SD_store,int nElt1, int i, int j,int origins,int skips,int xsrt,int xend, int xski);
__global__ void SD_calc(float *A,float *B,float *C,int SD_store,int nElt1, int i,int nstruct,int skips,int xsrt,int xend, int xski);
__global__ void SDreduce(float *a, float *b, int SD_store, int nElt1);
__global__ void SDreducef(float *a, float *b, int nElt1, int SD_store, int origins,int skips);
__global__ void SDreducef(float *a, float *b, int nElt1, int SD_store, int nstruct,int skips,int i);
__global__ void Xdist_calc(float *A,int *density,int nstruct,int natoms,int zsplit, float *ztick, int Dirn);
__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec);
__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float *latti);
__global__ void Hbondmat(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec);
__global__ void Hbondmat(float *x, float*y, int *a, int *b, float *c,float bondist, int nElt1, int nElt2, int nstruct, float *latti);
__global__ void Hbondmat(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid);
__global__ void Hbond_calc(float *A,float *Hbond_density,float *angle_density,float *distance_density, float *OHdensity, float *O_Hdensity, float *angdensity,int nstruct,int natoms,int zsplit, float *ztick, float angle, float max_O_O,int Dirn);
__global__ void Hbondmatmono(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec);
__global__ void Hbondmatmono(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float *latti);
__global__ void Hbondmatmono(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid);
__global__ void VAF_calc(float *A,float *B,float *C,int *D,int SD_store,int nElt1, int i, int j,int origins,int skips,int whichwater);
__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float rdf_max_rad, int whichwater);
__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater);
__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater,float mid,int dirn);
__global__ void distmatmono(float *x, float*y, float *a,int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float rdf_max_rad, int whichwater);
__global__ void distmatmono(float *x, float*y, float *a,int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater);
__global__ void distmatmono(float *x, float*y, float *a,int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater,float mid,int dirn);
__global__ void HBAF_calc(float *A,float *B,float *C,int SD_store,int initHBnum, int i, int j,int origins,int skips);
__global__ void SDreducefHB(float *a, float *b, int initHBnum, int SD_store, int origins,int skips);
__global__ void SDreduceHB(float *a, float *b, int SD_store, int initHBnum);
__global__ void density_filter(int *density, int *densityf, int mainsplit,int sub1split, int sub2split);
__global__ void Vdist_calc(float *A,float *B, float *density,float *xdensity,float *ydensity,float *zdensity,int nstruct,int natoms,int zsplit, float *ztick, int Dirn);
__global__ void SD_calc(float *A,float *B,float *C,int *D, int *E, int *F, int SD_store,int nElt1, int i, int j,int origins,int skips, int xsrt,int xend, int xski,int whichwater);
__global__ void SDreduce(float *a, float *b, int SD_store, int nElt1,int fairy);
__global__ void boxmat(float *x, int *a, int nElt1, int nstruct,float maxx, float minx, float maxy,float miny, float maxz, float minz);
__global__ void covermat(float *x, int *a, int nElt1, int nstruct, float minaz, float maxaz,int dirn);
__global__ void covermatmono(float *x, int *a, int nElt1, int nstruct, float minaz, float maxaz,int dirn);
__global__ void covermat(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti, int dirn, float surf,float mintop, float maxtop);
__global__ void covermatmono(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti, int dirn, float surf,float mintop, float maxtop);
__global__ void covermat(float *x, int *a, float minbondist, float maxbondist, int nElt1, int nstruct, int dirn, float surf,float maxtop,float mintop);
__global__ void SDreduce(float *a, float *b, int SD_store, int nElt1, int *MSDchaos, int fairy);
 __global__ void SD_calc(float *A,float *B,float *C,int *D, int *E, int *F, int SD_store,int nElt1, int i, int j,int origins,int skips, int xsrt,int xend, int xski,int whichwater,int chase);
 __global__ void Hbondmat(float *x, float*y, int*exch, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid);
 __global__ void Hbondmatmono(float *x, float*y, int*exch, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid);
__global__ void covermat(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti, float mintop, float maxtop,int dirn, float surf);
__global__ void covermatmono(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti , float mintop, float maxtop,int dirn, float surf);
__global__ void Xdist_calc(float *A, float *B,float *tetradensity, int *density,int nstruct,int natoms,int zsplit, float *ztick, int Dirn);
__global__ void covermatmono(float *x, float *y, int *a, float minbondist, float maxbondist, int index,float *latti, int nElt1, int nstruct);
__global__ void covermat(float *x, float *y, int *a, float minbondist, float maxbondist, int index,float *latti, int nElt1, int nstruct);
__global__ void Hbondmat(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid,int dirn);
__global__ void Hbondmatmono(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid,int dirn);
__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct,float *latti,int dirn);
__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct,float *latti,int dirn);
__global__ void distmatmono(float *x, float*y, float *a, int nElt1, int nElt2, int nstruct, float *latti, float cutoff);
__global__ void distmat(float *x, float*y, float *a, int nElt1, int nElt2, int nstruct, float *latti, float cutoff);
#endif