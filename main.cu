//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


#include "cudatools.cuh"

float Coulomb_cutoff;
float mintop,maxtop;
char c;
string Dummy;
string type;
string order;
int xsplit,ysplit,zsplit;
int skip,skips;
int color,symbol,switchint;
float size, bondist, Mbondist, minbondist,maxbondist;
//char Elt1[2],Elt2[2];
string Elt1,Elt2,dipole_orient,atop_orient;
string set_max_x,set_max_y,set_max_z;
float minx,maxx,minz,maxz,miny,maxy,set_minz,set_maxz,set_minx,set_maxx,set_miny,set_maxy;
int natoms,nstruct;
int *exch, *dev_exch;
float *A, *Aint, *B, *VEL, *METAL ; //Host
float *A_int, *B_int, *A_bulk, *B_bulk;
float *A_all, *B_all;
float *dev_A,*dev_B; 
float *Aavg;
int dircover;
float xvec,yvec,zvec,yxvec;
float *lattice, *dev_lattice;
//char felt[2];
string felt,choose_dirn;
int nElt1,nElt2;
string average;
string unwrap,choose_atoms;
string command,unwrapout;
int SIZER;
int maxneigh, minneigh;
float timestep, diffuse_time, total_time,bondistOH,bondistOM;
int num_bins, SD_store, origins;
float *SD, *SDsum;
float *dev_SD , *dev_SDsum; 
float *SD1, *SDsum1;
float *dev_SD1 , *dev_SDsum1;
float *SDavgf;
float *dev_SDavgf;
float *SDavg;
float *dev_SDavg;
float *SDavg1;
float *dev_SDavg1;
float MSDint;
int max_frames, *dev_MSDchaos, check_index;
string keyword,filter_density,HBlifetime,HB_criterium_set,density_plot;
string vacuum, transient_HB_definition,HB_lifestyle;
string veloc,diffusion_kernel,rdf_between,HB_crossing,strict_dirn,choose_atoms_from;
int count_metal;
string rdf_metal_exclude,diffuse_direction,HB_histograms,write_cube;
float metal_atom_volume, Hbond_angle_dev,bondist_int_bulk,bondist_int_bulk2;
float rdf_max_rad, max_O_O;
string cell_type, metal_species,scope,plot,msd_for,HB_for,diffco_for,rdf_for, vdos_for, include_transient_bonds,lifetime_for;
int nA_int, nB_int, nA_bulk, nB_bulk, *Aintnum, *Bintnum, *Abulknum, *Bbulknum;
int xsrt, xend, xski,start_frame;
float *velocity, *dev_velocity,*Vdensity, *dev_Vdensity;
float *dev_VEL;
int transtime;
string file_path;
string ext_charge;
float  *MSDtics;
int *MSDchaos, MSDsteps;
int *MSDcounter, *dev_MSDcounter;



int main(void)
{
  

  
 cout <<"#2015"<<endl;
 cout <<"################# TrajXCUDA #########################"<<endl;
 cout <<"########## Suresh Kondati Natarajan #################"<<endl;
 cout <<"##### Lehrstuhl fuer theoretische chemie ############"<<endl;
 cout <<"########## Ruhr Universitaet Bochum #################"<<endl<<endl;  
 cout <<"You can do one of the following operations:"<<endl<<endl;
 cout <<"1. Density profile - side view animation"<<endl;
 cout <<"2. Density profile - top view"<<endl;
 cout <<"3. Radial distribution function"<<endl;
 cout <<"4. Velocity profile - side view animation"<<endl;
 cout <<"5. Water fragments counter"<<endl;
 cout <<"6. Percent coverage of species at surface (need to be tested - orthorhombic only)"<<endl;
 cout <<"7. Diffusion coefficient"<<endl;
 cout <<"8. Vibrational density of states"<<endl;
 cout <<"9. Z-distribution"<<endl;
 cout <<"10. Hbonds"<<endl;
 cout <<"11. exchange"<<endl;
 cout <<"12. V-distribution"<<endl<<endl; 
 
 //char lElt1[2],lElt2[2];
 //string lElt1,lElt2;
 skip=1;
 SIZER=100;
 scope="bulk";
 ifstream elts;
 elts.open("input.dat",ios::in);
 while (!elts.eof())
 { 
 elts >> keyword;
 if (keyword.at(0) == '#') {elts.ignore(256,'\n');}
 else
 {
 if(keyword == "compute"){elts >> type;
  if(type=="density") {switchint=1;}
  else if(type=="density-top") {switchint=2;}
  else if(type=="rdf") {switchint=3;}
  else if(type=="densityvelocity"|| type =="densityvelocity-top") {switchint=4;} 
  else if(type=="O-fragments") {switchint=5;}
  else if(type=="percentage-coverage"){switchint=6;}
  else if(type=="diffusion-coefficient"){switchint=7;}
  else if(type=="VDOS"){switchint=8;}
  else if(type=="Z-dist"){switchint=9;}
  else if(type=="Hbonds"){switchint=10;}
  else if(type=="exchange"){switchint=11;}
  else if(type=="V-dist"){switchint=12;}
  else if(type=="orient"){switchint=13;}
  else if(type=="residenttime"){switchint=14;}
  else if(type=="Adatoms"){switchint=15;}
  else if(type=="tetrawater"){switchint=16;}
  else if(type=="coulomb"){switchint=17;}
  else {cout << type <<" is Not a valid value for "<<keyword << endl; exit(0);}}
 else if(keyword == "unwrap_traj"){elts >> unwrap;}
 else if(keyword == "element_1"){elts >> Elt1 ;}
 else if(keyword == "element_2"){elts >> Elt2 ;}
 else if(keyword == "surface_element"){elts >> Elt1;}
 else if(keyword == "adsorbate_element"){elts >> Elt2;}
 else if(keyword == "max_bond_distance"){elts >> bondist;}
 else if(keyword == "make_grid"){elts >> xsplit >> ysplit >> zsplit;}
 else if(keyword == "atoms_between"){elts >> mintop >> maxtop;}
 else if(keyword == "rdf_between"){elts >> rdf_between >> mintop >> maxtop;}
 else if(keyword == "vacuum_direction")
 {
 elts >> vacuum;
 if (vacuum.at(0) == 'z') { dircover=2; if(type == "density" || type =="densityvelocity") {order = "XYZ";} else if (type == "density-top"|| type =="densityvelocity-top") {order = "YZX";} }
 else if (vacuum.at(0) == 'y') { dircover=1; if(type == "density" || type =="densityvelocity") {order = "ZXY";} else if (type == "density-top"|| type =="densityvelocity-top") {order = "XYZ";} }
 else if (vacuum.at(0) == 'x') { dircover=0; if(type == "density" || type =="densityvelocity") {order = "YZX";} else if (type == "density-top"|| type =="densityvelocity-top") {order = "ZXY";} }
 else { cout << " enter x,y or z for 'vacuum_direction' keyword"<< endl;}
 }
 else if(keyword == "rdf_bins"){elts >> bin;}
 else if(keyword == "write_cube"){elts >> write_cube;}
 else if(keyword == "rdf_max_rad"){elts >> rdf_max_rad;} 
 else if(keyword == "equilibrium_metal_bond_distance"){elts >> Mbondist;}
 else if(keyword == "metal_coordination_number"){elts >> maxneigh >>minneigh;}
 else if(keyword == "skip_frames"){elts >> skip;}
 else if(keyword == "msd_skip_restarts" || keyword == "HB_skip_restarts" || keyword == "vdos_skip_restarts"){elts >> skips;}
 else if(keyword == "msd_timestep_fs" || keyword == "HB_timestep_fs"|| keyword == "vdos_timestep_fs"){elts >> timestep;}
 else if(keyword == "msd_lookback_ps" || keyword == "HB_lookback_ps" || keyword == "vdos_lookback_ps"){elts >> diffuse_time;}
 else if(keyword == "average_positions"){elts >> average;}
 else if(keyword == "xfarbe_symbol"){elts >> symbol;}
 else if(keyword == "xfarbe_size"){elts >> size;}
 else if(keyword == "xfarbe_color"){elts >> color;}
 else if(keyword == "read_velocities"){elts >> veloc;}
 else if(keyword == "HBlifetime"){elts >> HBlifetime;}
 else if(keyword == "scope"){elts >> scope;}
 else if(keyword == "choose_atoms"){elts >>Elt1 >>choose_dirn>> strict_dirn>>minbondist >> maxbondist>>metal_species; choose_atoms="yes";}
 else if(keyword == "metal_species"){elts >> metal_species;}
 else if(keyword == "rdf_metal_exclude"){elts >> rdf_metal_exclude >> metal_atom_volume;}
 else if(keyword == "set_maxz"){elts >> set_maxz;}
 else if(keyword == "set_minz"){elts >> set_minz;}
 else if(keyword == "set_maxx"){elts >> set_maxx;}
 else if(keyword == "set_minx"){elts >> set_minx;}
 else if(keyword == "set_maxy"){elts >> set_maxy;}
 else if(keyword == "set_miny"){elts >> set_miny;}
 else if(keyword == "check_index"){elts >> check_index;}
 else if(keyword == "choose_atoms_from"){elts >>Elt1 >>choose_dirn>> strict_dirn>>minbondist >> maxbondist>>check_index; choose_atoms_from="yes";}
 else if(keyword == "set_max_x"){elts >> set_max_x;}
 else if(keyword == "set_max_y"){elts >> set_max_y;}
 else if(keyword == "set_max_z"){elts >> set_max_z;}
 else if(keyword == "max_frames"){elts >> max_frames;}
 else if(keyword == "start_frame"){elts >> start_frame;}
 else if(keyword == "file_path"){elts >> file_path;}
 else if(keyword == "msd_for"){elts >> msd_for;}
 else if(keyword == "transient_HB_definition"){elts >> transient_HB_definition;}
 else if(keyword == "HB_for"){elts >> HB_for; }
 else if(keyword == "lifetime_for"){elts >> lifetime_for;}
 else if(keyword == "unwrapout"){elts >> unwrapout;}
 else if(keyword == "filter_density"){elts >> filter_density;}
 else if(keyword == "max_O_O"){elts >> max_O_O;}
 else if(keyword == "diffusion_kernel"){elts >> diffusion_kernel; if(diffusion_kernel=="z-binning") {elts >> MSDsteps;}}
 else if(keyword == "bondist_int_bulk"){elts >> bondist_int_bulk >> bondist_int_bulk2;}
 else if(keyword == "Hbond_angle_dev"){elts >> Hbond_angle_dev;}
 else if(keyword == "include_transient_bonds"){elts >> include_transient_bonds;}
 else if(keyword == "rdf_for"){elts >> rdf_for;}
 else if(keyword == "diffuse_direction") {elts >> diffuse_direction;}
 else if(keyword == "vdos_for"){elts >> vdos_for;}
 else if(keyword == "atop_orient"){elts >> atop_orient>>bondistOM;}
 else if(keyword == "dipole_orient"){elts >> dipole_orient>>bondistOH;}
 else if(keyword == "HB_crossing"){elts >> HB_crossing;}
 else if(keyword == "plot"){elts >> plot;}
 else if(keyword == "HB_lifestyle"){elts >> HB_lifestyle;  if(HB_lifestyle=="transient"){elts >> transtime;}}
 else if(keyword == "density_plot"){elts >> density_plot;}
 else if(keyword == "cell_type"){elts >> cell_type;}
 else if(keyword == "HB_histograms"){elts >> HB_histograms;}
 else if(keyword == "HB_criterium_set"){elts >> HB_criterium_set;}
 else if(keyword == "coulomb_cutoff"){elts >> Coulomb_cutoff;}
 else if(keyword == "use_external_charges"){elts >> ext_charge;}
 else {cout << "Keyword "<<keyword<<" not supported"<<endl;}
 }
 }
 elts.close();
 //elts.getline(type,SIZER);

//Reads the input file input.dat and based on the string in first line 
//decides what to do!
/* if(type=="density") {switchint=1;}
 else if(type=="density-top") {switchint=2;}
 else if(type=="rdf") {switchint=3;}
 else if(type=="densityvelocity") {switchint=4;} 
 else if(type=="O-fragments") {switchint=5;}
 else if(type=="percentage-coverage"){switchint=6;}
 else if(type=="diffusion-coefficient"){switchint=7;}
 else if(type=="VDOS"){switchint=8;}
 else {cout << "Not a valid keyword" << endl; exit(0);}
*/ 
 switch (switchint)
 {   case 1:     cout <<"Density profile- side view animation chosen:" <<endl <<endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >> zsplit;
		 //elts >> order;
		 //elts >> skip;
		 //elts >> average;
                 //elts.close();
		 read_inputs();
                 cout << nElt1 << " "<< Elt1 <<" atoms are there in each structure"<< endl;
                 if(average=="yes")
                 {
		   density_side_avg_pos();
                 }
                 else
                 {
                  compute_densitygrid();
                 }
		 break;
		 
     case 2:	 cout <<"Density profile- top view chosen:"<<endl<<endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >> zsplit;
		 //elts >> mintop >> maxtop;
		 //elts >> order;
		 //elts >> skip;
		 //elts >> average;    
                 //if(average=="yes")
                 //{
                 //  elts >> symbol;
                 //  elts >> size;
		 //  elts >> color;     
                 //}
                 if(dipole_orient == "yes")
		 {
		   Elt2="H";
		 }
                 cout << "Atom chosen: " << Elt1 <<endl;
                 //elts.close();
		 read_inputs();
		 cout <<nElt1 <<" "<<Elt1 << " atoms are there in each structure"<<endl;
                 if(average=="yes")
                 {
                  density_top_avg_pos();
                 }
                 else
                 {
                  compute_densitygrid();               
                 }
                 
                 free(Aavg);
		 break;
		 
     case 3:     cout<<"Radial distribution function chosen:"<<endl<<endl;
                 //elts >> Elt1 >> Elt2;
		 //elts >> bin;
		 //elts >> skip;
		 //elts.close();
		 read_inputs();
		 cout <<nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 cout <<nElt2<<" "<<Elt2<<" atoms are there in each structure"<<endl;
                 compute_rdf_frames();
		 //compute_rdf();
		 break;
		 
     case 4:	 cout<<"Velocity profile- side view animation chosen:"<<endl<<endl;
                 //elts >> Elt1;
                 //elts >> xsplit >> ysplit>> zsplit;
		 //elts >> order;
		 //elts >> skip;
		 //elts >> average;
                 //cout << "Atom chosen: "<<Elt1<<endl;
                 //elts.close(); 
		 read_inputs();
                 cout << nElt1 <<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 compute_densityvelocitygrid();
		 break;
		 
     case 5:     cout<<"Water fragments counter chosen:"<<endl<<endl;
                 Elt1="O";
                 Elt2="H";
                 //elts>>bondist;
                 //elts.close();
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 cout << nElt2<<" "<<Elt2<<" atoms are there in each structure"<<endl;
                 compute_Ofrag();
		 break;
		 
     case 6:     cout <<"Percentage coverage of "<<Elt2<<" on the chosen surface:"<<endl;;
                 //elts >> Elt1;
                 //elts >> Elt2;
		 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 cout << Elt2 << " atoms are chosen for analysis" << endl;
		 //elts >> Mbondist;
		 //cout << Mbondist<<endl;
		 //elts >> bondist;
		 //elts >> maxneigh;		 
                 //elts.close();
		 read_inputs();
                 cout<<nElt1 << " "<<Elt1<< " atoms are there in each structure" << endl;
                 cout<<nElt2 << " "<<Elt2<< " atoms are there in each structure" << endl;
                 compute_coverage();
		 break;
		 
     case 7:     cout << "You are computing diffusion coefficient:" << endl;
                 //elts >> Elt1;
		 //elts >> Elt2;
		 Elt2 = Elt1;
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 //elts >> timestep; //in fs
		 //elts >> diffuse_time; //in ps
		 //elts >> skips;
		 /*if(unwrap == "yes")
		 {
		 cout << "Unwrapping the trajectory"<<endl;  
		 if(cell_type =="orthorhombic"){system("unwrapper.x ortho");}
		 else if(cell_type =="monoclinic"){system("unwrapper.x mono");}
		 }*/
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 if(diffusion_kernel == "long_time") compute_diffusion1();
                 else if (diffusion_kernel == "binning") compute_diffusion_lifetime();
		 else if (diffusion_kernel == "z-binning") compute_diffusion_z();
		 //else if (diffusion_kernel == "binning") compute_diffusion();
                 cout << "MSD data given in (angstrom)^2 in msd.data"<<endl;
		 cout << "Diffusion coefficient data in (angstrom)^2/fs in diffco.data"<<endl;
		 cout << "However, the best way to find diffusion coefficient is"<<endl
		      << "to plot msd.data in gnuplot and fit a linear plot and "<<endl
		      << "take its slope"<<endl;
		 break;
		 
     case 8:     cout << "You are computing VDOS:" << endl;
                 //elts >> Elt1;
		 //elts >> Elt2;
		 Elt2 = Elt1;
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 //elts >> timestep; //in fs
		 //elts >> diffuse_time; //in ps
		 //elts >> skips;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 compute_VDOS();
		 break;
		 
     case 9:     cout << "You are computing Z-distribution:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 compute_Zdist();
		 break;
		 
     case 10:    cout << "You are computing Hbonds:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
		 cout << nElt2<<" "<<Elt2<<" atoms are there in each structure"<<endl;
		 if(HB_crossing == "yes" && HB_lifestyle!="SSP"&& HB_lifestyle!="switchbond") compute_Hbonds2();
		 else if(HB_crossing == "yes" && HB_lifestyle=="SSP") compute_HbondsSSP();
		 else if(HB_crossing == "yes" && HB_lifestyle=="switchbond") {cout<<"switchbond type chosen"<<endl;compute_HbondsSSP();}
		 else if (HB_crossing == "no") compute_Hbonds1();
		 break;
		 
      case 11:     cout <<"Exchange of "<<Elt1<<" atoms"<<endl;;
                 //elts >> Elt1;
                 //elts >> Elt2;
		 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 //cout << Elt2 << " atoms are chosen for analysis" << endl;
		 //elts >> Mbondist;
		 //cout << Mbondist<<endl;
		 //elts >> bondist;
		 //elts >> maxneigh;		 
                 //elts.close();
		 read_inputs();
                 cout<<nElt1 << " "<<Elt1<< " atoms are there in each structure" << endl;
                 //cout<<nElt2 << " "<<Elt2<< " atoms are there in each structure" << endl;
                 compute_exchange();
		 break;	 
		 
      case 12:     cout << "You are computing Velocity-distribution:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 compute_Vdist();
		 break;
		 
      case 13:     cout << "You are computing water orientation:" << endl;
                 Elt1="O";
                 Elt2="H";
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
                 compute_angledist();
		 break;
		 
      case 14:    cout << "You are computing Resident times:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
		 compute_resident_time();
		 break;	 
		 
      case 15:    cout << "You are computing Adatoms:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
                 cout << Elt1 << " atoms are chosen for analysis" << endl;
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
		 compute_adatoms();
		 break;	
       
       case 16:  cout << "You are computing Tetrahedral parameter of water:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
		 compute_tetrawater();
		 break;
		 
       case 17:  cout << "You are computing Coulomb charge:" << endl;
                 //elts >> Elt1;
		 //elts >> xsplit >> ysplit >>zsplit;
		 //elts >> order;
		 //elts >> skip;
                 //elts >> set_max >> set_maxz >> set_minz;
                 //elts.close();
		 read_inputs();
                 cout << nElt1<<" "<<Elt1<<" atoms are there in each structure"<<endl;
		 compute_coulomb_charge();
		 break;
}
 

cudaDeviceReset();
return 0;
}
