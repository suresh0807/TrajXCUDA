//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void check_together(float a, float b, float c, float d, float e, float f, float g, float h, float i,float latx,float laty,float latz,float latxy,float latxz,float latyz)
{

float chk1,chk2,chk3,chk4,chk5,chk6;

chk1= (a-d);
chk2= (b-e);
chk3= (c-f);
chk4= (a-g);
chk5= (b-h);
chk6= (c-i);


if( fabs(chk1) > latx/2.0 || fabs(chk4) > latx/2.0 || fabs(chk2) > laty/2.0 || fabs(chk5) > laty/2.0 ||\
 fabs(chk3) > latz/2.0 || fabs(chk6) >latz/2.0))
{
if(cell_type == "orthorhombic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d=d+latx; else d=d-latx}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) e=e+laty; else e=e-laty}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f=f+latz; else f=f-latz}
if(fabs(chk4) > latx/2.0) { if (chk4 >0) g=g+latx; else g=g-latx}
if(fabs(chk5) > laty/2.0) { if (chk5 >0) h=h+laty; else h=h-laty}
if(fabs(chk6) > latz/2.0) { if (chk6 >0) i=i+latz; else i=i-latz}
}

else if(cell_type == "monoclinic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d=d+latx; else d=d-latx}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) {e=e+laty; d=d+latxy;} else {e=e-laty; d=d-latxy;}}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f=f+latz; else f=f-latz}
if(fabs(chk4) > latx/2.0) { if (chk4 >0) g=g+latx; else g=g-latx}
if(fabs(chk5) > laty/2.0) { if (chk5 >0) {h=h+laty; g=g+latxy;} else {h=h-laty; g=g-latxy;}}
if(fabs(chk6) > latz/2.0) { if (chk6 >0) i=i+latz; else i=i-latz}
}
}


}


void check_together(float a, float b, float c, float d, float e, float f, float latx,float laty,float latz,float latxy,float latxz,float latyz)
{

float chk1,chk2,chk3;

chk1= (a-d);
chk2= (b-e);
chk3= (c-f);


if( fabs(chk1) > latx/2.0 || fabs(chk2) > laty/2.0 ||  fabs(chk3) > latz/2.0)
{
if(cell_type == "orthorhombic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d=d+latx; else d=d-latx}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) e=e+laty; else e=e-laty}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f=f+latz; else f=f-latz}
}

else if(cell_type == "monoclinic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d=d+latx; else d=d-latx}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) {e=e+laty; d=d+latxy;} else {e=e-laty; d=d-latxy;}}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f=f+latz; else f=f-latz}
}
}


}