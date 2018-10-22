//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"
//central atom then the other two
float angle(float x, float y, float z, float x1, float y1, float z1, float x2, float y2, float z2)
{
  //float pi = 3.14159;
  float AB[3], AC[3];
  float ABdotAC=0.0;
  float lenAB=0.0, lenAC=0.0;
  float angleret;
  AB[0] = (x1-x);
  AB[1] = (y1-y);
  AB[2] = (z1-z);
  AC[0] = (x2-x);
  AC[1] = (y2-y);
  AC[2] = (z2-z);
  for(int k=0;k<3;k++)
  {
  ABdotAC += AB[k]*AC[k];
  }
  for(int k=0;k<3;k++)
  {
    lenAB += square(AB[k]);
    lenAC += square(AC[k]);
  }
  lenAB=sqrt(lenAB);
  lenAC=sqrt(lenAC);
  if(ABdotAC/(lenAB*lenAC) > 1.0 || ABdotAC/(lenAB*lenAC) < -1.0)
  {
    angleret = 3.14159;
  }
  else
  {
  angleret= acos(ABdotAC/(lenAB*lenAC)) ;//* 180/pi;
  }
  return (angleret);
}

float angle(float x1, float y1, float z1, float x2, float y2, float z2,float x3, float y3, float z3, float x4, float y4, float z4)
{
  //float pi = 3.14159;
  float AB[3], AC[3];
  float ABdotAC=0.0;
  float lenAB=0.0, lenAC=0.0;
  float angleret;
  AB[0] = (x2-x1);
  AB[1] = (y2-y1);
  AB[2] = (z2-z1);
  AC[0] = (x4-x3);
  AC[1] = (y4-y3);
  AC[2] = (z4-z3);
  for(int k=0;k<3;k++)
  {
  ABdotAC += AB[k]*AC[k];
  }
  //ABdotAC=fabs(ABdotAC);
  for(int k=0;k<3;k++)
  {
    lenAB += square(AB[k]);
    lenAC += square(AC[k]);
  }
  lenAB=sqrt(lenAB);
  lenAC=sqrt(lenAC);
  if(ABdotAC/(lenAB*lenAC) > 1.0 || ABdotAC/(lenAB*lenAC) < -1.0)
  {
    angleret = 3.14159;
  }
  else
  {
  angleret= acos(ABdotAC/(lenAB*lenAC)) ;//* 180/pi;
  }
  return (angleret);
}