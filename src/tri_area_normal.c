typedef double scalarD;

void tri_area_normal_(double *out, const double *x) {
  //input
  scalarD ax = x[0];
  scalarD ay = x[1];
  scalarD az = x[2];
  scalarD bx = x[3];
  scalarD by = x[4];
  scalarD bz = x[5];
  scalarD cx = x[6];
  scalarD cy = x[7];
  scalarD cz = x[8];

  //temp
  scalarD tt1;
  scalarD tt2;
  scalarD tt3;
  scalarD tt4;
  scalarD tt5;
  scalarD tt6;

  tt1=bz-az;
  tt2=cy-by;
  tt3=by-ay;
  tt4=cz-bz;
  tt5=cx-bx;
  tt6=bx-ax;
  out[0]=(tt3*tt4-tt1*tt2)/2.0;
  out[1]=(tt1*tt5-tt6*tt4)/2.0;
  out[2]=(tt6*tt2-tt3*tt5)/2.0;      
}

void tri_area_normal_jac_(double *out, const double *x) {
  //input
  scalarD ax = x[0];
  scalarD ay = x[1];
  scalarD az = x[2];
  scalarD bx = x[3];
  scalarD by = x[4];
  scalarD bz = x[5];
  scalarD cx = x[6];
  scalarD cy = x[7];
  scalarD cz = x[8];

  //temp
  scalarD tt1;
  scalarD tt2;
  scalarD tt3;
  scalarD tt4;
  scalarD tt5;
  scalarD tt6;
  scalarD tt7;
  scalarD tt8;
  scalarD tt9;

  tt1=-bz;
  tt2=-cy;
  tt3=-cz;
  tt4=-bx;
  tt5=-by;
  tt6=-cx;
  tt7=-ay;
  tt8=-az;
  tt9=-ax;
  out[0]=0;
  out[1]=(cz+tt1)/2.0;
  out[2]=(tt2+by)/2.0;
  out[3]=(tt3+bz)/2.0;
  out[4]=0;
  out[5]=(cx+tt4)/2.0;
  out[6]=(cy+tt5)/2.0;
  out[7]=(tt6+bx)/2.0;
  out[8]=0;
  out[9]=0;
  out[10]=(tt3+az)/2.0;
  out[11]=(cy+tt7)/2.0;
  out[12]=(cz+tt8)/2.0;
  out[13]=0;
  out[14]=(tt6+ax)/2.0;
  out[15]=(tt2+ay)/2.0;
  out[16]=(cx+tt9)/2.0;
  out[17]=0;
  out[18]=0;
  out[19]=(bz+tt8)/2.0;
  out[20]=(tt5+ay)/2.0;
  out[21]=(tt1+az)/2.0;
  out[22]=0;
  out[23]=(bx+tt9)/2.0;
  out[24]=(by+tt7)/2.0;
  out[25]=(tt4+ax)/2.0;
  out[26]=0;
}
