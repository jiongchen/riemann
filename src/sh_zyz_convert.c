#include <math.h>
#include <minpack.h>

typedef double scalarD;

void zyz_to_sh(const double *zyz, double *out) {
  //input
  scalarD aa = zyz[0];
  scalarD bb = zyz[1];
  scalarD cc = zyz[2];

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
  scalarD tt10;
  scalarD tt11;
  scalarD tt12;
  scalarD tt13;
  scalarD tt14;
  scalarD tt15;
  scalarD tt16;
  scalarD tt17;
  scalarD tt18;
  scalarD tt19;
  scalarD tt20;
  scalarD tt21;
  scalarD tt22;
  scalarD tt23;
  scalarD tt24;
  scalarD tt25;
  scalarD tt26;
  scalarD tt27;
  scalarD tt28;
  scalarD tt29;
  scalarD tt30;
  scalarD tt31;
  scalarD tt32;
  scalarD tt33;
  scalarD tt34;
  scalarD tt35;
  scalarD tt36;
  scalarD tt37;
  scalarD tt38;
  scalarD tt39;
  scalarD tt40;
  scalarD tt41;

  tt1=sqrt(5);
  tt2=4*aa;
  tt3=sin(tt2);
  tt4=cos(bb);
  tt5=3*bb;
  tt6=cos(tt5);
  tt7=tt1*tt3*tt6/8.0+7.0*tt1*tt3*tt4/8.0;
  tt8=4*cc;
  tt9=cos(tt8);
  tt10=sqrt(7);
  tt11=cos(tt2);
  tt12=5.0*tt10*tt11/8.0+3.0*tt10/8.0;
  tt13=tt1*tt10/4.0-tt1*tt10*tt11/4.0;
  tt14=2*bb;
  tt15=cos(tt14);
  tt16=tt1*tt11/8.0+7.0*tt1/8.0;
  tt17=4*bb;
  tt18=cos(tt17);
  tt19=tt16*tt18/8.0-tt10*tt13*tt15/4.0+tt1*tt10*tt12/8.0;
  tt20=sin(tt8);
  tt21=sqrt(2);
  tt22=1/pow(tt21,7);
  tt23=sin(bb);
  tt24=sin(tt5);
  tt25=3*tt22*tt1*tt3*tt24+7*tt22*tt1*tt3*tt23;
  tt26=3*cc;
  tt27=cos(tt26);
  tt28=1/pow(tt21,3);
  tt29=sin(tt14);
  tt30=sin(tt17);
  tt31=tt28*tt16*tt30-tt28*tt10*tt13*tt29;
  tt32=sin(tt26);
  tt33=tt1*tt10*tt3*tt4/8.0-tt1*tt10*tt3*tt6/8.0;
  tt34=2*cc;
  tt35=cos(tt34);
  tt36=-tt10*tt16*tt18/4.0+tt13*tt15/2.0+tt1*tt12/4.0;
  tt37=sin(tt34);
  tt38=3*tt22*tt1*tt10*tt3*tt23-tt22*tt1*tt10*tt3*tt24;
  tt39=cos(cc);
  tt40=-tt28*tt10*tt16*tt30-tt28*tt13*tt29;
  tt41=sin(cc);
  out[0]=tt19*tt20+tt7*tt9;
  out[1]=tt31*tt32+tt25*tt27;
  out[2]=tt36*tt37+tt33*tt35;
  out[3]=tt40*tt41+tt38*tt39;
  out[4]=tt1*tt10*tt16*tt18/8.0+tt1*tt13*tt15/4.0+3.0*tt12/8.0;
  out[5]=tt40*tt39-tt38*tt41;
  out[6]=tt36*tt35-tt33*tt37;
  out[7]=tt31*tt27-tt25*tt32;
  out[8]=tt19*tt9-tt7*tt20;
}

void zyz_to_sh_jac(const double *zyz, double *out) {
  //input
  scalarD aa = zyz[0];
  scalarD bb = zyz[1];
  scalarD cc = zyz[2];

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
  scalarD tt10;
  scalarD tt11;
  scalarD tt12;
  scalarD tt13;
  scalarD tt14;
  scalarD tt15;
  scalarD tt16;
  scalarD tt17;
  scalarD tt18;
  scalarD tt19;
  scalarD tt20;
  scalarD tt21;
  scalarD tt22;
  scalarD tt23;
  scalarD tt24;
  scalarD tt25;
  scalarD tt26;
  scalarD tt27;
  scalarD tt28;
  scalarD tt29;
  scalarD tt30;
  scalarD tt31;
  scalarD tt32;
  scalarD tt33;
  scalarD tt34;
  scalarD tt35;
  scalarD tt36;
  scalarD tt37;
  scalarD tt38;
  scalarD tt39;
  scalarD tt40;
  scalarD tt41;
  scalarD tt42;
  scalarD tt43;
  scalarD tt44;
  scalarD tt45;
  scalarD tt46;
  scalarD tt47;
  scalarD tt48;
  scalarD tt49;
  scalarD tt50;
  scalarD tt51;
  scalarD tt52;
  scalarD tt53;
  scalarD tt54;
  scalarD tt55;
  scalarD tt56;
  scalarD tt57;
  scalarD tt58;
  scalarD tt59;
  scalarD tt60;

  tt1=sqrt(5);
  tt2=4*aa;
  tt3=cos(tt2);
  tt4=cos(bb);
  tt5=3*bb;
  tt6=cos(tt5);
  tt7=tt1*tt3*tt6/2.0+7.0*tt1*tt3*tt4/2.0;
  tt8=4*cc;
  tt9=cos(tt8);
  tt10=pow(tt1,3);
  tt11=sin(tt2);
  tt12=2*bb;
  tt13=cos(tt12);
  tt14=4*bb;
  tt15=cos(tt14);
  tt16=-tt1*tt11*tt15/16.0+(-7.0)*tt1*tt11*tt13/4.0+(-7.0)*tt10*tt11/16.0;
  tt17=sin(tt8);
  tt18=sqrt(2);
  tt19=1/pow(tt18,3);
  tt20=sin(bb);
  tt21=sin(tt5);
  tt22=3*tt19*tt1*tt3*tt21+7*tt19*tt1*tt3*tt20;
  tt23=3*cc;
  tt24=cos(tt23);
  tt25=sin(tt12);
  tt26=1/pow(tt18,5);
  tt27=sin(tt14);
  tt28=-tt26*tt1*tt11*tt27-7*tt19*tt1*tt11*tt25;
  tt29=sin(tt23);
  tt30=sqrt(7);
  tt31=tt1*tt30*tt3*tt4/2.0-tt1*tt30*tt3*tt6/2.0;
  tt32=2*cc;
  tt33=cos(tt32);
  tt34=tt1*tt30*tt11*tt15/8.0+tt1*tt30*tt11*tt13/2.0-tt10*tt30*tt11/8.0;
  tt35=sin(tt32);
  tt36=3*tt19*tt1*tt30*tt3*tt20-tt19*tt1*tt30*tt3*tt21;
  tt37=cos(cc);
  tt38=tt26*tt1*tt30*tt11*tt27-tt19*tt1*tt30*tt11*tt25;
  tt39=sin(cc);
  tt40=(-3.0)*tt1*tt11*tt21/8.0+(-7.0)*tt1*tt11*tt20/8.0;
  tt41=tt1*tt30/4.0-tt1*tt30*tt3/4.0;
  tt42=tt1*tt3/8.0+7.0*tt1/8.0;
  tt43=tt30*tt41*tt25/2.0-tt42*tt27/2.0;
  tt44=1/pow(tt18,7);
  tt45=9*tt44*tt1*tt11*tt6+7*tt44*tt1*tt11*tt4;
  tt46=1/tt18;
  tt47=tt18*tt42*tt15-tt46*tt30*tt41*tt13;
  tt48=3.0*tt1*tt30*tt11*tt21/8.0-tt1*tt30*tt11*tt20/8.0;
  tt49=tt30*tt42*tt27-tt41*tt25;
  tt50=3*tt44*tt1*tt30*tt11*tt4-3*tt44*tt1*tt30*tt11*tt6;
  tt51=-tt18*tt30*tt42*tt15-tt46*tt41*tt13;
  tt52=5.0*tt30*tt3/8.0+3.0*tt30/8.0;
  tt53=tt42*tt15/8.0-tt30*tt41*tt13/4.0+tt1*tt30*tt52/8.0;
  tt54=tt1*tt11*tt6/8.0+7.0*tt1*tt11*tt4/8.0;
  tt55=tt19*tt42*tt27-tt19*tt30*tt41*tt25;
  tt56=3*tt44*tt1*tt11*tt21+7*tt44*tt1*tt11*tt20;
  tt57=-tt30*tt42*tt15/4.0+tt41*tt13/2.0+tt1*tt52/4.0;
  tt58=tt1*tt30*tt11*tt4/8.0-tt1*tt30*tt11*tt6/8.0;
  tt59=-tt19*tt30*tt42*tt27-tt19*tt41*tt25;
  tt60=3*tt44*tt1*tt30*tt11*tt20-tt44*tt1*tt30*tt11*tt21;
  out[0]=tt16*tt17+tt7*tt9;
  out[1]=tt28*tt29+tt22*tt24;
  out[2]=tt34*tt35+tt31*tt33;
  out[3]=tt38*tt39+tt36*tt37;
  out[4]=(-5.0)*tt30*tt11*tt15/16.0+5.0*tt30*tt11*tt13/4.0+(-15.0)*tt30*tt11/16.0;
  out[5]=tt38*tt37-tt36*tt39;
  out[6]=tt34*tt33-tt31*tt35;
  out[7]=tt28*tt24-tt22*tt29;
  out[8]=tt16*tt9-tt7*tt17;
  out[9]=tt43*tt17+tt40*tt9;
  out[10]=tt47*tt29+tt45*tt24;
  out[11]=tt49*tt35+tt48*tt33;
  out[12]=tt51*tt39+tt50*tt37;
  out[13]=-tt1*tt30*tt42*tt27/2.0-tt1*tt41*tt25/2.0;
  out[14]=tt51*tt37-tt50*tt39;
  out[15]=tt49*tt33-tt48*tt35;
  out[16]=tt47*tt24-tt45*tt29;
  out[17]=tt43*tt9-tt40*tt17;
  out[18]=4*tt53*tt9-4*tt54*tt17;
  out[19]=3*tt55*tt24-3*tt56*tt29;
  out[20]=2*tt57*tt33-2*tt58*tt35;
  out[21]=tt59*tt37-tt60*tt39;
  out[22]=0;
  out[23]=-tt59*tt39-tt60*tt37;
  out[24]=-2*tt57*tt35-2*tt58*tt33;
  out[25]=-3*tt55*tt29-3*tt56*tt24;
  out[26]=-4*tt53*tt17-4*tt54*tt9;
}

static double *g_SH;

static void lmder_callback(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag) {
  if ( *iflag == 1 ) {
    zyz_to_sh(x, fvec);
    unsigned int i = 0;
    for (i = 0; i < 9; ++i)
      fvec[i] -= g_SH[i];
  } else if ( *iflag == 2 ) {
    zyz_to_sh_jac(x, fjac);
  }
}

int sh_to_zyz(const double *sh, double *zyz, const unsigned int maxits) {
  g_SH = sh;
  
  int m = 9, n = 3;
  double fvec[9], fjac[27];
  int ipvt[3];
  int info;
  double tol = 1e-8;

  int lwa = 9*3+5*3+9;
  double wa[9*3+5*3+9];

  unsigned int i = 0;
  for (i = 0; i < maxits; ++i) {
    lmder1_(lmder_callback, &m, &n, zyz, fvec, fjac, &m,
            &tol, &info, &ipvt[0], &wa[0], &lwa);
    if ( info < 4 )
      break;
  }  
  return info;
}
