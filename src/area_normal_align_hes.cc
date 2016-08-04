#include <algorithm>
#include <cmath>

namespace riemann {

static void area_normal_align_hes_d1(const double *x, const double eps, double *out) {
  const double ax = x[0], ay = x[1], az = x[2],
      bx = x[3], by = x[4], bz = x[5],
      cx = x[6], cy = x[7], cz = x[8];
  const double 
      tt1=bz-cz,
      tt2=tt1*tt1,
      tt3=az-bz,
      tt4=by-ay,
      tt5=tt4*cz+tt3*cy+ay*bz-az*by,
      tt6=tt5*tt5,
      tt7=sqrt(eps+tt6),
      tt8=1/(tt7*tt7*tt7),
      tt9=1/tt7,
      tt10=cy-by,
      tt11=0.5*tt10*tt1*tt9-0.5*tt10*tt1*tt6*tt8,
      tt12=cz-az,
      tt13=0.5*tt1*tt12*tt9-0.5*tt1*tt12*tt6*tt8,
      tt14=ay-cy,
      tt15=0.5*tt5*tt9,
      tt16=tt15+0.5*tt14*tt1*tt9-0.5*tt14*tt1*tt6*tt8,
      tt17=0.5*tt3*tt1*tt9-0.5*tt3*tt1*tt6*tt8,
      tt18=-0.5*tt5*tt9,
      tt19=tt18+0.5*tt4*tt1*tt9-0.5*tt4*tt1*tt6*tt8,
      tt20=tt10*tt10,
      tt21=tt18+0.5*tt10*tt12*tt9-0.5*tt10*tt12*tt6*tt8,
      tt22=0.5*tt14*tt10*tt9-0.5*tt14*tt10*tt6*tt8,
      tt23=tt15+0.5*tt3*tt10*tt9-0.5*tt3*tt10*tt6*tt8,
      tt24=0.5*tt4*tt10*tt9-0.5*tt4*tt10*tt6*tt8,
      tt25=tt12*tt12,
      tt26=0.5*tt14*tt12*tt9-0.5*tt14*tt12*tt6*tt8,
      tt27=0.5*tt3*tt12*tt9-0.5*tt3*tt12*tt6*tt8,
      tt28=tt15+0.5*tt4*tt12*tt9-0.5*tt4*tt12*tt6*tt8,
      tt29=tt14*tt14,
      tt30=tt18+0.5*tt3*tt14*tt9-0.5*tt3*tt14*tt6*tt8,
      tt31=0.5*tt4*tt14*tt9-0.5*tt4*tt14*tt6*tt8,
      tt32=tt3*tt3,
      tt33=0.5*tt4*tt3*tt9-0.5*tt4*tt3*tt6*tt8,
      tt34=tt4*tt4;
  const double rtn[81] = {0,0,0,0,0,0,0,0,0,0,0.5*tt2*tt9-0.5*tt2*tt6*tt8,tt11,0,tt13,tt16,0,tt17,tt19,0,tt11,0.5*tt20*tt9-0.5*tt20*tt6*tt8,0,tt21,tt22,0,tt23,tt24,0,0,0,0,0,0,0,0,0,0,tt13,tt21,0,0.5*tt25*tt9-0.5*tt25*tt6*tt8,tt26,0,tt27,tt28,0,tt16,tt22,0,tt26,0.5*tt29*tt9-0.5*tt29*tt6*tt8,0,tt30,tt31,0,0,0,0,0,0,0,0,0,0,tt17,tt23,0,tt27,tt30,0,0.5*tt32*tt9-0.5*tt32*tt6*tt8,tt33,0,tt19,tt24,0,tt28,tt31,0,tt33,0.5*tt34*tt9-0.5*tt34*tt6*tt8};
  std::copy(rtn, rtn+81, out);
}

static void area_normal_align_hes_d2(const double *x, const double eps, double *out) {
  const double ax = x[0], ay = x[1], az = x[2],
      bx = x[3], by = x[4], bz = x[5],
      cx = x[6], cy = x[7], cz = x[8];
  const double
      tt1=cz-bz,
      tt2=tt1*tt1,
      tt3=bz-az,
      tt4=ax-bx,
      tt5=tt4*cz+tt3*cx-ax*bz+az*bx,
      tt6=tt5*tt5,
      tt7=sqrt(eps+tt6),
      tt8=1/(tt7*tt7*tt7),
      tt9=1/tt7,
      tt10=bx-cx,
      tt11=0.5*tt10*tt1*tt9-0.5*tt10*tt1*tt6*tt8,
      tt12=az-cz,
      tt13=0.5*tt12*tt1*tt9-0.5*tt12*tt1*tt6*tt8,
      tt14=cx-ax,
      tt15=-0.5*tt5*tt9,
      tt16=tt15+0.5*tt14*tt1*tt9-0.5*tt14*tt1*tt6*tt8,
      tt17=0.5*tt3*tt1*tt9-0.5*tt3*tt1*tt6*tt8,
      tt18=0.5*tt5*tt9,
      tt19=tt18+0.5*tt4*tt1*tt9-0.5*tt4*tt1*tt6*tt8,
      tt20=tt10*tt10,
      tt21=tt18+0.5*tt10*tt12*tt9-0.5*tt10*tt12*tt6*tt8,
      tt22=0.5*tt10*tt14*tt9-0.5*tt10*tt14*tt6*tt8,
      tt23=tt15+0.5*tt3*tt10*tt9-0.5*tt3*tt10*tt6*tt8,
      tt24=0.5*tt4*tt10*tt9-0.5*tt4*tt10*tt6*tt8,
      tt25=tt12*tt12,
      tt26=0.5*tt14*tt12*tt9-0.5*tt14*tt12*tt6*tt8,
      tt27=0.5*tt3*tt12*tt9-0.5*tt3*tt12*tt6*tt8,
      tt28=tt15+0.5*tt4*tt12*tt9-0.5*tt4*tt12*tt6*tt8,
      tt29=tt14*tt14,
      tt30=tt18+0.5*tt3*tt14*tt9-0.5*tt3*tt14*tt6*tt8,
      tt31=0.5*tt4*tt14*tt9-0.5*tt4*tt14*tt6*tt8,
      tt32=tt3*tt3,
      tt33=0.5*tt4*tt3*tt9-0.5*tt4*tt3*tt6*tt8,
      tt34=tt4*tt4;
  const double rtn[81] = {0.5*tt2*tt9-0.5*tt2*tt6*tt8,0,tt11,tt13,0,tt16,tt17,0,tt19,0,0,0,0,0,0,0,0,0,tt11,0,0.5*tt20*tt9-0.5*tt20*tt6*tt8,tt21,0,tt22,tt23,0,tt24,tt13,0,tt21,0.5*tt25*tt9-0.5*tt25*tt6*tt8,0,tt26,tt27,0,tt28,0,0,0,0,0,0,0,0,0,tt16,0,tt22,tt26,0,0.5*tt29*tt9-0.5*tt29*tt6*tt8,tt30,0,tt31,tt17,0,tt23,tt27,0,tt30,0.5*tt32*tt9-0.5*tt32*tt6*tt8,0,tt33,0,0,0,0,0,0,0,0,0,tt19,0,tt24,tt28,0,tt31,tt33,0,0.5*tt34*tt9-0.5*tt34*tt6*tt8};
  std::copy(rtn, rtn+81, out);
}

static void area_normal_align_hes_d3(const double *x, const double eps, double *out) {
  const double ax = x[0], ay = x[1], az = x[2],
      bx = x[3], by = x[4], bz = x[5],
      cx = x[6], cy = x[7], cz = x[8];
  const double
      tt1=by-cy,
      tt2=tt1*tt1,
      tt3=ay-by,
      tt4=bx-ax,
      tt5=tt4*cy+tt3*cx+ax*by-ay*bx,
      tt6=tt5*tt5,
      tt7=sqrt(eps+tt6),
      tt8=1/(tt7*tt7*tt7),
      tt9=1/tt7,
      tt10=cx-bx,
      tt11=0.5*tt10*tt1*tt9-0.5*tt10*tt1*tt6*tt8,
      tt12=cy-ay,
      tt13=0.5*tt1*tt12*tt9-0.5*tt1*tt12*tt6*tt8,
      tt14=ax-cx,
      tt15=0.5*tt5*tt9,
      tt16=tt15+0.5*tt14*tt1*tt9-0.5*tt14*tt1*tt6*tt8,
      tt17=0.5*tt3*tt1*tt9-0.5*tt3*tt1*tt6*tt8,
      tt18=-0.5*tt5*tt9,
      tt19=tt18+0.5*tt4*tt1*tt9-0.5*tt4*tt1*tt6*tt8,
      tt20=tt10*tt10,
      tt21=tt18+0.5*tt10*tt12*tt9-0.5*tt10*tt12*tt6*tt8,
      tt22=0.5*tt14*tt10*tt9-0.5*tt14*tt10*tt6*tt8,
      tt23=tt15+0.5*tt3*tt10*tt9-0.5*tt3*tt10*tt6*tt8,
      tt24=0.5*tt4*tt10*tt9-0.5*tt4*tt10*tt6*tt8,
      tt25=tt12*tt12,
      tt26=0.5*tt14*tt12*tt9-0.5*tt14*tt12*tt6*tt8,
      tt27=0.5*tt3*tt12*tt9-0.5*tt3*tt12*tt6*tt8,
      tt28=tt15+0.5*tt4*tt12*tt9-0.5*tt4*tt12*tt6*tt8,
      tt29=tt14*tt14,
      tt30=tt18+0.5*tt3*tt14*tt9-0.5*tt3*tt14*tt6*tt8,
      tt31=0.5*tt4*tt14*tt9-0.5*tt4*tt14*tt6*tt8,
      tt32=tt3*tt3,
      tt33=0.5*tt4*tt3*tt9-0.5*tt4*tt3*tt6*tt8,
      tt34=tt4*tt4;
  const double rtn[81] = {0.5*tt2*tt9-0.5*tt2*tt6*tt8,tt11,0,tt13,tt16,0,tt17,tt19,0,tt11,0.5*tt20*tt9-0.5*tt20*tt6*tt8,0,tt21,tt22,0,tt23,tt24,0,0,0,0,0,0,0,0,0,0,tt13,tt21,0,0.5*tt25*tt9-0.5*tt25*tt6*tt8,tt26,0,tt27,tt28,0,tt16,tt22,0,tt26,0.5*tt29*tt9-0.5*tt29*tt6*tt8,0,tt30,tt31,0,0,0,0,0,0,0,0,0,0,tt17,tt23,0,tt27,tt30,0,0.5*tt32*tt9-0.5*tt32*tt6*tt8,tt33,0,tt19,tt24,0,tt28,tt31,0,tt33,0.5*tt34*tt9-0.5*tt34*tt6*tt8,0,0,0,0,0,0,0,0,0,0};
  std::copy(rtn, rtn+81, out);
}

void area_normal_align_hes(const double *x, const double eps, const int id, double *out) {
  switch ( id ) {
    case 0: area_normal_align_hes_d1(x, eps, out); break;
    case 1: area_normal_align_hes_d2(x, eps, out); break;
    case 2: area_normal_align_hes_d3(x, eps, out); break;
  }
}

}
