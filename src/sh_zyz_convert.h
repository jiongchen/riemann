#ifndef SH_ZYZ_CONVERT_H
#define SH_ZYZ_CONVERT_H

extern "C" {

  void sh_residual_(double *val, const double *zyz, const double *sh);
  void sh_residual_jac_(double *jac, const double *zyz, const double *sh);
  void sh_residual_hes_(double *hes, const double *zyz, const double *sh);

  void zyz_to_sh(const double *zyz, double *val);
  void zyz_to_sh_jac(const double *zyz, double *jac);
  void sh_to_zyz(const double *sh, double *zyz, const unsigned int maxits);
  
}

#endif
