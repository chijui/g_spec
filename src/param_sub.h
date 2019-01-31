#ifndef PARAM_SUB_H
#define PARAM_SUB_H

#ifndef DEF_H
#include "def.h"
#endif

#ifndef TYPES_H
#include "types.h"
#endif

#ifndef PARAM_CHECKER_H
#include "param_checker.h"
#endif

int gmx_parse_args(int *argc, char *argv[], fnm_param *fnm_param, spec_param *spec_param, w_param *w_param, POL_info *pol_info, real *delta, output_env_t *oenv);

int initialize_spec_param(spec_param *spec_param);
int initialize_w_param(w_param *w_param);
int set_w_param(w_param *w_param, spec_param *spec_param);
int set_waxis_nise(w_param *w_param);
int initialize_fnm_param(fnm_param *fnm_param);

int set_traj_param(traj_param *traj_param, spec_param *spec_param, w_param *w_param);
int set_BLAS_gemv_opt(BLAS_gemv_opt *gemv_opt);

#endif
