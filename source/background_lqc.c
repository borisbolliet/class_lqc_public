/** @file background_lqc.c Documented background_lqc module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background_lqc evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background_lqc, i.e. to integrate
 *    the background_lqc equations, and store all background_lqc quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background_lqc quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondence between redshift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background_lqc parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantities {B} can be found as a function of conformal
 * time by integrating the background_lqc equations.
 *
 * 3. some other quantities that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be inferred analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_lqc_functions() returns all background_lqc
 *    quantities {A} as a function of quantities {B}.
 *
 * - background_lqc_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_lqc_functions().
 *
 * - the result is stored in the form of a big table in the background_lqc
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantities {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most precise way is to call directly
 * background_lqc_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other quantity, we can call background_lqc_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_lqc_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background_lqc quantities, we should call first
 * background_lqc_tau_of_z() and then background_lqc_at_tau().
 *
 *
 * In order to save time, background_lqc_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantities, or rarely useful
 * quantities). Each line in the interpolation table is a vector whose
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantities like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantities needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_lqc_init() at the beginning
 * -# background_lqc_at_tau(), background_lqc_tau_of_z() at any later time
 * -# background_lqc_free() at the end, when no more calls to the previous functions are needed
 */

#include "background_lqc.h"

/**
 * Background_Lqc quantities at given conformal time tau.
 *
 * Evaluates all background_lqc quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background_lqc structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_lqc_at_tau(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;
// printf("future_branch = %d \n",pba->future_branch);


 if(pba->future_branch == 0){
  class_test(tau < pba->tau_table_pb[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table_pb[0]);

  class_test(tau > pba->tau_table_pb[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table_pb[pba->bt_size-1]);
}
else{
  class_test(tau < pba->tau_table_fb[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table_fb[0]);

  class_test(tau > pba->tau_table_fb[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table_fb[pba->bt_size-1]);

}
  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  
  /*
      if (pba->has_lqc == _FALSE_ && tau>0. && 1<0){


  class_test(tau < pba->tau_table_expanding[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table_expanding[0]);

  class_test(tau > pba->tau_table_expanding[pba->bt_expanding_phase_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table_expanding[pba->bt_expanding_phase_size-1]);


  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  //printf("okokoko");
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table_expanding,
                                        pba->bt_expanding_phase_size,
                                        pba->background_lqc_table_expanding,
                                        pba->d2background_lqc_dtau2_table_expanding,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table_expanding,
                                                        pba->bt_expanding_phase_size,
                                                        pba->background_lqc_table_expanding,
                                                        pba->d2background_lqc_dtau2_table_expanding,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;

  }*/
  
  
 if(pba->future_branch == 0){
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table_pb,
                                        pba->bt_size,
                                        pba->background_lqc_table_pb,
                                        pba->d2background_lqc_dtau2_table_pb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table_pb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_pb,
                                                        pba->d2background_lqc_dtau2_table_pb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  }
else{
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table_fb,
                                        pba->bt_size,
                                        pba->background_lqc_table_fb,
                                        pba->d2background_lqc_dtau2_table_fb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table_fb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_fb,
                                                        pba->d2background_lqc_dtau2_table_fb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }


}

        
  return _SUCCESS_;
}




int background_lqc_at_N(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;

if(pba->future_branch == 0){
  class_test(tau < pba->N_table_pb[0],
             pba->error_message,
             "out of range: N=%e < N_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->N_table_pb[0]);

  class_test(tau > pba->N_table_pb[pba->bt_size-1],
             pba->error_message,
             "out of range: N=%e > N_max=%e\n",tau,pba->N_table_pb[pba->bt_size-1]);
}
else{
  class_test(tau < pba->N_table_fb[0],
             pba->error_message,
             "out of range: N=%e < N_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->N_table_fb[0]);

  class_test(tau > pba->N_table_fb[pba->bt_size-1],
             pba->error_message,
             "out of range: N=%e > N_max=%e\n",tau,pba->N_table_fb[pba->bt_size-1]);

}
  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }


if(pba->future_branch == 0){
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->N_table_pb,
                                        pba->bt_size,
                                        pba->background_lqc_table_pb,
                                        pba->d2background_lqc_dN2_table_pb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->N_table_pb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_pb,
                                                        pba->d2background_lqc_dN2_table_pb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
}
else{
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->N_table_fb,
                                        pba->bt_size,
                                        pba->background_lqc_table_fb,
                                        pba->d2background_lqc_dN2_table_fb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->N_table_fb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_fb,
                                                        pba->d2background_lqc_dN2_table_fb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
}
  /*
              if (pba->has_lqc == _FALSE_ && tau>0.){




  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->N_table_expanding,
                                        pba->bt_expanding_phase_size,
                                        pba->background_lqc_table_expanding,
                                        pba->d2background_lqc_dN2_table_expanding,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);

    printf("ok");
  }
  if (intermode == pba->inter_closeby) {

    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->N_table_expanding,
                                                        pba->bt_expanding_phase_size,
                                                        pba->background_lqc_table_expanding,
                                                        pba->d2background_lqc_dN2_table_expanding,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }



  }*/
  return _SUCCESS_;
}

int background_lqc_tau_of_N(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  //- check that \f$ z \f$ is in the pre-computed range 
   /*class_test(z < pba->N_table[pba->bt_size-1],
             pba->error_message,
             "out of range: N=%e < N_min=%e\n",z,pba->N_table[pba->bt_size-1]);

  class_test(z > pba->N_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->N_table[0]);
  */
  //- interpolate from pre-computed table with array_interpolate()
if(pba->future_branch == 0) {
  class_call(array_interpolate_spline(
                                      pba->N_table_pb,
                                      pba->bt_size,
                                      pba->tau_table_pb,
                                      pba->d2tau_dN2_table_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
              if (pba->has_lqc == _FALSE_ && *tau>0.)
  class_call(array_interpolate_spline(
                                      pba->N_table_expanding_pb,
                                      pba->bt_expanding_phase_size,
                                      pba->tau_table_expanding_pb,
                                      pba->d2tau_dN2_table_expanding_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  }
else{
  class_call(array_interpolate_spline(
                                      pba->N_table_fb,
                                      pba->bt_size,
                                      pba->tau_table_fb,
                                      pba->d2tau_dN2_table_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
              if (pba->has_lqc == _FALSE_ && *tau>0.)
  class_call(array_interpolate_spline(
                                      pba->N_table_expanding_fb,
                                      pba->bt_expanding_phase_size,
                                      pba->tau_table_expanding_fb,
                                      pba->d2tau_dN2_table_expanding_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

}



  
  return _SUCCESS_;
}


int background_lqc_N_of_tau(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  //- check that \f$ z \f$ is in the pre-computed range 
  /*class_test(z < pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: t=%e < t_min=%e\n",z,pba->tau_table[pba->bt_size-1]);

  class_test(z > pba->tau_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->tau_table[0]);
  */
  //- interpolate from pre-computed table with array_interpolate() 
 if(pba->future_branch == 0){

 class_call(array_interpolate_spline(
                                      pba->tau_table_pb,
                                      pba->bt_size,
                                      pba->N_table_pb,
                                      pba->d2N_dtau2_table_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);
}
else{
 class_call(array_interpolate_spline(
                                      pba->tau_table_fb,
                                      pba->bt_size,
                                      pba->N_table_fb,
                                      pba->d2N_dtau2_table_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

}
  return _SUCCESS_;
}


///////////////////////////////////////////////////////////////////////////
//V. Sreenath 280617
/**
 * Background_Lqc quantities at given conformal time tau.
 *
 * Evaluates all background_lqc quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background_lqc structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_lqc_at_tau_pb(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;
// printf("future_branch = %d \n",pba->future_branch);


  class_test(tau < pba->tau_table_pb[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table_pb[0]);

  class_test(tau > pba->tau_table_pb[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table_pb[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  
 
  
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table_pb,
                                        pba->bt_size,
                                        pba->background_lqc_table_pb,
                                        pba->d2background_lqc_dtau2_table_pb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table_pb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_pb,
                                                        pba->d2background_lqc_dtau2_table_pb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }


        
  return _SUCCESS_;
}




int background_lqc_at_N_pb(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;

  class_test(tau < pba->N_table_pb[0],
             pba->error_message,
             "out of range: N=%e < N_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->N_table_pb[0]);

  class_test(tau > pba->N_table_pb[pba->bt_size-1],
             pba->error_message,
             "out of range: N=%e > N_max=%e\n",tau,pba->N_table_pb[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }


 if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->N_table_pb,
                                        pba->bt_size,
                                        pba->background_lqc_table_pb,
                                        pba->d2background_lqc_dN2_table_pb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->N_table_pb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_pb,
                                                        pba->d2background_lqc_dN2_table_pb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;

}

int background_lqc_tau_of_N_pb(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  class_call(array_interpolate_spline(
                                      pba->N_table_pb,
                                      pba->bt_size,
                                      pba->tau_table_pb,
                                      pba->d2tau_dN2_table_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
              if (pba->has_lqc == _FALSE_ && *tau>0.)
  class_call(array_interpolate_spline(
                                      pba->N_table_expanding_pb,
                                      pba->bt_expanding_phase_size,
                                      pba->tau_table_expanding_pb,
                                      pba->d2tau_dN2_table_expanding_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);



  
  return _SUCCESS_;
}


int background_lqc_N_of_tau_pb(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  //- check that \f$ z \f$ is in the pre-computed range 
  /*class_test(z < pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: t=%e < t_min=%e\n",z,pba->tau_table[pba->bt_size-1]);

  class_test(z > pba->tau_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->tau_table[0]);
  */
  //- interpolate from pre-computed table with array_interpolate() 

 class_call(array_interpolate_spline(
                                      pba->tau_table_pb,
                                      pba->bt_size,
                                      pba->N_table_pb,
                                      pba->d2N_dtau2_table_pb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

//V. Sreenath _pb over
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//V. Sreenath 280617
/**
 * Background_Lqc quantities at given conformal time tau.
 *
 * Evaluates all background_lqc quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background_lqc structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_lqc_at_tau_fb(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;
// printf("future_branch = %d \n",pba->future_branch);


  class_test(tau < pba->tau_table_fb[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table_fb[0]);

  class_test(tau > pba->tau_table_fb[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table_fb[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  
 
  
  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table_fb,
                                        pba->bt_size,
                                        pba->background_lqc_table_fb,
                                        pba->d2background_lqc_dtau2_table_fb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table_fb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_fb,
                                                        pba->d2background_lqc_dtau2_table_fb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }


        
  return _SUCCESS_;
}




int background_lqc_at_N_fb(
                      struct background_lqc *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  int pvecback_size;

  class_test(tau < pba->N_table_fb[0],
             pba->error_message,
             "out of range: N=%e < N_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->N_table_fb[0]);

  class_test(tau > pba->N_table_fb[pba->bt_size-1],
             pba->error_message,
             "out of range: N=%e > N_max=%e\n",tau,pba->N_table_fb[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }


 if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->N_table_fb,
                                        pba->bt_size,
                                        pba->background_lqc_table_fb,
                                        pba->d2background_lqc_dN2_table_fb,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->N_table_fb,
                                                        pba->bt_size,
                                                        pba->background_lqc_table_fb,
                                                        pba->d2background_lqc_dN2_table_fb,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;

}

int background_lqc_tau_of_N_fb(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  class_call(array_interpolate_spline(
                                      pba->N_table_fb,
                                      pba->bt_size,
                                      pba->tau_table_fb,
                                      pba->d2tau_dN2_table_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
              if (pba->has_lqc == _FALSE_ && *tau>0.)
  class_call(array_interpolate_spline(
                                      pba->N_table_expanding_fb,
                                      pba->bt_expanding_phase_size,
                                      pba->tau_table_expanding_fb,
                                      pba->d2tau_dN2_table_expanding_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);



  
  return _SUCCESS_;
}


int background_lqc_N_of_tau_fb(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  //- check that \f$ z \f$ is in the pre-computed range 
  /*class_test(z < pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: t=%e < t_min=%e\n",z,pba->tau_table[pba->bt_size-1]);

  class_test(z > pba->tau_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->tau_table[0]);
  */
  //- interpolate from pre-computed table with array_interpolate() 

 class_call(array_interpolate_spline(
                                      pba->tau_table_fb,
                                      pba->bt_size,
                                      pba->N_table_fb,
                                      pba->d2N_dtau2_table_fb,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

//V. Sreenath _fb over
///////////////////////////////////////////////////////////////////////////////



/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background_lqc structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */
/*
int background_lqc_tau_of_z(
                        struct background_lqc *pba,
                        double z,
                        double * tau
                        ) {


  // - define local variables 

  // necessary for calling array_interpolate(), but never used 
  int last_index;

  //- check that \f$ z \f$ is in the pre-computed range 
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  //- interpolate from pre-computed table with array_interpolate() 
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}
*/
/**
 * Background_Lqc quantities at given \f$ a \f$.
 *
 * Function evaluating all background_lqc quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background_lqc structure
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background_lqc quantities (assumed to be already allocated)
 * @return the error status
 */

int background_lqc_functions(
                         struct background_lqc *pba,
                         double * pvecback_B, /* Vector containing all {B} quantities. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;

  double rho_scf;
  double p_scf;
 
  /* scale factor relative to scale factor at the bounce */
  double a_bounce;
  
  /* scale factor */
  //double a;
  /* scalar field quantities */
  double phi, phi_prime;
  double x_scf, y_scf;

  /** - pass value of \f$ a\f$ to output */
  // pvecback[pba->index_bg_a] = a;

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [c^2/8\pi G] \f$, ie*/
 
  pvecback[pba->index_bg_a] = pow(pvecback_B[pba->index_bi_v]*pow(8.*_PI_,3./2.)/4.,1./3.);
  pvecback[pba->index_bg_N] = pvecback_B[pba->index_bi_N];
  pvecback[pba->index_bg_conformal_time] = pvecback_B[pba->index_bi_conformal_time];
  pvecback[pba->index_bg_cosmic_time] = pvecback_B[pba->index_bi_cosmic_time];
  pvecback[pba->index_bg_phi_scf] = pvecback_B[pba->index_bi_phi];
//220517  pvecback[pba->index_bg_phi_prime_scf] = pvecback_B[pba->index_bi_phi_dot];
  pvecback[pba->index_bg_Pphi_scf] = pow(8.*_PI_,3./2.)*
	pvecback_B[pba->index_bi_Pphi];
  pvecback[pba->index_bg_v] = pvecback_B[pba->index_bi_v];
  pvecback[pba->index_bg_h] = pvecback_B[pba->index_bi_h];
//Define potential of inflaton here
//Quadratic Potential
  pvecback[pba->index_bg_V] = pow(pba->m_scf_lqc*pvecback[pba->index_bg_phi_scf],2.)/2.;
  pvecback[pba->index_bg_V_phi] = pow(pba->m_scf_lqc,2.)*pvecback[pba->index_bg_phi_scf];
  pvecback[pba->index_bg_V_phiphi] = pow(pba->m_scf_lqc,2.);
  pvecback[pba->index_bg_V_phiphiphi] = 0.;


  pvecback[pba->index_bg_phi_prime_scf] = 4.*pvecback_B[pba->index_bi_Pphi]/pvecback[pba->index_bg_v];
  pvecback[pba->index_bg_x] = sqrt(pvecback[pba->index_bg_V]/pba->rho_bounce);
  pvecback[pba->index_bg_y] = pvecback[pba->index_bg_phi_prime_scf]/sqrt(2.*pba->rho_bounce);

  pvecback[pba->index_bg_pia] = 6.*pow(pvecback[pba->index_bg_a],2.)*pvecback[pba->index_bg_h];

/*  pvecback[pba->index_bg_pia] = 
	-sqrt( 6.*pow(pvecback[pba->index_bg_Pphi_scf],2.)/pow(pvecback[pba->index_bg_a],2.) 
	+ 12.*pow(pvecback[pba->index_bg_a],4.)*pvecback[pba->index_bg_V]);//Derived from GR constraint*/

/*  if(pba->future_branch == 0){
  pvecback[pba->index_bg_pia] = 
	-sqrt( 6.*pow(pvecback[pba->index_bg_Pphi_scf],2.)/pow(pvecback[pba->index_bg_a],2.) 
	+ 12.*pow(pvecback[pba->index_bg_a],4.)*pvecback[pba->index_bg_V]);//Derived from GR constraint
	}
else{
  pvecback[pba->index_bg_pia] = 
	-sqrt( 6.*pow(pvecback[pba->index_bg_Pphi_scf],2.)/pow(pvecback[pba->index_bg_a],2.) 
	+ 12.*pow(pvecback[pba->index_bg_a],4.)*pvecback[pba->index_bg_V]);//Derived from GR constraint
}*/

/*  if(pvecback[pba->index_bg_phi_scf]<0)
     pvecback[pba->index_bg_x] = -sqrt(pvecback[pba->index_bg_V]/pba->rho_bounce);*/ //V. Sreenath not needed anymore
   /** - compute each component's density and pressure */
 
  rho_tot = 0.;
  p_tot = 0.;
  
  rho_scf = pba->rho_bounce*(
            pvecback[pba->index_bg_x]*pvecback[pba->index_bg_x]
            +pvecback[pba->index_bg_y]*pvecback[pba->index_bg_y]
                             );
    
  p_scf = pba->rho_bounce*(
            pvecback[pba->index_bg_y]*pvecback[pba->index_bg_y]
            -pvecback[pba->index_bg_x]*pvecback[pba->index_bg_x]
                             );

  rho_tot += rho_scf;
  p_tot += p_scf;
  pvecback[pba->index_bg_rho_scf] =rho_scf;
  pvecback[pba->index_bg_p_scf] =p_scf;
 
//  pvecback[pba->index_bg_pia] = -6.*pow(pvecback[pba->index_bg_a],2.)*sqrt(pvecback[pba->index_bg_rho_scf]/3.);//V. Sreenath: redefining pia

  double r_lqc2;
 if (pba->has_lqc == _FALSE_){
   r_lqc2 = 0.;
 }
 else{
   r_lqc2 =
     pow(pvecback[pba->index_bg_x],2.)
     +pow(pvecback[pba->index_bg_y],2.);
 }

/*pvecback[pba->index_bg_H] =(sqrt(pba->rho_bounce*
                                            ( pow(pvecback[pba->index_bg_x],2.)
     +pow(pvecback[pba->index_bg_y],2.)))/sqrt(3.))*sqrt(1.-r_lqc2);

 
  if (pba->post_bounce == _FALSE_) pvecback[pba->index_bg_H] =
                                     -(sqrt(pba->rho_bounce*
                                            ( pow(pvecback[pba->index_bg_x],2.)
     +pow(pvecback[pba->index_bg_y],2.)))/sqrt(3.))*sqrt(1.-r_lqc2);


  if (r_lqc2 == 1.) pvecback[pba->index_bg_H] = 0.;
*///220517
pvecback[pba->index_bg_H] = -sin(2.*pba->l0_lqc*pvecback[pba->index_bg_h])/(2.*pba->l0_lqc);  
//   if (r_lqc2 == 1.) pvecback[pba->index_bg_H] = 0.;
 
  if (pba->remote_past == _TRUE_){


      pvecback[pba->index_bg_H] =
        pba->H_ini_remote_past
        *exp(-(3./2.)
             *fabs(pvecback_B[pba->index_bi_N]
                    -pba->N_ini_remote_past)
             ) ;
      pvecback[pba->index_bg_epsilon_H] =1.5;
      pvecback[pba->index_bg_x] = 0.;
      pvecback[pba->index_bg_y] =0.;
      pvecback[pba->index_bg_z_S] = 0.;
      pvecback[pba->index_bg_rho_scf] = 3.*pow(pvecback[pba->index_bg_H],2.);
  }
  else{
    
    if (pba->post_bounce == _FALSE_){
  pba->H_ini_remote_past = pvecback[pba->index_bg_H];
  pba->N_ini_remote_past = pvecback[pba->index_bg_N];
 }
    
pvecback[pba->index_bg_epsilon_H] =
  (3.*pow(pvecback[pba->index_bg_y],2.))
   /(pow(pvecback[pba->index_bg_x],2.)
     +pow(pvecback[pba->index_bg_y],2.));
 
  pvecback[pba->index_bg_z_S] =
    pvecback[pba->index_bg_a]
    *sqrt(6.)
    *pvecback[pba->index_bg_y]
    /sqrt(pow(pvecback[pba->index_bg_x],2.)
     +pow(pvecback[pba->index_bg_y],2.));
    }
  
  if (fabs(pvecback[pba->index_bg_y])==0.)
    pvecback[pba->index_bg_t] = _HUGE_;
  else  pvecback[pba->index_bg_t] =
          pvecback[pba->index_bg_x]
          /pvecback[pba->index_bg_y];
  
  double V_phiphi = pvecback[pba->index_bg_V_phiphi];
  double V_phi = pvecback[pba->index_bg_V_phi];
  double rho =  pvecback[pba->index_bg_rho_scf];
  double rho_B = pba->rho_bounce;
  double phi_dot = pvecback[pba->index_bg_phi_prime_scf];

/*    pvecback[pba->index_bg_H2_f_S_Q] =
    9./2.*pow(phi_dot,4.)*pow(pvecback[pba->index_bg_H]/rho,2.)
    -9.*phi_dot*phi_dot*pow(pvecback[pba->index_bg_H],2.)/rho
    -6.*pvecback[pba->index_bg_H]
    /rho
    *V_phi
    *phi_dot
    -V_phiphi;*/

//U1
/*    pvecback[pba->index_bg_H2_f_S_Q] = 
	-(-18.*pow(pvecback[pba->index_bg_Pphi_scf],4.)*pow(pvecback[pba->index_bg_a],-8.)*pow(pvecback[pba->index_bg_pia],-2.)
	+ 3.*pow(pvecback[pba->index_bg_Pphi_scf],2.)*pow(pvecback[pba->index_bg_a],-6.)
	- 12.*pvecback[pba->index_bg_Pphi_scf]*pvecback[pba->index_bg_V_phi]/(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_pia])
	+ pvecback[pba->index_bg_V_phiphi]);//This is equal to -U.
*/

//U2
	double Pphi_dot = -pow(8.*_PI_,3./2.)*pvecback[pba->index_bg_v]*pvecback[pba->index_bg_V_phi]/4.;
	pvecback[pba->index_bg_Pphidot_scf] = Pphi_dot;

	double pia_dot = 12.*pow(pvecback[pba->index_bg_a],2.)*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_h]
		+ 6.*pow(pvecback[pba->index_bg_a],2.)*(4.*pow(pvecback_B[pba->index_bi_Pphi]/pvecback[pba->index_bg_v],2.)
		+3.*pow(sin(pba->l0_lqc*pvecback[pba->index_bg_h])/pba->l0_lqc,2.)/2. 
		- pvecback[pba->index_bg_V]/2.);

/*	double pia_dot = ( 12.* pvecback[pba->index_bg_Pphi_scf]*Pphi_dot*pow(pvecback[pba->index_bg_a],-2.)
		- 12.*pow(pvecback[pba->index_bg_Pphi_scf]/pvecback[pba->index_bg_a],2.)*pvecback[pba->index_bg_H]
		+ 48.*pow(pvecback[pba->index_bg_a],4.)*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_V]
		+ 12.*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_V_phi]*pvecback[pba->index_bg_Pphi_scf])
		/(2.*pvecback[pba->index_bg_pia]);*/
	
	pvecback[pba->index_bg_piadot] = pia_dot;

	 pvecback[pba->index_bg_H2_f_S_Q] = 
		-(-9.*pow(pvecback[pba->index_bg_Pphi_scf],4.)*pow(pvecback[pba->index_bg_a],-8.)*pow(pvecback[pba->index_bg_pia],-2.)
		+ 3.*pow(pvecback[pba->index_bg_Pphi_scf],2.)*pow(pvecback[pba->index_bg_a],-6.)/2.
		- 6.*pvecback[pba->index_bg_Pphi_scf]*pvecback[pba->index_bg_V_phi]/(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_pia])
		+ pvecback[pba->index_bg_V_phiphi]
		+ 6.*pvecback[pba->index_bg_Pphi_scf]*Pphi_dot*pow(pvecback[pba->index_bg_a],-4.)/pvecback[pba->index_bg_pia]
		- 3.*pow(pvecback[pba->index_bg_Pphi_scf],2.)*pvecback[pba->index_bg_H]*pow(pvecback[pba->index_bg_a],-4.)/pvecback[pba->index_bg_pia]
		-  3.*pow(pvecback[pba->index_bg_Pphi_scf],2.)*pia_dot*pow(pvecback[pba->index_bg_a],-4.)*pow(pvecback[pba->index_bg_pia],-2.));

//U3
/*	double f = 3.*pow(pvecback[pba->index_bg_Pphi_scf],2.)/(pow(pvecback[pba->index_bg_Pphi_scf],2.)/2. + pow(pvecback[pba->index_bg_a],6.)*pvecback[pba->index_bg_V]);
	pvecback[pba->index_bg_H2_f_S_Q] = -(f*pvecback[pba->index_bg_V] - 2.*sqrt(f)*pvecback[pba->index_bg_V_phi] + pvecback[pba->index_bg_V_phiphi]);
*/

    pvecback[pba->index_bg_H2_f_S_Q_phys] = 
	-(-3./2.*pow(phi_dot,4.)/rho*(1. - rho/rho_B) + 3.*pow(phi_dot,2.) + 6.*V_phi*phi_dot*pvecback[pba->index_bg_H]/rho + V_phiphi);

//    pvecback[pba->index_bg_H2_f_S_Q] = pvecback[pba->index_bg_H2_f_S_Q_phys]; 

    pvecback[pba->index_bg_zT_primeprime_over_zT] =
        pvecback[pba->index_bg_a]
        *pvecback[pba->index_bg_a]
        *(pvecback[pba->index_bg_H]
        *pvecback[pba->index_bg_H]
        + (- pvecback[pba->index_bg_rho_scf]/6.*(1. - 4.*pvecback[pba->index_bg_rho_scf]/pba->rho_bounce) 
	   - pvecback[pba->index_bg_p_scf]/2.*(1. - 2.*pvecback[pba->index_bg_rho_scf]/pba->rho_bounce) ) 
	);

/*   pvecback[pba->index_bg_zT_primeprime_over_zT] =
        pvecback[pba->index_bg_a]
        *pvecback[pba->index_bg_a]
        *pvecback[pba->index_bg_H]
        *pvecback[pba->index_bg_H]
        *(2-pvecback[pba->index_bg_epsilon_H]);*/


    pvecback[pba->index_bg_zS_primeprime_over_zS] =
        +pvecback[pba->index_bg_a]
        *pvecback[pba->index_bg_a]
        *pvecback[pba->index_bg_H2_f_S_Q]
        +pvecback[pba->index_bg_zT_primeprime_over_zT];

  
  pvecback[pba->index_bg_z_T] = pvecback[pba->index_bg_a];
      
  pvecback[pba->index_bg_H_prime] = -pow(pvecback[pba->index_bg_H],2.)
                                    *pvecback[pba->index_bg_epsilon_H]
                                    *pvecback[pba->index_bg_a];
pvecback[pba->index_bg_epsilon_2] = 
  -6.
  -2.
  *pvecback[pba->index_bg_V_phi]
  /pvecback[pba->index_bg_H]
  /pvecback[pba->index_bg_phi_prime_scf]
  +2.*pvecback[pba->index_bg_epsilon_H];
   
  if(pba->remote_past == _TRUE_){
	pvecback[pba->index_bg_H2_f_S_Q] = 0.;
	pvecback[pba->index_bg_pia] = 0.;
	pvecback[pba->index_bg_piadot] = 0.;
   }
                    
  return _SUCCESS_;

}

/**
 * Initialize the background_lqc structure, and in particular the
 * background_lqc interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background_lqc structure
 * @return the error status
 */

int background_lqc_init(
                    struct precision * ppr,
                    struct background_lqc * pba
                    ) {

  /** Summary: */

  /** - define local variables */
  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  double Neff;
  int filenum=0;
  pba->background_lqc_verbose =10;

  /** - in verbose mode, provide some information */
  if (pba->background_lqc_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");
    printf("Mass of the Scalar Field: m = %e M_Pl\n",pba->m_scf_lqc);

 if (pba->has_lqc == _FALSE_){
   printf("GR Dynamics\n");
 }
 else{
   printf("LQC Dynamics\n");
   printf("l0 = %e\n",pba->l0_lqc);
   printf("rho_bounce = %e\n",pba->rho_bounce);
}
 
 }

  
  if (pba->has_contracting_phase==_TRUE_) {
  pba->post_bounce = _FALSE_;

     }


  class_call(background_lqc_indices(pba),
             pba->error_message,
             pba->error_message);

  class_call(background_lqc_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_lqc_init().
 *
 *
 * @param pba Input: pointer to background_lqc structure (to be freed)
 * @return the error status
 */

int background_lqc_free(
                    struct background_lqc *pba
                    ) {
  int err;

  free(pba->tau_table);
  free(pba->N_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->d2tau_dN2_table);
  free(pba->d2N_dtau2_table);
  free(pba->background_lqc_table);
  free(pba->d2background_lqc_dtau2_table);
  free(pba->d2background_lqc_dN2_table);
  
 free(pba->tau_table_expanding);
  free(pba->N_table_expanding);
  free(pba->d2tau_dN2_table_expanding);
  free(pba->d2N_dtau2_table_expanding);
  free(pba->background_lqc_table_expanding);
  free(pba->d2background_lqc_dtau2_table_expanding);
  free(pba->d2background_lqc_dN2_table_expanding);

  free(pba->tau_table_pb);
  free(pba->N_table_pb);
  free(pba->z_table_pb);
  free(pba->d2tau_dz2_table_pb);
  free(pba->d2tau_dN2_table_pb);
  free(pba->d2N_dtau2_table_pb);
  free(pba->background_lqc_table_pb);
  free(pba->d2background_lqc_dtau2_table_pb);
  free(pba->d2background_lqc_dN2_table_pb);
  
  free(pba->tau_table_expanding_pb);
  free(pba->N_table_expanding_pb);
  free(pba->d2tau_dN2_table_expanding_pb);
  free(pba->d2N_dtau2_table_expanding_pb);
  free(pba->background_lqc_table_expanding_pb);
  free(pba->d2background_lqc_dtau2_table_expanding_pb);
  free(pba->d2background_lqc_dN2_table_expanding_pb);

  free(pba->tau_table_fb);
  free(pba->N_table_fb);
  free(pba->z_table_fb);
  free(pba->d2tau_dz2_table_fb);
  free(pba->d2tau_dN2_table_fb);
  free(pba->d2N_dtau2_table_fb);
  free(pba->background_lqc_table_fb);
  free(pba->d2background_lqc_dtau2_table_fb);
  free(pba->d2background_lqc_dN2_table_fb);
  
  free(pba->tau_table_expanding_fb);
  free(pba->N_table_expanding_fb);
  free(pba->d2tau_dN2_table_expanding_fb);
  free(pba->d2N_dtau2_table_expanding_fb);
  free(pba->background_lqc_table_expanding_fb);
  free(pba->d2background_lqc_dtau2_table_expanding_fb);
  free(pba->d2background_lqc_dN2_table_expanding_fb);
 
 err = background_lqc_free_input(pba);

  return err;
}

/**
 * Free pointers inside background_lqc structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input: pointer to background_lqc structure
 * @return the error status
 */

int background_lqc_free_input(
                          struct background_lqc *pba
                          ) {

  int k;


  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background_lqc quantities.
 *
 * @param pba Input: pointer to background_lqc structure
 * @return the error status
 */

int background_lqc_indices(
                       struct background_lqc *pba
                       ) {


  int index_bg;
  int index_bi;


  pba->has_scf = _FALSE_;


  if (pba->Omega0_scf != 0.)
    pba->has_scf = _TRUE_;
  

  index_bg=0;

  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_x,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_y,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_mu,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_t_times_epsilon_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_t,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_V,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_V_phi,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_V_phiphi,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_V_phiphiphi,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_epsilon_H,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_epsilon_V,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_epsilon_2,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_f_S_Q,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H2_f_S_Q,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H2_f_S_Q_phys,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_f_T_Q,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_f_S_R,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_f_T_R,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_z_S,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_z_T,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_f_S,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_f_T,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_N,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_conformal_time,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_cosmic_time,_TRUE_,index_bg,1);
//  class_define_index(pba->index_bg_effective_potential,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_zS_primeprime_over_zS,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_zT_primeprime_over_zT,_TRUE_,index_bg,1);

  pba->bg_size_short = index_bg;
  class_define_index(pba->index_bg_phi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_phi_prime_scf,pba->has_scf,index_bg,1);

  class_define_index(pba->index_bg_rho_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_scf,pba->has_scf,index_bg,1);

  class_define_index(pba->index_bg_Pphi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_Pphidot_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_v,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_h,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_pia,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_piadot,pba->has_scf,index_bg,1);

  pba->bg_size_normal = index_bg;
  pba->bg_size = index_bg;

  index_bi=0;

  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_conformal_time,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_N,_TRUE_,index_bi,1);
//  class_define_index(pba->index_bi_x,pba->has_scf,index_bi,1);
//  class_define_index(pba->index_bi_y,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_phi,pba->has_scf,index_bi,1);
//  class_define_index(pba->index_bi_phi_dot,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_Pphi,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_v,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_h,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_H,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_cosmic_time,pba->has_scf,index_bi,1);

  pba->bi_B_size = index_bi;

  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  pba->bi_size = index_bi;

  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background_lqc integration requires index_bi_tau to be the last of all index_bi's");


  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 *  This function integrates the background_lqc over time, allocates and
 *  fills the background_lqc table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background_lqc structure
 */

int background_lqc_solve(
                     struct precision *ppr,
                     struct background_lqc *pba
                     ) {


  struct generic_integrator_workspace gi;
  struct background_lqc_parameters_and_workspace bpaw;
  growTable gTable;
  double * pData;
  void * memcopy_result;
  void * memcopy_result_expanding;
  double tau_start;
  double tau_end;
  int i;
  int j;
  int index_bi;
  double temp;
  int index_tau;
  double * pvecback_integration;
  double * pvecback;
  int last_index=0;


  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

/////////////
  pba->future_branch = 0;
  while(pba->future_branch <=1){
  printf("branch = %d\n",pba->future_branch);

  pba->post_bounce = _FALSE_;

  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  class_call(background_lqc_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);


  tau_end=pvecback_integration[pba->index_bi_tau];

  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  pba->bt_size=0;
  pba->bt_contracting_phase_size=0;
  pba->bt_remote_past_size=0;

   class_call(background_lqc_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);//180517

  if (pba->has_contracting_phase==_TRUE_) {
    short first_step = _TRUE_;
      double aH = HUGE;

      /*   while ((pow(pvecback_integration[pba->index_bi_x],2.)
            + pow(pvecback_integration[pba->index_bi_y],2.)) >
            (1./pow(pba->ALPHA_end*pba->Gamma,2)))*/    
      while ( pba->k_min_for_pk< pba->integration_depth_ini*fabs(aH) || (pow(pvecback[pba->index_bg_x],2.)
       + pow(pvecback[pba->index_bg_y],2.)) >
       (1./pow(pba->ALPHA_end*pba->Gamma,2)))      
    {
      if((pow(pvecback[pba->index_bg_x],2.)
            + pow(pvecback[pba->index_bg_y],2.)) >
           (1./pow(pba->ALPHA_end*pba->Gamma,2)))  
        pba->remote_past =_FALSE_;
        
      else  pba->remote_past =_TRUE_;
      /*
      printf("a=%e H=%e x=%e y=%e\n",pvecback_integration[pba->index_bi_a],
             pvecback[pba->index_bg_H],
             pvecback_integration[pba->index_bi_x],
             pvecback_integration[pba->index_bi_y]);

      */
      /*       
   if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>0.99) ppr->back_lqc_integration_stepsize = 1.e-4;
   else if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-3;
   else ppr->back_lqc_integration_stepsize = 1.e-2;

      */
        
/*   if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>.999) ppr->back_lqc_integration_stepsize = 1.e-6;
   else if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-4;
   else ppr->back_lqc_integration_stepsize = 1.e-2;*/
     
   if (pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2)>.999) ppr->back_lqc_integration_stepsize = 1.e-6;
   else if (pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-4;
   else ppr->back_lqc_integration_stepsize = 1.e-2;
   

      /*
      double timescale;
    timescale=  1.e-4*MIN(
                 (1./pba->m_scf_lqc),
                 fabs(1./pvecback[pba->index_bg_H]));
if (pvecback_integration[pba->index_bi_y]>0.)
  timescale = MIN(timescale,1.e-2/sqrt(pba->rho_bounce));
    ppr->back_lqc_integration_stepsize = 1.e-2*timescale;

    printf("timescale=%e\n",timescale);
      */
    tau_start = tau_end;
    
    class_call(background_lqc_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);
    tau_end = tau_start + ppr->back_lqc_integration_stepsize;
    aH = fabs(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a]);

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    if (first_step == _FALSE_){
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
  pba->bt_size++;
  pba->bt_contracting_phase_size++;
  if (pba->remote_past ==_TRUE_) pba->bt_remote_past_size++;


    }
    first_step = _FALSE_;
    /* -> perform one step */
    class_call(generic_integrator(background_lqc_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_lqc_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;


    
  }

    
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;
  pba->bt_contracting_phase_size++;


  }//has_contracting



  pba->remote_past =_FALSE_;
  pba->post_bounce = _TRUE_;
  pba->ALPHA_end = 10.;
  class_call(background_lqc_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);
  tau_end=pvecback_integration[pba->index_bi_tau];
   class_call(background_lqc_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);//180517

  /*
  while ((pow(pvecback_integration[pba->index_bi_x],2.)
         + pow(pvecback_integration[pba->index_bi_y],2.)) >
         (1./pow(pba->ALPHA_end*pba->Gamma,2))) {

        
   if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>.999) ppr->back_lqc_integration_stepsize = 1.e-6;
   else if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-4;
   else ppr->back_lqc_integration_stepsize = 1.e-2;*/
/*   while ((pow(pvecback_integration[pba->index_bi_x],2.)
         + pow(pvecback_integration[pba->index_bi_y],2.)) >
         (1./pow(pba->ALPHA_end*pba->Gamma,2))) {

        
   if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>.999) ppr->back_lqc_integration_stepsize = 1.e-6;
   else if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-4;
   else ppr->back_lqc_integration_stepsize = 1.e-2;*/
    
   while ((pow(pvecback[pba->index_bg_x],2.)
         + pow(pvecback[pba->index_bg_y],2.)) >
         (1./pow(pba->ALPHA_end*pba->Gamma,2))) {

        
   if (pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2)>.999) ppr->back_lqc_integration_stepsize = 1.e-6;
   else if (pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-4;
   else ppr->back_lqc_integration_stepsize = 1.e-2;

   /*
      if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>0.99) ppr->back_lqc_integration_stepsize = 1.e-4;
   else if (pow(pvecback_integration[pba->index_bi_x],2)+pow(pvecback_integration[pba->index_bi_y],2)>1./2.)
   ppr->back_lqc_integration_stepsize = 1.e-3;
   else ppr->back_lqc_integration_stepsize = 1.e-2;
   */
      /*
    double timescale;
    timescale=  1.e-2*MIN(
                 (1./pba->m_scf_lqc),
                 fabs(1./pvecback[pba->index_bg_H]));
    if(pvecback_integration[pba->index_bi_y]>0.)
    timescale = MIN(timescale,1.e-2/sqrt(pba->rho_bounce));
    ppr->back_lqc_integration_stepsize =1.e-2*timescale;

        
      
      printf("a=%e H=%e x=%e y=%e\n",pvecback_integration[pba->index_bi_a],
             pvecback[pba->index_bg_H],
             pvecback_integration[pba->index_bi_x],
             pvecback_integration[pba->index_bi_y]);
    printf("timescale=%e\n",timescale);
      */
      tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_lqc_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);
    tau_end = tau_start + ppr->back_lqc_integration_stepsize;


    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_lqc_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_lqc_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


 
  

  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  //Reverse the table for indices corresponding to contracting phase
  
    for (index_bi=0; index_bi<pba->bi_size; index_bi++){

      // for (index_tau=0; index_tau<pba->bt_contracting_phase_size; index_tau++){
      //printf("ok tau = %d\n",index_tau);

   j = pba->bt_contracting_phase_size-1;   // j will Point to last Element
   index_tau = 0;       // i will be pointing to first element
 
   while (index_tau < j) {
      temp = pData[index_tau*pba->bi_size+index_bi];
      pData[index_tau*pba->bi_size+index_bi] = pData[j*pba->bi_size+index_bi];
      pData[j*pba->bi_size+index_bi] = temp;
      index_tau++;             // increment i
      j--;          // decrement j

   }
   
   
   }

    //for (index_tau=0; index_tau<pba->bt_contracting_phase_size; index_tau++){
    for (index_tau=0; index_tau<pba->bt_size; index_tau++){
      pData[index_tau*pba->bi_size+pba->index_bi_tau]=-pData[index_tau*pba->bi_size+pba->index_bi_tau];
      pData[index_tau*pba->bi_size+pba->index_bi_tau]=pData[index_tau*pba->bi_size+pba->index_bi_cosmic_time];

    if(pba->future_branch == 0){  
      if (index_tau == 0){
        pba->tau_start_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_tau];
        pba->a_start_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_a];
        pba->N_start_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_N];

         }
            if (index_tau == pba->bt_remote_past_size-1){
        pba->tau_end_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_tau];
        pba->a_end_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_a];
        pba->N_end_remote_past =
          pData[index_tau*pba->bi_size+pba->index_bi_N];
         }}
  }

 if(pba->future_branch == 0){
 i = pba->bt_contracting_phase_size-1;
 int index_start_deflation;
 int index_end_deflation;
 while (pba->flag_end_deflation == _FALSE_) {


    pba->post_bounce = _FALSE_;
/*  if (pData[i*pba->bi_size+pba->index_bi_y]>0. &&
      pba->flag_end_deflation == _FALSE_ &&
      pba->flag_start_deflation == _TRUE_  &&
      pba->post_bounce == _FALSE_
      ){*/
  if (pData[i*pba->bi_size+pba->index_bi_Pphi]>0. &&
      pba->flag_end_deflation == _FALSE_ &&
      pba->flag_start_deflation == _TRUE_  &&
      pba->post_bounce == _FALSE_
      ){
    pba->tau_start_deflation =pData[i*pba->bi_size+pba->index_bi_tau];
    pba->N_start_deflation=pData[i*pba->bi_size+pba->index_bi_N];
    pba->a_start_deflation=pData[i*pba->bi_size+pba->index_bi_a];
    index_start_deflation = i;
    pba->flag_end_deflation = _TRUE_;
    }


/*  if (pData[i*pba->bi_size+pba->index_bi_y]<0. &&
      pba->flag_start_deflation == _FALSE_ &&
      pba->post_bounce == _FALSE_
      ){*/
  if (pData[i*pba->bi_size+pba->index_bi_Pphi]<0. &&
      pba->flag_start_deflation == _FALSE_ &&
      pba->post_bounce == _FALSE_
      ){

    pba->tau_end_deflation =pData[i*pba->bi_size+pba->index_bi_tau];
 pba->N_end_deflation=pData[i*pba->bi_size+pba->index_bi_N];
 pba->a_end_deflation=pData[i*pba->bi_size+pba->index_bi_a];
 index_end_deflation = i;
 pba->flag_start_deflation = _TRUE_;
    }
  i-=1;
  }

}
else{
  for (i=0; i < pba->bt_size; i++) {

    if (i<pba->bt_contracting_phase_size) pba->post_bounce = _FALSE_;
    else pba->post_bounce=_TRUE_;

    if (i<pba->bt_remote_past_size) pba->remote_past = _TRUE_;
    else pba->remote_past=_FALSE_;
    
    
/*  if (pData[i*pba->bi_size+pba->index_bi_y]<0. &&
      pba->flag_start_inflation == _FALSE_ &&
      pba->post_bounce == _TRUE_
      ){*/
  if (pData[i*pba->bi_size+pba->index_bi_Pphi]<0. &&
      pba->flag_start_inflation == _FALSE_ &&
      pba->post_bounce == _TRUE_
      ){
    pba->tau_start_inflation = pData[i*pba->bi_size+pba->index_bi_tau];
//pba->tau_table[i];
    pba->N_start_inflation =pData[i*pba->bi_size+pba->index_bi_N];
    pba->a_start_inflation =pData[i*pba->bi_size+pba->index_bi_a];
    pba->flag_start_inflation = _TRUE_;
    
    }




  
/*  if (pData[i*pba->bi_size+pba->index_bi_y]>0. &&
      pba->flag_end_inflation == _FALSE_ &&
      pba->flag_start_inflation == _TRUE_  &&
      pba->post_bounce == _TRUE_
      ){*/
  if (pData[i*pba->bi_size+pba->index_bi_Pphi]>0. &&
      pba->flag_end_inflation == _FALSE_ &&
      pba->flag_start_inflation == _TRUE_  &&
      pba->post_bounce == _TRUE_
      ){
/*    pba->tau_end_inflation =
       pba->tau_start_inflation
       +pData[i*pba->bi_size+pba->index_bi_N]
       -pba->N_start_inflation;*/
    pba->tau_end_inflation = pData[i*pba->bi_size+pba->index_bi_tau];
//pba->tau_table[i];

    pba->N_end = pData[i*pba->bi_size+pba->index_bi_N];
    pba->a_end_inflation = pData[i*pba->bi_size+pba->index_bi_a];
    pba->N_exit = pData[i*pba->bi_size+pba->index_bi_N]-pba->N_star;

    //printf("Inflation ends at N=%e\n", pData[i*pba->bi_size+pba->index_bi_N]);
    //printf("Inflation ends at t=%e\n", pData[i*pba->bi_size+pba->index_bi_tau]);
    pba->flag_end_inflation = _TRUE_;
    }

  }
}

  
    printf("contracting %d\n",pba->bt_contracting_phase_size);
    pba->bt_expanding_phase_size =
      pba->bt_size
      -pba->bt_contracting_phase_size-1;
    printf("expanding %d\n",pba->bt_expanding_phase_size);
  


  if(pba->future_branch == 0){
  class_alloc(pba->tau_table_pb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table_pb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table_pb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table_pb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table_pb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table_pb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table_pb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->tau_table_expanding_pb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table_expanding_pb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table_expanding_pb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table_expanding_pb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table_expanding_pb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table_expanding_pb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table_expanding_pb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);
 
  for (i=0; i < pba->bt_size; i++) {

    pba->tau_table_pb[i] = pData[i*pba->bi_size+pba->index_bi_tau];
    pba->N_table_pb[i] = pData[i*pba->bi_size+pba->index_bi_N];

    class_call(background_lqc_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

//-> write in the table 
     memcopy_result = memcpy(pba->background_lqc_table_pb + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_lqc_table_pb + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");
  }
  
  for (i=0; i < pba->bt_expanding_phase_size; i++) {

  pba->post_bounce=_TRUE_;
  pba->remote_past=_FALSE_;
   
  int j = i + pba->bt_contracting_phase_size;
  
    pba->tau_table_expanding_pb[i] = pData[j*pba->bi_size+pba->index_bi_tau];
    pba->N_table_expanding_pb[i] = pData[j*pba->bi_size+pba->index_bi_N];

    class_call(background_lqc_functions(pba,pData+j*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

 
  
//-> write in the table
/*
     memcopy_result_expanding = memcpy(pba->background_lqc_table_expanding + j*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result_expanding != pba->background_lqc_table_expanding + j*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");*/
  }

  // printf("tau_1=%e\n",pba->tau_table_expanding[1]);

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);



  class_call(array_spline_table_lines(pba->N_table_pb,
                                      pba->bt_size,
                                      pba->tau_table_pb,
                                      1,
                                      pba->d2tau_dN2_table_pb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
  class_call(array_spline_table_lines(pba->tau_table_pb,
                                      pba->bt_size,
                                      pba->N_table_pb,
                                      1,
                                      pba->d2N_dtau2_table_pb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table_pb,
                                      pba->bt_size,
                                      pba->background_lqc_table_pb,
                                      pba->bg_size,
                                      pba->d2background_lqc_dtau2_table_pb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);


  class_call(array_spline_table_lines(pba->N_table_pb,
                                      pba->bt_size,
                                      pba->background_lqc_table_pb,
                                      pba->bg_size,
                                      pba->d2background_lqc_dN2_table_pb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);
 }
  else{
  class_alloc(pba->tau_table_fb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table_fb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table_fb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table_fb,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table_fb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table_fb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table_fb,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->tau_table_expanding_fb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table_expanding_fb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table_expanding_fb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table_expanding_fb,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table_expanding_fb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table_expanding_fb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table_expanding_fb,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);


  for (i=0; i < pba->bt_size; i++) {

    pba->tau_table_fb[i] = pData[i*pba->bi_size+pba->index_bi_tau];
    pba->N_table_fb[i] = pData[i*pba->bi_size+pba->index_bi_N];

    class_call(background_lqc_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

//-> write in the table 
     memcopy_result = memcpy(pba->background_lqc_table_fb + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_lqc_table_fb + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");
  }
  
  for (i=0; i < pba->bt_expanding_phase_size; i++) {

  pba->post_bounce=_TRUE_;
  pba->remote_past=_FALSE_;
   
  int j = i + pba->bt_contracting_phase_size;
  
    pba->tau_table_expanding_fb[i] = pData[j*pba->bi_size+pba->index_bi_tau];
    pba->N_table_expanding_fb[i] = pData[j*pba->bi_size+pba->index_bi_N];

    class_call(background_lqc_functions(pba,pData+j*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

 
  
//-> write in the table
/*
     memcopy_result_expanding = memcpy(pba->background_lqc_table_expanding + j*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result_expanding != pba->background_lqc_table_expanding + j*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");*/
  }

  // printf("tau_1=%e\n",pba->tau_table_expanding[1]);

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);



  class_call(array_spline_table_lines(pba->N_table_fb,
                                      pba->bt_size,
                                      pba->tau_table_fb,
                                      1,
                                      pba->d2tau_dN2_table_fb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
  class_call(array_spline_table_lines(pba->tau_table_fb,
                                      pba->bt_size,
                                      pba->N_table_fb,
                                      1,
                                      pba->d2N_dtau2_table_fb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table_fb,
                                      pba->bt_size,
                                      pba->background_lqc_table_fb,
                                      pba->bg_size,
                                      pba->d2background_lqc_dtau2_table_fb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);


  class_call(array_spline_table_lines(pba->N_table_fb,
                                      pba->bt_size,
                                      pba->background_lqc_table_fb,
                                      pba->bg_size,
                                      pba->d2background_lqc_dN2_table_fb,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  }
 



  //////////EXPANDING SPLINE

  /*
  class_call(array_spline_table_lines(pba->N_table_expanding,
                                      pba->bt_expanding_phase_size,
                                      pba->tau_table_expanding,
                                      1,
                                      pba->d2tau_dN2_table_expanding,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  
  class_call(array_spline_table_lines(pba->tau_table_expanding,
                                      pba->bt_expanding_phase_size,
                                      pba->N_table_expanding,
                                      1,
                                      pba->d2N_dtau2_table_expanding,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table_expanding,
                                      pba->bt_expanding_phase_size,
                                      pba->background_lqc_table_expanding,
                                      pba->bg_size,
                                      pba->d2background_lqc_dtau2_table_expanding,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);


  class_call(array_spline_table_lines(pba->N_table_expanding,
                                      pba->bt_expanding_phase_size,
                                      pba->background_lqc_table_expanding,
                                      pba->bg_size,
                                      pba->d2background_lqc_dN2_table_expanding,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);


 */ 
pba->future_branch +=1;
}//loop over future_branch ends

  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);



  class_alloc(pba->tau_table_expanding,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->N_table_expanding,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2tau_dN2_table_expanding,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2N_dtau2_table_expanding,pba->bt_expanding_phase_size * sizeof(double),pba->error_message);
  class_alloc(pba->background_lqc_table_expanding,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_lqc_dtau2_table_expanding,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_lqc_dN2_table_expanding,pba->bt_expanding_phase_size * pba->bg_size * sizeof(double),pba->error_message);

  for(i=0;i<pba->bt_size;i++){
  
   if(i < pba->bt_contracting_phase_size-1) {
    pba->tau_table[i] = pba->tau_table_pb[i];
    pba->N_table[i] = pba->N_table_pb[i];

   memcopy_result = memcpy(pba->background_lqc_table + i*pba->bg_size,
        pba->background_lqc_table_pb + i*pba->bg_size,pba->bg_size*sizeof(double));
     class_test(memcopy_result != pba->background_lqc_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");
//   printf(" t = %e, a = %e \n",pba->background_lqc_table[i*pba->bg_size+pba->index_bg_cosmic_time],pba->background_lqc_table[i*pba->bg_size+pba->index_bg_a]);
   }
   else{
     pba->tau_table[i] = pba->tau_table_fb[i];
    pba->N_table[i] = pba->N_table_fb[i];

  memcopy_result = memcpy(pba->background_lqc_table + i*pba->bg_size,
        pba->background_lqc_table_fb + i*pba->bg_size,pba->bg_size*sizeof(double));
     class_test(memcopy_result != pba->background_lqc_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_lqc_table");
//   printf(" t = %e, a = %e \n",pba->background_lqc_table[i*pba->bg_size+pba->index_bg_cosmic_time],pba->background_lqc_table[i*pba->bg_size+pba->index_bg_a]);

   }
}

  //////END
  if (pba->background_lqc_verbose > 0) {
    printf("LQC background integration successfull\n");
    printf("Remote past starts at t = %e\n",pba->tau_start_remote_past);
    printf("Remote past ends at t = %e\n",pba->tau_end_remote_past);
    printf("Deflation starts at t=%e\n",
           pba->tau_start_deflation);
    printf("Deflation ends at t=%e\n",
           pba->tau_end_deflation);
    printf("Inflation starts at t=%e\n",
           pba->tau_start_inflation);
    printf("Inflation starts at N=%e\n",
           pba->N_start_inflation);
    printf("Inflation ends at t=%e\n",
           pba->tau_end_inflation);
     printf("Inflation ends at N=%e\n",
           pba->N_end);
 

    
  }
if (pba->has_perturbations_lqc == _FALSE_) 
  return _SUCCESS_;
 ////////////////////////////////////////////////////////////040617
  pba->future_branch = 1;
  double ti= 1.e5;
   double tf = pba->tau_table[pba->bt_size-1];
   double tm = (ti+tf)/2.;
  int index_back;
  double * pvecback2;
   class_alloc(pvecback2,pba->bg_size_normal*sizeof(double),pba->error_message);
   class_call(background_lqc_at_tau(pba,
                                 tm,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
   printf("Value of scalar field at the bounce = %e, corresponding value of potential = %e\n", pba->phi_ini,pow(pba->m_scf_lqc*pba->phi_ini,2.)/2.);
   while(fabs(pvecback2[pba->index_bg_epsilon_H] - 1.) > 1.e-3){
      printf(" ti = %e, tf = %e, delta_epsion_H = %e\n",ti, tf, fabs(pvecback2[pba->index_bg_epsilon_H] - 1.));
      if(pvecback2[pba->index_bg_epsilon_H]>1.)tf = tm;
      else ti = tm;
      tm = (ti+tf)/2.;
      class_alloc(pvecback2,pba->bg_size_normal*sizeof(double),pba->error_message);
      class_call(background_lqc_at_tau(pba,
                                 tm,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
   }
   printf(" Inflation ends at cosmic time = %e tpl, efolds = %e.\n",pvecback2[pba->index_bg_cosmic_time], pvecback2[pba->index_bg_N]);
   printf(" Value of scalar field at the end of inflation = %ecorresponding potential = %e.\n", pvecback2[pba->index_bg_phi_scf], pvecback2[pba->index_bg_V]);
   printf(" Value of slow roll parameter at the end of inflation = %e.\n",pvecback2[pba->index_bg_epsilon_H]);
/////////////////////
    class_call(background_lqc_at_tau(pba,
                                 0.,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
   printf("k_LQC at the bounce = %e,%e \n",sqrt(fabs(pvecback2[pba->index_bg_zT_primeprime_over_zT])) ,sqrt(fabs(-pba->rho_bounce*pow(pvecback2[pba->index_bg_y],2)*(1.-2.*(pow(pvecback2[pba->index_bg_x],2)+pow(pvecback2[pba->index_bg_y],2))))));
   printf("z''/z at the bounce = %e \n", pvecback2[pba->index_bg_zS_primeprime_over_zS]);

/////////////////////////////////////////////////////////////////

  /*find window of wavenumbers*/

  
  double tau_max=pba->tau_end_inflation;
    //pba->tau_start_inflation
    //+pba->N_end
    //-pba->N_start_inflation;
  double tau_min=pba->tau_start_inflation;

  //printf("taumin=%e, taumax=%e\n",tau_min,tau_max);
  
  double tau_half;
  double N_sample=_HUGE_;
  double error_N = .001;
//  int index_back;
//    double * pvecback2;

//    class_alloc(pvecback2,pba->bg_size_normal*sizeof(double),pba->error_message);


    /*
  while (fabs(pba->N_exit-N_sample)>error_N*pba->N_exit){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
  if (pvecback2[pba->index_bg_N]>pba->N_exit) tau_max = tau_half;
  else tau_min = tau_half;
  pba->k_star = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  N_sample = pvecback2[pba->index_bg_N];

 }

*/


  class_call(background_lqc_at_tau(pba,
                                 pba->tau_start_inflation,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);




   printf("As = %e\n",pba->A_s);

    double pk_s_slow_roll=HUGE;

    tau_max=pba->tau_end_inflation;
    //pba->tau_start_inflation;
   //+pba->N_end
   //-pba->N_start_inflation;
 tau_min=pba->tau_start_inflation;

    
  while (fabs(pba->A_s-pk_s_slow_roll)>error_N*pba->A_s){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);


 pk_s_slow_roll = 
  pvecback2[pba->index_bg_H]
   *pvecback2[pba->index_bg_H]
   /pvecback2[pba->index_bg_epsilon_H]
   /(8.*_PI_*_PI_);

   
  if ( pk_s_slow_roll<pba->A_s) tau_max = tau_half;
  else tau_min = tau_half;
  pba->k_star = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  pba->N_star = pba->N_end - pvecback2[pba->index_bg_N];
  pba->N_exit = pvecback2[pba->index_bg_N];
 }
    
  printf("Pivot scale: k_pivot = %e M_Pl\n",pba->k_star);
  printf("Observable Inflation, N_star =%e\n",pba->N_star);
  printf("Cosmic time at the hubble exit of pivot scale = %e\n", pvecback2[pba->index_bg_cosmic_time]);
  printf("Value of scalar field at the hubble exit of pivot scale = %e, corresponding value of potential = %e.\n", pvecback2[pba->index_bg_phi_scf],pvecback2[pba->index_bg_V]);
  
    double k_min_for_pk_s_slow_roll =1.e-3*pba->k_star;
    //1.e2*pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
    double k_max_for_pk_s_slow_roll = 1.e3*pba->k_star;
  printf("k_min_for_pk_s_slow_roll=%e, k_star=%e\n",k_min_for_pk_s_slow_roll,pba->k_star);

  
  FILE * pkfile;     
  FileName file_name;
  double k_per_decade_for_pk_slow_roll = 10.;
  int colnum = 1;


 double one_pk;
 double one_k;
 double error_aH = 0.001;
 double k;
 ///////////////
 k= 6.*pba->k_star;//V. Sreenath

 double aH_sample=_HUGE_;
   tau_max=pba->tau_end_inflation;
   tau_min=pba->tau_start_inflation;

   while (fabs(aH_sample-k)>
         error_aH*k){

    tau_half = (tau_max+tau_min)/2;

  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);

    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  if (aH_sample>k) tau_max = tau_half;
  else tau_min = tau_half;
 }

	double H =pvecback2[pba->index_bg_H];
  pk_s_slow_roll =
  pvecback2[pba->index_bg_H]
   *pvecback2[pba->index_bg_H]
   /pvecback2[pba->index_bg_epsilon_H]
   /(8.*_PI_*_PI_);

double ns;
double epsilon_1;
double delta_H;
double epsilon_2; 
double As;
 epsilon_1 = pvecback2[pba->index_bg_epsilon_H];
// pba->epsilon_2 = -6.-2.*pvecback2[pba->index_bg_mu]*pvecback2[pba->index_bg_t]+2.*pvecback2[pba->index_bg_epsilon_H];
  pba->epsilon_2 = -6. -2.*pvecback[pba->index_bg_V_phi]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_phi_prime_scf]
  + 2.*pvecback[pba->index_bg_epsilon_H];

epsilon_2 = pvecback2[pba->index_bg_epsilon_2];
delta_H = epsilon_1 - epsilon_2/2.;
ns = 1. - 4.*epsilon_1 + 2.*delta_H;
As =  pk_s_slow_roll;
 printf("As = %e; ns =%e; epsilon_1 =%e; delta_H=%e; e2 = %e \n", As, ns, epsilon_1, delta_H, epsilon_2);

 double k_3 = 1.e-1*pba->k_star; //comment for EL
 double pk_s_3 = As*pow((k_3/(6.*pba->k_star)),ns-1.);//comment for EL
 printf("pk_3=%e \n",pk_s_3);

  sprintf(file_name,"%s_%s",pba->root,"pk_slow_roll_s.dat");

  class_open(pkfile,file_name,"w",pba->error_message);

    fprintf(pkfile,"#Slow roll primordial power spectrum P(k) %s\n","scalars");
    fprintf(pkfile,"# for k=%g to %g,\n",
             k_min_for_pk_s_slow_roll,
             k_max_for_pk_s_slow_roll);
    fprintf(pkfile,"#");
    class_fprintf_columntitle(pkfile,"k",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"P",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"fNL",_TRUE_,colnum);
   fprintf(pkfile," ");
    fprintf(pkfile,"\n");



k = k_min_for_pk_s_slow_roll;
// k= pba->k_star;//V. Sreenath
while (k < k_max_for_pk_s_slow_roll) {
      one_k = k;
      class_fprintf_double(pkfile,one_k,_TRUE_);
	double k_1 = k;
	double k_2 = k;
//	double k_3 = k;// uncomment for EL


	 double aH_sample=_HUGE_;
	   tau_max=pba->tau_end_inflation;
	   tau_min=pba->tau_start_inflation;

	   while (fabs(aH_sample-k)>
	         error_aH*k){

	    tau_half = (tau_max+tau_min)/2;

		  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
	               pba->error_message,
        	       pba->error_message);

		    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
		  if (aH_sample>k) tau_max = tau_half;
		  else tau_min = tau_half;
			 }

	  pk_s_slow_roll =  pvecback2[pba->index_bg_H]*pvecback2[pba->index_bg_H]
			    /pvecback2[pba->index_bg_epsilon_H]/(8.*_PI_*_PI_);
//	 pk_s_slow_roll = As*pow((k/(6.*pba->k_star)),ns-1.);
         one_pk = pk_s_slow_roll;
         class_fprintf_double(pkfile,one_pk,_TRUE_);

        
//        double V_phi = pba->m_scf_lqc*pba->m_scf_lqc*pvecback2[pba->index_bg_phi_scf]; 
	double V_phi = pvecback2[pba->index_bg_V_phi]; 
	epsilon_1 = pvecback2[pba->index_bg_epsilon_H];
         pba->epsilon_2 = -6.-2.*pvecback2[pba->index_bg_mu]*pvecback2[pba->index_bg_t]+2.*pvecback2[pba->index_bg_epsilon_H];
//	pba->epsilon_2 = -6.-2.*V_phi/pvecback2[pba->index_bg_H]/pvecback2[pba->index_bg_phi_prime_scf]
//	+pvecback2[pba->index_bg_phi_prime_scf]*pvecback2[pba->index_bg_phi_prime_scf]/pvecback2[pba->index_bg_H]/pvecback2[pba->index_bg_H];
        epsilon_2 = pba->epsilon_2;
        delta_H = epsilon_1 - epsilon_2/2.;


      double pk_s_1=pk_s_slow_roll;
	double pk_s_2=pk_s_slow_roll;
//	double pk_s_3=pk_s_slow_roll;//uncomment for EL
	printf("pk=%e\n",pk_s_slow_roll);
	
	H = pvecback2[pba->index_bg_H];
	double H_1 = H;//pvecback2[pba->index_bg_H];	
	double H_2 = H;//pvecback2[pba->index_bg_H];
	double H_3 = H;//pvecback2[pba->index_bg_H];
	
	double epsilon_H_1 = epsilon_1;
	double epsilon_H_2 = epsilon_1;
	double epsilon_H_3 = epsilon_1;
	
	double eta_H_1 = delta_H + epsilon_1;
	double eta_H_2 = delta_H + epsilon_1;
	double eta_H_3 = delta_H + epsilon_1;
	
	double g_123 = 1./(k_1+k_2+k_3);
	double kt = k_1 + k_2 + k_3 ;
  

    double G_RRR =
    pow(2.*_PI_*_PI_,2.)
    *pow(k_1,-3.)*pk_s_1
    *pow(k_2,-3.)*pk_s_2
    *(0.5*(3.*epsilon_H_3
           -2.*eta_H_3
           +epsilon_H_3*(pow(k_1,2.)+pow(k_2,2.))/pow(k_3,2.))
      +4.*epsilon_H_3
      *pow(k_1*k_2,2.)
      /pow(k_3,3.)
      *g_123)
+  pow(2.*_PI_*_PI_,2.)
    *pow(k_1,-3.)*pk_s_1
    *pow(k_3,-3.)*pk_s_3
    *(0.5*(3.*epsilon_H_2
           -2.*eta_H_2
           +epsilon_H_2*(pow(k_1,2.)+pow(k_3,2.))/pow(k_2,2.))
      +4.*epsilon_H_2
      *pow(k_1*k_3,2.)
      /pow(k_2,3.)
      *g_123)
    +pow(2.*_PI_*_PI_,2.)
    *pow(k_3,-3.)*pk_s_3
    *pow(k_2,-3.)*pk_s_2
    *(0.5*(3.*epsilon_H_1
           -2.*eta_H_1
           +epsilon_H_1*(pow(k_3,2.)+pow(k_2,2.))/pow(k_1,2.))
      +4.*epsilon_H_1
      *pow(k_3*k_2,2.)
      /pow(k_1,3.)
      *g_123);
  double fNL =//G_RRR;
  -10./3.
  *pow(2.*_PI_,-4.)
  *pow(k_1*k_2*k_3,3.)
  *G_RRR
  *pow(pow(k_1,3.)*pk_s_2*pk_s_3
       +pow(k_2,3.)*pk_s_1*pk_s_3
       +pow(k_3,3.)*pk_s_2*pk_s_1,-1.);

  class_fprintf_double(pkfile,fNL,_TRUE_);
//  fprintf(pkfile,"\n");



  double G1 = pow(H_1,4.)/(16.*epsilon_1)*pow(k_1*k_2*k_3,-3.)*( (1+k_1/kt)*pow(k_2*k_3,2.)/kt +
              (1+k_2/kt)*pow(k_3*k_1,2.)/kt +(1+k_3/kt)*pow(k_1*k_2,2.)/kt );
  double G2 = pow(H_1,4.)/(16.*epsilon_1)*pow(k_1*k_2*k_3,-3.)*(-k_1*k_1 -k_2*k_2 -k_3*k_3)/2.
		*(-kt + (k_1*k_2 + k_2*k_3 + k_3*k_1)/kt + k_1*k_2*k_3/pow(kt,2.));
  double G3 = -pow(H_1,4.)/(16.*epsilon_1)*pow(k_1*k_2*k_3,-3.)*( 
		(k_3*k_3-k_1*k_1-k_2*k_2)*k_3*k_3/(2.*kt)*(2. + (k_1+k_2)/kt) 
		+ (k_1*k_1-k_2*k_2-k_3*k_3)*k_1*k_1/(2.*kt)*(2. + (k_2+k_3)/kt)
		+ (k_2*k_2-k_3*k_3-k_1*k_1)*k_2*k_2/(2.*kt)*(2. + (k_3+k_1)/kt) ); 
  double Gfr = 2*pow(_PI_,4.)*epsilon_2*pow(k_1*k_2*k_3,-3.)*( pow(k_1,3.)*pk_s_2*pk_s_3 
		+ pow(k_2,3.)*pk_s_3*pk_s_1 + pow(k_3,3.)*pk_s_1*pk_s_2 );
  G_RRR = G1 + G2 + G3 + Gfr;
  fNL = -10./3.
  *pow(2.*_PI_,-4.)
  *pow(k_1*k_2*k_3,3.)
  *G_RRR
  *pow(pow(k_1,3.)*pk_s_2*pk_s_3
       +pow(k_2,3.)*pk_s_1*pk_s_3
       +pow(k_3,3.)*pk_s_2*pk_s_1,-1.);
	
  class_fprintf_double(pkfile,fNL,_TRUE_);  
  fprintf(pkfile,"\n");
   k *= pow(10.,1./(k_per_decade_for_pk_slow_roll));
}

/*      //////////(SLOW-ROLL) Scalar Power Spectrum //V. Sreenath commented
//  sprintf(file_name,"%s_%s",pba->root,"pk_slow_roll_s.dat");

 // class_open(pkfile,file_name,"w",pba->error_message);

//    fprintf(pkfile,"#Slow roll primordial power spectrum P(k) %s\n","scalars");
 //   fprintf(pkfile,"# for k=%g to %g,\n",
   //          k_min_for_pk_s_slow_roll,
     //        k_max_for_pk_s_slow_roll);
//    fprintf(pkfile,"#");
//    class_fprintf_columntitle(pkfile,"k",_TRUE_,colnum);
//    class_fprintf_columntitle(pkfile,"P",_TRUE_,colnum);
//    class_fprintf_columntitle(pkfile,"fNL",_TRUE_,colnum);
//   fprintf(pkfile," ");
//    fprintf(pkfile,"\n");
//
//
//
//k = k_min_for_pk_s_slow_roll;
// k= pba->k_star;//V. Sreenath
//while (k < k_max_for_pk_s_slow_roll) {
//      one_k = k;
//      class_fprintf_double(pkfile,one_k,_TRUE_);

//
// double aH_sample=_HUGE_;
//  tau_max=pba->tau_end_inflation;
//   tau_min=pba->tau_start_inflation;
//
//   while (fabs(aH_sample-k)>
//         error_aH*k){
//
//    tau_half = (tau_max+tau_min)/2;
//    
//  class_call(background_lqc_at_tau(pba,
//                                 tau_half,
//                                 pba->normal_info,
//                                 pba->inter_normal,
//                                 &index_back,
//                                 pvecback2),
//               pba->error_message,
//               pba->error_message);
//  
//    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
//  if (aH_sample>k) tau_max = tau_half;
//  else tau_min = tau_half;
// }
//
      
// pk_s_slow_roll = 
//  pvecback2[pba->index_bg_H]
//   *pvecback2[pba->index_bg_H]
//   /pvecback2[pba->index_bg_epsilon_H]
//   /(8.*_PI_*_PI_);
//
// one_pk =  pk_s_slow_roll;
// 
//  class_fprintf_double(pkfile,one_pk,_TRUE_);
//
//
//
//  double k_1 = k;
//  double pk_s_1 =pk_s_slow_roll;
//  double epsilon_H_1 = pvecback2[pba->index_bg_epsilon_H];
//  double eta_H_1 = epsilon_H_1;
//
//
//pba->epsilon_2 = 
  -6.
  -2.*pvecback[pba->index_bg_mu]*pvecback[pba->index_bg_t]
  +2.*pvecback[pba->index_bg_epsilon_H];
 


  
//   k = pba->lambda_k_2*k;   //V. Sreenath
   //k = 3.*pba->k_star;  

   aH_sample=_HUGE_;
   tau_max=pba->tau_end_inflation;
   tau_min=pba->tau_start_inflation;

   while (fabs(aH_sample-k)>
         error_aH*k){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  if (aH_sample>k) tau_max = tau_half;
  else tau_min = tau_half;
 }

      
 pk_s_slow_roll = 
  pvecback2[pba->index_bg_H]
   *pvecback2[pba->index_bg_H]
   /pvecback2[pba->index_bg_epsilon_H]
   /(8.*_PI_*_PI_);



  double k_2 = k;
  double pk_s_2 =pk_s_slow_roll;
  double epsilon_H_2 =  epsilon_H_1;//pvecback2[pba->index_bg_epsilon_H];
  double eta_H_2 = epsilon_H_2; 

// k = pba->lambda_k_3*k/pba->lambda_k_2; //V. Sreenath

// k = 1.e-6*pba->k_star;
   aH_sample=_HUGE_;
   tau_max=pba->tau_end_inflation;
   tau_min=pba->tau_start_inflation;


  
  class_call(background_lqc_at_tau(pba,
                                 tau_min,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  if (pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H]>k) 
    printf("k=%e,aH=%e\n",k,pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H]);

  else {
   while (fabs(aH_sample-k)>
         error_aH*k){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  if (aH_sample>k) tau_max = tau_half;
  else tau_min = tau_half;
 }

      
 pk_s_slow_roll = 
  pvecback2[pba->index_bg_H]
   *pvecback2[pba->index_bg_H]
   /pvecback2[pba->index_bg_epsilon_H]
   /(8.*_PI_*_PI_);
 }


//  double k_3 = k;//V. Sreenath commented this
// double k_3 = 1e-6*pbs->k_star;//V. Sreenath
  double pk_s_3 =pk_s_slow_roll;
  double epsilon_H_3 =  epsilon_H_1;//pvecback2[pba->index_bg_epsilon_H];
  double eta_H_3 = epsilon_H_3;
    
  double g_123 = 1./(k_1+k_2+k_3);
  double G_RRR = 
    pow(2.*_PI_*_PI_,2.)
    *pow(k_1,-3.)*pk_s_1
    *pow(k_2,-3.)*pk_s_2
    *(0.5*(3.*epsilon_H_3
           -2.*eta_H_3
           +epsilon_H_3*(pow(k_1,2.)+pow(k_2,2.))/pow(k_3,2.))
      +4.*epsilon_H_3
      *pow(k_1*k_2,2.)
      /pow(k_3,3.)
      *g_123)
+  pow(2.*_PI_*_PI_,2.)
    *pow(k_1,-3.)*pk_s_1
    *pow(k_3,-3.)*pk_s_3
    *(0.5*(3.*epsilon_H_2
           -2.*eta_H_2
           +epsilon_H_2*(pow(k_1,2.)+pow(k_3,2.))/pow(k_2,2.))
      +4.*epsilon_H_2
      *pow(k_1*k_3,2.)
      /pow(k_2,3.)
      *g_123)
    +pow(2.*_PI_*_PI_,2.)
    *pow(k_3,-3.)*pk_s_3
    *pow(k_2,-3.)*pk_s_2
    *(0.5*(3.*epsilon_H_1
           -2.*eta_H_1
           +epsilon_H_1*(pow(k_3,2.)+pow(k_2,2.))/pow(k_1,2.))
      +4.*epsilon_H_1
      *pow(k_3*k_2,2.)
      /pow(k_1,3.)
      *g_123);

  
  double fNL =//G_RRR;
  10./3.
  *pow(2.*_PI_,-4.)
  *pow(k_1*k_2*k_3,3.)
  *G_RRR
  *pow(pow(k_1,3.)*pk_s_2*pk_s_3
       +pow(k_2,3.)*pk_s_1*pk_s_3
       +pow(k_3,3.)*pk_s_2*pk_s_1,-1.);
   

  class_fprintf_double(pkfile,fNL,_TRUE_);  
  fprintf(pkfile,"\n");


  k = k/pba->lambda_k_3;
  //k = k_1;
        k *= pow(10.,1./(k_per_decade_for_pk_slow_roll));

 }*/

    fprintf(pkfile,"\n");



    //////////(SLOW-ROLL) Tensor Power Spectrum
  sprintf(file_name,"%s_%s",pba->root,"pk_slow_roll_t.dat");
  class_open(pkfile,file_name,"w",pba->error_message);
  colnum = 1;

    fprintf(pkfile,"#Slow roll primordial power spectrum P(k) %s\n","tensors");
    fprintf(pkfile,"# for k=%g to %g,\n",
             k_min_for_pk_s_slow_roll,
             k_max_for_pk_s_slow_roll);
    fprintf(pkfile,"#");
    class_fprintf_columntitle(pkfile,"k",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"P",_TRUE_,colnum);
   fprintf(pkfile," ");



k = k_min_for_pk_s_slow_roll;
while (k < k_max_for_pk_s_slow_roll) {
      k *= pow(10.,1./(k_per_decade_for_pk_slow_roll));
      one_k = k;
      class_fprintf_double(pkfile,one_k,_TRUE_);


 double aH_sample=_HUGE_;

tau_max=pba->tau_end_inflation;
 //tau_max=pba->tau_start_inflation;
   //+pba->N_end
   //-pba->N_start_inflation;
  tau_min=pba->tau_start_inflation;

 
  while (fabs(aH_sample-k)>
         error_aH*k){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
  if (aH_sample>k) tau_max = tau_half;
  else tau_min = tau_half;
 }

      
 pk_s_slow_roll = 
  pvecback2[pba->index_bg_H]
   *pvecback2[pba->index_bg_H]
   *16.
   /(8.*_PI_*_PI_);

 one_pk =  pk_s_slow_roll;
 
  class_fprintf_double(pkfile,one_pk,_TRUE_);
  fprintf(pkfile,"\n");
 }
    fprintf(pkfile,"\n");

    ///Derived Parameters
    
  sprintf(file_name,"%s_%s",pba->root,"lqc_derived_parameters.dat");
  class_open(pkfile,file_name,"w",pba->error_message);
  colnum = 1;
  
    fprintf(pkfile,"#Derived Parameters");
  fprintf(pkfile,"\n");
    fprintf(pkfile,"#");
    
    class_fprintf_columntitle(pkfile,"k_star",_TRUE_,colnum);
    //class_fprintf_columntitle(pkfile,"N_star",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"k_min",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"k_max",_TRUE_,colnum);
    
    class_fprintf_columntitle(pkfile,"t_start_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"t_end_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"t_start_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"t_end_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"t_start_inflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"t_end_inflation",_TRUE_,colnum);
    
    class_fprintf_columntitle(pkfile,"a_start_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"a_end_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"a_start_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"a_end_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"a_start_inflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"a_end_inflation",_TRUE_,colnum);

     
    class_fprintf_columntitle(pkfile,"N_start_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"N_end_remote_past",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"N_start_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"N_end_deflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"N_start_inflation",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"N_end_inflation",_TRUE_,colnum);
    
    class_fprintf_columntitle(pkfile,"N_exit",_TRUE_,colnum);
   
  fprintf(pkfile," ");
  fprintf(pkfile,"\n");
  class_fprintf_double(pkfile,pba->k_star,_TRUE_);
  fprintf(pkfile," ");
  //class_fprintf_double(pkfile,pba->N_star,_TRUE_);
  //fprintf(pkfile,"\n");
  class_fprintf_double(pkfile,pba->alpha_min_for_pk_lqc*pba->k_star,_TRUE_);
  //fprintf(pkfile," ");
  class_fprintf_double(pkfile,pba->alpha_max_for_pk_lqc*pba->k_star,_TRUE_);
  //fprintf(pkfile," ");
  
  class_fprintf_double(pkfile,pba->tau_start_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->tau_end_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->tau_start_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->tau_end_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->tau_start_inflation,_TRUE_);
  class_fprintf_double(pkfile,pba->tau_end_inflation,_TRUE_);

  class_fprintf_double(pkfile,pba->a_start_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->a_end_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->a_start_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->a_end_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->a_start_inflation,_TRUE_);
  class_fprintf_double(pkfile,pba->a_end_inflation,_TRUE_);

  class_fprintf_double(pkfile,pba->N_start_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->N_end_remote_past,_TRUE_);
  class_fprintf_double(pkfile,pba->N_start_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->N_end_deflation,_TRUE_);
  class_fprintf_double(pkfile,pba->N_start_inflation,_TRUE_);
  class_fprintf_double(pkfile,pba->N_end,_TRUE_);
  
  class_fprintf_double(pkfile,pba->N_exit,_TRUE_);




  fprintf(pkfile,"\n");    

    

    
    fclose(pkfile);

 
    free(pvecback2);                        
  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background_lqc integrated variables.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background_lqc structure
 * @param pvecback             Input: vector of background_lqc quantities used as workspace
 * @param pvecback_integration Output: vector of background_lqc quantities to be integrated, returned with proper initial values
 * @return the error status
 */

int background_lqc_initial_conditions(
                                  struct precision *ppr,
                                  struct background_lqc *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {


     pvecback_integration[pba->index_bi_a] = pba->a_bounce;
     pvecback_integration[pba->index_bi_conformal_time] = 0.;
     pvecback_integration[pba->index_bi_cosmic_time] = 0.;
     pvecback_integration[pba->index_bi_N] = 0.;
     
     pvecback_integration[pba->index_bi_phi] = pba->phi_ini;
//220517     pvecback_integration[pba->index_bi_phi_dot] = pba->phi_dot_ini;
     pvecback_integration[pba->index_bi_Pphi] = sqrt(2.*pba->rho_bounce - pow(pba->m_scf_lqc*pba->phi_ini,2.))*(4./pow(8.*_PI_,3./2.))/4.;
     pvecback_integration[pba->index_bi_v] = 4./pow(8.*_PI_,3./2.);

//     pvecback_integration[pba->index_bi_h] = -_PI_/(2.*pba->l0_lqc);
//   printf("inside initial conditions.... branch = %d",pba->future_branch); 
     if(pba->future_branch == 0)
      pvecback_integration[pba->index_bi_h] = -_PI_/(2.*pba->l0_lqc);
     else 
      pvecback_integration[pba->index_bi_h] = -_PI_/(2.*pba->l0_lqc);

    
//     pvecback_integration[pba->index_bi_x] = pba->x_ini;

//     pvecback_integration[pba->index_bi_y] = sqrt(pba->rho_ini/pba->rho_bounce-pba->x_ini*pba->x_ini);
      
             
      /*
        printf("H init=%e\n",pvecback_integration[pba->index_bi_x]
               *pvecback_integration[pba->index_bi_x]
               +pvecback_integration[pba->index_bi_y]
               *pvecback_integration[pba->index_bi_y]);    

               */
 
  class_call(background_lqc_functions(pba, pvecback_integration, pba->normal_info, pvecback),
	     pba->error_message,
	     pba->error_message);
  
   
            if (pba->has_lqc == _TRUE_)
              pvecback_integration[pba->index_bi_H] = 0.;

            
    else pvecback_integration[pba->index_bi_H] =pvecback[pba->index_bg_H];

      
    pvecback_integration[pba->index_bi_tau] = 0.;

    
  return _SUCCESS_;

}

/**
 * Subroutine for formatting background_lqc output
 * 
 */

int background_lqc_output_titles(struct background_lqc * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){

  /** - Length of the column title should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[50];

  //class_store_columntitle(titles,"z",_TRUE_);
  //class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  //class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"time",_TRUE_);
  class_store_columntitle(titles,"scale factor",_TRUE_);
  class_store_columntitle(titles,"x",_TRUE_);
  class_store_columntitle(titles,"y",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"epsilon_H",_TRUE_);
  class_store_columntitle(titles,"epsilon_2",_TRUE_);
//  class_store_columntitle(titles,"mu",_TRUE_);
//  class_store_columntitle(titles,"t*epsilon_H",_TRUE_);
  class_store_columntitle(titles,"(.)rho_scf",_TRUE_);
  class_store_columntitle(titles,"(.)p_scf",_TRUE_);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
//  class_store_columntitle(titles,"f_S",pba->has_scf);
//  class_store_columntitle(titles,"f_T",pba->has_scf);
  class_store_columntitle(titles,"H2_f_S_Q",pba->has_scf);
  class_store_columntitle(titles,"z_S",pba->has_scf);
  class_store_columntitle(titles,"z_T",pba->has_scf);
  class_store_columntitle(titles,"N",pba->has_scf);
  class_store_columntitle(titles,"conformal time",pba->has_scf);
  class_store_columntitle(titles,"cosmic time",pba->has_scf);
//  class_store_columntitle(titles,"V_s",pba->has_scf);
  class_store_columntitle(titles,"zs_pp_over_zs",pba->has_scf);
  class_store_columntitle(titles,"zt_pp_over_zt",pba->has_scf);
  class_store_columntitle(titles,"v",_TRUE_);
  class_store_columntitle(titles,"h",_TRUE_);
  class_store_columntitle(titles,"phi",_TRUE_);
  class_store_columntitle(titles,"Pphi",_TRUE_);

  class_store_columntitle(titles,"pia=6a^2h",_TRUE_);
  class_store_columntitle(titles,"pia_c = -6a^2H",_TRUE_);

  class_store_columntitle(titles,"H2_f_S_Q_phys",pba->has_scf);
  class_store_columntitle(titles,"V(phi_scf)",pba->has_scf);
  class_store_columntitle(titles,"cosmic_time(Tpl)",pba->has_scf);
  class_store_columntitle(titles,"Pphidot",pba->has_scf);
  class_store_columntitle(titles,"piadot",_TRUE_);
  class_store_columntitle(titles,"Up",_TRUE_);
  class_store_columntitle(titles,"Uf",_TRUE_);
  class_store_columntitle(titles,"piap",_TRUE_);
  class_store_columntitle(titles,"piaf",_TRUE_);
  class_store_columntitle(titles,"piadotp",_TRUE_);
  class_store_columntitle(titles,"piadotf",_TRUE_);

  return _SUCCESS_;
}

int background_lqc_output_data(
                           struct background_lqc *pba,
                           int number_of_titles,
                           double *data){
  int index_tau, storeidx, n;
  double *dataptr, *pvecback, *pvecback_pb, *pvecback_fb;

  /** Stores quantities */
  for (index_tau=0; index_tau<pba->bt_size; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_lqc_table + index_tau*pba->bg_size;
    pvecback_pb = pba->background_lqc_table_pb + index_tau*pba->bg_size;
    pvecback_fb = pba->background_lqc_table_fb + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,pba->tau_table[index_tau],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_a],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_x],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_y],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_epsilon_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_epsilon_2],_TRUE_,storeidx);
//    class_store_double(dataptr,pvecback[pba->index_bg_mu],_TRUE_,storeidx);
//    class_store_double(dataptr,pvecback[pba->index_bg_t_times_epsilon_H],_TRUE_,storeidx);   
    class_store_double(dataptr,pvecback[pba->index_bg_rho_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_scf],pba->has_scf,storeidx);
//    class_store_double(dataptr,pvecback[pba->index_bg_f_S],pba->has_scf,storeidx);
//    class_store_double(dataptr,pvecback[pba->index_bg_f_T],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H2_f_S_Q],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_z_S],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_z_T],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_N],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_conformal_time],pba->has_scf,storeidx);  
    class_store_double(dataptr,pvecback[pba->index_bg_cosmic_time],pba->has_scf,storeidx);  
//            class_store_double(dataptr,pvecback[pba->index_bg_effective_potential],pba->has_scf,storeidx);  
    class_store_double(dataptr,pvecback[pba->index_bg_zS_primeprime_over_zS],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_zT_primeprime_over_zT],pba->has_scf,storeidx);  
    class_store_double(dataptr,pvecback[pba->index_bg_v],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_h],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_Pphi_scf],_TRUE_,storeidx);
    
    class_store_double(dataptr,pvecback[pba->index_bg_pia],_TRUE_,storeidx);
    class_store_double(dataptr,-6.*pow(pvecback[pba->index_bg_a],2.)*pvecback[pba->index_bg_H],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_H2_f_S_Q_phys],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_V],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_cosmic_time],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_Pphidot_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_piadot],pba->has_scf,storeidx);

    class_store_double(dataptr,pvecback_pb[pba->index_bg_H2_f_S_Q],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback_fb[pba->index_bg_H2_f_S_Q],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback_pb[pba->index_bg_pia],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback_fb[pba->index_bg_pia],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback_pb[pba->index_bg_piadot],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback_fb[pba->index_bg_piadot],pba->has_scf,storeidx);


  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background_lqc structure
 * and to a background_lqc vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_lqc_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_lqc_parameters_and_workspace * pbpaw;
  struct background_lqc * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_lqc_functions() */
  class_call(background_lqc_functions(pba, y, pba->normal_info, pvecback),
             pba->error_message,
             error_message);

 
   dy[pba->index_bi_H] = 0.; 
   dy[pba->index_bi_a] = 0.;
//   dy[pba->index_bi_phi_dot] = 0.;


   dy[pba->index_bi_conformal_time] = 1./(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);
   dy[pba->index_bi_N] = 1.; 
   dy[pba->index_bi_cosmic_time] = 1./pvecback[pba->index_bg_H];
   dy[pba->index_bi_v] = (-3.* y[pba->index_bi_v]/(2.*pba->l0_lqc)*sin(2.*pba->l0_lqc*y[pba->index_bi_h]))/pvecback[pba->index_bg_H];
   dy[pba->index_bi_h] = (4.*pow(y[pba->index_bi_Pphi]/y[pba->index_bi_v],2.) 
	- pvecback[pba->index_bg_V]/2. + 3.*pow(sin(pba->l0_lqc*y[pba->index_bi_h]),2.)/(2.*pow(pba->l0_lqc,2.)))
	/pvecback[pba->index_bg_H];
   dy[pba->index_bi_phi] = (y[pba->index_bi_Pphi]*4./y[pba->index_bi_v])/pvecback[pba->index_bg_H];
   dy[pba->index_bi_Pphi] = (-y[pba->index_bi_v]*pvecback[pba->index_bg_V_phi]/4.)/pvecback[pba->index_bg_H];

/*   printf("v=%e\n", y[pba->index_bi_v]);
   printf("h=%e\n", y[pba->index_bi_h]);
   printf("phi=%e\n", y[pba->index_bi_phi]);
   printf("Pphi=%e\n", y[pba->index_bi_Pphi]);
   printf("H=%e\n", pvecback[pba->index_bg_H]);
   printf("a=%e\n", pvecback[pba->index_bg_a]);*/

  if (pba->post_bounce == _FALSE_){

   dy[pba->index_bi_conformal_time] = 1./(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);
   dy[pba->index_bi_N] = -1.; 
   dy[pba->index_bi_cosmic_time] = 1./(pvecback[pba->index_bg_H]);

   dy[pba->index_bi_v] = (-3.* y[pba->index_bi_v]/(2.*pba->l0_lqc)*sin(2.*pba->l0_lqc*y[pba->index_bi_h]))/pvecback[pba->index_bg_H];
   dy[pba->index_bi_h] = (4.*pow(y[pba->index_bi_Pphi]/y[pba->index_bi_v],2.) 
	- pvecback[pba->index_bg_V]/2. + 3.*pow(sin(pba->l0_lqc*y[pba->index_bi_h]),2.)/(2.*pow(pba->l0_lqc,2.)))
	/pvecback[pba->index_bg_H];
   dy[pba->index_bi_phi] = (y[pba->index_bi_Pphi]*4./y[pba->index_bi_v])/pvecback[pba->index_bg_H];
   dy[pba->index_bi_Pphi] = (-y[pba->index_bi_v]*pvecback[pba->index_bg_V_phi]/4.)/pvecback[pba->index_bg_H];

/*  dy[pba->index_bi_x] = (pba->m_scf_lqc/pvecback[pba->index_bg_H])
                         *y[pba->index_bi_y];
    dy[pba->index_bi_y] = - 3.*y[pba->index_bi_y]
                          -(pba->m_scf_lqc/pvecback[pba->index_bg_H])
                          *y[pba->index_bi_x];*/
//220517    dy[pba->index_bi_phi] = y[pba->index_bi_phi_dot]/pvecback[pba->index_bg_H];
//220517   dy[pba->index_bi_phi_dot] = -3.*y[pba->index_bi_phi_dot] - pvecback[pba->index_bg_V_phi]/pvecback[pba->index_bg_H];
/*   printf("t=%e\n", y[pba->index_bi_cosmic_time]);
   printf("v=%e\n", y[pba->index_bi_v]);
   printf("h=%e\n", y[pba->index_bi_h]);
   printf("phi=%e\n", y[pba->index_bi_phi]);
   printf("Pphi=%e\n", y[pba->index_bi_Pphi]);
   printf("H=%e\n", pvecback[pba->index_bg_H]);
   printf("a=%e\n", pvecback[pba->index_bg_a]);
*/
    }


   if (pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2)>0.999) {

    
   dy[pba->index_bi_cosmic_time] = 1.;
   
//220517   if (pba->post_bounce == _FALSE_)dy[pba->index_bi_cosmic_time] = -dy[pba->index_bi_cosmic_time];

   dy[pba->index_bi_conformal_time] = 1./pvecback[pba->index_bg_a];             
   dy[pba->index_bi_N] = pvecback[pba->index_bg_H]; 

/*   dy[pba->index_bi_H] =
     -pba->rho_bounce
     *pow(y[pba->index_bi_y],2)
     *(1.-2.*(pow(y[pba->index_bi_x],2)+pow(y[pba->index_bi_y],2))); */
//220517   dy[pba->index_bi_H] = -pow(y[pba->index_bi_phi_dot],2)*(1.-2.*(pow(pvecback[pba->index_bg_x],2)+pow(pvecback[pba->index_bg_y],2))); 
   
/* if (pba->has_lqc == _FALSE_)
 dy[pba->index_bi_H] =
   -pba->rho_bounce
   *pow(y[pba->index_bi_y],2);*/
//220517   if (pba->has_lqc == _FALSE_)dy[pba->index_bi_H] = -pow(y[pba->index_bi_phi_dot],2);

//220517   dy[pba->index_bi_a] = pow(y[pba->index_bi_a],1)*y[pba->index_bi_H];


              
/*   dy[pba->index_bi_x] =
     pba->m_scf_lqc
     *y[pba->index_bi_y];
   
   dy[pba->index_bi_y] =
     -3.
     *y[pba->index_bi_y]
     *y[pba->index_bi_H]
     -pba->m_scf_lqc
     *y[pba->index_bi_x];*/
//220517   dy[pba->index_bi_phi] = y[pba->index_bi_phi_dot];
//220517   dy[pba->index_bi_phi_dot] = -3.*y[pba->index_bi_H]*y[pba->index_bi_phi_dot] - pvecback[pba->index_bg_V_phi];

   dy[pba->index_bi_v] = -3.* y[pba->index_bi_v]/(2.*pba->l0_lqc)*sin(2.*pba->l0_lqc*y[pba->index_bi_h]);
   dy[pba->index_bi_h] = 4.*pow(y[pba->index_bi_Pphi]/y[pba->index_bi_v],2.) 
        - pvecback[pba->index_bg_V]/2. + 3.*pow(sin(pba->l0_lqc*y[pba->index_bi_h]),2.)/(2.*pow(pba->l0_lqc,2.));
   dy[pba->index_bi_phi] = y[pba->index_bi_Pphi]*4./y[pba->index_bi_v];
   dy[pba->index_bi_Pphi] = -y[pba->index_bi_v]*pvecback[pba->index_bg_V_phi]/4.;

/*   printf("v=%e\n", y[pba->index_bi_v]);
   printf("h=%e\n", y[pba->index_bi_h]);
   printf("phi=%e\n", y[pba->index_bi_phi]);
   printf("Pphi=%e\n", y[pba->index_bi_Pphi]);
   printf("H=%e\n", pvecback[pba->index_bg_H]);*/

  if (pba->post_bounce == _FALSE_){
     
    //dy[pba->index_bi_a] =
    //            -pow(y[pba->index_bi_a],1.)
    //           *y[pba->index_bi_H];
   dy[pba->index_bi_cosmic_time] = -1.;
   dy[pba->index_bi_conformal_time] = -1./pvecback[pba->index_bg_a];// removed an overall - sign           
   dy[pba->index_bi_N] = pvecback[pba->index_bg_H]; //removed an overall - sign
       //if (pba->has_lqc == _FALSE_)
       //dy[pba->index_bi_N] =y[pba->index_bi_H];

              
/*   dy[pba->index_bi_x] = -pba->m_scf_lqc*y[pba->index_bi_y];
   dy[pba->index_bi_y] = -3.*y[pba->index_bi_y]*y[pba->index_bi_H]
                          +pba->m_scf_lqc*y[pba->index_bi_x];*/

//220517   dy[pba->index_bi_phi] = -y[pba->index_bi_phi_dot];
//220517   dy[pba->index_bi_phi_dot] = -3.*y[pba->index_bi_H]*y[pba->index_bi_phi_dot] + pvecback[pba->index_bg_V_phi];
   

/*    dy[pba->index_bi_x] = -(pvecback[pba->index_bg_V_phi]*y[pba->index_bi_y]*y[pba->index_bi_H]
        /(sqrt(2.*pvecback[pba->index_bg_V])));

    dy[pba->index_bi_y] = (- 3.*y[pba->index_bi_y]*y[pba->index_bi_H] + pvecback[pba->index_bg_V_phi]
        /(sqrt(2.*pba->rho_bounce)));*/
   dy[pba->index_bi_v] = 3.* y[pba->index_bi_v]/(2.*pba->l0_lqc)*sin(2.*pba->l0_lqc*y[pba->index_bi_h]);
   dy[pba->index_bi_h] = -(4.*pow(y[pba->index_bi_Pphi]/y[pba->index_bi_v],2.) 
	- pvecback[pba->index_bg_V]/2. + 3.*pow(sin(pba->l0_lqc*y[pba->index_bi_h]),2.)/(2.*pow(pba->l0_lqc,2.)));
   dy[pba->index_bi_phi] = -(y[pba->index_bi_Pphi]*4./y[pba->index_bi_v]);
   dy[pba->index_bi_Pphi] = -(-y[pba->index_bi_v]*pvecback[pba->index_bg_V_phi]/4.);

/*   printf("v=%e\n", y[pba->index_bi_v]);
   printf("h=%e\n", y[pba->index_bi_h]);
   printf("phi=%e\n", y[pba->index_bi_phi]);
   printf("Pphi=%e\n", y[pba->index_bi_Pphi]);
   printf("H=%e\n", pvecback[pba->index_bg_H]);
*/
//220517    if (pba->has_lqc == _FALSE_){
//220517       dy[pba->index_bi_H] = - dy[pba->index_bi_H];
/*       dy[pba->index_bi_y] = 3.*y[pba->index_bi_y]*y[pba->index_bi_H]
                          +pba->m_scf_lqc*y[pba->index_bi_x];*/

//220517   dy[pba->index_bi_phi_dot] = 3.*y[pba->index_bi_H]*y[pba->index_bi_phi_dot] + pvecback[pba->index_bg_V_phi];

/*       dy[pba->index_bi_y] = ( 3.*y[pba->index_bi_y] + pvecback[pba->index_bg_V_phi]
        /(sqrt(2.*pba->rho_bounce)*y[pba->index_bi_H]))*y[pba->index_bi_H];*/

//220517  }
   
    }

   
 }


  ////////////////////////////////////
  if (pba->remote_past == _TRUE_){
/*    dy[pba->index_bi_x] = 0.;
    dy[pba->index_bi_y] = 0.;*/
//220517    dy[pba->index_bi_phi] = 0.;
//220517    dy[pba->index_bi_phi_dot] = 0.;
   dy[pba->index_bi_v] = 0.;
   dy[pba->index_bi_h] = 0.;
   dy[pba->index_bi_phi] = 0.;
   dy[pba->index_bi_Pphi] = 0.;

  }
/*   printf("t=%e\n", y[pba->index_bi_cosmic_time]);
   printf("v=%e\n", y[pba->index_bi_v]);
   printf("h=%e\n", y[pba->index_bi_h]);
   printf("phi=%e\n", y[pba->index_bi_phi]);
   printf("Pphi=%e\n", y[pba->index_bi_Pphi]);
   printf("H=%e\n", pvecback[pba->index_bg_H]);
   printf("a=%e\n", pvecback[pba->index_bg_a]);*/

//   y[pba->index_bi_x] = sqrt(pvecback[pba->index_bg_V]/pba->rho_bounce); 
//   y[pba->index_bi_y] = y[pba->index_bi_phi_dot]/sqrt(2.*pba->rho_bounce);
  return _SUCCESS_;

}
