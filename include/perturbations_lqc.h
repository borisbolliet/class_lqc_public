/** @file perturbations_lqc.h Documented includes for perturbation module */

#ifndef __PERTURBATIONS_LQC__
#define __PERTURBATIONS_LQC__

#include "background_lqc.h"
#include "evolver_ndf15.h"
#include "evolver_rkck.h"

#define _scalars_lqc_ ((ppt->has_scalars == _TRUE_) && (index_md == ppt->index_md_scalars))
#define _tensors_lqc_ ((ppt->has_tensors == _TRUE_) && (index_md == ppt->index_md_tensors))

#define _set_source_lqc_(index) ppt->sources[index_md][index][index_tau * ppt->k_size[index_md]+index_k]
#define _set_pk_lqc_(index) ppt->pk_lqc[index_md][index][index_k]

#define _set_fields_at_tau_end_lqc_(index) ppt->fields_at_tau_end_lqc[index_md][index][index_k]

/**
 * List of coded gauges. More gauges can in principle be defined.
 */

//@{



/**
 * maximum number of k-values for perturbation output
 */
#define _MAX_NUMBER_OF_K_FILES_ 30

//@}



/**
 * Structure containing everything about perturbations_lqc that other
 * modules need to know, in particular tabled values of the source
 * functions \f$ S(k, \tau) \f$ for all requested modes
 * (scalar/vector/tensor), initial conditions, types (temperature,
 * E-polarization, B-polarization, lensing potential, etc), multipole
 * l and wavenumber k.
 *
 */

enum time_flags {efolds, cosmic_time};



struct perturbs_lqc
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  int ln_k_size;    /**< number ln(k) values */
  double * ln_k;    /**< list of ln(k) values ln_k[index_k] */
  double * ln_pk;

  
  double * bispectrum_s_equi; 
  double * bispectrum_t_equi; 


  int index_tau_pk;
  
  FileName root; /**< root for all file names */

  
  short has_perturbations_lqc; /**< do we need to compute perturbations_lqc at all ? */

  int index_tp_gw;  /**< index value for delta of scalar field */
  int index_tp_phi_scf;  /**< index value for delta of scalar field */
  int index_tp_pk_s_1;  /**< index value for delta of scalar field */
  int index_tp_pk_s_2;  /**< index value for delta of scalar field */
  int index_tp_pk_s_3;  /**< index value for delta of scalar field */


  int index_tp_bispectrum_s_123_111;
  int index_tp_bispectrum_s_123_110;
  int index_tp_bispectrum_s_123_101;
  int index_tp_bispectrum_s_123_011;
  int index_tp_bispectrum_s_123_100;
  int index_tp_bispectrum_s_123_010;
  int index_tp_bispectrum_s_123_001;
  int index_tp_bispectrum_s_123_000;
  int index_tp_bispectrum_s_prefactor_123_000;

  int index_tp_bispectrum_s_1d23_111;
  int index_tp_bispectrum_s_1d23_110;
  int index_tp_bispectrum_s_1d23_101;
  int index_tp_bispectrum_s_1d23_011;
  int index_tp_bispectrum_s_1d23_100;
  int index_tp_bispectrum_s_1d23_010;
  int index_tp_bispectrum_s_1d23_001;
  int index_tp_bispectrum_s_1d23_000;

   int index_tp_bispectrum_s_12d3_111;
  int index_tp_bispectrum_s_12d3_110;
  int index_tp_bispectrum_s_12d3_101;
  int index_tp_bispectrum_s_12d3_011;
  int index_tp_bispectrum_s_12d3_100;
  int index_tp_bispectrum_s_12d3_010;
  int index_tp_bispectrum_s_12d3_001;
  int index_tp_bispectrum_s_12d3_000;

   
  int index_tp_bispectrum_s_123d_111;
  int index_tp_bispectrum_s_123d_110;
  int index_tp_bispectrum_s_123d_101;
  int index_tp_bispectrum_s_123d_011;
  int index_tp_bispectrum_s_123d_100;
  int index_tp_bispectrum_s_123d_010;
  int index_tp_bispectrum_s_123d_001;
  int index_tp_bispectrum_s_123d_000;



  int index_tp_bispectrum_s_12d3d_111;
  int index_tp_bispectrum_s_12d3d_110;
  int index_tp_bispectrum_s_12d3d_101;
  int index_tp_bispectrum_s_12d3d_011;
  int index_tp_bispectrum_s_12d3d_100;
  int index_tp_bispectrum_s_12d3d_010;
  int index_tp_bispectrum_s_12d3d_001;
  int index_tp_bispectrum_s_12d3d_000;

  
  int index_tp_bispectrum_s_1d2d3_111;
  int index_tp_bispectrum_s_1d2d3_110;
  int index_tp_bispectrum_s_1d2d3_101;
  int index_tp_bispectrum_s_1d2d3_011;
  int index_tp_bispectrum_s_1d2d3_100;
  int index_tp_bispectrum_s_1d2d3_010;
  int index_tp_bispectrum_s_1d2d3_001;
  int index_tp_bispectrum_s_1d2d3_000;

  int index_tp_bispectrum_s_1d23d_111;
  int index_tp_bispectrum_s_1d23d_110;
  int index_tp_bispectrum_s_1d23d_101;
  int index_tp_bispectrum_s_1d23d_011;
  int index_tp_bispectrum_s_1d23d_100;
  int index_tp_bispectrum_s_1d23d_010;
  int index_tp_bispectrum_s_1d23d_001;
  int index_tp_bispectrum_s_1d23d_000;


    
 
  
      

  
  
  int index_tp_bispectrum_t_equi_111;
  int index_tp_bispectrum_t_equi_110;
  int index_tp_bispectrum_t_equi_100;
  int index_tp_bispectrum_t_equi_000;



  int index_tp_bispectrum_s_tau_end_111;
  int index_tp_bispectrum_s_tau_end_110;
  int index_tp_bispectrum_s_tau_end_101;
  int index_tp_bispectrum_s_tau_end_011;
  int index_tp_bispectrum_s_tau_end_100;
  int index_tp_bispectrum_s_tau_end_010;
  int index_tp_bispectrum_s_tau_end_001;
  int index_tp_bispectrum_s_tau_end_000;


  int index_tp_bispectrum_s_tau_end_field_redef;


  
  int index_tp_bispectrum_t_equi_tau_end_111;
  int index_tp_bispectrum_t_equi_tau_end_110;
  int index_tp_bispectrum_t_equi_tau_end_100;
  int index_tp_bispectrum_t_equi_tau_end_000;
  
  int index_tp_bispectrum_s_equi_integrand;
  int index_tp_bispectrum_t_equi_integrand;
  
  short has_source_gw;
  short has_source_phi_scf;
  
  short has_source_bispectrum_s_equi; 
  short has_source_bispectrum_t_equi; 

  
  int * tp_size; /**< number of types tp_size[index_md] included in computation for each mode */

  double gw_ini;
  double tau_ini;
  double tau_ini_for_k;
  double tau_end;
  double tol;
  double integration_depth_ini;
  double integration_depth_end;

  double lambda_k_2;
  double lambda_k_3;
  
  double k_min_for_pk;
  double k_max_for_pk;
  double delta_cut_off_bispectrum;
  

  short has_scalars; /**< do we need scalars? */
  short has_tensors; /**< do we need tensors? */

  short has_ad;      /**< do we need adiabatic mode? */

  short second_loop; /**< do we need scalars? */

  short has_pk_lqc;                /**< do we need lqc Fourier spectrum? */
  short tensors_only;                /**< do we need lqc Fourier spectrum? */
  short has_bispectrum_lqc;                /**< do we need lqc Fourier spectrum? */

  int k_output_values_num;       /**< Number of perturbation outputs (default=0) */
  
  int index_k_output;
  double k_output_values[_MAX_NUMBER_OF_K_FILES_];    /**< List of k values where perturbation output is requested. */
  int *index_k_output_values; /**< List of indices corresponding to k-values close to k_output_values for each mode. [index_md*k_output_values_num+ik]*/
  char scalar_titles[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for scalar perturbation output files. */
  char tensor_titles[_MAXTITLESTRINGLENGTH_]; /**< _DELIMITER_ separated string of titles for tensor perturbation output files. */
  int number_of_scalar_titles; /**< number of titles/columns in scalar perturbation output files */
  int number_of_tensor_titles; /**< number of titles/columns in tensor perturbation output files*/

  int store_perturbations;  /**< Do we want to store perturbations? */

  double * scalar_perturbations_lqc_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for scalars */
  double * tensor_perturbations_lqc_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of double pointers to perturbation output for tensors */
 int size_scalar_perturbation_lqc_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of scalar double pointers  */
 int size_tensor_perturbation_lqc_data[_MAX_NUMBER_OF_K_FILES_]; /**< Array of sizes of tensor double pointers  */

  //@}

  /** @name - useful flags inferred from the ones above */

  //@{

  //@}

  //@{


  //@}

  /** @name - indices running on modes (scalar, vector, tensor) */

  //@{

  int index_md_scalars; /**< index value for scalars */
  int index_md_tensors; /**< index value for vectors */

  int md_size; /**< number of modes included in computation */

  //@}


  //@{


  //@}

  /** @name - list of k values for each mode */

  //@{
int * x_size;  

  int * k_size;     /**< k_size[index_md] = total number of k
                       values, including those needed for P(k) but not
                       for \f$ C_l \f$'s */

  double ** k;      /**< k[index_md][index_k] = list of values */



  double ** table_x2;      /**< k[index_md][index_k] = list of values */
  double ** table_x3;      /**< k[index_md][index_k] = list of values */
  double ** table_fNL;      /**< k[index_md][index_k] = list of values */

  int index_x_2_x_3;
  
  double k_min;     /**< minimum value (over all modes) */
  double k_max;     /**< maximum value (over all modes) */

  double k_per_decade_for_pk_lqc;
  //@}

  /** @name - list of conformal time values in the source table
      (common to all modes and types) */

  //@{

  int tau_size;          /**< tau_size = number of values */

  double * tau_sampling; /**< tau_sampling[index_tau] = list of tau values */

  //@}

  /** @name - source functions interpolation table */

  //@{

  double *** sources; /**< Pointer towards the source interpolation table
                         sources[index_md]
                         [index_ic * ppt->tp_size[index_md] + index_type]
                         [index_tau * ppt->k_size + index_k] */

 double *** pk_lqc; /**< Pointer towards the power spectrum
                         sources[index_md]
                         [index_type]
                         [index_k] */

 double *** fields_at_tau_end_lqc; /**< Pointer towards the power spectrum
                         sources[index_md]
                         [index_type]
                         [index_k] */


  //@}

  /** @name - technical parameters */

  //@{

  short perturbations_lqc_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 */

struct perturb_lqc_vector
{
  int index_pt_phi_scf_re;  
  int index_pt_phi_prime_scf_re;  

  int index_pt_phi_scf_re_1;  
  int index_pt_phi_prime_scf_re_1;  
  int index_pt_phi_scf_re_2;  
  int index_pt_phi_prime_scf_re_2;  
  int index_pt_phi_scf_re_3;  
  int index_pt_phi_prime_scf_re_3;  

  
  int index_pt_gw_re;        
  int index_pt_gwdot_re;

  int index_pt_phi_scf_im;  
  int index_pt_phi_prime_scf_im;  

  int index_pt_phi_scf_im_1;  
  int index_pt_phi_prime_scf_im_1;  

  int index_pt_phi_scf_im_2;  
  int index_pt_phi_prime_scf_im_2;  

  int index_pt_phi_scf_im_3;  
  int index_pt_phi_prime_scf_im_3;  

  int index_pt_gw_im;        
  int index_pt_gwdot_im;

  
  int index_pt_bispectrum_t_000;
  int index_pt_bispectrum_t_100;
  int index_pt_bispectrum_t_110;
  int index_pt_bispectrum_t_111;


  int index_pt_bispectrum_s_prefactor_123_000;
  int index_pt_bispectrum_s_123_000;
  int index_pt_bispectrum_s_123_100;
  int index_pt_bispectrum_s_123_010;
  int index_pt_bispectrum_s_123_001;
  int index_pt_bispectrum_s_123_110;
  int index_pt_bispectrum_s_123_101;
  int index_pt_bispectrum_s_123_011;
  int index_pt_bispectrum_s_123_111;

    
  int index_pt_bispectrum_s_123d_000;
  int index_pt_bispectrum_s_123d_100;
  int index_pt_bispectrum_s_123d_010;
  int index_pt_bispectrum_s_123d_001;
  int index_pt_bispectrum_s_123d_110;
  int index_pt_bispectrum_s_123d_101;
  int index_pt_bispectrum_s_123d_011;
  int index_pt_bispectrum_s_123d_111;

  int index_pt_bispectrum_s_12d3d_000;
  int index_pt_bispectrum_s_12d3d_100;
  int index_pt_bispectrum_s_12d3d_010;
  int index_pt_bispectrum_s_12d3d_001;
  int index_pt_bispectrum_s_12d3d_110;
  int index_pt_bispectrum_s_12d3d_101;
  int index_pt_bispectrum_s_12d3d_011;
  int index_pt_bispectrum_s_12d3d_111;

    
  int index_pt_bispectrum_s_12d3_000;
  int index_pt_bispectrum_s_12d3_100;
  int index_pt_bispectrum_s_12d3_010;
  int index_pt_bispectrum_s_12d3_001;
  int index_pt_bispectrum_s_12d3_110;
  int index_pt_bispectrum_s_12d3_101;
  int index_pt_bispectrum_s_12d3_011;
  int index_pt_bispectrum_s_12d3_111;

  int index_pt_bispectrum_s_1d23d_000;
  int index_pt_bispectrum_s_1d23d_100;
  int index_pt_bispectrum_s_1d23d_010;
  int index_pt_bispectrum_s_1d23d_001;
  int index_pt_bispectrum_s_1d23d_110;
  int index_pt_bispectrum_s_1d23d_101;
  int index_pt_bispectrum_s_1d23d_011;
  int index_pt_bispectrum_s_1d23d_111;

     
    
  int index_pt_bispectrum_s_1d23_000;
  int index_pt_bispectrum_s_1d23_100;
  int index_pt_bispectrum_s_1d23_010;
  int index_pt_bispectrum_s_1d23_001;
  int index_pt_bispectrum_s_1d23_110;
  int index_pt_bispectrum_s_1d23_101;
  int index_pt_bispectrum_s_1d23_011;
  int index_pt_bispectrum_s_1d23_111;



  
  int index_pt_bispectrum_s_1d2d3_000;
  int index_pt_bispectrum_s_1d2d3_100;
  int index_pt_bispectrum_s_1d2d3_010;
  int index_pt_bispectrum_s_1d2d3_001;
  int index_pt_bispectrum_s_1d2d3_110;
  int index_pt_bispectrum_s_1d2d3_101;
  int index_pt_bispectrum_s_1d2d3_011;
  int index_pt_bispectrum_s_1d2d3_111;

     
      
  
  
  
  int pt_size;            /**< size of perturbation vector */

  double * y;             /**< vector of perturbations_lqc to be integrated */
  double * dy;            /**< time-derivative of the same vector */

  int * used_in_sources; /**< boolean array specifying which
                            perturbations_lqc enter in the calculation of
                            source functions */

};


/**
 * Workspace containing, among other things, the value at a given time
 * of all background/perturbed quantities, as well as their indices.
 * There will be one such structure created for each mode
 * (scalar/.../tensor) and each thread (in case of parallel computing)
 */

struct perturb_lqc_workspace
{

  /** @name - all possible useful indices for those metric
      perturbations_lqc which are not integrated over time, but just
      inferred from Einstein equations. "_mt_" stands for "metric".*/

  //@{

  int index_mt_phi_scf_prime_prime_re;    
  int index_mt_phi_scf_prime_prime_re_1;    
  int index_mt_phi_scf_prime_prime_re_2;    
  int index_mt_phi_scf_prime_prime_re_3;    
  int index_mt_gw_prime_prime_re;

  int index_mt_phi_scf_prime_prime_im;    
  int index_mt_phi_scf_prime_prime_im_1;    
  int index_mt_phi_scf_prime_prime_im_2;    
  int index_mt_phi_scf_prime_prime_im_3;    
  int index_mt_gw_prime_prime_im;
  
  int index_mt_phi_scf_prime_re;    
  int index_mt_phi_scf_prime_im;    

  int index_mt_phi_scf_prime_re_1;    
  int index_mt_phi_scf_prime_im_1;    

  int index_mt_phi_scf_prime_re_2;    
  int index_mt_phi_scf_prime_im_2;    

  int index_mt_phi_scf_prime_re_3;    
  int index_mt_phi_scf_prime_im_3;    

  int index_mt_gwdot_re;    
  int index_mt_gwdot_im;    

  int index_mt_phi_scf_abs;    
  int index_mt_gw_abs;


  int index_mt_gw_000;
  int index_mt_gw_100;
  int index_mt_gw_110;
  int index_mt_gw_111;


  int index_mt_phi_scf_re;    
  int index_mt_phi_scf_re_1;    
  int index_mt_phi_scf_re_2;    
  int index_mt_phi_scf_re_3;    
  int index_mt_gw_re;
  
  int index_mt_phi_scf_im;    
  int index_mt_phi_scf_im_1;    
  int index_mt_phi_scf_im_2;    
  int index_mt_phi_scf_im_3;    
  int index_mt_gw_im;

  int index_mt_pk_s;
  int index_mt_pk_s_1;
  int index_mt_pk_s_2;
  int index_mt_pk_s_3;
  int index_mt_pk_t;
  
  int index_mt_bispectrum_s_field_redef;


  int index_mt_R_s_re;    
  int index_mt_R_s_im;

  int index_mt_R_s_prime_re;    
  int index_mt_R_s_prime_im;


  int index_mt_v_s_re;    
  int index_mt_v_s_im;

  int index_mt_v_s_prime_re;    
  int index_mt_v_s_prime_im;


  int index_mt_v_t_re;    
  int index_mt_v_t_im;

  int index_mt_v_t_prime_re;    
  int index_mt_v_t_prime_im;


  
  int index_mt_phi_prefactor_123_000;
  int index_mt_phi_123_000;
  int index_mt_phi_123_100;
  int index_mt_phi_123_010;
  int index_mt_phi_123_001;
  int index_mt_phi_123_110;
  int index_mt_phi_123_101;
  int index_mt_phi_123_011;
  int index_mt_phi_123_111;


  
  int index_mt_phi_1d23_000;
  int index_mt_phi_1d23_100;
  int index_mt_phi_1d23_010;
  int index_mt_phi_1d23_001;
  int index_mt_phi_1d23_110;
  int index_mt_phi_1d23_101;
  int index_mt_phi_1d23_011;
  int index_mt_phi_1d23_111;



  
  int index_mt_phi_12d3_000;
  int index_mt_phi_12d3_100;
  int index_mt_phi_12d3_010;
  int index_mt_phi_12d3_001;
  int index_mt_phi_12d3_110;
  int index_mt_phi_12d3_101;
  int index_mt_phi_12d3_011;
  int index_mt_phi_12d3_111;


  
  int index_mt_phi_123d_000;
  int index_mt_phi_123d_100;
  int index_mt_phi_123d_010;
  int index_mt_phi_123d_001;
  int index_mt_phi_123d_110;
  int index_mt_phi_123d_101;
  int index_mt_phi_123d_011;
  int index_mt_phi_123d_111;


  
  int index_mt_phi_1d2d3_000;
  int index_mt_phi_1d2d3_100;
  int index_mt_phi_1d2d3_010;
  int index_mt_phi_1d2d3_001;
  int index_mt_phi_1d2d3_110;
  int index_mt_phi_1d2d3_101;
  int index_mt_phi_1d2d3_011;
  int index_mt_phi_1d2d3_111;


  int index_mt_phi_1d23d_000;
  int index_mt_phi_1d23d_100;
  int index_mt_phi_1d23d_010;
  int index_mt_phi_1d23d_001;
  int index_mt_phi_1d23d_110;
  int index_mt_phi_1d23d_101;
  int index_mt_phi_1d23d_011;
  int index_mt_phi_1d23d_111;


  int index_mt_phi_12d3d_000;
  int index_mt_phi_12d3d_100;
  int index_mt_phi_12d3d_010;
  int index_mt_phi_12d3d_001;
  int index_mt_phi_12d3d_110;
  int index_mt_phi_12d3d_101;
  int index_mt_phi_12d3d_011;
  int index_mt_phi_12d3d_111;


  
  int index_mt_bispectrum_t_000;
  int index_mt_bispectrum_t_100;
  int index_mt_bispectrum_t_110;
  int index_mt_bispectrum_t_111;


  int index_mt_bispectrum_s_000;
  int index_mt_bispectrum_s_100;
  int index_mt_bispectrum_s_110;
  int index_mt_bispectrum_s_111;



  int index_mt_bispectrum_s_prefactor_123_000;
  int index_mt_bispectrum_s_123_000;
  int index_mt_bispectrum_s_123_100;
  int index_mt_bispectrum_s_123_010;
  int index_mt_bispectrum_s_123_001;
  int index_mt_bispectrum_s_123_110;
  int index_mt_bispectrum_s_123_101;
  int index_mt_bispectrum_s_123_011;
  int index_mt_bispectrum_s_123_111;



  int index_mt_bispectrum_s_1d23_000;
  int index_mt_bispectrum_s_1d23_100;
  int index_mt_bispectrum_s_1d23_010;
  int index_mt_bispectrum_s_1d23_001;
  int index_mt_bispectrum_s_1d23_110;
  int index_mt_bispectrum_s_1d23_101;
  int index_mt_bispectrum_s_1d23_011;
  int index_mt_bispectrum_s_1d23_111;

  
    
  int index_mt_bispectrum_s_123d_000;
  int index_mt_bispectrum_s_123d_100;
  int index_mt_bispectrum_s_123d_010;
  int index_mt_bispectrum_s_123d_001;
  int index_mt_bispectrum_s_123d_110;
  int index_mt_bispectrum_s_123d_101;
  int index_mt_bispectrum_s_123d_011;
  int index_mt_bispectrum_s_123d_111;

  int index_mt_bispectrum_s_12d3d_000;
  int index_mt_bispectrum_s_12d3d_100;
  int index_mt_bispectrum_s_12d3d_010;
  int index_mt_bispectrum_s_12d3d_001;
  int index_mt_bispectrum_s_12d3d_110;
  int index_mt_bispectrum_s_12d3d_101;
  int index_mt_bispectrum_s_12d3d_011;
  int index_mt_bispectrum_s_12d3d_111;

    
  int index_mt_bispectrum_s_12d3_000;
  int index_mt_bispectrum_s_12d3_100;
  int index_mt_bispectrum_s_12d3_010;
  int index_mt_bispectrum_s_12d3_001;
  int index_mt_bispectrum_s_12d3_110;
  int index_mt_bispectrum_s_12d3_101;
  int index_mt_bispectrum_s_12d3_011;
  int index_mt_bispectrum_s_12d3_111;




    
  int index_mt_bispectrum_s_1d2d3_000;
  int index_mt_bispectrum_s_1d2d3_100;
  int index_mt_bispectrum_s_1d2d3_010;
  int index_mt_bispectrum_s_1d2d3_001;
  int index_mt_bispectrum_s_1d2d3_110;
  int index_mt_bispectrum_s_1d2d3_101;
  int index_mt_bispectrum_s_1d2d3_011;
  int index_mt_bispectrum_s_1d2d3_111;

  
  int index_mt_bispectrum_s_1d23d_000;
  int index_mt_bispectrum_s_1d23d_100;
  int index_mt_bispectrum_s_1d23d_010;
  int index_mt_bispectrum_s_1d23d_001;
  int index_mt_bispectrum_s_1d23d_110;
  int index_mt_bispectrum_s_1d23d_101;
  int index_mt_bispectrum_s_1d23d_011;
  int index_mt_bispectrum_s_1d23d_111;






  
  //int index_mt_R_t_re;    
  //int index_mt_R_t_im;
  

  int mt_size;                /**< size of metric perturbation vector */

  //@}

  /** @name - value at a given time of all background/perturbed
      quantities
  */

  //@{

  double * pvecback;          /**< background quantities */
  double * pvecmetric;        /**< metric quantities */
  struct perturb_lqc_vector * pv; /**< pointer to vector of integrated
                                 perturbations_lqc and their
                                 time-derivatives */

  
  FILE * perturb_lqc_output_file; /**< filepointer to output file*/
  int index_ikout;            /**< index for output k value (when k_output_values is set) */

  //@}

  /** @name - indices useful for searching background/thermo quantities in tables */

  //@{

  short inter_mode;	/**< flag defining the method used for interpolation background/thermo quantities tables */

  int last_index_back;   /**< the background interpolation function background_at_tau() keeps memory of the last point called through this index */

  //@}

  /** @name - approximations used at a given time */

  //@{


  int approx;     /**<approximation flags holding at a given time*/
  int approx_interval;     /**<approximation flags holding at a given time*/
  int num_approx_intervals;     /**<approximation flags holding at a given time*/
  double tau_ini_for_k;
  double N_exit_for_k;
  double t_exit_for_k;
  double a_exit_for_k;

};

/**
 * Structure pointing towards all what the function that perturb_lqc_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
 */

struct perturb_lqc_parameters_and_workspace {

  struct precision * ppr;         /**< pointer to the precision structure */
  struct background_lqc * pba;        /**< pointer to the background structure */
  struct perturbs_lqc * ppt;          /**< pointer to the precision structure */
  int index_md;                   /**< index of mode (scalar/.../vector/tensor) */
  int index_k;			          /**< index of wavenumber */
  double k;			              /**< current value of wavenumber in 1/Mpc */
  double k_2;
  double k_3;
  struct perturb_lqc_workspace * ppw; /**< workspace defined above */

};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int perturb_lqc_sources_at_tau(
                             struct perturbs_lqc * ppt,
                             int index_md,
                             int index_type,
                             double tau,
                             double * pvecsources
                             );

  int perturb_lqc_init(
                   struct precision * ppr,
                   struct background_lqc * pba,
                   struct perturbs_lqc * ppt
                   );

  int perturb_lqc_free(
                   struct perturbs_lqc * ppt
                   );

  int perturb_lqc_indices_of_perturbs_lqc(
                                  struct precision * ppr,
                                  struct background_lqc * pba,
                                  struct perturbs_lqc * ppt
                                  );

  int perturb_lqc_timesampling_for_sources(
                                       struct precision * ppr,
                                       struct background_lqc * pba,
                                       struct perturbs_lqc * ppt
                                       );
  int perturb_lqc_get_k_list(
                         struct precision * ppr,
                         struct background_lqc * pba,
                         struct perturbs_lqc * ppt
                         );

  int perturb_lqc_workspace_init(
                             struct precision * ppr,
                             struct background_lqc * pba,
                             struct perturbs_lqc * ppt,
                             int index_md,
                             struct perturb_lqc_workspace * ppw
                             );

  int perturb_lqc_workspace_free(
                             struct perturbs_lqc * ppt,
                             int index_md,
                             struct perturb_lqc_workspace * ppw
                             );

  int perturb_lqc_solve(
                    struct precision * ppr,
                    struct background_lqc * pba,
                    struct perturbs_lqc * ppt,
                    int index_md,
                    int index_k_1,
                    int index_k_2,
                    int index_k_3,
                    int index_x_2_x_3,
                    struct perturb_lqc_workspace * ppw
                    );


  int perturb_lqc_vector_init(
                          struct precision * ppr,
                          struct background_lqc * pba,
                          struct perturbs_lqc * ppt,
                          int index_md,
                          double k1,
                          double k2,
                          double k3,
                          double tau,
                          struct perturb_lqc_workspace * ppw,
                          int pa_old
                          );

  int perturb_lqc_vector_free(
                          struct perturb_lqc_vector * pv
                          );

  int perturb_lqc_initial_conditions(
                                 struct precision * ppr,
                                 struct background_lqc * pba,
                                 struct perturbs_lqc * ppt,
                                 int index_md,
                                 double k,
                                 double k2,
                                 double k3,
                                 double tau,
                                 struct perturb_lqc_workspace * ppw
                                 );

  int perturb_lqc_approximations(
                             struct precision * ppr,
                             struct background_lqc * pba,
                             struct perturbs_lqc * ppt,
                             int index_md,
                             double k,
                             double tau,
                             struct perturb_lqc_workspace * ppw
                             );

  int perturb_lqc_timescale(
                        double tau,
                        void * parameters_and_workspace,
                        double * timescale,
                        ErrorMsg error_message
                        );

  int perturb_lqc_einstein(
                       struct precision * ppr,
                       struct background_lqc * pba,
                       struct perturbs_lqc * ppt,
                       int index_md,
                       double k,
                       double k2,
                       double k3,
                       double tau,
                       double * y,
                       struct perturb_lqc_workspace * ppw
                       );

  int perturb_lqc_sources(
                      double tau,
                      double * pvecperturbations_lqc,
                      double * pvecderivs,
                      int index_tau,
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      );

  int perturb_lqc_print_variables(
                              double tau,
                              double * y,
                              double * dy,
                              void * parameters_and_workspace,
                              ErrorMsg error_message
                              );

  int perturb_lqc_derivs(
                     double tau,
                     double * y,
                     double * dy,
                     void * parameters_and_workspace,
                     ErrorMsg error_message
                     );



  int perturb_lqc_prepare_output_file(struct background_lqc * pba,
                                  struct perturbs_lqc * ppt,
                                  struct perturb_lqc_workspace * ppw,
                                  int index_ikout,
                                  int index_md);

  int perturb_lqc_prepare_output(struct background_lqc * pba,
                             struct perturbs_lqc * ppt);


  int output_one_line_of_pk_lqc(
                          FILE * pkfile,
                          double one_k,
                          double one_pk
                                );
int output_open_pk_file_lqc(
                        struct background_lqc * pba,
                        struct perturbs_lqc * ppt,
                        FILE * * pkfile,
                        FileName filename,
                        char * first_line
                            );

  
int output_pk_lqc(
              struct background_lqc * pba,
              struct perturbs_lqc * ppt,
              int index_md);

  
  int spectra_pk_lqc(
               struct background_lqc * pba,
               struct perturbs_lqc * ppt,
               int index_md
                     );

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
