/** @file background.h Documented includes for background module */

#ifndef __BACKGROUND_LQC__
#define __BACKGROUND_LQC__

#include "common.h"
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

enum spatial_curvature_lqc {flat_lqc,open_lqc,closed_lqc};

/**
 * All background parameters and evolution that other modules need to know.
 *
 * Once initialized by the backgound_init(), contains all necessary
 * information on the background evolution (except thermodynamics),
 * and in particular, a table of all background quantities as a
 * function of time and scale factor, used for interpolation in other
 * modules.
 */

struct background_lqc
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' structure)
   *
   * The background_lqc cosmological parameters listed here form a parameter
   * basis which is directly usable by the background_lqc module. Nothing
   * prevents from defining the input cosmological parameters
   * differently, and to pre-process them into this format, using the input
   * module (this might require iterative calls of background_lqc_init()
   * e.g. for dark energy or decaying dark matter). */

  //@{

  double H0; /**< \f$ H_0 \f$: Hubble parameter (in fact, [\f$H_0/c\f$]) in \f$ Mpc^{-1} \f$ */

  double x_ini;       /**< \f$ \phi(t_0) \f$: scalar field initial value */
  double phi_ini;       /**< \f$ \phi(t_0) \f$: scalar field initial value */
  double phi_dot_ini;       /**< \f$ \phi(t_0) \f$: scalar field initial value */

  double Omega0_scf; /**< \f$ \Omega_{0_k} \f$: scalar field omega */

  /** @name - related parameters */

  //@{
  double epsilon_2;

  double h; /**< reduced Hubble parameter */
  double age; /**< age in Gyears */

  //@}

  /** @name - other background_lqc parameters */

  //@{

  //double a_today; /**< scale factor today (arbitrary and irrelevant for most purposes) */
  double a_bounce;
  int flag_start_inflation;
  int flag_start_deflation;
  int flag_end_inflation;
  int flag_end_deflation;
  double tau_start_inflation;
  double N_start_inflation;
  double a_start_inflation;
  double tau_start_remote_past;
  double a_start_remote_past;
  double N_start_remote_past;
  double tau_start_deflation;
  double a_start_deflation;
  double N_start_deflation;
  double tau_end_inflation;
  double a_end_inflation;
  double tau_end_remote_past;
  double a_end_remote_past;
  double N_end_remote_past;
  double tau_end_deflation;
  double a_end_deflation;
  double N_end_deflation;
  double N_exit;
  double N_end;
  double tau_exit;
  double integration_depth_ini;
  double alpha_max_for_pk_lqc;
  double alpha_min_for_pk_lqc;
  //@}
  int has_perturbations_lqc;

  double H_ini_remote_past;
  double N_ini_remote_past;
  double rho_bounce;
  double l0_lqc;
  double Delta0_lqc;
  double rho_ini;
  double m_scf_lqc;
  double A_s; //amplitude of the scalar power spectrum at N_star;
  double k_star; //Planck pivot scale
  double N_star; //Number of efolds from horizon exit of the pivot scale to the end of inflation. 
  double phi_star; //potential at horizon exit of the pivot scale.
  double Gamma;
  double H_max;
  double ALPHA_end;
  double ALPHA_ini;
  
  double tau_end;

  int future_branch;
  /** @name - all indices for the vector of background_lqc (=bg) quantities stored in table */

  //@{
  FileName root; /**< root for all file names */

  int index_bg_a;             /**< scale factor */
  int index_bg_H;             /**< Hubble parameter in \f$Mpc^{-1}\f$ */
  int index_bg_v;
  int index_bg_h;
  int index_bg_pia;
  int index_bg_piadot;

  /* end of vector in short format, now quantities in normal format */


  int index_bg_phi_scf;       /**< scalar field value */
  int index_bg_phi_prime_scf; /**< scalar field derivative wrt conformal time */
  int index_bg_Pphi_scf;
  int index_bg_Pphidot_scf;
 
  int index_bg_rho_scf;       /**< scalar field energy density */
  int index_bg_p_scf;         /**< scalar field pressure */
  
  int index_bg_H_prime;

  int index_bg_x;
  int index_bg_y;
  
  int index_bg_mu;
  int index_bg_t_times_epsilon_H;
  int index_bg_t;

  int index_bg_V;
  int index_bg_V_phi;
  int index_bg_V_phiphi;
  int index_bg_V_phiphiphi;
  
  int index_bg_epsilon_H;
  int index_bg_effective_potential;
  int index_bg_epsilon_V;
  int index_bg_epsilon_2;

  int index_bg_f_S;
  int index_bg_f_T;
  double k_min_for_pk;
  double lambda_k_2;
  double lambda_k_3;
  
  int index_bg_f_S_Q;
  int index_bg_H2_f_S_Q;
  int index_bg_H2_f_S_Q_phys;
  int index_bg_zS_primeprime_over_zS;
  int index_bg_zT_primeprime_over_zT;
  int index_bg_f_S_R;

  int index_bg_f_T_Q;
  int index_bg_f_T_R;
  
  int index_bg_z_T;
  int index_bg_z_S;
 
  int index_bg_time;          /**< proper (cosmological) time in Mpc */

  int bg_size_short;  /**< size of background_lqc vector in the "short format" */
  int bg_size_normal; /**< size of background_lqc vector in the "normal format" */
  int bg_size;        /**< size of background_lqc vector in the "long format" */

  //@}

  /** @name - background_lqc interpolation tables */

  //@{

  int bt_size;               /**< number of lines (i.e. time-steps) in the array */
  int bt_contracting_phase_size; 
  int bt_remote_past_size;
  int bt_expanding_phase_size;
  
  double * tau_table;      
  double * N_table;       
  double * z_table;        
  double * background_lqc_table; 
  double * background_lqc_table_N; 
  double * d2tau_dz2_table;
  double * d2tau_dN2_table; 
  double * d2N_dtau2_table; 
  double * d2background_lqc_dtau2_table; 
  double * d2background_lqc_dN2_table;

  double * tau_table_expanding;      
  double * N_table_expanding;       
  double * background_lqc_table_expanding; 
  double * d2tau_dN2_table_expanding; 
  double * d2N_dtau2_table_expanding; 
  double * d2background_lqc_dtau2_table_expanding; 
  double * d2background_lqc_dN2_table_expanding;

//future branch (fb)
  int bt_fb_size;
  int bt_fb_contracting_phase_size; 
  int bt_fb_remote_past_size;
  int bt_fb_expanding_phase_size;
 
  double * tau_table_fb;      
  double * N_table_fb;       
  double * z_table_fb;        
  double * background_lqc_table_fb; 
  double * background_lqc_table_N_fb; 
  double * d2tau_dz2_table_fb;
  double * d2tau_dN2_table_fb; 
  double * d2N_dtau2_table_fb; 
  double * d2background_lqc_dtau2_table_fb; 
  double * d2background_lqc_dN2_table_fb;

  double * tau_table_expanding_fb;      
  double * N_table_expanding_fb;       
  double * background_lqc_table_expanding_fb; 
  double * d2tau_dN2_table_expanding_fb; 
  double * d2N_dtau2_table_expanding_fb; 
  double * d2background_lqc_dtau2_table_expanding_fb; 
  double * d2background_lqc_dN2_table_expanding_fb;

  
// past branch(pb)
  int bt_pb_size;
  int bt_pb_contracting_phase_size; 
  int bt_pb_remote_past_size;
  int bt_pb_expanding_phase_size;

  double * tau_table_pb;      
  double * N_table_pb;       
  double * z_table_pb;        
  double * background_lqc_table_pb; 
  double * background_lqc_table_N_pb; 
  double * d2tau_dz2_table_pb;
  double * d2tau_dN2_table_pb; 
  double * d2N_dtau2_table_pb; 
  double * d2background_lqc_dtau2_table_pb; 
  double * d2background_lqc_dN2_table_pb;

  double * tau_table_expanding_pb;      
  double * N_table_expanding_pb;       
  double * background_lqc_table_expanding_pb; 
  double * d2tau_dN2_table_expanding_pb; 
  double * d2N_dtau2_table_expanding_pb; 
  double * d2background_lqc_dtau2_table_expanding_pb; 
  double * d2background_lqc_dN2_table_expanding_pb;

  //@}


  /** @name - all indices for the vector of background_lqc quantities to be integrated (=bi)
   *
   * Most background_lqc quantities can be immediately inferred from the
   * scale factor. Only few of them require an integration with
   * respect to conformal time (in the minimal case, only one quantity needs to
   * be integrated with time: the scale factor, using the Friedmann
   * equation). These indices refer to the vector of
   * quantities to be integrated with time.
   * {B} quantities are needed by background_lqc_functions() while {C} quantities are not.
   */

  //@{

  int index_bi_a;       /**< {B} scale factor */
  int index_bi_conformal_time;       /**< {B} scale factor */
  int index_bg_conformal_time;       /**< {B} scale factor */
  int index_bi_cosmic_time;       /**< {B} scale factor */
  int index_bg_cosmic_time;       /**< {B} scale factor */
  int index_bi_N;       /**< {B} scale factor */
  int index_bg_N;       /**< {B} scale factor */

//  int index_bi_x;       /**< {B} scale factor */
//  int index_bi_y;       /**< {B} scale factor */ 
  int index_bi_phi;       /**< {B} scale factor */
//  int index_bi_phi_dot;       /**< {B} scale factor */
  int index_bi_Pphi;       /**< {B} scale factor */
  int index_bi_v;       /**< {B} scale factor */
  int index_bi_h;       /**< {B} scale factor */

  int index_bi_H;       /**< {B} scale factor */

  //int index_bi_time;    

  
  int index_bi_tau;     


  int bi_B_size;        /**< Number of {B} parameters */
  int bi_size;          /**< Number of {B}+{C} parameters */

  //@}

  /** @name - flags describing the absence or presence of cosmological
      ingredients
      *
      * having one of these flag set to zero allows to skip the
      * corresponding contributions, instead of adding null contributions.
      */


  //@{

  short has_scf;       /**< presence of a scalar field? */
  short has_lqc;       /**< GR or LQC*/
  short has_contracting_phase;       /**< presence of a scalar field? */
  short post_bounce;
  short remote_past;
  //@}

  

  /**
   *@name - some flags needed for calling background_lqc functions
   */

  //@{

  short short_info;  /**< flag for calling background_lqc_at_eta and return little information */
  short normal_info; /**< flag for calling background_lqc_at_eta and return medium information */
  short long_info;   /**< flag for calling background_lqc_at_eta and return all information */

  short inter_normal;  /**< flag for calling background_lqc_at_eta and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling background_lqc_at_eta and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short shooting_failed;  /**< flag is set to true if shooting failed. */

  ErrorMsg shooting_error; /**< Error message from shooting failed. */

  short background_lqc_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/**
 * temporary parameters and workspace passed to the background_lqc_derivs function
 */

struct background_lqc_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background_lqc * pba;

  /* workspace */
  double * pvecback;

};

/**
 * temporary parameters and workspace passed to phase space distribution function
 */


/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_lqc_at_tau(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);
  int background_lqc_at_N(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);
  int background_lqc_at_tau_pb(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);
  int background_lqc_at_N_pb(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);
  int background_lqc_at_tau_fb(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);
  int background_lqc_at_N_fb(
			struct background_lqc *pba,
			double tau,
			short return_format,
			short inter_mode,
			int * last_index,
			double * pvecback
			);

  int background_lqc_functions(
			   struct background_lqc *pba,
			   double * pvecback_B,
			   short return_format,
			   double * pvecback
			   );

  int background_lqc_tau_of_z(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );
  int background_lqc_tau_of_N(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );

  int background_lqc_N_of_tau(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );
  int background_lqc_tau_of_N_pb(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );

  int background_lqc_N_of_tau_pb(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );
  int background_lqc_tau_of_N_fb(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );

  int background_lqc_N_of_tau_fb(
			  struct background_lqc *pba,
			  double z,
			  double * tau
			  );

  
  int background_lqc_init(
		      struct precision *ppr,
		      struct background_lqc *pba
		      );

  int background_lqc_free(
		      struct background_lqc *pba
		      );

  int background_lqc_free_input(
                            struct background_lqc *pba
                            );

  int background_lqc_indices(
			 struct background_lqc *pba
			 );




  int background_lqc_solve(
		       struct precision *ppr,
		       struct background_lqc *pba
		       );

  int background_lqc_initial_conditions(
				    struct precision *ppr,
				    struct background_lqc *pba,
				    double * pvecback,
				    double * pvecback_integration
				    );

  int background_lqc_output_titles(struct background_lqc * pba,
                               char titles[_MAXTITLESTRINGLENGTH_]
                               );

  int background_lqc_output_data(
                           struct background_lqc *pba,
                           int number_of_titles,
                           double *data);

  int background_lqc_derivs(
			 double z,
			 double * y,
			 double * dy,
			 void * parameters_and_workspace,
			 ErrorMsg error_message
			 );

 
#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name Some conversion factors and fundamental constants needed by background_lqc module:
 */

//@{

#define _gamma_lqc_  2375.e-4/**< Babero-Immirzi parameter */

#define _Gyr_over_Mpc_ 3.06601394e2 /**< conversion factor from megaparsecs to gigayears
				         (c=1 units, Julian years of 365.25 days) */

#endif
/* @endcond */
