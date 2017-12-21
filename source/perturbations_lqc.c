/** @file perturbations_lqc.c Documented perturbation module
 *
 * Julien Lesgourgues, 23.09.2010
 *
 * Deals with the perturbation evolution.
 * This module has two purposes:
 *
 * - at the beginning; to initialize the perturbations_lqc, i.e. to
 * integrate the perturbation equations, and store temporarily the terms
 * contributing to the source functions as a function of conformal
 * time. Then, to perform a few manipulations of these terms in order to
 * infer the actual source functions \f$ S^{X} (k, \tau) \f$, and to
 * store them as a function of conformal time inside an interpolation
 * table.
 *
 * - at any time in the code; to evaluate the source functions at a
 * given conformal time (by interpolating within the interpolation
 * table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# perturb_lqc_init() at the beginning (but after background_lqc_init() and thermodynamics_init())
 * -# perturb_lqc_sources_at_tau() at any later time
 * -# perturb_lqc_free() at the end, when no more calls to perturb_lqc_sources_at_tau() are needed
 */

#include "perturbations_lqc.h"


/**
 * Source function \f$ S^{X} (k, \tau) \f$ at a given conformal time tau.
 *
 * Evaluate source functions at given conformal time tau by reading
 * the pre-computed table and interpolating.
 *
 * @param ppt        Input: pointer to perturbation structure containing interpolation tables
 * @param index_md   Input: index of requested mode
 * @param index_type Input: index of requested source function type
 * @param tau        Input: any value of conformal time
 * @param psource    Output: vector (already allocated) of source function as a function of k
 * @return the error status
 */

int perturb_lqc_sources_at_tau(
                           struct perturbs_lqc * ppt,
                           int index_md,
                           int index_type,
                           double tau,
                           double * psource
                           ) {


  /** Summary: */

  /** - interpolate in pre-computed table contained in ppt */
  class_call(array_interpolate_two_bis(ppt->tau_sampling,
                                       1,
                                       0,
                                       ppt->sources[index_md][index_type],
                                       ppt->k_size[index_md],
                                       ppt->tau_size,
                                       tau,
                                       psource,
                                       ppt->k_size[index_md],
                                       ppt->error_message),
             ppt->error_message,
             ppt->error_message);


  return _SUCCESS_;
}

/**
 * Initialize the perturbs_lqc structure, and in particular the table of source functions.
 *
 * Main steps:
 *
 * - given the values of the flags describing which kind of
 *   perturbations_lqc should be considered (modes: scalar/vector/tensor,
 *   initial conditions, type of source functions needed...),
 *   initialize indices and wavenumber list
 *
 * - define the time sampling for the output source functions
 *
 * - for each mode (scalar/vector/tensor): initialize the indices of
 *   relevant perturbations_lqc, integrate the differential system,
 *   compute and store the source functions.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background_lqc structure
 * @param ppt Output: Initialized perturbation structure
 * @return the error status
 */

int perturb_lqc_init(
                 struct precision * ppr,
                 struct background_lqc * pba,
                 struct perturbs_lqc * ppt
                 ) {

  /** Summary: */

  /** - define local variables */
    ppt->tol=ppr->tol_perturb_lqc_integration;

  /* running index for modes */
  int index_md;
  /* running index for initial conditions */
  /* running index for wavenumbers */
  int index_k;
  int index_k_1;
  int index_k_2;
  int index_k_3;
  int index_x_2_x_3 = 0;

                      int index_k_3_min,index_k_3_max;

  
  /* pointer to one struct perturb_lqc_workspace per thread (one if no openmp) */
  struct perturb_lqc_workspace ** pppw;
  /* number of threads (always one if no openmp) */
  int number_of_threads=1;
  /* index of the thread (always 0 if no openmp) */
  int thread=0;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" just after leaving the
     parallel region. */
  int abort;

  /* unsigned integer that will be set to the size of the workspace */
  size_t sz;





  
#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop, tspent;
#endif

  /** - perform preliminary checks */

  if (ppt->has_perturbations_lqc == _FALSE_) {
    if (ppt->perturbations_lqc_verbose > 0)
      printf("No LQC power spectrum requested. LQC module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ppt->perturbations_lqc_verbose > 0)
      printf("Computing perturbations LQC\n");
  }

 ppt->k_min_for_pk = pba->alpha_min_for_pk_lqc*pba->k_star;
 ppt->k_max_for_pk = pba->alpha_max_for_pk_lqc*pba->k_star;




   
  FILE * pkfile;     
  FileName file_name;
  int colnum = 1;


  //double one_fNL;
 //double one_x2;
 //double one_x3;
 //int index;
 //printf("x_size_end=%d",ppt->x_size[index_md]);

  sprintf(file_name,"%s_%s",pba->root,"bispectrum_s.dat");

  class_open(pkfile,file_name,"w",pba->error_message);

    fprintf(pkfile,"#Bispectrum  %s\n","scalars");
    fprintf(pkfile,"#");
    class_fprintf_columntitle(pkfile,"index_x2",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"index_x3",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"x2",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"x3",_TRUE_,colnum);
    class_fprintf_columntitle(pkfile,"fNL",_TRUE_,colnum);
   fprintf(pkfile," ");
     fprintf(pkfile,"\n");
     fprintf(pkfile,"\n");
     //fprintf(pkfile," ");
       fclose(pkfile);


     /*
     
    for (index=0; index< index_x_2_x_3; index++){

      one_x2 =ppt->table_x2[index_md][index];
      one_x3 =ppt->table_x3[index_md][index];
  one_fNL = ppt->table_fNL[index_md][index];;

  fprintf(pkfile," ");
  class_fprintf_double(pkfile,one_x2,_TRUE_);
  class_fprintf_double(pkfile,one_x3,_TRUE_);
  class_fprintf_double(pkfile,one_fNL,_TRUE_);
  fprintf(pkfile,"\n");

    }


       fclose(pkfile);

   */


 
/*find window of wavenumbers*/
 double tau_max=pba->tau_end_inflation;
    //pba->tau_start_inflation
    //            +pba->N_end
    //            -pba->N_start_inflation;
  double tau_min=pba->tau_start_inflation;
  double tau_half;
  double aH_sample=_HUGE_;
  double error_aH = .01;
  int index_back;
    double * pvecback2;

    class_alloc(pvecback2,pba->bg_size_normal*sizeof(double),pba->error_message);



    if (ppt->tau_ini == 0. /*&& ppt->tau_end == 0.*/){
     printf("Computing Inflationary Perturbations\n");

  class_call(background_lqc_at_tau(pba,
                                 tau_min,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
    aH_sample = pvecback2[pba->index_bg_a]*pvecback2[pba->index_bg_H];
    // printf("depth = %e \n",ppt->integration_depth_ini);
    // printf("depth = %e \n",ppt->integration_depth_end);

  if (ppt->k_min_for_pk/ppt->integration_depth_ini <  aH_sample){
    ppt->tau_ini = pba->tau_start_inflation;
  printf("Note that at the start of inflation k_min/aH = %e \n",
         ppt->k_min_for_pk/aH_sample);
  printf("If you want proper inflationary results,  increase the x_ini in order to achieve a longer inflationary phase.\n");

 }
else {
  while (fabs(aH_sample-ppt->k_min_for_pk
              /ppt->integration_depth_ini)>
         error_aH*ppt->k_min_for_pk
              /ppt->integration_depth_ini){

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
  if (aH_sample>ppt->k_min_for_pk/ppt->integration_depth_ini)
    tau_max = tau_half;
  else tau_min = tau_half;
  ppt->tau_ini = tau_half;
 }
  if (ppt->tau_ini == 0.)  ppt->tau_ini = pba->tau_start_inflation;
 }
  //ppt->tau_ini = pba->tau_start_inflation;
  
  
aH_sample=_HUGE_;
 tau_max=pba->tau_end_inflation;
  //pba->tau_start_inflation
  //              +pba->N_end
  //              -pba->N_start_inflation;
tau_min=pba->tau_start_inflation;
 
  while (fabs(aH_sample-ppt->k_max_for_pk*ppt->integration_depth_end)>
         error_aH*ppt->k_max_for_pk*ppt->integration_depth_end){

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
  if (aH_sample>ppt->k_max_for_pk*ppt->integration_depth_end) tau_max = tau_half;
  else tau_min = tau_half;
  ppt->tau_end = tau_half;
 }

 }
  
 if (ppt->tau_end == 0.){
   //printf("Computing tau_end\n");
aH_sample=_HUGE_;
 tau_max=pba->tau_end_inflation;
  //pba->tau_start_inflation
  //              +pba->N_end
  //               -pba->N_start_inflation;
tau_min=pba->tau_start_inflation;
 
  while (fabs(aH_sample-ppt->k_max_for_pk*ppt->integration_depth_end)>
         error_aH*ppt->k_max_for_pk
              *ppt->integration_depth_end){

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
  if (aH_sample>ppt->k_max_for_pk*ppt->integration_depth_end) tau_max = tau_half;
  else tau_min = tau_half;
  ppt->tau_end = tau_half;
 }

 
   }
 // else ppt->tau_end = pba->tau_end_inflation;
 //else ppt->tau_end = pba->tau_start_inflation;//???
//ppt->tau_end = -1.e-5;//pba->tau_start_inflation;

 //free(pvecback2);
printf("tau_end = %e\n",ppt->tau_end);


 
 
  



  /** - initialize all indices and lists in perturbs_lqc structure using perturb_lqc_indices_of_perturbs_lqc() */

  class_call(perturb_lqc_indices_of_perturbs_lqc(ppr,
                                         pba,
                                         ppt),
             ppt->error_message,
             ppt->error_message);
  


  /** - define the common time sampling for all sources using
      perturb_lqc_timesampling_for_sources() */
  
  class_call(perturb_lqc_timesampling_for_sources(ppr,
                                              pba,
                                              ppt),
             ppt->error_message,
             ppt->error_message);
  
  /** - if we want to store perturbations_lqc, write titles and allocate storage */
  class_call(perturb_lqc_prepare_output(pba,ppt),
             ppt->error_message,
             ppt->error_message);


  /** - create an array of workspaces in multi-thread case */





  
#ifdef _OPENMP

#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

  class_alloc(pppw,number_of_threads * sizeof(struct perturb_lqc_workspace *),ppt->error_message);

  /** - loop over modes (scalar, tensors, etc). For each mode: */

//  for (index_md = 0; index_md < ppt->md_size; index_md++) {
  for (index_md = 0; index_md <  1; index_md++) { //V.Sreenath 060617

    if (ppt->perturbations_lqc_verbose > 1)
      printf("Evolving mode %d/%d\n",index_md+1,ppt->md_size);
    //ppt->x_size[index_md] = 0;
    abort = _FALSE_;

    sz = sizeof(struct perturb_lqc_workspace);

#pragma omp parallel                                             \
  shared(pppw,ppr,pba,ppt,index_md,abort,number_of_threads)  \
  private(thread)                                                \
  num_threads(number_of_threads)

    {

#ifdef _OPENMP
      thread=omp_get_thread_num();
#endif

      /** - --> (a) create a workspace (one per thread in multi-thread case) */

      class_alloc_parallel(pppw[thread],sz,ppt->error_message);

      /** - --> (b) initialize indices of vectors of perturbations_lqc with perturb_lqc_indices_of_current_vectors() */

      class_call_parallel(perturb_lqc_workspace_init(ppr,
                                                 pba,
                                                 ppt,
                                                 index_md,
                                                 pppw[thread]),
                          ppt->error_message,
                          ppt->error_message);

    } /* end of parallel region */

    if (abort == _TRUE_) return _FAILURE_;

    /** - --> (c) loop over initial conditions and wavenumbers; for each of them, evolve perturbations_lqc and compute source functions with perturb_lqc_solve() */


        if (ppt->perturbations_lqc_verbose > 1)
          printf("evolving %d wavenumbers\n",ppt->k_size[index_md]);

      abort = _FALSE_;
      
      ppt->second_loop = 0;

      // index_k_1=ppt->k_size[index_md]-1;


      
#pragma omp parallel                                                    \
  shared(pppw,ppr,pba,ppt,index_md,abort,index_k_1,number_of_threads,index_x_2_x_3) \
  private(aH_sample,tau_min,tau_max,tau_half,index_k_3_min,index_k_3_max,index_k,index_k_3,index_k_2,thread,tstart,tstop,tspent) \
  num_threads(number_of_threads)

      {

#ifdef _OPENMP
        thread=omp_get_thread_num();
        tspent=0.;
#endif

#pragma omp for schedule (dynamic)

                for (index_k = ppt->k_size[index_md]-1; index_k >=0; index_k--) {

      
          if ((ppt->perturbations_lqc_verbose > 2) && (abort == _FALSE_)) {
            printf("evolving mode k=%e /Mpc  (%d/%d)",ppt->k[index_md][index_k],index_k+1,ppt->k_size[index_md]);
            printf("\n");
          }

 
                  /*for (index_k_2 = ppt->k_size[index_md]-1; index_k_2 >=0; index_k_2--) {

          
          if(ppt->k[index_md][index_k_2]>0.5*ppt->k[index_md][index_k_1]){
          index_k_3_max = ppt->k_size[index_md]-1;
          index_k_3_min = index_k_2;
              }

          else{
          index_k_3_max = ppt->k_size[index_md]-1;
          index_k_3_min = index_k_3_max-index_k_2;
              }


      
        for (index_k_3 = index_k_3_max; index_k_3 >=index_k_3_min; index_k_3--) {

          printf("k_2 k_3=%d %d\n",index_k_2,index_k_3);
                  */
        
          ppt->tau_ini_for_k = ppt->tau_ini;          
          //printf("\n");
            index_k_1 = index_k;
            index_k_2 = index_k;
            index_k_3 = index_k;
 
          
#ifdef _OPENMP
          tstart = omp_get_wtime();
#endif


          
          class_call_parallel(perturb_lqc_solve(ppr,
                                            pba,
                                            ppt,
                                            index_md,
                                            index_k_1,
                                            index_k_2,
                                            index_k_3,
                                                index_x_2_x_3,
                                            pppw[thread]),
                              ppt->error_message,
                              ppt->error_message);


          
#ifdef _OPENMP
          tstop = omp_get_wtime();

          tspent += tstop-tstart;
#endif

#pragma omp flush(abort)

          } /* end of loop over wavenumbers */
          // }/* end of loop over wavenumbers 2 */
          //}/* end of loop over wavenumbers 3 */




        
#ifdef _OPENMP
        if (ppt->perturbations_lqc_verbose>1)
          printf("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
                 __func__,tspent,omp_get_thread_num());
#endif

      } /* end of parallel region */

      //ppt->x_size[index_md]=index_x_2_x_3;
      if (abort == _TRUE_) return _FAILURE_;

    abort = _FALSE_;
  
#pragma omp parallel                                    \
  shared(pppw,ppt,index_md,abort,number_of_threads)     \
  private(thread)                                       \
  num_threads(number_of_threads)

    {

#ifdef _OPENMP
      thread=omp_get_thread_num();
#endif

      class_call_parallel(perturb_lqc_workspace_free(ppt,index_md,pppw[thread]),
                          ppt->error_message,
                          ppt->error_message);

    } /* end of parallel region */



    


    class_call(output_pk_lqc(pba,
                             ppt,
                             index_md),
               ppt->error_message,
               ppt->error_message);

    if (abort == _TRUE_) return _FAILURE_;
    
  } /* end loop over modes */

  free(pppw);
  free(pvecback2);

  
  return _SUCCESS_;
}

/**
 * Free all memory space allocated by perturb_lqc_init().
 *
 * To be called at the end of each run, only when no further calls to
 * perturb_lqc_sources_at_tau() are needed.
 *
 * @param ppt Input: perturbation structure to be freed
 * @return the error status
 */

int perturb_lqc_free(
                 struct perturbs_lqc * ppt
                 ) {

  int index_md;
  int index_type;
  int filenum;

  if (ppt->has_perturbations_lqc == _TRUE_) {

    for (index_md = 0; index_md < ppt->md_size; index_md++) {


        for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) 
          free(ppt->fields_at_tau_end_lqc[index_md][index_type]);
      free(ppt->fields_at_tau_end_lqc[index_md]);
      free(ppt->k[index_md]);
      /*free(ppt->table_x2[index_md]);
      free(ppt->table_x3[index_md]);
      free(ppt->table_fNL[index_md]);*/
    }

    free(ppt->tau_sampling);

    free(ppt->tp_size);

    free(ppt->k);
    /*free(ppt->table_x2);
    free(ppt->table_x3);
    free(ppt->table_fNL);*/
    free(ppt->ln_k);
    free(ppt->ln_pk);
    
    free(ppt->bispectrum_s_equi);
    free(ppt->bispectrum_t_equi);

    free(ppt->k_size);
    //free(ppt->x_size);

    //free(ppt->sources);

    /** Stuff related to perturbations_lqc output: */

    /** - Free non-NULL pointers */
    if (ppt->index_k_output_values != NULL)
      free(ppt->index_k_output_values);

    for (filenum = 0; filenum<_MAX_NUMBER_OF_K_FILES_; filenum++){
      if (ppt->scalar_perturbations_lqc_data[filenum] != NULL)
        free(ppt->scalar_perturbations_lqc_data[filenum]);

      if (ppt->tensor_perturbations_lqc_data[filenum] != NULL)
        free(ppt->tensor_perturbations_lqc_data[filenum]);
    }

  }

  return _SUCCESS_;

}

/**
 * Initialize all indices and allocate most arrays in perturbs_lqc structure.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */

int perturb_lqc_indices_of_perturbs_lqc(
                                struct precision * ppr,
                                struct background_lqc * pba,
                                struct perturbs_lqc * ppt
                                ) {

  /** Summary: */

  /** - define local variables */

  int index_type;
  int index_md;
  int index_type_common;
  int index_back;

  /** - count modes (scalar, vector, tensor) and assign corresponding indices */

  index_md = 0;
  class_define_index(ppt->index_md_scalars,ppt->has_scalars,index_md,1);
  class_define_index(ppt->index_md_tensors,ppt->has_tensors,index_md,1);
  ppt->md_size = index_md;

  class_test(index_md == 0,
             ppt->error_message,
             "you should have at least one out of {scalars, vectors, tensors} !!!");

  /** - allocate array of number of types for each mode, ppt->tp_size[index_md] */

  class_alloc(ppt->tp_size,ppt->md_size*sizeof(int),ppt->error_message);

  //class_alloc(ppt->sources,ppt->md_size * sizeof(double *),ppt->error_message);

  class_alloc(ppt->pk_lqc,ppt->md_size * sizeof(double *),ppt->error_message);
  class_alloc(ppt->fields_at_tau_end_lqc,ppt->md_size * sizeof(double *),ppt->error_message);


   ppt->has_source_phi_scf = _FALSE_;
   ppt->has_source_gw = _FALSE_;


  index_type = 0;
  
  index_type_common = index_type;


  /*window of wavenumbers*/
  
  double * pvecback;
  class_alloc(pvecback,pba->bg_size_short*sizeof(double),ppt->error_message);

  class_call(background_lqc_at_tau(pba,
                                 ppt->tau_ini,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback),
               pba->error_message,
               ppt->error_message);

   double a_prime_over_a_ini = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];


   
    class_call(background_lqc_at_tau(pba,
                                 ppt->tau_end,
                                 pba->short_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback),
               pba->error_message,
               ppt->error_message);


   double a_prime_over_a_end = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
  

  class_test(ppt->k_min_for_pk>ppt->k_max_for_pk,
             ppt->error_message,
             "you have to increase delta_tau perturbations");

    
     free(pvecback);

  /** - define k values with perturb_lqc_get_k_list() */

  class_call(perturb_lqc_get_k_list(ppr,
                                pba,
                                ppt),
             ppt->error_message,
             ppt->error_message);

  /** - loop over modes. Initialize flags and indices which are specific to each mode. */

  for (index_md = 0; index_md < ppt->md_size; index_md++) {

    /** - (a) scalars */

    if (_scalars_lqc_) {


      if ((ppt->has_pk_lqc == _TRUE_)) {
   ppt->has_source_phi_scf = _TRUE_;
      }

       
      if ((ppt->has_bispectrum_lqc == _TRUE_)) {
   ppt->has_source_bispectrum_s_equi = _TRUE_;
      }

       

      index_type = index_type_common;

      class_define_index(ppt->index_tp_phi_scf,
                         ppt->has_source_phi_scf,
                         index_type,1);


      class_define_index(ppt->index_tp_pk_s_1,
                         ppt->has_source_phi_scf,
                         index_type,1);


      class_define_index(ppt->index_tp_pk_s_2,
                         ppt->has_source_phi_scf,
                         index_type,1);


      class_define_index(ppt->index_tp_pk_s_3,
                         ppt->has_source_phi_scf,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_s_123_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_123_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
     class_define_index(ppt->index_tp_bispectrum_s_prefactor_123_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);





      class_define_index(ppt->index_tp_bispectrum_s_1d23_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_1d23_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);



      class_define_index(ppt->index_tp_bispectrum_s_12d3_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_12d3_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);




      class_define_index(ppt->index_tp_bispectrum_s_123d_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_123d_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_123d_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);






      class_define_index(ppt->index_tp_bispectrum_s_1d23d_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d23d_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);



      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_1d2d3_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);




      class_define_index(ppt->index_tp_bispectrum_s_12d3d_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);     
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      class_define_index(ppt->index_tp_bispectrum_s_12d3d_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);










      
      class_define_index(ppt->index_tp_bispectrum_s_tau_end_111,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      

      class_define_index(ppt->index_tp_bispectrum_s_tau_end_110,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);

      class_define_index(ppt->index_tp_bispectrum_s_tau_end_011,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      
      class_define_index(ppt->index_tp_bispectrum_s_tau_end_101,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
            

      class_define_index(ppt->index_tp_bispectrum_s_tau_end_001,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      
      class_define_index(ppt->index_tp_bispectrum_s_tau_end_010,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      
      class_define_index(ppt->index_tp_bispectrum_s_tau_end_100,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      

      class_define_index(ppt->index_tp_bispectrum_s_tau_end_000,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      

      class_define_index(ppt->index_tp_bispectrum_s_tau_end_field_redef,
                         ppt->has_source_bispectrum_s_equi,
                         index_type,1);
      

      
      ppt->tp_size[index_md] = index_type;


    }



    if (_tensors_lqc_) {

      if ((ppt->has_pk_lqc == _TRUE_)) {
   ppt->has_source_gw = _TRUE_;
      }


       
      if ((ppt->has_bispectrum_lqc == _TRUE_)) {
   ppt->has_source_bispectrum_t_equi = _TRUE_;
      }


      
      index_type = index_type_common;
      
      class_define_index(ppt->index_tp_gw,
                         ppt->has_source_gw,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_t_equi_111,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_t_equi_110,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);

  
      class_define_index(ppt->index_tp_bispectrum_t_equi_100,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_t_equi_000,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);
            
      class_define_index(ppt->index_tp_bispectrum_t_equi_tau_end_111,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_t_equi_tau_end_110,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);

  
      class_define_index(ppt->index_tp_bispectrum_t_equi_tau_end_100,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);


      class_define_index(ppt->index_tp_bispectrum_t_equi_tau_end_000,
                         ppt->has_source_bispectrum_t_equi,
                         index_type,1);
      


      
      ppt->tp_size[index_md] = index_type;

  
    }

      printf("Defined indices for mode %d... number of indices = %d\n", index_md, ppt->tp_size[index_md]);

   class_alloc(ppt->fields_at_tau_end_lqc[index_md],
               ppt->tp_size[index_md]*sizeof(double *),
                ppt->error_message);


    
  }

  return _SUCCESS_;

}

/**
 * Define time sampling for source functions.
 *
 * For each type, compute the list of values of tau at which sources
 * will be sampled.  Knowing the number of tau values, allocate all
 * arrays of source functions.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */

int perturb_lqc_timesampling_for_sources(
                                     struct precision * ppr,
                                     struct background_lqc * pba,
                                     struct perturbs_lqc * ppt
                                     ) {

  int index_md;
  int index_type;
  ppt->tau_size = 1;
  class_alloc(ppt->tau_sampling,ppt->tau_size*sizeof(double),ppt->error_message);
//  ppt->tau_sampling[0] = 0.;
  ppt->tau_sampling[0] = ppt->tau_end;

 printf("tau_end = %e\n",ppt->tau_end); 
  for (index_md = 0; index_md < ppt->md_size; index_md++) {
      for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++) {


        class_alloc(ppt->fields_at_tau_end_lqc[index_md][index_type],
                    ppt->k_size[index_md]*sizeof(double),
                    ppt->error_message);
        

      }
    }
  

  return _SUCCESS_;
}

/**
 * Define the number of comoving wavenumbers using the information
 * passed in the precision structure.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to perturbation structure
 * @return the error status
 */

int perturb_lqc_get_k_list(
                        struct precision * ppr,
                        struct background_lqc * pba,
                        struct perturbs_lqc * ppt
                        ) {
  int index_k, index_k_output, index_mode;
  double k,k_min=0.,k_rec,step,tau1;
  double * k_max_cmb;
  double * k_max_cl;
  double k_max=0.;
  double scale2;
  double *tmp_k_list;
  int newk_size, index_newk, add_k_output_value;

  /** Summary: */

  class_test(ppr->k_step_transition == 0.,
             ppt->error_message,
             "stop to avoid division by zero");

  /** - allocate arrays related to k list for each mode */

  class_alloc(ppt->k_size,
              ppt->md_size*sizeof(int),
              ppt->error_message);

  class_alloc(ppt->k,
              ppt->md_size*sizeof(double*),
              ppt->error_message);
  
  /** - scalar modes */

  if (ppt->has_scalars == _TRUE_) {


    k_min =  ppt->k_min_for_pk;
    k_max = ppt->k_max_for_pk;
   


    /* if K>0, the transfer function will be calculated for discrete
       integer values of nu=3,4,5,... where nu=sqrt(k2+(1+m)K) and
       m=0,1,2 for scalars/vectors/tensors. However we are free to
       define in the perturbation module some arbitrary values of k:
       later on, the transfer module will interpolate at values of k
       corresponding exactly to integer values of nu. Hence, apart
       from the value of k_min and the step size in the vicinity of
       k_min, we define exactly the same sampling in the three cases
       K=0, K<0, K>0 */

    /* allocate array with, for the moment, the largest possible size */
    class_alloc(ppt->k[ppt->index_md_scalars],
                (int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)*sizeof(double),ppt->error_message);

    index_k=0;
    k = k_min;
    ppt->k[ppt->index_md_scalars][index_k] = k;
    index_k++;

    /* values until k_max[ppt->index_md_scalars] */

    while (k < k_max) {
    //while (index_k < ppt->k_per_decade_for_pk_lqc) {

      k *= pow(10.,1./(ppt->k_per_decade_for_pk_lqc));
                     
      // k += (k_max-k_min)/(ppt->k_per_decade_for_pk_lqc-1.);
      ppt->k[ppt->index_md_scalars][index_k] = k;

      index_k++;


    }


    ppt->k_size[ppt->index_md_scalars] = index_k;

    class_realloc(ppt->k[ppt->index_md_scalars],
                  ppt->k[ppt->index_md_scalars],
                  ppt->k_size[ppt->index_md_scalars]*sizeof(double),
                  ppt->error_message);
  }

  /** - tensor modes */

  if (ppt->has_tensors == _TRUE_) {

  /* allocate array with, for the moment, the largest possible size */
    class_alloc(ppt->k[ppt->index_md_tensors],
                ((int)(ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
    *sizeof(double),ppt->error_message);


    /* allocate array with, for the moment, the largest possible size */
    /*
    class_alloc(ppt->table_x2[ppt->index_md_tensors],
                (int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
                *(int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
                *sizeof(double),ppt->error_message);
    class_alloc(ppt->table_x3[ppt->index_md_tensors],
                (int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
                                *(int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
*sizeof(double),ppt->error_message);
   class_alloc(ppt->table_fNL[ppt->index_md_tensors],
                (int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
                               *(int)((ppt->k_per_decade_for_pk_lqc*log(k_max/k_min)/log(10.))+3)
*sizeof(double),ppt->error_message);

*/

    
    /* first value */

    index_k=0;
    k = k_min;
    ppt->k[ppt->index_md_tensors][index_k] = k;
    index_k++;

    /* values until k_max[ppt->index_md_scalars] */

    while (k < k_max) {

      k *= pow(10.,1./(ppt->k_per_decade_for_pk_lqc));
                      

      ppt->k[ppt->index_md_tensors][index_k] = k;
      index_k++;
    }


    //ppt->x_size[ppt->index_md_tensors] = 0;

    ppt->k_size[ppt->index_md_tensors] = index_k;

    class_realloc(ppt->k[ppt->index_md_tensors],
                  ppt->k[ppt->index_md_tensors],
                  ppt->k_size[ppt->index_md_tensors]*sizeof(double),
                  ppt->error_message);
  }


  /** - If user asked for k_output_values, add those to all k lists: */
  if (ppt->k_output_values_num>0){
    /* Allocate storage */
    class_alloc(ppt->index_k_output_values,sizeof(double)*ppt->md_size*ppt->k_output_values_num,ppt->error_message);

    /** - --> Find indices in ppt->k[index_md] corresponding to 'k_output_values'.
        We are assuming that ppt->k is sorted and growing, and we have made sure
        that ppt->k_output_values is also sorted and growing.*/
    for (index_mode=0; index_mode<ppt->md_size; index_mode++){

      newk_size = ppt->k_size[index_mode]+ppt->k_output_values_num;

      class_alloc(tmp_k_list,sizeof(double)*newk_size,ppt->error_message);

      index_k=0;
      index_k_output=0;
      for (index_newk=0; index_newk<newk_size; index_newk++){
        /** - --> Decide if we should add k_output_value now. This has to be this complicated, since we
            can only compare the k-values when both indices are in range.*/
        if (index_k >= ppt->k_size[index_mode])
          add_k_output_value = _TRUE_;
        else if (index_k_output >= ppt->k_output_values_num)
          add_k_output_value = _FALSE_;
        else if (ppt->k_output_values[index_k_output] < ppt->k[index_mode][index_k])
          add_k_output_value = _TRUE_;
        else
          add_k_output_value = _FALSE_;

        if (add_k_output_value == _TRUE_){
          tmp_k_list[index_newk] = ppt->k_output_values[index_k_output];
          ppt->index_k_output_values[index_mode*ppt->k_output_values_num+index_k_output]=index_newk;
          index_k_output++;
        }
        else{
          tmp_k_list[index_newk] = ppt->k[index_mode][index_k];
          index_k++;
        }
      }

      free(ppt->k[index_mode]);
      ppt->k[index_mode] = tmp_k_list;
      ppt->k_size[index_mode] = newk_size;

    }
  }

  /* For testing, can be useful to print the k list in a file:

  FILE * out=fopen("output/k","w");

  for (index_k=0; index_k < ppt->k_size[0]; index_k++) {

    fprintf(out,"%e\n",ppt->k[0][index_k],pba->K);

  }
     fclose(out);
  */

  /** - finally, find the global k_min and k_max for the ensemble of all modes 9scalars, vectors, tensors) */

  ppt->k_min = _HUGE_;
  ppt->k_max = 0.;
  if (ppt->has_scalars == _TRUE_) {
    ppt->k_min = MIN(ppt->k_min,ppt->k[ppt->index_md_scalars][0]); /* first value, inferred from perturbations_lqc structure */
    ppt->k_max = MAX(ppt->k_max,ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1]); /* last value, inferred from perturbations_lqc structure */
  }

  if (ppt->has_tensors == _TRUE_) {
    ppt->k_min = MIN(ppt->k_min,ppt->k[ppt->index_md_tensors][0]); /* first value, inferred from perturbations_lqc structure */
    ppt->k_max = MAX(ppt->k_max,ppt->k[ppt->index_md_tensors][ppt->k_size[ppt->index_md_tensors]-1]); /* last value, inferred from perturbations_lqc structure */
  }

  return _SUCCESS_;

}

/**
 * Initialize a perturb_lqc_workspace structure. All fields are allocated
 * here, with the exception of the perturb_lqc_vector '-->pv' field, which
 * is allocated separately in perturb_lqc_vector_init. We allocate one
 * such perturb_lqc_workspace structure per thread and per mode
 * (scalar/../tensor). Then, for each thread, all initial conditions
 * and wavenumbers will use the same workspace.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input/Output: pointer to perturb_lqc_workspace structure which fields are allocated or filled here
 * @return the error status
 */

int perturb_lqc_workspace_init(
                           struct precision * ppr,
                           struct background_lqc * pba,
                           struct perturbs_lqc * ppt,
                           int index_md,
                           struct perturb_lqc_workspace * ppw
                           ) {

  /** Summary: */

  /** - define local variables */

  int index_mt=0;
 
  /** - define indices of metric perturbations_lqc obeying constraint
      equations (this can be done once and for all, because the
      vector of metric perturbations_lqc is the same whatever the
      approximation scheme, unlike the vector of quantities to
      be integrated, which is allocated separately in
      perturb_lqc_vector_init) */

  if (_scalars_lqc_) {

        class_define_index(ppw->index_mt_phi_scf_prime_prime_re,_TRUE_,index_mt,1); 
        class_define_index(ppw->index_mt_phi_scf_prime_prime_im,_TRUE_,index_mt,1);

        class_define_index(ppw->index_mt_phi_scf_prime_prime_re_1,_TRUE_,index_mt,1); 
        class_define_index(ppw->index_mt_phi_scf_prime_prime_im_1,_TRUE_,index_mt,1);

        class_define_index(ppw->index_mt_phi_scf_prime_prime_re_2,_TRUE_,index_mt,1); 
        class_define_index(ppw->index_mt_phi_scf_prime_prime_im_2,_TRUE_,index_mt,1);

        class_define_index(ppw->index_mt_phi_scf_prime_prime_re_3,_TRUE_,index_mt,1); 
        class_define_index(ppw->index_mt_phi_scf_prime_prime_im_3,_TRUE_,index_mt,1);

      class_define_index(ppw->index_mt_phi_scf_prime_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_prime_im,_TRUE_,index_mt,1);


      class_define_index(ppw->index_mt_phi_scf_prime_re_1,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_prime_im_1,_TRUE_,index_mt,1);


      class_define_index(ppw->index_mt_phi_scf_prime_re_2,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_prime_im_2,_TRUE_,index_mt,1);


      class_define_index(ppw->index_mt_phi_scf_prime_re_3,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_prime_im_3,_TRUE_,index_mt,1);



      
      class_define_index(ppw->index_mt_phi_scf_abs,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_pk_s,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_pk_s_1,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_pk_s_2,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_pk_s_3,_TRUE_,index_mt,1); 


      class_define_index(ppw->index_mt_phi_scf_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_im,_TRUE_,index_mt,1);

      
      class_define_index(ppw->index_mt_phi_scf_re_1,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_im_1,_TRUE_,index_mt,1);

      
      class_define_index(ppw->index_mt_phi_scf_re_2,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_im_2,_TRUE_,index_mt,1);

      
      class_define_index(ppw->index_mt_phi_scf_re_3,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_scf_im_3,_TRUE_,index_mt,1);

      
      class_define_index(ppw->index_mt_R_s_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_R_s_im,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_R_s_prime_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_R_s_prime_im,_TRUE_,index_mt,1);

      class_define_index(ppw->index_mt_v_s_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_s_im,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_s_prime_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_s_prime_im,_TRUE_,index_mt,1); 
      
      //1
      class_define_index(ppw->index_mt_phi_prefactor_123_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123_111,_TRUE_,index_mt,1);


      //2
      class_define_index(ppw->index_mt_phi_1d23_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23_111,_TRUE_,index_mt,1);

      //3
      class_define_index(ppw->index_mt_phi_12d3_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3_111,_TRUE_,index_mt,1);

      //4
      class_define_index(ppw->index_mt_phi_123d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_123d_111,_TRUE_,index_mt,1);

      //5
      class_define_index(ppw->index_mt_phi_1d2d3_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d2d3_111,_TRUE_,index_mt,1);

          
      //6
      class_define_index(ppw->index_mt_phi_1d23d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_1d23d_111,_TRUE_,index_mt,1);

      //7
      class_define_index(ppw->index_mt_phi_12d3d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_phi_12d3d_111,_TRUE_,index_mt,1);

   
      //1
      class_define_index(ppw->index_mt_bispectrum_s_prefactor_123_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123_111,_TRUE_,index_mt,1);


      //2
      class_define_index(ppw->index_mt_bispectrum_s_1d23_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23_111,_TRUE_,index_mt,1);


      //3
      class_define_index(ppw->index_mt_bispectrum_s_12d3_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3_111,_TRUE_,index_mt,1);


      //4
      class_define_index(ppw->index_mt_bispectrum_s_123d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_123d_111,_TRUE_,index_mt,1);

      //5
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d23d_111,_TRUE_,index_mt,1);

      //6
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_1d2d3_111,_TRUE_,index_mt,1);


      //7
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_010,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_001,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_101,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_011,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_s_12d3d_111,_TRUE_,index_mt,1);


     class_define_index(ppw->index_mt_bispectrum_s_field_redef,_TRUE_,index_mt,1);

      
    }


  if (_tensors_lqc_) {

    class_define_index(ppw->index_mt_gwdot_re,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gwdot_im,_TRUE_,index_mt,1);
    
    class_define_index(ppw->index_mt_gw_prime_prime_re,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_prime_prime_im,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_abs,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_pk_t,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_re,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_im,_TRUE_,index_mt,1);


      class_define_index(ppw->index_mt_v_t_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_t_im,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_t_prime_re,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_v_t_prime_im,_TRUE_,index_mt,1); 
      



    class_define_index(ppw->index_mt_gw_000,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_100,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_110,_TRUE_,index_mt,1);
    class_define_index(ppw->index_mt_gw_111,_TRUE_,index_mt,1);

      class_define_index(ppw->index_mt_bispectrum_t_000,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_t_100,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_t_110,_TRUE_,index_mt,1); 
      class_define_index(ppw->index_mt_bispectrum_t_111,_TRUE_,index_mt,1);

    
  }

  ppw->mt_size = index_mt;

  /** - allocate some workspace in which we will store temporarily the
      values of background, thermodynamics, metric and source
      quantities at a given time */

  class_alloc(ppw->pvecback,pba->bg_size_normal*sizeof(double),ppt->error_message);
  class_alloc(ppw->pvecmetric,ppw->mt_size*sizeof(double),ppt->error_message);

 

 
  return _SUCCESS_;
}

/**
 * Free the perturb_lqc_workspace structure (with the exception of the
 * perturb_lqc_vector '-->pv' field, which is freed separately in
 * perturb_lqc_vector_free).
 *
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input: pointer to perturb_lqc_workspace structure to be freed
 * @return the error status
 */

int perturb_lqc_workspace_free (
                            struct perturbs_lqc * ppt,
                            int index_md,
                            struct perturb_lqc_workspace * ppw
                            ) {
  free(ppw->pvecback);
  free(ppw->pvecmetric);

  //free(ppw->approx);


  free(ppw);

  return _SUCCESS_;
}

/**
 * Solve the perturbation evolution for a given mode, initial
 * condition and wavenumber, and compute the corresponding source
 * functions.
 *
 * For a given mode, initial condition and wavenumber, this function
 * finds the time ranges over which the perturbations_lqc can be described
 * within a given approximation. For each such range, it initializes
 * (or redistributes) perturbations_lqc using perturb_lqc_vector_init(), and
 * integrates over time. Whenever a "source sampling time" is passed,
 * the source terms are computed and stored in the source table using
 * perturb_lqc_sources().
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input/Output: pointer to the perturbation structure (output source functions S(k,tau) written here)
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param index_k    Input: index of wavenumber
 * @param ppw        Input: pointer to perturb_lqc_workspace structure containing index values and workspaces
 * @return the error status
 */

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
                  ) {

  /** Summary: */

  /** - define local variables */

  /* contains all fixed parameters, indices and workspaces used by the perturb_lqc_derivs function */
  struct perturb_lqc_parameters_and_workspace ppaw;

  /* conformal time */
  double tau,tau_lower,tau_upper,tau_mid;

  /* index running over time */
  int index_tau;

  /* number of values in the tau_sampling array that should be considered for a given mode */
  int tau_actual_size;

  /* running index over types (temperature, etc) */
  int index_type;

  /* Fourier mode */
  double k;
  double k_1;
  double k_2;
  double k_3;

  /* function pointer to ODE evolver and names of possible evolvers */

  extern int evolver_rk();
  extern int evolver_ndf15();
  int (*generic_evolver)();


  /* Related to the perturbation output */
  int (*perhaps_print_variables)();
  int index_ikout;
  int index_k;
  /** - initialize indices relevant for back/thermo tables search */
  ppw->last_index_back=0;
  ppw->inter_mode = pba->inter_normal;

  /** - get wavenumber value */
  k_1 = ppt->k[index_md][index_k_1];
  k_2 = ppt->lambda_k_2*ppt->k[index_md][index_k_2];
  //k_2 = 3.e0*pba->k_star;
  //k_3 = 3.e0*pba->k_star;
  // double theta = 5.*_PI_/5.;
  //k_3 = sqrt(k_1*k_1+k_2*k_2+2.*k_1*k_2*cos(theta));
  k_3 = ppt->lambda_k_3*ppt->k[index_md][index_k_3];
//  k_3 = 1.e-4*pba->k_star;//V. Sreenath

  

  
  index_k = MIN(index_k_1,index_k_2);
  index_k = MIN(index_k,index_k_3);

  k = ppt->k[index_md][index_k];



  

  
  class_test(k == 0.,
             ppt->error_message,
             "stop to avoid division by zero");


  tau_actual_size = ppt->tau_size;


  ppw->inter_mode = pba->inter_normal;

  /** - fill the structure containing all fixed parameters, indices
      and workspaces needed by perturb_lqc_derivs */

  ppaw.ppr = ppr;
  ppaw.pba = pba;
  ppaw.ppt = ppt;
  ppaw.index_md = index_md;
  ppaw.index_k = index_k_1;
  ppaw.k = k_1;
  ppaw.k_2 = k_2;
  ppaw.k_3 = k_3;
  ppaw.ppw = ppw;
  ppaw.ppw->inter_mode = pba->inter_closeby;
  ppaw.ppw->last_index_back = 0;

  /** - check whether we need to print perturbations_lqc to a file for this wavenumber */

  perhaps_print_variables = NULL;
  ppw->index_ikout = -1;
  for (index_ikout=0; index_ikout<ppt->k_output_values_num; index_ikout++){
    if (ppt->index_k_output_values[index_md*ppt->k_output_values_num+index_ikout] == index_k){
      ppw->index_ikout = index_ikout;
      perhaps_print_variables = perturb_lqc_print_variables;

    }
  }

    /** - --> (c) define the vector of perturbations_lqc to be integrated
        over. If the current interval starts from the initial time
        tau_ini, fill the vector with initial conditions for each
        mode. If it starts from an approximation switching point,
        redistribute correctly the perturbations_lqc from the previous to
        the new vector of perturbations_lqc. */

  //For each k, the evolution starts when k=integration_depth_ini*aH
  //Fisrt, we compute the value of tau_ini_for_k

 //find tau ini

  ppw->tau_ini_for_k = ppt->tau_ini_for_k;
  double tau_max=pba->tau_end_inflation;
    //pba->tau_start_inflation
    //            +pba->N_end
    //            -pba->N_start_inflation;
  double tau_min=pba->tau_start_inflation;
  double tau_half;
  double aH_sample=_HUGE_;
  double error_aH = .01;
  int index_back;
    double * pvecback2;
     class_alloc(pvecback2,pba->bg_size_normal*sizeof(double),pba->error_message);
            
      aH_sample =
      pba->a_bounce*sqrt(pba->rho_bounce/3.);


      if (ppw->tau_ini_for_k<0.){
      
    if (ppt->k[index_md][index_k]>ppt->integration_depth_ini*aH_sample){
           tau_max=ppt->tau_end;
           tau_min = ppw->tau_ini_for_k;           
         }
   else{
     // ppt->tau_ini_for_k=ppt->tau_ini;
           tau_min = ppt->tau_ini;
           tau_max=0.;
       
        }

    tau_min = ppw->tau_ini_for_k;
    //printf("Start tau_ini_for_k = %e \n",ppt->tau_ini_for_k);
    // printf("Start tau_min = %e \n",tau_min);
    //printf("Start tau_max = %e \n",tau_max);
 
  class_call(background_lqc_at_tau(pba,
                                 tau_min,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
    aH_sample =
      pvecback2[pba->index_bg_a]
      *sqrt(pvecback2[pba->index_bg_rho_scf]/3.);


    
    if (ppt->k[index_md][index_k]>ppt->integration_depth_ini*aH_sample){

      //printf("aH-k/K_ini = %e \n",aH_sample-ppt->k[index_md][index_k]
      //        /ppt->integration_depth_ini);

      
 while (fabs(aH_sample-ppt->k[index_md][index_k]
              /ppt->integration_depth_ini)>
        error_aH*ppt->k[index_md][index_k]
              /ppt->integration_depth_ini){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
   aH_sample =
      pvecback2[pba->index_bg_a]
      *sqrt(pvecback2[pba->index_bg_rho_scf]/3.);

   if (aH_sample>ppt->k[index_md][index_k]/ppt->integration_depth_ini)
    tau_max = tau_half;
  else tau_min = tau_half;
  ppw->tau_ini_for_k = tau_half;

  //printf("in while= %e \n",1.);

 }//end while

 //printf("Start perturbations tau_ini = %e \n",ppt->tau_ini_for_k);
 //    printf("error_aH = %e \n",error_aH);
 //    printf("aH-k/K_ini = %e \n",aH_sample-ppt->k[index_md][index_k]
 //            /ppt->integration_depth_ini);

 
    }//end if k>K_ini*aH_ini
    
    else { //case k<K_ini*aH_ini, so one has to go back in time to find
      // the time at which k=K_ini*aH_ini

      //printf("Start perturbations tau_ini = %e \n",ppt->tau_ini_for_k);

      tau_max = ppw->tau_ini_for_k;
      //tau_min = pba->tau_table[0];
      tau_min = pba->tau_end_deflation;
         
  class_call(background_lqc_at_tau(pba,
                                 tau_min,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
    aH_sample =
      pvecback2[pba->index_bg_a]
      *sqrt(pvecback2[pba->index_bg_rho_scf]/3.);


    if (ppt->k[index_md][index_k]>ppt->integration_depth_ini*aH_sample){

 while (fabs(aH_sample-ppt->k[index_md][index_k]
              /ppt->integration_depth_ini)>error_aH*ppt->k[index_md][index_k]
              /ppt->integration_depth_ini){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
   aH_sample =
      pvecback2[pba->index_bg_a]
      *sqrt(pvecback2[pba->index_bg_rho_scf]/3.);

   if (aH_sample>ppt->k[index_md][index_k]/ppt->integration_depth_ini)
    tau_max = tau_half;
  else tau_min = tau_half;
  ppw->tau_ini_for_k = tau_half;

 }//end while
    }//end if k>K_ini*aH_ini

    else {printf("Can not find a time at which  k>K_ini*aH_ini.\n");
      ppw->tau_ini_for_k = pba->tau_end_deflation;
     
    }
}
    ppw->tau_ini_for_k = pba->tau_end_deflation;

         printf("k/k_star=%e, tau_ini = %e\n",
             ppt->k[index_md][index_k]/pba->k_star,
             ppw->tau_ini_for_k);
   
      }//end if tau_ini<0
//    ppw->tau_ini_for_k = -5.e4;//1.e-5/sqrt(8.*_PI_);//V. Sreenath

      printf("k/k_star=%e, tau_ini = %e, tau_end=%e\n",
             ppt->k[index_md][index_k]/pba->k_star,
               ppw->tau_ini_for_k,
                              ppt->tau_end);
  

      
  int index_interval = 0;
  double * interval_limit;
  int interval_number=0;


  //if (ppt->tau_ini_for_k<2.e0*pba->tau_end_remote_past)
  //     ppt->tau_ini_for_k = 2.e0*pba->tau_end_remote_past;

  /*
if (index_md==ppt->index_md_scalars
    && ppw->tau_ini_for_k<pba->tau_end_remote_past){

 tau_max = pba->tau_table[0];
 tau_min = pba->tau_end_remote_past;
 class_call(background_lqc_at_tau(pba,
                                 tau_min,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);

  double effective_potential_over_omega2 =
    fabs(pvecback2[pba->index_bg_effective_potential]
    /(k*k
    +pvecback2[pba->index_bg_a]
    *pvecback2[pba->index_bg_a]
    *pba->m_scf_lqc
      *pba->m_scf_lqc));

  //printf("Mode k=%e, effective_potential_over_omega2 = %e \n",k,effective_potential_over_omega2);

  double scalar_limit = 0.001;
    if (effective_potential_over_omega2>scalar_limit){
      
      while (fabs(effective_potential_over_omega2-scalar_limit)/scalar_limit>0.001){

    tau_half = (tau_max+tau_min)/2;
    
  class_call(background_lqc_at_tau(pba,
                                 tau_half,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &index_back,
                                 pvecback2),
               pba->error_message,
               pba->error_message);
  
  effective_potential_over_omega2 = fabs(
    pvecback2[pba->index_bg_effective_potential]
    /(k*k
    +pvecback2[pba->index_bg_a]
    *pvecback2[pba->index_bg_a]
    *pba->m_scf_lqc
      *pba->m_scf_lqc));
  //printf("whilMode k=%e, effective_potential_over_omega2 = %e \n",k,effective_potential_over_omega2);

   if (effective_potential_over_omega2<scalar_limit)
    tau_max = tau_half;
  else tau_min = tau_half;
  ppw->tau_ini_for_k = tau_half;

 }//end while

      //printf("Mode k=%e, at t_ini, effective_potential_over_omega2 = %e \n",k,effective_potential_over_omega2);
      
    }//end if k>K_ini*aH_ini



  if (ppw->tau_ini_for_k<2.e0*pba->tau_end_remote_past)
       ppw->tau_ini_for_k = 10.e0*pba->tau_end_remote_past;

    
 }

  */

  // ppw->tau_ini_for_k = ppt->tau_ini;
  
if (ppw->tau_ini_for_k<pba->tau_end_deflation){
  interval_number=4;
 }
 else if (ppw->tau_ini_for_k<0.){
   interval_number=3;
  }
 else if (ppw->tau_ini_for_k<pba->tau_start_inflation){
   interval_number=2;
  }
 else{
   interval_number=1;
       }
//printf("Mode k=%e, interval_number = %d \n",k,interval_number);

  ppw->num_approx_intervals = interval_number;
 
  class_alloc(interval_limit,(interval_number+1)*sizeof(double),ppt->error_message);



  int * interval_approx;
 class_alloc(interval_approx,interval_number,ppt->error_message);

  
if (ppw->tau_ini_for_k<pba->tau_end_deflation){
  interval_limit[index_interval]=ppw->tau_ini_for_k;
  interval_approx[index_interval] = (int)efolds;
  index_interval += 1;
  interval_limit[index_interval]=pba->tau_end_deflation;
  interval_approx[index_interval] = (int)cosmic_time;
  index_interval += 1;
  interval_limit[index_interval]=0.;
  interval_approx[index_interval] = (int)cosmic_time;
  index_interval += 1;
  interval_limit[index_interval]=pba->tau_start_inflation;
  interval_approx[index_interval] = (int)efolds;
  index_interval += 1;
  interval_limit[index_interval]=ppt->tau_end;  
   }
else if (ppw->tau_ini_for_k< 0.){
  interval_limit[index_interval]=ppw->tau_ini_for_k;
  interval_approx[index_interval] = (int)cosmic_time;
  index_interval += 1;
  interval_limit[index_interval]=0.;
  interval_approx[index_interval] = (int)cosmic_time;
  index_interval += 1;
  interval_limit[index_interval]=pba->tau_start_inflation;
  interval_approx[index_interval] = (int)efolds;
  index_interval += 1;
  interval_limit[index_interval]=ppt->tau_end;  
   }
 else if (ppw->tau_ini_for_k<pba->tau_start_inflation){
  interval_limit[index_interval]=ppw->tau_ini_for_k;
  interval_approx[index_interval] = (int)cosmic_time;
  index_interval += 1;
  interval_limit[index_interval]=pba->tau_start_inflation;
  interval_approx[index_interval] = (int)efolds;
  index_interval += 1;
  interval_limit[index_interval]=ppt->tau_end;  
   }

 else {
  interval_limit[index_interval]=ppw->tau_ini_for_k;
  interval_approx[index_interval] = (int)efolds;
  index_interval += 1;
  interval_limit[index_interval]=ppt->tau_end;
   }



  
      
  tau_lower = ppw->tau_ini_for_k;
  tau_upper = ppt->tau_end;

  //printf("Mode k=%e, t_end = %e\n",k,tau_upper);


  
  double tau_ini = tau_lower;
  double tau_end = tau_upper;

  


  /** - loop over intervals over which approximation scheme is uniform. For each interval: */
  //index_interval = 0;
  //interval_number=1;

  //interval_limit[index_interval] = tau_ini;
  //interval_limit[index_interval+1] = tau_end;
  
  for (index_interval=0; index_interval<interval_number; index_interval++) {


    ppw->approx=interval_approx[index_interval];
    ppw->approx_interval = index_interval;


if(ppw->approx == (int)cosmic_time){

  
    tau_lower = interval_limit[index_interval];
    tau_upper = interval_limit[index_interval+1];
/*  if(tau_lower < 0.){
	 ppaw.pba->future_branch = 0;
	pba->future_branch = 0;
//	ppt->tau_sampling[0] = 0.;
 	}
  else {
	ppaw.pba->future_branch = 1;
	pba->future_branch = 1;
//	ppt->tau_sampling[0] = tau_end;
//	ppr->perturb_lqc_integration_stepsize = 1.e-2;
	}*/
//  printf("Mode k=%e, cosmic time starts at t = %e, ends at t = %e... future branch = %d\n",
//         k,
//          interval_limit[index_interval],
//         interval_limit[index_interval+1], pba->future_branch);
//printf("future branch = %d\n",pba->future_branch);
     }
else{

//    if(interval_limit[index_interval]<0.) pba->future_branch = 0;
//    else 
/*	ppaw.pba->future_branch = 1;
	pba->future_branch = 1;*/
    double N_min,N_max;
    double tau_of_N_min,tau_of_N_max;

    if(ppw->num_approx_intervals > 3 && ppw->approx_interval < 1){
    class_call(background_lqc_N_of_tau_pb(pba,
                               interval_limit[index_interval],
                               &N_min),
                   pba->error_message,
                   ppt->error_message);
    
    class_call(background_lqc_N_of_tau_pb(pba,
                               interval_limit[index_interval+1],
                               &N_max),
                   pba->error_message,
                   ppt->error_message);
  
//    printf("Mode k=%e, efolds starts at N = %e, ends at N = %e....future branch = %d\n",
//           k,
//          N_min,
//          N_max, pba->future_branch);

  class_call(background_lqc_tau_of_N_pb(pba,
                               N_max,
                               &tau_of_N_max),
                   pba->error_message,
                   ppt->error_message);
   class_call(background_lqc_tau_of_N_pb(pba,
                               N_min,
                               &tau_of_N_min),
                   pba->error_message,
                   ppt->error_message);
   }
   else{
    class_call(background_lqc_N_of_tau_fb(pba,
                               interval_limit[index_interval],
                               &N_min),
                   pba->error_message,
                   ppt->error_message);
    
    class_call(background_lqc_N_of_tau_fb(pba,
                               interval_limit[index_interval+1],
                               &N_max),
                   pba->error_message,
                   ppt->error_message);
  
//    printf("Mode k=%e, efolds starts at N = %e, ends at N = %e....future branch = %d\n",
//           k,
//          N_min,
//          N_max, pba->future_branch);

  class_call(background_lqc_tau_of_N_fb(pba,
                               N_max,
                               &tau_of_N_max),
                   pba->error_message,
                   ppt->error_message);
   class_call(background_lqc_tau_of_N_fb(pba,
                               N_min,
                               &tau_of_N_min),
                   pba->error_message,
                   ppt->error_message);


   }
   
    
/* printf("Mode k=%e, cosmic time starts at t = %e, ends at t = %e\n",
           k,
           tau_of_N_min,
           tau_of_N_max);
      
   printf("Mode k=%e, cosmic time starts at t = %e, ends at t = %e\n",
           k,
           interval_limit[index_interval],
           interval_limit[index_interval+1]);
   printf("future branch = %d\n",pba->future_branch);
*/ 
    tau_lower = N_min;
    tau_upper = N_max;
//    if(tau_upper < 0.) pba->future_branch = 0;
//    else pba->future_branch = 1;
}

//ppaw.pba = pba;
// printf("branch = %d",ppaw.pba->future_branch); 
    
    class_call(perturb_lqc_vector_init(ppr,
                                   pba,
                                   ppt,
                                   index_md,
                                   k_1,
                                   k_2,
                                   k_3,
                                   tau_lower,
                                   ppw,
                                   index_interval),
               ppt->error_message,
               ppt->error_message);

    /** - --> (d) integrate the perturbations_lqc over the current interval. */

    if(ppr->evolver == rk){
      generic_evolver = evolver_rk;
    }
    else{
      generic_evolver = evolver_ndf15;
    }
//       generic_evolver = evolver_ndf15;
   //printf("H=%e\n",tau_lower);
    //printf("H=%e\n",ppw->pvecback[pba->index_bg_H]);

    /*
    class_call(generic_evolver(perturb_lqc_derivs,
                               tau_lower,
                               tau_upper,
                               ppw->pv->y,
                               ppw->pv->used_in_sources,
                               ppw->pv->pt_size,
                               &ppaw,
                               ppr->tol_perturb_lqc_integration,
                               ppr->smallest_allowed_variation,
                               perturb_lqc_timescale,
                               ppr->perturb_lqc_integration_stepsize,
                               ppt->tau_sampling,
                               tau_actual_size,
                               perturb_lqc_sources,
                               perhaps_print_variables,
                               ppt->error_message),
               ppt->error_message,
               ppt->error_message);
    */


    generic_evolver(perturb_lqc_derivs,
                               tau_lower,
                               tau_upper,
                               ppw->pv->y,
                               ppw->pv->used_in_sources,
                               ppw->pv->pt_size,
                               &ppaw,
                               ppr->tol_perturb_lqc_integration,
                               ppr->smallest_allowed_variation,
                               perturb_lqc_timescale,
                               ppr->perturb_lqc_integration_stepsize,
                               ppt->tau_sampling,
                               tau_actual_size,
                               perturb_lqc_sources,
                               perhaps_print_variables,
                    ppt->error_message);


 }


//printf("integration completed... \n");
    
    //FILL THE ARRAYS WITH VALUES OF THE FIELDS
    //AT THE END OF INTEGRATION
    
   if (index_md==ppt->index_md_scalars){


     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_phi_scf) =
       ppw->pvecmetric[ppw->index_mt_pk_s];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_pk_s_1) =
       ppw->pvecmetric[ppw->index_mt_pk_s_1];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_pk_s_2) =
       ppw->pvecmetric[ppw->index_mt_pk_s_2];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_pk_s_3) =
       ppw->pvecmetric[ppw->index_mt_pk_s_3];

	printf("Ps(%e) = %e = %e\n", ppt->k[index_md][index_k], ppw->pvecmetric[ppw->index_mt_pk_s_3], ppt->fields_at_tau_end_lqc[index_md][ppt->index_tp_pk_s_1][index_k]);

    if (ppt->has_bispectrum_lqc == _TRUE_) {

double pk_s_1 = ppw->pvecmetric[ppw->index_mt_pk_s_1];
double pk_s_2 = ppw->pvecmetric[ppw->index_mt_pk_s_2];
double pk_s_3 = ppw->pvecmetric[ppw->index_mt_pk_s_3];

      
    double Z_123_000,Z_123_100,Z_123_010,Z_123_001,Z_123_110,Z_123_101,Z_123_011,Z_123_111;
    double Z_1d23_000,Z_1d23_100,Z_1d23_010,Z_1d23_001,Z_1d23_110,Z_1d23_101,Z_1d23_011,Z_1d23_111;
    double Z_12d3_000,Z_12d3_100,Z_12d3_010,Z_12d3_001,Z_12d3_110,Z_12d3_101,Z_12d3_011,Z_12d3_111;
    double Z_123d_000,Z_123d_100,Z_123d_010,Z_123d_001,Z_123d_110,Z_123d_101,Z_123d_011,Z_123d_111;
    double Z_12d3d_000,Z_12d3d_100,Z_12d3d_010,Z_12d3d_001,Z_12d3d_110,Z_12d3d_101,Z_12d3d_011,Z_12d3d_111;
    double Z_1d2d3_000,Z_1d2d3_100,Z_1d2d3_010,Z_1d2d3_001,Z_1d2d3_110,Z_1d2d3_101,Z_1d2d3_011,Z_1d2d3_111;
    double Z_1d23d_000,Z_1d23d_100,Z_1d23d_010,Z_1d23d_001,Z_1d23d_110,Z_1d23d_101,Z_1d23d_011,Z_1d23d_111;


        double Y_000,Y_100,Y_010,Y_001,Y_110,Y_101,Y_011,Y_111;

 
 double Q_s_1_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_1];
 double Q_s_1d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1];
 
 double Q_s_2_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_2];
 double Q_s_2d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2];

 double Q_s_3_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_3];
 double Q_s_3d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3];
 
 double Q_s_1_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_1];
 double Q_s_1d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1];

 double Q_s_2_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_2];
 double Q_s_2d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2];

 double Q_s_3_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_3];
 double Q_s_3d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3];

 //H over phi_dot
 double a_over_z_3 =
     pow(ppw->pvecback[pba->index_bg_H]
         /ppw->pvecback[pba->index_bg_phi_prime_scf],3.);


Y_111 = a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);
 
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_111) =
       a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);


Y_110 = a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);

      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_110) =
       a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);

Y_101 = a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);
      
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_101) =
       a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);


Y_011 =  a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);
   
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_011) =
        a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);


Y_100 = a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);

     
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_100) =
        a_over_z_3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);

 Y_010 =       a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


   
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_010) =
       a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);

Y_001 = a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);
  
      
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_001) =
       a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);


Y_000 = a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


      
       _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_000) =
       a_over_z_3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


double field_redef = -ppw->pvecmetric[ppw->index_mt_bispectrum_s_field_redef];
 printf("Mode k=%e, field_redef=%e.\n",k,field_redef);

       

       _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_tau_end_field_redef) =
      -ppw->pvecmetric[ppw->index_mt_bispectrum_s_field_redef];

       
       //1

Z_123_111 = ppw->pvecmetric[ppw->index_mt_phi_123_111];
       
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_111) =
        ppw->pvecmetric[ppw->index_mt_phi_123_111];
      
Z_123_110 = ppw->pvecmetric[ppw->index_mt_phi_123_110];
      
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_110) =
        ppw->pvecmetric[ppw->index_mt_phi_123_110];
     
Z_123_101 = ppw->pvecmetric[ppw->index_mt_phi_123_101];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_101) =
        ppw->pvecmetric[ppw->index_mt_phi_123_101];

Z_123_011 = ppw->pvecmetric[ppw->index_mt_phi_123_011];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_011) =
        ppw->pvecmetric[ppw->index_mt_phi_123_011];

 Z_123_100 = ppw->pvecmetric[ppw->index_mt_phi_123_100];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_100) =
        ppw->pvecmetric[ppw->index_mt_phi_123_100];

Z_123_000 = ppw->pvecmetric[ppw->index_mt_phi_123_000];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_000) =
        ppw->pvecmetric[ppw->index_mt_phi_123_000];
 
double prefactor_123_000 = ppw->pvecmetric[ppw->index_mt_phi_prefactor_123_000];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_prefactor_123_000) =
        ppw->pvecmetric[ppw->index_mt_phi_prefactor_123_000];
 
Z_123_010 = ppw->pvecmetric[ppw->index_mt_phi_123_010];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_010) =
        ppw->pvecmetric[ppw->index_mt_phi_123_010];

Z_123_001 = ppw->pvecmetric[ppw->index_mt_phi_123_001];

 _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123_001) =
        ppw->pvecmetric[ppw->index_mt_phi_123_001];
 
     //2

Z_1d23_111 = ppw->pvecmetric[ppw->index_mt_phi_1d23_111];

 
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_111) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_111];

Z_1d23_110 = ppw->pvecmetric[ppw->index_mt_phi_1d23_110];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_110) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_110];

Z_1d23_101 = ppw->pvecmetric[ppw->index_mt_phi_1d23_101];

      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_101) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_101];

Z_1d23_011 = ppw->pvecmetric[ppw->index_mt_phi_1d23_011];

    _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_011) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_011];

Z_1d23_100 = ppw->pvecmetric[ppw->index_mt_phi_1d23_100];

      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_100) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_100];

Z_1d23_000 = ppw->pvecmetric[ppw->index_mt_phi_1d23_000];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_000) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_000];

Z_1d23_010 = ppw->pvecmetric[ppw->index_mt_phi_1d23_010];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_010) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_010];

Z_1d23_001 = ppw->pvecmetric[ppw->index_mt_phi_1d23_001];

     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23_001) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23_001];
 
     //3

Z_12d3_111 = ppw->pvecmetric[ppw->index_mt_phi_12d3_111];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_111) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_111];

Z_12d3_110 = ppw->pvecmetric[ppw->index_mt_phi_12d3_110];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_110) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_110];

Z_12d3_101 = ppw->pvecmetric[ppw->index_mt_phi_12d3_101];
     
       _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_101) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_101];

Z_12d3_011 = ppw->pvecmetric[ppw->index_mt_phi_12d3_011];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_011) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_011];

Z_12d3_100 = ppw->pvecmetric[ppw->index_mt_phi_12d3_100];
     
       _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_100) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_100];

Z_12d3_000 = ppw->pvecmetric[ppw->index_mt_phi_12d3_000];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_000) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_000];

Z_12d3_010 = ppw->pvecmetric[ppw->index_mt_phi_12d3_010];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_010) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_010];

Z_12d3_001 = ppw->pvecmetric[ppw->index_mt_phi_12d3_001];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3_001) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3_001];
 
     //4

Z_123d_111 = ppw->pvecmetric[ppw->index_mt_phi_123d_111];
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_111) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_111];

Z_123d_110 = ppw->pvecmetric[ppw->index_mt_phi_123d_110];
      
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_110) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_110];

Z_123d_101 = ppw->pvecmetric[ppw->index_mt_phi_123d_101];
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_101) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_101];

Z_123d_011 = ppw->pvecmetric[ppw->index_mt_phi_123d_011];
      
    _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_011) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_011];

Z_123d_100 = ppw->pvecmetric[ppw->index_mt_phi_123d_100];
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_100) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_100];

Z_123d_000 = ppw->pvecmetric[ppw->index_mt_phi_123d_000];
      
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_000) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_000];

Z_123d_010 = ppw->pvecmetric[ppw->index_mt_phi_123d_010];
      
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_010) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_010];

Z_123d_001 = ppw->pvecmetric[ppw->index_mt_phi_123d_001];
      
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_123d_001) =
        ppw->pvecmetric[ppw->index_mt_phi_123d_001];
 
     //5

Z_1d2d3_111 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_111];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_111) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_111];

Z_1d2d3_110 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_110];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_110) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_110];

Z_1d2d3_101 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_101];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_101) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_101];

Z_1d2d3_011 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_011];
     
    _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_011) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_011];

Z_1d2d3_100 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_100];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_100) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_100];

Z_1d2d3_000 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_000];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_000) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_000];

Z_1d2d3_010 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_010];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_010) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_010];

Z_1d2d3_001 = ppw->pvecmetric[ppw->index_mt_phi_1d2d3_001];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d2d3_001) =
        ppw->pvecmetric[ppw->index_mt_phi_1d2d3_001];
 
     //6

Z_1d23d_111 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_111];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_111) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_111];

Z_1d23d_110 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_110];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_110) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_110];

Z_1d23d_101 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_101];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_101) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_101];

Z_1d23d_011 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_011];
     
    _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_011) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_011];

Z_1d23d_100 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_100];
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_100) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_100];

Z_1d23d_000 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_000];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_000) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_000];

Z_1d23d_010 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_010];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_010) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_010];

Z_1d23d_001 = ppw->pvecmetric[ppw->index_mt_phi_1d23d_001];
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_1d23d_001) =
        ppw->pvecmetric[ppw->index_mt_phi_1d23d_001];
 

     //7

Z_12d3d_111 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_111];     
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_111) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_111];

Z_12d3d_110 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_110];     
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_110) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_110];

Z_12d3d_101 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_101];     
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_101) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_101];

Z_12d3d_011 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_011];     
     
    _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_011) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_011];

Z_12d3d_100 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_100];     
     
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_100) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_100];

Z_12d3d_000 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_000];     
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_000) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_000];

Z_12d3d_010 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_010];     
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_010) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_010];

Z_12d3d_001 = ppw->pvecmetric[ppw->index_mt_phi_12d3d_001];     
     
     _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_s_12d3d_001) =
        ppw->pvecmetric[ppw->index_mt_phi_12d3d_001];



      //1
        double integral_123 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_123_111
            -Z_123_001
            -Z_123_010
            -Z_123_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_123_000
            -Z_123_110
            -Z_123_101
            -Z_123_011);


      //2
        double integral_1d23 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d23_111
            -Z_1d23_001
            -Z_1d23_010
            -Z_1d23_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d23_000
            -Z_1d23_110
            -Z_1d23_101
            -Z_1d23_011);



        //3
        double integral_12d3 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_12d3_111
            -Z_12d3_001
            -Z_12d3_010
            -Z_12d3_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_12d3_000
            -Z_12d3_110
            -Z_12d3_101
            -Z_12d3_011);
        //4
        double integral_123d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_123d_111
            -Z_123d_001
            -Z_123d_010
            -Z_123d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_123d_000
            -Z_123d_110
            -Z_123d_101
            -Z_123d_011);
        //5
        double integral_1d2d3 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d2d3_111
            -Z_1d2d3_001
            -Z_1d2d3_010
            -Z_1d2d3_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d2d3_000
            -Z_1d2d3_110
            -Z_1d2d3_101
            -Z_1d2d3_011);
        
        //6
        double integral_1d23d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d23d_111
            -Z_1d23d_001
            -Z_1d23d_010
            -Z_1d23d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d23d_000
            -Z_1d23d_110
            -Z_1d23d_101
            -Z_1d23d_011);
        
        //7
        double integral_12d3d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_12d3d_111
            -Z_12d3d_001
            -Z_12d3d_010
            -Z_12d3d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_12d3d_000
            -Z_12d3d_110
            -Z_12d3d_101
            -Z_12d3d_011);

        //field_redef
        /*   double field_redef =
         ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_field_redef]
               [index_k];
        */

        
  double      integral =
         integral_123 //1
         +integral_1d23 //2
         +integral_12d3 //3
         +integral_123d //4
         +integral_1d2d3 //5
         +integral_1d23d //6
          +integral_12d3d //7
        +field_redef; //redef
 printf("Mode k=%e, integral = %e.\n",k,integral);


double fNL =
  -10./3.
  *pow(2.*_PI_,-4.)
  *pow(k_1*k_2*k_3,3.)
  *integral
  *pow(pow(k_1,3.)*pk_s_2*pk_s_3
       +pow(k_2,3.)*pk_s_1*pk_s_3
       +pow(k_3,3.)*pk_s_2*pk_s_1,-1.);
  
//printf("IN SOLVE fNL=%e\n",fNL);

        
//printf("\n");
             /*
printf("IN SOLVE x_2=%e",
       k_2/k_1);
            printf("\n");

printf("IN SOLVE x_3=%e",
       k_3/k_1);
            printf("\n");
             */
  FILE * pkfile;
   FileName file_name;
  sprintf(file_name,"%s_%s",pba->root,"bispectrum_s.dat");

  class_open(pkfile,file_name,"a",pba->error_message);
  fprintf(pkfile," ");
  //fprintf(pkfile,"\n");

  class_fprintf_int(pkfile,index_k_2,_TRUE_);
  class_fprintf_int(pkfile,index_k_3,_TRUE_);
  class_fprintf_double(pkfile,k_2/k_1,_TRUE_);
  class_fprintf_double(pkfile,k_3/k_1,_TRUE_);
  class_fprintf_double(pkfile,fNL,_TRUE_);
       fprintf(pkfile,"\n");

       fclose(pkfile);


       /*

          ppt->table_x2[index_md][index_x_2_x_3] =
            ppt->k[index_md][index_k_2]/ppt->k[index_md][index_k_1];
           ppt->table_x3[index_md][index_x_2_x_3] =
            ppt->k[index_md][index_k_3]/ppt->k[index_md][index_k_1];
           ppt->table_fNL[index_md][index_x_2_x_3] = fNL;
       */
           
  index_x_2_x_3 += 1;
          
            
      
   }

   }//end if scalars

   
    else if (index_md==ppt->index_md_tensors){
      
       _set_fields_at_tau_end_lqc_(ppt->index_tp_gw) =
              ppw->pvecmetric[ppw->index_mt_pk_t];


    if (ppt->has_bispectrum_lqc == _TRUE_) {


      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_tau_end_111) =
         pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.);
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_tau_end_110) =
         pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.);

        
      
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_tau_end_100) =
         pow(ppw->pvecmetric[ppw->index_mt_gw_im],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],0.);

        
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_tau_end_000) =
         pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_im],0.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.)
        *pow(ppw->pvecmetric[ppw->index_mt_gw_re],1.);


      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_111) =
        ppw->pvecmetric[ppw->index_mt_gw_111];
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_110) =
        ppw->pvecmetric[ppw->index_mt_gw_110];
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_100) =
        ppw->pvecmetric[ppw->index_mt_gw_100];
      _set_fields_at_tau_end_lqc_(ppt->index_tp_bispectrum_t_equi_000) =
        ppw->pvecmetric[ppw->index_mt_gw_000];
      
      

      
   }

       
  
   }//end if tensors

// 	printf("Ps(%e) = %e\n", ppt->k[index_md][index_k], pk_s_1);
 
  /** - free quantities allocated at the beginning of the routine */

  class_call(perturb_lqc_vector_free(ppw->pv),
             ppt->error_message,
             ppt->error_message);


  free(interval_limit);
  free(interval_approx);
  free(pvecback2);

  return _SUCCESS_;
}

int perturb_lqc_prepare_output(struct background_lqc * pba,
			   struct perturbs_lqc * ppt){

  char tmp[40];

  ppt->scalar_titles[0]='\0';
  ppt->tensor_titles[0]='\0';


  if (ppt->k_output_values_num > 0) {

    /** Write titles for all perturbations_lqc that we would like to print/store. */
    if (ppt->has_scalars == _TRUE_){

      class_store_columntitle(ppt->scalar_titles,"t [Planck Sec]",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"a",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"phi_scf_re",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"phi_scf_im",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"phi_prime_scf_re",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"phi_prime_scf_im",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"phi_scf_abs",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"v_re",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"v_im",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"v_prime_re",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"v_prime_im",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Hubble",_TRUE_);
      
      class_store_columntitle(ppt->scalar_titles,"Z_123_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123_111",_TRUE_);

        
      class_store_columntitle(ppt->scalar_titles,"Z_123d_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123d_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123d_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_123d_111",_TRUE_);

        
      class_store_columntitle(ppt->scalar_titles,"Z_12d3d_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3d_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3d_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3d_111",_TRUE_);

      
      
      class_store_columntitle(ppt->scalar_titles,"Y_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Y_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Y_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Y_111",_TRUE_);


      class_store_columntitle(ppt->scalar_titles,"field_redef",_TRUE_);

      
      
      class_store_columntitle(ppt->scalar_titles,"N",_TRUE_);



      class_store_columntitle(ppt->scalar_titles,"Z_12d3_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_12d3_111",_TRUE_);

        
      class_store_columntitle(ppt->scalar_titles,"Z_1d23_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23_111",_TRUE_);

        
      class_store_columntitle(ppt->scalar_titles,"Z_1d23d_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23d_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23d_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d23d_111",_TRUE_);

      class_store_columntitle(ppt->scalar_titles,"Z_1d2d3_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d2d3_100",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d2d3_110",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"Z_1d2d3_111",_TRUE_);

      
      class_store_columntitle(ppt->scalar_titles,"prefactor_123_000",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"a_over_z",_TRUE_);

      class_store_columntitle(ppt->scalar_titles,"prefactor_123",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"prefactor_1d23",_TRUE_);
      class_store_columntitle(ppt->scalar_titles,"prefactor_1d2d3",_TRUE_);


      ppt->number_of_scalar_titles =
        get_number_of_titles(ppt->scalar_titles);

     
    }

    if (ppt->has_tensors == _TRUE_){

      class_store_columntitle(ppt->tensor_titles,"t [Planck Sec]",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"a",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"H (gw)_re",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"H (gw)_im",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"Hdot (gwdot)_re",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"Hdot (gwdot)_im",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"H (gw)_abs",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"v_re",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"v_im",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"v_prime_re",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"v_prime_im",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"Hubble",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"I_000",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"I_100",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"I_110",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"I_111",_TRUE_);
      class_store_columntitle(ppt->tensor_titles,"N",_TRUE_);

      ppt->number_of_tensor_titles =
        get_number_of_titles(ppt->tensor_titles);

    }

  }
  return _SUCCESS_;

}



/**
 * Initialize the field '-->pv' of a perturb_lqc_workspace structure, which
 * is a perturb_lqc_vector structure. This structure contains indices and
 * values of all quantities which need to be integrated with respect
 * to time (and only them: quantities fixed analytically or obeying
 * constraint equations are NOT included in this vector). This routine
 * distinguishes between two cases:
 *
 * --> the input pa_old is set to the NULL pointer:
 *
 * This happens when we start integrating over a new wavenumber and we
 * want to set initial conditions for the perturbations_lqc. Then, it is
 * assumed that ppw-->pv is not yet allocated. This routine allocates
 * it, defines all indices, and then fills the vector ppw-->pv-->y with
 * the initial conditions defined in perturb_lqc_initial_conditions.
 *
 * --> the input pa_old is not set to the NULL pointer and describes
 * some set of approximations:
 *
 * This happens when we need to change approximation scheme while
 * integrating over a given wavenumber. The new approximation
 * described by ppw-->pa is then different from pa_old. Then, this
 * routine allocates a new vector with a new size and new index
 * values; it fills this vector with initial conditions taken from the
 * previous vector passed as an input in ppw-->pv, and eventually with
 * some analytic approximations for the new variables appearing at
 * this time; then the new vector comes in replacement of the old one,
 * which is freed.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: workspace containing in input the approximation scheme, the background/thermodynamics/metric quantities, and eventually the previous vector y; and in output the new vector y.
 * @param pa_old     Input: NULL is we need to set y to initial conditions for a new wavenumber; points towards a perturb_lqc_approximations if we want to switch of approximation.
 * @return the error status
 */

int perturb_lqc_vector_init(
                        struct precision * ppr,
                        struct background_lqc * pba,
                        struct perturbs_lqc * ppt,
                        int index_md,
                        double k,
                        double k2,
                        double k3,
                        double tau,
                        struct perturb_lqc_workspace * ppw, /* ppw->pv unallocated if pa_old = NULL, allocated and filled otherwise */
                        int pa_old
                        ) {

  /** Summary: */

  /** - define local variables */

  struct perturb_lqc_vector * ppv;

  int index_pt;
  int l;
  int n_ncdm,index_q,ncdm_l_size;
  double rho_plus_p_ncdm,q,q2,epsilon,a,factor;

  /** - allocate a new perturb_lqc_vector structure to which ppw-->pv will point at the end of the routine */

  class_alloc(ppv,sizeof(struct perturb_lqc_vector),ppt->error_message);
  //ppv->l_max_ncdm = NULL;
  //ppv->q_size_ncdm = NULL;


  /** - define all indices in this new vector (depends on approximation scheme, described by the input structure ppw-->pa) */

  index_pt = 0;

  if (_scalars_lqc_) {

    /* scalar field */

    class_define_index(ppv->index_pt_phi_scf_re,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_re,pba->has_scf,index_pt,1);
    class_define_index(ppv->index_pt_phi_scf_im,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_im,pba->has_scf,index_pt,1);

    class_define_index(ppv->index_pt_phi_scf_re_1,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_re_1,pba->has_scf,index_pt,1);
    class_define_index(ppv->index_pt_phi_scf_im_1,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_im_1,pba->has_scf,index_pt,1);

    class_define_index(ppv->index_pt_phi_scf_re_2,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_re_2,pba->has_scf,index_pt,1);
    class_define_index(ppv->index_pt_phi_scf_im_2,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_im_2,pba->has_scf,index_pt,1);


    class_define_index(ppv->index_pt_phi_scf_re_3,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_re_3,pba->has_scf,index_pt,1);
    class_define_index(ppv->index_pt_phi_scf_im_3,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_phi_prime_scf_im_3,pba->has_scf,index_pt,1);

    //1
    class_define_index(ppv->index_pt_bispectrum_s_prefactor_123_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123_111,pba->has_scf,index_pt,1); 


    class_define_index(ppv->index_pt_bispectrum_s_123d_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_123d_111,pba->has_scf,index_pt,1); 

    class_define_index(ppv->index_pt_bispectrum_s_12d3d_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3d_111,pba->has_scf,index_pt,1); 


    class_define_index(ppv->index_pt_bispectrum_s_12d3_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_12d3_111,pba->has_scf,index_pt,1); 



    class_define_index(ppv->index_pt_bispectrum_s_1d23_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23_111,pba->has_scf,index_pt,1); 



    class_define_index(ppv->index_pt_bispectrum_s_1d23d_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d23d_111,pba->has_scf,index_pt,1); 



    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_010,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_001,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_101,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_011,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_s_1d2d3_111,pba->has_scf,index_pt,1); 



    
  }
  if (_tensors_lqc_) {


    class_define_index(ppv->index_pt_gw_re,_TRUE_,index_pt,1);    
    class_define_index(ppv->index_pt_gwdot_re,_TRUE_,index_pt,1);  
    class_define_index(ppv->index_pt_gw_im,_TRUE_,index_pt,1);    
    class_define_index(ppv->index_pt_gwdot_im,_TRUE_,index_pt,1);

        
    class_define_index(ppv->index_pt_bispectrum_t_000,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_t_100,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_t_110,pba->has_scf,index_pt,1); 
    class_define_index(ppv->index_pt_bispectrum_t_111,pba->has_scf,index_pt,1); 

    
  }

  ppv->pt_size = index_pt;

  /** - allocate vectors for storing the values of all these
      quantities and their time-derivatives at a given time */

  class_calloc(ppv->y,ppv->pt_size,sizeof(double),ppt->error_message);
  class_alloc(ppv->dy,ppv->pt_size*sizeof(double),ppt->error_message);
  class_alloc(ppv->used_in_sources,ppv->pt_size*sizeof(int),ppt->error_message);

  /** - specify which perturbations_lqc are needed in the evaluation of source terms */

  /* take none of them by default */
  for (index_pt=0; index_pt<ppv->pt_size; index_pt++)
    ppv->used_in_sources[index_pt] = _FALSE_;
 
  /** - case of setting initial conditions for a new wavenumber */

  if (pa_old == 0) {
	printf("pa_old = %d\n",pa_old);
    if (ppt->perturbations_lqc_verbose>2)
      // fprintf(stdout,"Mode k=%e: initializing vector at tau=%e\n",k,tau);
    /** - --> (b) let ppw-->pv points towards the perturb_lqc_vector structure
        that we just created */
    ppw->pv = ppv;

    /** - --> (c) fill the vector ppw-->pv-->y with appropriate initial conditions */

    class_call(perturb_lqc_initial_conditions(ppr,
                                          pba,
                                          ppt,
                                          index_md,
                                          k,
                                          k2,
                                          k3,
                                          tau,
                                          ppw),
               ppt->error_message,
               ppt->error_message);



  }

 /** - case of switching approximation while a wavenumber is being integrated */

  
  else {
    //fprintf(stdout,"Mode k=%e: switching time variable at tau=%e\n",k,tau);
	printf("pa_old = %d\n",pa_old);

  if (_scalars_lqc_) {
   
  ppv->y[ppv->index_pt_phi_scf_re] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_re];
  ppv->y[ppv->index_pt_phi_scf_im] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_im];
  
 ppv->y[ppv->index_pt_phi_prime_scf_re] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re];
  ppv->y[ppv->index_pt_phi_prime_scf_im] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im];
  
  ppv->y[ppv->index_pt_phi_scf_re_1] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_re_1];
  ppv->y[ppv->index_pt_phi_scf_im_1] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_im_1];
  
 ppv->y[ppv->index_pt_phi_prime_scf_re_1] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1];
  ppv->y[ppv->index_pt_phi_prime_scf_im_1] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1];
  
  ppv->y[ppv->index_pt_phi_scf_re_2] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_re_2];
  ppv->y[ppv->index_pt_phi_scf_im_2] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_im_2];
  
 ppv->y[ppv->index_pt_phi_prime_scf_re_2] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2];
  ppv->y[ppv->index_pt_phi_prime_scf_im_2] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2];
  
  ppv->y[ppv->index_pt_phi_scf_re_3] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_re_3];
  ppv->y[ppv->index_pt_phi_scf_im_3] =
ppw->pv->y[ppw->pv->index_pt_phi_scf_im_3];
  
 ppv->y[ppv->index_pt_phi_prime_scf_re_3] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3];
  ppv->y[ppv->index_pt_phi_prime_scf_im_3] =
        ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3];

//printf("");
/*
/////////////////////////////////////////////////////////////////
//This part will have to be modified if an interval in cosmic_time 
//follows an interval in efolds. (20-06-17)
ppv->y[ppv->index_pt_phi_prime_scf_re] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re];
// *(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im];
// *(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_1];
// *(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_1];
// *(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_2];
// *(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_2];
// *(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_3];
// *(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_3];
// *(ppw->pvecback[pba->index_bg_H]);




////////////////////////////////////////////////////////////////  
*/
if(ppw->approx == (int)efolds){
  
ppv->y[ppv->index_pt_phi_prime_scf_re] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re]
*(1./ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im]
*(1./ppw->pvecback[pba->index_bg_H]);

ppv->y[ppv->index_pt_phi_prime_scf_re_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_1]
*(1./ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_1]
*(1./ppw->pvecback[pba->index_bg_H]);

ppv->y[ppv->index_pt_phi_prime_scf_re_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_2]
*(1./ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_2]
*(1./ppw->pvecback[pba->index_bg_H]);

ppv->y[ppv->index_pt_phi_prime_scf_re_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_3]
*(1./ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_3]
*(1./ppw->pvecback[pba->index_bg_H]);

// printf("Mode k=%e: resetting IC in efolds\n",k);

 
 }
/*   
 else {

ppv->y[ppv->index_pt_phi_prime_scf_re] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re]
*(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im]
*(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_1]
*(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_1] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_1]
*(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_2]
*(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_2] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_2]
*(ppw->pvecback[pba->index_bg_H]);


ppv->y[ppv->index_pt_phi_prime_scf_re_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_3]
*(ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_phi_prime_scf_im_3] =
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_3]
*(ppw->pvecback[pba->index_bg_H]);


 
// printf("Mode k=%e: resetting IC in cosmic time\n",k);

 }
  
 */
/* printf( "Q_re( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppv->y[ppv->index_pt_phi_scf_re]);
 printf( "Q_im( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppv->y[ppv->index_pt_phi_scf_im]);
 printf( "Q'_re( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppv->y[ppv->index_pt_phi_prime_scf_re]);
 printf( "Q'_im( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppv->y[ppv->index_pt_phi_prime_scf_im]);
*/
//1
  ppv->y[ppv->index_pt_bispectrum_s_prefactor_123_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_prefactor_123_000];
  ppv->y[ppv->index_pt_bispectrum_s_123_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_000];
  ppv->y[ppv->index_pt_bispectrum_s_123_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_100];
  ppv->y[ppv->index_pt_bispectrum_s_123_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_010];
  ppv->y[ppv->index_pt_bispectrum_s_123_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_001];
  ppv->y[ppv->index_pt_bispectrum_s_123_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_110];
  ppv->y[ppv->index_pt_bispectrum_s_123_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_101];
  ppv->y[ppv->index_pt_bispectrum_s_123_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_011];
  ppv->y[ppv->index_pt_bispectrum_s_123_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_111];



  //2
  ppv->y[ppv->index_pt_bispectrum_s_123d_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_000];
  ppv->y[ppv->index_pt_bispectrum_s_123d_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_100];
  ppv->y[ppv->index_pt_bispectrum_s_123d_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_010];
  ppv->y[ppv->index_pt_bispectrum_s_123d_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_001];
  ppv->y[ppv->index_pt_bispectrum_s_123d_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_110];
  ppv->y[ppv->index_pt_bispectrum_s_123d_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_101];
  ppv->y[ppv->index_pt_bispectrum_s_123d_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_011];
  ppv->y[ppv->index_pt_bispectrum_s_123d_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_111];



  //3
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_000];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_100];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_010];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_001];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_110];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_101];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_011];
  ppv->y[ppv->index_pt_bispectrum_s_12d3d_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_111];

  //4
  ppv->y[ppv->index_pt_bispectrum_s_12d3_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_000];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_100];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_010];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_001];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_110];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_101];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_011];
  ppv->y[ppv->index_pt_bispectrum_s_12d3_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_111];

  //5
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_000];
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_100];
 ppv->y[ppv->index_pt_bispectrum_s_1d23d_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_010];
 ppv->y[ppv->index_pt_bispectrum_s_1d23d_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_001];
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_110];
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_101];
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_011];
  ppv->y[ppv->index_pt_bispectrum_s_1d23d_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_111];

  //6
  ppv->y[ppv->index_pt_bispectrum_s_1d23_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_000];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_100];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_010];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_001];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_110];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_101];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_011];
  ppv->y[ppv->index_pt_bispectrum_s_1d23_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_111];


  //7
  ppv->y[ppv->index_pt_bispectrum_s_1d2d3_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_000];
  ppv->y[ppv->index_pt_bispectrum_s_1d2d3_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_100];
  ppv->y[ppv->index_pt_bispectrum_s_1d2d3_010] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_010];
   ppv->y[ppv->index_pt_bispectrum_s_1d2d3_001] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_001];
   ppv->y[ppv->index_pt_bispectrum_s_1d2d3_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_110];
   ppv->y[ppv->index_pt_bispectrum_s_1d2d3_101] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_101];
    ppv->y[ppv->index_pt_bispectrum_s_1d2d3_011] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_011];
   ppv->y[ppv->index_pt_bispectrum_s_1d2d3_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_111];

  
 }

  if (_tensors_lqc_) {


 
  ppv->y[ppv->index_pt_gw_re] =
ppw->pv->y[ppw->pv->index_pt_gw_re];
  ppv->y[ppv->index_pt_gw_im] =
ppw->pv->y[ppw->pv->index_pt_gw_im];
  
 ppv->y[ppv->index_pt_gwdot_re] =
        ppw->pvecmetric[ppw->index_mt_gwdot_re];
  ppv->y[ppv->index_pt_gwdot_im] =
        ppw->pvecmetric[ppw->index_mt_gwdot_im];
  
if(ppw->approx == (int)efolds){
ppv->y[ppv->index_pt_gwdot_re] =
  ppw->pv->y[ppw->pv->index_pt_gwdot_re]
*(1./ppw->pvecback[pba->index_bg_H]);
ppv->y[ppv->index_pt_gwdot_im] =
  ppw->pv->y[ppw->pv->index_pt_gwdot_im]
*(1./ppw->pvecback[pba->index_bg_H]);


//printf("Mode k=%e: resetting IC in efolds\n",k,tau);
//    printf("H =%e: resetting IC in efolds\n",ppw->pvecback[pba->index_bg_H]);
 }
 
 else {
   ppv->y[ppv->index_pt_gwdot_re] =
  ppw->pv->y[ppw->pv->index_pt_gwdot_re]
*(ppw->pvecback[pba->index_bg_H]);
   
ppv->y[ppv->index_pt_gwdot_im] =
  ppw->pv->y[ppw->pv->index_pt_gwdot_im]
*(ppw->pvecback[pba->index_bg_H]);
 }
  
  

  ppv->y[ppv->index_pt_bispectrum_t_000] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_000];
  ppv->y[ppv->index_pt_bispectrum_t_100] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_100];
  ppv->y[ppv->index_pt_bispectrum_t_110] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_110];
  ppv->y[ppv->index_pt_bispectrum_t_111] =
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_111];  
    
    }



    /** - --> (d) free the previous vector of perturbations */

    class_call(perturb_lqc_vector_free(ppw->pv),
               ppt->error_message,
               ppt->error_message);

    /** - --> (e) let ppw-->pv points towards the perturb_vector structure
        that we just created */

    ppw->pv = ppv;


  
     }
  //printf("Mode k=%e, approx = %d\n",k,ppw->approx);


  
  return _SUCCESS_;
}

/**
 * Free the perturb_lqc_vector structure.
 *
 * @param pv        Input: pointer to perturb_lqc_vector structure to be freed
 * @return the error status
 */

int perturb_lqc_vector_free(
                        struct perturb_lqc_vector * pv
                        ) {

  free(pv->y);
  free(pv->dy);
  free(pv->used_in_sources);
  free(pv);

  return _SUCCESS_;
}

/**
 * For each mode, wavenumber and initial condition, this function
 * initializes in the vector all values of perturbed variables (in a
 * given gauge). It is assumed here that all values have previously been
 * set to zero, only non-zero values are set here.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md   Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: workspace containing in input the approximation scheme, the background/thermodynamics/metric quantities, and eventually the previous vector y; and in output the new vector y.
 * @return the error status
 */

int perturb_lqc_initial_conditions(struct precision * ppr,
                               struct background_lqc * pba,
                               struct perturbs_lqc * ppt,
                               int index_md,
                               double k,
                               double k2,
                               double k3,
                               double tau,
                               struct perturb_lqc_workspace * ppw
                               ) {
  /** Summary: */

  /** --> Declare local variables */

 if(ppw->approx == (int)efolds) {
   //  printf("Mode k=%e, initial conditions in the efolds interval.\n",k);
  double tau_of_N;
if(ppw->num_approx_intervals>3 && ppw->approx_interval < 1){
class_call(background_lqc_tau_of_N_pb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


    class_call(background_lqc_at_tau_pb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
 }
else{
class_call(background_lqc_tau_of_N_fb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


    class_call(background_lqc_at_tau_fb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);

}
/*
    class_call(background_lqc_at_N(pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
*/
 }

    else {

      //printf("Mode k=%e, initial conditions in the cosmic time interval.\n",k);
    if(ppw->num_approx_intervals >2 && ppw->approx_interval <1){
 
    class_call(background_lqc_at_tau_pb(pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
    }
    else{
    class_call(background_lqc_at_tau_fb(pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);

   }

   }
  double a,a_prime_over_a;



    a = ppw->pvecback[pba->index_bg_a];

    a_prime_over_a = ppw->pvecback[pba->index_bg_H]*a;
   
    //printf("N=%e\n", ppw->pvecback[pba->index_bg_N]);
    // printf("a=%e\n", ppw->pvecback[pba->index_bg_a]);
    // printf("H=%e\n", ppw->pvecback[pba->index_bg_H]);

  if (_scalars_lqc_) {

        ppw->pv->y[ppw->pv->index_pt_phi_scf_re] =
          (1./a)
          *pow(2.*k,-0.5);
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re] =
         -(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(2.*k,-0.5);

        ppw->pv->y[ppw->pv->index_pt_phi_scf_im] = 0.;
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im] =
          -(1./a_prime_over_a)
          *(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(k/2.,0.5);


        ////////k1
        ppw->pv->y[ppw->pv->index_pt_phi_scf_re_1] =
          (1./a)
          *pow(2.*k,-0.5);
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_1] =
         -(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(2.*k,-0.5);

        ppw->pv->y[ppw->pv->index_pt_phi_scf_im_1] = 0.;
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_1] =
          -(1./a_prime_over_a)
          *(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(k/2.,0.5);

        
        //////k2
            ppw->pv->y[ppw->pv->index_pt_phi_scf_re_2] =
          (1./a)
          *pow(2.*k2,-0.5);
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_2] =
         -(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(2.*k2,-0.5);

        ppw->pv->y[ppw->pv->index_pt_phi_scf_im_2] = 0.;
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_2] =
          -(1./a_prime_over_a)
          *(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(k2/2.,0.5);

        
        //////k3
            ppw->pv->y[ppw->pv->index_pt_phi_scf_re_3] =
          (1./a)
          *pow(2.*k3,-0.5);
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_3] =
         -(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(2.*k3,-0.5);

        ppw->pv->y[ppw->pv->index_pt_phi_scf_im_3] = 0.;
        
        ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_3] =
          -(1./a_prime_over_a)
          *(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(k3/2.,0.5);

/*  printf( "Q_re( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppw->pv->y[ppw->pv->index_pt_phi_scf_re]);
 printf( "Q_im( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppw->pv->y[ppw->pv->index_pt_phi_scf_im]);
 printf( "Q'_re( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re]);
 printf( "Q'_im( %e ) = %e\n", ppw->pvecback[pba->index_bg_cosmic_time], ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im]);
  */ 


        //1
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_prefactor_123_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123_111] =0.;

  //2
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_123d_111] =0.;

  //3
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3_111] =0.;


  //4
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23_111] =0.;

  //5

  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d2d3_111] =0.;


  //6

  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_1d23d_111] =0.;

  //7
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_000] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_100] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_010] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_001] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_110] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_101] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_011] =0.;
  ppw->pv->y[ppw->pv->index_pt_bispectrum_s_12d3d_111] =0.;







        
        
if(ppw->approx == (int)efolds){
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re] *=
    (1./ppw->pvecback[pba->index_bg_H]);      
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im] *=
    (1./ppw->pvecback[pba->index_bg_H]);


  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_1] *=
    (1./ppw->pvecback[pba->index_bg_H]);      
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_1] *=
    (1./ppw->pvecback[pba->index_bg_H]);

  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_2] *=
    (1./ppw->pvecback[pba->index_bg_H]);      
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_2] *=
    (1./ppw->pvecback[pba->index_bg_H]);

  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re_3] *=
    (1./ppw->pvecback[pba->index_bg_H]);      
  ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_im_3] *=
    (1./ppw->pvecback[pba->index_bg_H]);


  
  //printf("Mode k=%e: initializing vector with phi_prime_scf_re =%e\n",
  //      k,ppw->pv->y[ppw->pv->index_pt_phi_prime_scf_re]);

         }
        
    }
  /** --> For tensors */

  if (_tensors_lqc_) {
    

        ppw->pv->y[ppw->pv->index_pt_gw_re] =
          (1./a)
          *pow(2.*k,-0.5);

        ppw->pv->y[ppw->pv->index_pt_gwdot_re] =
          -(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(2.*k,-0.5);

        ppw->pv->y[ppw->pv->index_pt_gw_im] = 0.;
        ppw->pv->y[ppw->pv->index_pt_gwdot_im] =
          -(1./a_prime_over_a)
          *(1./a)
          *ppw->pvecback[pba->index_bg_H]
          *pow(k/2.,0.5);



        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_000] = 0.;
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_100] = 0.;
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_110] = 0.;
        ppw->pv->y[ppw->pv->index_pt_bispectrum_t_111] = 0.;



if(ppw->approx == (int)efolds){
  ppw->pv->y[ppw->pv->index_pt_gwdot_re] *=
    (1./ppw->pvecback[pba->index_bg_H]);      
  ppw->pv->y[ppw->pv->index_pt_gwdot_im] *=
    (1./ppw->pvecback[pba->index_bg_H]);

  // printf("Mode k=%e: initializing vector with gwdot_re =%e\n",
  //      k,ppw->pv->y[ppw->pv->index_pt_gwdot_re]);

         }
              
        

  }

  return _SUCCESS_;
}


/**
 * Compute typical timescale over which the perturbation equations
 * vary. Some integrators (e.g. Runge-Kunta) benefit from calling this
 * routine at each step in order to adapt the next step.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know the content of this pointer.
 *   error_message passed in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param parameters_and_workspace Input: fixed parameters (e.g. indices), workspace, approximation used, etc.
 * @param timescale                Output: perturbation variation timescale (given the approximation used)
 * @param error_message            Output: error message
 */

int perturb_lqc_timescale(
                      double tau,
                      void * parameters_and_workspace,
                      double * timescale,
                      ErrorMsg error_message
                      ) {

  double tau_k;
  double tau_h;
  double a;

  /* various pointers allowing to extract the fields of the
     parameter_and_workspace input structure */
  struct perturb_lqc_parameters_and_workspace * pppaw;
  struct background_lqc * pba;
  struct perturbs_lqc * ppt;
  struct perturb_lqc_workspace * ppw;
  double * pvecback;

  pppaw = parameters_and_workspace;
  pba = pppaw->pba;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;
  pvecback = ppw->pvecback;


  /** - compute Fourier mode time scale = \f$ \tau_k = 1/k \f$ */

  class_test(pppaw->k == 0.,
             ppt->error_message,
             "stop to avoid division by zero");

  
       if (ppw->approx == (int)efolds){
 double tau_of_N;
/*class_call(background_lqc_tau_of_N(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


    class_call(background_lqc_at_tau(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
*/
/*    
         class_call(background_lqc_at_N(pba,tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), pvecback),
             pba->error_message,
             error_message);*/
//Modified on 17/07/17

if(ppw->num_approx_intervals > 3 && ppw->approx_interval < 1){
class_call(background_lqc_tau_of_N_pb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


    class_call(background_lqc_at_tau_pb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 pvecback),
               pba->error_message,
               ppt->error_message);


}
else{
class_call(background_lqc_tau_of_N_fb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


    class_call(background_lqc_at_tau_fb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 pvecback),
               pba->error_message,
               ppt->error_message);

}
         }
       else{

/*         class_call(background_lqc_at_tau(pba,tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), pvecback),
             pba->error_message,
             error_message); */
//Modified on 17/07/17

        if(ppw->num_approx_intervals > 2 && ppw->approx_interval < 1){

         class_call(background_lqc_at_tau_pb(pba,tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), pvecback),
             pba->error_message,
             error_message); 
	}
	else{
         class_call(background_lqc_at_tau_fb(pba,tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), pvecback),
             pba->error_message,
             error_message);  }
	}	

   a = pvecback[pba->index_bg_a];
   //printf("a=%e\n", pvecback[pba->index_bg_a]);

//printf("pia = %e, pia = %e \n",pvecback[pba->index_bg_pia],ppw->pvecback[pba->index_bg_pia]);
   
   tau_h = 1./sqrt(pvecback[pba->index_bg_rho_scf]/3.);
   tau_k = a/pppaw->k;

       *timescale = MIN(fabs(tau_h),fabs(tau_k));
       if(ppw->approx == (int)efolds){
       *timescale *= sqrt(pvecback[pba->index_bg_rho_scf]/3.);
       }
       
  return _SUCCESS_;
}


int perturb_lqc_einstein(
                     struct precision * ppr,
                     struct background_lqc * pba,
                     struct perturbs_lqc * ppt,
                     int index_md,
                     double k,
                     double k_2,
                     double k_3,
                     double tau,
                     double * y,
                     struct perturb_lqc_workspace * ppw
                     ) {


  double k2,k2_1,k2_2,k2_3,a,a2,a_prime_over_a;
  double K2;
  double H;
  double k_1 = k;
  H = ppw->pvecback[pba->index_bg_H];
 

  

  k2 = k*k;

  k2_1 = k2; 
  k2_2 = k_2*k_2; 
  k2_3 = k_3*k_3;

  
  a = ppw->pvecback[pba->index_bg_a];
  a2 = a * a;
  a_prime_over_a = ppw->pvecback[pba->index_bg_H]*a;
  //K2 = k2/pow( a_prime_over_a,2.);


  double K = fabs(k/a/pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.));
  double K_1 = fabs(k_1/a/pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.));
  double K_2 = fabs(k_2/a/pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.));
  double K_3 = fabs(k_3/a/pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.));
  
  double cut_off_bispectrum;
  double cut_off_bispectrum_1;
  double cut_off_bispectrum_2;
  double cut_off_bispectrum_3;
 if (_scalars_lqc_)  {
   
/*double sqrt_z_prime_prime_over_z =
  a*pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.)
  *sqrt(fabs(ppw->pvecback[pba->index_bg_f_S]));*/
double sqrt_z_prime_prime_over_z = 
//sqrt(fabs(ppw->pvecback[pba->index_bg_zT_primeprime_over_zT]));
//a*pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1./2.);
	sqrt(fabs(ppw->pvecback[pba->index_bg_zS_primeprime_over_zS]));
// if(ppw->approx == (int)efolds){
// cut_off_bispectrum = exp(-K/ppt->delta_cut_off_bispectrum);
   cut_off_bispectrum_1 = exp(-k_1/ppt->delta_cut_off_bispectrum/sqrt_z_prime_prime_over_z/3.);
   cut_off_bispectrum_2 = exp(-k_2/ppt->delta_cut_off_bispectrum/sqrt_z_prime_prime_over_z/3.);
   cut_off_bispectrum_3 = exp(-k_3/ppt->delta_cut_off_bispectrum/sqrt_z_prime_prime_over_z/3.);
   cut_off_bispectrum =  cut_off_bispectrum_1*cut_off_bispectrum_2*cut_off_bispectrum_3;
//    }
//  else cut_off_bispectrum = 1.;
   cut_off_bispectrum = exp(-(K_1 + K_2 + K_3)/3./ppt->delta_cut_off_bispectrum);

 }

   
 else cut_off_bispectrum = exp(-K/ppt->delta_cut_off_bispectrum);
    
    double prefactor_bispectrum_s =
      (3.*pow(ppw->pvecback[pba->index_bg_a],5.)/(k*k))
                                   *(4.*pba->rho_bounce*pba->rho_bounce)
                                   *pow(ppw->pvecback[pba->index_bg_y],1.)
                                   *cut_off_bispectrum;
    
    double prefactor_bispectrum_t = (1./4.)*a*cut_off_bispectrum;
    
  
  double prefactor;
//  printf("U = %e\n", ppw->pvecback[pba->index_bg_a]);
  if (_scalars_lqc_) {
    prefactor =  prefactor_bispectrum_s;

    //pvecmetric always contains the cosmic time quantities

 ppw->pvecmetric[ppw->index_mt_phi_scf_re] =  y[ppw->pv->index_pt_phi_scf_re];
 ppw->pvecmetric[ppw->index_mt_phi_scf_im] =  y[ppw->pv->index_pt_phi_scf_im];

    
 ppw->pvecmetric[ppw->index_mt_phi_scf_re_1] =  y[ppw->pv->index_pt_phi_scf_re_1];
 ppw->pvecmetric[ppw->index_mt_phi_scf_im_1] =  y[ppw->pv->index_pt_phi_scf_im_1];

    
 ppw->pvecmetric[ppw->index_mt_phi_scf_re_2] =  y[ppw->pv->index_pt_phi_scf_re_2];
 ppw->pvecmetric[ppw->index_mt_phi_scf_im_2] =  y[ppw->pv->index_pt_phi_scf_im_2];

    
 ppw->pvecmetric[ppw->index_mt_phi_scf_re_3] =  y[ppw->pv->index_pt_phi_scf_re_3];
 ppw->pvecmetric[ppw->index_mt_phi_scf_im_3] =  y[ppw->pv->index_pt_phi_scf_im_3];

    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re] =
      y[ppw->pv->index_pt_phi_prime_scf_re];
      
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im] =
      y[ppw->pv->index_pt_phi_prime_scf_im];
    

    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1] =
      y[ppw->pv->index_pt_phi_prime_scf_re_1];
      
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1] =
      y[ppw->pv->index_pt_phi_prime_scf_im_1];
    

    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2] =
      y[ppw->pv->index_pt_phi_prime_scf_re_2];
      
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2] =
      y[ppw->pv->index_pt_phi_prime_scf_im_2];
    

    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3] =
      y[ppw->pv->index_pt_phi_prime_scf_re_3];
      
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3] =
      y[ppw->pv->index_pt_phi_prime_scf_im_3];
    

    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re]/a
      -(k2-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_re]/a/a;



    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_1] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1]/a
      -(k2_1-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_re_1]/a/a;



    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_2] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2]/a
      -(k2_2-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_re_2]/a/a;



    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_3] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3]/a
      -(k2_3-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_re_3]/a/a;

    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im]/a
      -(k2-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_im]/a/a;

    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_1] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1]/a
      -(k2_1-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_im_1]/a/a;

    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_2] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2]/a
      -(k2_2-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_im_2]/a/a;

    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_3] =
      -3.*a_prime_over_a
      *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3]/a
      -(k2_3-
        a
        *a
        *ppw->pvecback[pba->index_bg_H2_f_S_Q])
      *ppw->pvecmetric[ppw->index_mt_phi_scf_im_3]/a/a;

    
if(ppw->approx == (int)efolds){

ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re] *=
  (ppw->pvecback[pba->index_bg_H]);
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im] *=
  (ppw->pvecback[pba->index_bg_H]);

   
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1] *=
  (ppw->pvecback[pba->index_bg_H]);
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1] *=
  (ppw->pvecback[pba->index_bg_H]);

   
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2] *=
  (ppw->pvecback[pba->index_bg_H]);
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2] *=
  (ppw->pvecback[pba->index_bg_H]);

   
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3] *=
  (ppw->pvecback[pba->index_bg_H]);
ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3] *=
  (ppw->pvecback[pba->index_bg_H]);

   
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im]
      -(K*K-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im];




    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_1] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_1]
      -(K_1*K_1-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_1];


    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_2] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_2]
      -(K_2*K_2-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_2];

              

    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_3] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_3]
      -(K_3*K_3-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_3];

       
    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re]
      -(K*K-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re];


    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_1] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_1]
      -(K_1*K_1-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_1];


    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_2] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_2]
      -(K_2*K_2-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_2];


    
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_3] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_3]
      -(K_3*K_3-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_3];




    
if(ppw->pvecback[pba->index_bg_H]<0.){
  
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im]
      -(K*K-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im];

       
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re]
      -(K*K-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re];




    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_1] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_1]
      -(K_1*K_1-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_1];

       
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_1] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_1]
      -(K_1*K_1-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_1];


    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_2] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_2]
      -(K_2*K_2-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_2];

       
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_2] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_2]
      -(K_2*K_2-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_2];



    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_3] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im_3]
      -(K_3*K_3-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_im_3];

       
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_3] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re_3]
      -(K_3*K_3-
        ppw->pvecback[pba->index_bg_H2_f_S_Q]
        /pow(ppw->pvecback[pba->index_bg_rho_scf]/3.,1.))
      *y[ppw->pv->index_pt_phi_scf_re_3];

    

    
    /*
    if(a*a*pba->m_scf_lqc*pba->m_scf_lqc>10.*k*k){

  
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_im] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_im];

       
    ppw->pvecmetric[ppw->index_mt_phi_scf_prime_prime_re] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_phi_prime_scf_re];

      }*/
      
    }
    }
    
    

 ppw->pvecmetric[ppw->index_mt_phi_scf_abs] =
   pow(pow(ppw->pvecmetric[ppw->index_mt_phi_scf_re],2.)
       +pow(ppw->pvecmetric[ppw->index_mt_phi_scf_im],2.),1./2.);



 double R_s_re_1,R_s_re_2,R_s_re_3,R_s_im_1,R_s_im_2,R_s_im_3;
 
      
     if ( ppw->pvecback[pba->index_bg_y] != 0)
       {
  double One_over_z_S =   (1./ (ppw->pvecback[pba->index_bg_a]
                                   *sqrt(6.)
                                       *ppw->pvecback[pba->index_bg_y]))
           *sqrt((pow(ppw->pvecback[pba->index_bg_x],2.)
                  +pow(ppw->pvecback[pba->index_bg_y],2.))
                 *(1./*-(pow(ppw->pvecback[pba->index_bg_x],2.)
                       +pow(ppw->pvecback[pba->index_bg_y],2.))*/));
  double One_over_Hz_S = (1./ (ppw->pvecback[pba->index_bg_a]
                                   *sqrt(2.*pba->rho_bounce)
                               *ppw->pvecback[pba->index_bg_y]));


 ppw->pvecmetric[ppw->index_mt_R_s_re] =
   ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_re]; 

 ppw->pvecmetric[ppw->index_mt_R_s_im] =
   ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_im];


R_s_re_1 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_re_1]; 


R_s_re_2 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_re_2]; 


R_s_re_3 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_re_3]; 


R_s_im_1 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_im_1]; 


R_s_im_2 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_im_2]; 


R_s_im_3 =
  ppw->pvecback[pba->index_bg_a]
   *One_over_z_S
   *ppw->pvecmetric[ppw->index_mt_phi_scf_im_3]; 


 
 ppw->pvecmetric[ppw->index_mt_R_s_prime_re] =
   (2.-ppw->pvecback[pba->index_bg_epsilon_H]
    +ppw->pvecback[pba->index_bg_mu]
    *ppw->pvecback[pba->index_bg_x]
    /ppw->pvecback[pba->index_bg_y])
   *ppw->pvecmetric[ppw->index_mt_R_s_re]
   +ppw->pvecback[pba->index_bg_a]
   *(One_over_z_S*ppw->pvecmetric[ppw->index_mt_phi_scf_re]
     +One_over_Hz_S*y[ppw->pv->index_pt_phi_prime_scf_re]);

   
 ppw->pvecmetric[ppw->index_mt_R_s_prime_im] =
   (2.-ppw->pvecback[pba->index_bg_epsilon_H]
    +ppw->pvecback[pba->index_bg_mu]
    *ppw->pvecback[pba->index_bg_x]
    /ppw->pvecback[pba->index_bg_y])
   *ppw->pvecmetric[ppw->index_mt_R_s_im]
   +ppw->pvecback[pba->index_bg_a]
   *(One_over_z_S*ppw->pvecmetric[ppw->index_mt_phi_scf_im]
     +One_over_Hz_S*y[ppw->pv->index_pt_phi_prime_scf_im]);
       }
     
else
{
 ppw->pvecmetric[ppw->index_mt_R_s_re] = 0.;
 ppw->pvecmetric[ppw->index_mt_R_s_im] = 0.;
 ppw->pvecmetric[ppw->index_mt_R_s_prime_re] = 0.;
 ppw->pvecmetric[ppw->index_mt_R_s_prime_im] = 0.;
 }

     

 /////////////////////////
     /*
     double k_1 = k_1;
     double k_2 = k_2;
     double k_3 = k_3;
     */
     double k_1_dot_k_2 =
       (k_3*k_3
       -k_1*k_1
        -k_2*k_2)
       /2.;
     double k_1_dot_k_3 =
       (k_2*k_2
       -k_1*k_1
        -k_3*k_3)
       /2.;
     double k_2_dot_k_3 =
       (k_1*k_1
       -k_2*k_2
        -k_3*k_3)
       /2.;

     double h12 = k_1_dot_k_2/k_1/k_1;
     double h21 = k_1_dot_k_2/k_2/k_2;
     double h13 = k_1_dot_k_3/k_1/k_1;
     double h31 = k_1_dot_k_3/k_3/k_3;
     double h23 = k_2_dot_k_3/k_2/k_2;
     double h32 = k_2_dot_k_3/k_3/k_3;

     double h12_s = (h12+h21)/2.;
     double h13_s = (h13+h31)/2.;
     double h23_s = (h23+h32)/2.;


     //V. Sreenath
     double phi_dot =  ppw->pvecback[pba->index_bg_phi_prime_scf];
     double V_phi =   ppw->pvecback[pba->index_bg_V_phi];
     double V_phiphi =   ppw->pvecback[pba->index_bg_V_phiphi];
     double V_phiphiphi =   ppw->pvecback[pba->index_bg_V_phiphiphi];
     double Pphi =   ppw->pvecback[pba->index_bg_Pphi_scf];
     double PIa =   ppw->pvecback[pba->index_bg_pia];

     double rho = (pow(ppw->pvecback[pba->index_bg_x],2.)
              +pow(ppw->pvecback[pba->index_bg_y],2.))*pba->rho_bounce;

/*     double z = //ppw->pvecback[pba->index_bg_z_S];
      sqrt(6.)*a
      *ppw->pvecback[pba->index_bg_y]
       /sqrt((pow(ppw->pvecback[pba->index_bg_x],2.)
              +pow(ppw->pvecback[pba->index_bg_y],2.)));

             double z_dot =
             (-2.*z
             -ppw->pvecback[pba->index_bg_mu]
              *ppw->pvecback[pba->index_bg_x]
              *sqrt(6.)*a
       /sqrt((pow(ppw->pvecback[pba->index_bg_x],2.)
              +pow(ppw->pvecback[pba->index_bg_y],2.)))
             +ppw->pvecback[pba->index_bg_epsilon_H]*z)
             *ppw->pvecback[pba->index_bg_H];

     double z_tilde = z/a;
     double z_tilde_dot = -(3. - ppw->pvecback[pba->index_bg_epsilon_H])*phi_dot
                - V_phi*ppw->pvecback[pba->index_bg_H]/(rho/3.);*/

/*     double z_tilde_dot = 
             z_dot/a
             -ppw->pvecback[pba->index_bg_H]
             *z/a;*/

//     double phi_dot = ppw->pvecback[pba->index_bg_phi_prime_scf];

     double t123 = (k_1*k_1+k_2*k_2+k_3*k_3)/(a*a);

     double g23 =
       k_2_dot_k_3*k_2_dot_k_3
       /k_2/k_2/k_3/k_3-1.;
     
       double g13 =
       k_1_dot_k_3*k_1_dot_k_3
       /k_1/k_1/k_3/k_3-1.;

   double g12 =
       k_1_dot_k_2*k_1_dot_k_2
       /k_2/k_2/k_1/k_1-1.;
    
       
/*   double prefactor_123 = 
     -(-t123
      +9.*phi_dot*phi_dot
     -(3./2.)*z_tilde*z_tilde*phi_dot*phi_dot
     -6.*pba->m_scf_lqc*pba->m_scf_lqc
     -(1./2.)*(g12+g13+g23)*z_tilde_dot*z_tilde_dot*z_tilde//V. Sreenath
     +3.*phi_dot*z_tilde_dot
     )
   *(1./4.)*a*a*a*z_tilde
   *cut_off_bispectrum;*/

   double prefactor_123 = -(-( (243.*pow(Pphi,7.)*pow(a,-8.)*pow(PIa,-5.)/2. 
	- 81.*pow(Pphi,5.)*pow(a*a*PIa,-3.)/2. + 27.*pow(Pphi,3.)*pow(a,-4.)/PIa/8. 
	+ 81.*pow(Pphi/PIa,4.)*V_phi/a - 27.*a*pow(Pphi/PIa,2.)*V_phi/2 + 27.*pow(a*a/PIa,3.)*Pphi*V_phi*V_phi/2.)
	*2.*(-g12 -g23 -g13)
	+ 3.*a*a*Pphi/PIa/2.*(-t123) + 6.*9.*a*pow(Pphi/PIa,2.)*V_phi/2. //010717
	- 6.*3.*a*a*Pphi*V_phiphi/PIa/2. + 6.*a*a*a*V_phiphiphi/6.)
	 *cut_off_bispectrum);
	
 
/*   double prefactor_123_000 = 
        -(-t123
      +9.*phi_dot*phi_dot
     -(3./2.)*z_tilde*z_tilde*phi_dot*phi_dot
     -6.*pba->m_scf_lqc*pba->m_scf_lqc
     -(1./2.)*(g12+g13+g23)*z_tilde_dot*z_tilde_dot*z_tilde//V. Sreenath
     +3.*phi_dot*z_tilde_dot
     )
   *(1./4.)*a*a*a*z_tilde
        *cut_off_bispectrum;
   //prefactor_123 = a*a*cut_off_bispectrum;*/
   double prefactor_123_000 = -(-( (243.*pow(Pphi,7.)*pow(a,-8.)*pow(PIa,-5.)/2. 
        - 81.*pow(Pphi,5.)*pow(a*a*PIa,-3.)/2. + 27.*pow(Pphi,3.)*pow(a,-4.)/PIa/8.
        + 81.*pow(Pphi/PIa,4.)*V_phi/a - 27.*a*pow(Pphi/PIa,2.)*V_phi/2 + 27.*pow(a*a/PIa,3.)*Pphi*V_phi*V_phi/2.)
        *2.*(-g12 -g23 -g13)
        + 3.*a*a*Pphi/PIa/2.*(-t123) + 6.*9.*a*pow(Pphi/PIa,2.)*V_phi/2.  //010717
        - 6.*3.*a*a*Pphi*V_phiphi/PIa/2. + 6.*a*a*a*V_phiphiphi/6.)
         *cut_off_bispectrum);


/* double prefactor_1d23 =
   -(-8.*z_tilde_dot*h23_s
    +2.*z_tilde*z_tilde*phi_dot
    +z_tilde*z_tilde*z_tilde_dot*(g12+g13)
    )*(1./8.)*a*a*a
   *cut_off_bispectrum;*/
   double prefactor_1d23 = -( -a*a*a*(
	( 81.*pow(Pphi,5.)*pow(a,-7.)*pow(PIa,-4.) - 27.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ 27.*pow(Pphi,2.)*V_phi*pow(PIa,-3.))*(-g12 -g13) - 2.*9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ ( -3.*Pphi*pow(a,-3.)/2. + 9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2) + 3.*a*a*V_phi/PIa )*h23_s)
	*cut_off_bispectrum);
	


/* double prefactor_12d3 =
   -(-8.*z_tilde_dot*h13_s
    +2.*z_tilde*z_tilde*phi_dot
    +z_tilde*z_tilde*z_tilde_dot*(g23+g12)
    )*(1./8.)*a*a*a
   *cut_off_bispectrum;*/
   double prefactor_12d3 = -( -a*a*a*(
	( 81.*pow(Pphi,5.)*pow(a,-7.)*pow(PIa,-4.) - 27.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ 27.*pow(Pphi,2.)*V_phi*pow(PIa,-3.))*(-g12 -g23) - 2.*9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ ( -3.*Pphi*pow(a,-3.)/2. + 9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2) + 3.*a*a*V_phi/PIa )*h13_s)
	*cut_off_bispectrum);



/* double prefactor_123d =
   -(-8.*z_tilde_dot*h12_s
    +2.*z_tilde*z_tilde*phi_dot
    +z_tilde*z_tilde*z_tilde_dot*(g23+g13)
    )*(1./8.)*a*a*a
   *cut_off_bispectrum;*/
   double prefactor_123d = -( -a*a*a*(
	( 81.*pow(Pphi,5.)*pow(a,-7.)*pow(PIa,-4.) - 27.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ 27.*pow(Pphi,2.)*V_phi*pow(PIa,-3.))*(-g13 -g23) - 2.*9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
	+ ( -3.*Pphi*pow(a,-3.)/2. + 9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2) + 3.*a*a*V_phi/PIa )*h12_s)
	*cut_off_bispectrum);




  
/* double prefactor_12d3d =
   -(4.*z_tilde*(h31+h21-1.)
    -g23*z_tilde*z_tilde*z_tilde
    )*(1./8.)*a*a*a*cut_off_bispectrum;*/
   double prefactor_12d3d = -( -pow(a,6.)*(
	27.*pow(Pphi,3.)*pow(a*a*PIa,-3)*(-g23) - 3.*Pphi*pow(a,-4.)/PIa
	+ 3.*Pphi*pow(a,-4.)/PIa*(h12 + h13))
	*cut_off_bispectrum);


/* double prefactor_1d23d =
   -(4.*z_tilde*(h32+h12-1.)
    -g13*z_tilde*z_tilde*z_tilde
    )*(1./8.)*a*a*a*cut_off_bispectrum;*/
   double prefactor_1d23d = -( -pow(a,6.)*(
	27.*pow(Pphi,3.)*pow(a*a*PIa,-3)*(-g13) - 3.*Pphi*pow(a,-4.)/PIa
	+ 3.*Pphi*pow(a,-4.)/PIa*(h12 + h23))
	*cut_off_bispectrum);


/* double prefactor_1d2d3 =
   -(4.*z_tilde*(h23+h13-1.)
    -g12*z_tilde*z_tilde*z_tilde
    )*(1./8.)*a*a*a*cut_off_bispectrum;*/
   double prefactor_1d2d3 = -( -pow(a,6.)*(
	27.*pow(Pphi,3.)*pow(a*a*PIa,-3)*(-g12) - 3.*Pphi*pow(a,-4.)/PIa
	+ 3.*Pphi*pow(a,-4.)/PIa*(h13 + h23))
	*cut_off_bispectrum);


 
if(ppw->approx == (int)efolds){
 prefactor_123_000 = prefactor_123_000/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_123 = prefactor_123/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_1d23 = prefactor_1d23/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_12d3 = prefactor_12d3/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_123d = prefactor_123d/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_1d2d3 = prefactor_1d2d3/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_1d23d = prefactor_1d23d/fabs(ppw->pvecback[pba->index_bg_H]);
 prefactor_12d3d = prefactor_12d3d/fabs(ppw->pvecback[pba->index_bg_H]);
}
 

 double Q_s_1_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_1];
 double Q_s_1d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1];
 double Q_s_2_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_2];
 double Q_s_2d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2];
 double Q_s_3_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_3];
 double Q_s_3d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3];
 
 double Q_s_1_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_1];
 double Q_s_1d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1];
 double Q_s_2_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_2];
 double Q_s_2d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2];
 double Q_s_3_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_3];
 double Q_s_3d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3];


 double Q_s_1_abs =pow(pow(Q_s_1_im,2.)+pow(Q_s_1_re,2.),0.5);
 double Q_s_2_abs =pow(pow(Q_s_2_im,2.)+pow(Q_s_2_re,2.),0.5);
 double Q_s_3_abs =pow(pow(Q_s_3_im,2.)+pow(Q_s_3_re,2.),0.5);


/* double V_phi =
   pba->m_scf_lqc
   *pba->m_scf_lqc
   *ppw->pvecback[pba->index_bg_phi_scf];*/
 
 double epsilon_2 =
   -6.
   -2.*V_phi
   /ppw->pvecback[pba->index_bg_H]
   /ppw->pvecback[pba->index_bg_phi_prime_scf]
   +ppw->pvecback[pba->index_bg_phi_prime_scf]
   *ppw->pvecback[pba->index_bg_phi_prime_scf]
   /ppw->pvecback[pba->index_bg_H]
   /ppw->pvecback[pba->index_bg_H];

  ppw->pvecmetric[ppw->index_mt_bispectrum_s_field_redef]=
   epsilon_2/2.
   *(pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
         *Q_s_1_abs,2.)
     *pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
         *Q_s_2_abs,2.)
    +pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
         *Q_s_1_abs,2.)
     *pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
         *Q_s_3_abs,2.)
     +pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
         *Q_s_2_abs,2.)
     *pow(ppw->pvecback[pba->index_bg_H]
        /ppw->pvecback[pba->index_bg_phi_prime_scf]
          *Q_s_3_abs,2.));
   

 
 //1
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_prefactor_123_000]=
   prefactor_123_000;

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_000]=
 prefactor_123
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_100]=
 prefactor_123
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_010]=
 prefactor_123
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_001]=
 prefactor_123
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_110]=
 prefactor_123
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_011]=
 prefactor_123
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_101]=
 prefactor_123
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_111]=
 prefactor_123
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);
 
 //2
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_000]=
 prefactor_123d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_100]=
 prefactor_123d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_010]=
 prefactor_123d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_001]=
 prefactor_123d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_110]=
 prefactor_123d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_011]=
 prefactor_123d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_101]=
 prefactor_123d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_111]=
 prefactor_123d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,0.);
 
 //3 
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_000]=
 prefactor_12d3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_100]=
 prefactor_12d3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_010]=
 prefactor_12d3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_001]=
 prefactor_12d3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,0.);
 

 
  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_110]=
 prefactor_12d3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_011]=
 prefactor_12d3
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_101]=
 prefactor_12d3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_111]=
 prefactor_12d3
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,0.);
 
 //4
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_000]=
 prefactor_1d23
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_100]=
 prefactor_1d23
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_010]=
 prefactor_1d23
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_001]=
 prefactor_1d23
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_110]=
 prefactor_1d23
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_011]=
 prefactor_1d23
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_101]=
 prefactor_1d23
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_111]=
 prefactor_1d23
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);


 //5
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_000]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_100]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_010]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_001]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_110]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_011]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_101]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_111]=
 prefactor_1d2d3
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3_re,0.);
 
 //6
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_000]=
 prefactor_1d23d
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_100]=
 prefactor_1d23d
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_010]=
 prefactor_1d23d
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_001]=
 prefactor_1d23d
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_110]=
 prefactor_1d23d
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_011]=
 prefactor_1d23d
       *pow(Q_s_1d_im,0.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1d_re,1.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_101]=
 prefactor_1d23d
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_111]=
 prefactor_1d23d
       *pow(Q_s_1d_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1d_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3d_re,0.);
 

 //7
 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_000]=
 prefactor_12d3d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_100]=
 prefactor_12d3d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3d_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_010]=
 prefactor_12d3d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3d_re,1.);


  ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_001]=
 prefactor_12d3d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3d_re,0.);
 

  

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_110]=
 prefactor_12d3d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3d_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3d_re,1.);


 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_011]=
 prefactor_12d3d
       *pow(Q_s_1_im,0.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_101]=
 prefactor_12d3d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,0.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,1.)
       *pow(Q_s_3d_re,0.);

 ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_111]=
 prefactor_12d3d
       *pow(Q_s_1_im,1.)
       *pow(Q_s_2d_im,1.)
       *pow(Q_s_3d_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2d_re,0.)
       *pow(Q_s_3d_re,0.);
 
 
 


 
 ////////////////////////
 
   ppw->pvecmetric[ppw->index_mt_v_s_im] =
     a*ppw->pvecmetric[ppw->index_mt_phi_scf_im];
   
   ppw->pvecmetric[ppw->index_mt_v_s_re] =
     a*ppw->pvecmetric[ppw->index_mt_phi_scf_re];
   
   ppw->pvecmetric[ppw->index_mt_v_s_prime_im] =
          a*ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im]
          +H*ppw->pvecmetric[ppw->index_mt_v_s_im];
   
   ppw->pvecmetric[ppw->index_mt_v_s_prime_re] =
          a*ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re]
          +H*ppw->pvecmetric[ppw->index_mt_v_s_re];

   
 ppw->pvecmetric[ppw->index_mt_pk_s] =
   pow(k,3.)*(
   pow(ppw->pvecmetric[ppw->index_mt_R_s_im],2.)
   +pow(ppw->pvecmetric[ppw->index_mt_R_s_re],2.))
   *(1./(2.*_PI_*_PI_));

 /*
 ppw->pvecmetric[ppw->index_mt_pk_s] =
   pow(k,3.)*(
   pow(ppw->pvecmetric[ppw->index_mt_phi_scf_im],2.)
   +pow(ppw->pvecmetric[ppw->index_mt_phi_scf_re],2.))
   *(1./(2.*_PI_*_PI_));
   */

 ppw->pvecmetric[ppw->index_mt_pk_s_1] =
   pow(k_1,3.)*(
   pow(R_s_im_1,2.)
   +pow(R_s_re_1,2.))
   *(1./(2.*_PI_*_PI_));


 ppw->pvecmetric[ppw->index_mt_pk_s_2] =
   pow(k_2,3.)*(
   pow(R_s_im_2,2.)
   +pow(R_s_re_2,2.))
   *(1./(2.*_PI_*_PI_));


 ppw->pvecmetric[ppw->index_mt_pk_s_3] =
   pow(k_3,3.)*(
   pow(R_s_im_3,2.)
   +pow(R_s_re_3,2.))
   *(1./(2.*_PI_*_PI_));




    //Store the integrand: products of the real and imaginary parts
    //of the fields
 //1
 ppw->pvecmetric[ppw->index_mt_phi_prefactor_123_000] = y[ppw->pv->index_pt_bispectrum_s_prefactor_123_000];
 ppw->pvecmetric[ppw->index_mt_phi_123_000] = y[ppw->pv->index_pt_bispectrum_s_123_000];
 ppw->pvecmetric[ppw->index_mt_phi_123_100] = y[ppw->pv->index_pt_bispectrum_s_123_100];
 ppw->pvecmetric[ppw->index_mt_phi_123_010] = y[ppw->pv->index_pt_bispectrum_s_123_010];
 ppw->pvecmetric[ppw->index_mt_phi_123_001] = y[ppw->pv->index_pt_bispectrum_s_123_001];
 ppw->pvecmetric[ppw->index_mt_phi_123_110] = y[ppw->pv->index_pt_bispectrum_s_123_110];
 ppw->pvecmetric[ppw->index_mt_phi_123_101] = y[ppw->pv->index_pt_bispectrum_s_123_101];
 ppw->pvecmetric[ppw->index_mt_phi_123_011] = y[ppw->pv->index_pt_bispectrum_s_123_011];
 ppw->pvecmetric[ppw->index_mt_phi_123_111] = y[ppw->pv->index_pt_bispectrum_s_123_111];
  
 //2
 ppw->pvecmetric[ppw->index_mt_phi_123d_000] = y[ppw->pv->index_pt_bispectrum_s_123d_000];
 ppw->pvecmetric[ppw->index_mt_phi_123d_100] = y[ppw->pv->index_pt_bispectrum_s_123d_100];
 ppw->pvecmetric[ppw->index_mt_phi_123d_010] = y[ppw->pv->index_pt_bispectrum_s_123d_010];
 ppw->pvecmetric[ppw->index_mt_phi_123d_001] = y[ppw->pv->index_pt_bispectrum_s_123d_001];
 ppw->pvecmetric[ppw->index_mt_phi_123d_110] = y[ppw->pv->index_pt_bispectrum_s_123d_110];
 ppw->pvecmetric[ppw->index_mt_phi_123d_101] = y[ppw->pv->index_pt_bispectrum_s_123d_101];
 ppw->pvecmetric[ppw->index_mt_phi_123d_011] = y[ppw->pv->index_pt_bispectrum_s_123d_011];
 ppw->pvecmetric[ppw->index_mt_phi_123d_111] = y[ppw->pv->index_pt_bispectrum_s_123d_111];
  
 //3
 ppw->pvecmetric[ppw->index_mt_phi_12d3_000] = y[ppw->pv->index_pt_bispectrum_s_12d3_000];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_100] = y[ppw->pv->index_pt_bispectrum_s_12d3_100];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_010] = y[ppw->pv->index_pt_bispectrum_s_12d3_010];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_001] = y[ppw->pv->index_pt_bispectrum_s_12d3_001];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_110] = y[ppw->pv->index_pt_bispectrum_s_12d3_110];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_101] = y[ppw->pv->index_pt_bispectrum_s_12d3_101];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_011] = y[ppw->pv->index_pt_bispectrum_s_12d3_011];
 ppw->pvecmetric[ppw->index_mt_phi_12d3_111] = y[ppw->pv->index_pt_bispectrum_s_12d3_111];
  
 //4
 ppw->pvecmetric[ppw->index_mt_phi_1d23_000] = y[ppw->pv->index_pt_bispectrum_s_1d23_000];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_100] = y[ppw->pv->index_pt_bispectrum_s_1d23_100];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_010] = y[ppw->pv->index_pt_bispectrum_s_1d23_010];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_001] = y[ppw->pv->index_pt_bispectrum_s_1d23_001];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_110] = y[ppw->pv->index_pt_bispectrum_s_1d23_110];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_101] = y[ppw->pv->index_pt_bispectrum_s_1d23_101];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_011] = y[ppw->pv->index_pt_bispectrum_s_1d23_011];
 ppw->pvecmetric[ppw->index_mt_phi_1d23_111] = y[ppw->pv->index_pt_bispectrum_s_1d23_111];
  
 //5
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_000] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_000];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_100] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_100];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_010] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_010];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_001] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_001];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_110] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_110];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_101] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_101];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_011] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_011];
 ppw->pvecmetric[ppw->index_mt_phi_1d2d3_111] = y[ppw->pv->index_pt_bispectrum_s_1d2d3_111];
  
 //6
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_000] = y[ppw->pv->index_pt_bispectrum_s_1d23d_000];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_100] = y[ppw->pv->index_pt_bispectrum_s_1d23d_100];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_010] = y[ppw->pv->index_pt_bispectrum_s_1d23d_010];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_001] = y[ppw->pv->index_pt_bispectrum_s_1d23d_001];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_110] = y[ppw->pv->index_pt_bispectrum_s_1d23d_110];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_101] = y[ppw->pv->index_pt_bispectrum_s_1d23d_101];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_011] = y[ppw->pv->index_pt_bispectrum_s_1d23d_011];
 ppw->pvecmetric[ppw->index_mt_phi_1d23d_111] = y[ppw->pv->index_pt_bispectrum_s_1d23d_111];
  
 //7
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_000] = y[ppw->pv->index_pt_bispectrum_s_12d3d_000];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_100] = y[ppw->pv->index_pt_bispectrum_s_12d3d_100];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_010] = y[ppw->pv->index_pt_bispectrum_s_12d3d_010];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_001] = y[ppw->pv->index_pt_bispectrum_s_12d3d_001];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_110] = y[ppw->pv->index_pt_bispectrum_s_12d3d_110];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_101] = y[ppw->pv->index_pt_bispectrum_s_12d3d_101];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_011] = y[ppw->pv->index_pt_bispectrum_s_12d3d_011];
 ppw->pvecmetric[ppw->index_mt_phi_12d3d_111] = y[ppw->pv->index_pt_bispectrum_s_12d3d_111];
  
 
  }



  if (_tensors_lqc_) {
    prefactor =  prefactor_bispectrum_t;

    
    //pvecmetric always contains the cosmic time quantities


 ppw->pvecmetric[ppw->index_mt_gw_re] = y[ppw->pv->index_pt_gw_re];
 ppw->pvecmetric[ppw->index_mt_gw_im] = y[ppw->pv->index_pt_gw_im];

 ppw->pvecmetric[ppw->index_mt_gwdot_re] = y[ppw->pv->index_pt_gwdot_re];
 ppw->pvecmetric[ppw->index_mt_gwdot_im] = y[ppw->pv->index_pt_gwdot_im];


    
ppw->pvecmetric[ppw->index_mt_gw_prime_prime_re] =
  -3.*a_prime_over_a
  *ppw->pvecmetric[ppw->index_mt_gwdot_re]/a
  -(k2
    -a_prime_over_a
    *a_prime_over_a
    *ppw->pvecback[pba->index_bg_f_T_Q])
  *ppw->pvecmetric[ppw->index_mt_gw_re]/a/a;

ppw->pvecmetric[ppw->index_mt_gw_prime_prime_im] =
  -3.*a_prime_over_a
  *ppw->pvecmetric[ppw->index_mt_gwdot_im]/a
  -(k2
    -a_prime_over_a
    *a_prime_over_a
    *ppw->pvecback[pba->index_bg_f_T_Q])
  *ppw->pvecmetric[ppw->index_mt_gw_im]/a/a;


if(ppw->approx == (int)efolds){
ppw->pvecmetric[ppw->index_mt_gwdot_re] *=
  (ppw->pvecback[pba->index_bg_H]);
ppw->pvecmetric[ppw->index_mt_gwdot_im] *=
  (ppw->pvecback[pba->index_bg_H]);

 prefactor =   prefactor/fabs(ppw->pvecback[pba->index_bg_H]);
  
    ppw->pvecmetric[ppw->index_mt_gw_prime_prime_re] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_gwdot_re]
      -K*K*y[ppw->pv->index_pt_gw_re];
 
    ppw->pvecmetric[ppw->index_mt_gw_prime_prime_im] =
      -(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_gwdot_im]
      -K*K*y[ppw->pv->index_pt_gw_im];

    if(ppw->pvecback[pba->index_bg_H]<0.){
      
  ppw->pvecmetric[ppw->index_mt_gw_prime_prime_re] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_gwdot_re]
      -K*K*y[ppw->pv->index_pt_gw_re];
 
    ppw->pvecmetric[ppw->index_mt_gw_prime_prime_im] =
      +(3.-ppw->pvecback[pba->index_bg_epsilon_H])
      *y[ppw->pv->index_pt_gwdot_im]
      -K*K*y[ppw->pv->index_pt_gw_im];

     }
    }

 ppw->pvecmetric[ppw->index_mt_gw_abs] =
   pow(pow(y[ppw->pv->index_pt_gw_re],2.)
   +pow(y[ppw->pv->index_pt_gw_im],2.),1./2.);

 
 ppw->pvecmetric[ppw->index_mt_pk_t] =
   (8./(2.*_PI_*_PI_))
   *a2*pow(k,3.)
   *pow(ppw->pvecmetric[ppw->index_mt_gw_abs],2.)
   /pow(ppw->pvecback[pba->index_bg_z_T],2.);



    //Store the integrand: products of the real and imaginary parts
    //of the fields
 
 ppw->pvecmetric[ppw->index_mt_gw_000] = y[ppw->pv->index_pt_bispectrum_t_000];
 ppw->pvecmetric[ppw->index_mt_gw_100] = y[ppw->pv->index_pt_bispectrum_t_100];
 ppw->pvecmetric[ppw->index_mt_gw_110] = y[ppw->pv->index_pt_bispectrum_t_110];
 ppw->pvecmetric[ppw->index_mt_gw_111] = y[ppw->pv->index_pt_bispectrum_t_111];
  
   ppw->pvecmetric[ppw->index_mt_v_t_im] =
     a*ppw->pvecmetric[ppw->index_mt_gw_im];
   
   ppw->pvecmetric[ppw->index_mt_v_t_re] =
     a*ppw->pvecmetric[ppw->index_mt_gw_re];

   

   ppw->pvecmetric[ppw->index_mt_v_t_prime_im] =
          a*ppw->pvecmetric[ppw->index_mt_gwdot_im]
          +H*ppw->pvecmetric[ppw->index_mt_v_t_im];
   
   ppw->pvecmetric[ppw->index_mt_v_t_prime_re] =
     a*ppw->pvecmetric[ppw->index_mt_gwdot_re]
          +H*ppw->pvecmetric[ppw->index_mt_v_t_re];



ppw->pvecmetric[ppw->index_mt_bispectrum_t_111] =
     prefactor*pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_re],0.)
      *pow(y[ppw->pv->index_pt_gw_re],0.)
      *pow(y[ppw->pv->index_pt_gw_re],0.);
 
ppw->pvecmetric[ppw->index_mt_bispectrum_t_110]=
      prefactor*pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_re],1.)
      *pow(y[ppw->pv->index_pt_gw_re],0.)
      *pow(y[ppw->pv->index_pt_gw_re],0.);      


ppw->pvecmetric[ppw->index_mt_bispectrum_t_100]= 
      prefactor*pow(y[ppw->pv->index_pt_gw_im],1.)
      *pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_re],1.)
      *pow(y[ppw->pv->index_pt_gw_re],1.)
        *pow(y[ppw->pv->index_pt_gw_re],0.);

ppw->pvecmetric[ppw->index_mt_bispectrum_t_000]=
      prefactor*pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_im],0.)
      *pow(y[ppw->pv->index_pt_gw_re],1.)
      *pow(y[ppw->pv->index_pt_gw_re],1.)
      *pow(y[ppw->pv->index_pt_gw_re],1.);
 
  }


    
  return _SUCCESS_;

}


/**
 * Compute the source functions 
 */
 
int perturb_lqc_sources(
                    double tau,
                    double * y,
                    double * dy,
                    int index_tau,
                    void * parameters_and_workspace,
                    ErrorMsg error_message
                    ) {
  return _SUCCESS_;

}




int perturb_lqc_print_variables(double tau,
                            double * y,
                            double * dy,
                            void * parameters_and_workspace,
                            ErrorMsg error_message
                            ) {

  struct perturb_lqc_parameters_and_workspace * pppaw;
  /** Summary: */

  /** - define local variables */
  double k;
  int index_md;
  //struct precision * ppr;
  struct background_lqc * pba;
  struct perturbs_lqc * ppt;
  struct perturb_lqc_workspace * ppw;
  double * pvecback;
  double * pvecmetric;

  double a,a2,H,N;
  int idx,index_q, storeidx;
  double *dataptr;


  /** - rename structure fields (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;
  k = pppaw->k;
  index_md = pppaw->index_md;
  //ppr = pppaw->ppr;
  pba = pppaw->pba;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;
  pvecback = ppw->pvecback;
  pvecmetric = ppw->pvecmetric;

  a = pvecback[pba->index_bg_a];
  a2 = a*a;
  H = pvecback[pba->index_bg_H];
  N = pvecback[pba->index_bg_N];

  double k_1 = k;
  double k_2 = ppt->lambda_k_2*k_1;
  double k_3 = ppt->lambda_k_3*k_1;
  double Pphi = pvecback[pba->index_bg_Pphi_scf];
  double PIa =  pvecback[pba->index_bg_pia];
  double V_phi =  pvecback[pba->index_bg_V_phi];
  double V_phiphi =  pvecback[pba->index_bg_V_phiphi];
  double V_phiphiphi =  pvecback[pba->index_bg_V_phiphiphi];
  double k_1_dot_k_2 = (k_3*k_3-k_1*k_1-k_2*k_2)/2.;
  double k_1_dot_k_3 = (k_2*k_2-k_1*k_1-k_3*k_3)/2.;
  double k_2_dot_k_3 = (k_1*k_1-k_2*k_2-k_3*k_3)/2.;
  double h12 = k_1_dot_k_2/k_1/k_1;
  double h21 = k_1_dot_k_2/k_2/k_2;
  double h13 = k_1_dot_k_3/k_1/k_1;
  double h31 = k_1_dot_k_3/k_3/k_3;
  double h23 = k_2_dot_k_3/k_2/k_2;
  double h32 = k_2_dot_k_3/k_3/k_3;
  double h12_s = (h12+h21)/2.;
  double h13_s = (h13+h31)/2.;
  double h23_s = (h23+h32)/2.;
  double t123 = (k_1*k_1+k_2*k_2+k_3*k_3)/(a*a);
  double g23 = k_2_dot_k_3*k_2_dot_k_3/k_2/k_2/k_3/k_3-1.;
  double g13 = k_1_dot_k_3*k_1_dot_k_3/k_1/k_1/k_3/k_3-1.;
  double g12 = k_1_dot_k_2*k_1_dot_k_2/k_2/k_2/k_1/k_1-1.;
  double  cut_off_bispectrum = exp(-k_1/a/sqrt(pvecback[pba->index_bg_rho_scf]/3.)/ppt->delta_cut_off_bispectrum);






/** - for scalar modes */



  
  if (_scalars_lqc_) {



    //    fprintf(ppw->perturb_lqc_output_file," ");
    /** - --> Handle (re-)allocation */
    if (ppt->scalar_perturbations_lqc_data[ppw->index_ikout] == NULL){
      class_alloc(ppt->scalar_perturbations_lqc_data[ppw->index_ikout],
                  sizeof(double)*ppt->number_of_scalar_titles,
                  error_message);
      ppt->size_scalar_perturbation_lqc_data[ppw->index_ikout] = 0;
    }
    else{
      ppt->scalar_perturbations_lqc_data[ppw->index_ikout] =
        realloc(ppt->scalar_perturbations_lqc_data[ppw->index_ikout],
                sizeof(double)*(ppt->size_scalar_perturbation_lqc_data[ppw->index_ikout]+ppt->number_of_scalar_titles));
    }
    storeidx = 0;
    dataptr = ppt->scalar_perturbations_lqc_data[ppw->index_ikout]+
      ppt->size_scalar_perturbation_lqc_data[ppw->index_ikout];
    ppt->size_scalar_perturbation_lqc_data[ppw->index_ikout] += ppt->number_of_scalar_titles;


    class_store_double(dataptr, pvecback[pba->index_bg_cosmic_time], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[pba->index_bg_H2_f_S_Q], _TRUE_, storeidx);//This was 'a'.
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_scf_re], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_scf_im], _TRUE_, storeidx);
    class_store_double(dataptr,
                       pvecback[pba->index_bg_a]
                       *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re], _TRUE_, storeidx);
    class_store_double(dataptr,
                       pvecback[pba->index_bg_a]
                       *ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im], _TRUE_, storeidx);
    class_store_double(dataptr, /*pow(k,3.)*(
   pow(ppw->pvecmetric[ppw->index_mt_phi_scf_im],2.)
   +pow(ppw->pvecmetric[ppw->index_mt_phi_scf_re],2.))
   *(1./(2.*_PI_*_PI_))*/
                       
    pow(k,3.)*(
   pow(ppw->pvecmetric[ppw->index_mt_R_s_im],2.)
   +pow(ppw->pvecmetric[ppw->index_mt_R_s_re],2.))
   *(1./(2.*_PI_*_PI_))
                       
, _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_s_re], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_s_im], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_s_prime_re], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_s_prime_im], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[pba->index_bg_H], _TRUE_, storeidx);

    
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123_111], _TRUE_, storeidx);

     
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123d_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123d_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123d_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_123d_111], _TRUE_, storeidx);

    
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3d_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3d_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3d_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3d_111], _TRUE_, storeidx);


    //printf("IN FILE: Y_100=%e.\n",ppw->pvecmetric[ppw->index_mt_phi_12d3d_110]);



        double Y_000,Y_100,Y_010,Y_001,Y_110,Y_101,Y_011,Y_111;

 
 double Q_s_1_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_1];
 double Q_s_1d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_1];
 
 double Q_s_2_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_2];
 double Q_s_2d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_2];

 double Q_s_3_im = ppw->pvecmetric[ppw->index_mt_phi_scf_im_3];
 double Q_s_3d_im = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_im_3];
 
 double Q_s_1_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_1];
 double Q_s_1d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_1];

 double Q_s_2_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_2];
 double Q_s_2d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_2];

 double Q_s_3_re = ppw->pvecmetric[ppw->index_mt_phi_scf_re_3];
 double Q_s_3d_re = ppw->pvecmetric[ppw->index_mt_phi_scf_prime_re_3];

 //H over phi_dot
 double a_over_z_3 =
     pow(ppw->pvecback[pba->index_bg_H]
         /ppw->pvecback[pba->index_bg_phi_prime_scf],3.);
 double a_over_z =
     pow(ppw->pvecback[pba->index_bg_H]
         /ppw->pvecback[pba->index_bg_phi_prime_scf],1.);





Y_111 = 
       pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,1.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,0.);
 

Y_110 =
       pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,1.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,0.)
       *pow(Q_s_3_re,1.);


Y_100 =
       pow(Q_s_1_im,1.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,0.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);

     
Y_000 =
       pow(Q_s_1_im,0.)
       *pow(Q_s_2_im,0.)
       *pow(Q_s_3_im,0.)
       *pow(Q_s_1_re,1.)
       *pow(Q_s_2_re,1.)
       *pow(Q_s_3_re,1.);




    class_store_double(dataptr,Y_000, _TRUE_, storeidx);
    class_store_double(dataptr,Y_100, _TRUE_, storeidx);
    class_store_double(dataptr,Y_110, _TRUE_, storeidx);
    class_store_double(dataptr,Y_111, _TRUE_, storeidx);

double field_redef = -ppw->pvecmetric[ppw->index_mt_bispectrum_s_field_redef];

    class_store_double(dataptr,field_redef, _TRUE_, storeidx);

    

    class_store_double(dataptr, pvecback[pba->index_bg_N], _TRUE_, storeidx);
    //fprintf(ppw->perturb_lqc_output_file,"\n");



    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_12d3_111], _TRUE_, storeidx);


    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23_111], _TRUE_, storeidx);

    
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23d_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23d_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23d_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d23d_111], _TRUE_, storeidx);


    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d2d3_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d2d3_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d2d3_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_1d2d3_111], _TRUE_, storeidx);

    /*
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_phi_prefactor_123_000], _TRUE_, storeidx);
    */

    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_bispectrum_s_prefactor_123_000], _TRUE_, storeidx);


    class_store_double(dataptr, a_over_z, _TRUE_, storeidx);

  double prefactor_123 = -(-( (243.*pow(Pphi,7.)*pow(a,-8.)*pow(PIa,-5.)/2.
	            - 81.*pow(Pphi,5.)*pow(a*a*PIa,-3.)/2. + 27.*pow(Pphi,3.)*pow(a,-4.)/PIa/8.
 	            + 81.*pow(Pphi/PIa,4.)*V_phi/a - 27.*a*pow(Pphi/PIa,2.)*V_phi/2 + 27.*pow(a*a/PIa,3.)*Pphi*V_phi*V_phi/2.)
	            *2.*(-g12 -g23 -g13)
	             + 3.*a*a*Pphi/PIa/2.*(-t123) + 6.*9.*a*pow(Pphi/PIa,2.)*V_phi/2.  //010717
	             - 6.*3.*a*a*Pphi*V_phiphi/PIa/2. + 6.*a*a*a*V_phiphiphi/6.)
	              *cut_off_bispectrum);

  double prefactor_1d23 = -( -a*a*a*(
		        ( 81.*pow(Pphi,5.)*pow(a,-7.)*pow(PIa,-4.) - 27.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
		          + 27.*pow(Pphi,2.)*V_phi*pow(PIa,-3.))*(-g12 -g13) - 2.*9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2.)/2.
			+ ( -3.*Pphi*pow(a,-3.)/2. + 9.*pow(Pphi,3.)*pow(a,-5.)*pow(PIa,-2) + 3.*a*a*V_phi/PIa )*h23_s)
		        *cut_off_bispectrum);

  double prefactor_1d2d3 = -( -pow(a,6.)*(
		           27.*pow(Pphi,3.)*pow(a*a*PIa,-3)*(-g12) - 3.*Pphi*pow(a,-4.)/PIa
		           + 3.*Pphi*pow(a,-4.)/PIa*(h13 + h23))
		           *cut_off_bispectrum);

  class_store_double(dataptr, prefactor_123, _TRUE_, storeidx);
  class_store_double(dataptr, prefactor_1d23, _TRUE_, storeidx);
  class_store_double(dataptr, prefactor_1d2d3, _TRUE_, storeidx);
    
  }
  /** - for tensor modes: */

  if (_tensors_lqc_) {

    /** - --> Handle (re-)allocation */
    if (ppt->tensor_perturbations_lqc_data[ppw->index_ikout] == NULL){
      class_alloc(ppt->tensor_perturbations_lqc_data[ppw->index_ikout],
                  sizeof(double)*ppt->number_of_tensor_titles,
                  error_message);
      ppt->size_tensor_perturbation_lqc_data[ppw->index_ikout] = 0;
    }
    else{
      ppt->tensor_perturbations_lqc_data[ppw->index_ikout] =
        realloc(ppt->tensor_perturbations_lqc_data[ppw->index_ikout],
                sizeof(double)*(ppt->size_tensor_perturbation_lqc_data[ppw->index_ikout]+ppt->number_of_tensor_titles));
    }
    storeidx = 0;
    dataptr = ppt->tensor_perturbations_lqc_data[ppw->index_ikout]+
      ppt->size_tensor_perturbation_lqc_data[ppw->index_ikout];
    ppt->size_tensor_perturbation_lqc_data[ppw->index_ikout] += ppt->number_of_tensor_titles;

    class_store_double(dataptr, pvecback[pba->index_bg_cosmic_time], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[pba->index_bg_a], _TRUE_, storeidx);
    class_store_double(dataptr,  ppw->pvecmetric[ppw->index_mt_gw_re], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_gw_im], _TRUE_, storeidx);
    class_store_double(dataptr,
                       pvecback[pba->index_bg_a]
                       *ppw->pvecmetric[ppw->index_mt_gwdot_re], _TRUE_, storeidx);
    class_store_double(dataptr,
                       pvecback[pba->index_bg_a]
                       *ppw->pvecmetric[ppw->index_mt_gwdot_im], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_pk_t], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_t_re], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_t_im], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_t_prime_re], _TRUE_, storeidx);
    class_store_double(dataptr,ppw->pvecmetric[ppw->index_mt_v_t_prime_im], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[pba->index_bg_H], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_gw_000], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_gw_100], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_gw_110], _TRUE_, storeidx);
    class_store_double(dataptr, ppw->pvecmetric[ppw->index_mt_gw_111], _TRUE_, storeidx);
    class_store_double(dataptr, pvecback[pba->index_bg_N], _TRUE_, storeidx);

  }



    
  return _SUCCESS_;

  }



int perturb_lqc_derivs(double tau,
                   double * y,
                   double * dy,
                   void * parameters_and_workspace,
                   ErrorMsg error_message
                   ) {
  
  /* scale factor and other background quantities */
  double a,a2,a_prime_over_a,R;

  /* short-cut names for the fields of the input structure */
  struct perturb_lqc_parameters_and_workspace * pppaw;
  double k,k2,k_2,k_3;
  int index_md;
  struct precision * ppr;
  struct background_lqc * pba;
  struct perturbs_lqc * ppt;
  struct perturb_lqc_workspace * ppw;
  double * pvecback;
  double * pvecmetric;
  struct perturb_lqc_vector * pv;


  /** - rename the fields of the input structure (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;

  k = pppaw->k;
  k_2 =  pppaw->k_2;
  k_3 =  pppaw->k_3;
  int index_k = pppaw->index_k;
  
  k2=k*k;
  index_md = pppaw->index_md;
  ppr = pppaw->ppr;
  pba = pppaw->pba;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;

  pvecback = ppw->pvecback;
  pvecmetric = ppw->pvecmetric;
  pv = ppw->pv;

  /** - get background/thermo quantities in this point */

  
  //if(ppw->approx == (int)efolds)
      //printf("Mode k=%e, efolds.\n",k);
    // else printf("Mode k=%e, cosmic time.\n",k);
  
  if(ppw->approx == (int)efolds){
    double tau_of_N;
if(ppw->num_approx_intervals > 3 && ppw->approx_interval < 1){
class_call(background_lqc_tau_of_N_pb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);

 
    class_call(background_lqc_at_tau_pb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
 
/*
  class_call(background_lqc_at_N(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_closeby,
                               &(ppw->last_index_back),
                               pvecback),
             pba->error_message,
             error_message);
*/
class_call(background_lqc_N_of_tau_pb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);
    
}
else{
class_call(background_lqc_tau_of_N_fb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);

 
    class_call(background_lqc_at_tau_fb(pba,
                                 tau_of_N,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
               pba->error_message,
               ppt->error_message);
 
/*
  class_call(background_lqc_at_N(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_closeby,
                               &(ppw->last_index_back),
                               pvecback),
             pba->error_message,
             error_message);
*/
class_call(background_lqc_N_of_tau_fb(pba,
                               tau,
                               &tau_of_N),
                   pba->error_message,
                   ppt->error_message);


}
  
 class_call(perturb_lqc_einstein(ppr,
                              pba,
                              ppt,
                              index_md,
                              k,
                              k_2,
                              k_3,
                              tau_of_N,
                              y,
                              ppw),
             ppt->error_message,
             error_message);


 }
 else {  
  if(ppw->num_approx_intervals > 2 && ppw->approx_interval < 1){
  class_call(background_lqc_at_tau_pb(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_closeby,
                               &(ppw->last_index_back),
                               pvecback),
             pba->error_message,
             error_message);
  }
  else{
  class_call(background_lqc_at_tau_fb(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_closeby,
                               &(ppw->last_index_back),
                               pvecback),
             pba->error_message,
             error_message);

  }

  /** - get metric perturbations_lqc with perturb_lqc_einstein() */
  class_call(perturb_lqc_einstein(ppr,
                              pba,
                              ppt,
                              index_md,
                              k,
                              k_2,
                              k_3,
                              tau,
                              y,
                              ppw),
             ppt->error_message,
             error_message);
}
  /** - compute related background quantities */

  a = pvecback[pba->index_bg_a];
  a2 = a*a;
  a_prime_over_a = pvecback[pba->index_bg_H] * a;

// printf("derivs at t = %e...future_branch = %d\n",tau, pba->future_branch);
//if(k>10.)printf("%e %e \n",pvecback[pba->index_bg_cosmic_time],pvecback[pba->index_bg_zS_primeprime_over_zS]);
  
  if (_scalars_lqc_) {
 
      dy[pv->index_pt_phi_scf_re] =
        y[ppw->pv->index_pt_phi_prime_scf_re];
      dy[pv->index_pt_phi_scf_im] =
        y[ppw->pv->index_pt_phi_prime_scf_im];
      dy[pv->index_pt_phi_prime_scf_re] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_re];
      dy[pv->index_pt_phi_prime_scf_im] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_im];
   

      
      
 
      dy[pv->index_pt_phi_scf_re_1] =
        y[ppw->pv->index_pt_phi_prime_scf_re_1];
      dy[pv->index_pt_phi_scf_im_1] =
        y[ppw->pv->index_pt_phi_prime_scf_im_1];
      dy[pv->index_pt_phi_prime_scf_re_1] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_1];
      dy[pv->index_pt_phi_prime_scf_im_1] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_1];
   

      
      
 
      dy[pv->index_pt_phi_scf_re_2] =
        y[ppw->pv->index_pt_phi_prime_scf_re_2];
      dy[pv->index_pt_phi_scf_im_2] =
        y[ppw->pv->index_pt_phi_prime_scf_im_2];
      dy[pv->index_pt_phi_prime_scf_re_2] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_2];
      dy[pv->index_pt_phi_prime_scf_im_2] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_2];
   

      
      
 
      dy[pv->index_pt_phi_scf_re_3] =
        y[ppw->pv->index_pt_phi_prime_scf_re_3];
      dy[pv->index_pt_phi_scf_im_3] =
        y[ppw->pv->index_pt_phi_prime_scf_im_3];
      dy[pv->index_pt_phi_prime_scf_re_3] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_re_3];
      dy[pv->index_pt_phi_prime_scf_im_3] =
        pvecmetric[ppw->index_mt_phi_scf_prime_prime_im_3];
   

      
      
                if (ppt->has_bispectrum_lqc == _TRUE_) {

       if (ppt->tensors_only == _FALSE_) {
         //1

     dy[pv->index_pt_bispectrum_s_prefactor_123_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_prefactor_123_000];
         
     dy[pv->index_pt_bispectrum_s_123_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_000];
     dy[pv->index_pt_bispectrum_s_123_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_100];
      dy[pv->index_pt_bispectrum_s_123_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_010];
      dy[pv->index_pt_bispectrum_s_123_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_001];
     dy[pv->index_pt_bispectrum_s_123_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_110];
     dy[pv->index_pt_bispectrum_s_123_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_101];
       dy[pv->index_pt_bispectrum_s_123_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_011];
       dy[pv->index_pt_bispectrum_s_123_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123_111];
       /* dy[pv->index_pt_bispectrum_s_123_000] =0.;
       dy[pv->index_pt_bispectrum_s_123_100] =0.;
       dy[pv->index_pt_bispectrum_s_123_010] =0.;
       dy[pv->index_pt_bispectrum_s_123_001] =0.;
       dy[pv->index_pt_bispectrum_s_123_110] =0.;
       dy[pv->index_pt_bispectrum_s_123_101] =0.;
       dy[pv->index_pt_bispectrum_s_123_011] =0.;
       dy[pv->index_pt_bispectrum_s_123_111] =0.;*/         
       //2
     dy[pv->index_pt_bispectrum_s_123d_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_000];
     dy[pv->index_pt_bispectrum_s_123d_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_100];
     dy[pv->index_pt_bispectrum_s_123d_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_010];
     dy[pv->index_pt_bispectrum_s_123d_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_001];
     dy[pv->index_pt_bispectrum_s_123d_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_110];
     dy[pv->index_pt_bispectrum_s_123d_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_101];
     dy[pv->index_pt_bispectrum_s_123d_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_011];
     dy[pv->index_pt_bispectrum_s_123d_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_123d_111];
            
            
     //3
     dy[pv->index_pt_bispectrum_s_12d3d_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_000];
     dy[pv->index_pt_bispectrum_s_12d3d_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_100];
     dy[pv->index_pt_bispectrum_s_12d3d_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_010];
     dy[pv->index_pt_bispectrum_s_12d3d_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_001];
     dy[pv->index_pt_bispectrum_s_12d3d_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_110];
     dy[pv->index_pt_bispectrum_s_12d3d_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_101];
     dy[pv->index_pt_bispectrum_s_12d3d_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_011];
     dy[pv->index_pt_bispectrum_s_12d3d_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3d_111];
            
     //4
     dy[pv->index_pt_bispectrum_s_12d3_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_000];
     dy[pv->index_pt_bispectrum_s_12d3_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_100];
     dy[pv->index_pt_bispectrum_s_12d3_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_010];
     dy[pv->index_pt_bispectrum_s_12d3_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_001];
     dy[pv->index_pt_bispectrum_s_12d3_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_110];
     dy[pv->index_pt_bispectrum_s_12d3_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_101];
     dy[pv->index_pt_bispectrum_s_12d3_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_011];
     dy[pv->index_pt_bispectrum_s_12d3_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_12d3_111];
            
     //5
     dy[pv->index_pt_bispectrum_s_1d2d3_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_000];
     dy[pv->index_pt_bispectrum_s_1d2d3_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_100];
     dy[pv->index_pt_bispectrum_s_1d2d3_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_010];
       dy[pv->index_pt_bispectrum_s_1d2d3_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_001];
       dy[pv->index_pt_bispectrum_s_1d2d3_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_110];
       dy[pv->index_pt_bispectrum_s_1d2d3_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_101];
       dy[pv->index_pt_bispectrum_s_1d2d3_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_011];
     dy[pv->index_pt_bispectrum_s_1d2d3_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d2d3_111];
            

     //6
     dy[pv->index_pt_bispectrum_s_1d23_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_000];
     dy[pv->index_pt_bispectrum_s_1d23_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_100];
     dy[pv->index_pt_bispectrum_s_1d23_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_010];
     dy[pv->index_pt_bispectrum_s_1d23_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_001];
     dy[pv->index_pt_bispectrum_s_1d23_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_110];
     dy[pv->index_pt_bispectrum_s_1d23_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_101];
     dy[pv->index_pt_bispectrum_s_1d23_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_011];
     dy[pv->index_pt_bispectrum_s_1d23_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23_111];
            


     //7
     dy[pv->index_pt_bispectrum_s_1d23d_000] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_000];
     dy[pv->index_pt_bispectrum_s_1d23d_100] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_100];
      dy[pv->index_pt_bispectrum_s_1d23d_010] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_010];
      dy[pv->index_pt_bispectrum_s_1d23d_001] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_001];
     dy[pv->index_pt_bispectrum_s_1d23d_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_110];
     dy[pv->index_pt_bispectrum_s_1d23d_101] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_101];
      dy[pv->index_pt_bispectrum_s_1d23d_011] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_011];
      dy[pv->index_pt_bispectrum_s_1d23d_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_s_1d23d_111];
            




     
                      }
       else {

       dy[pv->index_pt_bispectrum_s_prefactor_123_000] =0.;
       dy[pv->index_pt_bispectrum_s_123_000] =0.;
       dy[pv->index_pt_bispectrum_s_123_100] =0.;
       dy[pv->index_pt_bispectrum_s_123_010] =0.;
       dy[pv->index_pt_bispectrum_s_123_001] =0.;
       dy[pv->index_pt_bispectrum_s_123_110] =0.;
       dy[pv->index_pt_bispectrum_s_123_101] =0.;
       dy[pv->index_pt_bispectrum_s_123_011] =0.;
       dy[pv->index_pt_bispectrum_s_123_111] =0.;
       
       dy[pv->index_pt_bispectrum_s_123d_000] =0.;
       dy[pv->index_pt_bispectrum_s_123d_100] =0.;
       dy[pv->index_pt_bispectrum_s_123d_010] =0.;
       dy[pv->index_pt_bispectrum_s_123d_001] =0.;
       dy[pv->index_pt_bispectrum_s_123d_110] =0.;
       dy[pv->index_pt_bispectrum_s_123d_101] =0.;
       dy[pv->index_pt_bispectrum_s_123d_011] =0.;
       dy[pv->index_pt_bispectrum_s_123d_111] =0.;
       
       dy[pv->index_pt_bispectrum_s_12d3d_000] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_100] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_010] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_001] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_110] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_101] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_011] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_111] =0.;
      
       dy[pv->index_pt_bispectrum_s_12d3_000] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_100] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_010] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_001] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_110] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_101] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_011] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_111] =0.;

     
       dy[pv->index_pt_bispectrum_s_1d23_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_111] =0.;

       dy[pv->index_pt_bispectrum_s_1d23d_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_111] =0.;

 
     
       dy[pv->index_pt_bispectrum_s_1d2d3_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_111] =0.;
             
       
          }

 }
      else {

       dy[pv->index_pt_bispectrum_s_prefactor_123_000] =0.;
       dy[pv->index_pt_bispectrum_s_123_000] =0.;
       dy[pv->index_pt_bispectrum_s_123_100] =0.;
       dy[pv->index_pt_bispectrum_s_123_010] =0.;
       dy[pv->index_pt_bispectrum_s_123_001] =0.;
       dy[pv->index_pt_bispectrum_s_123_110] =0.;
       dy[pv->index_pt_bispectrum_s_123_101] =0.;
       dy[pv->index_pt_bispectrum_s_123_011] =0.;
       dy[pv->index_pt_bispectrum_s_123_111] =0.;
       
       dy[pv->index_pt_bispectrum_s_123d_000] =0.;
       dy[pv->index_pt_bispectrum_s_123d_100] =0.;
       dy[pv->index_pt_bispectrum_s_123d_010] =0.;
       dy[pv->index_pt_bispectrum_s_123d_001] =0.;
       dy[pv->index_pt_bispectrum_s_123d_110] =0.;
       dy[pv->index_pt_bispectrum_s_123d_101] =0.;
       dy[pv->index_pt_bispectrum_s_123d_011] =0.;
       dy[pv->index_pt_bispectrum_s_123d_111] =0.;
       
       dy[pv->index_pt_bispectrum_s_12d3d_000] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_100] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_010] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_001] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_110] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_101] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_011] =0.;
       dy[pv->index_pt_bispectrum_s_12d3d_111] =0.;
      
       dy[pv->index_pt_bispectrum_s_12d3_000] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_100] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_010] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_001] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_110] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_101] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_011] =0.;
       dy[pv->index_pt_bispectrum_s_12d3_111] =0.;

     
       dy[pv->index_pt_bispectrum_s_1d23_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d23_111] =0.;

       dy[pv->index_pt_bispectrum_s_1d23d_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d23d_111] =0.;

 
     
       dy[pv->index_pt_bispectrum_s_1d2d3_000] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_100] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_010] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_001] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_110] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_101] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_011] =0.;
       dy[pv->index_pt_bispectrum_s_1d2d3_111] =0.;
             
      }
  }
  

 if (_tensors_lqc_) {
    ppr->tol_perturb_lqc_integration= ppt->tol;

    dy[pv->index_pt_gw_re] =
      y[pv->index_pt_gwdot_re];
    dy[pv->index_pt_gw_im] =
      y[pv->index_pt_gwdot_im];
    dy[pv->index_pt_gwdot_re] = pvecmetric[ppw->index_mt_gw_prime_prime_re];
    dy[pv->index_pt_gwdot_im] = pvecmetric[ppw->index_mt_gw_prime_prime_im];


                if (ppt->has_bispectrum_lqc == _TRUE_) {

     dy[pv->index_pt_bispectrum_t_111] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_t_111];

     dy[pv->index_pt_bispectrum_t_110] =
       ppw->pvecmetric[ppw->index_mt_bispectrum_t_110];


      dy[pv->index_pt_bispectrum_t_100] =
        ppw->pvecmetric[ppw->index_mt_bispectrum_t_100];


      dy[pv->index_pt_bispectrum_t_000] =
        ppw->pvecmetric[ppw->index_mt_bispectrum_t_000];

      
  }
else {
  dy[pv->index_pt_bispectrum_t_111] =0.;
  dy[pv->index_pt_bispectrum_t_110] =0.;
  dy[pv->index_pt_bispectrum_t_100] =0.;
  dy[pv->index_pt_bispectrum_t_000] =0.;
           }       
  }

  return _SUCCESS_;
}

int spectra_pk_lqc(
               struct background_lqc * pba,
               struct perturbs_lqc * ppt,
               int index_md
               ) {

  /** Summary: */
    ppt->ln_k_size = ppt->k_size[index_md];

  class_alloc(ppt->ln_k,
              sizeof(double)*ppt->ln_k_size,
              ppt->error_message);


  /** - allocate and fill array of \f$P(k,\tau)\f$ values */

  class_alloc(ppt->ln_pk,
              sizeof(double)*ppt->ln_k_size,
              ppt->error_message);

  class_alloc(ppt->bispectrum_s_equi,
              sizeof(double)*ppt->ln_k_size,
              ppt->error_message);


  
  class_alloc(ppt->bispectrum_t_equi,
              sizeof(double)*ppt->ln_k_size,
              ppt->error_message);





  
  /** - define local variables */

    int index_tau;

    int index_type;
    int index_k;
    double k;

    double Y_000,Y_100,Y_010,Y_001,Y_110,Y_101,Y_011,Y_111;
    double Z_000,Z_100,Z_010,Z_001,Z_110,Z_101,Z_011,Z_111;

    double Z_123_000,Z_123_100,Z_123_010,Z_123_001,Z_123_110,Z_123_101,Z_123_011,Z_123_111;
    double Z_1d23_000,Z_1d23_100,Z_1d23_010,Z_1d23_001,Z_1d23_110,Z_1d23_101,Z_1d23_011,Z_1d23_111;
    double Z_12d3_000,Z_12d3_100,Z_12d3_010,Z_12d3_001,Z_12d3_110,Z_12d3_101,Z_12d3_011,Z_12d3_111;
    double Z_123d_000,Z_123d_100,Z_123d_010,Z_123d_001,Z_123d_110,Z_123d_101,Z_123d_011,Z_123d_111;
    double Z_12d3d_000,Z_12d3d_100,Z_12d3d_010,Z_12d3d_001,Z_12d3d_110,Z_12d3d_101,Z_12d3d_011,Z_12d3d_111;
    double Z_1d2d3_000,Z_1d2d3_100,Z_1d2d3_010,Z_1d2d3_001,Z_1d2d3_110,Z_1d2d3_101,Z_1d2d3_011,Z_1d2d3_111;
    double Z_1d23d_000,Z_1d23d_100,Z_1d23d_010,Z_1d23d_001,Z_1d23d_110,Z_1d23d_101,Z_1d23d_011,Z_1d23d_111;

 
      double integral;

      double integrand_tau_plus;
      double integrand_tau_minus;
      
      int ik;
      int counter_integrand_output;
      counter_integrand_output = 0;

      double K;
      double a_prime_over_a;
      double  cut_off_bispectrum;
      double prefactor_bispectrum_s;
      double prefactor_bispectrum_t;
      double prefactor;



    for (index_k=0; index_k<ppt->k_size[index_md]; index_k++) {


     k = ppt->k[index_md][index_k];
     
    ppt->index_k_output=0;
      
      integral = 0.;
      
      if ((ppt->has_bispectrum_lqc == _TRUE_)) {
        
        
if (index_md==ppt->index_md_scalars)  {

  //collect the 'monomial' integrands at the
  //end of integration
  
        Y_000 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_000]
               [index_k];
        Y_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_100]
               [index_k];
        Y_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_010]
               [index_k];
        Y_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_001]
               [index_k];
        Y_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_110]
               [index_k];
        Y_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_101]
               [index_k];
        Y_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_011]
               [index_k];
        Y_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_111]
               [index_k];


  //collect the 'monomial' integrands at the
  //end of integration
  
        
        Z_123_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_000]
               [index_k];
        double prefactor_123_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_prefactor_123_000]
               [index_k];
        Z_123_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_100]
               [index_k];
        Z_123_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_010]
               [index_k];
        Z_123_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_001]
               [index_k];
        Z_123_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_110]
               [index_k];
        Z_123_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_101]
               [index_k];
        Z_123_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_011]
               [index_k];
        Z_123_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123_111]
               [index_k];
        


        
        Z_1d23_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_000]
               [index_k];
        Z_1d23_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_100]
               [index_k];
        Z_1d23_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_010]
               [index_k];
        Z_1d23_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_001]
               [index_k];
        Z_1d23_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_110]
               [index_k];
        Z_1d23_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_101]
               [index_k];
        Z_1d23_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_011]
               [index_k];
        Z_1d23_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23_111]
               [index_k];
        

        
        Z_12d3_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_000]
               [index_k];
        Z_12d3_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_100]
               [index_k];
        Z_12d3_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_010]
               [index_k];
        Z_12d3_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_001]
               [index_k];
        Z_12d3_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_110]
               [index_k];
        Z_12d3_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_101]
               [index_k];
        Z_12d3_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_011]
               [index_k];
        Z_12d3_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3_111]
               [index_k];
        

        
        Z_123d_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_000]
               [index_k];
        Z_123d_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_100]
               [index_k];
        Z_123d_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_010]
               [index_k];
        Z_123d_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_001]
               [index_k];
        Z_123d_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_110]
               [index_k];
        Z_123d_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_101]
               [index_k];
        Z_123d_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_011]
               [index_k];
        Z_123d_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_123d_111]
               [index_k];
        



        
        Z_1d2d3_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_000]
               [index_k];
        Z_1d2d3_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_100]
               [index_k];
        Z_1d2d3_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_010]
               [index_k];
        Z_1d2d3_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_001]
               [index_k];
        Z_1d2d3_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_110]
               [index_k];
        Z_1d2d3_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_101]
               [index_k];
        Z_1d2d3_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_011]
               [index_k];
        Z_1d2d3_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d2d3_111]
               [index_k];
        


        
        Z_1d23d_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_000]
               [index_k];
        Z_1d23d_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_100]
               [index_k];
        Z_1d23d_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_010]
               [index_k];
        Z_1d23d_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_001]
               [index_k];
        Z_1d23d_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_110]
               [index_k];
        Z_1d23d_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_101]
               [index_k];
        Z_1d23d_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_011]
               [index_k];
        Z_1d23d_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_1d23d_111]
               [index_k];
        


        
        Z_12d3d_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_000]
               [index_k];
        Z_12d3d_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_100]
               [index_k];
        Z_12d3d_010 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_010]
               [index_k];
        Z_12d3d_001 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_001]
               [index_k];
        Z_12d3d_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_110]
               [index_k];
        Z_12d3d_101 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_101]
               [index_k];
        Z_12d3d_011 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_011]
               [index_k];
        Z_12d3d_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_12d3d_111]
               [index_k];
        



        
                                      }


if (index_md==ppt->index_md_tensors)  {

  //collect the values of the fields at the
  //end of integration
  
        Y_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_tau_end_000]
               [index_k];
        Y_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_tau_end_100]
               [index_k];
        Y_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_tau_end_110]
               [index_k];
        Y_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_tau_end_111]
               [index_k];

  //collect the 'monomial' integrands at the
  //end of integration
  
        Z_000 = ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_000]
               [index_k];
        Z_100 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_100]
               [index_k];
        Z_110 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_110]
               [index_k];
        Z_111 =  ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_t_equi_111]
               [index_k];

        

                              }

//Interaction term:
      if (index_md==ppt->index_md_tensors)
        integral = //pow(ppt->k[index_md][index_k],9./2.)*Y_000;
      2.*(Y_000-3.*Y_110)*(Z_111-3.*Z_100)
         -2.*(Y_111-3.*Y_100)*(Z_000-3.*Z_110);
          //2.*(Y_000-3.*Y_110)*exp(-3.*ppt->k[index_md][index_k]/sqrt(2.39))/2.
          //-2.*(Y_111-3.*Y_100)
          // *exp(-3.*ppt->k[index_md][index_k]/sqrt(2.39))
          // /pow(ppt->k[index_md][index_k],3.);
      
      else{

        //1
        double integral_123 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(
            Z_123_111
            -Z_123_001
            -Z_123_010
            -Z_123_100
          )
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_123_000
            -Z_123_110
            -Z_123_101
           -Z_123_011);


        //printf("\nintegral123=%e, efolds.\n",Z_123_001);
        //printf("integral123=%e, efolds.\n",Z_123_010);
        //printf("\nintegral123=%e, efolds.\n",Z_123_100);
        //printf("integral123=%e, efolds.\n",Z_123_000);
        // printf("integral123=%e, efolds.\n",Z_123_110);
        //printf("integral123=%e, efolds.\n",Z_123_101);
        //printf("integral123=%e, efolds.\n",Z_123_011);
        //printf("integral123=%e, efolds.\n",Z_123_111);
        
        // printf("\nintegral123=%e, efolds.\n",1.
        //       -(Y_000-3.*Y_011/*-Y_101-Y_110*/)*(
        //     Z_123_111
             //-Z_123_001
             //-Z_123_010
        //      -3.*Z_123_100)
        //       /((Y_111-3.*Y_100/*-Y_010-Y_001*/)*(
        //     Z_123_000
             //-Z_123_110
             //-Z_123_101
        //    -3.*Z_123_011)));
        /*printf("integral123=%e, efolds.\n",((Y_111-Y_100-Y_010-Y_001)*(Z_123_000
            -Z_123_110
            -Z_123_101
            -Z_123_011)));*/

        //2
        double integral_1d23 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d23_111
            -Z_1d23_001
            -Z_1d23_010
            -Z_1d23_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d23_000
            -Z_1d23_110
            -Z_1d23_101
            -Z_1d23_011);
        /*printf("integral1d23=%e, efolds.\n\n\n",1.-(Y_000-Y_011-Y_101-Y_110)*(Z_1d23_111
            -Z_1d23_001
            -Z_1d23_010
                                                                               -Z_1d23_100)/((Y_111-Y_100-Y_010-Y_001)*(Z_1d23_000
            -Z_1d23_110
            -Z_1d23_101
            -Z_1d23_011)));*/
                /*     printf("integral1d23=%e, efolds.\n\n\n\n",(Y_111-Y_100-Y_010-Y_001)*(Z_1d23_000
            -Z_1d23_110
            -Z_1d23_101
            -Z_1d23_011));*/

        //3
        double integral_12d3 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_12d3_111
            -Z_12d3_001
            -Z_12d3_010
            -Z_12d3_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_12d3_000
            -Z_12d3_110
            -Z_12d3_101
            -Z_12d3_011);
        //4
        double integral_123d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_123d_111
            -Z_123d_001
            -Z_123d_010
            -Z_123d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_123d_000
            -Z_123d_110
            -Z_123d_101
            -Z_123d_011);
        //5
        double integral_1d2d3 =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d2d3_111
            -Z_1d2d3_001
            -Z_1d2d3_010
            -Z_1d2d3_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d2d3_000
            -Z_1d2d3_110
            -Z_1d2d3_101
            -Z_1d2d3_011);
        
        //6
        double integral_1d23d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(Z_1d23d_111
            -Z_1d23d_001
            -Z_1d23d_010
            -Z_1d23d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
          *(Z_1d23d_000
            -Z_1d23d_110
            -Z_1d23d_101
            -Z_1d23d_011);
        
        //7
        double integral_12d3d =
          2.*(Y_000-Y_011-Y_101-Y_110)
          *(
            Z_12d3d_111
          -Z_12d3d_001
          -Z_12d3d_010
          -Z_12d3d_100)
          -2.*(Y_111-Y_100-Y_010-Y_001)
           *(Z_12d3d_000
           -Z_12d3d_110
           -Z_12d3d_101
          -Z_12d3d_011);


        //printf("IN TABLE: Y_110=%e.\n",Z_12d3d_110);

        
        //field_redef
       double field_redef =
         ppt->fields_at_tau_end_lqc
               [index_md]
               [ppt->index_tp_bispectrum_s_tau_end_field_redef]
               [index_k];


       
        
        integral =
          integral_123 //1
          +integral_1d23 //2
          +integral_12d3 //3
          +integral_123d //4
          +integral_1d2d3 //5
          +integral_1d23d //6
          +integral_12d3d //7
          +field_redef; //redef

        /*
        printf("integral_123=%e, efolds.\n",integral_123);
        printf("integral_1d23=%e, efolds.\n",integral_1d23);
        printf("integral_12d3=%e, efolds.\n",integral_12d3);
        printf("Z_123d_100=%e, efolds.\n",Z_123d_001);
        printf("integral_1d2d3=%e, efolds.\n",integral_1d2d3);
        printf("integral_12d3d=%e, efolds.\n",integral_12d3d);
        printf("integral_1d23d=%e, efolds.\n",integral_1d23d);
        printf("field_redef=%e, efolds.\n",field_redef);

        printf("integral=%e, efolds.\n",integral);
        */

      
 }
                           

double pk_s_1 =ppt->fields_at_tau_end_lqc[index_md]
                                         [ppt->index_tp_pk_s_1]
                                         [index_k];
double pk_s_2 = ppt->fields_at_tau_end_lqc[index_md]
                                          [ppt->index_tp_pk_s_2]
                                          [index_k];
double pk_s_3 = ppt->fields_at_tau_end_lqc[index_md]
                                          [ppt->index_tp_pk_s_3]
                                          [index_k];



double pk_t_1 =ppt->fields_at_tau_end_lqc[index_md]
                                         [ppt->index_tp_gw]
                                         [index_k];
double pk_t_2 = ppt->fields_at_tau_end_lqc[index_md]
                                          [ppt->index_tp_gw]
                                          [index_k];
double pk_t_3 = ppt->fields_at_tau_end_lqc[index_md]
                                          [ppt->index_tp_gw]
                                          [index_k];



 
 double k_1 = ppt->k[index_md][index_k];
 double k_2 = ppt->lambda_k_2*k_1;
 double k_3 = ppt->lambda_k_3*k_1;
//double k_3 = 1.e-4*pba->k_star;//V. Sreenath

 
//k_2 = 3.e0*pba->k_star;
//k_3 = 3.e0*pba->k_star;
  // double theta = 5.*_PI_/5.;
  //k_3 = sqrt(k_1*k_1+k_2*k_2+2.*k_1*k_2*cos(theta));
  //k_3 = ppt->lambda_k_3*ppt->k[index_md][index_k];

  
 
 
 double integral_s= -10./3.
  *pow(2.*_PI_,-4.)
  *pow(k_1*k_2*k_3,3.)
  *integral
  *pow(pow(k_1,3.)*pk_s_2*pk_s_3
       +pow(k_2,3.)*pk_s_1*pk_s_3
       +pow(k_3,3.)*pk_s_2*pk_s_1,-1.);//V. Sreenath's comment fNL

//double integral_s= pow(k_1*k_3,3.)*integral; //V. Sreenath Bispectrum

double integral_t =
  8.
  *4.*4.
  *pow(2.*_PI_*_PI_,-2.)
  *pow(k_1*k_2*k_3,3.)
  *(pow(k_1,2.)+pow(k_2,2.)+pow(k_3,2.))
  *integral
  *pow(pow(k_1,3.)*pk_t_2*pk_t_3
       +pow(k_2,3.)*pk_t_1*pk_t_3
       +pow(k_3,3.)*pk_t_2*pk_t_1,-1.)
  /2.;

printf(" k = %e, integral = %e, fnl = %e \n", k_1, integral, integral_s);
 
 //integral = (Z_123_000-3.*Z_123_011);
 double epsilon_2 = pba->epsilon_2;

 /*integral =  epsilon_2/(4.*pk_s_2)
   *pow(2.*_PI_,-3.)
   *pow(k_1,3.)
   *integral;
   */
      
     if (index_md==ppt->index_md_scalars) 
       ppt->bispectrum_s_equi[index_k] = integral_s;
     if (index_md==ppt->index_md_tensors)
       ppt->bispectrum_t_equi[index_k] = integral_t;
         integral_t;
         (Z_111-3*Z_100);
         integral_t/(
  8.
  *4.*4.
  *pow(2.*_PI_*_PI_,-2.)
  *pow(k_1*k_2*k_3,3.)
  *(pow(k_1,2.)+pow(k_2,2.)+pow(k_3,2.))
   *(Y_000-3.*Y_110)
  *exp(-3.*k_1/sqrt(2.39))/2.
  *pow(pow(k_1,3.)*pk_t_2*pk_t_3
       +pow(k_2,3.)*pk_t_1*pk_t_3
       +pow(k_3,3.)*pk_t_2*pk_t_1,-1.)
  /2.);

        

    }//end loop over k 
    }//end has bispeectrum

    //Now collect the primordial power spectra in ppt->ln_pk
  double source_ic1;
    for (index_k=0; index_k<ppt->ln_k_size; index_k++) {
      
    ppt->ln_k[index_k]=log(ppt->k[index_md][index_k]);

    if (index_md==ppt->index_md_scalars) source_ic1 =
                                           ppt->fields_at_tau_end_lqc
                                           [index_md]
                                           [ppt->index_tp_phi_scf]
                                           [index_k];

        if (index_md==ppt->index_md_tensors) source_ic1 =
                                               ppt->fields_at_tau_end_lqc
                                               [index_md]
                                               [ppt->index_tp_gw]
                                               [index_k];
       ppt->ln_pk[index_k] = log(source_ic1);

  }



  return _SUCCESS_;
}



int output_pk_lqc(
              struct background_lqc * pba,
              struct perturbs_lqc * ppt,
              int index_md
              ) {

  /** Summary: */


  FILE * out;     
  double * pk_tot;
  int index_k;

  FileName file_name;
  char first_line[_LINE_LENGTH_MAX_];


    spectra_pk_lqc(pba,ppt,index_md);

    if (index_md==ppt->index_md_scalars){
    sprintf(file_name,"%s_%s",ppt->root,"pk_lqc_s.dat");
    class_call(output_open_pk_file_lqc(pba,
                                   ppt,
                                   &out,
                                   file_name,
                                   "scalars"),
               ppt->error_message,
               ppt->error_message);   }

    
    else if (index_md==ppt->index_md_tensors){
    sprintf(file_name,"%s_%s",ppt->root,"pk_lqc_t.dat");

    class_call(output_open_pk_file_lqc(pba,
                                   ppt,
                                   &out,
                                   file_name,
                                   "tensors"),
               ppt->error_message,
               ppt->error_message);  }

    
    class_alloc(pk_tot,
                ppt->ln_k_size*sizeof(double),
                ppt->error_message);



    //POWER SPECTRUM
    
      for (index_k=0; index_k<ppt->ln_k_size; index_k++) {

          pk_tot[index_k] = exp(ppt->ln_pk[index_k]);

      class_call(output_one_line_of_pk_lqc(out,
                                       exp(ppt->ln_k[index_k]),
                                       pk_tot[index_k]),
                 ppt->error_message,
                 ppt->error_message);
    }


    //SCALAR BISPECTRUM


    if (index_md==ppt->index_md_scalars){
    sprintf(file_name,"%s_%s",ppt->root,"bispectrum_equi_lqc_s.dat");
    class_call(output_open_pk_file_lqc(pba,
                                   ppt,
                                   &out,
                                   file_name,
                                   "scalars"),
               ppt->error_message,
               ppt->error_message);   }


      //TENSOR  BISPECTRUM

    
    else if (index_md==ppt->index_md_tensors){
    sprintf(file_name,"%s_%s",ppt->root,"bispectrum_equi_lqc_t.dat");

    class_call(output_open_pk_file_lqc(pba,
                                   ppt,
                                   &out,
                                   file_name,
                                   "tensors"),
               ppt->error_message,
               ppt->error_message);  }

    
    class_alloc(pk_tot,
                ppt->ln_k_size*sizeof(double),
                ppt->error_message);



    
      for (index_k=0; index_k<ppt->ln_k_size; index_k++) {




        
if (index_md==ppt->index_md_scalars)
          pk_tot[index_k] = ppt->bispectrum_s_equi[index_k];
        
else if (index_md==ppt->index_md_tensors)
          pk_tot[index_k] = ppt->bispectrum_t_equi[index_k];

//printf("Mode b=%e, efolds.\n",pk_tot[index_k]);
      class_call(output_one_line_of_pk_lqc(out,
                                       exp(ppt->ln_k[index_k]),
                                       pk_tot[index_k]),
                 ppt->error_message,
                 ppt->error_message);
    }


      
    free(pk_tot);
    fclose(out);

    


  return _SUCCESS_;

}

int output_one_line_of_pk_lqc(
                          FILE * pkfile,
                          double one_k,
                          double one_pk
                          ) {

  fprintf(pkfile," ");
  class_fprintf_double(pkfile,one_k,_TRUE_);
  class_fprintf_double(pkfile,one_pk,_TRUE_);
  fprintf(pkfile,"\n");

  return _SUCCESS_;

}


int output_open_pk_file_lqc(
                        struct background_lqc * pba,
                        struct perturbs_lqc * ppt,
                        FILE * * pkfile,
                        FileName filename,
                        char * first_line
                            ) {

  int colnum = 1;

  class_open(*pkfile,filename,"w",ppt->error_message);

    fprintf(*pkfile,"# LQC primordial power spectrum P(k) %s\n",first_line);
    fprintf(*pkfile,"# for k=%g to %g,\n",
            exp(ppt->ln_k[0]),
            exp(ppt->ln_k[ppt->ln_k_size-1]));
    fprintf(*pkfile,"# number of wavenumbers equal to %d\n",ppt->ln_k_size);

    fprintf(*pkfile,"#");
    class_fprintf_columntitle(*pkfile,"k",_TRUE_,colnum);
    class_fprintf_columntitle(*pkfile,"P",_TRUE_,colnum);

    fprintf(*pkfile,"\n");
  return _SUCCESS_;
}

