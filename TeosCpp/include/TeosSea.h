#ifndef TEOSSEA_H
#define TEOSSEA_H


/******************************************
  "TeosSea.h" Version 1.06
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is a modification of the
  TEOS-10 file "gswteos-10.h"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  OR converted from .m Matlab files in
  the Teos-10 Matlab version at
  http://www.teos-10.org
  ---------------------------------------
  All copyrights and all license issues
  are the same as for the C version
*******************************************/

/** Base-class object header */
#include "TeosBase.h"

using namespace std;

/***************************
  TEOS C++ CLASS
  ------------------------
  TEOS descriptions of the
  methods are in the .cpp
  file.
***************************/
class TeosSea:  public TeosBase
{
	public:
		TeosSea();
		virtual ~TeosSea();

		/** 'A' */
		double gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p);
		double gsw_alpha_on_beta(double sa, double ct, double p);
		double gsw_alpha_wrt_t_exact(double sa, double t, double p);

		/** 'B' */
		double gsw_beta_const_t_exact(double sa, double t, double p);

		/** 'C' */
		double gsw_cabbeling(double sa, double ct, double p);
		double gsw_cabbeling_ct_exact(double SA, double CT, double p);
		double gsw_c_from_sp(double sp, double t, double p);
		double gsw_cp_t_exact(double sa, double t, double p);
		double gsw_ct_from_entropy(double sa, double entropy);
      void gsw_ct_from_rho_exact(double rho, double SA, double p, double &CT,
											double &CT_multiple);

		void   gsw_ct_first_derivatives(double sa, double pt,
		         double *ct_sa, double *ct_pt);

		void   gsw_ct_second_derivatives(double sa, double pt, double *ct_sa_sa,
		         double *ct_sa_pt, double *ct_pt_pt);

		/** 'E' */
		void gsw_enthalpy_first_derivatives_wrt_t_exact(double SA, double t, double p,
				double &h_SA_wrt_t,
				double &h_T_wrt_t,
				double &h_P_wrt_t);

		void   gsw_enthalpy_second_derivatives(double sa, double ct, double p,
		         double *h_sa_sa, double *h_sa_ct, double *h_ct_ct);

		double gsw_entropy_from_pt(double sa, double pt);
		double gsw_entropy_from_t(double sa, double t, double p);
		/** 'G' */
		double *gsw_geo_strf_dyn_height(double *sa, double *ct, double *p,
		         double p_ref, int n_levels, double *dyn_height);

		int    gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p,
		         double p_ref, int n_levels, double *dyn_height,
		         double max_dp_i, int interp_method);

      /** I */
		void gsw_internal_energy_first_derivatives(double SA, double CT, double p,
				double &u_SA, double &u_CT,
				double &u_P);

		void gsw_internal_energy_first_derivatives_ct_exact(double SA, double CT,
				double p, double &u_SA,
				double &u_CT, double &u_P);

      void gsw_internal_energy_second_derivatives(double SA, double CT, double p,
		                                             double &u_SA_SA, double &u_SA_CT,
		                                             double &u_CT_CT, double &u_SA_P,
                                                   double &u_CT_P);

      double gsw_isopycnal_slope_ratio(double SA, double CT, double p,
                                       double p_ref=0.0);

		/** 'L' */
		double gsw_latentheat_evap_ct(double sa, double ct);
		double gsw_latentheat_evap_t(double sa, double t);
		int    gsw_linear_interp_SA_CT_for_dh(double *sa, double *ct,
		         double *p, int nz, double *p_i, int n_i, double *sa_i,
		         double *ct_i);

      /** 'N' */
      double gsw_ntp_pt_vs_CT_ratio(double SA, double CT, double p);

		/** 'P' */
		void   gsw_p_sequence(double p1, double p2, double max_dp_i,
		         double *pseq, int *nps);

		double gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref);
		void   gsw_pt_first_derivatives(double sa, double ct, double *pt_sa, double *pt_ct);
		double gsw_pt_from_entropy(double sa, double entropy);
		void   gsw_pt_second_derivatives(double sa, double ct, double *pt_sa_sa,
		         double *pt_sa_ct, double *pt_ct_ct);

		/** 'R' */
		int    gsw_refine_grid_for_dh(double *p, double p_ref, int nz, double dp,
		         double *p_i, int ni_max, int *p_indices, int *p_ref_ind_ptr);

		void gsw_rho_alpha_beta_ct_exact(double SA, double CT, double p,
													double &rho_CT_exact, double &alpha_CT_exact,
													double &beta_CT_exact);

      double gsw_rho_ct_exact(double SA, double CT, double p);

		void   gsw_rho_second_derivatives(double sa, double ct, double p,
		         double *rho_sa_sa, double *rho_sa_ct, double *rho_ct_ct,
		         double *rho_sa_p, double *rho_ct_p);

		void   gsw_rho_first_derivatives(double sa, double ct, double p,
		         double *drho_dsa, double *drho_dct, double *drho_dp);

		void   gsw_rho_first_derivatives_wrt_enthalpy(double sa, double ct,
		         double p, double *rho_sa, double *rho_h);

		void   gsw_rho_second_derivatives_wrt_enthalpy(double sa, double ct,
		         double p, double *rho_sa_sa, double *rho_sa_h, double *rho_h_h);

		double gsw_rho_t_exact(double sa, double t, double p);
		void   gsw_rr68_interp_sa_ct(double *sa, double *ct, double *p,
		         int mp, double *p_i, int mp_i, double *sa_i, double *ct_i);

		/** 'S' */
		double gsw_sa_from_rho(double rho, double ct, double p);
		double gsw_sigma0(double sa, double ct);
		double gsw_sigma0_ct_exact(double SA, double CT);
		double gsw_sigma0_pt0_exact(double SA, double pt0);
		double gsw_sigma1(double sa, double ct);
		double gsw_sigma1_ct_exact(double SA, double CT);
		double gsw_sigma2(double sa, double ct);
		double gsw_sigma2_ct_exact(double SA, double CT);
		double gsw_sigma3(double sa, double ct);
		double gsw_sigma3_ct_exact(double SA, double CT);
		double gsw_sigma4(double sa, double ct);
		double gsw_sigma4_ct_exact(double SA, double CT);
		double gsw_sound_speed(double sa, double ct, double p);
		double gsw_sound_speed_ct_exact(double SA, double CT, double p);
		double gsw_sound_speed_t_exact(double sa, double t, double p);
		double gsw_sp_from_sk(double sk);
		double gsw_sp_from_sr(double sr);
		double gsw_sp_from_sstar(double sstar, double p,double lon,double lat);
		double gsw_sp_salinometer(double rt, double t);
		void   gsw_specvol_alpha_beta(double sa, double ct, double p,
		         double *specvol, double *alpha, double *beta);

		double gsw_specvol_anom_standard(double sa, double ct, double p);
		void   gsw_specvol_first_derivatives(double sa, double ct, double p,
		         double *v_sa, double *v_ct, double *v_p);

		void   gsw_specvol_first_derivatives_wrt_enthalpy(double sa, double ct,
		         double p, double *v_sa, double *v_h);

		void   gsw_specvol_second_derivatives(double sa, double ct, double p,
		         double *v_sa_sa, double *v_sa_ct, double *v_ct_ct,
		         double *v_sa_p, double *v_ct_p);

		void   gsw_specvol_second_derivatives_wrt_enthalpy(double sa,
		         double ct, double p, double *v_sa_sa, double *v_sa_h, double *v_h_h);

		double gsw_spiciness0(double sa, double ct);
		double gsw_spiciness1(double sa, double ct);
		double gsw_spiciness2(double sa, double ct);

		/** 'T' */
		double gsw_t_from_entropy(double SA, double entropy, double p);
      void gsw_t_from_rho_exact(double rho, double SA, double p, double &t,
										  double &t_multiple);

		double gsw_thermobaric(double sa, double ct, double p);
		double gsw_thermobaric_ct_exact(double SA, double CT, double p);
		void   gsw_turner_rsubrho(double *sa, double *ct, double *p,
		         int nz, double *tu, double *rsubrho, double *p_mid);

		/** 'U' */
		double *gsw_util_interp1q_int(int nx, double *x, int *iy,
		         int nxi, double *x_i, double *y_i);

		double *gsw_util_linear_interp(int nx, double *x, int ny,
		         double *y, int nxi, double *x_i, double *y_i);

		int    gsw_util_pchip_interp(double *x, double *y, int n,
		         double *xi, double *yi, int ni);

};
#endif
