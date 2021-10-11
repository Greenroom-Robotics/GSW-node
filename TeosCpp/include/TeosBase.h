#ifndef TEOSBASE_H
#define TEOSBASE_H

/*****************************************
  "TeosBase.h" Version 1.06
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is a modification of the
  TEOS-10 file "gswteos-10.h"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  This version also contains conversion
  of .m Matlab version files into C++.
  ---------------------------------------
  All copyrights and all license issues
  are the same as for the C version and
  the Matlab version as expressed on the
  TEOS-10 website http://www.teos-10.org
-----------------------------------------
  "TeosBase" is the base or ancestor
  class for the "TeosSea", "TeosIce".
*****************************************/

/** std headers */
#include <string.h>
#include <math.h>
#include <complex>
#include <list>
#include <errno.h>

using namespace std;

/****************************
  The "gsw_structs.h" file
  has been converted to C++
  encapsulated format by
  placing data within
  structs.
*****************************/
#include "gsw_structs.h"

/************************
  Header file to the
  SAAR object for access
  to SAAR data arrays.
************************/
#include "SaarDataHandler.h"

/** error return values from various gsw methods */
const double cppGSW_INVALID_VALUE = 9e15;
const double cppGSW_ERROR_LIMIT = 1e10;
const double rtNaN = NAN;
const double rtNanF = NAN;

/** used to list/sort arrays */
struct DI
{
	double d;
	int i;
};

/********************************
  TEOS-10 data in the "TeosBase"
  class is handled by an object
  of the "SaarDataHandler" class
  and various C++ structs.
  ------------------------------
  Descendant classes have full
  access to the "gsw_structs"
  and "gsw_saar_data"
*********************************/

class TeosBase
{
	public:
		TeosBase();
		virtual ~TeosBase();

		/** these complex numbers are initialized in the constructor */
		complex<double> t1;
		complex<double> t2;
		complex<double> r1;
		complex<double> r20;
		complex<double> r21;
		complex<double> r22;

		/** Specialized data handler */
		SaarDataHandler *pSaarData;

		/********************************
		  struct data used by "TeosSea"
		  & "TeosIce" objects are
		  initialized in the
		  "TeosBase" constructor.
		********************************/
		GSW_TEOS10_CONSTANTS gtc;
		GSW_SPECVOL_COEFFICIENTS gsvco;
		GSW_SP_COEFFICIENTS gspc;
		GSW_BALTIC_DATA gbalt;
		GSW_SAAR_DATA gsaar;
		GSW_FREEZING_POLY_COEFFICIENTS gfpc;
		GSW_GIBBS_ICE_COEFFICIENTS gic;

		bool statusOk;

		/*************************************
		  Previous "C functions" have become
		  "C++ methods" that are used by both
		  the "TeosSea" & "TeosIce" objects
		  via the "TeosBase" base class or
		  their own methods.
		  ----------------------------------
		  Their C library names as well as
		  functionality, as described in
		  the "TeosSea.cpp" and "TeosIce.cpp"
		  files (descriptions by TEOS org),
		  have been preserved.
		**************************************/

		/**************************************
		   The following function/method names
		   are alphabetized beginning with the
		   first letter following 'gsw_'
		***************************************/

		/** A */
		double gsw_Abs_Pressure_from_p(double p);
		void gsw_add_barrier(double *input_data, double lon, double lat, double long_grid,
									double lat_grid, double dlong_grid,
									double dlat_grid, double *output_data);

		void gsw_add_mean(double *data_in, double *data_out);
		double gsw_adiabatic_lapse_rate_from_t(double SA, double t, double p);
		double gsw_adiabatic_lapse_rate_t_exact(double SA, double t, double p);
		double gsw_alpha(double sa, double ct, double p);
		double gsw_alpha_CT_exact(double SA, double CT, double p);
		double gsw_alpha_on_beta_CT_exact(double SA, double CT, double p);
		double gsw_alpha_wrt_CT_t_exact(double SA, double t, double p);
		double gsw_alpha_wrt_pt_t_exact(double SA, double t, double p);
		double gsw_Arsol(double SA, double CT, double p, double lon, double lat);
		double gsw_Arsol_SP_pt(double SP, double pt);
		double gsw_atomic_weight(void);


		/** B */
		double gsw_beta(double sa, double ct, double p);
		double gsw_beta_const_ct_t_exact(double SA, double t, double p);
		double gsw_beta_const_pt_t_exact(double SA, double t, double p);


		/** C */
		double gsw_chem_potential_relative_t_exact(double SA, double t, double p);
		double gsw_chem_potential_salt_t_exact(double SA, double t, double p);
		double gsw_chem_potential_water_t_exact(double sa, double t, double p);
		void gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p,
				                                    double *ct_sa_wrt_t,
				                                    double *ct_t_wrt_t,
				                                    double *ct_p_wrt_t);

		double gsw_ct_freezing(double sa, double p, double saturation_fraction=0.0);
		double gsw_CT_freezing(double SA, double saturation_fraction=0.0);
		double gsw_ct_from_enthalpy(double sa, double h, double p);
		double gsw_ct_from_enthalpy_exact(double sa, double h, double p);
		double gsw_ct_from_pt(double sa, double pt);
		double gsw_ct_from_t(double sa, double t, double p);
		double gsw_ct_maxdensity(double sa, double p);
		double gsw_CT_maxdensity_exact(double SA, double p);

		/** D */
		double gsw_deltasa_atlas(double p, double lon, double lat);
		double gsw_depth_from_z(double z);
		double gsw_dilution_coefficient_t_exact(double sa, double t, double p);
		double gsw_dynamic_enthalpy(double sa, double ct, double p);
		double gsw_dynamic_enthalpy_CT_exact(double SA, double CT, double p);
		double gsw_dynamic_enthalpy_t_exact(double SA, double t, double p);


		/** E */
		double gsw_enthalpy(double sa, double ct, double p);
		double gsw_enthalpy_ct_exact(double sa, double ct, double p);
		double gsw_enthalpy_diff_CT_exact(double SA, double CT, double p_shallow, double p_deep);
		void gsw_enthalpy_first_derivatives(double sa, double ct, double p,
														double *h_sa, double *h_ct);

		void gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct,
				double p, double *h_sa, double *h_ct);

		void gsw_enthalpy_second_derivatives_ct_exact(double sa, double ct, double p,
				double &h_sa_sa, double &h_sa_ct, double &h_ct_ct);

		double gsw_enthalpy_sso_0(double p);
		double gsw_enthalpy_t_exact(double sa, double t, double p);

		void gsw_entropy_first_derivatives(double sa, double ct, double *eta_sa,
													  double *eta_ct);

		double gsw_entropy_from_ct(double sa, double ct);
		double gsw_entropy_part(double sa, double t, double p);
		double gsw_entropy_part_zerop(double sa, double pt0);
		void gsw_entropy_second_derivatives(double sa, double ct, double *eta_sa_sa,
														double *eta_sa_ct, double *eta_ct_ct);

      /** F */
		double gsw_f(double lat);
		double gsw_fdelta(double p, double lon, double lat);

      /** G */
		double gsw_gibbs(int ns,int nt,int np, double sa, double t, double p);
		double gsw_Gibbs_energy_t_exact(double SA, double t, double p);
		double gsw_gibbs_ice(int nt,int np, double t, double p);
		double gsw_gibbs_pt0_pt0(double sa, double pt0);
		double gsw_grav(double lat, double p);

		/** H */
		double gsw_Helmholtz_energy_t_exact(double SA, double t, double p);
		double gsw_Hesol(double SA, double CT, double p, double lon, double lat);
		double gsw_Hesol_SP_pt(double SP, double pt);
		double gsw_hill_ratio_at_sp2(double t);
		void gsw_Hill_ratio_at_SP2(const double t_data[],const int t_size[2],
											double Hill_ratio_data[],int Hill_ratio_size[2]);


		/** I */
		double gsw_internal_energy(double sa, double ct, double p);
		double gsw_internal_energy_CT_exact(double SA, double CT, double p);

		void gsw_internal_energy_second_derivatives(double SA, double CT, double p,
				double &u_SA_SA, double &u_SA_CT,
				double &u_CT_CT, double &u_SA_P,
				double &u_CT_P);

		double gsw_internal_energy_t_exact(double SA, double t, double p);
		double gsw_isochoric_heat_cap_t_exact(double SA, double t, double p);

		/** K */
		double gsw_kappa(double sa, double ct, double p);
		double gsw_kappa_const_t_exact(double SA, double t, double p);
		double gsw_kappa_CT_exact(double SA, double CT, double p);
		double gsw_kappa_t_exact(double sa, double t, double p);
		double gsw_Krsol(double SA, double CT, double p, double lon, double lat);
		double gsw_Krsol_SP_pt(double SP, double pt);

		/** L */
		void gsw_linear_interp_sa_ct(double *sa, double *ct, double *p,int np,
											  double *p_i,int npi, double *sa_i, double *ct_i);

		int gsw_linear_interp_sa_ct_for_dh(double *sa, double *ct, double *p,
													  int nz, double *p_i,int n_i, double *sa_i, double *ct_i);

      /** M */
		double gsw_molality_from_SA(double SA);


      /** N */
		double gsw_N2sol(double SA, double CT, double p, double lon, double lat);
		double gsw_N2sol_SP_pt(double SP, double pt);
		double gsw_N2Osol(double SA, double CT, double p, double lon, double lat);
		double gsw_N2Osol_SP_pt(double SP, double pt);
		double gsw_Nsquared_lowerlimit(double p, double b_long, double lat);

      /** O */
		double gsw_osmotic_coefficient_t_exact(double SA, double t, double p);
		double gsw_osmotic_pressure_t_exact(double SA, double t, double pw);
		double gsw_o2sol(double sa, double ct, double p, double lon, double lat);
		double gsw_o2sol_sp_pt(double sp, double pt);

		/** P */
		double gsw_p_from_Abs_Pressure(double Absolute_Pressure);
		double gsw_p_from_z(double z, double lat, double geo_strf_dyn_height=0.0,
								  double sea_surface_geopotential=0.0);
		int gsw_pchip_derivs(double *x, double *y,int n, double *d);
		double gsw_pchip_edge_case(double h0, double h1, double m0, double m1);
		double gsw_pot_enthalpy_from_pt(double SA, double pt);
		double gsw_pot_enthalpy_from_pt_ice(double pt0_ice);
		double gsw_pot_enthalpy_from_specvol_ice_poly(double specvol_ice, double p);

		void gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa, double p,
				double *pot_enthalpy_ice_freezing_sa, double *pot_enthalpy_ice_freezing_p);

		double gsw_pot_enthalpy_ice_freezing_poly(double sa, double p);
		double gsw_pt_from_ct(double sa, double ct);
		double gsw_pt_from_t(double sa, double t, double p, double p_ref=0.0);
		double gsw_pt0_from_t(double sa, double t, double p);

		/** R */
		double gsw_R_from_SP(double SP, double t, double p);
		double gsw_rho(double sa, double ct, double p);
		void gsw_rho_alpha_beta(double sa, double ct, double p, double *rho,
										double *alpha, double *beta);

		void gsw_rho_alpha_beta_indexed(double SA, double CT, double p, double &rho,
												  double &alpha, double &beta);

		void gsw_rho_first_derivatives_CT_exact(double SA, double CT, double p,
															 double &rho_SA, double &rho_CT,
															 double &rho_P);

		void gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(double SA, double CT,
				double p, double &rho_SA,
				double &rho_h);

		void gsw_rho_second_derivatives_CT_exact(double SA, double CT, double p,
				double *rho_SA_SA, double *rho_SA_CT, double *rho_CT_CT,
				double *rho_SA_P, double *rho_CT_P);

		void gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(double SA, double CT,
				double p, double &rho_SA_SA, double &rho_SA_h, double &rho_h_h);

		void gsw_rr68_interp_section(int sectnum, double *sa, double *ct, double *p,
											  int mp,int nsect, double *ip_sect,int *ip_isect, double *p_i,
											  double *sa_i, double *ct_i);

		bool gsw_rtIsNaN(double val);
		bool gsw_rtIsNaNF(double val);

      /** S */
		double gsw_saar(double p, double lon, double lat);
		double gsw_SA_freezing_from_CT_poly(double CT, double p,
														double saturation_fraction=0.0);

		double gsw_sa_from_sp(double sp, double p, double lon, double lat);
		double gsw_sa_from_sp_baltic(double sp, double lon, double lat);
		double gsw_sa_from_sstar(double sstar, double p, double lon, double lat);

		int gsw_sgn(double x);
		double gsw_sp_from_c(double c, double t, double p);
		double gsw_SP_from_R(double R, double t, double p);
		double gsw_sp_from_sa(double sa, double p, double lon, double lat);
		double gsw_sp_from_sa_baltic(double sa, double lon, double lat);
		double gsw_specvol(double sa, double ct, double p);

		void gsw_specvol_alpha_beta_CT_exact(double SA, double CT, double p,
														 double &specvol_CT_exact,
														 double &alpha_CT_exact,
														 double &beta_CT_exact);

		double gsw_specvol_anom(double SA, double CT, double p, double SA_ref,
										double CT_ref);

		double gsw_specvol_anom_CT_exact(double SA, double CT, double p, double SA_ref,
													double CT_ref);

		double gsw_specvol_anom_standard_CT_exact(double SA, double CT, double p);
		double gsw_specvol_anom_standard_t_exact(double SA, double t, double p);
		double gsw_specvol_anom_t_exact(double SA, double t, double p);
		double gsw_specvol_CT_exact(double SA, double CT, double p);
		double gsw_specvol_diff(double SA, double CT, double p_shallow, double p_deep);

		void gsw_specvol_first_derivatives_CT_exact(double SA, double CT, double p,
				double &v_SA, double &v_CT, double &v_P);

		void gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(double SA, double CT,
				double p, double &v_SA,
				double &v_h);

		double gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(double pot_enthalpy_ice,
				double p);

		double gsw_specvol_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice, double p);
		void gsw_specvol_p_parts(double SA, double CT, double &vp0, double &vp1,
										 double &vp2, double &vp3, double &vp4, double &vp5,
										 double &vp6);

		void gsw_specvol_second_derivatives_CT_exact(double SA, double CT, double p,
				double &v_SA_SA, double &v_SA_CT, double &v_CT_CT, double &v_SA_P,
				double &v_CT_P);

		double gsw_specvol_t_exact(double sa, double t, double p);
		double gsw_specvol_sso_0(double p);
		double gsw_sr_from_sp(double sp);
		double gsw_SSO(void);
		double gsw_sstar_from_sa(double sa, double p, double lon, double lat);
		double gsw_sstar_from_sp(double sp, double p, double lon, double lat);
		double gsw_sum(double *x,int n);

		/** T */
		double gsw_t_deriv_chem_potential_water_t_exact(double sa, double t, double p);
		double gsw_t_freezing(double sa, double p, double saturation_fraction=0.0);
		double gsw_t_from_ct(double sa, double ct, double p);
		double gsw_t_from_pt0(double SA, double pt0, double p);
		double gsw_t_maxdensity_exact(double SA, double p);
		double gsw_t90_from_t48(double t48);
		double gsw_t90_from_t68(double t68);

		/** U */
		int gsw_util_indx(double *x,int n, double z);
		void gsw_util_sort_dbl(double *dArray,int nx,int *iArray,bool ASC=false);
		double gsw_util_xinterp1(double *x, double *y,int n, double x0);


		/** W */
		void gsw_weighted_nanmean(double data, double weights, double &weighted_mean,
										  double &weighted_weights);

		/** Z */
		double gsw_z_from_depth(double depth);
		double gsw_z_from_p(double p, double lat, double geo_strf_dyn_height=0.0,
								  double sea_surface_geopotential=0.0);

	protected:
		const unsigned badBaseAlloc = 9999;

		double gsw_gsvco_a(double xs, double ys, double z);
		double gsw_gsvco_b(double xs, double ys, double z);
		double gsw_gsvco_c(double xs, double ys, double z);
		double gsw_gsvco_v(double xs, double ys, double z);

		double gsw_max_d(double a, double b);
		double gsw_min_d(double a, double b);
		double gsw_get_pot_ent(double x2, double x, double y);

		unsigned gsw_max_u(unsigned a,unsigned b);
		unsigned gsw_min_u(unsigned a,unsigned b);

		int gsw_max_i(int a,int b);
		int gsw_min_i(int a,int b);

		int gsw_util_indxPref(int n, double z);

	private:
		double gsw_get_dyn_enth_prt(double xs, double ys, double z);

		double gsw_get_g03plusg08(double z, double y, double x2, double x);
		double gsw_get_g03plusg08_h(double z, double y, double x2, double x);

		double gsw_get_v_ct_part(double xs, double ys, double z);

		double gsw_gibbs_g03plusg08_a(double sa, double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_b(double sa, double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_c(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_d(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_e(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_f(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_g(double z, double y, double x2, double x);

		double gsw_gibbs_g08_a(double sa, double z, double y, double x);
		double gsw_gibbs_g08_b(double z, double y, double x);
		double gsw_gibbs_g08_c(double sa, double z, double y, double x);
		double gsw_g08_d(double x2, double z, double y, double x);
};
#endif

