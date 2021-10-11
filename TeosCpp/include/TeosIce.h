#ifndef TEOSICE_H
#define TEOSICE_H


/******************************************
  "TeosIce.h" Version 1.06
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

/****************************
  "TeosIce" C++ Class
  -------------------
  The methods of this class
  focus more on ice related
  issues (the "TeosSea"
  class focuses more on
  non-frozen seawater).
  ---------------------
  TEOS descriptions of the
  various methods are in
  the relevant ".cpp" file.
****************************/
class TeosIce:  public TeosBase
{
public:
    TeosIce();
    virtual ~TeosIce();

   /************ 'A' ***********/
   double gsw_adiabatic_lapse_rate_ice(double t, double p);
   double gsw_alpha_wrt_t_ice(double t, double p);

   /************ 'C' ***********/
   double gsw_chem_potential_water_ice(double t, double p);
   double gsw_cp_ice(double t, double p);
   double gsw_ct_freezing_exact(double sa, double p, double saturation_fraction);
   void gsw_ct_freezing_first_derivatives(double sa, double p, double saturation_fraction,
															double *ctfreezing_sa, double *ctfreezing_p);

   void gsw_ct_freezing_first_derivatives_poly(double sa, double p,
				double saturation_fraction, double *ctfreezing_sa, double *ctfreezing_p);

   double gsw_ct_freezing_poly(double sa, double p, double saturation_fraction=0.0);
   void gsw_ct_from_rho(double rho, double sa, double p, double *ct, double *ct_multiple);

   /************ 'D' ***********/
   double gsw_deltasa_from_sp(double sp, double p, double lon, double lat);

   /************ 'E' ***********/
   double gsw_enthalpy_diff(double sa, double ct, double p_shallow, double p_deep);
   double gsw_enthalpy_ice(double t, double p);
   double gsw_entropy_ice(double t, double p);

   /************ 'F' **********/
   void gsw_frazil_properties(double sa_bulk, double h_bulk, double p,
											double *sa_final, double *ct_final, double *w_ih_final);

   void gsw_frazil_properties_potential(double sa_bulk, double h_pot_bulk,
														 double p, double *sa_final,
														 double *ct_final, double *w_ih_final);

   void gsw_frazil_properties_potential_poly(double sa_bulk, double h_pot_bulk,
				double p, double *sa_final, double *ct_final, double *w_ih_final);

   void gsw_frazil_ratios_adiabatic(double sa, double p, double w_ih,
													double *dsa_dct_frazil,
													double *dsa_dp_frazil, double *dct_dp_frazil);

   void gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih,
														  double *dsa_dct_frazil,
														  double *dsa_dp_frazil,
														  double *dct_dp_frazil);

    /************ 'G' ***********/
    double *gsw_geo_strf_dyn_height_pc(double *sa, double *ct, double *delta_p,
            int n_levels, double *geo_strf_dyn_height_pc, double *p_mid);

   double gsw_gibbs_ice_part_t(double t, double p);
   double gsw_gibbs_ice_pt0(double pt0);
   double gsw_gibbs_ice_pt0_pt0(double pt0);

    /************ 'H' ***********/
    double gsw_helmholtz_energy_ice(double t, double p);

    /************ 'I' ***********/
    void   gsw_ice_fraction_to_freeze_seawater(double sa, double ct, double p,
            double t_ih, double *sa_freeze, double *ct_freeze, double *w_ih);

    double gsw_internal_energy_ice(double t, double p);
    void   gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p,
            double p_ref, int nz, double *ipv_vs_fnsquared_ratio, double *p_mid);

    /************ 'K' ***********/
    double gsw_kappa_const_t_ice(double t, double p);
    double gsw_kappa_ice(double t, double p);

    /************ 'L' ***********/
    double gsw_latentheat_melting(double sa, double p);

    /************ 'M' ***********/
    double gsw_melting_ice_equilibrium_sa_ct_ratio(double sa, double p);
    double gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa, double p);
    void   gsw_melting_ice_into_seawater(double sa, double ct, double p,
            double w_ih, double t_ih, double *sa_final, double *ct_final,
            double *w_ih_final);

    double gsw_melting_ice_sa_ct_ratio(double sa, double ct, double p,
            double t_ih);

    double gsw_melting_ice_sa_ct_ratio_poly(double sa, double ct, double p,
            double t_ih);

    double gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa, double p);
    double gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa, double p);
    void   gsw_melting_seaice_into_seawater(double sa, double ct, double p,
            double w_seaice, double sa_seaice, double t_seaice, double *sa_final,
            double *ct_final);

    double gsw_melting_seaice_sa_ct_ratio(double sa, double ct, double p,
            double sa_seaice, double t_seaice);

    double gsw_melting_seaice_sa_ct_ratio_poly(double sa, double ct,
            double p, double sa_seaice, double t_seaice);

    /************ 'N' ***********/
    void   gsw_nsquared(double *sa, double *ct, double *p, double *lat,
            int nz, double *n2, double *p_mid);

    /************ 'P' ***********/
    double gsw_pot_enthalpy_ice_freezing(double sa, double p);
    void gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa, double p,
				double *pot_enthalpy_ice_freezing_sa, double *pot_enthalpy_ice_freezing_p);

    double gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice);
    double gsw_pot_enthalpy_from_specvol_ice(double specvol_ice, double p);
    double gsw_pressure_coefficient_ice(double t, double p);
    double gsw_pressure_freezing_ct(double sa, double ct,
            double saturation_fraction);

    double gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice);
    double gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice);
    double gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice);
    double gsw_pt_from_t_ice(double t, double p, double p_ref);
    double gsw_pt0_cold_ice_poly(double pot_enthalpy_ice);


    /************ 'R' ***********/
    double gsw_rho_ice(double t, double p);

    /************ 'S' ***********/
    double gsw_sa_freezing_estimate(double p, double saturation_fraction,
            double *ct, double *t);

    double gsw_sa_freezing_from_ct(double ct, double p,
            double saturation_fraction);

    double gsw_sa_freezing_from_ct_poly(double ct, double p,
            double saturation_fraction);

    double gsw_sa_freezing_from_t(double t, double p,
            double saturation_fraction);

    double gsw_sa_freezing_from_t_poly(double t, double p,
            double saturation_fraction);

    int    gsw_sa_p_inrange(double sa, double p);
    void   gsw_seaice_fraction_to_freeze_seawater(double sa, double ct,
            double p, double sa_seaice, double t_seaice, double *sa_freeze,
            double *ct_freeze, double *w_seaice);

   double gsw_specvol_from_pot_enthalpy_ice(double pot_enthalpy_ice, double p);
   double gsw_specvol_ice(double t, double p);

   double gsw_sound_speed_ice(double t, double p);

   /************ 'T' ***********/
   double gsw_pt0_from_t_ice(double t, double p);
   double gsw_t_freezing_exact(double sa, double p, double saturation_fraction);
   void gsw_t_freezing_first_derivatives(double sa, double p,
														  double saturation_fraction,
														  double *tfreezing_sa, double *tfreezing_p);

   void gsw_t_freezing_first_derivatives_poly(double sa, double p,
				double saturation_fraction, double *tfreezing_sa, double *tfreezing_p);

   double gsw_t_freezing_poly(double sa, double p, double saturation_fraction=0.0);
   double gsw_t_freezing_poly(double sa, double p, double saturation_fraction,int polynomial);
   double gsw_t_from_pt0_ice(double pt0_ice, double p);
   double gsw_t_from_rho_ice(double rho_ice, double p);
};
#endif
