import { Buffer } from "@mcesystems/nbind/dist/shim";

export class NBindBase { free?(): void }

export class TeosBase extends NBindBase {
	/** TeosBase(); */
	constructor();

	/** float64_t gsw_Abs_Pressure_from_p(float64_t); */
	gsw_Abs_Pressure_from_p(p0: number): number;

	/** void gsw_add_barrier(float64_t *, float64_t, float64_t, float64_t, float64_t, float64_t, float64_t, float64_t *); */
	gsw_add_barrier(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** void gsw_add_mean(float64_t *, float64_t *); */
	gsw_add_mean(p0: number, p1: number): void;

	/** float64_t gsw_adiabatic_lapse_rate_from_t(float64_t, float64_t, float64_t); */
	gsw_adiabatic_lapse_rate_from_t(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_adiabatic_lapse_rate_t_exact(float64_t, float64_t, float64_t); */
	gsw_adiabatic_lapse_rate_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha(float64_t, float64_t, float64_t); */
	gsw_alpha(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_CT_exact(float64_t, float64_t, float64_t); */
	gsw_alpha_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_on_beta_CT_exact(float64_t, float64_t, float64_t); */
	gsw_alpha_on_beta_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_wrt_CT_t_exact(float64_t, float64_t, float64_t); */
	gsw_alpha_wrt_CT_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_wrt_pt_t_exact(float64_t, float64_t, float64_t); */
	gsw_alpha_wrt_pt_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_Arsol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_Arsol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_Arsol_SP_pt(float64_t, float64_t); */
	gsw_Arsol_SP_pt(p0: number, p1: number): number;

	/** float64_t gsw_atomic_weight(); */
	gsw_atomic_weight(): number;

	/** float64_t gsw_beta(float64_t, float64_t, float64_t); */
	gsw_beta(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_beta_const_ct_t_exact(float64_t, float64_t, float64_t); */
	gsw_beta_const_ct_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_beta_const_pt_t_exact(float64_t, float64_t, float64_t); */
	gsw_beta_const_pt_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_chem_potential_relative_t_exact(float64_t, float64_t, float64_t); */
	gsw_chem_potential_relative_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_chem_potential_salt_t_exact(float64_t, float64_t, float64_t); */
	gsw_chem_potential_salt_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_chem_potential_water_t_exact(float64_t, float64_t, float64_t); */
	gsw_chem_potential_water_t_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_ct_first_derivatives_wrt_t_exact(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_ct_first_derivatives_wrt_t_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_ct_freezing(float64_t, float64_t, float64_t); */
	gsw_ct_freezing(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_ct_from_enthalpy(float64_t, float64_t, float64_t); */
	gsw_ct_from_enthalpy(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_ct_from_enthalpy_exact(float64_t, float64_t, float64_t); */
	gsw_ct_from_enthalpy_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_ct_from_pt(float64_t, float64_t); */
	gsw_ct_from_pt(p0: number, p1: number): number;

	/** float64_t gsw_ct_from_t(float64_t, float64_t, float64_t); */
	gsw_ct_from_t(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_ct_maxdensity(float64_t, float64_t); */
	gsw_ct_maxdensity(p0: number, p1: number): number;

	/** float64_t gsw_CT_maxdensity_exact(float64_t, float64_t); */
	gsw_CT_maxdensity_exact(p0: number, p1: number): number;

	/** float64_t gsw_deltasa_atlas(float64_t, float64_t, float64_t); */
	gsw_deltasa_atlas(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_depth_from_z(float64_t); */
	gsw_depth_from_z(p0: number): number;

	/** float64_t gsw_dilution_coefficient_t_exact(float64_t, float64_t, float64_t); */
	gsw_dilution_coefficient_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_dynamic_enthalpy(float64_t, float64_t, float64_t); */
	gsw_dynamic_enthalpy(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_dynamic_enthalpy_CT_exact(float64_t, float64_t, float64_t); */
	gsw_dynamic_enthalpy_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_dynamic_enthalpy_t_exact(float64_t, float64_t, float64_t); */
	gsw_dynamic_enthalpy_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_enthalpy(float64_t, float64_t, float64_t); */
	gsw_enthalpy(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_enthalpy_ct_exact(float64_t, float64_t, float64_t); */
	gsw_enthalpy_ct_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_enthalpy_diff_CT_exact(float64_t, float64_t, float64_t, float64_t); */
	gsw_enthalpy_diff_CT_exact(p0: number, p1: number, p2: number, p3: number): number;

	/** void gsw_enthalpy_first_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_enthalpy_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_enthalpy_first_derivatives_ct_exact(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_enthalpy_first_derivatives_ct_exact(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_enthalpy_second_derivatives_ct_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_enthalpy_second_derivatives_ct_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_enthalpy_sso_0(float64_t); */
	gsw_enthalpy_sso_0(p0: number): number;

	/** float64_t gsw_enthalpy_t_exact(float64_t, float64_t, float64_t); */
	gsw_enthalpy_t_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_entropy_first_derivatives(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_entropy_first_derivatives(p0: number, p1: number, p2: number, p3: number): void;

	/** float64_t gsw_entropy_from_ct(float64_t, float64_t); */
	gsw_entropy_from_ct(p0: number, p1: number): number;

	/** float64_t gsw_entropy_part(float64_t, float64_t, float64_t); */
	gsw_entropy_part(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_entropy_part_zerop(float64_t, float64_t); */
	gsw_entropy_part_zerop(p0: number, p1: number): number;

	/** void gsw_entropy_second_derivatives(float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_entropy_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_f(float64_t); */
	gsw_f(p0: number): number;

	/** float64_t gsw_fdelta(float64_t, float64_t, float64_t); */
	gsw_fdelta(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_gibbs(int32_t, int32_t, int32_t, float64_t, float64_t, float64_t); */
	gsw_gibbs(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): number;

	/** float64_t gsw_Gibbs_energy_t_exact(float64_t, float64_t, float64_t); */
	gsw_Gibbs_energy_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_gibbs_ice(int32_t, int32_t, float64_t, float64_t); */
	gsw_gibbs_ice(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_gibbs_pt0_pt0(float64_t, float64_t); */
	gsw_gibbs_pt0_pt0(p0: number, p1: number): number;

	/** float64_t gsw_grav(float64_t, float64_t); */
	gsw_grav(p0: number, p1: number): number;

	/** float64_t gsw_Helmholtz_energy_t_exact(float64_t, float64_t, float64_t); */
	gsw_Helmholtz_energy_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_Hesol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_Hesol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_Hesol_SP_pt(float64_t, float64_t); */
	gsw_Hesol_SP_pt(p0: number, p1: number): number;

	/** float64_t gsw_hill_ratio_at_sp2(float64_t); */
	gsw_hill_ratio_at_sp2(p0: number): number;

	/** float64_t gsw_internal_energy(float64_t, float64_t, float64_t); */
	gsw_internal_energy(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_internal_energy_CT_exact(float64_t, float64_t, float64_t); */
	gsw_internal_energy_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_internal_energy_t_exact(float64_t, float64_t, float64_t); */
	gsw_internal_energy_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_isochoric_heat_cap_t_exact(float64_t, float64_t, float64_t); */
	gsw_isochoric_heat_cap_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_kappa(float64_t, float64_t, float64_t); */
	gsw_kappa(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_kappa_const_t_exact(float64_t, float64_t, float64_t); */
	gsw_kappa_const_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_kappa_CT_exact(float64_t, float64_t, float64_t); */
	gsw_kappa_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_kappa_t_exact(float64_t, float64_t, float64_t); */
	gsw_kappa_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_Krsol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_Krsol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_Krsol_SP_pt(float64_t, float64_t); */
	gsw_Krsol_SP_pt(p0: number, p1: number): number;

	/** void gsw_linear_interp_sa_ct(float64_t *, float64_t *, float64_t *, int32_t, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_linear_interp_sa_ct(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** int32_t gsw_linear_interp_sa_ct_for_dh(float64_t *, float64_t *, float64_t *, int32_t, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_linear_interp_sa_ct_for_dh(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): number;

	/** float64_t gsw_molality_from_SA(float64_t); */
	gsw_molality_from_SA(p0: number): number;

	/** float64_t gsw_N2sol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_N2sol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_N2sol_SP_pt(float64_t, float64_t); */
	gsw_N2sol_SP_pt(p0: number, p1: number): number;

	/** float64_t gsw_N2Osol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_N2Osol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_N2Osol_SP_pt(float64_t, float64_t); */
	gsw_N2Osol_SP_pt(p0: number, p1: number): number;

	/** float64_t gsw_Nsquared_lowerlimit(float64_t, float64_t, float64_t); */
	gsw_Nsquared_lowerlimit(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_osmotic_coefficient_t_exact(float64_t, float64_t, float64_t); */
	gsw_osmotic_coefficient_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_osmotic_pressure_t_exact(float64_t, float64_t, float64_t); */
	gsw_osmotic_pressure_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_o2sol(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_o2sol(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_o2sol_sp_pt(float64_t, float64_t); */
	gsw_o2sol_sp_pt(p0: number, p1: number): number;

	/** float64_t gsw_p_from_Abs_Pressure(float64_t); */
	gsw_p_from_Abs_Pressure(p0: number): number;

	/** float64_t gsw_p_from_z(float64_t, float64_t, float64_t, float64_t); */
	gsw_p_from_z(p0: number, p1: number, p2: number, p3: number): number;

	/** int32_t gsw_pchip_derivs(float64_t *, float64_t *, int32_t, float64_t *); */
	gsw_pchip_derivs(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_pchip_edge_case(float64_t, float64_t, float64_t, float64_t); */
	gsw_pchip_edge_case(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_pot_enthalpy_from_pt(float64_t, float64_t); */
	gsw_pot_enthalpy_from_pt(p0: number, p1: number): number;

	/** float64_t gsw_pot_enthalpy_from_pt_ice(float64_t); */
	gsw_pot_enthalpy_from_pt_ice(p0: number): number;

	/** float64_t gsw_pot_enthalpy_from_specvol_ice_poly(float64_t, float64_t); */
	gsw_pot_enthalpy_from_specvol_ice_poly(p0: number, p1: number): number;

	/** void gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(p0: number, p1: number, p2: number, p3: number): void;

	/** float64_t gsw_pot_enthalpy_ice_freezing_poly(float64_t, float64_t); */
	gsw_pot_enthalpy_ice_freezing_poly(p0: number, p1: number): number;

	/** float64_t gsw_pt_from_ct(float64_t, float64_t); */
	gsw_pt_from_ct(p0: number, p1: number): number;

	/** float64_t gsw_pt_from_t(float64_t, float64_t, float64_t, float64_t); */
	gsw_pt_from_t(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_pt0_from_t(float64_t, float64_t, float64_t); */
	gsw_pt0_from_t(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_R_from_SP(float64_t, float64_t, float64_t); */
	gsw_R_from_SP(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_rho(float64_t, float64_t, float64_t); */
	gsw_rho(p0: number, p1: number, p2: number): number;

	/** void gsw_rho_alpha_beta(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_rho_alpha_beta(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_rho_alpha_beta_indexed(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_rho_alpha_beta_indexed(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_rho_first_derivatives_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_rho_first_derivatives_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &); */
	gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_rho_second_derivatives_CT_exact(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *, float64_t *, float64_t *); */
	gsw_rho_second_derivatives_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** void gsw_rr68_interp_section(int32_t, float64_t *, float64_t *, float64_t *, int32_t, int32_t, float64_t *, int32_t *, float64_t *, float64_t *, float64_t *); */
	gsw_rr68_interp_section(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number, p8: number, p9: number, p10: number): void;

	/** bool gsw_rtIsNaN(float64_t); */
	gsw_rtIsNaN(p0: number): boolean;

	/** bool gsw_rtIsNaNF(float64_t); */
	gsw_rtIsNaNF(p0: number): boolean;

	/** float64_t gsw_saar(float64_t, float64_t, float64_t); */
	gsw_saar(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sa_from_sp(float64_t, float64_t, float64_t, float64_t); */
	gsw_sa_from_sp(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sa_from_sp_baltic(float64_t, float64_t, float64_t); */
	gsw_sa_from_sp_baltic(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sa_from_sstar(float64_t, float64_t, float64_t, float64_t); */
	gsw_sa_from_sstar(p0: number, p1: number, p2: number, p3: number): number;

	/** int32_t gsw_sgn(float64_t); */
	gsw_sgn(p0: number): number;

	/** float64_t gsw_sp_from_c(float64_t, float64_t, float64_t); */
	gsw_sp_from_c(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_SP_from_R(float64_t, float64_t, float64_t); */
	gsw_SP_from_R(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sp_from_sa(float64_t, float64_t, float64_t, float64_t); */
	gsw_sp_from_sa(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sp_from_sa_baltic(float64_t, float64_t, float64_t); */
	gsw_sp_from_sa_baltic(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol(float64_t, float64_t, float64_t); */
	gsw_specvol(p0: number, p1: number, p2: number): number;

	/** void gsw_specvol_alpha_beta_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_specvol_alpha_beta_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_specvol_anom(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_specvol_anom(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_specvol_anom_CT_exact(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_specvol_anom_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_specvol_anom_standard_CT_exact(float64_t, float64_t, float64_t); */
	gsw_specvol_anom_standard_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol_anom_standard_t_exact(float64_t, float64_t, float64_t); */
	gsw_specvol_anom_standard_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol_anom_t_exact(float64_t, float64_t, float64_t); */
	gsw_specvol_anom_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol_CT_exact(float64_t, float64_t, float64_t); */
	gsw_specvol_CT_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol_diff(float64_t, float64_t, float64_t, float64_t); */
	gsw_specvol_diff(p0: number, p1: number, p2: number, p3: number): number;

	/** void gsw_specvol_first_derivatives_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_specvol_first_derivatives_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &); */
	gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(float64_t, float64_t); */
	gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(p0: number, p1: number): number;

	/** float64_t gsw_specvol_from_pot_enthalpy_ice_poly(float64_t, float64_t); */
	gsw_specvol_from_pot_enthalpy_ice_poly(p0: number, p1: number): number;

	/** void gsw_specvol_p_parts(float64_t, float64_t, float64_t &, float64_t &, float64_t &, float64_t &, float64_t &, float64_t &, float64_t &); */
	gsw_specvol_p_parts(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number, p8: number): void;

	/** void gsw_specvol_second_derivatives_CT_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &, float64_t &, float64_t &); */
	gsw_specvol_second_derivatives_CT_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_specvol_t_exact(float64_t, float64_t, float64_t); */
	gsw_specvol_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_specvol_sso_0(float64_t); */
	gsw_specvol_sso_0(p0: number): number;

	/** float64_t gsw_sr_from_sp(float64_t); */
	gsw_sr_from_sp(p0: number): number;

	/** float64_t gsw_sstar_from_sa(float64_t, float64_t, float64_t, float64_t); */
	gsw_sstar_from_sa(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sstar_from_sp(float64_t, float64_t, float64_t, float64_t); */
	gsw_sstar_from_sp(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sum(float64_t *, int32_t); */
	gsw_sum(p0: number, p1: number): number;

	/** float64_t gsw_t_deriv_chem_potential_water_t_exact(float64_t, float64_t, float64_t); */
	gsw_t_deriv_chem_potential_water_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_t_freezing(float64_t, float64_t, float64_t); */
	gsw_t_freezing(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_t_from_ct(float64_t, float64_t, float64_t); */
	gsw_t_from_ct(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_t_from_pt0(float64_t, float64_t, float64_t); */
	gsw_t_from_pt0(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_t_maxdensity_exact(float64_t, float64_t); */
	gsw_t_maxdensity_exact(p0: number, p1: number): number;

	/** float64_t gsw_t90_from_t48(float64_t); */
	gsw_t90_from_t48(p0: number): number;

	/** float64_t gsw_t90_from_t68(float64_t); */
	gsw_t90_from_t68(p0: number): number;

	/** int32_t gsw_util_indx(float64_t *, int32_t, float64_t); */
	gsw_util_indx(p0: number, p1: number, p2: number): number;

	/** void gsw_util_sort_dbl(float64_t *, int32_t, int32_t *, bool); */
	gsw_util_sort_dbl(p0: number, p1: number, p2: number, p3: boolean): void;

	/** float64_t gsw_util_xinterp1(float64_t *, float64_t *, int32_t, float64_t); */
	gsw_util_xinterp1(p0: number, p1: number, p2: number, p3: number): number;

	/** void gsw_weighted_nanmean(float64_t, float64_t, float64_t &, float64_t &); */
	gsw_weighted_nanmean(p0: number, p1: number, p2: number, p3: number): void;

	/** float64_t gsw_z_from_depth(float64_t); */
	gsw_z_from_depth(p0: number): number;

	/** float64_t gsw_z_from_p(float64_t, float64_t, float64_t, float64_t); */
	gsw_z_from_p(p0: number, p1: number, p2: number, p3: number): number;
}

export class TeosIce extends NBindBase {
	/** TeosIce(); */
	constructor();

	/** float64_t gsw_adiabatic_lapse_rate_ice(float64_t, float64_t); */
	gsw_adiabatic_lapse_rate_ice(p0: number, p1: number): number;

	/** float64_t gsw_alpha_wrt_t_ice(float64_t, float64_t); */
	gsw_alpha_wrt_t_ice(p0: number, p1: number): number;

	/** float64_t gsw_chem_potential_water_ice(float64_t, float64_t); */
	gsw_chem_potential_water_ice(p0: number, p1: number): number;

	/** float64_t gsw_cp_ice(float64_t, float64_t); */
	gsw_cp_ice(p0: number, p1: number): number;

	/** float64_t gsw_ct_freezing_exact(float64_t, float64_t, float64_t); */
	gsw_ct_freezing_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_ct_freezing_first_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_ct_freezing_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_ct_freezing_first_derivatives_poly(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_ct_freezing_first_derivatives_poly(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_ct_freezing_poly(float64_t, float64_t, float64_t); */
	gsw_ct_freezing_poly(p0: number, p1: number, p2: number): number;

	/** void gsw_ct_from_rho(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_ct_from_rho(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_deltasa_from_sp(float64_t, float64_t, float64_t, float64_t); */
	gsw_deltasa_from_sp(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_enthalpy_diff(float64_t, float64_t, float64_t, float64_t); */
	gsw_enthalpy_diff(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_enthalpy_ice(float64_t, float64_t); */
	gsw_enthalpy_ice(p0: number, p1: number): number;

	/** float64_t gsw_entropy_ice(float64_t, float64_t); */
	gsw_entropy_ice(p0: number, p1: number): number;

	/** void gsw_frazil_properties(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_frazil_properties(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_frazil_properties_potential(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_frazil_properties_potential(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_frazil_properties_potential_poly(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_frazil_properties_potential_poly(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_frazil_ratios_adiabatic(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_frazil_ratios_adiabatic(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_frazil_ratios_adiabatic_poly(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_frazil_ratios_adiabatic_poly(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t * gsw_geo_strf_dyn_height_pc(float64_t *, float64_t *, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_geo_strf_dyn_height_pc(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): number;

	/** float64_t gsw_gibbs_ice_part_t(float64_t, float64_t); */
	gsw_gibbs_ice_part_t(p0: number, p1: number): number;

	/** float64_t gsw_gibbs_ice_pt0(float64_t); */
	gsw_gibbs_ice_pt0(p0: number): number;

	/** float64_t gsw_gibbs_ice_pt0_pt0(float64_t); */
	gsw_gibbs_ice_pt0_pt0(p0: number): number;

	/** float64_t gsw_helmholtz_energy_ice(float64_t, float64_t); */
	gsw_helmholtz_energy_ice(p0: number, p1: number): number;

	/** void gsw_ice_fraction_to_freeze_seawater(float64_t, float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_ice_fraction_to_freeze_seawater(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number): void;

	/** float64_t gsw_internal_energy_ice(float64_t, float64_t); */
	gsw_internal_energy_ice(p0: number, p1: number): number;

	/** void gsw_ipv_vs_fnsquared_ratio(float64_t *, float64_t *, float64_t *, float64_t, int32_t, float64_t *, float64_t *); */
	gsw_ipv_vs_fnsquared_ratio(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number): void;

	/** float64_t gsw_kappa_const_t_ice(float64_t, float64_t); */
	gsw_kappa_const_t_ice(p0: number, p1: number): number;

	/** float64_t gsw_kappa_ice(float64_t, float64_t); */
	gsw_kappa_ice(p0: number, p1: number): number;

	/** float64_t gsw_latentheat_melting(float64_t, float64_t); */
	gsw_latentheat_melting(p0: number, p1: number): number;

	/** float64_t gsw_melting_ice_equilibrium_sa_ct_ratio(float64_t, float64_t); */
	gsw_melting_ice_equilibrium_sa_ct_ratio(p0: number, p1: number): number;

	/** float64_t gsw_melting_ice_equilibrium_sa_ct_ratio_poly(float64_t, float64_t); */
	gsw_melting_ice_equilibrium_sa_ct_ratio_poly(p0: number, p1: number): number;

	/** void gsw_melting_ice_into_seawater(float64_t, float64_t, float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_melting_ice_into_seawater(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_melting_ice_sa_ct_ratio(float64_t, float64_t, float64_t, float64_t); */
	gsw_melting_ice_sa_ct_ratio(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_melting_ice_sa_ct_ratio_poly(float64_t, float64_t, float64_t, float64_t); */
	gsw_melting_ice_sa_ct_ratio_poly(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_melting_seaice_equilibrium_sa_ct_ratio(float64_t, float64_t); */
	gsw_melting_seaice_equilibrium_sa_ct_ratio(p0: number, p1: number): number;

	/** float64_t gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(float64_t, float64_t); */
	gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(p0: number, p1: number): number;

	/** void gsw_melting_seaice_into_seawater(float64_t, float64_t, float64_t, float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_melting_seaice_into_seawater(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_melting_seaice_sa_ct_ratio(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_melting_seaice_sa_ct_ratio(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** float64_t gsw_melting_seaice_sa_ct_ratio_poly(float64_t, float64_t, float64_t, float64_t, float64_t); */
	gsw_melting_seaice_sa_ct_ratio_poly(p0: number, p1: number, p2: number, p3: number, p4: number): number;

	/** void gsw_nsquared(float64_t *, float64_t *, float64_t *, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_nsquared(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number): void;

	/** float64_t gsw_pot_enthalpy_ice_freezing(float64_t, float64_t); */
	gsw_pot_enthalpy_ice_freezing(p0: number, p1: number): number;

	/** void gsw_pot_enthalpy_ice_freezing_first_derivatives(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_pot_enthalpy_ice_freezing_first_derivatives(p0: number, p1: number, p2: number, p3: number): void;

	/** float64_t gsw_pot_enthalpy_from_pt_ice_poly(float64_t); */
	gsw_pot_enthalpy_from_pt_ice_poly(p0: number): number;

	/** float64_t gsw_pot_enthalpy_from_specvol_ice(float64_t, float64_t); */
	gsw_pot_enthalpy_from_specvol_ice(p0: number, p1: number): number;

	/** float64_t gsw_pressure_freezing_ct(float64_t, float64_t, float64_t); */
	gsw_pressure_freezing_ct(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_pt_from_pot_enthalpy_ice(float64_t); */
	gsw_pt_from_pot_enthalpy_ice(p0: number): number;

	/** float64_t gsw_pt_from_pot_enthalpy_ice_poly(float64_t); */
	gsw_pt_from_pot_enthalpy_ice_poly(p0: number): number;

	/** float64_t gsw_pt_from_pot_enthalpy_ice_poly_dh(float64_t); */
	gsw_pt_from_pot_enthalpy_ice_poly_dh(p0: number): number;

	/** float64_t gsw_pt_from_t_ice(float64_t, float64_t, float64_t); */
	gsw_pt_from_t_ice(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_pt0_cold_ice_poly(float64_t); */
	gsw_pt0_cold_ice_poly(p0: number): number;

	/** float64_t gsw_rho_ice(float64_t, float64_t); */
	gsw_rho_ice(p0: number, p1: number): number;

	/** float64_t gsw_sa_freezing_estimate(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_sa_freezing_estimate(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sa_freezing_from_ct(float64_t, float64_t, float64_t); */
	gsw_sa_freezing_from_ct(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sa_freezing_from_ct_poly(float64_t, float64_t, float64_t); */
	gsw_sa_freezing_from_ct_poly(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sa_freezing_from_t(float64_t, float64_t, float64_t); */
	gsw_sa_freezing_from_t(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sa_freezing_from_t_poly(float64_t, float64_t, float64_t); */
	gsw_sa_freezing_from_t_poly(p0: number, p1: number, p2: number): number;

	/** int32_t gsw_sa_p_inrange(float64_t, float64_t); */
	gsw_sa_p_inrange(p0: number, p1: number): number;

	/** void gsw_seaice_fraction_to_freeze_seawater(float64_t, float64_t, float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_seaice_fraction_to_freeze_seawater(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_specvol_from_pot_enthalpy_ice(float64_t, float64_t); */
	gsw_specvol_from_pot_enthalpy_ice(p0: number, p1: number): number;

	/** float64_t gsw_specvol_ice(float64_t, float64_t); */
	gsw_specvol_ice(p0: number, p1: number): number;

	/** float64_t gsw_sound_speed_ice(float64_t, float64_t); */
	gsw_sound_speed_ice(p0: number, p1: number): number;

	/** float64_t gsw_pt0_from_t_ice(float64_t, float64_t); */
	gsw_pt0_from_t_ice(p0: number, p1: number): number;

	/** float64_t gsw_t_freezing_exact(float64_t, float64_t, float64_t); */
	gsw_t_freezing_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_t_freezing_first_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_t_freezing_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_t_freezing_first_derivatives_poly(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_t_freezing_first_derivatives_poly(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_t_from_pt0_ice(float64_t, float64_t); */
	gsw_t_from_pt0_ice(p0: number, p1: number): number;

	/** float64_t gsw_t_from_rho_ice(float64_t, float64_t); */
	gsw_t_from_rho_ice(p0: number, p1: number): number;
}

export class TeosSea extends NBindBase {
	/** TeosSea(); */
	constructor();

	/** float64_t gsw_adiabatic_lapse_rate_from_ct(float64_t, float64_t, float64_t); */
	gsw_adiabatic_lapse_rate_from_ct(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_on_beta(float64_t, float64_t, float64_t); */
	gsw_alpha_on_beta(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_alpha_wrt_t_exact(float64_t, float64_t, float64_t); */
	gsw_alpha_wrt_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_beta_const_t_exact(float64_t, float64_t, float64_t); */
	gsw_beta_const_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_cabbeling(float64_t, float64_t, float64_t); */
	gsw_cabbeling(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_cabbeling_ct_exact(float64_t, float64_t, float64_t); */
	gsw_cabbeling_ct_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_c_from_sp(float64_t, float64_t, float64_t); */
	gsw_c_from_sp(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_cp_t_exact(float64_t, float64_t, float64_t); */
	gsw_cp_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_ct_from_entropy(float64_t, float64_t); */
	gsw_ct_from_entropy(p0: number, p1: number): number;

	/** void gsw_ct_from_rho_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &); */
	gsw_ct_from_rho_exact(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_ct_first_derivatives(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_ct_first_derivatives(p0: number, p1: number, p2: number, p3: number): void;

	/** void gsw_ct_second_derivatives(float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_ct_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_enthalpy_first_derivatives_wrt_t_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_enthalpy_first_derivatives_wrt_t_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_enthalpy_second_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_enthalpy_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_entropy_from_pt(float64_t, float64_t); */
	gsw_entropy_from_pt(p0: number, p1: number): number;

	/** float64_t gsw_entropy_from_t(float64_t, float64_t, float64_t); */
	gsw_entropy_from_t(p0: number, p1: number, p2: number): number;

	/** void gsw_internal_energy_first_derivatives(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_internal_energy_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_internal_energy_first_derivatives_ct_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_internal_energy_first_derivatives_ct_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_internal_energy_second_derivatives(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &, float64_t &, float64_t &); */
	gsw_internal_energy_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_isopycnal_slope_ratio(float64_t, float64_t, float64_t, float64_t); */
	gsw_isopycnal_slope_ratio(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_latentheat_evap_ct(float64_t, float64_t); */
	gsw_latentheat_evap_ct(p0: number, p1: number): number;

	/** float64_t gsw_latentheat_evap_t(float64_t, float64_t); */
	gsw_latentheat_evap_t(p0: number, p1: number): number;

	/** float64_t gsw_ntp_pt_vs_CT_ratio(float64_t, float64_t, float64_t); */
	gsw_ntp_pt_vs_CT_ratio(p0: number, p1: number, p2: number): number;

	/** void gsw_p_sequence(float64_t, float64_t, float64_t, float64_t *, int32_t *); */
	gsw_p_sequence(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_pot_rho_t_exact(float64_t, float64_t, float64_t, float64_t); */
	gsw_pot_rho_t_exact(p0: number, p1: number, p2: number, p3: number): number;

	/** void gsw_pt_first_derivatives(float64_t, float64_t, float64_t *, float64_t *); */
	gsw_pt_first_derivatives(p0: number, p1: number, p2: number, p3: number): void;

	/** float64_t gsw_pt_from_entropy(float64_t, float64_t); */
	gsw_pt_from_entropy(p0: number, p1: number): number;

	/** void gsw_pt_second_derivatives(float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_pt_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** int32_t gsw_refine_grid_for_dh(float64_t *, float64_t, int32_t, float64_t, float64_t *, int32_t, int32_t *, int32_t *); */
	gsw_refine_grid_for_dh(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): number;

	/** void gsw_rho_alpha_beta_ct_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &, float64_t &); */
	gsw_rho_alpha_beta_ct_exact(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_rho_ct_exact(float64_t, float64_t, float64_t); */
	gsw_rho_ct_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_rho_second_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *, float64_t *, float64_t *); */
	gsw_rho_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** void gsw_rho_first_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_rho_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_rho_first_derivatives_wrt_enthalpy(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_rho_first_derivatives_wrt_enthalpy(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_rho_second_derivatives_wrt_enthalpy(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_rho_second_derivatives_wrt_enthalpy(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_rho_t_exact(float64_t, float64_t, float64_t); */
	gsw_rho_t_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_rr68_interp_sa_ct(float64_t *, float64_t *, float64_t *, int32_t, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_rr68_interp_sa_ct(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** float64_t gsw_sa_from_rho(float64_t, float64_t, float64_t); */
	gsw_sa_from_rho(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sigma0(float64_t, float64_t); */
	gsw_sigma0(p0: number, p1: number): number;

	/** float64_t gsw_sigma0_ct_exact(float64_t, float64_t); */
	gsw_sigma0_ct_exact(p0: number, p1: number): number;

	/** float64_t gsw_sigma0_pt0_exact(float64_t, float64_t); */
	gsw_sigma0_pt0_exact(p0: number, p1: number): number;

	/** float64_t gsw_sigma1(float64_t, float64_t); */
	gsw_sigma1(p0: number, p1: number): number;

	/** float64_t gsw_sigma1_ct_exact(float64_t, float64_t); */
	gsw_sigma1_ct_exact(p0: number, p1: number): number;

	/** float64_t gsw_sigma2(float64_t, float64_t); */
	gsw_sigma2(p0: number, p1: number): number;

	/** float64_t gsw_sigma2_ct_exact(float64_t, float64_t); */
	gsw_sigma2_ct_exact(p0: number, p1: number): number;

	/** float64_t gsw_sigma3(float64_t, float64_t); */
	gsw_sigma3(p0: number, p1: number): number;

	/** float64_t gsw_sigma3_ct_exact(float64_t, float64_t); */
	gsw_sigma3_ct_exact(p0: number, p1: number): number;

	/** float64_t gsw_sigma4(float64_t, float64_t); */
	gsw_sigma4(p0: number, p1: number): number;

	/** float64_t gsw_sigma4_ct_exact(float64_t, float64_t); */
	gsw_sigma4_ct_exact(p0: number, p1: number): number;

	/** float64_t gsw_sound_speed(float64_t, float64_t, float64_t); */
	gsw_sound_speed(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sound_speed_ct_exact(float64_t, float64_t, float64_t); */
	gsw_sound_speed_ct_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sound_speed_t_exact(float64_t, float64_t, float64_t); */
	gsw_sound_speed_t_exact(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_sp_from_sk(float64_t); */
	gsw_sp_from_sk(p0: number): number;

	/** float64_t gsw_sp_from_sr(float64_t); */
	gsw_sp_from_sr(p0: number): number;

	/** float64_t gsw_sp_from_sstar(float64_t, float64_t, float64_t, float64_t); */
	gsw_sp_from_sstar(p0: number, p1: number, p2: number, p3: number): number;

	/** float64_t gsw_sp_salinometer(float64_t, float64_t); */
	gsw_sp_salinometer(p0: number, p1: number): number;

	/** void gsw_specvol_alpha_beta(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_specvol_alpha_beta(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_specvol_anom_standard(float64_t, float64_t, float64_t); */
	gsw_specvol_anom_standard(p0: number, p1: number, p2: number): number;

	/** void gsw_specvol_first_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_specvol_first_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** void gsw_specvol_first_derivatives_wrt_enthalpy(float64_t, float64_t, float64_t, float64_t *, float64_t *); */
	gsw_specvol_first_derivatives_wrt_enthalpy(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** void gsw_specvol_second_derivatives(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *, float64_t *, float64_t *); */
	gsw_specvol_second_derivatives(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number, p7: number): void;

	/** void gsw_specvol_second_derivatives_wrt_enthalpy(float64_t, float64_t, float64_t, float64_t *, float64_t *, float64_t *); */
	gsw_specvol_second_derivatives_wrt_enthalpy(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): void;

	/** float64_t gsw_spiciness0(float64_t, float64_t); */
	gsw_spiciness0(p0: number, p1: number): number;

	/** float64_t gsw_spiciness1(float64_t, float64_t); */
	gsw_spiciness1(p0: number, p1: number): number;

	/** float64_t gsw_spiciness2(float64_t, float64_t); */
	gsw_spiciness2(p0: number, p1: number): number;

	/** float64_t gsw_t_from_entropy(float64_t, float64_t, float64_t); */
	gsw_t_from_entropy(p0: number, p1: number, p2: number): number;

	/** void gsw_t_from_rho_exact(float64_t, float64_t, float64_t, float64_t &, float64_t &); */
	gsw_t_from_rho_exact(p0: number, p1: number, p2: number, p3: number, p4: number): void;

	/** float64_t gsw_thermobaric(float64_t, float64_t, float64_t); */
	gsw_thermobaric(p0: number, p1: number, p2: number): number;

	/** float64_t gsw_thermobaric_ct_exact(float64_t, float64_t, float64_t); */
	gsw_thermobaric_ct_exact(p0: number, p1: number, p2: number): number;

	/** void gsw_turner_rsubrho(float64_t *, float64_t *, float64_t *, int32_t, float64_t *, float64_t *, float64_t *); */
	gsw_turner_rsubrho(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number): void;

	/** float64_t * gsw_util_interp1q_int(int32_t, float64_t *, int32_t *, int32_t, float64_t *, float64_t *); */
	gsw_util_interp1q_int(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): number;

	/** float64_t * gsw_util_linear_interp(int32_t, float64_t *, int32_t, float64_t *, int32_t, float64_t *, float64_t *); */
	gsw_util_linear_interp(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number, p6: number): number;

	/** int32_t gsw_util_pchip_interp(float64_t *, float64_t *, int32_t, float64_t *, float64_t *, int32_t); */
	gsw_util_pchip_interp(p0: number, p1: number, p2: number, p3: number, p4: number, p5: number): number;
}
