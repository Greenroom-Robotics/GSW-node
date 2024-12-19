#include "TeosSea.h"

/*****************************************
  "TeosSea.cpp" Version 2.0
  by Randall Kent Whited
  rkwhited@gmail.com
  ---------------------------------------
  All copyrights and all license issues
  are the same as for previous versions
******************************************/

/*****************************************
  Most C function calls in the original
  TEOS-10 "gsw_oceanographic_toolbox.c"
  and "gsw_saar.c" source code file that
  have been used are in the base-class
  "TeosBase.cpp" or in the descendant
  classes "TeosIce.cpp" or "TeosSea.cpp".
  ---------------------------------------
  Methods used by both "TeosIce",
  "TeosSea" and "Teos_Matlab" are
  implemented in the "TeosBase" class.
*****************************************/

/** descendant class constructor : base class constructor*/
TeosSea::TeosSea(): TeosBase()
{
	try
	{
		statusOk = true;
	}
	catch(const unsigned e)
	{
		statusOk = !(e == badBaseAlloc);
	}
	catch (...)
	{
		statusOk = false;
	}
}

/** descendant class destructor */
TeosSea::~TeosSea()
{
}

/**************************************************************************
% gsw_adiabatic_lapse_rate_from_CT                     adiabatic lapse rate
%==========================================================================
%
% USAGE:
%  gsw_adiabatic_lapse_rate_from_CT(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the adiabatic lapse rate of sea water from Conservative
%  Temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  adiabatic_lapse_rate  =  adiabatic lapse rate                   [ K/Pa ]
%    Note.  The output is in unit of degress Celsius per Pa,
%      (or equivilently K/Pa) not in units of K/dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.22.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_adiabatic_lapse_rate_from_ct(double sa,double ct,double p)
{
	int n0=0, n1=1, n2=2;
	double pt0, pr0=0.0, t;

	pt0 = gsw_pt_from_ct(sa,ct);
	t = gsw_pt_from_t(sa,pt0,pr0,p);

	return (-gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n2,n0,sa,t,p));
}

/***************************************************************************
% gsw_alpha_on_beta                           alpha/beta (75-term equation)
%==========================================================================
%
% USAGE:
%    gsw_alpha_on_beta(SA,CT,p)
%
% DESCRIPTION:
%  Calculates alpha divided by beta, where alpha is the thermal expansion
%  coefficient and beta is the saline contraction coefficient of seawater
%  from Absolute Salinity and Conservative Temperature.  This function uses
%  the computationally-efficient expression for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  alpha_on_beta  =  thermal expansion coefficient with respect to
%                    Conservative Temperature divided by the saline
%                    contraction coefficient at constant Conservative
%                    Temperature                           [ kg g^-1 K^-1 ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2014: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_alpha_on_beta(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double xs, ys, z, v_ct_part, v_sa_part;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v_ct_part = gsw_gsvco_a(xs, ys, z);
	v_sa_part = gsw_gsvco_b(xs, ys, z);

	return (-(v_ct_part*xs)/(20.0*gtc.gsw_sfac*v_sa_part));
}

/**************************************************************************
% gsw_alpha_wrt_t_exact                       thermal expansion coefficient
%                                       with respect to in-situ temperature
%==========================================================================
%
% USAGE:
%  gsw_alpha_wrt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to
%  in-situ temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  alpha_wrt_t_exact  =  thermal expansion coefficient              [ 1/K ]
%                        with respect to in-situ temperature
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.18.1) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_alpha_wrt_t_exact(double sa,double t,double p)
{
	return (gsw_gibbs(0,1,1,sa,t,p)/gsw_gibbs(0,0,1,sa,t,p));
}

/***************************************************************************
% gsw_beta_const_t_exact                     saline contraction coefficient
%                                           at constant in-situ temperature
%==========================================================================
%
% USAGE:
%  gsw_beta_const_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the saline (i.e. haline) contraction coefficient of seawater
%  at constant in-situ temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  beta_const_t_exact  =  saline contraction coefficient           [ kg/g ]
%                         at constant in-situ temperature
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.19.1) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_beta_const_t_exact(double sa,double t,double p)
{
	return (-gsw_gibbs(1,0,1,sa,t,p)/gsw_gibbs(0,0,1,sa,t,p));
}

/***************************************************************************
% gsw_C_from_SP                                        conductivity from SP
%==========================================================================
%
% USAGE:
%  gsw_C_from_SP(SP,t,p)
%
% DESCRIPTION:
%  Calculates conductivity, C, from (SP,t,p) using PSS-78 in the range
%  2 < SP < 42.  If the input Practical Salinity is less than 2 then a
%  modified form of the Hill et al. (1986) fomula is used for Practical
%  Salinity.  The modification of the Hill et al. (1986) expression is to
%  ensure that it is exactly consistent with PSS-78 at SP = 2.
%
%  The conductivity ratio returned by this function is consistent with the
%  input value of Practical Salinity, SP, to 2x10^-14 psu over the full
%  range of input parameters (from pure fresh water up to SP = 42 psu).
%  This error of 2x10^-14 psu is machine precision at typical seawater
%  salinities.  This accuracy is achieved by having four different
%  polynomials for the starting value of Rtx (the square root of Rt) in
%  four different ranges of SP, and by using one and a half iterations of
%  a computationally efficient modified Newton-Raphson technique (McDougall
%  and Wotherspoon, 2013) to find the root of the equation.
%
%  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
%  Salinity in terms of the conductivity ratio, R, without actually
%  specifying the value of C(35,15,0) (which we currently take to be
%  42.9140 mS/cm).
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  C  =  conductivity                                             [ mS/cm ]
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Hill, K.D., T.M. Dauphinee and D.J. Woods, 1986: The extension of the
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   OE-11, 1, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix E of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Unesco, 1983: Algorithms for computation of fundamental properties of
%   seawater. Unesco Technical Papers in Marine Science, 44, 53 pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_c_from_sp(double sp,double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SP_COEFFICIENTS use gspc */
	double  p0 = 4.577801212923119e-3,      p1 = 1.924049429136640e-1,
			  p2 = 2.183871685127932e-5,      p3 = -7.292156330457999e-3,
			  p4 = 1.568129536470258e-4,      p5 = -1.478995271680869e-6,
			  p6 = 9.086442524716395e-4,      p7 = -1.949560839540487e-5,
			  p8 = -3.223058111118377e-6,     p9 = 1.175871639741131e-7,
			  p10 = -7.522895856600089e-5,    p11 = -2.254458513439107e-6,
			  p12 = 6.179992190192848e-7,     p13 = 1.005054226996868e-8,
			  p14 = -1.923745566122602e-9,    p15 = 2.259550611212616e-6,
			  p16 = 1.631749165091437e-7,     p17 = -5.931857989915256e-9,
			  p18 = -4.693392029005252e-9,    p19 = 2.571854839274148e-10,
			  p20 = 4.198786822861038e-12,
			  q0 = 5.540896868127855e-5,      q1 = 2.015419291097848e-1,
			  q2 = -1.445310045430192e-5,     q3 = -1.567047628411722e-2,
			  q4 = 2.464756294660119e-4,      q5 = -2.575458304732166e-7,
			  q6 = 5.071449842454419e-3,      q7 = 9.081985795339206e-5,
			  q8 = -3.635420818812898e-6,     q9 = 2.249490528450555e-8,
			  q10 = -1.143810377431888e-3,    q11 = 2.066112484281530e-5,
			  q12 = 7.482907137737503e-7,     q13 = 4.019321577844724e-8,
			  q14 = -5.755568141370501e-10,   q15 = 1.120748754429459e-4,
			  q16 = -2.420274029674485e-6,    q17 = -4.774829347564670e-8,
			  q18 = -4.279037686797859e-9,    q19 = -2.045829202713288e-10,
			  q20 = 5.025109163112005e-12,
			  s0 = 3.432285006604888e-3,      s1 = 1.672940491817403e-1,
			  s2 = 2.640304401023995e-5,      s3 = 1.082267090441036e-1,
			  s4 = -6.296778883666940e-5,     s5 = -4.542775152303671e-7,
			  s6 = -1.859711038699727e-1,     s7 = 7.659006320303959e-4,
			  s8 = -4.794661268817618e-7,     s9 = 8.093368602891911e-9,
			  s10 = 1.001140606840692e-1,     s11 = -1.038712945546608e-3,
			  s12 = -6.227915160991074e-6,    s13 = 2.798564479737090e-8,
			  s14 = -1.343623657549961e-10,   s15 = 1.024345179842964e-2,
			  s16 = 4.981135430579384e-4,     s17 = 4.466087528793912e-6,
			  s18 = 1.960872795577774e-8,     s19 = -2.723159418888634e-10,
			  s20 = 1.122200786423241e-12,
			  u0 = 5.180529787390576e-3,      u1 = 1.052097167201052e-3,
			  u2 = 3.666193708310848e-5,      u3 = 7.112223828976632e0,
			  u4 = -3.631366777096209e-4,     u5 = -7.336295318742821e-7,
			  u6 = -1.576886793288888e+2,     u7 = -1.840239113483083e-3,
			  u8 = 8.624279120240952e-6,      u9 = 1.233529799729501e-8,
			  u10 = 1.826482800939545e+3,     u11 = 1.633903983457674e-1,
			  u12 = -9.201096427222349e-5,    u13 = -9.187900959754842e-8,
			  u14 = -1.442010369809705e-10,   u15 = -8.542357182595853e+3,
			  u16 = -1.408635241899082e0,     u17 = 1.660164829963661e-4,
			  u18 = 6.797409608973845e-7,     u19 = 3.345074990451475e-10,
			  u20 = 8.285687652694768e-13;

	double t68, ft68, x, rtx=0.0, dsp_drtx, sqrty, part1, part2, hill_ratio,
								sp_est, rtx_old, rt, aa, bb, cc, dd, ee, ra,r, rt_lc, rtxm,
								sp_hill_raw;

	t68 = t*1.00024e0;
	ft68 = (t68 - 15e0)/(1e0 + gspc.k*(t68 - 15e0));
	x = sqrt(sp);
	/**
	-------------------------------------------------------------------------
	   Finding the starting value of Rtx, the square root of Rt, using four
	   different polynomials of SP and t68.
	--------------------------------------------------------------------------
	*/
	if (sp >= 9.0)
	{
		rtx = p0
				+ x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6
											 + p11*t68 + x*(p10 + p16*t68 + x*p15))))
				+ t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x
									 + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x
											  + p20*t68))));
	}
	else if (sp >= 0.25 && sp < 9.0)
	{
		rtx = q0
				+ x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6
											 + q11*t68 + x*(q10 + q16*t68 + x*q15))))
				+ t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x
									 + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x
											  + q20*t68))));
	}
	else if (sp >= 0.003 && sp < 0.25)
	{
		rtx = s0
				+ x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6
											 + s11*t68 + x*(s10 + s16*t68 + x*s15))))
				+ t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x
									 + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x
											  + s20*t68))));
	}
	else if (sp < 0.003)
	{
		rtx = u0
				+ x*(u1 + u4*t68 + x*(u3 + u7*t68 + x*(u6
											 + u11*t68 + x*(u10 + u16*t68 + x*u15))))
				+ t68*(u2+ t68*(u5 + x*x*(u12 + x*u17) + u8*x
									 + t68*(u9 + x*(u13 + x*u18)+ t68*(u14 + u19*x
											  + u20*t68))));
	}
	/**
	   -------------------------------------------------------------------------
	   Finding the starting value of dSP_dRtx, the derivative of SP with respect
	   to Rtx.
	   -------------------------------------------------------------------------
	*/
	dsp_drtx =  gspc.a1
					+ (2e0*gspc.a2 + (3e0*gspc.a3
											+ (4e0*gspc.a4 + 5e0*gspc.a5*rtx)*rtx)*rtx)*rtx
					+ ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3 + (4e0*gspc.b4
											 + 5e0*gspc.b5*rtx)*rtx)*rtx)*rtx);
	if (sp < 2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10.0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		dsp_drtx = dsp_drtx
					  + gspc.a0*800e0*rtx*(1.5e0 + 2e0*x)/(part1*part1)
					  + gspc.b0*ft68*(10e0 + sqrty*(20e0
											+ 30e0*sqrty))/(part2*part2);
		dsp_drtx = hill_ratio*dsp_drtx;
	}
	/**
	---------------------------------------------------------------------------
	   One iteration through the modified Newton-Raphson method (McDougall and
	   Wotherspoon, 2012) achieves an error in Practical Salinity of about
	   10^-12 for all combinations of the inputs.  One and a half iterations of
	   the modified Newton-Raphson method achevies a maximum error in terms of
	   Practical Salinity of better than 2x10^-14 everywhere.
	   We recommend one and a half iterations of the modified Newton-Raphson
	   method.
	   Begin the modified Newton-Raphson method.
	----------------------------------------------------------------------------
	*/
	sp_est = gspc.a0
				+ (gspc.a1 + (gspc.a2 + (gspc.a3
												 + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
				+ ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3
														+ (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
	if (sp_est <  2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		sp_est = hill_ratio*sp_hill_raw;
	}
	rtx_old = rtx;
	rtx = rtx_old - (sp_est - sp)/dsp_drtx;
	rtxm = 0.5e0*(rtx + rtx_old); /** This mean value of Rtx, Rtxm, is the
                  value of Rtx at which the derivative dSP_dRtx is evaluated. */
	dsp_drtx = gspc.a1
				  + (2e0*gspc.a2 + (3e0*gspc.a3 + (4e0*gspc.a4
										  + 5e0*gspc.a5*rtxm)*rtxm)*rtxm)*rtxm
				  + ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3
											+ (4e0*gspc.b4 + 5e0*gspc.b5*rtxm)*rtxm)*rtxm)*rtxm);
	if (sp_est <  2.0)
	{
		x = 400e0*(rtxm*rtxm);
		sqrty = 10e0*rtxm;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		dsp_drtx = dsp_drtx
					  + gspc.a0*800e0*rtxm*(1.5e0 + 2e0*x)/(part1*part1)
					  + gspc.b0*ft68*(10e0 + sqrty*(20e0
											+ 30e0*sqrty))/(part2*part2);
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		dsp_drtx = hill_ratio*dsp_drtx;
	}
	/**
	------------------------------------------------------------------------
	   The line below is where Rtx is updated at the end of the one full
	   iteration of the modified Newton-Raphson technique.
	-------------------------------------------------------------------------
	*/
	rtx = rtx_old - (sp_est - sp)/dsp_drtx;
	/**
	------------------------------------------------------------------------
	   Now we do another half iteration of the modified Newton-Raphson
	   technique, making a total of one and a half modified N-R iterations.
	-------------------------------------------------------------------------
	*/
	sp_est = gspc.a0
				+ (gspc.a1 + (gspc.a2 + (gspc.a3
            + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
				+ ft68*(gspc.b0 + (gspc.b1 + (gspc.b2 + (gspc.b3
            + (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);

	if (sp_est <  2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		sp_est = hill_ratio*sp_hill_raw;
	}
	rtx = rtx - (sp_est - sp)/dsp_drtx;
	/**
	----------------------------------------------------------------------------
	   Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
	----------------------------------------------------------------------------
	*/
	rt = rtx*rtx;
	aa = gspc.d3 + gspc.d4*t68;
	bb = 1e0 + t68*(gspc.d1 + gspc.d2*t68);
	cc = p*(gspc.e1 + p*(gspc.e2 + gspc.e3*p));
	/** rt_lc (i.e. rt_lower_case) corresponds to rt as defined in
	   the UNESCO 44 (1983) routines. */
	rt_lc = gspc.c0
			  + (gspc.c1 + (gspc.c2 + (gspc.c3
												+ gspc.c4*t68)*t68)*t68)*t68;
	dd = bb - aa*rt_lc*rt;
	ee = rt_lc*rt*aa*(bb + cc);
	ra = sqrt(dd*dd + 4e0*ee) - dd;
	r = 0.5e0*ra/aa;
	/**
	 The dimensionless conductivity ratio, R, is the conductivity input, C,
	   divided by the present estimate of C(SP=35, t_68=15, p=0) which is
	   42.9140 mS/cm (=4.29140 S/m^).
	*/
	return (gtc.gsw_c3515*r);
}

/***************************************************************************
% gsw_cp_t_exact                                     isobaric heat capacity
%==========================================================================
%
% USAGE:
%  gsw_cp_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the isobaric heat capacity of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  cp_t_exact  =  heat capacity of seawater                    [ J/(kg*K) ]
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_cp_t_exact(double sa,double t,double p)
{
	return (-(t+273.15e0)*gsw_gibbs(0,2,0,sa,t,p));
}

/***************************************************************************
% gsw_CT_from_entropy                         Conservative Temperature with
%                                                          entropy as input
% =========================================================================
%
% USAGE:
%   gsw_CT_from_entropy(SA,entropy)
%
% DESCRIPTION:
%  Calculates Conservative Temperature with entropy as an input variable.
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                   [ deg C ]
%
% OUTPUT:
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_ct_from_entropy(double sa,double entropy)
{
	double pt = gsw_pt_from_entropy(sa,entropy);

	return (gsw_ct_from_pt(sa,pt));
}

/***************************************************************************
% gsw_CT_first_derivatives    first derivatives of Conservative Temperature
%==========================================================================
%
% USAGE:
%  gsw_CT_first_derivatives(SA,pt)
%
% DESCRIPTION:
%  Calculates the following two derivatives of Conservative Temperature
%  (1) CT_SA, the derivative with respect to Absolute Salinity at
%      constant potential temperature (with pr = 0 dbar), and
%   2) CT_pt, the derivative with respect to potential temperature
%      (the regular potential temperature which is referenced to 0 dbar)
%      at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%         (whose reference pressure is 0 dbar)
%
% OUTPUT:
%  CT_SA  =  The derivative of Conservative Temperature with respect to
%            Absolute Salinity at constant potential temperature
%            (the regular potential temperature which has reference
%            sea pressure of 0 dbar).
%            The CT_SA output has units of:                     [ K/(g/kg)]
%  CT_pt  =  The derivative of Conservative Temperature with respect to
%            potential temperature (the regular one with pr = 0 dbar)
%            at constant SA. CT_pt is dimensionless.           [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.3a,b) and (A.15.8) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_ct_first_derivatives(double sa,double pt,double *ct_sa,
                                          double *ct_pt)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double abs_pt, g_sa_mod, g_sa_t_mod, x, y_pt;

	abs_pt = gtc.gsw_t0 + pt ;
	if (ct_pt != NULL)
		*ct_pt = -(abs_pt*gsw_gibbs_pt0_pt0(sa,pt))/gtc.gsw_cp0;

	if (ct_sa == NULL)
		return;

	x = sqrt(gtc.gsw_sfac*sa);
	y_pt = 0.025*pt;
	g_sa_t_mod = 1187.3715515697959
					 + x*(-1480.222530425046
                  + x*(2175.341332000392 + x*(-980.14153344888
                  + 220.542973797483*x) + y_pt*(-548.4580073635929
                  + y_pt*(592.4012338275047 + y_pt*(-274.2361238716608
                  + 49.9394019139016*y_pt)))) + y_pt*(-258.3988055868252
                  + y_pt*(-90.2046337756875 + y_pt*10.50720794170734)))
					   + y_pt*(3520.125411988816 + y_pt*(-1351.605895580406
                  + y_pt*(731.4083582010072 + y_pt*(-216.60324087531103
										  + 25.56203650166196*y_pt))));

	g_sa_t_mod = 0.5*gtc.gsw_sfac*0.025*g_sa_t_mod;

	g_sa_mod = 8645.36753595126
				   + x*(-7296.43987145382
               + x*(8103.20462414788 + y_pt*(2175.341332000392
               + y_pt*(-274.2290036817964 + y_pt*(197.4670779425016
               + y_pt*(-68.5590309679152 + 9.98788038278032*y_pt))))
               + x*(-5458.34205214835 - 980.14153344888*y_pt
               + x*(2247.60742726704 - 340.1237483177863*x
               + 220.542973797483*y_pt))) + y_pt*(-1480.222530425046
               + y_pt*(-129.1994027934126 + y_pt*(-30.0682112585625
               + y_pt*(2.626801985426835 ))))) + y_pt*(1187.3715515697959
               + y_pt*(1760.062705994408 + y_pt*(-450.535298526802
               + y_pt*(182.8520895502518 + y_pt*(-43.3206481750622
               + 4.26033941694366*y_pt)))));

	g_sa_mod = 0.5*gtc.gsw_sfac*g_sa_mod;

	*ct_sa = (g_sa_mod - abs_pt*g_sa_t_mod)/gtc.gsw_cp0;
}

/***************************************************************************
% gsw_entropy_from_pt                          specific entropy of seawater
%==========================================================================
%
% USAGE:
%  gsw_entropy_from_pt(SA,pt)
%
% DESCRIPTION:
%  Calculates specific entropy of seawater as a function of potential
%  temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%
% OUTPUT:
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.10 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_entropy_from_pt(double sa,double pt)
{
	return (-gsw_gibbs(0,1,0,sa,pt,0.0));
}

/***************************************************************************
% gsw_enthalpy_second_derivatives            second derivatives of enthalpy
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_enthalpy_second_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of specific
%  enthalpy (h),using the computationally-efficient expression for
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).
%   (1) h_SA_SA, second-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) h_SA_CT, second-order derivative with respect to SA & CT at
%       constant p.
%   (3) h_CT_CT, second-order derivative with respect to CT at constant SA
%       and p.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  h_SA_SA  =  The second derivative of specific enthalpy with respect to
%              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
%  h_SA_CT  =  The second derivative of specific enthalpy with respect to
%              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
%  h_CT_CT  =  The second derivative of specific enthalpy with respect to
%              CT at constant SA and p.                      [ J/(kg K^2) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic
%   variable for evaluating heat content and heat fluxes. Journal of
%   Physical Oceanography, 33, 945-963.
%    See Eqns. (18) and (22)
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_enthalpy_second_derivatives(double sa,double ct,double p,
                                                double *h_sa_sa,double *h_sa_ct,
                                                double *h_ct_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double dynamic_h_ct_ct_part, dynamic_h_sa_ct_part, dynamic_h_sa_sa_part, xs, xs2, ys, z;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;

	if (h_sa_sa != NULL)
	{
		xs2 = pow(xs,2);
		dynamic_h_sa_sa_part =
			z*(-gsvco.h101 + xs2*(3.0*gsvco.h301 + xs*(8.0*gsvco.h401
			+ xs*(15.0*gsvco.h501 + 24.0*gsvco.h601*xs))) + ys*(- gsvco.h111
			+ xs2*(3.0*gsvco.h311 + xs*(8.0*gsvco.h411 + 15.0*gsvco.h511*xs)) + ys*(-gsvco.h121
			+ xs2*(3.0*gsvco.h321 + 8.0*gsvco.h421*xs) + ys*(-gsvco.h131 + 3.0*gsvco.h331*xs2
			+ ys*(-gsvco.h141 - gsvco.h151*ys)))) + z*(-gsvco.h102 + xs2*(3.0*gsvco.h302
			+ xs*(8.0*gsvco.h402 + 15.0*gsvco.h502*xs)) + ys*(-gsvco.h112 + xs2*(3.0*gsvco.h312
			+ 8.0*gsvco.h412*xs) + ys*(-gsvco.h122 + 3.0*gsvco.h322*xs2
			+ ys*(-gsvco.h132 - gsvco.h142*ys ))) + z*(xs2*(8.0*gsvco.h403*xs + 3.0*gsvco.h313*ys)
			+ z*(-gsvco.h103 + 3.0*gsvco.h303*xs2 + ys*(-gsvco.h113 + ys*(-gsvco.h123 - gsvco.h133*ys))
			+ z*(-gsvco.h104 - gsvco.h114*ys - gsvco.h105*z)))));

		*h_sa_sa = 1e8*0.25*gtc.gsw_sfac*gtc.gsw_sfac*dynamic_h_sa_sa_part/pow(xs,3);
	}

	if (h_sa_ct != NULL)
	{
		dynamic_h_sa_ct_part =
			z*(gsvco.h111 + xs*(2.0*gsvco.h211 + xs*(3.0*gsvco.h311
			+ xs*(4.0*gsvco.h411 + 5.0*gsvco.h511*xs))) + ys*(2.0*gsvco.h121
			+ xs*(4.0*gsvco.h221 + xs*(6.0*gsvco.h321 + 8.0*gsvco.h421*xs))
			+ ys*(3.0*gsvco.h131 + xs*(6.0*gsvco.h231 + 9.0*gsvco.h331*xs)
			+ ys*(4.0*gsvco.h141 + 8.0*gsvco.h241*xs + 5.0*gsvco.h151*ys ))) + z*(gsvco.h112
			+ xs*(2.0*gsvco.h212 + xs*(3.0*gsvco.h312 + 4.0*gsvco.h412*xs))
			+ ys*(2.0*gsvco.h122 + xs*(4.0*gsvco.h222 + 6.0*gsvco.h322*xs)
			+ ys*(3.0*gsvco.h132 + 6.0*gsvco.h232*xs + 4.0*gsvco.h142*ys)) + z*(gsvco.h113
			+ xs*(2.0*gsvco.h213 + 3.0*gsvco.h313*xs) + ys*(2.0*gsvco.h123
			+ 4.0*gsvco.h223*xs + 3.0*gsvco.h133*ys) + gsvco.h114*z)));

		*h_sa_ct = 1e8*0.025*0.5*gtc.gsw_sfac*dynamic_h_sa_ct_part/xs;
	}

	if (h_ct_ct != NULL)
	{
		dynamic_h_ct_ct_part =
			z*(2.0*gsvco.h021 + xs*(2.0*gsvco.h121 + xs*(2.0*gsvco.h221
			+ xs*(2.0*gsvco.h321 + 2.0*gsvco.h421*xs))) + ys*(6.0*gsvco.h031
			+ xs*(6.0*gsvco.h131 + xs*(6.0*gsvco.h231 + 6.0*gsvco.h331*xs))
			+ ys*(12.0*gsvco.h041 + xs*(12.0*gsvco.h141 + 12.0*gsvco.h241*xs)
			+ ys*(20.0*gsvco.h051 + 20.0*gsvco.h151*xs + 30.0*gsvco.h061*ys)))
			+ z*(2.0*gsvco.h022 + xs*(2.0*gsvco.h122 + xs*(2.0*gsvco.h222
			+ 2.0*gsvco.h322*xs)) + ys*(6.0*gsvco.h032 + xs*(6.0*gsvco.h132
			+ 6.0*gsvco.h232*xs) + ys*(12.0*gsvco.h042 + 12.0*gsvco.h142*xs
			+ 20.0*gsvco.h052*ys)) + z*(2.0*gsvco.h023 + xs*(2.0*gsvco.h123
			+ 2.0*gsvco.h223*xs) + ys*(6.0*gsvco.h133*xs + 6.0*gsvco.h033
			+ 12.0*gsvco.h043*ys) + 2.0*gsvco.h024*z)));

		*h_ct_ct = 1e8*6.25e-4*dynamic_h_ct_ct_part;
	}
}

/**************************************************************************
% gsw_entropy_from_t                           specific entropy of seawater
%==========================================================================
%
% USAGE:
%  gsw_entropy_from_t(SA,t,p)
%
% DESCRIPTION:
%  Calculates specific entropy of seawater from in-situ temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_entropy_from_t(double sa,double t,double p)
{
	return (-gsw_gibbs(0,1,0,sa,t,p));
}

/***************************************************************************
% gsw_latentheat_evap_CT                         latent heat of evaporation
%==========================================================================
%
% USAGE:
%  gsw_latentheat_evap_CT(SA,CT)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of evaporation at p = 0 (the
%  surface).  It is defined as a function of Absolute Salinity, SA, and
%  Conservative Temperature, CT, and is valid in the ranges
%  0 < SA < 42 g/kg and 0 < CT < 40 deg C.  The errors range between
%  -0.4 and 0.6 J/kg.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  latentheat_evap = latent heat of evaporation                    [ J/kg ]
%
% AUTHOR:
%  Paul Barker, Trevor McDougall & Rainer Feistel      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.39 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_latentheat_evap_ct(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  c0  =  2.499065844825125e6, c1  = -1.544590633515099e-1,
			  c2  = -9.096800915831875e4, c3  =  1.665513670736000e2,
			  c4  =  4.589984751248335e1, c5  =  1.894281502222415e1,
			  c6  =  1.192559661490269e3, c7  = -6.631757848479068e3,
			  c8  = -1.104989199195898e2, c9  = -1.207006482532330e3,
			  c10 = -3.148710097513822e3, c11 =  7.437431482069087e2,
			  c12 =  2.519335841663499e3, c13 =  1.186568375570869e1,
			  c14 =  5.731307337366114e2, c15 =  1.213387273240204e3,
			  c16 =  1.062383995581363e3, c17 = -6.399956483223386e2,
			  c18 = -1.541083032068263e3, c19 =  8.460780175632090e1,
			  c20 = -3.233571307223379e2, c21 = -2.031538422351553e2,
			  c22 =  4.351585544019463e1, c23 = -8.062279018001309e2,
			  c24 =  7.510134932437941e2, c25 =  1.797443329095446e2,
			  c26 = -2.389853928747630e1, c27 =  1.021046205356775e2;

	double x, y;

	x = sqrt(gtc.gsw_sfac*sa);
	y = ct/40.0;

	return (c0 + x*(c1 + c4*y + x*(c3
				+ y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))
				+ x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))
			   + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)
				+ y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y))))));
}

/***************************************************************************
% gsw_latentheat_evap_t                          latent heat of evaporation
%==========================================================================
%
% USAGE:
%  gsw_latentheat_evap_t(SA,t)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of evaporation at p = 0 (the
%  surface).  It is defined as a function of Absolute Salinity, SA, and
%  in-situ temperature, t, and is valid in the ranges 0 < SA < 40 g/kg
%  and 0 < CT < 42 deg C. The errors range between -0.4 and 0.6 J/kg.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% OUTPUT:
%  latentheat_evap = latent heat of evaporation                    [ J/kg ]
%
% AUTHOR:
%  Paul Barker, Trevor McDougall & Rainer Feistel      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.39 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_latentheat_evap_t(double sa,double t)
{
	double  ct = gsw_ct_from_pt(sa,t);

	return (gsw_latentheat_evap_ct(sa,ct));
}

/**************************************************************************
% gsw_pot_rho_t_exact                                     potential density
%==========================================================================
%
% USAGE:
%  gsw_pot_rho_t_exact(SA,t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential density of seawater.  Note. This function outputs
%  potential density, not potential density anomaly; that is, 1000 kg/m^3
%  is not subtracted.
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  t      =  in-situ temperature (ITS-90)                         [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref  =  reference pressure                                    [ dbar ]
%            ( i.e. reference absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_rho_t_exact  =  potential density (not potential density anomaly)
%                                                                [ kg/m^3 ]
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.4 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_pot_rho_t_exact(double sa,double t,double p,double p_ref)
{
	double pt = gsw_pt_from_t(sa,t,p,p_ref);

	return (gsw_rho_t_exact(sa,pt,p_ref));
}

/**************************************************************************
function pt = gsw_pt_from_entropy(SA,entropy)

% gsw_pt_from_entropy                          potential temperature with a
%                                       reference sea pressure of zero dbar
%                                                  as a function of entropy
% =========================================================================
%
% USAGE:
%  gsw_pt_from_entropy(SA,entropy)
%
% DESCRIPTION:
%  Calculates potential temperature with reference pressure p_ref = 0 dbar
%  and with entropy as an input variable.
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                   [ deg C ]
%
% OUTPUT:
%  pt   =  potential temperature                                  [ deg C ]
%          with reference sea pressure (p_ref) = 0 dbar.
%  Note. The reference sea pressure of the output, pt, is zero dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_pt_from_entropy(double sa,double entropy)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int number_of_iterations;
	double c, dentropy, dentropy_dt, ent_sa, part1, part2, pt, ptm, pt_old;

	/**Find the initial value of pt*/
	part1 = 1.0 - sa/gtc.gsw_sso;
	part2 = 1.0 - 0.05*part1;
	ent_sa = (gtc.gsw_cp0/gtc.gsw_t0)*part1*(1.0 - 1.01*part1);
	c = (entropy - ent_sa)*(part2/gtc.gsw_cp0);
	pt = gtc.gsw_t0*(exp(c) - 1.0);
	dentropy_dt = gtc.gsw_cp0/((gtc.gsw_t0 + pt)*part2);
	for (number_of_iterations = 1; number_of_iterations <= 2;
			number_of_iterations++)
	{
		pt_old = pt;
		dentropy = gsw_entropy_from_pt(sa,pt_old) - entropy;
		pt = pt_old - dentropy/dentropy_dt;
		ptm = 0.5*(pt + pt_old);
		dentropy_dt = -gsw_gibbs_pt0_pt0(sa,ptm);
		pt = pt_old - dentropy/dentropy_dt;
	}
	/**
	  Maximum error of 2.2x10^-6 degrees C for one iteration.
	  Maximum error is 1.4x10^-14 degrees C for two iterations
	  (two iterations is the default, "for Number_of_iterations = 1:2").
	*/
	return (pt);
}

/***************************************************************************
% gsw_rho_first_derivatives_wrt_enthalpy                  first derivatives
%                                  specific volume with respect to enthalpy
% =========================================================================
%
% USAGE:
%  gsw_rho_first_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of specific
%  volume (v),
%   (1) rho_SA, first-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) rho_h, first-order derivative with respect to SA & CT at
%       constant p.
%
%  Note that this function uses the using the computationally-efficient
%  expression for specific volume (Roquet et al., 2015).  There is an
%  alternative to calling this function, namely
%  gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses
%  the full Gibbs function (IOC et al., 2010).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_SA =  The first derivative of rho with respect to
%              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
%  rho_h  =  The first derivative of rho with respect to
%              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/

void TeosSea::gsw_rho_first_derivatives_wrt_enthalpy(double sa,double ct,
                                                      double p,double *rho_sa,
                                                      double *rho_h)
{
	double rec_v2, v_h=0.0, v_sa;

	if ((rho_sa != NULL) && (rho_h != NULL))
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,&v_h);
	}
	else if (rho_sa != NULL)
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,NULL);
	}
	else if (rho_h != NULL)
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,NULL,&v_h);
	}

	rec_v2 = pow(1.0/gsw_specvol(sa,ct,p), 2);

	if (rho_sa != NULL) *rho_sa = -v_sa*rec_v2;

	if (rho_h != NULL) *rho_h = -v_h*rec_v2;
}

/***************************************************************************
% gsw_rho_t_exact                                       density of seawater
%==========================================================================
%
% USAGE:
%  gsw_rho_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates in-situ density of seawater from Absolute Salinity and
%  in-situ temperature.  Note that the output, rho, is density,
%  not density anomaly; that is, 1000 kg/m^3 is not subracted from it.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_t_exact  =  in-situ density (not density anomaly)         [ kg/m^3 ]
%
% AUTHOR:
%  Paul Barker, David Jackett & Trevor McDougall       [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.8 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_rho_t_exact(double sa,double t,double p)
{
	return (1.0/gsw_gibbs(0,0,1,sa,t,p));
}

/**************************************************************************
% gsw_SA_from_rho                            Absolute Salinity from density
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_SA_from_rho(rho,CT,p)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of a seawater sample, for given values
%  of its density, Conservative Temperature and sea pressure (in dbar).
%  This function uses the computationally-efficient 75-term expression for
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  SA  =  Absolute Salinity.                                       [ g/kg ]
%   Note. This is expressed on the Reference-Composition Salinity
%     Scale of Millero et al. (2008).
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.5 of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008:
%   The composition of Standard Seawater and the definition of the
%   Reference-Composition Salinity Scale. Deep-Sea Res. I, 55, 50-72.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sa_from_rho(double rho,double ct,double p)
{
	int no_iter;
	double sa, v_lab, v_0, v_50, v_sa, sa_old, delta_v, sa_mean;

	v_lab = 1.0/rho;
	v_0 = gsw_specvol(0.0,ct,p);
	v_50 = gsw_specvol(50.0,ct,p);
	sa = 50.0*(v_lab - v_0)/(v_50 - v_0);
	if (sa < 0.0 || sa > 50.0)
		return (cppGSW_INVALID_VALUE);
	v_sa = (v_50 - v_0)/50.0;
	for (no_iter=1; no_iter <= 2; no_iter++)
	{
		sa_old = sa;
		delta_v = gsw_specvol(sa_old,ct,p) - v_lab;
		sa = sa_old - delta_v/v_sa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_specvol_first_derivatives(sa_mean,ct,p,&v_sa,NULL,NULL);
		sa = sa_old - delta_v/v_sa;

		if (sa < 0.0 || sa > 50.0)
			return (cppGSW_INVALID_VALUE);
	}
	return (sa);
}

/***************************************************************************
% gsw_sigma0                       potential density anomaly with reference
%                                 sea pressure of 0 dbar (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sigma0(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 0 dbar,
%  this being this particular potential density minus 1000 kg/m^3.  This
%  function has inputs of Absolute Salinity and Conservative Temperature.
%  This function uses the computationally-efficient expression for
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma0  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 0 dbar,
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma0(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double vp0, xs, ys;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	vp0 = gsvco.v000
			+ xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030
			+ xs*(gsvco.v040 + xs*(gsvco.v050 + gsvco.v060*xs)))))
			+ ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120
			+ xs*(gsvco.v130 + xs*(gsvco.v140 + gsvco.v150*xs))))
			+ ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220
			+ xs*(gsvco.v230 + gsvco.v240*xs))) + ys*(gsvco.v300
			+ xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs))
			+ ys*(gsvco.v400 + xs*(gsvco.v410 + gsvco.v420*xs)
			+ ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys)))));

	return (1.0/vp0 - 1000.0);
}

/***************************************************************************
% gsw_sound_speed                            sound speed (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sound_speed(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the speed of sound in seawater.  This function has inputs of
%  Absolute Salinity and Conservative Temperature.  This function uses the
%  computationally-efficient expression for specific volume in terms of SA,
%  CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  sound_speed  =  speed of sound in seawater                       [ m/s ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.17.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sound_speed(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double v, v_p, xs, ys, z;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	v_p = gsw_gsvco_c(xs, ys, z);

	return (10000.0*sqrt(-v*v/v_p));
}

/***************************************************************************
% gsw_sound_speed_t_exact                                       sound speed
%==========================================================================
%
% USAGE:
%  sound_speed_t_exact = gsw_sound_speed_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the speed of sound in seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  sound_speed_t_exact  =  speed of sound in seawater               [ m/s ]
%
% AUTHOR:
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.17.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sound_speed_t_exact(double sa,double t,double p)
{
	double g_tt, g_tp;

	g_tt = gsw_gibbs(0,2,0,sa,t,p);
	g_tp = gsw_gibbs(0,1,1,sa,t,p);

	return (gsw_gibbs(0,0,1,sa,t,p) *
			  sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(0,0,2,sa,t,p))));
}

/***************************************************************************
% gsw_SP_from_SK                   Practical Salinity from Knudsen Salinity
%==========================================================================
%
% USAGE:
%  gsw_SP_from_SK(SK)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Knudsen Salinity.
%
% INPUT:
%  SK  =  Knudsen Salinity                        [parts per thousand, ppt]
%
% OUTPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%      See Appendix A.3 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sp_from_sk(double sk)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double gsw_sp_from_sk_value = (sk - 0.03e0)*(gtc.gsw_soncl/1.805e0);

	if (gsw_sp_from_sk_value < 0e0)
		gsw_sp_from_sk_value = cppGSW_INVALID_VALUE;

	return (gsw_sp_from_sk_value);
}

/**************************************************************************
% gsw_SP_from_SR                 Practical Salinity from Reference Salinity
%==========================================================================
%
% USAGE:
%  gsw_SP_from_SR(SR)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Reference Salinity.
%
% INPUT:
%  SR  =  Reference Salinity                                       [ g/kg ]
%
% OUTPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sp_from_sr(double sr)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return(sr/gtc.gsw_ups);
}

/***************************************************************************
% gsw_SP_salinometer       Practical Salinity from a laboratory salinometer
%==========================================================================
%
% USAGE:
%  gsw_SP_salinometer(Rt,t)
%
% DESCRIPTION:
%  Calculates Practical Salinity SP from a salinometer, primarily using the
%  PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity
%  is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
%  produces a Practical Salinity that is less than 2 then the Practical
%  Salinity is recalculated with a modified form of the Hill et al. (1986)
%  formula.  The modification of the Hill et al. (1986) expression is to
%  ensure that it is exactly consistent with PSS-78 at SP = 2.
%
%  A laboratory salinometer has the ratio of conductivities, Rt, as an
%  output, and the present function uses this conductivity ratio and the
%  temperature t of the salinometer bath as the two input variables.
%
% INPUT:
%  Rt  =  C(SP,t_68,0)/C(SP=35,t_68,0)                         [ unitless ]
%  t   =  temperature of the bath of the salinometer,
%         measured on the ITS-90 scale (ITS-90)                   [ deg C ]
%
% OUTPUT:
%  SP  =  Practical Salinity on the PSS-78 scale               [ unitless ]
%
% AUTHOR:
%  Paul Barker, Trevor McDougall and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation of
%   fundamental properties of seawater. Unesco Tech. Pap. in Mar. Sci., 44,
%   53 pp.
%
%  Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   11, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See appendix E of this TEOS-10 Manual, and in particular,
%     Eqns. (E.2.1) and (E.2.6).
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sp_salinometer(double rt,double t)
{
	/** for GSW_SP_COEFFICIENTS use gspc */
	double t68, ft68, rtx, sp, hill_ratio, x, sqrty, part1, part2, sp_hill_raw;

	if (rt < 0)
	{
		return NAN;
	}
	t68 = t*1.00024;
	ft68 = (t68 - 15)/(1 + gspc.k*(t68 - 15));
	rtx = sqrt(rt);
	sp = gspc.a0
      + (gspc.a1 + (gspc.a2 + (gspc.a3
      + (gspc.a4 + gspc.a5 * rtx) * rtx) * rtx) * rtx) * rtx
      + ft68 * (gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3
		+ (gspc.b4 + gspc.b5 * rtx) * rtx) * rtx) * rtx) * rtx);

	/**
	   The following section of the code is designed for SP < 2 based on the
	   Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
	   exactly equal to the PSS-78 algorithm at SP = 2.
	*/
	if (sp < 2)
	{
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		x = 400e0*rt;
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp - gspc.a0/part1 - gspc.b0*ft68/part2;
		sp = hill_ratio*sp_hill_raw;
	}

	return sp;
}

/***************************************************************************
% gsw_specvol_anom_standard                    specific volume anomaly with
%                         reference of SA = SSO & CT = 0 (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_specvol_anom_standard(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, Conservative
%  Temperature and pressure. It uses the computationally-efficient
%  expression for specific volume as a function of SA, CT and p (Roquet
%  et al., 2015).  The reference value to which the anomally is calculated
%  has an Absolute Salinity of SSO and Conservative Temperature equal to
%  0 degress C.
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol_anom  =  specific volume anomaly                      [ m^3/kg ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_specvol_anom_standard(double sa,double ct,double p)
{
	return (gsw_specvol(sa,ct,p) - gsw_specvol_sso_0(p));
}

/***************************************************************************
% gsw_specvol_first_derivatives                     first order derivatives
%                                     of specific volume (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_specvol_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) v_CT, first-order derivative with respect to CT at
%       constant SA & p.
%   (3) v_P, first-order derivative with respect to P at constant SA
%       and CT.
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_specvol_first_derivatives_CT_exact(SA,CT,p) which uses the full
%  Gibbs function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2010).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  v_SA  =  The first derivative of specific volume with respect to
%           Absolute Salinity at constant CT & p.     [ (m^3/kg)(g/kg)^-1 ]
%  v_CT  =  The first derivative of specific volume with respect to
%           CT at constant SA and p.                         [ m^3/(K kg) ]
%  v_P   =  The first derivative of specific volume with respect to
%           P at constant SA and CT.                        [ m^3/(Pa kg) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_specvol_first_derivatives(double sa,double ct,double p,
                                             double *v_sa,double *v_ct,double *v_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double  v_ct_part, v_p_part, v_sa_part, xs, ys, z;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;

	if (v_sa != NULL)
	{
		v_sa_part = gsw_gsvco_b(xs, ys, z);
		*v_sa = 0.5*gtc.gsw_sfac*v_sa_part/xs;
	}

	if (v_ct != NULL)
	{
		v_ct_part =gsw_gsvco_a(xs, ys, z);
		*v_ct = 0.025*v_ct_part;
	}

	if (v_p != NULL)
	{
		v_p_part = gsw_gsvco_c(xs, ys, z);
		*v_p = 1e-8*v_p_part;
	}
}

/***************************************************************************
% gsw_specvol_first_derivatives_wrt_enthalpy        first order derivatives
%                               of specific volume with respect to enthalpy
% =========================================================================
%
% USAGE:
%  gsw_specvol_first_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of specific
%  volume (v),
%   (1) v_SA_wrt_h, first-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) v_h, first-order derivative with respect to h at
%       constant SA & p.
%
%  Note that this function uses the using the computationally-efficient
%  75 term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses
%  the full Gibbs function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2010).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  v_SA_wrt_h  =  The first derivative of specific volume with respect to
%              Absolute Salinity at constant CT & p.
%                                           [ (m^3/kg)(g/kg)^-1 (J/kg)^-1 ]
%  v_h  =  The first derivative of specific volume with respect to
%              SA and CT at constant p.               [ (m^3/kg)(J/kg)^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_specvol_first_derivatives_wrt_enthalpy(double sa,double ct,
                                                            double p,double *v_sa,
                                                            double *v_h)
{
	double h_ct=1.0, h_sa, rec_h_ct, vct_ct, vct_sa;

	if (v_sa != NULL)
	{
		gsw_specvol_first_derivatives(sa,ct,p,&vct_sa,&vct_ct,NULL);
		gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);
	}
	else if (v_h != NULL)
	{
		gsw_specvol_first_derivatives(sa,ct,p,NULL,&vct_ct,NULL);
		gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);
	}

	rec_h_ct = 1.0/h_ct;

	if (v_sa != NULL)
		*v_sa = vct_sa - (vct_ct*h_sa)*rec_h_ct;

	if (v_h != NULL)
		*v_h = vct_ct*rec_h_ct;
	return;
}

/***************************************************************************
% gsw_specvol_second_derivatives                   second order derivatives
%                                     of specific volume (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_specvol_second_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of specific
%  volume (v),
%   (1) v_SA_SA, second-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) v_SA_CT, second-order derivative with respect to SA & CT at
%       constant p.
%   (3) v_CT_CT, second-order derivative with respect to CT at constant SA
%       and p.
%   (4) v_SA_P, second-order derivative with respect to SA & P at
%       constant CT.
%   (5) v_CT_P, second-order derivative with respect to CT & P at
%       constant SA.
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_specvol_second_derivatives_CT_exact(SA,CT,p) which uses the full
%  Gibbs function (IOC et al., 2010).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  v_SA_SA  =  The second derivative of specific volume with respect to
%              Absolute Salinity at constant CT & p.  [ (m^3/kg)(g/kg)^-2 ]
%  v_SA_CT  =  The second derivative of specific volume with respect to
%              SA and CT at constant p.           [ (m^3/kg)(g/kg)^-1 K^-1]
%  v_CT_CT  =  The second derivative of specific volume with respect to
%              CT at constant SA and p.                  [ (m^3/kg) K^-2) ]
%  v_SA_P  =  The second derivative of specific volume with respect to
%              SA and P at constant CT.                  [ (m^3/kg) Pa^-1 ]
%  v_CT_P  =  The second derivative of specific volume with respect to
%              CT and P at constant SA.             [ (m^3/kg) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_specvol_second_derivatives(double sa,double ct,double p,
		                                          double *v_sa_sa,double *v_sa_ct,
		                                          double *v_ct_ct,double *v_sa_p,
		                                          double *v_ct_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double v_ct_ct_part, v_ct_p_part, v_sa_ct_part, v_sa_p_part,
			 v_sa_sa_part, xs, xs2, ys, z;

	xs2 = gtc.gsw_sfac*sa + gtc.offset;
	xs = sqrt(xs2);
	ys = ct*0.025;
	z = p*1e-4;

	if (v_sa_sa != NULL)
	{
		v_sa_sa_part = (-gsvco.b000
                     + xs2*(gsvco.b200 + xs*(2.0*gsvco.b300 + xs*(3.0*gsvco.b400
							+ 4.0*gsvco.b500*xs))) + ys*(-gsvco.b010 + xs2*(gsvco.b210 + xs*(2.0*gsvco.b310
							+ 3.0*gsvco.b410*xs)) + ys*(-gsvco.b020 + xs2*(gsvco.b220 + 2.0*gsvco.b320*xs)
							+ ys*(-gsvco.b030 + gsvco.b230*xs2 + ys*(-gsvco.b040 - gsvco.b050*ys)))) + z*(-gsvco.b001
							+ xs2*(gsvco.b201 + xs*(2.0*gsvco.b301 + 3.0*gsvco.b401*xs)) + ys*(-gsvco.b011
							+ xs2*(gsvco.b211 + 2.0*gsvco.b311*xs) + ys*(-gsvco.b021 + gsvco.b221*xs2
							+ ys*(-gsvco.b031 - gsvco.b041*ys))) + z*(-gsvco.b002 + xs2*(gsvco.b202 + 2.0*gsvco.b302*xs)
							+ ys*(-gsvco.b012 + gsvco.b212*xs2 + ys*(-gsvco.b022 - gsvco.b032*ys))
							+ z*(-gsvco.b003 - gsvco.b013*ys - gsvco.b004*z))))/xs2;

		*v_sa_sa = 0.25*gtc.gsw_sfac*gtc.gsw_sfac*v_sa_sa_part/xs; /** "/x2 added" */
	}
	if (v_sa_ct != NULL)
	{
		v_sa_ct_part = (gsvco.b010
							+ xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310 + gsvco.b410*xs)))
							+ ys*(2.0*(gsvco.b020 + xs*(gsvco.b120 + xs*(gsvco.b220 + gsvco.b320*xs)))
							+ ys*(3.0*(gsvco.b030 + xs*(gsvco.b130 + gsvco.b230*xs)) + ys*(4.0*(gsvco.b040
							+ gsvco.b140*xs) + 5.0*gsvco.b050*ys))) + z*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211
							+ gsvco.b311*xs)) + ys*(2.0*(gsvco.b021 + xs*(gsvco.b121 + gsvco.b221*xs))
							+ ys*(3.0*(gsvco.b031 + gsvco.b131*xs) + 4.0*gsvco.b041*ys)) + z*(gsvco.b012
							+ xs*(gsvco.b112 + gsvco.b212*xs) + ys*(2.0*(gsvco.b022 + gsvco.b122*xs)
							+ 3.0*gsvco.b032*ys) + gsvco.b013*z)))/xs;

		*v_sa_ct = 0.025*0.5*gtc.gsw_sfac*v_sa_ct_part;
	}
	if (v_ct_ct != NULL)
	{
		v_ct_ct_part = gsvco.a010
							+ xs*(gsvco.a110 + xs*(gsvco.a210 + xs*(gsvco.a310 + gsvco.a410*xs)))
							+ ys*(2.0*(gsvco.a020 + xs*(gsvco.a120 + xs*(gsvco.a220 + gsvco.a320*xs)))
							+ ys*(3.0*(gsvco.a030 + xs*(gsvco.a130 + gsvco.a230*xs)) + ys*(4.0*(gsvco.a040
							+ gsvco.a140*xs) + 5.0*gsvco.a050*ys))) + z*(gsvco.a011 + xs*(gsvco.a111 + xs*(gsvco.a211
							+ gsvco.a311*xs)) + ys*(2.0*(gsvco.a021 + xs*(gsvco.a121 + gsvco.a221*xs))
							+ ys*(3.0*(gsvco.a031 + gsvco.a131*xs) + 4.0*gsvco.a041*ys)) + z*(gsvco.a012
							+ xs*(gsvco.a112 + gsvco.a212*xs) + ys*(2.0*(gsvco.a022 + gsvco.a122*xs)
							+ 3.0*gsvco.a032*ys) + gsvco.a013*z));

		*v_ct_ct = 0.025*0.025*v_ct_ct_part;
	}
	if (v_sa_p != NULL)
	{
		v_sa_p_part = gsvco.b001
						+ xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301
						+ gsvco.b401*xs))) + ys*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211
						+ gsvco.b311*xs)) + ys*(gsvco.b021 + xs*(gsvco.b121 + gsvco.b221*xs)
						+ ys*(gsvco.b031 + gsvco.b131*xs + gsvco.b041*ys))) + z*(2.0*(gsvco.b002 + xs*(gsvco.b102
						+ xs*(gsvco.b202 + gsvco.b302*xs)) + ys*(gsvco.b012 + xs*(gsvco.b112
						+ gsvco.b212*xs) + ys*(gsvco.b022 + gsvco.b122*xs + gsvco.b032*ys)))
						+ z*(3.0*(gsvco.b003 + gsvco.b103*xs + gsvco.b013*ys) + 4.0*gsvco.b004*z));

		*v_sa_p = 1e-8*0.5*gtc.gsw_sfac*v_sa_p_part/xs;
	}
	if (v_ct_p != NULL)
	{
		v_ct_p_part = gsvco.a001
						   + xs*(gsvco.a101 + xs*(gsvco.a201 + xs*(gsvco.a301
							+ gsvco.a401*xs))) + ys*(gsvco.a011
							+ xs*(gsvco.a111 + xs*(gsvco.a211 + gsvco.a311*xs)) + ys*(gsvco.a021
							+ xs*(gsvco.a121 + gsvco.a221*xs) + ys*(gsvco.a031 + gsvco.a131*xs
							+ gsvco.a041*ys))) + z*(2.0*(gsvco.a002 + xs*(gsvco.a102
							+ xs*(gsvco.a202 + gsvco.a302*xs)) + ys*(gsvco.a012 + xs*(gsvco.a112 + gsvco.a212*xs)
							+ ys*(gsvco.a022 + gsvco.a122*xs + gsvco.a032*ys))) + z*(3.0*(gsvco.a003
							+ gsvco.a103*xs + gsvco.a013*ys) + 4.0*gsvco.a004*z));

		*v_ct_p = 1e-8*0.025*v_ct_p_part;
	}
}

/***************************************************************************
% gsw_specvol_second_derivatives_wrt_enthalpy      second order derivatives
%                               of volume specific with respect to enthalpy
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_specvol_second_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v) with respect to enthalpy,
%   (1) v_SA_SA_wrt_h, second-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) v_SA_h, second-order derivative with respect to SA & h at
%       constant p.
%   (3) v_h_h, second-order derivative with respect to h at
%       constant SA & p.
%
%  Note that this function uses the using the computationally-efficient
%  75 term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses
%  the full Gibbs function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2010).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  v_SA_SA_wrt_h = The second-order derivative of specific volume with
%                  respect to Absolute Salinity at constant h & p.
%                                           [ (m^3/kg)(g/kg)^-2 (J/kg)^-1 ]
%  v_SA_h  = The second-order derivative of specific volume with respect to
%            SA and h at constant p.        [ (m^3/kg)(g/kg)^-1 (J/kg)^-1 ]
%  v_h_h   = The second-order derivative with respect to h at
%            constant SA & p.                         [ (m^3/kg)(J/kg)^-2 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_specvol_second_derivatives_wrt_enthalpy(double sa,double ct,
                                                            double p,
		                                                      double *v_sa_sa,
		                                                      double *v_sa_h,
		                                                      double *v_h_h)
{
	double h_ct, h_ct_ct, h_sa, h_sa_ct, h_sa_sa, rec_h_ct, v_h_h_part,
			  rec_h_ct2, v_ct, vct_ct_ct, vct_sa_ct, vct_sa_sa, v_sa_h_part;

	gsw_specvol_first_derivatives(sa,ct,p,NULL, &v_ct, NULL);

	if ((v_sa_sa != NULL) || (v_sa_h != NULL))
		gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);
	else
		gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

	if (v_sa_sa != NULL)
		gsw_specvol_second_derivatives(sa,ct,p,&vct_sa_sa,&vct_sa_ct,
												 &vct_ct_ct, NULL, NULL);
	else if (v_sa_h != NULL)
		gsw_specvol_second_derivatives(sa,ct,p,NULL,&vct_sa_ct,&vct_ct_ct,
												 NULL, NULL);
	else
		gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,&vct_ct_ct,
												 NULL, NULL);

	if (v_sa_sa != NULL)
		gsw_enthalpy_second_derivatives(sa,ct,p,&h_sa_sa,&h_sa_ct,&h_ct_ct);
	else if (v_sa_h != NULL)
		gsw_enthalpy_second_derivatives(sa,ct,p,NULL,&h_sa_ct,&h_ct_ct);
	else
		gsw_enthalpy_second_derivatives(sa,ct,p,NULL,NULL,&h_ct_ct);

	rec_h_ct = 1.0/h_ct;
	rec_h_ct2 = rec_h_ct*rec_h_ct;

	v_h_h_part = (vct_ct_ct*h_ct - h_ct_ct*v_ct)*(rec_h_ct2*rec_h_ct);

	if (v_h_h != NULL) *v_h_h = v_h_h_part;

	if ((v_sa_sa != NULL) || (v_sa_h != NULL))
	{
		v_sa_h_part = (vct_sa_ct*h_ct - v_ct*h_sa_ct)*rec_h_ct2
						  - h_sa*v_h_h_part;

		if (v_sa_h != NULL) *v_sa_h = v_sa_h_part;

		if (v_sa_sa != NULL)
			*v_sa_sa = vct_sa_sa - (h_ct*(vct_sa_ct*h_sa
            - v_ct*h_sa_sa) + v_ct*h_sa*h_sa_ct)*rec_h_ct2
            - h_sa*v_sa_h_part;
	}
}

/***************************************************************************
% gsw_Turner_Rsubrho              Turner angle & Rsubrho (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_Turner_Rsubrho(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the Turner angle and the Rsubrho as a function of pressure
%  down a vertical water column.  These quantities express the relative
%  contributions of the vertical gradients of Conservative Temperature
%  and Absolute Salinity to the vertical stability (the square of the
%  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at
%  the mid pressure between the individual data points in the vertical.
%  This function uses computationally-efficient 75-term expression for
%  specific volume in terms of SA, CT and p (Roquet et al., 2015).  Note
%  that in the double-diffusive literature, papers concerned with the
%  "diffusive" form of double-diffusive convection often define the
%  stability ratio as the reciprocal of what is defined here as the
%  stability ratio.
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  Tu       =  Turner angle, on the same (M-1)xN grid as p_mid.
%              Turner angle has units of:           [ degrees of rotation ]
%  Rsubrho  =  Stability Ratio, on the same (M-1)xN grid as p_mid.
%              Rsubrho is dimensionless.                       [ unitless ]
%  p_mid    =  mid pressure between the indivual points of the p grid.
%              That is, p_mid is on a (M-1)xN grid in the vertical.
%              p_mid has units of:                                 [ dbar ]
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.15.1) and (3.16.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_turner_rsubrho(double *sa,double *ct,double *p,int nz,
                                    double *tu,double *rsubrho,double *p_mid)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int k;
	double dsa, sa_mid, dct, ct_mid, alpha_mid, beta_mid;

	if (nz < 2)
   {
      return;
   }

	for (k = 0; k < nz-1; k++)
	{
		dsa         = (sa[k] - sa[k+1]);
		sa_mid      = 0.5e0*(sa[k] + sa[k+1]);
		dct         = (ct[k] - ct[k+1]);
		ct_mid      = 0.5e0*(ct[k] + ct[k+1]);
		p_mid[k]    = 0.5e0*(p[k] + p[k+1]);

		gsw_specvol_alpha_beta(sa_mid,ct_mid,p_mid[k],NULL,&alpha_mid,
									  &beta_mid);

		tu[k] = gtc.rad2deg*atan2((alpha_mid*dct + beta_mid*dsa),
										  (alpha_mid*dct - beta_mid*dsa));

		if (dsa == 0.0)
			rsubrho[k] = cppGSW_INVALID_VALUE;
		else
			rsubrho[k] = (alpha_mid*dct)/(beta_mid*dsa);
	}
}

/***************************************************************************
% gsw_ntp_pt_vs_CT_ratio                    ratio of gradients of potential
%                             temperature and Conservative Temperature in a
%                            neutral tangent plane (in a locally-referenced
%                              potential density surface)(75-term equation)
% =========================================================================
%
% USAGE: gsw_ntp_pt_vs_CT_ratio(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the ratio of the two-dimensional gradient of potential
%  temperature versus that of Conservative Temperature, CT, along the
%  neutral tangent plane.  The potential temperature is the regular one
%  which has a reference sea pressure of 0 dbar.  Part of the calculation
%  uses the computationally-efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  ntp_pt_vs_CT_ratio  =  The ratio of the spatial gradient of
%                         potential temperature versus that of
%                         Conservative Temperature in the
%                         neutral tangent plane (ntp).         [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (A.14.5) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_ntp_pt_vs_CT_ratio(double SA, double CT, double p)
{
   double dummy, alpha, beta, pt_SA, pt_CT;

   gsw_specvol_alpha_beta(SA,CT,p, &dummy, &alpha, &beta);

   gsw_pt_first_derivatives(SA, CT, &pt_SA, &pt_CT);

   return pt_CT + pt_SA * (alpha / beta);
}

/***************************************************************************
% gsw_rho_alpha_beta_CT_exact            in-situ density, thermal expansion
%                                  & saline contraction coefficient from CT
%==========================================================================
%
% USAGE: gsw_rho_alpha_beta_CT_exact(SA,CT,p,
                                       *rho_CT_exact,
                                       *alpha_CT_exact,
                                       *beta_CT_exact)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_alpha_beta(SA,CT,p), which uses the computationally-efficient
%  75-term expression for density in terms of SA, CT and p (Roquet et al.,
%  2015)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_CT_exact    =  in-situ density                            [ kg/m^3 ]
%  alpha_CT_exact  =  thermal expansion coefficient                 [ 1/K ]
%                     with respect to Conservative Temperature
%  beta_CT_exact   =  saline contraction coefficient               [ kg/g ]
%                     at constant Conservative Temperature
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See sections (2.8), (2.18) and (2.19) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_rho_alpha_beta_ct_exact(double SA, double CT, double p,
		                                    double &rho_CT_exact,
		                                    double &alpha_CT_exact,
		                                    double &beta_CT_exact)
{
	double t = gsw_t_from_ct(SA, CT, p);

	rho_CT_exact = gsw_rho_t_exact(SA, t, p);
	alpha_CT_exact = gsw_alpha_wrt_CT_t_exact(SA, t, p);
	beta_CT_exact = gsw_beta_const_ct_t_exact(SA, t, p);
}

/***************************************************************************
% gsw_thermobaric_CT_exact                          thermobaric coefficient
%==========================================================================
%
% USAGE: gsw_thermobaric_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates the thermobaric
%  coefficient with the full TEOS-10 Gibbs function expression for density.
%  This function uses finite differences to calculate the temperature and
%  pressure derivatives.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_thermobaric(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  thermobaric_CT_exact  =  thermobaric coefficient with       [ 1/(K Pa) ]
%                           respect to Conservative Temperature.
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_thermobaric_ct_exact(double SA, double CT, double p)
{
	double v_SA_P, v_CT_P, v_SA, v_CT, dummy, dummy1, dummy2, dummy3, rho;
	double sa = 0.0;

	if (SA > 0.0) sa = SA;

	rho = gsw_rho_ct_exact(sa, CT, p);

	gsw_specvol_first_derivatives_CT_exact(sa, CT, p, v_SA, v_CT, dummy);

	gsw_specvol_second_derivatives_CT_exact(sa, CT, p, dummy1,
	                                          dummy2, dummy3, v_SA_P,
	                                          v_CT_P);

	return rho * (v_CT_P - (v_CT / v_SA) * v_SA_P);
}

/***************************************************************************
% gsw_internal_energy_second_derivatives              second derivatives of
%                                       specific interal energy of seawater
%                                                        (75-term equation)
%==========================================================================
%
% USAGE: gsw_internal_energy_second_derivatives(SA,CT,p,
                                                &u_SA_SA,
                                                &u_SA_CT,
                                                &u_CT_CT,
                                                &u_SA_P,
                                                &u_CT_P)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of
%  internal energy,
%  (1) u_SA_SA, second order derivative with respect to Absolute Salinity
%      at constant CT & p.
%  (2) u_SA_CT, second order derivative with respect to SA & CT at
%      constant p.
%  (3) u_CT_CT, second order derivative with respect to CT at constant
%      SA & p.
%  (4) u_SA_P, second-order derivative with respect to SA & P at
%      constant CT.
%  (5) u_CT_P, second-order derivative with respect to CT & P at
%      constant SA.
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_internal_energy_second_derivatives_CT_exact(SA,CT,p) which uses the
%  full Gibbs function (IOC et al., 2010).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  u_SA_SA =  The second derivative of internal energy with respect to
%             Absolute Salinity at constant CT & p.     [ (J/kg)(g/kg)^-2 ]
%  u_SA_CT =  The second derivative of internal energy with respect to
%             SA & CT at constant p.                [ (J/kg)(g/kg)^-1 K^-1]
%  u_CT_CT =  The second derivative of internal energy with respect to
%             CT at constant SA and p.                      [ (J/kg) K^-2 ]
%  u_SA_P  =  The second derivative of internal energy with respect to
%             SA & P at constant CT.              [ (J/kg)(g/kg)^-1 Pa^-1 ]
%  u_CT_P  =  The second derivative of internal energy with respect to
%             CT & P at constant SA.                  [ (J/kg) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_internal_energy_second_derivatives(double SA, double CT,
                                                      double p,
		                                                double &u_SA_SA,
		                                                double &u_SA_CT,
		                                                double &u_CT_CT,
		                                                double &u_SA_P,
		                                                double &u_CT_P)
{
	double *h_SA_SA=new double, *h_SA_CT=new double, *h_CT_CT=new double;
	double *v_SA_SA=new double, *v_SA_CT=new double, *v_CT_CT=new double;
	double *v_SA_P=new double, *v_CT_P=new double;
	double P, sa = 0.0;

	if (SA > 0.0) sa = SA;

	P = (gtc.db2pa * p + gtc.gsw_p0);

	gsw_enthalpy_second_derivatives(sa,CT,p, h_SA_SA, h_SA_CT, h_CT_CT);
	gsw_specvol_second_derivatives(sa,CT,p,v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P);

	u_SA_SA = *h_SA_SA - P * *v_SA_SA;

	u_SA_CT = *h_SA_CT - P * *v_SA_CT;

	u_CT_CT = *h_CT_CT - P * *v_CT_CT;

	u_SA_P = -P * *v_SA_P;

	u_CT_P = -P * *v_CT_P;

	/** clean up */
	delete h_SA_SA;
	delete h_SA_CT;
	delete h_CT_CT;
	delete v_SA_SA;
	delete v_SA_CT;
	delete v_CT_CT;
	delete v_SA_P;
	delete v_CT_P;
}

/***************************************************************************
% gsw_sound_speed_ct_exact                                      sound speed
%==========================================================================
%
% USAGE: gsw_sound_speed_ct_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the speed of sound in seawater from Absolute Salinity and
%  Conservative Temperature and pressure.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sound_speed(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  sound_speed_CT_exact  =  speed of sound in seawater              [ m/s ]
%
% AUTHOR:
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.17.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sound_speed_ct_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA,CT,p);

	return gsw_sound_speed_t_exact(SA,t,p);
}

/***************************************************************************
% gsw_internal_energy_first_derivatives_CT_exact       first derivatives of
%                                       specific interal energy of seawater
%==========================================================================
%
% USAGE: gsw_internal_energy_first_derivatives_CT_exact(SA,CT,p,
                                                         &u_SA,
                                                         &u_CT,
                                                         &u_P)
%
% DESCRIPTION:
%  Calculates the first order derivates of specific internal energy of
%  seawater.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_internal_energy_first_derivatives(SA,CT,p), which uses the
%  computationally efficient polynomial for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  u_SA = The first derivative of internal energy with respect to
%           Absolute Salinity at constant CT & p.
%                                          [ (J/kg)(g/kg)^-1 ] i.e. [ J/g ]
%  u_CT = The first derivative of internal energy with respect to
%           Conservative Temperature at constant SA & p.    [ (J/kg) K^-1 ]
%  u_P = The first derivative of internal energy with respect to
%           pressure at constant SA & CT.                  [ (J/kg) Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_internal_energy_first_derivatives_ct_exact(double SA,
		double CT,
		double p,
		double &u_SA,
		double &u_CT,
		double &u_P)
{
	double v_SA, v_CT, v_P, v, h_SA, h_CT, P, sa = 0.0;

	if (SA > 0.0) sa = SA;

	P = gtc.db2pa * p + gtc.gsw_p0;

	gsw_enthalpy_first_derivatives_ct_exact(sa, CT, p, &h_SA, &h_CT);
	v = gsw_specvol_CT_exact(sa,CT,p);
	gsw_specvol_first_derivatives_CT_exact(sa, CT, p, v_SA, v_CT, v_P);

	u_SA = h_SA - P * v_SA;

	u_CT = h_CT - P * v_CT;

	u_P = v - P * v_P;
}

/***************************************************************************
% gsw_t_from_entropy                                    in-situ temperature
%                                                  as a function of entropy
% =========================================================================
%
% USAGE: gsw_t_from_entropy(SA,entropy,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature with entropy as an input variable.
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_t_from_entropy(double SA, double entropy, double p)
{
	double p0 = 0.0, pt = gsw_pt_from_entropy(SA,entropy);

	/**
	   Note that pt is potential temperature
	   with a reference pressure of zero.
	*/

	return gsw_pt_from_t(SA, pt, p0, p);
}

/**************************************************************************
% gsw_sigma0_pt0_exact                           potential density anomaly,
%                                 being potential density minus 1000 kg/m^3
%==========================================================================
%
% USAGE: gsw_sigma0_pt0_exact(SA,pt0)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference sea pressure of
%  zero (0) dbar.  The temperature input to this function is potential
%  temperature referenced to zero dbar.
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  pt0  =  potential temperature with respect to a
%          reference sea pressure of 0 dbar (ITS-90)              [ deg C ]
%
% OUTPUT:
%  sigma0_pt0_exact  =  potential density anomaly with           [ kg/m^3 ]
%                       respect to a reference pressure of 0 dbar,
%                       that is, potential density minus 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.6.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma0_pt0_exact(double SA, double pt0)
{
	double sa = 0.0;

	if (SA > 0.0) sa = SA;

	return gsw_rho_t_exact(sa, pt0, 0.0) - 1000.0;
}

/***************************************************************************
% gsw_CT_from_rho_exact               Conservative Temperature from density
% =========================================================================
%
% USAGE: gsw_CT_from_rho_exact(rho,SA,p, *CT, *CT_multiple)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_CT_from_rho(rho,SA,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is |density|, not |density anomaly|.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  CT  =  Conservative Temperature                                [ deg C ]
%  CT_multiple  =  Conservative Temperature                       [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both
%      valid solutions.  To see this second solution the user must call the
%      programme with two outputs (i.e. [CT,CT_multiple]), if there is only
%      one possible solution and the programme has been called with two
%      outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_ct_from_rho_exact(double rho, double SA, double p,
													 double &CT,
													 double &CT_multiple)
{
	double t, t_multiple;

	gsw_t_from_rho_exact(rho, SA, p, t, t_multiple);

	CT = gsw_ct_from_t(SA, t, p);
	CT_multiple = gsw_ct_from_t(SA, t_multiple, p);
}

/***************************************************************************
% gsw_cabbeling_CT_exact                              cabbeling coefficient
%==========================================================================
%
% USAGE: gsw_cabbeling_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the cabbeling coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates the cabbeling
%  coefficient with the full TEOS-10 Gibbs function expression for specific
%  volume.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_cabbeling(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  cabbeling  =  cabbeling coefficient with respect to            [ 1/K^2 ]
%                Conservative Temperature.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.9.2) and (P.4) of this TEOS-10 manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_cabbeling_ct_exact(double SA, double CT, double p)
{
	double v_SA, v_CT, dummy1, rho, alpha_CT;
	double v_SA_SA, v_SA_CT, v_CT_CT, dummy2, dummy3;
	double alpha_SA, beta_SA, alpha_on_beta;

	gsw_specvol_first_derivatives_CT_exact(SA,CT,p, v_SA, v_CT, dummy1);

	gsw_specvol_second_derivatives_CT_exact(SA,CT,p,v_SA_SA, v_SA_CT, v_CT_CT, dummy2, dummy3);

	rho = gsw_rho_ct_exact(SA,CT,p);

	alpha_CT = rho * (v_CT_CT - rho * (v_CT * v_CT));

	alpha_SA = rho * (v_SA_CT - rho * v_SA * v_CT);

	beta_SA = -rho *(v_SA_SA - rho * (v_SA * v_SA));

	alpha_on_beta = gsw_alpha_on_beta_CT_exact(SA,CT,p);

	return alpha_CT + alpha_on_beta * (2.0 * alpha_SA - alpha_on_beta * beta_SA);
}

/***************************************************************************
% gsw_sigma3_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 3000 dbar
%==========================================================================
%
% USAGE: gsw_sigma3_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 3000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma3(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma3_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 3000 dbar,
%                      that is, this potential density - 1000 kg m^-3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma3_ct_exact(double SA, double CT)
{
	double t = gsw_t_from_ct(SA, CT, 3000.0);

	return gsw_rho_t_exact(SA, t, 3000.0) - 1000.0;
}

/***************************************************************************
% gsw_rho_CT_exact                                          in-situ density
%==========================================================================
%
% USAGE: gsw_rho_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density from Absolute Salinity and Conservative
%  Temperature.
%
%  Note that potential density with respect to reference pressure, p_ref,
%  is obtained by calling this function with the pressure argument being
%  p_ref (i.e. "gsw_rho_CT_exact(SA,CT,p_ref)").
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_rho(SA,CT,p),
%  which uses the computationally efficient 75-term expression for density
%  in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_CT_exact  =  in-situ density                              [ kg/m^3 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.8.2) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_rho_ct_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA, CT, p);

	return gsw_rho_t_exact(SA, t, p);
}

/***************************************************************************
% gsw_isopycnal_slope_ratio               ratio of the slopes of isopycnals
%                                      on the SA-CT diagram for p and p_ref
%                                                        (75-term equation)
% =========================================================================
%
% USAGE: gsw_isopycnal_slope_ratio(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the ratio of alpha/beta at pressure, p, to that at reference
%  pressure, p_ref.  This function uses the computationally-efficient
%  75-term expression for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  pr  =  reference pressure                                       [ dbar ]
%         ( i.e. absolute reference pressure - 10.1325 dbar )
%
% OUTPUT:
%  isopycnal_slope_ratio
%               =  The ratio of alpha/beta evaluated at        [ unitless ]
%                  pressure, p, to that at reference pressure, p_ref.
%
% AUTHOR:
%  Trevor McDougall, Paul Barker & David Jackett       [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (3.17.2) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_isopycnal_slope_ratio(double SA, double CT, double p, double p_ref)
{
   double *dummy1=new double, *alpha=new double, *beta=new double;
   double *dummy2=new double, *alpha_pref=new double, *beta_pref=new double;

   gsw_specvol_alpha_beta(SA, CT, p, dummy1, alpha, beta);
   gsw_specvol_alpha_beta(SA, CT, p_ref, dummy2, alpha_pref, beta_pref);

/**
%--------------------------------------------------------------------------
% This function calculates isopycnal_slope_ratio using the computationally
% efficient 75-term expression for specific volume as a function of SA, CT
% and p.  If one wanted to compute this with the full TEOS-10 Gibbs
% function expression for specific volume, the following lines of code will
% enable this.
%
%     t = gsw_pt_from_CT(SA,CT,p);
%     alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%     beta = gsw_beta_const_CT_t_exact(SA,t,p);
%     tr = gsw_pt_from_t(SA,pt,p_ref0,p_ref);
%     alpha_pref = gsw_alpha_wrt_CT_t_exact(SA,tr,p_ref);
%     beta_pref = gsw_beta_const_CT_t_exact(SA,tr,p_ref);
%
%--------------This is the end of the alternative code---------------------
*/
   double v1 = *alpha, v2 = *beta_pref, v3 = *alpha_pref, v4 = *beta;
   double slope_ratio = v1 * v2 / v3 * v4;

   /** clean up */
   delete dummy1;
   delete alpha;
   delete beta;
   delete dummy2;
   delete alpha_pref;
   delete beta_pref;

   return slope_ratio;
}

/**************************************************************************
% gsw_spiciness0                                    spiciness at p = 0 dbar
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_spiciness0(SA,CT)
%
% DESCRIPTION:
%  Calculates spiciness from Absolute Salinity and Conservative
%  Temperature at a pressure of 0 dbar, as described by McDougall and
%  Krzysik (2015).  This routine is based on the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  spiciness0  =  spiciness referenced to a pressure of 0 dbar,
%                 i.e. the surface                               [ kg/m^3 ]
%
% AUTHOR:
%  Oliver Krzysik and Trevor McDougall                 [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_spiciness0(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double xs, ys, spiciness0;
	double  s01 = -9.22982898371678e1,      s02 = -1.35727873628866e1,
			  s03 =  1.87353650994010e1,      s04 = -1.61360047373455e1,
			  s05 =  3.76112762286425e1,      s06 = -4.27086671461257e1,
			  s07 =  2.00820111041594e1,      s08 =  2.87969717584045e2,
			  s09 =  1.13747111959674e1,      s10 =  6.07377192990680e1,
			  s11 = -7.37514033570187e1,      s12 = -7.51171878953574e1,
			  s13 =  1.63310989721504e2,      s14 = -8.83222751638095e1,
			  s15 = -6.41725302237048e2,      s16 =  2.79732530789261e1,
			  s17 = -2.49466901993728e2,      s18 =  3.26691295035416e2,
			  s19 =  2.66389243708181e1,      s20 = -2.93170905757579e2,
			  s21 =  1.76053907144524e2,      s22 =  8.27634318120224e2,
			  s23 = -7.02156220126926e1,      s24 =  3.82973336590803e2,
			  s25 = -5.06206828083959e2,      s26 =  6.69626565169529e1,
			  s27 =  3.02851235050766e2,      s28 = -1.96345285604621e2,
			  s29 = -5.74040806713526e2,      s30 =  7.03285905478333e1,
			  s31 = -2.97870298879716e2,      s32 =  3.88340373735118e2,
			  s33 = -8.29188936089122e1,      s34 = -1.87602137195354e2,
			  s35 =  1.27096944425793e2,      s36 =  2.11671167892147e2,
			  s37 = -3.15140919876285e1,      s38 =  1.16458864953602e2,
			  s39 = -1.50029730802344e2,      s40 =  3.76293848660589e1,
			  s41 =  6.47247424373200e1,      s42 = -4.47159994408867e1,
			  s43 = -3.23533339449055e1,      s44 =  5.30648562097667,
			  s45 = -1.82051249177948e1,      s46 =  2.33184351090495e1,
			  s47 = -6.22909903460368,        s48 = -9.55975464301446,
			  s49 =  6.61877073960113;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	spiciness0 = s01
                  +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					   +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
						+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
						+xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
						+xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
						+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
						+xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));

	return (spiciness0);
}

/***************************************************************************
% gsw_spiciness1                                 spiciness at p = 1000 dbar
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_spiciness1(SA,CT,p)
%
% DESCRIPTION:
%  Calculates spiciness from Absolute Salinity and Conservative
%  Temperature at a pressure of 1000 dbar, as described by McDougall and
%  Krzysik (2015).  This routine is based on the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  spiciness1  =  spiciness referenced to a pressure of 1000 dbar
%                                                                [ kg/m^3 ]
%
% AUTHOR:
%  Oliver Krzysik and Trevor McDougall                 [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_spiciness1(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
   double xs, ys, spiciness1;
	double  s01 = -9.19874584868912e1,      s02 = -1.33517268529408e1,
			  s03 =  2.18352211648107e1,      s04 = -2.01491744114173e1,
			  s05 =  3.70004204355132e1,      s06 = -3.78831543226261e1,
			  s07 =  1.76337834294554e1,      s08 =  2.87838842773396e2,
			  s09 =  2.14531420554522e1,      s10 =  3.14679705198796e1,
			  s11 = -4.04398864750692e1,      s12 = -7.70796428950487e1,
			  s13 =  1.36783833820955e2,      s14 = -7.36834317044850e1,
			  s15 = -6.41753415180701e2,      s16 =  1.33701981685590,
			  s17 = -1.75289327948412e2,      s18 =  2.42666160657536e2,
			  s19 =  3.17062400799114e1,      s20 = -2.28131490440865e2,
			  s21 =  1.39564245068468e2,      s22 =  8.27747934506435e2,
			  s23 = -3.50901590694775e1,      s24 =  2.87473907262029e2,
			  s25 = -4.00227341144928e2,      s26 =  6.48307189919433e1,
			  s27 =  2.16433334701578e2,      s28 = -1.48273032774305e2,
			  s29 = -5.74545648799754e2,      s30 =  4.50446431127421e1,
			  s31 = -2.30714981343772e2,      s32 =  3.15958389253065e2,
			  s33 = -8.60635313930106e1,      s34 = -1.22978455069097e2,
			  s35 =  9.18287282626261e1,      s36 =  2.12120473062203e2,
			  s37 = -2.21528216973820e1,      s38 =  9.19013417923270e1,
			  s39 = -1.24400776026014e2,      s40 =  4.08512871163839e1,
			  s41 =  3.91127352213516e1,      s42 = -3.10508021853093e1,
			  s43 = -3.24790035899152e1,      s44 =  3.91029016556786,
			  s45 = -1.45362719385412e1,      s46 =  1.96136194246355e1,
			  s47 = -7.06035474689088,        s48 = -5.36884688614009,
			  s49 =  4.43247303092448;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;

	spiciness1 = s01
					   +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					   +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
                  +xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
						+xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
						+xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
						+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
						+xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));

	return (spiciness1);
}

/***************************************************************************
% gsw_spiciness2                                 spiciness at p = 2000 dbar
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_spiciness2(SA,CT,p)
%
% DESCRIPTION:
%  Calculates spiciness from Absolute Salinity and Conservative
%  Temperature at a pressure of 2000 dbar, as described by McDougall and
%  Krzysik (2015).  This routine is based on the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  spiciness2  =  spiciness referenced to a pressure of 2000 dbar
%                                                                [ kg/m^3 ]
%
% AUTHOR:
%  Oliver Krzysik and Trevor McDougall                 [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_spiciness2(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
   double xs, ys, spiciness2;
	double  s01 = -9.17327320732265e1,      s02 = -1.31200235147912e1,
			  s03 =  2.49574345782503e1,      s04 = -2.41678075247398e1,
			  s05 =  3.61654631402053e1,      s06 = -3.22582164667710e1,
			  s07 =  1.45092623982509e1,      s08 =  2.87776645983195e2,
			  s09 =  3.13902307672447e1,      s10 =  1.69777467534459,
			  s11 = -5.69630115740438,        s12 = -7.97586359017987e1,
			  s13 =  1.07507460387751e2,      s14 = -5.58234404964787e1,
			  s15 = -6.41708068766557e2,      s16 = -2.53494801286161e1,
			  s17 = -9.86755437385364e1,      s18 =  1.52406930795842e2,
			  s19 =  4.23888258264105e1,      s20 = -1.60118811141438e2,
			  s21 =  9.67497898053989e1,      s22 =  8.27674355478637e2,
			  s23 =  5.27561234412133e-1,     s24 =  1.87440206992396e2,
			  s25 = -2.83295392345171e2,      s26 =  5.14485994597635e1,
			  s27 =  1.29975755062696e2,      s28 = -9.36526588377456e1,
			  s29 = -5.74911728972948e2,      s30 =  1.91175851862772e1,
			  s31 = -1.59347231968841e2,      s32 =  2.33884725744938e2,
			  s33 = -7.87744010546157e1,      s34 = -6.04757235443685e1,
			  s35 =  5.27869695599657e1,      s36 =  2.12517758478878e2,
			  s37 = -1.24351794740528e1,      s38 =  6.53904308937490e1,
			  s39 = -9.44804080763788e1,      s40 =  3.93874257887364e1,
			  s41 =  1.49425448888996e1,      s42 = -1.62350721656367e1,
			  s43 = -3.25936844276669e1,      s44 =  2.44035700301595,
			  s45 = -1.05079633683795e1,      s46 =  1.51515796259082e1,
			  s47 = -7.06609886460683,        s48 = -1.48043337052968,
			  s49 =  2.10066653978515;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;

	spiciness2 = s01
					   +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					   +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
						+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
						+xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
						+xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
						+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
						+xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));

	return (spiciness2);
}

/***************************************************************************
% gsw_CT_second_derivatives                           second derivatives of
%                                                  Conservative Temperature
%==========================================================================
%
% USAGE:
%  gsw_CT_second_derivatives(SA,pt)
%
% DESCRIPTION:
%  Calculates the following three, second-order derivatives of Conservative
%  Temperature
%   (1) CT_SA_SA, the second derivative with respect to Absolute Salinity
%       at constant potential temperature (with p_ref = 0 dbar),
%   (2) CT_SA_pt, the derivative with respect to potential temperature
%       (the regular potential temperature which is referenced to 0 dbar)
%       and Absolute Salinity, and
%   (3) CT_pt_pt, the second derivative with respect to potential
%       temperature (the regular potential temperature which is referenced
%       to 0 dbar) at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%         (whose reference pressure is 0 dbar)
%
% OUTPUT:
%  CT_SA_SA  =  The second derivative of Conservative Temperature with
%               respect to Absolute Salinity at constant potential
%               temperature (the regular potential temperature which
%               has reference sea pressure of 0 dbar).
%               CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  CT_SA_pt  =  The derivative of Conservative Temperature with
%               respect to potential temperature (the regular one with
%               p_ref = 0 dbar) and Absolute Salinity.
%               CT_SA_pt has units of:                        [ 1/(g/kg) ]
%  CT_pt_pt  =  The second derivative of Conservative Temperature with
%               respect to potential temperature (the regular one with
%               p_ref = 0 dbar) at constant SA.
%               CT_pt_pt has units of:                              [ 1/K ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (1st September, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix A.12 of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_ct_second_derivatives(double sa,double pt,double *ct_sa_sa,
                                          double *ct_sa_pt,double *ct_pt_pt)
{
	double  ct_pt_l, ct_pt_u, ct_sa_l, ct_sa_u, pt_l, pt_u, sa_l, sa_u,
			  dsa = 1e-3, dpt = 1e-2;

	if ((ct_sa_sa != NULL))
	{
		if ((sa_l = sa - dsa) < 0.0)
      {
         sa_l = 0.0;
      }
		sa_u = sa + dsa;
		gsw_ct_first_derivatives(sa_l,pt,&ct_sa_l,NULL);
		gsw_ct_first_derivatives(sa_u,pt,&ct_sa_u,NULL);
		*ct_sa_sa = (ct_sa_u - ct_sa_l)/(sa_u - sa_l);
	}

	if ((ct_sa_pt != NULL) || (ct_pt_pt != NULL))
	{
		pt_l = pt - dpt;
		pt_u = pt + dpt;

		if ((ct_sa_pt != NULL) && (ct_pt_pt != NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,&ct_pt_l);
			gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,&ct_pt_u);
			*ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
			*ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);
		}
		else if ((ct_sa_pt != NULL) && (ct_pt_pt == NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,NULL);
			gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,NULL);
			*ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
		}
		else if ((ct_sa_pt == NULL) && (ct_pt_pt != NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,NULL,&ct_pt_l);
			gsw_ct_first_derivatives(sa,pt_u,NULL,&ct_pt_u);
			*ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);
		}
	}
}

/***************************************************************************
% gsw_geo_strf_dyn_height                            dynamic height anomaly
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_geo_strf_dyn_height(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume
%  anomaly from the pressure p of the bottle to the reference pressure
%  p_ref.
%
%  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
%  to a given reference pressure.  This is the geostrophic streamfunction
%  for the difference between the horizontal velocity at the pressure
%  concerned, p, and the horizontal velocity at p_ref.  Dynamic height
%  anomaly is the geostrophic streamfunction in an isobaric surface.  The
%  reference values used for the specific volume anomaly are
%  SSO = 35.16504 g/kg and CT = 0 deg C.  This function calculates
%  specific volume anomaly using the computationally efficient
%  expression for specific volume of Roquet et al. (2015).
%
%  This function evaluates the pressure integral of specific volume using
%  SA and CT interpolated using the MRST-PCHIP method of Barker and
%  McDougall (2020).  This "curve fitting" method uses a Piecewise Cubic
%  Hermite Interpolating Polynomial to produce a smooth curve with minimal
%  artificial watermasses between the observed data points.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  For dynamical oceanography we may
%  take the 75-term rational function expression for specific volume as
%  essentially reflecting the full accuracy of TEOS-10.  The GSW library
%  function "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to
%  test if some of one's data lies outside this "funnel".
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref =  reference pressure                                     [ dbar ]
%           ( i.e. reference absolute pressure - 10.1325 dbar )
%
%
% OUTPUT:
%  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
%   Note. If p_ref exceeds the pressure of the deepest bottle on a
%     vertical profile, the dynamic height anomaly for each bottle
%     on the whole vertical profile is returned as NaN.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2020: Two interpolation methods using
%   multiply-rotated piecewise cubic hermite interpolating polynomials.
%   J. Atmosph. Ocean. Tech., 37, pp. 605-619.
%   doi: 10.1175/JTECH-D-19-0211.1.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) and section 3.27 of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double *TeosSea::gsw_geo_strf_dyn_height(double *sa,double *ct,double *p,
                                          double p_ref,int n_levels,
                                          double *dyn_height)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int m_levels = (n_levels <= 0) ? 1 : n_levels,
		 p_cnt, top_pad, i, nz, ibottle, ipref,
		 np_max, np, ibpr=0, *iidata=NULL;

	double dp_min, dp_max, p_min, p_max, max_dp_i,
			 *b=NULL, *b_av=NULL, *dp=NULL, *dp_i=NULL, *sa_i=NULL, *ct_i=NULL,
			  *p_i=NULL, *geo_strf_dyn_height0=NULL;
	/**
	----------------------------------------------------------------------------
	   This max_dp_i is the limit we choose for the evaluation of specific
	   volume in the pressure integration.  That is, the vertical integration
	   of specific volume with respect to pressure is perfomed with the pressure
	   increment being no more than max_dp_i (the default value being 1 dbar).
	----------------------------------------------------------------------------
	*/
	max_dp_i = 1.0;
	if ((nz = m_levels) <= 1)
   {
      return (NULL);
   }

	dp = new double[nz];
	dp_min = 11000.0;
	dp_max = -11000.0;

	for (i=0; i<nz-1; i++)
	{
		if ((dp[i] = p[i+1] - p[i]) < dp_min)
      {
         dp_min = dp[i];
      }

		if (dp[i] > dp_max)
      {
         dp_max = dp[i];
      }
	}

	if (dp_min <= 0.0)
	{
		/** pressure must be monotonic */
		if (dp) delete []dp;

		return (NULL);
	}

	p_min = p[0];
	p_max = p[nz-1];

	if (p_ref > p_max)
	{
		/**the reference pressure p_ref is deeper than all bottles*/
		if (dp) delete []dp;

		return (NULL);
	}

	/** Determine if there is a "bottle" at exactly p_ref */
	ipref = -1;
	for (ibottle = 0; ibottle < nz; ibottle++)
	{
		if (p[ibottle] == p_ref)
		{
			ipref = ibottle;
			break;
		}
	}

	if ((dp_max <= max_dp_i) && (p[0] == 0.0) && (ipref >= 0))
	{
		/**
		  vertical resolution is good (bottle gap is no larger than max_dp_i)
		  & the vertical profile begins at the surface (i.e. at p = 0 dbar)
		  & the profile contains a "bottle" at exactly p_ref.
		*/
		b = new double[3*nz];
		b_av = b+nz;
		geo_strf_dyn_height0 = b_av+nz;

		for (i=0; i<nz; i++)
		{
			b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
			if (i > 0)
         {
            b_av[i-1] = 0.5*(b[i] + b[i-1]);
         }
		}
		/**
		  "geo_strf_dyn_height0" is the dynamic height anomaly with respect
		  to p_ref = 0 (the surface).
		*/
		geo_strf_dyn_height0[0] = 0.0;

		for (i=1; i<nz; i++)
      {
         geo_strf_dyn_height0[i] = b_av[i]*dp[i]*gtc.db2pa;
      }

		for (i=1; i<nz; i++) /** cumulative sum */
		{
			geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
											  - geo_strf_dyn_height0[i];
      }

		for (i=0; i<nz; i++)
      {
         dyn_height[i] = geo_strf_dyn_height0[i]
								 - geo_strf_dyn_height0[ipref];
      }

		if (b) delete []b;
	}
	else
	{
		/**
		  Test if there are vertical gaps between adjacent "bottles" which are
		  greater than max_dp_i, and that there is a "bottle" exactly at the
		  reference pressure.
		*/
		iidata = new int[nz+1];
		if ((dp_max <= max_dp_i) && (ipref >= 0))
		{
			/**
			  Vertical resolution is already good (no larger than max_dp_i), and
			  there is a "bottle" at exactly p_ref.
			*/
			sa_i = new double[2*(nz+1)];
			ct_i = sa_i+nz+1;
			p_i = new double[nz+1];

			if (p_min > 0.0)
			{
				/**
				  resolution is fine and there is a bottle at p_ref, but
				  there is not a bottle at p = 0. So add an extra bottle.
				*/
				for (i=0; i<nz; i++)
				{
					sa_i[i+1]       = sa[i];
					ct_i[i+1]       = ct[i];
					p_i[i+1]        = p[i];
				}

				sa_i[0] = sa[0];
				ct_i[0] = ct[0];
				p_i[0] = 0.0;
				ibpr = ipref+1;
				p_cnt = nz+1;

				for (i=0; i<p_cnt; i++)
            {
               iidata[i] = i;
            }
			}
			else
			{
				/**
				  resolution is fine, there is a bottle at p_ref, and
				  there is a bottle at p = 0
				*/
				memmove(sa_i, sa, nz*sizeof (double));
				memmove(ct_i, ct, nz*sizeof (double));
				memmove(p_i, p, nz*sizeof (double));
				ibpr = ipref;

				for (i=0; i<nz; i++)
				{
				   iidata[i] = i;
            }

				p_cnt = nz;
			}
		}
		else
		{
			/**
			  interpolation is needed.
			*/
			np_max = 2*rint(p[nz-1]/max_dp_i+0.5);
			if (p_i) delete []p_i;
			p_i = new double[np_max];
			/** sa_i is allocated below, when its size is known */
			if (p_min > 0.0)
			{
				/**
				  there is not a bottle at p = 0.
				*/
				if (p_ref < p_min)
				{
					/**
					  p_ref is shallower than the minimum bottle pressure.
					*/
					p_i[0] = 0.0;
					gsw_p_sequence(p_i[0],p_ref,max_dp_i, p_i+1,&np);
					ibpr = p_cnt = np;
					p_cnt++;
					gsw_p_sequence(p_ref,p_min,max_dp_i, p_i+p_cnt,&np);
					p_cnt += np;
					top_pad = p_cnt;
				}
				else
				{
					/**
					  p_ref is deeper than the minimum bottle pressure.
					*/
					p_i[0] = 0.0;
					p_i[1] = p_min;
					top_pad = 2;
					p_cnt = 2;
				}
			}
			else
			{
				/**
				  there is a bottle at p = 0.
				*/
				p_i[0] = p_min;
				top_pad = 1;
				p_cnt = 1;
			}

			for (ibottle=0; ibottle < nz-1; ibottle++)
			{
				iidata[ibottle] = p_cnt-1;
				if (p[ibottle] == p_ref) ibpr = p_cnt-1;
				if (p[ibottle] < p_ref && p[ibottle+1] > p_ref)
				{
					/**
					  ... reference pressure is spanned by bottle pairs -
					  need to include p_ref as an interpolated pressure.
					*/
					gsw_p_sequence(p[ibottle],p_ref,max_dp_i, p_i+p_cnt,&np);
					p_cnt += np;
					ibpr = p_cnt-1;
					gsw_p_sequence(p_ref,p[ibottle+1],max_dp_i,p_i+p_cnt,&np);
					p_cnt += np;
				}
				else
				{
					/**
					  ... reference pressure is not spanned by bottle pairs.
					*/
					gsw_p_sequence(p[ibottle],p[ibottle+1],max_dp_i,
										p_i+p_cnt,&np);
					p_cnt += np;
				}
			}

			iidata[nz-1] = p_cnt-1;

			if (p[nz-1] == p_ref)
			{
			   ibpr = p_cnt-1;
         }

			sa_i = new double[2*p_cnt];
			ct_i = sa_i+p_cnt;

			if (top_pad > 1)
			{
				gsw_linear_interp_sa_ct(sa,ct,p,nz,
												p_i,top_pad-1,sa_i,ct_i);
			}

			gsw_rr68_interp_sa_ct(sa,ct,p,nz,p_i+top_pad-1,p_cnt-top_pad+1,
										 sa_i+top_pad-1,ct_i+top_pad-1);
		}

		b = new double[4*p_cnt];
		b_av = b+p_cnt;
		dp_i = b_av+p_cnt;
		geo_strf_dyn_height0 = dp_i+p_cnt;

		for (i=0; i<p_cnt; i++)
		{
			b[i] = gsw_specvol_anom_standard(sa_i[i],ct_i[i],p_i[i]);
			if (i > 0)
			{
				dp_i[i-1] = p_i[i]-p_i[i-1];
				b_av[i-1] = 0.5*(b[i] + b[i-1]);
			}
		}
		/**
		  "geo_strf_dyn_height0" is the dynamic height anomaly with respect
		  to p_ref = 0 (the surface).
		*/
		geo_strf_dyn_height0[0] = 0.0;
		for (i=1; i<p_cnt; i++)
      {
         geo_strf_dyn_height0[i] = b_av[i-1]*dp_i[i-1];
      }

		for (i=1; i<p_cnt; i++) /** cumulative sum */
      {
         geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
											  - geo_strf_dyn_height0[i];
      }

		for (i=0; i<nz; i++)
		{
		   dyn_height[i] = (geo_strf_dyn_height0[iidata[i]]
								  - geo_strf_dyn_height0[ibpr])*gtc.db2pa;
      }

		if (b) delete []b;
		if (iidata) delete []iidata;
		if (sa_i) delete []sa_i;
		if (p_i)	delete []p_i;
	}

	if (dp) delete []dp;

	return (dyn_height);
}

/****************************************************************
=======
method: gsw_p_sequence(p1, p2, max_dp_i, *pseq, *nps)
=======
******************************************************************/
void TeosSea::gsw_p_sequence(double p1,double p2,double max_dp_i,
                              double *pseq,int *nps)
{
	double dp, pstep;
	int n, i;

	dp = p2 - p1;
	n = ceil(dp/max_dp_i);
	pstep = dp/n;

	if (nps != NULL)
	{
	   *nps = n;

   }

	/**
	  Generate the sequence ensuring that the value of p2 is exact to
	  avoid round-off issues, ie. dont do "pseq = p1+pstep*(i+1)".
	*/
	for (i=0; i<n; i++)
   {
      pseq[i] = p2-pstep*(n-1-i);
   }
}

/**************************************************************************
% gsw_pt_second_derivatives     second derivatives of potential temperature
% =========================================================================
%
% USAGE:
%  [pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of potential
%  temperature (the regular potential temperature which has a reference
%  sea pressure of 0 dbar),
%   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity
%       at constant Conservative Temperature,
%   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
%       and Absolute Salinity, and
%   (3) pt_CT_CT, the second derivative with respect to Conservative
%       Temperature at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt_SA_SA  =  The second derivative of potential temperature (the
%               regular potential temperature which has reference sea
%               pressure of 0 dbar) with respect to Absolute Salinity
%               at constant Conservative Temperature.
%               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  pt_SA_CT  =  The derivative of potential temperature with respect
%               to Absolute Salinity and Conservative Temperature.
%               pt_SA_CT has units of:                         [ 1/(g/kg) ]
%  pt_CT_CT  =  The second derivative of potential temperature (the
%               regular one with p_ref = 0 dbar) with respect to
%               Conservative Temperature at constant SA.
%               pt_CT_CT has units of:                              [ 1/K ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.9) and (A.12.10) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_pt_second_derivatives(double sa,double ct,double *pt_sa_sa,
                                          double *pt_sa_ct,double *pt_ct_ct)
{
	double ct_l, ct_u, pt_ct_l, pt_ct_u, pt_sa_l, pt_sa_u, sa_l, sa_u,
			 dct = 1e-2, dsa = 1e-3;

	if (pt_sa_sa != NULL)
	{
		if ((sa_l = sa - dsa) < 0.0)
			sa_l = 0.0;
		sa_u = sa + dsa;
		gsw_pt_first_derivatives(sa_l,ct,&pt_sa_l,NULL);
		gsw_pt_first_derivatives(sa_u,ct,&pt_sa_u,NULL);
		*pt_sa_sa = (pt_sa_u - pt_sa_l)/(sa_u - sa_l);
	}

	if (pt_sa_ct != NULL || pt_ct_ct != NULL)
	{
		ct_l = ct - dct;
		ct_u = ct + dct;

		if ((pt_sa_ct != NULL) && (pt_ct_ct != NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,&pt_ct_l);
			gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,&pt_ct_u);
			*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
			*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
		}
		else if ((pt_sa_ct != NULL) && (pt_ct_ct == NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,NULL);
			gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,NULL);
			*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
		}
		else if ((pt_sa_ct == NULL) && (pt_ct_ct != NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,NULL,&pt_ct_l);
			gsw_pt_first_derivatives(sa,ct_u,NULL,&pt_ct_u);
			*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
		}
	}
}

/***************************************************************************
=======
method:
gsw_refine_grid_for_dh(*p,p_ref,nz,dp,*p_i,ni_max,*p_indices,*p_ref_ind_ptr)
======
****************************************************************************/
int TeosSea::gsw_refine_grid_for_dh(double *p,double p_ref,int nz,double dp,
                                    double *p_i,int ni_max,int *p_indices,
                                    int *p_ref_ind_ptr)
{
	int i, iuniform, iorig;
	double p_next;
	/** Dont add a new point if it is within p_tol of an original. */
	double p_tol = 0.001 * dp;

	p_i[0] = p[0];
	p_indices[0] = 0;
	*p_ref_ind_ptr = -1;  /** initialize to a flag value */

	if (p_ref <= p[0] + p_tol)
	{
		*p_ref_ind_ptr = 0;
	}

	for (i=1, iuniform=1, iorig=1; i<ni_max && iorig<nz; i++)
	{
		/** Candidate insertion based on uniform grid: */
		p_next = p[0] + dp * iuniform;
		/** See if we need to insert p_ref: */
		if (*p_ref_ind_ptr == -1 && p_ref <= p_next && p_ref <= p[iorig])
		{
			p_i[i] = p_ref;
			*p_ref_ind_ptr = i;

			if (p_ref == p[iorig])
			{
				p_indices[iorig] = i;
				iorig++;
			}

			if (p_ref > p_next - p_tol)
			{
				iuniform++;
			}
			continue;
		}

		/** We did not insert p_ref, so insert either p_next or p[iorig]. */
		if (p_next < p[iorig] - p_tol)
		{
			p_i[i] = p_next;
			iuniform++;
		}
		else
		{
			p_i[i] = p[iorig];
			p_indices[iorig] = i;
			/** Skip this p_next if it is close to the point we just added. */
			if (p_next < p[iorig] + p_tol)
			{
				iuniform++;
			}
			iorig++;
		}
	}
	if (i == ni_max)
	{
		return (-1);  /** error  */
	}

	return (i);  /** number of elements in p_i */
}

/***************************************************************************
% gsw_rho_first_derivatives                SA, CT and p partial derivatives
%                                             of density (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_rho_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the three (3) partial derivatives of in-situ density with
%  respect to Absolute Salinity, Conservative Temperature and pressure.
%  Note that the pressure derivative is done with respect to pressure in
%  Pa, not dbar.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_SA  =  partial derivative of density           [ (kg/m^3)(g/kg)^-1 ]
%                 with respect to Absolute Salinity
%  rho_CT  =  partial derivative of density                  [ kg/(m^3 K) ]
%                 with respect to Conservative Temperature
%  rho_P   =  partial derivative of density                 [ kg/(m^3 Pa) ]
%                 with respect to pressure in Pa
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_rho_first_derivatives(double sa,double ct,double p,
                                          double *drho_dsa,double *drho_dct,
                                          double *drho_dp)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rho2, v_ct, v_p, v_sa, xs, ys, z, v;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	rho2 = pow(1.0/v, 2.0);

	if (drho_dsa != NULL)
	{
		v_sa = gsw_gsvco_b(xs, ys, z);
		*drho_dsa = -rho2*0.5*gtc.gsw_sfac*v_sa/xs;
	}

	if (drho_dct != NULL)
	{
		v_ct = gsw_gsvco_a(xs, ys, z);
		*drho_dct = -rho2*0.025*v_ct;
	}

	if (drho_dp != NULL)
	{
		v_p = gsw_gsvco_c(xs, ys, z);
		*drho_dp = 1e-4*gtc.pa2db*-rho2*v_p;
	}

	return;
}

/***************************************************************************
% gsw_sigma0_ct_exact                        potential density anomaly with
%                                         reference sea pressure of 0 dbar.
%==========================================================================
%
% USAGE: gsw_sigma0_ct_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 0 dbar,
%  this being this particular potential density minus 1000 kg/m^3.  This
%  function has inputs of Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma0(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma0_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 0 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma0_ct_exact(double SA, double CT)
{
   double pt0 = gsw_pt_from_ct(SA, CT);

   return gsw_sigma0_pt0_exact(SA, pt0);
}

/***************************************************************************
% gsw_rho_second_derivatives_wrt_enthalpy                second derivatives
%                        of rho with respect to enthalpy (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_rho_second_derivatives_wrt_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of rho with
%  respect to enthalpy,
%   (1) rho_SA_SA, second-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) rho_SA_h, second-order derivative with respect to SA & h at
%       constant p.
%   (3) rho_h_h, second-order derivative with respect to h at
%       constant SA & p.
%
%  Note that this function uses the using the computationally-efficient
%  expression for specific volume (Roquet et al., 2015).  There is an
%  alternative to calling this function, namely
%  gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses
%  the full Gibbs function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_SA_SA = The second-order derivative of rho with respect to
%              Absolute Salinity at constant h & p.     [ J/(kg (g/kg)^2) ]
%  rho_SA_h  = The second-order derivative of rho with respect to
%              SA and h at constant p.                   [ J/(kg K(g/kg)) ]
%  rho_h_h   = The second-order derivative of rho with respect to h at
%              constant SA & p
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_rho_second_derivatives_wrt_enthalpy(double sa,double ct,
                                                      double p,
                                                      double *rho_sa_sa,
                                                      double *rho_sa_h,
                                                      double *rho_h_h)
{
	double  rec_v, rec_v2, rec_v3, v_h, v_h_h, v_sa, v_sa_h, v_sa_sa,
			  *pv_sa=NULL, *pv_h=NULL, *pv_sa_sa=NULL, *pv_sa_h=NULL, *pv_h_h=NULL;

	pv_sa   = ((rho_sa_sa != NULL) || (rho_sa_h != NULL)) ? &v_sa : NULL;
	pv_h    = ((rho_sa_h != NULL) || (rho_h_h != NULL)) ?  &v_h : NULL;

	gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,pv_sa,pv_h);

	pv_sa_sa = ((rho_sa_sa != NULL)) ? &v_sa_sa : NULL;
	pv_sa_h  = ((rho_sa_h != NULL)) ? &v_sa_h : NULL;
	pv_h_h   = ((rho_h_h != NULL)) ? &v_h_h : NULL;

	gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,pv_sa_sa,pv_sa_h,pv_h_h);

	rec_v = 1.0/gsw_specvol(sa,ct,p);
	rec_v2 = rec_v*rec_v;
	rec_v3 = rec_v2*rec_v;

	if (rho_sa_sa != NULL)
		*rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;

	if (rho_sa_h != NULL)
		*rho_sa_h = -v_sa_h*rec_v2 + 2.0*v_sa*v_h*rec_v3;

	if (rho_h_h != NULL)
		*rho_h_h = -v_h_h*rec_v2 + 2.0*v_h*v_h*rec_v3;

}

/**************************************************************************
==========================================================================
method: gsw_util_interp1q_int (x, iy, x_i)
==========================================================================
  Returns the value of the 1-D method iy (integer) at the points of column
  vector x_i using linear interpolation. The vector x specifies the
  coordinates of the underlying interval.
  result(y_i)
==========================================================================
*************************************************************************/
double *TeosSea::gsw_util_interp1q_int(int nx,double *x,int *iy,int nxi,
                                          double *x_i,double *y_i)
{
	char *in_rng=NULL;
	int *j=NULL, *k=NULL, *r=NULL, *jrev=NULL, *ki=NULL, imax_x, imin_x, i, n, m, ii;
	double *xi=NULL, *xxi=NULL, u, max_x, min_x;

	if (nx <= 0 || nxi <= 0)
	{
	   return (NULL);
   }

	min_x = max_x = x[0];
	imin_x = imax_x = 0;

	for (i=0; i<nx; i++)
	{
		if (x[i] < min_x)
		{
			min_x = x[i];
			imin_x = i;
		}
		else if (x[i] > max_x)
		{
			max_x = x[i];
			imax_x = i;
		}
	}

	in_rng = new char[nxi];
	memset(in_rng, 0, nxi*sizeof (char));

	for (i=n=0; i<nxi; i++)
	{
		if (x_i[i] <= min_x)
		{
			y_i[i] = iy[imin_x];
		}
		else if (x_i[i] >= max_x)
		{
			y_i[i] = iy[imax_x];
		}
		else
		{
			in_rng[i] = 1;
			n++;
		}
	}

	if (n==0)
	{
		if (in_rng) delete []in_rng;
		return (y_i);
	}

	xi = new double[n];
	k = new int[3*n];
	ki = k+n;
	r = ki+n;
	m  = nx + n;
	xxi = new double[m];
	j = new int[2*m];
	jrev = j+m;
	ii = 0;

	for (i = 0; i<nxi; i++)
	{
		if (in_rng[i])
		{
			xi[ii] = x_i[i];
			ki[ii] = i;
			ii++;
		}
	}

	if (in_rng) delete []in_rng;
	/**
	 Note that the following operations on the index
	 vectors jrev and r depend on the sort utility
	 gsw_util_sort_dbl() consistently ordering the
	 sorting indexes either in ascending or descending
	 sequence for replicate values in the real vector.
	*/
	gsw_util_sort_dbl(xi, n, k);
	for (i = 0; i<nx; i++)
	{
	   xxi[i] = x[i];
   }

	for (i = 0; i<n; i++)
	{
		xxi[nx+i] = xi[k[i]];
   }

	gsw_util_sort_dbl(xxi, nx+n, j);

	for (i = 0; i<nx+n; i++)
	{
		jrev[j[i]] = i;
   }

	for (i = 0; i<n; i++)
	{
		r[k[i]] = jrev[nx+i] - i-1;
   }

	for (i = 0; i<n; i++)
	{
		u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
		y_i[ki[i]] = iy[r[i]] + (iy[r[i]+1]-iy[r[i]])*u;
	}

	if (j) delete []j;
	if (xxi) delete []xxi;
	if (k) delete []k;
	if (xi) delete []xi;

	return (y_i);
}

/**************************************************************************
==========================================================================
method: gsw_util_linear_interp (x, y, x_i)
==========================================================================
  Returns the values of the methods y{ny} at the points of column
  vector x_i using linear interpolation. The vector x specifies the
  coordinates of the underlying interval, and the matrix y specifies
  the method values at each x coordinate. Note that y has dimensions
  nx x ny and y_i has dimensions nxi x ny.
  This method was adapted from Matlabs interp1q.
  result(y_i)
==========================================================================
*************************************************************************/
double *TeosSea::gsw_util_linear_interp(int nx,double *x,int ny,double *y,int nxi,double *x_i,double *y_i)
{
	char *in_rng=NULL;
	int *j=NULL, *k=NULL, *r=NULL, *jrev=NULL, *ki=NULL,
		  imax_x, imin_x, i, n, m, ii, jy, jy0, jyi0, r0;
	double *xi=NULL, *xxi=NULL, u, max_x, min_x;

	if (nx <= 0 || nxi <= 0 || ny <= 0)
		return (NULL);
	min_x = max_x = x[0];
	imin_x = imax_x = 0;
	for (i=0; i<nx; i++)
	{
		if (x[i] < min_x)
		{
			min_x = x[i];
			imin_x = i;
		}
		else if (x[i] > max_x)
		{
			max_x = x[i];
			imax_x = i;
		}
	}
	in_rng = new char[nxi];
	memset(in_rng, 0, nxi*sizeof (char));
	for (i=n=0; i<nxi; i++)
	{
		if (x_i[i] <= min_x)
		{
			for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
				y_i[jyi0+i] = y[jy0+imin_x];
		}
		else if (x_i[i] >= max_x)
		{
			for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
				y_i[jyi0+i] = y[jy0+imax_x];
		}
		else
		{
			in_rng[i] = 1;
			n++;
		}
	}
	if (n==0)
	{
		if (in_rng) delete []in_rng;
		return (y_i);
	}
	xi = new double[n];
	k = new int[3*n];
	ki = k+n;
	r = ki+n;
	m  = nx + n;
	xxi = new double[m];
	j = new int[2*m];
	jrev = j+m;
	ii = 0;
	for (i = 0; i<nxi; i++)
	{
		if (in_rng[i])
		{
			xi[ii] = x_i[i];
			ki[ii] = i;
			ii++;
		}
	}
	if (in_rng) delete []in_rng;
	/**
	   This algorithm mimics the Matlab interp1q method.
	   An explaination of this algorithm:
	   We have points we are interpolating from (x) and
	   points that we are interpolating to (xi).  We
	   sort the interpolating from points, concatenate
	   them with the interpolating to points and sort the result.
	   We then construct index r, the interpolation index in x for
	   each point in xi.
	   Note that the following operations on the index
	   vectors jrev and r depend on the sort utility
	   gsw_util_sort_dbl() consistently ordering the
	   sorting indexes either in ascending or descending
	   sequence for replicate values in the real vector.
	*/
	gsw_util_sort_dbl(xi, n, k);
	memmove(xxi, x, nx*sizeof (double));
	memmove(xxi+nx, xi, n*sizeof (double));
	gsw_util_sort_dbl(xxi, m, j);
	for (i = 0; i<m; i++)
		jrev[j[i]] = i;
	for (i = 0; i<n; i++)
		r[k[i]] = jrev[nx+i] - i - 1;
	/** this is now the interpolation index in x for a point in xi */
	for (jy=jy0=jyi0=0; jy < ny; jy++, jy0+=nx, jyi0+=nxi)
	{
		for (i = 0; i<n; i++)
		{
			u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
			r0 = jy0+r[i];
			y_i[jyi0+ki[i]] = y[r0] + (y[r0+1]-y[r0])*u;
		}
	}
	if (j) delete []j;
	if (xxi) delete []xxi;
	if (k) delete []k;
	if (xi) delete []xi;

	return (y_i);
}

/*************************************************************************
=================================================
method: gsw_util_pchip_interp(*x,*y,n,*xi,*yi,ni)
=================================================
   Piecewise-Hermite algorithm from
   https://en.wikipedia.org/wiki/Cubic_Hermite_spline
   Extrapolation to points outside the range is done by setting those
   points to the corresponding end values.
   The input x must be monotonically increasing; the interpolation points,
   xi, may be in any order, but the algorithm will be faster if they are
   monotonic, increasing or decreasing.
   Returns 0 on success, 1 if it fails because there are fewer than 2 points,
   2 if it fails because x is not increasing.
   Consistent with other GSW-C code at present, the memory allocations
   are assumed to succeed.
*************************************************************************/
int TeosSea::gsw_util_pchip_interp(double *x,double *y,int n,double *xi,
                                    double *yi,int ni)
{
	double *d=NULL;
	double t, tt, ttt, xx, dx;
	int i, j0, j1, err;
	double h00, h10, h01, h11;

	if (n<2)
	{
		return 1;
	}
	d = new double[n];
	err = gsw_pchip_derivs(x, y, n, d);
	if (err)
	{
		if (d) delete []d;
		return 2;
	}
	j0 = 0;
	for (i=0; i<ni; i++)
	{
		xx = xi[i];
		/** Linear search is appropriate and probably optimal for the
		   expected primary use case of interpolation to a finer grid.
		   It is inefficient but still methodal in the worst case of
		   randomly distributed xi.
		*/
		while (xx < x[j0] && j0 > 0)
		{
			j0--;
		}
		while (xx > x[j0+1] && j0 < n - 2)
		{
			j0++;
		}
		j1 = j0 + 1;
		if (xx >= x[j0] && xx <= x[j1])
		{
			dx = x[j1] - x[j0];
			t = (xx - x[j0]) / dx;
			tt = t * t;
			ttt = tt * t;
			/** Using intermediate variables for readability. */
			h00 = (2*ttt - 3*tt + 1);
			h10 =  (ttt - 2*tt + t);
			h01 = (-2*ttt + 3*tt);
			h11 = (ttt - tt);
			yi[i] = y[j0] * h00 + d[j0] * dx * h10 +
					  y[j1] * h01 + d[j1] * dx * h11;
		}
		else
		{
			/** extrapolate with constant end values */
			yi[i] = (xx < x[0]) ? y[0] : y[n-1];
		}
	}
	if (d)
	{
	   delete []d;
   }

	return 0;
}

/***************************************************************************
% gsw_specvol_alpha_beta                 specific volume, thermal expansion
%                       & saline contraction coefficient (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_specvol_alpha_beta(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.  This function uses the
%  computationally-efficient expression for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol =  specific volume                                      [ m/kg ]
%  alpha   =  thermal expansion coefficient                         [ 1/K ]
%             with respect to Conservative Temperature
%  beta    =  saline (i.e. haline) contraction                     [ kg/g ]
%             coefficient at constant Conservative Temperature
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_specvol_alpha_beta(double sa,double ct,double p,
                                       double *specvol,double *alpha,
                                       double *beta)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double  v, v_ct, v_sa_part, xs, ys, z;

	xs = sqrt(gtc.gsw_sfac * sa + gtc.offset);
	ys = ct * 0.025;
	z = p * gtc.rec_db2pa;
	v = gsw_gsvco_v(xs, ys, z);

	if (specvol != NULL)
		*specvol = v;

	if (alpha != NULL)
	{
		v_ct = gsw_gsvco_a(xs, ys, z);
		*alpha = 0.025*v_ct/v;
	}

	if (beta != NULL)
	{
		v_sa_part = gsw_gsvco_b(xs, ys, z);
		*beta = -v_sa_part*0.5*gtc.gsw_sfac/(v*xs);
	}
}

/***************************************************************************
% gsw_SP_from_Sstar              Practical Salinity from Preformed Salinity
%==========================================================================
%
% USAGE:
%  gsw_SP_from_Sstar(Sstar,p,long,lat)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Preformed Salinity.
%
% INPUT:
%  Sstar  =  Preformed Salinity                                    [ g/kg ]
%  p      =  sea pressure                                          [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  long   =  longitude in decimal degrees                    [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat    =  latitude in decimal degrees north              [ -90 ... +90 ]
%
% OUTPUT:
%  SP        =  Practical Salinity  (PSS-78)                   [ unitless ]
%  in_ocean  =  0, if long and lat are a long way from the ocean
%            =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several
%    hundred kilometres inland from the coast.
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and
%   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
%   Ocean Science, 8, 1123-1134.
%   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sp_from_sstar(double sstar,double p,double lon,double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double saar, sp_baltic;

	/**
	**  In the Baltic Sea, SA = Sstar.
	*/
	sp_baltic = gsw_sp_from_sa_baltic(sstar,lon,lat);

	if (sp_baltic < cppGSW_ERROR_LIMIT)
	{
	   return (sp_baltic);
   }

	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
	{
	   return (saar);
   }

	return ((sstar/gtc.gsw_ups)/(1.0 - 0.35e0*saar));
}

/***************************************************************************
% gsw_cabbeling                                       cabbeling coefficient
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_cabbeling(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the cabbeling coefficient of seawater with respect to
%  Conservative Temperature.  This function uses the computationally-
%  efficient expression for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  cabbeling  =  cabbeling coefficient with respect to            [ 1/K^2 ]
%                Conservative Temperature.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.9.2) and (P.4) of this TEOS-10 manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_cabbeling(double sa,double ct,double p)
{
	double alpha_ct, alpha_on_beta, alpha_sa, beta_sa, rho,
			 v_sa, v_ct, v_sa_sa, v_sa_ct, v_ct_ct;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct, NULL);
	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
											 NULL, NULL);

	rho = gsw_rho(sa,ct,p);
	alpha_ct = rho*(v_ct_ct - rho*v_ct*v_ct);
	alpha_sa = rho*(v_sa_ct - rho*v_sa*v_ct);
	beta_sa = -rho*(v_sa_sa - rho*v_sa*v_sa);
	alpha_on_beta = gsw_alpha_on_beta(sa,ct,p);

	return (alpha_ct +
			  alpha_on_beta*(2.0*alpha_sa - alpha_on_beta*beta_sa));
}

/***************************************************************************
% gsw_thermobaric                thermobaric coefficient (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_thermobaric(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine is based on the
%  computationally-efficient expression for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  thermobaric  =  thermobaric coefficient with                [ 1/(K Pa) ]
%                  respect to Conservative Temperature.
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_thermobaric(double sa,double ct,double p)
{
	double v_ct, v_ct_p, v_sa, v_sa_p;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,NULL);
	gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,NULL,&v_sa_p,&v_ct_p);

	return (gsw_rho(sa,ct,p)*(v_ct_p - (v_ct/v_sa)*v_sa_p));
}

/***************************************************************************
% gsw_rho_second_derivatives                             second derivatives
%                                                 of rho (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_rho_second_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of rho,
%   (1) rho_SA_SA, second-order derivative with respect to Absolute
%       Salinity at constant CT & p.
%   (2) rho_SA_CT, second-order derivative with respect to SA & CT at
%       constant p.
%   (3) rho_CT_CT, second-order derivative with respect to CT at
%       constant SA & p.
%   (4) rho_SA_P, second-order derivative with respect to SA & P at
%       constant CT.
%   (5) rho_CT_P, second-order derivative with respect to CT & P at
%       constant SA.
%
%  Note that this function uses the using the computationally-efficient
%  expression for specific volume (Roquet et al., 2015).  There is an
%  alternative to calling this function, namely
%  gsw_rho_second_derivatives_CT_exact(SA,CT,p) which uses the full Gibbs
%  function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_SA_SA = The second-order derivative of rho with respect to
%              Absolute Salinity at constant CT & p.  [ (kg/m^3)(g/kg)^-2 ]
%  rho_SA_CT = The second-order derivative of rho with respect to
%              SA and CT at constant p.          [ (kg/m^3)(g/kg)^-1 K^-1 ]
%  rho_CT_CT = The second-order derivative of rho with respect to CT at
%              constant SA & p                            [ (kg/m^3) K^-2 ]
%  rho_SA_P  = The second-order derivative with respect to SA & P at
%              constant CT.                     [ (kg/m^3)(g/kg)^-1 Pa^-1 ]
%  rho_CT_P  = The second-order derivative with respect to CT & P at
%              constant SA.                         [ (kg/m^3) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_rho_second_derivatives(double sa,double ct,double p,
                                          double *rho_sa_sa,double *rho_sa_ct,
                                          double *rho_ct_ct,double *rho_sa_p,
                                          double *rho_ct_p)
{
	double rec_v, rec_v2, rec_v3, v_ct, v_ct_ct, v_ct_p, v_p, v_sa,
			 v_sa_ct, v_sa_p, v_sa_sa;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,&v_p);
	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
											 &v_sa_p,&v_ct_p);

	rec_v = 1.0 / gsw_specvol(sa,ct,p);
	rec_v2 = pow(rec_v, 2);
	rec_v3 = rec_v2*rec_v;

	*rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;
	*rho_sa_ct = -v_sa_ct*rec_v2 + 2.0*v_sa*v_ct*rec_v3;
	*rho_ct_ct = -v_ct_ct*rec_v2 + 2.0*v_ct*v_ct*rec_v3;
	*rho_sa_p = -v_sa_p*rec_v2 + 2.0*v_sa*v_p*rec_v3;
	*rho_ct_p = -v_ct_p*rec_v2 + 2.0*v_ct*v_p*rec_v3;
}

/***************************************************************************
% gsw_pt_first_derivatives       first derivatives of potential temperature
% =========================================================================
%
% USAGE:
%  [pt_SA, pt_CT] = gsw_pt_first_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following two partial derivatives of potential
%  temperature (the regular potential temperature whose reference sea
%  pressure is 0 dbar)
%  (1) pt_SA, the derivative with respect to Absolute Salinity at
%       constant Conservative Temperature, and
%  (2) pt_CT, the derivative with respect to Conservative Temperature at
%       constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  pt_SA =  The derivative of potential temperature with respect to
%           Absolute Salinity at constant Conservative Temperature.
%                                                               [ K/(g/kg)]
%  pt_CT =  The derivative of potential temperature with respect to
%           Conservative Temperature at constant Absolute Salinity.
%           pt_CT is dimensionless.                            [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.6), (A.12.3), (P.6) and (P.8) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_pt_first_derivatives(double sa,double ct,double *pt_sa,
                                          double *pt_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double abs_pt, ct_pt, ct_sa, pt, pr0 = 0.0;

	pt = gsw_pt_from_ct(sa,ct);
	abs_pt = (gtc.gsw_t0 + pt);
	ct_pt = -(abs_pt*gsw_gibbs(0,2,0,sa,pt,pr0))/gtc.gsw_cp0;

	if (pt_sa != NULL)
	{
		ct_sa = (gsw_gibbs(1,0,0,sa,pt,pr0) -
					abs_pt*gsw_gibbs(1,1,0,sa,pt,pr0))/gtc.gsw_cp0;

		*pt_sa = -ct_sa/ct_pt;
	}

	if (pt_ct != NULL)
   {
      *pt_ct = 1.0/ct_pt;
   }
}

/***************************************************************************
% gsw_sigma2_ct_exact                        potential density anomaly with
%                                       reference sea pressure of 2000 dbar
%==========================================================================
%
% USAGE: gsw_sigma2_ct_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 2000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma2(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma2_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 2000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma2_ct_exact(double SA, double CT)
{
   double t = gsw_t_from_ct(SA, CT, 2000.0);

   return gsw_rho_t_exact(SA, t, 2000.0) - 1000.0;
}

/**************************************************************************
% gsw_sigma4_ct_exact                        potential density anomaly with
%                                       reference sea pressure of 4000 dbar
%==========================================================================
%
% USAGE: gsw_sigma4_ct_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 4000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma4(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma4_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 4000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosSea::gsw_sigma4_ct_exact(double SA, double CT)
{
	double t = gsw_t_from_ct(SA, CT, 4000.0);

	return gsw_rho_t_exact(SA, t, 4000.0) - 1000.0;
}

/***************************************************************************
% gsw_sigma1_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 1000 dbar
%==========================================================================
%
% USAGE: gsw_sigma1_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 1000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma1(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma1_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 1000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma1_ct_exact(double SA, double CT)
{
   double t = gsw_t_from_ct(SA, CT, 1000.0);

   return gsw_rho_t_exact(SA, t, 1000.0) - 1000.0;
}

/***************************************************************************
% gsw_t_from_rho_exact                     in-situ temperature from density
% =========================================================================
%
% USAGE: gsw_t_from_rho_exact(rho,SA,p, *t, *t_multiple)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg m^-3 subtracted from it.
%     That is, it is |density|, not |density anomaly|.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  t  =  in-situ temperature  (ITS-90)                            [ deg C ]
%  t_multiple  =  in-situ temperature  (ITS-90)                   [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both
%      valid solutions.  To see this second solution the user must call the
%      programme with two outputs (i.e. [t,t_multiple]), if there is only
%      one possible solution and the programme has been called with two
%      outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_t_from_rho_exact(double rho,double SA,double p,double &t,
													double &t_multiple)
{
	double sa, P, rho_40;

	/** check range, if OOR set to NAN */
	if (SA < 0.0 || SA > 42.0)
	{
	   sa = NAN;
   }
   else
   {
      sa = SA;
   }

	if (p < -1.5 || p > 12000.0) P = NAN;

	rho_40 = gsw_rho_t_exact(sa, 40.0 * sa, p);
	if (rho - rho_40 < 0.0) rho_40 = NAN;

	if (sa == NAN || P == NAN || rho_40 == NAN)
	{
		t = t_multiple = NAN;

		return;
	}

	double t_max_rho, max_rho, v_t_t, rho_t_t, deriv_lin, discrim;
	double a, b, c, v_t, rho_t, t_old, t_upper, t_freezing;

	t_max_rho = gsw_t_maxdensity_exact(sa, p);
	max_rho = gsw_rho_t_exact(sa, t_max_rho, p);
	v_t_t = gsw_gibbs(0, 2, 1, sa, t_max_rho, p);
	rho_t_t = -v_t_t * rho * rho;   //This is approximate

	deriv_lin = rho_40 - max_rho / 40.0 - t_max_rho;
	discrim = (deriv_lin * deriv_lin) - 4.0 * 0.5 * rho_t_t * (max_rho - rho);

	if (discrim < 0.0) discrim = 0.0;

	t = t_max_rho + 2.0 * (max_rho - rho / -deriv_lin + sqrt(discrim));

	/**
	   Having found this initial value of t, begin the iterative solution
	   for the t_upper part, the solution warmer than the TMD
	*/
	for (unsigned iter = 0; iter < 3; iter++)
	{
		v_t = gsw_gibbs(0,1,1,sa,t,p);
		rho_t = -v_t * rho * rho;        // This is approximate

		v_t_t = gsw_gibbs(0,2,1,sa,t,p);
		rho_t_t = -v_t_t *rho *rho;      //This is approximate

		b = rho_t;
		a = 0.5 *rho_t_t;
		c = gsw_rho_t_exact(sa,t,p) - rho;
		discrim = b * b - 4.0 * a * c;
		if (discrim < 0.0) discrim = 0.0;
		t_old = t;
		t = t_old + (2.0 * c / -b + sqrt(discrim));
	} /** for */

	t_upper = t;

	/** Now start the t_multiple part, the solution cooler than the TMD */
	t = 2.0 * t_max_rho - t_upper;

	for (unsigned iter = 0; iter < 6; iter++)
	{
		v_t = gsw_gibbs(0,1,1,sa,t,p);
		rho_t = -v_t * rho * rho;           //This is approximate

		v_t_t = gsw_gibbs(0,2,1,sa,t,p);
		rho_t_t = -v_t_t * rho * rho;       //This is approximate

		b = rho_t;
		a = 0.5 * rho_t_t;
		c = gsw_rho_t_exact(sa,t,p) - rho;
		discrim = b * b - 4.0 * a * c;
		if (discrim < 0.0) discrim = 0.0;
		t_old = t;
		t = t_old + (2.0 * c / -b - sqrt(discrim)); //Note the sign change of the sqrt term
	} /** for */

	t_multiple = t;

	/**This assumes that the seawater is unsaturated with air */
	t_freezing = gsw_t_freezing(sa, p, 0.0);

	/** Set values outside the relevant bounds to NaNs */
	if (t_upper < t_freezing || t_upper < t_max_rho) t_upper = NAN;
	if (t_multiple > t_max_rho || t_multiple > t_max_rho) t_multiple = NAN;

	t = t_upper;
}

/***************************************************************************
% gsw_enthalpy_first_derivatives_wrt_t_exact           first derivatives of
%                                                                  enthalpy
%==========================================================================
%
% USAGE: gsw_enthalpy_first_derivatives_wrt_t_exact(SA,t,p,
                                                      &h_SA_wrt_t,
                                                      &h_T_wrt_t,
                                                      &h_P_wrt_t)
%
% DESCRIPTION:
%  Calculates the following three derivatives of specific enthalpy, h.
%  These derivatives are done with respect to in-situ temperature t (in the
%  case of h_T_wrt_t) or at constant in-situ tempertature (in the cases of
%  h_SA_wrt_t and h_P_wrt_t).
%   (1) h_SA_wrt_t, the derivative with respect to Absolute Salinity at
%       constant t and p.
%   (2) h_T_wrt_t, derivative with respect to in-situ temperature t at
%       constant SA and p.
%   (3) h_P_wrt_t, derivative with respect to pressure P (in Pa) at constant
%       SA and t.  This output has the same dimensions as specific volume,
%       but is not equal to specific volume.
%
%  Note that this function uses the full Gibbs function.  This function
%  avoids the Nan that would exist in h_sub_SA at SA=0 if it were
%  evaluated in the straightforward way from the gibbs function (as in the
%  commented line 111 of the code below).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  h_SA_wrt_t =  The first derivative of specific enthalpy with respect to
%                Absolute Salinity at constant t and p.
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_T_wrt_t  =  The first derivative of specific enthalpy with respect to
%                in-situ temperature, t, at constant SA and p. [ J/(kg K) ]
%
%  h_P_wrt_t  =  The first derivative of specific enthalpy with respect to
%                pressure P (in Pa) at constant SA and t.        [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_enthalpy_first_derivatives_wrt_t_exact(double SA, double t,
                                                         double p,
                                                         double &h_SA_wrt_t,
                                                         double &h_T_wrt_t,
                                                         double &h_P_wrt_t)
{
   double gibbs_SA_T, g08_SA_T, gibbs_SA, g08_SA, x2, x, y, z, sa = 0.0;

   if (SA > 0.0) sa = SA;

   h_T_wrt_t = gsw_cp_t_exact(sa, t, p);
   h_P_wrt_t = gsw_specvol_t_exact(sa, t, p) - (gtc.gsw_t0 + t) * gsw_gibbs(0,1,1, sa, t, p);

   x2 = gtc.gsw_sfac * sa;
   x = sqrt(x2);
   y = t * 0.025;
   z = p * gtc.rec_db2pa; //1e-4; Note.The input pressure (p) is sea pressure in units of dbar.

   g08_SA = 8645.36753595126 + z *(-6620.98308089678 +
        z *(769.588305957198 + z *(-193.0648640214916 + (31.6816345533648 - 5.24960313181984 *z) *z))) +
        x *(-7296.43987145382 + x *(8103.20462414788 +
        y *(2175.341332000392 + y *(-274.2290036817964 +
        y *(197.4670779425016 + y *(-68.5590309679152 + 9.98788038278032 *y))) - 90.6734234051316 *z) +
        x *(-5458.34205214835 - 980.14153344888 *y +
        x *(2247.60742726704 - 340.1237483177863 *x + 220.542973797483 *y) + 180.142097805543 *z) +
        z *(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644 *z) *z)) +
        z *(598.378809221703 + z *(-156.8822727844005 + (204.1334828179377 - 10.23755797323846 *z) *z)) +
        y *(-1480.222530425046 + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
        y *(-129.1994027934126 + z *(1149.174198007428 + z *(-162.5751787551336 + 76.9195462169742 *z)) +
        y *(-30.0682112585625 - 1380.9597954037708 *z + y *(2.626801985426835 + 703.695562834065 *z))))) +
        y *(1187.3715515697959 + z *(1458.233059470092 +
        z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
        y *(1760.062705994408 + y *(-450.535298526802 +
        y *(182.8520895502518 + y *(-43.3206481750622 + 4.26033941694366 *y) +
        z *(-595.457483974374 + (149.452282277512 - 72.9745838003176 *z) *z)) +
        z *(1388.489628266536 + z *(-409.779283929806 + (227.123395681188 - 22.2565468652826 *z) *z))) +
        z *(-1721.528607567954 + z *(674.819060538734 +
        z *(-356.629112415276 + (88.4080716616 - 15.84003094423364 *z) *z)))));

   gibbs_SA = 0.5 * gtc.gsw_sfac * g08_SA;

   g08_SA_T = 1187.3715515697959 + z *(1458.233059470092 +
        z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
        x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
        y *(-548.4580073635929 + y *(592.4012338275047 + y *(-274.2361238716608 + 49.9394019139016 *y))) -
        90.6734234051316 *z) + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
        y *(-258.3988055868252 + z *(2298.348396014856 + z *(-325.1503575102672 + 153.8390924339484 *z)) +
        y *(-90.2046337756875 - 4142.8793862113125 *z + y *(10.50720794170734 + 2814.78225133626 *z)))) +
        y *(3520.125411988816 + y *(-1351.605895580406 +
        y *(731.4083582010072 + y *(-216.60324087531103 + 25.56203650166196 *y) +
        z *(-2381.829935897496 + (597.809129110048 - 291.8983352012704 *z) *z)) +
        z *(4165.4688847996085 + z *(-1229.337851789418 + (681.370187043564 - 66.7696405958478 *z) *z))) +
        z *(-3443.057215135908 + z *(1349.638121077468 +
        z *(-713.258224830552 + (176.8161433232 - 31.68006188846728 *z) *z))));

   gibbs_SA_T = 0.5 * gtc.gsw_sfac *0.025 * g08_SA_T;

   h_SA_wrt_t = gibbs_SA - (gtc.gsw_t0 + t) * gibbs_SA_T;
}

/***************************************************************************
% gsw_internal_energy_first_derivatives       first derivatives of specific
%                             interal energy of seawater (75-term equation)
%==========================================================================
%
% USAGE: gsw_internal_energy_first_derivatives(SA,CT,p,&u_SA,
                                                &u_CT, &u_P)
%
% DESCRIPTION:
%  Calculates the first order derivates of specific internal energy of
%  seawater using the computationally-efficient expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  u_SA = The first derivative of internal energy with respect to
%           Absolute Salinity at constant CT & p.
%                                          [ (J/kg)(g/kg)^-1 ] i.e. [ J/g ]
%  u_CT = The first derivative of internal energy with respect to
%           Conservative Temperature at constant SA & p.    [ (J/kg) K^-1 ]
%  u_P = The first derivative of internal energy with respect to
%           pressure at constant SA & CT.                  [ (J/kg) Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_internal_energy_first_derivatives(double SA, double CT,
                                                      double p,double &u_SA,
                                                      double &u_CT,
                                                      double &u_P)
{
	double v, P, h_SA, h_CT, v_SA, v_CT, v_P, sa = 0.0;

	if (SA > 0.0) sa = SA;

	P = (gtc.db2pa * p + gtc.gsw_p0);

	gsw_enthalpy_first_derivatives(sa, CT, p, &h_SA, &h_CT);
	v = gsw_specvol(sa, CT, p);
	gsw_specvol_first_derivatives(sa, CT, p, &v_SA, &v_CT, &v_P);

	u_SA = h_SA - (P * v_SA);

	u_CT = h_CT - (P * v_CT);

	u_P = v - (P * v_P);
}

/***************************************************************************
% gsw_sigma1                       potential density anomaly with reference
%                              sea pressure of 1000 dbar (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sigma1(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 1000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma1  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 1000 dbar,
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma1(double sa,double ct)
{
	return (gsw_rho(sa,ct,1000.0) - 1000.0);
}

/***************************************************************************
% gsw_sigma2                       potential density anomaly with reference
%                              sea pressure of 2000 dbar (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sigma2(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 2000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  Temperature.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma2  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 2000 dbar,
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma2(double sa,double ct)
{
	return (gsw_rho(sa,ct,2000.0) - 1000.0);
}

/***************************************************************************
% gsw_sigma3                       potential density anomaly with reference
%                              sea pressure of 3000 dbar (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sigma3(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 3000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  Temperature.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma3  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 3000 dbar,
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma3(double sa,double ct)
{
	return (gsw_rho(sa,ct,3000.0) - 1000.0);
}

/***************************************************************************
% gsw_sigma4                       potential density anomaly with reference
%                              sea pressure of 4000 dbar (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_sigma4(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 4000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  Temperature.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OUTPUT:
%  sigma4  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 4000 dbar,
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosSea::gsw_sigma4(double sa,double ct)
{
	return (gsw_rho(sa,ct,4000.0) - 1000.0);
}

/***************************************************************************
% gsw_rr68_interp_SA_CT              Reiniger and Ross (1968) interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  gsw_rr68_interp_SA_CT(SA,CT,p,p_i)
%
% DESCRIPTION:
%  Interpolate Absolute Salinity and Conservative Temperature values to
%  arbitrary pressures using the Reiniger and Ross (1968) interpolation
%  scheme.
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%  Any interpolated bottles that have pressures shallower than the
%  shallowest observed bottle are set equal to the shallowest observed
%  bottle.
%
% INPUT:
%  SA   =  Absolute Salinity                                  [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
%  p    =  sea pressure                                       [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  pressures to interpolate to.
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i.
%  CT_i = interpolated CT values at pressures p_i.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall             [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% References
%  Reiniger, R.F. and C.K. Ross, 1968: A method of interpolation with
%   application to oceanographic data. Deep-Sea Res., 15, 185-193.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosSea::gsw_rr68_interp_sa_ct(double *sa,double *ct,double *p,int mp,
                                       double *p_i,int mp_i,double *sa_i,
                                       double *ct_i)
{
	int i, j, nshallow, ncentral, ndeep, *ip=NULL, *ip_i=NULL,
         *ip_ishallow=NULL, *ip_icentral=NULL, *ip_ideep=NULL;

	char *shallow=NULL, *central=NULL, *deep=NULL;

	double *ip_shallow=NULL, *ip_central=NULL, *ip_deep=NULL, *dp=NULL, *p_ii=NULL;

	if (mp < 4)
	{
		/** need at least four bottles to perform this interpolation */
		ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
		return;
	}

	dp = new double[mp];

	for (i=1; i<mp; i++)
	{
		if ((dp[i-1] = (p[i] - p[i-1])) <= 0.0)
		{
			if (dp) delete []dp;
			ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
			return;
		}
	}

	shallow = new char[3*mp_i];
	central = shallow+mp_i;
	deep = central+mp_i;
	nshallow=ncentral=ndeep=0;
	memset(shallow, 0, 3*mp_i*sizeof (char));

	for (i=0; i<mp_i; i++)
	{
		if (p_i[i] >= p[0] && p_i[i] <= p[1])
		{
			nshallow++;
			shallow[i] = 1;
		}
		if (p_i[i] >= p[1] && p_i[i] <= p[mp-2])
		{
			ncentral++;
			central[i] = 1;
		}
		if (p_i[i] >= p[mp-2] && p_i[i] <= p[mp-1])
		{
			ndeep++;
			deep[i] = 1;
		}
	}

	if ((nshallow == 0) || (ncentral == 0) || (ndeep == 0))
	{
		if (shallow) delete []shallow;
		if (dp) delete []dp;
		ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;

		return;
	}

	ip = new int[mp+mp_i];
	ip_i = ip+mp;

	for (i=0; i<mp; i++)
	{
		ip[i] = i;
   }

	for (i=0; i<mp_i; i++)
	{
		ip_i[i] = i;
   }

	ip_ishallow = new int[nshallow+ncentral+ndeep];
	ip_icentral = ip_ishallow+nshallow;
	ip_ideep = ip_icentral+ncentral;
	ip_shallow = new double[2*(nshallow+ncentral+ndeep)];
	ip_central = ip_shallow+nshallow;
	ip_deep = ip_central+ncentral;
	p_ii = ip_deep+ndeep;

	/**
	  Calculate the 2 outer extrapolated values
	  and the inner interpolated values
	*/
	for (i=j=0; i<mp_i; i++)
	{
		if (central[i])
		{
			ip_icentral[j] = ip_i[i];
			j++;
		}
	}

	for (i=0; i<ncentral; i++)
	{
		p_ii[i] = p_i[ip_icentral[i]];
   }

	gsw_util_interp1q_int(mp,p,ip,ncentral,p_ii,ip_central);
	gsw_rr68_interp_section(0,sa,ct,p,mp,ncentral,ip_central,ip_icentral,
									p_i,sa_i,ct_i);

	for (i=j=0; i<mp_i; i++)
	{
		if (shallow[i])
		{
			ip_ishallow[j] = ip_i[i];
			j++;
		}
	}

	for (i=0; i<nshallow; i++)
	{
		p_ii[i] = p_i[ip_ishallow[i]];
   }

	gsw_util_interp1q_int(mp,p,ip,nshallow,p_ii,ip_shallow);
	gsw_rr68_interp_section(-1,sa,ct,p,mp,nshallow,ip_shallow,ip_ishallow,
									p_i,sa_i,ct_i);

	for (i=j=0; i<mp_i; i++)
	{
		if (deep[i])
		{
			ip_ideep[j] = ip_i[i];
			j++;
		}
	}

	for (i=0; i<ndeep; i++)
	{
		p_ii[i] = p_i[ip_ideep[i]];
   }

	gsw_util_interp1q_int(mp,p,ip,ndeep,p_ii,ip_deep);
	gsw_rr68_interp_section(1,sa,ct,p,mp,ndeep,ip_deep,ip_ideep,p_i,sa_i,ct_i);

	/**
	  Insert any observed bottles that are
	  at the required interpolated pressures
	*/
	for (i=0; i<mp_i; i++)
	{
		for (j=0; j<mp; j++)
		{
			if (p_i[i] == p[j])
			{
				sa_i[i] = sa[j];
				ct_i[i] = ct[j];
			}
		}
	}

	/** clean up */
	if (ip_shallow) delete []ip_shallow;
	if (ip_ishallow) delete []ip_ishallow;
	if (ip) delete []ip;
	if (shallow) delete []shallow;
	if (dp) delete []dp;
}
