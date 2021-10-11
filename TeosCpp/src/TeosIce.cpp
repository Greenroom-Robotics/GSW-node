#include "TeosIce.h"

/******************************************
  TeosIce.cpp      Version 1.06
  by Randall Kent Whited
  rkwhited@gmail.com
  ---------------------------------------
  All copyrights and all license issues
  are the same as previous versions
  ---------------------------------------
  OR converted from .m Matlab files in
  the Teos-10 Matlab version at
  http://www.teos-10.org
*******************************************/

/** class constructor : base class constructor */
TeosIce::TeosIce(): TeosBase()
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
TeosIce::~TeosIce()
{
}

/***************************************************************************
% gsw_pot_enthalpy_from_specvol_ice                 potential enthalpy from
%                                                    specific volume of ice
%==========================================================================
%
% USAGE: gsw_pot_enthalpy_from_specvol_ice(specvol_ice,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from the specific volume
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
% INPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pot_enthalpy_from_specvol_ice(double specvol_ice, double p)
{
	double rho_ice = 1.0 / specvol_ice;
	double t_ice = gsw_t_from_rho_ice(rho_ice,p);
	double pt0_ice = gsw_pt0_from_t_ice(t_ice,p);

	return gsw_pot_enthalpy_from_pt_ice(pt0_ice);
}

/***************************************************************************
% gsw_t_from_rho_ice                        temperature from density of ice
% =========================================================================
%
% USAGE: gsw_t_from_rho_ice(rho_ice,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of ice, for given values of its
%  density and sea pressure (in dbar).
%
% INPUT:
%  rho_ice =  density of ice                                     [ kg/m^3 ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  t_ice  =  in-situ temperature of ice                           [ deg C ]
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
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_t_from_rho_ice(double rho_ice, double p)
{
   double v_t_ice, t_ice, v_freezing, t_freezing, v_173_15,v_93_15,v_63_15;
   double t_ice_mean, delta_v_ice, t_ice_old, specvol_ice = 1.0 / rho_ice;

   v_173_15 = gsw_specvol_ice(-100.0 * 1.0, p);
   t_freezing = gsw_t_freezing(0.0, p, 0.0);
   v_freezing = gsw_specvol_ice(t_freezing, p);

   if (specvol_ice < v_173_15)
   {
      v_93_15 = gsw_specvol_ice(-180.0 * 1.0,p);
      v_63_15 = gsw_specvol_ice(-210.0 * 1.0,p);

      t_ice = -180.0 + 30.0 * (specvol_ice - v_93_15) / (v_93_15 - v_63_15);  // initial estimate.
      v_t_ice = (v_63_15 - v_93_15) / (-30.0); //initial estimate of v_t_ice, the t derivative of v

      t_ice = -100.0 + 80.0 * (specvol_ice - v_173_15) / (v_173_15 - v_93_15);  // initial estimate.
      v_t_ice = (v_93_15 - v_173_15) / (-80.0); //initial estimate of v_t_ice, the t derivative of v

      t_ice =  (100.0 - t_freezing * (specvol_ice -
                  v_freezing) / (v_freezing - v_173_15));  // initial estimate.

      v_t_ice = (v_freezing - v_173_15) /
                  (-100.0 - t_freezing); //initial estimate of v_t_ice, the t derivative of v
   }
   else
   {
      t_ice = (100.0 - t_freezing) * (specvol_ice -
                  v_freezing) / (v_freezing - v_173_15);  //% the initial estimate of t_ice

      v_t_ice = (v_freezing - v_173_15) /
                  (-100.0 - t_freezing); //%initial estimate of v_t_ice, the t derivative of v
   }

   /**
      --------------------------------------------------------------------------
            Begin the modified Newton-Raphson iterative procedure
      --------------------------------------------------------------------------
   */

   for (unsigned iter = 0; iter < 3; iter++) //Number_of_iterations = 1:3
   {
      t_ice_old = t_ice;
      delta_v_ice = gsw_specvol_ice(t_ice_old,p) - specvol_ice;
      t_ice = t_ice_old - delta_v_ice / v_t_ice ;

      // this is half way through the modified
      // N-R method (McDougall and Wotherspoon, 2012)

      t_ice_mean = 0.5 * (t_ice + t_ice_old);
      v_t_ice = gsw_gibbs_ice(1,1,t_ice_mean,p);
      t_ice = t_ice_old - delta_v_ice / v_t_ice;
   }

   /**
      After 3 iterations of this modified Newton-Raphson iteration,
      the error in t_ice is no larger than 2.3x10^-12 deg C, which
      is machine precision for this calculation.
   */
   return t_ice;
}

/***************************************************************************
% gsw_specvol_ice                                    specific volume of ice
%==========================================================================
%
% USAGE:
%  specvol_ice = gsw_specvol_ice(t,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_specvol_ice(double t,double p)
{
	return (gsw_gibbs_ice(0,1,t,p));
}

/***************************************************************************
% gsw_pot_enthalpy_ice_freezing_poly              potential enthalpy of ice
%                                                 at which seawater freezes
%==========================================================================
%
% USAGE:
%  gsw_pot_enthalpy_ice_freezing(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice at which seawater freezes.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing
%                              of seawater                         [ J/kg ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pot_enthalpy_ice_freezing(double sa,double p)
{
	double pt0_ice, t_freezing = gsw_t_freezing(sa,p,0.0);

	pt0_ice = gsw_pt0_from_t_ice(t_freezing,p);

	return (gsw_pot_enthalpy_from_pt_ice(pt0_ice));
}

/***************************************************************************
% gsw_pt0_from_t_ice                           potential temperature of ice
%                                       with a reference pressure of 0 dbar
% =========================================================================
%
% USAGE:
%  gsw_pt0_from_t_ice(t,p)
%
% DESCRIPTION:
%  Calculates potential temperature of ice Ih with a reference pressure of
%  0 dbar, from in-situ temperature, t.
%
% INPUT:
%  t   =  in-situ temperature  (ITS-90)                           [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pt0_ice  =  potential temperature of ice Ih with reference pressure of
%              zero dbar (ITS-90)                                 [ deg C ]
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
%    See appendix I of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pt0_from_t_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int number_of_iterations;
	double dentropy, dentropy_dt, pt0_ice,
			 pt0_ice_old, ptm_ice, true_entropy,
			 /**This is the starting polynomial for pt0 of ice Ih.*/
			 s1 = -2.256611570832386e-4,
			 s2 = -6.045305921314694e-7,
			 s3 =  5.546699019612661e-9,
			 s4 =  1.795030639186685e-11,
			 s5 =  1.292346094030742e-9,
			 p1 = -2.259745637898635e-4,
			 p2 =  1.486236778150360e-9,
			 p3 =  6.257869607978536e-12,
			 p4 = -5.253795281359302e-7,
			 p5 =  6.752596995671330e-9,
			 p6 =  2.082992190070936e-11,
			 q1 = -5.849191185294459e-15,
			 q2 =  9.330347971181604e-11,
			 q3 =  3.415888886921213e-13,
			 q4 =  1.064901553161811e-12,
			 q5 = -1.454060359158787e-10,
			 q6 = -5.323461372791532e-13;

	true_entropy = -gsw_gibbs_ice_part_t(t,p);
	if (t < -45.0 && t > -273.0)
	{
		pt0_ice = t + p*(p1 + p*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));
		if (pt0_ice < -gtc.gsw_t0) pt0_ice = -gtc.gsw_t0;
		/**
		   we add 0.05d0 to the initial estimate of pt0_ice at
		   temps less than -273 to ensure that it is never less than -273.15.
		*/
		if (pt0_ice < -273.0) pt0_ice = pt0_ice + 0.05;
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

		for (number_of_iterations = 1; number_of_iterations <= 3;
				number_of_iterations++)
		{
			pt0_ice_old = pt0_ice;
			dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
			ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
			dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		}
	}
	else
	{
		pt0_ice = t + p*(s1 + t*(s2 + t*(s3 + t*s4)) + s5*p);
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);
		pt0_ice_old = pt0_ice;
		dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
		pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
		pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	}

	if (pt0_ice < -273.0)
	{
		pt0_ice = t + p*(q1 + p*(q2 + q3*t) + t*(q4 + t*(q5 + q6*t)));
		/**
		   add 0.01d0 to the initial estimate of pt_ice used in the
		   derivative to ensure that it is never less than -273.15d0
		   because the derivative approaches zero at absolute zero.
		*/
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice+0.01);
		for (number_of_iterations = 1; number_of_iterations <= 3;
				number_of_iterations++)
		{
			pt0_ice_old = pt0_ice;
			dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
			ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
			/**
			   add 0.01d0 to the estimate of ptm_ice for temperatures less
			   than -273 to ensure that they are never less than -273.15d0
			   because the derivative approaches zero at absolute zero and
			   the addition of 0.01d0 degrees c ensures that when we divide
			   by the derivatve in the modified newton routine the method
			   does not blow up.
			*/
			ptm_ice = ptm_ice + 0.01;
			dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		}
	}
	/**
	   For temperatures less than -273.1 degsC the maximum error is less
	   than 2x10^-7 degsC. For temperatures between -273.1 and 273 the
	   maximum error is less than 8x10^-8 degsC, and for temperatures
	   greater than -273 degsC the   maximum error is 1.5x10^-12 degsC.
	   These errors are over the whole ocean depths with p varying between
	   0 and 10,000 dbar, while the in-situ temperature varied independently
	   between -273.15 and +2 degsC.
	*/
	return (pt0_ice);
}



/***************************************************************************
% gsw_gibbs_ice_part_t        part of the derivative of Gibbs energy of ice
% =========================================================================
%
% USAGE:
%  gsw_gibbs_ice_part_t(nt,np,t,p)
%
% DESCRIPTION:
%  part of the the first temperature derivative of Gibbs energy of ice
%  that is the output is gibbs_ice(1,0,t,p) + S0
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%
% OUTPUT:
%  gibbs_ice_part_t = part of temperature derivative       [ J kg^-1 K^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IAPWS, 2009: Revised Release on the Equation of State 2006 for H2O Ice
%   Ih. The International Association for the Properties of Water and
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_gibbs_ice_part_t(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  dzi, tau;
	complex<double> g, tau_t1, tau_t2, r2;

	tau = (t + gtc.gsw_t0)*gic.rec_tt;
	dzi = gtc.db2pa*p*gic.rec_pt;
	tau_t1 = tau/t1;
	tau_t2 = tau/t2;
	r2 = r20 + dzi*(r21 + r22*dzi);

	g = r1*(log((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
		 + r2*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	return (real(g));
}

/***************************************************************************
% gsw_gibbs_ice_part_pt0      part of the derivative of Gibbs energy of ice
% =========================================================================
%
% USAGE:
%  gsw_gibbs_ice_pt0(pt0)
%
% DESCRIPTION:
%  part of the the first temperature derivative of Gibbs energy of ice
%  that is the outout is "gibbs_ice(1,0,pt0,0) + s0"
%
% INPUT:
%  pt0  =  potential temperature with reference sea pressure of zero dbar
%                                                                 [ deg C ]
%
% OUTPUT:
%  gibbs_ice_part_pt0 = part of temperature derivative     [ J kg^-1 K^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IAPWS, 2009: Revised Release on the Equation of State 2006 for H2O Ice
%   Ih. The International Association for the Properties of Water and
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_gibbs_ice_pt0(double pt0)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  tau;
	complex<double> g, tau_t1, tau_t2;

	tau = (pt0 + gtc.gsw_t0)*gic.rec_tt;
	tau_t1 = tau/t1;
	tau_t2 = tau/t2;

	g = r1*(log((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
		 + r20*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	return (real(g));
}

/***************************************************************************
% gsw_gibbs_ice_pt0_pt0       part of the derivative of Gibbs energy of ice
% =========================================================================
%
% USAGE:
%  gsw_gibbs_ice_pt0_pt0(pt0)
%
% DESCRIPTION:
%  The second temperature derivative of Gibbs energy of ice at the
%  potential temperature with reference sea pressure of zero dbar.  That is
%  the output is gibbs_ice(2,0,pt0,0).
%
% INPUT:
%  pt0  =  potential temperature with reference sea pressure of zero dbar
%                                                                 [ deg C ]
%
% OUTPUT:
%  gibbs_ice_pt0_pt0 = temperature second derivative at pt0
%                                                          [ J kg^-1 K^-2 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IAPWS, 2009: Revised Release on the Equation of State 2006 for H2O Ice
%   Ih. The International Association for the Properties of Water and
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_gibbs_ice_pt0_pt0(double pt0)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  tau;
	complex<double> g;

	tau = (pt0 + gtc.gsw_t0)*gic.rec_tt;

	g = r1*(1.0/(t1 - tau) + 1.0/(t1 + tau) - 2.0/t1)
		 + r20*(1.0/(t2 - tau) + 1.0/(t2 + tau) - 2.0/t2);

	return (gic.rec_tt*real(g));
}


/***************************************************************************
% gsw_entropy_ice                                   specific entropy of ice
%==========================================================================
%
% USAGE:
%  gsw_entropy_ice(t,p)
%
% DESCRIPTION:
%  Calculates specific entropy of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  ice_entropy  =  specific entropy of ice                 [ J kg^-1 K^-1 ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_entropy_ice(double t,double p)
{
	return (-gsw_gibbs_ice(1,0,t,p));
}

/***************************************************************************
% gsw_enthalpy_ice                                 specific enthalpy of ice
%==========================================================================
%
% USAGE:
%  gsw_enthalpy_ice(t,p)
%
% DESCRIPTION:
%  Calculates the specific enthalpy of ice (h_Ih).
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  enthalpy_ice  =  specific enthalpy of ice                       [ J/kg ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_enthalpy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double tau, dzi, g0;
	complex<double> r2, sqtau_t1, sqtau_t2, g;

	tau = (t + gtc.gsw_t0)*gic.rec_tt;
	dzi = gtc.db2pa*p*gic.rec_pt;
	g0 = gic.g00 + dzi*(gic.g01 + dzi*(gic.g02 + dzi*(gic.g03 + gic.g04*dzi)));
	r2 = r20 + dzi*(r21 + r22*dzi);
	sqtau_t1 = (tau*tau)/(t1*t1);
	sqtau_t2 = (tau*tau)/(t2*t2);
	g = r1*t1*(log(1.0 - sqtau_t1) + sqtau_t1)
		 + r2*t2*(log(1.0 - sqtau_t2) + sqtau_t2);

	return (g0 + gic.tt*real(g));
}

/**************************************************************************
% gsw_latentheat_melting                             latent heat of melting
%==========================================================================
%
% USAGE:
%  gsw_latentheat_melting(SA,p)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of melting.  It is defined in terms
%  of Absolute Salinity, SA, and sea pressure, p, and is valid in the
%  ranges 0 < SA < 42 g kg^-1 and 0 < p < 10,000 dbar.  This is based on
%  the IAPWS Releases IAPWS-09 (for pure water), IAPWS-08 (for the saline
%  compoonent of seawater and IAPWS-06 for ice Ih.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  latentheat_melting  =  latent heat of melting                   [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
% IAPWS, 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic
%  Properties of Seawater. The International Association for the Properties
%  of Water and Steam. Berlin, Germany, September 2008.  This Release is
%  known as IAPWS-09.
%
% IAPWS, 2009a: Revised Release on the Equation of State 2006 for H2O Ice
%  Ih. The International Association for the Properties of Water and Steam.
%  Doorwerth, The Netherlands, September 2009. This Release is known as
%  IAPWS-06
%
% IAPWS, 2009b: Supplementary Release on a Computationally Efficient
%  Thermodynamic Formulation for Liquid Water for Oceanographic Use. The
%  International Association for the Properties of Water and Steam.
%  Doorwerth, The Netherlands, September 2009.  This Release is known as
%  IAPWS-09.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.34 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_latentheat_melting(double sa,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double tf = gsw_t_freezing(sa,p,0.0);

	return (1000.0*(gsw_chem_potential_water_t_exact(sa,tf,p)
						   - (gtc.gsw_t0 + tf)
						   *gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p))
			            - gsw_enthalpy_ice(tf,p));
}


/***************************************************************************
% gsw_t_freezing_first_derivatives         first derivatives of the in-situ
%                                     temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  gsw_t_freezing_first_derivatives(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the in-situ temperature at which
%  seawater freezes with respect to Absolute Salinity SA and pressure P (in
%  Pa).  These expressions come from differentiating the expression that
%  defines the freezing temperature, namely the equality between the
%  chemical potentials of water in seawater and in ice.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  tfreezing_SA = the derivative of the in-situ freezing temperature
%                 (ITS-90) with respect to Absolute Salinity at fixed
%                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  tfreezing_P  = the derivative of the in-situ freezing temperature
%                 (ITS-90) with respect to pressure (in Pa) at fixed
%                 Absolute Salinity                                [ K/Pa ]
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_t_freezing_first_derivatives(double sa,double p,
                                                double saturation_fraction,
                                                double *tfreezing_sa,
                                                double *tfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double rec_denom,tf,g_per_kg = 1000.0;

	tf = gsw_t_freezing(sa,p,saturation_fraction);
	rec_denom = 1.0/
					(g_per_kg*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p)
					 + gsw_entropy_ice(tf,p));
	if (tfreezing_sa != NULL)
		*tfreezing_sa =
			gsw_dilution_coefficient_t_exact(sa,tf,p)*rec_denom
			+ saturation_fraction*(1e-3)/(2.0*gtc.gsw_sso);
	if (tfreezing_p != NULL)
		*tfreezing_p =
			-(gsw_specvol_t_exact(sa,tf,p) - sa*gsw_gibbs(1,0,1,sa,tf,p)
			  - gsw_specvol_ice(tf,p))*rec_denom;
}

/***************************************************************************
% gsw_t_freezing_first_derivatives_poly            first derivatives of the
%                             in-situ temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  gsw_t_freezing_first_derivatives_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the frist derivatives of the in-situ temperature at which
%  seawater freezes with respect to Absolute Salinity SA and pressure P (in
%  Pa).  These expressions come from differentiating the expression that
%  defines the freezing temperature, namely the equality between the
%  chemical potentials of water in seawater and in ice.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  tfreezing_SA = the derivative of the in-situ freezing temperature
%                 (ITS-90) with respect to Absolute Salinity at fixed
%                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  tfreezing_P  = the derivative of the in-situ freezing temperature
%                 (ITS-90) with respect to pressure (in Pa) at fixed
%                 Absolute Salinity                                [ K/Pa ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_t_freezing_first_derivatives_poly(double sa,double p,
                                                      double saturation_fraction,
                                                      double *tfreezing_sa,
                                                      double *tfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, c = 1e-3/(2.0*gtc.gsw_sso);

	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	if (tfreezing_sa != NULL)
	{
		*tfreezing_sa = (gfpc.t1
		                  + x*(1.5*gfpc.t2 + x*(2.0*gfpc.t3 + x*(2.5*gfpc.t4 + x*(3.0*gfpc.t5
								+ 3.5*gfpc.t6*x)))) + p_r*(gfpc.t10 + x*(1.5*gfpc.t11 + x*(2.0*gfpc.t13
								+ x*(2.5*gfpc.t16 + x*(3.0*gfpc.t19 + 3.5*gfpc.t22*x))))
								+ p_r*(gfpc.t12 + x*(1.5*gfpc.t14 + x*(2.0*gfpc.t17 + 2.5*gfpc.t20*x))
								+ p_r*(gfpc.t15 + x*(1.5*gfpc.t18 + 2.0*gfpc.t21*x)))))*1e-2
			               + saturation_fraction*c;
	}

	if (tfreezing_p != NULL)
	{
		*tfreezing_p = (gfpc.t7
		                  + sa_r*(gfpc.t10 + x*(gfpc.t11 + x*(gfpc.t13 + x*(gfpc.t16
		                  + x*(gfpc.t19 + gfpc.t22*x))))) + p_r*(2.0*gfpc.t8 + sa_r*(2.0*gfpc.t12
		                  + x*(2.0*gfpc.t14 + x*(2.0*gfpc.t17 + 2.0*gfpc.t20*x)))
		                  + p_r*(3.0*gfpc.t9 + sa_r*(3.0*gfpc.t15 + x*(3.0*gfpc.t18
								+ 3.0*gfpc.t21*x)))))*1e-8;
	}
}


/***************************************************************************
% gsw_t_freezing_poly         in-situ temperature at which seawater freezes
%                                                                    (poly)
%==========================================================================
%
% USAGE:
%  gsw_t_freezing_poly(SA,p,saturation_fraction,polynomial)
%
% DESCRIPTION:
%  Calculates the in-situ temperature at which seawater freezes from a
%  comptationally efficient polynomial.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
%               (ITS-90)
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rainer Feistal    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_t_freezing_poly(double sa,double p,
                                       double saturation_fraction,
                                       int polynomial)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, ctf, sfrac, return_value;
	int direct_poly=polynomial;

	if (direct_poly)
	{
		sfrac = saturation_fraction;
		ctf = gsw_ct_freezing_poly(sa,p,sfrac);
		return_value = gsw_t_from_ct(sa,ctf,p);
	}
	else
	{
		/** Alternative calculation ... */
		sa_r = sa*1e-2;
		x = sqrt(sa_r);
		p_r = p*1e-4;
		return_value = gfpc.t0
							+ sa_r*(gfpc.t1 + x*(gfpc.t2 + x*(gfpc.t3 + x*(gfpc.t4
							+ x*(gfpc.t5 + gfpc.t6*x))))) + p_r*(gfpc.t7 + p_r*(gfpc.t8 + gfpc.t9*p_r))
							+ sa_r*p_r*(gfpc.t10 + p_r*(gfpc.t12 + p_r*(gfpc.t15 + gfpc.t21*sa_r))
							+ sa_r*(gfpc.t13 + gfpc.t17*p_r + gfpc.t19*sa_r)
							+ x*(gfpc.t11 + p_r*(gfpc.t14 + gfpc.t18*p_r) + sa_r*(gfpc.t16
							+ gfpc.t20*p_r + gfpc.t22*sa_r)));

		/** Adjust for the effects of dissolved air */
		return_value = return_value -
							saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gtc.gsw_sso));
	}

	return (return_value);
}


/***************************************************************************
% gsw_t_freezing_poly         in-situ temperature at which seawater freezes
%                                                                    (poly)
%==========================================================================
%
% USAGE:
%  gsw_t_freezing_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the in-situ temperature at which seawater freezes from a
%  comptationally efficient polynomial.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
%               (ITS-90)
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rainer Feistal    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_t_freezing_poly(double sa,double p,
                                       double saturation_fraction)
{
	double ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);

	return gsw_t_from_ct(sa,ctf,p);
}

/***************************************************************************
% gsw_frazil_properties_potential_poly      Absolute Salinity, Conservative
%                               Temperature, and the ice mass fraction when
%              at thermodynamic equilibrium between seawater and ice (poly)
%==========================================================================
%
% USAGE:
%   gsw_frazil_properties_potential_poly(SA_bulk,h_pot_bulk,p)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), w_Ih_eq, which results from given values of the bulk
%  Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
%  occuring at pressure p.  The final equilibrium values of Absolute
%  Salinity, SA_eq, and Conservative Temperature, CT_eq, of the
%  interstitial seawater phase are also returned.  This code assumes that
%  there is no dissolved air in the seawater (that is, saturation_fraction
%  is assumed to be zero thoughout the code).
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
%  present in the final state.  In this case the final state consists of
%  only seawater rather than being an equlibrium mixture of seawater and
%  ice which occurs when w_Ih_final is positive.  Note that when
%  w_Ih_final = 0, the final seawater is not at the freezing temperature.
%
%  Note that this code uses the polynomial forms of CT_freezing and
%  pot_enthalpy_ice_freezing.  This code is intended to be used in ocean
%  models where the model prognostic variables are SA_bulk and h_pot_bulk.
%
% INPUT:
%  SA_bulk     =  bulk Absolute Salinity of the seawater and ice mixture
%                                                                  [ g/kg ]
%  h_pot_bulk  =  bulk potential enthalpy of the seawater and ice mixture
%                                                                  [ J/kg ]
%  p           =  sea pressure                                  [ dbar ]
%                  ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state,
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at
%                 thermodynamic equilibrium.  If this ice mass fraction is
%                 zero there is no ice in the final state which consists
%                 only of seawater which is warmer than the freezing
%                 temperature.                                   [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Anonymous, 2014: Modelling the interaction between seawater and frazil
%   ice.  Manuscript, March 2015.  See Eqns. (8)-(15) of this  manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_frazil_properties_potential_poly(double sa_bulk,double h_pot_bulk,
                                                      double p,double *sa_final,
                                                      double *ct_final,double *w_ih_final)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int     iterations, max_iterations;

	double  ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly, dpot_h_ihf_dsa,
			  func, func0, h_pot_ihf, sa, w_ih_old, w_ih, x, xa, y, z;

	double  f01 = -9.041191886754806e-1,
			  f02 =  4.169608567309818e-2,
			  f03 = -9.325971761333677e-3,
			  f04 =  4.699055851002199e-2,
			  f05 = -3.086923404061666e-2,
			  f06 =  1.057761186019000e-2,
			  f07 = -7.349302346007727e-2,
			  f08 =  1.444842576424337e-1,
			  f09 = -1.408425967872030e-1,
			  f10 =  1.070981398567760e-1,
			  f11 = -1.768451760854797e-2,
			  f12 = -4.013688314067293e-1,
			  f13 =  7.209753205388577e-1,
			  f14 = -1.807444462285120e-1,
			  f15 =  1.362305015808993e-1,
			  f16 = -9.500974920072897e-1,
			  f17 =  1.192134856624248,
			  f18 = -9.191161283559850e-2,
			  f19 = -1.008594411490973,
			  f20 =  8.020279271484482e-1,
			  f21 = -3.930534388853466e-1,
			  f22 = -2.026853316399942e-2,
			  f23 = -2.722731069001690e-2,
			  f24 =  5.032098120548072e-2,
			  f25 = -2.354888890484222e-2,
			  f26 = -2.454090179215001e-2,
			  f27 =  4.125987229048937e-2,
			  f28 = -3.533404753585094e-2,
			  f29 =  3.766063025852511e-2,
			  f30 = -3.358409746243470e-2,
			  f31 = -2.242158862056258e-2,
			  f32 =  2.102254738058931e-2,
			  f33 = -3.048635435546108e-2,
			  f34 = -1.996293091714222e-2,
			  f35 =  2.577703068234217e-2,
			  f36 = -1.292053030649309e-2,
			  g01 =  3.332286683867741e5,
			  g02 =  1.416532517833479e4,
			  g03 = -1.021129089258645e4,
			  g04 =  2.356370992641009e4,
			  g05 = -8.483432350173174e3,
			  g06 =  2.279927781684362e4,
			  g07 =  1.506238790315354e4,
			  g08 =  4.194030718568807e3,
			  g09 = -3.146939594885272e5,
			  g10 = -7.549939721380912e4,
			  g11 =  2.790535212869292e6,
			  g12 =  1.078851928118102e5,
			  g13 = -1.062493860205067e7,
			  g14 =  2.082909703458225e7,
			  g15 = -2.046810820868635e7,
			  g16 =  8.039606992745191e6,
			  g17 = -2.023984705844567e4,
			  g18 =  2.871769638352535e4,
			  g19 = -1.444841553038544e4,
			  g20 =  2.261532522236573e4,
			  g21 = -2.090579366221046e4,
			  g22 = -1.128417003723530e4,
			  g23 =  3.222965226084112e3,
			  g24 = -1.226388046175992e4,
			  g25 =  1.506847628109789e4,
			  g26 = -4.584670946447444e4,
			  g27 =  1.596119496322347e4,
			  g28 = -6.338852410446789e4,
			  g29 =  8.951570926106525e4,
			  saturation_fraction = 0.0;
	/**
	    -----------------------------------------------------------------------
	      Finding func0.  This is the value of the function, func, that would
	      result in the output w_Ih_final being exactly zero.
	    -----------------------------------------------------------------------
	*/
	func0 = h_pot_bulk - gtc.gsw_cp0
			  *gsw_ct_freezing_poly(sa_bulk,p,saturation_fraction);
	/**
	    -----------------------------------------------------------------------
	      Setting the three outputs for data points that have func0 non-negative
	    -----------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		/**
		      When func0 is zero or positive then the final answer will contain
		      no frazil ice; that is, it will be pure seawater that is warmer
		      han the freezing temperature.  If func0 >= 0 we do not need to go
		      through the modified Newton-Raphson procedure and we can simply
		      write down the answer, as in the following 4 lines of code.
		*/
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	    Begin finding the solution for data points that have func0 < 0, so that
	    the output will have a positive ice mass fraction w_Ih_final.
	    -----------------------------------------------------------------------
	*/
	/**Evalaute a polynomial for w_Ih in terms of SA_bulk, func0 and p*/
	x = sa_bulk*1e-2;
	y = func0/3e5;
	z = p*1e-4;
	w_ih = y*(f01
	         + x*(f02 + x*(f03 + x*(f04 + x*(f05 + f06*x))))
            + y*(f07 + x*(f08 + x*(f09 + x*(f10 + f11*x))) + y*(f12 + x*(f13
				+ x*(f14 + f15*x)) + y*(f16 + x*(f17 + f18*x) + y*(f19 + f20*x
				+ f21*y)))) + z*(f22 + x*(f23 + x*(f24 + f25*x)) + y*(x*(f26
				+ f27*x) + y*(f28 + f29*x + f30*y)) + z*(f31 + x*(f32 + f33*x)
				+ y*(f34 + f35*x + f36*y))));

	if (w_ih > 0.9)
	{
		/**
		    The ice mass fraction out of this code is restricted to be
		    less than 0.9.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	     The initial guess at the absolute salinity of the interstitial
	     seawater
	*/
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Doing a Newton step with a separate polynomial estimate of the mean
	      derivative dfunc_dw_Ih_mean_poly.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
	func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
	xa = sa*1e-2;
	dfunc_dw_ih_mean_poly = g01
	                        + xa*(g02 + xa*(g03 + xa*(g04 + g05*xa)))
									+ w_ih*(xa*(g06 + xa*(g07 + g08*xa)) + w_ih*(xa*(g09 + g10*xa)
									+ w_ih*xa*(g11 + g12*xa + w_ih*(g13 + w_ih*(g14 + w_ih*(g15
									+ g16*w_ih)))))) + z*(g17 + xa*(g18 + g19*xa) + w_ih*(g20
									+ w_ih*(g21 + g22*w_ih) + xa*(g23 + g24*xa*w_ih))
									+ z*(g25 + xa*(g26 + g27*xa) + w_ih*(g28 + g29*w_ih)));

	w_ih_old = w_ih;
	w_ih = w_ih_old - func/dfunc_dw_ih_mean_poly;
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
	      fed into Newton s Method.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,&ctf_sa,
														NULL);
	gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa,p,
			&dpot_h_ihf_dsa, NULL);
	dfunc_dw_ih = gtc.gsw_cp0*ctf - h_pot_ihf -
					  sa*(gtc.gsw_cp0*ctf_sa + w_ih*dpot_h_ihf_dsa/(1.0 - w_ih));

	if (w_ih >= 0.0 && w_ih <= 0.20 && sa > 15.0
			&& sa < 60.0 && p <= 3000.0)
	{
		max_iterations = 1;
	}
	else if (w_ih >= 0.0 && w_ih <= 0.85 && sa > 0.0
				&& sa < 120.0 && p <= 3500.0)
	{
		max_iterations = 2;
	}
	else
	{
		max_iterations = 3;
	}

	for (iterations = 1; iterations <= max_iterations; iterations++)
	{
		if (iterations > 1)
		{
			/**On the first iteration ctf and h_pot_ihf are both known*/
			ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
			h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
		}
		/**This is the function, func, whose zero we seek ...*/
		func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */

	if (w_ih < 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
	}
	else
	{
		*sa_final = sa;
		*ct_final = gsw_ct_freezing_poly(sa,p,saturation_fraction);
		*w_ih_final = w_ih;
	}
}




/***************************************************************************
% gsw_ct_from_rho                     Conservative Temperature from density
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  gsw_ct_from_rho(rho,SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar),
%  using the computationally-efficient expression for specific volume in
%  terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  CT  =  Conservative Temperature  (ITS-90)                      [ deg C ]
%  CT_multiple  =  Conservative Temperature  (ITS-90)             [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      Conservative Temperatures for a single density.  This programme will
%      output both valid solutions.  To see this second solution the user
%      must call the programme with two outputs (i.e. [CT,CT_multiple]), if
%      there is only one possible solution and the programme has been
%      called with two outputs the second variable will be set to NaN.
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
void TeosIce::gsw_ct_from_rho(double rho,double sa,double p,double *ct,
                                 double *ct_multiple)
{
	int number_of_iterations;
	double a, alpha_freezing, alpha_mean, b, c, ct_a, ct_b, ct_diff,
			 ct_freezing, ct_max_rho, ct_mean, ct_old,
			 delta_ct, delta_v, factor, factorqa, factorqb,
			 rho_40, rho_extreme, rho_freezing, rho_max, rho_mean,
			 rho_old, sqrt_disc, top, v_ct, v_lab;
	/**
	   alpha_limit is the positive value of the thermal expansion coefficient
	   which is used at the freezing temperature to distinguish between
	   salty and fresh water.
	*/
	double alpha_limit = 1e-5;
	/**
	   rec_half_rho_TT is a constant representing the reciprocal of half the
	   second derivative of density with respect to temperature near the
	   temperature of maximum density.
	*/
	double rec_half_rho_tt = -110.0;
	rho_40 = gsw_rho(sa,40.0,p);
	if (rho < rho_40)
	{
		*ct = cppGSW_INVALID_VALUE;
		if (ct_multiple != NULL) *ct_multiple = *ct;
		return;
	}
	ct_max_rho = gsw_ct_maxdensity(sa,p);
	rho_max = gsw_rho(sa,ct_max_rho,p);
	rho_extreme = rho_max;
	/** Assumes that the seawater is always unsaturated with air */
	ct_freezing = gsw_ct_freezing_poly(sa,p,0.0);
	gsw_rho_alpha_beta(sa,ct_freezing,p,&rho_freezing,&alpha_freezing,NULL);
	/** reset the extreme values */
	if (ct_freezing > ct_max_rho)
	{
		rho_extreme = rho_freezing;
	}
	if (rho > rho_extreme)
	{
		*ct = cppGSW_INVALID_VALUE;
		if (ct_multiple != NULL) *ct_multiple = *ct;
		return;
	}
	if (alpha_freezing > alpha_limit)
	{
		ct_diff = 40.0 - ct_freezing;
		top = rho_40 - rho_freezing + rho_freezing*alpha_freezing*ct_diff;
		a = top/(ct_diff*ct_diff);
		b = -rho_freezing*alpha_freezing;
		c = rho_freezing - rho;
		sqrt_disc = sqrt(b*b - 4*a*c);
		*ct = ct_freezing + 0.5*(-b - sqrt_disc)/a;
	}
	else
	{
		ct_diff = 40.0 - ct_max_rho;
		factor = (rho_max - rho)/(rho_max - rho_40);
		delta_ct = ct_diff*sqrt(factor);
		if (delta_ct > 5.0)
		{
			*ct = ct_max_rho + delta_ct;
		}
		else
		{
			/** Set the initial value of the quadratic solution roots. */
			ct_a = ct_max_rho + sqrt(rec_half_rho_tt*(rho - rho_max));
			for (number_of_iterations = 1; number_of_iterations <= 7;
					number_of_iterations++)
			{
				ct_old = ct_a;
				rho_old = gsw_rho(sa,ct_old,p);
				factorqa = (rho_max - rho)/(rho_max - rho_old);
				ct_a = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqa);
			}
			if ((ct_freezing - ct_a) < 0.0)
			{
				*ct = cppGSW_INVALID_VALUE;
				if (ct_multiple != NULL) *ct_multiple = *ct;
				return;
			}
			*ct = ct_a;
			if (ct_multiple == NULL)
			{
				return;
			}
			/** Set the initial value of the quadratic solution roots. */
			ct_b = ct_max_rho - sqrt(rec_half_rho_tt*(rho - rho_max));
			for (number_of_iterations = 1; number_of_iterations <= 7;
					number_of_iterations++)
			{
				ct_old = ct_b;
				rho_old = gsw_rho(sa,ct_old,p);
				factorqb = (rho_max - rho)/(rho_max - rho_old);
				ct_b = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqb);
			}
			/**
			   After seven iterations of this quadratic iterative procedure,
			   the error in rho is no larger than 4.6x10^-13 kg/m^3.
			*/
			if ((ct_freezing - ct_b) < 0.0)
			{
				*ct = cppGSW_INVALID_VALUE;
				*ct_multiple = *ct;
				return;
			}
			*ct_multiple = ct_b;
			return;
		} /** if delta ct > 5 else */
	} /** if alpha_freezing else */
	/** Begin the modified Newton-Raphson iterative function */
	v_lab = 1.0/rho;
	gsw_rho_alpha_beta(sa,*ct,p,&rho_mean,&alpha_mean,NULL);
	v_ct = alpha_mean/rho_mean;
	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		ct_old = *ct;
		delta_v = gsw_specvol(sa,ct_old,p) - v_lab;
		*ct = ct_old - delta_v/v_ct;
		ct_mean = 0.5*(*ct + ct_old);
		gsw_rho_alpha_beta(sa,ct_mean,p,&rho_mean,&alpha_mean,NULL);
		v_ct = alpha_mean/rho_mean;
		*ct = ct_old - delta_v/v_ct ;
	}
	/**
	   After three iterations of this modified Newton-Raphson iteration,
	   the error in rho is no larger than 1.6x10^-12 kg/m^3.
	*/
	if (ct_multiple != NULL)
	{
		*ct_multiple = cppGSW_INVALID_VALUE;
	}
	return;
}

/***************************************************************************
% gsw_ct_freezing_first_derivatives       first derivatives of Conservative
%                                     Temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  gsw_ct_freezing_first_derivatives(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the Conservative Temperature at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, completely unsaturated)
%
% OUTPUT:
%  CTfreezing_SA = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to Absolute Salinity at
%                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  CTfreezing_P  = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to pressure (in Pa) at
%                  fixed Absolute Salinity                         [ K/Pa ]
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_ct_freezing_first_derivatives(double sa,double p,
                                                   double saturation_fraction,
                                                   double *ctfreezing_sa,
                                                   double *ctfreezing_p)
{
	double tf_sa, tf_p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t, tf;

	tf = gsw_t_freezing(sa,p,saturation_fraction);

	if (ctfreezing_sa != NULL && ctfreezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa,&tf_p);

		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 &ct_sa_wrt_t,&ct_t_wrt_t,&ct_p_wrt_t);

		*ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;
		*ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;
	}
	else if (ctfreezing_sa != NULL && ctfreezing_p == NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa, NULL);

		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 &ct_sa_wrt_t,&ct_t_wrt_t,NULL);

		*ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;
	}
	else if (ctfreezing_sa == NULL && ctfreezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													NULL, &tf_p);

		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 NULL,&ct_t_wrt_t,&ct_p_wrt_t);

		*ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;
	}
}

/***************************************************************************
% gsw_ct_freezing_first_derivatives_poly               first derivatives of
%                 Conservative Temperature at which seawater freezes (poly)
%==========================================================================
%
% USAGE:
%  gsw_ct_freezing_first_derivatives_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the Conservative Temperature at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa) of the comptationally efficient polynomial fit of the
%  freezing temperature (McDougall et al., 2014).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  CTfreezing_SA = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to Absolute Salinity at
%                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  CTfreezing_P  = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to pressure (in Pa) at
%                  fixed Absolute Salinity                         [ K/Pa ]
%
% AUTHOR:
%  Trevor McDougall, Paul Barker  [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_ct_freezing_first_derivatives_poly(double sa,double p,
                                                      double saturation_fraction,
                                                      double *ctfreezing_sa,
                                                      double *ctfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x,
			 d = -gfpc.a - gfpc.a*gfpc.b - 2.4*gfpc.b/gtc.gsw_sso,
			 e = 2.0*gfpc.a*gfpc.b/gtc.gsw_sso;

	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	if (ctfreezing_sa != NULL)
	{
		*ctfreezing_sa = (gfpc.c1
		                  + x*(1.5*gfpc.c2 + x*(2.0*gfpc.c3 + x*(2.5*gfpc.c4 + x*(3.0*gfpc.c5
								+ 3.5*gfpc.c6*x)))) + p_r*(gfpc.c10 + x*(1.5*gfpc.c11 + x*(2.0*gfpc.c13
								+ x*(2.5*gfpc.c16 + x*(3.0*gfpc.c19 + 3.5*gfpc.c22*x))))
								+ p_r*(gfpc.c12 + x*(1.5*gfpc.c14 + x*(2.0*gfpc.c17 + 2.5*gfpc.c20*x))
								+ p_r*(gfpc.c15 + x*(1.5*gfpc.c18
								+ 2.0*gfpc.c21*x)))))*1e-2 - saturation_fraction*1e-3*(d - sa*e);
	}

	if (ctfreezing_p != NULL)
	{
		*ctfreezing_p = (gfpc.c7
		                  + sa_r*(gfpc.c10 + x*(gfpc.c11 + x*(gfpc.c13 + x*(gfpc.c16
								+ x*(gfpc.c19 + gfpc.c22*x))))) + p_r*(2.0*gfpc.c8
								+ sa_r*(2.0*gfpc.c12 + x*(2.0*gfpc.c14 + x*(2.0*gfpc.c17
								+ 2.0*gfpc.c20*x))) + p_r*(3.0*gfpc.c9 + sa_r*(3.0*gfpc.c15
								+ x*(3.0*gfpc.c18 + 3.0*gfpc.c21*x)))))*1e-8;
	}
}

/***************************************************************************
% gsw_CT_freezing_poly                    Conservative Temperature at which
%                                                   seawater freezes (poly)
%==========================================================================
%
% USAGE:
%  gsw_CT_freezing_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature at which seawater freezes.
%  The error of this fit ranges between -5e-4 K and 6e-4 K when compared
%  with the Conservative Temperature calculated from the exact in-situ
%  freezing temperature which is found by a Newton-Raphson iteration of the
%  equality of the chemical potentials of water in seawater and in ice.
%  Note that the Conservative Temperature freezing temperature can be found
%  by this exact method using the function gsw_CT_freezing.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, completely unsaturated)
%
% OUTPUT:
%  CT_freezing = Conservative Temperature at freezing of seawater [ deg C ]
%                That is, the freezing temperature expressed in
%                terms of Conservative Temperature (ITS-90).
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rainer Feistal    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_ct_freezing_poly(double sa,double p,double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, return_value;

	sa_r = sa*1.0e-2;
	x = sqrt(sa_r);
	p_r = p*1.0e-4;

	return_value = gfpc.c0
						+ sa_r*(gfpc.c1 + x*(gfpc.c2 + x*(gfpc.c3 + x*(gfpc.c4 + x*(gfpc.c5 + gfpc.c6*x)))))
						+ p_r*(gfpc.c7 + p_r*(gfpc.c8 + gfpc.c9*p_r)) + sa_r*p_r*(gfpc.c10 + p_r*(gfpc.c12
                  + p_r*(gfpc.c15 + gfpc.c21*sa_r)) + sa_r*(gfpc.c13 + gfpc.c17*p_r + gfpc.c19*sa_r)
                  + x*(gfpc.c11 + p_r*(gfpc.c14 + gfpc.c18*p_r) + sa_r*(gfpc.c16 + gfpc.c20*p_r
                  + gfpc.c22*sa_r)));

	/** Adjust for the effects of dissolved air */
	return_value = return_value - saturation_fraction*
						(1e-3)*(2.4 - gfpc.a*sa)*(1.0 + gfpc.b*(1.0 - sa/gtc.gsw_sso));

	return (return_value);
}


/**************************************************************************
% gsw_frazil_properties_potential           Absolute Salinity, Conservative
%                               Temperature, and the ice mass fraction when
%                     at thermodynamic equilibrium between seawater and ice
%==========================================================================
%
% USAGE:
%  gsw_frazil_properties_potential(SA_bulk,h_pot_bulk,p)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), w_Ih_eq, which results from given values of the bulk
%  Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
%  occuring at pressure p.  The final equilibrium values of Absolute
%  Salinity, SA_eq, and Conservative Temperature, CT_eq, of the
%  interstitial seawater phase are also returned.  This code assumes that
%  there is no dissolved air in the seawater (that is, saturation_fraction
%  is assumed to be zero thoughout the code).
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
%  present in the final state.  In this case the final state consists of
%  only seawater rather than being an equlibrium mixture of seawater and
%  ice which occurs when w_Ih_final is positive.  Note that when
%  w_Ih_final = 0, the final seawater is not at the freezing temperature.
%
%  Note that this code uses the exact forms of CT_freezing and
%  pot_enthalpy_ice_freezing.
%
% INPUT:
%  SA_bulk     =  bulk Absolute Salinity of the seawater and ice mixture
%                                                                  [ g/kg ]
%  h_pot_bulk  =  bulk potential enthalpy of the seawater and ice mixture
%                                                                  [ J/kg ]
%  p           =  sea pressure                                  [ dbar ]
%                  ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state,
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at
%                 thermodynamic equilibrium.  If this ice mass fraction is
%                 zero there is no ice in the final state which consists
%                 only of seawater which is warmer than the freezing
%                 temperature.                                   [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Anonymous, 2014: Modelling the interaction between seawater and frazil
%   ice.  Manuscript, March 2015.  See Eqns. (8)-(15) of this  manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_frazil_properties_potential(double sa_bulk,double h_pot_bulk,
                                                double p,double *sa_final,
                                                double *ct_final,double *w_ih_final)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int iterations, max_iterations;
	double ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
			 dpot_h_ihf_dsa, func, func0, h_pot_ihf, sa, w_ih_old, w_ih,
			 x, xa, y, z;

	double f01 = -9.041191886754806e-1,
			 f02 =  4.169608567309818e-2,
			 f03 = -9.325971761333677e-3,
			 f04 =  4.699055851002199e-2,
			 f05 = -3.086923404061666e-2,
			 f06 =  1.057761186019000e-2,
			 f07 = -7.349302346007727e-2,
			 f08 =  1.444842576424337e-1,
			 f09 = -1.408425967872030e-1,
			 f10 =  1.070981398567760e-1,
			 f11 = -1.768451760854797e-2,
			 f12 = -4.013688314067293e-1,
			 f13 =  7.209753205388577e-1,
			 f14 = -1.807444462285120e-1,
			 f15 =  1.362305015808993e-1,
			 f16 = -9.500974920072897e-1,
			 f17 =  1.192134856624248,
			 f18 = -9.191161283559850e-2,
			 f19 = -1.008594411490973,
			 f20 =  8.020279271484482e-1,
			 f21 = -3.930534388853466e-1,
			 f22 = -2.026853316399942e-2,
			 f23 = -2.722731069001690e-2,
			 f24 =  5.032098120548072e-2,
			 f25 = -2.354888890484222e-2,
			 f26 = -2.454090179215001e-2,
			 f27 =  4.125987229048937e-2,
			 f28 = -3.533404753585094e-2,
			 f29 =  3.766063025852511e-2,
			 f30 = -3.358409746243470e-2,
			 f31 = -2.242158862056258e-2,
			 f32 =  2.102254738058931e-2,
			 f33 = -3.048635435546108e-2,
			 f34 = -1.996293091714222e-2,
			 f35 =  2.577703068234217e-2,
			 f36 = -1.292053030649309e-2,
			 g01 =  3.332286683867741e5,
			 g02 =  1.416532517833479e4,
			 g03 = -1.021129089258645e4,
			 g04 =  2.356370992641009e4,
			 g05 = -8.483432350173174e3,
			 g06 =  2.279927781684362e4,
			 g07 =  1.506238790315354e4,
			 g08 =  4.194030718568807e3,
			 g09 = -3.146939594885272e5,
			 g10 = -7.549939721380912e4,
			 g11 =  2.790535212869292e6,
			 g12 =  1.078851928118102e5,
			 g13 = -1.062493860205067e7,
			 g14 =  2.082909703458225e7,
			 g15 = -2.046810820868635e7,
			 g16 =  8.039606992745191e6,
			 g17 = -2.023984705844567e4,
			 g18 =  2.871769638352535e4,
			 g19 = -1.444841553038544e4,
			 g20 =  2.261532522236573e4,
			 g21 = -2.090579366221046e4,
			 g22 = -1.128417003723530e4,
			 g23 =  3.222965226084112e3,
			 g24 = -1.226388046175992e4,
			 g25 =  1.506847628109789e4,
			 g26 = -4.584670946447444e4,
			 g27 =  1.596119496322347e4,
			 g28 = -6.338852410446789e4,
			 g29 =  8.951570926106525e4,
			 saturation_fraction = 0.0;
	/**
	-----------------------------------------------------------------------
	   Finding func0.  This is the value of the function, func, that would
	   result in the output w_Ih_final being exactly zero.
	-----------------------------------------------------------------------
	*/
	func0 = h_pot_bulk - gtc.gsw_cp0
			  *gsw_ct_freezing(sa_bulk,p,saturation_fraction);
	/**
	-------------------------------------------------------------------------
	   Setting the three outputs for data points that have func0 non-negative
	-------------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		/**
		    When func0 is zero or positive then the final answer will contain
		    no frazil ice; that is, it will be pure seawater that is warmer
		    than the freezing temperature. If func0 >= 0 we do not need to go
		    through the modified Newton-Raphson procedure and we can simply
		    write down the answer, as in the following 4 lines of code.
		*/
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	    Begin finding the solution for data points that have func0 < 0, so that
	    the output will have a positive ice mass fraction w_Ih_final.
	    -----------------------------------------------------------------------
	*/
	/**Evalaute a polynomial for w_Ih in terms of SA_bulk, func0 and p*/
	x = sa_bulk*1e-2;
	y = func0/3e5;
	z = p*1e-4;
	w_ih = y*(f01
	         + x*(f02 + x*(f03 + x*(f04 + x*(f05 + f06*x))))
            + y*(f07 + x*(f08 + x*(f09 + x*(f10 + f11*x))) + y*(f12 + x*(f13
				+ x*(f14 + f15*x)) + y*(f16 + x*(f17 + f18*x) + y*(f19 + f20*x
				+ f21*y)))) + z*(f22 + x*(f23 + x*(f24 + f25*x)) + y*(x*(f26 + f27*x)
				+ y*(f28 + f29*x + f30*y)) + z*(f31 + x*(f32 + f33*x)
				+ y*(f34 + f35*x + f36*y))));

	if (w_ih > 0.9)
	{
		/**
		   The ice mass fraction out of this code is restricted to be
		   less than 0.9.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	    The initial guess at the absolute salinity of the interstitial seawater
	*/
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Doing a Newton step with a separate polynomial estimate of the mean
	      derivative dfunc_dw_Ih_mean_poly.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
	func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
	xa = sa*1e-2;
	dfunc_dw_ih_mean_poly = g01
	                        + xa*(g02 + xa*(g03 + xa*(g04 + g05*xa)))
									+ w_ih*(xa*(g06 + xa*(g07 + g08*xa)) + w_ih*(xa*(g09 + g10*xa)
									+ w_ih*xa*(g11 + g12*xa + w_ih*(g13 + w_ih*(g14 + w_ih*(g15
									+ g16*w_ih)))))) + z*(g17 + xa*(g18 + g19*xa) + w_ih*(g20
									+ w_ih*(g21 + g22*w_ih) + xa*(g23 + g24*xa*w_ih))
									+ z*(g25 + xa*(g26 + g27*xa) + w_ih*(g28 + g29*w_ih)));

	w_ih_old = w_ih;
	w_ih = w_ih_old - func/dfunc_dw_ih_mean_poly;
	sa = sa_bulk/(1.0 - w_ih);
	/**
	  -----------------------------------------------------------------------
	     Calculating the estimate of the derivative of func, dfunc_dw_Ih, to
	     be fed into Newton s Method.
	  -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
	gsw_pot_enthalpy_ice_freezing_first_derivatives(sa,p,&dpot_h_ihf_dsa, NULL);
	dfunc_dw_ih = gtc.gsw_cp0*ctf - h_pot_ihf - sa*
					  (gtc.gsw_cp0*ctf_sa + w_ih*dpot_h_ihf_dsa/(1.0 - w_ih));

	if (w_ih >= 0.0 && w_ih <= 0.20 && sa > 15.0
			&& sa < 60.0 && p <= 3000.0)
	{
		max_iterations = 1;
	}
	else if (w_ih >= 0.0 && w_ih <= 0.85 && sa > 0.0
				&& sa < 120.0 && p <= 3500.0)
	{
		max_iterations = 2;
	}
	else
	{
		max_iterations = 3;
	}
	for (iterations = 1; iterations <= max_iterations; iterations++)
	{
		if (iterations > 1)
		{
			/**On the first iteration ctf and h_pot_ihf are both known*/
			ctf = gsw_ct_freezing(sa,p,saturation_fraction);
			h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
		}
		/**This is the method, func, whose zero we seek ...*/
		func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */
	if (w_ih < 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
	}
	else
	{
		*sa_final = sa;
		*ct_final = gsw_ct_freezing(sa,p,saturation_fraction);
		*w_ih_final = w_ih;
	}
}

/***************************************************************************
% gsw_pot_enthalpy_ice_freezing_first_derivatives         first derivatives
%                        of the potential enthalpy of ice at the conditions
%                   where ice and seawater are in thermodynamic equilibrium
%==========================================================================
%
% USAGE:
%  gsw_pot_enthalpy_ice_freezing_first_derivatives(SA,p)
%
% DESCRIPTION:
%  Calculates the first derivatives of the potential enthalpy of ice at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_enthalpy_ice_freezing_SA = the derivative of the potential enthalpy
%                  of ice at freezing (ITS-90) with respect to Absolute
%                  salinity at fixed pressure  [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
%                  of ice at freezing (ITS-90) with respect to pressure
%                  (in Pa) at fixed Absolute Salinity              [ K/Pa ]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa,double p,
                                                               double *pot_enthalpy_ice_freezing_sa,
                                                               double *pot_enthalpy_ice_freezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  cp_ihf, pt_icef, ratio_temp, tf, tf_p, tf_sa, saturation_fraction = 0.0;

	tf = gsw_t_freezing(sa,p,saturation_fraction);
	pt_icef = gsw_pt0_from_t_ice(tf,p);
	ratio_temp = (gtc.gsw_t0 + pt_icef)/(gtc.gsw_t0 + tf);
	cp_ihf = gsw_cp_ice(tf,p);

	if ((pot_enthalpy_ice_freezing_sa != NULL) &&
			(pot_enthalpy_ice_freezing_p != NULL))
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa,&tf_p);
	}
	else if (pot_enthalpy_ice_freezing_sa != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa, NULL);
	}
	else if (pot_enthalpy_ice_freezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													NULL,&tf_p);
	}
	if (pot_enthalpy_ice_freezing_sa != NULL)
		*pot_enthalpy_ice_freezing_sa = ratio_temp*cp_ihf*tf_sa;

	if (pot_enthalpy_ice_freezing_p != NULL)
		*pot_enthalpy_ice_freezing_p = ratio_temp*cp_ihf*tf_p
												 - (gtc.gsw_t0 + pt_icef)*gsw_gibbs_ice(1,1,tf,p);
}

/***************************************************************************
% gsw_frazil_properties        Absolute Salinity, Conservative Temperature,
%                                       and the ice mass fraction from bulk
%                                       Absolute Salinity and bulk enthalpy
%==========================================================================
%
% USAGE:
%  gsw_frazil_properties(SA_bulk,h_bulk,p)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), w_Ih_final, which results from given values of the bulk
%  Absolute Salinity, SA_bulk, bulk enthalpy, h_bulk, occuring at pressure
%  p.  The final values of Absolute Salinity, SA_final, and Conservative
%  Temperature, CT_final, of the interstitial seawater phase are also
%  returned.  This code assumes that there is no dissolved air in the
%  seawater (that is, saturation_fraction is assumed to be zero
%  throughout the code).
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
%  present in the final state.  In this case the final state consists of
%  only seawater rather than being an equlibrium mixture of seawater and
%  ice which occurs when w_Ih_final is positive.  Note that when
%  w_Ih_final = 0, the final seawater is not at the freezing temperature.
%
%  Note that there is another GSW code,
%  gsw_frazil_properties_potential_poly(SA_bulk,h_pot_bulk,p) which
%  treats potential enthalpy as the conservative variable, while, in
%  contrast, the present code treats in situ enthalpy as the conservative
%  variable during the interaction of seawater and ice Ih.
%
% INPUT:
%  SA_bulk =  bulk Absolute Salinity of the seawater and ice mixture
%                                                                  [ g/kg ]
%  h_bulk  =  bulk enthalpy of the seawater and ice mixture        [ J/kg ]
%  p       =  sea pressure                                         [ dbar ]
%             ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state,
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at
%                 thermodynamic equilibrium.  If this ice mass fraction is
%                 zero there is no ice in the final state which consists
%                 only of seawater which is warmer than the freezing
%                 temperature.                                   [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Anonymous, 2014: Modelling the interaction between seawater and frazil
%   ice.  Manuscript, March 2015.  See Eqns. (8) - (15) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_frazil_properties(double sa_bulk,double h_bulk,double p,
                                       double *sa_final,double *ct_final,
                                       double *w_ih_final)
{
	int number_of_iterations;
	double cp_ih, ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
			 func, func0, hf, h_hat_ct, h_hat_sa,
			 h_ihf, sa, tf_sa, tf, w_ih_mean, w_ih_old, w_ih,
			 /**
			    Throughout this code seawater is taken to contain
			    no dissolved air.
			 */
			 saturation_fraction = 0.0,
			 num_f = 5.0e-2, num_f2 = 6.9e-7, num_p = 2.21;
	/**
	    ---------------
	      Finding func0
	    --------------
	*/
	ctf = gsw_ct_freezing(sa_bulk,p,saturation_fraction);
	func0 = h_bulk - gsw_enthalpy_ct_exact(sa_bulk,ctf,p);
	/**
	    -----------------------------------------------------------------------
	      When func0 is zero or positive we can immediately calculate the three
	      outputs, as the bulk enthalpy, h_bulk, is too large to allow any ice
	      at thermodynamic equilibrium. The result will be (warm) seawater with
	      no frazil ice being present. The three outputs can be set and the rest
	      of this code does not need to be performed.
	    -----------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = gsw_ct_from_enthalpy_exact(sa_bulk,h_bulk,p);
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	      Begin to find the solution for those data points that have func0 < 0,
	      implying that the output will be a positive ice mass fraction
	      w_Ih_final.
	      Do a quasi-Newton step with a separate polynomial estimate of the
	      derivative of func with respect to the ice mass fraction.  This
	      section of the code delivers initial values of both w_Ih and SA to
	      the rest of the more formal modified Newtons Method approach of
	      McDougall and Wotherspoon (2014).
	    -----------------------------------------------------------------------
	*/
	dfunc_dw_ih_mean_poly = 3.347814e+05
									- num_f*func0*(1.0 + num_f2*func0) - num_p*p;
	w_ih = gsw_min_d(-func0/dfunc_dw_ih_mean_poly, 0.95);
	sa = sa_bulk/(1.0 - w_ih);
	if (sa < 0.0 || sa > 120.0)
	{
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	      Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
	      fed into the iterative Newton s Method.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	hf = gsw_enthalpy_ct_exact(sa,ctf,p);
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	h_ihf = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p, &h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,NULL);
	dfunc_dw_ih = hf - h_ihf - sa*(h_hat_sa + h_hat_ct*ctf_sa +
											 w_ih*cp_ih*tf_sa/(1.0 - w_ih));
	/**
	    -----------------------------------------------------------------------
	     Enter the main McDougall-Wotherspoon (2014) modified Newton-Raphson
	     loop
	    -----------------------------------------------------------------------
	*/
	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		if (number_of_iterations > 1)
		{
			/** on the first iteration these values are already known */
			ctf = gsw_ct_freezing(sa,p,saturation_fraction);
			hf = gsw_enthalpy_ct_exact(sa,ctf,p);
			tf = gsw_t_freezing(sa,p,saturation_fraction);
			h_ihf = gsw_enthalpy_ice(tf,p);
		}
		func = h_bulk - (1.0 - w_ih)*hf - w_ih*h_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		w_ih_mean = 0.5*(w_ih + w_ih_old);
		if (w_ih_mean > 0.9)
		{
			/**This ensures that the mass fraction of ice never exceeds 0.9*/
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih_mean);
		ctf = gsw_ct_freezing(sa,p,saturation_fraction);
		hf = gsw_enthalpy_ct_exact(sa,ctf,p);
		tf = gsw_t_freezing(sa,p,saturation_fraction);
		h_ihf = gsw_enthalpy_ice(tf,p);
		cp_ih = gsw_cp_ice(tf,p);
		gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p, &h_hat_sa, &h_hat_ct);
		gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa, NULL);
		dfunc_dw_ih = hf - h_ihf - sa*(h_hat_sa + h_hat_ct*ctf_sa
												 + w_ih_mean*cp_ih*tf_sa/
												 (1.0 - w_ih_mean));
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			/**This ensures that the mass fraction of ice never exceeds 0.9*/
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */
	*sa_final = sa;
	*ct_final = gsw_ct_freezing(sa,p,saturation_fraction);
	*w_ih_final = w_ih;
	if (*w_ih_final < 0.0)
	{
		/**
		      This will only trap cases that are smaller than zero by just
		      machine precision
		*/
		*sa_final = sa_bulk;
		*ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_bulk,p);
		*w_ih_final = 0.0;
	}
}


/***************************************************************************
% gsw_cp_ice                                  isobaric heat capacity of ice
%==========================================================================
%
% USAGE:
%  gsw_cp_ice(t,p)
%
% DESCRIPTION:
%  Calculates the isobaric heat capacity of seawater.
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  cp_ice  =  heat capacity of ice                           [J kg^-1 K^-1]
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
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_cp_ice(double t, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (-(t + gtc.gsw_t0)*gsw_gibbs_ice(2,0,t,p));
}


/***************************************************************************
% gsw_frazil_ratios_adiabatic_poly      ratios of SA, CT and P changes when
%                           frazil ice forms due to the adiabatic change in
%                          pressure of a mixture of seawater and ice (poly)
%==========================================================================
%
% USAGE:
%  gsw_frazil_ratios_adiabatic_poly(SA,p,w_Ih)
%
% DESCRIPTION:
%  Calculates the ratios of SA, CT and P changes when frazil ice forms or
%  melts in response to an adiabatic change in pressure of a mixture of
%  seawater and frazil ice crystals.
%
%  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
%  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
%  whereas dCT/dSA would then be infinite.
%
%  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
%  derivatives with the pressure measured in Pa not dbar.
%
%  This function uses the computationally-efficient expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015) and the polynomial
%  expression for freezing temperature based on Conservative Temperature
%  (McDougall et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
%           sum of the masses of ice and seawater.  That is, the mass of
%           ice divided by the mass of the final mixed fluid.
%           w_Ih must be between 0 and 1.                      [ unitless ]
%
% OUTPUT:
%  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity
%                    to that of Conservative Temperature       [ g/(kg K) ]
%  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity
%                    to that of pressure (in Pa)              [ g/(kg Pa) ]
%  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
%                    to that of pressure (in Pa)                   [ K/Pa ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqns. (47), (48) and (49) of this manuscript.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_frazil_ratios_adiabatic_poly(double sa,double p,double w_ih,
                                                double *dsa_dct_frazil,
                                                double *dsa_dp_frazil,
                                                double *dct_dp_frazil)
{
	double bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
			 rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
			 ctf, ctf_sa, ctf_p, saturation_fraction = 0.0;

	tf = gsw_t_freezing_poly(sa,p,saturation_fraction);
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gamma_ih = gsw_adiabatic_lapse_rate_ice(tf,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction, &tf_sa,&tf_p);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction, &ctf_sa,&ctf_p);
	wcp = cp_ih*w_ih/(1.0 - w_ih);
	part = (tf_p - gamma_ih)/ctf_p;
	bracket1 = h_hat_ct + wcp*part;
	bracket2 = h - h_ih - sa*(h_hat_sa + wcp*(tf_sa - part*ctf_sa));
	rec_bracket3 = 1.0/(h - h_ih - sa*(h_hat_sa + h_hat_ct*ctf_sa + wcp*tf_sa));

	*dsa_dct_frazil = sa*(bracket1/bracket2);
	*dsa_dp_frazil = sa*ctf_p*bracket1*rec_bracket3;
	*dct_dp_frazil = ctf_p*bracket2*rec_bracket3;
}

/***************************************************************************
% gsw_frazil_ratios_adiabatic    ratios of SA, CT and P changes when frazil
%                                  ice forms due to the adiabatic change in
%                                 pressure of a mixture of seawater and ice
%==========================================================================
%
% USAGE:
%  gsw_frazil_ratios_adiabatic(SA,p,w_Ih)
%
% DESCRIPTION:
%  Calculates the ratios of SA, CT and P changes when frazil ice forms or
%  melts in response to an adiabatic change in pressure of a mixture of
%  seawater and frazil ice crystals.
%
%  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
%  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
%  whereas dCT/dSA would then be infinite.
%
%  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
%  derivatives with the pressure measured in Pa not dbar.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
%           sum of the masses of ice and seawater.  That is, the mass of
%           ice divided by the mass of the final mixed fluid.
%           w_Ih must be between 0 and 1.                      [ unitless ]
%
% OUTPUT:
%  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity
%                    to that of Conservative Temperature       [ g/(kg K) ]
%  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity
%                    to that of pressure (in Pa)              [ g/(kg Pa) ]
%  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
%                    to that of pressure (in Pa)                   [ K/Pa ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqns. (47), (48) and (49) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_frazil_ratios_adiabatic(double sa,double p,double w_ih,
                                             double *dsa_dct_frazil,
                                             double *dsa_dp_frazil,
                                             double *dct_dp_frazil)
{
	double bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
			 rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
			 ctf, ctf_sa, ctf_p, saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	h = gsw_enthalpy_ct_exact(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gamma_ih = gsw_adiabatic_lapse_rate_ice(tf,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,&tf_p);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction, &ctf_sa,&ctf_p);
	wcp = cp_ih*w_ih/(1.0 - w_ih);
	part = (tf_p - gamma_ih)/ctf_p;
	bracket1 = h_hat_ct + wcp*part;
	bracket2 = h - h_ih - sa*(h_hat_sa + wcp*(tf_sa - part*ctf_sa));
	rec_bracket3 = 1.0/(h - h_ih - sa*(h_hat_sa + h_hat_ct*ctf_sa + wcp*tf_sa));

	*dsa_dct_frazil = sa*(bracket1/bracket2);
	*dsa_dp_frazil = sa*ctf_p*bracket1*rec_bracket3;
	*dct_dp_frazil = ctf_p*bracket2*rec_bracket3;
}

/***************************************************************************
% gsw_adiabatic_lapse_rate_ice                         adiabatic lapse rate
%==========================================================================
%
% USAGE:
%  adiabatic_lapse_rate_ice = gsw_adiabatic_lapse_rate_ice(t,p)
%
% DESCRIPTION:
%  Calculates the adiabatic lapse rate of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  adiabatic_lapse_rate_ice  =  adiabatic lapse rate               [ K/Pa ]
%    Note.  The output is in unit of degress Celsius per Pa,
%      (or equivilently K/Pa) not in units of K/dbar.
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_adiabatic_lapse_rate_ice(double t,double p)
{
	return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(2,0,t,p));
}


/***************************************************************************
% gsw_alpha_wrt_t_ice                  thermal expansion coefficient of ice
%                                       with respect to in-situ temperature
%==========================================================================
%
% USAGE:
%  gsw_alpha_wrt_t_ice(t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of ice with respect to
%  in-situ temperature.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  alpha_wrt_t_ice  =  thermal expansion coefficient of ice with respect
%                      to in-situ temperature                       [ 1/K ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.18.1) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_alpha_wrt_t_ice(double t, double p)
{
	return (gsw_gibbs_ice(1,1,t,p) / gsw_gibbs_ice(0,1,t,p));
}
/***************************************************************************
% gsw_chem_potential_water_ice           chemical potential of water in ice
%==========================================================================
%
% USAGE:
%  gsw_chem_potential_water_ice(t,p)
%
% DESCRIPTION:
%  Calculates the chemical potential of water in ice from in-situ
%  temperature and pressure.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  chem_potential_water_ice  =  chemical potential of ice          [ J/kg ]
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
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_chem_potential_water_ice(double t,double p)
{
	return (gsw_gibbs_ice(0,0,t,p));
}


/***************************************************************************
% gsw_enthalpy_diff            difference of enthalpy between two pressures
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_enthalpy_diff(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This
%  function uses the computationally-efficient expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).  The output
%  (enthalpy_diff) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA         =  Absolute Salinity                                 [ g/kg ]
%  CT         =  Conservative Temperature (ITS-90)                [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar )
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  enthalpy_diff  =  difference of specific enthalpy            [ J/kg ]
%                       (deep minus shallow)
%
% AUTHOR:
%  Trevor McDougall & Paul Barker.                     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual.
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_enthalpy_diff(double sa,double ct,
                                    double p_shallow,
                                    double p_deep)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double  dynamic_enthalpy_shallow, dynamic_enthalpy_deep,
			  part_1, part_2, part_3, part_4, part_5, xs, ys,
			  z_deep, z_shallow;

	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z_shallow = p_shallow*1e-4;
	z_deep = p_deep*1e-4;

	part_1 = gsvco.h001
				+ xs*(gsvco.h101 + xs*(gsvco.h201 + xs*(gsvco.h301 + xs*(gsvco.h401
				+ xs*(gsvco.h501 + gsvco.h601*xs))))) + ys*(gsvco.h011
				+ xs*(gsvco.h111 + xs*(gsvco.h211+ xs*(gsvco.h311
				+ xs*(gsvco.h411 + gsvco.h511*xs)))) + ys*(gsvco.h021 + xs*(gsvco.h121
				+ xs*(gsvco.h221 + xs*(gsvco.h321 + gsvco.h421*xs))) + ys*(gsvco.h031
				+ xs*(gsvco.h131 + xs*(gsvco.h231 + gsvco.h331*xs)) + ys*(gsvco.h041
				+ xs*(gsvco.h141 + gsvco.h241*xs) + ys*(gsvco.h051 + gsvco.h151*xs
				+ gsvco.h061*ys)))));

	part_2 = gsvco.h002
				   + xs*( gsvco.h102 + xs*( gsvco.h202 + xs*( gsvco.h302
					+ xs*( gsvco.h402 + gsvco.h502*xs))))
				   + ys*( gsvco.h012 + xs*( gsvco.h112 + xs*( gsvco.h212
					+ xs*( gsvco.h312 + gsvco.h412*xs))) + ys*( gsvco.h022
					+ xs*( gsvco.h122 + xs*( gsvco.h222 + gsvco.h322*xs))
					+ ys*( gsvco.h032 + xs*( gsvco.h132 + gsvco.h232*xs)
					+ ys*( gsvco.h042 + gsvco.h142*xs + gsvco.h052*ys))));

	part_3 = gsvco.h003
				   + xs*(gsvco.h103 + xs*(gsvco.h203 + xs*(gsvco.h303
					+ gsvco.h403*xs))) + ys*(gsvco.h013
					+ xs*(gsvco.h113 + xs*(gsvco.h213 + gsvco.h313*xs))
					+ ys*(gsvco.h023 + xs*(gsvco.h123 + gsvco.h223*xs)
					+ ys*(gsvco.h033 + gsvco.h133*xs + gsvco.h043*ys)));

	part_4 = gsvco.h004
	            + xs*(gsvco.h104 + gsvco.h204*xs)
				   + ys*(gsvco.h014 + gsvco.h114*xs + gsvco.h024*ys);

	part_5 = gsvco.h005 + gsvco.h105*xs + gsvco.h015*ys;

	dynamic_enthalpy_shallow = z_shallow*(part_1
	                              + z_shallow*(part_2 + z_shallow*(part_3
	                              + z_shallow*(part_4 + z_shallow*(part_5
											+ z_shallow*(gsvco.h006 + gsvco.h007*z_shallow))))));

	dynamic_enthalpy_deep = z_deep*(part_1
	                           + z_deep*(part_2 + z_deep*(part_3
										+ z_deep*(part_4 + z_deep*(part_5
										+ z_deep*(gsvco.h006 + gsvco.h007 * z_deep))))));

	return ((dynamic_enthalpy_deep - dynamic_enthalpy_shallow)*gtc.db2pa*1e4);
}

/***************************************************************************
% gsw_Helmholtz_energy_ice                          Helmholtz energy of ice
%==========================================================================
%
% USAGE:
%  gsw_Helmholtz_energy_ice(t,p)
%
% DESCRIPTION:
%  Calculates the Helmholtz energy of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  Helmholtz_energy_ice  =  Helmholtz energy of ice                [ J/kg ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_helmholtz_energy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gsw_gibbs_ice(0,0,t,p) - (gtc.db2pa*p
	         + gtc.gsw_p0) * gsw_gibbs_ice(0,1,t,p));
}

/***************************************************************************
% gsw_ice_fraction_to_freeze_seawater         ice mass fraction, which when
%                                  melted into seawater, brings the diluted
%                                      seawater to the freezing temperature
%==========================================================================
%
% USAGE:
%  gsw_ice_fraction_to_freeze_seawater(SA,CT,p,t_Ih)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), which, when melted into seawater having (SA,CT,p) causes
%  the final dilute seawater to be at the freezing temperature.  The other
%  outputs are the Absolute Salinity and Conservative Temperature of the
%  final diluted seawater.
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  t_Ih =  in-situ temperature of the ice at pressure p (ITS-90)  [ deg C ]
%
% OUTPUT:
%  SA_freeze = Absolute Salinity of seawater after the mass fraction of
%              ice, ice_fraction, at temperature t_Ih has melted into the
%              original seawater, and the final mixture is at the freezing
%              temperature of seawater.                            [ g/kg ]
%
%  CT_freeze = Conservative Temperature of seawater after the mass
%              fraction, w_Ih, of ice at temperature t_Ih has melted into
%              the original seawater, and the final mixture is at the
%              freezing temperature of seawater.                  [ deg C ]
%
%  w_Ih      = mass fraction of ice, having in-situ temperature t_Ih,
%              which, when melted into seawater at (SA,CT,p) leads to the
%              final diluted seawater being at the freezing temperature.
%              This output must be between 0 and 1.              [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (9) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_ice_fraction_to_freeze_seawater(double sa,double ct,
                                                   double p,double t_ih,
                                                   double *sa_freeze,
                                                   double *ct_freeze,
                                                   double *w_ih)
{
	int no_iter;
	double ctf, ctf_mean, ctf_old, ctf_plus1, ctf_zero,
			 dfunc_dsaf, func, func_plus1, func_zero, h, h_ih,
			 saf, saf_mean, saf_old, tf, h_hat_sa, h_hat_ct, ctf_sa;

	double  sa0 = 0.0, saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);

	if (ct < ctf)
	{
		/**The seawater ct input is below the freezing temp*/
		*sa_freeze = cppGSW_INVALID_VALUE;
		*ct_freeze = *sa_freeze;
		*w_ih = *sa_freeze;
		return;
	}

	tf = gsw_t_freezing(sa0,p,saturation_fraction);

	if (t_ih > tf)
	{
		/**The input, t_Ih, exceeds the freezing temperature at sa = 0*/
		*sa_freeze = cppGSW_INVALID_VALUE;
		*ct_freeze = *sa_freeze;
		*w_ih = *sa_freeze;
		return;
	}

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);
	ctf_zero = gsw_ct_freezing(sa0,p,saturation_fraction);
	func_zero = sa*(gsw_enthalpy_ct_exact(sa0,ctf_zero,p) - h_ih);
	ctf_plus1 = gsw_ct_freezing(sa+1.0,p,saturation_fraction);
	func_plus1 = sa*(gsw_enthalpy_ct_exact(sa+1.0,ctf_plus1,p) - h) - (h - h_ih);
	saf = -(sa+1.0)*func_zero/(func_plus1 - func_zero);   /**initial guess*/
	ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	gsw_enthalpy_first_derivatives_ct_exact(saf,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(saf,p,1.0,&ctf_sa,NULL);
	dfunc_dsaf = sa*(h_hat_sa + h_hat_ct*ctf_sa) - (h - h_ih);

	for (no_iter = 1; no_iter <= 2; no_iter++)
	{
		saf_old = saf;
		ctf_old = ctf;
		func = sa*(gsw_enthalpy_ct_exact(saf_old,ctf_old,p) - h)
				 - (saf_old - sa)*(h - h_ih);

		saf = saf_old - func/dfunc_dsaf;
		saf_mean = 0.5*(saf + saf_old);
		ctf_mean = gsw_ct_freezing(saf_mean,p,saturation_fraction);
		gsw_enthalpy_first_derivatives_ct_exact(saf_mean,ctf_mean,p,
															 &h_hat_sa, &h_hat_ct);

		gsw_ct_freezing_first_derivatives(saf_mean,p,saturation_fraction,
													 &ctf_sa, NULL);

		dfunc_dsaf = sa*(h_hat_sa + h_hat_ct*ctf_sa) - (h - h_ih);
		saf = saf_old - func/dfunc_dsaf;
		ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	}
	/**
	    After these 2 iterations of this modified Newton-Raphson method, the
	    error in SA_freeze is less than 1.3d0x10^-13 g/kg, in CT_freeze is
	    less than   4x10^-13 deg C and in w_Ih is less than 3.8d0x10^-15
	    which represent machine precision for these calculations.
	*/
	*sa_freeze = saf;
	*ct_freeze = ctf;
	*w_ih = (h - gsw_enthalpy_ct_exact(*sa_freeze,*ct_freeze,p))/(h - h_ih);
}


/***************************************************************************
% gsw_internal_energy_ice                   specific internal energy of ice
%==========================================================================
%
% USAGE:
%  gsw_internal_energy_ice(t,p)
%
% DESCRIPTION:
%  Calculates the specific internal energy of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  internal_energy_ice  =  specific internal energy (u)              [J/kg]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_internal_energy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gsw_gibbs_ice(0,0,t,p) - (gtc.gsw_t0
	         + t) * gsw_gibbs_ice(1,0,t,p) - (gtc.db2pa*p
	         + gtc.gsw_p0)*gsw_gibbs_ice(0,1,t,p));
}


/***************************************************************************
% gsw_kappa_const_t_ice                          isothermal compressibility
%==========================================================================
%
% USAGE:
%  gsw_kappa_const_t_ice(t,p)
%
% DESCRIPTION:
%  Calculates isothermal compressibility of ice.
%  Note. This is the compressibility of ice AT CONSTANT IN-SITU
%    TEMPERATURE
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  kappa_const_t_ice  =  isothermal compressibility                [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_kappa_const_t_ice(double t,double p)
{
	return (-gsw_gibbs_ice(0,2,t,p) / gsw_gibbs_ice(0,1,t,p));
}


/***************************************************************************
% gsw_kappa_ice                                  isentropic compressibility
%==========================================================================
%
% USAGE:
%  gsw_kappa_ice(t,p)
%
% DESCRIPTION:
%  Calculates the isentropic compressibility of ice.
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  kappa_ice  =  isentropic compressibility                        [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_kappa_ice(double t,double p)
{
	double gi_tp, gi_tt;

	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);

	return ((gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p)) /
			  (gsw_gibbs_ice(0,1,t,p)*gi_tt));
}


/***************************************************************************
% gsw_melting_ice_equilibrium_SA_CT_ratio    ratio of SA to CT changes when
%                             ice melts into a large mass of seawater, with
%                              both the seawater and ice temperatures being
%                      almost equal to the equilibrium freezing temperature
%==========================================================================
%
% USAGE:
%  gsw_melting_ice_equilibrium_SA_CT_ratio(SA,p)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when ice melts into seawater
%  with both the seawater and the seaice temperatures being almost equal to
%  the equilibrium freezing temperature.  It is assumed that a small mass
%  of ice melts into an infinite mass of seawater.  If indeed the
%  temperature of the seawater and the ice were both equal to the freezing
%  temperature, then no melting or freezing would occur; an imbalance
%  between these three temperatures is needed for freezing or melting to
%  occur (the three temperatures being (1) the seawater temperature,
%  (2) the ice temperature, and (3) the freezing temperature.
%
%  The output, melting_ice_equilibrium_SA_CT_ratio, is dSA/dCT rather than
%  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is zero
%  whereas dCT/dSA would be infinite.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  melting_ice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
%                                changes when ice melts into seawater, with
%                                the seawater and seaice being close to the
%                                freezing temperature.         [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (16) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_ice_equilibrium_sa_ct_ratio(double sa,double p)
{
	double ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing(sa,p,saturation_fraction);
	h = gsw_enthalpy_ct_exact(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);

	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	/**note that h_hat_ct is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa * h_hat_ct / (h - h_ih - sa * h_hat_sa));
}


/***************************************************************************
% gsw_melting_ice_SA_CT_ratio                ratio of SA to CT changes when
%                                                   ice melts into seawater
%==========================================================================
%
% USAGE:
%  gsw_melting_ice_SA_CT_ratio(SA,CT,p,t_Ih)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when ice melts into seawater.
%  It is assumed that a small mass of ice melts into an infinite mass of
%  seawater.  Because of the infinite mass of seawater, the ice will always
%  melt.
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
%  This is done so that when SA = 0, the output, dSA/dCT is zero whereas
%  dCT/dSA would be infinite.
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure at which the melting occurs                [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
%
% OUTPUT:
%  melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
%                            into a large mass of seawater
%                                                          [ g kg^-1 K^-1 ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (13) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_ice_sa_ct_ratio(double sa,double ct,double p,
                                             double t_ih)
{
	double  ctf, h, h_ih, tf, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temperature*/
		return (cppGSW_INVALID_VALUE);
	}
	tf = gsw_t_freezing(0.0,p,saturation_fraction);
	if (t_ih > tf)
	{
		/**t_ih exceeds the freezing temperature at sa = 0*/
		return (cppGSW_INVALID_VALUE);
	}
	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);

	gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,&h_hat_sa,&h_hat_ct);
	/**Note that h_hat_CT is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa * h_hat_ct / (h - h_ih - sa * h_hat_sa));
}



/***************************************************************************
% gsw_melting_seaice_equilibrium_SA_CT_ratio      ratio of SA to CT changes
%                         when sea ice melts into a large mass of seawater,
%                     with both the seawater and sea ice temperatures being
%                      almost equal to the equilibrium freezing temperature
%==========================================================================
%
% USAGE:
%  gsw_melting_seaice_equilibrium_SA_CT_ratio(SA,p)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into
%  seawater with both the seawater and the sea ice temperatures being
%  almost equal to the equilibrium freezing temperature.  It is assumed
%  that a small mass of seaice melts into an infinite mass of seawater.  If
%  indeed the temperature of the seawater and the sea ice were both equal
%  to the freezing temperature, then no melting or freezing would occur; an
%  imbalance between these three temperatures is needed for freezing or
%  melting to occur (the three temperatures being (1) the seawater
%  temperature, (2) the sea ice temperature, and (3) the freezing
%  temperature.
%
%  Note that the output of this function, dSA/dCT is independent of the
%  sea ice salinity, SA_seaice.  That is, the output applies equally to
%  pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result
%  is proven in McDougall et al. (2014).
%
%  The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather
%  than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is
%  zero whereas dCT/dSA would be infinite.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
%                            changes when sea ice melts into seawater, with
%                            the seawater and sea ice being close to the
%                            freezing temperature.             [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (29) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa,
                                                            double p)
{
	double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing(sa,p,saturation_fraction);
	h = gsw_enthalpy_ct_exact(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);

	return (sa * h_hat_ct / (h - h_ih - sa * h_hat_sa));
}


/****************************************************************************
% gsw_melting_seaice_equilibrium_SA_CT_ratio_poly
%   ratio of SA to CT changes when sea ice melts into a large mass of
%   seawater, with both the seawater and sea ice temperatures being
%   almost equal to the equilibrium freezing temperature  (poly)
%==========================================================================
%
% USAGE:
%  gsw_melting_seaice_equilibrium_SA_CT_ratio(SA,p)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into
%  seawater with both the seawater and the sea ice temperatures being
%  almost equal to the equilibrium freezing temperature.  It is assumed
%  that a small mass of seaice melts into an infinite mass of seawater.  If
%  indeed the temperature of the seawater and the sea ice were both equal
%  to the freezing temperature, then no melting or freezing would occur; an
%  imbalance between these three temperatures is needed for freezing or
%  melting to occur (the three temperatures being (1) the seawater
%  temperature, (2) the sea ice temperature, and (3) the freezing
%  temperature.
%
%  Note that the output of this function, dSA/dCT is independent of the
%  sea ice salinity, SA_seaice.  That is, the output applies equally to
%  pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result
%  is proven in McDougall et al. (2014).
%
%  The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather
%  than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is
%  zero whereas dCT/dSA would be infinite.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
% OUTPUT:
%  melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
%                            changes when sea ice melts into seawater, with
%                            the seawater and sea ice being close to the
%                            freezing temperature.             [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (29) of this manuscript.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa,
                                                                  double p)
{
	double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing_poly(sa,p,saturation_fraction);
	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);

	return (sa * h_hat_ct / (h - h_ih - sa  *h_hat_sa));
}



/****************************************************************************
% gsw_melting_seaice_into_seawater             the resulting SA and CT when
%                                           sea ice is melted into seawater
%==========================================================================
%
% USAGE:
%   gsw_melting_seaice_into_seawater(SA,CT,p,w_seaice,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity and Conservative Temperature that
%  results when a given mass of sea ice (or ice) melts and is mixed into a
%  known mass of seawater (whose properties are (SA,CT,p)).
%
%  If the ice contains no salt (e.g. if it is of glacial origin), then the
%  input 'SA_seaice' should be set to zero.
%
%  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
%  and 12 g/kg of salt (defined as the mass of salt divided by the mass of
%  ice Ih plus brine) and this programme returns NaN's if the input
%  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
%  usually this would imply that the pressure p should be zero, as sea ice
%  only occurs near the sea surface.  The code does not impose that p = 0
%  if SA_seaice is non-zero.  Rather, this is left to the user.
%
%  The Absolute Salinity, SA_brine, of the brine trapped in little pockets
%  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
%  surrounds these pockets.  As the sea ice temperature, t_seaice, may be
%  less than the freezing temperature, SA_brine is usually greater than the
%  Absolute Salinity of the seawater at the time and place when and where
%  the sea ice was formed.  So usually SA_brine will be larger than SA.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  w_seaice  =  mass fraction of sea ice, that is the mass of sea ice
%               divided by the sum of the masses of sea ice and seawater.
%               That is, the mass of sea ice divided by the mass of the
%               final mixed fluid.  w_seaice must be between 0 and 1.
%                                                              [ unitless ]
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of
%               salt in sea ice, expressed in g of salt per kg of sea ice.
%                                                                  [ g/kg ]
%  t_seaice  =  the in-situ temperature of the sea ice (or ice) (ITS-90)
%                                                                 [ deg C ]
%
% OUTPUT:
%  SA_final  =  Absolute Salinity of the mixture of the melted sea ice
%               (or ice) and the orignal seawater                  [ g/kg ]
%  CT_final  =  Conservative Temperature of the mixture of the melted
%               sea ice (or ice) and the orignal seawater         [ deg C ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%     Eqns. (8) and (9) are the simplifications when SA_seaice = 0.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_melting_seaice_into_seawater(double sa,double ct,double p,
                                                double w_seaice,double sa_seaice,
                                                double t_seaice,double *sa_final,
                                                double *ct_final)
{
	double  ctf, h, h_brine, h_final, h_ih, sa_brine, tf_sa_seaice;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**The seawater ct input is below the freezing temp*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		return;
	}

	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;

	if (t_seaice > tf_sa_seaice)
	{
		/**
		  The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
		  some ice Ih in the sea ice. Without this buffer, that is if t_seaice
		  is allowed to be exactly equal to tf_sa_seaice, the seaice is
		  actually 100% brine at Absolute Salinity of SA_seaice.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		return;
	}

	sa_brine = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);

	if (sa_brine >= cppGSW_ERROR_LIMIT)
	{
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		return;
	}

	h_brine = gsw_enthalpy_t_exact(sa_brine,t_seaice,p);
	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	h_final = h - w_seaice*(h - h_ih - (h_brine - h_ih)*sa_seaice/sa_brine);
	*sa_final = sa - w_seaice*(sa - sa_seaice);
	/**
	ctf = gsw_ct_freezing(sa_final,p,saturation_fraction)
	if (h_final .lt. gsw_enthalpy_ct_exact(sa_final,ctf,p)) then
	       Melting this much seaice is not possible as it would result in
	       frozen seawater
	     sa_final = gsw_error_code(4,func_name)
	     ct_final = sa_final
	     return
	end if
	*/
	*ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_final,p);

	if (*ct_final > cppGSW_ERROR_LIMIT)
	{
		*sa_final = *ct_final;
	}
}

/***************************************************************************
% gsw_melting_seaice_SA_CT_ratio_poly        ratio of SA to CT changes when
%                                        sea ice melts into seawater (poly)
%==========================================================================
%
% USAGE:
%  gsw_melting_seaice_SA_CT_ratio_poly(SA,CT,p,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into
%  seawater.  It is assumed that a small mass of sea ice melts into an
%  infinite mass of seawater.  Because of the infinite mass of seawater,
%  the sea ice will always melt.
%
%  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
%  and 12 g/kg of salt (defined as the mass of salt divided by the mass of
%  ice Ih plus brine) and this programme returns NaN's if the input
%  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
%  usually this would imply that the pressure p should be zero, as sea ice
%  only occurs near the sea surface.  The code does not impose that p = 0
%  if SA_seaice is non-zero.  Rather, this is left to the user.
%
%  The Absolute Salinity, SA_brine, of the brine trapped in little pockets
%  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
%  surrounds these pockets.  As the seaice temperature, t_seaice, may be
%  less than the freezing temperature, SA_brine is usually greater than the
%  Absolute Salinity of the seawater at the time and place when and where
%  the sea ice was formed.  So usually SA_brine will be larger than SA.
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
%  This is done so that when (SA - seaice_SA) = 0, the output, dSA/dCT is
%  zero whereas dCT/dSA would be infinite.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction
%               of salt in sea ice expressed in g of salt per kg of
%               sea ice                                            [ g/kg ]
%  t_seaice =   the in-situ temperature of the sea ice (ITS-90)   [ deg C ]
%
% OUTPUT:
%  melting_seaice_SA_CT_ratio = the ratio dSA/dCT of SA to CT changes when
%                sea ice melts into a large mass of seawater   [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (31) of this manuscript.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_seaice_sa_ct_ratio_poly(double sa,double ct,
                                                      double p,double sa_seaice,
                                                      double t_seaice)
{
	double  ctf, delsa, h, h_brine, h_ih, sa_brine, ct_brine,
			  tf_sa_seaice, h_hat_sa, h_hat_ct, saturation_fraction = 0.0;

	if (sa_seaice < 0.0 || sa_seaice > 15.0)
	{
		return (NAN);
	}
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temp*/
		return (cppGSW_INVALID_VALUE);
	}
	tf_sa_seaice = gsw_t_freezing_poly(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)
	{
		/**t_seaice exceeds the freezing sa*/
		return (cppGSW_INVALID_VALUE);
	}
	/**
	-----------------------------------------------------------------------
	The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
	some ice Ih in the sea ice.  Without this buffer, that is if t_seaice
	is allowed to be exactly equal to tf_sa_seaice, the sea ice is actually
	100% brine at Absolute Salinity of SA_seaice.
	-----------------------------------------------------------------------
	*/
	h = gsw_enthalpy(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives(sa,ct,p,&h_hat_sa,&h_hat_ct);
	sa_brine = gsw_sa_freezing_from_t_poly(t_seaice,p,saturation_fraction);
	if (sa_brine > cppGSW_ERROR_LIMIT)
	{
		return (cppGSW_INVALID_VALUE);
	}
	ct_brine = gsw_ct_from_t(sa_brine,t_seaice,p);
	h_brine = gsw_enthalpy(sa_brine,ct_brine,p);
	delsa = sa - sa_seaice;

	return (h_hat_ct * delsa / (h - h_ih
	         - delsa * h_hat_sa
	         - sa_seaice * (h_brine
	         - h_ih) / sa_brine));
}

/***************************************************************************
% gsw_pressure_freezing_CT             pressure of seawater at the freezing
%                                                               temperature
%==========================================================================
%
% USAGE:
%  gsw_pressure_freezing_CT(SA,CT,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the pressure (in dbar) of seawater at the freezing
%  temperature.  That is, the output is the pressure at which seawater,
%  with Absolute Salinity SA, Conservative Temperature CT, and with
%  saturation_fraction of dissolved air, freezes.  If the input values are
%  such that there is no value of pressure in the range between 0 dbar and
%  10,000 dbar for which seawater is at the freezing temperature, the
%  output, pressure_freezing, is put equal to NaN.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  pressure_freezing = sea pressure at which the seawater freezes  [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
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
%    See section 3.33 of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pressure_freezing_ct(double sa,double ct,
                                          double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int     i_iter, number_of_iterations = 3;

	double  ct_freezing_p0, ct_freezing_p10000, dctf_dp, f,
			  pf, pfm, pf_old, ctfreezing_p;

	/**
	  rec_Pa2dbar is to have dCTf_dp in units of K/dbar rather than K/Pa
	  Find CT > CT_freezing_p0.  If this is the case, the input CT value
	  represent seawater that will not be frozen at any positive p.
	*/
	ct_freezing_p0 = gsw_ct_freezing(sa,0.0,saturation_fraction);
	if (ct > ct_freezing_p0)
	{
		return (cppGSW_INVALID_VALUE);
	}
	/**
	  Find CT < CT_freezing_p10000.  If this is the case, the input CT value
	  represent seawater that is frozen even at p = 10,000 dbar.
	*/
	ct_freezing_p10000 = gsw_ct_freezing(sa,1e4,saturation_fraction);
	if (ct < ct_freezing_p10000)
	{
		return (cppGSW_INVALID_VALUE);
	}
	/**
	  This is the initial (linear) guess of the freezing pressure, in dbar.
	*/
	pf = gtc.rec_pa2db*(ct_freezing_p0 - ct)/
		  (ct_freezing_p0 - ct_freezing_p10000);
	gsw_ct_freezing_first_derivatives(sa,pf,saturation_fraction,
												 NULL,&ctfreezing_p);
	dctf_dp = gtc.rec_pa2db*ctfreezing_p;
	/**
	   this dctf_dp is the initial value of the partial derivative of
	   ct_freezing with respect to pressure (in dbar) at fixed sa,
	   assuming that the saturation_fraction is zero.
	*/
	for (i_iter = 0; i_iter < number_of_iterations; i_iter++)
	{
		pf_old = pf;
		f = gsw_ct_freezing(sa,pf_old,saturation_fraction) - ct;
		pf = pf_old - f/dctf_dp;
		pfm = 0.5*(pf + pf_old);
		gsw_ct_freezing_first_derivatives(sa,pfm,saturation_fraction,
													 NULL,&ctfreezing_p);
		dctf_dp = gtc.rec_pa2db*ctfreezing_p;
		pf = pf_old - f/dctf_dp;
	}

	if (gsw_sa_p_inrange(sa,pf))
	{
	   return (pf);
   }

	return (cppGSW_INVALID_VALUE);
}

/**************************************************************************
==========================================================================
  method: gsw_pt_from_pot_enthalpy_ice_poly_dh (pot_enthalpy_ice)
==========================================================================
   This method supplements 'gsw_pt_from_pot_enthalpy_ice_poly'
-------------------------------------------------------------------------
*****************************************************************************/
double TeosIce::gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice)
{
	double q1 = 2.594351081876611e-3,
			 p2 = 3.530155620427630e-8,
			 p3 = 2.330421169287162e-13,
			 p4 = 8.139369017110120e-19,
			 p5 = 1.610007265856420e-24,
			 p6 = 1.707103685781641e-30,
			 p7 = 7.658041152250651e-37;

	return (q1
	         + pot_enthalpy_ice*(p2
	         + pot_enthalpy_ice*(p3
				+ pot_enthalpy_ice*(p4
				+ pot_enthalpy_ice*(p5
				+ pot_enthalpy_ice*(p6
            + pot_enthalpy_ice*p7))))));
}


/***************************************************************************
% gsw_pt_from_pot_enthalpy_ice_poly        potential temperature refered to
%                         the surface from potential enthalpy of ice (poly)
%==========================================================================
%
% USAGE:
%  gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates the potential temperature of ice (whose reference sea
%  pressure is zero dbar) from the potential enthalpy of ice.  This is a
%  compuationally efficient polynomial fit to the potential enthalpy of
%  ice.
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% OUTPUT:
%  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
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
double TeosIce::gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice)
{
	double q0 = 2.533588268773218e2,
			 q1 = 2.594351081876611e-3,
			 q2 = 1.765077810213815e-8,
			 q3 = 7.768070564290540e-14,
			 q4 = 2.034842254277530e-19,
			 q5 = 3.220014531712841e-25,
			 q6 = 2.845172809636068e-31,
			 q7 = 1.094005878892950e-37;
	/**
	   The error of this fit ranges between -5e-5 and 2e-4 deg C over the potential
	   temperature range of -100 to 2 deg C, or the potential enthalpy range of
	   -5.7 x 10^5 to -3.3 x 10^5 J/kg.
	*/
	return (q0
			   + pot_enthalpy_ice*(q1
			   + pot_enthalpy_ice*(q2
			   + pot_enthalpy_ice*(q3
            + pot_enthalpy_ice*(q4
            + pot_enthalpy_ice*(q5
            + pot_enthalpy_ice*(q6
            + pot_enthalpy_ice*q7)))))));
}




/***************************************************************************
% gsw_pt_from_t_ice                            potential temperature of ice
% =========================================================================
%
% USAGE:
%  gsw_pt_from_t_ice(t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential temperature of ice Ih with the general reference
%  pressure, p_ref, from in-situ temperature, t.
%
%  A faster gsw routine exists if p_ref is indeed zero dbar.  This routine
%  is "gsw_pt0_from_t_ice(t,p)".
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  p_ref  =  reference pressure                                    [ dbar ]
%  (If reference pressure is not given then it is assumed that reference
%   pressure is zero).
%
%
% OUTPUT:
%  pt_ice  =  potential temperature of ice Ih with reference pressure,
%             p_ref, on the ITS-90 temperature scale              [ deg C ]
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
%    See appendix I of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pt_from_t_ice(double t,double p,double p_ref)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int number_of_iterations;
	double dentropy, dentropy_dt, dp,
			 pt_ice, pt_ice_old, ptm_ice, true_entropy,
			 p1 = -2.259745637898635e-4,
			 p2 =  1.486236778150360e-9,
			 p3 =  6.257869607978536e-12,
			 p4 = -5.253795281359302e-7,
			 p5 =  6.752596995671330e-9,
			 p6 =  2.082992190070936e-11,
			 q1 = -5.849191185294459e-15,
			 q2 =  9.330347971181604e-11,
			 q3 =  3.415888886921213e-13,
			 q4 =  1.064901553161811e-12,
			 q5 = -1.454060359158787e-10,
			 q6 = -5.323461372791532e-13;

	/**This is the starting polynomial for pt of ice Ih.*/
	dp = p - p_ref;
	pt_ice = t + dp*(p1 + (p + p_ref)*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));
	if (pt_ice < -gtc.gsw_t0) pt_ice = -gtc.gsw_t0;
	if (pt_ice < -273.0) pt_ice = pt_ice + 0.05;
	/**
	  we add 0.05 to the initial estimate of pt_ice at temps less than
	  -273 to ensure that it is never less than -273.15.
	*/
	dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice,p_ref);
	true_entropy = -gsw_gibbs_ice_part_t(t,p);

	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		pt_ice_old = pt_ice;
		dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref) - true_entropy;
		pt_ice = pt_ice_old - dentropy/dentropy_dt;
		ptm_ice = 0.5*(pt_ice + pt_ice_old);
		dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
		pt_ice = pt_ice_old - dentropy/dentropy_dt;
	}

	if (pt_ice < -273.0)
	{
		pt_ice = t + (p - p_ref)*(q1 + (p + p_ref)*(q2 + q3*t)
										  + t*(q4 + t*(q5 + q6*t)));
		dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice+0.01,p_ref);
		/**
		   we add 0.01 to the initial estimate of pt_ice used in the derivative
		   to ensure that it is never less than -273.15 because the derivative
		   approaches zero at absolute zero.
		*/
		for (number_of_iterations = 1; number_of_iterations <= 3;
				number_of_iterations++)
		{
			pt_ice_old = pt_ice;
			dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref)
						  - true_entropy;
			pt_ice = pt_ice_old - dentropy/dentropy_dt;
			ptm_ice = 0.5*(pt_ice + pt_ice_old);
			ptm_ice = ptm_ice + 0.01;
			dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
			pt_ice = pt_ice_old - dentropy/dentropy_dt;
		}
	}
	/**
	  For temperatures less than -273.1 degsC the maximum error is less than
	  2x10^-7 degsC. For temperatures between -273.1 and 273 the maximum error
	  is less than 8x10^-8 degsC, and for temperatures greater than -273 degsC the
	  maximum error is 1.5x10^-12 degsC.  These errors are over the whole
	  ocean depths with both p and pref varying independently between 0 and
	  10,000 dbar, while the in-situ temperature varied independently between
	  -273.15 and +2 degsC.
	*/
	return (pt_ice);
}

/***************************************************************************
% gsw_rho_ice                                                density of ice
%==========================================================================
%
% USAGE:
%  gsw_rho_ice(t,p)
%
% DESCRIPTION:
%  Calculates in-situ density of ice from in-situ temperature and pressure.
%  Note that the output, rho_ice, is density, not density anomaly;  that
%  is, 1000 kg/m^3 is not subracted from it.
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  rho_ice  =  in-situ density of ice (not density anomaly)      [ kg/m^3 ]
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
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_rho_ice(double t,double p)
{
	return (1.0/gsw_gibbs_ice(0,1,t,p));
}

/***************************************************************************
% gsw_SA_freezing_from_t               Absolute Salinity of seawater at the
%                                                      freezing temperature
%==========================================================================
%
% USAGE:
%  gsw_SA_freezing_from_t(t,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.
%  That is, the output is the Absolute Salinity of seawater, with the
%  fraction saturation_fraction of dissolved air, that is in equilibrium
%  with ice at in-situ temperature t and pressure p.  If the input values
%  are such that there is no positive value of Absolute Salinity for which
%  seawater is frozen, the output, SA_freezing, is set to NaN.
%
% INPUT:
%  t  =  in-situ Temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  SA_freezing =  Absolute Salinity of seawater when it freezes, for
%                 given input values of in situ temperature, pressure and
%                 air saturation fraction.                         [ g/kg ]
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
%    See section 3.33 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_sa_freezing_from_t(double t,double p,
                                          double saturation_fraction)
{
	int i_iter, number_of_iterations = 2;
	double f, sa, sa_mean, sa_old, t_freezing_zero_sa, tfreezing_sa;
	/**
	  This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	  differently in calculating the initial values of both SA and dCT_dSA.
	*/
	double  sa_cut_off = 2.5;
	/**
	  Find t > t_freezing_zero_SA.  If this is the case, the input values
	  represent seawater that is not frozen (at any positive SA).
	*/
	t_freezing_zero_sa = gsw_t_freezing(0.0,p,saturation_fraction);
	if (t > t_freezing_zero_sa)
		return (cppGSW_INVALID_VALUE);
	/**
	  This is the inital guess of SA using a purpose-built
	  polynomial in CT and p
	*/
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,NULL,&t);
	if (sa < -sa_cut_off)
		return (cppGSW_INVALID_VALUE);
	/**
	  Form the first estimate of tfreezing_SA, the derivative of
	  CT_freezing with respect to SA at fixed p.
	*/
	sa = gsw_max_d(sa,0.0);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
												&tfreezing_sa,NULL);
	/**     For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	  with one based on (t_freezing_zero_SA - t).
	*/
	if (fabs(sa) < sa_cut_off)
		sa = (t - t_freezing_zero_sa)/tfreezing_sa;
	/**
	-----------------------------------------------------------------------
	  Begin the modified Newton-Raphson method to find the root of
	  f = (t_freezing - t) = 0 for SA.
	-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
	{
		sa_old = sa;
		f = gsw_t_freezing(sa_old,p,saturation_fraction) - t;
		sa = sa_old - f/tfreezing_sa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_t_freezing_first_derivatives(sa_mean,p,saturation_fraction,
													&tfreezing_sa, NULL);
		sa = sa_old - f/tfreezing_sa;
	}

	if (gsw_sa_p_inrange(sa,p))
		return (sa);

	return (cppGSW_INVALID_VALUE);
}

/***************************************************************************
% gsw_SA_freezing_from_CT              Absolute Salinity of seawater at the
%                                                            freezing point
%==========================================================================
%
% USAGE:
%  gsw_SA_freezing_from_CT(CT,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.
%  That is, the output is the Absolute Salinity of seawater, with
%  Conservative Temperature CT, pressure p and the fraction
%  saturation_fraction of dissolved air, that is in equilibrium
%  with ice at the same in situ temperature and pressure.  If the input
%  values are such that there is no positive value of Absolute Salinity for
%  which seawater is frozen, the output, SA_freezing, is made a NaN.
%
% INPUT:
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction  =  the saturation fraction of dissolved air in
%                          seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  SA_freezing =  Absolute Salinity of seawater when it freezes, for
%                 given input values of its Conservative Temperature,
%                 pressure and air saturation fraction.            [ g/kg ]
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
%    See section 3.33 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_sa_freezing_from_ct(double ct,double p,
                                          double saturation_fraction)
{
	int i_iter, number_of_iterations = 3;
	double ct_freezing_zero_sa, f, ctfreezing_sa, sa, sa_mean, sa_old;
	/**
	  This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	  differently in calculating the initial values of both SA and dCT_dSA.
	*/
	double sa_cut_off = 2.5;
	/**
	  Find CT > CT_freezing_zero_SA.  If this is the case, the input values
	  represent seawater that is not frozen (for any positive SA).
	*/
	ct_freezing_zero_sa = gsw_ct_freezing(0.0,p,saturation_fraction);
	if (ct > ct_freezing_zero_sa)
		return (cppGSW_INVALID_VALUE);
	/**Form the first estimate of SA from a polynomial in CT and p*/
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
	if (sa < -sa_cut_off)
		return (cppGSW_INVALID_VALUE);
	/**
	  Form the first estimate of CTfreezing_SA,
	  the derivative of CT_freezing with respect to SA at fixed p.
	*/
	sa = gsw_max_d(sa,0.0);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,
												 &ctfreezing_sa, NULL);
	/**
	  For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	  with one based on (CT_freezing_zero_SA - CT).
	*/
	if (fabs(sa) < sa_cut_off)
		sa = (ct - ct_freezing_zero_sa)/ctfreezing_sa;
	/**
	------------------------------------------------------
	  Begin the modified Newton-Raphson method to solve
	  f = (CT_freezing - CT) = 0 for SA.
	------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
	{
		sa_old = sa;
		f = gsw_ct_freezing(sa,p,saturation_fraction) - ct;
		sa = sa_old - f/ctfreezing_sa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_ct_freezing_first_derivatives(sa_mean,p,saturation_fraction,
													 &ctfreezing_sa, NULL);
		sa = sa_old - f/ctfreezing_sa;
	}

	if (gsw_sa_p_inrange(sa,p))
		return (sa);

	return (cppGSW_INVALID_VALUE);
}

/***************************************************************************
% gsw_SA_freezing_from_t_poly          Absolute Salinity of seawater at the
%                                                     freezing point (poly)
%==========================================================================
%
% USAGE:
%  gsw_SA_freezing_from_t_poly(t,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.
%  That is, the output is the Absolute Salinity of seawater, with the
%  fraction saturation_fraction of dissolved air, that is in equilibrium
%  with ice at in-situ temperature t and pressure p.  If the input values
%  are such that there is no positive value of Absolute Salinity for which
%  seawater is frozen, the output, SA_freezing, is put equal to NaN.
%
% INPUT:
%  t  =  in-situ Temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  SA_freezing  =  Absolute Salinity of seawater when it freezes, for
%                  given input values of in situ temperature, pressure and
%                  air saturation fraction.                        [ g/kg ]
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
%    See section 3.33 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_sa_freezing_from_t_poly(double t,double p,
                                             double saturation_fraction)
{
	int i_iter, number_of_iterations = 5;
	double dt_dsa, sa, sa_old, sa_mean, t_freezing, t_freezing_zero_sa;
	/**
	  This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	  differently in calculating the initial values of both SA and dCT_dSA.
	*/
	double  sa_cut_off = 2.5;
	/**
	  Find t > t_freezing_zero_SA.  If this is the case, the input values
	  represent seawater that is not frozen (at any positive SA).
	*/
	t_freezing_zero_sa = gsw_t_freezing_poly(0.0,p,saturation_fraction);
	if (t > t_freezing_zero_sa)
		return (cppGSW_INVALID_VALUE);
	/**
	  This is the inital guess of SA using a purpose-built
	  polynomial in CT and p
	*/
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,NULL,&t);
	if (sa < -sa_cut_off)
		return (cppGSW_INVALID_VALUE);
	/**
	   Form the first estimate of dt_dSA, the derivative of t with respect
	   to SA at fixed p.
	*/
	sa = gsw_max_d(sa,0.0);
	gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,
													  &dt_dsa, NULL);
	/**
	  For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	  with one based on (t_freezing_zero_SA - t).
	*/
	if (fabs(sa) < sa_cut_off)
		sa = (t - t_freezing_zero_sa)/dt_dsa;
	/**
	-----------------------------------------------------------------------
	  Begin the modified Newton-Raphson method to find the root of
	  t_freezing = t for SA.
	-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
	{
		sa_old = sa;
		t_freezing = gsw_t_freezing_poly(sa_old,p,saturation_fraction);
		sa = sa_old - (t_freezing - t)/dt_dsa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_t_freezing_first_derivatives_poly(sa_mean,p,saturation_fraction,
														  &dt_dsa, NULL);
		sa = sa_old - (t_freezing - t)/dt_dsa;
	}

	if (gsw_sa_p_inrange(sa,p))
		return (sa);

	return (cppGSW_INVALID_VALUE);
}

/**************************************************************************
==========================================================================
  method: gsw_sa_freezing_estimate (p, saturation_fraction,*ct,*t)
==========================================================================
         an estimate of SA from a polynomial in CT and p
--------------------------------------------------------------------------
****************************************************************************/
double TeosIce::gsw_sa_freezing_estimate(double p,double saturation_fraction,
                                          double *ct,double *t)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double ctx, ctsat, sa,
			 /**note that aa = 0.502500117621d0/35.16504*/
			 aa = 0.014289763856964,
			 bb = 0.057000649899720,
			 p0  =  2.570124672768757e-1,
			 p1  = -1.917742353032266e1,
			 p2  = -1.413382858617969e-2,
			 p3  = -5.427484830917552e-1,
			 p4  = -4.126621135193472e-4,
			 p5  = -4.176407833276121e-7,
			 p6  =  4.688217641883641e-5,
			 p7  = -3.039808885885726e-8,
			 p8  = -4.990118091261456e-11,
			 p9  = -9.733920711119464e-9,
			 p10 = -7.723324202726337e-12,
			 p11 =  7.121854166249257e-16,
			 p12 =  1.256474634100811e-12,
			 p13 =  2.105103897918125e-15,
			 p14 =  8.663811778227171e-19;

	/** A very rough estimate of sa to get the saturated ct */
	if (ct != NULL)
	{
		sa = gsw_max_d(-(*ct + 9e-4*p)/0.06, 0.0);
		ctx = *ct;
	}
	else if (t != NULL)
	{
		sa = gsw_max_d(-(*t + 9e-4*p)/0.06, 0.0);
		ctx = gsw_ct_from_t(sa,*t,p);
	}
	else
	{
		return (0.0);
	}
	/**
	  CTsat is the estimated value of CT if the seawater were saturated with
	  dissolved air, recognizing that it actually has the air fraction
	  saturation_fraction; see McDougall, Barker and Feistel, 2014).
	*/
	ctsat = ctx - (1.0-saturation_fraction)*
			  (1e-3)*(2.4-aa*sa)*(1.0+bb*(1.0-sa/gtc.gsw_sso));

	return (p0 + p*(p2 + p4*ctsat + p*(p5 + ctsat*(p7 + p9*ctsat)
            + p*(p8  + ctsat*(p10 + p12*ctsat) + p*(p11 + p13*ctsat + p14*p))))
			   + ctsat*(p1 + ctsat*(p3 + p6*p)));
}

/***************************************************************************
% gsw_pot_enthalpy_from_pt_ice_poly          potential enthalpy of ice from
%                              potential temperature refered to the surface
%==========================================================================
%
% USAGE:
%  gsw_pot_enthalpy_from_pt_ice_poly(pt0_ice)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from potential temperature of
%  ice (whose reference sea pressure is zero dbar).  This is a
%  compuationally efficient polynomial fit to the potential enthalpy of
%  ice.
%
% INPUT:
%  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
%
% OUTPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
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
double TeosIce::gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice)
{
	int iteration;
	double df_dt, f, pot_enthalpy_ice, pot_enthalpy_ice_mid,
	         pot_enthalpy_ice_old,
			   p0 = -3.333601570157700e5,
			   p1 =  2.096693916810367e3,
			   p2 =  3.687110754043292,
			   p3 =  4.559401565980682e-4,
			   p4 = -2.516011957758120e-6,
			   p5 = -1.040364574632784e-8,
			   p6 = -1.701786588412454e-10,
			   p7 = -7.667191301635057e-13;

	/**initial estimate of the potential enthalpy.*/
	pot_enthalpy_ice = p0
	                     + pt0_ice*(p1
	                     + pt0_ice*(p2
	                     + pt0_ice*(p3
								+ pt0_ice*(p4
								+ pt0_ice*(p5
								+ pt0_ice*(p6
								+ pt0_ice*p7))))));

	df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice);

	for (iteration = 1; iteration <= 5; iteration++)
	{
		pot_enthalpy_ice_old = pot_enthalpy_ice;

		f = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old) - pt0_ice;
		pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
		pot_enthalpy_ice_mid = 0.5*(pot_enthalpy_ice+pot_enthalpy_ice_old);
		df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice_mid);
		pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
	}
	/**
	The error of this fit ranges between -6e-3 and 6e-3 J/kg over the
	potential temperature range of -100 to 2 deg C, or the potential
	enthalpy range of -5.7 x 10^5 to -3.3 x 10^5 J/kg.
	*/
	return (pot_enthalpy_ice);
}

/**************************************************************************
==========================================================================
method: gsw_ct_freezing_exact(sa,p,saturation_fraction)
==========================================================================
   Calculates the Conservative Temperature at which seawater freezes.  The
   Conservative Temperature freezing point is calculated from the exact
   in-situ freezing temperature which is found by a modified Newton-Raphson
   iteration (McDougall and Wotherspoon, 2013) of the equality of the
   chemical potentials of water in seawater and in ice.
   An alternative GSW function, gsw_CT_freezing_poly, it is based on a
   computationally-efficient polynomial, and is accurate to within -5e-4 K
   and 6e-4 K, when compared with this function.
   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   CT_freezing : Conservative Temperature at freezing of seawater [ deg C ]
---------------------------------------------------------------------------
***************************************************************************/
double TeosIce::gsw_ct_freezing_exact(double sa,double p,double saturation_fraction)
{
	double t_freezing = gsw_t_freezing_exact(sa,p,saturation_fraction);

	return (gsw_ct_from_t(sa,t_freezing,p));
}


/***************************************************************************
% gsw_geo_strf_dyn_height_pc                     dynamic height anomaly for
%                            piecewise constant profiles (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_geo_strf_dyn_height_pc(SA,CT,delta_p)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume
%  anomaly from the the sea surface pressure (0 Pa) to the pressure p.
%  This function, gsw_geo_strf_dyn_height_pc, is to used when the
%  Absolute Salinity and Conservative Temperature are piecewise constant in
%  the vertical over sucessive pressure intervals of delta_p (such as in
%  a forward "z-coordinate" ocean model).  "geo_strf_dyn_height_pc" is
%  the dynamic height anomaly with respect to the sea surface.  That is,
%  "geo_strf_dyn_height_pc" is the geostrophic streamfunction for the
%  difference between the horizontal velocity at the pressure concerned, p,
%  and the horizontal velocity at the sea surface.  Dynamic height anomaly
%  is the geostrophic streamfunction in an isobaric surface.  The reference
%  values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
%  and CT = 0 deg C.  The output values of geo_strf_dyn_height_pc are
%  given at the mid-point pressures, p_mid, of each layer in which SA and
%  CT are vertically piecewice constant (pc).  This function calculates
%  enthalpy using the computationally-efficient 75-term expression for
%  specific volume of Roquet et al., (2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  CT       =  Conservative Temperature (ITS-90)                  [ deg C ]
%  delta_p  =  difference in sea pressure between the deep and     [ dbar ]
%              shallow extents of each layer in which SA and CT
%              are vertically constant. delta_p must be positive.
%
%  Note. sea pressure is absolute pressure minus 10.1325 dbar.
%
% OUTPUT:
%  geo_strf_dyn_height_pc =  dynamic height anomaly             [ m^2/s^2 ]
%  p_mid                  =  mid-point pressure in each layer      [ dbar ]
%
% AUTHOR:
%  Trevor McDougall and Claire Roberts-Thomson         [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J. and A. Klocker, 2010: An approximate geostrophic
%   streamfunction for use in density surfaces. Ocean Modelling, 32,
%   105-117.
%    See section 8 of this paper.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double *TeosIce::gsw_geo_strf_dyn_height_pc(double *sa,double *ct,
                                             double *delta_p,int n_levels,
                                             double *geo_strf_dyn_height_pc,
                                             double *p_mid)
{
    int i, np;
    double *delta_h=NULL, delta_h_half, dyn_height_deep=0.0,
      *p_deep=NULL, *p_shallow=NULL;

    for (i=0; i<n_levels; i++)
    {
      if (delta_p[i] < 0.0)
      {
         return (NULL);
      }
    }
    np = n_levels;
    delta_h = new double[3*np];
    p_deep = delta_h+np;
    p_shallow = p_deep+np;

    for (i=0; i<np; i++)
    {
      p_deep[i] = (i==0)? delta_p[0] : p_deep[i-1] + delta_p[i];
      p_shallow[i] = p_deep[i] - delta_p[i];
      delta_h[i] = gsw_enthalpy_diff(sa[i],ct[i],p_shallow[i],p_deep[i]);
    }

    for (i=0; i<np; i++)
    {
      dyn_height_deep = dyn_height_deep - delta_h[i];
      /** This is Phi minus Phi_0 of Eqn. (3.32.2) of IOC et al. (2010).*/
      p_mid[i] = 0.5*(p_shallow[i]  + p_deep[i]);
      delta_h_half = gsw_enthalpy_diff(sa[i],ct[i],p_mid[i],p_deep[i]);
      geo_strf_dyn_height_pc[i] = gsw_enthalpy_sso_0(p_mid[i]) +
               dyn_height_deep + delta_h_half;
    }

    if (delta_h) delete []delta_h;

    return (geo_strf_dyn_height_pc);
}

/***************************************************************************
% gsw_IPV_vs_fNsquared_ratio              ratio of the vertical gradient of
%                       potential density (with reference pressure, p_ref),
%                            to the vertical gradient of locally-referenced
%                                      potential density (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_IPV_vs_fNsquared_ratio(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the ratio of the vertical gradient of potential density to
%  the vertical gradient of locally-referenced potential density.  This
%  ratio is also the ratio of the planetary Isopycnal Potential Vorticity
%  (IPV) to f times N^2, hence the name for this variable,
%  IPV_vs_fNsquared_ratio (see Eqn. (3.20.17) of IOC et al. (2010)).
%  The reference sea pressure, p_ref, of the potential density surface must
%  have a constant value.
%
%  IPV_vs_fNsquared_ratio is evaluated at the mid pressure between the
%  individual data points in the vertical.  This function uses the
%  computationally-efficient 75-term expression for specific volume in
%  terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA    = Absolute Salinity                                       [ g/kg ]
%  CT    = Conservative Temperature (ITS-90)                      [ deg C ]
%  p     = sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref = reference sea pressure of the potential density surface
%         ( i.e. absolute reference pressure - 10.1325 dbar )      [ dbar ]
%
% OUTPUT:
%  IPV_vs_fNsquared_ratio = The ratio of the vertical gradient of
%          potential density referenced to p_ref, to the vertical
%          gradient of locally-referenced potential density.  It is
%          output on the same vertical (M-1)xN grid as p_mid.
%          IPV_vs_fNsquared_ratio is dimensionless.            [ unitless ]
%  p_mid = mid pressure between the individual points of the p grid.
%          That is, p_mid is on a (M-1)xN grid.                    [ dbar ]
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
%    See Eqn. (3.20.5) of this TEOS-10 Manual.
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
void TeosIce::gsw_ipv_vs_fnsquared_ratio(double *sa,double *ct,double *p,
                                          double p_ref,int nz,
                                          double *ipv_vs_fnsquared_ratio,
                                          double *p_mid)
{
    int     k;
    double  dsa, sa_mid, dct, ct_mid;
    double  alpha_mid, beta_mid;
    double  alpha_pref, beta_pref, numerator, denominator;

    if (nz < 2)
    {
      *p_mid = *ipv_vs_fnsquared_ratio = cppGSW_INVALID_VALUE;
      return;
    }

    for (k = 0; k < nz-1; k++)
    {
      dsa = (sa[k] - sa[k+1]);
      dct = (ct[k] - ct[k+1]);
      sa_mid = 0.5*(sa[k] + sa[k+1]);
      ct_mid = 0.5*(ct[k] + ct[k+1]);
      p_mid[k] = 0.5*(p[k] + p[k+1]);
      alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
      beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);
      alpha_pref = gsw_alpha(sa_mid,ct_mid,p_ref);
      beta_pref = gsw_beta(sa_mid,ct_mid,p_ref);
      numerator = dct*alpha_pref - dsa*beta_pref;
      denominator = dct*alpha_mid - dsa*beta_mid;

      if (denominator == 0.0)
      {
         ipv_vs_fnsquared_ratio[k] = cppGSW_INVALID_VALUE;
      }
      else
      {
         ipv_vs_fnsquared_ratio[k] = numerator/denominator;
      }
   }
}

/***************************************************************************
% gsw_melting_ice_equilibrium_SA_CT_ratio_poly    ratio of SA to CT changes
%                        when ice melts into a large mass of seawater, with
%                       both the seawater and ice temperatures being almost
%                      equal to the equilibrium freezing temperature (poly)
%==========================================================================
%
% USAGE:
%  gsw_melting_ice_equilibrium_SA_CT_ratio_poly(SA,p)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when ice melts into seawater
%  with both the seawater and the seaice temperatures being almost equal to
%  the equilibrium freezing temperature.  It is assumed that a small mass
%  of ice melts into an infinite mass of seawater.  If indeed the
%  temperature of the seawater and the ice were both equal to the freezing
%  temperature, then no melting or freezing would occur; an imbalance
%  between these three temperatures is needed for freezing or melting to
%  occur (the three temperatures being (1) the seawater temperature,
%  (2) the ice temperature, and (3) the freezing temperature.
%
%  The output, melting_ice_equilibrium_SA_CT_ratio, is dSA/dCT rather than
%  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is zero
%  whereas dCT/dSA would be infinite.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  melting_ice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
%                                changes when ice melts into seawater, with
%                                the seawater and seaice being close to the
%                                freezing temperature.         [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (16) of this manuscript.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa,double p)
{
	double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing_poly(sa,p,saturation_fraction);
	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	/**note that h_hat_ct is equal to cp0*(273.15 + t)/(273.15 + pt0)*/
	return (sa*h_hat_ct / (h - h_ih - sa*h_hat_sa));
}

/**************************************************************************
% gsw_melting_ice_into_seawater             Absolute Salinity, Conservative
%                  Temperature and final ice mass fraction when ice of mass
%                fraction w_Ih and temperature t_Ih is melted into seawater
%==========================================================================
%
% USAGE:
%  gsw_melting_ice_into_seawater(SA,CT,p,w_Ih,t_Ih)
%
% DESCRIPTION:
%  Calculates the final Absolute Salinity, final Conservative Temperature
%  and final ice mass fraction that results when a given mass fraction of
%  ice melts and is mixed into seawater whose properties are (SA,CT,p).
%  This code takes the seawater to contain no dissolved air.
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
%  present in the final state.  In this case the final state consists of
%  only seawater rather than being an equlibrium mixture of seawater and
%  ice which occurs when w_Ih_final is positive.  Note that when
%  w_Ih_final = 0, the final seawater is not at the freezing temperature.
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure at which the melting occurs                [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  w_Ih =  mass fraction of ice, that is the mass of ice divided by the
%          sum of the masses of ice and seawater.  That is, the mass of
%          ice divided by the mass of the final mixed fluid.
%          w_Ih must be between 0 and 1.                       [ unitless ]
%  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
%
%  SA, CT, w_Ih and t_Ih must all have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where where SA, CT,
%  w_Ih and t_Ih are MxN.
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state,
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at
%                 thermodynamic equilibrium.  If this ice mass fraction is
%                 zero there is no ice in the final state which consists
%                 only of seawater which is warmer than the freezing
%                 temperature.                                   [unitless]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_melting_ice_into_seawater(double sa,double ct,double p,
                                             double w_ih,double t_ih,
                                             double *sa_final,
                                             double *ct_final,
                                             double *w_ih_final)
{
    double  ctf, h_bulk, sa_bulk, tf_ih;
    double  saturation_fraction = 0.0;

    ctf = gsw_ct_freezing(sa,p,saturation_fraction);

    if (ct < ctf)
    {
      /**The seawater ct input is below the freezing temp*/
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;

      return;
    }

    tf_ih = gsw_t_freezing(0.0,p,saturation_fraction) - 1e-6;

    if (t_ih > tf_ih)
    {
      /**
         t_ih input exceeds the freezing temp.
         The 1e-6 C buffer in the allowable
         t_Ih is to ensure that there is some ice Ih in the sea ice.
      */
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;

      return;
    }

    sa_bulk = (1.0 - w_ih)*sa;
    h_bulk = (1.0 - w_ih)*gsw_enthalpy_ct_exact(sa,ct,p)
           + w_ih*gsw_enthalpy_ice(t_ih,p);

    gsw_frazil_properties(sa_bulk,h_bulk,p,sa_final,ct_final,w_ih_final);

    if (*sa_final > cppGSW_ERROR_LIMIT)
    {
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;
    }

   return;
}

/**************************************************************************
% gsw_Nsquared             buoyancy (Brunt-Vaisala) frequency squared (N^2)
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  gsw_Nsquared(SA,CT,p,lat)
%
% DESCRIPTION:
%  Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala
%  frequency squared) at the mid pressure from the equation,
%
%           2      2     beta.dSA - alpha.dCT
%         N   =  g  . -------------------------
%                         specvol_local.dP
%
%  The pressure increment, dP, in the above formula is in Pa, so that it is
%  10^4 times the pressure increment dp in dbar.
%
%  Note. This routine uses specvol from "gsw_specvol", which is based on
%  the computationally efficient expression for specific volume in terms of
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
% OPTIONAL:
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ]
%  Note. If lat is not supplied, a default gravitational acceleration
%    of 9.7963 m/s^2 (Griffies, 2004) will be applied.
%
% OUTPUT:
%  N2     =  Brunt-Vaisala Frequency squared  (M-1xN)             [ 1/s^2 ]
%  p_mid  =  Mid pressure between p grid      (M-1xN)              [ dbar ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (30th June, 2020)
%
% REFERENCES:
%  Griffies, S.M., 2004: Fundamentals of Ocean Climate Models. Princeton,
%   NJ: Princeton University Press, 518 pp + xxxiv.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.10 and Eqn. (3.10.2) of this TEOS-10 Manual.
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
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_nsquared(double *sa,double *ct,double *p,double *lat,
                           int nz,double *n2,double *p_mid)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */
   int     k;
   double  p_grav, n_grav, grav_local, dsa, sa_mid, dct, ct_mid,
      dp, rho_mid, alpha_mid, beta_mid;

   if (nz < 2)
   {
      return;
   }

   p_grav  = gsw_grav(lat[0],p[0]);

   for (k = 0; k < nz-1; k++)
   {
       n_grav = gsw_grav(lat[k+1],p[k+1]);
       grav_local = 0.5*(p_grav + n_grav);
       dsa = (sa[k+1] - sa[k]);
       sa_mid = 0.5*(sa[k] + sa[k+1]);
       dct = (ct[k+1] - ct[k]);
       ct_mid      = 0.5*(ct[k] + ct[k+1]);
       dp = (p[k+1] - p[k]);
       p_mid[k] = 0.5*(p[k] + p[k+1]);
       rho_mid = gsw_rho(sa_mid,ct_mid,p_mid[k]);
       alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
       beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);
       n2[k] = (grav_local*grav_local)*(rho_mid/(gtc.db2pa*dp))*
           (beta_mid*dsa - alpha_mid*dct);

       p_grav = n_grav;
   }
}

/***************************************************************************
% gsw_seaice_fraction_to_freeze_seawater  sea ice mass fraction, which when
%                                          melted into seawater, brings the
%                                      seawater to the freezing temperature
%==========================================================================
%
% USAGE:
%  gsw_seaice_fraction_to_freeze_seawater(SA,CT,p,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the mass fraction of sea ice (mass of sea ice divided by mass
%  of sea ice plus seawater), which, when melted into seawater having the
%  properties (SA,CT,p) causes the final seawater to be at the freezing
%  temperature.  The other outputs are the Absolute Salinity and
%  Conservative Temperature of the final seawater.
%
% INPUT:
%  SA        =  Absolute Salinity of seawater                      [ g/kg ]
%  CT        =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
%  p         =  sea pressure                                       [ dbar ]
%            ( i.e. absolute pressure - 10.1325 dbar )
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of
%               salt in sea ice, expressed in g of salt per kg of sea ice.
%                                                                  [ g/kg ]
%  t_seaice  =  in-situ temperature of the sea ice at pressure p (ITS-90)
%                                                                 [ deg C ]
%
% OUTPUT:
%  SA_freeze  =  Absolute Salinity of seawater after the mass fraction of
%                sea ice, w_seaice, at temperature t_seaice has melted into
%                the original seawater, and the final mixture is at the
%                freezing temperature of seawater.                 [ g/kg ]
%
%  CT_freeze  =  Conservative Temperature of seawater after the mass
%                fraction, w_seaice, of sea ice at temperature t_seaice has
%                melted into the original seawater, and the final mixture
%                is at the freezing temperature of seawater.      [ deg C ]
%
%  w_seaice   =  mass fraction of sea ice, at SA_seaice and t_seaice,
%                which, when melted into seawater at (SA,CT,p) leads to the
%                final mixed seawater being at the freezing temperature.
%                This output is between 0 and 1.                 [unitless]
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
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall T.J. and S.J. Wotherspoon, 2013: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (23) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
void TeosIce::gsw_seaice_fraction_to_freeze_seawater(double sa,double ct,
                                                      double p,double sa_seaice,
                                                      double t_seaice,double *sa_freeze,
                                                      double *ct_freeze,double *w_seaice)
{
	int number_of_iterations;
	double ctf, ctf_mean, ctf_old, ctf_plus1, ctf_zero,
			 dfunc_dsaf, func, func_plus1, func_zero, h, h_brine,
			 h_ih, sa_freezing, saf, saf_mean, saf_old,
			 salt_ratio, tf_sa_seaice, h_hat_sa, h_hat_ct, ctf_sa,
			 sa0 = 0.0, saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);

	if (ct < ctf)
	{
		/**The seawater ct input is below the freezing temp*/
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}
	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)
	{
		/**
		The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
		some ice Ih in the sea ice.   Without this buffer, that is if t_seaice
		is allowed to be exactly equal to tf_sa_seaice, the sea ice is
		actually 100% brine at Absolute Salinity of SA_seaice.
		*/
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}
	sa_freezing = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);
	if (sa_freezing > cppGSW_ERROR_LIMIT)
	{
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}
	h_brine = gsw_enthalpy_t_exact(sa_freezing,t_seaice,p);
	salt_ratio = sa_seaice/sa_freezing;
	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	ctf_plus1 = gsw_ct_freezing(sa+1.0,p,saturation_fraction);
	func_plus1 = (sa - sa_seaice)
					 *(gsw_enthalpy_ct_exact(sa+1.0,ctf_plus1,p)
						- h) - (h - h_ih) + salt_ratio*(h_brine - h_ih);
	ctf_zero = gsw_ct_freezing(sa0,p,saturation_fraction);
	func_zero = (sa - sa_seaice)
					*(gsw_enthalpy_ct_exact(sa0,ctf_zero,p) - h)
					+ sa*((h - h_ih) - salt_ratio*(h_brine - h_ih));
	saf = -(sa+1.0)*func_zero/(func_plus1 - func_zero);
	/**initial guess of saf*/
	ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	gsw_enthalpy_first_derivatives_ct_exact(saf,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(saf,p,saturation_fraction,
												 &ctf_sa, NULL);
	dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
					 - (h - h_ih) + salt_ratio*(h_brine - h_ih);
	for (number_of_iterations = 1; number_of_iterations <= 4;
			number_of_iterations++)
	{
		saf_old = saf;
		ctf_old = ctf;
		func = (sa - sa_seaice)
				 *(gsw_enthalpy_ct_exact(saf_old,ctf_old,p) - h)
				 - (saf_old - sa)*((h - h_ih) - salt_ratio*(h_brine - h_ih));
		saf = saf_old - func/dfunc_dsaf;
		saf_mean = 0.5*(saf + saf_old);
		ctf_mean = gsw_ct_freezing(saf_mean,p,saturation_fraction);
		gsw_enthalpy_first_derivatives_ct_exact(saf_mean,ctf_mean,p,
															 &h_hat_sa,&h_hat_ct);
		gsw_ct_freezing_first_derivatives(saf_mean,p,saturation_fraction,
													 &ctf_sa, NULL);
		dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
						 - (h - h_ih) + salt_ratio*(h_brine - h_ih);
		saf = saf_old - func/dfunc_dsaf;
		ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	}
	/**
	 After these 4 iterations of this modified Newton-Raphson method, the
	  errors in SA_freeze is less than 1.5x10^-12 g/kg, in CT_freeze is less than
	  2x10^-13 deg C and in w_seaice is less than 2.8x10^-13 which represent machine
	  precision for these calculations.
	*/
	*sa_freeze = saf;
	*ct_freeze = ctf;
	*w_seaice = (h - gsw_enthalpy_ct_exact(*sa_freeze,*ct_freeze,p)) /
					(h - h_ih - salt_ratio*(h_brine - h_ih));
}

/***************************************************************************
% gsw_sound_speed_ice                               compression sound speed
%==========================================================================
%
% USAGE:
%  gsw_sound_speed_ice(t,p)
%
% DESCRIPTION:
%  Calculates the compression speed of sound in ice.
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  sound_speed_ice  =  compression speed of sound in ice            [ m/s ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_sound_speed_ice(double t,double p)
{
	double gi_tp, gi_tt;

	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);

	return (gsw_gibbs_ice(0,1,t,p) *
			  sqrt(gi_tt/(gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))));
}

/***************************************************************************
% gsw_t_from_pt0_ice                             in-situ temperature of ice
% =========================================================================
%
% USAGE:
%  gsw_t_from_pt0_ice(pt0_ice,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature from the potential temperature of ice Ih
%  with reference pressure, p_ref, of 0 dbar (the surface), and the
%  in-situ pressure.
%
% INPUT:
%  pt0_ice  =  potential temperature of ice Ih with reference pressure of
%              zero dbar (ITS-90)                                 [ deg C ]
%  p        =  sea pressure                                        [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
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
%    See appendix I of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_t_from_pt0_ice(double pt0_ice,double p)
{
	double  p0 = 0.0;
	return (gsw_pt_from_t_ice(pt0_ice,p0,p));
}

/***************************************************************************
% gsw_melting_ice_SA_CT_ratio_poly           ratio of SA to CT changes when
%                                            ice melts into seawater (poly)
%==========================================================================
%
% USAGE:
%  gsw_melting_ice_SA_CT_ratio_poly(SA,CT,p,t_Ih)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when ice melts into seawater.
%  It is assumed that a small mass of ice melts into an infinite mass of
%  seawater.  Because of the infinite mass of seawater, the ice will always
%  melt.
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
%  This is done so that when SA = 0, the output, dSA/dCT is zero whereas
%  dCT/dSA would be infinite.
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA   =  Absolute Salinity of seawater                           [ g/kg ]
%  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
%  p    =  sea pressure at which the melting occurs                [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
%
% OUTPUT:
%  melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
%                            into a large mass of seawater
%                                                          [ g kg^-1 K^-1 ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (13) of this manuscript.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_ice_sa_ct_ratio_poly(double sa,double ct,
                                                   double p,double t_ih)
{
	double  ctf, h, h_ih, tf, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temperature*/
		return (cppGSW_INVALID_VALUE);
	}

	tf = gsw_t_freezing_poly(0.0,p,saturation_fraction);

	if (t_ih > tf)
	{
		/**t_ih exceeds the freezing temperature at sa = 0*/
		return (cppGSW_INVALID_VALUE);
	}

	h = gsw_enthalpy(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);
	gsw_enthalpy_first_derivatives(sa,ct,p,&h_hat_sa,&h_hat_ct);
	/**Note that h_hat_CT is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct/(h - h_ih - sa*h_hat_sa));
}

/***************************************************************************
% gsw_melting_seaice_SA_CT_ratio             ratio of SA to CT changes when
%                                               sea ice melts into seawater
%==========================================================================
%
% USAGE:
%  gsw_melting_seaice_SA_CT_ratio(SA,CT,p,SA_seaice,t_seaice)
%
% DESCRIPTION:
%  Calculates the ratio of SA to CT changes when sea ice melts into
%  seawater.  It is assumed that a small mass of sea ice melts into an
%  infinite mass of seawater.  Because of the infinite mass of seawater,
%  the sea ice will always melt.
%
%  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
%  and 12 g/kg of salt (defined as the mass of salt divided by the mass of
%  ice Ih plus brine) and this programme returns NaN's if the input
%  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
%  usually this would imply that the pressure p should be zero, as sea ice
%  only occurs near the sea surface.  The code does not impose that p = 0
%  if SA_seaice is non-zero.  Rather, this is left to the user.
%
%  The Absolute Salinity, SA_brine, of the brine trapped in little pockets
%  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
%  surrounds these pockets.  As the seaice temperature, t_seaice, may be
%  less than the freezing temperature, SA_brine is usually greater than the
%  Absolute Salinity of the seawater at the time and place when and where
%  the sea ice was formed.  So usually SA_brine will be larger than SA.
%
%  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
%  This is done so that when (SA - seaice_SA) = 0, the output, dSA/dCT is
%  zero whereas dCT/dSA would be infinite.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
%  p   =  sea pressure at which the melting occurs                 [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction
%               of salt in sea ice expressed in g of salt per kg of
%               sea ice                                            [ g/kg ]
%  t_seaice =   the in-situ temperature of the sea ice (ITS-90)   [ deg C ]
%
% OUTPUT:
%  melting_seaice_SA_CT_ratio = the ratio dSA/dCT of SA to CT changes when
%                sea ice melts into a large mass of seawater   [ g/(kg K) ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%    See Eqn. (28) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_melting_seaice_sa_ct_ratio(double sa,double ct,double p,
                                                double sa_seaice,double t_seaice)
{
	double  ctf, delsa, h, h_brine, h_ih, sa_brine,
			  tf_sa_seaice, h_hat_sa, h_hat_ct,
			  saturation_fraction = 0.0;

	if (sa_seaice < 0.0 || sa_seaice > 15.0)
	{
		return (cppGSW_INVALID_VALUE);
	}

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);

   /*the seawater ct input is below the freezing temp*/
	if (ct < ctf)
	{
		return (cppGSW_INVALID_VALUE);
	}

	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)
	{
	   /*t_seaice exceeds the freezing sa*/
		return (cppGSW_INVALID_VALUE);
	}

	/**
	------------------------------------------------------------------------
	The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
	some ice Ih in the sea ice.  Without this buffer, that is if t_seaice
	is allowed to be exactly equal to tf_sa_seaice, the sea ice is actually
	100% brine at Absolute Salinity of SA_seaice.
	------------------------------------------------------------------------
	 */
	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,&h_hat_sa,&h_hat_ct);
	sa_brine = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);
	if (sa_brine > cppGSW_ERROR_LIMIT)
	{
		return (cppGSW_INVALID_VALUE);
	}

	h_brine = gsw_enthalpy_t_exact(sa_brine,t_seaice,p);
	delsa = sa - sa_seaice;

	return (h_hat_ct*delsa /
			  (h - h_ih - delsa*h_hat_sa - sa_seaice*(h_brine - h_ih)/sa_brine));
}

/**************************************************************************
% gsw_deltaSA_from_SP                             Absolute Salinity Anomaly
%                                                   from Practical Salinity
%==========================================================================
%
% USAGE:
%  deltaSA = gsw_deltaSA_from_SP(SP,p,long,lat)
%
% DESCRIPTION:
%  Calculates Absolute Salinity Anomaly from Practical Salinity.  Since SP
%  is non-negative by definition, this function changes any negative input
%  values of SP to be zero.
%
% INPUT:
%  SP   =  Practical Salinity  (PSS-78)                        [ unitless ]
%  p    =  sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 ... +90 ]
%
%  p, lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SP is MxN.
%
% OUTPUT:
%  deltaSA  =  Absolute Salinity Anomaly                           [ g/kg ]
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
%    See section 2.5 and appendices A.4 and A.5 of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and
%   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
%   Ocean Science, 8, 1117-1128.
%   http://www.ocean-sci.net/8/1117/2012/os-8-1117-2012.pdf
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_deltasa_from_sp(double sp,double p,double lon,double lat)
{
	double  res = gsw_sa_from_sp(sp,p,lon,lat) - gsw_sr_from_sp(sp);

	if (res > cppGSW_ERROR_LIMIT)
	{
		res = cppGSW_INVALID_VALUE;
	}

	return (res);
}

/***************************************************************************
% gsw_SA_freezing_from_CT_poly             Absolute Salinity of seawater at
%                                                 the freezing point (poly)
%==========================================================================
%
% USAGE:
%  gsw_SA_freezing_from_CT_poly(CT,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.
%  That is, the output is the Absolute Salinity of seawater, with the
%  fraction saturation_fraction of dissolved air, that is in equilibrium
%  with ice at Conservative Temperature CT and pressure p.  If the input
%  values are such that there is no positive value of Absolute Salinity for
%  which seawater is frozen, the output, SA_freezing, is put equal to NaN.
%
% INPUT:
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction  =  the saturation fraction of dissolved air in
%                          seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default
%    is 0, air free)
%
% OUTPUT:
%  SA_freezing =  Absolute Salinity of seawater when it freezes, for
%                 given input values of Conservative Temperature
%                 pressure and air saturation fraction.            [ g/kg ]
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
%    See section 3.33 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_sa_freezing_from_ct_poly(double ct,double p,
                                                double saturation_fraction)
{
	int i_iter, number_of_iterations = 2;
	double ct_freezing, ct_freezing_zero_sa, dct_dsa, sa, sa_old, sa_mean;
	/**
	  This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	  differently in calculating the initial values of both SA and dCT_dSA.
	*/
	double sa_cut_off = 2.5;
	/**
	  Find CT > CT_freezing_zero_SA.  If this is the case, the input values
	  represent seawater that is not frozen (at any positive SA).
	*/
	ct_freezing_zero_sa = gsw_ct_freezing_poly(0.0,p,saturation_fraction);
	if (ct > ct_freezing_zero_sa)
		return (cppGSW_INVALID_VALUE);
	/**Form the first estimate of SA from a polynomial in CT and p */
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
	if (sa < -sa_cut_off)
		return (cppGSW_INVALID_VALUE);
	/**
	  Form the first estimate of dCT_dSA, the derivative of CT with respect
	  to SA at fixed p.
	*/
	sa = gsw_max_d(sa,0.0);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,
														&dct_dsa, NULL);
	/**
	  For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	  with one based on (CT_freezing_zero_SA - CT).
	*/
	if (fabs(sa) < sa_cut_off)
		sa = (ct - ct_freezing_zero_sa)/dct_dsa;
	/**
	-----------------------------------------------------------------------
	  Begin the modified Newton-Raphson method to solve the root of
	  CT_freezing = CT for SA.
	-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
	{
		sa_old = sa;
		ct_freezing = gsw_ct_freezing_poly(sa_old,p,saturation_fraction);
		sa = sa_old - (ct_freezing - ct)/dct_dsa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_ct_freezing_first_derivatives_poly(sa_mean,p,
															saturation_fraction, &dct_dsa, NULL);
		sa = sa_old - (ct_freezing - ct)/dct_dsa;
	}
	if (gsw_sa_p_inrange(sa,p))
   {
      return (sa);
   }

	return (cppGSW_INVALID_VALUE);
}

/***************************************************************************
% gsw_specvol_from_pot_enthalpy_ice                    specific volume from
%                                                 potential enthalpy of ice
%==========================================================================
%
% USAGE:
%  specvol_ice = gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice from the potential enthalpy
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_specvol_from_pot_enthalpy_ice(double pot_enthalpy_ice, double p)
{
   double pt0_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice);

   double t_ice = gsw_t_from_pt0_ice(pt0_ice, p);

   return gsw_specvol_ice(t_ice, p);
}

/***************************************************************************
% gsw_pt_from_pot_enthalpy_ice      potential temperature of ice refered to
%                                the surface from potential enthalpy of ice
%==========================================================================
%
% USAGE:
%  gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice)
%
% DESCRIPTION:
%  Calculates the potential temperature of ice from the potential enthalpy
%  of ice.  The reference sea pressure of both the potential temperature
%  and the potential enthalpy is zero dbar.
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% OUTPUT:
%  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
***************************************************************************/
double TeosIce::gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice)
{
	int iteration;
	double df_dt, f, mod_pot_enthalpy_ice, pt0_cold_ice, recip_df_dt,
			 pt0_cold_ice_old, pt0_ice, pt0_ice_old, ptm_cold_ice, ptm_ice;
	double  h00 = -6.320202333358860e5, /*gsw_enthalpy_ice(-gtc.gsw_t0,0)*/
			  p0 = 0.0;

	mod_pot_enthalpy_ice = gsw_max_d(pot_enthalpy_ice,h00);
	if (mod_pot_enthalpy_ice >= -5.1e5)
	{
		/**
		   For input potential enthalpies greater than -5.1e-5, the above part of
		   the code gives the output potential temperature of ice accurate to
		   1e-13 degrees C.
		*/
		pt0_ice = gsw_pt_from_pot_enthalpy_ice_poly(mod_pot_enthalpy_ice);
		/**
		The variable "df_dt" below is the derivative of the above polynomial
		with respect to pot_enthalpy_ice.  This is the initial value of the
		derivative of the method f.
		*/
		recip_df_dt =
			gsw_pt_from_pot_enthalpy_ice_poly_dh(mod_pot_enthalpy_ice);
		pt0_ice_old = pt0_ice;
		f = gsw_pot_enthalpy_from_pt_ice(pt0_ice_old)
			 - mod_pot_enthalpy_ice;
		pt0_ice = pt0_ice_old - f*recip_df_dt;
		ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
		recip_df_dt = 1.0/gsw_cp_ice(ptm_ice,p0);
		pt0_ice = pt0_ice_old - f*recip_df_dt;
	}
	else
	{
		/**
		   For  pot_enthalpy_ice < -5.1e5 (or pt0_ice less than about -100 deg c)
		   these temperatures are less than those found in nature on planet earth
		*/
		pt0_cold_ice = gsw_pt0_cold_ice_poly(mod_pot_enthalpy_ice);
		df_dt = gsw_cp_ice(pt0_cold_ice+0.02,p0);
		/**         the heat capacity, cp, is
		  evaluated at 0.02 c greater than usual in order to avoid stability
		  issues and to ensure convergence near zero absolute temperature.
		*/
		for (iteration = 1; iteration <= 6; iteration++)
		{
			pt0_cold_ice_old = pt0_cold_ice;
			f = gsw_pot_enthalpy_from_pt_ice(pt0_cold_ice_old)
				 - mod_pot_enthalpy_ice;
			pt0_cold_ice = pt0_cold_ice_old - f/df_dt;
			ptm_cold_ice = 0.5*(pt0_cold_ice + pt0_cold_ice_old);
			df_dt = gsw_cp_ice(ptm_cold_ice+0.02,p0);
			/**note the extra 0.02 c here as well*/
			pt0_cold_ice = pt0_cold_ice_old - f/df_dt;
		}
		pt0_ice = pt0_cold_ice;
	}
	/**
	The potential temerature has a maximum error as listed in the table below.
	   potential temerature error (deg C)  |  @ potential temerature (deg C)
	--------------------------------------|---------------------------------
	       0.012       |     -273.15 to -273.12
	          4 x 10^-4          |     -232.12 to -273.0
	         2.5 x 10^-6         |     -273
	          7 x 10^-9          |     -272
	        3.7 x 10^-10         |     -270
	          6 x 10^-11         |     -268
	         2.5 x 10^11         |     -266
	         3 x 10^-12          |     -260
	         7 x 10^-13          |     -250
	        2.2 x 10^-13         |     -220
	        1.7 x 10^-13         |    >= -160
	  Note.  The above errors in each temperature range are machine precisions
	  for this calculation.
	*/
	return (pt0_ice);
}

/**************************************************************************
==========================================================================
method: gsw_t_freezing_exact (sa,p,saturation_fraction)
==========================================================================
   Calculates the in-situ temperature at which seawater freezes. The
   in-situ temperature freezing point is calculated from the exact
   in-situ freezing temperature which is found by a modified Newton-Raphson
   iteration (McDougall and Wotherspoon, 2013) of the equality of the
   chemical potentials of water in seawater and in ice.
   An alternative GSW function, gsw_t_freezing_poly, it is based on a
   computationally-efficient polynomial, and is accurate to within -5e-4 K
   and 6e-4 K, when compared with this function.
   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   (i.e., saturation_fraction must be between 0 and 1, and the default
     is 1, completely saturated)
   t_freezing : in-situ temperature at which seawater freezes.    [ deg C ]
---------------------------------------------------------------------------
***************************************************************************/
double TeosIce::gsw_t_freezing_exact(double sa,double p,double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double df_dt, tf, tfm, tf_old, f, return_value;

	int polynomial=1;
	/** The initial value of t_freezing_exact (for air-free seawater) */
	tf = gsw_t_freezing_poly(sa,p,saturation_fraction,polynomial);
	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p) -
	gsw_gibbs_ice(1,0,tf,p);
	/**
	  df_dt here is the initial value of the derivative of the function f whose
	  zero (f = 0) we are finding (see Eqn. (3.33.2) of IOC et al (2010)).
	*/
	tf_old = tf;
	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
	gsw_gibbs_ice(0,0,tf_old,p);
	tf = tf_old - f/df_dt;
	tfm = 0.5*(tf + tf_old);
	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tfm,p) -
	gsw_gibbs_ice(1,0,tfm,p);
	tf = tf_old - f/df_dt;
	tf_old = tf;
	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
	gsw_gibbs_ice(0,0,tf_old,p);
	tf = tf_old - f/df_dt;
	/** Adjust for the effects of dissolved air */
	return_value = tf - saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gtc.gsw_sso));

	return (return_value);
}

/**************************************************************************
==========================================================================
  method: gsw_sa_p_inrange (sa, p)
==========================================================================
   Check for any values that are out of the TEOS-10 range ...
   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
---------------------------------------------------------------------------
***************************************************************************/
int TeosIce::gsw_sa_p_inrange(double sa,double p)
{
	if (p > 10000.0 || sa > 120.0 ||
			(p + sa*71.428571428571402) > 13571.42857142857)
	{
		return (0);
	}

	return (1);
}

/**************************************************************************
==========================================================================
  method: gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
==========================================================================
   Calculates an initial estimate of pt0_ice when it is less than about
   -100 deg C.
   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   pt0_cold_ice_poly  :  initial estimate of potential temperatur
          of very cold ice in dgress C (not K)     [ deg C ]
--------------------------------------------------------------------------
*************************************************************************/
double TeosIce::gsw_pt0_cold_ice_poly(double pot_enthalpy_ice)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  log_abs_theta0, log_h_diff,
			  /**h00 = gsw_enthalpy_ice(-gtc.gsw_t0,0)*/
			  h00 = -6.320202333358860e5,
			  s0 =  1.493103204647916,
			  s1 =  2.372788609320607e-1,
			  s2 = -2.014996002119374e-3,
			  s3 =  2.640600197732682e-6,
			  s4 =  3.134706016844293e-5,
			  s5 =  2.733592344937913e-6,
			  s6 =  4.726828010223258e-8,
			  s7 = -2.735193883189589e-9,
			  s8 = -8.547714991377670e-11;

	log_h_diff = log(pot_enthalpy_ice - h00);

	log_abs_theta0 = s0
	                  + log_h_diff*(s1
	                  + log_h_diff*(s2
	                  + log_h_diff*(s3
							+ log_h_diff*(s4
							+ log_h_diff*(s5
							+ log_h_diff*(s6
							+ log_h_diff*(s7
							+ log_h_diff*s8)))))));

	return (exp(log_abs_theta0) - gtc.gsw_t0);
}
