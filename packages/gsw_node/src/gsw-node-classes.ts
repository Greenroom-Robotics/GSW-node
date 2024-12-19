export declare class TeosBase {
  public constructor();

  public gsw_z_from_p(
    p: number,
    lat: number,
    geo_strf_dyn_height: number,
    sea_surface_geopotential: number,
  ): number;

  public gsw_sp_from_c(c: number, t: number, p: number): number;

  public gsw_sa_from_sp(
    sp: number,
    p: number,
    lon: number,
    lat: number,
  ): number;

  public gsw_ct_from_t(sa: number, t: number, p: number): number;

  public gsw_depth_from_z(z: number): number;
}

export declare class TeosIce extends TeosBase {
  public constructor();

  public gsw_cp_ice(t: number, p: number): number;
}

export declare class TeosSea extends TeosBase {
  public constructor();

  public gsw_c_from_sp(sp: number, t: number, p: number): number;

  public gsw_sound_speed(sa: number, ct: number, p: number): number;
}
