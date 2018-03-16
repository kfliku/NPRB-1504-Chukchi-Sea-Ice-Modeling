/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Central Channel simulation
*/

/* COAWST stuff */
#undef ROMS_MODEL
#undef WRF_MODEL
#undef MCT_LIB
#undef ATM2OCN_FLUXES  /* not sure about this with ice */
#undef NO_LBC_ATT

/* ROMS stuff */
#undef NO_HIS
#undef HDF5
#undef DEFLATE
#define PERFECT_RESTART
#define USE_NETCDF4
#define MPI
#define DOUBLE_PRECISION
#define NONLINEAR
#define POWER_LAW
#define K_GSCHEME
#define VAR_RHO_2D
#define VISC_GRID
#define DIAGNOSTICS_TS
/* general */

#define SOLVE3D
#ifdef SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif

#define ASSUMED_SHAPE

/* boundary condition */

#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC
#define ANA_TOBC
#define FSOBC_REDUCED

/* mixing */

# undef LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
# endif

# undef GLS_MIXING
# define MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#endif


/* ice */

#ifdef SOLVE3D
# undef CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  define SNOW_FROM_RAIN
# endif

# define  ICE_MODEL
# define NO_SNOW
# ifdef ICE_MODEL
#  define ANA_ICE
#  define ANA_HIOBC
#  define ANA_AIOBC
#  undef SNOWFALL
#  define  ICE_THERMO
#  define  ICE_MK
#  define   ICE_MOMENTUM
#  define   ICE_MOM_BULK
#  define   ICE_EVP
#  define  ICE_STRENGTH_QUAD
#  define  ICE_ADVECT   /* Note that we need these two for the */
#  define  ICE_SMOLAR   /* timestepping to work correctly.     */
#  define  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  undef ICE_CONVSNOW
#  undef  MELT_PONDS
#  define ICE_DIAGS
#  define ICE_I_O
# endif
#endif

/* output stuff */

#undef NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#undef AVERAGES

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#define UV_ADV
#define UV_COR
#define UV_U3HADVECTION
#define UV_C4VADVECTION
#undef UV_SADVECTION

#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif

#define UV_QDRAG
#define UV_VIS2
#define TS_DIF2
#define MIX_GEO_TS
#define MIX_S_UV

/* surface forcing */

#ifdef SOLVE3D
# ifndef ATM2OCN_FLUXES
#  define CORE_FORCING
#  define BULK_FLUXES
#  define CCSM_FLUXES
#  undef ARCTIC_MERRA_HACK
# endif
# if defined BULK_FLUXES
#  undef LONGWAVE_OUT
#  define SOLAR_SOURCE
#  define EMINUSP
#  define ANA_LRFLUX
#  define ANA_SRFLUX
#  undef ANA_ALBEDO
#  define ALBEDO
#  undef ANA_SNOW
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92  /* for ice */
#  undef ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE  /* for both */
#  define LONGWAVE
#  define ANA_PAIR
#  define ANA_TAIR
#  define ANA_HUMIDITY
#  define ANA_WINDS
#  define ANA_CLOUD
#  define ANA_RAIN
# endif
#endif


/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

