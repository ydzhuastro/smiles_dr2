README file for the Systematic Mid-Infrared Instrument Legacy Extragalactic Survey (SMILES)
MAST webpage: https://archive.stsci.edu/hlsp/smiles/
Refer to this HLSP with DOI: http://dx.doi.org/10.17909/et3f-zd57

## Contributor
Stacey Alberts, salberts@arizona.edu, Steward Observatory, University of Arizona
Yongda Zhu, yongdaz@arizona.edu, Steward Observatory, University of Arizona

## INTRODUCTION

The Systematic Mid-Infrared Instrument Legacy Extragalactic Survey (SMILES) (PID 1207) is a Cycle 1 imaging program that obtained MIRI imaging in eight photometric bands (F560W, F770W, F1000W, F1280W, F1500W, F1800W, F2100W, F2550W) plus NIRSpec spectroscopic follow-up using the Multi-Shutter Array (MSA).  SMILES covers an area of 34.5 arcmin2 in the GOODS-S/HUDF field.  MIRI imaging was obtained via ~650-2200 s exposures, reaching 5sig point source sensitivities of 0.2-17 uJy (25.7 - 20.8 AB) in 65% encircled energy apertures (see Alberts et al., 2024 for details).  Data release 1 (DR1) includes science-ready MIRI mosaics as well as a catalog contains 3,096 sources with photometric measurements in all bands. Data Release 2 (DR2) presents medium-resolution JWST/NIRSpec spectroscopy of 166 objects spanning 0 < z < 7.5, targeting a mix of star-forming galaxies, quiescent systems, and AGN with a focus on galaxies at cosmic noon (z ∼ 1–3). Observations were obtained with the G140M/F100LP and G235M/F170LP grating/filter pairs. DR2 includes calibrated 2D and 1D spectra, redshift measurements, emission-line catalogs from GELATO and pPXF, and ancillary SED-based properties (see Zhu et al. 2025 for details). 

==========================================================================================

## Data Products

%%%%%%%%%%%%%%%%%%%%%%%%%%
Data Release 2 (August 2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%

SMILES DR2 presents medium-resolution JWST/NIRSpec spectroscopy of 166 galaxies spanning 0 < z < 7.5 in the GOODS-South field, observed with the Micro-Shutter Array (MSA) using the G140M/F100LP and G235M/F170LP grating/filter pairs. The sample includes star-forming galaxies, quiescent systems, and AGN, with an emphasis on galaxies at cosmic noon (z ~ 1–3). DR2 includes calibrated 2D and 1D spectra, a spectroscopic redshift catalog, emission line measurements, and stellar population properties derived from SED fitting. 

The NIRSpec files follow the naming convention:

hlsp_smiles_jwst_nirspec_goodss<NNNNNNNN>_<suffix>_v1.0_<type>.fits

where:

<NNNNNNNN>     8-digit NIRCam source ID. Omitted for the catalog file.

<suffix>       Describes the spectral configuration:
               - For 2D spectra: <filter>_<grating>_nod<n>, where:
                   - <filter>  = f100lp or f170lp
                   - <grating> = g140m or g235m
                   - <n>       = 2 or 3, indicating the nodding strategy used
               - For 1D spectra: 'multi', indicating that the 1D spectrum combines both gratings
               - For the catalog: also 'multi'

<type>         Data product type:
               - 'spec2d'   : Rectified 2D spectrum for one filter/grating/nod configuration
               - 'spec1d'   : Stitched 1D spectrum (with slitloss correction) across G140M and G235M based on appropriote nodding configuration
               - 'catalog'  : Redshift, emission-line, and SED catalog for all sources

---------

# 2D Spectra

Each file contains the nod-specific, calibrated 2D spectrum with the following structure:

No.   Name       	Type          	Description
0     PRIMARY    	PrimaryHDU    	Metadata and HLSP keywords
1     SCI        	ImageHDU      	2D flux image [MJy]
2     ERR        	ImageHDU      	2D uncertainty image [MJy]
3     EXTRACTED  	BinTableHDU   	1D spectrum from x1d extraction

The EXTRACTED extension contains:

-----------------------
Description of columns:
-----------------------

'WAVELENGTH'     Wavelength in microns (μm)  
'FLUX'           1D extracted flux density in Jy  
'FLUX_ERROR'     Uncertainty on the extracted flux in Jy  

---------

# 1D Spectra

These files contain stitched 1D spectra combining G140M and G235M observations. The extraction uses nod2 spectra for extended sources and nod3 for compact sources (as classified in the DR2 catalog). Each file includes both the original and slitloss-corrected versions of the spectrum (the latter for reference only — see Zhu et al. 2025 for details).

No.   Name           	Type          	Description
0     PRIMARY        	PrimaryHDU    	Metadata and extraction info
1     SPEC           	BinTableHDU   	Stitched 1D spectrum
2     SPEC_SLITLOSS  	BinTableHDU   	Slitloss-corrected 1D spectrum (for reference use only)

Each table includes:

-----------------------
Description of columns:
-----------------------

'WAVELENGTH'     Wavelength in Angstrom  
'FLUX'           Flux density in erg/cm^2/s/Angstrom  
'FLUX_ERROR'     1-sigma uncertainty in flux (erg/cm^2/s/Angstrom)  

---------

# Spectroscopic Catalog

The SMILES DR2 spectroscopic catalog provides derived quantities for 166 galaxies observed with JWST/NIRSpec, including spectroscopic redshifts, emission- and absorption-line properties, and stellar population parameters from SED fitting. Each source is assigned a unique NIRCam ID that links the catalog to its associated 1D/2D spectra files and imaging cutouts. Source classifications (e.g., AGN, extended source) are also included to support filtering by scientific use case.

The catalog is a multi-extension FITS file with the following structure:

No.	Name           	Type          	Description
0   PRIMARY        	PrimaryHDU    	Header only
1   Z_SPEC         	BinTableHDU   	Main redshift catalog
2   GELATO_LINES   	BinTableHDU   	Emission-line fits for star-forming and AGN sources
3   PPXF_LINES     	BinTableHDU   	Absorption-line fits for quiescent sources
4   SED_FITS       	BinTableHDU   	Stellar population properties from SED fitting
5   AGN_SED_FITS   	BinTableHDU   	Placeholder for AGN-specific SED fits (to be updated)

EXT 1 contains the following:

-----------------------
Description of columns:
-----------------------

'NIRCAM_ID'	        Unique source ID matched to NIRCam catalog  
'MIRI_ID'	        Source ID in SMILES DR1 MIRI catalog; -1 if no match  
'RA'		        Right Ascension (J2000; degrees)  
'DEC'		        Declination (J2000; degrees)  
'Z_PHOTO'	        Photometric redshift from ancillary data  
'Z_SPEC'	        Best available spectroscopic redshift  
'UNCERTAIN_Z_SPEC'	Flag (1 = redshift uncertain; 0 = reliable)  
'AGN'			    Flag (1 = AGN; 0 = not AGN)  
'EXTENDED'		    Flag (1 = spatially extended; 0 = compact)  
'KEYWORDS'		    Science label and selection tags (e.g., FRESCO-overdensity, asagao_specz_miri)


EXT 2 contains provides Gaussian emission-line fitting results from the GELATO code (Zhu et al. 2025, in prep) for 152 galaxies. Fluxes, rest-frame equivalent widths (REWs), and line-based redshifts are reported for each detected transition, along with 1σ uncertainties. Lines are grouped by physical origin: AGN narrow lines, Balmer and Paschen series, [O III] outflows, and star formation tracers.

-----------------------
Description of columns:
-----------------------

'NIRCAM_ID'               Unique source ID matched to NIRCam catalog  
'SSP_Redshift'            Best-fit redshift from stellar population fitting (float64)  

For each line <LINE>, the following columns are included:

'<LINE>_Flux'             Line flux in erg/cm^2/s  
'<LINE>_Flux_err'         1-sigma uncertainty in flux  
'<LINE>_REW'              Rest-frame equivalent width in Angstrom  
'<LINE>_REW_err'          1-sigma uncertainty in REW  
'<LINE>_Redshift'         Redshift based on the line center  
'<LINE>_Redshift_err'     1-sigma uncertainty in redshift  

Example <LINE> values include:

AGN narrow lines:  
AGN_[NeV]_3426.85, AGN_[NeIII]_3869.86, AGN_[OIII]_5008.24, AGN_[NII]_6585.27, AGN_[SII]_6718.29  

Balmer and Paschen recombination lines:  
Balmer_HI_4862.68, Balmer_HI_Broad_6564.61, Paschen_HI_10052.6  

Outflow components (fitted independently with offset Gaussians):  
Outflow_[OIII]_Outflow_4960.295, Outflow_[OIII]_Outflow_5008.24  

Star-formation lines:  
SF_[OII]_3728.48, SF_[OI]_6302.046, SF_[SIII]_9533.2  

Helium:  
Helium_HeI_10030.48

If a line is undetected or unconstrained in a given object, all related columns will be marked as null.

Note: For quescent galaxies, refer to the pPXF fits in EXT 3.


EXT 3 contains pPXF line fitting for 16 quiescent galaxies:

-----------------------
Description of columns:
-----------------------

'NIRCAM_ID'         Unique source ID matched to NIRCam catalog
'z_ppxf'            Redshift from pPXF fit
'Hb'                H-beta 4862.68 flux (erg/cm^2/s)
'Hb_err'            1-sigma uncertainty on H-beta flux
'Ha'                H-alpha 6564.61 flux (erg/cm^2/s)
'Ha_err'            1-sigma uncertainty on H-alpha flux
'OII_3727'          [O II] 3727.09 flux (erg/cm^2/s)
'OII_3727_err'      1-sigma uncertainty on [O II] 3727.09
'OII_3729'          [O II] 3729.88 flux (erg/cm^2/s)
'OII_3729_err'      1-sigma uncertainty on [O II] 3729.88
'OI_6302'           [O I] 6302.046 flux (erg/cm^2/s)
'OI_6302_err'       1-sigma uncertainty on [O I] 6302.046
'OI_6365'           [O I] 6365.536 flux (erg/cm^2/s)
'OI_6365_err'       1-sigma uncertainty on [O I] 6365.536
'OIII_4960'         [O III] 4960.295 flux (erg/cm^2/s)
'OIII_4960_err'     1-sigma uncertainty on [O III] 4960.295
'OIII_5008'         [O III] 5008.24 flux (erg/cm^2/s)
'OIII_5008_err'     1-sigma uncertainty on [O III] 5008.24
'NII_6549'          [N II] 6549.86 flux (erg/cm^2/s)
'NII_6549_err'      1-sigma uncertainty on [N II] 6549.86
'NII_6585'          [N II] 6585.27 flux (erg/cm^2/s)
'NII_6585_err'      1-sigma uncertainty on [N II] 6585.27
'SII_6718'          [S II] 6718.29 flux (erg/cm^2/s)
'SII_6718_err'      1-sigma uncertainty on [S II] 6718.29
'SII_6732'          [S II] 6732.67 flux (erg/cm^2/s)
'SII_6732_err'      1-sigma uncertainty on [S II] 6732.67



EXT 4 contains the results of spectral energy distribution (SED) fitting for 151 SMILES NIRSpec sources using the JADES NIRCam photometry.

-----------------------
Description of columns:
-----------------------

'NIRCAM_ID'         Unique source ID matched to NIRCam catalog  
'logmass'           Stellar mass in log(Msun)  
'logmass_e'         16th percentile of stellar mass  
'logmass_E'         84th percentile of stellar mass  
'logzsol'           Stellar metallicity in log(Z/Zsun)  
'logzsol_e'         16th percentile of stellar metallicity  
'logzsol_E'         84th percentile of stellar metallicity  
'dust2'             V-band optical depth for diffuse dust  
'dust2_e'           16th percentile of dust2  
'dust2_E'           84th percentile of dust2  
'dust_index'        Power-law slope of the attenuation curve  
'dust_index_e'      16th percentile of dust_index  
'dust_index_E'      84th percentile of dust_index  
'dust1_fraction'    Ratio of birth-cloud to diffuse dust attenuation  
'dust1_fraction_e'  16th percentile of dust1_fraction  
'dust1_fraction_E'  84th percentile of dust1_fraction  
'log_fagn'          AGN luminosity fraction in log scale  
'log_fagn_e'        16th percentile of log_fagn  
'log_fagn_E'        84th percentile of log_fagn  
'log_agn_tau'       Optical depth of AGN torus (log scale)  
'log_agn_tau_e'     16th percentile of log_agn_tau  
'log_agn_tau_E'     84th percentile of log_agn_tau  
'gas_logz'          Nebular gas-phase metallicity in log(Z/Zsun)  
'gas_logz_e'        16th percentile of gas_logz  
'gas_logz_E'        84th percentile of gas_logz  
'duste_qpah'        PAH mass fraction in dust emission model  
'duste_qpah_e'      16th percentile of duste_qpah  
'duste_qpah_E'      84th percentile of duste_qpah  
'duste_umin'        Minimum radiation field in dust emission model  
'duste_umin_e'      16th percentile of duste_umin  
'duste_umin_E'      84th percentile of duste_umin  
'log_duste_gamma'   Log of fraction of dust heated by strong radiation fields  
'log_duste_gamma_e' 16th percentile of log_duste_gamma  
'log_duste_gamma_E' 84th percentile of log_duste_gamma  
'mfrac'             Ratio of surviving stellar mass to formed mass
'SFR'               Star formation rate (Msun/yr)  
'SFR_e'             16th percentile of SFR  
'SFR_E'             84th percentile of SFR

Note: AGN-related properties (log_fagn, log_agn_tau) may be unreliable for some sources. See EXT 5 for specialized AGN SED fits.

EXT 5 contains ... [TODO]



%%%%%%%%%%%%%%%%%%%%%%%%%%
Data Release 1 (June 2024)
%%%%%%%%%%%%%%%%%%%%%%%%%%

The SMILES data products are named according to the following convention:

hlsp_smiles_jwst_<instrument>_goodss_<filter>_<version>_<data file type>.fits

where:

<inst> is the instrument, here 'miri'
<filter> is the filter or 'multi' for the photometric catalog
<version> is the version number ("v1.0").

Data file types:

_drz.fits		Image mosaic
_catalog.fits	Photometric catalog

---------

# Mosaics

MIRI image mosaics are presented with a pixel scale of 0.06" and have the following structure:

No.	Name  		Type  	 	          Format
  0  	PRIMARY   	 PrimaryHDU 		– 	 
  1  	SCI       	 ImageHDU    		float32     basic image 
  2  	ERR       	 ImageHDU    		float32     flux uncertainty
  3  	EXP       	 ImageHDU   		float32     exposure time
  4  	WHT       	 ImageHDU     		float32     weight
  
----------------------

# Photometric catalog

The SMILES DR1 photometric catalog contains 3,096 sources.  Measurements are provided in all bands for sources detected at ≥4sigma in F560W or F770W, with false positives identified and removed using ultra-deep NIRCam imaging from JADES (see Alberts et al., 2024 for details).  Flux densities and their associated uncertainties are provided in 5 circular apertures (CIRC3-7, r = 0.25, 0.3, 0.35, 0.5, 0.6”) and 2.5x scaled Kron apertures. 

Fluxes and uncertainties are in units of nano-Janskys.

The photometric catalog is a fits file with the following structure:

No.	Name  	   		Type  	             
  0  PRIMARY       	  PrimaryHDU	 
  1  FILTERS   	      BinTableHDU       
  2  SIZE      	      BinTableHDU    
  3  CIRC      	      BinTableHDU
  4  KRON       	  BinTableHDU  
  
EXT 1 contains information about the filters used and their corresponding Point Spread Functions (PSFs), including the FWHM and aperture corrections for the circular apertures used.

EXT 2 contains the following:

-----------------------
Description of columns:
-----------------------

'ID'			Observation ID
'RA'			Source RA
'DEC'			Source Dec
'NPIX_DET'		Number of pixels used for detection
'X'				Pixel location in x direction
'Y'				Pixel location in y direction
'XC'			Pixel location of x centroid in detection image
'YC'			Pixel location of x centroid in detection image
'BBOX_XMIN'		The boundaries of the area used for detection
'BBOX_XMAX'		"
'BBOX_YMIN'		"
'BBOX_YMAX'		"
'R_KRON'		Kron radius (2.5x scaled)
'PA'			Position angle of Kron aperture
'Q'				Axial ratio of Kron aperture
'A'				Semi-major axis
'B'				Semi-minor axis

EXT 3 contains the circular aperture photometry:

-----------------------
Description of columns:
-----------------------

'ID'					Observation ID
'RA'					Source RA
'DEC'					Source Dec
<filter>_<aperture>		Flux densities in a given filter for a given circular aperture
<filter>_<aperture>_en	Flux uncertainties in a given filter for a given circular aperture

Flux densities and uncertainties are given in nJy.  The apertures are CIRC3, CIRC4, CIRC5, CIRC6, CIRC7 with radii of r=0.25, 0.3, 0.35, 0.5, 0.6 arcseconds.

EXT 3 contains the Kron aperture photometry:

-----------------------
Description of columns:
-----------------------

'ID'					Observation ID
'RA'					Source RA
'DEC'					Source Dec
<filter>_KRON			Flux densities in a given filter
<filter>_KRON_en		Flux uncertainties in a given filter

Flux densities and uncertainties are given in nJy.  The Kron apertures are 2.5x scaled.
