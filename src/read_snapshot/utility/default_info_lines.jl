"""
    default_info_lines

Storage Array as a fallback if no `INFO` block is present.
"""
const default_info_lines = [ InfoLine("POS",  Float32, 3,  [1, 1, 1, 1, 1, 1]), # Positions (internal units)
                             InfoLine("VEL",  Float32, 3,  [1, 1, 1, 1, 1, 1]), # Velocities (internal units - not v_comoving but v_com*sqrt(a))
                             InfoLine("ID",   UInt32,  1,  [1, 1, 1, 1, 1, 1]), # Particle ID
                             InfoLine("MASS", Float32, 1,  [1, 1, 1, 1, 1, 1]), # Mass of Particle (internal units)
                             InfoLine("U",    Float32, 1,  [1, 0, 0, 0, 0, 0]), # Internal Energy of gas particles (internal units)
                             InfoLine("RHO",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # Density (internal units)
                             InfoLine("NE",   Float32, 1,  [1, 0, 0, 0, 0, 0]), # Number density of free electrons
                             InfoLine("NH",   Float32, 1,  [1, 0, 0, 0, 0, 0]), # Number density of neutral hydrogen
                             InfoLine("HSML", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SFR",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # star formation rate (internal units)
                             InfoLine("AGE",  Float32, 1,  [0, 0, 0, 0, 1, 0]), # Expansion factor at which star (or BH) is born
                             InfoLine("BHMA", Float32, 1,  [0, 0, 0, 0, 0, 1]), # True blackhole mass (MASS contains the dynamical mass !!!)
                             InfoLine("BHMD", Float32, 1,  [0, 0, 0, 0, 0, 1]), # blackhole accretion rate

                             InfoLine("POT",  Float32, 1,  [1, 1, 1, 1, 1, 1]), # Potential
                             InfoLine("TSTP", Float32, 1,  [1, 1, 1, 1, 1, 1]), # timestep
                             InfoLine("ACCE", Float32, 3,  [1, 1, 1, 1, 1, 1]), # acceleration
                              
                             InfoLine("BFLD", Float32, 3,  [1, 0, 0, 0, 0, 0]), # magnetic field
                             InfoLine("BFSM", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Smoothed magnetic field
                             InfoLine("DBDT", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Rate of change of magnetic field
                             InfoLine("DIVB", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Divergence of magnetic field
                             InfoLine("RDIB", Float32, 1,  [1, 0, 0, 0, 0, 0]), # relative Divergence of magnetic field
                             InfoLine("ROTB", Float32, 3,  [1, 0, 0, 0, 0, 0]), # rot of magnetic field
                             InfoLine("SRTB", Float32, 3,  [1, 0, 0, 0, 0, 0]), # smoothed rot of magnetic field 

                             InfoLine("Zs",   Float32, 8,  [1, 0, 0, 0, 1, 0]), # or 11?  # Mass of metals
                             InfoLine("Z",    Float32, 1,  [1, 0, 0, 0, 1, 0]), # Metallicity
                             InfoLine("iM",   Float32, 1,  [0, 0, 0, 0, 1, 0]), # initial mass of star particle
                             InfoLine("HSMS", Float32, 1,  [0, 0, 0, 0, 1, 0]), # Chemical spreading (e.g. Smoothing) length of stars

                             InfoLine("TIPS", Float32, 9,  [1, 1, 1, 1, 1, 1]), # 3x3 configuration-space tidal tensor that is driving the GDE
                             InfoLine("DIPS", Float32, 36, [1, 1, 1, 1, 1, 1]), # full 6D phase-space distortion tensor from GDE integration
                             InfoLine("CACO", Float32, 1,  [1, 1, 1, 1, 1, 1]), # caustic counter
                             InfoLine("FLDE", Float32, 1,  [1, 1, 1, 1, 1, 1]), # physical NON-CUTOFF corrected stream determinant = 1.0/normed stream density * 1.0/initial stream density
                             InfoLine("STDE", Float32, 1,  [1, 1, 1, 1, 1, 1]), # physical CUTOFF corrected stream density = normed stream density * initial stream density 
                             InfoLine("PSDE", Float32, 1,  [1, 1, 1, 1, 1, 1]), # determinant of phase-space distortion tensor -> should be 1 due to Liouville theorem OR pre shock density
                             
                             InfoLine("ANRA", Float32, 3,  [1, 1, 1, 1, 1, 1]), # annihilation radiation: time integrated stream density in physical units
                             InfoLine("LACA", Float32, 20, [1, 1, 1, 1, 1, 1]), # extensive information on the last caustic the particle has passed
                             InfoLine("SHIN", Float32, 3,  [1, 1, 1, 1, 1, 1]), 
                             InfoLine("SHOR", Float32, 9,  [1, 1, 1, 1, 1, 1]), # initial orientation of the CDM sheet where the particle started
                             InfoLine("INDE", Float32, 1,  [1, 1, 1, 1, 1, 1]), # initial stream density in physical units
                             InfoLine("HII",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # ionized hydrogen abundance
                             InfoLine("HeI",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # neutral Helium
                             InfoLine("HeII", Float32, 1,  [1, 0, 0, 0, 0, 0]), # ionized Heluum
                             InfoLine("He3",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # double ionised Helium 
                             InfoLine("H2I",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # H2 molecule
                             InfoLine("H2II", Float32, 1,  [1, 0, 0, 0, 0, 0]), # ionised H2 molecule
                             InfoLine("HM" ,  Float32, 1,  [1, 0, 0, 0, 0, 0]), # H minus
                             InfoLine("HD" ,  Float32, 1,  [1, 0, 0, 0, 0, 0]), # HD
                             InfoLine("DI" ,  Float32, 1,  [1, 0, 0, 0, 0, 0]), # deuterium
                             InfoLine("DII",  Float32, 1,  [1, 0, 0, 0, 0, 0]), # deuteriumII
                             InfoLine("HeHp", Float32, 1,  [1, 0, 0, 0, 0, 0]), # HeH+
                             #InfoLine("ARHO", Float32, 1,  [0, 0, 0, 0, 0, 1]),
                             #InfoLine("HSMB", Float32, 1,  [0, 0, 0, 0, 0, 1]),
                             #InfoLine("HSBP", Float32, 1,  [0, 0, 1, 1, 0, 0]),
                             # 
                             InfoLine("ACRS", Float32, 1,  [0, 0, 0, 0, 1, 0]), # accreation length for star particles
                             InfoLine("BHPC", Float32, 1,  [0, 0, 0, 0, 0, 1]), # progenitor count blackholes
                             InfoLine("ACRB", Float32, 1,  [0, 0, 0, 0, 0, 1]), # Black hole sphere of influence (e.g. smoothing) of Black Holes
                             # IPOT # gravitational potential
                             InfoLine("PHID", Float32, 1,  [1, 0, 0, 0, 0, 0]), # dPotential/dt
                             InfoLine("ABVC", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Artificial bulk viscosity constant for gas particle
                             # ACVC # artificial conductivity of particle
                             # AMDC # artificial magnetic dissipation of particle
                             # VTRB # Local Turbulence estimation from SF
                             # LTRB # Local Coherence lenght estimation from SF
                             # ADYN # alfa^2 Dynamo
                             # EDYN # Eta in Dynamo
                             # PHI # divBcleaning fuction of particle
                             # XPHI # Cold fraction in SF
                             # GPHI # divBcleaning fuction of particle
                             # ROTB # rot of magnetic field
                             # SRTB # smoothed rot of magnetic field
                             # EULA # EulerPotentialA
                             # EULB # EulerPotentialB
                             # COOR # Cooling rate
                             # CONR # current heating/cooling due to thermal conduction
                             # DENN # density normalization factor
                             # EGYP # EnergyReservoirForFeeback
                             # EGYC # EnergyReservoirForColdPhase
                             # CRC0
                             # CRQ0
                             # CRP0
                             # CRE0
                             # CRn0
                             # CRco
                             # CRdi
                             # BHMA
                             # ACRB
                             # BHMD
                             # MHPC
                             # BHMD
                             # BHMI
                             # BHMR
                             # VSIG
                             InfoLine("SHNR", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("MACH", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHSP", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHCP", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHRU", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHFL", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHVU", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             InfoLine("SHVD", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Smoothing length of gas particles
                             # SHRH
                             # SHPU
                             # SHPD
                             # SHVU
                             # SHVD
                             # MALF
                             # SHAU
                             # SHAD
                             # SHOB
                             # SHCP
                             # SHNR
                             # DTED
                             # PSCS
                             # PSDE
                             # PSEN
                             # PSXC
                             # DJMP
                             # EJMP
                             # CRDE
                             # TIPS
                             # DIPS
                             # CACO
                             # FLDE
                             # STDE
                             # SOMA
                             # PSDE
                             # ANRA
                             # LACA
                             # SHOR
                             # INDE
                             # TEMP
                             # XNUC
                             # P
                             # RADG
                             # RADA
                             # ET
                             # PTSU
                             # DMNB
                             # NTSC
                             # SHSM
                             # SRHO
                             # SVEL
                             # DMHS
                             # DMDE
                             # DMVD
                             # Zs
                             # CII
                             # CIa
                             # CAGB
                             # ZAge
                             # ZAlv
                             # iM
                             # HOTT
                             # TEMP
                             # TRCK
                             # ZsMT
                             # ZSMT
                             # MHOT
                             # MCLD
                             # MMOL
                             # EHOT
                             # MSF
                             # MFST
                             # NMF
                             # EOUT
                             # CLCK
                             # Egy0
                             # TDYN
                             # EREC
                             # GRAD
                             # TCOO
                             # EKRC
                             # TSMP
                             # CRpN
                             # CReN
                             # CRpS
                             # CReS
                             # CRpC
                             # CReC
                             # CRpP
                             # CReP
                             # RHOO
                             # SYNE
                             # AGSH
                             # AGSD
                             # AGSZ
                             # AGSO
                             # AGSC
                             # AGSN
                             # MGPH
                             # MGGP
                             # MGAC
                             # GRDU
                             # GRDP
                             # MGVX
                             # MGYY
                             # MGVZ
                             # MGRH
                             # MGU
                             # QP
                             # QDP
                             # RHOS
                             # SMST
                             # NEIS
                             InfoLine("VRMS", Float32, 1,  [1, 0, 0, 0, 0, 0]), # RMS velocity within kernel for gas particles
                             InfoLine("VBLK", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Mean bulk velocity within smoothing kernel
                             # VTAN # turbulent velocity around mean, tangential part
                             # VRAD # turbulent velocity around mean, radial part
                             # VROT # Velocity Curl
                             # VORT # Vorticity
                             InfoLine("TNGB", UInt32,  1,  [1, 0, 0, 0, 1, 1]), # true number of neighbors (gas particles)
                             # NGB
                             InfoLine("CLDX", Float32, 1,  [1, 0, 0, 0, 0, 0]), # cold fraction from Springel & Hernquist model
                             InfoLine("HOTT", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Hot gas temperature
                             InfoLine("TEMP", Float32, 1,  [1, 0, 0, 0, 0, 0]), # mean gas temperature (use this) [K]
                             InfoLine("DETI", Float32, 1,  [1, 0, 0, 0, 0, 0]),
                             # ENDT # DtEntropy
                             # STRD # Stressdiag
                             # STRO # Stressoffdiag
                             # STRB # Stressbulk
                             # SHCO # Shearcoeff
                             # DPP # Reacceleration Coefficient
                             # VDIV # Divergence of Vel

                             # MHIX # Mainhalo
                             
                             # subfind
                             InfoLine("GLEN", UInt32,  1,  [1, 0, 0, 0, 0, 0]), # Total number of particles in the FoF
                             InfoLine("GOFF", UInt32,  1,  [1, 0, 0, 0, 0, 0]), # Offset of the corresponding ID particles
                             InfoLine("MTOT", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Total mass of Halo [10^10 h-1 M_☉]
                             InfoLine("GPOS", Float32, 3,  [1, 0, 0, 0, 0, 0]), # Position [h-1 kpc]
                             InfoLine("MTOP", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_vir [10^10 h-1 M_☉]
                             InfoLine("RTOP", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_vir [h-1 kpc]
                             InfoLine("MVIR", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_vir [10^10 h-1 M_☉]
                             InfoLine("RVIR", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_vir [h-1 kpc]
                             InfoLine("MMEA", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_200 mean [10^10 h-1 M_☉]
                             InfoLine("RMEA", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_200 wrt the mean density of the Universe [h-1 kpc]
                             InfoLine("M200", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_200 mean [10^10 h-1 M_☉]
                             InfoLine("R200", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_200 wrt the mean density of the Universe [h-1  kpc]
                             InfoLine("MCRI", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_200 crit [10^10 h-1 M_☉]
                             InfoLine("RCRI", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_200 wrt the critical density of the Universe [h-1 kpc]
                             InfoLine("M500", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_500 mean [10^10 h-1 M_☉]
                             InfoLine("R500", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_500 mean [h-1 kpc]
                             InfoLine("M5CC", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_500 crit [h-1 kpc]
                             InfoLine("R5CC", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_500 crit [10^10 h-1 M_☉]
                             InfoLine("M25K", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass inside r_2500 crit [10^10 h-1 M_☉]
                             InfoLine("R25K", Float32, 1,  [1, 0, 0, 0, 0, 0]), # r_2500 crit [h-1 kpc]
                             InfoLine("MGAS", Float32, 6,  [1, 0, 0, 0, 0, 0]), # Mass of gas in the six radii RVIR/RMEA/R500/RCRI/R5CC/R25K [10^10 h-1 M_☉]
                             InfoLine("MSTR", Float32, 6,  [1, 0, 0, 0, 0, 0]), # Mass of stars in the six radii RVIR/RMEA/R500/RCRI/R5CC/R25K [10^10 h-1 M_☉]
                             InfoLine("TGAS", Float32, 6,  [1, 0, 0, 0, 0, 0]), # Gas temperature in the six radii RVIR/RMEA/R500/RCRI/R5CC/R25K [KeV]
                             InfoLine("LGAS", Float32, 6,  [1, 0, 0, 0, 0, 0]), # X-ray bolometric luminosity in the six radii RVIR/RMEA/R500/RCRI/R5CC/R25K [10^44 erg/s]
                             InfoLine("YGAS", Float32, 6,  [1, 0, 0, 0, 0, 0]), # Ygas, unitless, in the six radii RVIR/RMEA/R500/RCRI/R5CC/R25K
                             InfoLine("NCON", UInt32,  1,  [1, 0, 0, 0, 0, 0]), # Number of contaminting particles in the FoF
                             InfoLine("MCON", Float32, 1,  [1, 0, 0, 0, 0, 0]), # Mass of contaminating particles in the FoF
                             InfoLine("BGPO", Float32, 3,  [0, 0, 0, 1, 0, 0]), # Position of Big Groups (M>10^14)
                             InfoLine("BGMA", Float32, 1,  [0, 0, 0, 1, 0, 0]), # Mass of Big Groups (M>10^14)
                             InfoLine("BGRA", Float32, 1,  [0, 0, 0, 1, 0, 0]), # Radius of Big Groups (M>10^14)
                             InfoLine("NSUB", UInt32,  1,  [1, 0, 0, 0, 0, 0]), # Number of SubHalos in FoF-group
                             InfoLine("FSUB", UInt32,  1,  [1, 0, 0, 0, 0, 0]), # Index of first SubHalo in FoF-group
                             InfoLine("SLEN", UInt32,  1,  [0, 1, 0, 0, 0, 0]), # Number of particles in the SubHalo
                             InfoLine("SOFF", UInt32,  1,  [0, 1, 0, 0, 0, 0]), # Offset for the particles IDs in the PID array
                             InfoLine("SSUB", UInt32,  1,  [0, 1, 0, 0, 0, 0]), # Sub parent halo index (one level above of the SubHalo). Relative index, total index=fsub+ssub
                             InfoLine("MSUB", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Total mass of the SubHalo
                             InfoLine("SPOS", Float32, 3,  [0, 1, 0, 0, 0, 0]), # Positions of SubHalos
                             InfoLine("SVEL", Float32, 3,  [0, 1, 0, 0, 0, 0]), # Velocities of SubHalos (physical units!)
                             InfoLine("SCM",  Float32, 3,  [0, 1, 0, 0, 0, 0]), # Position of the centre of Mass
                             InfoLine("SPIN", Float32, 3,  [0, 1, 0, 0, 0, 0]), # SubHalo spin vector
                             InfoLine("DSUB", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Subhalo internal velocity dispersion
                             InfoLine("VMAX", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Maximum velocity
                             InfoLine("RMAX", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Radius corresponding to the maximum velocity
                             InfoLine("RHMS", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Half-mass radius of SubHalos (DM)
                             InfoLine("SHMR", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Stellar half-mass radius
                             InfoLine("REFF", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Stellar half-mass radius (already contains the scale factor!)
                             InfoLine("MBID", UInt64,  1,  [0, 1, 0, 0, 0, 0]), # ID of the Most bound particle
                             InfoLine("GRNR", UInt32,  1,  [0, 1, 0, 0, 0, 0]), # Parent FoF index
                             InfoLine("SMST", Float32, 6,  [0, 1, 0, 0, 0, 0]), # SubHalo mass table (for each particle type)
                             InfoLine("SLUM", Float32, 12, [0, 1, 0, 0, 0, 0]), # Rest-frame magnitudes in different bands(without dust correction)
                             InfoLine("SLAT", Float32, 12, [0, 1, 0, 0, 0, 0]), # Rest frame magnitudes in different bands(including dust correction)
                             InfoLine("SLOB", Float32, 12, [0, 1, 0, 0, 0, 0]), # Observer frame magnitudes (including dust correction)
                             InfoLine("DUST", Float32, 11, [0, 1, 0, 0, 0, 0]), # Dust profiles in 11 bins (should double check the units)
                             InfoLine("SAGE", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Mean stellar age (in expansion factor parameter)
                             InfoLine("SZ",   Float32, 1,  [0, 1, 0, 0, 0, 0]), # Mean metallicity
                             InfoLine("SSFR", Float32, 1,  [0, 1, 0, 0, 0, 0]), # Sfr [msun/yr]
                             InfoLine("SANG", Float32, 1,  [0, 1, 0, 0, 0, 0]), # SubHalo specific angular momentum in physical units [kpc km/s]
                             InfoLine("PID ", UInt64,  1,  [0, 0, 1, 0, 0, 0]), # Arrays of IDs of all the particles in the FoFs

                             InfoLine("SLg",  Float32, 1,  [0, 0, 0, 0, 1, 0]),
                             InfoLine("SUB",  Float32, 1,  [0, 0, 0, 0, 1, 0]),
                             InfoLine("IDU",  UInt32,  1,  [1, 1, 1, 1, 1, 1])
                             ]