module canopy_drydep_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        ZK, FCH, TEMPA, PRESSA, &
        RELHUMA, FSUN, PPFD_SUN, PPFD_SHADE, UBAR, &
        SRAD, RA, DEP_IND, DEP_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy gas dry deposition based on Zhang et al. (2003)

! Revision History:
!     Prototype 02/25 by PCC, based on Zhang et al. (2003) algorithms as implemented
!                             in ACCESS (Saylor 2013)
! Citation:
! Zhang, L., Brook, J. R., and Vet, R.: A revised parameterization for gaseous dry
! deposition in air-quality models, Atmos. Chem. Phys., 3, 2067–2082,
! https://doi.org/10.5194/acp-3-2067-2003, 2003.
!
! Saylor, R. D.: The Atmospheric Chemistry and Canopy Exchange Simulation System (ACCESS):
! model description and application to a temperate deciduous forest canopy,
! Atmos. Chem. Phys., 13, 693–715, https://doi.org/10.5194/acp-13-693-2013, 2013.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Feb 2025 P.C. Campbell: Initial Zhang et al. gas dry deposition
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,rs_zhang_gas,EffHenrysLawCoeff,&
            ReactivityParam,rbl,rcl,rml

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: ZK(:)              ! Model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                ! Canopy height (m)
        REAL(RK),    INTENT( IN )       :: FSUN(:)            ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)        ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)      ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TEMPA(:)           ! Ambient Temperature profile in canopy [K]
        REAL(RK),    INTENT( IN )       :: PRESSA(:)          ! Ambient Pressure profile in canopy [mb]
        REAL(RK),    INTENT( IN )       :: RELHUMA(:)         ! Ambient Relative Humidity profile in canopy [%]
        REAL(RK),    INTENT( IN )       :: UBAR(:)            ! Mean above/in-canopy wind speed [m/s]
        REAL(RK),    INTENT( IN )       :: SRAD               ! Incoming solar irradiation top of canopy (W/m^2)
        REAL(RK),    INTENT( IN )       :: RA                 ! Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech, set in constants)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT(:)         ! Output canopy layer gas dry deposition rate for each DEP_IND

! Local Variables
        REAL(RK) ::  PPFD(SIZE(ZK))                           ! PPFD ave sun and shade (umol/m2 s)
        REAL(RK) ::  mdiffl(SIZE(ZK))                         ! Molecular diffusivity for species l based on DEP_IND [cm2/s]
        REAL(RK) ::  rs(SIZE(ZK))                             ! Stomatal resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rb(SIZE(ZK))                             ! Boundary layer resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rc(SIZE(ZK))                             ! Cuticular resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rm(SIZE(ZK))                             ! Mesophyll resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  hstarl                                   ! effective Henry's law coefficient based on DEP_IND (M/atm)
        REAL(RK) ::  f01                                      ! reactivity parameter based on DEP_IND (0-1)

        REAL(RK) :: rnum,rden,rlx,vdlx
        INTEGER i

        PPFD = (PPFD_SUN*FSUN) + (PPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

!! Calculate molecular diffusivity (cm^2/s) and resistances (cm/s) of species l using from DEP_IND
        hstarl  = EffHenrysLawCoeff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        f01     = ReactivityParam(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        do i=1, SIZE(ZK)
            if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           ! above ground level and at/below canopy top
                mdiffl(i)  = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPA(i),PRESSA(i))
                rs(i)      = rs_zhang_gas(mdiffl(i),TEMPA(i),PRESSA(i),PPFD(i),SRAD,RELHUMA(i)) !stomatal resistance (s/cm)
                rb(i)      = rbl(mdiffl(i), UBAR(i)*100.0_rk) !leaf boundary layer resistance (s/cm)
                rc(i)      = rcl(hstarl, f01)                              !leaf cuticular resistance (s/cm)
                rm(i)      = rml(hstarl, f01)                              !leaf mesophyll resistance (s/cm)
                rnum = rc(i) * (rs(i) + rm(i))
                rden = rc(i) + 2.0_rk * (rs(i) + rm(i))
                rlx   = rb(i) + (rnum/rden) + RA                           !surface+boundary+aerodynamic resistances (s/cm)
                vdlx  = 1.0_rk/rlx
                dep_out(i) = vdlx                                          !calculate deposition velocity (cm/s)
            else
                rb(i) = 0.0_rk
                rc(i) = 0.0_rk
                rm(i) = 0.0_rk
                rs(i) = 0.0_rk
                dep_out(i) = 0.0_rk
            endif
        end do

    END SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG


    SUBROUTINE CANOPY_GAS_DRYDEP_SOIL( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        TEMPSOIL, PRESSA, UBAR, SOCAT, SOTYP, DSOIL, STHETA, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,SoilResist,SoilRbg

        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: TEMPSOIL           ! soil temperature in topsoil (K)
        REAL(RK),    INTENT( IN )       :: PRESSA             ! Ambient Pressure just above surface [mb]
        REAL(RK),    INTENT( IN )       :: UBAR               ! Mean wind speed just above surface [m/s]
        INTEGER,     INTENT( IN )       :: SOCAT              ! input soil category datset used
        INTEGER,     INTENT( IN )       :: SOTYP              ! input soil type integer associated with soilcat
        REAL(RK),    INTENT( IN )       :: DSOIL              ! depth of topsoil (cm)
        REAL(RK),    INTENT( IN )       :: STHETA             ! volumetric soil water content in topsoil(m^3/m^3)
        REAL(RK),    INTENT( IN )       :: RA                 ! Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            ! Output soil layer gas dry deposition rate for each DEP_IND

        real(rk)                        :: mdiffl             ! molecular diffusivity (cm^2/s)
        real(rk)                        :: rsoill             ! resistance to diffusion thru soil pore space for chemical species (s/cm)
        real(rk)                        :: rbg                ! ground boundary layer resistance (s/cm)


        mdiffl = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPSOIL,PRESSA)  !Use soil temperature and pressure

        rsoill = SoilResist(mdiffl,SOCAT,SOTYP,DSOIL,STHETA) !Depends on soil type, depth, and moisture

        rbg = SoilRbg(UBAR*100.0_rk) !convert wind to cm/s   !Rbg(ground boundary layer resistance, s/cm)
        !Rbg is invariant to species not layers
        !Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0
        DEP_OUT = 1.0_rk/(rbg+rsoill+RA)                         !deposition velocity to ground surface under canopy or outside of
        !canopy, e.g., barren land (cm/s)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_SOIL

    SUBROUTINE CANOPY_GAS_DRYDEP_SNOW( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: ReactivityParamHNO3, SoilRbg

        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: UBAR               ! Mean wind speed just above surface [m/s]
        REAL(RK),    INTENT( IN )       :: RA                 ! Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            ! Output soil layer gas dry deposition rate for each DEP_IND

        real(rk), parameter             :: ar_0   = 8.0       ! used to scale other species to HNO3 (dimensionless)
        real(rk)                        :: ar_l               ! reactivity denominator relative to HNO3 for each species (dimensionless)
        real(rk), parameter             :: rsnow0 = 100.0     ! resistance to deposition to snow (cm/s) based on Helmig et al.
        real(rk)                        :: rsnowl             ! resistance to diffusion thru snow space for chemical species (s/cm)
        real(rk)                        :: rbg                ! ground boundary layer resistance (s/cm)

        ar_l = ReactivityParamHNO3(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        rsnowl = rsnow0 * (ar_0/ar_l)  !Based on CMAQv5.3.1 formulation scaled to reactivity relative to HNO3

        rbg = SoilRbg(UBAR*100.0_rk) !convert wind to cm/s   !Rbg(ground boundary layer resistance, s/cm)
        !Rbg is invariant to species not layers
        !Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0

        DEP_OUT = 1.0_rk/(rbg+rsnowl+RA)                               !deposition velocity to ground surface under canopy or outside of
        !canopy, e.g., barren land (cm/s) when snow cover is present

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_SNOW

    SUBROUTINE CANOPY_GAS_DRYDEP_URBAN( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, TEMP, GAMMA_BUILD, &
        RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, pi, rgasuniv                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolarMassGas, SoilRbg

        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: UBAR               ! Mean wind speed just above surface [m/s]
        REAL(RK),    INTENT( IN )       :: TEMP               ! Mean temperature just above surface [K]
        REAL(RK),    INTENT( IN )       :: GAMMA_BUILD        ! Reaction probability with building type (dimensionless)
        ! Default NL is average of range in gamma from as low as 10−8 for
        ! glass and metal to 10−4 for activated carbon and brick.
        ! =5.0D-5.  Reference (Shen and Gao, 2018;
        ! https://doi.org/10.1016/j.buildenv.2018.02.046)
        REAL(RK),    INTENT( IN )       :: RA                 ! Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            ! Output soil layer gas dry deposition rate for each DEP_IND

        real(rk)                        :: mmg_l              ! molar mass for each gas species (kg/mol)
        real(rk)                        :: cave_l             ! Maxwell-Boltzmann average speed of gas distribution (m/s)
        real(rk)                        :: rurbanl            ! resistance to diffusion thru snow space for chemical species (s/m)
        real(rk)                        :: rbg                ! ground boundary layer resistance (s/cm)

        !Get molar mass for each gas species (kg/mol)
        mmg_l = MolarMassGas(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        !Based on Maxwell-Boltzmann distribution, average speed of gas distribution (cm/s)
        cave_l = sqrt((8.0_rk*rgasuniv*TEMP)/(pi*mmg_l))*100.0_rk

        !Based on Shen and Gao (2018), Eq. (2):  https://doi.org/10.1016/j.buildenv.2018.02.046)
        rurbanl = 4.0_rk/(GAMMA_BUILD*cave_l) !already converted in to units of s/cm from cave_l

        rbg = SoilRbg(UBAR*100.0_rk) !convert wind to cm/s   !Rbg(ground boundary layer resistance, s/cm)
        !Rbg is invariant to species not layers
        !Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0

        !deposition velocity to urban surfaces (cm/s)
        DEP_OUT = 1.0_rk/(rbg+rurbanl+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_URBAN

    SUBROUTINE CANOPY_GAS_DRYDEP_WATER( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, TEMP2, QV2, &
        USTAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, cpd, lv0, dlvdt, stdtemp      !constants for canopy models
        use canopy_utils_mod,  ONLY: EffHenrysLawCoeff, LeBasMVGas, WaterRbw

        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT     ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT     ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: TEMP2               ! Mean temperature just above surface [K]
        REAL(RK),    INTENT( IN )       :: QV2                 ! Mean mixing ratio just above surface [kg/kg]
        REAL(RK),    INTENT( IN )       :: USTAR               ! Friction velocity at surface (m/s)
        REAL(RK),    INTENT( IN )       :: RA                 ! Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            ! Output soil layer gas dry deposition rate for each DEP_IND

        real(rk)                        :: hstarl             ! effective Henry's law coefficient based on DEP_IND (M/atm)
        real(rk)                        :: ctemp2             ! Mean temperature just above surface [C]
        real(rk)                        :: lv                 ! latent heat of vaporization [ J/kg ]
        real(rk)                        :: cp_air             ! specific heat of moist air [J/kg-K]
        real(rk)                        :: tw                 ! wet bulb temperature [K]

        real(rk)                        :: lebas_l            ! Le Bas molar volumes are from the Schroeder additive method [cm3/mol ]
        real(rk)                        :: dw                 ! diffusivity of water
        real(rk)                        :: dw25               ! diffusivity of water at 298.15 K
        real(rk)                        :: kvisw              ! kinematic viscosity of water [cm^2/s]
        real(rk)                        :: scw_pr_23          ! (scw/pr)**2/3
        real(rk)                        :: rbw                ! water boundary layer resistance (s/cm)
        real(rk)                        :: rwaterl            ! water surface resistance (s/cm)

        real(rk), Parameter             :: pr         = 0.709 ! Prandtl Number [dim'less]
        real(rk), Parameter             :: rt25inK    = 1.0_rk/(stdtemp + 25.0_rk) ! 298.15K = 25C
        real(rk), Parameter             :: twothirds  = 2.0_rk / 3.0_rk
        real(rk), Parameter             :: d3         = 1.38564e-2 ! scaling parameter used to estimate the friction velocity in surface waters from
        ! the atmospheric friction velocity to a value following Slinn et al. (1978)
        ! and Fairall et al. (2007)

        hstarl  = EffHenrysLawCoeff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        ! Calculate the water surface film temperature: wet bulb temperature.
        ! Wet bulb temperature based on eqn in Fritschen and Gay (1979).

        ctemp2 = TEMP2 - stdtemp                                      ! [C]
        lv     = lv0 - dlvdt * ctemp2                                 ! [J/kg]
        cp_air = cpd * ( 1.0_rk + 0.84_rk * QV2 )                     ! [J/kg-K]
        tw     = ( ( 4.71e4 * cp_air / lv ) - 0.870_rk ) + stdtemp    ! [K]

        ! Make Henry's Law constant non-dimensional.

        hstarl  = hstarl * 0.08205_rk * tw

        !Get Le Bas molar volumes for each gas species [[cm3/mol ]
        lebas_l = LeBasMVGas(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        ! from Hayduk and Laudie
        dw25 = 13.26e-5 / ( 0.8904_rk**1.14_rk * lebas_l**0.589_rk )
        kvisw = 0.017_rk * EXP( -0.025_rk * ( tw - stdtemp ) )
        dw    = dw25 * ( tw * rt25inK ) * ( 0.009025_rk / kvisw )
        scw_pr_23 = ( ( kvisw / dw ) / pr ) ** twothirds

        rwaterl = scw_pr_23 / ( hstarl * d3 * USTAR*100.0_rk ) ! Resistance to water surface, (s/cm)

        rbw = WaterRbw(d3*USTAR*100.0_rk) !convert ustar to cm/s and scale by d3  !Rbw(water boundary layer resistance, s/cm)
        !Rbw is invariant to species not layers

        !deposition velocity to water surfaces (cm/s)
        DEP_OUT = 1.0_rk/(rbw+rwaterl+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_WATER

end module canopy_drydep_mod
