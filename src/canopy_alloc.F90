
SUBROUTINE canopy_alloc

!-------------------------------------------------------------------------------
! Name:     Allocate Arrays for Canopy Inputs and Outputs
! Purpose:  Allocate arrays for Canopy Inputs and Outputs
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_canopts_mod
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_canvars_mod

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Allocate canopy input variables
!-------------------------------------------------------------------------------

    if(.not.allocated(variables))     allocate(variables(nlat*nlon))
    if(.not.allocated(variables_2d))  allocate(variables_2d(nlon,nlat))
    if (var3d_opt .eq. 1) then
        if(.not.allocated(variables_3d))  allocate(variables_3d(nlon,nlat,var3d_set))
        if(.not.allocated(variables_1d))  allocate(variables_1d(var3d_set))
        if(.not.allocated(variables_can)) allocate(variables_can(nlat*nlon))
        if(.not.allocated(pavdref))       allocate(pavdref(var3d_set))
        if(.not.allocated(levref))        allocate(levref(var3d_set))
        if(.not.allocated(pavd_arr))      allocate(pavd_arr(var3d_set))
        if(.not.allocated(lev_arr))       allocate(lev_arr(var3d_set))
    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Internal Canopy Distribution Variables
!-------------------------------------------------------------------------------

    if(.not.allocated(zhc))                allocate(zhc(modlays))
    if(.not.allocated(fafraczInt))         allocate(fafraczInt(modlays))
    if(.not.allocated(fsun))               allocate(fsun(modlays))
    if(.not.allocated(tleaf_sun))          allocate(tleaf_sun(modlays))
    if(.not.allocated(tleaf_shade))        allocate(tleaf_shade(modlays))
    if(.not.allocated(tleaf_ave))          allocate(tleaf_ave(modlays))
    if(.not.allocated(ppfd_sun))           allocate(ppfd_sun(modlays))
    if(.not.allocated(ppfd_shade))         allocate(ppfd_shade(modlays))
    if(.not.allocated(ppfd_ave))           allocate(ppfd_ave(modlays))
    if(.not.allocated(lad))                allocate(lad(nlat*nlon,modlays))
    if(.not.allocated(lad_3d))             allocate(lad_3d(nlon,nlat,modlays))
    if(.not.allocated(zo_h))               allocate(zo_h(nlat*nlon))
    if(.not.allocated(zo_h_2d))            allocate(zo_h_2d(nlon,nlat))
    if(.not.allocated(d_h))                allocate(d_h(nlat*nlon))
    if(.not.allocated(d_h_2d))             allocate(d_h_2d(nlon,nlat))
    if(.not.allocated(tka))                allocate(tka(nlat*nlon,modlays))
    if(.not.allocated(tka_3d))             allocate(tka_3d(nlon,nlat,modlays))
    if(.not.allocated(pressa))             allocate(pressa(nlat*nlon,modlays))
    if(.not.allocated(pressa_3d))          allocate(pressa_3d(nlon,nlat,modlays))
    if(.not.allocated(relhuma))            allocate(relhuma(nlat*nlon,modlays))
    if(.not.allocated(relhuma_3d))         allocate(relhuma_3d(nlon,nlat,modlays))
    if(.not.allocated(spechuma))           allocate(spechuma(nlat*nlon,modlays))
    if(.not.allocated(spechuma_3d))        allocate(spechuma_3d(nlon,nlat,modlays))

    if (hist_opt .eq. 1) then
        if(.not.allocated(tleaf_sun24_tmp))    allocate(tleaf_sun24_tmp(25,nlat*nlon,modlays))
        if(.not.allocated(tleaf_shade24_tmp))  allocate(tleaf_shade24_tmp(25,nlat*nlon,modlays))
        if(.not.allocated(tleaf_ave24_tmp))    allocate(tleaf_ave24_tmp(25,nlat*nlon,modlays))
        if(.not.allocated(ppfd_sun24_tmp))     allocate(ppfd_sun24_tmp(25,nlat*nlon,modlays))
        if(.not.allocated(ppfd_shade24_tmp))   allocate(ppfd_shade24_tmp(25,nlat*nlon,modlays))
        if(.not.allocated(tleaf_sun240_tmp))   allocate(tleaf_sun240_tmp(241,nlat*nlon,modlays))
        if(.not.allocated(tleaf_shade240_tmp)) allocate(tleaf_shade240_tmp(241,nlat*nlon,modlays))
        if(.not.allocated(tleaf_ave240_tmp))   allocate(tleaf_ave240_tmp(241,nlat*nlon,modlays))
        if(.not.allocated(ppfd_sun240_tmp))    allocate(ppfd_sun240_tmp(241,nlat*nlon,modlays))
        if(.not.allocated(ppfd_shade240_tmp))  allocate(ppfd_shade240_tmp(241,nlat*nlon,modlays))
        if(.not.allocated(tmp2mref_tmp))       allocate(tmp2mref_tmp(25,nlat*nlon))
        if(.not.allocated(ubzref_tmp))         allocate(ubzref_tmp(25,nlat*nlon))
    end if

    if(.not.allocated(tleaf_sun24))        allocate(tleaf_sun24(nlat*nlon,modlays))
    if(.not.allocated(tleaf_shade24))      allocate(tleaf_shade24(nlat*nlon,modlays))
    if(.not.allocated(tleaf_ave24))        allocate(tleaf_ave24(nlat*nlon,modlays))
    if(.not.allocated(ppfd_sun24))         allocate(ppfd_sun24(nlat*nlon,modlays))
    if(.not.allocated(ppfd_shade24))       allocate(ppfd_shade24(nlat*nlon,modlays))
    if(.not.allocated(tleaf_sun240))       allocate(tleaf_sun240(nlat*nlon,modlays))
    if(.not.allocated(tleaf_shade240))     allocate(tleaf_shade240(nlat*nlon,modlays))
    if(.not.allocated(tleaf_ave240))       allocate(tleaf_ave240(nlat*nlon,modlays))
    if(.not.allocated(ppfd_sun240))        allocate(ppfd_sun240(nlat*nlon,modlays))
    if(.not.allocated(ppfd_shade240))      allocate(ppfd_shade240(nlat*nlon,modlays))

    if(.not.allocated(daily_maxt2m))       allocate(daily_maxt2m(nlat*nlon))
    if(.not.allocated(daily_mint2m))       allocate(daily_mint2m(nlat*nlon))
    if(.not.allocated(daily_maxws10m))     allocate(daily_maxws10m(nlat*nlon))

    if (hist_opt .eq. 1) then
        if(.not.allocated(tleaf_sun24_tmp_3d))    allocate(tleaf_sun24_tmp_3d(25,nlon,nlat,modlays))
        if(.not.allocated(tleaf_shade24_tmp_3d))  allocate(tleaf_shade24_tmp_3d(25,nlon,nlat,modlays))
        if(.not.allocated(tleaf_ave24_tmp_3d))    allocate(tleaf_ave24_tmp_3d(25,nlon,nlat,modlays))
        if(.not.allocated(ppfd_sun24_tmp_3d))     allocate(ppfd_sun24_tmp_3d(25,nlon,nlat,modlays))
        if(.not.allocated(ppfd_shade24_tmp_3d))   allocate(ppfd_shade24_tmp_3d(25,nlon,nlat,modlays))
        if(.not.allocated(tleaf_sun240_tmp_3d))   allocate(tleaf_sun240_tmp_3d(241,nlon,nlat,modlays))
        if(.not.allocated(tleaf_shade240_tmp_3d)) allocate(tleaf_shade240_tmp_3d(241,nlon,nlat,modlays))
        if(.not.allocated(tleaf_ave240_tmp_3d))   allocate(tleaf_ave240_tmp_3d(241,nlon,nlat,modlays))
        if(.not.allocated(ppfd_sun240_tmp_3d))    allocate(ppfd_sun240_tmp_3d(241,nlon,nlat,modlays))
        if(.not.allocated(ppfd_shade240_tmp_3d))  allocate(ppfd_shade240_tmp_3d(241,nlon,nlat,modlays))
        if(.not.allocated(tmp2mref_tmp_3d))       allocate(tmp2mref_tmp_3d(25,nlon,nlat))
        if(.not.allocated(ubzref_tmp_3d))         allocate(ubzref_tmp_3d(25,nlon,nlat))
    end if

    if(.not.allocated(tleaf_sun24_3d))        allocate(tleaf_sun24_3d(nlon,nlat,modlays))
    if(.not.allocated(tleaf_shade24_3d))      allocate(tleaf_shade24_3d(nlon,nlat,modlays))
    if(.not.allocated(tleaf_ave24_3d))        allocate(tleaf_ave24_3d(nlon,nlat,modlays))
    if(.not.allocated(ppfd_sun24_3d))         allocate(ppfd_sun24_3d(nlon,nlat,modlays))
    if(.not.allocated(ppfd_shade24_3d))       allocate(ppfd_shade24_3d(nlon,nlat,modlays))
    if(.not.allocated(tleaf_sun240_3d))       allocate(tleaf_sun240_3d(nlon,nlat,modlays))
    if(.not.allocated(tleaf_shade240_3d))     allocate(tleaf_shade240_3d(nlon,nlat,modlays))
    if(.not.allocated(tleaf_ave240_3d))       allocate(tleaf_ave240_3d(nlon,nlat,modlays))
    if(.not.allocated(ppfd_sun240_3d))        allocate(ppfd_sun240_3d(nlon,nlat,modlays))
    if(.not.allocated(ppfd_shade240_3d))      allocate(ppfd_shade240_3d(nlon,nlat,modlays))

    if(.not.allocated(daily_maxt2m_2d))       allocate(daily_maxt2m_2d(nlon,nlat))
    if(.not.allocated(daily_mint2m_2d))       allocate(daily_mint2m_2d(nlon,nlat))
    if(.not.allocated(daily_maxws10m_2d))     allocate(daily_maxws10m_2d(nlon,nlat))

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Wind Outputs
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        write(*,*)  'Canopy wind and/or WAF option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(canBOT))        allocate(canBOT(modlays))
        if(.not.allocated(canTOP))        allocate(canTOP(modlays))
        if(.not.allocated(canWIND))       allocate(canWIND(nlat*nlon,modlays))
        if(.not.allocated(canWIND_3d))    allocate(canWIND_3d(nlon,nlat,modlays))
        if(.not.allocated(dx))            allocate(dx(nlat*nlon))
        if(.not.allocated(dx_2d))         allocate(dx_2d(nlon,nlat))
        if(.not.allocated(waf))           allocate(waf(nlat*nlon))
        if(.not.allocated(waf_2d))        allocate(waf_2d(nlon,nlat))
        if(.not.allocated(flameh))        allocate(flameh(nlat*nlon))
        if(.not.allocated(flameh_2d))     allocate(flameh_2d(nlon,nlat))
    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Diffusivity Profile Outputs
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        write(*,*)  'Canopy eddy Kz option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(Kz))            allocate(Kz(nlat*nlon,modlays))
        if(.not.allocated(Kz_3d))         allocate(Kz_3d(nlon,nlat,modlays))
    end if


!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Photolysis Correction Factor Outputs
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        write(*,*)  'Canopy photolysis option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(rjcf))            allocate(rjcf(nlat*nlon,modlays))
        if(.not.allocated(rjcf_3d))         allocate(rjcf_3d(nlon,nlat,modlays))
    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Biogenic Emissions Outputs
!-------------------------------------------------------------------------------

    if (ifcanbio) then
        write(*,*)  'Canopy biogenic emissions option selected'
        write(*,*)  '-------------------------------'
        if (biospec_opt == 0 .or. biospec_opt == 1) then
            if(.not.allocated(emi_isop))         allocate(emi_isop(nlat*nlon,modlays))
            if(.not.allocated(emi_isop_3d))      allocate(emi_isop_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 2) then
            if(.not.allocated(emi_myrc))         allocate(emi_myrc(nlat*nlon,modlays))
            if(.not.allocated(emi_myrc_3d))      allocate(emi_myrc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 3) then
            if(.not.allocated(emi_sabi))         allocate(emi_sabi(nlat*nlon,modlays))
            if(.not.allocated(emi_sabi_3d))      allocate(emi_sabi_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 4) then
            if(.not.allocated(emi_limo))         allocate(emi_limo(nlat*nlon,modlays))
            if(.not.allocated(emi_limo_3d))      allocate(emi_limo_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 5) then
            if(.not.allocated(emi_care))         allocate(emi_care(nlat*nlon,modlays))
            if(.not.allocated(emi_care_3d))      allocate(emi_care_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 6) then
            if(.not.allocated(emi_ocim))         allocate(emi_ocim(nlat*nlon,modlays))
            if(.not.allocated(emi_ocim_3d))      allocate(emi_ocim_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 7) then
            if(.not.allocated(emi_bpin))         allocate(emi_bpin(nlat*nlon,modlays))
            if(.not.allocated(emi_bpin_3d))      allocate(emi_bpin_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 8) then
            if(.not.allocated(emi_apin))         allocate(emi_apin(nlat*nlon,modlays))
            if(.not.allocated(emi_apin_3d))      allocate(emi_apin_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 9) then
            if(.not.allocated(emi_mono))         allocate(emi_mono(nlat*nlon,modlays))
            if(.not.allocated(emi_mono_3d))      allocate(emi_mono_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 10) then
            if(.not.allocated(emi_farn))         allocate(emi_farn(nlat*nlon,modlays))
            if(.not.allocated(emi_farn_3d))      allocate(emi_farn_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 11) then
            if(.not.allocated(emi_cary))         allocate(emi_cary(nlat*nlon,modlays))
            if(.not.allocated(emi_cary_3d))      allocate(emi_cary_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 12) then
            if(.not.allocated(emi_sesq))         allocate(emi_sesq(nlat*nlon,modlays))
            if(.not.allocated(emi_sesq_3d))      allocate(emi_sesq_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 13) then
            if(.not.allocated(emi_mbol))         allocate(emi_mbol(nlat*nlon,modlays))
            if(.not.allocated(emi_mbol_3d))      allocate(emi_mbol_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 14) then
            if(.not.allocated(emi_meth))         allocate(emi_meth(nlat*nlon,modlays))
            if(.not.allocated(emi_meth_3d))      allocate(emi_meth_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 15) then
            if(.not.allocated(emi_acet))         allocate(emi_acet(nlat*nlon,modlays))
            if(.not.allocated(emi_acet_3d))      allocate(emi_acet_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 16) then
            if(.not.allocated(emi_co))           allocate(emi_co(nlat*nlon,modlays))
            if(.not.allocated(emi_co_3d))        allocate(emi_co_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 17) then
            if(.not.allocated(emi_bvoc))         allocate(emi_bvoc(nlat*nlon,modlays))
            if(.not.allocated(emi_bvoc_3d))      allocate(emi_bvoc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 18) then
            if(.not.allocated(emi_svoc))         allocate(emi_svoc(nlat*nlon,modlays))
            if(.not.allocated(emi_svoc_3d))      allocate(emi_svoc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 19) then
            if(.not.allocated(emi_ovoc))         allocate(emi_ovoc(nlat*nlon,modlays))
            if(.not.allocated(emi_ovoc_3d))      allocate(emi_ovoc_3d(nlon,nlat,modlays))
        end if
    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Gas Dry Deposition Outputs
!-------------------------------------------------------------------------------

    if (ifcanddepgas) then
        write(*,*)  'Canopy gas dry deposition option selected'
        write(*,*)  '-------------------------------'
        if (chemmechgas_opt == 0) then !RACM2 --> 31 species
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                if(.not.allocated(ddep_no))              allocate(ddep_no(nlat*nlon,modlays))
                if(.not.allocated(ddep_no_3d))           allocate(ddep_no_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                if(.not.allocated(ddep_no2))              allocate(ddep_no2(nlat*nlon,modlays))
                if(.not.allocated(ddep_no2_3d))           allocate(ddep_no2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                if(.not.allocated(ddep_o3))              allocate(ddep_o3(nlat*nlon,modlays))
                if(.not.allocated(ddep_o3_3d))           allocate(ddep_o3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                if(.not.allocated(ddep_hono))              allocate(ddep_hono(nlat*nlon,modlays))
                if(.not.allocated(ddep_hono_3d))           allocate(ddep_hono_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                if(.not.allocated(ddep_hno4))              allocate(ddep_hno4(nlat*nlon,modlays))
                if(.not.allocated(ddep_hno4_3d))           allocate(ddep_hno4_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                if(.not.allocated(ddep_hno3))              allocate(ddep_hno3(nlat*nlon,modlays))
                if(.not.allocated(ddep_hno3_3d))           allocate(ddep_hno3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                if(.not.allocated(ddep_n2o5))              allocate(ddep_n2o5(nlat*nlon,modlays))
                if(.not.allocated(ddep_n2o5_3d))           allocate(ddep_n2o5_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                if(.not.allocated(ddep_co))              allocate(ddep_co(nlat*nlon,modlays))
                if(.not.allocated(ddep_co_3d))           allocate(ddep_co_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                if(.not.allocated(ddep_h2o2))              allocate(ddep_h2o2(nlat*nlon,modlays))
                if(.not.allocated(ddep_h2o2_3d))           allocate(ddep_h2o2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                if(.not.allocated(ddep_ch4))              allocate(ddep_ch4(nlat*nlon,modlays))
                if(.not.allocated(ddep_ch4_3d))           allocate(ddep_ch4_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                if(.not.allocated(ddep_mo2))              allocate(ddep_mo2(nlat*nlon,modlays))
                if(.not.allocated(ddep_mo2_3d))           allocate(ddep_mo2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                if(.not.allocated(ddep_op1))              allocate(ddep_op1(nlat*nlon,modlays))
                if(.not.allocated(ddep_op1_3d))           allocate(ddep_op1_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                if(.not.allocated(ddep_moh))              allocate(ddep_moh(nlat*nlon,modlays))
                if(.not.allocated(ddep_moh_3d))           allocate(ddep_moh_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                if(.not.allocated(ddep_no3))              allocate(ddep_no3(nlat*nlon,modlays))
                if(.not.allocated(ddep_no3_3d))           allocate(ddep_no3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                if(.not.allocated(ddep_o3p))              allocate(ddep_o3p(nlat*nlon,modlays))
                if(.not.allocated(ddep_o3p_3d))           allocate(ddep_o3p_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                if(.not.allocated(ddep_o1d))              allocate(ddep_o1d(nlat*nlon,modlays))
                if(.not.allocated(ddep_o1d_3d))           allocate(ddep_o1d_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                if(.not.allocated(ddep_ho))              allocate(ddep_ho(nlat*nlon,modlays))
                if(.not.allocated(ddep_ho_3d))           allocate(ddep_ho_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                if(.not.allocated(ddep_ho2))              allocate(ddep_ho2(nlat*nlon,modlays))
                if(.not.allocated(ddep_ho2_3d))           allocate(ddep_ho2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                if(.not.allocated(ddep_ora1))              allocate(ddep_ora1(nlat*nlon,modlays))
                if(.not.allocated(ddep_ora1_3d))           allocate(ddep_ora1_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                if(.not.allocated(ddep_hac))              allocate(ddep_hac(nlat*nlon,modlays))
                if(.not.allocated(ddep_hac_3d))           allocate(ddep_hac_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                if(.not.allocated(ddep_paa))              allocate(ddep_paa(nlat*nlon,modlays))
                if(.not.allocated(ddep_paa_3d))           allocate(ddep_paa_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                if(.not.allocated(ddep_dhmob))              allocate(ddep_dhmob(nlat*nlon,modlays))
                if(.not.allocated(ddep_dhmob_3d))           allocate(ddep_dhmob_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                if(.not.allocated(ddep_hpald))              allocate(ddep_hpald(nlat*nlon,modlays))
                if(.not.allocated(ddep_hpald_3d))           allocate(ddep_hpald_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                if(.not.allocated(ddep_ishp))              allocate(ddep_ishp(nlat*nlon,modlays))
                if(.not.allocated(ddep_ishp_3d))           allocate(ddep_ishp_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                if(.not.allocated(ddep_iepox))              allocate(ddep_iepox(nlat*nlon,modlays))
                if(.not.allocated(ddep_iepox_3d))           allocate(ddep_iepox_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                if(.not.allocated(ddep_propnn))              allocate(ddep_propnn(nlat*nlon,modlays))
                if(.not.allocated(ddep_propnn_3d))           allocate(ddep_propnn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                if(.not.allocated(ddep_isopnb))              allocate(ddep_isopnb(nlat*nlon,modlays))
                if(.not.allocated(ddep_isopnb_3d))           allocate(ddep_isopnb_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                if(.not.allocated(ddep_isopnd))              allocate(ddep_isopnd(nlat*nlon,modlays))
                if(.not.allocated(ddep_isopnd_3d))           allocate(ddep_isopnd_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                if(.not.allocated(ddep_macrn))              allocate(ddep_macrn(nlat*nlon,modlays))
                if(.not.allocated(ddep_macrn_3d))           allocate(ddep_macrn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                if(.not.allocated(ddep_mvkn))              allocate(ddep_mvkn(nlat*nlon,modlays))
                if(.not.allocated(ddep_mvkn_3d))           allocate(ddep_mvkn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                if(.not.allocated(ddep_isnp))              allocate(ddep_isnp(nlat*nlon,modlays))
                if(.not.allocated(ddep_isnp_3d))           allocate(ddep_isnp_3d(nlon,nlat,modlays))
            end if
        else
            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
            write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
            call exit(2)
        end if
    end if

END SUBROUTINE canopy_alloc
