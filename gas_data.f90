!!********************************** stepsize parameter (by feng fan)******
      module gas_data

      real, parameter :: max_stepsize=5    !unit sec     1<=max_stepsize<=60
      real, parameter :: min_stepsize=5    !unit sec     1<=min_stepsize<=5  

      real, parameter  :: cri_dydt=7e6  ! unit molec/(cc*sec)

      real, parameter :: growth_coeff=1.2 ! growth_coeff must be >=1.
      real, parameter :: decay_coeff=1.8 ! decay_coeff must be >=1.

      integer, parameter :: &
          nx=148, &
          ny = 160,  &
          nzz = 20, &
          npt = nx*ny*nzz

      integer, parameter :: VLEN=16 ! vector length

!--------------------------------------------------------------------------

      parameter(ngas_com = 34,  & !li jie has changed 32 to 34
                ngas_urb = 21,  & !li jie has changed 19 to 21
                ngas_bio = 10,  & !li jie has changed 5 to 10
                ngas_mar = 11)

      parameter(ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar)

      parameter(nmax = ngas_max)

!------------------------------------------------------------------------
!--------------------------------------------------------------------------

      parameter(nperox = 10)	! total number of alkylperoxy radicals
      parameter(nphoto = 20)	! total number of photolyzing species

      

      parameter(nrxn_het = 28, & !li jie has changed 7 to 28
                nrxn_com = 76, & !feng fan has changed 74 to 76 for no3 sinks with vocs
                nrxn_urb = 49, & !feng fan has changed 44 to 49 for no3 sinks with vocs
                nrxn_bio = 20, & !li jie has changed 16 to 20
                nrxn_mar = 35)

 

      parameter(nreg1 = ngas_com, &
                nreg2 = ngas_com + ngas_urb, &
                nreg3 = ngas_com + ngas_urb + ngas_bio, &
                nreg4 = ngas_com + ngas_mar, &
                nreg5 = ngas_com + ngas_urb + ngas_mar, &
                nreg6 = ngas_com + ngas_urb + ngas_bio + ngas_mar)

!------------------------------------------------------------------------


! global species indices

      common /globalgas/ &
       kh2so4,      khno3,       khcl,        knh3,        kno, &
       kno2,        kno3,        kn2o5,       khono,       khno4, &
       ko3,         ko1d,        ko3p,        koh,         kho2, &
       kh2o2,       kco,         kso2,        kch4,        kc2h6, &          
       kch3o2,      kethp,       khcho,       kch3oh,      kanol, &
       kch3ooh,     kethooh,     kald2,       khcooh,      krcooh, &
       kc2o3,       kpan,        kpar,        kaone,       kmgly, &
       keth,        kolet,       kolei,       ktol,        kxyl, &
       kcres,       kto2,        kcro,        kopen,       konit, &
       krooh,       kro2,        kano2,       knap,        kxo2, &
       kxpar,       kisop,       kisoprd,     kisopp,      kisopn, &
       kisopo2, &
       kdms,        kmsa,        kdmso,       kdmso2,      kch3so2h, &
       kch3sch2oo,  kch3so2,     kch3so3,     kch3so2ch2oo,kch3so2oo, &
       ksulfhox, &
       kterp,       ksv1,        ksv2,        ksv3,        ksv4, &  
       ksv5,        ksv6,        kdso4,       kdno3  ! kterp,ksv1-6,kdso4,kdno3 are added by li jie

!-------------------------------------------------------------------------

      common /constants/ pi, avogad, deg2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef KNL_OPT

   !dir$ attributes align:64 :: cnn
      !common /globalconcentration/ cnn(nmax)
      !common/emissions/ emission(nmax), emit(nmax)

      common/gasrateconst/cnn(VLEN,nmax),emission(VLEN,nmax), emit(VLEN,nmax), &
                          rk_bio(VLEN,nrxn_bio),rk_mar(VLEN,nrxn_mar), &
                          !rk_het(nrxn_het),
                          rk_param(VLEN,nperox),&
                          !r_urb(nrxn_urb), &
                          

      !common/gasratescom/ r_urb(nrxn_urb), &
                          r_bio(VLEN,nrxn_bio),r_mar(VLEN,nrxn_mar), &
                          !r_het(nrxn_het), &
                          !r1_urb(nrxn_urb), &
                          r1_bio(VLEN,nrxn_bio),r1_mar(VLEN,nrxn_mar), &
                          !r1_het(nrxn_het), &
                          !r2_urb(nrxn_urb), &
                          r2_bio(VLEN,nrxn_bio),r2_mar(VLEN,nrxn_mar), &
                          !r2_het(nrxn_het), &

      !common/gaspd/ &
                    !p_urb(ngas_max),rl_urb(ngas_max), &
                    !p_bio(VLEN,ngas_max),rl_bio(VLEN,ngas_max), &
                    !p_mar(VLEN,ngas_max),rl_mar(VLEN,ngas_max), &
                    !p_het(ngas_max),rl_het(ngas_max), &
                    total_p(VLEN,ngas_max),total_l(VLEN,ngas_max), &
       tcur_sec(VLEN), tmid_sec(VLEN), told_sec(VLEN), dt_sec(VLEN),&
       rlon(VLEN), rlat(VLEN), zalt_m(VLEN),    cos_sza(VLEN), &
                    msolar, mphoto, &
      !common /timenposition/ &
       tbeg_dd,   tbeg_hh,   tbeg_mm,   tbeg_ss, &
       trun_dd,   trun_hh,   trun_mm,   trun_ss, dt_min, &
       !tbeg_plus_trun_sec  !!notice this (by feng fan)


      !common /photolysis/ &
#endif
#ifdef VEC_OPT
      common /localgas/ &
!       ih2so4,      ihno3,       ihcl,        inh3,        ino, &
!       ino2,        ino3,        in2o5,       ihono,       ihno4, &
!       io3,         io1d,        io3p,        ioh,         iho2, &
!       ih2o2,       ico,         iso2,        ich4,        ic2h6, &           
!       ich3o2,      iethp,       ihcho,       ich3oh,      ianol, &
!       ich3ooh,     iethooh,     iald2,       ihcooh,      ircooh, &
!       ic2o3,       ipan,        idso4,       idno3, &       
!       ipar,        iaone,       imgly, &
!       ieth,        iolet,       iolei,       itol,        ixyl, &
!       icres,       ito2,        icro,        iopen,       ionit, &
!       irooh,       iro2,        iano2,       inap,        ixo2, &
!       ixpar,       isv1, isv2, &

!       iisop,       iisoprd,     iisopp,      iisopn, &
!       iisopo2,     iterp,       isv3,        isv4,        isv5,        isv6 &

       idms,        imsa,        idmso,       idmso2,      ich3so2h, &
       ich3sch2oo,  ich3so2,     ich3so3,     ich3so2ch2oo,ich3so2oo, isulfhox
#endif         
              ! kterp,ksv1-6,kdso4,kdno3 are added by li jie
      type cbmztype
        ! Padding to align 64, by Junmin, 10/25/2016, may cause 693 Error
        !real :: rk_com(nrxn_com+4),rk_urb(nrxn_urb+15),rk_het(nrxn_het+4)
        real :: rk_com(VLEN,nrxn_com),rk_urb(VLEN,nrxn_urb),rk_het(VLEN,nrxn_het)
        real :: aperox(nperox,nperox),bperox(nperox,nperox) ! constants to promote
        real :: rk_photo(VLEN,nphoto)

        real :: cair_mlc(VLEN)
        real :: te(VLEN)
        real :: h2o(VLEN)
        real :: h2(VLEN)
        real :: o2(VLEN)
        real :: pr_atm(VLEN), rh(VLEN)

        integer :: iregime(VLEN)
        logical :: pmask(VLEN) ! indicate early completion, true=>false one-way change
        logical :: bmask(VLEN) ! current control flow status
        
        real :: p_het(VLEN,76), &! help to reduce memset overheads
                p_com(VLEN,76), &
                rl_com(VLEN,76), &
                r2_com(VLEN,76), &
                p_urb(VLEN,76), &
                rl_het(VLEN,76), &
                rl_urb(VLEN,76), &
                p_bio(VLEN,76), &
                rl_bio(VLEN,76), &
                p_mar(VLEN,76), &
                rl_mar(VLEN,76)
      end type cbmztype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef OMP_OPT

!$omp threadprivate(/gasrateconst/)
#endif


!      common/het/ npcasp(15)


!------------------------------------------------------------------------


      common /gasperoxyradicals/ &
       jch3o2,      jethp,       jro2,        jc2o3,       jano2, &
       jnap,        jisopp,      jisopn,      jisopo2,     jxo2

      common /photolyzingspecies/ &
       jphoto_no2,    jphoto_no3,   jphoto_hono,   jphoto_hno3, & 
       jphoto_hno4,   jphoto_n2o5,  jphoto_o3a,    jphoto_o3b, &    
       jphoto_h2o2,   jphoto_hchoa, jphoto_hchob,  jphoto_ch3ooh, &
       jphoto_ethooh, jphoto_ald2,  jphoto_aone,   jphoto_mgly, &
       jphoto_open,   jphoto_rooh,  jphoto_onit,   jphoto_isoprd

!!*****************************************************************
!!*****************************************************************

!!   place various constant or initial values into common


!!-------------------------------------------------
!!     define fundamental constants...
      data pi		/3.141592654/
      data avogad	/6.02217e+23/
      data deg2rad	/0.017453293/

!!--------------------------------
!! define species indices

!! species in inorgani!! chemistry
      data kh2so4       /  1/
      data khno3        /  2/
      data khcl         /  3/
      data knh3         /  4/
      data kno          /  5/
      data kno2         /  6/
      data kno3         /  7/
      data kn2o5        /  8/
      data khono        /  9/
      data khno4        / 10/
      data ko3          / 11/
      data ko1d         / 12/
      data ko3p         / 13/
      data koh          / 14/
      data kho2         / 15/
      data kh2o2        / 16/
      data kco          / 17/
      data kso2         / 18/

!! species in methane, ethane, formaldehyde chemistry
      data kch4         / 19/
      data kc2h6        / 20/
      data kch3o2       / 21/
      data kethp        / 22/
      data khcho        / 23/
      data kch3oh       / 24/
      data kanol	/ 25/
      data kch3ooh      / 26/
      data kethooh	/ 27/
      data kald2        / 28/
      data khcooh	/ 29/
      data krcooh	/ 30/
      data kc2o3	/ 31/
      data kpan		/ 32/
 
!! species in hc1 mechanism. initialize indices to zero
      data kpar         / 33/
      data kaone        / 34/
      data kmgly        / 35/
      data keth         / 36/
      data kolet        / 37/
      data kolei        / 38/
      data ktol         / 39/
      data kxyl         / 40/
      data kcres        / 41/
      data kto2         / 42/
      data kcro         / 43/
      data kopen        / 44/
      data konit        / 45/
      data krooh	/ 46/
      data kro2         / 47/
      data kano2	/ 48/
      data knap		/ 49/
      data kxo2		/ 50/
      data kxpar	/ 51/

!! species in hc2 mechanism. initialize indices to zero
      data kisop	/ 52/
      data kisoprd	/ 53/
      data kisopp	/ 54/
      data kisopn	/ 55/
      data kisopo2	/ 56/

!! species in dms mechanism. initialize indices to zero
      data kdms         / 57/
      data kmsa         / 58/
      data kdmso        / 59/
      data kdmso2       / 60/
      data kch3so2h     / 61/
      data kch3sch2oo   / 62/
      data kch3so2      / 63/
      data kch3so3      / 64/
      data kch3so2oo    / 65/
      data kch3so2ch2oo / 66/
      data ksulfhox     / 67/
!! species in soa formation 2 species in urban, 5 species in biogenic
      data kterp        / 68/ 
      data ksv1         / 69/
      data ksv2         / 70/
      data ksv3         / 71/
      data ksv4         / 72/
      data ksv5         / 73/
      data ksv6         / 74/
!! species in aso4 and ano3 from heteorogeneous chemistry
      data kdso4        / 75/
      data kdno3        / 76/

!! regime-dependent chemistry definitions

      data iregime	/  1/

!!     gas

!! species in common chemistry

!! alkylperoxy radical indices for parameterized permutation reactions
      data jch3o2	/  1/
      data jethp	/  2/
      data jro2		/  3/
      data jc2o3	/  4/
      data jano2	/  5/
      data jnap		/  6/
      data jisopp	/  7/
      data jisopn	/  8/
      data jisopo2	/  9/
      data jxo2		/ 10/

!! photolyzing species indices
      data jphoto_no2	/  1/
      data jphoto_no3	/  2/
      data jphoto_hono	/  3/
      data jphoto_hno3	/  4/
      data jphoto_hno4	/  5/
      data jphoto_n2o5  /  6/
      data jphoto_o3a	/  7/
      data jphoto_o3b	/  8/
      data jphoto_h2o2	/  9/
      data jphoto_hchoa	/ 10/
      data jphoto_hchob	/ 11/
      data jphoto_ch3ooh/ 12/
      data jphoto_ethooh/ 13/
      data jphoto_ald2	/ 14/
      data jphoto_aone	/ 15/
      data jphoto_mgly	/ 16/
      data jphoto_open	/ 17/
      data jphoto_rooh	/ 18/
      data jphoto_onit	/ 19/
      data jphoto_isoprd/ 20/

      end module gas_data

