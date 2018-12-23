!!*********************************************************************
!!*** revised by feng fan, feb,11th,2015 
!!*** mbe solver, 76 species
!!*** modified backward euler (mbe)
!!*********************************************************************

      module mod_cbmz

#ifdef KNL_OPT
      integer,parameter :: &
      ih2so4	= 1, &
      ihno3	= 2, &
      ihcl	= 3, &
      inh3	= 4, &
      ino	= 5, &
      ino2	= 6, &
      ino3	= 7, &
      in2o5	= 8, &
      ihono	= 9, &
      ihno4	= 10, &
      io3	= 11, &
      io1d	= 12, &
      io3p	= 13, &
      ioh	= 14, &
      iho2	= 15, &
      ih2o2	= 16, &
      ico	= 17, &
      iso2	= 18, &
      ich4	= 19, &
      ic2h6	= 20, &
      ich3o2	= 21, &
      iethp	= 22, &
      ihcho	= 23, &
      ich3oh	= 24, &
      ianol	= 25, &
      ich3ooh	= 26, &
      iethooh	= 27, &
      iald2	= 28, &
      ihcooh	= 29, &
      ircooh	= 30, &
      ic2o3	= 31, &
      ipan	= 32

#endif
#ifdef VEC_OPT
      !dir$ attributes align:64 :: rkindex,arrp1,arrp2,arnm,ark0i,cnnindex,hetindex,sindex,sindex2
      integer,parameter,dimension(10) :: cnnindex = &
        (/21,22,47,31,48,49,54,55,56,50/)
      real,parameter,dimension(48) :: arrp1 = &
        (/3.2e-11,1.8e-11,8.0e-12,6.5e-12,2.0e-12,1.2e-13,1.6e-12,1.1e-14, &
        5.5e-12,1.8e-11,1.3e-12,4.8e-11,2.9e-12,3.5e-12,1.5e-11,4.5e-14, &
        8.5e-13,2.8e-14,1.5e-17,6.7e-12,3.4e-13,3.8e-12,3.8e-12,3.0e-12, &
        2.6e-12,3.8e-13,7.5e-13,7.0e-12,5.6e-12,1.4e-12,5.3e-12,4.5e-13, &
        9.4e-13,4.8e26,3.7e26,1.1e28, &
        6.0e-34,7.2e-15,4.1e-16,1.9e-33,2.3e-13,1.7e-33,1.4e-21,8.18e-23,0.,0.,0.,0./)
      
      real,parameter,dimension(48) :: arrp2 = &
        (/70.0,110.0,-2060.0,-120.0,-1400.0,-2450.0,-940.0,-500.0, &
        -2000.0,-390.0,380.0,250.0,-160.0,250.0,170.0,-1260.0, &
        -2450.0,-1575.0,-492.0,-600.0,-1900.0,200.0,200.0,280.0, &
        365.0,800.0,700.0,-235.0,270.0,-1900.0,360.0,1000.0,-2650.0, &
        -10900.0,-11000.0,-14000.0, &
        0.0,785.0,1440.0,725.0,600.0,1000.0,2200.0,0.0, 0.,0.,0.,0./)
      !real,parameter,dimension(8) :: ark0 = &
      !  (/9.0e-32,9.0e-32,7.0e-31,2.5e-30,1.8e-31,2.2e-30,3.0e-31,9.7e-29/)
      !real,parameter,dimension(8) :: arnn = &
      !  (/2.0,1.5,2.6,4.4,3.2,3.9,3.3,5.6/)
      !real,parameter,dimension(8) :: arki = &
      !  (/2.2e-11,3.0e-11,3.6e-11,1.6e-11,4.7e-12,1.5e-12,1.5e-12,9.3e-12/)
      !real,parameter,dimension(8) :: armm = &
      !  (/0.0,0.0,0.1,1.7,1.4,0.7,0.0,1.5/)

      real, parameter,dimension(16) :: arnm = &
        (/2.0,1.5,2.6,4.4,3.2,3.9,3.3,5.6,0.0,0.0,0.1,1.7,1.4,0.7,0.0,1.5/)
      real,parameter,dimension(16) :: ark0i = &
        (/9.0e-32,9.0e-32,7.0e-31,2.5e-30,1.8e-31,2.2e-30,3.0e-31,9.7e-29,2.2e-11,3.0e-11,3.6e-11,1.6e-11,4.7e-12,1.5e-12,1.5e-12,9.3e-12/)

      !real,parameter,dimension(8) :: expp1 = &
      !  (/6.0e-34,7.2e-15,4.1e-16,1.9e-33,2.3e-13,1.7e-33,1.4e-21,8.18e-23/)
      !real,parameter,dimension(8) :: expp2 = &
      !  (/0.0,785.0,1440.0,725.0,600.0,1000.0,2200.0,0.0/)
      integer,parameter,dimension(28) :: hetindex = &
        (/8,6,7,15,23,14,11,6,2,8,11,2,6,7,&
        8,14,15,16,18,48,24,23,8,7,15,18,7,2/)
        ! Fixed by Junmin, 10/21/2016
      integer,parameter,dimension(76) :: sindex = &
        (/6,7,9,2,10,8,11,11,16,& ! 1 - 9
        12,12,12,13,13,13,13,13,11,11,11,11, & ! 10-21
        14,14,14,14,14,14,14,14,14, & ! 22-30
        15,15,15,15,15,10, & ! 31-36
        7,7,7,7,7,8,8,17,18,19,20,24,23,23,23,23, & ! 37-52
        26,27,26,27,21,22,21,22,21,22,21,22,25, & ! 53-65
        28,28,28,31,32,31,31,31,31,7,7/) ! 66-76 
      integer,parameter,dimension(76) :: sindex2 = &
        (/33,33,33,33,33,33,33,33,33, & ! 1-9
        33,33,33,33, & ! 10-13
        11,6,6,5,5,6,14,15,33,5,6,7,9,2,10,15,16,15,15,5,6,6,33, & ! 14-36
        5,6,6,7,15,33,33,14,14,14,14,14,33,33,14,7,33,33,14,14,5,5, & ! 37-58
        7,7,15,15,33,33,14,33,14,7,6,33,5,7,15,33,20,24/) 
      integer,parameter,dimension(44) :: rkindex = &
        (/10,11,14,15,18,19,20,21,22,26,28,29,30,33,37,38,40, &
        46,47,48,52,55,56,57,58,61,62,65,67,68,71,73,76,36,43,70,16,17,23,24,34,39,45,69/)

#endif

      contains

      subroutine  cbmz(myid,cbmzobj,cppb,factcld)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      !implicit real(a-h,o-z), integer(i-n) 
      integer :: myid
!dir$ assume_aligned cppb:64, factcld:64, cbmzobj:64
      type(cbmztype), target :: cbmzobj
      real,dimension(VLEN,76) :: cppb
      real,dimension(VLEN) :: factcld

      real :: tlocal,tdec,sidec,codec,tloc,thou
      real, dimension(VLEN) :: rlat1,rlat2
      real, dimension(:,:), pointer :: rk_com,rk_urb
      real :: tbeg_plus_trun_sec,tbeg_sec,trun_sec
      integer :: it
      integer :: i2 ! loop idx
      rk_com => cbmzobj%rk_com
      rk_urb => cbmzobj%rk_urb



      
      
      

      tbeg_sec = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss
      trun_sec = ((trun_dd*24.+trun_hh)*60.+trun_mm)*60.+trun_ss
      dt_sec(:)   = 60.*dt_min	! time step [seconds]
      tcur_sec(:) = tbeg_sec	! initialize current time to tbeg_sec
      !tcur_min = tcur_sec/60.
      !tcur_hrs = tcur_min/60.
      !it = 0
      rlon(:) = rlon(:)*deg2rad
      rlat(:) = rlat(:)*deg2rad       

      cbmzobj%rk_com(:,12) = 2.2e-10
      rk_com(:,25) = 2.2e-11
      rk_com(:,35) = 5.0e-16
      rk_com(:,41) = 3.5e-12
      rk_com(:,42) = 2.0e-21
      rk_com(:,51) = 1.0e-11
      rk_com(:,59) = 1.1e-12
      rk_com(:,60) = 2.5e-12
      rk_com(:,72) = 4.0e-12
      rk_com(:,75) = 1.0e-17   ! ethane + no3

      call setaircomposition(cbmzobj)
      call loadperoxyparameters(cbmzobj)! aperox and bperox  

      dt_sec(:)= min_stepsize   !!!!!! notice this, by feng fan
      tbeg_plus_trun_sec=tbeg_sec+trun_sec


      ! Need only one call at the end of subroutine
      !call printresult(cbmzobj,cppb)                       

       it = 1
!!------------------------------------------------------------------

!!----------------- main time-loop begins...  ----------------------

     rlat1(:) = sin(rlat(:))
     rlat2(:) = cos(rlat(:))

     !TODO tcur_sec set mask
     do while (tcur_sec(1) <tbeg_plus_trun_sec)

       it = it+1 

       call updatetime
      if(msolar.eq.1)then
      !!dir$ simd private(tlocal,tdec,sidec,codec,tloc,thou)
      do i2=1,VLEN
       tlocal=tmid_sec(i2)
       tdec=0.4092797*sin(1.992385e-7*tlocal)                           
       sidec=sin(tdec)                                                  
       codec=cos(tdec)                                                  
       tloc=7.272205e-5*tlocal                                           
       thou=cos(rlon(i2)+tloc)                                      
       cos_sza(i2)=rlat1(i2)*sidec+rlat2(i2)*codec*thou         
      enddo !i2
      endif

      !call updatemetfields     !at present, "call updatemetfields" is of no use.

      !call updateemissions     !at present, "call updateemissions" is of no use.

#ifdef KNL_OPT
      call integratechemistry(cbmzobj,factcld)   !!! solver is in this
#else
      call integratechemistry(factcld)   !!! solver is in this
      ! xinweix no use of result
      call domassbalance       
#endif

      !call printresult(cbmzobj,cppb)                       

    end do

    call printresult(cbmzobj,cppb)                       

  return
  end            

!!**********************************************************************
!!*********************  functions *************************************


!!------------------------------------------------------
!!------------------------------------------------------

#ifdef KNL_OPT
      real function arr(te,aa,bb)
      use gas_data
      implicit none
      real :: te,aa,bb
      arr = aa*exp(bb/te)
      return
      end
#else
      function arr(aa,bb)
      use gas_data
      arr = aa*exp(bb/te)
      return
      end
#endif


!!------------------------------------------------------
!!------------------------------------------------------
	integer function nbllen( str )

!!   returns the position of the last non-blank character in str

	character*(*) str

	j = len(str)

	if (j .gt. 0) then
1000	    if (str(j:j) .eq. ' ') then
		j = j - 1
		if (j .gt. 0) goto 1000
	    end if
	end if
	nbllen = max0( j, 0 )

	return
	end

!!------------------------------------------------------
!!------------------------------------------------------
      real function troe(cair_mlc,te,rk0,rnn,rki,rmm)
      implicit none
      real :: cair_mlc,te,rk0,rnn,rki,rmm
      real :: expo
      rk0 = rk0*cair_mlc*(te/300.)**(-rnn)
      rki = rki*(te/300.)**(-rmm)
      expo= 1./(1. + (alog10(rk0/rki))**2)
      troe  = (rk0*rki/(rk0+rki))*.6**expo
      return
      end


!!------------------------------------------------------
!!------------------------------------------------------
!!--  function fuchs(kn, a)   !!!(used to compute rk_  )!!!!

      real function fuchs(kn, a)
      implicit none
      real :: kn,a

      fuchs = 0.75*a*(1. + kn)/(kn**2 + kn + 0.283*kn*a + 0.75*a)

      return
      end

!!********************************************************************
!! function watervapor
!! purpose: calculates concentration of h2o using the method given in 
!!          seinfeld's book, pg. 181
!!---------------------------------------------------------------------

      real function watervapor(rh, cair_mlc, te, pr_atm)
      implicit none
      real :: rh,cair_mlc,te,pr_atm
      real :: t_steam,pr_std,a,arg,pr_h2o

      t_steam = 373.15			! steam temperature  [k]
      pr_std   = 1.0			! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)  			! [atm]
      watervapor = rh*(pr_h2o/pr_atm)*cair_mlc/100.	! [molec/cc]

      return
      end







!!**********************************************************************
!!*********************  subroutines ***********************************

#ifndef KNL_OPT
      subroutine readinputfile(lin)
      use gas_data

      character*40 dword

      read(lin,*)dword

!!----------- begin time from 12:00 (noon) march 21 [min]
      read(lin,*)dword, tbeg_dd, tbeg_hh, tbeg_mm, tbeg_ss, dword
      read(lin,*)dword, trun_dd, trun_hh, trun_mm, trun_ss, dword
      read(lin,*)dword, dt_min,dword ! transport time-step [min]
      read(lin,*)dword, rlon,  dword ! longitude [deg]
      read(lin,*)dword, rlat,  dword ! latitude [deg]
      read(lin,*)dword, zalt_m,dword ! altitude  above mean sea level [m]
      read(lin,*)dword, rh,    dword ! relative humidity [%]
      read(lin,*)dword, te,    dword ! temperature [k]
      read(lin,*)dword, pr_atm,dword ! pressure [atm]
      read(lin,*)dword, msolar,dword ! msolar flag
      read(lin,*)dword, mphoto,dword ! mphoto flag
      read(lin,*)dword, iprint,dword ! freq of output

!!----------- read index, species name, initial conc, and emissions
      read(lin,*)dword
      read(lin,*)dword

      read(lin,*)dword ! gas

      do i=1, 76
        read(lin,*)k, species(k), cnn(k), emission(k)   ! k?
      
      enddo

      write(6,*)'finished reading all inputs...'
      return
      end
#endif



!!**********************************************************************
!!***********************************************************************
      subroutine integratechemistry(cbmzobj,factcld)
      !use gas_data, only : told_sec,cbmztype
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
!dir$ assume_aligned factcld:64
      real, dimension(VLEN) :: factcld
!dir$ attributes align:64 :: t_in
      real, dimension(VLEN) :: t_in


      t_in = told_sec

      call gaschemistry(cbmzobj,t_in,factcld)     !!! odesolver(ntot,stot,t_in) is in this

      return
      end


!!***********************************************************************
!!********************************************************************
      subroutine gaschemistry(cbmzobj,t_in,factcld)
      !use gas_data, only: cbmztype
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
!dir$ assume_aligned t_in:64, factcld:64
      real, dimension(VLEN) :: t_in, factcld

!dir$ attributes align:64 :: stot, ntot
      real, dimension(VLEN,76) :: stot !local species array
      integer, dimension(VLEN) :: ntot

      call selectgasregime(cbmzobj,ntot)! selects iregime and calculates ntot
      call peroxyrateconstants(cbmzobj)
      call gasrateconstants(cbmzobj,factcld)
      call setgasindices(cbmzobj) ! set gas indices for selected iregime 
      call mapgasspecies(cbmzobj,stot,0)! map cnn into stot for selected iregime
      call odesolver(cbmzobj,ntot,stot,t_in)
      call mapgasspecies(cbmzobj,stot,1)! map stot back into cnn


      
      return
      end


#ifndef KNL_OPT
!!**********************************************************************
!!**********************************************************************
      subroutine updateemissions
      use gas_data

!! update emission fluxes here

      return
      end



!!**********************************************************************
!!***********************************************************************
      subroutine updatemetfields
      use gas_data

!! update temperature, pressure, cloud flag, etc. here

      return
      end
#endif





!!**********************************************************************
!! all variables are array(VLEN), by Junmin, 01/18/2017
!!***********************************************************************
      subroutine updatetime

      use gas_data

     
      told_sec = tcur_sec

      tcur_sec = tcur_sec + dt_sec
#ifndef KNL_OPT
      tcur_min = tcur_sec/60.
      tcur_hrs = tcur_min/60.
#endif
      tmid_sec = told_sec + 0.5*dt_sec

      return
      end




!!*****************************************************************
!!     subroutine for printing output at iprint time steps
!!***************************************************************** 
      subroutine printresult(cbmzobj,cppb)
      !use gas_data, only : cbmztype,cnn
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(VLEN,76) :: cppb
!dir$ assume_aligned cppb:64
      real, dimension(:), pointer :: cair_mlc
      integer :: i, i2
      cair_mlc => cbmzobj%cair_mlc
      !ppb => cbmzobj%ppb


!=================================================================
!
!     converting (molecules/cc) to (ppb)
!-----------------------------------------------------------------
#ifdef KNL_OPT
      !dir$ ivdep
#endif
      do i=1,76
        cppb(:,i) = (cnn(:,i)/cair_mlc(:))*1.e+9
      enddo

!=================================================================
!
!     gas-phase species
!-----------------------------------------------------------------
!                                                            comm by lijie
!      if(it.eq.0)write(20,770)

!      write(20,7)tcur_hrs,cppb(ko3),cppb(kno),cppb(kno2),
!     &                    cppb(khno3),cppb(khono),cppb(kpan),
!     &                    cppb(konit),cppb(kh2o2),cppb(kolet),
!     &                    cppb(kolei),cppb(kpar),cppb(kald2),
!     &                    cppb(kisop)

!
!===================================================================
!
!     formats
!-------------------------------------------------------------------

!7     format(f6.1, 13(2x,f7.3)) ! original

!!!!!----- comm by feng fan ----------------------------------
!!7       format(f6.1, 13(2x,f10.3)) !lijie modify
!!770   format('  time     o3       no       no2      hno3' &
!!            '     hono     pan     onit      h2o2    olet' &
!!            '     olei     par      ald2     isop')
!!!!!----------------------------------------------------------

      return
      end


!!***********************************************************************
!!***********************************************************************
#ifdef KNL_OPT
      subroutine setaircomposition(cbmzobj)
      !use gas_data, only: cbmztype,avogad,cnn,emission
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:), pointer :: cair_mlc
      integer :: l, i2
      cair_mlc => cbmzobj%cair_mlc

!!-----------------------------------------------------------------------
!! set bulk air composition in (molec/cc)
    !dir$ ivdep
    do i2=1,VLEN
      cair_mlc(i2) = avogad*cbmzobj%pr_atm(i2)/(82.056*cbmzobj%te(i2))	! air conc [molec/cc]
!!      write(6,*)' air concentraton = ', cair_mlc, ' molec/cc'
      cbmzobj%o2(i2)       = 0.21*cair_mlc(i2)
      cbmzobj%h2(i2)       = 0.58*1.e-6*cair_mlc(i2)
      cbmzobj%h2o(i2)      = watervapor(cbmzobj%rh(i2), cair_mlc(i2), cbmzobj%te(i2), cbmzobj%pr_atm(i2))
    enddo !i2

!!      write(6,*)'h2o = ', h2o

!!-----------------------------------------------------------------------
!! conversion factor for converting [ppb] to [molec/cc]
!!
      !cbmzobj%ppb = 1.e+9
!!
!!-------------------------------------------------------------
!! converting gas-phase con!! from [ppb] to [molec/cc]
      do l=1, 76
        cnn(:,l) = cnn(:,l)*cair_mlc(:)/1.0e+9
      enddo

!! convert from ppb/hr to molec/cc/s
      do l=1, 76
        emission(:,l) = emission(:,l)*cair_mlc(:)/1.0e+9/3600.
      enddo

      return
      end            
#else
      subroutine setaircomposition
      use gas_data

!!-----------------------------------------------------------------------
!! set bulk air composition in (molec/cc)
      cair_mlc = avogad*pr_atm/(82.056*te)	! air conc [molec/cc]
!!      write(6,*)' air concentraton = ', cair_mlc, ' molec/cc'
      o2       = 0.21*cair_mlc
      h2       = 0.58*1.e-6*cair_mlc
      h2o      = watervapor(rh, cair_mlc, te, pr_atm)

!!      write(6,*)'h2o = ', h2o

!!-----------------------------------------------------------------------
!! conversion factor for converting [ppb] to [molec/cc]
!!
      ppb = 1.e+9
!!
!!-------------------------------------------------------------
!! converting gas-phase con!! from [ppb] to [molec/cc]
      do l=1, 76
        cnn(l) = cnn(l)*cair_mlc/ppb
      enddo

!! convert from ppb/hr to molec/cc/s
      do l=1, 76
        emission(l) = emission(l)*cair_mlc/ppb/3600.
      enddo

      return
      end            

#endif







!!************************************************************************
!! subroutine peroxyrateconstants: calculates parameterized thermal rate 
!!                     constants for the alkylperoxy radical permutation 
!!                     reactions for the entire mechanism.
!! nomenclature:
!! rk_param  = parameterized reaction rate constants (1/s)
!! rk_perox  = individual permutation reaction rate constants (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! 
!! author: rahul a. zaveri
!! date  : june 1998
!!
!!------------------------------------------------------------------------
      !use gas_data,only:nperox,rk_param,te,cnn,aperox,bperox
      subroutine peroxyrateconstants(cbmzobj)
      !use gas_data,only:nperox,cbmztype,rk_param,cnn
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:), pointer :: te
      real, dimension(:,:), pointer :: aperox, bperox
      integer :: i,j
      aperox => cbmzobj%aperox
      bperox => cbmzobj%bperox
      te => cbmzobj%te


!! initialize to zero
      do i = 1, nperox
      rk_param(:,i) = 0.0
      enddo


      do j = 1, nperox
      do i = 1, nperox
      rk_param(:,i) = rk_param(:,i) + aperox(i,j)*exp(bperox(i,j)/te(:))*cnn(:,cnnindex(j))

        enddo
      enddo
      return
      end



!!************************************************************************
!! subroutine loadperoxyparameters: loads thermal rate coefficients 
!!                                  for peroxy-peroxy permutation reactions 
!! nomenclature:
!! aperox  = pre-exponential factor (mole!!c!!s)
!! bperox  = activation energy (-e/r)  (k)
!! 
!! author: rahul a. zaveri
!! date  : june 1998
!!------------------------------------------------------------------------
      subroutine loadperoxyparameters(cbmzobj)
      use gas_data
      type(cbmztype), target :: cbmzobj
      real, dimension(:,:), pointer :: aperox, bperox
      aperox => cbmzobj%aperox
      bperox => cbmzobj%bperox

      aperox(jch3o2,jch3o2)   = 2.5e-13
      aperox(jethp,jethp)     = 6.8e-14
      aperox(jro2,jro2)       = 5.3e-16
      aperox(jc2o3,jc2o3)     = 2.9e-12
      aperox(jano2,jano2)     = 8.0e-12
      aperox(jnap,jnap)       = 1.0e-12
      aperox(jisopp,jisopp)   = 3.1e-14
      aperox(jisopn,jisopn)   = 3.1e-14
      aperox(jisopo2,jisopo2) = 3.1e-14
      aperox(jxo2,jxo2)       = 3.1e-14

      bperox(jch3o2,jch3o2)   = 190.
      bperox(jethp,jethp)     = 0.0
      bperox(jro2,jro2)       = 1980.
      bperox(jc2o3,jc2o3)     = 500.
      bperox(jano2,jano2)     = 0.0
      bperox(jnap,jnap)       = 0.0
      bperox(jisopp,jisopp)   = 1000.
      bperox(jisopn,jisopn)   = 1000.
      bperox(jisopo2,jisopo2) = 1000.
      bperox(jxo2,jxo2)       = 1000.

      do i = 1, nperox
      do j = 1, nperox
        if(i.ne.j)then
          aperox(i,j) = 2.0*sqrt(aperox(i,i)*aperox(j,j))
          bperox(i,j) = 0.5*(bperox(i,i) + bperox(j,j))
        endif
      enddo
      enddo

!! except for
      aperox(jc2o3,jch3o2) = 1.3e-12
      aperox(jch3o2,jc2o3) = 1.3e-12
      bperox(jc2o3,jch3o2) = 640.
      bperox(jch3o2,jc2o3) = 640.

      return
      end










!!***********************************************************************
!! subroutine photoconstants_solar: calculates photochemical rate constants (1/s)
!!
!! input: cos_sza (cosine of solar zenith angle from solarzenithangle.f)
!!        zalt_m (altitude above sea level in meters)
!!
!!------------------------------------------------------------------------
 
      subroutine photoconstants_solar(cbmzobj,factcld)
      use gas_data
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:) :: factcld
      real, parameter :: sza_cut = 89.0 ! cutoff solar zenith angle
      real, parameter :: cos_sza_cut = 0.017452406 ! cos of sza_cut
      integer :: jphoto, i2
      real, dimension(:,:), pointer :: rk_photo
      logical, dimension(VLEN) :: tmask_cond_1,fmask_cond_1
      logical :: has_t, has_f
      rk_photo => cbmzobj%rk_photo

     ! factcld = 1.0
#ifdef KNL_OPT
!      !dir$ simd reduction(.and.:has_t,has_f)
      do i2 = 1, VLEN
         tmask_cond_1(i2) = (cos_sza(i2) .ge. cos_sza_cut) .and. cbmzobj%pmask(i2)
         fmask_cond_1(i2) = .not.(cos_sza(i2) .ge. cos_sza_cut) .and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(tmask_cond_1)
      has_f = any(fmask_cond_1)

      !sequential blocks
      ! code on true branch can be executed without mask, only if the rk_photo
      ! is updated later correctly in false branch.
      if (has_t) then
         !cbmzobj%bmask(:) = tmask_cond_1(:)
         if(mphoto.eq.1)then
           call photoparam1(cbmzobj)
         elseif(mphoto.eq.2)then
           call photoparam2(cbmzobj)
         endif

         !! apply cloudiness correction factor
         do jphoto = 1, nphoto
         do i2 = 1, VLEN
           !if(tmask_cond_1(i2)) then
             rk_photo(i2,jphoto) = rk_photo(i2,jphoto)*factcld(i2)
           !endif !if tmask
         enddo !i2
         enddo
      endif
      if (has_f) then
         do jphoto = 1, nphoto
         do i2 = 1, VLEN
           if(fmask_cond_1(i2)) then
             rk_photo(i2,jphoto) = 0.0
           endif !if fmask
         enddo !i2
         enddo
      endif
#else
      if(cos_sza .ge. cos_sza_cut)then	! daytime

         if(mphoto.eq.1)then
           call photoparam1(cbmzobj)
         elseif(mphoto.eq.2)then
           call photoparam2(cbmzobj)
         endif

!! apply cloudiness correction factor
         do jphoto = 1, nphoto
           rk_photo(jphoto) = rk_photo(jphoto)*factcld
         enddo

      else				! nighttime


         do jphoto = 1, nphoto
           rk_photo(jphoto) = 0.0
         enddo

      endif
#endif
      return
      end












#if 0
!!*****************************************************************
!!      subroutine solar - calculates cosine of zenith angle
!!                         for use in photochemical rate coefficient
!!                         calculations.
!!
!!      nomenclature:
!!  
!!      tmid_sec = time in seconds from greenich noon march 21
!!
!!      cos_za = cosine of zenith angle
!!
!!      rlon   = local longitude, w hemis. values are negative (radians)
!!                   
!!      rlat   = local latitude, s hemis. values are negative (radians)
!!
!!*****************************************************************
 
       subroutine solarzenithangle
       !use gas_data,only:tmid_sec,cos_sza,rlon,rlat
       use gas_data
       implicit none
       real :: tlocal,tdec,sidec,codec,tloc,thou
 
       tlocal=tmid_sec                                                  
       tdec=0.4092797*sin(1.992385e-7*tlocal)                           
       sidec=sin(tdec)                                                  
       codec=cos(tdec)                                                  
       tloc=7.272205e-5*tlocal                                           
       thou=cos(rlon+tloc)                                      
       cos_sza=sin(rlat)*sidec+cos(rlat)*codec*thou         
 
       return
       end
#endif




!!******************************************************************
!!******************************************************************
      subroutine domassbalance
#ifdef KNL_OPT
      !use gas_data,only:cnn,tni,tsi,tcli,dn,ds,dcl,kdms,ksulfhox
      !use gas_data,only:kso2,kh2so4,kno,khno4,kpan,konit,kisopn
      !use gas_data,only:khcl,kn2o5,knap,it
      use gas_data
      implicit none
      integer :: i
      real :: tots,totn,totcl,tnid,tsid,tclid
#else
      use gas_data
      implicit real(a-h,o-z), integer(i-n) 
#endif


!!=================================================================
!! mass balances
!!=================================================================
!!     initialize...
      tots = 0.0
      totn = 0.0
      totcl= 0.0
!!--------------------------------------------------
!! sulfur balance

#ifndef KNL_OPT
      do i=57,67
      tots = tots + cnn(i)
      enddo
      tots = tots + cnn(18) + cnn(1)
 
!!--------------------------------------------------
!! nitrogen balance

      do i=5,10
      totn = totn + cnn(i)
      enddo
      totn = totn + cnn(8)+cnn(32)+cnn(45)+cnn(49)+cnn(55)

!!--------------------------------------------------
!! chlorine balance

      totcl = cnn(3)

!#else
      do i=kdms,ksulfhox
      tots = tots + cnn(i)
      enddo
      tots = tots + cnn(kso2) + cnn(kh2so4)
 
!!--------------------------------------------------
!! nitrogen balance

      do i=kno,khno4
      totn = totn + cnn(i)
      enddo
      totn = totn + cnn(kn2o5)+cnn(kpan)+cnn(konit)+cnn(knap)+cnn(kisopn)

!!--------------------------------------------------
!! chlorine balance

      totcl = cnn(khcl)

!!--------------------------------------------------
!!  initial total elements (n,s,cl)
      if(it.eq.0) then

          if(totn.gt.1.e-30)then
          tni = totn
          tnid= totn
          else
          tni = totn
          tnid= 1.0
          endif

          if(tots.gt.1.e-30)then
          tsi = tots
          tsid= tots
          else
          tsi = tots
          tsid= 1.0
          endif

          if(totcl.gt.1.e-30)then
          tcli = totcl
          tclid= totcl
          else
          tcli = totcl
          tclid= 1.0
          endif

      end if


!!  calculate percent deviation in elemental mass balance
      dn = 100.*(totn-tni)/tnid
      ds = 100.*(tots-tsi)/tsid
      dcl= 100.*(totcl-tcli)/tclid
#endif
 
      return
      end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!************************************************************************
      subroutine photoconstants_fixed(cbmzobj)
      use gas_data
      type(cbmztype), target :: cbmzobj
      real, dimension(:,:), pointer :: rk_photo
      rk_photo => cbmzobj%rk_photo

      rk_photo(:,1)   = 0.0	! no2 + hv    --> no + o(:,3p)
      rk_photo(:,2)   = 0.0	! no3 + hv    --> .89no2 + .89o(:,3p) + .11no
      rk_photo(:,3)   = 0.0	! hono + hv   --> oh + no
      rk_photo(:,4)   = 0.0	! hno3 + hv   --> oh + no2
      rk_photo(:,5)   = 0.0	! hno4 + hv   --> ho2 + no2
      rk_photo(:,6)   = 0.0	! n2o5 + hv   --> no2 + no3
      rk_photo(:,7)   = 0.0	! o3 + hv     --> o(:,3p)
      rk_photo(:,8)   = 0.0	! o3 + hv     --> o(:,1d)
      rk_photo(:,9)   = 0.0	! h2o2 + hv   --> 2oh
      rk_photo(:,10)   = 0.0	! hcho + hv   --> 2ho2 + co
      rk_photo(:,11)  = 0.0	! hcho + hv   --> co
      rk_photo(:,12)  = 0.0	! ch3ooh + hv --> hcho + ho2 + oh
      rk_photo(:,13)  = 0.0	! ethooh + hv --> ald2 + ho2 + oh
      rk_photo(:,14)  = 0.0	! ald2 + hv   --> 
      rk_photo(:,15)  = 0.0	! aone + hv   --> 
      rk_photo(:,16)  = 0.0	! mgly + hv   --> 
      rk_photo(:,17)  = 0.0	! open + hv   --> 
      rk_photo(:,18)  = 0.0	! rooh + hv   -->
      rk_photo(:,19)  = 0.0	! onit + hv   --> 	
      rk_photo(:,20)  = 0.0	! isoprd + hv -->

      return
      end






!!***********************************************************************
!! subroutine photoparam1: calculates photochemical rate constants (1/s)
!!
!! input: cos_sza (cosine of solar zenith angle from solarzenithangle.f)
!!        zalt_m (altitude above sea level in meters)
!!
!!------------------------------------------------------------------------
      subroutine photoparam1(cbmzobj)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      type(cbmztype), target :: cbmzobj
      real, dimension(:,:), pointer :: rk_photo
      integer :: i2
      rk_photo => cbmzobj%rk_photo

    do i2=1,VLEN
      z = min(1.e-4*zalt_m(i2), 1.1)	! zalt in meters

!! no2 + hv --> no + o3p
      alpha0=2.10223e-2-1.6e-3*(z-1.15)**2.
      alpha1=-12.2e-2+alpha0**.5
      beta0=1.258-.16*(z-2.15)**2.
      beta1=-1.2+beta0**.5
      rk_photo(i2,jphoto_no2)    = alpha1*exp(beta1/cos_sza(i2))

!! no3 + hv --> .89 no2 + .89 o3p + .11 no
      rk_photo(i2,jphoto_no3)    = 23.8*rk_photo(i2,jphoto_no2)

!! hono + hv --> oh + no
      rk_photo(i2,jphoto_hono)   = 0.197*rk_photo(i2,jphoto_no2)

!! hno3 + hv --> oh + no2
      rk_photo(i2,jphoto_hno3)   = 3.3e-5*rk_photo(i2,jphoto_no2)

!! hno4 + hv --> ho2 + no2
      rk_photo(i2,jphoto_hno4)   = 5.9e-4*rk_photo(i2,jphoto_no2)

!! n2o5 + hv --> no2 + no3
      rk_photo(i2,jphoto_n2o5)   = 0.0*rk_photo(i2,jphoto_no2)

!! o3 + hv --> o3p
      rk_photo(i2,jphoto_o3a)    = 0.053*rk_photo(i2,jphoto_no2)

!! o3 + hv --> o1d
      alpha0=2.69924e-7-4.0e-8*(z-.375)**2.
      alpha7=-3.25e-4+sqrt(alpha0)
      beta0=4.173-.64*(z-2.0)**2.
      beta7=-3.2+sqrt(beta0)
      rk_photo(i2,jphoto_o3b)    = alpha7*exp(beta7/cos_sza(i2))

!! h2o2 + hv --> 2 oh
      alpha0=2.540534e-9-4.e-11*(z-0.75)**2.
      alpha8=-3.e-5+sqrt(alpha0)
      beta0=.539284-.16*(z-1.75)**2.
      beta8=-1+sqrt(beta0)
      rk_photo(i2,jphoto_h2o2)   = alpha8*exp(beta8/cos_sza(i2))

!! hcho + hv ---> 2 ho2 + co
      alpha0=7.1747e-9-1.6e-9*(z-.75)**2.
      alpha49=-4.e-5+sqrt(alpha0)
      beta0=.7631144-.16*(z-2.)**2.
      beta49=-1.2+sqrt(beta0)
      rk_photo(i2,jphoto_hchoa)  = alpha49*exp(beta49/cos_sza(i2))

!! hcho + hv ---> co
      alpha0=1.693813e-7-4.e-8*(z-.875)**2.
      alpha50=-2.5e-4+sqrt(alpha0)
      beta0=.7631144-.16*(z-1.875)**2.
      beta50=-1.1+sqrt(beta0)
      rk_photo(i2,jphoto_hchob)  = alpha50*exp(beta50/cos_sza(i2))

      rk_photo(i2,jphoto_ch3ooh) = 0.7   *rk_photo(i2,jphoto_h2o2)
      rk_photo(i2,jphoto_ethooh) = 0.7   *rk_photo(i2,jphoto_h2o2)
      rk_photo(i2,jphoto_ald2)   = 4.6e-4*rk_photo(i2,jphoto_no2)
      rk_photo(i2,jphoto_aone)   = 7.8e-5*rk_photo(i2,jphoto_no2)
      rk_photo(i2,jphoto_mgly)   = 9.64  *rk_photo(i2,jphoto_hchoa)
      rk_photo(i2,jphoto_open)   = 9.04  *rk_photo(i2,jphoto_hchoa)
      rk_photo(i2,jphoto_rooh)   = 0.7   *rk_photo(i2,jphoto_h2o2)
      rk_photo(i2,jphoto_onit)   = 1.0e-4*rk_photo(i2,jphoto_no2)
      rk_photo(i2,jphoto_isoprd) = .025  *rk_photo(i2,jphoto_hchob)
    enddo !i2

!!check if rk_photo is negative   (by feng fan)

      do jdum=1,nphoto
        do i2=1,VLEN
         if(rk_photo(i2,jdum)<0)then
           write(*,*)"subroutine photoparam1 rk_photo(",i2,jdum,")<0"
           stop
          end if
        enddo  !i2
      end do 

    
      return
      end




!!***********************************************************************
!! subroutine photoparam2: calculates photochemical rate constants (1/s)
!!
!! input: cos_sza (cosine of solar zenith angle from solarzenithangle.f)
!!        zalt_m (altitude above sea level in meters)
!!
!!------------------------------------------------------------------------

      subroutine photoparam2(cbmzobj)
      !use gas_data,only : cbmztype,cos_sza,zalt_m
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      !implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:,:), pointer :: rk_photo
      real kr(3,nphoto)	! kr(level, species #):

!#ifdef KNL_OPTBAD
      !dir$ attributes align:64 :: krp0,krp1,krp2,krp3,krp4,krp5
      real,parameter,dimension(16,3) :: krp0 = &
        (/-1.0184e-3,4.3945e-3,-1.7863e-4,1.9592e-8,2.1392e-8,0.,0.,-3.2168e-8,-2.1672e-7,-5.4160e-8,-2.6855e-6,0.,0.0,0.0,0.0,0.0,&
        -1.3136e-3,1.59e-2,-2.6339e-4,2.2106e-8,-1.0672e-7,0.,0.,1.6295e-7,-4.73e-7,-7.4493e-7,-5.1681e-6,0.,0.0,0.0,0.0,0.0,&
        -1.3748e-3,2.80132e-2,-3.1944e-4,1.9176e-8,-3.2131e-7,0.,0.,1.6295e-7,-7.6642e-7,-1.7563e-6,-7.9124e-6,0.,0.0,0.0,0.0,0.0/)
      real,parameter,dimension(16,3) :: krp1 = &
        (/1.8542e-2,0.556446,3.2272e-3,-2.8147e-7,-2.0854e-7,0.,0.,2.4588e-6,4.0070e-6,9.4694e-7,4.9102e-5,0.,0.0,0.0,0.0,0.0,&
        2.4948e-2,0.54202,4.6997e-3,-3.4422e-7,1.1698e-6,0.,0.,4.9940e-7,7.4881e-6,8.7149e-6,8.4398e-5,0.,0.0,0.0,0.0,0.0,&
        2.9757e-2,0.51381,6.0983e-3,-3.4083e-7,3.6898e-6,0.,0.,4.9940e-7,1.1717e-5,2.0714e-5,1.258e-4,0.,0.0,0.0,0.0,0.0/)
      real,parameter,dimension(16,3) :: krp2 = &
        (/-9.5368e-3,-0.71996,-8.5989e-4,1.3533e-6,1.6131e-5,0.,0.,-2.6274e-5,1.1189e-5,6.4697e-5,5.7086e-5,0.,0.0,0.0,0.0,0.0,&
        -1.9513e-2,-0.72079,-2.9408e-3,1.8449e-6,1.9044e-5,0.,0.,-2.9055e-5,9.7183e-6,7.1885e-5,2.6478e-5,0.,0.0,0.0,0.0,0.0,&
        -2.8355e-2,-0.68839,-5.2694e-3,2.1560e-6,1.8481e-5,0.,0.,-2.9055e-5,5.3611e-6,6.5668e-5,-2.8767e-5,0.,0.0,0.0,0.0,0.0/)
      real,parameter,dimension(16,3) :: krp3 = &
        (/1.8165e-3,0.34253,-1.8987e-4,-4.2010e-7,-7.2297e-6,0.,0.,1.2005e-4,-6.4306e-6,-3.2594e-5,-4.3040e-5,0.,0.0,0.0,0.0,0.0,&
        6.611e-3,0.34898,7.4996e-4,-6.799e-7,-9.4072e-6,0.,0.,1.8187e-4,-6.4955e-6,-3.9526e-5,-3.4452e-5,0.,0.0,0.0,0.0,0.0,&
        1.1168e-2,0.33448,1.9111e-3,-8.7941e-7,-9.826e-6,0.,0.,1.8187e-4,-4.9358e-6,-3.9386e-5,-1.0505e-5,0.,0.0,0.0,0.0,0.0/)
      real,parameter,dimension(16,3) :: krp4 = &
        (/0.,0.,0.,0.,0.,0.,0.,-7.6216e-5,0.,0.,0.,0.,0.0,0.0,0.0,0.0,&
        0.,0.,0.,0.,0.,0.,0.,-1.5627e-4,0.,0.,0.,0.,0.0,0.0,0.0,0.0,&
        0.,0.,0.,0.,0.,0.,0.,-1.5627e-4,0.,0.,0.,0.,0.0,0.0,0.0,0.0/)
      real,parameter,dimension(16,3) :: krp5 = &
        (/0.,0.,0.,0.,0.,0.,0.,1.4628e-5,0.,0.,0.,0.,0.0,0.0,0.0,0.0,&
        0.,0.,0.,0.,0.,0.,0.,4.5975e-5,0.,0.,0.,0.,0.0,0.0,0.0,0.0,&
        0.,0.,0.,0.,0.,0.,0.,4.5975e-5,0.,0.,0.,0.,0.0,0.0,0.0,0.0/)

      real :: kr2(VLEN,12,3)
      real, dimension(VLEN) :: cz,cz2,cz3,cz4,cz5
      real alpha,a,a1
      integer :: i,j,k,km1
      integer :: i2
      rk_photo => cbmzobj%rk_photo
      cz(:) = cos_sza(:)
      cz2(:) = cz(:)*cz(:)
      cz3(:) = cz(:)*cz2(:)
      cz4(:) = cz2(:)*cz2(:)
      cz5(:) = cz2(:)*cz3(:)

      do j = 1,3
        do i = 1,12
          kr2(:,i,j) = krp0(i,j)+cz(:)*krp1(i,j)+cz2(:)*krp2(i,j)+cz3(:)*krp3(i,j)
          if (i .eq. 8) then
            kr2(:,i,j) = kr2(:,i,j)+krp4(i,j)*cz4(:)+krp5(i,j)*cz5(:)
          endif
          kr2(:,i,j) = max(kr2(:,i,j),0.0)
        end do
      end do

      !!dir$ simd private(alpha,k,km1,a,a1)
      do i2 = 1, VLEN
if (zalt_m(i2)<0)then
        zalt_m(i2) = 1
end if

      if (  zalt_m(i2) .lt.  4000 )then
        alpha = (4000. - zalt_m(i2))/4000
        k = 2
        km1 = 1
      else if (zalt_m(i2) .ge. 8000) then
        alpha = 1
        k = 3
        km1 = 3
      else if (zalt_m(i2) .lt. 8000. .and. zalt_m(i2) .ge. 4000.) then
        alpha = (8000. - zalt_m(i2))/4000.
        k = 3
        km1 = 2
      end if

      a = alpha
      a1= 1. - alpha
      
      do i=1,12
        rk_photo(i2,i) = a1*kr2(i2,i,k) + a*kr2(i2,i,km1)
      end do
    enddo !i2

      do i2 = 1, VLEN
!! n2o5 6, 1
      !rk_photo(i2,jphoto_n2o5)   = 0.0   *rk_photo(i2,jphoto_no2)
      rk_photo(i2,6)   = 0.0 
!! o3a 7, 1      
      rk_photo(i2,7)   = 0.053*rk_photo(i2,1)
!! ald2 14, 1
      rk_photo(i2,14)   = 4.6e-4*rk_photo(i2,1)
!! aone 15, 1
      rk_photo(i2,15)   = 7.8e-5*rk_photo(i2,1)
!! onit 19, 1
      rk_photo(i2,19)   = 1.0e-4*rk_photo(i2,1)
!! ch3ooh 12, 9
      rk_photo(i2,12) = 0.7   *rk_photo(i2,9)
!! ethooh 13, 9
      rk_photo(i2,13) = 0.7   *rk_photo(i2,9)
!! rooh 18, 9
      rk_photo(i2,18)   = 0.7   *rk_photo(i2,9)
!! mgly 16, 10
      rk_photo(i2,16)   = 9.64  *rk_photo(i2,10)
!! open 17, 10
      rk_photo(i2,17)   = 9.04  *rk_photo(i2,10)
!! isoprd 20, 11
      rk_photo(i2,20) = .025  *rk_photo(i2,11)
      enddo !i2

!!check if rk_photo is negative   (by feng fan)

#ifndef KNL_OPT
      do jdum=1,nphoto
      do i2 = 1, VLEN
         if(rk_photo(i2,jdum)<0)then
           write(*,*)"subroutine photoparam2 rk_photo(",i2,jdum,")<0"
           stop
          end if
      enddo !i2
      end do    
#endif

      return
      end





             









!!************************************************************************
!! subroutine gasrateconstants_bio: generates thermal rate coefficients 
!!                   for the selected mechanism
!! nomenclature:
!! rk_bio    = reaction rate constants for hc2 mechanism    (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine gasrateconstants_bio(cbmzobj)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:), pointer :: te

      real, parameter, dimension(9) :: biopar1 = &
        (/2.6e-11,1.2e-14,3.0e-12,1.7e-13,1.7e-13,1.7e-13,6.712e-11,7.44e-17,6.642e-12/)
      real, parameter, dimension(9) :: biopar2 = &
        (/409.,-2013.,-446.,1300.,1300.,1300.,-449.,821.,-175./)
      integer, parameter, dimension(9) :: bioindex = &
        (/1,2,3,11,12,13,18,19,20/)
      integer :: i
      integer :: i2

      te => cbmzobj%te

      !!dir$ simd private(i)
      !!dir$ unroll(9) ? Check
      do i = 1,9
      do i2=1,VLEN
      !if (cbmzobj%bmask(i2)) then
        rk_bio(i2,bioindex(i)) = biopar1(i)*exp(biopar2(i)/te(i2))
      !endif !if iregime
      enddo !i2
      end do
      !dir$ simd 
      do i2=1,VLEN
      !rk_bio(1)  = arr(te,2.6e-11, 409.)
      !rk_bio(2)  = arr(te,1.2e-14, -2013.)
      !rk_bio(3)  = arr(te,3.0e-12, -446.)
      rk_bio(i2,4)  = cbmzobj%rk_photo(i2,20)
      rk_bio(i2,5)  = 3.3e-11
      rk_bio(i2,6)  = 7.0e-18
      rk_bio(i2,7)  = 1.0e-15
      rk_bio(i2,8)  = 4.0e-12
      rk_bio(i2,9)  = 4.0e-12
      rk_bio(i2,10) = 4.0e-12
      !rk_bio(11) = arr(te,1.7e-13, 1300.)
      !rk_bio(12) = arr(te,1.7e-13, 1300.)
      !rk_bio(13) = arr(te,1.7e-13, 1300.)
      rk_bio(i2,14) = rk_param(i2,7)
      rk_bio(i2,15) = rk_param(i2,8)
      rk_bio(i2,16) = rk_param(i2,9)
      rk_bio(i2,17) = 3.563e-11
      !rk_bio(18) = arr(te,6.712e-11,-449.)
      !rk_bio(19) = arr(te,7.44e-17,821.)
      !rk_bio(20) = arr(te,6.642e-12,-175.)
      !endif !if iregime
      enddo !i2

      return
      end





!!************************************************************************
!! subroutine gasrateconstants_com: generates thermal rate coefficients 
!!                   for the selected mechanism
!! nomenclature:
!! rk_com    = reaction rate constants for common mechanism (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! iregime = selected mechanism for the current chemical regime (1-6) 
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!------------------------------------------------------------------------
      subroutine gasrateconstants_com(cbmzobj)
      !use gas_data,only:te,cair_mlc,rk_com,76
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real :: rkexp,rk3,rk0,rki,expo,rk_2ho2
      integer :: i
      integer :: i2

!dir$ attributes align:64 :: rk,rk_com3
      real rk(VLEN,16),rk_com3(VLEN,76)

      real, dimension(:,:), pointer :: rk_com, rk_photo
      real, dimension(:), pointer :: te, cair_mlc
      rk_com => cbmzobj%rk_com
      rk_photo => cbmzobj%rk_photo
      te => cbmzobj%te
      cair_mlc => cbmzobj%cair_mlc

      !arrp1(18) = arrp1(18)*(te**0.667)
      !arrp1(19) = arrp1(19)*(te**2) 
      !TODO: backup finished rk_com and eliminate mask
      !!$omp simd private(i,rkexp)
      do i = 1, 8
      do i2=1,VLEN
        rkexp = exp(-arnm(i)*alog(te(i2)/300.))
        rk(i2,i) = ark0i(i)*cair_mlc(i2)*rkexp
      enddo !i2
      end do
      do i = 9, 16
      do i2=1,VLEN
        rkexp = exp(-arnm(i)*alog(te(i2)/300.))
        rk(i2,i) = ark0i(i)*rkexp
      enddo !i2
      end do
      !!$omp simd
      do i=1,48
      do i2=1,VLEN
        rk_com3(i2,i) = arrp1(i)*exp(arrp2(i)/te(i2))
      enddo !i2
      enddo
      !$omp simd private(rk_2ho2)
      do i2=1,VLEN
      !if(cbmzobj%bmask(i2)) then
      rk_com3(i2,18) = rk_com3(i2,18)*(te(i2)**0.667)
      rk_com3(i2,19) = rk_com3(i2,19)*(te(i2)**2)
      rk3 = rk_com3(i2,40)*cair_mlc(i2)
      rk_com(i2,27) = rk_com3(i2,38) + rk3/(1.0+rk3/rk_com3(i2,39))

      rk_2ho2 = rk_com3(i2,41)+rk_com3(i2,42)*cair_mlc(i2)
      rk_com(i2,31) = rk_2ho2
      rk_com(i2,32) = rk_2ho2*rk_com3(i2,43)
      !endif ! bmask
      enddo !i2

      do i=37,44
      !$omp simd private(rk0,rki,expo)
      do i2=1,VLEN
        rk0 = rk(i2,i-36)
        rki = rk(i2,i-28)
        expo = 1./(1. + (alog10(rk0/rki))**2)
        rk_com3(i2,i) = (rk0*rki/(rk0+rki))*exp(expo*alog(.6))
      enddo !i2
      enddo

      !$omp simd private(i)
      do i2=1,VLEN
      !if(cbmzobj%bmask(i2)) then
      !!!!dir$ vector nontemporal
      do i=1,44
        rk_com(i2,rkindex(i)) = rk_com3(i2,i)
      enddo
      rk_com(i2,1) = rk_photo(i2,1)                                                         
      rk_com(i2,2) = rk_photo(i2,2)                                                         
      rk_com(i2,3) = rk_photo(i2,3)                                                        
      rk_com(i2,4) = rk_photo(i2,4)                                                        
      rk_com(i2,5) = rk_photo(i2,5)                                                       
      rk_com(i2,6) = rk_photo(i2,6) ! 0.0                                                    
      rk_com(i2,7) = rk_photo(i2,7)                                                         
      rk_com(i2,8) = rk_photo(i2,8)                                                         
      rk_com(i2,9) = rk_photo(i2,9)

      rk_com(i2,13) = cair_mlc(i2)*6.e-34*exp(-2.3*alog(te(i2)/300.))                                                        
      rk_com(i2,36) = rk_com(i2,34)*rk_com(i2,36)
      rk_com(i2,43) = rk_com(i2,39)*rk_com(i2,43)
      rk_com(i2,44) = 1.5e-13 * (1.+8.18e-23*te(i2)*cair_mlc(i2))       
            
      !rk_com(i2,46) = te(i2)**.667*ARR(te(i2),2.8e-14, -1575.)                                
      !rk_com(i2,47) = te(i2)**2*ARR(te(i2),1.5e-17, -492.) 

      rk_com(i2,49) = rk_photo(i2,10)                                                      
      rk_com(i2,50) = rk_photo(i2,11)                                                      
      rk_com(i2,53) = rk_photo(i2,12)                                                    
      rk_com(i2,54) = rk_photo(i2,13)                                                    
      rk_com(i2,63) = rk_param(i2,1)                                                 
      rk_com(i2,64) = rk_param(i2,2)                                                  
      rk_com(i2,66) = rk_photo(i2,14)                                               

      rk_com(i2,70) = rk_com(i2,69)*rk_com(i2,70)
      rk_com(i2,74) = rk_param(i2,4)  
      !endif !if bmask
      enddo !i2

      return
      end





!!!! to comment all lines in this subroutine and add a new fron lijie
!!!             for heterogeneous chemistry  !!!
!!************************************************************************
!! subroutine gasrateconstants_het: generates thermal rate coefficients 
!!                   for the selected mechanism
!! 
!! rk_het    = reaction rate constants for heterogeneous rxns (1/s)
!! 
!! 
!! revised by li jie
!!------------------------------------------------------------------------
      subroutine gasrateconstants_het
      use gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! to comment all lines in this subroutine and add a new form (by lijie)

#ifndef VEC_OPT
         rk_het(1) = rk_het(1)  !1  n2o5 + aso4 -> 2hno3 1/s
         rk_het(2) = rk_het(2)  !2  no2  + bc   -> 0.5hono+0.5hno3
         rk_het(3) = rk_het(3)  !3  no3  + aso4 -> hno3
         rk_het(4) = rk_het(4)  !4  ho2  + aso4 -> 0.5h2o2
         rk_het(5) = rk_het(5)  !5  hcho + aso4 -> products
         rk_het(6) = rk_het(6)  !6  oh   + aso4 -> products
         rk_het(7) = rk_het(7)  !7  o3   + bc   -> products
         rk_het(8) = rk_het(8)  !8  no2  + bc   -> hono
         rk_het(9) = rk_het(9)  !9  hno3 + bc   -> no2
         rk_het(10) = rk_het(10)    ! 10  n2o5 + bc   -> 2hno3
         rk_het(11) = rk_het(11)    ! 11  o3   + dust -> products
         rk_het(12) = rk_het(12)    ! 12  hno3 + dust -> ano3 + products
         rk_het(13) = rk_het(13)    ! 13  no2  + dust -> 0.5hono + 0.5hno3
         rk_het(14) = rk_het(14)    ! 14  no3  + dust -> hno3
         rk_het(15) = rk_het(15)    ! 15  n2o5 + dust -> 2hno3
         rk_het(16) = rk_het(16)    ! 16  oh   + dust -> products
         rk_het(17) = rk_het(17)    ! 17  ho2  + dust -> 0.5h2o2
         rk_het(18) = rk_het(18)    ! 18  h2o2 + dust -> products
         rk_het(19) = rk_het(19)    ! 19  so2  + dust -> aso4 
         rk_het(20) = rk_het(20)    ! 20  ch3cooh + dust -> products
         rk_het(21) = rk_het(21)    ! 21  ch3oh   + dust -> products
         rk_het(22) = rk_het(22)    ! 22  hcho    + dust -> products
         rk_het(23) = rk_het(23)    ! 23  n2o5 + ssa  -> 2hno3
         rk_het(24) = rk_het(24)    ! 24  no3  + ssa  -> hno3
         rk_het(25) = rk_het(24)    ! 25  ho2  + ssa  -> 0.5hono
         rk_het(26) = rk_het(26)    ! 26  so2  + ssa  -> aso4
         rk_het(27) = rk_het(27)    ! 27  no3  + ssa  -> ano3
         rk_het(28) = rk_het(28)    ! 28  hno3 + ssa  -> ano3
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the following values are temperary values only used for test. the above 28 lines are written by li jie
!         rk_het(1) = 1.1451681e-27  
!         rk_het(2) = 9.2583208e-29
!         rk_het(3) = 2.3924208e-27
!         rk_het(4) = 2.7327207e-25
!         rk_het(5) = 2.5221690e-26
!         rk_het(6) = 3.0459142e-25
!         rk_het(7) = 5.6911365e-30
!         rk_het(8) = 3.0552454e-28
!         rk_het(9) = 6.3289386e-27
!         rk_het(10) = 3.0211289e-27
!         rk_het(11) = 0.0
!         rk_het(12) = 0.0
!         rk_het(13) = 0.0
!         rk_het(14) = 0.0
!         rk_het(15) = 0.0
!         rk_het(16) = 0.0
!         rk_het(17) = 0.0
!         rk_het(18) = 0.0
!         rk_het(19) = 0.0
!         rk_het(20) = 0.0
!         rk_het(21) = 0.0
!         rk_het(22) = 0.0
!         rk_het(23) = 0.0
!         rk_het(24) = 0.0
!         rk_het(25) = 0.0
!         rk_het(26) = 0.0
!         rk_het(27) = 0.0
!         rk_het(28) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return
      end





!!************************************************************************
!! subroutine gasrateconstants_mar: generates thermal rate coefficients 
!!                   for the selected mechanism
!! nomenclature:
!! rk_mar    = reaction rate constants for marine mechanism (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine gasrateconstants_mar(cbmzobj)
      !use gas_data, only : cbmztype, rk_mar
      use gas_data
      implicit none
      type(cbmztype), target :: cbmzobj
      real :: rk_tot_num, rk_tot_den, rk_tot, babs, badd
      real, dimension(:), pointer :: te
      integer :: i2
      te => cbmzobj%te

!! abstraction reaction
!! hynes et al. (1986)
      !dir$ simd private(rk_tot_num, rk_tot_den, rk_tot, babs, badd)
      do i2=1,VLEN
      !if (cbmzobj%bmask(i2)) then
      rk_tot_num =       te(i2) * exp(-234./te(i2)) +  &
                   8.46e-10 * exp(7230./te(i2)) +  &
                   2.68e-10 * exp(7810./te(i2))
      rk_tot_den = 1.04e+11 * te(i2) + 88.1 * exp(7460./te(i2))
      rk_tot	 = rk_tot_num/rk_tot_den

      rk_mar(i2,1)   = 9.60e-12 * exp(-234./te(i2)) ! ch3sch3 + oh --> ch3sch2
      babs       = rk_mar(i2,1)/rk_tot
      badd	 = 1. - babs
      rk_mar(i2,2)   = 1.40e-13 * exp(500./te(i2))  ! ch3sch3 + no3 --> 
      rk_mar(i2,3)   = 1.26e-11 * exp(409./te(i2))  ! ch3sch3 + o3p --> 

!! addition reaction
      rk_mar(i2,4)   = badd*rk_tot		     ! ch3sch3 + oh --> ch3s(oh)ch3
      rk_mar(i2,5)   = 8.0e-12
      rk_mar(i2,6)   = 1.8e-13
      rk_mar(i2,7)   = 5.8e-11
      rk_mar(i2,8)   = 1.0e-14
      rk_mar(i2,9)   = 5.0e-12
      rk_mar(i2,10)  = 1.8e-13
      rk_mar(i2,11)  = 1.0e-15
      rk_mar(i2,12)  = 1.0e-13
      rk_mar(i2,13)  = 1.0e-15
      rk_mar(i2,14)  = 1.6e-11
      rk_mar(i2,15)  = 1.0e-13
      rk_mar(i2,16)  = 2.5e+13 * exp(-8686./te(i2)) ! ch3so2 --> so2 + ch3o2
      rk_mar(i2,17)  = 1.0e-14
      rk_mar(i2,18)  = 5.0e-15
      rk_mar(i2,19)  = 2.5e-13
      rk_mar(i2,20)  = 2.5e-13
      rk_mar(i2,21)  = 5.0e-11
      rk_mar(i2,22)  = 2.6e-18
      rk_mar(i2,23)  = 3.3
      rk_mar(i2,24)  = 1.0e-11
      rk_mar(i2,25)  = 5.5e-12
      rk_mar(i2,26)  = 2.0e+17 * exp(-12626./te(i2)) ! ch3so3 --> h2so4 + ch3o2
      rk_mar(i2,27)  = 3.0e-15
      rk_mar(i2,28)  = 3.0e-15
      rk_mar(i2,29)  = 5.0e-11
      rk_mar(i2,30)  = 1.6e-15

      rk_mar(i2,31)  = 2.5e-13	! ch3sch2oo + ch3so2 --> ch3so3 + ch3so2
      rk_mar(i2,32)  = 8.6e-14	! 2ch3sch2oo --> .15mtf + 1.85ch3so2

!! dry deposition
      rk_mar(i2,33)  = 0.0 ! 2.0e-5	! 1/s
      rk_mar(i2,34)  = 0.0 ! 2.0e-5
      rk_mar(i2,35)  = 0.0 ! 2.0e-5
      !endif !if iregime
      enddo !i2

      return
      end
!!************************************************************************
!! subroutine gasrateconstants_urb: generates thermal rate coefficients 
!!                   for the selected mechanism
!! nomenclature:
!! rk_urb    = reaction rate constants for hc1 mechanism    (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine gasrateconstants_urb(cbmzobj)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      real, dimension(:,:), pointer :: rk_urb, rk_photo
      real, dimension(:), pointer :: cair_mlc, te
      type(cbmztype), target :: cbmzobj
#ifdef VEC_OPT
      real :: rk0,rnn,rki,rmm
      integer :: i
      integer :: i2
      !real :: rk_urb3(18)
      real, parameter,dimension(18) :: p1 = &
        (/5.3e-18,1.4e-12,1.2e-14,4.2e-15,8.9e-16,5.8e-12,2.9e-11,3.1e-13, &
          2.1e-12,1.7e-11,5.4e-17,3.8e-12,1.6e-11,1.7e-13,1.2e-13,1.7e-13, &
          1.7e-13,3.3e-12/)
      real, parameter,dimension(18) :: p2 = &
        (/-230.,-1900.,-2630.,-1800.,-392.,478.,255.,-1010.,322.,116., &
          -500.,200.,-540.,1300.,1300.,1300.,1300.,-2880./)
      integer, parameter, dimension(18) :: urbindex = &
        (/3,6,7,9,10,11,12,13,15,16,23,25,26,36,37,38,39,46/)
#endif
      rk_urb => cbmzobj%rk_urb
      rk_photo => cbmzobj%rk_photo
      cair_mlc => cbmzobj%cair_mlc
      te => cbmzobj%te

      !dir$ simd private(i)
      do i2 = 1,VLEN
      if (cbmzobj%bmask(i2)) then
      !!dir$ simd
      !!dir$ unroll(18) ? Check
      do i = 1,18
        rk_urb(i2,urbindex(i)) = p1(i)*exp(p2(i)/te(i2))
      enddo
      rk_urb(i2,3) = rk_urb(i2,3)*te(i2)*te(i2)
      rk_urb(i2,1) = 8.1e-13                                                           
      rk_urb(i2,5) = 1.7e-11                                                           
      rk_urb(i2,14) = 2.5e-12                                                          
      rk_urb(i2,17) = 8.1e-12                                                          
      rk_urb(i2,18) = 4.1e-11                                                          
      rk_urb(i2,19) = 2.2e-11                                                          
      rk_urb(i2,20) = 1.4e-11                                                          
      rk_urb(i2,21) = 3.0e-11                                                          
      rk_urb(i2,28) = 4.0e-12                                                          
      rk_urb(i2,29) = 4.0e-12                                                          
      rk_urb(i2,30) = 4.0e-12                                                          
      rk_urb(i2,31) = 4.0e-12                                                          
      rk_urb(i2,32) = 2.5e-12                                                          
      rk_urb(i2,33) = 1.2e-12                                                          
      rk_urb(i2,34) = 4.0e-12                                                          
      rk_urb(i2,35) = 2.5e-12                                                          
      rk_urb(i2,44) = 1.0e-11   
      rk_urb(i2,45) = 7.0e-17   ! par+no3
      rk_urb(i2,47) = 7.0e-17   ! tol+no3
      rk_urb(i2,48) = 3.0e-17   ! anoe+no3
      rk_urb(i2,49) = 5.0e-16   ! xyl+no3



      rk_urb(i2,2) = rk_photo(i2,15)                                                
      rk_urb(i2,4) = rk_photo(i2,16)                                                  
      rk_urb(i2,22) = rk_photo(i2,17)                                                 
      rk_urb(i2,24) = rk_photo(i2,18)                                                    
      rk_urb(i2,27) = rk_photo(i2,19)                                               
      rk_urb(i2,40) = rk_param(i2,3)                                                   
      rk_urb(i2,41) = rk_param(i2,5)                                                  
      rk_urb(i2,42) = rk_param(i2,6)                                                   
      rk_urb(i2,43) = rk_param(i2,10)                                                   

      rk0 = 1.0e-28
      rnn = 0.8
      rki = 8.8e-12
      rmm = 0.0                                        
      rk_urb(i2,8) = troe(cair_mlc(i2),te(i2),rk0,rnn,rki,rmm)    
                                 
      endif !if bmask
      enddo !i2

!      rk_urb(3) = te**2*arr(te,5.3e-18, -230.)                                         
!      rk_urb(6) = arr(te,1.4e-12, -1900.)                                              
!      rk_urb(7) = arr(te,1.2e-14, -2630.)      
!      rk_urb(9) = arr(te,4.2e-15, -1800.)                                              
!      rk_urb(10) = arr(te,8.9e-16, -392.)                                              
!      rk_urb(11) = arr(te,5.8e-12, 478.)                                               
!      rk_urb(12) = arr(te,2.9e-11, 255.)                                               
!      rk_urb(13) = arr(te,3.1e-13, -1010.)                                             
!      rk_urb(15) = arr(te,2.1e-12, 322.)                                               
!      rk_urb(16) = arr(te,1.7e-11, 116.)                                               
!      rk_urb(23) = arr(te,5.4e-17, -500.)                                              
!      rk_urb(25) = arr(te,3.8e-12, 200.)                                               
!      rk_urb(26) = arr(te,1.6e-11, -540.)                                              
!      rk_urb(36) = arr(te,1.7e-13, 1300.)                                              
!      rk_urb(37) = arr(te,1.2e-13, 1300.)                                              
!      rk_urb(38) = arr(te,1.7e-13, 1300.)                                              
!      rk_urb(39) = arr(te,1.7e-13, 1300.)                                              
!      rk_urb(46) = arr(te,3.3e-12,-2880.) ! eth+no3=no2+xo2+2hcho
!
! move from rk_com to rk_urb by feng fan


 

      return
      end





!!!-----------------------------------------------------------------
!!!------------------------------------------------------------------
      subroutine gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      use gas_data !,only:cbmztype
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      real :: s(VLEN,76),rk_s
      !implicit none
      type(cbmztype), target :: cbmzobj

      real, intent(out) :: r_com(VLEN,76),r1_com(VLEN,76),r2_com(VLEN,76)
!dir$ assume_aligned s:64,r_com:64,r1_com:64,r2_com:64
      real, pointer :: rk_com(:,:), cair_mlc(:),h2o(:),h2(:)
      integer :: i, i2
      rk_com => cbmzobj%rk_com
      cair_mlc => cbmzobj%cair_mlc
      h2o => cbmzobj%h2o
      h2 => cbmzobj%h2

!#ifdef VEC_OPTBAD

!dir$ assume_aligned r_com:64,r1_com:64,r2_com:64,s:64
      do i = 1,76
      !$omp simd private(rk_s)
    do i2 =1,VLEN
    !if (cbmzobj%pmask(i2)) then
        rk_s = rk_com(i2,i)*s(i2,sindex(i))
        if (sindex2(i) == 33 ) then
          r_com(i2,i) = rk_s
          r1_com(i2,i) = rk_com(i2,i)
        else 
          r_com(i2,i) = rk_s*s(i2,sindex2(i))
          r1_com(i2,i) = rk_s
          r2_com(i2,i) = rk_com(i2,i)*s(i2,sindex2(i))
        endif
    !endif !if pmask
    enddo !i2
      enddo
!by WH 2016/11/15
    !$omp simd
    do i2 =1,VLEN
    !if (cbmzobj%pmask(i2)) then
      r2_com(i2,32) = 0.
      r2_com(i2,31) = 0.
      r2_com(i2,40) = 0.
      !r_com(i2,32)  = 1.!rk_com(i2,32)*s(i2,15)*s(i2,15) 
      r_com(i2,10) = r_com(i2,10)*.21*cair_mlc(i2) 
      r1_com(i2,10) = rk_com(i2,10)*.21*cair_mlc(i2)
      r_com(i2,11) = r_com(i2,11)*.79*cair_mlc(i2)
      r1_com(i2,11) = rk_com(i2,11)*.79*cair_mlc(i2) 
      r_com(i2,12) = r_com(i2,12)*h2o(i2)
      r1_com(i2,12) = rk_com(i2,12)*h2o(i2)  
      r_com(i2,13) = r_com(i2,13)*.21*cair_mlc(i2)
      r1_com(i2,13) = rk_com(i2,13)*.21*cair_mlc(i2)  
      r_com(i2,22) = r_com(i2,22)*h2(i2) 
      r1_com(i2,22) = rk_com(i2,22)*h2(i2) 
      r_com(i2,32) = r_com(i2,32)*h2o(i2)
      r1_com(i2,32) = r1_com(i2,32)*h2o(i2)  
      r_com(i2,42) = r_com(i2,42)*h2o(i2)
      r1_com(i2,42) = r1_com(i2,42)*h2o(i2)   
    !endif !if pmask
    enddo !i2


      return
      end


!!!-----------------------------------------------------------------
!!!------------------------------------------------------------------
      subroutine gasrates_bio(s)
      use gas_data
!dir$ assume_aligned s:64
      real :: s(VLEN,76)
      !subroutine gasrates_bio(cbmzobj,s,r_bio,r1_bio,r2_bio)
      !use gas_data
      !implicit none
      !type(cbmztype), target :: cbmzobj
      !real :: s(76),r_bio(80),r1_bio(80),r2_bio(80)
      integer :: iisop,iisoprd,iisopp,iisopn,iisopo2,iterp
      integer :: i2
      iisop = 56
      iisoprd = 57
      iisopp = 58
      iisopn = 59
      iisopo2= 60
      iterp = 61
!      integer, parameter, dimension(20) :: bindex1 = &
!        (/56,56,56,57,57,57,57,58,59,60,58,59,60,58,59,60,61,61,61,61/)
!      integer, parameter, dimension(20) :: bindex2 = &
!        (/14,11, 7,-1,14,11, 7, 5, 5, 5,15,15,15,-1,-1,-1,12,14,11, 7/)
!      integer :: i
!      real :: res
!      do i = 1, 20
!        res = rk_bio(i)*s(bindex1(i))
!        if (bindex2(i) == -1) then
!          r_bio(i) = res
!          r1_bio(i) = rk_bio(i)
!        else
!          r_bio(i) = res*s(bindex2(i))
!          r1_bio(i) = res
!          r2_bio(i) = rk_bio(i)*s(bindex2(i))
!        endif 
!      end do  

    !dir$ simd
    do i2 =1,VLEN
    !if (cbmzobj%bmask(i2)) then
! 56,14
      r_bio(i2,1)  = rk_bio(i2,1)*s(i2,iisop)*s(i2,ioh)
      r1_bio(i2,1)  = rk_bio(i2,1)*s(i2,iisop)
      r2_bio(i2,1)  = rk_bio(i2,1)*s(i2,ioh)

! 56,11
      r_bio(i2,2)  = rk_bio(i2,2)*s(i2,iisop)*s(i2,io3)
      r1_bio(i2,2)  = rk_bio(i2,2)*s(i2,iisop)
      r2_bio(i2,2)  = rk_bio(i2,2)*s(i2,io3)

! 56,7
      r_bio(i2,3)  = rk_bio(i2,3)*s(i2,iisop)*s(i2,ino3)
      r1_bio(i2,3)  = rk_bio(i2,3)*s(i2,iisop)
      r2_bio(i2,3)  = rk_bio(i2,3)*s(i2,ino3)

! 57
      r_bio(i2,4)  = rk_bio(i2,4) *s(i2,iisoprd)
      r1_bio(i2,4)  = rk_bio(i2,4)

! 57,14
      r_bio(i2,5)  = rk_bio(i2,5)*s(i2,iisoprd)*s(i2,ioh)
      r1_bio(i2,5)  = rk_bio(i2,5)*s(i2,iisoprd)
      r2_bio(i2,5)  = rk_bio(i2,5)*s(i2,ioh)

! 57,11
      r_bio(i2,6)  = rk_bio(i2,6)*s(i2,iisoprd)*s(i2,io3)
      r1_bio(i2,6)  = rk_bio(i2,6)*s(i2,iisoprd)
      r2_bio(i2,6)  = rk_bio(i2,6)*s(i2,io3)

! 57,7
      r_bio(i2,7)  = rk_bio(i2,7)*s(i2,iisoprd)*s(i2,ino3)
      r1_bio(i2,7)  = rk_bio(i2,7)*s(i2,iisoprd)
      r2_bio(i2,7)  = rk_bio(i2,7)*s(i2,ino3)

! 58,5
      r_bio(i2,8)  = rk_bio(i2,8)*s(i2,iisopp)*s(i2,ino)
      r1_bio(i2,8)  = rk_bio(i2,8)*s(i2,iisopp)
      r2_bio(i2,8)  = rk_bio(i2,8)*s(i2,ino)

! 59,5
      r_bio(i2,9)  = rk_bio(i2,9)*s(i2,iisopn)*s(i2,ino)
      r1_bio(i2,9)  = rk_bio(i2,9)*s(i2,iisopn)
      r2_bio(i2,9)  = rk_bio(i2,9)*s(i2,ino)

! 60 ,5
      r_bio(i2,10) = rk_bio(i2,10)*s(i2,iisopo2)*s(i2,ino)
      r1_bio(i2,10) = rk_bio(i2,10)*s(i2,iisopo2)
      r2_bio(i2,10) = rk_bio(i2,10)*s(i2,ino)
! 58,15
      r_bio(i2,11) = rk_bio(i2,11)*s(i2,iisopp)*s(i2,iho2)
      r1_bio(i2,11) = rk_bio(i2,11)*s(i2,iisopp)
      r2_bio(i2,11) = rk_bio(i2,11)*s(i2,iho2)
! 59,15
      r_bio(i2,12) = rk_bio(i2,12)*s(i2,iisopn)*s(i2,iho2)
      r1_bio(i2,12) = rk_bio(i2,12)*s(i2,iisopn)
      r2_bio(i2,12) = rk_bio(i2,12)*s(i2,iho2)
! 60,15
      r_bio(i2,13) = rk_bio(i2,13)*s(i2,iisopo2)*s(i2,iho2)
      r1_bio(i2,13) = rk_bio(i2,13)*s(i2,iisopo2)
      r2_bio(i2,13) = rk_bio(i2,13)*s(i2,iho2)
! 58
      r_bio(i2,14) = rk_bio(i2,14)*s(i2,iisopp)
      r1_bio(i2,14) = rk_bio(i2,14)
! 59
      r_bio(i2,15) = rk_bio(i2,15)*s(i2,iisopn)
      r1_bio(i2,15) = rk_bio(i2,15)
! 60
      r_bio(i2,16) = rk_bio(i2,16)*s(i2,iisopo2)
      r1_bio(i2,16) = rk_bio(i2,16)

! li jie added r_bio(i2,17)(i2,18)(i2,19)(i2,20)
! 61,12
      r_bio(i2,17) = rk_bio(i2,17)*s(i2,iterp)*s(i2,io1d)
      r1_bio(i2,17) = rk_bio(i2,17)*s(i2,iterp)
      r2_bio(i2,17) = rk_bio(i2,17)*s(i2,io1d)
! 61,14
      r_bio(i2,18) = rk_bio(i2,18)*s(i2,iterp)*s(i2,ioh)
      r1_bio(i2,18) = rk_bio(i2,18)*s(i2,iterp)
      r2_bio(i2,18) = rk_bio(i2,18)*s(i2,ioh)
! 61,11
      r_bio(i2,19) = rk_bio(i2,19)*s(i2,iterp)*s(i2,io3)
      r1_bio(i2,19) = rk_bio(i2,19)*s(i2,iterp)
      r2_bio(i2,19) = rk_bio(i2,19)*s(i2,io3)
! 61,7
      r_bio(i2,20) = rk_bio(i2,20)*s(i2,iterp)*s(i2,ino3)
      r1_bio(i2,20) = rk_bio(i2,20)*s(i2,iterp)
      r2_bio(i2,20) = rk_bio(i2,20)*s(i2,ino3)
    !endif !if bmask
    enddo !i2

      return
      end










!!!-----------------------------------------------------------------
!!!------------------------------------------------------------------
      subroutine gasrates_het(cbmzobj,s,r_het)

      use gas_data
      type(cbmztype),target :: cbmzobj
      real :: s(VLEN,76), r_het(VLEN,76)
!dir$ assume_aligned s:64,r_het:64

      real, dimension(:,:), pointer :: rk_het
      rk_het => cbmzobj%rk_het


!! heterogeneous chemistry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! to comment all lines in this subroutine and add a new form (by li jie)

      do i=1,28
       do i2=1,VLEN
        !if (cbmzobj%pmask(i2)) then
        r_het(:,i) = rk_het(:,i)*s(:,hetindex(i))
        !r1_het(i) = rk_het(i)
        !endif ! if pmask
       enddo !i2
      end do

      return
      end

!!!-----------------------------------------------------------------
!!!------------------------------------------------------------------

      subroutine gasrates_mar(cbmzobj,s)
      use gas_data
      dimension s(VLEN,76)
!dir$ assume_aligned s:64
      type(cbmztype), target :: cbmzobj
      integer :: i2

      !dir$ simd
      do i2=1,VLEN
      !if (cbmzobj%bmask(i2)) then
      r_mar(i2,1) = rk_mar(i2,1)*s(i2,idms)*s(i2,ioh)
      r1_mar(i2,1) = rk_mar(i2,1)*s(i2,idms)
      r2_mar(i2,1) = rk_mar(i2,1)*s(i2,ioh)

      r_mar(i2,2) = rk_mar(i2,2)*s(i2,idms)*s(i2,ino3)
      r1_mar(i2,2) = rk_mar(i2,2)*s(i2,idms)
      r2_mar(i2,2) = rk_mar(i2,2)*s(i2,ino3)

      r_mar(i2,3) = rk_mar(i2,3)*s(i2,idms)*s(i2,io3p)
      r1_mar(i2,3) = rk_mar(i2,3)*s(i2,idms)
      r2_mar(i2,3) = rk_mar(i2,3)*s(i2,io3p)

      r_mar(i2,4) = rk_mar(i2,4)*s(i2,idms)*s(i2,ioh)
      r1_mar(i2,4) = rk_mar(i2,4)*s(i2,idms)
      r2_mar(i2,4) = rk_mar(i2,4)*s(i2,ioh)

      r_mar(i2,5) = rk_mar(i2,5)*s(i2,ich3sch2oo)*s(i2,ino)
      r1_mar(i2,5) = rk_mar(i2,5)*s(i2,ich3sch2oo)
      r2_mar(i2,5) = rk_mar(i2,5)*s(i2,ino)

      r_mar(i2,6) = rk_mar(i2,6)*s(i2,ich3sch2oo)*s(i2,ich3o2)
      r1_mar(i2,6) = rk_mar(i2,6)*s(i2,ich3sch2oo)
      r2_mar(i2,6) = rk_mar(i2,6)*s(i2,ich3o2)

      r_mar(i2,7) = rk_mar(i2,7)*s(i2,idmso)*s(i2,ioh)
      r1_mar(i2,7)  = rk_mar(i2,7)*s(i2,idmso)
      r2_mar(i2,7) = rk_mar(i2,7)*s(i2,ioh)

      r_mar(i2,8) = rk_mar(i2,8)*s(i2,idmso2)*s(i2,ioh)
      r1_mar(i2,8) = rk_mar(i2,8)*s(i2,idmso2)
      r2_mar(i2,8) = rk_mar(i2,8)*s(i2,ioh)

      r_mar(i2,9) = rk_mar(i2,9)*s(i2,ich3so2ch2oo)*s(i2,ino)
      r1_mar(i2,9) = rk_mar(i2,9)*s(i2,ich3so2ch2oo)
      r2_mar(i2,9) = rk_mar(i2,9)*s(i2,ino)

      r_mar(i2,10) = rk_mar(i2,10)*s(i2,ich3so2ch2oo)*s(i2,ich3o2)
      r1_mar(i2,10) = rk_mar(i2,10)*s(i2,ich3so2ch2oo)
      r2_mar(i2,10) = rk_mar(i2,10)*s(i2,ich3o2)

      r_mar(i2,11) = rk_mar(i2,11)*s(i2,ich3so2h)*s(i2,iho2)
      r1_mar(i2,11) = rk_mar(i2,11)*s(i2,ich3so2h)
      r2_mar(i2,11) = rk_mar(i2,11)*s(i2,iho2)

      r_mar(i2,12) = rk_mar(i2,12)*s(i2,ich3so2h)*s(i2,ino3)
      r1_mar(i2,12) = rk_mar(i2,12)*s(i2,ich3so2h)
      r2_mar(i2,12) = rk_mar(i2,12)*s(i2,ino3)

      r_mar(i2,13) = rk_mar(i2,13)*s(i2,ich3so2h)*s(i2,ich3o2)
      r1_mar(i2,13) = rk_mar(i2,13)*s(i2,ich3so2h)
      r2_mar(i2,13) = rk_mar(i2,13)*s(i2,ich3o2)

      r_mar(i2,14) = rk_mar(i2,14)*s(i2,ich3so2h)*s(i2,ioh)
      r1_mar(i2,14) = rk_mar(i2,14)*s(i2,ich3so2h)
      r2_mar(i2,14) = rk_mar(i2,14)*s(i2,ioh)

      r_mar(i2,15) = rk_mar(i2,15)*s(i2,ich3so2h)*s(i2,ich3so3)
      r1_mar(i2,15) = rk_mar(i2,15)*s(i2,ich3so2h)
      r2_mar(i2,15) = rk_mar(i2,15)*s(i2,ich3so3)

      r_mar(i2,16) = rk_mar(i2,16)*s(i2,ich3so2)
      r1_mar(i2,16) = rk_mar(i2,16) 
      

      r_mar(i2,17) = rk_mar(i2,17)*s(i2,ich3so2)*s(i2,ino2)
      r1_mar(i2,17) = rk_mar(i2,17)*s(i2,ich3so2)
      r2_mar(i2,17) = rk_mar(i2,17)*s(i2,ino2)

      r_mar(i2,18) = rk_mar(i2,18)*s(i2,ich3so2)*s(i2,io3)
      r1_mar(i2,18) = rk_mar(i2,18)*s(i2,ich3so2)
      r2_mar(i2,18) = rk_mar(i2,18)*s(i2,io3)

      r_mar(i2,19) = rk_mar(i2,19)*s(i2,ich3so2)*s(i2,iho2)
      r1_mar(i2,19) = rk_mar(i2,19)*s(i2,ich3so2)
      r2_mar(i2,19) = rk_mar(i2,19)*s(i2,iho2)

      r_mar(i2,20) = rk_mar(i2,20)*s(i2,ich3so2)*s(i2,ich3o2)
      r1_mar(i2,20) = rk_mar(i2,20)*s(i2,ich3so2)
      r2_mar(i2,20) = rk_mar(i2,20)*s(i2,ich3o2)

      r_mar(i2,21) = rk_mar(i2,21)*s(i2,ich3so2)*s(i2,ioh)
      r1_mar(i2,21) = rk_mar(i2,21)*s(i2,ich3so2)
      r2_mar(i2,21) = rk_mar(i2,21)*s(i2,ioh)

      r_mar(i2,22) = rk_mar(i2,22)*s(i2,ich3so2)*cbmzobj%o2(i2)
      r1_mar(i2,22) = rk_mar(i2,22)*cbmzobj%o2(i2)

      r_mar(i2,23) = rk_mar(i2,23)*s(i2,ich3so2oo)
      r1_mar(i2,23) = rk_mar(i2,23)  

      r_mar(i2,24) = rk_mar(i2,24)*s(i2,ich3so2oo)*s(i2,ino)
      r1_mar(i2,24) = rk_mar(i2,24)*s(i2,ich3so2oo)
      r2_mar(i2,24) = rk_mar(i2,24)*s(i2,ino)

      r_mar(i2,25) = rk_mar(i2,25)*s(i2,ich3so2oo)*s(i2,ich3o2)
      r1_mar(i2,25) = rk_mar(i2,25)*s(i2,ich3so2oo)
      r2_mar(i2,25) = rk_mar(i2,25)*s(i2,ich3o2)

      r_mar(i2,26) = rk_mar(i2,26)*s(i2,ich3so3)
      r1_mar(i2,26) = rk_mar(i2,26)  

      r_mar(i2,27) = rk_mar(i2,27)*s(i2,ich3so3)*s(i2,ino2)
      r1_mar(i2,27) = rk_mar(i2,27)*s(i2,ich3so3)
      r2_mar(i2,27) = rk_mar(i2,27)*s(i2,ino2)

      r_mar(i2,28) = rk_mar(i2,28)*s(i2,ich3so3)*s(i2,ino)
      r1_mar(i2,28) = rk_mar(i2,28)*s(i2,ich3so3)
      r2_mar(i2,28) = rk_mar(i2,28)*s(i2,ino)

      r_mar(i2,29) = rk_mar(i2,29)*s(i2,ich3so3)*s(i2,iho2)
      r1_mar(i2,29) = rk_mar(i2,29)*s(i2,ich3so3)
      r2_mar(i2,29) = rk_mar(i2,29)*s(i2,iho2)

      r_mar(i2,30) = rk_mar(i2,30)*s(i2,ich3so3)*s(i2,ihcho)
      r1_mar(i2,30) = rk_mar(i2,30)*s(i2,ich3so3)
      r2_mar(i2,30) = rk_mar(i2,30)*s(i2,ihcho)

      r_mar(i2,31) = rk_mar(i2,31)*s(i2,ich3sch2oo)*s(i2,ich3so2)
      r1_mar(i2,31) = rk_mar(i2,31)*s(i2,ich3sch2oo)
      r2_mar(i2,31) = rk_mar(i2,31)*s(i2,ich3so2)

      r_mar(i2,32) = rk_mar(i2,32)*s(i2,ich3sch2oo)*s(i2,ich3sch2oo)
      r1_mar(i2,32) = rk_mar(i2,32)*s(i2,ich3sch2oo)
      

      r_mar(i2,33) = rk_mar(i2,33)*s(i2,iso2)
      r1_mar(i2,33) = rk_mar(i2,33) 

      r_mar(i2,34) = rk_mar(i2,34)*s(i2,idmso)
      r1_mar(i2,34) = rk_mar(i2,34) 

      r_mar(i2,35) = rk_mar(i2,35)*s(i2,idmso2)
      r1_mar(i2,35) = rk_mar(i2,35)
      !endif !if bmask
      enddo !i2

      return
      end
!!!-----------------------------------------------------------------
!!!------------------------------------------------------------------
      subroutine gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
!dir$ assume_aligned s:64,r_urb:64,r1_urb:64,r2_urb
      implicit none
      type(cbmztype), target :: cbmzobj

      real :: s(VLEN,76),r_urb(VLEN,76),r1_urb(VLEN,76),r2_urb(VLEN,76)
!dir$ attributes align:64 :: uindex1,uindex2
      real, parameter, dimension(49) :: uindex1 = &
        (/35,36,36,37,37,37,38,38,39,40,39,40,39,40,41,42,44,43,43,45, &
          46,46,46,48,48,47,47,49,50,51,52,49,50,51,52,49,50,51,52,49, &
          50,51,52,35, 7, 7, 7, 7, 7/)
      real, parameter, dimension(49) :: uindex2 = &
        (/14,-1,14,-1,14, 7,11,14,11,11,14,14, 7, 7,14,14, 5,14, 7, 6, &
          14,-1,11,-1,14,14,-1, 5, 5, 5, 5, 7, 7, 7, 7,15,15,15,15,-1, &
          -1,-1,-1,53,35,38,41,36,42/)
      real :: res_urb
      integer :: i, i2
#ifdef KNL_OPT
      real, pointer :: rk_urb(:,:)
      rk_urb => cbmzobj%rk_urb
!dir$ assume_aligned r_urb:64,r1_urb:64,r2_urb:64,s:64
!!dir$ simd
#endif

      do i=1,49
       if (uindex2(i) == -1) then
         !dir$ simd private(res_urb)
         do i2=1,VLEN
         !if (cbmzobj%bmask(i2)) then
            res_urb = rk_urb(i2,i)*s(i2,uindex1(i))
            r_urb(i2,i) = res_urb
            r1_urb(i2,i) = rk_urb(i2,i)
         !endif !if bmask
         enddo !i2
       else
         !dir$ simd private(res_urb)
         do i2=1,VLEN
          res_urb = rk_urb(i2,i)*s(i2,uindex1(i))
          r_urb(i2,i) = res_urb*s(i2,uindex2(i))
          r1_urb(i2,i) = res_urb
          r2_urb(i2,i) = rk_urb(i2,i)*s(i2,uindex2(i))
         enddo !i2
       endif
      enddo

      return
      end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!*********************************************************************
!!*********************************************************************
#ifdef KNL_OPT
      !!!dir$ attributes forceinline :: ode_com
      subroutine ode_com(cbmzobj,r_com,r1_com,r2_com,p_com,rl_com)
      !use gas_data, only : 76,76
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
!dir$ assume_aligned cppb:64,r_com:64,r1_com:64,r2_com:64,p_com
      implicit none
      type(cbmztype), target :: cbmzobj
      real :: r_com(VLEN,76),r1_com(VLEN,76),r2_com(VLEN,76)
      real :: p_com(VLEN,76),rl_com(VLEN,76)
      integer :: i2

#else
      subroutine ode_com
      use gas_data
#endif

    !dir$ simd
    do i2 = 1, VLEN
    !if (cbmzobj%pmask(i2)) then
      p_com(i2,ih2so4)= r_com(i2,45)                                                  
      p_com(i2,ihno3)= r_com(i2,24)+.3*r_com(i2,41)+r_com(i2,42)+r_com(i2,42)+r_com(i2,52)+r_com(i2,68)
      !p_com(i2,ihcl)= 0.0                                                              
      !p_com(i2,inh3)= 0.0                                                              
      p_com(i2,ino)= r_com(i2,1)+0.11*r_com(i2,2)+r_com(i2,3)+r_com(i2,15)+r_com(i2,38)  
      p_com(i2,ino2)= 0.89*r_com(i2,2)+r_com(i2,4)+r_com(i2,5) &
             +r_com(i2,6)+r_com(i2,17)+r_com(i2,18)     &                        
             +r_com(i2,25)+r_com(i2,26)+r_com(i2,28)    &
             +r_com(i2,33)+r_com(i2,36)+r_com(i2,37)    &                          
             +r_com(i2,37)+r_com(i2,38)+r_com(i2,40)+r_com(i2,40)+.7*r_com(i2,41)   &                            
             +r_com(i2,43)+r_com(i2,57)+r_com(i2,58)+r_com(i2,59)                &
             +r_com(i2,60)+r_com(i2,70)+r_com(i2,71)+r_com(i2,72)
      p_com(i2,ino3)= r_com(i2,6)+r_com(i2,16)+r_com(i2,19)+r_com(i2,27)+r_com(i2,43) 
      p_com(i2,in2o5)= r_com(i2,39)                                                           
      p_com(i2,ihono)= r_com(i2,23)+r_com(i2,35)                                                     
      p_com(i2,ihno4)= r_com(i2,34)                                                           
      p_com(i2,io3)= r_com(i2,13)+.4*r_com(i2,73) 
      p_com(i2,io1d)= r_com(i2,8)                                                             
      p_com(i2,io3p)= r_com(i2,1)+0.89*r_com(i2,2)+r_com(i2,7)        &
                  +r_com(i2,10)+r_com(i2,11)
      p_com(i2,ioh)= r_com(i2,3)+r_com(i2,4)+2*r_com(i2,9)   &
            +2*r_com(i2,12)+r_com(i2,21)+r_com(i2,33)     &                         
            +.7*r_com(i2,41)+r_com(i2,53)+r_com(i2,54)    &
            +.3*r_com(i2,55)+.5*r_com(i2,56) 
      p_com(i2,iho2)= r_com(i2,5)+r_com(i2,20)+r_com(i2,22)        &
             +r_com(i2,25)+r_com(i2,30)+r_com(i2,36)            &                   
             +r_com(i2,44)+r_com(i2,45)+r_com(i2,48)            &
             +2*r_com(i2,49)+r_com(i2,51)                    &              
             +r_com(i2,52)+r_com(i2,53)+r_com(i2,54)            &
             +r_com(i2,57)+r_com(i2,58)+r_com(i2,59)            &                  
             +r_com(i2,60)+.32*r_com(i2,63)+.6*r_com(i2,64)     &
             +r_com(i2,65)+r_com(i2,66) 
      p_com(i2,ih2o2)= r_com(i2,31)+r_com(i2,32)                                                     
      p_com(i2,ico)= r_com(i2,49)+r_com(i2,50)+r_com(i2,51)     &
                 +r_com(i2,52)+r_com(i2,66)
      !p_com(i2,iso2)= 0.0                                                              
      !p_com(i2,ich4)= 0.0                                                              
      p_com(i2,ic2h6)= .2*r_com(i2,64)                                                        
      p_com(i2,ich3o2)= r_com(i2,46)+.7*r_com(i2,55)+r_com(i2,66)   &
               +r_com(i2,71)+r_com(i2,72)                     &              
               +r_com(i2,74) 
      p_com(i2,iethp)= r_com(i2,47)+.5*r_com(i2,56) 
      p_com(i2,ihcho)= r_com(i2,48)+r_com(i2,53)+.3*r_com(i2,55)   &
              +r_com(i2,57)+r_com(i2,59)                     &           
              +.66*r_com(i2,63)
      p_com(i2,ich3oh)= .34*r_com(i2,63)                                                      
      !p_com(i2,ianol)= 0.0                                                             
      p_com(i2,ich3ooh)= r_com(i2,61)                                                         
      p_com(i2,iethooh)= r_com(i2,62)                                                         
      p_com(i2,iald2)= r_com(i2,54)+.5*r_com(i2,56)+r_com(i2,58) &
              +r_com(i2,60)+.8*r_com(i2,64)                &             
              +r_com(i2,65)   
      !p_com(i2,ihcooh)= 0.0                                                            
      p_com(i2,ircooh)= .4*r_com(i2,73)                                                       
      p_com(i2,ic2o3)= r_com(i2,67)+r_com(i2,68)+r_com(i2,70) 
      p_com(i2,ipan)= r_com(i2,69)                                                            
      !p_com(i2,idso4)= 0.0                                         
      !p_com(i2,idno3)= 0.0                                                            
!
      !rl_com(i2,ih2so4)= 0.0                                                            
      rl_com(i2,ihno3)= r1_com(i2,4)+r1_com(i2,27)   
      !rl_com(i2,ihcl)= 0.0                                                              
      !rl_com(i2,inh3)= 0.0                                                              
      rl_com(i2,ino)= r1_com(i2,17)+r1_com(i2,18)+r1_com(i2,23)+r1_com(i2,33)+r1_com(i2,37)  &
                  +r1_com(i2,57)+r1_com(i2,58)+r1_com(i2,71)                                               
      rl_com(i2,ino2)= r1_com(i2,1)+r1_com(i2,15)+r1_com(i2,16)+r1_com(i2,19)+r1_com(i2,24)  &
                 +r1_com(i2,34)+r1_com(i2,35)+r1_com(i2,38)+r1_com(i2,39)+r1_com(i2,69) 
      rl_com(i2,ino3)= r1_com(i2,2)+r1_com(i2,25)+r2_com(i2,37)+r2_com(i2,38)    &
             +r2_com(i2,39)+r1_com(i2,40)                            &       
             +r1_com(i2,40)+r2_com(i2,41)+r1_com(i2,52)+r1_com(i2,59)        &
             +r1_com(i2,60)+r1_com(i2,68)                            &                      
             +r1_com(i2,72)    &                                                        
             +r2_com(i2,75)+r2_com(i2,76) !lijie add no3+vocs on 20130327  
      rl_com(i2,in2o5)= r1_com(i2,6)+r1_com(i2,42)+r1_com(i2,43)                                              
      rl_com(i2,ihono)= r1_com(i2,3)+r1_com(i2,26)                                                      
      rl_com(i2,ihno4)= r1_com(i2,5)+r1_com(i2,28)+r1_com(i2,36)                                              
      rl_com(i2,io3)= r1_com(i2,7)+r1_com(i2,8)+r1_com(i2,14)+r2_com(i2,18)  &
                 +r2_com(i2,19)+r2_com(i2,20)+r2_com(i2,21)                                                
      rl_com(i2,io1d)= r1_com(i2,10)+r1_com(i2,11)+r1_com(i2,12)                                              
      rl_com(i2,io3p)= r1_com(i2,13)+r2_com(i2,14)+r2_com(i2,15)+r2_com(i2,16)+r2_com(i2,17)
      rl_com(i2,ioh)= r1_com(i2,20)+r1_com(i2,22)+r2_com(i2,23)        &
            +r2_com(i2,24)+r2_com(i2,25)+r2_com(i2,26)             &                  
            +r2_com(i2,27)+r2_com(i2,28)+r2_com(i2,29)+r2_com(i2,30)   &
            +r1_com(i2,44)+r1_com(i2,45)                       &        
            +r1_com(i2,46)+r1_com(i2,47)+r1_com(i2,48)             &
            +r1_com(i2,51)+r1_com(i2,55)+r1_com(i2,56)             &                  
            +r1_com(i2,65)+r1_com(i2,67)                                                       
      rl_com(i2,iho2)= r1_com(i2,21)+r1_com(i2,29)+r1_com(i2,31)    &
             +r1_com(i2,31)+r1_com(i2,32)+r1_com(i2,32)         &                     
             +r2_com(i2,33)+r2_com(i2,34)+r2_com(i2,35)         &
             +r1_com(i2,41)+r1_com(i2,61)+r1_com(i2,62)         &                     
             +r1_com(i2,73)                                                            
      rl_com(i2,ih2o2)= r1_com(i2,9)+r1_com(i2,30)                                                      
      rl_com(i2,ico)= r2_com(i2,44)                                                             
      rl_com(i2,iso2)= r2_com(i2,45)                                                            
      rl_com(i2,ich4)= r2_com(i2,46)                                                            
      rl_com(i2,ic2h6)= r2_com(i2,47)                                                           
      rl_com(i2,ich3o2)= r2_com(i2,57)+r2_com(i2,59)+r2_com(i2,61)+r1_com(i2,63) 
      rl_com(i2,iethp)= r2_com(i2,58)+r2_com(i2,60)+r2_com(i2,62)+r1_com(i2,64) 
      rl_com(i2,ihcho)= r1_com(i2,49)+r1_com(i2,50)+r2_com(i2,51)+r2_com(i2,52)       
      rl_com(i2,ich3oh)= r2_com(i2,48)                                                          
      rl_com(i2,ianol)= r2_com(i2,65)                                                           
      rl_com(i2,ich3ooh)= r1_com(i2,53)+r2_com(i2,55)                                                   
      rl_com(i2,iethooh)= r1_com(i2,54)+r2_com(i2,56)                                                   
      rl_com(i2,iald2)= r1_com(i2,66)+r2_com(i2,67)+r2_com(i2,68)                                           
      !rl_com(i2,ihcooh)= 0.0                                                            
      !rl_com(i2,ircooh)= 0.0                                                            
      rl_com(i2,ic2o3)= r2_com(i2,69)+r2_com(i2,71)+r2_com(i2,72) &
             +r2_com(i2,73)+r1_com(i2,74)    
      rl_com(i2,ipan)= r1_com(i2,70)  
      !rl_com(i2,idso4)= 0.0   
      !rl_com(i2,idno3)= 0.0  
    !endif ! if pmask
    enddo !i2
      return
      end


!!*********************************************************************
!!*********************************************************************
#ifdef KNL_OPT
      !subroutine ode_bio(r_bio,r1_bio,r2_bio,p_bio,rl_bio)
      subroutine ode_bio(cbmzobj)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      !real :: r_bio(80),r1_bio(80),r2_bio(80)
      !real :: p_bio(80),rl_bio(80)
      integer :: i2
      real, dimension(:,:), pointer :: p_bio,rl_bio
   
      integer :: ih2so4,ihno3,ihcl,inh3,ino,ino2,ino3,in2o5,ihono,ihno4
      integer :: io3,io1d,io3p,ioh,iho2,ih2o2,ico,iso2,ich4,ic2h6,ich3o2
      integer :: iethp,ihcho,ich3oh,ianol,ich3ooh,iethooh,iald2,ihcooh,ircooh
      integer :: ic2o3,ipan,ipar,iaone,imgly,ieth,iolet,iolei,itol,ixyl
      integer :: icres,ito2,icro,iopen,ionit,irooh,iro2,iano2,inap,ixo2
      integer :: ixpar,isv1,isv2,iisop,iisoprd,iisopp,iisopn,iisopo2,iterp,isv3,isv4,isv5,isv6

      p_bio => cbmzobj%p_bio
      rl_bio => cbmzobj%rl_bio
      ih2so4	= 1
      ihno3	= 2
      ihcl	= 3
      inh3	= 4
      ino	= 5
      ino2	= 6
      ino3	= 7
      in2o5	= 8
      ihono	= 9
      ihno4	= 10
      io3	= 11
      io1d	= 12
      io3p	= 13
      ioh	= 14
      iho2	= 15
      ih2o2	= 16
      ico	= 17
      iso2	= 18
      ich4	= 19
      ic2h6	= 20
      ich3o2	= 21
      iethp	= 22
      ihcho	= 23
      ich3oh	= 24
      ianol	= 25
      ich3ooh	= 26
      iethooh	= 27
      iald2	= 28
      ihcooh	= 29
      ircooh	= 30
      ic2o3	= 31
      ipan	= 32

      ipar	= 34 + 1
      iaone	= 34 + 2
      imgly	= 34 + 3
      ieth	= 34 + 4
      iolet	= 34 + 5
      iolei	= 34 + 6
      itol	= 34 + 7
      ixyl	= 34 + 8
      icres	= 34 + 9
      ito2	= 34 + 10
      icro	= 34 + 11
      iopen	= 34 + 12
      ionit	= 34 + 13
      irooh	= 34 + 14
      iro2	= 34 + 15
      iano2	= 34 + 16
      inap	= 34 + 17
      ixo2	= 34 + 18
      ixpar	= 34 + 19
      isv1      = 34 + 20
      isv2      = 34 + 21

      iisop	= 55 + 1 !56
      iisoprd	= 55 + 2 !57
      iisopp	= 55 + 3 !58
      iisopn	= 55 + 4 !59
      iisopo2	= 55 + 5 !60
      iterp     = 55 + 6 !61
      isv3      = 55 + 7 !62
      isv4      = 55 + 8 !63
      isv5      = 55 + 9 !64
      isv6      = 55 + 10 !65
 
#else
      subroutine ode_bio
      use gas_data
#endif
    do i2=1, VLEN
      p_bio(i2,ino)= 0.0
      rl_bio(i2,ino)= r1_bio(i2,8)+r1_bio(i2,9)+r1_bio(i2,10)

      p_bio(i2,ino2)= .91*r_bio(i2,8)+1.2*r_bio(i2,9)+r_bio(i2,10)+0.47*r_bio(i2,20) !(i2,20)added by li jie
      rl_bio(i2,ino2)= 0.0

      p_bio(i2,ino3)= 0.0
      rl_bio(i2,ino3)= r1_bio(i2,3)+r1_bio(i2,7)+r1_bio(i2,20) !(i2,20)added by li jie

      p_bio(i2,ihno3)= .07*r_bio(i2,7)
      rl_bio(i2,ihno3)= 0.0

      p_bio(i2,io3)= 0.0
      rl_bio(i2,io3)= r1_bio(i2,2)+r1_bio(i2,6)+r1_bio(i2,19) !(i2,19)added by li jie

      p_bio(i2,ioh)= .27*r_bio(i2,2)+.27*r_bio(i2,6)+0.57*r_bio(i2,19) !(i2,19)added by li jie
      rl_bio(i2,ioh)= r1_bio(i2,1)+r1_bio(i2,5)+r1_bio(i2,18)     !(i2,18)added by li jie

      p_bio(i2,iho2)= .07*r_bio(i2,2)+.33*r_bio(i2,4)+.1*r_bio(i2,6)+.93*r_bio(i2,7)+.91*r_bio(i2,8)+.8*r_bio(i2,9)+r_bio(i2,10) &
                   +.75*r_bio(i2,18)+.07*r_bio(i2,19)+.28*r_bio(i2,20) !(i2,18)(i2,19)(i2,20)added by li jie

      rl_bio(i2,iho2)= r1_bio(i2,11)+r1_bio(i2,12)+r1_bio(i2,13)

      p_bio(i2,ih2o2)= 0.0
      rl_bio(i2,ih2o2)= 0.0

      p_bio(i2,ico)= .07*r_bio(i2,2)+.33*r_bio(i2,4)+.16*r_bio(i2,6)+.64*r_bio(i2,7)+.59*r_bio(i2,10) &
                  +0.001*r_bio(i2,19) !(i2,19)added by li jie

      rl_bio(i2,ico)= 0.0

      p_bio(i2,ihcho)= .6*r_bio(i2,2)+.2*r_bio(i2,4)+.15*r_bio(i2,6)+.28*r_bio(i2,7)+.63*r_bio(i2,8)+.25*r_bio(i2,10) &
                    +.28*r_bio(i2,18)+.24*r_bio(i2,19) !(i2,18)(i2,19)added by li jie

      rl_bio(i2,ihcho)= 0.0


      p_bio(i2,iald2)= .15*r_bio(i2,2)+.07*r_bio(i2,4)+.02*r_bio(i2,6)+.28*r_bio(i2,7)+.8*r_bio(i2,9)+.55*r_bio(i2,10)+r_bio(i2,15)+.5*r_bio(i2,16) &
                    +.15*r_bio(i2,17)+.47*r_bio(i2,18)+.21*r_bio(i2,19)+0.47*r_bio(i2,20)  !(i2,17)(i2,18)(i2,19)(i2,20)added by li jie

      rl_bio(i2,iald2)= 0.0


      p_bio(i2,ipar)= 1.86*r_bio(i2,7)+.18*r_bio(i2,8)+1.6*r_bio(i2,9)+2*r_bio(i2,12)+2*r_bio(i2,15) &
                   +0.51*r_bio(i2,17)+7.*r_bio(i2,19) !(i2,17)(i2,19)added by li jie
      rl_bio(i2,ipar)= 0.0

      p_bio(i2,iaone)= .03*r_bio(i2,4)+.09*r_bio(i2,6)+.63*r_bio(i2,10)+.5*r_bio(i2,16)
      rl_bio(i2,iaone)= 0.0

      p_bio(i2,imgly)= .85*r_bio(i2,6)+.34*r_bio(i2,10)
      rl_bio(i2,imgly)= 0.0

      p_bio(i2,ionit)= .93*r_bio(i2,7)+.09*r_bio(i2,8)+.8*r_bio(i2,9)+r_bio(i2,12)+r_bio(i2,15) &
                    +0.53*r_bio(i2,20) !(i2,20)added by li jie
      rl_bio(i2,ionit)= 0.0

      p_bio(i2,ircooh)= .39*r_bio(i2,2)+.46*r_bio(i2,6)
      rl_bio(i2,ircooh)= 0.0

      p_bio(i2,irooh)= r_bio(i2,11)+r_bio(i2,13)
      rl_bio(i2,irooh)= 0.0

      p_bio(i2,ich3o2)= .7*r_bio(i2,4)+.05*r_bio(i2,6)
      rl_bio(i2,ich3o2)= 0.0

      p_bio(i2,ic2o3)= .2*r_bio(i2,2)+.97*r_bio(i2,4)+.5*r_bio(i2,5)+.11*r_bio(i2,6)+.07*r_bio(i2,7)
      rl_bio(i2,ic2o3)= 0.0

      p_bio(i2,ixo2)= .08*r_bio(i2,1)+.2*r_bio(i2,2)+.2*r_bio(i2,5)+.07*r_bio(i2,6)+.93*r_bio(i2,7) &
                   +1.25*r_bio(i2,18)+.76*r_bio(i2,19)+.76*r_bio(i2,20) !(i2,18)(i2,19)(i2,20)added by li jie
      rl_bio(i2,ixo2)= 0.0

      p_bio(i2,iisop)= 0.0
      rl_bio(i2,iisop)= r2_bio(i2,1)+r2_bio(i2,2)+r2_bio(i2,3)

      p_bio(i2,iisoprd)= .65*r_bio(i2,2)+.91*r_bio(i2,8)+.2*r_bio(i2,9)+r_bio(i2,14)
      rl_bio(i2,iisoprd)= r1_bio(i2,4)+r2_bio(i2,5)+r2_bio(i2,6)+r2_bio(i2,7)

      p_bio(i2,iisopp)= r_bio(i2,1)
      rl_bio(i2,iisopp)= r2_bio(i2,8)+r2_bio(i2,11)+r1_bio(i2,14)

      p_bio(i2,iisopn)= r_bio(i2,3)
      rl_bio(i2,iisopn)= r2_bio(i2,9)+r2_bio(i2,12)+r1_bio(i2,15)

      p_bio(i2,iisopo2)= .5*r_bio(i2,5)
      rl_bio(i2,iisopo2)= r2_bio(i2,10)+r2_bio(i2,13)+r1_bio(i2,16)

!!! the following p_bio(i2,io1d),rl_bio(i2,io1d),sv3,sv4,sv5,sv6,terp are added by li jie 

      p_bio(i2,io1d)= 0.0
      rl_bio(i2,io1d)= r1_bio(i2,17)

      p_bio(i2,isv3) = .110681*r_bio(i2,17)+.063304*r_bio(i2,18)+.10998*r_bio(i2,19)+.092287*r_bio(i2,20)
      rl_bio(i2,isv3) = 0.0

      p_bio(i2,isv4) = .4172*r_bio(i2,17)+.324763*r_bio(i2,18)+.389277*r_bio(i2,19)+.395069*r_bio(i2,20)
      rl_bio(i2,isv4) = 0.0

      p_bio(i2,isv5) = .232*r_bio(i2,1)
      rl_bio(i2,isv5) = 0.0

      p_bio(i2,isv6) = .0288*r_bio(i2,1)
      rl_bio(i2,isv6) = 0.0

      p_bio(i2,iterp) = 0.0
      rl_bio(i2,iterp) = r2_bio(i2,17)+r2_bio(i2,18)+r2_bio(i2,19)+r2_bio(i2,20)
    enddo !i2

      return
      end









!!*********************************************************************
!!*********************************************************************
      subroutine ode_het(r_het,r1_het,p_het,rl_het)
      use gas_data, only:VLEN
!dir$ assume_aligned r_het:64,r1_het:64,p_het:64,rl_het:64
      !implicit none
      real :: r_het(VLEN,76), r1_het(VLEN,28), p_het(VLEN,76), rl_het(VLEN,76)
      integer :: i2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! to comment all lines in this subroutine and add a new form (by lijie)



!       p_het(in2o5) = 0.0
!       p_het(ino3)  = 0.0
!       p_het(iho2)  = 0.0
!       p_het(ihcho) = 0.0
!       p_het(ioh)   = 0.0
!       p_het(io3)   = 0.0
!       rl_het(ihono) = 0.0
!       p_het(iso2) = 0.0
!       p_het(irooh) = 0.0
!       p_het(ich3oh) = 0.0
!       rl_het(idso4) = 0.0
!       rl_het(idno3) = 0.0

    !dir$ simd
    do i2 = 1, VLEN
    !if (cbmzobj%pmask(i2)) then
       rl_het(i2,8) = r1_het(i2,1) + r1_het(i2,10) + r1_het(i2,15) + r1_het(i2,23)

       p_het(i2,6)  = r_het(i2,9)
       rl_het(i2,6)  = r1_het(i2,8) + r1_het(i2,2) +r1_het(i2,13)

       rl_het(i2,7)  = r1_het(i2,3) + r1_het(i2,14) + r1_het(i2,24) + r1_het(i2,27)

       rl_het(i2,15)  = r1_het(i2,4) + r1_het(i2,17) + r1_het(i2,25)

       rl_het(i2,23) = r1_het(i2,5) + r1_het(i2,22)

       rl_het(i2,14)   = r1_het(i2,6) + r1_het(i2,16)

       rl_het(i2,11)   = r1_het(i2,7) + r1_het(i2,11)


       p_het(i2,2) = 2.*r_het(i2,1) + 0.5*r_het(i2,2) + r_het(i2,3) + 2.* r_het(i2,10)  &
                    + 0.5*r_het(i2,13) + r_het(i2,14) + 2.*r_het(i2,15)       &
                     + 2.*r_het(i2,23) + r_het(i2,24)

       rl_het(i2,2) = r1_het(i2,9) + r1_het(i2,12) +r1_het(i2,28)


       p_het(i2,9) = 0.5*r_het(i2,2) + r_het(i2,8) + 0.5*r_het(i2,13) + 0.5 * r_het(i2,25)
                      

       p_het(i2,16) = 0.5 * r_het(i2,4) + 0.5 * r_het(i2,17)
       rl_het(i2,16) = r1_het(i2,18)

       rl_het(i2,18) = r1_het(i2,19) + r1_het(i2,26)

       rl_het(i2,48) = r1_het(i2,20)

       rl_het(i2,24) = r1_het(i2,21)

       p_het(i2,33) = r_het(i2,19) + r_het(i2,26)

       p_het(i2,34) = r_het(i2,12) + r_het(i2,27)+r_het(i2,28)

    !endif ! if pmask
    enddo !i2
            
      return
      end





!!*********************************************************************
!!*********************************************************************
#ifdef KNL_OPT
      subroutine ode_mar(cbmzobj)
      use gas_data
!dir$ assume_aligned cnn:64,emission:64,emit:64,rk_bio:64,rk_mar:64,&
        rk_param:64,r_bio:64,r_mar:64,r1_bio:64,r1_mar:64,&
        r2_bio:64,r2_mar:64,p_bio:64,rl_bio:64,&
        p_mar:64,rl_mar:64,total_p:64,total_l:64,&
        tcur_sec:64,tmid_sec:64,told_sec:64,dt_sec:64,&
        rlon:64,rlat:64,zalt_m:64,cos_sza:64
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(VLEN) :: a, b
      integer :: i2
      integer :: ino,ino2,ino3,ihono,ihno3,io3,io3p,ioh,iho2,ih2o2,iso2
      integer :: ih2so4,ich3o2,ich3ooh,ihcho
      real, dimension(:,:), pointer :: p_mar,rl_mar
      p_mar => cbmzobj%p_mar
      rl_mar => cbmzobj%rl_mar
   
      a(:) = 5.e+5/(5.e+5 + cbmzobj%o2(:)*3.e-12)
      b(:) = 1.5e+7/(1.5e+7 + cbmzobj%o2(:)*1.2e-12)
      ih2so4	= 1
      ino	= 5
      ino2	= 6
      ino3	= 7
      ihono	= 9
      ihno3	= 2
      io3	= 11
      io3p	= 13
      ioh	= 14
      iho2	= 15
      ih2o2	= 16
      iso2	= 18
      ich3o2	= 21
      ich3ooh	= 26
      ihcho	= 23
#else
      subroutine ode_mar
      use gas_data
      a = 5.e+5/(5.e+5 + o2*3.e-12)
      b = 1.5e+7/(1.5e+7 + o2*1.2e-12)
#endif

    !dir$ simd
    do i2 = 1, VLEN
    !if (cbmzobj%bmask(i2)) then
      p_mar(i2,ino)= r_mar(i2,17) 
  
      rl_mar(i2,ino)= r1_mar(i2,5)+r1_mar(i2,9)+r1_mar(i2,24)+r1_mar(i2,28)



      p_mar(i2,ino2)= r_mar(i2,5)+r_mar(i2,9)+r_mar(i2,24)   

      rl_mar(i2,ino2)= r1_mar(i2,17)+r1_mar(i2,27)



      p_mar(i2,ino3)=   0.0
      rl_mar(i2,ino3)= r1_mar(i2,2)+r1_mar(i2,12)   



      p_mar(i2,ihono)= r_mar(i2,28)   
      rl_mar(i2,ihono)=   0.0



      p_mar(i2,ihno3)= r_mar(i2,2)+r_mar(i2,12)+r_mar(i2,27)   
      rl_mar(i2,ihno3)=   0.0



      p_mar(i2,io3)=   0.0
      rl_mar(i2,io3)= r1_mar(i2,18)   



      p_mar(i2,io3p)=   0.0
      rl_mar(i2,io3p)= r1_mar(i2,3)    



      p_mar(i2,ioh)= r_mar(i2,19)  
 
      rl_mar(i2,ioh)= r1_mar(i2,1)+r1_mar(i2,4)+r1_mar(i2,7)+r1_mar(i2,8)+  &
                 r1_mar(i2,14)+r1_mar(i2,21)   



      p_mar(i2,iho2)= (1.-a(i2))*r_mar(i2,4)+r_mar(i2,6)+(1.-b(i2))*r_mar(i2,7)+  &
                  r_mar(i2,10)+r_mar(i2,20)+r_mar(i2,25)+r_mar(i2,30)   

      rl_mar(i2,iho2)= r1_mar(i2,11)+r1_mar(i2,19)+r1_mar(i2,29)   



      p_mar(i2,ih2o2)= r_mar(i2,11)   
      rl_mar(i2,ih2o2)=   0.0



      p_mar(i2,iso2)= r_mar(i2,16)   
      rl_mar(i2,iso2)= r1_mar(i2,33)



      p_mar(i2,ih2so4)= r_mar(i2,26)   
      rl_mar(i2,ih2so4)=   0.0



      p_mar(i2,ich3o2)= r_mar(i2,3)+a(i2)*r_mar(i2,4)+b(i2)*r_mar(i2,7)+r_mar(i2,16)+   &
                   r_mar(i2,26)   

      rl_mar(i2,ich3o2)= r1_mar(i2,6)+r1_mar(i2,10)+r1_mar(i2,13)+r1_mar(i2,20)+    &
                    r1_mar(i2,25)   



      p_mar(i2,ich3ooh)= r_mar(i2,13)   
      rl_mar(i2,ich3ooh)=   0.0



      p_mar(i2,ihcho)= r_mar(i2,5)+r_mar(i2,6)+r_mar(i2,6)+r_mar(i2,9)+     &
                   r_mar(i2,10)+r_mar(i2,10)+r_mar(i2,20)+r_mar(i2,25) 
  
      rl_mar(i2,ihcho)= r1_mar(i2,30)   



      p_mar(i2,idms)= 0.0
      rl_mar(i2,idms)= r2_mar(i2,1)+r2_mar(i2,2)+r2_mar(i2,3)+r2_mar(i2,4)    



      p_mar(i2,imsa)= r_mar(i2,15)+r_mar(i2,21)+r_mar(i2,27)+r_mar(i2,28)+     &
                  r_mar(i2,29)+r_mar(i2,30)   
      rl_mar(i2,imsa)=   0.0



      p_mar(i2,idmso)= (1.-a(i2))*r_mar(i2,4)   
 
      rl_mar(i2,idmso)= r2_mar(i2,7) + r1_mar(i2,34)



      p_mar(i2,idmso2)= (1.-b(i2))*r_mar(i2,7)   
 
      rl_mar(i2,idmso2)= r2_mar(i2,8) + r1_mar(i2,35)



      p_mar(i2,ich3so2h)= b(i2)*r_mar(i2,7)   
 
      rl_mar(i2,ich3so2h)= r2_mar(i2,11)+r2_mar(i2,12)+r2_mar(i2,13)+r2_mar(i2,14)+    &
                      r2_mar(i2,15)   



      p_mar(i2,ich3sch2oo)= r_mar(i2,1)+r_mar(i2,2)

      rl_mar(i2,ich3sch2oo)= r2_mar(i2,5)+r2_mar(i2,6)+r2_mar(i2,31)+2.*r1_mar(i2,32)



      p_mar(i2,ich3so2)= r_mar(i2,3)+a(i2)*r_mar(i2,4)+r_mar(i2,5)+r_mar(i2,6)+     &
                      r_mar(i2,9)+r_mar(i2,10)+r_mar(i2,11)+r_mar(i2,12)+    &
                      r_mar(i2,13)+r_mar(i2,14)+r_mar(i2,15)+r_mar(i2,23)+   &
                      r_mar(i2,31)+1.85*r_mar(i2,32)

      rl_mar(i2,ich3so2)= r1_mar(i2,16)+r2_mar(i2,17)+r2_mar(i2,18)+r2_mar(i2,19)+  &
                     r2_mar(i2,20)+r2_mar(i2,21)+r1_mar(i2,22)+r1_mar(i2,31)



      p_mar(i2,ich3so3)= r_mar(i2,17)+r_mar(i2,18)+r_mar(i2,19)+r_mar(i2,20)+  &
                      r_mar(i2,24)+r_mar(i2,25)+r_mar(i2,31)

      rl_mar(i2,ich3so3)= r1_mar(i2,15)+r1_mar(i2,26)+r2_mar(i2,27)+r2_mar(i2,28)+  &
                      r2_mar(i2,29)+r2_mar(i2,30)  

 

      p_mar(i2,ich3so2oo)= r_mar(i2,22) 
  
      rl_mar(i2,ich3so2oo)= r1_mar(i2,23)+r2_mar(i2,24)+r2_mar(i2,25) 

  

      p_mar(i2,ich3so2ch2oo)= r_mar(i2,8) 
   
      rl_mar(i2,ich3so2ch2oo)= r2_mar(i2,9)+r2_mar(i2,10)



      p_mar(i2,isulfhox)= 0.15*r_mar(i2,32)
      rl_mar(i2,isulfhox)= 0.0

    !endif ! if pmask
    enddo !i2

      return
      end




!!*********************************************************************
!!*********************************************************************
#ifdef KNL_OPT
      subroutine ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
      use gas_data, only:VLEN
!dir$ assume_aligned r_urb:64,r1_urb:64,r2_urb:64,p_urb:64
      !implicit none
      real :: r_urb(VLEN,76),r1_urb(VLEN,76),r2_urb(VLEN,76)
      real :: p_urb(VLEN,76),rl_urb(VLEN,76)
      integer :: i2

      ipar	= 34 + 1
      iaone	= 34 + 2
      imgly	= 34 + 3
      ieth	= 34 + 4
      iolet	= 34 + 5
      iolei	= 34 + 6
      itol	= 34 + 7
      ixyl	= 34 + 8
      icres	= 34 + 9
      ito2	= 34 + 10
      icro	= 34 + 11
      iopen	= 34 + 12
      ionit	= 34 + 13
      irooh	= 34 + 14
      iro2	= 34 + 15
      iano2	= 34 + 16
      inap	= 34 + 17
      ixo2	= 34 + 18
      ixpar	= 34 + 19
      isv1      = 34 + 20
      isv2      = 34 + 21

#else
      subroutine ode_urb
      use gas_data
#endif

    !dir$ simd
    do i2 = 1, VLEN
    !if (cbmzobj%bmask(i2)) then
      p_urb(i2,2)= r_urb(i2,6)+r_urb(i2,19)                                                      
      rl_urb(i2,2)= 0.0                                                             

      p_urb(i2,5)= 0.0                                                               
      rl_urb(i2,5)= r1_urb(i2,17)+r1_urb(i2,28)+r1_urb(i2,29)+r1_urb(i2,30)+r1_urb(i2,31)    

                                 

      p_urb(i2,6)= .95*r_urb(i2,17)+r_urb(i2,27)+.84*r_urb(i2,28)+r_urb(i2,29)  &                                 
             +1.5*r_urb(i2,30)+r_urb(i2,31)+r_urb(i2,32)+r_urb(i2,33)           &
             +1.5*r_urb(i2,34)+r_urb(i2,35)+.5*r_urb(i2,42)                  &
             +r_urb(i2,46)  !move from p_com(i2,ino2) by feng fan on 20140822  
                                            
      rl_urb(i2,6)= r1_urb(i2,20)                                                            



      p_urb(i2,7)= 0.0                                                              
      rl_urb(i2,7)= r1_urb(i2,6)+r1_urb(i2,13)+r1_urb(i2,14)+r1_urb(i2,19)     &
             +r1_urb(i2,32)+r1_urb(i2,33)+r1_urb(i2,34)+r1_urb(i2,35)         &
             +r2_urb(i2,45)+r2_urb(i2,46)+r2_urb(i2,47)+r2_urb(i2,48)+r2_urb(i2,49)  !move from rl_com(i2,ino3) by feng fan on 20140822

                                                    

      p_urb(i2,11)= 0.0                                                               
      rl_urb(i2,11)= r1_urb(i2,7)+r1_urb(i2,9)+r1_urb(i2,10)+r1_urb(i2,23)    
                                         

      p_urb(i2,14)= .12*r_urb(i2,7)+.33*r_urb(i2,9)+.6*r_urb(i2,10)      &
            +.08*r_urb(i2,23)+r_urb(i2,24)+.23*r_urb(i2,25) 
                                                  
      rl_urb(i2,14)= r1_urb(i2,1)+r1_urb(i2,3)+r1_urb(i2,5)                      &
            +r1_urb(i2,8)+r1_urb(i2,11)+r1_urb(i2,12)+r1_urb(i2,15)               &              
            +r1_urb(i2,16)+r1_urb(i2,18)+r1_urb(i2,21)+r1_urb(i2,25)+r1_urb(i2,26)                                     



      p_urb(i2,15)= r_urb(i2,4)+.22*r_urb(i2,7)+r_urb(i2,8)                 &
             +.26*r_urb(i2,9)+.22*r_urb(i2,10)                          &  
             +r_urb(i2,11)+r_urb(i2,12)+.2*r_urb(i2,15)                    &
             +.55*r_urb(i2,16)+.95*r_urb(i2,17)                         &
             +.6*r_urb(i2,18)+2*r_urb(i2,21)+r_urb(i2,22)+.76*r_urb(i2,23)    &                            
             +.9*r_urb(i2,24)+.9*r_urb(i2,27)+.76*r_urb(i2,28)             &                
             +.5*r_urb(i2,30)+.9*r_urb(i2,32)+.5*r_urb(i2,34)              &
             +.54*r_urb(i2,40)    
                                  
      rl_urb(i2,15)= r1_urb(i2,36)+r1_urb(i2,37)+r1_urb(i2,38)+r1_urb(i2,39)                                          



      p_urb(i2,17)= r_urb(i2,4)+r_urb(i2,6)+.24*r_urb(i2,7)       &
            +.31*r_urb(i2,9)+.3*r_urb(i2,10)                 &             
            +2*r_urb(i2,21)+r_urb(i2,22)+.69*r_urb(i2,23)                                           
      rl_urb(i2,17)= 0.0                                                               



      p_urb(i2,21)= r_urb(i2,2)+.07*r_urb(i2,9)+.1*r_urb(i2,10)                                         
      rl_urb(i2,21)= 0.0                                                            



      p_urb(i2,22)= .06*r_urb(i2,9)+.05*r_urb(i2,10)+.1*r_urb(i2,24)    &
              +.1*r_urb(i2,27)                                    &
              +.08*r_urb(i2,28)+.1*r_urb(i2,32)+.06*r_urb(i2,40)                                    
      rl_urb(i2,22)= 0.0                                                             



      p_urb(i2,23)= r_urb(i2,7)+1.56*r_urb(i2,8)+.57*r_urb(i2,9)        &
              +r_urb(i2,11)+r_urb(i2,21)                             &
              +.7*r_urb(i2,23)+r_urb(i2,29)+.5*r_urb(i2,30)             &
              +r_urb(i2,33)+.5*r_urb(i2,34)                          &
              +.7*r_urb(i2,41)+.5*r_urb(i2,42)                       &
              + 2.0*r_urb(i2,46) !move from p_com(i2,ihcho) by feng fan on 20140822 
                                             
      rl_urb(i2,23)= 0.0                                                             



      p_urb(i2,24)= .03*r_urb(i2,9)+.04*r_urb(i2,10)                                             
      rl_urb(i2,24)= 0.0                                                            



      p_urb(i2,28)= .22*r_urb(i2,8)+.47*r_urb(i2,9)+1.03*r_urb(i2,10)    &
              +r_urb(i2,11)                                        &
              +1.77*r_urb(i2,12)+.03*r_urb(i2,23)+.3*r_urb(i2,24)        &
              +.04*r_urb(i2,25)                                    &
              +.3*r_urb(i2,27)+.25*r_urb(i2,28)+.5*r_urb(i2,30)          &
              +.3*r_urb(i2,32)                                     &
              +.5*r_urb(i2,34)+.21*r_urb(i2,40)+.5*r_urb(i2,42)                                     
      rl_urb(i2,28)= 0.0                                                             



      p_urb(i2,29)= .52*r_urb(i2,7)+.22*r_urb(i2,9)                                              
      rl_urb(i2,29)= 0.0                                                            



      p_urb(i2,30)= .09*r_urb(i2,9)+.16*r_urb(i2,10)                                             
      rl_urb(i2,30)= 0.0                                                            



      p_urb(i2,31)= r_urb(i2,2)+r_urb(i2,4)+r_urb(i2,5)+r_urb(i2,6)         &
              +.13*r_urb(i2,9)+.19*r_urb(i2,10)                       &   
              +r_urb(i2,21)+r_urb(i2,22)+.62*r_urb(i2,23)                &
              +r_urb(i2,29)+r_urb(i2,33)                              & 
              +.7*r_urb(i2,41)                                                        
      rl_urb(i2,31)= 0.0                                                             



      p_urb(i2,35)= 1.1*r_urb(i2,16)    
                                                    
      rl_urb(i2,35)= r2_urb(i2,1)+r2_urb(i2,44)                                                       



      p_urb(i2,36)= .07*r_urb(i2,10)+.23*r_urb(i2,12)             &
              +.74*r_urb(i2,24)+.74*r_urb(i2,27)                  &       
              +.62*r_urb(i2,28)+.74*r_urb(i2,32)                  &
              +.57*r_urb(i2,40)+.15*r_urb(i2,41) 
                        
      rl_urb(i2,36)= r1_urb(i2,2)+r2_urb(i2,3)                                                       



      p_urb(i2,37)= .04*r_urb(i2,9)+.07*r_urb(i2,10)     &
              +.8*r_urb(i2,16)+.2*r_urb(i2,23)           &                 
              +.19*r_urb(i2,25)+.15*r_urb(i2,41)  
                                           
      rl_urb(i2,37)= r1_urb(i2,4)+r2_urb(i2,5)+r2_urb(i2,6)  

                                                

      p_urb(i2,38)= 0.0                                                              
      rl_urb(i2,38)= r2_urb(i2,7)+r2_urb(i2,8)   

                                                     
 
      p_urb(i2,39)= 0.0                                                             
      rl_urb(i2,39)= r2_urb(i2,9)+r2_urb(i2,11)+r2_urb(i2,13)   

                                             

      p_urb(i2,40)= 0.0                                                             
      rl_urb(i2,40)= r2_urb(i2,10)+r2_urb(i2,12)+r2_urb(i2,14)   

                                            

      p_urb(i2,41)= 0.0                                                              
      rl_urb(i2,41)= r2_urb(i2,15)   

                                                         

      p_urb(i2,42)= 0.0                                                              
      rl_urb(i2,42)= r2_urb(i2,16)   

                                                         

      p_urb(i2,43)= .12*r_urb(i2,15)+.05*r_urb(i2,16)  
                                           
      rl_urb(i2,43)= r2_urb(i2,18)+r2_urb(i2,19)  

                                                   

      p_urb(i2,44)= .8*r_urb(i2,15)+.45*r_urb(i2,16) 
                                              
      rl_urb(i2,44)= r2_urb(i2,17) 

                                                           

      p_urb(i2,45)= .4*r_urb(i2,18)+r_urb(i2,19)   
                                                
      rl_urb(i2,45)= r2_urb(i2,20)  

                                                          

      p_urb(i2,46)= .95*r_urb(i2,17)+.3*r_urb(i2,18)  
                                            
      rl_urb(i2,46)= r2_urb(i2,21)+r1_urb(i2,22)+r2_urb(i2,23) 

                                              

      p_urb(i2,47)= .05*r_urb(i2,17)+r_urb(i2,20)    &
              +.16*r_urb(i2,28)+.5*r_urb(i2,30)      &                        
              +.5*r_urb(i2,34)+r_urb(i2,38)+.5*r_urb(i2,42)   
                                      
      rl_urb(i2,47)= r2_urb(i2,26)+r1_urb(i2,27)  
                        
                           

      p_urb(i2,48)= r_urb(i2,36)+r_urb(i2,37)  
                                                   
      rl_urb(i2,48)= r1_urb(i2,24)+r2_urb(i2,25)  

                                                   

      p_urb(i2,49)= r_urb(i2,1)+.03*r_urb(i2,9)+.09*r_urb(i2,10) &
             +.77*r_urb(i2,25)  
                              
      rl_urb(i2,49)= r2_urb(i2,28)+r2_urb(i2,32)+r2_urb(i2,36)+r1_urb(i2,40)  

                                        

      p_urb(i2,50)= r_urb(i2,3)+.11*r_urb(i2,10)  
                                                
      rl_urb(i2,50)= r2_urb(i2,29)+r2_urb(i2,33)+r2_urb(i2,37)+r1_urb(i2,41)  

                                       

      p_urb(i2,51)= r_urb(i2,13)+r_urb(i2,14)+r_urb(i2,26) 
                                               
      rl_urb(i2,51)= r2_urb(i2,30)+r2_urb(i2,34)+r2_urb(i2,38)+r1_urb(i2,42) 

                                         

      p_urb(i2,52)= r_urb(i2,5)+r_urb(i2,8)+r_urb(i2,11)+r_urb(i2,12)          &
             +.08*r_urb(i2,15)                                       &
             +.5*r_urb(i2,16)+.6*r_urb(i2,18)+r_urb(i2,21)+.03*r_urb(i2,23)   &                             
             +.4*r_urb(i2,24)+.41*r_urb(i2,27)+.34*r_urb(i2,28)            &
             +.4*r_urb(i2,32)+.24*r_urb(i2,40)                          &
             +r_urb(i2,46) ! move from p_com to p_urb by feng fan on 20140822
                                                    
      rl_urb(i2,52)= r2_urb(i2,31)+r2_urb(i2,35)+r2_urb(i2,39)+r1_urb(i2,43) 

                                         

      p_urb(i2,53)= 1.06*r_urb(i2,9)+2.26*r_urb(i2,10)                  &
              +r_urb(i2,11)+2.23*r_urb(i2,12)                           &
              +1.98*r_urb(i2,24)+.42*r_urb(i2,25)+1.98*r_urb(i2,27)        &                         
              +1.68*r_urb(i2,28)+r_urb(i2,30)+1.98*r_urb(i2,32)+r_urb(i2,34)  &                             
              +1.25*r_urb(i2,40)+r_urb(i2,42) 
                                               
      rl_urb(i2,53)= r1_urb(i2,44)     


!!! the following sv1,sv2 added by li jie

      p_urb(i2,54) = r_urb(i2,15)*0.071+r_urb(i2,16)*0.038+r_urb(i2,1)*0.0718+r_urb(i2,18)*0.05
      rl_urb(i2,54) = 0.0

      p_urb(i2,55) = r_urb(i2,15)*0.138+r_urb(i2,16)*0.167
      rl_urb(i2,55) = 0.0     
     
    !endif ! if pmask
    enddo !i2

      return
      end



!!*********************************************************************
!! subroutine mapgas_com: maps cnn to and fro stot for the common 
!!                        gas-phase mechanism.
!!
!! nomenclature:
!! cnn       = full species concentration array.
!! stot      = subset of cnn. species concentration array to be supplied to
!!             lsodes. length of stot depends on the selected mechanism
!! iregime   = selected chemical regime (1-6)
!! imap      = 0 : map cnn to stot
!!           = 1 : map stot to cnn
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!------------------------------------------------------------------------
      subroutine mapgas_com(stot,imap)
      use gas_data
      dimension stot(VLEN,76)
!dir$ assume_aligned stot:64
#ifndef VEC_OPT
      emit(ih2so4)	=emission(kh2so4)
      emit(ihno3)	=emission(khno3)
      emit(ihcl)	=emission(khcl)
      emit(inh3)	=emission(knh3)
      emit(ino)		=emission(kno)
      emit(ino2)	=emission(kno2)
      emit(ino3)	=emission(kno3)
      emit(in2o5)	=emission(kn2o5)
      emit(ihono)	=emission(khono)
      emit(ihno4)	=emission(khno4)
      emit(io3)		=emission(ko3)
      emit(io1d)	=emission(ko1d)
      emit(io3p)	=emission(ko3p)
      emit(ioh)		=emission(koh)
      emit(iho2)	=emission(kho2)
      emit(ih2o2)	=emission(kh2o2)
      emit(ico)		=emission(kco)
      emit(iso2)	=emission(kso2)
      emit(ich4)	=emission(kch4)
      emit(ic2h6)	=emission(kc2h6)
      emit(ich3o2)	=emission(kch3o2)
      emit(iethp)	=emission(kethp)
      emit(ihcho)	=emission(khcho)
      emit(ich3oh)	=emission(kch3oh)
      emit(ianol)	=emission(kanol)
      emit(ich3ooh)	=emission(kch3ooh)
      emit(iethooh)	=emission(kethooh)
      emit(iald2)	=emission(kald2)
      emit(ihcooh)	=emission(khcooh)
      emit(ircooh)	=emission(krcooh)
      emit(ic2o3)	=emission(kc2o3)
      emit(ipan)	=emission(kpan)
      emit(idso4)       =emission(kdso4)
      emit(idno3)       =emission(kdno3)
#else
      do i=1,32
        emit(:,i) = emission(:,i)
      enddo
      emit(:,33) = emission(:,75)
      emit(:,34) = emission(:,76)
#endif

      if(imap.eq.0)then    ! map cnn into stot
#ifdef VEC_OPT
      do i=1,32
        stot(:,i) = cnn(:,i)
      enddo
      stot(:,33) = cnn(:,75)
      stot(:,34) = cnn(:,76)
#else
      stot(ih2so4)=cnn(kh2so4)
      stot(ihno3)	=cnn(khno3)
      stot(ihcl)	=cnn(khcl)
      stot(inh3)	=cnn(knh3)
      stot(ino)		=cnn(kno)
      stot(ino2)	=cnn(kno2)
      stot(ino3)	=cnn(kno3)
      stot(in2o5)	=cnn(kn2o5)
      stot(ihono)	=cnn(khono)
      stot(ihno4)	=cnn(khno4)
      stot(io3)		=cnn(ko3)
      stot(io1d)	=cnn(ko1d)
      stot(io3p)	=cnn(ko3p)
      stot(ioh)		=cnn(koh)
      stot(iho2)	=cnn(kho2)
      stot(ih2o2)	=cnn(kh2o2)
      stot(ico)		=cnn(kco)
      stot(iso2)	=cnn(kso2)
      stot(ich4)	=cnn(kch4)
      stot(ic2h6)	=cnn(kc2h6)
      stot(ich3o2)	=cnn(kch3o2)
      stot(iethp)	=cnn(kethp)
      stot(ihcho)	=cnn(khcho)
      stot(ich3oh)	=cnn(kch3oh)
      stot(ianol)	=cnn(kanol)
      stot(ich3ooh)	=cnn(kch3ooh)
      stot(iethooh)	=cnn(kethooh)
      stot(iald2)	=cnn(kald2)
      stot(ihcooh)	=cnn(khcooh)
      stot(ircooh)	=cnn(krcooh)
      stot(ic2o3)	=cnn(kc2o3)
      stot(ipan)	=cnn(kpan)
      stot(idso4)       =cnn(kdso4)
      stot(idno3)       =cnn(kdno3)
#endif
      else                 ! map stot back into cnn
#ifdef VEC_OPT
      do i=1,32
        cnn(:,i) = stot(:,i)
      enddo
      cnn(:,75) = stot(:,33)
      cnn(:,76) = stot(:,34)
#else
      cnn(kh2so4)	=stot(ih2so4)
      cnn(khno3)	=stot(ihno3)
      cnn(khcl)		=stot(ihcl)
      cnn(knh3)		=stot(inh3)
      cnn(kno)		=stot(ino)
      cnn(kno2)		=stot(ino2)
      cnn(kno3)		=stot(ino3)
      cnn(kn2o5)	=stot(in2o5)
      cnn(khono)	=stot(ihono)
      cnn(khno4)	=stot(ihno4)
      cnn(ko3)		=stot(io3)
      cnn(ko1d)		=stot(io1d)
      cnn(ko3p)		=stot(io3p)
      cnn(koh)		=stot(ioh)
      cnn(kho2)		=stot(iho2)
      cnn(kh2o2)	=stot(ih2o2)
      cnn(kco)		=stot(ico)
      cnn(kso2)		=stot(iso2)
      cnn(kch4)		=stot(ich4)
      cnn(kc2h6)	=stot(ic2h6)
      cnn(kch3o2)	=stot(ich3o2)
      cnn(kethp)	=stot(iethp)
      cnn(khcho)	=stot(ihcho)
      cnn(kch3oh)	=stot(ich3oh)
      cnn(kanol)	=stot(ianol)
      cnn(kch3ooh)	=stot(ich3ooh)
      cnn(kethooh)	=stot(iethooh)
      cnn(kald2)	=stot(iald2)
      cnn(khcooh)	=stot(ihcooh)
      cnn(krcooh)	=stot(ircooh)
      cnn(kc2o3)	=stot(ic2o3)
      cnn(kpan)		=stot(ipan)
      cnn(kdso4)        =stot(idso4)
      cnn(kdno3)        =stot(idno3)
#endif
      endif

      return
      end





!!***********************************************************************
!! subroutine mapgas_bio: maps cnn to and fro stot for the biogeni!! 
!!                        gas-phase mechanism.
!!
!! nomenclature:
!! cnn       = full species concentration array.
!! stot      = subset of cnn. species concentration array to be supplied to
!!             lsodes. length of stot depends on the selected mechanism
!! iregime   = selected chemical regime (1-6)
!! imap      = 0 : map cnn to stot
!!           = 1 : map stot to cnn
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!-------------------------------------------------------------------------
      subroutine mapgas_bio(stot,imap)
      use gas_data
      dimension stot(VLEN,76)
!dir$ assume_aligned stot:64
#ifdef VEC_OPT
      emit(:,56)	=emission(:,52)
      emit(:,57)	=emission(:,53)
      emit(:,58)	=emission(:,54)
      emit(:,59)	=emission(:,55)
      emit(:,60)	=emission(:,56)
      emit(:,61)       =emission(:,68)
      emit(:,62)        =emission(:,71)
      emit(:,63)        =emission(:,72)
      emit(:,64)        =emission(:,73)
      emit(:,65)        =emission(:,74)

      if(imap.eq.0)then    ! map cnn into stot
      stot(:,56)	=cnn(:,52)
      stot(:,57)	=cnn(:,53)
      stot(:,58)	=cnn(:,54)
      stot(:,59)	=cnn(:,55)
      stot(:,60)	=cnn(:,56)
      stot(:,61)       =cnn(:,68)
      stot(:,62)        =cnn(:,71)
      stot(:,63)        =cnn(:,72)
      stot(:,64)        =cnn(:,73)
      stot(:,65)        =cnn(:,74)


      else                 ! map stot back into cnn
      cnn(:,52)	=stot(:,56)
      cnn(:,53)	=stot(:,57)
      cnn(:,54)	=stot(:,58)
      cnn(:,55)	=stot(:,59)
      cnn(:,56)	=stot(:,60)
      cnn(:,68)        =stot(:,61)
      cnn(:,71)         =stot(:,62)
      cnn(:,72)         =stot(:,63)
      cnn(:,73)         =stot(:,64)
      cnn(:,74)         =stot(:,65) 

      endif
#else
      emit(iisop)	=emission(kisop)
      emit(iisoprd)	=emission(kisoprd)
      emit(iisopp)	=emission(kisopp)
      emit(iisopn)	=emission(kisopn)
      emit(iisopo2)	=emission(kisopo2)
      emit(iterp)       =emission(kterp)
      emit(isv3)        =emission(ksv3)
      emit(isv4)        =emission(ksv4)
      emit(isv5)        =emission(ksv5)
      emit(isv6)        =emission(ksv6)

      if(imap.eq.0)then    ! map cnn into stot
      stot(iisop)	=cnn(kisop)
      stot(iisoprd)	=cnn(kisoprd)
      stot(iisopp)	=cnn(kisopp)
      stot(iisopn)	=cnn(kisopn)
      stot(iisopo2)	=cnn(kisopo2)
      stot(iterp)       =cnn(kterp)
      stot(isv3)        =cnn(ksv3)
      stot(isv4)        =cnn(ksv4)
      stot(isv5)        =cnn(ksv5)
      stot(isv6)        =cnn(ksv6)


      else                 ! map stot back into cnn
      cnn(kisop)	=stot(iisop)
      cnn(kisoprd)	=stot(iisoprd)
      cnn(kisopp)	=stot(iisopp)
      cnn(kisopn)	=stot(iisopn)
      cnn(kisopo2)	=stot(iisopo2)
      cnn(kterp)        =stot(iterp)
      cnn(ksv3)         =stot(isv3)
      cnn(ksv4)         =stot(isv4)
      cnn(ksv5)         =stot(isv5)
      cnn(ksv6)         =stot(isv6) 

      endif
#endif
      return
      end










!!*********************************************************************
!! subroutine mapgas_mar: maps cnn to and fro stot for the marine
!!                        gas-phase mechanism.
!!
!! nomenclature:
!! cnn       = full species concentration array.
!! stot      = subset of cnn. species concentration array to be supplied to
!!             lsodes. length of stot depends on the selected mechanism
!! iregime   = selected chemical regime (1-6)
!! imap      = 0 : map cnn to stot
!!           = 1 : map stot to cnn
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!------------------------------------------------------------------------
      subroutine mapgas_mar(stot,imap)
      use gas_data
      dimension stot(VLEN,76)
!dir$ assume_aligned stot:64
#ifdef VEC_OPT
      emit(:,idms)	=emission(:,57)
      emit(:,imsa)	=emission(:,58)
      emit(:,idmso)	=emission(:,59)
      emit(:,idmso2)	=emission(:,60)
      emit(:,ich3so2h)	=emission(:,61)
      emit(:,ich3sch2oo)	=emission(:,62)
      emit(:,ich3so2)	=emission(:,63)
      emit(:,ich3so3)	=emission(:,64)
      emit(:,ich3so2oo)	=emission(:,65)
      emit(:,ich3so2ch2oo)=emission(:,66)
      emit(:,isulfhox)	=emission(:,67)

      if(imap.eq.0)then    ! map cnn into stot
      stot(:,idms)	=cnn(:,57)
      stot(:,imsa)	=cnn(:,58)
      stot(:,idmso)	=cnn(:,59)
      stot(:,idmso2)	=cnn(:,60)
      stot(:,ich3so2h)	=cnn(:,61)
      stot(:,ich3sch2oo)	=cnn(:,62)
      stot(:,ich3so2)	=cnn(:,63)
      stot(:,ich3so3)	=cnn(:,64)
      stot(:,ich3so2oo)	=cnn(:,65)
      stot(:,ich3so2ch2oo)=cnn(:,66)
      stot(:,isulfhox)	=cnn(:,67)


      else                 ! map stot back into cnn
      cnn(:,57)		=stot(:,idms)
      cnn(:,58)		=stot(:,imsa)
      cnn(:,59)	=stot(:,idmso)
      cnn(:,60)	=stot(:,idmso2)
      cnn(:,61)	=stot(:,ich3so2h)
      cnn(:,62)	=stot(:,ich3sch2oo)
      cnn(:,63)	=stot(:,ich3so2)
      cnn(:,64)	=stot(:,ich3so3)
      cnn(:,65)	=stot(:,ich3so2oo)
      cnn(:,66)	=stot(:,ich3so2ch2oo)
      cnn(:,67)	=stot(:,isulfhox)
      endif
#else
      emit(idms)	=emission(kdms)
      emit(imsa)	=emission(kmsa)
      emit(idmso)	=emission(kdmso)
      emit(idmso2)	=emission(kdmso2)
      emit(ich3so2h)	=emission(kch3so2h)
      emit(ich3sch2oo)	=emission(kch3sch2oo)
      emit(ich3so2)	=emission(kch3so2)
      emit(ich3so3)	=emission(kch3so3)
      emit(ich3so2oo)	=emission(kch3so2oo)
      emit(ich3so2ch2oo)=emission(kch3so2ch2oo)
      emit(isulfhox)	=emission(ksulfhox)

      if(imap.eq.0)then    ! map cnn into stot
      stot(idms)	=cnn(kdms)
      stot(imsa)	=cnn(kmsa)
      stot(idmso)	=cnn(kdmso)
      stot(idmso2)	=cnn(kdmso2)
      stot(ich3so2h)	=cnn(kch3so2h)
      stot(ich3sch2oo)	=cnn(kch3sch2oo)
      stot(ich3so2)	=cnn(kch3so2)
      stot(ich3so3)	=cnn(kch3so3)
      stot(ich3so2oo)	=cnn(kch3so2oo)
      stot(ich3so2ch2oo)=cnn(kch3so2ch2oo)
      stot(isulfhox)	=cnn(ksulfhox)


      else                 ! map stot back into cnn
      cnn(kdms)		=stot(idms)
      cnn(kmsa)		=stot(imsa)
      cnn(kdmso)	=stot(idmso)
      cnn(kdmso2)	=stot(idmso2)
      cnn(kch3so2h)	=stot(ich3so2h)
      cnn(kch3sch2oo)	=stot(ich3sch2oo)
      cnn(kch3so2)	=stot(ich3so2)
      cnn(kch3so3)	=stot(ich3so3)
      cnn(kch3so2oo)	=stot(ich3so2oo)
      cnn(kch3so2ch2oo)	=stot(ich3so2ch2oo)
      cnn(ksulfhox)	=stot(isulfhox)
      endif
#endif
      return
      end









      subroutine mapgas_all(stot,imap)
      use gas_data
      implicit none
      real, dimension(VLEN,76) :: stot
      integer :: imap
      integer :: i
!dir$ assume_aligned stot:64
      do i=1,32
        emit(:,i) = emission(:,i)
      enddo
      emit(:,33) = emission(:,75)
      emit(:,34) = emission(:,76)
      do i = 33, 51
        emit(:,i+2) = emission(:,i)
      end do
      emit(:,54)        =emission(:,69)
      emit(:,55)        =emission(:,70)
      emit(:,56)	=emission(:,52)
      emit(:,57)	=emission(:,53)
      emit(:,58)	=emission(:,54)
      emit(:,59)	=emission(:,55)
      emit(:,60)	=emission(:,56)
      emit(:,61)       =emission(:,68)
      emit(:,62)        =emission(:,71)
      emit(:,63)        =emission(:,72)
      emit(:,64)        =emission(:,73)
      emit(:,65)        =emission(:,74)
      emit(:,idms)	=emission(:,57)
      emit(:,imsa)	=emission(:,58)
      emit(:,idmso)	=emission(:,59)
      emit(:,idmso2)	=emission(:,60)
      emit(:,ich3so2h)	=emission(:,61)
      emit(:,ich3sch2oo)	=emission(:,62)
      emit(:,ich3so2)	=emission(:,63)
      emit(:,ich3so3)	=emission(:,64)
      emit(:,ich3so2oo)	=emission(:,65)
      emit(:,ich3so2ch2oo)=emission(:,66)
      emit(:,isulfhox)	=emission(:,67)

      if(imap.eq.0)then    ! map cnn into stot
      do i=1,32
        stot(:,i) = cnn(:,i)
      enddo
      stot(:,33) = cnn(:,75)
      stot(:,34) = cnn(:,76)
      do i = 33,51
        stot(:,i+2) = cnn(:,i)
      enddo
      stot(:,54)        =cnn(:,69)
      stot(:,55)        =cnn(:,70)
      stot(:,56)	=cnn(:,52)
      stot(:,57)	=cnn(:,53)
      stot(:,58)	=cnn(:,54)
      stot(:,59)	=cnn(:,55)
      stot(:,60)	=cnn(:,56)
      stot(:,61)       =cnn(:,68)
      stot(:,62)        =cnn(:,71)
      stot(:,63)        =cnn(:,72)
      stot(:,64)        =cnn(:,73)
      stot(:,65)        =cnn(:,74)
      stot(:,idms)	=cnn(:,57)
      stot(:,imsa)	=cnn(:,58)
      stot(:,idmso)	=cnn(:,59)
      stot(:,idmso2)	=cnn(:,60)
      stot(:,ich3so2h)	=cnn(:,61)
      stot(:,ich3sch2oo)	=cnn(:,62)
      stot(:,ich3so2)	=cnn(:,63)
      stot(:,ich3so3)	=cnn(:,64)
      stot(:,ich3so2oo)	=cnn(:,65)
      stot(:,ich3so2ch2oo)=cnn(:,66)
      stot(:,isulfhox)	=cnn(:,67)

 
      else                 ! map stot back into cnn
      do i=1,32
        cnn(:,i) = stot(:,i)
      enddo
      cnn(:,75) = stot(:,33)
      cnn(:,76) = stot(:,34)
      do i = 33,51
        cnn(:,i) = stot(:,i+2)
      enddo
      cnn(:,69)         =stot(:,54)
      cnn(:,70)         =stot(:,55)
      cnn(:,52)	=stot(:,56)
      cnn(:,53)	=stot(:,57)
      cnn(:,54)	=stot(:,58)
      cnn(:,55)	=stot(:,59)
      cnn(:,56)	=stot(:,60)
      cnn(:,68)        =stot(:,61)
      cnn(:,71)         =stot(:,62)
      cnn(:,72)         =stot(:,63)
      cnn(:,73)         =stot(:,64)
      cnn(:,74)         =stot(:,65) 
      cnn(:,57)		=stot(:,idms)
      cnn(:,58)		=stot(:,imsa)
      cnn(:,59)	=stot(:,idmso)
      cnn(:,60)	=stot(:,idmso2)
      cnn(:,61)	=stot(:,ich3so2h)
      cnn(:,62)	=stot(:,ich3sch2oo)
      cnn(:,63)	=stot(:,ich3so2)
      cnn(:,64)	=stot(:,ich3so3)
      cnn(:,65)	=stot(:,ich3so2oo)
      cnn(:,66)	=stot(:,ich3so2ch2oo)
      cnn(:,67)	=stot(:,isulfhox)

      endif
      return
      end


!!*********************************************************************
!! subroutine mapgas_urb: maps cnn to and fro stot for the urban 
!!                        gas-phase mechanism.
!!
!! nomenclature:
!! cnn       = full species concentration array.
!! stot      = subset of cnn. species concentration array to be supplied to
!!             lsodes. length of stot depends on the selected mechanism
!! iregime   = selected chemical regime (1-6)
!! imap      = 0 : map cnn to stot
!!           = 1 : map stot to cnn
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!------------------------------------------------------------------------
      subroutine mapgas_urb(stot,imap)
      use gas_data
      dimension stot(VLEN,76)
!dir$ assume_aligned stot:64
#ifdef VEC_OPT
      do i = 33, 51
        emit(:,i+2) = emission(:,i)
      end do
      emit(:,54)        =emission(:,69)
      emit(:,55)        =emission(:,70)
#else
      emit(ipar)	=emission(kpar)
      emit(iaone)	=emission(kaone)
      emit(imgly)	=emission(kmgly)
      emit(ieth)	=emission(keth)
      emit(iolet)	=emission(kolet)
      emit(iolei)	=emission(kolei)
      emit(itol)	=emission(ktol)
      emit(ixyl)	=emission(kxyl)
      emit(icres)	=emission(kcres)
      emit(ito2)	=emission(kto2)
      emit(icro)	=emission(kcro)
      emit(iopen)	=emission(kopen)
      emit(ionit)	=emission(konit)
      emit(irooh)	=emission(krooh)
      emit(iro2)	=emission(kro2)
      emit(iano2)	=emission(kano2)
      emit(inap)	=emission(knap)
      emit(ixo2)	=emission(kxo2)
      emit(ixpar)	=emission(kxpar)
      emit(isv1)        =emission(ksv1)
      emit(isv2)        =emission(ksv2)
#endif
      if(imap.eq.0)then    ! map cnn into stot
#ifdef KNL_OPT
      do i = 33,51
        stot(:,i+2) = cnn(:,i)
      enddo
      stot(:,54)        =cnn(:,69)
      stot(:,55)        =cnn(:,70)
#else
      stot(ipar)	=cnn(kpar)
      stot(iaone)	=cnn(kaone)
      stot(imgly)	=cnn(kmgly)
      stot(ieth)	=cnn(keth)
      stot(iolet)	=cnn(kolet)
      stot(iolei)	=cnn(kolei)
      stot(itol)	=cnn(ktol)
      stot(ixyl)	=cnn(kxyl)
      stot(icres)	=cnn(kcres)
      stot(ito2)	=cnn(kto2)
      stot(icro)	=cnn(kcro)
      stot(iopen)	=cnn(kopen)
      stot(ionit)	=cnn(konit)
      stot(irooh)	=cnn(krooh)
      stot(iro2)	=cnn(kro2)
      stot(iano2)	=cnn(kano2)
      stot(inap)	=cnn(knap)
      stot(ixo2)	=cnn(kxo2)
      stot(ixpar)	=cnn(kxpar)
      stot(isv1)        =cnn(ksv1)
      stot(isv2)        =cnn(ksv2)
#endif
 
      else                 ! map stot back into cnn
#ifdef VEC_OPT
      do i = 33,51
        cnn(:,i) = stot(:,i+2)
      enddo
      cnn(:,69)         =stot(:,54)
      cnn(:,70)         =stot(:,55)
#else
      cnn(kpar)		=stot(ipar)
      cnn(kaone)	=stot(iaone)
      cnn(kmgly)	=stot(imgly)
      cnn(keth)		=stot(ieth)
      cnn(kolet)	=stot(iolet)
      cnn(kolei)	=stot(iolei)
      cnn(ktol)		=stot(itol)
      cnn(kxyl)		=stot(ixyl)
      cnn(kcres)	=stot(icres)
      cnn(kto2)		=stot(ito2)
      cnn(kcro)		=stot(icro)
      cnn(kopen)	=stot(iopen)
      cnn(konit)	=stot(ionit)
      cnn(krooh)	=stot(irooh)
      cnn(kro2)		=stot(iro2)
      cnn(kano2)	=stot(iano2)
      cnn(knap)		=stot(inap)
      cnn(kxo2)		=stot(ixo2)
      cnn(kxpar)	=stot(ixpar)
      cnn(ksv1)         =stot(isv1)
      cnn(ksv2)         =stot(isv2)
#endif 
      endif
      return
      end







!!*********************************************************************
!! subroutine setgas_com: sets up gas-phase species indices for 
!! the selected mechanism.
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine setgas_com(ilast)
      use gas_data
#ifdef VEC_OPT
      integer :: ilast
      ilast  = 34
#else
      ih2so4	= 1
      ihno3	= 2
      ihcl	= 3
      inh3	= 4
      ino	= 5
      ino2	= 6
      ino3	= 7
      in2o5	= 8
      ihono	= 9
      ihno4	= 10
      io3	= 11
      io1d	= 12
      io3p	= 13
      ioh	= 14
      iho2	= 15
      ih2o2	= 16
      ico	= 17
      iso2	= 18
      ich4	= 19
      ic2h6	= 20
      ich3o2	= 21
      iethp	= 22
      ihcho	= 23
      ich3oh	= 24
      ianol	= 25
      ich3ooh	= 26
      iethooh	= 27
      iald2	= 28
      ihcooh	= 29
      ircooh	= 30
      ic2o3	= 31
      ipan	= 32
      idso4     = 33
      idno3     = 34

      ilast	= idno3
#endif
 
      return
      end






!!*********************************************************************
!! subroutine setgas_bio: sets up gas-phase species indices for 
!! the selected mechanism.
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine setgas_bio(ilast)
#ifdef KNL_OPT
#else
      use gas_data
#endif

      iisop	= ilast + 1 !56
      iisoprd	= ilast + 2 !57
      iisopp	= ilast + 3 !58
      iisopn	= ilast + 4 !59
      iisopo2	= ilast + 5 !60
      iterp     = ilast + 6 !61
      isv3      = ilast + 7 !62
      isv4      = ilast + 8 !63
      isv5      = ilast + 9 !64
      isv6      = ilast + 10 !65
      
      ilast	= isv6 ! 34+21+10=65

      return
      end








!!*********************************************************************
!! subroutine setgas_mar: sets up gas-phase species indices for 
!! the selected mechanism.
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine setgas_mar(ilast)
      use gas_data

      idms         = ilast + 1
      imsa         = ilast + 2
      idmso        = ilast + 3
      idmso2       = ilast + 4
      ich3so2h     = ilast + 5
      ich3sch2oo   = ilast + 6
      ich3so2      = ilast + 7
      ich3so3      = ilast + 8
      ich3so2oo    = ilast + 9
      ich3so2ch2oo = ilast + 10
      isulfhox     = ilast + 11

      ilast	   = isulfhox

      return
      end












!!*********************************************************************
!! subroutine setgas_urb: sets up gas-phase species indices for 
!! the selected mechanism.
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine setgas_urb(ilast)
#ifdef KNL_OPTXX
      ilast = 34+21
#else
      ipar	= ilast + 1
      iaone	= ilast + 2
      imgly	= ilast + 3
      ieth	= ilast + 4
      iolet	= ilast + 5
      iolei	= ilast + 6
      itol	= ilast + 7
      ixyl	= ilast + 8
      icres	= ilast + 9
      ito2	= ilast + 10
      icro	= ilast + 11
      iopen	= ilast + 12
      ionit	= ilast + 13
      irooh	= ilast + 14
      iro2	= ilast + 15
      iano2	= ilast + 16
      inap	= ilast + 17
      ixo2	= ilast + 18
      ixpar	= ilast + 19
      isv1      = ilast + 20
      isv2      = ilast + 21
      
      ilast	= isv2
#endif
      return
      end


!!************************************************************************
!! subroutine gasrateconstants: generates thermal rate coefficients 
!!                   for the selected mechanism
!! nomenclature:
!! rk_com    = reaction rate constants for common mechanism (mole!!c!!s)
!! rk_urb    = reaction rate constants for hc1 mechanism    (mole!!c!!s)
!! rk_bio    = reaction rate constants for hc2 mechanism    (mole!!c!!s)
!! rk_mar    = reaction rate constants for marine mechanism (mole!!c!!s)
!! te        = ambient atmospheri!! temperature (k)
!! iregime = selected mechanism for the current chemical regime (1-6) 
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine gasrateconstants(cbmzobj,factcld)
      !use gas_data,only:iregime,msolar
      use gas_data
      implicit none
      real, dimension(VLEN) :: factcld
      type(cbmztype) cbmzobj
      integer :: i2
      logical :: has_t

      if(msolar.eq.1)then
       !call solarzenithangle ! calculates cos_sza
       call photoconstants_solar(cbmzobj,factcld)! natural diurnal variation
      elseif(msolar.eq.2)then
       call photoconstants_fixed(cbmzobj) ! artificial as in a smog chamber
      endif

      call gasrateconstants_het

#ifdef KNL_OPT
      cbmzobj%bmask(:) = cbmzobj%pmask(:)
      call gasrateconstants_com(cbmzobj)

!      !dir$ simd reduction(.or.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ne.1 .and. cbmzobj%iregime(i2).ne.4 .and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) &
      call gasrateconstants_urb(cbmzobj)

      !!dir$ simd reduction(.and.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = (cbmzobj%iregime(i2).eq.3 .or. cbmzobj%iregime(i2).eq.6).and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) &
      call gasrateconstants_bio(cbmzobj)

!      has_t = .false.
!      !dir$ simd reduction(.or.:has_t)
!      do i2 = 1, VLEN
!        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ge.4.and. cbmzobj%pmask(i2)
!        has_t = has_t .or. cbmzobj%bmask(i2)
!      enddo !i2
      if (has_t) &
      call gasrateconstants_mar(cbmzobj)
#else
      goto (1,2,3,4,5,6), cbmzobj%iregime

1     call gasrateconstants_com(cbmzobj)
      return

2     call gasrateconstants_com(cbmzobj)
      call gasrateconstants_urb(cbmzobj)
      return

3     call gasrateconstants_com(cbmzobj)
      call gasrateconstants_urb(cbmzobj)
      call gasrateconstants_bio(cbmzobj)
      return

4     call gasrateconstants_com(cbmzobj)
      call gasrateconstants_mar(cbmzobj)
      return
 
5     call gasrateconstants_com(cbmzobj)
      call gasrateconstants_urb(cbmzobj)
      call gasrateconstants_mar(cbmzobj)
      return

6     call gasrateconstants_com(cbmzobj)
      call gasrateconstants_urb(cbmzobj)
      call gasrateconstants_bio(cbmzobj)
      call gasrateconstants_mar(cbmzobj)
      return
#endif

      end




!!************************************************************************
!! subroutine gasrate: calculates reaction rates for the selected mechanism
!!
!! nomenclature:
!! r_com, r_urb, r_bio, r_mar = reaction rates (molec/cc/sec)
!! rk_com,rk_urb,rk_bio,rk_mar= rate constants in appropriate units
!! s                          = species concentrations (molec/cc)
!! o2                         = oxygen concentration   (molec/cc)
!! cair_mlc!! (used for m)      = air concentration      (molec/cc)
!! h2o                        = water vapor            (molec/cc)
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!-------------------------------------------------------------------------
      !subroutine gasrates(cbmzobj,s,r_com,r1_com,r2_com,r_urb,r1_urb,r2_urb,r_het,r_bio,r1_bio,r2_bio)
      subroutine gasrates(cbmzobj,s,r_com,r1_com,r2_com,r_urb,r1_urb,r2_urb,r_het)
      !use gas_data, only:cbmztype,76,76
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
      real, intent(out) :: r_com(VLEN,76),r1_com(VLEN,76),r2_com(VLEN,76)
      real :: s(VLEN,76),r_urb(VLEN,76),r1_urb(VLEN,76),r2_urb(VLEN,76),r_het(VLEN,76)
      !real :: r_bio(80),r1_bio(80),r2_bio(80)
!dir$ assume_aligned s:64,r_com:64,r1_com:64,r2_com:64,r_urb:64,r1_urb:64,r2_urb:64,r_het:64
      integer :: i2
      logical :: has_t

#ifdef KNL_OPT
      cbmzobj%bmask(:) = cbmzobj%pmask(:)
      call gasrates_het(cbmzobj,s,r_het)
      call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)

!      !dir$ simd reduction(.or.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ne.1 .and. cbmzobj%iregime(i2).ne.4 .and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) &
      call gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)

!      !dir$ simd reduction(.or.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = (cbmzobj%iregime(i2).eq.3 .or. cbmzobj%iregime(i2).eq.6).and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) &
      call gasrates_bio(s)

!      !dir$ simd reduction(.or.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ge.4.and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) &
      call gasrates_mar(cbmzobj,s)

#else
      goto (1,2,3,4,5,6), cbmzobj%iregime
1     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      return

2     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      call gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)
      return

3     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      call gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)
      !call gasrates_bio(cbmzobj,s,r_bio,r1_bio,r2_bio)
      call gasrates_bio(s)
      return

4     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      call gasrates_mar(s)
      return
 
5     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      call gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)
      call gasrates_mar(s)
      return

6     call gasrates_com(cbmzobj,s,r_com,r1_com,r2_com)
      call gasrates_urb(cbmzobj,s,r_urb,r1_urb,r2_urb)
      !call gasrates_bio(cbmzobj,s,r_bio,r1_bio,r2_bio)
      call gasrates_bio(s)
      call gasrates_mar(s)
      return
#endif
      end




!!*********************************************************************
!! subroutine mapgasspecies: maps cnn to and fro stot for the selected 
!!                           gas-phase mechanism.
!!
!! nomenclature:
!! cnn       = full species concentration array.
!! stot      = subset of cnn. species concentration array to be supplied to
!!             lsodes. length of stot depends on the selected mechanism
!! iregime   = selected chemical regime (1-6)
!! imap      = 0 : map cnn to stot
!!           = 1 : map stot to cnn
!! 
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!------------------------------------------------------------------------
      subroutine mapgasspecies(cbmzobj,stot,imap)
      !use gas_data,only:cbmztype
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
      real, dimension(VLEN,76) :: stot
!dir$ assume_aligned stot:64
      integer :: imap

#ifdef KNL_OPT
      call mapgas_all(stot,imap)

#else
      goto (1,2,3,4,5,6), cbmzobj%iregime
1     call mapgas_com(stot,imap)
      return


2     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      return


3     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      return


4     call mapgas_com(stot,imap)
      call mapgas_mar(stot,imap)
      return

 
5     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_mar(stot,imap)
      return


6     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      call mapgas_mar(stot,imap)
      return
#endif
      end

!!*********************************************************************
!! subroutine setgasindices: sets up gas-phase species indices for 
!! the selected mechanism.
!!
!! input: iregime    = 1     : com
!!                   = 2     : com + urb
!!                   = 3     : com + urb + bio
!!                   = 4     : com + mar
!!                   = 5     : com + urb + mar
!!                   = 6     : com + urb + bio + mar
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!------------------------------------------------------------------------
      subroutine setgasindices(cbmzobj)
      !use gas_data, only:cbmztype
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
      integer :: ilast
      ilast = 0
#ifdef KNL_OPT
      call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      call setgas_mar(ilast)

#else
      goto (1,2,3,4,5,6), cbmzobj%iregime


1     call setgas_com(ilast)
      return


2     call setgas_com(ilast)
      call setgas_urb(ilast)
      return


3     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      return


4     call setgas_com(ilast)
      call setgas_mar(ilast)
      return


5     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_mar(ilast)
      return


6     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      call setgas_mar(ilast)
      return
#endif

      end



!!************************************************************************
!! subroutine selectgasregime: selects an optimum combination of gas-phase
!!                             mechanisms 
!!
!! input : cnn       = full species concentrations array (molec/cc)
!!
!! output: iregime   = 1     : com
!!                   = 2     : com + urb
!!                   = 3     : com + urb + bio
!!                   = 4     : com + mar
!!                   = 5     : com + urb + mar
!!                   = 6     : com + urb + bio + mar    
!!         ngas      = number of gas-phase species in the selected mechanism
!!
!! author: rahul a. zaveri
!! date  : february 1996
!!
!!---------------------------------------------------------------------
#ifdef KNL_OPT
      subroutine selectgasregime(cbmzobj,ntot)
      !use gas_data,only:cbmztype,cnn,nreg1,nreg2,nreg3,nreg4,nreg5,nreg6
      !use gas_data,only:emission
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
      integer, dimension(VLEN) :: ntot
!dir$ assume_aligned ntot:64

      integer :: k
      integer, dimension(VLEN) :: m_com,m_urb,m_bio,m_mar
      real :: cutoff_conc
      integer :: i2
#else
      subroutine selectgasregime(ntot)
      use gas_data
#endif

      cutoff_conc = 5.e+6     ! [molec/cc]


!!---------------------------------------------
!! initialize regime flags to zero...
      m_com(:) = 1 ! 1 (always)
      m_urb(:) = 0	! 0 or 1
      m_bio(:) = 0	! 0 or 2
      m_mar(:) = 0	! 0 or 3


!! decide mechanism flags...
#ifdef VEC_OPT
!bad for perfomance ,BY WH
      do i2 = 1,VLEN
      do k = 33,67
        if ((cnn(i2,k) .gt. cutoff_conc)  .or.  (emission(i2,k) .gt. 0.0)) then
          if (k .lt. 52) then
            m_urb(i2)=1
          else if (k .gt. 56) then
            m_mar(i2)=3
          else
            m_bio(i2)=2
          endif
        endif
      enddo
      enddo !i2
#else
      do k = kpar, kxpar
      if( (cnn(k) .gt. cutoff_conc)  .or.  (emission(k) .gt. 0.0)    ) then 
       m_urb=1
      end if
      enddo

      do k = kisop, kisopo2
      if( (cnn(k) .gt. cutoff_conc)  .or.  (emission(k) .gt. 0.0)    ) then 
       m_bio=2
      end if
      enddo

      do k = kdms, ksulfhox
      if( (cnn(k) .gt. cutoff_conc)  .or.  (emission(k) .gt. 0.0)     ) then 
       m_mar=3
      end if
      enddo
#endif

    do i2=1,VLEN
      cbmzobj%iregime(i2) = m_com(i2) + m_urb(i2)*((2-m_bio(i2))/2) + m_bio(i2) + m_mar(i2)
      select case (cbmzobj%iregime(i2))
        case (1)
          ntot(i2) = 34 
        case (2)
          ntot(i2) = 55 
        case (3)
          ntot(i2) = 65 
        case (4)
          ntot(i2) = 45 
        case (5)
          ntot(i2) = 66 
        case (6)
          ntot(i2) = 76 
      end select
    enddo !i2
      end




!!*************************************************************************
!! subroutine ode_gas
!!
!! purpose: computes p and l.  dy/dt = p-ly.
!!          calls ode_com, ode_urb, ode_bio, ode_mar depending on the
!!          chemical regime (iregime)
!!
!! author: fan feng
!!
!!--------------------------------------------------------------------------

      subroutine ode_gas(cbmzobj,ntot,y,t_in)
      use gas_data
      implicit none
      type(cbmztype), target :: cbmzobj
      integer, dimension(VLEN) :: ntot
      real, dimension(VLEN,76) :: y
      real, dimension(VLEN) :: t_in
!dir$ assume_aligned ntot:64,y:64,t_in:64
      integer :: i, i2
      logical :: has_t

#ifdef KNL_OPT
!dir$ attributes align:64 :: r_com,r1_com,r2_com,p_com,rl_com
!dir$ attributes align:64 :: r_urb,r1_urb,r2_urb,p_urb,rl_urb
!dir$ attributes align:64 :: r_het,p_het,rl_het
!!!dir$ attributes align:64 :: r_bio,r1_bio,r2_bio,p_bio,rl_bio
      real :: r_com(VLEN,76),r1_com(VLEN,76)
      real :: r_urb(VLEN,76),r1_urb(VLEN,76),r2_urb(VLEN,76)

      real :: r_het(VLEN,76)
      real, dimension(:,:), pointer :: p_het,p_com,rl_com,r2_com,p_urb,rl_het,rl_urb,p_bio,rl_bio,p_mar,rl_mar
      !!real :: r_bio(80),r1_bio(80),r2_bio(80)
      !!real :: p_bio(80),rl_bio(80)

      p_het => cbmzobj%p_het
      p_com => cbmzobj%p_com
      rl_com => cbmzobj%rl_com
      r2_com => cbmzobj%r2_com
      p_urb => cbmzobj%p_urb
      rl_het => cbmzobj%rl_het
      rl_urb => cbmzobj%rl_urb
      p_bio => cbmzobj%p_bio
      rl_bio => cbmzobj%rl_bio
      p_mar => cbmzobj%p_mar
      rl_mar => cbmzobj%rl_mar
#endif
#ifndef KNL_OPT
      do i=1,76
        r_com(i) = 0.
      enddo

      do i=1,nrxn_urb
        r_urb(i) = 0.
      enddo

      do i=1,nrxn_bio
        r_bio(i) = 0.
      enddo

      do i=1,nrxn_mar
        r_mar(i) = 0.
      enddo

      do i=1,nrxn_het
        r_het(i) = 0.
      enddo
#endif

#ifdef KNL_OPT
      !call gasrates(cbmzobj,y,r_com,r1_com,r2_com,r_urb,r1_urb,r2_urb,r_het,r_bio,r1_bio,r2_bio)
      call gasrates(cbmzobj,y,r_com,r1_com,r2_com,r_urb,r1_urb,r2_urb,r_het)
#else
      call gasrates(y)
#endif

#ifndef KNL_OPT
      do i=1,76
        p_com(i) = 0.
        p_urb(i) = 0.
        p_bio(i) = 0.
        p_mar(i) = 0.
        p_het(i) = 0.
        total_p(i)= 0.

        rl_com(i) = 0.
        rl_urb(i) = 0.
        rl_bio(i) = 0.
        rl_mar(i) = 0.
        rl_het(i) = 0.
        total_l(i)= 0.
      enddo
#endif
   
#ifdef KNL_OPT
      ! use pmask
      call ode_com(cbmzobj,r_com,r1_com,r2_com,p_com,rl_com)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)
      do i=1,76
        do i2=1,VLEN
         if (cbmzobj%pmask(i2)) then
         total_p(i2,i)=     (p_com(i2,i)+p_het(i2,i)+emit(i2,i))
         total_l(i2,i)=     (rl_com(i2,i)+rl_het(i2,i))
         endif ! pmask
        enddo
      enddo

!      !dir$ simd reduction(.or.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ne.1 .and. cbmzobj%iregime(i2).ne.4 .and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) then
        call ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
        do i=1,76
        do i2=1,VLEN
         if (cbmzobj%bmask(i2)) then
           total_p(i2,i)= total_p(i2,i) +     (p_urb(i2,i))
           total_l(i2,i)= total_l(i2,i) +     (rl_urb(i2,i))
         endif ! bmask
        enddo
        enddo
      endif !has_t

!      !dir$ simd reduction(.and.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = (cbmzobj%iregime(i2).eq.3 .or. cbmzobj%iregime(i2).eq.6).and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) then
        call ode_bio(cbmzobj)
        do i=1,76
        do i2=1,VLEN
         if (cbmzobj%bmask(i2)) then
           total_p(i2,i)= total_p(i2,i) +     (p_bio(i2,i))
           total_l(i2,i)= total_l(i2,i) +     (rl_bio(i2,i))
         endif ! bmask
        enddo
        enddo
      endif !has_t

!      !dir$ simd reduction(.and.:has_t)
      do i2 = 1, VLEN
        cbmzobj%bmask(i2) = cbmzobj%iregime(i2).ge.4.and. cbmzobj%pmask(i2)
      enddo !i2
      has_t = any(cbmzobj%bmask)
      if (has_t) then
        call ode_mar(cbmzobj)
        do i=1,76
        do i2=1,VLEN
         if (cbmzobj%bmask(i2)) then
           total_p(i2,i)= total_p(i2,i) +     (p_mar(i2,i))
           total_l(i2,i)= total_l(i2,i) +     (rl_mar(i2,i))
         endif ! bmask
        enddo
        enddo
      endif !has_t

#else
      goto (1,2,3,4,5,6), cbmzobj%iregime

1     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)
      do i=1,34
         total_p(i)= real( dble(p_com(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_het(i)) ) 
      enddo

      return

!!---------------------------------------------------------
2     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)

      do i=1,55
         total_p(i)= real( dble(p_com(i)+p_urb(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_urb(i)+rl_het(i)) ) 
      enddo

      return

!!---------------------------------------------------------
3     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)
      !call ode_bio(r_bio,r1_bio,r2_bio,p_bio,rl_bio)
      call ode_bio

      do i=1,65
         total_p(i)= real( dble(p_com(i)+p_urb(i)+p_bio(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_urb(i)+rl_bio(i)+rl_het(i)) ) 
      enddo

      return

!!---------------------------------------------------------
4     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_mar(cbmzobj)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)

      do i=1,45
         total_p(i)= real( dble(p_com(i)+p_mar(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_mar(i)+rl_het(i)) ) 
      enddo

      return

!!---------------------------------------------------------
5     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
      call ode_mar(cbmzobj)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)

      do i=1,66
         total_p(i)= real( dble(p_com(i)+p_urb(i)+p_mar(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_urb(i)+rl_mar(i)+rl_het(i)) ) 
      enddo

      return

!!---------------------------------------------------------
6     call ode_com(r_com,r1_com,r2_com,p_com,rl_com)
      call ode_urb(r_urb,r1_urb,r2_urb,p_urb,rl_urb)
      !call ode_bio(r_bio,r1_bio,r2_bio,p_bio,rl_bio)
      call ode_bio
      call ode_mar(cbmzobj)
      call ode_het(r_het,cbmzobj%rk_het,p_het,rl_het)

      do i=1,76
         total_p(i)= real( dble(p_com(i)+p_urb(i)+p_bio(i)+p_mar(i)+p_het(i)) )+emit(i)
         total_l(i)= real( dble(rl_com(i)+rl_urb(i)+rl_bio(i)+rl_mar(i)+rl_het(i)) ) 
      enddo

      return
#endif

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!modified backward euler (mbe)
!!!---------------      subroutine mbe_solver(odefunc,neq,y,t)        ----
!!! renew dt_sec, y 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mbe_solver(cbmzobj,neq,y,t)
!! neq: number of equations

      !external odefunc
      !use gas_data,only:total_p,total_l,dt_sec,cri_dydt,decay_coeff
      !use gas_data,only:growth_coeff,max_stepsize,min_stepsize

      use gas_data
      implicit none
!dir$ assume_aligned neq:64,y:64,t:64
      integer, dimension(VLEN), intent(in) :: neq
      real, dimension(VLEN,76), intent(inout) :: y
      real, dimension(VLEN), intent(in) :: t
      type(cbmztype), intent(in) :: cbmzobj
      real :: dydt

      integer :: j, max_neq
      logical, dimension(VLEN) :: flag_change_fast_y
      integer :: i2


  call ode_gas(cbmzobj,neq,y,t)

  !do i2=1,VLEN
  !if (cbmzobj%pmask(i2)) then

 flag_change_fast_y=.false.
 max_neq = maxval(neq(:))
#ifdef KNL_OPT
  !dir$ assume_aligned y:64
  !!dir$ assume_aligned total_p:64
  !!dir$ assume_aligned total_l:64
  !!dir$ simd
#endif
  do j=1,max_neq
   do i2=1,VLEN
    if (j .le. neq(i2)) then

 !!---- p,l,y should be nonnegative ------------------------------

 !!--in fact, with any positive dt_sec, modified backward euler (mbe) can guarantee nonnegativity of y ---------
  y(i2,j)=(y(i2,j)+total_p(i2,j)*dt_sec(i2))/(1+total_l(i2,j)*dt_sec(i2)) !mbe
  dydt=total_p(i2,j)-total_l(i2,j)*y(i2,j) 
  if (dydt > abs(cri_dydt)) then
    flag_change_fast_y(i2)=.true.
  end if ! end of if
    
    endif ! j .le. neq(i2)
   enddo !i2
  end do !end of do j=1,neq
 


!!------------------  prepare for stepsize changing --------------------------------------------

 


!!------------------  change stepsize ---------------------------------------------------------
   do i2=1,VLEN

  if (flag_change_fast_y(i2) .eqv. .true.) then  !there exists value that changes fast .
      dt_sec(i2)=dt_sec(i2)/decay_coeff 
  else
    dt_sec(i2)=dt_sec(i2)*growth_coeff 
  end if

  
  if (dt_sec(i2)>max_stepsize) then
  dt_sec(i2)=max_stepsize 
  end if

  if (dt_sec(i2)<min_stepsize) then
  dt_sec(i2)=min_stepsize 
  end if
   enddo !i2

!endif ! if pmask
!enddo !i2

return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  subroutine odesolver(ntot,stot,t_in)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine odesolver(cbmzobj,ntot,stot,t_in)
      use gas_data
      implicit none
      type(cbmztype) cbmzobj
!dir$ assume_aligned ntot:64,stot:64,t_in:64
      integer, dimension(VLEN) :: ntot
      real, dimension(VLEN,76) :: stot! local species array
      real, dimension(VLEN) :: t_in

      !external ode_gas


      call mbe_solver(cbmzobj,ntot,stot,t_in)


      return
      end

    subroutine chemprod(cbmzobj,ig,i,j,k,delta,dt,ne)  ! calculate the production of species
    use gas_data
    implicit none
    type(cbmztype), target :: cbmzobj
    real, dimension(:,:), pointer :: rk_urb,rk_com
    real, dimension(:), pointer :: cair_mlc
 !include 'chm1.inc'
 !include 'gas1.inc'
        integer :: ig,i,j,k,ne
        integer :: i2
        real, dimension(VLEN) :: delta
!dir$ assume_aligned delta:64
        real :: delta0,dt
!ozone and its the production occurs in daytime         
!  F(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
    delta(:)=0.0
    rk_com => cbmzobj%rk_com
    rk_urb => cbmzobj%rk_urb
    cair_mlc => cbmzobj%cair_mlc
   IF(ig==11) THEN
    !dir$ simd private(delta0)
    do i2 = 1, VLEN
    if (cbmzobj%pmask(i2) .and. &
        cbmzobj%rk_photo(i2,1)>0.0) then
     delta0=60.*dt*(rk_com(i2,33)*cnn(i2,kno)*cnn(i2,kho2)& !k4[NO][HO2] in molec/cc
              +rk_com(i2,57)*cnn(i2,kch3o2)*cnn(i2,kno)& !k5[NO][CH3O2]) 
              +rk_urb(i2,17)*cnn(i2,kto2)*cnn(i2,kno)&   !k[TO2][NO]
              +rk_com(i2,58)*cnn(i2,kethp)*cnn(i2,kno)&  !K[ETHP][NO] 
              +rk_urb(i2,28)*cnn(i2,kro2)*cnn(i2,kno)&   !K[RO2][NO]
              +rk_com(i2,71)*cnn(i2,kc2o3)*cnn(i2,kno)&  !K[c2o3][no] 
              +rk_urb(i2,29)*cnn(i2,kano2)*cnn(i2,kno)&  !k[ano2][no]
              +rk_urb(i2,30)*cnn(i2,knap)*cnn(i2,kno)&  !k[nap][no]
              +rk_urb(i2,31)*cnn(i2,kxo2)*cnn(i2,kno)&  !k[xo2][no]
              +rk_bio(i2,8)*cnn(i2,kisopp)*cnn(i2,kno)& !k[isopp][no]
              +rk_bio(i2,9)*cnn(i2,kisopn)*cnn(i2,kno)& !k[isopn][no]
              +rk_bio(i2,10)*cnn(i2,kisopo2)*cnn(i2,kno)) !k[isopo2][no]

!-----------convert molec/cc to ppbv-----
    delta(i2) = delta0/cair_mlc(i2)*1.e+9     
    endif
    enddo !i2
   ENDIF
!      if(i==46.and.j==30.and.k==3.and.ig==11.and.ne==1)   print*,delta,'555'

    return
    end
   
#ifdef KNL_OPT
    subroutine chemprodloss(cbmzobj,myid,delta1,delta2,delta3,delta4,dt,i,j,k,ne) 
    use gas_data
    implicit none
    type(cbmztype), target :: cbmzobj
    real, pointer :: rk_urb(:,:), rk_com(:,:), cair_mlc(:),h2o(:)
    real, dimension(VLEN) :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11
     
#else
    subroutine chemprodloss(myid,delta1,delta2,delta3,delta4,dt,i,j,k,ne) 
    use gas_data
#endif
 !include 'chm1.inc'
 !include 'gas1.inc'
      integer myid,i,j,k,ne
      real, dimension(VLEN) :: delta1,delta2,delta3,delta4
!dir$ assume_aligned delta1:64,delta2:64,delta3:64,delta4:64
      real :: dt
!ozone production and loss,and loss due to nox,radicals from radical-radical reactions
!F(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
!L(o3)=k3[O1D][H2O]+k8[O3][HO2]+k7[O3][OH]+(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
!Ln(radical)=k[no2][oh]
!Lr(ridical)=k[RO2][HO2]+2k[HO2]^2
! Ln/(Ln+Lr)<0.4 is nox limited, >0,6 VOCs limited
! refrence :Kleinman, L. I., and Coauthors, 1997: Dependence of O3 production on
! NO and hydrocarbons in the troposphere. Geophys. Res. Lett., 24, 2299-2302.
#ifdef KNL_OPT
    rk_com => cbmzobj%rk_com
    rk_urb => cbmzobj%rk_urb
    cair_mlc => cbmzobj%cair_mlc
    h2o => cbmzobj%h2o
delta1=60.*dt*(rk_com(:,33)*cnn(:,kno)*cnn(:,kho2)& !k4[NO][HO2] in molec/cc
            +rk_com(:,57)*cnn(:,kch3o2)*cnn(:,kno)& !k5[NO][CH3O2])
            +rk_urb(:,17)*cnn(:,kto2)*cnn(:,kno)&   !k[TO2][NO]
            +rk_com(:,58)*cnn(:,kethp)*cnn(:,kno)&  !K[ETHP][NO]
            +rk_urb(:,28)*cnn(:,kro2)*cnn(:,kno)&   !K[RO2][NO]
            +rk_com(:,71)*cnn(:,kc2o3)*cnn(:,kno)&  !K[c2o3][no]
            +rk_urb(:,29)*cnn(:,kano2)*cnn(:,kno)&  !k[ano2][no]
            +rk_urb(:,30)*cnn(:,knap)*cnn(:,kno)&   !k[nap][no]
            +rk_urb(:,31)*cnn(:,kxo2)*cnn(:,kno)&   !k[xo2][no]
            +rk_bio(:,8)*cnn(:,kisopp)*cnn(:,kno)&  !k[isopp][no]
            +rk_bio(:,9)*cnn(:,kisopn)*cnn(:,kno)&  !k[isopn][no]
            +rk_bio(:,10)*cnn(:,kisopo2)*cnn(:,kno)&!k[isopo2][no]
            )/cair_mlc(:)*1.e+9

!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta1,'111'

      delta2= 60.*dt*(2.2e-10*cnn(:,ko1d)*h2o(:)& ! k3[O1D][h2o]
            +rk_com(:,21)*cnn(:,ko3)*cnn(:,kho2)&  ! k8[O3][HO2]
            +rk_com(:,20)*cnn(:,ko3)*cnn(:,koh)&   ! k7[O3][OH]
            +(rk_com(:,18)*cnn(:,kno)*cnn(:,ko3)*& !k10[NO][O3]
             rk_com(:,24)*cnn(:,kno2)*cnn(:,koh)/& !*k12[NO2][OH]
            (rk_com(:,1)*cnn(:,kno2)+rk_com(:,24)*cnn(:,kno2)*cnn(:,koh))&! k11[NO2]+k12[NO2][OH]
            )&      !(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
            )/cair_mlc(:)*1.e+9            
            
!             if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta2,'2222'

     delta3 = 60.*dt*(rk_com(:,24)*cnn(:,kno2)*cnn(:,koh))/cair_mlc(:)*1.e+9 !k[no2][oh]
!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta3,'333'
     
     delta4 = 60.*dt*(rk_com(:,61)*cnn(:,kch3o2)*cnn(:,kho2)& !k[HO2][ch3o2]
             +rk_com(:,62)*cnn(:,kethp)*cnn(:,kho2)&!k[HO2][ETHP]
             +rk_urb(:,36)*cnn(:,kro2)*cnn(:,kho2)&!k[HO2][RO2]
             +rk_com(:,73)*cnn(:,kc2o3)*cnn(:,kho2)&!k[HO2][C2O3]
             +rk_urb(:,37)*cnn(:,kano2)*cnn(:,kho2)&!k[HO2][ANO2]
             +rk_urb(:,38)*cnn(:,knap)*cnn(:,kho2)&!k[HO2][NAP]
             +rk_bio(:,11)*cnn(:,kisopp)*cnn(:,kho2)&!k[HO2][ISOPP]
             +rk_bio(:,12)*cnn(:,kisopn)*cnn(:,kho2)&!k[HO2][ISOPN]
             +rk_bio(:,13)*cnn(:,kisopo2)*cnn(:,kho2)&!k[HO2][ISOPO2]
             +rk_urb(:,39)*cnn(:,kxo2)*cnn(:,kho2)&!k[HO2][XO2]
             +2.*rk_com(:,31)*cnn(:,kho2)*cnn(:,kho2)& !2k[HO2]^2
             )/cair_mlc(:)*1.e+9

       d1=rk_com(:,61)*cnn(:,kch3o2)*cnn(:,kho2)
       d2=rk_com(:,62)*cnn(:,kethp)*cnn(:,kho2)
       d3=rk_urb(:,36)*cnn(:,kro2)*cnn(:,kho2)
       d4=rk_com(:,73)*cnn(:,kc2o3)*cnn(:,kho2)
       d5=rk_urb(:,37)*cnn(:,kano2)*cnn(:,kho2)
       d6=rk_urb(:,38)*cnn(:,knap)*cnn(:,kho2)
       d7=rk_bio(:,11)*cnn(:,kisopp)*cnn(:,kho2)
       d8=rk_bio(:,12)*cnn(:,kisopn)*cnn(:,kho2)
       d9=rk_bio(:,13)*cnn(:,kisopo2)*cnn(:,kho2)
       d10=rk_urb(:,39)*cnn(:,kxo2)*cnn(:,kho2)
       d11=2.*rk_com(:,31)*cnn(:,kho2)*cnn(:,kho2)
#else
delta1=60.*dt*(rk_com(33)*cnn(kno)*cnn(kho2)& !k4[NO][HO2] in molec/cc
            +rk_com(57)*cnn(kch3o2)*cnn(kno)& !k5[NO][CH3O2])
            +rk_urb(17)*cnn(kto2)*cnn(kno)&   !k[TO2][NO]
            +rk_com(58)*cnn(kethp)*cnn(kno)&  !K[ETHP][NO]
            +rk_urb(28)*cnn(kro2)*cnn(kno)&   !K[RO2][NO]
            +rk_com(71)*cnn(kc2o3)*cnn(kno)&  !K[c2o3][no]
            +rk_urb(29)*cnn(kano2)*cnn(kno)&  !k[ano2][no]
            +rk_urb(30)*cnn(knap)*cnn(kno)&   !k[nap][no]
            +rk_urb(31)*cnn(kxo2)*cnn(kno)&   !k[xo2][no]
            +rk_bio(8)*cnn(kisopp)*cnn(kno)&  !k[isopp][no]
            +rk_bio(9)*cnn(kisopn)*cnn(kno)&  !k[isopn][no]
            +rk_bio(10)*cnn(kisopo2)*cnn(kno)&!k[isopo2][no]
            )/cair_mlc*1.e+9

!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta1,'111'

      delta2= 60.*dt*(rk_com(12)*cnn(ko1d)*h2o& ! k3[O1D][H2O]
            +rk_com(21)*cnn(ko3)*cnn(kho2)&  ! k8[O3][HO2]
            +rk_com(20)*cnn(ko3)*cnn(koh)&   ! k7[O3][OH]
            +(rk_com(18)*cnn(kno)*cnn(ko3)*& !k10[NO][O3]
             rk_com(24)*cnn(kno2)*cnn(koh)/& !*k12[NO2][OH]
            (rk_com(1)*cnn(kno2)+rk_com(24)*cnn(kno2)*cnn(koh))&! k11[NO2]+k12[NO2][OH]
            )&      !(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
            )/cair_mlc*1.e+9            
            
!             if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta2,'2222'

     delta3 = 60.*dt*(rk_com(24)*cnn(kno2)*cnn(koh))/cair_mlc*1.e+9 !k[no2][oh]
!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta3,'333'
     
     delta4 = 60.*dt*(rk_com(61)*cnn(kch3o2)*cnn(kho2)& !k[HO2][ch3o2]
             +rk_com(62)*cnn(kethp)*cnn(kho2)&!k[HO2][ETHP]
             +rk_urb(36)*cnn(kro2)*cnn(kho2)&!k[HO2][RO2]
             +rk_com(73)*cnn(kc2o3)*cnn(kho2)&!k[HO2][C2O3]
             +rk_urb(37)*cnn(kano2)*cnn(kho2)&!k[HO2][ANO2]
             +rk_urb(38)*cnn(knap)*cnn(kho2)&!k[HO2][NAP]
             +rk_bio(11)*cnn(kisopp)*cnn(kho2)&!k[HO2][ISOPP]
             +rk_bio(12)*cnn(kisopn)*cnn(kho2)&!k[HO2][ISOPN]
             +rk_bio(13)*cnn(kisopo2)*cnn(kho2)&!k[HO2][ISOPO2]
             +rk_urb(39)*cnn(kxo2)*cnn(kho2)&!k[HO2][XO2]
             +2.*rk_com(31)*cnn(kho2)*cnn(kho2)& !2k[HO2]^2
             )/cair_mlc*1.e+9

       d1=rk_com(61)*cnn(kch3o2)*cnn(kho2)
       d2=rk_com(62)*cnn(kethp)*cnn(kho2)
       d3=rk_urb(36)*cnn(kro2)*cnn(kho2)
       d4=rk_com(73)*cnn(kc2o3)*cnn(kho2)
       d5=rk_urb(37)*cnn(kano2)*cnn(kho2)
       d6=rk_urb(38)*cnn(knap)*cnn(kho2)
       d7=rk_bio(11)*cnn(kisopp)*cnn(kho2)
       d8=rk_bio(12)*cnn(kisopn)*cnn(kho2)
       d9=rk_bio(13)*cnn(kisopo2)*cnn(kho2)
       d10=rk_urb(39)*cnn(kxo2)*cnn(kho2)
       d11=2.*rk_com(31)*cnn(kho2)*cnn(kho2)
#endif

     !       if(i==46.and.j==30.and.k==3.and.ne==1)&
     !        print*,d1,d2,d3,d4, '111'
     !        print*,d5,d6,d7,d8, '222'
     !        print*,d9,d10,d11 ,'333' 
     !         print*,delta3,delta4,60.*dt*d11/cair_mlc*1.e+9,60.*dt*(d1+d2+d3+d4)/cair_mlc*1.e+9,60.*dt*(d5+d6+d7+d8)/cair_mlc*1.e+9,&
     !                  60.*dt*(d9+d10)/cair_mlc*1.e+9
     return
    end

      subroutine chemope(cbmzobj,myid, ope, dt)
      use gas_data
      implicit none
      type(cbmztype), target :: cbmzobj
      real, dimension(:,:), pointer :: rk_urb,rk_com
      real, dimension(:), pointer :: cair_mlc,h2o
!c  to calculate the net ope (ozone production efficiency)for o3      
!c     ope = (p(o3)-p(o3))/p(noz) 
   
 !include 'chm1.inc'
 !include 'gas1.inc'
      integer myid
      real, dimension(VLEN) :: ope
!dir$ assume_aligned ope:64
      real    dt,po3,lo3,pnoz
      integer i2

!ozone production and loss,and loss due to nox,radicals from radical-radical reactions
!p(o3)=k4[no][ho2]+k5[no][ch3o2]+k6[ro2][no]
!l(o3)=k3[o1d][h2o]+k8[o3][ho2]+k7[o3][oh]+(k10[no][o3]*k12[no2][oh]/{k11[no2]+k12[no2][oh]})
      
    rk_com => cbmzobj%rk_com
    rk_urb => cbmzobj%rk_urb
    cair_mlc => cbmzobj%cair_mlc
    h2o => cbmzobj%h2o
    do i2=1,VLEN
      po3  = 60.*dt*(rk_com(i2,33)*cnn(i2,kno)*cnn(i2,kho2)& !k4[no][ho2] inmolec/cc
            +rk_com(i2,57)*cnn(i2,kch3o2)*cnn(i2,kno)& !k5[no][ch3o2])
            +rk_urb(i2,17)*cnn(i2,kto2)*cnn(i2,kno)&   !k[to2][no]
            +rk_com(i2,58)*cnn(i2,kethp)*cnn(i2,kno)&  !k[ethp][no]
            +rk_urb(i2,28)*cnn(i2,kro2)*cnn(i2,kno)&   !k[ro2][no]
            +rk_com(i2,71)*cnn(i2,kc2o3)*cnn(i2,kno)&  !k[c2o3][no]
            +rk_urb(i2,29)*cnn(i2,kano2)*cnn(i2,kno)&  !k[ano2][no]
            +rk_urb(i2,30)*cnn(i2,knap)*cnn(i2,kno)&   !k[nap][no]
            +rk_urb(i2,31)*cnn(i2,kxo2)*cnn(i2,kno)&   !k[xo2][no]
            +rk_bio(i2,8)*cnn(i2,kisopp)*cnn(i2,kno)&  !k[isopp][no]
            +rk_bio(i2,9)*cnn(i2,kisopn)*cnn(i2,kno)&  !k[isopn][no]
            +rk_bio(i2,10)*cnn(i2,kisopo2)*cnn(i2,kno)&!k[isopo2][no]
            )/cair_mlc(i2)*1.e+9

 
      lo3 = 60.*dt*(2.2e-10*cnn(i2,ko1d)*h2o(i2)& ! k3[o1d][h2o(i2)]
            +rk_com(i2,21)*cnn(i2,ko3)*cnn(i2,kho2)&  ! k8[o3][ho2]
            +rk_com(i2,20)*cnn(i2,ko3)*cnn(i2,koh)&   ! k7[o3][oh]
            +(rk_com(i2,18)*cnn(i2,kno)*cnn(i2,ko3)*& !k10[no][o3]
             rk_com(i2,24)*cnn(i2,kno2)*cnn(i2,koh)/& !*k12[no2][oh]
            (rk_com(i2,1)*cnn(i2,kno2)+rk_com(i2,24)*cnn(i2,kno2)*cnn(i2,koh))& !k11[no2]+k12[no2][oh]
            )&      !(i2,k10[no][o3]*k12[no2][oh]/{k11[no2]+k12[no2][oh]})
            )/cair_mlc(i2)*1.e+9

! ccccccccccc   noz production
! oh + no2       = hno3    
! no3 + no2      = n2o5      
! c2o3 + no2     = pan           
! cro + no2      = onit


      pnoz =  60.*dt*(rk_com(i2,24)*cnn(i2,koh)*cnn(i2,kno2) & ! oh +no2
                  + rk_com(i2,39)*cnn(i2,kno3)*cnn(i2,kno2) & ! no3+no2
                  + rk_com(i2,69)*cnn(i2,kc2o3)*cnn(i2,kno2)& ! c2o3+no2
                  + rk_urb(i2,20)*cnn(i2,kcro)*cnn(i2,kno2) &  ! cro +no2 
                  ) / cair_mlc(i2)*1.e+9


      if ( (po3-lo3) .ge. 1.e-20 .and. pnoz.ge.1.e-20) then
       ope(i2) = (po3 - lo3) / pnoz
      else
       ope(i2) = -1.e20
      endif
    enddo !i2

      return 
      endsubroutine
      
      end module mod_cbmz
