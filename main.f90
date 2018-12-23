!------------------------------------------------------
! naqpms        version 3.0
!
! lijie  from 12, 4, 2006 version 2 in cbm-z
! to fix the confusions and problems among different versions
!
! zifa   on oct.1,2006
! to add the source map factor based on the version in 2005/01
! the following revised
!     to track any number species/any sources
! lijie  add cbmz
! lijie  modify cbmz source mark for accurate calculating ozone ! production,gross,net and destruction
! lijie  add convection , vertical diffusion and direction, dry deposition for  gases
! lijie  set a input file in emit direction for automatical different emissions  files
! lijie  add a program global.tar.gz for set up global boundary conditions
! lijie  add isorropia for inorganic aerosols
! lijie  add soa formation schedume
! lijie  add a module for aerosol optical parameters
! lijie  modify aod
! lijie  add aqeous chemistry and new wet depostion
! lijie  add ope and visibility calculation
! lijie  remove tuv for fast calculation
! lijie  fix bugs for source marking when deltac or concentration is very small (close to 0)
! lijie  add dust module from luo geatm and modify it
! lijie  add seasalt emissions module
! chen huansheng change naqpms to global version
! chen huagsheng add hg module

 program main
#ifdef OMP_OPT
 use omp_lib
#endif
 use mod_cbmz
 use gas_data
 implicit none
 include 'mpif.h'

      integer :: lun_input,lun_output1,nchar
      real :: t1,t2,Elapsed_time
      integer :: comp_stime, comp_etime, clock_rate, clock_max

      !dir$ attributes align:64 :: cbmzobj
      type(cbmztype) cbmzobj
      !dir$ attributes align:64 :: cppb
      real ::  cppb(VLEN,ngas_max)
      !dir$ attributes align:64 :: fcld
      real,dimension(VLEN)   :: fcld          ! the coffi of cloud to photo

      character*40 inputfilename, outputfilename, dumword

      real, dimension(:,:), allocatable :: ggas
      real, dimension(:), allocatable :: grlat,grlon,gte,gRH,gpr_atm,gzalt_m
      integer :: i,i2,ig,ilen
      logical :: linit

      integer :: ierr, myid, numprocs
      integer :: mynpt,myst
      
  call mpi_init( ierr )
  call mpi_comm_rank( mpi_comm_world, myid, ierr )
  call mpi_comm_size( mpi_comm_world, numprocs, ierr )
  mynpt = (npt-1)/numprocs+1
  myst = myid*mynpt+1
  if (myid.eq.(numprocs-1)) then
     mynpt = npt - (numprocs-1)*mynpt
  endif
  !write(6,*),myid,mynpt,myst

      allocate(ggas(ngas_max,mynpt))
      allocate(grlat(mynpt),grlon(mynpt),gte(mynpt),&
               gRH(mynpt),gpr_atm(mynpt),gzalt_m(mynpt))
 

  if (myid.eq.0) then
      write(6,*)'   '
      write(6,*)'******************************************************'
      write(6,*)'             Optimized CBM-Z'
      write(6,*)'   '
      write(6,*)'******************************************************'
      write(6,*)'   '
      write(6,*)'   '
  endif


!!----- input ----------------------
      lun_input = 10
      !write(6,*)'Enter input filename. Example: *.input'
      !read(5,*)inputfilename
      inputfilename = 'cases.input'

  if (myid.eq.0) then
!!----- output ----------------------
      lun_output1 = 20
      nchar = nbllen(inputfilename) - 6
      outputfilename = inputfilename(1:nchar)//'.output'
      open(lun_output1, file = outputfilename)
  endif

      open(lun_input, file = inputfilename)
        call ReadInputFile(mynpt,myst,lun_input,ggas,grlat,grlon,gte,gRH,gpr_atm,gzalt_m)
      close(lun_input)

  call mpi_barrier( mpi_comm_world, ierr )
  if (myid.eq.0) then
      write(6,*)'Finished reading all inputs...'
      call system_clock(comp_stime,clock_rate,clock_max)
      call cpu_time(t1) ! time counting
  endif

      linit = .false.
      fcld(:) = 1.0
      emission(:,:) = 0.

#ifdef OMP_OPT
!$omp parallel default(none) &
!$omp shared(comp_stime,clock_rate,clock_max) &
!$omp shared(ggas,grlat,grlon,gte,gRH,gpr_atm,gzalt_m) &
!$omp shared(fcld,mynpt) &

!$omp private(i,i2,ig,ilen,linit) &
!$omp copyin(msolar,mphoto,dt_min) &
!$omp copyin(tbeg_dd,tbeg_hh,tbeg_mm,tbeg_ss,trun_dd,trun_hh,trun_mm) &
#ifdef KNL_OPT
!$omp private(cbmzobj) &
#endif
!$omp private(cppb)

!$omp do schedule(dynamic)
  do i=1,mynpt,VLEN
#else
gas_chemistry:   do i=1,mynpt,VLEN
#endif
      ilen = min(mynpt-i+1, VLEN)
      if (.not.linit) then
          cbmzobj%p_het(:,:) = 0.           !init only once
          cbmzobj%p_com(:,:) = 0.           !init only once
          cbmzobj%rl_com(:,:) = 0.           !init only once
          cbmzobj%r2_com(:,:) = 0.           !init only once
          cbmzobj%p_urb(:,:) = 0.           !init only once
          cbmzobj%rl_het(:,:) = 0.           !init only once
          cbmzobj%rl_urb(:,:) = 0.           !init only once
          cbmzobj%p_bio(:,:) = 0.           !init only once
          cbmzobj%rl_bio(:,:) = 0.           !init only once
          cbmzobj%p_mar(:,:) = 0.           !init only once
          cbmzobj%rl_mar(:,:) = 0.           !init only once
          linit = .true.
      endif
      do ig=1,ngas_max
        do i2 = 1, ilen
          cnn(i2,ig) = ggas(ig,i+i2-1)
        enddo !ig
      enddo !ig

      cbmzobj%pmask(:)  =  .false.             ! enabled flag
      do i2 = 1, ilen
        rlon(i2)   =  grlon(i+i2-1)       !rlat the box lat
        rlat(i2)   =  grlat(i+i2-1)       ! as the above ,but lon
        zalt_m(i2) =  gzalt_m(i+i2-1)     ! the altitude(asl) of box(m)
        cbmzobj%pmask(i2)  =  .true.             ! enabled flag
        cbmzobj%rh(i2)     =  gRH(i+i2-1)        ! the rh
        cbmzobj%te(i2)     =  gte(i+i2-1)        ! the temp
        cbmzobj%pr_atm(i2) =  gpr_atm(i+i2-1)  ! the pressure but as the atm
      enddo !i2

      call cbmz(0,cbmzobj,cppb,fcld)  !! need to modify

#ifdef OMP_OPT
  enddo !i
!$omp end do
!$omp end parallel
#else
  enddo gas_chemistry !i

#endif

  !write(6,*)'Finished.',myid
  !call mpi_barrier( mpi_comm_world, ierr )
  if (myid.eq.0) then
    call cpu_time(t2)

    call system_clock(comp_etime,clock_rate,clock_max)
  endif

  call mpi_finalize(ierr)
!!------------------------------------------------------------------
    deallocate(ggas)
    deallocate(grlat,grlon,gte,gRH,gpr_atm,gzalt_m)

  if (myid.eq.0) then
      write(6,*)'   '
      write(6,*)'   '
      write(6,*)'         End of Simulation'
      write(6,*)'   '
      write(6,*)'   '
      write(6,*)'******************************************************'

      Elapsed_time=t2-t1
      write(6,*)'Elapsed_time=',Elapsed_time
      print *, 'Elapsed wall time:', real(comp_etime-comp_stime)/real(clock_rate)
  endif
      stop
      end            

!!**********************************************************************
!!*********************  Subroutines ***********************************


      subroutine ReadInputFile(mynpt,myst,lin,ggas,grlat,grlon,gte,gRH,gpr_atm,gzalt_m)
      use gas_data

      integer :: mynpt,myst
      real, dimension(ngas_max,mynpt), intent(out) :: ggas
      real, dimension(mynpt), intent(out) :: grlat,grlon,gte,gRH,gpr_atm,gzalt_m
      character*40 dword
      integer :: ig,i,j,k, x,y,z
      real :: dummy

      read(lin,*)dword ! PARAMETERS

!!----------- begin time from 12:00 (noon) March 21 [min]
      read(lin,*)dword, tbeg_dd, tbeg_hh, tbeg_mm, tbeg_ss, dword
      read(lin,*)dword, trun_dd, trun_hh, trun_mm, trun_ss, dword
      read(lin,*)dword, dt_min,dword ! transport time-step [min]
      read(lin,*)dword, dummy,  dword ! longitude [deg]
      read(lin,*)dword, dummy,  dword ! latitude [deg]
      read(lin,*)dword, dummy,dword ! altitude  above mean sea level [m]
      read(lin,*)dword, dummy,    dword ! relative humidity [%]
      read(lin,*)dword, dummy,    dword ! temperature [K]
      read(lin,*)dword, dummy,dword ! pressure [atm]
      read(lin,*)dword, msolar,dword ! msolar flag
      read(lin,*)dword, mphoto,dword ! mphoto flag
      read(lin,*)dword, iprint,dword ! freq of output

!!----------- read lat,lon,te,RH and gas
      read(lin,*)dword ! GAS

      !!---- skip (myst-1) space points
      do i=1,myst-1
      do ig=1,9+67
        read(lin,*)
      enddo
      enddo  !skip
      do i=1,mynpt
      read(lin,*)dword, x ! X index
      read(lin,*)dword, y ! Y index
      read(lin,*)dword, z ! Z index
      !if (i.ne.x .or. j.ne.y .or. k.ne.z) then
      !  write(6,*)'inputs error at: ',i,j,k,x,y,z
      !endif
      !write(6,*)'X:',x,' Y:',y,' Z:',z
      read(lin,*)dword, grlat(i) ! LATITUDE
      !write(6,*)'rlat:',grlat(i)
      read(lin,*)dword, grlon(i) ! LONGITUDE
      !write(6,*)'rlon:',grlon(i)
      read(lin,*)gte(i) ! temperature [K]
      !write(6,*)'te:',gte(i)
      read(lin,*)gRH(i) ! relative humidity [%]
      !write(6,*)'RH:',gRH(i)
      read(lin,*)gpr_atm(i) ! pressure [atm]
      !write(6,*)'pr_atm:',gpr_atm(i)
      read(lin,*)gzalt_m(i) ! pressure [atm]
      !write(6,*)'zalt_m:',gzalt_m(i)
      do ig=1, 67!ngas_max
        read(lin,*)ggas(ig,i)
        !write(6,*)'gas(',ig,'):',ggas(ig,i)
      enddo
      ggas(68:ngas_max,i) = 0.
      enddo ! i

      !write(6,*)'Finished reading all inputs...'
      return
      end

