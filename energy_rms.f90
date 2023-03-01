program energy_rms    ! By Yalong Cong on 02/11/2020
implicit none
  
  integer*4 :: n
  integer*8 :: p,l
  integer*8 :: i,j,ist
  integer*8 :: natom
  integer*8 :: ntype
  integer*8 :: nstep
  integer*8 :: nsnapshot
  integer*8 :: iacj,ic
  integer*8 :: snapshot
  integer*8 :: snapshot_start
  integer*8 :: snapshot_end
  integer*8 :: snapshot_out
  integer*8 :: snapshot_step
  integer*8 :: natom_start_pro
  integer*8 :: natom_start_lig
  integer*8 :: natom_end_pro
  integer*8 :: natom_end_lig
  integer*8 :: date_time(8)                                !the time
  integer*8 :: ifbox                                       !the box in crdfile, 0:NO; 1:YES
  integer*8,allocatable :: iac(:),ico(:)
  real*8,allocatable :: box(:)
  real*8,allocatable :: ele(:)
  real*8,allocatable :: vdw(:)
  real*8,allocatable :: ave(:)
  real*8,allocatable :: rms(:)
  real*8,allocatable :: charge(:)
  real*8,allocatable :: coord(:,:)
  real*8,allocatable :: cn1(:),cn2(:)
  real*8 :: cut
  real*8 :: dielc
  real*8 :: dist,dist2,dist6,dist12
 
  character(len=160) :: infile
  character(len=160) :: outfile
  character(len=160) :: topfile
  character(len=160) :: crdfile
  character(len=160) :: energyfile
  character(len=160) :: totenergyfile
  character(len=160) :: fline
  character(len=160) :: date,time,zone,correntworkdirectory !the time and directory

  namelist /energy_rms_cntrl/ cut,dielc,ifbox,energyfile,totenergyfile,&
                              snapshot_start, snapshot_end,&
                              snapshot_step,  snapshot_out,&
                              natom_start_pro,natom_end_pro,&
                              natom_start_lig,natom_end_lig

  cut=999.0
  dielc=1.0
  ifbox=1
  natom_start_pro=1
  i=iargc()
  if(i.eq.0)then
    write(*,*)"energy_rmd  -h : for help"
    stop
  endif
  n=1
  do while(.true.)
    if(n.gt.i)exit
    call getarg(n,fline)
    if(fline.eq."-i")then
      call getarg(n+1,infile)
      n=n+2
    elseif(fline.eq."-o")then
      call getarg(n+1,outfile)
      n=n+2
    elseif(fline.eq."-p")then
      call getarg(n+1,topfile)
      n=n+2
    elseif(fline.eq."-c")then
      call getarg(n+1,crdfile)
      n=n+2
    elseif(fline.eq."-h")then
      write(*,*)"energy_rmd  -h : for help"
      write(*,*)"usage: energy_rms  -i infile -o outfile -p topfile -c crdfile"
      stop
    else
      write(*,*)'Error unknown flag: '//fline//''
      write(*,*)"energy_rmd  -h : for help"    
      stop
    endif

  enddo

  open(10,file=infile)
    read(10,nml=energy_rms_cntrl)
  close(10)
  allocate(ele(snapshot_start:snapshot_end))
  allocate(vdw(snapshot_start:snapshot_end))
  ele(snapshot_start:snapshot_end)=0.0
  vdw(snapshot_start:snapshot_end)=0.0

  open(20,file=outfile)
    call date_and_time(date,time,zone,date_time)
    call getcwd(correntworkdirectory)
    write(20,*)
    write(20,*)
    write(20,'(A41)')"#----------------------------------------"
    write(20,'(A41)')"#               ENERGY_RMS               "
    write(20,'(A41)')"#----------------------------------------"
    write(20,*)
    write(20,*)
    write(20,'(A41)')"#This program can pick out the trajectory"
    write(20,'(A41)')"#of convergence based on gas-phase energy"
    write(20,'(A15)')"#(ELE and VDW)."
    write(20,'(A29)')"#By Yalong Cong on 02/11/2020"
    write(20,*)
    write(20,*)
    write(20,'(A8,A2,A1,A2,A1,2A4,A2,A1,A2,A1,A2)')&
    "#Run on: ",date(7:8),"/",date(5:6),"/",date(1:4)," at ",time(1:2),":",time(3:4),":",time(5:6)
    write(20,'(A27)')"#Corrent Working Directory:"
    write(20,'(A1,T2,A80)')"#",correntworkdirectory
    write(20,*)
    write(20,*)
    write(20,'(A10)')"parameters"
    write(20,'(2(A10,F10.2))')"cut:",cut, "dielc:",dielc
 
    open(11,file=topfile)
      do while(.true.)
        read(11,'(A80)',iostat=ist)fline
        if(ist/=0)exit
        if(fline(1:14).eq.'%FLAG POINTERS')then
          read(11,*)
          read(11,'(2I8)')natom,ntype
          allocate(charge(natom))
          allocate(coord(natom,3))
          allocate(box(3))
          allocate(iac(natom))
          allocate(ico(ntype*ntype))
          allocate(cn1(-1:ntype*(ntype+1)/2))
          allocate(cn2(-1:ntype*(ntype+1)/2))
        endif
        if(fline(1:12).eq.'%FLAG CHARGE')then
          read(11,*)
          read(11,'(5E16.8)')(charge(i),i=1,natom)
        endif
        if(fline(2:21).eq.'FLAG ATOM_TYPE_INDEX')then
          read(11,*)
          read(11,'(10I8)')(iac(i),i=1,natom)
        endif
        if(fline(2:26).eq.'FLAG NONBONDED_PARM_INDEX')then
          read(11,*)
          read(11,'(10I8)')(ico(i),i=1,ntype*ntype)
        endif
        if(fline(2:25).eq.'FLAG LENNARD_JONES_ACOEF')then
          read(11,*)
          read(11,'(5E16.8)')(cn1(i),i=1,ntype*(ntype+1)/2)
        endif
        if(fline(2:25).eq.'FLAG LENNARD_JONES_BCOEF')then
          read(11,*)
          read(11,'(5E16.8)')(cn2(i),i=1,ntype*(ntype+1)/2)
        endif
      enddo
    close(11)
 
    open(12,file=crdfile)
    open(21,file=totenergyfile)
      write(21,'(4A16)')"Snapshot","ELE","VDW","ELE_VDW"
      snapshot=0
      read(12,*)
      do while(.true.)
        read(12,'(10F8.3)',iostat=ist)((coord(n,i),i=1,3),n=1,natom)
        if(ist/=0)exit
        if(ifbox.eq.1)then
          read(12,'(3F8.3)',iostat=ist)(box(i),i=1,3)
          if(ist/=0)exit
        endif
        snapshot=snapshot+1
        if(snapshot.ge.snapshot_start .and. snapshot.le.snapshot_end)then
          do p=natom_start_pro,natom_end_pro
            do l=natom_start_lig,natom_end_lig
              dist2=(coord(p,1)-coord(l,1))**2+(coord(p,2)-coord(l,2))**2+(coord(p,3)-coord(l,3))**2
              dist=sqrt(dist2)
              if(dist.le.cut)then
                dist6=dist2**3
                dist12=dist2**6
                ele(snapshot)=ele(snapshot)+charge(p)*charge(l)/(dist*dielc)
                iacj=ntype*(iac(p)-1)
                ic=ico(iacj+iac(l))
                vdw(snapshot)=vdw(snapshot)+(cn1(ic)/dist12-cn2(ic)/dist6)
              endif
            enddo
          enddo
          write(21,'(I16,3F16.4)')snapshot,ele(snapshot),vdw(snapshot),ele(snapshot)+vdw(snapshot)
        elseif(snapshot.gt.snapshot_end)then
          exit
        endif 
      enddo
    close(21)
    close(12)

    write(20,*)
    write(20,'(A48)')"================================================"
    write(20,'(3A16)')"NSTEP","AVERAGE","RMS"
    write(20,'(A48)')"------------------------------------------------"
    nsnapshot=snapshot_end-snapshot_start+1
    nstep=(nsnapshot-snapshot_out)/snapshot_step+1
    allocate(ave(nstep))
    allocate(rms(nstep))
    ave(nstep)=0.0
    rms(nstep)=0.0
    do i=1,nstep
      do j=snapshot_start+(i-1)*snapshot_step,snapshot_start+(i-1)*snapshot_step+snapshot_out-1
        ave(i)=ave(i)+(ele(j)+vdw(j))
      enddo
      ave(i)=ave(i)/snapshot_out
      do j=snapshot_start+(i-1)*snapshot_step,snapshot_start+(i-1)*snapshot_step+snapshot_out-1
        rms(i)=rms(i)+(ele(j)+vdw(j)-ave(i))**2
      enddo
      rms(i)=sqrt(rms(i)/snapshot_out)
      write(20,'(I16,2F16.4)')i,ave(i),rms(i)
    enddo
    write(20,'(A48)')"================================================"
    i=minloc(rms,1)
    write(20,*)
    write(20,'(A8,I8,2F16.4)')"MIX:",i,ave(i),rms(i)
    write(20,'(A8,I8,I16)')"Snap:",snapshot_start+(i-1)*snapshot_step,snapshot_start+(i-1)*snapshot_step+snapshot_out-1
    
    open(22,file=energyfile)
      write(22,'(4A16)')"Snapshot","ELE","VDW","ELE_VDW"
      do snapshot=snapshot_start+(i-1)*snapshot_step,snapshot_start+(i-1)*snapshot_step+snapshot_out-1
        write(22,'(I16,3F16.4)')snapshot,ele(snapshot),vdw(snapshot),ele(snapshot)+vdw(snapshot) 
      enddo
    close(22)

    write(20,*)
    write(20,'(A36,T37,A30)')"#Energy of corresponding snap is to ",energyfile
    call date_and_time(date,time,zone,date_time)
    write(20,'(A8,A2,A1,A2,A1,2A4,A2,A1,A2,A1,A2)')&
   "#End on: ",date(7:8),"/",date(5:6),"/",date(1:4)," at ",time(1:2),":",time(3:4),":",time(5:6)

  close(20)

end program energy_rms

