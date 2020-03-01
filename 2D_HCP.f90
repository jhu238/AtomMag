Program atomic_spin_dynamic
    implicit none  
    !!!!!!!!!!previous module parameter!!!!!!!!!!
    !Lattice
    integer Natom, Natom0      !Number of total atoms,number of remaining atoms
    integer Nx,Ny,Nz     !Number of grids at  x, y, z direction
    real(kind=8) ax,ay,az     !lattice constant at x, y, z direction
    real(kind=8),allocatable,dimension(:,:) :: atomposition,atompositionold !atom position
    integer,allocatable,dimension(:) :: oMag    !mark whether an atom has magnetization
    integer flagLattice  !1-FCC 2-SC 3-HCP 4-2D Hexagonal Latiice
    integer flagShape    !1-normal square   2-circular
    real(kind=8) center(3),radius(3)              !center and radius of the circular
    integer,allocatable,dimension(:,:) :: neighboratom     !number of neighboratom
    real(kind=8) distance   !distance between two atoms
    real(kind=8) dcriteria  !criteria for determining neighbors
    real(kind=8),allocatable,dimension(:,:)::distanceneighbor     !the minimum distance of neighbors
    
    !parameter
    real,parameter :: mu0=1.2566370614E-6
    real,parameter :: pi=3.141592653
    real,parameter :: kB=1.38E-23

    !Material
    real(kind=8) Ms        !saturation magnetization
    real(kind=8) J         !exchange costant     
    real(kind=8) D         !DMI constant
    integer flagDMI        !1-interfacial 2-bulk
    real(kind=8) Ku        !uniaxial anisoropy constant
    !real(kind=8) Kc        !cubic anisotropy constant
    real(kind=8) kaxis(3)  !anisotropy axis
    real(kind=8) B0(3)     !external field
    real(kind=8) gamma     !gyromagnetic ratio
    real(kind=8) alpha     !damping constant
    !real(kind=8) B_const   !B0 for unitless

    !thermal
    ! real(kind=8) temperature   !temperature for thermal effect
    ! real(kind=8) constThermal   !thermal_constant
    ! real(kind=8) randomnumber1(3),randomnumber2(3) !randomnumber for thermal field
    ! real(kind=8) randomnumber_gaussion(3)   !random number of gaussin distribution

    !SOT
    real(kind=8) betaFL,betaDL            !strength of field-like torque and damping-like torque
    real(kind=8) SPorientation(3)             !spin orientation of spin-orbit-torque
    real(kind=8) FLterm(3), DLterm(3)      !m cross product P, m cross product (m cross product P)
    
    !M,E&H
    real(kind=8),allocatable,dimension(:,:)::magnt     !spin polarization
    real(kind=8),allocatable,dimension(:,:)::Beff,Bexc,Bani,BDMI,Bext,Bstray!,Bthermal    !effective field
    real(kind=8),allocatable,dimension(:)::Eeff,Eexc,Eani,EDMI,Eext,Estray       !Hamiltonia
    integer flagspin          !1-up 2-random   3-for test
    
    !time
    integer stepBegin, stepEnd,istep,dstep
    real(kind=8) dt
    
    !solver for RK4
    integer flagSolver     !1-RK4 2-Heun with projection
    real(kind=8),allocatable,dimension(:,:,:)::dm
    real(kind=8),allocatable,dimension(:,:)::magnt_old 
    
    !for output
    real(kind=8) avmagnt(3)
    real(kind=8) avBeff(3),avBexc(3),avBani(3),avBDMI(3),avBext(3),avBstray(3)
    real(kind=8) avEeff,avEexc,avEani,avEDMI,avEext,avEstray
    integer flagoutputB    !0-not outputB 1-outputB
    integer filenumber
    character(len=8) filename

    logical:: file_exists

    integer ii,jj,kk,n0,nn
    integer i,k
    real(kind=8) product

    real(kind=8) randn(3)    !random number array
    real(kind=8) magnttotal   !for normalization
    real(kind=8) centerx,centery
    real(kind=8) rij(3), rij_total
    real(kind=8) product1, product2, product3

    !for input
    logical:: flagreadin  !flag for input: 0- no input 1-input
    integer Natomtemp     !number of atoms(input)
    real(kind=8),allocatable,dimension(:,:)::magnttemp

    !!!!main variables!!!!!!!!!!!!!!!!!!!!
    real:: start,finish
    integer nstep
    
    call cpu_time(start)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!Init!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!1.Initparameter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Lattice related
    Nx=100;Ny=1;Nz=1
    ax=2.715E-10;ay=2.715E-10;az=4.08E-10
    flagLattice=2    !1-FCC 2-SC 3-HCP 4-2D hexagonal
    flagShape=1      !1-normal square 2-circular
    
    !material related
    Ms=2.7822E-23
    J=9.408E-22
    D=2.496E-22
    flagDMI=1       !1-interfacial 2-bulk
    Ku=6.56e-23       
    !Kc=0.d0
    kaxis(1)=0.d0;   kaxis(2)=0.d0;       kaxis(3)=1.d0
    B0(1)=0.d0;      B0(2)=0.d0;          B0(3)=2.d0; 
    ! temperature=0.d0;
    betaFL=0.d0;
    betaDL=0.d0;
    SPorientation(1)=0.d0;  SPorientation(2)=0.d0;    SPorientation(3)=0.d0
    gamma=1.76E11
    alpha=0.017
    
    !initial spin
    flagspin=3
    !B_const=1.d0
    
    !time
    stepBegin=0; stepEnd=100; istep=10;dstep=10000
    dt=1E-20
    !eff_dt=gamma/(1+alpha**2)*B_const*dt
    
    !solver
    flagSolver=1

    !output
    flagoutputB=0    

    !input
    flagreadin=0
    Natomtemp=0
    
    !read parameter.in
    Inquire(File='parameter.in',EXIST=file_exists)
    if(file_exists) then
        open(unit=4,file='parameter.in')
        read(4,*),Nx,Ny,Nz
        read(4,*),ax,ay,az
        read(4,*),flagLattice
        read(4,*),flagShape
        read(4,*)
        read(4,*),Ms
        read(4,*),J
        read(4,*),D
        read(4,*),flagDMI
        read(4,*),Ku
        read(4,*),kaxis(1),kaxis(2),kaxis(3)
        read(4,*),B0(1),B0(2),B0(3)
        ! read(4,*),temperature
        read(4,*),betaFL
        read(4,*),betaDL
        read(4,*),SPorientation(1),SPorientation(2),SPorientation(3)
        read(4,*),gamma
        read(4,*),alpha
        read(4,*)
        read(4,*),flagspin
        read(4,*)
        read(4,*),stepBegin
        read(4,*),stepEnd
        read(4,*),istep
        read(4,*),dstep
        read(4,*),dt
        read(4,*)
        read(4,*),flagSolver
        read(4,*)
        read(4,*),flagoutputB
        read(4,*)
        read(4,*),flagreadin
        read(4,*),Natomtemp
        close(4)
    end if

    !thermal constant
    ! constThermal=sqrt(2*alpha*kB*temperature/(gamma*Ms*dt))
    call random_seed()

    !!!!!!!!2.InitLattice!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(flagreadin==0) then
        n0=1
        !1-FCC 2-SC 3-HCP 4-2D Hexagonal
        if(flagLattice==1) then
            Natom=Nx*Ny*Nz*4
            allocate(atomposition(Natom,3));  atomposition=0.d0
            do ii=1,Nx
            do jj=1,Ny
            do kk=1,Nz
                atomposition(n0,1)=(ii-1)*ax;             atomposition(n0,2)=(jj-1)*ay;                 atomposition(n0,3)=(kk-1)*az                  !atom 1 in unit cell
                atomposition(n0+1,1)=(ii-0.5)*ax;         atomposition(n0+1,2)=(jj-0.5)*ay;             atomposition(n0+1,3)=(kk-1)*az                !atom 2 in unit cell
                atomposition(n0+2,1)=(ii-0.5)*ax;         atomposition(n0+2,2)=(jj-1)*ay;               atomposition(n0+2,3)=(kk-0.5)*az            !atom 3 in unit cell
                atomposition(n0+3,1)=(ii-1)*ax;           atomposition(n0+3,2)=(jj-0.5)*ay;             atomposition(n0+3,3)=(kk-0.5)*az            !atom 4 in unit cell
                n0=n0+4
            enddo
            enddo
            enddo
            dcriteria=ax                     !for FCC, the nearest neighbor distance must be lower than ax and we could assume ax=ay=az, the nearest distance is sqrt(2)/2*a
        endif
    
        if(flagLattice==2) then
            Natom=Nx*Ny*Nz
            allocate(atomposition(Natom,3));  atomposition=0.d0
            do ii=1,Nx
            do jj=1,Ny
            do kk=1,Nz
                atomposition(n0,1)=(ii-1)*ax;       atomposition(n0,2)=(jj-1)*ay;               atomposition(n0,3)=(kk-1)*az;
                n0=n0+1
            enddo
            enddo
            enddo
        endif

        if(flagLattice==3) then
            Natom=Nx*Ny*Nz*3
            allocate(atomposition(Natom,3));  atomposition=0.d0
            do ii=1,Nx
            do jj=1,Ny
            do kk=1,Nz
                atomposition(n0,1)=(ii-1)*ax;              atomposition(n0,2)=(jj-1)*sqrt(3.0d0)*ay;                         atomposition(n0,3)=(kk-1)*az
                atomposition(n0+1,1)=(ii-1)*ax+0.5*ax;     atomposition(n0+1,2)=(jj-1)*sqrt(3.0d0)*ay+0.5d0*sqrt(3.0d0)*ay;      atomposition(n0,3)=(kk-1)*az
                atomposition(n0+2,1)=(ii-1)*ax+0.5*ax;     atomposition(n0+1,2)=(jj-1)*sqrt(3.0d0)*ay+(sqrt(3.0d0)/6.0d0)*ay;     atomposition(n0,3)=(kk-1)*az+0.5d0*az
                n0=n0+3
            enddo
            enddo
            enddo
            dcriteria=1.5*ax               !for HCP, ax=ay, az=1.63ax, the nearest distance is exatcly ax
        endif

        if(flagLattice==4) then
            Natom=Nx*Ny*Nz*2
            allocate(atomposition(Natom,3));  atomposition=0.d0
            do ii=1,Nx      
            do jj=1,Ny
            do kk=1,Nz
                atomposition(n0,1)=(ii-1)*ax;              atomposition(n0,2)=(jj-1)*sqrt(3.0d0)*ay;                         atomposition(n0,3)=(kk-1)*az
                atomposition(n0+1,1)=(ii-1)*ax+0.5*ax;     atomposition(n0+1,2)=(jj-1)*sqrt(3.0d0)*ay+0.5d0*sqrt(3.0d0)*ay;      atomposition(n0,3)=(kk-1)*az
                n0=n0+2
            enddo
            enddo
            enddo
            dcriteria=1.5*ax               !for 2D hexagonal, ax=ay, the nearest distance is exactly ax
        endif
    
        allocate(atompositionold(Natom,3));  atompositionold=0.d0
        allocate(oMag(Natom));               oMag=0
     
        !!!!!!!!!!!!3.InitShape!!!!!!!!!!!!!!!!!!!!!!!
        if(flagShape==1) then
            Natom0=0
            do n0=1,Natom
                oMag(n0)=1
                Natom0=Natom0+1
            enddo
        endif
        if(flagshape==2) then
            center=0.5d0*(maxval(atomposition,1)+minval(atomposition,1))
            radius=0.5d0*(maxval(atomposition,1)-minval(atomposition,1))
        
            Natom0=0
            do n0=1,Natom
                if(((atomposition(n0,1)-center(1))**2+(atomposition(n0,2)-center(2))**2)<=(radius(1)**2)) then
                    oMag(n0)=1
                    Natom0=Natom0+1
                endif
            enddo
        
            atompositionold=atomposition
            deallocate(atomposition)
            allocate(atomposition(Natom0,3));   atomposition=0.d0;
            i=1
            do n0=1,Natom
                if(oMag(n0)==1) then
                    atomposition(i,:)=atompositionold(n0,:)
                    i=i+1
                endif
            enddo
        endif
 
        deallocate(atompositionold,oMag)
        Natom=Natom0
    endif

    if(flagreadin==1) then
        Natom=Natomtemp
        allocate(atomposition(Natom,3))
        allocate(magnttemp(Natom,3))
        Inquire(File='magnt.in',EXIST=file_exists)
        if(file_exists) then
             open(unit=2,file='magnt.in')
             do n0=1,Natom
                read(2,*),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),magnttemp(n0,1),magnttemp(n0,2),magnttemp(n0,3)
             enddo
        endif
        Inquire(File='atomposition.in',EXIST=file_exists)
        if(file_exists) then
             open(unit=3,file='atomposition.in')
             do n0=1,Natom
                read(3,*),i,atomposition(n0,1),atomposition(n0,2),atomposition(n0,3)
             enddo
        endif
        dcriteria=1.5*ax
    endif

    open(unit=5,file='atomposition.dat')
    do n0=1,Natom
        write(5,'(i10,3ES16.7)') n0,atomposition(n0,1),atomposition(n0,2),atomposition(n0,3)
    enddo
    close(5)

    !!allocation of arrays
    
    !M,E&H
    allocate(magnt(Natom,3));           magnt=0.d0
    allocate(Beff(Natom,3));            Beff=0.d0
    allocate(Bexc(Natom,3));            Bexc=0.d0
    allocate(Bani(Natom,3));            Bani=0.d0
    allocate(BDMI(Natom,3));            BDMI=0.d0
    allocate(Bext(Natom,3));            Bext=0.d0
    allocate(Bstray(Natom,3));          Bstray=0.d0
    !allocate(Bthermal(Natom,3));        Bthermal=0.d0
    allocate(Eeff(Natom));              Eeff=0.d0
    allocate(Eexc(Natom));              Eexc=0.d0
    allocate(Eani(Natom));              Eani=0.d0
    allocate(EDMI(Natom));              EDMI=0.d0
    allocate(Eext(Natom));              Eext=0.d0
    allocate(Estray(Natom));            Estray=0.d0
    
    !solver
    allocate(dm(Natom,3,4));       dm=0.d0
    allocate(magnt_old(Natom,3)); magnt_old=0.d0
    
    !output
    filenumber=17
    
    
    !!!!!!!!3.InitM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(flagreadin==0) then
        !flagspin:1-up 2-random
        if(flagspin==1) then
            do n0=1,Natom
                magnt(n0,1)=0.d0;            magnt(n0,2)=0.d0;            magnt(n0,3)=1.d0
            enddo
        endif
    
        if(flagspin==2) then   !call random seed
            do n0=1,Natom
                call RANDOM_NUMBER(randn);
                magnt(n0,1)=2.d0*randn(1)-1.d0
                magnt(n0,2)=2.d0*randn(2)-1.d0
                magnt(n0,3)=2.d0*randn(3)-1.d0
                magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                magnt(n0,1)=magnt(n0,1)/magnttotal
                magnt(n0,2)=magnt(n0,2)/magnttotal
                magnt(n0,3)=magnt(n0,3)/magnttotal
            enddo
        endif
    
        if(flagspin==3) then
            do n0=1,Natom
                if(n0<44) then
                    magnt(n0,1)=0.d0;       magnt(n0,2)=0.1;           magnt(n0,3)=0.9
                    magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                    magnt(n0,1)=magnt(n0,1)/magnttotal
                    magnt(n0,2)=magnt(n0,2)/magnttotal
                    magnt(n0,3)=magnt(n0,3)/magnttotal
                else if((n0<=57).and.(n0>=44)) then
                    magnt(n0,1)=0.d0;       magnt(n0,2)=0.1;           magnt(n0,3)=-0.9
                    magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                    magnt(n0,1)=magnt(n0,1)/magnttotal
                    magnt(n0,2)=magnt(n0,2)/magnttotal
                    magnt(n0,3)=magnt(n0,3)/magnttotal
                else
                    magnt(n0,1)=0.d0;       magnt(n0,2)=0.1;           magnt(n0,3)=0.9
                    magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                    magnt(n0,1)=magnt(n0,1)/magnttotal
                    magnt(n0,2)=magnt(n0,2)/magnttotal
                    magnt(n0,3)=magnt(n0,3)/magnttotal
                endif
            enddo
        endif
    
        if(flagspin==4) then
            centerx=0.5d0*(atomposition(1,1)+atomposition((Nx-1)*Ny*Nz+1,1))
            centery=0.5d0*(atomposition(1,2)+atomposition((Ny-1)*Nz+1,2))
            do n0=1,Natom
                if((atomposition(n0,1)-centerx)**2+(atomposition(n0,2)-centery)**2<(1E-9)**2) then
                    magnt(n0,1)=0.d0;       magnt(n0,2)=0.d0;           magnt(n0,3)=-1.0
                else
                    magnt(n0,1)=0.d0;       magnt(n0,2)=0.d0;           magnt(n0,3)=1.0
                endif
            enddo
        endif
    endif
        
    if(flagreadin==1) then
        magnt=magnttemp
    endif
    !!!!!!!!4.InitNeighbors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!FCC or HCP
    if((flagLattice==1).or.(flagLattice==3)) then
        allocate(neighboratom(Natom,12));            neighboratom=0.d0       !12 nearest neighbors in FCC
        allocate(distanceneighbor(Natom,12));        distanceneighbor=0.d0
        do n0=1,Natom
            distanceneighbor(n0,1:12)=dcriteria      !given an initial large value of distance
            neighboratom(n0,1:12)=0                           !given an initial value of neighboring atoms
            do nn=1,Natom
                if(nn/=n0) then
                    distance=sqrt((atomposition(nn,1)-atomposition(n0,1))**2+(atomposition(nn,2)-atomposition(n0,2))**2+(atomposition(nn,3)-atomposition(n0,3))**2)
                    if(distance<distanceneighbor(n0,1)) then
                        do i=12,2,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,1)=distance
                        neighboratom(n0,1)=nn
                    elseif(distance<distanceneighbor(n0,2)) then
                        do i=12,3,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,2)=distance
                        neighboratom(n0,2)=nn
                    elseif(distance<distanceneighbor(n0,3)) then
                        do i=12,4,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,3)=distance
                        neighboratom(n0,3)=nn
                    elseif(distance<distanceneighbor(n0,4)) then
                        do i=12,5,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,4)=distance
                        neighboratom(n0,4)=nn
                    elseif(distance<distanceneighbor(n0,5)) then
                        do i=12,6,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,5)=distance
                        neighboratom(n0,5)=nn
                    elseif(distance<distanceneighbor(n0,6)) then
                        do i=12,7,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,6)=distance
                        neighboratom(n0,6)=nn
                    elseif(distance<distanceneighbor(n0,7)) then
                        do i=12,8,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,7)=distance
                        neighboratom(n0,7)=nn
                    elseif(distance<distanceneighbor(n0,8)) then
                        do i=12,9,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,8)=distance
                        neighboratom(n0,8)=nn
                    elseif(distance<distanceneighbor(n0,9)) then
                        do i=12,10,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,9)=distance
                        neighboratom(n0,9)=nn
                    elseif(distance<distanceneighbor(n0,10)) then
                        do i=12,11,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                        enddo
                        distanceneighbor(n0,10)=distance
                        neighboratom(n0,10)=nn
                    elseif(distance<distanceneighbor(n0,11)) then
                        distanceneighbor(n0,12)=distanceneighbor(n0,11)
                        neighboratom(n0,12)=neighboratom(n0,11)
                        distanceneighbor(n0,11)=distance
                        neighboratom(n0,11)=nn
                    elseif(distance<distanceneighbor(n0,12)) then
                        distanceneighbor(n0,12)=distance
                        neighboratom(n0,12)=nn
                    endif
                endif
            enddo
        enddo
    endif

    !SC
    if(flagLattice==2) then
        allocate(neighboratom(Natom,6));       neighboratom=0.d0       !6 nearest neighbors in SC
        allocate(distanceneighbor(Natom,6));   distanceneighbor=0.d0
        do n0=1,Natom
            if((n0-1)>=1) then
                neighboratom(n0,1)=n0-1        !z-
            endif
            if((n0+1)<=Natom) then
                neighboratom(n0,2)=n0+1        !z+
            endif
            !y neighbor
            if((n0-Nz)>=1) then
                neighboratom(n0,3)=n0-Nz        !y- 
            endif
            if((n0+Nz)<=Natom) then
                neighboratom(n0,4)=n0+Nz        !y+
            endif
            !x neighbor
            if((n0-Nz*Ny)>=1) then
                neighboratom(n0,5)=n0-Nz*Ny    !x-
            endif
            if((n0+Nz*Ny)<=Natom) then
                neighboratom(n0,6)=n0+Nz*Ny    !x+
            endif
        enddo
    endif

    !2D Hexagonal
    if(flagLattice==4) then
        allocate(neighboratom(Natom,6));       neighboratom=0.d0       !6 nearest neighbors in SC
        allocate(distanceneighbor(Natom,6));   distanceneighbor=0.d0
        do n0=1,Natom
           distanceneighbor(n0,1:6)=dcriteria      !given an initial large value of distance
           neighboratom(n0,1:6)=0                           !given an initial value of neighboring atoms
           do nn=1,Natom
               if(nn/=n0) then
                   distance=sqrt((atomposition(nn,1)-atomposition(n0,1))**2+(atomposition(nn,2)-atomposition(n0,2))**2+(atomposition(nn,3)-atomposition(n0,3))**2)
                   if(distance<distanceneighbor(n0,1)) then
                       do i=6,2,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                       enddo
                       distanceneighbor(n0,1)=distance
                       neighboratom(n0,1)=nn
                   else if(distance<distanceneighbor(n0,2)) then
                       do i=6,3,-1
                          distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                          neighboratom(n0,i)=neighboratom(n0,i-1)
                       enddo
                       distanceneighbor(n0,2)=distance
                       neighboratom(n0,2)=nn
                   else if(distance<distanceneighbor(n0,3)) then
                       do i=6,4,-1
                          distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                          neighboratom(n0,i)=neighboratom(n0,i-1) 
                       enddo
                       distanceneighbor(n0,3)=distance             
                       neighboratom(n0,3)=nn
                   else if(distance<distanceneighbor(n0,4)) then
                       do i=6,5,-1
                           distanceneighbor(n0,i)=distanceneighbor(n0,i-1)
                           neighboratom(n0,i)=neighboratom(n0,i-1)
                       enddo
                       distanceneighbor(n0,4)=distance
                       neighboratom(n0,4)=nn
                   else if(distance<distanceneighbor(n0,5)) then
                       distanceneighbor(n0,6)=distanceneighbor(n0,5)
                       neighboratom(n0,6)=neighboratom(n0,5)
                       distanceneighbor(n0,5)=distance
                       neighboratom(n0,5)=nn
                   else if(distance<distanceneighbor(n0,6)) then
                       distanceneighbor(n0,6)=distance
                       neighboratom(n0,6)=nn
                   end if
               endif
           enddo
        enddo
    endif
    
    !check the neighbor situation of 2D Hexagonal
    open(unit=7,file='neighbor.dat')
    do n0=1,Natom
        write(7,'(7i10)') n0,neighboratom(n0,1),neighboratom(n0,2),neighboratom(n0,3),neighboratom(n0,4),neighboratom(n0,5),neighboratom(n0,6)
    enddo
    close(7)
        
    !!!!!!!!5.calculateBinit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!calculateBexcint!!!!!!!!!!!!!
    do n0=1,Natom
        Bexc(n0,1)=0;    Bexc(n0,2)=0;    Bexc(n0,3)=0
        if((flagLattice==1).or.(flagLattice==3)) then
            do nn=1,12
                if(neighboratom(n0,nn)/=0) then
                    Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                    Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                    Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                endif
            enddo
        endif
        if(flagLattice==2) then      !only for SC
            if(Nz/=1) then
                do nn=1,2
                    if(neighboratom(n0,nn)/=0) then
                        Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                        Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                        Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                    endif
                enddo
            endif
            if(Ny/=1) then
                do nn=3,4
                    if(neighboratom(n0,nn)/=0) then
                        Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                        Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                        Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                    endif
                enddo
            endif
            if(Nx/=1) then
                do nn=5,6
                    if(neighboratom(n0,nn)/=0) then
                        Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                        Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                        Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                    endif
                enddo
            endif
        endif
        if(flagLattice==4) then
            do nn=1,6
                if(neighboratom(n0,nn)/=0) then
                    Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                    Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                    Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                endif
            enddo
        endif
        Eexc(n0)=-0.5d0*Ms*(magnt(n0,1)*Bexc(n0,1)+magnt(n0,2)*Bexc(n0,2)+magnt(n0,3)*Bexc(n0,3))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!calculateBaniint!!!!!!!!!!!!!!!!
    do n0=1,Natom
        product=0
        do i=1,3
            product=product+magnt(n0,i)*kaxis(i)
        enddo
        Bani(n0,1)=(2*Ku/Ms)*product*kaxis(1)
        Bani(n0,2)=(2*Ku/Ms)*product*kaxis(2)
        Bani(n0,3)=(2*Ku/Ms)*product*kaxis(3)
        Eani(n0)=-0.5d0*Ms*(magnt(n0,1)*Bani(n0,1)+magnt(n0,2)*Bani(n0,2)+magnt(n0,3)*Bani(n0,3))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!calculateBDMIint!!!!!!!!!!!!!!!!!!!!
    do n0=1,Natom
        BDMI(n0,1)=0;      BDMI(n0,2)=0;         BDMI(n0,3)=0
        if(flagDMI==1) then
            if((flagLattice==1).or.(flagLattice==3)) then
                do nn=1,12
                    if(neighboratom(n0,nn)/=0) then
                        rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                        rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                        rij(3)=0                        
                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                        rij(1)=rij(1)/rij_total
                        rij(2)=rij(2)/rij_total
                        rij(3)=rij(3)/rij_total
                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))                               
                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                    endif
                enddo
            endif
            if(flagLattice==2) then
                if(Nz/=1) then
                    do nn=1,2
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                            rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                            rij(3)=0                        
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
                endif
                if(Ny/=1) then
                    do nn=3,4
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                            rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                            rij(3)=0  
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
               endif
               if(Nx/=1) then
                    do nn=5,6
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                            rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                            rij(3)=0 
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
               endif
            endif
            if(flagLattice==4) then
                do nn=1,6
                    if(neighboratom(n0,nn)/=0) then
                        rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                        rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                        rij(3)=0                       
                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                        rij(1)=rij(1)/rij_total
                        rij(2)=rij(2)/rij_total
                        rij(3)=rij(3)/rij_total
                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                    endif
                enddo
            endif
        endif
        if(flagDMI==2) then
            if((flagLattice==1).or.(flagLattice==3)) then
                do nn=1,12
                    if(neighboratom(n0,nn)/=0) then
                        rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                        rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                        rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                        rij(1)=rij(1)/rij_total
                        rij(2)=rij(2)/rij_total
                        rij(3)=rij(3)/rij_total
                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                    endif
                enddo
            endif
            if(flagLattice==2) then
                if(Nz/=1) then
                    do nn=1,2
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                            rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                            rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
                endif
                if(Ny/=1) then
                    do nn=3,4
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                            rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                            rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
                endif
                if(Nx/=1) then
                    do nn=5,6
                        if(neighboratom(n0,nn)/=0) then
                            rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                            rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                            rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                            BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                            BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                        endif
                    enddo
                endif
            endif
            if(flagLattice==4) then
                do nn=1,6
                    if(neighboratom(n0,nn)/=0) then
                        rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                        rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                        rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                        rij(1)=rij(1)/rij_total
                        rij(2)=rij(2)/rij_total
                        rij(3)=rij(3)/rij_total
                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                    endif
                enddo
            endif
        endif
        EDMI(n0)=-0.5d0*Ms*(magnt(n0,1)*BDMI(n0,1)+magnt(n0,2)*BDMI(n0,2)+magnt(n0,3)*BDMI(n0,3))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateBextint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n0=1,Natom
        Bext(n0,1)=B0(1)
        Bext(n0,2)=B0(2)
        Bext(n0,3)=B0(3)
        Eext(n0)=-Ms*(magnt(n0,1)*Bext(n0,1)+magnt(n0,2)*Bext(n0,2)+magnt(n0,3)*Bext(n0,3))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateBstrayint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n0=1,Natom
        Bstray(n0,1)=0;      Bstray(n0,2)=0;         Bstray(n0,3)=0
        do i=1,Natom
            if(i/=n0) then
                rij(1)=atomposition(i,1)-atomposition(n0,1)
                rij(2)=atomposition(i,2)-atomposition(n0,2)
                rij(3)=atomposition(i,3)-atomposition(n0,3)
                rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                rij(1)=rij(1)/rij_total
                rij(2)=rij(2)/rij_total
                rij(3)=rij(3)/rij_total
                product=magnt(i,1)*rij(1)+magnt(i,2)*rij(2)+magnt(i,3)*rij(3)
                Bstray(n0,1)=Bstray(n0,1)+(mu0*Ms/(4*pi))*(3*product*rij(1)-magnt(i,1))/rij_total**3
                Bstray(n0,2)=Bstray(n0,2)+(mu0*Ms/(4*pi))*(3*product*rij(2)-magnt(i,2))/rij_total**3
                Bstray(n0,3)=Bstray(n0,3)+(mu0*Ms/(4*pi))*(3*product*rij(3)-magnt(i,3))/rij_total**3
            endif
        enddo
        Estray(n0)=-0.5d0*Ms*(magnt(n0,1)*Bstray(n0,1)+magnt(n0,2)*Bstray(n0,2)+magnt(n0,3)*Bstray(n0,3))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateBthermalint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do n0=1,Natom
        ! call random_number(randomnumber1)
        ! call random_number(randomnumber2)
        ! randomnumber_gaussion(1)=sqrt(-2*log(randomnumber1(1))*cos(2*pi*randomnumber2(1)))
        ! randomnumber_gaussion(2)=sqrt(-2*log(randomnumber1(2))*cos(2*pi*randomnumber2(2)))
        ! randomnumber_gaussion(3)=sqrt(-2*log(randomnumber1(3))*cos(2*pi*randomnumber2(3)))
        ! Bthermal(n0,1)=constThermal*randomnumber_gaussion(1)
        ! Bthermal(n0,2)=constThermal*randomnumber_gaussion(2)
        ! Bthermal(n0,3)=constThermal*randomnumber_gaussion(3)
    ! enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateBint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n0=1,Natom
        Beff(n0,1)=Bexc(n0,1)+Bani(n0,1)+BDMI(n0,1)+Bext(n0,1)+Bstray(n0,1)!+Bthermal(n0,1)
        Beff(n0,2)=Bexc(n0,2)+Bani(n0,2)+BDMI(n0,2)+Bext(n0,2)+Bstray(n0,2)!+Bthermal(n0,2)
        Beff(n0,3)=Bexc(n0,3)+Bani(n0,3)+BDMI(n0,3)+Bext(n0,3)+Bstray(n0,3)!+Bthermal(n0,3)
        Eeff(n0)=Eexc(n0)+Eani(n0)+EDMI(n0)+Eext(n0)+Estray(n0)
    enddo
    
    write(*,*),"Finish initial"
    
    open(unit=8,file='avMagnt.dat')
    open(unit=16,file='avEnergy.dat')
   
    if(flagSolver==1) then
    !remember to add Bthermal when reopen thermal field
    !$acc data copy(magnt,Beff,Bexc,Bani,BDMI,Bext,Bstray,Eeff,Eexc,Eani,EDMI,Eext,Estray,dm,magnt_old) copyin(kaxis,neighboratom,atomposition,b0,SPorientation,FLterm,DLterm)  
        do nstep=stepBegin, stepEnd
            if((nstep==stepBegin).or.(mod(nstep,istep)==0)) then
            !$acc update host(magnt,Eeff,Eexc,Eani,EDMI,Eext,Estray)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!outputM!!!!!!!!!!!!!!!!!!!!!!!
                avmagnt(:)=0    
                do n0=1,Natom
                    avmagnt(1)=avmagnt(1)+magnt(n0,1)
                    avmagnt(2)=avmagnt(2)+magnt(n0,2)
                    avmagnt(3)=avmagnt(3)+magnt(n0,3)
                enddo
                avmagnt(:)=avmagnt(:)/Natom
                write(8,'(i10,3ES16.7)'),nstep, avmagnt(1), avmagnt(2), avmagnt(3)
                !write(*,*),nstep,avmagnt(1),avmagnt(2),avmagnt(3)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!outputE!!!!!!!!!!!!!!!!!!!!!!!!
                avEeff=0;      avEexc=0;        avEani=0;           avEDMI=0;             avEext=0;             avEstray=0;         
                do n0=1,Natom
                    avEeff=avEeff+Eeff(n0)
                    avEexc=avEexc+Eexc(n0)
                    avEani=avEani+Eani(n0)
                    avEDMI=avEDMI+EDMI(n0)
                    avEext=avEext+Eext(n0)
                    avEstray=avEstray+Estray(n0)
                enddo
                
                avEeff=avEeff/Natom
                avEexc=avEexc/Natom
                avEani=avEani/Natom
                avEDMI=avEDMI/Natom
                avEext=avEext/Natom
                avEstray=avEstray/Natom
                
                write(16,'(i10,6ES16.7)'),nstep,avEeff,avEexc,avEani,avEDMI,avEext,avEstray
                !write(*,'(i10,6ES16.7)'),nstep,avEeff,avEexc,avEani,avEDMI,avEext,avEstray
                
                !write(*,*),nstep,avmagnt(1),avmagnt(2),avmagnt(3),avEeff,avEexc,avEani,avEDMI,avEExt,avEstray
            endif
            if((nstep==stepBegin).or.(mod(nstep,dstep)==0)) then
           !$acc update host(magnt)
           !!!!!!!!!!!!!!!!!!!!!!!!!outputdistriM!!!!!!!!!!!!!!!!!!!!
                write(filename,"(I8.8)") nstep
                open(unit=filenumber,file='magnt.'//filename//'.dat')
                do n0=1,Natom
                    write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),magnt(n0,1),magnt(n0,2),magnt(n0,3)
                enddo
                close(filenumber)
                
                filenumber=filenumber+1
            endif
            if(flagoutputB==1) then
                if((nstep==stepBegin).or.(mod(nstep,dstep)==0)) then
                    !remember to add Bthermal when reopen thermal field
                    !$acc update host(Beff,Bexc,Bani,BDMI,Bext,Bstray)
                    !!!!!!!!!!!!!!!!!!!!!!!!!outputdistriBexc!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='Bexc.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Bexc(n0,1),Bexc(n0,2),Bexc(n0,3)
                    enddo
                    close(filenumber)
                    
                    filenumber=filenumber+1
                    !!!!!!!!!!!!!!!!!!!!!!!!outputBani!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='Bani.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Bani(n0,1),Bani(n0,2),Bani(n0,3)
                    enddo
                    close(filenumber)
                    
                    filenumber=filenumber+1
                    !!!!!!!!!!!!!!!!!!!!!outputBDMI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='BDMI.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),BDMI(n0,1),BDMI(n0,2),BDMI(n0,3)
                    enddo
                    close(filenumber)
                    
                    filenumber=filenumber+1
                    !!!!!!!!!!!!!!!!!!!!!!outputBext!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='Bext.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Bext(n0,1),Bext(n0,2),Bext(n0,3)
                    enddo
                    close(filenumber)
                    
                    filenumber=filenumber+1
                    !!!!!!!!!!!!!!!!!!!!!!outputBstray!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='Bstray.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Bstray(n0,1),Bstray(n0,2),Bstray(n0,3)
                    enddo
                    close(filenumber)

                    filenumber=filenumber+1
                    !!!!!!!!!!!!!!!!!!!!!!outputBthermal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! write(filename,"(I8.8)") nstep
                    ! open(unit=filenumber,file='Bthermal.'//filename//'.dat')
                    ! do n0=1,Natom
                        ! write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Bthermal(n0,1),Bthermal(n0,2),Bthermal(n0,3)
                    ! enddo
                    ! close(filenumber)

                    !!!!!!!!!!!!!!!!!!!!!!outputBeff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    write(filename,"(I8.8)") nstep
                    open(unit=filenumber,file='Beff.'//filename//'.dat')
                    do n0=1,Natom
                        write(filenumber,'(3ES9.2,3ES16.7)'),atomposition(n0,1),atomposition(n0,2),atomposition(n0,3),Beff(n0,1),Beff(n0,2),Beff(n0,3)
                    enddo
                    close(filenumber)
                    
                    filenumber=filenumber+1

                endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!RK4solve!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc kernels
            do n0=1,Natom
                magnt_old(n0,1)=magnt(n0,1)
                magnt_old(n0,2)=magnt(n0,2)
                magnt_old(n0,3)=magnt(n0,3)
            enddo
            !$acc end kernels
            
            !RK4: 1st
            !$acc parallel loop private(product1, product2, product3,FLterm,DLterm)
            do n0=1,Natom
                product1=magnt(n0,2)*Beff(n0,3)-magnt(n0,3)*Beff(n0,2)
                product2=magnt(n0,3)*Beff(n0,1)-magnt(n0,1)*Beff(n0,3)
                product3=magnt(n0,1)*Beff(n0,2)-magnt(n0,2)*Beff(n0,1)
        
                FLterm(1)=magnt(n0,2)*SPorientation(3)-magnt(n0,3)*SPorientation(2)
                FLterm(2)=magnt(n0,3)*SPorientation(1)-magnt(n0,1)*SPorientation(3)
                FLterm(3)=magnt(n0,1)*SPorientation(2)-magnt(n0,2)*SPorientation(1)
                DLterm(1)=magnt(n0,2)*FLterm(3)-magnt(n0,3)*FLterm(2)
                DLterm(2)=magnt(n0,3)*FLterm(1)-magnt(n0,1)*FLterm(3)
                DLterm(3)=magnt(n0,1)*FLterm(2)-magnt(n0,2)*FLterm(1)

                dm(n0,1,1)=(dt*gamma/(1+alpha**2))*(-product1-alpha*(magnt(n0,2)*product3-magnt(n0,3)*product2))
                dm(n0,2,1)=(dt*gamma/(1+alpha**2))*(-product2-alpha*(magnt(n0,3)*product1-magnt(n0,1)*product3))
                dm(n0,3,1)=(dt*gamma/(1+alpha**2))*(-product3-alpha*(magnt(n0,1)*product2-magnt(n0,2)*product1))

                dm(n0,1,1)=dm(n0,1,1)+dt*betaFL*FLterm(1)+dt*betaDL*DLterm(1)
                dm(n0,2,1)=dm(n0,2,1)+dt*betaFL*FLterm(2)+dt*betaDL*DLterm(2)
                dm(n0,3,1)=dm(n0,3,1)+dt*betaFL*FLterm(3)+dt*betaDL*DLterm(3)
            enddo
            
            do i=2,4
                !$acc parallel loop private(magnttotal) 
                do n0=1,Natom
                    if(i==4) then
                        magnt(n0,1)=magnt_old(n0,1)+dm(n0,1,i-1)
                        magnt(n0,2)=magnt_old(n0,2)+dm(n0,2,i-1)
                        magnt(n0,3)=magnt_old(n0,3)+dm(n0,3,i-1)
                    else
                        magnt(n0,1)=magnt_old(n0,1)+0.5d0*dm(n0,1,i-1)
                        magnt(n0,2)=magnt_old(n0,2)+0.5d0*dm(n0,2,i-1)
                        magnt(n0,3)=magnt_old(n0,3)+0.5d0*dm(n0,3,i-1)
                    endif
                    magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                    magnt(n0,1)=magnt(n0,1)/magnttotal
                    magnt(n0,2)=magnt(n0,2)/magnttotal
                    magnt(n0,3)=magnt(n0,3)/magnttotal
                enddo
                !!!!!!!!!!!!!!!!!!!!!!!calculateB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!1.calculateBexc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc parallel loop private(nn)
                do n0=1,Natom
                    Bexc(n0,1)=0;    Bexc(n0,2)=0;    Bexc(n0,3)=0;
                    if((flagLattice==1).or.(flagLattice==3)) then
                        do nn=1,12
                            if(neighboratom(n0,nn)/=0) then
                                Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                            endif
                        enddo
                    endif
                    if(flagLattice==2) then      !only for SC
                        if(Nz/=1) then
                            do nn=1,2
                                if(neighboratom(n0,nn)/=0) then
                                    Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                    Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                    Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                                endif
                            enddo
                        endif
                        if(Ny/=1) then
                            do nn=3,4
                                if(neighboratom(n0,nn)/=0) then
                                    Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                    Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                    Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                                endif
                            enddo
                        endif
                        if(Nx/=1) then
                            do nn=5,6
                                if(neighboratom(n0,nn)/=0) then
                                    Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                    Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                    Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                                endif
                            enddo
                        endif
                    endif
                    if(flagLattice==4) then
                        do nn=1,6
                            if(neighboratom(n0,nn)/=0) then
                                Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                            endif
                        enddo
                    endif
                    Eexc(n0)=-0.5d0*Ms*(magnt(n0,1)*Bexc(n0,1)+magnt(n0,2)*Bexc(n0,2)+magnt(n0,3)*Bexc(n0,3))
                enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!2.calculateBani!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc kernels
                do n0=1,Natom
                    product=0
                    do k=1,3
                        product=product+magnt(n0,k)*kaxis(k)
                    enddo
                    Bani(n0,1)=(2*Ku/Ms)*product*kaxis(1)
                    Bani(n0,2)=(2*Ku/Ms)*product*kaxis(2)
                    Bani(n0,3)=(2*Ku/Ms)*product*kaxis(3)
                    Eani(n0)=-0.5d0*Ms*(magnt(n0,1)*Bani(n0,1)+magnt(n0,2)*Bani(n0,2)+magnt(n0,3)*Bani(n0,3))
                enddo
                !$acc end kernels
                !!!!!!!!!!!!!!!!!!!!!!!!!!3.calculateBDMI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc parallel loop private(rij,rij_total,nn)
                do n0=1,Natom
                    BDMI(n0,1)=0;      BDMI(n0,2)=0;         BDMI(n0,3)=0
                    if(flagDMI==1) then
                        if((flagLattice==1).or.(flagLattice==3)) then
                            do nn=1,12
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                    rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                    rij(3)=0                       
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(flagLattice==2) then
                            if(Nz/=1) then
                                do nn=1,2
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                        rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                        rij(3)=0                      
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                            if(Ny/=1) then
                                do nn=3,4
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                        rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                        rij(3)=0 
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                            if(Nx/=1) then
                                do nn=5,6
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                        rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                        rij(3)=0  
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                        endif
                        if(flagLattice==4) then
                            do nn=1,6
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                    rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                    rij(3)=0                        
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                    endif
                    if(flagDMI==2) then
                        if((flagLattice==1).or.(flagLattice==3)) then
                            do nn=1,12
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                    rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                    rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(flagLattice==2) then
                            if(Nz/=1) then
                                do nn=1,2
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                        rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                        rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                            if(Ny/=1) then
                                do nn=3,4
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                        rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                        rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                            if(Nx/=1) then
                                do nn=5,6
                                    if(neighboratom(n0,nn)/=0) then
                                        rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                        rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                        rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                        rij(1)=rij(1)/rij_total
                                        rij(2)=rij(2)/rij_total
                                        rij(3)=rij(3)/rij_total     
                                        BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                        BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                        BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                    endif
                                enddo
                            endif
                        endif
                        if(flagLattice==4) then
                            do nn=1,6
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                    rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                    rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                    endif
                    EDMI(n0)=-0.5d0*Ms*(magnt(n0,1)*BDMI(n0,1)+magnt(n0,2)*BDMI(n0,2)+magnt(n0,3)*BDMI(n0,3))
                enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4.calculateBext!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc kernels
                do n0=1,Natom
                    Bext(n0,1)=B0(1)
                    Bext(n0,2)=B0(2)
                    Bext(n0,3)=B0(3)
                    Eext(n0)=-Ms*(magnt(n0,1)*Bext(n0,1)+magnt(n0,2)*Bext(n0,2)+magnt(n0,3)*Bext(n0,3))
                enddo
                !$acc end kernels
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5.calculateBstray!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc parallel loop private(rij,rij_total,product,i)
                do n0=1,Natom
                    Bstray(n0,1)=0;      Bstray(n0,2)=0;         Bstray(n0,3)=0
                    do k=1,Natom
                        if(k/=n0) then
                            rij(1)=atomposition(k,1)-atomposition(n0,1)
                            rij(2)=atomposition(k,2)-atomposition(n0,2)
                            rij(3)=atomposition(k,3)-atomposition(n0,3)
                            rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                            rij(1)=rij(1)/rij_total
                            rij(2)=rij(2)/rij_total
                            rij(3)=rij(3)/rij_total
                            product=magnt(k,1)*rij(1)+magnt(k,2)*rij(2)+magnt(k,3)*rij(3)
                            Bstray(n0,1)=Bstray(n0,1)+(mu0*Ms/(4*pi))*(3*product*rij(1)-magnt(k,1))/rij_total**3
                            Bstray(n0,2)=Bstray(n0,2)+(mu0*Ms/(4*pi))*(3*product*rij(2)-magnt(k,2))/rij_total**3
                            Bstray(n0,3)=Bstray(n0,3)+(mu0*Ms/(4*pi))*(3*product*rij(3)-magnt(k,3))/rij_total**3
                        endif
                    enddo
                    Estray(n0)=-0.5d0*Ms*(magnt(n0,1)*Bstray(n0,1)+magnt(n0,2)*Bstray(n0,2)+magnt(n0,3)*Bstray(n0,3))
                enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !$acc kernels
                do n0=1,Natom
                    Beff(n0,1)=Bexc(n0,1)+Bani(n0,1)+BDMI(n0,1)+Bext(n0,1)+Bstray(n0,1)!+Bthermal(n0,1)
                    Beff(n0,2)=Bexc(n0,2)+Bani(n0,2)+BDMI(n0,2)+Bext(n0,2)+Bstray(n0,2)!+Bthermal(n0,2)
                    Beff(n0,3)=Bexc(n0,3)+Bani(n0,3)+BDMI(n0,3)+Bext(n0,3)+Bstray(n0,3)!+Bthermal(n0,3)
                    Eeff(n0)=Eexc(n0)+Eani(n0)+EDMI(n0)+Eext(n0)+Estray(n0)
                enddo
                !$acc end kernels

                !$acc parallel loop private(product1,product2,product3,FLterm,DLterm)
                do n0=1,Natom
                    product1=magnt(n0,2)*Beff(n0,3)-magnt(n0,3)*Beff(n0,2)
                    product2=magnt(n0,3)*Beff(n0,1)-magnt(n0,1)*Beff(n0,3)
                    product3=magnt(n0,1)*Beff(n0,2)-magnt(n0,2)*Beff(n0,1)
                    
                    FLterm(1)=magnt(n0,2)*SPorientation(3)-magnt(n0,3)*SPorientation(2)
                    FLterm(2)=magnt(n0,3)*SPorientation(1)-magnt(n0,1)*SPorientation(3)
                    FLterm(3)=magnt(n0,1)*SPorientation(2)-magnt(n0,2)*SPorientation(1)
                    DLterm(1)=magnt(n0,2)*FLterm(3)-magnt(n0,3)*FLterm(2)
                    DLterm(2)=magnt(n0,3)*FLterm(1)-magnt(n0,1)*FLterm(3)
                    DLterm(3)=magnt(n0,1)*FLterm(2)-magnt(n0,2)*FLterm(1)
                    
                    dm(n0,1,i)=(dt*gamma/(1+alpha**2))*(-product1-alpha*(magnt(n0,2)*product3-magnt(n0,3)*product2))
                    dm(n0,2,i)=(dt*gamma/(1+alpha**2))*(-product2-alpha*(magnt(n0,3)*product1-magnt(n0,1)*product3))
                    dm(n0,3,i)=(dt*gamma/(1+alpha**2))*(-product3-alpha*(magnt(n0,1)*product2-magnt(n0,2)*product1))
                    
                    dm(n0,1,i)=dm(n0,1,i)+dt*betaFL*FLterm(1)+dt*betaDL*DLterm(1)
                    dm(n0,2,i)=dm(n0,2,i)+dt*betaFL*FLterm(2)+dt*betaDL*DLterm(2)
                    dm(n0,3,i)=dm(n0,3,i)+dt*betaFL*FLterm(3)+dt*betaDL*DLterm(3)
                enddo
            enddo
            
            !$acc parallel loop private(magnttotal)
            do n0=1,Natom
                magnt(n0,1)=magnt_old(n0,1)+(dm(n0,1,1)+2*dm(n0,1,2)+2*dm(n0,1,3)+dm(n0,1,4))/6.d0
                magnt(n0,2)=magnt_old(n0,2)+(dm(n0,2,1)+2*dm(n0,2,2)+2*dm(n0,2,3)+dm(n0,2,4))/6.d0
                magnt(n0,3)=magnt_old(n0,3)+(dm(n0,3,1)+2*dm(n0,3,2)+2*dm(n0,3,3)+dm(n0,3,4))/6.d0
                magnttotal=sqrt(magnt(n0,1)**2+magnt(n0,2)**2+magnt(n0,3)**2)
                magnt(n0,1)=magnt(n0,1)/magnttotal
                magnt(n0,2)=magnt(n0,2)/magnttotal
                magnt(n0,3)=magnt(n0,3)/magnttotal
            enddo

            !!!!!!!!!!!!!!!!!!!!!!!calculateB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!1.calculateBexc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc parallel loop private(nn)
            do n0=1,Natom
                Bexc(n0,1)=0;    Bexc(n0,2)=0;    Bexc(n0,3)=0;
                if((flagLattice==1).or.(flagLattice==3)) then
                    do nn=1,12
                        if(neighboratom(n0,nn)/=0) then
                            Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                            Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                            Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                        endif
                    enddo
                endif
                if(flagLattice==2) then      !only for SC
                    if(Nz/=1) then
                        do nn=1,2
                            if(neighboratom(n0,nn)/=0) then
                                Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                            endif
                        enddo
                    endif
                    if(Ny/=1) then
                        do nn=3,4
                            if(neighboratom(n0,nn)/=0) then
                                Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                            endif
                        enddo
                    endif
                    if(Nx/=1) then
                        do nn=5,6
                            if(neighboratom(n0,nn)/=0) then
                                Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                                Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                                Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                            endif
                        enddo
                    endif
                endif
                if(flagLattice==4) then
                    do nn=1,6
                        if(neighboratom(n0,nn)/=0) then
                            Bexc(n0,1)=Bexc(n0,1)+(J/Ms)*magnt(neighboratom(n0,nn),1)
                            Bexc(n0,2)=Bexc(n0,2)+(J/Ms)*magnt(neighboratom(n0,nn),2)
                            Bexc(n0,3)=Bexc(n0,3)+(J/Ms)*magnt(neighboratom(n0,nn),3)
                        endif
                    enddo
                endif
                Eexc(n0)=-0.5d0*Ms*(magnt(n0,1)*Bexc(n0,1)+magnt(n0,2)*Bexc(n0,2)+magnt(n0,3)*Bexc(n0,3))
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!2.calculateBani!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc kernels
            do n0=1,Natom
                product=0
                do i=1,3
                    product=product+magnt(n0,i)*kaxis(i)
                enddo
                Bani(n0,1)=(2*Ku/Ms)*product*kaxis(1)
                Bani(n0,2)=(2*Ku/Ms)*product*kaxis(2)
                Bani(n0,3)=(2*Ku/Ms)*product*kaxis(3)
                Eani(n0)=-0.5d0*Ms*(magnt(n0,1)*Bani(n0,1)+magnt(n0,2)*Bani(n0,2)+magnt(n0,3)*Bani(n0,3))
            enddo
            !$acc end kernels
            !!!!!!!!!!!!!!!!!!!!!!!!!!3.calculateBDMI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc parallel loop private(rij,rij_total,nn)
            do n0=1,Natom
                BDMI(n0,1)=0;      BDMI(n0,2)=0;         BDMI(n0,3)=0
                if(flagDMI==1) then
                    if((flagLattice==1).or.(flagLattice==3)) then
                        do nn=1,12
                            if(neighboratom(n0,nn)/=0) then
                                rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                rij(3)=0                       
                                rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                rij(1)=rij(1)/rij_total
                                rij(2)=rij(2)/rij_total
                                rij(3)=rij(3)/rij_total
                                BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                            endif
                        enddo
                    endif
                    if(flagLattice==2) then
                        if(Nz/=1) then
                            do nn=1,2
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                    rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                    rij(3)=0                     
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(Ny/=1) then
                            do nn=3,4
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                    rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                    rij(3)=0
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(Nx/=1) then
                            do nn=5,6
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                    rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                    rij(3)=0
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                    endif
                    if(flagLattice==4) then
                        do nn=1,6
                            if(neighboratom(n0,nn)/=0) then
                                rij(1)=-(atomposition(neighboratom(n0,nn),2)-atomposition(n0,2))
                                rij(2)=(atomposition(neighboratom(n0,nn),1)-atomposition(n0,1))
                                rij(3)=0                     
                                rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                rij(1)=rij(1)/rij_total
                                rij(2)=rij(2)/rij_total
                                rij(3)=rij(3)/rij_total
                                BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                            endif
                        enddo
                    endif
                endif
                if(flagDMI==2) then
                    if((flagLattice==1).or.(flagLattice==3)) then
                        do nn=1,12
                            if(neighboratom(n0,nn)/=0) then
                                rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                rij(1)=rij(1)/rij_total
                                rij(2)=rij(2)/rij_total
                                rij(3)=rij(3)/rij_total
                                BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                            endif
                        enddo
                    endif
                    if(flagLattice==2) then
                        if(Nz/=1) then
                            do nn=1,2
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                    rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                    rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(Ny/=1) then
                            do nn=3,4
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                    rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                    rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                        if(Nx/=1) then
                            do nn=5,6
                                if(neighboratom(n0,nn)/=0) then
                                    rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                    rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                    rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                    rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                    rij(1)=rij(1)/rij_total
                                    rij(2)=rij(2)/rij_total
                                    rij(3)=rij(3)/rij_total
                                    BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                    BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                    BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                                endif
                            enddo
                        endif
                    endif
                    if(flagLattice==4) then
                        do nn=1,6
                            if(neighboratom(n0,nn)/=0) then
                                rij(1)=atomposition(neighboratom(n0,nn),1)-atomposition(n0,1)
                                rij(2)=atomposition(neighboratom(n0,nn),2)-atomposition(n0,2)
                                rij(3)=atomposition(neighboratom(n0,nn),3)-atomposition(n0,3)
                                rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                                rij(1)=rij(1)/rij_total
                                rij(2)=rij(2)/rij_total
                                rij(3)=rij(3)/rij_total
                                BDMI(n0,1)=BDMI(n0,1)+(D/Ms)*(rij(2)*magnt(neighboratom(n0,nn),3)-rij(3)*magnt(neighboratom(n0,nn),2))
                                BDMI(n0,2)=BDMI(n0,2)+(D/Ms)*(rij(3)*magnt(neighboratom(n0,nn),1)-rij(1)*magnt(neighboratom(n0,nn),3))
                                BDMI(n0,3)=BDMI(n0,3)+(D/Ms)*(rij(1)*magnt(neighboratom(n0,nn),2)-rij(2)*magnt(neighboratom(n0,nn),1))
                            endif
                        enddo
                    endif
                endif
                EDMI(n0)=-0.5d0*Ms*(magnt(n0,1)*BDMI(n0,1)+magnt(n0,2)*BDMI(n0,2)+magnt(n0,3)*BDMI(n0,3))
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4.calculateBext!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc kernels
            do n0=1,Natom
                Bext(n0,1)=B0(1)
                Bext(n0,2)=B0(2)
                Bext(n0,3)=B0(3)
                Eext(n0)=-Ms*(magnt(n0,1)*Bext(n0,1)+magnt(n0,2)*Bext(n0,2)+magnt(n0,3)*Bext(n0,3))
            enddo
            !$acc end kernels
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5.calculateBstray!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc parallel loop private(rij,rij_total,product,i)
            do n0=1,Natom
                Bstray(n0,1)=0;      Bstray(n0,2)=0;         Bstray(n0,3)=0
                do i=1,Natom
                    if(i/=n0) then
                        rij(1)=atomposition(i,1)-atomposition(n0,1)
                        rij(2)=atomposition(i,2)-atomposition(n0,2)
                        rij(3)=atomposition(i,3)-atomposition(n0,3)
                        rij_total=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
                        rij(1)=rij(1)/rij_total
                        rij(2)=rij(2)/rij_total
                        rij(3)=rij(3)/rij_total
                        product=magnt(i,1)*rij(1)+magnt(i,2)*rij(2)+magnt(i,3)*rij(3)
                        Bstray(n0,1)=Bstray(n0,1)+(mu0*Ms/(4*pi))*(3*product*rij(1)-magnt(i,1))/rij_total**3
                        Bstray(n0,2)=Bstray(n0,2)+(mu0*Ms/(4*pi))*(3*product*rij(2)-magnt(i,2))/rij_total**3
                        Bstray(n0,3)=Bstray(n0,3)+(mu0*Ms/(4*pi))*(3*product*rij(3)-magnt(i,3))/rij_total**3
                    endif
                enddo
                Estray(n0)=-0.5d0*Ms*(magnt(n0,1)*Bstray(n0,1)+magnt(n0,2)*Bstray(n0,2)+magnt(n0,3)*Bstray(n0,3))
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculateB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$acc kernels
            do n0=1,Natom
                Beff(n0,1)=Bexc(n0,1)+Bani(n0,1)+BDMI(n0,1)+Bext(n0,1)+Bstray(n0,1)!+Bthermal(n0,1)
                Beff(n0,2)=Bexc(n0,2)+Bani(n0,2)+BDMI(n0,2)+Bext(n0,2)+Bstray(n0,2)!+Bthermal(n0,2)
                Beff(n0,3)=Bexc(n0,3)+Bani(n0,3)+BDMI(n0,3)+Bext(n0,3)+Bstray(n0,3)!+Bthermal(n0,3)
                Eeff(n0)=Eexc(n0)+Eani(n0)+EDMI(n0)+Eext(n0)+Estray(n0)
            enddo
            !$acc end kernels	
        enddo
        !$acc end data
    endif
    

    close(8)
    close(16)
    
    deallocate(atomposition)
    deallocate(magnt)          
    deallocate(Beff)           
    deallocate(Bexc)    
    deallocate(Bani)    
    deallocate(Bext)
    deallocate(Bstray)
    deallocate(Eeff)      
    deallocate(Eexc)      
    deallocate(Eani)      
    deallocate(Eext)
    deallocate(Estray)
    deallocate(dm)
    deallocate(magnt_old)
    
    call cpu_time(finish)

    write(*,'("Time=",f20.3," seconds.")') finish-start
   
    
end program
    
    
    
