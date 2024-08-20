program vmc
    implicit none
    ! defining our types
    integer, parameter :: r8=selected_real_kind(15,9)
    integer, parameter :: i4=selected_int_kind(9)

    ! defining our parameters
    real(kind=r8), parameter :: pi=acos(-1.0_r8) ! pi
    real(kind=r8), parameter :: a0=0.529 ! bohr radius
    real(kind=r8), parameter :: hc=1973 ! hbar times c
    real(kind=r8), parameter :: mec2=0.5110e6 ! electron reast mass

    ! defining our arrays
    real(kind=r8), allocatable, dimension(:,:) :: rdist

    ! defining our variables
    real(kind=r8) :: delta, acc, w1, w2, ene, sigma, eps
    integer(kind=i4) :: i, nconf, dim, lpot

    eps=hc**2.0_r8/(mec2*a0**2.0)
    ! reading the input values
    open(unit=100,file='vmc-input.txt')
    read(100,*) nconf   ! number of points that we want
    read(100,*) dim     ! dimension of the problem
    read(100,*) delta   ! delta, if this is negative we do not knew the good one, and will just make trys
    read(100,*) lpot    ! pot type 1-> square well, 2-> 2d harmonic oscilator
    read(100,*) w1, w2  ! first and second variational parameters, for lpot=2 they are the gaussian width
    close(100)
    ! allocating the memory
    allocate(rdist(nconf,dim))
    ! if delta is negative we'll search for thwe right one
    acc=0.0_r8
    if(delta.lt.0.0_r8)then
        do i=1,100
            read(*,*) delta
            if(delta.lt.0.0)stop 'process stoped, i hope youve found the rigt delta'
            call gives_dist(rdist,nconf,delta,acc,dim,lpot,w1,w2)
            print*,'the acceptance percentage is = ',acc
        enddo
    else ! if we already know the good delta, we just return an aleatory distribution
        call gives_dist(rdist,nconf,delta,acc,dim,lpot,w1,w2)
        print*,'the acceptance percentage is = ',acc
        call stat_ana(rdist,nconf,dim,lpot,w1,w2,ene,sigma)
    endif

    ! geting out with the results
    open(unit=200,file='vmc-energy.txt',access='append')
    open(unit=300,file='vmc-test-function.txt')
    if(lpot.eq.1)then
        write(200,*) 'delta','&','energy ','&',' sigma ','&','acc'
        write(200,*) delta,'&',ene,'&',sigma,'&',acc
    elseif(lpot.eq.2)then
        write(200,*) 'energy ','&',' sigma ','&','acc','&''w1','&','w2'
        write(200,*) ene,'&',sigma,'&',acc,'&',w1,'&',w2
    else
        write(200,*) '$w$','&','\varepsilon','&',' sigma_V ','&',' $E_V$','&',' \sigma_V'
        write(200,*) w1,'&',ene,'&',sigma,'&',ene*eps,'&',sigma*eps
        print*,ene*eps
    endif
!    do i=1,nconf
!        write(300,*) rdist(i,:),phit(rdist(i,:),dim,lpot,w1,w2) 
!    enddo
    close(200)
    close(300)
    print*,0 ! if everything went right
    contains
    ! function with the test funtions
    function phit(r,dim,lpot,w1,w2)
        real(kind=r8) :: phit, r(1,dim), w1, w2, r1, r2
        integer(kind=i4) :: dim, lpot
        selectcase(lpot)
        case(1) ! the square well
            if((r(1,1).le.1.0_r8).and.(r(1,1).ge.0.0_r8))then
                phit=r(1,1)*(1.0_r8-r(1,1))
            else
                phit=0.0_r8
            endif
        case(2) ! the 2d harmonic oscilator
            phit=exp(-w1*r(1,1)**2.0_r8-w2*r(1,2)**2.0_r8)
        case(3) ! the helium
            r1=sqrt(r(1,1)**2.0_r8+r(1,2)**2.0_r8+r(1,3)**2.0_r8)
            r2=sqrt(r(1,4)**2.0_r8+r(1,5)**2.0_r8+r(1,6)**2.0_r8)
            phit=exp(-w1*(r1+r2))
        case default
            stop 'no valid test function was selected'
        endselect
    endfunction phit

    ! funtion with the hamiltonian operator action on the test function
    function Hphit(r,dim,lpot,w1,w2)
        real(kind=r8) :: Hphit, r(1,dim), w1, w2, r1, r2, r12
        integer(kind=i4) :: dim, lpot
        selectcase(lpot)
        case(1) ! the square well
            if((r(1,1).le.1.0_r8).and.(r(1,1).ge.0.0_r8))then
                Hphit=2.0_r8
            else
                Hphit=0.0_r8
            endif
        case(2) ! the 2d harmonic oscilator
            Hphit=((0.5_r8-2.0_r8*w1*w1)*r(1,1)**2.0_r8+(2.0_r8-2.0_r8*w2*w2)*r(1,2)**2.0_r8+w1+w2)*&
            &exp(-w1*r(1,1)**2.0_r8-w2*r(1,2)**2.0_r8)
        case(3) ! the helium
            r1=sqrt(r(1,1)**2.0_r8+r(1,2)**2.0_r8+r(1,3)**2.0_r8)
            r2=sqrt(r(1,4)**2.0_r8+r(1,5)**2.0_r8+r(1,6)**2.0_r8)
            r12=sqrt((r(1,1)-r(1,4))**2.0_r8+(r(1,2)-r(1,5))**2.0_r8+(r(1,3)-r(1,6))**2.0_r8)
            Hphit=(-0.5_r8*(2.0_r8*w1**2.0_r8-2.0_r8*w1/r1-2.0_r8*w1/r2)&
            &-2.0_r8/r1-2.0_r8/r2+1.0_r8/r12)*exp(-w1*(r1+r2))
        case default
            stop 'no valid test function was selected'
        endselect
    endfunction Hphit

    ! subroututine that give us a good distribution
    subroutine gives_dist(r,nconf,delta,acc,dim,lpot,w1,w2)
        real(kind=r8) :: r(nconf,dim), delta, acc, rnew(1,dim)
        real(kind=r8) :: epsilon, zeta, w1, w2
        integer(kind=i4) :: i, j, nconf, dim, lpot

        rdist=0.0_r8
        ! initializing the process
        ! we gonna do the metropoles algorhitm a 100 to make sure 
        ! that the initial configuration does not intenterfer in the result
        do i=1,dim
            call random_number(r(1,i))
        enddo
        do i=1,400
            call random_number(zeta)
            do j=1,dim
                call random_number(epsilon)
                rnew(1,j)=r(1,j)+(epsilon-0.5_r8)*delta
            enddo
            !print*,rnew(1,1),phit(rnew(1,1)),r(1,1),phit(r(1,1)),zeta,phiT(rnew(1,1))/phit(r(1,1))
            if(phiT(rnew(1,:),dim,lpot,w1,w2)**2.0_r8/phit(r(1,:),dim,lpot,w1,w2)**2.0_r8.gt.zeta)then
                r(1,1)=rnew(1,1)
            endif
        enddo

        ! we've the initial configuration, let's get the reast of it
        acc=0
        do i=1,nconf-1
            call random_number(zeta)
            do j=1,dim
                call random_number(epsilon)
                rnew(1,j)=r(i,j)+(epsilon-0.5_r8)*delta
            enddo
            if(phit(rnew(1,:),dim,lpot,w1,w2)**2.0_r8/phit(r(i,:),dim,lpot,w1,w2)**2.0_r8.gt.zeta)then
                r(i+1,:)=rnew(1,:)
                acc=acc+1
            else
                r(i+1,:)=r(i,:)
            endif
        enddo
        acc=acc/nconf
        !print*,r
    endsubroutine gives_dist

    ! subroutine that calculates the media, standard deviation and all that kind of shit
    subroutine stat_ana(r,nconf,dim,lpot,w1,w2,enev,sigv)
        real(kind=r8) :: r(nconf,dim), enev, sigv, w1, w2
        integer(kind=i4) :: i, nconf, dim, lpot

        ! calculating the energy
        enev=0.0_r8
        do i=1,nconf
            enev=enev+Hphit(r(i,:),dim,lpot,w1,w2)/phit(r(i,:),dim,lpot,w1,w2)
        enddo
        enev=enev/nconf

        ! calculating the standard deviation
        sigv=0.0_r8
        do i=1,nconf
            sigv=sigv+(Hphit(r(i,:),dim,lpot,w1,w2)/phit(r(i,:),dim,lpot,w1,w2)-enev)**2.0_r8
        enddo
        sigv=sqrt(sigv/(nconf-1))
        sigv=sigv/sqrt(real(nconf))
        print*,'energy =',enev
        print*,'sigma =',sigv
    endsubroutine stat_ana
endprogram vmc