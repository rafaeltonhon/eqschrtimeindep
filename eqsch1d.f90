program eqsch1d
    implicit none
    ! declaring our types
    integer, parameter :: r8=selected_real_kind(15,9)
    integer, parameter :: i4=selected_int_kind(9)

    ! declaring our parameters
    real(kind=r8), parameter :: pi=acos(-1.0_r8)

    ! declaring our vectors
    real(kind=r8), allocatable, dimension(:) :: psi, x

    ! declaring our variables
    real(kind=r8) :: xin, xout, energy, denergy, limene, dx, l, xmeet, psi0, psip0
    integer(kind=i4) :: np, lint, lpot, ldir, i, n

    open(unit=100,file='eqsch1d-input.txt',status='old')
    read(100,*) np ! number of points
    read(100,*) xin ! x in
    read(100,*) xout ! x out
    read(100,*) energy ! initial energy
    read(100,*) lpot ! lpot
    read(100,*) lint ! lint
    read(100,*) l ! the width of the square well
    read(100,*) denergy ! denergy
    read(100,*) limene ! limit energy
    read(100,*) xmeet ! x of the meet point
    read(100,*) psi0,psip0 ! psi0 and psi prime 0
    read(100,*) n
    close(100)
    ldir=1
    dx=(xout-xin)/np
    allocate(psi(0:np),x(0:np))
    call scheqsolver(psi,x,np,xin,xout,energy,lint,ldir,lpot,l,denergy,limene,xmeet,psi0,psip0)
    ! writring the wave function
    open(unit=200,file='eqsch1d-psi.txt')
    do i=0,np
            write(200,*) x(i), psi(i)
    enddo
    close(200)

    open(unit=300,file='eqsch1d-energy.txt',access='append')
    if(lpot.eq.2)then
        write(300,*) n,'&',energy,'&',pi*pi*n*n/8.0_r8,'&',energy-pi*pi*n*n/8.0_r8
    else
        write(300,*) energy
    endif
    close(300)
    deallocate(psi,x)

    contains
    ! function that give us the potentials in the game
    real(kind=r8) function v(x,l,lpot)
        real(kind=r8) :: x, l
        integer(kind=i4) :: lpot
        selectcase(lpot)
        case(1) ! the free particle
            v=0.0_r8
        case(2) ! the square well
            if((x.gt.l).or.(x.lt.-l))then
                v=10000000.0_r8!infinity
            else
                v=0.0_r8
            endif
            !v=0.5_r8*x**2.0_r8 ! uncoment this line to solve the harmonic oscilator potential
        case(3) ! the lenard-jones potential
            v=4.0_r8*l*(x**(-12.0_r8)-x**(-6.0_r8))
        case default
            stop "no valid potential was select"
        endselect
    endfunction v

    ! funtion that integrate the vector by the simpson's 1/3 rule
    real(kind=r8) function simpson(psi,dx,np)
        real(kind=r8) :: psi(0:np), dx, sum
        integer(kind=i4) :: i, np
        sum=0.0_r8
        do i=1,np-1,2
            sum=sum+psi(i-1)+4.0_r8*psi(i)+psi(i+1)
        enddo
        simpson=dx*sum/3.0_r8
    endfunction simpson

    ! function that give us the derive given the two poins of the vector
    function psiprime(psir,psil,h)
        real(kind=r8) :: h, psil, psir, psiprime
        psiprime=(psir-psil)/(2.0_r8*h)
    endfunction psiprime

    ! function g of the numerov method
    function gnum(x,l,lpot,e)
        real(kind=r8) :: x, l, e, gnum
        integer(kind=i4) :: lpot
        gnum=2.0_r8*(e-v(x,l,lpot))
    endfunction gnum

    ! subroutine that does the integration of the schrodinger equation
    subroutine scheqint(psi,x,np,xout,xin,energy,lint,ldir,lpot,l,xmeet,lmeet)
        real(kind=r8), dimension(0:np) , intent(inout) :: psi, x
        real(kind=r8), intent(in) :: xin, xout, energy, l, xmeet
        real(kind=r8) :: dx, epsmeet, gp1, gm1, gc
        integer(kind=i4), intent(in) :: lpot, ldir, lint, np
        integer(kind=i4) :: i, j, lmeet
        dx=abs(xout-xin)/np
        epsmeet=dx/10.0_r8
        selectcase(lint)
        case(1) ! the simetric derivative
            if(ldir.eq.1)then ! integrating to the right
                do i=1,np-1
                    x(i+1)=x(i)+dx
                    psi(i+1)=2.0_r8*psi(i)-psi(i-1)-2.0_r8*(energy-v(x(i),l,lpot))*psi(i)*dx**2.0_r8
                    ! if we're on the meet point, we return this i
                    if(abs(xmeet-x(i)).lt.epsmeet)then
                        lmeet=i
                        !xnew=x(i)
                    endif
                enddo
            else
                do j=1,np-1
                    i=np-j
                    x(i-1)=x(i)-dx
                    psi(i-1)=2.0_r8*psi(i)-psi(i+1)-2.0_r8*(energy-v(x(i),l,lpot))*psi(i)*dx**2.0_r8
                    ! if we're on the meet point, we return this i
                    if(abs(xmeet-x(i)).lt.epsmeet)then
                        lmeet=i
                        !xnew=x(i)
                    endif
                enddo
            endif
        case(2) ! the numerov method
            if(ldir.eq.1)then ! integrating to the right
                do i=1,np-1
                    x(i+1)=x(i)+dx
                    gp1=1.0_r8+dx**2.0_r8*gnum(x(i+1),l,lpot,energy)/12.0_r8
                    gc=1.0_r8-5.0_r8*dx**2.0_r8*gnum(x(i),l,lpot,energy)/12.0_r8
                    gm1=1.0_r8+dx**2.0_r8*gnum(x(i-1),l,lpot,energy)/12.0_r8

                    psi(i+1)=2.0_r8*psi(i)*gc-psi(i-1)*gm1
                    psi(i+1)=psi(i+1)/gp1
                    ! if we're on the meet point, we return this i
                    if(abs(xmeet-x(i)).lt.epsmeet)then
                        lmeet=i
                        !xnew=x(i)
                    endif
                enddo
            else
                do j=1,np-1
                    i=np-j
                    x(i-1)=x(i)-dx
                    gp1=1.0_r8+dx**2.0_r8*gnum(x(i+1),l,lpot,energy)/12.0_r8
                    gc=1.0_r8-5.0_r8*dx**2.0_r8*gnum(x(i),l,lpot,energy)/12.0_r8
                    gm1=1.0_r8+dx**2.0_r8*gnum(x(i-1),l,lpot,energy)/12.0_r8

                    psi(i-1)=2.0_r8*psi(i)*gc-psi(i+1)*gp1
                    psi(i-1)=psi(i-1)/gm1
                    
                    ! if we're on the meet point, we return this i
                    if(abs(xmeet-x(i)).lt.epsmeet)then
                        lmeet=i
                        !xnew=x(i)
                    endif
                enddo
            endif
        case default 
            stop 'no valid integrator selected'
        endselect
    endsubroutine scheqint

    ! subroutnie that solves the schrodinger equation
    subroutine scheqsolver(psi,x,np,xin,xout,energy,lint,ldir,lpot,l,denergy,limene,xmeet,psi0,psip0)
        real(kind=r8), dimension(0:np), intent(inout) :: x, psi
        real(kind=r8), intent(in) :: xin, xout, l, limene, xmeet, psi0,psip0
        real(kind=r8), intent(inout) :: denergy, energy
        real(kind=r8) :: dx, psir(0:np), psil(0:np), xr(0:np), xl(0:np), a
        real(kind=r8) :: psiprimel, psiprimer, deltaprime, psit, k, c
        integer(kind=i4), intent(in) :: lint, lpot, ldir, np
        integer(kind=i4) :: sig, lsig, il, ir
        psi=0.0_r8
        x=0.0_r8
        dx=(xout-xin)/(np-1)
        selectcase(lpot) ! what problem are we solving
        case(1) ! the free particle
            dx=2.0_r8*dx
            psil=0.0_r8
            psir=0.0_r8
            xr=0.0_r8
            xl=0.0_r8

            ! incitializig the function
            ! for the right
            psir(0)=psi0
            psir(1)=psi0+psip0*dx
            xr(1)=xin+dx

            ! for the left
            psil(np)=psi0
            psil(np-1)=psi0-psip0*dx
            xl(np-1)=xin-dx

            ! integrating the function
            call scheqint(psir,xr,np,2.0_r8*xout,0.0_r8,energy,lint,ldir,lpot,l,xmeet,il)
            call scheqint(psil,xl,np,2.0_r8*xout,0.0_r8,energy,lint,-ldir,lpot,l,xmeet,il)

            ! new we join the two funtions in the desired interval (i know its not the easier way, but i saw that we want the integration from zero right now)
            do i=0,np
                if(i.lt.np/2)then
                    psi(i)=psil(np/2+i)
                    x(i)=xl(np/2+i)
                    !print*,i,xl(i)
                else
                    psi(i)=psir(i-np/2)
                    x(i)=xr(i-np/2)
                endif
            enddo
            ! normalizing
            psi=psi/simpson(psi**2.0_r8,dx,np)**0.5_r8
            ! writing the difference between the numerical and theorical wave function
            k=sqrt(2.0_r8*energy)
            c=((2.0_r8*k*l-sin(2.0_r8*k*l))*psip0*psip0+(k*k*sin(2.0_r8*k*l)+2.0_r8*l*k*k*k)*psi0*psi0)/(2.0_r8*k*k*k)
            c=c**0.50_r8
            do i=0,np
                psit=(psi0*cos(k*x(i))+psip0*sin(k*x(i))/k)/c
                write(500,*) x(i), psi(i), psit, psi(i)-psit
            enddo
        case(2) ! the square well by the shotting and macthing method
            sig=1
            lsig=1
            x=0.0_r8
            x(1)=dx
            do i=1,100
                psi=0.0_r8
                psi(0)=psi0
                psi(1)=psi0+psip0*dx
                ! new we integrate this function with the given energy
                call scheqint(psi,x,np,xout,xin,energy,lint,ldir,lpot,l,xmeet,il)

                ! new we check the siginal of diverge of psi in the origin
                if(psi(np).gt.0.0_r8)sig=1
                if(psi(np).lt.0.0_r8)sig=-1

                if((i.gt.1).and.(sig*lsig.lt.0))denergy=-denergy/2.0_r8
                energy=energy+denergy
                lsig=sig
                if(abs(denergy).lt.limene)exit
            enddo
            psi=psi/sqrt(2.0_r8*simpson(psi**2.0_r8,dx,np)) ! the sqrt2 is because i just found half of the wave function
            print*,energy
        case(3) ! the lenard-jones potential with the math method
            psil=0.0_r8
            psir=0.0_r8
            xr=0.0_r8
            xl=0.0_r8
            do i=1,100
                ! initializing the funtions
                psir(0)=0.0_r8
                psir(1)=dx
                xr(0)=xin
                xr(1)=xin+dx

                psil(np)=0.0_r8
                psil(np-1)=-dx
                xl(np)=xout
                xl(np-1)=xout-dx

                ! integrating both
                call scheqint(psir,xr,np,xout,xin,energy,lint,ldir,lpot,l,xmeet,ir)
                call scheqint(psil,xl,np,xout,xin,energy,lint,-ldir,lpot,l,xmeet,il)
                !print*,ir,il,xr(ir),xl(il)
                
                ! new that we've the left and right wavefunctions we reescalete and check the derivative
                psil=psir(ir)*psil/psil(il)
                psiprimer=psiprime(psir(ir+1),psir(ir-1),dx)
                psiprimel=psiprime(psil(il+1),psil(il-1),dx)
                deltaprime=(psiprimer-psiprimel)
                !print*,psiprimer,psiprimel
                if(deltaprime.lt.0.0_r8)sig=-1
                if(deltaprime.gt.0.0_r8)sig=1
                if((i.gt.1).and.(sig*lsig.lt.0))denergy=-denergy/2.0_r8
                lsig=sig
                energy=energy+denergy
                if(abs(deltaprime).lt.limene)exit
            enddo

            ! new that we've the energy we marry the two wave functions to find the correct psi
            do i=0,np	
                if(i.le.ir)then
                    x(i)=xr(i)
                    psi(i)=psir(i)
                elseif(i.gt.ir)then
                    x(i)=xl(i)
                    psi(i)=psil(i)
                endif
                !if(i.lt.ir+200)write(900,*) xr(i),psir(i)
                !if(i.gt.il-200)write(800,*) xl(i),psil(i)
            enddo

            ! normalizar
            psi=psi/simpson(psi**2.0_r8,dx,np)**0.5_r8
            a=simpson(psi**2.0_r8,dx,np)**0.5_r8
            !psir=psir/a!simpson(psi**2.0_r8,dx,np)**0.5_r8
            !psil=psil/simpson(psi**2.0_r8,dx,np)**0.5_r8
            do i=0,np	
                !if(i.lt.ir+200)write(900,*) xr(i),psir(i)
                if(x(i).lt.xmeet) write(800,*) x(i),psi(i)!/a
                if(x(i).gt.xmeet) write(900,*) x(i),psi(i)!/a
                !if(i.gt.il-200)write(800,*) xl(i),psil(i)
            enddo
            print*,energy,'fuck you'
        case default
            stop 'no valid potential selected'
        endselect
    endsubroutine scheqsolver
endprogram eqsch1d
