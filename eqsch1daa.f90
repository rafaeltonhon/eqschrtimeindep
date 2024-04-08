program eqsch1daa
    implicit none

    ! Declarando os tipos
    integer, parameter :: i4 = selected_int_kind(8)
    integer, parameter :: r8 = selected_real_kind(8,15)

    ! Declarando os vetores
    real(kind=r8), allocatable, dimension(:) :: D, E, W, work, x
    real(kind=r8), allocatable, dimension(:,:) :: Z
    integer(kind=i4), allocatable, dimension(:) :: isuppz, iwork

    ! Declarando as variáveis
    real(kind=r8) :: vl, vu, dx, l, lesq, ldir, sigma, epsilon, k
    integer(kind=i4) :: i, il, iu, m, ldz, info, tipo_pot
    integer(kind=i4) :: nzc, lwork, liwork, n, j
    character(len=1) :: jobz, range
    character(len=20) :: nome_pot
    logical :: tryrac
    
    ! Lendo as variáveis de entrada
    open(unit=1,file="eqsch1daa1.txt",status="old")
    read(1,*) n
    read(1,*) iu
    read(1,*) l
    read(1,*) ldir
    read(1,*) lesq
    read(1,*) nome_pot
    read(1,*) sigma
    read(1,*) epsilon
    read(1,*) k
    close(1)

    ! Selecionando o tipo de potencial
    select case(nome_pot)
    case("Caixa")
        tipo_pot = 1
    case("Harmonico")
        tipo_pot = 2
    case("Lennard-Jonnes")
        tipo_pot = 3
    case("Culomb")
        tipo_pot = 4
    case default
        write(*,*) "Nenhum potencial válido foi selecionado!"
        stop
    end select

    ! Declarando as variáveis da função que fará o calculo dos auto-valores e auto-vetores
    jobz = "v"
    range = "i"
    il = 1
    m = iu
    ldz = n
    nzc = iu - il + 1
    tryrac = .true.
    lwork = 18*n
    liwork = 10*n
    dx = (ldir - lesq)/(n+1)

    allocate(Z(ldz,m))
    allocate(isuppz(2*m))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(D(1:n))
    allocate(E(1:n))
    allocate(W(1:n))
    allocate(x(0:n+1))

    ! Preenchendo o vetor posição
    x(0) = -lesq
    x(n+1) = +ldir
    do i=1,n
        x(i) = lesq+i*dx
    enddo
    ! Preenchendo as diagonais da Hamiltoniana
    do i=1,n
            D(i) = 1.0_r8/dx**2 + potencial(x(i), l, sigma, epsilon, k, tipo_pot)
            E(i) = -0.5_r8/dx**2
    enddo

    call dstemr(jobz, range, n, D, E, vl, vu, il, iu, m, W, Z, ldz, &
    & nzc, isuppz, tryrac, work, lwork, iwork, liwork, info)
    
    ! Normalizando os auto-vetores
    !if (tipo_pot.Eq.4)then
    !    do i=1,n
    !        do j=1,m
    !            Z(i,j)=Z(i,j)/x(i)
    !        enddo
    !    enddo
    !endif
    call normaliza(Z, m, n, dx/2.0_r8)

    ! Saindo com os auto-valores normalizados
    open(unit=2,file="eqsch1daa2.txt")
    open(unit=3,file="eqsch1daa3.txt")
    do i=1,m
        write(2,10) i,"=",W(i)
        10 format(i4,a1,f15.5)
    enddo
    write(3,*) -lesq, 0.0_r8
    do i=1,n
        if(tipo_pot.eq.4)then
            write(3,*) x(i), (Z(i,j) , j = 1,m)
        else
            write(3,*) x(i), (Z(i,j) , j = 1,m)
        endif
    enddo
    write(3,*) ldir, 0.0_r8
    close(2)
    close(3)

    ! Liberando a memória utilizada
    deallocate(Z)
    deallocate(isuppz)
    deallocate(work)
    deallocate(iwork)
    deallocate(D)
    deallocate(E)
    deallocate(W)
    deallocate(x)
    write(*,*) info ! Apenas afirmando se tudo correu bem ou não
    contains

    ! Função potnecial
    function potencial(x, l, sigma, epsilon, k, tipo_pot)
        real(kind=r8) :: x, l, k, potencial, sigma, epsilon
        integer(kind=i4) :: tipo_pot
        select case(tipo_pot)
        case(1)
            if((-l<=x).and.(x<=l))then
                potencial = 0.0_r8
            else
                potencial = 1.0E+15
            endif
            return
        case(2)
            potencial = 0.5_r8*x**2
            return
        case(3)
            potencial = 4.0_r8*epsilon*((sigma/x)**12-(sigma/x)**6)
            return
        case(4)
            potencial = 0.5_r8*k*(k+1)/x**2-1/x
            return
        case default
            stop
        end select
    end function potencial

    ! Subrotina de normalização
    subroutine normaliza(vetor, colunas, pontos, h)
        real(kind=r8), dimension(pontos,colunas), intent(inout) :: vetor
        real(kind=r8), intent(in) :: h
        real(kind=r8) :: prob
        integer(kind=i4), intent(in) :: pontos, colunas
        integer(kind=i4) :: i, j

        ! Fazemos a normalização para cada coluna
        do i=1,colunas
            prob = 0.0_r8
            do j=2,pontos-1
                prob=prob+(h/3.0_r8)*(vetor(j+1,i)**2+4.0_r8*vetor(j,i)**2+vetor(j-1,i)**2)
            enddo

            do j=1,pontos
                vetor(j,i)=vetor(j,i)**2/prob
            enddo
        enddo
    end subroutine normaliza
end program eqsch1daa