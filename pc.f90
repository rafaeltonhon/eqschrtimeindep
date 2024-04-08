program PC    ! Partíclua confinada        
    ! Declaração das variáveis
    implicit none
    real * 8, allocatable, dimension(:) :: Psi
    real * 8 :: Energia, FatorDeExplosao, dx, de
    real * 8 :: Probabilidade, UltimaEnergia
    integer * 8 Divergiu, UltimaDiv
    integer * 8 :: N, i
    character * 1 :: Paridade
    
    ! Recebendo as condições iniciais
    write(*,*) "Entrar com o N :"
    read(*,*) N
    write(*,*) "Entre com a energia"
    read(*,*) Energia
    write(*,*) "Qual a paridade: [p] par [i] para impar :"
    read(*,*) Paridade

    ! Iniciando as variáveis globais
    dx = 1.3d0/N
    FatorDeExplosao = 1.1d0
    de = 0.5d0
    Divergiu = 0
    UltimaDiv = 0

    allocate(psi(-N:N))
    psi=0.0d0
    if(Paridade.Eq.'p')then   
        psi(0) = 1.0d0      ! Inicialização de função par.
        psi(1) = 1.0d0
    else if(Paridade.Eq.'i')then
        psi(0) = 0.0d0
        psi(1) = -dx
    else
        write(*,*) "Erro!!! Tente novamente."
    endif

    ! Vamos buscar a energia correta
    do while(abs(dE).Ge.1.0e-5)
        do i=2,N            ! Integrando a primeira metade da função de onda.
            psi(i)=EvoluiPsi(psi(i-1), psi(i-2), i*dx, Energia, dx, 0)
            if(abs(Psi(i)).Ge.FatorDeExplosao)then
                goto 10
            endif
        enddo
        10 continue

        ! Atualizando a divergência atual.
        if(psi(i).Gt.0.0d0)then
            Divergiu = 1    ! Divergiu para o positivo, energia menor do que deveria.
        else
            Divergiu = -1   ! Divergiu para o negativo, energia maior do que deveria.
        endif

        ! Verificando se a divergência está no mesmo sentido anterior.
        ! Se tiver mudado de direção passamos da energia real e temos que inverter o incremento.
        if(Divergiu*UltimaDiv.Lt.0)then
            dE = -dE/2.0d0
        endif
        UltimaDiv = Divergiu
        Energia=Energia+dE
        ! Deslocamos o código caso seja detectada a entrada num loop infinito.
        if(UltimaEnergia.Eq.Energia)then
            goto 20
        endif
        UltimaEnergia=Energia
    enddo
    20 continue

    ! Integrando a segunda metade da função de onda.
    if(Paridade.Eq.'p')then
        psi(-1) = psi(0)
    else if(Paridade.Eq.'i')then
        psi(-1) = -dx
    endif
    do i=2,N    ! Integrando a segunda metade da função de onda.
        psi(-i) = EvoluiPsi(psi(1-i), psi(2-i), -i*dx, Energia, dx, 1)
        if(abs(Psi(-i)).Ge.FatorDeExplosao)then
            go to 30
        endif
    enddo
    30 continue

    psi=psi**2 ! Temos agora a densidade de probabilidade.
    Probabilidade = 0.0d0 ! Vamos fazer a normalização.
    do i=1,N-1
        if(i*dx.Lt.1.0d0)then ! A normalização é feita apenas dentro do poço.
            Probabilidade=Probabilidade+Integral(psi(i+1),psi(i),psi(i-1),dx/2.0d0)
        endif
    enddo
    Probabilidade=2.0d0*Probabilidade
    psi=psi/Probabilidade ! Temos agora a função de onda normalizada.
    FatorDeExplosao = FatorDeExplosao**2/Probabilidade ! Normalizando o fator de explosão

    ! Saindo com os dados para o arquivo
    open(unit = 1, file = "pc1.txt")
    open(unit = 2, file = "pc2.txt", access = 'append')
    do i=-N, N
        if(Psi(i).Lt.FatorDeExplosao)then
            write(1,*) i*dx, psi(i)
        endif
    enddo

    write(2,100) "Energia = ",Energia,"Pontos utilizados = ",N
    100 format(a,f10.6,a25,i10)
    close(1)
    close(2)

    contains
    ! Função potencial.
    ! Inicialmente apenas a caixa rígida será implementada.
    function Potencial(Posicao)
        real * 8 :: Potencial, Posicao
        !if((Posicao.Lt.1.0d0).And.(Posicao.Gt.-1.0d0))then ! -1 < x < 1
            Potencial = 0.0d0
        !else
            Potencial = 0,5_r8*Posicao**2.0_r8!1.0E+12
        !endif
    end function Potencial

    ! Função de integração numérica da equação diferencial.
    ! Vamos utilizar o método de Numerov para a integração numérica da EDO
    function EvoluiPsi(Yatual, Yantes, Posicao, Energia, Incremento, Direcao)
        real * 8 :: EvoluiPsi, Incremento, Posicao
        real * 8 :: Yatual, Yantes, Energia
        real * 8 :: gN, gNmenos1, gNmais1
        integer * 4 :: Direcao

        ! Vamos calcular os fatores multip´licativos do método de Numerov a parte.
        gNmais1 = 1.0d0 + G(Energia, Posicao+Incremento)*Incremento**2/12.0d0
        gN = 1.0d0 - 5.0d0*G(Energia, Posicao)*Incremento**2/12.0d0
        gNmenos1 = 1.0d0 + G(Energia, Posicao - Incremento)*Incremento**2/12.0d0

        ! Fazendo o cálculo do valor posterior
        if(Direcao.Eq.0)then    ! integração da direita para esquerda
            EvoluiPsi = (2.0d0*gN*Yatual - Gnmenos1*Yantes) / gNmais1
        else                    ! Integração da esquerda para direita
            EvoluiPsi = (2.0d0*gN*Yatual - GNmais1*Yantes) / gNmenos1
        end if
    end function EvoluiPsi

    ! Função G do numerov
    function G(Energia, Posicao)
        real * 8 :: G, Energia, Posicao
        G = 2.0d0*(Energia - Potencial(Posicao))
    end function G

    ! Integral pela regra de Simpson
    function Integral(fp1,f0,fm1,h)
        real*8::fp1,f0,fm1,h
        real*8::integral
        Integral=(h/3.0d0)*(fp1+4.0d0*f0+fm1)
    end function Integral
end program PC