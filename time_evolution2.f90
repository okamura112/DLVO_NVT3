program main
  use functions
  use fcc_module_v2
  implicit none

  integer, parameter :: N=1000, Nstep=50000, ndiv=451, iteration=100
  !integer, parameter :: N=1000, Nstep=4000000, ndiv=451, iteration=1
  real(wp), parameter :: Q=1, deltat=0.002_wp, g = 3*N, rmax=4.5_wp, r_cutoff=2.5_wp
  real(wp) :: KBT0, momentum(3, N), momentum_av(3), V, s, eta, PV, position(3, N), Ftotal(3, N),&
  & zeta, zeta0, time, L, Pex, P0, Tinstant, Pinstant, deltar, rho, integrand(4, ndiv),r_array(ndiv),&
  & g_array(ndiv), error_integral, theta(4), dtheta(4), r, gtgt(ndiv), alpha, L0, M, position0(3, N), momentum0(3, N), v_cutoff
  integer(i4b) :: i, imin, j
  integer seedsize
  integer,allocatable:: seed(:)
  logical, parameter :: readfile = .true.
  M=0.001
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  ! seed = [  1632087667, -1853336134,   564385543,   814436085,  -340992037,   361556921,  -209566829,  1696273288]!test用
  ! call random_seed(put=seed)!test用
  open(11, file = 'position, momentum, M=0.001, Nstep=50000, iteration=100, alpha=0.001,&
  &r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(12, file = 'L, zeta, eta, s, PV, V, M=0.001, Nstep=50000, iteration=100, alpha=0.001,&
  &r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(13, file = 'time, Pinstant, Tinstant, V, PV, eta, M=0.001, Nstep=50000, iteration=100, alpha=0.001,&
  &r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(14, file = 'hamiltonian, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(15, file = 'g_of_r, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(16, file = 'vpair, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(17, file = 'nabra_theta_u, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(18, file = 'bcc_theta24.txt')
  open(19, file = 'integrate_gr')
  open(20, file = 'theta, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(21, file = 'randp_NVT_KBT1.5_N1000.txt', action = 'read')
  !open(21, file = 'random_position1000.txt', action = 'read')
  !open(22, file = 'error_integral, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  !open(23, file = 'error_values.dat')
  open(25, file = 'r, integrand, M=0.001, Nstep=50000, iteration=100, alpha=0.001, r_cutoff=2.5, theta42(16)2, KBT0=1.5,1(3)2.dat')
  open(101, file = 'setting.dat')
  write(101, *) seed
  imin = 102
  alpha = 0.001_wp

  !theta = [0.867*0.95, 177.2*0.95, 0.145*0.95, 0.578*0.95]!rcut=4.5でBCCが得られる値からずらす
  theta = [0.867*0.9, 177.2*1.1, 0.168*1.1, 0.591*1.1]!rcut=2.5でBCCが得られる値からずらす
  !theta = [0.824, 168.3, 0.150, 0.556]
  !theta = [0.910, 168.3, 0.150, 0.556]
  !theta = [0.867, 177.2, 0.145, 0.550]!rcut=4.5でBCCが得られる値
  !theta = [0.856, 161.28, 0.139, 0.578]!rcut=4.5でBCCが得られる値
  !theta = [0.889, 161.28, 0.142, 0.563]!rcut=4.5でBCCが得られる値

  v_cutoff = vpair(r_cutoff, theta)
  deltar = rmax/(ndiv-1)
  zeta = 1
  eta = 0
  s = 1
  PV = 0
  write(20, '(i8, 4g26.16)') 0, theta
  do i = 1, 200
    r = theta(4)+r_cutoff*i/200
    write(16, *) r, vpair_cut(r, theta, r_cutoff, v_cutoff)
    write(17, *) r, nabla_theta_u(r, theta, r_cutoff)
  end do
  write(16, '(/)')
  if (readfile) then
    call read_position_momentum(N, position0, momentum0, L0, zeta0, eta , s, PV, V, 21)
    L=L0
    position = position0
    call random_number(momentum0)
    momentum_av = sum(momentum0, dim = 2)/N
    do i = 1, N
      momentum0(:, i) = momentum0(:, i) - momentum_av
    end do
    momentum = momentum0
    zeta = zeta0
  else
    call random_number(momentum)
    momentum_av = sum(momentum, dim = 2)/N
    do i = 1, N
      momentum(:, i) = momentum(:, i) - momentum_av
    end do
    position = get_fcc_coordinates(N)*1.2_wp
    rho = 1.4_wp
    L0 = (N/rho)**(1.0_wp/3)
    L = L0
    call get_Ftotal_Pexcess(L,N,position, Ftotal, Pex, theta, r_cutoff)
    ! do i = 1, N
    !   write(999, *) i, norm2(Ftotal(:, i))
    ! end do
    zeta0 = 1
    zeta = zeta0
    V=L**3
  end if
  !試しに0.95倍
  ! L0 = 0.95*L0
  ! L = L0
  ! position0 = 0.95*position0
  ! position = position0
  !ここまで
  call write_L_zeta_eta_s_PV_V(L, zeta, eta , s, PV, V, 12)

  call write_position_momentum(N, position, momentum, L, 11)
  call read_gtgt(ndiv, gtgt, 18)
  do i = 1, ndiv
    if(gtgt(i)>0.01) then
      imin = i
      exit
    end if
  end do
  KBT0=1.0
  g_array = g_of_r(position, deltar, n, ndiv, L, V)
  do i = 1, ndiv
    write(15,*) (i-0.5)*deltar, g_array(i)
  end do
  write(15,'(/)')
  ! do i = 1, ndiv
  !   write(23, *) i, gtgt(i), g_array(i), (g_array(i) - gtgt(i))**2
  ! end do
  ! write(23,'(/)')
  do j = 1, iteration
    ! if((j==2).or.(j==3)) then
    !    M = 10*M
    ! end if
    position = position0
    momentum = momentum0
    ! L=L0
    ! V=L**3
    zeta = zeta0
    do i = 1, Nstep
      call time_evolutionNVT3(momentum, zeta, deltat)
      call get_Ftotal_Pexcess(L,N,position, Ftotal, Pex, theta, r_cutoff)
      call time_evolutionNVT2(momentum, Ftotal, deltat)
      call time_evolutionNVT1(position, momentum, zeta, deltat, Q, g, KBT0, L)
      call get_Ftotal_Pexcess(L,N,position, Ftotal, Pex, theta, r_cutoff)
      call time_evolutionNVT2(momentum, Ftotal, deltat)
      call time_evolutionNVT3(momentum, zeta, deltat)

      if ((i==1).or.(mod(i, Nstep/200)==0)) then
        Tinstant = instantaneous_temprature(momentum, g)
        Pinstant = N*Tinstant/V + Pex
        time = i * deltat
        write(13,'(6g24.14)') time, Pinstant, Tinstant, V, PV, eta
        write(14,'(3g18.8)') time, hamiltonian(momentum, N, position, L, Q, zeta, g, KBT0, s, theta, r_cutoff, v_cutoff),&
        &total_potential(N, position, L, theta, r_cutoff, v_cutoff)
      end if
    end do
    write(13, '(/)')
    write(14, '(/)')

    g_array = g_of_r(position, deltar, n, ndiv, L, V)
    if((j == 1).or.(mod(j, max(1, int(iteration/5)))==0))then
      call write_L_zeta_eta_s_PV_V(L, zeta, eta , s, PV, V, 12)
      call write_position_momentum(N, position, momentum, L, 11)
    end if
    if((j == 1).or.(mod(j, max(1, int(iteration/20)))==0))then
      do i = 1, ndiv
        write(15,*) (i-0.5)*deltar, g_array(i)
      end do
      write(15,'(/)')


    end if
    do i = 1, ndiv
      write(23, *) i, gtgt(i), g_array(i), (g_array(i) - gtgt(i))**2
    end do
    write(23,'(/)')
    !gtgt = 0
    !write(19, *) integrate_gr(r, ndiv, deltar, gtgt, position, N, L, V, theta)
!write(997,*) deltar, r_array(i), integrand(3, ndiv), gtgt(ndiv), nabla_theta_u(r_array(i), theta)
    !write(*, *) 'test', deltar, N, L, V, rmax
    ! do i = 1, ndiv
    !   if(g_array(i)>0.01) then
    !     imin = min(i, imin_tgt)
    !     exit
    !   end if
    ! end do
    dtheta = integrate_gr(ndiv, deltar, gtgt, g_array, theta, r_cutoff, imin)
    !write(*, *)  dtheta
    theta = theta + alpha*dtheta
    v_cutoff = vpair(r_cutoff, theta)
    if((j == 1).or.(mod(j, max(1, int(iteration/100)))==0))then
      error_integral = calculate_error_integral(ndiv, deltar, gtgt, g_array, r_cutoff)
      !write(22, '(g24.16)') error_integral
      if (error_integral < 0.05_wp) then
        call write_L_zeta_eta_s_PV_V(L, zeta, eta , s, PV, V, 12)
        call write_position_momentum(N, position, momentum, L, 11)
        exit
      end if
      write(20, '(i8, 9g26.16)') j, theta, dtheta, error_integral
      do i = imin, ndiv
        r_array(i) = (i - 1) * deltar
        integrand(:, i) = r_array(i)**2 * (g_array(i) - gtgt(i)) * nabla_theta_u(r_array(i), theta, r_cutoff)
        write(25, '(5g24.16)') r_array(i), integrand(:, i)
      end do
      write(25, '(/)')

      do i = 1, 200
        r = theta(4)+5._wp*i/200
        write(16, *) r, vpair_cut(r, theta, r_cutoff, v_cutoff)
        write(17, *) r, nabla_theta_u(r, theta, r_cutoff)
      end do
      write(16, '(/)')
    end if
  end do

end program main
