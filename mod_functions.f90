module functions
  use constants
  use mod_periodic_boundary
  implicit none
  
  contains
  
  function fpair(r_vec, theta, r_cutoff)
    real(wp)::r_vec(3), fpair(3), theta(4), r, r_cutoff
    r = norm2(r_vec)
    if(r>r_cutoff) then
      fpair = 0
    else
      fpair = -(theta(1)/theta(4)/(r/theta(4)-1)**2 &
      &- theta(2)/(theta(3)*theta(4))*exp(-(r/theta(4)-1)/theta(3)))*(r_vec/r)
    end if
  end function fpair
  
  function function_Ftotal(L, N, position, theta, r_cutoff) result(Ftotal)
    integer(i4b)::i,j,N
    real(wp)::L,Ftotal(3,N),relative_r(3),ftmp(3),position(3,N), theta(4), r_cutoff
    Ftotal=0
    do i = 1,N-1
      do j = i+1,N
        relative_r=relative_vector(L, position(:, i), position(:, j))
        ftmp=fpair(relative_r, theta, r_cutoff)
        Ftotal(:,j)=Ftotal(:,j)+ftmp
        Ftotal(:,i)=Ftotal(:,i)-ftmp
      end do
    end do
  end function function_Ftotal
  
  subroutine get_Ftotal_Pexcess(L,N,position, Ftotal, Pex, theta, r_cutoff)
    integer(i4b)::i,j,N
    real(wp)::L,Ftotal(3,N),relative_r(3),ftmp(3),position(3,N), Pex, theta(4), r_cutoff
    Ftotal=0
    Pex=0
    do i = 1,N-1
      do j = i+1,N
        relative_r=relative_vector(L, position(:, i), position(:, j))
        ftmp=fpair(relative_r, theta, r_cutoff)
        Ftotal(:,j)=Ftotal(:,j)+ftmp
        Ftotal(:,i)=Ftotal(:,i)-ftmp
        Pex=Pex+dot_product(relative_r, ftmp)
      end do
    end do
    Pex=Pex/(3*L**3)
  end subroutine get_Ftotal_Pexcess
  
  function instantaneous_temprature(momentum, g) result(KBT)
    real(wp) :: momentum(:, :), g, KBT
    KBT = sum(momentum**2)/g
  end function instantaneous_temprature
  
  !(9.54)
  subroutine time_evolutionNVT3(momentum, zeta, deltat)
    real(wp) :: momentum(:, :), zeta, deltat
    momentum=momentum*exp(-zeta*deltat/2)
  end subroutine time_evolutionNVT3
  
   !(9.55)
  subroutine time_evolutionNVT2(momentum, Ftotal, deltat)
    real(wp) :: momentum(:, :), Ftotal(:, :), deltat
    momentum=momentum+Ftotal*deltat/2
  end subroutine time_evolutionNVT2
  
   !(9.56)
  subroutine time_evolutionNVT1(position, momentum, zeta, deltat, Q, g, KBT0, L)
    real(wp) :: position(:, :), momentum(:, :), zeta, deltat, Q, g, KBT0, L
    position=position+momentum*deltat
    where (position < 0)
      position=position+L
    else where (position >= L)
      position=position-L
    end where
    zeta = zeta+(1/Q)*(sum(momentum**2)-g*kBT0)*deltat
  end subroutine time_evolutionNVT1
  
  
  subroutine write_position_momentum(N, position, momentum_tilde, L, file_number)
    integer file_number, N, i
    !real(wp), intent(in) :: position(3, N), momentum_tilde(3, N), L
    real(wp) :: position(3, N), momentum_tilde(3, N), L
    do i = 1, N
      write(file_number, '(6f24.16)') position(:,i), momentum_tilde(:,i)/L
    end do
    write(file_number, '(/)')
  end subroutine write_position_momentum
  
  subroutine write_L_zeta_eta_s_PV_V(L, zeta, eta , s, PV, V, file_number)
    integer file_number
    real(wp) :: L, zeta, eta , s, PV, V
    write(file_number, '(6f24.16)') L, zeta, eta , s, PV, V
  end subroutine write_L_zeta_eta_s_PV_V
  
  
  subroutine read_position_momentum(N, position, momentum_tilde, L, zeta, eta , s, PV, V, file_number)
    integer file_number, N, i
    real(wp) :: position(3, N), momentum_tilde(3, N), momentum(3, N), L, zeta, eta , s, PV, V
    do i = 1, N
      read(file_number, *) position(:,i), momentum(:,i)
    end do
    read(file_number, *) L, zeta, eta , s, PV, V
    momentum_tilde = L*momentum
  end subroutine read_position_momentum
  
  function vpair(r, theta)
    real(wp)::r, vpair, theta(4)
    vpair = -theta(1)/(r/theta(4)-1) + theta(2)*exp(-(r/theta(4)-1)/theta(3))
  end function vpair
  
  function vpair_cut(r, theta, r_cutoff, v_cutoff)
    real(wp)::r, vpair_cut, theta(4), r_cutoff, v_cutoff
    if(r>r_cutoff) then
      vpair_cut = 0
    else
      vpair_cut = -theta(1)/(r/theta(4)-1) + theta(2)*exp(-(r/theta(4)-1)/theta(3)) - v_cutoff
    end if
  end function vpair_cut
  
  function total_potential(N, position, L, theta, r_cutoff, v_cutoff)
    real(wp) :: position(3, N), L, total_potential, relative_r(3), r, theta(4), r_cutoff, v_cutoff
    integer N, i1, i2
    total_potential = 0
    do i2 = 1, N-1
      do i1 = i2+1, N
        relative_r=relative_vector(L, position(:, i1), position(:, i2))
        r = norm2(relative_r)
        total_potential = total_potential + vpair_cut(r, theta, r_cutoff, v_cutoff)
      end do
    end do
  end function total_potential
  
  function hamiltonian(momentum, N, position, L, Q, zeta, g, KBT0, s, theta, r_cutoff, v_cutoff)
    real(wp) :: hamiltonian, momentum(3, N), position(3, N), L, Q, zeta, g, KBT0, s, theta(4), r_cutoff, v_cutoff
    integer N
    hamiltonian = sum(momentum**2)/2+total_potential(N, position, L, theta, r_cutoff, v_cutoff)&
    &+(Q*zeta**2)/2+g*KBT0*log(s)
  end function hamiltonian
  
  function g_of_r(position, deltar, N, ndiv, L, V)
    real(wp) :: g_of_r(ndiv), ng(ndiv), r, deltar, relative_r(3), position(3, N), L, V
    integer N, i1, i2, ndiv, ir
    ng = 0
    do i1 = 1, N-1
      do i2 = i1+1, N
        relative_r=relative_vector(L, position(:, i1), position(:, i2))
        r = norm2(relative_r)
        ir = int(r/deltar)+1
        if(ir<=ndiv) ng(ir) = ng(ir)+1
      end do
    end do
    do i1 = 1, ndiv
      r = i1*deltar
      g_of_r(i1) = ng(i1)/(((4*pi)*(r**2+r*deltar+(deltar**2/3))*deltar)*(0.5*N)*(N/V))
    end do
  end function g_of_r
  
  function nabla_theta_u(r, theta, r_cutoff)
    real(wp) :: nabla_theta_u(4), r, theta(4), r_cutoff
    if(r>r_cutoff) then
      nabla_theta_u = 0
    else
      nabla_theta_u(1) = -1/(r/theta(4)-1)
      nabla_theta_u(2) = exp(-(r/theta(4)-1)/theta(3))
      nabla_theta_u(3) = theta(2)*(r/theta(4)-1)*exp(-(r/theta(4)-1)/theta(3))/theta(3)**2
      nabla_theta_u(4) = -r*theta(1)/(r-theta(4))**2 &
      &+ r*theta(2)/(theta(3)*theta(4)**2)*exp((theta(4)-r)/(theta(3)*theta(4)))
    end if
  end function nabla_theta_u
  
  
  function integrate(array, ndiv, deltar)
    real(wp) :: integrate, array(ndiv), deltar
    integer ndiv
    integrate = deltar*sum(array)
  end function integrate
  
  subroutine read_gtgt(ndiv, gtgt, file_number)
    integer :: file_number, ndiv, i
    real(wp) :: y, gtgt(ndiv)
    !open(18, file = 'bcc_tsujiData.txt')
    do i = 1, ndiv
      read(file_number, *) y, gtgt(i)
    end do
  end subroutine read_gtgt
  
  function integrate_gr(ndiv, deltar, gtgt, g_array, theta, r_cutoff , imin)
    integer :: ndiv, i, imin
    real(wp) :: integrate_gr(4), deltar, gtgt(ndiv), g_array(ndiv), theta(4), r_array(ndiv),&
    & integrand(4, ndiv), r_cutoff
    do i = imin, ndiv 
      r_array(i) = (i - 1) * deltar
       integrand(:, i) = r_array(i)**2 * (g_array(i) - gtgt(i)) * nabla_theta_u(r_array(i), theta, r_cutoff)
      end do 
    do i = 1, 4
      integrate_gr(i) = integrate(integrand(i, imin:ndiv), ndiv-imin+1, deltar)
    end do
  end function integrate_gr

  function calculate_error_integral(ndiv, deltar, gtgt, g_array, r_cutoff) result(error_integral)
    integer :: ndiv, i
    real(wp) :: deltar, gtgt(ndiv), g_array(ndiv), error_integral, r_cutoff
    real(wp) :: error_integrand(ndiv)
    do i = 1, ndiv
      if ((i-1)*deltar <= r_cutoff) then
        error_integrand(i) = (g_array(i) - gtgt(i))**2
      else
        error_integrand(i) = 0
      end if
    end do
    error_integral = integrate(error_integrand, ndiv, deltar)
  end function calculate_error_integral
  
  ! 与えられた g_in(1:ndiv_in) （第i成分は位置 r = i * delta_in における動径分布関数の値）に対して、
  ! 横軸方向に scale 倍された動径分布関数の値を保持する配列 g_out(1:ndiv_out) （第i成分は位置 r = i * delta_out
  ! における動径分布関数の値）を求めるサブルーチン
  subroutine scale_g_of_r(ndiv_in, deltar_in, g_in, scale, ndiv_out, deltar_out, g_out)
    integer, intent(in) :: ndiv_in, ndiv_out
    real(wp), intent(in) :: deltar_in, g_in(ndiv_in), scale, deltar_out
    real(wp), intent(out) :: g_out(ndiv_out)
    integer :: i1, i_r_in
    real(wp) :: lambda, r_in(ndiv_in), r_out
    do i1 = 1, ndiv_in
      r_in(i1) = i1 * deltar_in * scale
    end do
    do i1 = 1, ndiv_out
      r_out = i1 * deltar_out
      i_r_in = int(r_out / (deltar_in * scale)) + 1
      if ((i_r_in >= 1).and.(i_r_in <= ndiv_in)) then
        lambda = (r_out - r_in(i_r_in - 1)) / (r_in(i_r_in) - r_in(i_r_in - 1))
        g_out(i1) = (1 - lambda) * g_in(i_r_in) + lambda * g_in(i_r_in + 1)
      elseif (i_r_in == 0) then
        g_out(i1) = 0
      else
        g_out(i1:ndiv_out) = 1
        exit
      end if
    end do
  end subroutine scale_g_of_r
  
  end module functions