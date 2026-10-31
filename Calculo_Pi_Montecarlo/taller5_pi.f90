! Taller 5 -- Convergencia del estimador de pi por acierto y error.
! Salidas: convergencia.dat (N, media, desviación) y pendiente.dat (pendiente m).
program taller5_pi
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: M = 20                  ! corridas por cada N
  integer :: k, j, n_sem, u_dat, n_puntos
  integer(8) :: i, n_total, n_aciertos
  integer, allocatable :: semilla(:)
  real(dp) :: u, v, x, y, pi_est, suma, suma2, media, desv
  real(dp) :: X_k, Y_k, sx, sy, sxx, sxy, pendiente

  call random_seed(size=n_sem)
  allocate(semilla(n_sem))

  sx = 0.0_dp; sy = 0.0_dp; sxx = 0.0_dp; sxy = 0.0_dp; n_puntos = 0   ! sumas del ajuste
  open(newunit=u_dat,file='convergencia.dat',status='replace')
  do k = 2, 7                                   ! N = 10^2, ..., 10^7
     n_total = 10_8**k
     suma = 0.0_dp; suma2 = 0.0_dp              ! sumas sobre las M corridas
     do j = 1, M
        semilla = 1000*k + j                    ! una semilla distinta por corrida
        call random_seed(put=semilla)
        pi_est = -1.0_dp                        ! valor provisional

        ! TODO A1: obtenga pi_est con n_total puntos uniformes en el cuadrado (sección 2).
        ! Inicializamos el contador de aciertos dentro del círculo unitario
        n_aciertos = 0_8
        do i = 1_8, n_total
           ! Generamos dos números aleatorios uniformes en [0,1)
           call random_number(u)
           call random_number(v)
           ! Mapeamos el intervalo [0,1) a [-1,1) para las coordenadas (x,y)
           x = 2.0_dp*u - 1.0_dp
           y = 2.0_dp*v - 1.0_dp
           ! Verificamos la condición del círculo x^2 + y^2 < 1 (área = pi)
           if (x*x + y*y < 1.0_dp) then
              n_aciertos = n_aciertos + 1_8
           end if
        end do
        ! Estimación de pi = 4 * (puntos dentro) / (puntos totales)
        ! Usamos real(..., dp) para evitar truncamiento por división entera
        pi_est = 4.0_dp * real(n_aciertos, dp) / real(n_total, dp)

        ! TODO A2: acumule lo necessary para la media y la desviación de las M corridas.
        ! Acumulamos la suma lineal y la suma de cuadrados de pi_est
        suma = suma + pi_est
        suma2 = suma2 + pi_est*pi_est
     end do

     media = -1.0_dp; desv = -1.0_dp            ! valores provisionales
     ! TODO A3: calcule la media y la desviación estándar muestral de las M estimaciones.
     ! Calculamos la media muestral de las M corridas
     media = suma / real(M, dp)
     ! Calculamos la desviación estándar muestral s con el factor de corrección M/(M-1)
     desv = sqrt((real(M, dp)/(real(M, dp) - 1.0_dp)) * ((suma2 / real(M, dp)) - media*media))

     write(u_dat,'(i10,2es25.16e3)') n_total,media,desv
     print '(i10,2f12.6)', n_total,media,desv

     ! TODO A4: acumule lo necesario para ajustar log10(desv) contra log10(N).
     ! Variables para la regresión lineal en escala logarítmica: X_k = log10(N), Y_k = log10(s)
     X_k = log10(real(n_total, dp))
     Y_k = log10(desv)
     ! Acumulamos las sumatorios sum(X), sum(Y), sum(X^2) y sum(X*Y)
     sx = sx + X_k
     sy = sy + Y_k
     sxx = sxx + X_k*X_k
     sxy = sxy + X_k*Y_k
     n_puntos = n_puntos + 1
  end do
  close(u_dat)

  pendiente = -1.0_dp                           ! valor provisional
  ! TODO A5: calcule la pendiente del ajuste por mínimos cuadrados.
  ! Fórmula de la pendiente por mínimos cuadrados: m = (n*sxy - sx*sy) / (n*sxx - sx^2)
  pendiente = (real(n_puntos, dp)*sxy - sx*sy) / (real(n_puntos, dp)*sxx - sx*sx)

  open(newunit=u_dat,file='pendiente.dat',status='replace')
  write(u_dat,'(es25.16e3)') pendiente
  close(u_dat)
  print '(a,f9.4)', 'Pendiente de log10(s) contra log10(N): ', pendiente
end program taller5_pi
