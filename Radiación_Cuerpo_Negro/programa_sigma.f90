! ==============================================================================
! Programa: programa_sigma
! Asignatura: Física Computacional (106018C) - Taller 8 (Parte C)
! Descripción: Calcula la integral adimensional de Planck I y la constante de
!              Stefan-Boltzmann sigma mediante integración de Romberg en dos
!              regímenes: (1) con cambio de variable x = t/(1 - t) en [0, 1] y
!              (2) directamente sobre [0, 100] sin cambio de variable. Exporta
!              los resultados a resultado_sigma.dat y los compara con pi^4/15.
! Autor: Erick Alexander Muñoz Orozco
! Fecha: 16/11/2026
! ==============================================================================
program programa_sigma
  use modulo_romberg
  implicit none

  ! ----------------------------------------------------------------------------
  ! Constantes físicas exactas del Sistema Internacional (SI) y tolerancia
  ! ----------------------------------------------------------------------------
  real(dp), parameter :: h   = 6.62607015e-34_dp         ! Constante de Planck [J s]
  real(dp), parameter :: c   = 299792458.0_dp            ! Rapidez de la luz [m/s]
  real(dp), parameter :: kB  = 1.380649e-23_dp           ! Constante de Boltzmann [J/K]
  real(dp), parameter :: pi  = acos(-1.0_dp)             ! Número pi en doble precisión
  real(dp), parameter :: tol = 1.0e-8_dp                 ! Tolerancia de integración

  ! ----------------------------------------------------------------------------
  ! Variables de resultados y unidades de archivo
  ! ----------------------------------------------------------------------------
  real(dp) :: pref_sigma, I_exacta, sigma_exacta
  real(dp) :: I_cambio, sigma_cambio
  real(dp) :: I_directa, sigma_directa
  integer  :: k_cambio, eval_cambio
  integer  :: k_directa, eval_directa
  integer  :: u_out

  ! Prefactor físico de la constante de Stefan-Boltzmann: (2 * pi * kB^4) / (h^3 * c^2)
  pref_sigma   = (2.0_dp * pi * (kB**4)) / ((h**3) * (c**2))
  I_exacta     = (pi**4) / 15.0_dp
  sigma_exacta = pref_sigma * I_exacta

  ! ----------------------------------------------------------------------------
  ! Bloque 1: Integración con cambio de variable en t in [0, 1] (modo = 1)
  ! ----------------------------------------------------------------------------
  call integrar_romberg(0.0_dp, 1.0_dp, tol, 1, I_cambio, k_cambio, eval_cambio)
  sigma_cambio = pref_sigma * I_cambio

  ! ----------------------------------------------------------------------------
  ! Bloque 2: Integración directa sin cambio de variable en x in [0, 100] (modo = 2)
  ! ----------------------------------------------------------------------------
  call integrar_romberg(0.0_dp, 100.0_dp, tol, 2, I_directa, k_directa, eval_directa)
  sigma_directa = pref_sigma * I_directa

  ! ----------------------------------------------------------------------------
  ! Bloque 3: Exportación de resultados a resultado_sigma.dat
  !           Columnas: k_final, evaluaciones, integral, sigma
  ! ----------------------------------------------------------------------------
  open(newunit=u_out, file="resultado_sigma.dat", status="replace", action="write")
  write(u_out, '(i6, i10, 2es24.16)') k_cambio,  eval_cambio,  I_cambio,  sigma_cambio
  write(u_out, '(i6, i10, 2es24.16)') k_directa, eval_directa, I_directa, sigma_directa
  close(u_out)

  ! ----------------------------------------------------------------------------
  ! Bloque 4: Comparación con el valor analítico exacto en consola
  ! ----------------------------------------------------------------------------
  print *, "========================================================================="
  print *, " RESULTADOS PARTE C - CONSTANTE DE STEFAN-BOLTZMANN (ROMBERG)"
  print *, "========================================================================="
  print '(a, es24.16)', " Valor exacto I = pi^4/15      = ", I_exacta
  print '(a, es24.16)', " Valor exacto sigma [W/m^2K^4] = ", sigma_exacta
  print *, "-------------------------------------------------------------------------"
  print '(a, i4, a, i8)', " [1] Con cambio [0, 1]   -> k = ", k_cambio,  ", eval = ", eval_cambio
  print '(a, es24.16, a, es12.4)', "     I_num = ", I_cambio,  " | Error abs = ", abs(I_cambio - I_exacta)
  print '(a, es24.16)',            "     sigma = ", sigma_cambio
  print *, "-------------------------------------------------------------------------"
  print '(a, i4, a, i8)', " [2] Sin cambio [0, 100] -> k = ", k_directa, ", eval = ", eval_directa
  print '(a, es24.16, a, es12.4)', "     I_num = ", I_directa, " | Error abs = ", abs(I_directa - I_exacta)
  print '(a, es24.16)',            "     sigma = ", sigma_directa
  print *, "========================================================================="

end program programa_sigma
