! ==============================================================================
! Módulo: modulo_romberg
! Asignatura: Física Computacional (106018C) - Taller 8 (Parte C)
! Descripción: Implementa la integración numérica mediante el método de
!              Romberg (regla del trapecio compuesto + extrapolación de
!              Richardson) en doble precisión, junto con los integrandos de
!              la ley de Stefan-Boltzmann con y sin cambio de variable.
! Autor: Erick Alexander Muñoz Orozco
! Fecha: 16/11/2026
! ==============================================================================
module modulo_romberg
  implicit none

  ! ----------------------------------------------------------------------------
  ! Parámetros globales de precisión y control de niveles de Romberg
  ! ----------------------------------------------------------------------------
  integer, parameter :: dp    = kind(1.0d0) ! Doble precisión (64 bits)
  integer, parameter :: k_max = 20            ! Nivel máximo de refinamiento (2^20 + 1 pts)

contains

  ! ----------------------------------------------------------------------------
  ! Función: integrando_sin_cambio
  ! Entrada:  x (variable adimensional en [0, inf)) [real(dp)]
  ! Salida:   f(x) = x^3 / (exp(x) - 1) [real(dp)]
  ! Descripción: Evalúa el integrando de Stefan-Boltzmann en su variable
  !              original, garantizando el límite exacto f(0) = 0 y evitando
  !              desbordamiento exponencial para x muy grande.
  ! ----------------------------------------------------------------------------
  pure function integrando_sin_cambio(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res

    if (x <= 0.0_dp .or. x > 700.0_dp) then
      res = 0.0_dp
    else
      res = (x**3) / (exp(x) - 1.0_dp)
    end if
  end function integrando_sin_cambio

  ! ----------------------------------------------------------------------------
  ! Función: integrando_con_cambio
  ! Entrada:  t (variable transformada en el intervalo finito [0, 1]) [real(dp)]
  ! Salida:   g(t) = [x(t)^3 / (exp(x(t)) - 1)] * [1 / (1 - t)^2] [real(dp)]
  ! Descripción: Mapea x = t / (1 - t) sobre [0, 1] e impone los límites
  !              analíticos g(0) = 0 y g(1) = 0 en los extremos del intervalo.
  ! ----------------------------------------------------------------------------
  pure function integrando_con_cambio(t) result(res)
    real(dp), intent(in) :: t
    real(dp)             :: res
    real(dp)             :: x, one_minus_t

    ! Límites exactos en los extremos t = 0 y t = 1
    if (t <= 0.0_dp .or. t >= 1.0_dp) then
      res = 0.0_dp
      return
    end if

    one_minus_t = 1.0_dp - t
    x           = t / one_minus_t

    ! Evita overflow en exp(x) cuando t es extremadamente cercano a 1
    if (x > 700.0_dp) then
      res = 0.0_dp
    else
      res = ((x**3) / (exp(x) - 1.0_dp)) / (one_minus_t**2)
    end if
  end function integrando_con_cambio

  ! ----------------------------------------------------------------------------
  ! Función: evaluar_integrando
  ! Entradas: u (punto de evaluación) [real(dp)]
  !           modo (1 = con cambio de variable en [0,1];
  !                 2 = sin cambio de variable en [0,100]) [integer]
  ! Salida:   valor del integrando seleccionado en u [real(dp)]
  ! ----------------------------------------------------------------------------
  pure function evaluar_integrando(u, modo) result(res)
    real(dp), intent(in) :: u
    integer, intent(in)  :: modo
    real(dp)             :: res

    if (modo == 1) then
      res = integrando_con_cambio(u)
    else
      res = integrando_sin_cambio(u)
    end if
  end function evaluar_integrando

  ! ----------------------------------------------------------------------------
  ! Subrutina: integrar_romberg
  ! Entradas: a, b (límites de integración) [real(dp)]
  !           tol (tolerancia de convergencia) [real(dp)]
  !           modo (selector de integrando: 1 con cambio, 2 sin cambio) [integer]
  ! Salidas:  integral (valor estimado R(k,k)) [real(dp)]
  !           k_final (nivel final alcanzado) [integer]
  !           n_eval (número total de evaluaciones de función: 2^k_final + 1) [integer]
  ! ----------------------------------------------------------------------------
  subroutine integrar_romberg(a, b, tol, modo, integral, k_final, n_eval)
    real(dp), intent(in)  :: a, b, tol
    integer, intent(in)   :: modo
    real(dp), intent(out) :: integral
    integer, intent(out)  :: k_final, n_eval

    ! Tabla de Romberg y variables auxiliares
    real(dp) :: R(0:k_max, 0:k_max)
    real(dp) :: h, suma, x_i, factor4
    integer  :: k, j, i, n_sub

    R = 0.0_dp

    ! Nivel k = 0: Regla del trapecio simple con 1 subintervalo (2 evaluaciones)
    h       = b - a
    R(0, 0) = 0.5_dp * h * (evaluar_integrando(a, modo) + evaluar_integrando(b, modo))

    ! Refinamiento sucesivo de niveles k = 1, ..., k_max
    do k = 1, k_max
      h     = 0.5_dp * h
      n_sub = 2**(k - 1)
      suma  = 0.0_dp

      ! Suma únicamente los puntos nuevos impares para reutilizar evaluaciones previas
      do i = 1, n_sub
        x_i  = a + real(2*i - 1, dp) * h
        suma = suma + evaluar_integrando(x_i, modo)
      end do

      ! Trapecio compuesto en el nivel k
      R(k, 0) = 0.5_dp * R(k - 1, 0) + h * suma

      ! Extrapolación de Richardson para eliminar términos de error O(h^(2j))
      factor4 = 1.0_dp
      do j = 1, k
        factor4 = 4.0_dp * factor4
        R(k, j) = R(k, j - 1) + (R(k, j - 1) - R(k - 1, j - 1)) / (factor4 - 1.0_dp)
      end do

      ! Criterio de parada por convergencia diagonal
      if (abs(R(k, k) - R(k - 1, k - 1)) < tol) then
        k_final  = k
        n_eval   = 2**k + 1
        integral = R(k, k)
        return
      end if
    end do

    ! Si alcanza el límite máximo de niveles sin satisfacer tol
    k_final  = k_max
    n_eval   = 2**k_max + 1
    integral = R(k_max, k_max)
    print *, "Advertencia (Romberg): Se alcanzo el nivel maximo k_max sin converger."
  end subroutine integrar_romberg

end module modulo_romberg
