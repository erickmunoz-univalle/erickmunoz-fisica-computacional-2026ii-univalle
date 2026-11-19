! ==============================================================================
! Módulo: modulo_ceros
! Asignatura: Física Computacional (106018C) - Taller 8
! Descripción: Implementa la función trascendente asociada a la ley de
!              desplazamiento de Wien, su primera derivada y tres métodos
!              numéricos de búsqueda de raíces (Bisección, Newton-Raphson y
!              Secante) en doble precisión.
! Autor: Erick Alexander Muñoz Orozco
! Fecha: 14/11/2026
! ==============================================================================
module modulo_ceros
  implicit none

  ! ----------------------------------------------------------------------------
  ! Parámetros globales de precisión y control iterativo
  ! ----------------------------------------------------------------------------
  integer, parameter :: dp       = kind(1.0d0) ! Doble precisión (64 bits)
  integer, parameter :: max_iter = 100           ! Límite máximo de iteraciones

contains

  ! ----------------------------------------------------------------------------
  ! Función: f_wien
  ! Entradas: x (variable adimensional hc / (lambda * k_B * T)) [real(dp)]
  !           n (exponente del espectro: n = 5 en lambda, n = 3 en nu) [integer]
  ! Salida:   valor de F(x) = x - n * (1 - exp(-x)) [real(dp)]
  ! ----------------------------------------------------------------------------
  pure function f_wien(x, n) result(res)
    real(dp), intent(in) :: x
    integer, intent(in)  :: n
    real(dp)             :: res

    res = x - real(n, dp) * (1.0_dp - exp(-x))
  end function f_wien

  ! ----------------------------------------------------------------------------
  ! Función: df_wien
  ! Entradas: x (variable adimensional) [real(dp)]
  !           n (exponente del espectro) [integer]
  ! Salida:   primera derivada F'(x) = 1 - n * exp(-x) [real(dp)]
  ! ----------------------------------------------------------------------------
  pure function df_wien(x, n) result(res)
    real(dp), intent(in) :: x
    integer, intent(in)  :: n
    real(dp)             :: res

    res = 1.0_dp - real(n, dp) * exp(-x)
  end function df_wien

  ! ----------------------------------------------------------------------------
  ! Subrutina: biseccion
  ! Entradas: a_in, b_in (extremos del intervalo inicial [a, b]) [real(dp)]
  !           n (parámetro de la función de Wien) [integer]
  !           tol (tolerancia de convergencia epsilon) [real(dp)]
  ! Salidas:  raiz (aproximación al cero de F(x)) [real(dp)]
  !           iter (número de iteraciones realizadas) [integer]
  !           historial (vector con las estimaciones c_k de cada paso) [real(dp)]
  ! Descripción: Encuentra una raíz en [a, b] mediante bisección sucesiva.
  !              Verifica previamente el cambio de signo (Teorema de Bolzano)
  !              y detiene la búsqueda cuando b - c <= tol o tras 100 pasos.
  ! ----------------------------------------------------------------------------
  subroutine biseccion(a_in, b_in, n, tol, raiz, iter, historial)
    real(dp), intent(in)  :: a_in, b_in, tol
    integer, intent(in)   :: n
    real(dp), intent(out) :: raiz
    integer, intent(out)  :: iter
    real(dp), intent(out) :: historial(max_iter)

    ! Variables locales del intervalo y evaluaciones
    real(dp) :: a, b, c, fa, fb, fc

    ! Inicialización del estado
    a         = a_in
    b         = b_in
    fa        = f_wien(a, n)
    fb        = f_wien(b, n)
    iter      = 0
    historial = 0.0_dp

    ! Comprobación de cambio de signo antes de iterar
    if (fa * fb >= 0.0_dp) then
      print *, "Error (Biseccion): No hay cambio de signo en el intervalo inicial."
      raiz = 0.0_dp
      return
    end if

    ! Ciclo iterativo de reducción del intervalo
    do while (iter < max_iter)
      iter            = iter + 1
      c               = 0.5_dp * (a + b)
      historial(iter) = c
      fc              = f_wien(c, n)

      ! Criterio de parada por semi-ancho del intervalo o raíz exacta
      if (b - c <= tol .or. fc == 0.0_dp) then
        raiz = c
        return
      end if

      ! Selección del subintervalo que conserva el cambio de signo
      if (fa * fc < 0.0_dp) then
        b  = c
        fb = fc
      else
        a  = c
        fa = fc
      end if
    end do

    raiz = c
    print *, "Advertencia (Biseccion): Se alcanzo el maximo de 100 iteraciones sin converger."
  end subroutine biseccion

  ! ----------------------------------------------------------------------------
  ! Subrutina: newton
  ! Entradas: x0 (estimación inicial) [real(dp)]
  !           n (parámetro de la función de Wien) [integer]
  !           tol (tolerancia de convergencia epsilon) [real(dp)]
  ! Salidas:  raiz (aproximación al cero de F(x)) [real(dp)]
  !           iter (número de iteraciones realizadas) [integer]
  !           historial (vector con las estimaciones x_k de cada paso) [real(dp)]
  ! Descripción: Aplica el método de Newton-Raphson x_{k+1} = x_k - F(x_k)/F'(x_k).
  !              Se detiene con mensaje si F'(x_k) es nula o si |x_{k+1} - x_k| < tol.
  ! ----------------------------------------------------------------------------
  subroutine newton(x0, n, tol, raiz, iter, historial)
    real(dp), intent(in)  :: x0, tol
    integer, intent(in)   :: n
    real(dp), intent(out) :: raiz
    integer, intent(out)  :: iter
    real(dp), intent(out) :: historial(max_iter)

    ! Variables locales de iteración y evaluación
    real(dp) :: x_curr, x_next, fx, dfx

    ! Inicialización del estado
    x_curr    = x0
    iter      = 0
    historial = 0.0_dp

    ! Ciclo iterativo de aproximación tangencial
    do while (iter < max_iter)
      fx  = f_wien(x_curr, n)
      dfx = df_wien(x_curr, n)

      ! Protección contra división por cero (derivada nula)
      if (abs(dfx) <= tiny(1.0_dp)) then
        print *, "Error (Newton): Derivada nula. Deteniendo ejecucion."
        raiz = x_curr
        return
      end if

      iter            = iter + 1
      x_next          = x_curr - fx / dfx
      historial(iter) = x_next

      ! Criterio de parada por diferencia entre iteraciones consecutivas
      if (abs(x_next - x_curr) < tol) then
        raiz = x_next
        return
      end if

      x_curr = x_next
    end do

    raiz = x_curr
    print *, "Advertencia (Newton): Se alcanzo el maximo de 100 iteraciones sin converger."
  end subroutine newton

  ! ----------------------------------------------------------------------------
  ! Subrutina: secante
  ! Entradas: x0_in, x1_in (dos puntos iniciales distintos) [real(dp)]
  !           n (parámetro de la función de Wien) [integer]
  !           tol (tolerancia de convergencia epsilon) [real(dp)]
  ! Salidas:  raiz (aproximación al cero de F(x)) [real(dp)]
  !           iter (número de iteraciones realizadas) [integer]
  !           historial (vector con las estimaciones x_{k+1} de cada paso) [real(dp)]
  ! Descripción: Aproxima la derivada mediante diferencias finitas entre dos
  !              iteraciones previas. Se detiene con mensaje si f(x1) - f(x0) = 0
  !              o cuando |x_{k+1} - x_k| < tol.
  ! ----------------------------------------------------------------------------
  subroutine secante(x0_in, x1_in, n, tol, raiz, iter, historial)
    real(dp), intent(in)  :: x0_in, x1_in, tol
    integer, intent(in)   :: n
    real(dp), intent(out) :: raiz
    integer, intent(out)  :: iter
    real(dp), intent(out) :: historial(max_iter)

    ! Variables locales de puntos, evaluaciones y denominador
    real(dp) :: x0, x1, x2, f0, f1, denom

    ! Inicialización del estado
    x0        = x0_in
    x1        = x1_in
    f0        = f_wien(x0, n)
    f1        = f_wien(x1, n)
    iter      = 0
    historial = 0.0_dp

    ! Ciclo iterativo de aproximación secante
    do while (iter < max_iter)
      denom = f1 - f0

      ! Protección contra denominador nulo
      if (abs(denom) <= tiny(1.0_dp)) then
        print *, "Error (Secante): Denominador nulo (f(x1) - f(x0) = 0). Deteniendo."
        raiz = x1
        return
      end if

      iter            = iter + 1
      x2              = x1 - f1 * (x1 - x0) / denom
      historial(iter) = x2

      ! Criterio de parada por diferencia entre iteraciones consecutivas
      if (abs(x2 - x1) < tol) then
        raiz = x2
        return
      end if

      ! Desplazamiento de puntos para la siguiente iteración
      x0 = x1
      f0 = f1
      x1 = x2
      f1 = f_wien(x1, n)
    end do

    raiz = x1
    print *, "Advertencia (Secante): Se alcanzo el maximo de 100 iteraciones sin converger."
  end subroutine secante

end module modulo_ceros
