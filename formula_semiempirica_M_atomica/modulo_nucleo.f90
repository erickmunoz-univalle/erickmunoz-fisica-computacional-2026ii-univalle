! ==============================================================================
! Módulo: modulo_nucleo
! Asignatura: Física Computacional (106018C) - Taller 6
! Descripción: Implementa las constantes y los procedimientos para la fórmula
!              semiempírica de masa atómica (Modelo de la gota líquida de
!              Bethe y Weizsäcker).
! Autor: Erick Alexander Muñoz Orozco
! Fecha: 02/11/2026
! ==============================================================================
module modulo_nucleo
  implicit none

  ! Definiendo doble precisión como dp
  integer, parameter :: dp = kind(1.0d0)

  ! ----------------------------------------------------------------------------
  ! Constantes del modelo
  ! ----------------------------------------------------------------------------
  real(dp), parameter :: a_V = 15.8_dp     ! Coeficiente de volumen [MeV]
  real(dp), parameter :: a_S = 18.3_dp     ! Coeficiente de superficie [MeV]
  real(dp), parameter :: a_C = 0.71_dp     ! Coeficiente de Coulomb [MeV]
  real(dp), parameter :: a_A = 92.7_dp     ! Coeficiente de asimetría [MeV]
  real(dp), parameter :: M_n = 939.57_dp   ! Masa del neutrón libre [MeV]
  real(dp), parameter :: M_H = 938.78_dp   ! Masa del átomo de hidrógeno [MeV]

contains

  ! ----------------------------------------------------------------------------
  ! Función: e_volumen
  ! Entrada: a_num (número de nucleones A)
  ! Salida:  término de volumen por nucleón [MeV] (-a_V)
  ! ----------------------------------------------------------------------------
  real(dp) function e_volumen(a_num)
    integer, intent(in) :: a_num
    real(dp) :: a_real
    a_real = real(a_num, dp) ! Evita advertencias de argumento sin usar
    e_volumen = -a_V
  end function e_volumen

  ! ----------------------------------------------------------------------------
  ! Función: e_superficie
  ! Entrada: a_num (número de nucleones A)
  ! Salida:  término de superficie por nucleón [MeV] (a_S * A^(-1/3))
  ! ----------------------------------------------------------------------------
  real(dp) function e_superficie(a_num)
    integer, intent(in) :: a_num
    real(dp) :: a_real
    a_real = real(a_num, dp)
    e_superficie = a_S * (a_real**(-1.0_dp / 3.0_dp))
  end function e_superficie

  ! ----------------------------------------------------------------------------
  ! Función: e_coulomb
  ! Entradas: z_num (protones Z), a_num (nucleones A)
  ! Salida:   término de Coulomb por nucleón [MeV] (a_C * Z^2 * A^(-4/3))
  ! ----------------------------------------------------------------------------
  real(dp) function e_coulomb(z_num, a_num)
    integer, intent(in) :: z_num, a_num
    real(dp) :: z_real, a_real
    z_real = real(z_num, dp)
    a_real = real(a_num, dp)
    e_coulomb = a_C * (z_real**2) * (a_real**(-4.0_dp / 3.0_dp))
  end function e_coulomb

  ! ----------------------------------------------------------------------------
  ! Función: e_asimetria
  ! Entradas: z_num (protones Z), a_num (nucleones A)
  ! Salida:   término de asimetría por nucleón [MeV] (a_A * (Z - A/2)^2 / A^2)
  ! ----------------------------------------------------------------------------
  real(dp) function e_asimetria(z_num, a_num)
    integer, intent(in) :: z_num, a_num
    real(dp) :: z_real, a_real
    z_real = real(z_num, dp)
    a_real = real(a_num, dp)
    e_asimetria = a_A * ((z_real - a_real / 2.0_dp)**2) / (a_real**2)
  end function e_asimetria

  ! ----------------------------------------------------------------------------
  ! Función: E_B
  ! Entradas: z_num (protones Z), a_num (nucleones A)
  ! Salida:   energía de enlace por nucleón total [MeV] (suma de los 4 términos)
  ! ----------------------------------------------------------------------------
  real(dp) function E_B(z_num, a_num)
    integer, intent(in) :: z_num, a_num
    E_B = e_volumen(a_num) + e_superficie(a_num) + &
          e_coulomb(z_num, a_num) + e_asimetria(z_num, a_num)
  end function E_B

  ! ----------------------------------------------------------------------------
  ! Función: M
  ! Entradas: z_num (protones Z), a_num (nucleones A)
  ! Salida:   masa atómica aproximada en equivalente de energía [MeV]
  !           M = Z*M_H + (A - Z)*M_n + A*E_B(Z,A)
  ! ----------------------------------------------------------------------------
  real(dp) function M(z_num, a_num)
    integer, intent(in) :: z_num, a_num
    real(dp) :: z_real, a_real
    z_real = real(z_num, dp)
    a_real = real(a_num, dp)
    M = z_real * M_H + (a_real - z_real) * M_n + a_real * E_B(z_num, a_num)
  end function M

  ! ----------------------------------------------------------------------------
  ! Subrutina: Z_estable
  ! Entrada:  a_num (número de nucleones A) [intent(in)]
  ! Salidas:  z_optimo (Z de máximo enlace), eb_min (E_B más negativo) [intent(out)]
  ! Descripción: Recorre Z de 1 a A para encontrar el núcleo con menor E_B
  !              (más ligado por nucleón). Devuelve dos resultados.
  ! ----------------------------------------------------------------------------
  subroutine Z_estable(a_num, z_optimo, eb_min)
    integer, intent(in)  :: a_num
    integer, intent(out) :: z_optimo
    real(dp), intent(out):: eb_min

    ! Variables locales
    integer  :: z_cand
    real(dp) :: eb_cand

    ! Inicialización con la primera alternativa Z = 1
    z_optimo = 1
    eb_min   = E_B(1, a_num)

    ! Búsqueda del mínimo E_B en el intervalo Z = 2, ..., A
    do z_cand = 2, a_num
       eb_cand = E_B(z_cand, a_num)
       if (eb_cand < eb_min) then
          eb_min   = eb_cand
          z_optimo = z_cand
       end if
    end do
  end subroutine Z_estable

end module modulo_nucleo
