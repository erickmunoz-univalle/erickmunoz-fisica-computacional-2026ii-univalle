! ==============================================================================
! Programa Principal: programa_nucleo
! Asignatura: Física Computacional (106018C) - Taller 6
! Descripción: Calcula los términos de la fórmula semiempírica para 8 núclidos
!              de prueba mediante un bucle directo, genera la curva de máximo
!              enlace (valle.dat) e identifica el núcleo más ligado (maximo.dat).
! Autor: Erick Alexander Muñoz Orozco
! Fecha: 02/11/2026
! ==============================================================================
program programa_nucleo
  use modulo_nucleo, only: dp, e_volumen, e_superficie, e_coulomb, &
                           e_asimetria, E_B, M, Z_estable
  implicit none

  ! Declaración de variables locales
  integer :: i, a_loop, z_opt, u_nuc, u_valle, u_max
  real(dp):: eb_opt, eb_global_min
  integer :: a_global_min, z_global_min

  ! Variables temporales para la lectura dentro del bucle de los 8 núclidos
  integer :: z_i, a_i

  ! Arreglos con los 8 núclidos de referencia en el orden del enunciado:
  ! 1) Pu-244, 2) Ge-70, 3) Hf-179, 4) Ar-38, 5) Fe-56, 6) Ca-40, 7) U-235, 8) U-238
  integer, dimension(8), parameter :: Z_ref = [94, 32, 72, 18, 26, 20, 92, 92]
  integer, dimension(8), parameter :: A_ref = [244, 70, 179, 38, 56, 40, 235, 238]

  ! ----------------------------------------------------------------------------
  ! 1. Generación de nucleos.dat con un solo bucle DO
  ! ----------------------------------------------------------------------------
  open(newunit=u_nuc, file='nucleos.dat', status='replace', action='write')

  do i = 1, 8
     z_i = Z_ref(i)
     a_i = A_ref(i)

     ! Formato: Z, A, Evol, Esup, Ecoul, Easim, EB, M
     write(u_nuc, '(i4, 1x, i4, 1x, 5(f16.8, 1x), f16.4)') &
          z_i, a_i, &
          e_volumen(a_i), &
          e_superficie(a_i), &
          e_coulomb(z_i, a_i), &
          e_asimetria(z_i, a_i), &
          E_B(z_i, a_i), &
          M(z_i, a_i)
  end do

  close(u_nuc)

  ! ----------------------------------------------------------------------------
  ! 2. Generación de valle.dat (A = 4 hasta 260) e identificación de maximo.dat
  ! ----------------------------------------------------------------------------
  open(newunit=u_valle, file='valle.dat', status='replace', action='write')

  ! Inicializar el seguimiento del mínimo global con A = 4
  call Z_estable(4, z_opt, eb_opt)
  a_global_min  = 4
  z_global_min  = z_opt
  eb_global_min = eb_opt

  ! Escribir primera fila (A = 4)
  write(u_valle, '(i5, 1x, i5, 1x, f16.8)') 4, z_opt, eb_opt

  ! Recorrer desde A = 5 hasta 260
  do a_loop = 5, 260
     call Z_estable(a_loop, z_opt, eb_opt)
     write(u_valle, '(i5, 1x, i5, 1x, f16.8)') a_loop, z_opt, eb_opt

     ! Actualizar si se encuentra un núcleo más ligado (E_B más negativo)
     if (eb_opt < eb_global_min) then
        eb_global_min = eb_opt
        a_global_min  = a_loop
        z_global_min  = z_opt
     end if
  end do

  close(u_valle)

  ! ----------------------------------------------------------------------------
  ! 3. Escritura del núcleo más ligado de toda la curva en maximo.dat
  ! ----------------------------------------------------------------------------
  open(newunit=u_max, file='maximo.dat', status='replace', action='write')
  write(u_max, '(i5, 1x, i5, 1x, f16.8)') a_global_min, z_global_min, eb_global_min
  close(u_max)

end program programa_nucleo
