program ajuste_pendulo
  use iso_fortran_env, only: real64, iostat_end ! Precisión y código de fin de archivo
  implicit none                                  ! Evita variables implícitas

  ! Constantes
  real(real64), parameter :: PI = acos(-1.0_real64)

  ! Variables de archivo y control
  integer :: unidad_entrada, unidad_salida, estado
  integer :: id, n_osc, n_datos

  ! Variables de lectura por fila
  real(real64) :: longitud_cm, angulo_deg, tiempo_s
  real(real64) :: x, y, T_periodo

  ! Acumuladores y parámetros del ajuste lineal
  real(real64) :: sx, sy, sxx, sxy
  real(real64) :: den, a, b, g
  real(real64) :: y_ajustada, sse, sst, r2

  ! 1. Inicialización de acumuladores
  n_datos = 0
  sx  = 0.0_real64
  sy  = 0.0_real64
  sxx = 0.0_real64
  sxy = 0.0_real64

  ! 2. Apertura del archivo de datos limpios
  open(newunit=unidad_entrada, file='pendulo_limpio.dat', status='old', &
       action='read', iostat=estado)

  if (estado /= 0) error stop 'No fue posible abrir el archivo pendulo_limpio.dat'

  ! 3. Primera pasada: Lectura de datos y cálculo de acumuladores
  do
     read(unidad_entrada, *, iostat=estado) id, longitud_cm, angulo_deg, n_osc, tiempo_s

     if (estado == iostat_end) exit  ! Fin normal del archivo
     if (estado /= 0) error stop 'Hay una fila malformada en pendulo_limpio.dat'

     ! Conversión de unidades
     x = longitud_cm / 100.0_real64                   ! Longitud L en metros
     T_periodo = tiempo_s / real(n_osc, real64)       ! Período T en segundos
     y = T_periodo**2                                 ! y = T^2 en s^2

     ! Acumuladores
     n_datos = n_datos + 1
     sx  = sx  + x      ! Suma de x
     sy  = sy  + y      ! Suma de y
     sxx = sxx + x*x    ! Suma de x^2
     sxy = sxy + x*y    ! Suma de x*y
  end do

  ! Comprobación de mínimo 2 datos
  if (n_datos < 2) error stop 'Se requieren al menos 2 mediciones para realizar el ajuste'

  ! 4. Cálculo de pendiente (a), intercepto (b) y gravedad (g)
  den = real(n_datos, real64)*sxx - sx*sx

  if (abs(den) <= tiny(den)) error stop 'No se puede calcular la pendiente'

  a = (real(n_datos, real64)*sxy - sx*sy) / den ! Pendiente
  b = (sy - a*sx) / real(n_datos, real64)       ! Intercepto
  g = (4.0_real64 * PI**2) / a                  ! Gravedad

  ! 5. Segunda pasada: cálculo de residuos, SSE y variación total SST
  ! Inicialización de sse y  sst
  sse = 0.0_real64
  sst = 0.0_real64

  ! Volver al inicio del archivo
  rewind(unidad_entrada)

  ! Lectura de datos y cálculo de SSE y SST
  do
     read(unidad_entrada, *, iostat=estado) id, longitud_cm, angulo_deg, n_osc, tiempo_s

     if (estado == iostat_end) exit

     x = longitud_cm / 100.0_real64
     T_periodo = tiempo_s / real(n_osc, real64)
     y = T_periodo**2

     y_ajustada = a*x + b               ! Predicción de la recta
     sse = sse + (y - y_ajustada)**2    ! Suma de residuos al cuadrado
     sst = sst + (y - sy / real(n_datos, real64))**2  ! Suma de variación total alrededor de la media
  end do
  ! Cierre del archivo
  close(unidad_entrada)

  ! R^2 solo está definido aquí si los valores de y no son todos iguales
  if (sst > 0.0_real64) then
     r2 = 1.0_real64 - (sse / sst)
  else
     r2 = 0.0_real64
  end if

  ! 6. Impresión de resultados en pantalla
  print '(A)', '--------------------------------------------------'
  print '(A)', '   RESULTADOS DEL AJUSTE LINEAL (FORTRAN 2008)    '
  print '(A)', '--------------------------------------------------'
  print '(A, I5)',     'Numero de datos (N)    : ', n_datos
  print '(A, F12.6)', 'Pendiente a (s^2/m)    : ', a
  print '(A, F12.6)', 'Intercepto b (s^2)     : ', b
  print '(A, F12.6)', 'Coeficiente R^2        : ', r2
  print '(A, F12.6)', 'Gravedad g (m/s^2)     : ', g
  print '(A)', '--------------------------------------------------'

  ! 7. Escritura en resultados_ajuste.dat para verificacion PASS/FAIL
  open(newunit=unidad_salida, file='resultados_ajuste.dat', &
     status='replace', action='write')
  write(unidad_salida, *) n_datos, a, b, r2, g
  close(unidad_salida)

end program ajuste_pendulo
