program taller1_programa_hmax
  !definiendo variables y parámetros
  implicit none
  real:: v0, angulo_grados, angulo_rad, h_max
  real, parameter:: g = 9.8, pi = 4.0 * atan(1.0)

  !leyendo entradas de v0 y angulo_grados
  print *, "Ingrese la rapidez inicial del proyectil: "
  read *, v0
  print *, "Ingrese el ángulo respecto a la horizontal: "
  read *, angulo_grados

  !convirtiendo angulo_grados a radianes
  angulo_rad = angulo_grados * PI / 180

  !calculando h_max
  h_max = (v0**2 * (SIN(angulo_rad))**2) / (2*g)

  !imprimiendo resultados
  print *, ""
  print *, "===================== RESULTADOS ===================="
  print *, "rapidez inicial :", v0, "m/s "
  print *, "ángulo :", angulo_grados, "grados"
  print *, "Altura máxima del proyectil :", h_max, "metros"
  print *, "====================================================="
  print *, ""
end program taller1_programa_hmax
