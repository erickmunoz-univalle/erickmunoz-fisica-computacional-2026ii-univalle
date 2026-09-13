program taller1_pruebas_programa_hmax
!definiendo variables y parámetros
  implicit none
  real:: v0, angulo_grados, angulo_rad, h_max
  real, parameter:: g = 9.8, pi = 4.0 * atan(1.0)

  !variables para las pruebas
  real:: a1 = 90, a2 = 45 , a3 = 30 !ángulos
  real:: v1 = 9.8, v2 = 20, v3 = 20 !rapidez inicial

  !caso 1
  angulo_grados = a1
  v0 = v1
  !convirtiendo angulo_grados a radianes
  angulo_rad = angulo_grados * PI / 180
  !calculando h_max
  h_max = (v0**2 * (SIN(angulo_rad))**2) / (2*g)
  if (abs(h_max-4.9000)<0.0001) then !usamos cierto rango de tolerancia dado que los valores no serían exactamente iguales bit a bit
    print *, "caso 1: PASS"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 4.9000
  else
    print *, "caso 1: FAIL"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 4.9000
  end if
  print *, "Error relativo: ", abs(h_max-4.9000)/4.9000*100, "%" !calculamos el error relativo

  print *, ""

  !caso 2
  angulo_grados = a2
  v0 = v2
  !convirtiendo angulo_grados a radianes
  angulo_rad = angulo_grados * PI / 180
  !calculando h_max
  h_max = (v0**2 * (SIN(angulo_rad))**2) / (2*g)
  if (abs(h_max-10.2041)<0.0001) then !rango de tolerancia
    print *, "caso 2: PASS"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 10.2041
  else
    print *, "caso 2: FAIL"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 10.2041
  end if
  print *, "Error relativo: ", abs(h_max-10.2041)/10.2041*100, "%" !error relativo

  print *, ""

  !caso 3
  angulo_grados = a3
  v0 = v3
  !convirtiendo angulo_grados a radianes
  angulo_rad = angulo_grados * PI / 180
  !calculando h_max
  h_max = (v0**2 * (SIN(angulo_rad))**2) / (2*g)
  if (abs(h_max-5.1020)<0.0001) then !rango de tolerancia
    print *, "caso 3: PASS"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 5.1020
  else
    print *, "caso 3: FAIL"
    print *, "valor obtenido :", h_max
    print *, "valor esperado :", 5.1020
  end if
  print *, "Error relativo: ", abs(h_max-5.1020)/5.1020*100, "%" !error relativo

end program taller1_pruebas_programa_hmax
