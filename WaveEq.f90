!************************************
!Wave equation 1+1
!Pablo Villanueva Domingo
!Marzo 2016
!************************************

program WaveEq

implicit none

  real dt, dx, phi_in, psi_in, alpha, beta !Incrementos, funciones iniciales, funciones de la métrica
  integer i
  integer, parameter :: x_num=40 !Número de componentes de los vectores usados
  real t, t_max
  real, dimension(x_num) :: pi1, pi2, psi1, psi2, phi1, phi2 !Vectores de las variables pi, psi y phi, para instante inicial y final

!Incrementos de paso para t y x
dt=0.01
dx=0.2

!Tiempo incial y final
t=0
t_max=6

open (unit=6, file='wave.dat')

!Condiciones iniciales (ecs. 8-10 del pdf)
do i=1, x_num

  phi1(i) =  phi_in(i*dx,x_num*dx/2)
  psi1(i) =  psi_in(i*dx,x_num*dx/2)
  pi1(i)  =  0

end do

!Evolución
do while (t<= t_max)

	!Ecuación de evolución para pi (ec. 5 del pdf)
	call wave( x_num, t, dx, dt/(2*dx), pi1, pi2, psi1 )  
	
	!Ecuación de evolución para psi (ec. 6 del pdf)
	call wave( x_num, t, dx, dt/(2*dx), psi1, psi2, pi1 ) 
	
	!Evolución para phi, a partir de pi y psi (ec. 7 del pdf)

	do i=2, x_num-1

		phi2(i) =  phi1(i) + dt * ( alpha(i*dx)*pi1(i) + beta(i*dx)*psi1(i) )

	end do

    phi2(1)=phi2(2)
    phi2(x_num)=phi2(x_num-1)

	!Escribe en archivo los resultados
	do i=1, x_num

		write (6,*) t, i*dx, phi2(i)

	end do

	!Realiza el paso de tiempo
	pi1=pi2
	psi1=psi2
	phi1=phi2

	t=t+dt

end do

close(6)

end program WaveEq

!---------------------------------------------------
!Subrutina que calcula la evolución de las variables
!---------------------------------------------------
subroutine wave ( x_num, t, dx, c, u1, u2, v1 )

  implicit none

  integer x_num, i
  real c, t, dx, alpha, beta
  real, dimension(x_num) :: u1, u2, v1

do i=2, x_num-1

  !Evolución de variable u en el esquema FCTS
  u2(i) =  u1(i) + c * ( alpha((i+1)*dx)*v1(i+1) - alpha((i-1)*dx)*v1(i-1) &
  & + beta((i+1)*dx)*u1(i+1) - beta((i-1)*dx)*u1(i-1))

end do

u2(1)=u2(2)
u2(x_num)=u2(x_num-1)

end subroutine wave

!--------------------------------------
!Funciones de las condiciones iniciales
!--------------------------------------
function phi_in(x,x0)

	real x, phi_in, x0

	phi_in=exp(-(x-x0)**2)

return
end function phi_in

function psi_in(x,x0)

	real x, psi_in, x0

	psi_in=-2*(x-x0)*exp(-(x-x0)**2)

return
end function psi_in

!------------------------------------------------------------
!Funciones alpha y beta de la métrica
!(Por simplicidad, se han escogido independientes del tiempo)
!------------------------------------------------------------
function alpha(x)

	real x, alpha

	alpha=1

return
end function alpha

function beta(x)

	real x, beta

	beta=0

return
end function beta