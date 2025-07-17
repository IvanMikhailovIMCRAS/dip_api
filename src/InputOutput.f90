!######################################################################################
module InputOutput
!######################################################################################
use CommonParam
use ErrorList 

Contains

!**************************************************************************************
subroutine get_input()	
!**************************************************************************************
! считывание входных параметров из внешнего файла INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer(4) :: N1, N2, N_s, n_layer, nfree
	real(8) :: sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
	common /inter_int/ N1, N2, N_s, n_layer, nfree
	common /inter_real / sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
	integer(4) ioer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	open(n_input,file=name_input,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(0)
	
	read(n_input,*,iostat=ioer) n_layer ! Число слоёв (расстояние между стенками)
	if (ioer.ne.0.or.n_layer.lt.1) call ERROR(13)
	
	read(n_input,*,iostat=ioer) N1 ! Степень полимеризации привитых полимеров
	if (ioer.ne.0.or.N1.lt.1) call ERROR(1)
	
	read(n_input,*,iostat=ioer) N2 ! Степень полимеризации привитых полимеров
	if (ioer.ne.0.or.N2.lt.1) call ERROR(1)
	
	read(n_input,*,iostat=ioer) sigma1 ! Плотность прививки к одной из стенок
	if (ioer.ne.0.or.sigma1.gt.1.0.or.sigma1.lt.0.0) call ERROR(2)
	
	read(n_input,*,iostat=ioer) sigma2 ! Плотность прививки к одной из стенок
	if (ioer.ne.0.or.sigma2.gt.1.0.or.sigma2.lt.0.0) call ERROR(2)
	
	read(n_input,*,iostat=ioer) p1   ! Длина сегмента Куна
	if (ioer.ne.0.or.p1.lt.1.0) call ERROR(3)
	
	read(n_input,*,iostat=ioer) p2   ! Длина сегмента Куна
	if (ioer.ne.0.or.p2.lt.1.0) call ERROR(3)
	
	read(n_input,*,iostat=ioer) ps   ! Длина сегмента Куна
	if (ioer.ne.0.or.ps.lt.1.0) call ERROR(3)
	
	read(n_input,*,iostat=ioer) chi1 ! Параметр Флори-Хаггинса (полимер-растворитель)
	if (ioer.ne.0) call ERROR(4)
	
	read(n_input,*,iostat=ioer) chi2 ! Параметр Флори-Хаггинса (полимер-растворитель)
	if (ioer.ne.0) call ERROR(4)
	
	read(n_input,*,iostat=ioer) tau1 ! Коэффициент пропорциональный дипольному моменту
	if (ioer.ne.0) call ERROR(5)
	
	read(n_input,*,iostat=ioer) tau2 ! Коэффициент пропорциональный дипольному моменту
	if (ioer.ne.0) call ERROR(5)
	
	read(n_input,*,iostat=ioer) tau_s ! Коэффициент пропорциональный дипольному моменту
	if (ioer.ne.0) call ERROR(5)
	
	read(n_input,*,iostat=ioer) N_s ! Степень полимеризации молекул растворителя
	if (ioer.ne.0.or.N_s.lt.0) call ERROR(6)
	
	read(n_input,*,iostat=ioer) eta ! Начальная скорость градиентного спуска
	if (ioer.ne.0.or.eta.gt.0.5) call ERROR(9)
	
	read(n_input,*,iostat=ioer) nfree ! Максимальное количество "холостых шагов" спуска
	if (ioer.ne.0.or.nfree.lt.0) call ERROR(10)
			 
	close(n_input)
	
	if ((sigma1*dble(N1)+sigma2*dble(N2)).gt.dble(n_layer)) call ERROR(14)
	
	return
end subroutine get_input
!**************************************************************************************

!**************************************************************************************
subroutine print_profile(n_layer, alpha, phi1, phi2, dist1, dist2)
!**************************************************************************************
	implicit none
	integer(4), intent(in) :: n_layer
	real(8), intent(in) :: phi1(0:n_layer+1), phi2(0:n_layer+1)
	real(8), intent(in) :: dist1(0:n_layer+1), dist2(0:n_layer+1)  
	real(8), intent(in) :: alpha(0:n_layer+1)
	integer(4) z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_pro,file=name_pro)
	! print title
	write(n_pro,'(7x,a,9x,a,4(18x,a))') &
	&   				'z', 'alpha', 'phi1', 'phi2', 'end1', 'end2'
	do z = 1, n_layer
		write(n_pro,'(I8,5(1x,E20.10))') &
	&			 z, alpha(z), phi1(z), phi2(z), dist1(z), dist2(z)
	enddo
	
	close(n_pro)
	
end subroutine print_profile
!**************************************************************************************

!******************************************************************************
subroutine print_spins(n_layer,phi1,phi2)
!******************************************************************************
	implicit none
	integer(4), intent(in) :: n_layer
	real(8), intent(in) :: phi1(0:n_layer+1,-1:1), phi2(0:n_layer+1,-1:1) 
	integer(4) z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(n_spn,file=name_spn)
	! print title
	write(n_spn,'(7x,a,6(6x,a,7x))') 'z', 'phi1(-1) ', ' phi1(0) ', 'phi1(+1) ', &
	&									'phi2(-1)', 'phi2(0) ', 'phi2(+1)' 
	do z = 1, n_layer
		write(n_spn,'(I8,6(1x,E20.10))') z, phi1(z,-1), phi1(z,0), phi1(z,+1), &
		&                                   phi2(z,-1), phi2(z,0), phi2(z,+1)  
	enddo
	
	close(n_spn)
	
end subroutine print_spins
!******************************************************************************


!######################################################################################
end module InputOutput
!######################################################################################
