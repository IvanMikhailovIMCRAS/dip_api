!######################################################################################
program DIP
!######################################################################################
	use CommonParam
	use InputOutput
	implicit none
	integer(4) :: N1, N2, N_s, n_layer, nfree
	real(8) :: sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
	common /inter_int/ N1, N2, N_s, n_layer, nfree
	common /inter_real / sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
	open(n_info, file=name_info)
	
	call get_input() ! считываем входные данные из файла INPUT (see InputOutput.f90)
	
	call main()      ! выполняем расчёты (see OnePoint.f90 and Lib.f90)
	
	write(n_info,*) 'Awesome! That worked out!'
	close(n_info)
	stop
end program DIP
!######################################################################################

!**************************************************************************************
subroutine main()
!**************************************************************************************
	use CommonParam
	use OnePoint
	implicit none
	integer(4) :: N1, N2, N_s, n_layer, nfree
	real(8) :: sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
	real(8) :: F ! free energy
	common /inter_int/ N1, N2, N_s, n_layer, nfree
	common /inter_real / sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s, eta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
	!!!! расчёт для одной точки с конкретными параметрами
	!!!! при желании здесь удобно организовать цикл по какому-либо параметру
	
	call calc_one_point(N1, N2, sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, &
           & tau_s, N_s, eta, nfree, n_layer, F)
	
 
return
end subroutine main
!**************************************************************************************
