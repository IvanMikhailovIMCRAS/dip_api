!######################################################################################
module OnePoint
!######################################################################################
Contains

!***************************************************************************************
subroutine calc_one_point(N1, N2, sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, & 
         & tau_s, N_s, eta, nfree, n_layer, F)	
!***************************************************************************************
	use CommonParam
	use ErrorList
	use Lib
	use InputOutput
	implicit none
	integer(4), intent(in) :: N1, N2, N_s, nfree, n_layer
	real(8), intent(in) :: sigma1, sigma2, p1, p2, ps, chi1, chi2, tau1, tau2, tau_s
	real(8), intent(inout) :: eta
	real(8), intent(out) :: F
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8)  phi_solv1(0:n_layer+1), phi_solv2(0:n_layer+1), alpha(0:n_layer+1)
	real(8)  phi_pol1(0:n_layer+1), phi_pol2(0:n_layer+1)
	real(8)  phi1(0:n_layer+1,-1:1), phi2(0:n_layer+1,-1:1)
	real(8)  phis1(0:n_layer+1,-1:1), phis2(0:n_layer+1,-1:1)
	real(8)  WB1(0:n_layer+1,-1:1), WB2(0:n_layer+1,-1:1)
	real(8)  WBs1(0:n_layer+1,-1:1), WBs2(0:n_layer+1,-1:1) 
	
	real(8)  U1(0:n_layer+1,-1:1), U2(0:n_layer+1,-1:1)
	real(8)  UF1(0:n_layer+1,-1:1), UF2(0:n_layer+1,-1:1)
	real(8)  Us1(0:n_layer+1,-1:1), Us2(0:n_layer+1,-1:1)
	real(8)  W1(0:n_layer+1,-1:1), W2(0:n_layer+1,-1:1)
	real(8)  Ws1(0:n_layer+1,-1:1), Ws2(0:n_layer+1,-1:1)
	
	real(8)  s1_1(0:n_layer+1), s1_2(0:n_layer+1)
	real(8)  s1_s1(0:n_layer+1), s1_s2(0:n_layer+1) 
	real(8)  Gf1(0:n_layer+1,0:N1,-1:1), Gb1(0:n_layer+1,1:N1,-1:1)
	real(8)  Gf2(0:n_layer+1,0:N2,-1:1), Gb2(0:n_layer+1,1:N2,-1:1) 
	real(8)  Sf1(0:n_layer+1,1:N_s,-1:1), Sb1(0:n_layer+1,1:N_s,-1:1)
	real(8)  Sf2(0:n_layer+1,1:N_s,-1:1), Sb2(0:n_layer+1,1:N_s,-1:1)
	real(8)  dist1(0:n_layer+1), dist2(0:n_layer+1)
	real(8)  dists1(0:n_layer+1), dists2(0:n_layer+1) 
	real(8)  mp1(0:n_layer+1), mp2(0:n_layer+1), mps1(0:n_layer+1), mps2(0:n_layer+1) 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8) lambda1(-1:1), lambda2(-1:1), lambda_s(-1:1)
	real(8) deviation, deviation_old
	real(8) part_fun_pol1, part_fun_pol2, part_fun_solv1, part_fun_solv2, theta1, theta2
	integer(4) z, s
	integer(4) iter, counter_free
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

	theta1 = dble(N1) * sigma1 ! нормировка профиля объёмной доли левой щётки
	theta2 = dble(N2) * sigma2 ! нормировка профиля объёмной доли правой щётки
    
    !!!! zero initialization of arrays                                       
	phi_pol1(:) = 0.0; phi_pol2(:) = 0.0; phi_solv1(:) = 0.0; phi_solv2(:) = 0.0;
	alpha(:) = 0.0
	phi1(:,:) = 0.0; phi2(:,:) = 0.0; phis1(:,:) = 0.0; phis2(:,:) = 0.0; 
	WB1(:,:) = 0.0; WB2(:,:) = 0.0; WBs1(:,:) = 0.0; WBs2(:,:) = 0.0;
	U1(:,:) = 0.0; U2(:,:) = 0.0; W1(:,:) = 0.0; W2(:,:) = 0.0
	Us1(:,:) = 0.0; Us2(:,:) = 0.0; Ws1(:,:) = 0.0; Ws2(:,:) = 0.0
	s1_1(:) = 0.0; s1_2(:) = 0.0; s1_s1(:) = 0.0; s1_s2(:) = 0.0;
	Gf1(:,:,:) = 0.0; Gb1(:,:,:) = 0.0; Gf2(:,:,:) = 0.0; Gb2(:,:,:) = 0.0 
	Sf1(:,:,:) = 0.0; Sb1(:,:,:) = 0.0; Sf2(:,:,:) = 0.0; Sb2(:,:,:) = 0.0;
	dist1(:) = 0.0; dist2(:) = 0.0; dists1(:) = 0.0; dists2(:) = 0.0
	mp1(:) = 0.0; mp2(:) = 0.0; mps1(:) = 0.0; mps2(:) = 0.0;
	
! рассчитываем вероятности появления дискретных конфигураций угла (see CommonParam.f90)
	call calc_trans_probabilities(p1,lambda1)
	call calc_trans_probabilities(p2,lambda2)
	call calc_trans_probabilities(ps,lambda_s)
	
	iter = 0 ! cчётчик итераций
	counter_free = 0 ! счётчик "холостых" шагов градиентного спуска
	deviation = dble(n_layer)*max_number ! затравочное отклонение для входа в цикл
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ENGINE of PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do while (deviation.gt.tolerance)
		iter = iter + 1
		if (iter.gt.n_max_iter) call ERROR(11)
		!!!! old potentials <<< new potentials
		W1(:,:) = U1(:,:); W2(:,:) = U2(:,:); Ws1(:,:) = Us1(:,:); Ws2(:,:) = Us2(:,:)  	
		!!!! Boltzmann weights for polymers:
		! polymer tethered to LEFT wall:
		WB1(1:n_layer,-1) = exp(-alpha(1:n_layer)-U1(1:n_layer,-1))
		WB1(1:n_layer, 0) = exp(-alpha(1:n_layer)-U1(1:n_layer, 0))
		WB1(1:n_layer,+1) = exp(-alpha(1:n_layer)-U1(1:n_layer,+1))
		! polymer tethered to RIGHT wall:
		WB2(1:n_layer,-1) = exp(-alpha(1:n_layer)-U2(1:n_layer,-1))
		WB2(1:n_layer, 0) = exp(-alpha(1:n_layer)-U2(1:n_layer, 0))
		WB2(1:n_layer,+1) = exp(-alpha(1:n_layer)-U2(1:n_layer,+1))
		!!!! Boltzmann weights for solvent:
		WBs1(1:n_layer,-1) = exp(-alpha(1:n_layer)-Us1(1:n_layer,-1))
		WBs1(1:n_layer, 0) = exp(-alpha(1:n_layer)-Us1(1:n_layer, 0))
		WBs1(1:n_layer,+1) = exp(-alpha(1:n_layer)-Us1(1:n_layer,+1))
		!!!!!
		WBs2(1:n_layer,-1) = exp(-alpha(1:n_layer)-Us2(1:n_layer,-1)) ! +++++++++++++
		WBs2(1:n_layer, 0) = exp(-alpha(1:n_layer)-Us2(1:n_layer, 0))
		WBs2(1:n_layer,+1) = exp(-alpha(1:n_layer)-Us2(1:n_layer,+1))
				
		!!!! Initial guess for polymer propagators:
		!! polymer tethered to LEFT wall:
		! forward propagators
		Gf1(0,0,+1) = WB1(1,+1) 
		Gf1(0,0, 0) = 0.0
		Gf1(0,0,-1) = 0.0 
		! back propagators
		Gb1(1:n_layer,N1,+1) = WB1(1:n_layer,-1) 
		Gb1(1:n_layer,N1, 0) = WB1(1:n_layer, 0) 
		Gb1(1:n_layer,N1,-1) = WB1(1:n_layer,+1) 
		!! polymer tethered to RIGHT wall:
		! forward propagators
		Gf2(n_layer+1,0,+1) = WB2(n_layer,+1) 
		Gf2(n_layer+1,0, 0) = 0.0
		Gf2(n_layer+1,0,-1) = 0.0 
		! back propagators
		Gb2(1:n_layer,N2,+1) = WB2(1:n_layer,-1) 
		Gb2(1:n_layer,N2, 0) = WB2(1:n_layer, 0) 
		Gb2(1:n_layer,N2,-1) = WB2(1:n_layer,+1) 
		
		!!!! Initial guess for solvent propagators:
		! forward propagators
		Sf1(1:n_layer,1,+1) = WBs1(1:n_layer,+1) 
		Sf1(1:n_layer,1, 0) = WBs1(1:n_layer, 0) 
		Sf1(1:n_layer,1,-1) = WBs1(1:n_layer,-1) 
		! back propagators
		Sb1(1:n_layer,N_s,+1) = WBs1(1:n_layer,-1) 
		Sb1(1:n_layer,N_s, 0) = WBs1(1:n_layer, 0) 
		Sb1(1:n_layer,N_s,-1) = WBs1(1:n_layer,+1) 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! forward propagators
		Sf2(1:n_layer,1,+1) = WBs2(1:n_layer,+1) 
		Sf2(1:n_layer,1, 0) = WBs2(1:n_layer, 0) 
		Sf2(1:n_layer,1,-1) = WBs2(1:n_layer,-1) 
		! back propagators
		Sb2(1:n_layer,N_s,+1) = WBs2(1:n_layer,-1) 
		Sb2(1:n_layer,N_s, 0) = WBs2(1:n_layer, 0) 
		Sb2(1:n_layer,N_s,-1) = WBs2(1:n_layer,+1)  
		
		!!!! "Fill up" the propagators 
		call forward_propagate(N1, n_layer, lambda1, WB1, Gf1, 0)
		call back_propagate(N1, n_layer, lambda1, WB1, Gb1)
		!!!!
		call forward_propagateR(N2, n_layer, lambda2, WB2, Gf2, 0)
		call back_propagateR(N2, n_layer, lambda2, WB2, Gb2)
		!!!!
		call forward_propagate(N_s, n_layer, lambda_s, WBs1, Sf1, 1)
		call back_propagate(N_s, n_layer, lambda_s, WBs1, Sb1)
		!!!!
		call forward_propagate(N_s, n_layer, lambda_s, WBs2, Sf2, 1)
		call back_propagate(N_s, n_layer, lambda_s, WBs2, Sb2)
		
		!!!! Calculate density profiles $phi$ and partition functions $part_fun$
		call phi_calc(n_layer, N1, Gf1, Gb1, WB1, phi1, part_fun_pol1, 0)
		call phi_calc(n_layer, N2, Gf2, Gb2, WB2, phi2, part_fun_pol2, 0)
		call phi_calc(n_layer, N_s, Sf1, Sb1, WBs1, phis1, part_fun_solv1, 1)
		call phi_calc(n_layer, N_s, Sf2, Sb2, WBs2, phis2, part_fun_solv2, 1)
		!!!! 
		phi1(:,-1) = phi1(:,-1) / part_fun_pol1 * theta1
		phi1(:, 0) = phi1(:, 0) / part_fun_pol1 * theta1
		phi1(:,+1) = phi1(:,+1) / part_fun_pol1 * theta1
		phi_pol1(:) = phi1(:,-1) + phi1(:,0) + phi1(:,+1)
		!!!! 
		phi2(:,-1) = phi2(:,-1) / part_fun_pol2 * theta2
		phi2(:, 0) = phi2(:, 0) / part_fun_pol2 * theta2
		phi2(:,+1) = phi2(:,+1) / part_fun_pol2 * theta2
		
		phi_pol2(:) = phi2(:,-1) + phi2(:,0) + phi2(:,+1)
		!!!!
		if (int(sigma1*dble(N1)+sigma2*dble(N2)).eq.n_layer.or.N_s.eq.0) then
			phis1(:,-1) = 0.0
			phis1(:, 0) = 0.0
			phis1(:,+1) = 0.0
			phis2(:,-1) = 0.0 
			phis2(:, 0) = 0.0
			phis2(:,+1) = 0.0
		else
			phis1(:,-1) = phis1(:,-1) / dble(N_s) / 6.0 / 2.0
			phis1(:, 0) = phis1(:, 0) / dble(N_s) / 6.0 / 2.0
			phis1(:,+1) = phis1(:,+1) / dble(N_s) / 6.0 / 2.0
			phis2(:,-1) = phis2(:,-1) / dble(N_s) / 6.0 / 2.0
			phis2(:, 0) = phis2(:, 0) / dble(N_s) / 6.0 / 2.0
			phis2(:,+1) = phis2(:,+1) / dble(N_s) / 6.0 / 2.0
		endif
		phi_solv1(:) = phis1(:,-1) + phis1(:,0) + phis1(:,+1)
		phi_solv2(:) = phis2(:,-1) + phis2(:,0) + phis2(:,+1)
		
		
		!!!! Вычисляем параметр порядка ориентации сегментов s_1
		call calc_orientation(n_layer, phi1, s1_1)
		call calc_orientation(n_layer, phi2, s1_2)
		call calc_orientation(n_layer, phis1, s1_s1)
		call calc_orientation(n_layer, phis2, s1_s2)

			
		!!!! Вычисляем потенциалы взаимодействий
		U1(:,:) = 0.0; U2(:,:) = 0.0; Us1(:,:) = 0.0; Us2(:,:) = 0.0
		!! взаимодействия Флори-Хаггинса (полимер-растворитель)
		call calc_U_Flory(n_layer, chi1, phi_solv1, 0.0_8, U1)
		call calc_U_Flory(n_layer, chi2, phi_solv1, 0.0_8, U2)
		call calc_U_Flory(n_layer, chi1, phi_solv2, 0.0_8, U1)
		call calc_U_Flory(n_layer, chi2, phi_solv2, 0.0_8, U2)
		
		!write(*,*) sum(U1(:,:)), sum(U2(:,:))
		
		call calc_U_Flory(n_layer, chi1, phi_pol1, 1.0_8, Us1)
		call calc_U_Flory(n_layer, chi1, phi_pol1, 1.0_8, Us2)
		call calc_U_Flory(n_layer, chi2, phi_pol2, 1.0_8, Us1)
		call calc_U_Flory(n_layer, chi2, phi_pol2, 1.0_8, Us2)
		!! диполь-дипольные взаимодействия
		call calc_dip_dip_interaction(n_layer, s1_1, s1_2, s1_s1, s1_s2, tau1, tau2, tau_s, U1, U2, Us1, Us2)
		
		!!!! градиентный спуск по потенциалам с помощью
        !!!! exponentially weighted moving average (EWMA) 											
		U1(:,0) = eta*U1(:,0) +  (1.0 - eta)*W1(:,0)
		U1(:,1) = eta*U1(:,1) +  (1.0 - eta)*W1(:,1)
		U1(:,-1) = eta*U1(:,-1) +  (1.0 - eta)*W1(:,-1)
		!!!!
		U2(:,0) = eta*U2(:,0) +  (1.0 - eta)*W2(:,0)
		U2(:,1) = eta*U2(:,1) +  (1.0 - eta)*W2(:,1)
		U2(:,-1) = eta*U2(:,-1) +  (1.0 - eta)*W2(:,-1)
		!!!!	
		Us1(:,0) = eta*Us1(:,0) +  (1.0 - eta)*Ws1(:,0)
		Us1(:,1) = eta*Us1(:,1) +  (1.0 - eta)*Ws1(:,1)
		Us1(:,-1) = eta*Us1(:,-1) +  (1.0 - eta)*Ws1(:,-1)
		!!!!	
		Us2(:,0) = eta*Us2(:,0) +  (1.0 - eta)*Ws2(:,0)
		Us2(:,1) = eta*Us2(:,1) +  (1.0 - eta)*Ws2(:,1)
		Us2(:,-1) = eta*Us2(:,-1) +  (1.0 - eta)*Ws2(:,-1)
				
		!!!! градиентный спуск по полю Лагранжа (стохастический градиентный спуск)
		! https://ruder.io/optimizing-gradient-descent/ - какие бывают модификации
		alpha(1:n_layer) = alpha(1:n_layer) +  &
		&	eta*(phi_pol1(1:n_layer) + phi_pol2(1:n_layer) + &
		&   phi_solv1(1:n_layer) + phi_solv2(1:n_layer) -1.0) 
		
		!!!! рассчитываем целевую функцию (отклонение dev = sum of the loss function)			
		deviation_old = deviation
		deviation = 0.0
		do z = 1, n_layer
		deviation = deviation + (phi_pol1(z) + phi_pol2(z) + phi_solv1(z) + phi_solv2(z) - 1.0)**2  + &
		&      (U1(z,-1)-W1(z,-1))**2 + (U1(z,0)-W1(z,0))**2 + (U1(z,1)-W1(z,1))**2 + &
		&      (U2(z,-1)-W2(z,-1))**2 + (U2(z,0)-W2(z,0))**2 + (U2(z,1)-W2(z,1))**2 + &
		&      (Us1(z,-1)-Ws1(z,-1))**2 + (Us1(z,0)-Ws1(z,0))**2 + (Us1(z,1)-Ws1(z,1))**2 + &
		&      (Us2(z,-1)-Ws2(z,-1))**2 + (Us2(z,0)-Ws2(z,0))**2 + (Us2(z,1)-Ws2(z,1))**2 
		
		enddo
		deviation = sqrt(deviation)
				
		!write(*,'(a,I9,1x,a,E16.10,1x,a,E16.10)') &
		!&		'iter = ', iter, 'dev = ', deviation, 'eta = ', eta
				
		!!!! автоматизация подбора шага (скорости) сходимости, eta
		if (deviation.ge.deviation_old) then
        	counter_free = counter_free + 1
        endif
        
        if (counter_free.eq.nfree.or.deviation.ne.deviation) then
            call check_deviation(deviation, deviation_old, eta, iter)
            write(n_info,'(a,I9,1x,a,E16.10,1x,a,E16.10)') &
		&		'iter = ', iter, 'dev = ', deviation, 'eta = ', eta
            alpha(:) = 0.0
            U1(:,:) = 0.0
            U2(:,:) = 0.0
            Us1(:,:) = 0.0
            Us2(:,:) = 0.0
            counter_free = 0
        endif
        
                       
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!! КОНЕЦ ОСНОВНОГО ЦИКЛА !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	write(n_info,'(a,I9,1x,a,E16.10,1x,a,E16.10)') &
		&		'iter = ', iter, 'dev = ', deviation, 'eta = ', eta
	
		
	!!!! перерасчитываем стат.суммы по методу Леермакерса
	part_fun_pol1 = (Gb1(1,1,-1) + 4.0*Gb1(1,1,0) + Gb1(1,1,+1)) / 6.0
	part_fun_pol2 = (Gb2(n_layer,1,-1) + 4.0*Gb2(n_layer,1,0) + Gb2(n_layer,1,+1)) / 6.0
	
	!!!!!!!!!!!! РАССЧИТЫВАЕМ СВОБОДНУЮ ЭНЕРГИЮ (F) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call calc_free_energy(n_layer, N1, N2, N_s, sigma1, sigma2, phi1, phi2, phis1, phis2, s1_1, s1_2,  &
&      s1_s1, s1_s2, U1, U2, Us1, Us2, alpha, tau1, tau2, tau_s, chi1, chi2, part_fun_pol1, part_fun_pol2, F)

	call end_point(n_layer, N1, Gf1, Gb1, WB1, dist1, 0)
	call end_point(n_layer, N2, Gf2, Gb2, WB2, dist2, 0)
	
	call print_profile(n_layer, alpha, phi_pol1, phi_pol2, dist1, dist2)
	!call print_profile(n_layer, alpha, phi_pol1, phi_pol2, &
	!& alpha(:)+U1(:,-1)/6.0+U1(:,0)*4.0/6.0+U1(:,+1)/6.0, alpha(:)+U2(:,-1)/6.0+U2(:,0)*4.0/6.0+U2(:,+1)/6.0)

!	write(*,'(a,1x,I6,5x,a,1x,E16.10)') 'n_layer:', n_layer, 'Free_Energy:', F
	
		
return
end subroutine calc_one_point
!**************************************************************************************

!######################################################################################
end module OnePoint
!######################################################################################
		
