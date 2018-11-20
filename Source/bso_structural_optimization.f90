module struct_param_optimization

   use bayes_pest_control
   use model_input_output
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use jupiter_input_data_support
   use nelder_mead_simplex_routines
   
  contains
  
 subroutine ABIC_struct_param_optim (d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,it_bga,n)
  
   implicit none
   
   type(kernel_XQR),     intent(in)     :: d_XQR
   type(cv_observ),      intent(in)     :: cv_OBS
   type(d_observ),       intent(in)     :: d_OBS
   type(d_algorithmic),  intent(inout)  :: d_A 
   type(cv_algorithmic), intent(inout)  :: cv_A
   type(d_prior_mean),   intent(in)     :: d_PM
   type(cv_param),       intent(in)     :: cv_PAR 
   type(d_param),        intent(inout)  :: d_PAR        
   type(cv_struct),      intent(inout)  :: cv_S 
   type(d_struct),       intent(inout)  :: d_S
   type(Q0_compr),       intent(in)     :: Q0_All(:)
   type(cv_prior_mean),  intent(in)     :: cv_PM
   character (len=ERRORWIDTH)           :: retmsg
   
   integer                           :: it_bga, i, j, k
   integer                           :: nQ0 = 0 !Dimension of Q0_All(:) [0] in case of no compression [cv_PAR%p] in case of compression
   double precision                  :: stc = 1. !Parameter that control values of step. step = stc*initial values
   integer ( kind = 4 )              n       !Indicates the number of pars to estimate.
   integer ( kind = 4 )              icount  !Number of function evaluation used (output)
   integer ( kind = 4 )              ifault  !Error indicator (output)
   integer ( kind = 4 )              konvge
   integer ( kind = 4 )              numres  !Number of restarts (output)
   real    ( kind = 8 )              reqmin
   real    ( kind = 8 )              start(n)
   real    ( kind = 8 )              step(n)
   real    ( kind = 8 )              xmin(n)
   real    ( kind = 8 )              ynewlo
   real    ( kind = 8 ), external :: ABIC
   
  !***********************************************************************************************************************************
  !start is the starting point for the iteration, and must contain the power transformed values if power transformation is required.
  !step determines size and shape of initial simplex and must reflect the units of the variables. step = stc * initial values. 
  !In case of power transformation step is calculated to be the same of the not tansformed case.
  !*********************************************************************************************************************************** 
  !!!
  !!!if ((maxval(cv_S%struct_par_opt).eq.0)) then            !Only if theta opimization is not required. Then they are not recalculated in the ABIC function!  
  !!!    call bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All) 
  !!!endif
  !!!
  start = d_S%struct_par_opt_vec  
  step  = stc * d_S%struct_par_opt_vec
  
  !Loop to power transform where necessary
  if ((maxval(cv_S%trans_theta).eq.1).or.(d_S%trans_sig.eq.1)) then !This means we have at least a theta or sig to power transform
    k = 0
    if ((maxval(cv_S%struct_par_opt).eq.1)) then   
      do i = 1, cv_PAR%p
        if (cv_S%struct_par_opt(i).eq.1)  then    
          do j = 1,cv_S%num_theta_type (i)
            k = k + 1
            if (cv_S%trans_theta(i).eq.1) then
              start(k) = cv_S%alpha_trans(i)*((d_S%struct_par_opt_vec(k)**(1./cv_S%alpha_trans(i)))-1) !Forward-trans
              step(k) = cv_S%alpha_trans(i)*((((stc+1)*d_S%struct_par_opt_vec(k))**(1./cv_S%alpha_trans(i)))-1) - &
                        cv_S%alpha_trans(i)*((d_S%struct_par_opt_vec(k)**(1./cv_S%alpha_trans(i)))-1)
            endif
          enddo
        endif
      enddo
    endif
    
    if ((d_S%sig_opt.eq.1).and.(d_S%trans_sig.eq.1)) then 
      k = k + 1
      start(k) = d_S%alpha_trans_sig *((d_S%struct_par_opt_vec(k)**(1./d_S%alpha_trans_sig))-1) !Forward-trans
      step(k) = d_S%alpha_trans_sig *((((stc+1)*d_S%struct_par_opt_vec(k))**(1./d_S%alpha_trans_sig))-1) - &
                d_S%alpha_trans_sig *((d_S%struct_par_opt_vec(k)**(1./d_S%alpha_trans_sig))-1)
    endif
  endif 
  
  !***********************************************************************************************************************************
  !At this point start and step are ready and power transformed where and if required
  
   konvge = 1
   reqmin = cv_A%structural_conv*1e-3 !1.0D-2 !Need to be tested
   if (cv_A%Q_compression_flag.ne.0) nQ0 = cv_PAR%p
        
   call nelmin_sp(ABIC,n,start,xmin,ynewlo,reqmin,step,konvge,cv_A%it_max_structural,icount,numres,ifault, & 
        & d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S, d_S, d_PM, cv_PAR,cv_PM,nQ0)
        
   if (ifault.eq.2) then !Maximum number of iterations has exceeded --> Warning               
       write(retmsg,10) it_bga   
10     format('Warning: Maximum number of iterations exceeded in structural parameter optimization procedure during ebaPEST iteration',i4, & 
                 & '. The vector that gives the minimum ABIC during the procedure will be considered.') 
       call utl_writmess(6,retmsg) 
   endif
   
   cv_S%str_obj_fun = ynewlo ! Minimum value of the objective function. It is the value of ABIC
   
   !****************************************************************************************************************************************
   !------ Here we reform d_S%theta and d_S%sig overwriting the elements optimized by Nelder-Mead and leaving unchanged the others. --------
   !---------- The values that minimize the structural parameter objective function are also assigned to d_S%struct_par_opt_vec ------------
   !--------------------- In case of power transformation, the structural parameters are back-transformed. --------------------------------- 
   !------------------- At the end, d_S%struct_par_opt_vec, d_S%theta and d_S%sig are in the physical space. -------------------------------     
   !Form d_S%struct_par_opt_vec
   d_S%struct_par_opt_vec = xmin !May include sigma
   !Reform d_S%theta and d_S%sig
   k = 0
   if ((maxval(cv_S%struct_par_opt).eq.1)) then   
     do i = 1, cv_PAR%p
       if (cv_S%struct_par_opt(i).eq.1) then
         if (cv_S%trans_theta(i).eq.1) then
           do j = 1,cv_S%num_theta_type (i)
              k = k + 1
              d_S%struct_par_opt_vec(k) = ((d_S%struct_par_opt_vec(k) + cv_S%alpha_trans(i)) / (cv_S%alpha_trans(i)))**(cv_S%alpha_trans(i)) !Back-trans
              d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
           enddo
         else
           do j = 1,cv_S%num_theta_type (i)
             k = k + 1
             d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
           enddo
         endif  
       endif
     enddo
   endif
   
   if (d_S%sig_opt.eq.1) then
     k=k+1
     if (d_S%trans_sig.eq.1) then 
       d_S%struct_par_opt_vec(k) = ((d_S%struct_par_opt_vec(k) + d_S%alpha_trans_sig) / (d_S%alpha_trans_sig))**(d_S%alpha_trans_sig) !Back-trans
       d_S%sig = d_S%struct_par_opt_vec(k)
     else
       d_S%sig = d_S%struct_par_opt_vec(k)
     endif
   endif
   !Note: from here d_S%theta and d_S%sig contain the theta and sigma values optimized after the minimization of the structural parameter obj. function.
   !It's not true that the CURRENT Qss, Qsy, HQsy, Qyy and other variables, that depend on theta and sigma, are calculated with these optimized parameters.
   !This because in Nelder Mead the values that minimize the function can occur not during the last performed iteration.    
   !*******************************************************************************************************************************************     
   
 end subroutine ABIC_struct_param_optim


subroutine beg_str_object_fun(cv_OBS,d_OBS,d_A,cv_S,d_PM,cv_PAR,cv_PM)


   implicit none
   
   type(cv_observ),      intent(in)     :: cv_OBS
   type(d_observ),       intent(in)     :: d_OBS
   type(d_algorithmic),  intent(inout)  :: d_A 
   type(d_prior_mean),   intent(in)     :: d_PM
   type(cv_param),       intent(in)     :: cv_PAR 
   type(cv_struct),      intent(inout)  :: cv_S 
   type(cv_prior_mean),  intent(in)     :: cv_PM
   
   integer                              :: errcode,i, curr_struct 
   double precision, allocatable        :: z(:), pivot(:)
   double precision                     :: lndetGyy, ztiGyyz
   double precision, allocatable        :: UinvGyy(:,:) ! used as both U and InvGyy
   double precision, allocatable        :: Gyy(:,:)
   double precision, allocatable        :: HXB(:), TMPV(:)  
   double precision, allocatable        :: HXQbb(:,:)
   double precision, allocatable        :: OMEGA(:,:)
   
   
   allocate (pivot(cv_OBS%nobs))
   allocate (z(cv_OBS%nobs))
   allocate (UinvGyy(cv_OBS%nobs,cv_OBS%nobs))
   allocate (Gyy(cv_OBS%nobs,cv_OBS%nobs))
   
   
   !----- intitialize variables
     lndetGyy  = 0.D0
     ztiGyyz   = 0.D0
     UinvGyy   = UNINIT_REAL   ! matrix
     Gyy       = UNINIT_REAL   ! matrix
  
   ! At this point we need HXB, HXQbb, OMEGA only if we have prior information about beta. 
   if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta   
      !Calculate HXB
      allocate(HXB(cv_OBS%nobs))
      HXB = UNINIT_REAL ! -- matrix (nobs)
      call DGEMV('n',cv_OBS%nobs, cv_PAR%p, 1.D0, d_A%HX, cv_OBS%nobs, &
              d_PM%beta_0, 1, 0.D0, HXB,1)
      !Form the linearization-corrected residuals and deallocate HXB no more necessary.
      z = d_OBS%obs - d_OBS%h + d_A%Hsold - HXB
      if (allocated(HXB))  deallocate(HXB)
      !Form HXQbb
      allocate(HXQbb(cv_OBS%nobs,cv_PAR%p))
      call dgemm('n', 'n', cv_OBS%nobs, cv_PAR%p,  cv_PAR%p, 1.D0, d_A%HX, &
              cv_OBS%nobs, d_PM%Qbb, cv_PAR%p, 0.D0, HXQbb, cv_OBS%nobs)
      !Form OMEGA = HXQbb x (HX)' and deallocate HXQbb no more necessary.
      allocate(OMEGA(cv_OBS%nobs,cv_OBS%nobs))
      call dgemm('n', 't', cv_OBS%nobs, cv_OBS%nobs,  cv_PAR%p, 1.D0, HXQbb, &
              cv_OBS%nobs, d_A%HX,  cv_OBS%nobs, 0.D0, OMEGA, cv_OBS%nobs)
      if (allocated(HXQbb))  deallocate(HXQbb)
      !Form Gyy = Qyy + OMEGA and deallocate OMEGA no more necessary.
      Gyy = OMEGA + d_A%Qyy
      if (allocated(OMEGA))  deallocate(OMEGA)
   else !------> we don't have prior information about beta
     !Form the linearization residuals z  
     z = d_OBS%obs - d_OBS%h + d_A%Hsold
     !Form Gyy = Qyy 
     Gyy = d_A%Qyy
   endif  
   !***********************************************************************************************************
   
   !*******************************************************************************   
   !--- Calculate the determinant term of the str. pars objective function --------
   !----------------------------- 0.5ln(det(Gyy)) ---------------------------------
   !-- First perform LU decomposition on Gyy 
   UinvGyy = Gyy !-nobs x nobs --- note that this is used as U in this context
   call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode)
   if (allocated(pivot))  deallocate(pivot)
   do i = 1,cv_OBS%nobs                             
     lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i)))  
   enddo                                            
   lndetGyy = 0.5 * lndetGyy
   !*******************************************************************************
   
   !*******************************************************************************   
   !------------ Calculate the misfit term of the objective function --------------
   !------------------------ 0.5(z' x invGyy x z) ---------------------------------
   !Calculate the inverse of Gyy
   UinvGyy = Gyy   ! nobs x nobs, re-use UinvGyy, now as InvGyy
   call INVGM(cv_OBS%nobs,UinvGyy)
   !Form inv(Gyy)*z
   allocate(TMPV(cv_OBS%nobs))
     call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
            z, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
   !Multiply z' * TMPV and 0.5 
   call DGEMV('t',cv_OBS%nobs, 1, 5.0D-1, z, cv_OBS%nobs, &
          TMPV, 1, 0.D0, ztiGyyz,1)
   if (allocated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
   !*******************************************************************************
   
   if (allocated(z))       deallocate(z)
   if (allocated(UinvGyy)) deallocate(UinvGyy)
   if (allocated(Gyy))     deallocate(Gyy)
   
   !****************************************************************************** 
   !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
   cv_S%str_obj_fun = lndetGyy + ztiGyyz
   !******************************************************************************
   
end subroutine beg_str_object_fun





subroutine beg_ABIC(d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM)

   use bayes_pest_control
   use model_input_output
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use jupiter_input_data_support
   use param_trans
   
   implicit none
   
   type(kernel_XQR),     intent(in)      :: d_XQR
   type(cv_observ),      intent(in)     :: cv_OBS
   type(d_observ),       intent(in)     :: d_OBS
   type(d_algorithmic),  intent(inout)  :: d_A 
   type(cv_algorithmic), intent(inout)  :: cv_A
   type(d_prior_mean),   intent(in)     :: d_PM
   type(cv_param),       intent(in)     :: cv_PAR 
   type(d_param),        intent(inout)  :: d_PAR        
   type(cv_struct),      intent(inout)     :: cv_S 
   type(d_struct),       intent(inout)  :: d_S
   type(Q0_compr),       intent(in)     :: Q0_All(:)
   type(cv_prior_mean),  intent(in)     :: cv_PM
   
   integer                              :: errcode, i
   double precision, allocatable        :: pivot(:)
   double precision                     :: lndetGyy, ztiGyyz
   double precision, allocatable        :: UinvGyy(:,:) ! used as both U and InvGyy
   double precision, allocatable        :: TMPV(:)  
   
   allocate (pivot(cv_OBS%nobs))
   allocate (UinvGyy(cv_OBS%nobs,cv_OBS%nobs))

      
   !----- intitialize variables
     errcode   = UNINIT_INT
     lndetGyy  = 0.D0
     ztiGyyz   = 0.D0
     UinvGyy   = UNINIT_REAL   ! matrix

   
   !!!!!!****************************************************************************************************************************
   !!!!!! At this point d_S%theta, d_S%sig and d_S%struct_par_opt_vec are ready to be used in the calculation of the obj func.
   !!!!!!---------------------------------- d_S%struct_par_opt_vec is in the physical space. ----------------------------------------
   !!!!!!****************************************************************************************************************************
   !!!!!
   !!!!!!We need to recalculate Qss Qsy HQsy only if at least one theta optimization is required. Otherwise, if the optimization
   !!!!!!call bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All) 
   !!!!!!We need to recalculate Qyy.
   !!!!!!call bmo_form_Qyy(d_XQR, d_S, cv_OBS, d_A)
     
     
   !*******************************************************************************   
   !--- Calculate the determinant term of the ABIC --------------------------------
   !------------------------------- ln(det(Gyy)) ----------------------------------
   !-- First perform LU decomposition on Gyy 
   UinvGyy = d_A%Qyy !-nobs x nobs --- note that this is used as U in this context
   call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode)
   if (allocated(pivot))       deallocate(pivot)
   do i = 1,cv_OBS%nobs                             
     lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i)))  
   enddo                                            
   !*******************************************************************************
   
   !*******************************************************************************   
   !------------ Calculate the misfit term of the ABIC -------------------------
   !------------------------ (z' x invGyy x z) ---------------------------------
   !Calculate the inverse of Gyy
   UinvGyy = d_A%Qyy   ! nobs x nobs, re-use UinvGyy, now as InvGyy
   call INVGM(cv_OBS%nobs,UinvGyy)
   !Form inv(Gyy)*z
   allocate(TMPV(cv_OBS%nobs))
     call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
            d_OBS%linobs, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
   !Multiply z' * TMPV  
   call DGEMV('t',cv_OBS%nobs, 1, 1.D0, d_OBS%linobs, cv_OBS%nobs, &
          TMPV, 1, 0.D0, ztiGyyz,1)
   if (allocated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
   !*******************************************************************************

   
   if (allocated(UinvGyy)) deallocate(UinvGyy)
 
   
   !****************************************************************************** 
   !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
   cv_S%str_obj_fun = cv_OBS%nobs*log(2.*3.14159265) + lndetGyy + ztiGyyz + 2*cv_S%num_theta_opt 
   !******************************************************************************
end subroutine beg_ABIC






end module struct_param_optimization


!****************************************************************************************************************************
!---------------- External function that contains the structural parameters objective function to minimize ------------------
!****************************************************************************************************************************

real (kind = 8) function ABIC(str_par_opt_vec,d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,n,nQ0)

   use bayes_pest_control
   use model_input_output
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use jupiter_input_data_support
   
   implicit none
   
   integer                              :: nQ0, n
   type(kernel_XQR),     intent(in)      :: d_XQR
   type(cv_observ),      intent(in)     :: cv_OBS
   type(d_observ),       intent(in)     :: d_OBS
   type(d_algorithmic),  intent(inout)  :: d_A 
   type(cv_algorithmic), intent(inout)  :: cv_A
   type(d_prior_mean),   intent(in)     :: d_PM
   type(cv_param),       intent(in)     :: cv_PAR 
   type(d_param),        intent(inout)  :: d_PAR        
   type(cv_struct),      intent(in)     :: cv_S 
   type(d_struct),       intent(inout)  :: d_S
   type(Q0_compr),       intent(in)     :: Q0_All(nQ0)
   type(cv_prior_mean),  intent(in)     :: cv_PM
   
   integer                              :: errcode, i, j, k, munit,kk,kj
   double precision                     :: str_par_opt_vec (n)!This must be the pars vector to be optimized for. Must be the first argument in ABIC (NelMead requires this) 
   double precision, allocatable        :: pivot(:)
   double precision                     :: lndetGyy, ztiGyyz
   double precision, allocatable        :: UinvGyy(:,:) ! used as both U and InvGyy
   double precision, allocatable        :: TMPV(:)  
   
   allocate (pivot(cv_OBS%nobs))
   allocate (UinvGyy(cv_OBS%nobs,cv_OBS%nobs))

      
   !----- intitialize variables
     errcode   = UNINIT_INT
     lndetGyy  = 0.D0
     ztiGyyz   = 0.D0
     UinvGyy   = UNINIT_REAL   ! matrix

     
     
   !******************************************************************************************************
   !First we need to reform d_S%theta and d_S%sig overwriting the elements that must be optimized for and  
   !leaving unchanged the others. These elements are into str_par_opt_vec (passed by NelMead). 
   !d_S%theta and d_S%sig are used during the matrix operations.(Qss Qsy HQsy Qyy)
   !*****************************************************************************************************
   !In case of power transformation str_par_opt_vec may contains values in the estimation space. The value
   !will be back-transformed and assigned to d_S%struct_par_opt_vec.  
   !*****************************************************************************************************
   d_S%struct_par_opt_vec = str_par_opt_vec
   k = 0
   if ((maxval(cv_S%struct_par_opt).eq.1)) then   
     do i = 1, cv_PAR%p
       if (cv_S%struct_par_opt(i).eq.1) then
         if (cv_S%trans_theta(i).eq.1) then
           do j = 1,cv_S%num_theta_type (i)
              k = k + 1
              d_S%struct_par_opt_vec(k) = ((str_par_opt_vec(k) + cv_S%alpha_trans(i)) / (cv_S%alpha_trans(i)))**(cv_S%alpha_trans(i))
              d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
           enddo
         else
           do j = 1,cv_S%num_theta_type (i)
             k = k + 1
             d_S%theta(i,j) = d_S%struct_par_opt_vec(k)
           enddo
         endif  
       endif
     enddo
   endif
   
   if (d_S%sig_opt.eq.1) then
     k=k+1
     if (d_S%trans_sig.eq.1) then 
       d_S%struct_par_opt_vec(k) = ((str_par_opt_vec(k) + d_S%alpha_trans_sig) / (d_S%alpha_trans_sig))**(d_S%alpha_trans_sig)
       d_S%sig = d_S%struct_par_opt_vec(k)
     else
       d_S%sig = d_S%struct_par_opt_vec(k)
     endif
   endif
   !****************************************************************************************************************************
   !****************************************************************************************************************************
   
   !****************************************************************************************************************************
   ! At this point d_S%theta, d_S%sig and d_S%struct_par_opt_vec are ready to be used in the calculation of the obj func.
   !---------------------------------- d_S%struct_par_opt_vec is in the physical space. ----------------------------------------
   !****************************************************************************************************************************
   !We need to recalculate Qss Qsy HQsy only if at least one theta optimization is required. Otherwise, if the optimization
   !is just for sig, we need to recalculate just Qyy
   if ((maxval(cv_S%struct_par_opt).eq.1)) then            !Only if theta opimization is required.  
      call bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All) 
   endif
   !Here Qss Qsy HQsy are ready, recalculated if necessary. We need to recalculate Qyy.
   call bmo_form_Qyy(d_XQR, d_S, cv_OBS, d_A)
   
  
   !*******************************************************************************   
   !--- Calculate the determinant term of the ABIC --------------------------------
   !------------------------------- ln(det(Gyy)) ----------------------------------
   !-- First perform LU decomposition on Gyy 
   UinvGyy = d_A%Qyy !-nobs x nobs --- note that this is used as U in this context
   
       !!!!!! munit = utl_nextunit()
       !!!!!!open (munit,file='fileqyy.txt')
       !!!!!!do kk=1,289
       !!!!!!write (munit,"(97ES14.7)") (d_A%Qyy(kk,kj), kj=1,97)
       !!!!!!enddo
       !!!!!!close(munit)
       
       call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode)
   if (allocated(pivot))       deallocate(pivot)
   do i = 1,cv_OBS%nobs                             
     lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i)))  
   enddo                                            
   !*******************************************************************************
   
   !*******************************************************************************   
   !------------ Calculate the misfit term of the ABIC -------------------------
   !------------------------ (z' x invGyy x z) ---------------------------------
   !Calculate the inverse of Gyy
   UinvGyy = d_A%Qyy   ! nobs x nobs, re-use UinvGyy, now as InvGyy
   call INVGM(cv_OBS%nobs,UinvGyy)
   !Form inv(Gyy)*z
   allocate(TMPV(cv_OBS%nobs))
     call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
            d_OBS%linobs, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
   !Multiply z' * TMPV  
   call DGEMV('t',cv_OBS%nobs, 1, 1.D0, d_OBS%linobs, cv_OBS%nobs, &
          TMPV, 1, 0.D0, ztiGyyz,1)
   if (allocated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
   !*******************************************************************************

   
   if (allocated(UinvGyy)) deallocate(UinvGyy)
 
   
   !****************************************************************************** 
   !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
   ABIC = cv_OBS%nobs*log(2.*3.14159265) + lndetGyy + ztiGyyz + 2*cv_S%num_theta_opt 
   !******************************************************************************
   
return
end function ABIC