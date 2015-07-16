subroutine condition_kalman(mean_d,mean_obs,var_d,var_obs,m,dim_obs,n_u,&
                dim_t,no_pars,cor_len,input_dim,lambda_dim,hyperparam,&
                design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod
IMPLICIT NONE

type(linear_model_data) :: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: mean_d(m*n_u,dim_t,3),var_d(m*n_u,m*n_u,dim_t,3)
real :: mean_obs(dim_obs*n_u,dim_t,2),var_obs(dim_obs*n_u,dim_obs*n_u,dim_t,2)
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: design_data(n_u*dim_obs,dim_t),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)

shallow_water%m=m
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%n_used=n_u
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim      
call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)
shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)
if (dim_obs>1) then 
        shallow_water%k_level=hyperparam(lambda_dim+2*input_dim+1)
endif

call condition_kalman_sub(shallow_water)
mean_d=shallow_water%states_means
mean_obs=shallow_water%states_means_obs
var_d=shallow_water%states_variances
var_obs=shallow_water%states_variances_obs

END subroutine

subroutine evaluate_kalman(emulated_output,mean_d,mean_obs,var_d,var_obs,m,&
        dim_obs,n_u,dim_t,no_pars,cor_len,input_dim,lambda_dim,hyperparam,param,&
        design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod

IMPLICIT NONE

type(linear_model_data) :: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: mean_d(m*n_u,dim_t,3),var_d(m*n_u,m*n_u,dim_t,3)
real :: mean_obs(dim_obs*n_u,dim_t,2),var_obs(dim_obs*n_u,dim_obs*n_u,dim_t,2)
real :: emulated_output(dim_obs,dim_obs,dim_t,2)
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: param(no_pars)
real :: design_data(n_u*dim_obs,dim_t),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)


 
shallow_water%m=m
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%n_used=n_u
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim      
call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)
shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)
if (dim_obs>1) then 
        shallow_water%k_level=hyperparam(lambda_dim+2*input_dim+1)
endif
shallow_water%parameters(n_u+1,:)=param


!Allocate some initial vectors.

IF (ALLOCATED(shallow_water%states_means)) THEN
DEALLOCATE(shallow_water%states_means)
ENDIF
ALLOCATE(shallow_water%states_means(m*n_u,dim_t,3)) 
shallow_water%states_means=0  

IF (ALLOCATED(shallow_water%states_variances)) THEN
DEALLOCATE(shallow_water%states_variances)
ENDIF
ALLOCATE(shallow_water%states_variances(m*n_u,m*n_u,dim_t,3))
shallow_water%states_variances=0

IF (ALLOCATED(shallow_water%states_means_obs)) THEN
DEALLOCATE(shallow_water%states_means_obs)
ENDIF
ALLOCATE(shallow_water%states_means_obs(dim_obs*n_u,dim_t,2))
shallow_water%states_means_obs=0

IF (ALLOCATED(shallow_water%states_variances_obs)) THEN
DEALLOCATE(shallow_water%states_variances_obs)
ENDIF
ALLOCATE(shallow_water%states_variances_obs(dim_obs*n_u,dim_obs*n_u,dim_t,2))
shallow_water%states_variances_obs=0



shallow_water%states_means=mean_d
shallow_water%states_means_obs=mean_obs
shallow_water%states_variances=var_d
shallow_water%states_variances_obs=var_obs
call evaluate_kalman_sub(shallow_water)
emulated_output(1,:,:,1)=shallow_water%e_obs_total(1,:,:,1)
emulated_output(:,:,:,2)=shallow_water%var_e_obs_total(1,:,:,:)

END subroutine

subroutine condition_nonkalman(z_prime,m,dim_obs,n_u,dim_t,no_pars,cor_len,input_dim,&
        lambda_dim,hyperparam,design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: z_prime(n_u*dim_t*m)
real :: design_data(n_u,dim_t*dim_obs),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)

shallow_water%m=m
shallow_water%n_used=n_u
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim      

call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)

shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)
if (dim_obs>1) then 
        shallow_water%k_level=hyperparam(lambda_dim+2*input_dim+1)
endif

call write_g_to_file(shallow_water)
call set_z(shallow_water)
call set_sigma_tilde(shallow_water)
call set_sigma(shallow_water)
call set_z_prime(shallow_water)
z_prime=shallow_water%z_prime

END subroutine condition_nonkalman

subroutine evaluate_nonkalman(emulated_output,z_prime,m,dim_obs,&
                n_u,dim_t,no_pars,cor_len,input_dim,&
                lambda_dim,hyperparam,param,design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: z_prime(m*n_u*dim_t)
real :: param(no_pars)
real :: emulated_output(dim_obs*dim_t)
real :: design_data(n_u,dim_t*dim_obs),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)

shallow_water%m=m
shallow_water%n_used=n_u
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim  

 
call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)

shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)
if (dim_obs>1) then 
        shallow_water%k_level=hyperparam(lambda_dim+2*input_dim+1)
endif
shallow_water%parameters(n_u+1,:)=param
shallow_water%z_prime=z_prime

call set_y(shallow_water)
emulated_output=shallow_water%y
END subroutine evaluate_nonkalman

subroutine evaluate_nonkalman_variance(variance_output,m,dim_obs,&
                n_u,dim_t,no_pars,cor_len,input_dim,&
                lambda_dim,hyperparam,param,design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: param(no_pars)
real :: variance_output(dim_obs*dim_t,dim_obs*dim_t)
real :: design_data(n_u,dim_t*dim_obs),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)

shallow_water%m=m
shallow_water%n_used=n_u
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim  
 
call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)

shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)
if (dim_obs>1) then 
        shallow_water%k_level=hyperparam(lambda_dim+2*input_dim+1)
endif
shallow_water%parameters(n_u+1,:)=param

call write_g_to_file_var(shallow_water)
call write_g_to_file(shallow_water)
call set_z(shallow_water)
call set_sigma_tilde(shallow_water)
call set_sigma(shallow_water)
variance_output=calc_variance(shallow_water)
END subroutine evaluate_nonkalman_variance


subroutine ll_nonkalman(ll,m,dim_obs,n_u,dim_t,no_pars,cor_len,input_dim,&
        lambda_dim,hyperparam,design_data,design_pars,rain,pars_physical,v_ini,e_ini)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: shallow_water

integer :: m,n_u,dim_t,dim_obs,no_pars,input_dim,lambda_dim
real :: hyperparam(2*input_dim+lambda_dim+dim_obs-1),cor_len,v_ini,e_ini
real :: design_data(n_u,dim_t),design_pars(n_u,no_pars)
real :: rain(dim_t)
real :: pars_physical(5)
real :: ll

shallow_water%m=m
shallow_water%n_used=n_u
shallow_water%dim_obs=dim_obs
shallow_water%t_max=dim_t
shallow_water%no_of_pars=no_pars
shallow_water%cor_factor=cor_len      
shallow_water%V_ini_inp=v_ini      
shallow_water%e_ini_inp=e_ini      
shallow_water%input_dim=input_dim      
shallow_water%lambda_dim=lambda_dim      

call read_data(shallow_water)
call resample_pars(shallow_water)
shallow_water%output=design_data
shallow_water%parameters(1:n_u,:)=design_pars
shallow_water%input=rain
shallow_water%parameters_physical=pars_physical
call allocate_initial(shallow_water)

shallow_water%k_lam=hyperparam(1:lambda_dim)
shallow_water%k_loss=hyperparam(lambda_dim+1:lambda_dim+input_dim)
shallow_water%k_delay=hyperparam(lambda_dim+input_dim+1&
        :lambda_dim+2*input_dim)

call write_g_to_file(shallow_water)
call set_z(shallow_water)
call set_sigma_tilde(shallow_water)
call set_sigma(shallow_water)
ll=log_likelihood(shallow_water)

END subroutine ll_nonkalman



subroutine condition_kalman_sub(this)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: this
!indices, r for row, c for column
integer :: i,j,ir,ic
integer :: m,dim_obs,dim_t,n_used,n_pars,n_test

! design data variables
real,ALLOCATABLE :: var_V_D(:,:), var_V_1(:,:), H_D(:,:),&
  A(:,:), F_D(:,:),&
  g_D(:),pars_des(:,:),observations_des(:,:)
real :: beta(this%no_of_pars)

dim_obs=this%dim_obs
dim_t=this%t_max
n_used=this%n_used
n_test=this%n_test_sets
n_pars=this%no_of_pars
m=this%m
beta=1.0/(getdelta(this)*this%cor_factor_multi)

  IF (ALLOCATED(observations_des)) THEN
    DEALLOCATE(observations_des)
  ENDIF
  ALLOCATE(observations_des(n_used*dim_obs,dim_t))
  observations_des=this%output(1:(n_used*dim_obs),:)



  IF (ALLOCATED(pars_des)) THEN
    DEALLOCATE(pars_des)
  ENDIF
  ALLOCATE(pars_des(1:n_used,n_pars))
  pars_des=this%parameters(1:n_used,:)

  IF (ALLOCATED(var_V_D)) THEN
    DEALLOCATE(var_V_D)
  ENDIF
  ALLOCATE(var_V_D(m*n_used,m*n_used))
  var_V_D=0

  IF (ALLOCATED(var_V_1)) THEN
    DEALLOCATE(var_V_1)
  ENDIF
  ALLOCATE(var_V_1(m*n_used,m*n_used))
  var_V_1=0

  do ir=1,n_used
    do ic=1,n_used
      var_V_D(m*(ir-1)+1:m*ir,m*(ic-1)+1:m*ic) = rho(this%Sigma_ini,beta,this%gamma,pars_des(ir,:),&
        pars_des(ic,:))
      var_V_1(m*(ir-1)+1:m*ir,m*(ic-1)+1:m*ic) = rho(this%V_ini,beta,this%gamma,pars_des(ir,:),&
        pars_des(ic,:))
    end do
  end do 

  IF (ALLOCATED(H_D)) THEN
    DEALLOCATE(H_D)
  ENDIF
  ALLOCATE(H_D(dim_obs*n_used,m*n_used))
  H_D=0
  do i=1,n_used
    H_D((i-1)*dim_obs+1:i*dim_obs,(i-1)*m+1:i*m)=getH(this,i)
  end do

  IF (ALLOCATED(this%states_means)) THEN
   DEALLOCATE(this%states_means)
  ENDIF
  ALLOCATE(this%states_means(m*n_used,dim_t,3)) 
  this%states_means=0  

  IF (ALLOCATED(this%states_variances)) THEN
    DEALLOCATE(this%states_variances)
  ENDIF
  ALLOCATE(this%states_variances(m*n_used,m*n_used,dim_t,3))
  this%states_variances=0

  IF (ALLOCATED(this%states_means_obs)) THEN
    DEALLOCATE(this%states_means_obs)
  ENDIF
  ALLOCATE(this%states_means_obs(dim_obs*n_used,dim_t,2))
  this%states_means_obs=0

  IF (ALLOCATED(this%states_variances_obs)) THEN
    DEALLOCATE(this%states_variances_obs)
  ENDIF
  ALLOCATE(this%states_variances_obs(dim_obs*n_used,dim_obs*n_used,dim_t,2))
  this%states_variances_obs=0

  IF (ALLOCATED(A)) THEN
    DEALLOCATE(A)
  ENDIF
  ALLOCATE(A(m*n_used,dim_obs*n_used))
  A=0

  ! initial 
  this%states_means(:,1,1) = rep_array(this%E_ini,n_used)
    this%states_means(:,1,3) = rep_array(this%E_ini,n_used)
    this%states_variances(:,:,1,1) = var_V_1

  !observation prediction
  this%states_means_obs(:,1,1) = MATMUL(H_d,this%states_means(:,1,1))
  this%states_means_obs(:,1,2) = MATMUL(H_d,this%states_means(:,1,3))
  this%states_variances_obs(:,:,1,1) = MATMUL(MATMUL(H_d,this%states_variances(:,:,1,1)),&
    TRANSPOSE(H_d))

  !observation prediction inverse

  this%states_variances(:,:,1,2) = inv_mat(this%states_variances(:,:,1,1))   
  this%states_variances_obs(:,:,1,2) = inv_mat(this%states_variances_obs(:,:,1,1))


  !filtering
  A = matmul(matmul(this%states_variances(:,:,1,1),transpose(H_d)),&
    this%states_variances_obs(:,:,1,2))
  this%states_means(:,1,2) = this%states_means(:,1,1) + matmul(A, observations_des(:,1)-&
    this%states_means_obs(:,1,1))
  this%states_variances(:,:,1,3) = this%states_variances(:,:,1,1) - matmul(matmul(A,H_d),&
    this%states_variances(:,:,1,1))

  IF (ALLOCATED(F_D)) THEN
    DEALLOCATE(F_D)
  ENDIF
  ALLOCATE(F_D(m*n_used,m*n_used))
  F_D=0
    do i=1,n_used
      F_D((i-1)*m+1:i*m, (i-1)*m+1:i*m)=getFexp(this,i)
    end do

  do j=2,dim_t

    IF (ALLOCATED(g_D)) THEN
      DEALLOCATE(g_D)
    ENDIF
    ALLOCATE(g_D(m*n_used))

    do i=1,n_used
      g_D((i-1)*m+1:i*m)=getk(this,i,j,F_D((i-1)*m+1:i*m, (i-1)*m+1:i*m))
    end do

    IF (ALLOCATED(H_D)) THEN
      DEALLOCATE(H_D)
    ENDIF
    ALLOCATE(H_D(dim_obs*n_used,m*n_used))
    H_D=0

    do i=1,n_used
      H_D((i-1)*dim_obs+1:i*dim_obs,(i-1)*m+1:i*m)=getH(this,i)
    end do

    this%states_means(:,j,1) = matmul(F_D(:,:),this%states_means(:,j-1,2)) + g_D
      this%states_means(:,j,3) = matmul(F_D(:,:),this%states_means(:,j-1,3)) + g_D
    this%states_variances(:,:,j,1) = matmul(matmul(F_D(:,:),this%states_variances(:,:,j-1,3)),&
      transpose(F_D(:,:))) + var_V_D

    !observation prediction
    this%states_means_obs(:,j,1) = MATMUL(H_d,this%states_means(:,j,1))
      this%states_means_obs(:,j,2) = MATMUL(H_d,this%states_means(:,j,3))
    this%states_variances_obs(:,:,j,1) = MATMUL(MATMUL(H_d,this%states_variances(:,:,j,1)),&
      TRANSPOSE(H_d))


  !     !observation prediction inverse
    
    this%states_variances_obs(:,:,j,2) = inv_mat(this%states_variances_obs(:,:,j,1))

      this%states_variances(:,:,j,2) = inv_mat(this%states_variances(:,:,j,1))   

  !     !filtering

    A = matmul(matmul(this%states_variances(:,:,j,1),transpose(H_d)),&
      this%states_variances_obs(:,:,j,2))
    this%states_means(:,j,2) = this%states_means(:,j,1) + matmul(A,(observations_des(:,j)-&
      this%states_means_obs(:,j,1)))
     this%states_variances(:,:,j,3) = this%states_variances(:,:,j,1) - matmul(matmul(A,H_d),&
       this%states_variances(:,:,j,1))

  end do

! write(*,*) this%states_variances(2,8,500,2)

END subroutine condition_kalman_sub

subroutine evaluate_kalman_sub(this)
USE lin_mod
USE service_functions
IMPLICIT NONE
type(linear_model_data):: this
integer :: i,j,i_e,k
integer :: m,dim_obs,dim_t,n_used,n_pars,n_test

! design data variables

real,ALLOCATABLE :: H_D(:,:),&
  A(:,:), F_D(:,:),states_means_smooth(:,:),&
  pars_des(:,:),pars_est(:,:),observations_des(:,:)
real :: beta(this%no_of_pars)

! emulation data variables----------------------------------------------------------

real,ALLOCATABLE :: states_means_e(:,:,:), states_variances_e(:,:,:,:),&
  states_covariances_e(:,:,:,:),pars_e(:),states_means_e_smooth(:,:)
real,ALLOCATABLE :: F_E(:,:,:),Cov_V_DE(:,:),g_E(:),&
  Var_V_E(:,:),d(:,:)

! smoothing variabless
real,ALLOCATABLE :: E_DE(:),E_DE_allY(:),E_DE_prevY(:),M1(:,:),M2(:,:),M3(:,:),Mm(:,:)

! smoothing variables for variance
real,ALLOCATABLE :: var_DE(:,:),var_DE_prevY(:,:),var_DE_allY(:,:),&
  states_variances_smooth(:,:,:),states_covariances_smooth(:,:,:),&
  states_variances_e_smooth(:,:,:)

! observations
real,ALLOCATABLE :: e_obs(:,:,:), H_e(:,:)

!estimation variables, will be allocated only if it is ON!

dim_obs=this%dim_obs
dim_t=this%t_max
n_used=this%n_used
n_test=this%n_test_sets
n_pars=this%no_of_pars
m=this%m
beta=1.0/(getdelta(this)*this%cor_factor_multi)


!I may not do the conditioning, but i still need some variables from that part
  IF (ALLOCATED(observations_des)) THEN
      DEALLOCATE(observations_des)
    ENDIF
    ALLOCATE(observations_des(n_used*dim_obs,dim_t))
    observations_des=this%output(1:(n_used*dim_obs),:)
    IF (ALLOCATED(pars_des)) THEN
      DEALLOCATE(pars_des)
    ENDIF
    ALLOCATE(pars_des(1:n_used,n_pars))
    pars_des=this%parameters(1:n_used,:)
    IF (ALLOCATED(H_D)) THEN
      DEALLOCATE(H_D)
    ENDIF
    ALLOCATE(H_D(dim_obs*n_used,m*n_used))
    H_D=0
    do i=1,n_used
      H_D((i-1)*dim_obs+1:i*dim_obs,(i-1)*m+1:i*m)=getH(this,i)
    end do
    IF (ALLOCATED(F_D)) THEN
    DEALLOCATE(F_D)
   ENDIF
    ALLOCATE(F_D(m*n_used,m*n_used))
    F_D=0
      do i=1,n_used
        F_D((i-1)*m+1:i*m, (i-1)*m+1:i*m)=getFexp(this,i)
      end do
!end of part


IF (ALLOCATED(pars_est)) THEN
  DEALLOCATE(pars_est)
ENDIF
ALLOCATE(pars_est(n_test,n_pars))
pars_est=this%parameters(n_used+1:n_used+n_test,:)

!allocate file for writing all the results
IF (ALLOCATED(this%e_obs_total)) THEN
  DEALLOCATE(this%e_obs_total)
ENDIF
ALLOCATE(this%e_obs_total(n_test,dim_obs,dim_t,3))
this%e_obs_total=0

IF (ALLOCATED(this%var_e_obs_total)) THEN
  DEALLOCATE(this%var_e_obs_total)
ENDIF
ALLOCATE(this%var_e_obs_total(n_test,dim_obs,dim_obs,dim_t))
this%var_e_obs_total=0

do i_e=1,n_test

  IF (ALLOCATED(pars_e)) THEN
    DEALLOCATE(pars_e)
  ENDIF
  ALLOCATE(pars_e(n_pars))
  pars_e=0

  pars_e=pars_est(i_e,:)

  IF (ALLOCATED(states_means_e)) THEN
    DEALLOCATE(states_means_e)
  ENDIF
  ALLOCATE(states_means_e(m,dim_t,3))
  states_means_e=0

  IF (ALLOCATED(states_variances_e)) THEN
    DEALLOCATE(states_variances_e)
  ENDIF
  ALLOCATE(states_variances_e(m,m,dim_t,2))
  states_variances_e=0

  IF (ALLOCATED(states_covariances_e)) THEN
    DEALLOCATE(states_covariances_e)
  ENDIF
  ALLOCATE(states_covariances_e(m*n_used,m,dim_t,2))
  states_covariances_e=0

  IF (ALLOCATED(H_D)) THEN
    DEALLOCATE(H_D)
    ENDIF
    ALLOCATE(H_D(dim_obs*n_used,m*n_used))
    H_D=0
    do i=1,n_used
      H_D((i-1)*dim_obs+1:i*dim_obs,(i-1)*m+1:i*m)=getH(this,i)
  end do

  ! initial 
  states_means_e(:,1,1) = this%E_ini
  states_means_e(:,1,3) = this%E_ini
  states_variances_e(:,:,1,1) = this%V_ini
  do i=1,n_used
    states_covariances_e((i-1)*m+1:i*m,:,1,1) = rho(this%V_ini,beta,this%gamma,pars_e,pars_des(i,:))
  end do
  !     !filtering

  IF (ALLOCATED(A)) THEN
    DEALLOCATE(A)
  ENDIF
  ALLOCATE(A(m,dim_obs*n_used))
  A=0

  A = matmul(matmul(transpose(states_covariances_e(:,:,1,1)),transpose(H_d)),&
    this%states_variances_obs(:,:,1,2))
   states_means_e(:,1,2) = states_means_e(:,1,1) + matmul(A,(observations_des(:,1)-&
    this%states_means_obs(:,1,1)))
   states_variances_e(:,:,1,2) = states_variances_e(:,:,1,1) - matmul(matmul(A,H_d),&
     states_covariances_e(:,:,1,1))
   states_covariances_e(:,:,1,2)=states_covariances_e(:,:,1,1) - matmul(this%states_variances(:,:,1,1),&
    matmul(transpose(H_d),matmul(this%states_variances_obs(:,:,1,2),matmul(H_d,states_covariances_e(:,:,1,1)))))

  IF (ALLOCATED(F_E)) THEN
    DEALLOCATE(F_E)
  ENDIF
  ALLOCATE(F_E(dim_t,m,m))
  F_E=0
  F_E(1,:,:)=getFexp(this,n_used+i_e)

  do j=2,dim_t
    F_E(j,:,:)=F_E(1,:,:)

    IF (ALLOCATED(g_E)) THEN
      DEALLOCATE(g_E)
    ENDIF
    ALLOCATE(g_E(m))
    g_E(:)=getk(this,n_used+i_e,j,F_E(j,:,:))

    IF (ALLOCATED(Cov_V_DE)) THEN
      DEALLOCATE(Cov_V_DE)
    ENDIF
    ALLOCATE(Cov_V_DE(m*n_used,m))
    Cov_V_DE=0
    do k=1,n_used
      Cov_V_DE((k-1)*m+1:k*m,:)=rho(this%Sigma_ini,beta,this%gamma,pars_des(k,:),pars_e)
    end do

    IF (ALLOCATED(Var_V_E)) THEN
      DEALLOCATE(Var_V_E)
    ENDIF
    ALLOCATE(Var_V_E(m,m))
    Var_V_E=0
    Var_V_E=rho(this%Sigma_ini,beta,this%gamma,pars_e,pars_e)

    states_means_e(:,j,1) = matmul(F_E(j,:,:),states_means_e(:,j-1,2)) + g_E 
      states_means_e(:,j,3) = matmul(F_E(j,:,:),states_means_e(:,j-1,3)) + g_E 
    states_variances_e(:,:,j,1) = matmul(F_E(j,:,:),matmul(states_variances_e(:,:,j-1,2),transpose(F_E(j,:,:))))+&
      Var_V_E
    states_covariances_e(:,:,j,1) = matmul(F_D(:,:),matmul(states_covariances_e(:,:,j-1,2),&
      transpose(F_E(j,:,:)))) + Cov_V_DE

  !filtering of emulator data
    A = matmul(matmul(transpose(states_covariances_e(:,:,j,1)),transpose(H_d)),&
    this%states_variances_obs(:,:,j,2))
   states_means_e(:,j,2) = states_means_e(:,j,1) + matmul(A,(observations_des(:,j)-&
    this%states_means_obs(:,j,1)))
   states_variances_e(:,:,j,2) = states_variances_e(:,:,j,1) - matmul(matmul(A,H_d),&
     states_covariances_e(:,:,j,1))
   states_covariances_e(:,:,j,2)=states_covariances_e(:,:,j,1) - matmul(this%states_variances(:,:,j,1),&
    matmul(transpose(H_d),matmul(this%states_variances_obs(:,:,j,2),matmul(H_d,states_covariances_e(:,:,j,1)))))

  end do

  !smoothing


  IF (ALLOCATED(states_means_e_smooth)) THEN
  DEALLOCATE(states_means_e_smooth)
  ENDIF
  ALLOCATE(states_means_e_smooth(m,dim_t))
  states_means_e_smooth=0

  IF (ALLOCATED(states_means_smooth)) THEN
  DEALLOCATE(states_means_smooth)
  ENDIF
  ALLOCATE(states_means_smooth(m*n_used,dim_t))
  states_means_smooth=0

    IF (ALLOCATED(states_variances_smooth)) THEN
      DEALLOCATE(states_variances_smooth)
    ENDIF
    ALLOCATE(states_variances_smooth(m*n_used,m*n_used,dim_t))
    states_variances_smooth=0

    IF (ALLOCATED(states_variances_e_smooth)) THEN
      DEALLOCATE(states_variances_e_smooth)
    ENDIF
    ALLOCATE(states_variances_e_smooth(m,m,dim_t))
    states_variances_e_smooth=0

    IF (ALLOCATED(states_covariances_smooth)) THEN
      DEALLOCATE(states_covariances_smooth)
    ENDIF
    ALLOCATE(states_covariances_smooth(m*n_used,m,dim_t))
    states_covariances_smooth=0

  states_means_e_smooth(:,dim_t)=states_means_e(:,dim_t,2)
  states_means_smooth(:,dim_t)=this%states_means(:,dim_t,2)
    states_variances_e_smooth(:,:,dim_t)=states_variances_e(:,:,dim_t,2)
    states_variances_smooth(:,:,dim_t)=this%states_variances(:,:,dim_t,2)
    states_covariances_smooth(:,:,dim_t)=states_covariances_e(:,:,dim_t,2)

  IF (ALLOCATED(M1)) THEN
  DEALLOCATE(M1)
  ENDIF
  ALLOCATE(M1(m*(n_used+1),m*(n_used+1)))
  M1=0

  IF (ALLOCATED(M2)) THEN
  DEALLOCATE(M2)
  ENDIF
  ALLOCATE(M2(m*(n_used+1),m*(n_used+1)))
  M2=0

  IF (ALLOCATED(M3)) THEN
  DEALLOCATE(M3)
  ENDIF
  ALLOCATE(M3(m*(n_used+1),m*(n_used+1)))
  M3=0

  IF (ALLOCATED(mm)) THEN
  DEALLOCATE(mm)
  ENDIF
  ALLOCATE(mm(m*(n_used+1),m*(n_used+1)))
  mm=0

  DO j=dim_t-1,1,-1
    M1(1:m*n_used,1:m*n_used) = this%states_variances(:,:,dim_t,3)
    M1(1:m*n_used,m*n_used+1:m*(n_used+1)) = states_covariances_e(:,:,dim_t,2)
    M1(m*n_used+1:m*(n_used+1),1:m*n_used) = transpose(states_covariances_e(:,:,dim_t,2))
    M1(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1)) = states_variances_e(:,:,dim_t,2)

    M2(1:m*n_used,1:m*n_used) = F_D(:,:)
!   M2(1:m*n_used,m*n_used+1:m*(n_used+1)) = states_covariances_e(:,:,dim_t,2)
!   M2(m*n_used+1:m*(n_used+1),1:m*n_used) = transpose(states_covariances_e(:,:,dim_t,2))
    M2(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1)) = F_E(j+1,:,:)

    IF (ALLOCATED(a)) THEN
      DEALLOCATE(a)
    ENDIF
    ALLOCATE(a(m*n_used,m*n_used))
    a=0

    IF (ALLOCATED(d)) THEN
      DEALLOCATE(d)
    ENDIF
    ALLOCATE(d(m,m))
    d=0

    IF (ALLOCATED(E_DE)) THEN
      DEALLOCATE(E_DE)
    ENDIF
    ALLOCATE(E_DE(m*(n_used+1)))
    E_DE=0

    IF (ALLOCATED(E_DE_prevY)) THEN
      DEALLOCATE(E_DE_prevY)
    ENDIF
    ALLOCATE(E_DE_prevY(m*(n_used+1)))
    E_DE_prevY=0

    IF (ALLOCATED(E_DE_allY)) THEN
      DEALLOCATE(E_DE_allY)
    ENDIF
    ALLOCATE(E_DE_allY(m*(n_used+1)))
    E_DE_allY=0


      IF (ALLOCATED(var_DE)) THEN
        DEALLOCATE(var_DE)
      ENDIF
      ALLOCATE(var_DE(m*(n_used+1),m*(n_used+1)))
      var_DE=0

      IF (ALLOCATED(var_DE_prevY)) THEN
        DEALLOCATE(var_DE_prevY)
      ENDIF
      ALLOCATE(var_DE_prevY(m*(n_used+1),m*(n_used+1)))
      var_DE_prevY=0

      IF (ALLOCATED(var_DE_allY)) THEN
        DEALLOCATE(var_DE_allY)
      ENDIF
      ALLOCATE(var_DE_allY(m*(n_used+1),m*(n_used+1)))
      var_DE_allY=0

    d=inv_mat(states_variances_e(:,:,j+1,1)-matmul(matmul(transpose(states_covariances_e(:,:,j+1,1)),&
      this%states_variances(:,:,j+1,2)),states_covariances_e(:,:,j+1,1)))

    a=this%states_variances(:,:,j+1,2)+matmul(matmul(matmul(matmul(this%states_variances(:,:,j+1,2),&
      states_covariances_e(:,:,j+1,1)),d),transpose(states_covariances_e(:,:,j+1,1))),this%states_variances(:,:,j+1,2))

    M3(1:m*n_used,1:m*n_used) = a
    M3(1:m*n_used,m*n_used+1:m*(n_used+1)) = -matmul(matmul(this%states_variances(:,:,j+1,2),&
      states_covariances_e(:,:,j+1,1)),d)
    M3(m*n_used+1:m*(n_used+1),1:m*n_used) = -matmul(matmul(d,&
      transpose(states_covariances_e(:,:,j+1,1))),this%states_variances(:,:,j+1,2))
    M3(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1)) = d


    Mm=matmul(matmul(M1,transpose(M2)),M3)


    E_DE(1:m*n_used)=this%states_means(:,j,2)
    E_DE(m*n_used+1:m*(n_used+1))=states_means_e(:,j,2)
    E_DE_allY(1:m*n_used)=states_means_smooth(:,j+1)
    E_DE_allY(m*n_used+1:m*(n_used+1))=states_means_e_smooth(:,j+1)
    E_DE_prevY(1:m*n_used)=this%states_means(:,j+1,1)
    E_DE_prevY(m*n_used+1:m*(n_used+1))=states_means_e(:,j+1,1)
    E_DE=E_DE+matmul(Mm,E_DE_allY- E_DE_prevY )
    states_means_smooth(:,j)=E_DE(1:m*n_used)
    states_means_e_smooth(:,j)=E_DE(m*n_used+1:m*(n_used+1))
      var_DE(1:m*n_used,1:m*n_used)=states_variances_smooth(:,:,j)
      var_DE(m*n_used+1:m*(n_used+1),1:m*n_used)=transpose(states_covariances_smooth(:,:,j))
      var_DE(1:m*n_used,m*n_used+1:m*(n_used+1))=states_covariances_smooth(:,:,j)
      var_DE(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1))=states_variances_e_smooth(:,:,j)
      var_DE_prevY(1:m*n_used,1:m*n_used)=this%states_variances(:,:,j+1,1)
      var_DE_prevY(m*n_used+1:m*(n_used+1),1:m*n_used)=transpose(states_covariances_e(:,:,j+1,1))
      var_DE_prevY(1:m*n_used,m*n_used+1:m*(n_used+1))=states_covariances_e(:,:,j+1,1)
      var_DE_prevY(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1))=states_variances_e(:,:,j+1,1)
      var_DE_allY(1:m*n_used,1:m*n_used)=states_variances_smooth(:,:,j+1)
      var_DE_allY(m*n_used+1:m*(n_used+1),1:m*n_used)=transpose(states_covariances_smooth(:,:,j+1))
      var_DE_allY(1:m*n_used,m*n_used+1:m*(n_used+1))=states_covariances_smooth(:,:,j+1)
      var_DE_allY(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1))=states_variances_e_smooth(:,:,j+1)
      ! var_DE=var_DE + matmul(matmul(Mm,var_DE_prevY- var_DE_allY),transpose(Mm))
      var_DE=var_DE - matmul(matmul(Mm,var_DE_prevY- var_DE_allY),transpose(Mm))
      ! pozor tady

      states_variances_e_smooth(:,:,j)=var_DE(m*n_used+1:m*(n_used+1),m*n_used+1:m*(n_used+1))
      states_variances_smooth(:,:,j)=var_DE(1:m*n_used,1:m*n_used)
      states_covariances_smooth(:,:,j)=var_DE(1:m*n_used,m*n_used+1:m*(n_used+1))
  end do
  !observation

  IF (ALLOCATED(H_e)) THEN
    DEALLOCATE(H_e)
  ENDIF
  ALLOCATE(H_e(dim_obs,m))
    H_e=getH(this,i_e+n_used)

!e obs contains also the variance
IF (ALLOCATED(e_obs)) THEN
  DEALLOCATE(e_obs)
ENDIF
ALLOCATE(e_obs(dim_obs,dim_t,3))
e_obs=0
do j=1,dim_t
  e_obs(:,j,1)=matmul(H_e,states_means_e_smooth(:,j))
end do
do i=1,dim_t
  this%e_obs_total(i_e,:,i,1)=e_obs(:,i,1)
end do
  e_obs(:,:,2)=matmul(H_e,states_means_e(:,:,3))
  do i=1,dim_t
    this%e_obs_total(i_e,:,i,2)=e_obs(:,i,2)
     this%var_e_obs_total(i_e,:,:,i)=&
       matmul(H_e,matmul(states_variances_e_smooth(:,:,i),transpose(H_e)))
  end do

  do i=1,dim_t
    this%e_obs_total(i_e,:,i,3)=e_obs(:,i,3)
  end do

end do


END subroutine evaluate_kalman_sub


