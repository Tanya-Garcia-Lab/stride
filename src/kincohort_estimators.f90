subroutine kincohort_estimators(n,q,x,delta,t,qvs,p,m,r,lim,rat,setting,d,H0,&
     num_estimators,boot,bootvar,usetruth,useOLS,useWLS,useEFF,useNPMLEs,hts0,Fest,var_est,eflag)
  implicit none
  character(len=10),intent(in) :: setting
  logical,intent(in) :: bootvar,H0,usetruth,useOLS,useWLS,useEFF,useNPMLEs
  integer,intent(in) :: n,p,m,boot,num_estimators
  integer,dimension(m),intent(in) :: r
  double precision,intent(in) :: t,lim,rat,d
  double precision,dimension(n), intent(in) :: x,delta
  double precision,dimension(p,m),intent(in) :: qvs
  double precision,dimension(p,n),intent(in) :: q
  double precision,dimension(p,num_estimators),intent(out) :: Fest
  double precision,dimension(p,p,num_estimators),intent(out) :: var_est  
  integer,intent(out) :: eflag

  double precision,dimension(n) :: ws,wn,gtl,gt,st,ht
  double precision,dimension(m,n) :: xsub,dsub,sts,gts,hts
  double precision,dimension(m),intent(out) :: hts0
  double precision,dimension(m+1) :: pQ
  double precision,dimension(p,p,n) :: qq
  double precision,dimension(p,p) :: a1,K
  integer,dimension(n) :: id1
  integer,dimension(m,n) :: id
  integer,dimension(m) :: rt
  integer,dimension(m+1):: rc
 
  call hatG(n,x,delta,gt)
  call hatGl(n,x,delta,gtl)
  call hatS(n,x,delta,st,ht)
  call subinfo(p,n,m,qvs,q,x,delta,r,xsub,dsub,sts,gts,hts,ht,id1,id)
  !print*,'hi there'
  if(useOLS) then
     call ols(p,n,t,q,x,delta,gtl,st,Fest(:,1),var_est(:,:,1),ws,qq,a1)
     call olsaug(p,n,m,t,q,x,delta,gtl,st,id,r,qq,a1,Fest(:,2),var_est(:,:,2))
  end if
  if(useWLS) then
     call wls(p,n,t,q,x,delta,ws,gtl,st,Fest(:,4),var_est(:,:,4),wn,qq,a1)
     call wlsaug(p,n,m,t,q,x,delta,gtl,st,id,r,wn,qq,a1,Fest(:,5),var_est(:,:,5))
  end if
  if(useEFF) then
	call effest(p,n,m,t,qvs,q,x,delta,gtl,st,r,K,Fest(:,7),var_est(:,:,7),eflag)
     if (eflag.ne.1) then
        goto 20
     end if
    call effestaug(p,n,m,t,qvs,q,x,delta,gtl,st,r,id,K,Fest(:,8),var_est(:,:,8))
  end if
  call imput(p,n,m,bootvar,boot,t,q,x,delta,gt,id,r,lim,rat,setting,d,H0,usetruth,useOLS,useWLS,useEFF,&
       qvs,Fest(:,3),var_est(:,:,3),Fest(:,6),&
       var_est(:,:,6),Fest(:,9),var_est(:,:,9))
  if(useNPMLES) then
     call bootnpmles(bootvar,boot,n,p,t,m,qvs,q,x,delta,r,hts0,Fest(:,10),var_est(:,:,10),&
          Fest(:,11),var_est(:,:,11))
     !print*,'hts0out=',hts0
  end if
20  return  
end subroutine kincohort_estimators


! set values: qvs, lim, t,rat
subroutine Rsetvalues(setting,p,m,qvs,lim,rat)
  implicit none
  character(len=10),intent(in) :: setting
  integer, intent(in) :: p,m
  double precision,dimension(p,m),intent(out) :: qvs
  double precision, intent(out) :: lim,  rat
  double precision :: tmp,delta
  integer :: i
  select case(setting)
  case('mypenedata')
     ! mimic penetrance data
     lim=100.
     qvs(1,1)=0.
     qvs(1,2)=1.
     qvs(1,3)=0.5
     qvs(2,:)=1.-qvs(1,:)
  case('myrealdata')
     ! mimic huntington disease data
     lim = 100
     ! beginning: mimic the huntington's disease qvs 
     qvs(1,1)=0. 
     qvs(1,2)=0.5 
     qvs(1,3)=0.97 
     qvs(1,4)=0.25 
     qvs(1,m-1)=1. 
     qvs(1,m)=0.75 
     qvs(2,:)=1.-qvs(1,:) 
     ! end: mimic the hungtinton's disease qvs
  case default
     ! setting for cox and non-cox data
     lim = 10
     ! beginning:the wang and robinowitz set up of qvs 
     tmp=0.2 
     qvs(1,1)=1 
     qvs(1,2)=(1.+tmp)/2 
     qvs(1,3)=tmp 
     qvs(1,4)=1+tmp*tmp/4.-tmp/2.-3./4. 
     !  qvs(1,5)=-tmp*tmp/4.+3./4.*tmp+1./2. 
     !  qvs(1,6)=-tmp*tmp/4.+tmp 
     qvs(2,:)=1.-qvs(1,:) 
     ! end: the wang and robinowitz set up of qvs 
     ! begining:mimic the kin-cohort qvs 
     !  qvs(1,1)=0.8 
     !  qvs(1,2)=0.4 
     !  qvs(1,3)=0.1 
     !  qvs(2,:)=1.-qvs(1,:) 
     ! end:mimic the kin-cohort qvs 
  end select
  rat = 4.
  return
end subroutine Rsetvalues

! set values: qvs, lim, t,rat
subroutine setvalues(setting,fullrange,numt,tval,&
     p,m,qvs,lim,t,rat,timeval)
  implicit none
  character(len=10),intent(in) :: setting
  logical,intent(in) :: fullrange
  double precision,intent(in) :: tval
  integer, intent(in) :: p,m,numt
  double precision,dimension(p,m),intent(out) :: qvs
  double precision,dimension(numt),intent(out) :: timeval
  double precision, intent(out) :: lim, t, rat
  double precision :: tmp,delta
  integer :: i
  select case(setting)
  case('myrealdata')
     ! mimic huntington disease data
     lim = 100
     ! beginning: mimic the huntington's disease qvs 
     qvs(1,1)=0. 
     qvs(1,2)=0.5 
     qvs(1,3)=0.97 
     qvs(1,4)=0.25 
     qvs(1,m-1)=1. 
     qvs(1,m)=0.75 
     qvs(2,:)=1.-qvs(1,:) 
     ! end: mimic the hungtinton's disease qvs
     t=70
     delta =1.0
  case default
     ! setting for cox and non-cox data
     lim = 10
     ! beginning:the wang and robinowitz set up of qvs 
     tmp=0.2 
     qvs(1,1)=1 
     qvs(1,2)=(1.+tmp)/2 
     qvs(1,3)=tmp 
     qvs(1,4)=1+tmp*tmp/4.-tmp/2.-3./4. 
     !  qvs(1,5)=-tmp*tmp/4.+3./4.*tmp+1./2. 
     !  qvs(1,6)=-tmp*tmp/4.+tmp 
     qvs(2,:)=1.-qvs(1,:) 
     ! end: the wang and robinowitz set up of qvs 
     ! begining:mimic the kin-cohort qvs 
     !  qvs(1,1)=0.8 
     !  qvs(1,2)=0.4 
     !  qvs(1,3)=0.1 
     !  qvs(2,:)=1.-qvs(1,:) 
     ! end:mimic the kin-cohort qvs 
     delta=0.1
     t=2.5
  end select
  rat = 4.
  if(fullrange) then
     select case(setting)
     case('myrealdata')
        do i =1,numt
           ! original
           !timeval(i) =2*tval
           ! testing, 07/07/2012
           timeval(i) = 2*i
        end do
     case default
        do i=1,numt
           !timeval(i) = 2*i
           timeval(i)= delta*i 
           ! for getting extra calculations, to get endpoints
           !timeval(1) = tval
           !timeval(2) = 2.*tval
        end do
     end select
  else
     timeval(1)=t
  end if
  return
end subroutine setvalues






! subroutine to create rc and rt values
!subroutine rcrtvalues(iseed,n,m,rc,rt)
!  !use m_inssor  !! need this to run in gfortran (uncommented to run in R)
!  !use IFPORT
!  implicit none
!  integer,intent(in) :: iseed,n,m
!  integer,dimension(m),intent(out) :: rt
!  integer,dimension(m+1),intent(out) :: rc
!  double precision,dimension(m+1) :: pQ
!  integer :: i,j
!  double precision :: tmp
! 
!  ! initialize random number generator at iseed
!  tmp=rand(iseed)
!
! ! print*,'iseed=',iseed, 'tmp=',tmp
!  do i=2,m
!     rc(i)=int(rand(0)*(n-m))
!  end do
!  rc(1)=0
!  rc(m+1)=n
!  pQ=rc
!  j= m+1
!  call inssor(pQ)
!  do i=2,m
!     rc(i)=int(pQ(i))+i-1
!  end do
!  do i=1,m
!     rt(i)=rc(i+1)-rc(i)
!  end do
!  return
!end subroutine rcrtvalues




! organize data
subroutine orgdata(n,p,m,q,x,delta,r,qvs)
  implicit none
  integer, intent(in) :: n,p,m
  double precision,dimension(p,n),intent(in) :: q
  double precision,dimension(n),intent(in) :: x,delta
  integer,dimension(m),intent(out) :: r
  double precision,dimension(n) :: fail,uniqueID,gender,origin,t
  double precision,dimension(p,m),intent(out) :: qvs
  integer :: i,j,kflag,m1
  double precision :: tol
  tol=1e-6
  kflag=2
  !sort data in increasing order
  do i=1,p
     t=x
     call dsort(t,q(i,:),n,kflag)
  end do
  call dsort(x,delta,n,kflag)
  ! get the new r,m,qvs
  r=0
  m1=0
  qvs=0
  i=1
1 if (i.le.n) then
     j=1
2    if (j.le.m1) then
        if (sum(abs(qvs(:,j)-q(:,i))).lt.tol) then
           r(j)=r(j)+1
!           print*,'j=',j,'r=',r
           i=i+1
           goto 1
        else
           j=j+1
           goto 2
        end if
     else
        m1=m1+1
 !       print*,'m1=',m1
        r(m1)=1
        qvs(:,m1)=q(:,i)
        i=i+1
        goto 1
     end if
  end if
  return
end subroutine orgdata




! 2/19/2013: CDF function for new cure rate model
! for 0 <= t <=100
function Ct1(a,t)
  implicit none
  double precision,intent(in) :: a,t
  double precision :: Ct1
  Ct1=a / (1.+exp(-(t-80.)/5.))
  return
end function Ct1

! 2/19/2013: CDF function for new cure rate model
! for 100 <= t <=300
function Ct2(a,t)
  implicit none
  double precision, intent(in) :: a,t
  double precision :: Ct2, Ct1
  double precision :: alpha,beta,ttmp
  ttmp=100.
  beta = (1.-Ct1(a,ttmp))/200.
  alpha = 1.-300.*beta
  Ct2 = alpha + beta *t
  return
end function Ct2

! 2/19/2013: pdf function for new cure rate model
! for 0 <= t <=100
function Ct1pdf(a,t)
  implicit none
  double precision, intent(in) :: a,t
  double precision :: Ct1pdf,Ct1
  Ct1pdf = (a/5.) * exp(-(t-80.)/5.)/ &
       (1.+exp(-(t-80.)/5.))**2.
  return
end function Ct1pdf

! 2/19/2013: CDF function for new cure rate model
! for 100 <= t <=300
function Ct2pdf(a,t)
  implicit none
  double precision, intent(in) :: a,t
  double precision :: Ct2pdf, Ct1
  double precision :: alpha,beta,ttmp
  ttmp=100.
  beta = (1.-Ct1(a,ttmp))/200.
  Ct2pdf = beta
  return
end function Ct2pdf

! 2/19/2013: inverse function for new cure rate model
! for 0 <= t <=100
function Ct1inv(a,u)
  implicit none
  double precision,intent(in) :: a,u
  double precision :: Ct1inv
  Ct1inv = -5.*log(a/u-1.)+80.
  return
end function Ct1inv

! 2/19/2013: inverse function for new cure rate model
! for 100 <= t <=300
function Ct2inv(a,u)
  implicit none
  double precision, intent(in) :: a,u
  double precision :: Ct2inv, Ct1
  double precision :: alpha,beta,ttmp
  ttmp=100.
  beta = (1.-Ct1(a,ttmp))/200.
  alpha = 1.-300.*beta
  Ct2inv = (u-alpha)/beta
  return
end function Ct2inv


!the function Ut={1-exp(-t/2)} / {1-exp(-5d/2)} 
function Ut(d,t)
  implicit none
  double precision,intent(in) :: d,t
  double precision :: Ut
  Ut=(1.-exp(-t/2.))/(1.-exp(-5.*d/2.))
  return
end function Ut

!the function Utpdf=0.5 * exp(-t/2)} / {1-exp(-5d/2)}
function Utpdf(d,t)
  implicit none
  double precision,intent(in) :: d,t
  double precision :: Utpdf
  Utpdf=0.5 * exp(-t/2.)/(1.-exp(-5.*d/2.))
  return
end function Utpdf

!the function Utinv=-2 * log(1-a*(1-exp(-5d/2)))
function Utinv(d,a)
  implicit none
  double precision,intent(in) :: d,a
  double precision :: Utinv, tmp
  tmp = 1.-exp(-5.*d/2.)
  Utinv=-2. * log( 1.-a*tmp)
  return
end function Utinv


!the function Vt={1-exp(-t/4)} / {1-exp(-5d/4)}
function Vt(d,t)
  implicit none
  double precision,intent(in) :: d,t
  double precision :: Vt
  Vt=(1.-exp(-t/4.))/(1.-exp(-5.*d/4.))
  return
end function Vt

!the function Vtpdf=0.25 * exp(-t/4)} / {1-exp(-5d/4)}
function Vtpdf(d,t)
  implicit none
  double precision,intent(in) :: d,t
  double precision :: Vtpdf
  Vtpdf=0.25 * exp(-t/4.)/(1.-exp(-5.*d/4.))
  return
end function Vtpdf

!the function Vtinv=-4 * log(1-a*(1-exp(-5d/4)))
function Vtinv(d,a)
  implicit none
  double precision,intent(in) :: d,a
  double precision :: Vtinv, tmp
  tmp = 1.-exp(-5.*d/4.)
  Vtinv=-4. * log( 1.-a*tmp)
  return
end function Vtinv

! the function Et=(1-exp(-t/b))/(1-exp(-10./b))
function Et(b,t)
  implicit none
  double precision,intent(in) :: b,t
  double precision :: Et
  Et=(1.-exp(-t/b))/ (1-exp(-10./b))
  return
end function Et

! the function Etpdf=(1/b)*exp(-t/b))/(1-exp(-10./b))
function Etpdf(b,t)
  implicit none
  double precision,intent(in) :: b,t
  double precision :: Etpdf
  Etpdf=exp(-t/b)/ (1-exp(-10./b))/b
  return
end function Etpdf

! the function Etinv=-b*log(1-a*(1-exp(-10/b)))
function Etinv(b,a)
  implicit none
  double precision,intent(in) :: b,a
  double precision :: Etinv,tmp
  tmp=1.-exp(-10./b)
  Etinv=-b*log(1-a*tmp)
  return
end function Etinv


 
 
! generate random censor time 
subroutine genc(c,setting,d,H0,censorrate) 
  !use IFPORT
  character(len=10),intent(in) :: setting
  logical,intent(in) :: H0
  integer,intent(in) :: censorrate
  double precision,intent(in) :: d
  double precision,intent(out) :: c 
  double precision :: a,b,tau 
  b=rand(0) 
  tau=100.
  ! censoring for Cox PH setup 
  !  c=b*tau 

  if(censorrate.eq.0) then
     ! 0% censoring 
     c=-100 
  elseif(censorrate.eq.20) then
     ! 20% censoring 
     if (b.le.0.1) then 
        c=-100 
     else 
        b=rand(0) 
        select case(setting)
        case('newmodcure')
           !print*,'newmodcure'
           a=650.
        case('newcoxcase')
           a=10.
        case('regcoxcase')
           a=10.
        case('noncoxcase')
           a=8.
        case('regcoxcure')
           a=5.
        case('noncoxcure')
           a=3.
        case('regpowcase')
           if(d.le.1.25) then
              a=7.*d
           else if((d.ge.1.26).and.(d.le.1.45)) then
              a=6*d
           else if((d.ge.1.49).and.(d.le.1.55)) then
              a=5*d   
           else
              a=4.*d
           end if
        case('regpowcure')
           if(d.le.1.12) then
              a=3.5*d
           else if((d.ge.1.15).and.(d.le.1.35))then
              a=3.2*d
           else if((d.ge.1.37).and.(d.le.1.5)) then
              a=2.8*d
           else
              a=2.*d
           end if
        end select
        call trueGinvU(a,b,c) 
     end if
  elseif(censorrate.eq.40) then
     ! 40% censoring 
     if (b.le.0.1) then 
        c=-100 
     else 
        b=rand(0)
        select case(setting)
        case('newmodcure')
           a=350.
        case('newcoxcase')
           !print*,'newcoxcase'
           a=4.
        case('regcoxcase')
           a=4.
        case('noncoxcase')
           a=3.
        case default
           a=3.
        end select
        call trueGinvU(a,b,c) 
        ! a=3.0
        ! call trueGinvE(a,b,c) 
     end if
  elseif(censorrate.eq.50) then
     ! 50% censoring 
     if (b.le.0.1) then 
        c=-100 
     else 
        b=rand(0) 
        select case(setting)
        case('mypenedata')
           if(H0) then
              a=150.   
           else
              a=210.
           end if
        case('newmodcure')
           a=3.
        case('newcoxcase')
           a=3.
        case('regcoxcase')
           a=3.
        case('noncoxcase')
           a=2.4
        case('regcoxcure')
           a=1.1
        case('noncoxcure')
           a=0.7
        case('regpowcase')
           if(d.le.1.15) then
              a=2.5*d
           else if((d.ge.1.17).and.(d.le.1.29)) then
              a=2.3*d
           else if((d.ge.1.3).and.(d.le.1.35)) then
              a=2.*d
           else if((d.ge.1.4).and.(d.le.1.45)) then
              a=1.9*d
           else if((d.ge.1.5).and.(d.le.1.9)) then
              a=1.5*d
           else
              a=1.2*d
           end if
        case('regpowcure')
           if(d.le.1.05) then
              a=1*d
           else if((d.ge.0.9).and.(d.le.1.35)) then
              a=0.75*d
           else if((d.ge.1.4).and.(d.le.1.45)) then
              a=0.6*d
           else if((d.ge.1.49).and.(d.le.1.95)) then
              a=0.55*d
           else
              a=0.45*d
           end if
        end select
        call trueGinvU(a,b,c) 
        ! a=3.0 
        ! call trueGinvE(a,b,c) 
     end if
  elseif(censorrate.eq.65) then
     ! 65% censoring 
     b=rand(0) 
     select case(setting)
     case('myrealdata')
        c=b*tau
     case('regcoxcase')
        ! cox data 
        a=4.5 
     case('noncoxcase')
        ! non-cox data 
        a=2.
     end select
     call trueGinvU(a,b,c) 
  elseif(censorrate.eq.90) then
     ! 90% censoring 
     if (b.le.0.04) then 
        c=-100 
     else 
        b=rand(0) 
        a=0.2 
        call trueGinvE(a,b,c) 
     end if
  else   
     ! 98% censoring 
     if (b.le.0.01) then 
        c=-100 
     else 
        b=rand(0) 
        a=0.02 
        call trueGinvE(a,b,c) 
     end if
  end if
  return 
end subroutine genc 
 
!true G^{-1}(t) 
subroutine trueGinvU(a,t,ginvt) 
  double precision, intent(in) :: a,t 
  double precision, intent(out) :: ginvt 
!uniform  
  ginvt=(1.-t)*a 
  return 
end subroutine trueGinvU 
 
!true G^{-1}(t) 
subroutine trueGinvE(a,t,ginvt) 
  double precision, intent(in) :: a,t 
  double precision, intent(out) :: ginvt 
! exponential 
  ginvt=-log(t)*a 
  return 
end subroutine trueGinvE 
 
!Kaplan-Meier estimator of G(t) at x_i's 
subroutine hatG(n,x,delta,gt) 
  integer,intent(in) :: n 
  double precision,dimension(n),intent(in) :: x,delta 
  double precision,dimension(n),intent(out) :: gt 
  integer :: i 
  gt(1)=1.-(1.-delta(1))/n 
  do i=2,n 
     gt(i)=gt(i-1)*(1.-(1.-delta(i))/(n-i+1)) 
  end do 
  gt(n)=gt(n-1) !avoid 0 
  return 
end subroutine hatG 
 
!Kaplan-Meier estimator of G(t-) at x_i's,  
subroutine hatGl(n,x,delta,gt) 
  integer,intent(in) :: n 
  double precision,dimension(n),intent(in) :: x,delta 
  double precision,dimension(n),intent(out) :: gt 
  integer :: i 
  gt(1)=1. 
  gt(2)=1.-(1.-delta(1))/n 
  do i=2,(n-1) 
     gt(i+1)=gt(i)*(1.-(1.-delta(i))/(n-i+1)) 
  end do 
  return 
end subroutine hatGl 
 
!Kaplan-Meier estimator of event process S(t) and 1-S(t) at x_i's 
subroutine hatS(n,x,delta,st,ht) 
  integer,intent(in) :: n 
  double precision,dimension(n),intent(in) :: x,delta 
  double precision,dimension(n),intent(out) :: st,ht 
  integer :: i 
  st(1)=1.-delta(1)/n 
  do i=2,n 
     st(i)=st(i-1)*(1.-delta(i)/(n-i+1)) 
  end do 
  ht=1-st 
  st(n)=st(n-1) !avoid 0 
  ht(n)=ht(n-1) 
  return 
end subroutine hatS 
 
 
 
!get m subsample information 
subroutine subinfo(p,n,m,qvs,q,x,delta,r,xsub,dsub,sts,gts,hts,ht,id,ida) 
  implicit none
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(p,n),intent(in) :: q  
  double precision,dimension(n),intent(in) :: x,delta 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(m,n),intent(out) :: xsub,dsub,sts,gts,hts 
  double precision,dimension(n),intent(out) :: ht 
  integer,dimension(n),intent(out) :: id 
  integer,dimension(m,n),intent(out) :: ida 
  double precision,dimension(p) :: ys 
  double precision :: tol 
  integer :: i,j,cnt 
  tol=1e-6 
  xsub=0 
  dsub=0 
  do j=1,m 
     cnt=1 
     do i=1,n 
        if (sum(abs(q(:,i)-qvs(:,j))).le.tol) then 
           xsub(j,cnt)=x(i) 
           dsub(j,cnt)=delta(i) 
           id(i)=j 
           ida(j,cnt)=i 
           cnt=cnt+1 
        end if 
     end do 
     call hatS(r(j),xsub(j,:),dsub(j,:),sts(j,:),hts(j,:)) 
     !print*,'j=',j,'xsub=',xsub(j,:),'sts=',sts(j,:)
     call hatG(r(j),xsub(j,:),dsub(j,:),gts(j,:)) 
     do i=1,r(j) 
        ht(ida(j,i))=hts(j,i) 
     end do 
  end do 
  return  
end subroutine subinfo 
 
 
!get m subsample information at t, transform K-M to get estimate 
subroutine newsubtinfo(n,m,p,t,r,qvs,xsub,dsub,hts,ht0,Fa) 
  integer,intent(in) :: n,m,p 
  double precision,intent(in) :: t 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(m,n),intent(in) :: xsub,dsub,hts 
  double precision,dimension(m),intent(out) :: ht0 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p) :: va 
  double precision,dimension(p) :: ys 
  double precision,dimension(m) :: sig2,w,tmp 
  double precision,dimension(p,p) :: a1,A,Ainv,B,vv 
  double precision :: tol 
  integer :: i,j 
  tol=1e-6 
!  print*,'r=',r
  do j=1,m 
     ht0(j)=0 
     sig2(j)=0 
     i=1 
     !ht0(j)=hts(j,i)
     !ht0(j)=hts(j,r(j))
5    if (i.le.r(j)) then 
        if (t.ge.xsub(j,i)) then 
           ht0(j)=hts(j,i) 
           if (dsub(j,i).gt.0.5) then 
              sig2(j)=sig2(j)+1./(max((r(j)-i)*(r(j)-i+1.),tol)) 
           end if 
        end if 
        i=i+1 
        goto 5 
     end if 
   !  print*,'t=',t,'j=',j,'ht0=',ht0(j),'max x=',maxval(xsub(j,:))
     sig2(j)=sig2(j)*(1-ht0(j))**2. 
  end do 
  w=r  ! codes wls (weighted least squares), w=r 
!  w=1.  ! codes ols, equivalent to NPMLE 1 
  A=0 
  B=0 
  ys=0 
  do j=1,m 
     call mul(p,qvs(:,j),qvs(:,j),a1) 
     A=A+w(j)*a1 
     B=B+w(j)*sig2(j)*a1*w(j) 
     ys=ys+w(j)*ht0(j)*qvs(:,j)      ! ht0(j) is 1-\hat S(t)
  end do 
  !print*,'ht0=',ht0
  call inv(p,A,Ainv) 
  Fa=matmul(Ainv,ys) 
  va=matmul(matmul(Ainv,B),transpose(Ainv)) 
  return  
end subroutine newsubtinfo 
 
! generate bootstrap data 
subroutine genbootdata(n,p,om,oq,ox,odelta,qb,xb,deltab,r,m,qvs) 
  !use IFPORT
  implicit none 
  integer,intent(in) :: p,n,om 
  double precision,dimension(p,n),intent(in) :: oq 
  double precision,dimension(n),intent(in) :: ox,odelta 
  double precision,dimension(p,n),intent(out) :: qb 
  double precision,dimension(n),intent(out) :: xb,deltab 
  double precision,dimension(p,om),intent(out) :: qvs 
  integer,dimension(om),intent(out) :: r 
  integer,intent(out) :: m 
  double precision,dimension(n) :: t 
  integer :: kflag,i,ind,j 
  double precision :: a,tol 
  tol=1e-6 
  kflag=2 
  do i=1,n 
     ind=nint(rand(0)*n+0.5) 
     qb(:,i)=oq(:,ind) 
     xb(i)=ox(ind) 
     deltab(i)=odelta(ind) 
  end do
  !sort them in increasing order 
  do i=1,p 
     t=xb 
     call dsort(t,qb(i,:),n,kflag) 
  end do
  call dsort(xb,deltab,n,kflag) 
  ! get the new r,m,qvs 
  r=0 
  m=0 
  qvs=0 
  i=1 
1 if (i.le.n) then 
     j=1 
2    if (j.le.m) then 
        if (sum(abs(qvs(:,j)-qb(:,i))).lt.tol) then 
           r(j)=r(j)+1 
           i=i+1 
           goto 1 
        else 
           j=j+1 
           goto 2 
        end if
     else 
        m=m+1 
        r(m)=1 
        qvs(:,m)=qb(:,i) 
        i=i+1 
        goto 1 
     end if
  end if
  return 
end subroutine genbootdata


! calculates type I and type II NPMLE and variance via bootstrap
subroutine bootnpmles(bootvar,boot,n,p,t,m,qvs,q,x,delta,r,hts0out,F1,va1,F2,va2) 
  implicit none 
  logical,intent(in) :: bootvar
  integer,intent(in) :: n,p,m,boot 
  double precision,intent(in) :: t 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(p,n),intent(in) :: q  
  double precision,dimension(n),intent(in) :: x,delta 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(out) :: F1,F2
  double precision,dimension(p,p),intent(out) :: va1,va2  
  double precision,dimension(p,n) :: qb
  double precision,dimension(n) :: xb,deltab
  double precision,dimension(m,n) :: xsub,dsub,sts,gts,hts 
  double precision,dimension(n) :: ht 
  double precision,dimension(m),intent(out) :: hts0out 
  double precision,dimension(m) :: hts0
  integer,dimension(n) :: id1 
  integer,dimension(m,n) :: id 
  integer,dimension(m) :: r1 
  double precision,dimension(p,m) :: qvs1 
  double precision,dimension(p,boot) :: Ft1boot,Ft2boot 
  integer :: i,m1 
  ! type I NPMLE
  call subinfo(p,n,m,qvs,q,x,delta,r,xsub,dsub,sts,gts,hts,ht,id1,id) 
  call newsubtinfo(n,m,p,t,r,qvs,xsub,dsub,hts,hts0out,F1) 
  ! type II NPMLE
  call npmle2(p,n,q,x,delta,t,F2) 
  if(bootvar) then
     i=1 
10   if(i.le.boot) then 
        call genbootdata(n,p,m,q,x,delta,qb,xb,deltab,r1,m1,qvs1) 
        ! type I NPMLE
        call subinfo(p,n,m1,qvs1,qb,xb,deltab,r1,xsub,dsub,sts,gts,hts,ht,id1,id) 
        call newsubtinfo(n,m1,p,t,r1,qvs1,xsub,dsub,hts,hts0,Ft1boot(:,i)) 
        ! type II NPMLE
        call npmle2(p,n,qb,xb,deltab,t,Ft2boot(:,i)) 
        i=i+1 
        goto 10 
     end if
     call sampcov(p,boot,Ft1boot,va1)
     call sampcov(p,boot,Ft2boot,va2)
  else
     va1=0.
     va2=0.
  end if
  return 
end subroutine bootnpmles
 
 
subroutine ols(p,n,t,q,x,delta,gt,st,Fa,va,ws,qq,a1) 
  integer,intent(in) :: p,n 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va,a1 
  double precision,dimension(n),intent(out) :: ws 
  double precision,dimension(p,p,n),intent(out) :: qq 
  double precision,dimension(p,n) :: psi,phi 
  double precision,dimension(p,p) :: a,a2 
  double precision :: tmp,indicate 
  integer :: i 
!phi=(A'A)^{-1}{q_iI(s_i<t)-q_iq_i'F(t)} 
  do i=1,n 
     call mul(p,q(:,i),q(:,i),qq(:,:,i)) 
  end do 
  a=sum(qq,3) 
  a=a/n 
  call inv(p,a,a1) 
  psi=0 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t) 
        Fa=Fa+psi(:,i)/gt(i) 
        a=a+qq(:,:,i)/gt(i) 
     end if 
  end do 
  Fa=Fa/n 
  a=a/n 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
!pass out the initial weights and form the phi's 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     ws(i)=1./(tmp*(1.-tmp)) 
     if (delta(i).gt.0.5) then 
        phi(:,i)=psi(:,i)-q(:,i)*tmp 
     end if 
  end do 
  phi=matmul(a1,phi) 
  call getvar(p,n,phi,gt,st,delta,va) 
  return 
end subroutine ols 
 
!augmented ols 
subroutine olsaug(p,n,m,t,q,x,delta,gt,st,id,r,qq,a1,Fa,va) 
  implicit none 
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,p,n),intent(in) :: qq 
  double precision,dimension(p,p),intent(in) :: a1 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: psi,h,Bh,phi 
  double precision,dimension(p,p,n) :: g,Bf 
  double precision,dimension(p,p) :: a,a2 
  double precision :: indicate,tol,tmp 
  integer :: i 
  tol=1e-6 
  psi=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t) 
     end if 
  end do 
  call gethx(p,n,m,delta,psi,gt,r,id,h) 
  call B2hfu(p,n,m,x,psi,delta,gt,qq,r,id,Bh,Bf) 
  Bh=h-Bh 
  Bf=qq-Bf 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i)/gt(i) 
        a=a+qq(:,:,i)/gt(i) 
     else 
        Fa=Fa+Bh(:,i)/gt(i) 
        a=a+Bf(:,:,i)/gt(i) 
     end if 
  end do 
  Fa=Fa/n 
  a=a/n 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  phi=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        call inprod(p,q(:,i),Fa,tmp) 
        phi(:,i)=psi(:,i)-q(:,i)*tmp 
     end if 
  end do 
  phi=matmul(a1,phi) 
  call getvaraug(p,n,m,phi,gt,st,x,delta,id,r,va) 
  return 
end subroutine olsaug 
 
 
subroutine wls(p,n,t,q,x,delta,ws,gt,st,Fa,va,w,qq,a1) 
  integer,intent(in) :: p,n 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: ws 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va,a1 
  double precision,dimension(p,p,n),intent(out) :: qq 
  double precision,dimension(n),intent(out) :: w 
  double precision,dimension(p,n) :: psi,phi 
  double precision,dimension(p,p) :: a,a2 
  double precision,dimension(n) :: wold 
  double precision :: tmp,tol,indicate,rep 
  integer :: i,j 
!phi=(A'WA)^{-1}q_iw(i){I(s_i<t)-q_i'F(t)} 
  tol=1e-6 
  rep=1 
  w=ws 
20 do i=1,n 
     call mul(p,q(:,i),q(:,i),qq(:,:,i)) 
     qq(:,:,i)=qq(:,:,i)*w(i) 
  end do 
  a=sum(qq,3) 
  a=a/n 
  call inv(p,a,a1) 
  psi=0 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*w(i)*indicate(x(i),t) 
        Fa=Fa+psi(:,i)/gt(i) 
        a=a+qq(:,:,i)/gt(i) 
     end if 
  end do 
  Fa=Fa/n 
  a=a/n 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  do i=1,p 
     Fa(i)=max(min(Fa(i),1.-tol),tol) 
  end do 
  wold=w 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     w(i)=1./(tmp*(1.-tmp)) 
  end do 
  if ((sum(abs(w-wold))/n.gt.tol).and.(rep.lt.10)) then 
     rep=rep+1 
     goto 20 
  end if 
! we settle with this w 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     if (delta(i).gt.0.5) then 
        phi(:,i)=psi(:,i)-q(:,i)*tmp*w(i) 
     end if 
  end do 
  phi=matmul(a1,phi) 
  call getvar(p,n,phi,gt,st,delta,va) 
  return 
end subroutine wls 
 
!augmented wls 
subroutine wlsaug(p,n,m,t,q,x,delta,gt,st,id,r,w,qq,a1,Fa,va) 
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt,st,w 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,p,n),intent(in) :: qq 
  double precision,dimension(p,p),intent(in) :: a1 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: psi,h,Bh,phi 
  double precision,dimension(p,p,n) :: g,Bf 
  double precision,dimension(p,p) :: a,a2 
  double precision :: indicate,tol,tmp 
  integer :: i 
  tol=1e-6 
  psi=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t)*w(i) 
     end if 
  end do 
  call gethx(p,n,m,delta,psi,gt,r,id,h) 
  call B2hfu(p,n,m,x,psi,delta,gt,qq,r,id,Bh,Bf) 
  Bh=h-Bh 
  Bf=qq-Bf 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i)/gt(i) 
        a=a+qq(:,:,i)/gt(i) 
     else 
        Fa=Fa+Bh(:,i)/gt(i) 
        a=a+Bf(:,:,i)/gt(i) 
     end if 
  end do 
  Fa=Fa/n 
  a=a/n 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  phi=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        call inprod(p,q(:,i),Fa,tmp) 
        phi(:,i)=psi(:,i)-q(:,i)*tmp*w(i) 
     end if 
  end do 
  phi=matmul(a1,phi) 
  call getvaraug(p,n,m,phi,gt,st,x,delta,id,r,va) 
  return 
end subroutine wlsaug 
 
 
 
 
!the imputation method on ols, wls and teffest, using bootstrap 
subroutine imput(p,n,m,bootvar,boot,t,q,x,delta,gt,id,r,lim,rat,setting,d,H0,usetruth,&
     useOLS,useWLS,useEFF,qvs,Fa,Va,Fb,Vb,Fc,Vc) 
  implicit none 
  character(len=10),intent(in) :: setting
  logical,intent(in) :: H0,usetruth,useOLS,useWLS,useEFF
  integer,intent(in) :: boot,p,n,m
  logical,intent(in) :: bootvar
  double precision,intent(in) :: t,lim,rat,d
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(p),intent(out) :: Fa,Fb,Fc 
  double precision,dimension(p,p),intent(out) :: Va,Vb,Vc 
  double precision,dimension(p,n) :: qb
  double precision,dimension(n) :: xb,deltab
  double precision,dimension(p,n) :: psi,Bp1 
  double precision,dimension(n) :: gt1,w,ht 
  double precision,dimension(p,m) :: qvs1 
  integer,dimension(m) :: r1 
  double precision,dimension(p,p) :: a,a2,K 
  double precision,dimension(m,n) :: xsub,dsub,sts,gts,hts 
  integer,dimension(m,n) :: id1 
  integer,dimension(n) :: id2 
  double precision,dimension(p,boot) :: Fa1,Fb1,Fc1 
  double precision,dimension(p) :: tmpa,tmpb,tmpc 
  double precision :: indicate 
  integer :: i,m1,eflag 
  if(useOLS) then
     call getolsimputFa(p,n,m,t,q,x,delta,gt,id,r,Fa,w)
  end if
  if(useWLS) then 
     call getwlsimputFa(p,n,m,t,q,x,delta,gt,id,r,w,Fb)
  end if
  if(useEFF) then
     call geteffestimputFa(p,n,m,t,qvs,q,x,delta,gt,r,id,Fc,eflag) 
     !  print*,'eflag=',eflag
     !the above eflag is guaranteed to be 1
  end if
  if(bootvar) then
!     print*,'bootvar=',bootvar, 'boot=',boot
     i=1 
10   if (i.le.boot) then   
        call genbootdata(n,p,m,q,x,delta,qb,xb,deltab,r1,m1,qvs1) 
        call hatG(n,xb,deltab,gt1) 
        call subinfo(p,n,m1,qvs1,qb,xb,deltab,r1,xsub,dsub,sts,gts,hts,ht,id2,id1) 
        if(useEFF) then 
           call geteffestimputFa(p,n,m1,t,qvs1,qb,xb,deltab,gt1,r1,id1,Fc1(:,i),eflag) 
           if (eflag.ne.1) then 
              goto 10 
           end if
        end if
        if(useOLS) then
           call getolsimputFa(p,n,m1,t,qb,xb,deltab,gt1,id1,r1,Fa1(:,i),w)
        end if
        if(useWLS) then
           call getwlsimputFa(p,n,m1,t,qb,xb,deltab,gt1,id1,r1,w,Fb1(:,i)) 
           !         print*,'Fa1=',Fa1
        end if
        i=i+1 
        goto 10 
     end if
     tmpa=sum(Fa1,2)/boot 
     tmpb=sum(Fb1,2)/boot 
     tmpc=sum(Fc1,2)/boot 
     do i=1,boot 
        Fa1(:,i)=Fa1(:,i)-tmpa 
        Fb1(:,i)=Fb1(:,i)-tmpb 
        Fc1(:,i)=Fc1(:,i)-tmpc 
     end do
     if(useOLS) then
        Va=matmul(Fa1,transpose(Fa1))/boot 
     else
        Va=0.
     end if
     if(useWLS) then
        Vb=matmul(Fb1,transpose(Fb1))/boot
     else
        Vb=0.
     end if
     if(useEFF) then
        Vc=matmul(Fc1,transpose(Fc1))/boot
     else 
        Vc=0.
     end if
  else
     Va=0.
     Vb=0.
     Vc=0.
  end if
  return 
end subroutine imput
 
 
 
 
 
 
subroutine getolsimputFa(p,n,m,t,q,x,delta,gt,id,r,Fa,ws) 
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(n),intent(out) :: ws 
  double precision,dimension(p,n) :: psi,Bp1 
  double precision,dimension(p,p) :: a,a2 
  double precision,dimension(p,p,n) :: qq 
  double precision :: indicate,tmp 
  integer :: i 
  do i=1,n 
     call mul(p,q(:,i),q(:,i),qq(:,:,i)) 
  end do 
  psi=0 
  Fa=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t) 
     end if 
  end do 
  call gethx(p,n,m,delta,psi,gt,r,id,Bp1) 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i) 
     else 
        Fa=Fa+Bp1(:,i) 
     end if 
  end do 
  a=sum(qq,3) 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
!pass out the initial weights 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     ws(i)=1./(tmp*(1.-tmp)) 
  end do 
  return 
end subroutine getolsimputFa 
 
 
subroutine getwlsimputFa(p,n,m,t,q,x,delta,gt,id,r,w,Fa) 
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt,w 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p,n) :: qq 
  double precision,dimension(p,n) :: psi,Bp1,hs 
  double precision,dimension(p,p) :: a,a2 
  double precision :: indicate,tmp 
  integer :: i 
  do i=1,n 
     call mul(p,q(:,i),q(:,i),qq(:,:,i)) 
     qq(:,:,i)=qq(:,:,i)*w(i) 
  end do 
  psi=0 
  Fa=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t)*w(i) 
     end if 
  end do 
  call gethx(p,n,m,delta,psi,gt,r,id,Bp1) 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i) 
     else 
        Fa=Fa+Bp1(:,i) 
     end if 
  end do 
  a=sum(qq,3) 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  return 
end subroutine getwlsimputFa 


subroutine geteffestimputFa(p,n,m,t,qvs,q,x,delta,gt,r,id,Fa,eflag)
  integer, intent(in):: p,n,m
  double precision,dimension(p,n),intent(in) :: q
  double precision,dimension(n),intent(in) :: x,delta,gt
  double precision,intent(in) :: t
  double precision,dimension(p,m),intent(in) :: qvs
  integer,dimension(m),intent(in) :: r
  integer,dimension(m,n),intent(in) :: id
  double precision,dimension(p),intent(out) :: Fa
  integer,intent(out) :: eflag
  double precision,dimension(p,p) :: A1,K1,K2,K21,K
  double precision,dimension(m) :: uf
  double precision,dimension(p,n) :: psi,phi,Bp1
  double precision :: s0,step,tol,indicate,a,h
  integer :: i,j,iter,pts,n1
  eflag=1
  n1=0
  tol=1e-6
  h=(maxval(x)-minval(abs(x)))*(1.*n/m)**(-0.2)*2.5
  K1=0
  K2=0
  pts=10000
  step=(maxval(x)-minval(abs(x))+2.*h*n/m)/pts
!  step=(maxval(x)-minval(abs(x)))/pts
!  step=maxval(x)*3./pts
!  do iter=0,pts-1
  do iter=0,pts
    ! s0=minval(abs(x))+(iter+0.1)*step
     s0=minval(abs(x))-h*n/m+(iter+0.5)*step
    ! s0=(iter+0.5)*step
     call getAinv(p,n,m,h,s0,r,q,x,delta,gt,qvs,A1,uf)
     if (minval(uf).gt.tol) then
        K1=K1+step*A1*indicate(s0,t)
        K2=K2+step*A1
     end if
  end do
!  print*,'K1imp=',K1,'K2imp=',K2
  call inv(p,K2,K21)
  K=matmul(K1,K21)
  do i=1,n
     if (delta(i).gt.0.5) then
        K1=0
        !fix = .true.
        call getAinv(p,n,m,h,x(i),r,q,x,delta,gt,qvs,A1,uf)
        if (minval(uf).gt.tol) then
           K1=A1*indicate(x(i),t)-matmul(K,A1)
           j=1
30         if (sum(abs(q(:,i)-qvs(:,j))).le.tol) then
              psi(:,i)=matmul(K1,q(:,i))/uf(j)+sum(K,2)
           else
              j=j+1
              goto 30
           end if
           n1=n1+1
        end if
     else
        n1=n1+1
     end if
  end do
  if (n1.ne.n) then
     print*,'in geteffestimput, failed to obtain n=n1, n1=',n1
    ! write(22,*) 'in geteffestimput, failed to obtain n=n1, n1=',n1
     eflag=-1
     goto 10
  end if
  call gethx(p,n,m,delta,psi,gt,r,id,Bp1)
  Fa=0
  do i=1,n
     if (delta(i).gt.0.5) then
        Fa=Fa+psi(:,i)
     else
        Fa=Fa+Bp1(:,i)
     end if
  end do
  Fa=Fa/n
10  return
end subroutine geteffestimputFa



 
 
 
! calculate estimation variance of ipw estimator 
subroutine getvar(p,n,phi,gt,st,delta,va) 
  integer,intent(in) :: p,n 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: gt,st,delta 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: Bphi 
  double precision,dimension(p,p,n) :: Bphi2 
  double precision,dimension(p,p) :: a 
  call Bp(p,n,phi,gt,st,delta,Bphi,Bphi2) 
  va=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        call mul(p,phi(:,i),phi(:,i),a) 
        va=va+a/gt(i) 
     else 
        call mul(p,Bphi(:,i),Bphi(:,i),a) 
        va=va+(Bphi2(:,:,i)-a)/(gt(i)**2.) 
     end if 
  end do 
  va=va/n/n 
  return 
end subroutine getvar 
 
! calculate estimation variance of augmented estimator 
subroutine getvaraug(p,n,m,phi,gt,st,x,delta,id,r,va) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: gt,st,x,delta 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: f1 
  double precision,dimension(p,p,n) :: f2 
  double precision,dimension(p,p) :: a 
  integer :: i 
  call getBfhs(p,n,m,x,delta,phi,gt,st,r,id,f1,f2) 
  va=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        call mul(p,phi(:,i),phi(:,i),a) 
        va=va+a/gt(i) 
     else 
        call mul(p,f1(:,i),f1(:,i),a) 
        va=va+(f2(:,:,i)-a)/(gt(i)**2.) 
     end if 
  end do 
  va=va/n/n 
  return 
end subroutine getvaraug 
 
! the h function values h(q_i,x_i) 
! using h2 is better than using sts. 
subroutine gethx(p,n,m,delta,phi,gt,r,id,hs) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: delta,gt 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,n),intent(out) :: hs 
  double precision,dimension(p,n) :: h1 
  double precision,dimension(n) :: h2 
  double precision,dimension(p) :: ft,fu,f1,f2 
  double precision :: b1,b2 
  integer:: i,j,k1,k2,k3 
  hs=0 
  h2=0 
  do k1=1,m 
     do k2=1,r(k1) 
        j=id(k1,k2) 
        if (delta(j).lt.0.5) then 
           do k3=k2+1,r(k1) 
              i=id(k1,k3) 
              if (delta(i).gt.0.5) then 
                 hs(:,j)=hs(:,j)+phi(:,i)/gt(i) 
                 h2(j)=h2(j)+1./gt(i) 
              end if 
           end do 
! or we can calculate h2(j)=r(k1)-sum(delta_i I(k3<=k2)) 
           if (h2(j).gt.0) then 
              hs(:,j)=hs(:,j)/h2(j) 
           end if 
        end if 
     end do 
  end do   
  return 
end subroutine gethx 
 
! B2 function for arbitrary h, function h is always a function of u 
subroutine B2hfu(p,n,m,x,phi,delta,gt,qq,r,id,bh,bf) 
 ! implicit none
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  double precision,dimension(p,p,n),intent(in) :: qq 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,n),intent(out) :: bh 
  double precision,dimension(p,p,n),intent(out) :: bf 
  double precision,dimension(p,n) :: hs,ans1,ans2 
  integer:: j 
  bh=0 
  bf=0 
  do j=1,n 
     if (delta(j).lt.0.5) then 
        call geths(p,n,m,j,x(j),x,delta,phi,gt,r,id,hs) 
        bh(:,j)=sum(hs(:,j:n),2) 
        bf(:,:,j)=sum(qq(:,:,j:n),3) 
        bh(:,j)=bh(:,j)/(n+1-j) 
        bf(:,:,j)=bf(:,:,j)/(n+1-j) 
     end if 
  end do 
  return 
end subroutine B2hfu 
 
 
! the h function values h(q_i,u) 
subroutine geths(p,n,m,i0,u,x,delta,phi,gt,r,id,hs) 
  integer,intent(in) :: p,n,m,i0 
  double precision,intent(in) :: u 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,n),intent(out) :: hs 
  double precision,dimension(p,n) :: h1 
  double precision,dimension(n) :: h2 
  double precision,dimension(p) :: ft,fu,f1,f2 
  double precision :: b1,b2 
  integer :: i,j,k1,k2,k3 
  hs=0 
  h2=0 
  do k1=1,m 
     do k2=1,r(k1) 
        j=id(k1,k2) 
        if (j.ge.i0) then 
           do k3=1,r(k1) 
              i=id(k1,k3) 
              if ((delta(i).gt.0.5).and.(x(i).ge.u)) then 
                 hs(:,j)=hs(:,j)+phi(:,i)/gt(i) 
                 h2(j)=h2(j)+1./gt(i) 
              end if 
           end do 
           if (h2(j).gt.0) then 
              hs(:,j)=hs(:,j)/h2(j) 
           end if 
        end if 
     end do 
  end do 
  return 
end subroutine geths 
 
 
!subroutine to calculate B(phi-h,x_j) and B((phi-h)^2,x_j)'s 
subroutine getBfhs(p,n,m,x,delta,phi,gt,st,r,id,f1,f2) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,dimension(p,n),intent(in) :: phi 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,n),intent(out) :: f1 
  double precision,dimension(p,p,n),intent(out) :: f2 
  double precision,dimension(p,n) :: hs 
  double precision,dimension(p,p) :: a 
  integer :: i,j 
  f1=0 
  f2=0 
  do j=1,n 
     if (delta(j).lt.0.5) then 
        call geths(p,n,m,j,x(j),x,delta,phi,gt,r,id,hs) 
        hs=phi-hs 
        do i=j,n 
           if (delta(i).gt.0.5) then 
              f1(:,j)=f1(:,j)+hs(:,i)/gt(i) 
              call mul(p,hs(:,i),hs(:,i),a) 
              f2(:,:,j)=f2(:,:,j)+a/gt(i) 
           end if 
        end do 
        f1(:,j)=f1(:,j)/n/st(j) 
        f2(:,:,j)=f2(:,:,j)/n/st(j) 
     end if 
  end do 
  return 
end subroutine getBfhs 
 
! B1 function on function h and hh', h does not depend on u 
subroutine Bp(p,n,h,gt,st,delta,bh1,bh2) 
  integer,intent(in) :: p,n 
  double precision,dimension(p,n),intent(in) :: h 
  double precision,dimension(n),intent(in) :: gt,st,delta 
  double precision,dimension(p,n),intent(out) :: bh1 
  double precision,dimension(p,p,n),intent(out) :: bh2 
  double precision,dimension(p,p) :: a 
  integer:: i,j 
  bh1=0 
  bh2=0 
  do j=1,n 
     if (delta(j).lt.0.5) then 
        do i=j,n 
           if (delta(i).gt.0.5) then 
              bh1(:,j)=bh1(:,j)+h(:,i)/gt(i) 
              call mul(p,h(:,i),h(:,i),a) 
              bh2(:,:,j)=bh2(:,:,j)+a/gt(i) 
           end if 
        end do 
        bh1(:,j)=bh1(:,j)/n/st(j) 
        bh2(:,:,j)=bh2(:,:,j)/n/st(j) 
     end if 
  end do 
  return 
end subroutine Bp 
 
 
!type I NPMLE 
subroutine npmle1(p,n,m,t,qvs,xsub,dsub,gts,sts,r,Fa,va) 
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(m,n),intent(in) :: xsub,dsub 
  double precision,dimension(m,n),intent(in) :: sts,gts 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,p) :: a,a1 
  double precision,dimension(p,m) :: a2 
  double precision,dimension(m) :: stt 
  double precision,dimension(n) :: phi,Bphi,Bphi2 
  double precision,dimension(m,m) :: vs 
  double precision :: tol,indicate 
  integer :: i,j,cnt,d 
  tol=1e-6 
  stt=1 
  do j=1,m 
     do i=1,r(j) 
        if (xsub(j,i).le.t) then 
           stt(j)=sts(j,i) 
        end if 
     end do 
  end do 
  a=matmul(qvs,transpose(qvs)) 
  call inv(p,a,a1) 
  a2=matmul(a1,qvs) 
  stt=1-stt 
  Fa=matmul(a2,stt) 
!now calculate variance 
  d=1 
  vs=0 
  do j=1,m 
     do i=1,r(j) 
        phi(i)=indicate(xsub(j,i),t)-stt(j) 
     end do 
     call Bp(d,r(j),phi,gts(j,:),sts(j,:),dsub(j,:),Bphi,Bphi2)  
     do i=1,r(j) 
        if (dsub(j,i).gt.0.5) then 
           vs(j,j)=vs(j,j)+phi(i)**2/gts(j,i) 
        else 
           vs(j,j)=vs(j,j)+(Bphi2(i)-Bphi(i)**2.)/(gts(j,i)**2.) 
        end if 
     end do 
     vs(j,j)=vs(j,j)/r(j)/r(j) 
  end do 
  va=matmul(matmul(a2,vs),transpose(a2)) 
  return  
end subroutine npmle1 


! the type II npmle estimator 
subroutine npmle2(p,n,q,x,delta,t,Fa) 
  integer,intent(in) :: p,n 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta 
  double precision,intent(in) :: t 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,n) :: fpdf,F 
  double precision,dimension(p,n) :: S,c,Sold 
  integer :: j,i 
  double precision :: tol,tmp,indicate 
  tol=1e-6 
  fpdf=1.0/n 
  do i=1,n 
     S(:,i)=1.*(1.+n-i)/n 
  end do 
  Sold=0 
10 if (maxval(abs(Sold-S)/n).gt.tol) then 
     Sold=S 
     do i=1,n 
        if (delta(i).gt.0.5) then 
           call inprod(p,q(:,i),fpdf(:,i),tmp) 
           do j=1,p 
              c(j,i)=q(j,i)*fpdf(j,i)/tmp 
           end do 
        else 
          call inprod(p,q(:,i),S(:,i),tmp) 
           do j=1,p 
              c(j,i)=q(j,i)*S(j,i)/tmp 
           end do 
        end if 
     end do 
     S=1 
     if (delta(1).gt.0.5) then 
        do j=1,p 
           tmp=sum(c(j,:)) 
           if (c(j,1).gt.tol) then  
              S(j,1)=(1.-c(j,1)/tmp) 
           end if 
        end do 
     end if 
     do i=2,n 
        if (delta(i).gt.0.5) then 
           do j=1,p 
              tmp=sum(c(j,i:n)) 
              if (c(j,i).gt.tol) then 
                 S(j,i)=S(j,i-1)*(1.-c(j,i)/tmp) 
              else 
                 S(j,i)=S(j,i-1) 
              end if 
           end do 
        else 
           S(:,i)=S(:,i-1) 
        end if 
        fpdf(:,i)=S(:,i-1)-S(:,i) 
     end do 
     fpdf(:,1)=1-S(:,1) 
     goto 10 
  end if 
  F=1-S 
  do i=1,n 
     if (x(i).le.t) then 
        Fa=F(:,i) 
     end if 
  end do 
  return 
end subroutine npmle2 
 
! invert a matrix 
subroutine inv(lb,A,A1) 
  integer,intent(in) :: lb 
  double precision,dimension(lb,lb),intent(in) :: A  
  double precision,dimension(lb,lb),intent(out) :: A1  
  integer :: i 
  integer,dimension(lb) :: ipvt 
  double precision,dimension(lb,lb) :: B 
  double precision :: cond 
  double precision,dimension(lb) :: wv 
  B=A 
  A1=0 
  do i=1,lb 
     A1(i,i)=1 
  end do 
  call decomp(lb,lb,B,cond,ipvt,wv) 
  do i=1,lb 
     call solvels(lb,lb,B,A1(:,i),ipvt) 
  end do 
  return 
end subroutine inv 
 
! multiply vectors: y=st' 
subroutine mul(lb,s,t,y) 
  integer,intent(in) :: lb 
  double precision,dimension(lb),intent(in) :: s,t 
  double precision,dimension(lb,lb),intent(out) :: y 
  integer :: i,j 
  do i=1,lb 
     do j=1,lb 
        y(i,j)=s(i)*t(j) 
     end do 
  end do 
  return 
end subroutine mul 
 
! inner product of two vectors: y=st' 
subroutine inprod(l,s,t,y) 
  integer,intent(in) :: l 
  double precision,dimension(l),intent(in) :: s,t 
  double precision,intent(out) :: y 
  integer :: i 
  y=0 
  do i=1,l 
     y=y+s(i)*t(i) 
  end do 
  return 
end subroutine inprod 
 
!the indicate function I(s<=t) 
function indicate(s,t) 
  double precision,intent(in) :: s,t 
  double precision :: indicate 
  double precision:: tmp 
  indicate=0 
  if (s.le.t) then 
     indicate=1 
  end if 
end function indicate 
 
 
! sample covariance matrix 
subroutine sampcov(lb,n,x,covmat) 
  implicit none 
  integer, intent (in) :: lb,n 
  double precision,dimension(lb,n),intent(in) :: x 
  double precision,dimension(lb,lb),intent(out) :: covmat 
  double precision :: mean1,mean2 
  integer :: i,j 
  covmat=0. 
  do i=1,lb 
     do j=1,lb 
        mean1=sum(x(i,:))/n 
        mean2=sum(x(j,:))/n 
        covmat(i,j)=sum((x(i,:)-mean1) * (x(j,:)-mean2))/(n-1.) 
       ! print*,'i=',i,'j=',j,'sum(x(i,:))/n=',sum(x(i,:))/n        
     end do 
  end do 
  return 
end subroutine sampcov 
 
 
 
 
 
 
 
 
 
!!!!!!!!!!!!!!!The remaining code is not used, but save for future use 
! calculate estimation variance of augmented estimator, expand everything, slow 
subroutine getvaraug1(p,n,m,phi,gt,st,x,delta,id,r,va) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,n),intent(in) :: phi 
  double precision,dimension(n),intent(in) :: gt,st,x,delta 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: h1,f1 
  double precision,dimension(p,p,n) :: Bphi2,h2,h3,h4,h5,f2 
  double precision,dimension(p,p) :: a,a2,a3,b1,b2,b3 
  double precision,dimension(p,m) :: Fam 
  double precision,dimension(p) :: Fam0 
  integer :: i 
  call getBfs(p,n,phi,gt,st,x,delta,f1,f2) 
  call getBhs(p,n,m,x,delta,phi,f1,gt,r,id,h1,h2,h3) 
  call getBphih(p,n,m,x,delta,gt,st,phi,r,id,h4) 
  call getphiBh(p,n,delta,gt,st,phi,h1,h5) 
  va=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        call mul(p,phi(:,i),phi(:,i),a) 
        va=va+a/gt(i) 
     else 
        call mul(p,f1(:,i),f1(:,i),a) 
        call mul(p,h1(:,i),h1(:,i),a2) 
        call mul(p,f1(:,i),h1(:,i),a3) 
        b1=f2(:,:,i)-a+h2(:,:,i)-a2-a3-transpose(a3) 
        b2=-h4(:,:,i)-transpose(h4(:,:,i))+h5(:,:,i)+transpose(h5(:,:,i)) 
        b3=h3(:,:,i)+transpose(h3(:,:,i)) 
        va=va+(b1+b2+b3)/(gt(i)**2.) 
     end if 
  end do 
  va=va/n/n 
  return 
end subroutine getvaraug1 
 
subroutine getBfs(p,n,h,gt,st,x,delta,f1,f2) 
  integer,intent(in) :: p,n 
  double precision,dimension(p,n),intent(in) :: h 
  double precision,dimension(n),intent(in) :: gt,st,x,delta 
  double precision,dimension(p,n),intent(out) :: f1 
  double precision,dimension(p,p,n),intent(out) :: f2 
  double precision,dimension(p,p) :: a 
  integer:: i,j 
  f1=0 
  f2=0 
  do j=1,n 
     do i=j,n 
        if (delta(i).gt.0.5) then 
           f1(:,j)=f1(:,j)+h(:,i)/gt(i) 
           call mul(p,h(:,i),h(:,i),a) 
           f2(:,:,j)=f2(:,:,j)+a/gt(i) 
        end if 
     end do 
     f1(:,j)=f1(:,j)/n/st(j) 
     f2(:,:,j)=f2(:,:,j)/n/st(j) 
  end do 
  return 
end subroutine getBfs 
 
subroutine getBhs(p,n,m,x,delta,phi,f1,gt,r,id,h1,h2,h3) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  double precision,dimension(p,n),intent(in) :: f1,phi 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,n),intent(out) :: h1 
  double precision,dimension(p,p,n),intent(out) :: h2,h3 
  double precision,dimension(p,n) :: hs 
  double precision,dimension(p,p) :: a 
  integer :: i,j 
  h1=0 
  h2=0 
  h3=0 
  do j=1,n 
     call geths(p,n,m,j,x(j),x,delta,phi,gt,r,id,hs) 
     do i=j,n 
        call mul(p,hs(:,i),hs(:,i),a) 
        h2(:,:,j)=h2(:,:,j)+a 
        call mul(p,hs(:,i),f1(:,i),a) 
        h3(:,:,j)=h3(:,:,j)+a 
     end do 
     h1(:,j)=sum(hs(:,j:n),2) 
     h1(:,j)=h1(:,j)/(n+1-j) 
     h2(:,:,j)=h2(:,:,j)/(n+1-j) 
     h3(:,:,j)=h3(:,:,j)/(n+1-j) 
  end do 
  return 
end subroutine getBhs 
 
subroutine getBphih(p,n,m,x,delta,gt,st,phi,r,id,h) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,dimension(p,n),intent(in) :: phi 
  integer,dimension(m),intent(in) :: r 
  integer,dimension(m,n),intent(in) :: id 
  double precision,dimension(p,p,n),intent(out) :: h 
  double precision,dimension(p,n) :: hs 
  double precision,dimension(p,p) :: a 
  integer :: i,j 
  h=0 
  do j=1,n 
     call geths(p,n,m,j,x(j),x,delta,phi,gt,r,id,hs) 
     do i=j,n 
        if (delta(i).gt.0.5) then 
           call mul(p,phi(:,i),hs(:,i),a) 
           h(:,:,j)=h(:,:,j)+a/gt(i) 
        end if 
     end do 
     h(:,:,j)=h(:,:,j)/n/st(j) 
  end do 
  return 
end subroutine getBphih 
 
subroutine getphiBh(p,n,delta,gt,st,phi,h1,h) 
  integer,intent(in) :: p,n 
  double precision,dimension(n),intent(in) :: delta,gt,st 
  double precision,dimension(p,n),intent(in) :: phi,h1 
  double precision,dimension(p,p,n),intent(out) :: h 
  double precision,dimension(p,n) :: hs 
  double precision,dimension(p,p) :: a 
  integer :: i,j 
  hs=0 
  do j=1,n 
     do i=j,n 
        if (delta(i).gt.0.5) then 
           hs(:,j)=hs(:,j)+phi(:,i)/gt(i) 
        end if 
     end do 
     hs(:,j)=hs(:,j)/n/st(j) 
     call mul(p,hs(:,j),h1(:,j),h(:,:,j)) 
  end do 
  return 
end subroutine getphiBh 
 
 
! the efficient estimator 
subroutine effest(p,n,m,t,qvs,q,x,delta,gt,st,r,K,Fa,va,eflag) 
  integer, intent(in):: p,n,m 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,intent(in) :: t 
  double precision,dimension(p,m),intent(in) :: qvs 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va,K 
  integer, intent(out) :: eflag
  double precision,dimension(p,p) :: A1,K1,K2,K21 
  double precision,dimension(m) :: uf  
  double precision,dimension(p,n) :: psi,phi 
  double precision :: h,s0,step,tol,indicate,a 
  integer :: i,j,iter,pts,n1 
  eflag=1
  n1=0 
  tol=1e-6 
  Fa=0 
  a=0
  h=(maxval(x)-minval(abs(x)))*(1.*n/m)**(-0.2)*2.5 
  K1=0 
  K2=0 
  pts=10000 
  step=(maxval(x)-minval(abs(x))+2.*h*n/m)/pts 
  do iter=0,pts 
     s0=minval(abs(x))-h*n/m+(iter+0.5)*step 
     call getAinv(p,n,m,h,s0,r,q,x,delta,gt,qvs,A1,uf) 
     if (minval(uf).gt.tol) then  
        K1=K1+step*A1*indicate(s0,t) 
        K2=K2+step*A1 
     end if 
  end do 
  call inv(p,K2,K21) 
  K=matmul(K1,K21) 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        K1=0 
        call getAinv(p,n,m,h,x(i),r,q,x,delta,gt,qvs,A1,uf) 
        if (minval(uf).gt.tol) then                
           K1=A1*indicate(x(i),t)-matmul(K,A1) 
           j=1 
30         if (sum(abs(q(:,i)-qvs(:,j))).le.tol) then 
              psi(:,i)=matmul(K1,q(:,i))/uf(j)+sum(K,2) 
           else 
              j=j+1 
              goto 30 
           end if 
           Fa=Fa+psi(:,i)/gt(i)
           a=a+1/gt(i)
           n1=n1+1 
        end if 
     else 
        n1=n1+1 
     end if 
  end do 
  if (n1.ne.n) then 
     print*,'in effest, failed to obtain n=n1, n1=',n1 
    ! write(22,*) 'in effest, failed to obtain n=n1, n1=',n1 
     eflag=-1
     goto 10
  end if 
  Fa=Fa/a
  do i=1,n
     phi(:,i)=psi(:,i)-Fa
  end do
  call getvar(p,n,phi,gt,st,delta,va) 
10  return 
end subroutine effest 
 


! the efficient estimator using true ft with augmentation
subroutine effestaug(p,n,m,t,qvs,q,x,delta,gt,st,r,id,K,Fa,va)
  integer, intent(in):: p,n,m
  double precision,dimension(p,n),intent(in) :: q
  double precision,dimension(n),intent(in) :: x,delta,gt,st
  double precision,intent(in) :: t
  double precision,dimension(p,m),intent(in) :: qvs
  integer,dimension(m),intent(in) :: r
  integer,dimension(m,n) :: id
  double precision,dimension(p,p),intent(in) :: K
  double precision,dimension(p),intent(out) :: Fa
  double precision,dimension(p,p),intent(out) :: va
  double precision,dimension(p,p) :: A1,K1,K2,K21
  double precision,dimension(m) :: uf
  double precision,dimension(p,n) :: psi,phi,h,Bh
  double precision,dimension(p,p,n) :: Bf,qq
  double precision :: s0,step,tol,indicate,a,h0
  integer :: i,j,iter,pts,n1
  do i=1,n
     call mul(p,q(:,i),q(:,i),qq(:,:,i))
  end do
  n1=0
  tol=1e-6
  h0=(maxval(x)-minval(abs(x)))*(1.*n/m)**(-0.2)*2.5
  do i=1,n
     if (delta(i).gt.0.5) then
        K1=0
        call getAinv(p,n,m,h0,x(i),r,q,x,delta,gt,qvs,A1,uf)
        if (minval(uf).gt.tol) then
           K1=A1*indicate(x(i),t)-matmul(K,A1)
           j=1
30         if (sum(abs(q(:,i)-qvs(:,j))).le.tol) then
              psi(:,i)=matmul(K1,q(:,i))/uf(j)+sum(K,2)
           else
              j=j+1
              goto 30
           end if
        end if
     end if
  end do
  call gethx(p,n,m,delta,psi,gt,r,id,h)
  call B2hfu(p,n,m,x,psi,delta,gt,qq,r,id,Bh,Bf)
  Bh=h-Bh
  Fa=0
  a=0
  do i=1,n
     if (delta(i).gt.0.5) then
        Fa=Fa+psi(:,i)/gt(i)
        a=a+1/gt(i)
     else
        Fa=Fa+Bh(:,i)/gt(i)
     end if
  end do
  Fa=Fa/a
  do i=1,n
     if (delta(i).gt.0.5) then
        phi(:,i)=psi(:,i)-Fa
     end if
  end do
  call getvaraug(p,n,m,phi,gt,st,x,delta,id,r,va)
  return
end subroutine effestaug



 

!subroutine to calculate A^{-1}(s,q^T\hat f) and  \hat f(s) 
subroutine getAinv(p,n,m,h,s0,r,q,x,delta,gt,qvs,A1,uf) 
  integer, intent(in):: p,n,m 
  double precision,dimension(n),intent(in) :: x,delta,gt 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(p,n),intent(in) :: q 
  integer,dimension(m),intent(in) :: r 
  double precision,intent(in) :: h,s0 
  double precision,dimension(p,p),intent(out) :: A1 
  double precision,dimension(m),intent(out) :: uf  
  double precision,dimension(p,p) :: A,At 
  double precision :: tmp,quaker,step,tol 
  integer :: i,j  
  tol=1e-6 
!calculate uf 
  uf=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        j=1 
30      if (sum(abs(q(:,i)-qvs(:,j))).le.tol) then 
           uf(j)=uf(j)+quaker(h,x(i),s0)/gt(i) 
           !uf(j)=uf(j)+quaker(h,x(i),s0) 
        else 
           j=j+1 
           goto 30 
        end if 
     end if 
  end do 
  do j=1,m 
     uf(j)=uf(j)/r(j) 
  end do 
!calculate A matrix 
  A=0 
  do j=1,m 
     call mul(p,qvs(:,j),qvs(:,j),At) 
     A=A+At/uf(j)*r(j)/n 
  end do 
  call inv(p,A,A1) 
  return 
end subroutine getAinv 
 
 
 
    
!quatic kernel 
function quaker(h,xi,x) 
  double precision,intent(in) :: h,xi,x 
  double precision :: quaker 
  double precision:: tmp 
  tmp=((xi-x)/h)**2. 
  if (tmp.ge.1) then 
     quaker=0   
  else 
     quaker=(1.-tmp)/h 
  end if 
end function quaker 
 
!true G(t) 
subroutine trueG(t,g) 
  double precision, intent(in) :: t 
  double precision, intent(out) :: g 
  g=1.-t/10. 
  return 
end subroutine trueG 
 
!true G(t) and mixed infinity 
subroutine truemixG(n,x,g) 
  integer,intent(in) :: n 
  double precision,dimension(n),intent(in) :: x 
  double precision,dimension(n),intent(out) :: g 
  g=1.-x*.09 
  return 
end subroutine truemixG 
 
!the imputation method on ols, using asymptotic variance estimation, code not right. 
subroutine olsimputasymp(p,n,m,t,qvs,q,x,delta,gt,st,xsub,dsub,gts,sts,id,rt,r,qq,a1,&
     lim,rat,setting,d,H0,Fa,va)
  implicit none
  character(len=10),intent(in) :: setting
  logical,intent(in) :: H0
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: t,d 
  double precision,intent(in) :: lim,rat
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(p,m),intent(in) :: qvs  
  double precision,dimension(n),intent(in) :: x,delta,gt,st 
  double precision,dimension(m,n),intent(in) :: xsub,dsub,gts,sts 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: rt,r 
  double precision,dimension(p,p,n),intent(in) :: qq 
  double precision,dimension(p,p),intent(in) :: a1 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: psi,Bp1,hs,phi 
  double precision,dimension(p,p) :: a,a2 
  double precision :: indicate,tol,tmp 
  integer :: i 
  tol=1e-6 
  psi=0 
  Fa=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t) 
     end if 
  end do 
!the true heff 
!  call olsheffx(p,n,q,x,t,qq,Bp1,lim,rat,setting,d,H0) 
!the estimated heff 
  call gethx(p,n,m,delta,psi,gt,r,id,Bp1) 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i) 
     else 
        Fa=Fa+Bp1(:,i) 
     end if 
  end do 
  a=sum(qq,3) 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     if (delta(i).gt.0.5) then 
        phi(:,i)=psi(:,i)-q(:,i)*tmp 
     else 
        phi(:,i)=Bp1(:,i)-matmul(qq(:,:,i),Fa) 
     end if 
  end do 
  phi=matmul(a1,phi) 
! asymptotic variance, not correct.   
 call getvarimput(p,n,m,qvs,phi,gt,st,q,x,delta,r,Fa,va) 
  return 
end subroutine olsimputasymp 
 
!the imputation method on wls, using asymptotic variance estimation, code not right. 
subroutine wlsimputasymp(p,n,m,t,qvs,q,x,delta,gt,st,xsub,dsub,gts,sts,id,rt,r,w,qq,a1,&
     lim,rat,setting,d,H0,Fa,va)
  implicit none
  character(len=10),intent(in) :: setting
  logical,intent(in) :: H0
  integer,intent(in) :: p,n,m 
  double precision,intent(in) :: lim,rat,d
  double precision,intent(in) :: t 
  double precision,dimension(p,n),intent(in) :: q 
  double precision,dimension(p,m),intent(in) :: qvs  
  double precision,dimension(n),intent(in) :: x,delta,gt,st,w 
  double precision,dimension(m,n),intent(in) :: xsub,dsub,gts,sts 
  integer,dimension(m,n),intent(in) :: id 
  integer,dimension(m),intent(in) :: rt,r 
  double precision,dimension(p,p,n),intent(in) :: qq 
  double precision,dimension(p,p),intent(in) :: a1 
  double precision,dimension(p),intent(out) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,n) :: psi,Bp1,hs,phi 
  double precision,dimension(p,p) :: a,a2 
  double precision :: indicate,tol,tmp 
  integer :: i 
  tol=1e-6 
  psi=0 
  Fa=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        psi(:,i)=q(:,i)*indicate(x(i),t)*w(i) 
     end if 
  end do 
!the true heff 
!  call wlsheffx(p,n,q,x,w,t,qq,Bp1,lim,rat,setting,d,H0) 
!the estimated heff 
  call gethx(p,n,m,delta,psi,gt,r,id,Bp1) 
  Fa=0 
  a=0 
  do i=1,n 
     if (delta(i).gt.0.5) then 
        Fa=Fa+psi(:,i) 
     else 
        Fa=Fa+Bp1(:,i) 
     end if 
  end do 
  a=sum(qq,3) 
  call inv(p,a,a2) 
  Fa=matmul(a2,Fa) 
  do i=1,n 
     call inprod(p,q(:,i),Fa,tmp) 
     if (delta(i).gt.0.5) then 
        phi(:,i)=psi(:,i)-q(:,i)*tmp*w(i) 
     else 
        phi(:,i)=Bp1(:,i)-matmul(qq(:,:,i),Fa) 
     end if 
  end do 
  phi=matmul(a1,phi) 
 call getvarimput(p,n,m,qvs,phi,gt,st,q,x,delta,r,Fa,va) 
  return 
end subroutine wlsimputasymp 
 
 
! calculate estimation variance of imputation estimator 
subroutine getvarimput(p,n,m,qvs,phi,gt,st,q,x,delta,r,Fa,va) 
  integer,intent(in) :: p,n,m 
  double precision,dimension(p,m),intent(in) :: qvs 
  double precision,dimension(p,n),intent(in) :: phi,q 
  double precision,dimension(n),intent(in) :: gt,st,x,delta 
  integer,dimension(m),intent(in) :: r 
  double precision,dimension(p),intent(in) :: Fa 
  double precision,dimension(p,p),intent(out) :: va 
  double precision,dimension(p,p,n) :: Bphi2 
  double precision,dimension(p,p) :: a 
  double precision,dimension(p,m) :: Fam 
  integer :: i,j 
  va=0 
  do i=1,n 
     call mul(p,phi(:,i),phi(:,i),a) 
     va=va+a 
  end do 
  va=va/n/n 
  return 
end subroutine getvarimput 
 
