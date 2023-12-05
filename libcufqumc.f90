! MIT License
! 
! Copyright (c) 2023 Shinsuke Ideguchi
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mcmc_mod
  use cudafor
  implicit none
  integer         , parameter :: Nsrctype = 13
  integer         , parameter :: Npl(Nsrctype) = (/3,4,4,5,5,6,4,5,5,6,5,5,6/)
  double precision, parameter :: ls  = 2.99792458d8    ! light speed [m/s]
  double precision, parameter :: pi  = acos(-1.d0)
  double precision, parameter :: ipi = 1.d0/pi
  double precision, parameter :: isqrtpi = 1.d0/sqrt(pi)
contains

  attributes(global) subroutine calc_chisq_array(Ndata,Nchan,lamsq,q,u,q_err,u_err,mask_chan,src_list,Npar,Nbeta,p,chisq_array)
      implicit none
      integer          ,intent(in)    :: Ndata, Nchan, Npar, Nbeta, src_list(Nsrctype), mask_chan(Nchan,Ndata)
      double precision ,intent(in)    :: lamsq(Nchan), q(Nchan,Ndata), u(Nchan,Ndata), q_err(Nchan,Ndata), u_err(Nchan,Ndata)
      double precision ,intent(in)    :: p(Npar,Nbeta,Ndata)
      double precision ,intent(inout) :: chisq_array(Nchan,Nbeta,Ndata)
      integer                         :: i, ichan, ibeta, idata, k, isrc, s
      double precision                :: p1, p2, p3, p4, p5
      double precision                :: fact1, fact2
      complex(kind(0d0))              :: PI_sum
      double precision                :: tmpq, tmpu, sinc

      i = (blockIdx%x-1)*blockDim%x+threadIdx%x
      if(i<=Nchan*Nbeta*Ndata) then
        idata = int((i-1)/(Nchan*Nbeta))+1
        ibeta = mod(int((i-1)/Nchan),Nbeta)+1
        ichan = mod(i-1,Nchan)+1

        if (mask_chan(ichan,idata) == 0) then
          PI_sum = 0.d0
          s = 0

          ! 01: delta function
          k = 1
          do isrc = 1, src_list(k)
            p1 = p(Npl(k)*(isrc-1)+1+s,ibeta,idata)
            p2 = p(Npl(k)*(isrc-1)+2+s,ibeta,idata)
            p3 = p(Npl(k)*(isrc-1)+3+s,ibeta,idata)
            fact1 = p2
            fact2 = 2.d0*(p1*lamsq(ichan)+p3)
            PI_sum = PI_sum + fact1*cmplx(cos(fact2),sin(fact2))
          end do
          s = s + src_list(k)*Npl(k)

          ! 03: Gaussian function
          k = 3
          do isrc = 1, src_list(k)
            p1 = p(Npl(k)*(isrc-1)+1+s,ibeta,idata)
            p2 = p(Npl(k)*(isrc-1)+2+s,ibeta,idata)
            p3 = p(Npl(k)*(isrc-1)+3+s,ibeta,idata)
            p4 = p(Npl(k)*(isrc-1)+4+s,ibeta,idata)
            fact1 = p2*exp(-2.d0*p4**2*lamsq(ichan)**2)
            fact2 = 2.d0*(p1*lamsq(ichan)+p3)
            PI_sum = PI_sum + fact1*cmplx(cos(fact2),sin(fact2))
          end do
          s = s + src_list(k)*Npl(k)

          ! 07: tophat function
          k = 7
          do isrc = 1, src_list(k)
            p1 = p(Npl(k)*(isrc-1)+1+s,ibeta,idata)
            p2 = p(Npl(k)*(isrc-1)+2+s,ibeta,idata)
            p3 = p(Npl(k)*(isrc-1)+3+s,ibeta,idata)
            p4 = p(Npl(k)*(isrc-1)+4+s,ibeta,idata)
            fact1 = p2/(4.d0*p4*lamsq(ichan))*sin(2.d0*p4*lamsq(ichan))
            fact2 = 2.d0*(p1*lamsq(ichan)+p3)
            PI_sum = PI_sum + fact1*cmplx(cos(fact2),sin(fact2))
          end do
          s = s + src_list(k)*Npl(k)

          ! 11: osul17 model
          k = 11
          do isrc = 1, src_list(k)
            p1 = p(Npl(k)*(isrc-1)+1+s,ibeta,idata)
            p2 = p(Npl(k)*(isrc-1)+2+s,ibeta,idata)
            p3 = p(Npl(k)*(isrc-1)+3+s,ibeta,idata)
            p4 = p(Npl(k)*(isrc-1)+4+s,ibeta,idata)
            p5 = p(Npl(k)*(isrc-1)+5+s,ibeta,idata)
            fact1 = 2.d0*(p1*lamsq(ichan)+p3)
            fact2 = 2.d0*p5*lamsq(ichan)
            sinc = 1.d0
            if (fact2 /= 0.d0) then
              sinc = sin(fact2)/(2.d0*fact2)
            end if
            tmpq = cos(fact1)*isqrtpi
            tmpu = sin(fact1)*isqrtpi
            PI_sum = PI_sum + p2*exp(-2.d0*p4**2*lamsq(ichan)**2)*sinc*cmplx(tmpq,tmpu)
          end do
          s = s + src_list(k)*Npl(k)

          chisq_array(ichan,ibeta,idata) = (real(PI_sum)-q(ichan,idata))**2/q_err(ichan,idata)**2 + (aimag(PI_sum)-u(ichan,idata))**2/u_err(ichan,idata)**2
        else
          chisq_array(ichan,ibeta,idata) = 0.d0
        end if
      endif
  end subroutine calc_chisq_array

end module mcmc_mod


program main
  use cudafor
  use curand
  use mcmc_mod
  implicit none

  ! data
  integer           , managed     :: Ndata, Nchan
  integer           , allocatable :: mask_chan(:,:)
  double precision  , allocatable :: lamsq(:), q(:,:), u(:,:), q_err(:,:), u_err(:,:)
  logical           , allocatable :: mask_cand(:,:)
  character(len=128), allocatable :: data_list(:)

  ! beta
  integer           , managed     :: Nbeta
  double precision  , allocatable :: beta(:,:)

  ! model
  integer           , managed     :: Npar
  double precision  , allocatable :: par(:,:,:,:)
  integer           , managed     :: src_list(Nsrctype)
  integer           , allocatable :: peri(:)

  ! mcmc
  integer                         :: Nburnin, Nsample, Nadjust
  double precision                :: logu, logr, factor, sign
  double precision  , allocatable :: chain(:,:,:), chisq(:,:), cand(:,:,:), chisq_old(:,:)
  double precision  , allocatable :: par_chain(:,:,:), chisq_chain(:,:,:)

  ! cuda
  integer                         :: blockDim, gridDim, err, Nthreads
  double precision  , allocatable :: chisq_array(:,:,:)

  integer                         :: Nexchevry, Nexch
  double precision  , allocatable :: dum_chain(:)

  ! 
  integer                         :: i, n, idata, iiter, iiterB, ipar, ibeta, ichan, iexch, NiterB

  double precision                :: width_max, widthinc, adjust_width_fact, val
  double precision                :: ratio
  double precision                :: dum1, dum5(5)
  character(len=128)              :: dirname

  integer           , allocatable :: Nacpt(:,:,:), Nstep(:,:,:)

  ! random number
  integer                         :: Nrnd, seed, rnd_i
  type(curandGenerator)           :: gen
  double precision  , allocatable :: rnd(:)

  !!! getarg
  call getarg(1,dirname)

  !!! read files

  ! data list
  Ndata = 0
  open(10,file=trim(dirname)//'datalist.txt',status='old')
  do
      read(10,*,end=10)
      Ndata = Ndata + 1
  end do
10  rewind(10)
  allocate( data_list (Ndata) )
  do idata = 1, Ndata
      read(10,*) data_list(idata)
  end do
  close(10)

  ! freq
  Nchan = 0
  open(10,file=trim(dirname)//'freq.txt',status='old')
  do
      read(10,*,end=11)
      Nchan = Nchan + 1
  end do
11  rewind(10)
  allocate( lamsq (Nchan) )
  do ichan = 1, Nchan
      read(10,*) dum1
      lamsq(ichan) = ls**2/dble(dum1)**2
  end do
  close(10)

  ! data
  allocate(         q (Nchan,Ndata) )
  allocate(         u (Nchan,Ndata) )
  allocate(     q_err (Nchan,Ndata) )
  allocate(     u_err (Nchan,Ndata) )
  allocate( mask_chan (Nchan,Ndata) )
  do idata = 1, Ndata
    open(10,file=trim(dirname)//trim(data_list(idata))//'data.txt',status='old')
    do ichan = 1, Nchan
        read(10,*) dum5
        q(ichan,idata)         = dum5(1)
        u(ichan,idata)         = dum5(2)
        q_err(ichan,idata)     = dum5(3)
        u_err(ichan,idata)     = dum5(4)
        mask_chan(ichan,idata) = dum5(5)
    end do
    close(10)
  end do

  ! MCMC parameters
  open(10,file=trim(dirname)//'mcmc.txt',status='old')
  read(10,*) Nbeta
  read(10,*) Nburnin
  read(10,*) Nsample
  read(10,*) Nadjust
  read(10,*) Nexchevry
  read(10,*) Nexch
  read(10,*) Npar
  read(10,*) Nthreads
  read(10,*) seed
  close(10)

  ! beta
  allocate( beta (Nbeta,Ndata) )
  do idata = 1, Ndata
    open(10,file=trim(dirname)//trim(data_list(idata))//'beta.txt',status='old')
    do ibeta = 1, Nbeta
        read(10,*) beta(ibeta,idata)
    end do
    close(10)
  end do

  ! src_list
  open(10,file=trim(dirname)//'src_list.txt',status='old')
  do i = 1, Nsrctype
    read(10,*) src_list(i)
  end do
  close(10)

  ! peri
  allocate( peri (Npar) )
  open(10,file=trim(dirname)//'peri.txt',status='old')
  do ipar = 1, Npar
    read(10,*) peri(ipar)
  end do
  close(10)

  ! par
  allocate( par (4,Npar,Nbeta,Ndata) )
  do idata = 1, Ndata
    open(10,file=trim(dirname)//trim(data_list(idata))//'par.dat',form='unformatted')
    read(10) par(:,:,:,idata)
    close(10)
  end do

  ! allocate
  allocate( Nacpt (Npar,Nbeta,Ndata) )
  allocate( Nstep (Npar,Nbeta,Ndata) )
  allocate( mask_cand (Nbeta,Ndata) )
  allocate( cand (Npar,Nbeta,Ndata) )
  allocate( dum_chain (Npar+1) )
  allocate( chisq (Nbeta,Ndata) )
  allocate( chisq_old (Nbeta,Ndata) )
  allocate( chisq_array(Nchan,Nbeta,Ndata) )
  allocate( chain (Npar+1,Nbeta,Ndata) )

  ! gridDim, blockDim
  blockDim = Nthreads
  gridDim  = (Nchan*Nbeta*Ndata+blockDim-1)/blockDim

  ! preparation of random number
  err = curandCreateGeneratorHost(gen,CURAND_RNG_PSEUDO_DEFAULT)
  err = curandSetPseudoRandomGeneratorSeed(gen,seed)
  ! random number
  Nrnd = 2**30
  allocate( rnd (Nrnd) )
  rnd_i = 0
  err = curandGenerate(gen, rnd, Nrnd)



  ! start MCMC

  ! burnin
  if (Nburnin > 0) then

    val = 3.d-1
    widthinc = 1.2d0
    NiterB = int(Nburnin/Nadjust)
    ! allocate
    allocate(   par_chain (Npar,Nburnin,Ndata) )
    allocate( chisq_chain (Nbeta,Nburnin,Ndata) )
    ! calculate chain(1)
    chain(1:Npar,:,:) = par(1,:,:,:)
    call calc_chisq_array<<<gridDim,blockDim>>>(Ndata,Nchan,lamsq,q,u,q_err,u_err,mask_chan,src_list,Npar,Nbeta,chain(1:Npar,:,:),chisq_array)
    chain(Npar+1,:,:) = sum(chisq_array,dim=1)

    do iiterB = 1, NiterB
      ! Initialization
      Nacpt = 0
      Nstep = 0
      do iiter = 1, Nadjust

        n = int((iiterB-1)*Nadjust) + iiter
        do ipar = 1, Npar

          ! candidate
          cand(:,:,:)    = chain(1:Npar,:,:)
          chisq_old(:,:) = chain(Npar+1,:,:)

          ! choice of parameter to move
          i = rand_int_unif(1,Npar)

          ! work
          do idata = 1, Ndata

            ! move parameter p(Npar,Nbeta)
            do ibeta = 1, Nbeta
              cand(i,ibeta,idata) = rand_dble_norm(cand(i,ibeta,idata),par(2,i,ibeta,idata))
              if(peri(i).eq.1) then
                mask_cand(ibeta,idata) = .true.
                factor                 = int(abs(cand(i,ibeta,idata))*ipi+5.d-1)
                sign                   = abs(cand(i,ibeta,idata))/cand(i,ibeta,idata)
                cand(i,ibeta,idata)    = cand(i,ibeta,idata) - pi*sign*factor
              else
                if(cand(i,ibeta,idata)>par(3,i,ibeta,idata) .and. cand(i,ibeta,idata)<par(4,i,ibeta,idata)) then
                  mask_cand(ibeta,idata) = .true.
                else
                  mask_cand(ibeta,idata) = .false.
                end if
              end if
            end do
          end do
          ! call kernel: calculation of chisq array
          call calc_chisq_array<<<gridDim,blockDim>>>(Ndata,Nchan,lamsq,q,u,q_err,u_err,mask_chan,src_list,Npar,Nbeta,cand,chisq_array)
          chisq = sum(chisq_array,dim=1)

          ! adoption/rejection of the candidate parameter
          do idata = 1, Ndata
            do ibeta = 1, Nbeta
              if (mask_cand(ibeta,idata)) then
                logr = -5.d-1 * beta(ibeta,idata) * (chisq(ibeta,idata) - chisq_old(ibeta,idata))

                ! accept the candidate parameter with the probability equal to the Metropolis ratio
                logu = log(rand_num())
                Nstep(i,ibeta,idata) = Nstep(i,ibeta,idata) + 1
                if ( (logu < logr) .or. (logr >= 0) ) then
                  Nacpt(i,ibeta,idata) = Nacpt(i,ibeta,idata) + 1
                  chain(1:Npar,ibeta,idata) = cand(:,ibeta,idata)
                  chain(Npar+1,ibeta,idata) = chisq(ibeta,idata)
                end if
              end if
            end do
          end do
        end do

        ! exchange
        if(Nbeta>1 .and. n>1 .and. mod(n+1,Nexchevry)==0) then

          do idata = 1, Ndata
            do iexch = 1, Nbeta * Nexch

              i = rand_int_unif(1,Nbeta-1)

              logr = 5.d-1 * (beta(i+1,idata) - beta(i,idata)) &
              * (chain(Npar+1,i+1,idata) - chain(Npar+1,i,idata))

              logu = log(rand_num())

              if ( (logu < logr) .or. (logr >= 0) ) then
                dum_chain(:)       = chain(:,i+1,idata)
                chain(:,i+1,idata) = chain(:,i,idata)
                chain(:,i,idata)   = dum_chain(:)
              end if
            end do
          end do
        end if

        par_chain(:,n,:) = chain(1:Npar,Nbeta,:)
        chisq_chain(:,n,:) = chain(Npar+1,:,:)
      end do

      ! adjust step width
      do idata = 1, Ndata
        do ibeta = 1, Nbeta
          do ipar = 1, Npar
            ratio = dble(Nacpt(ipar,ibeta,idata))/dble(Nstep(ipar,ibeta,idata)) ! calculate acceptance ratio
            width_max = (par(4,ipar,ibeta,idata)-par(3,ipar,ibeta,idata))*2.d-1
            adjust_width_fact = abs(ratio-val) + widthinc
            if (ratio>val .and. par(2,ipar,ibeta,idata)<width_max) then
              par(2,ipar,ibeta,idata) = par(2,ipar,ibeta,idata)*adjust_width_fact
            else
              par(2,ipar,ibeta,idata) = par(2,ipar,ibeta,idata)/adjust_width_fact
            end if
          end do
        end do
      end do

    end do

    !!! output
    ! do idata = 1, Ndata
    !   open(10,file=trim(dirname)//trim(data_list(idata))//'chainB.dat',form='unformatted')
    !   write(10) par_chain(:,:,idata)
    !   close(10)
    ! end do

    ! do idata = 1, Ndata
    !   open(10,file=trim(dirname)//trim(data_list(idata))//'chisqB.dat',form='unformatted')
    !   write(10) chisq_chain(:,:,idata)
    !   close(10)
    ! end do

    !!! initial values for sampling step
    par(1,:,:,:) = chain(1:Npar,:,:)
    do idata = 1, Ndata
      open(10,file=trim(dirname)//trim(data_list(idata))//'par.dat',form='unformatted')
      write(10) par(:,:,:,idata)
      close(10)
    end do

    ! deallocate
    deallocate(par_chain,chisq_chain)
  end if


  ! sampling
  if (Nsample > 0) then

    ! allocate
    allocate(   par_chain (Npar,Nsample,Ndata) )
    allocate( chisq_chain (Nbeta,Nsample,Ndata) )
    ! calculate chain(1)
    chain(1:Npar,:,:) = par(1,:,:,:)
    call calc_chisq_array<<<gridDim,blockDim>>>(Ndata,Nchan,lamsq,q,u,q_err,u_err,mask_chan,src_list,Npar,Nbeta,chain(1:Npar,:,:),chisq_array)
    chain(Npar+1,:,:) = sum(chisq_array,dim=1)

    do iiter = 1, Nsample

      do ipar = 1, Npar

        ! candidate
        cand(:,:,:) = chain(1:Npar,:,:)
        chisq_old(:,:) = chain(Npar+1,:,:)

        ! choice of parameter to move
        i = rand_int_unif(1,Npar)

        ! work
        do idata = 1, Ndata

          ! move parameter p(Npar,Nbeta)
          do ibeta = 1, Nbeta
            cand(i,ibeta,idata) = rand_dble_norm(cand(i,ibeta,idata),par(2,i,ibeta,idata))
            if(peri(i).eq.1) then
              mask_cand(ibeta,idata) = .true.
              factor                 = int(abs(cand(i,ibeta,idata))*ipi+5.d-1)
              sign                   = abs(cand(i,ibeta,idata))/cand(i,ibeta,idata)
              cand(i,ibeta,idata)    = cand(i,ibeta,idata) - pi*sign*factor
            else
              if(cand(i,ibeta,idata).gt.par(3,i,ibeta,idata) .and. cand(i,ibeta,idata).lt.par(4,i,ibeta,idata)) then
                mask_cand(ibeta,idata) = .true.
              else
                mask_cand(ibeta,idata) = .false.
              end if
            end if
          end do
        end do
        ! call kernel: calculation of chisq array
        call calc_chisq_array<<<gridDim,blockDim>>>(Ndata,Nchan,lamsq,q,u,q_err,u_err,mask_chan,src_list,Npar,Nbeta,cand,chisq_array)
        chisq = sum(chisq_array,dim=1)

        ! adoption/rejection of the candidate parameter
        do idata = 1, Ndata
          do ibeta = 1, Nbeta
            if (mask_cand(ibeta,idata)) then
              logr = -5.d-1 * beta(ibeta,idata) * (chisq(ibeta,idata) - chisq_old(ibeta,idata))

              ! accept the cand parameter with the probability equal to the Metropolis ratio
              logu = log(rand_num())
              if ( (logu < logr) .or. (logr >= 0) ) then
                chain(1:Npar,ibeta,idata) = cand(:,ibeta,idata)
                chain(Npar+1,ibeta,idata) = chisq(ibeta,idata)
              end if
            end if
          end do
        end do
      end do

      ! exchange
      if(Nbeta>1 .and. iiter>1 .and. mod(iiter+1,Nexchevry)==0) then

        do idata = 1, Ndata
          do iexch = 1, Nbeta * Nexch

            i = rand_int_unif(1,Nbeta-1)

            logr = 5.d-1 * (beta(i+1,idata) - beta(i,idata)) &
            * (chain(Npar+1,i+1,idata) - chain(Npar+1,i,idata))

            logu = log(rand_num())

            if ( (logu < logr) .or. (logr >= 0) ) then
              dum_chain(:)       = chain(:,i+1,idata)
              chain(:,i+1,idata) = chain(:,i,idata)
              chain(:,i,idata)   = dum_chain(:)
            end if
          end do
        end do
      end if

      par_chain(:,iiter,:) = chain(1:Npar,Nbeta,:)
      chisq_chain(:,iiter,:) = chain(Npar+1,:,:)
    end do

    ! output
    do idata = 1, Ndata
      open(10,file=trim(dirname)//trim(data_list(idata))//'chainS.dat',form='unformatted')
      write(10) par_chain(:,:,idata)
      close(10)
    end do

    do idata = 1, Ndata
      open(10,file=trim(dirname)//trim(data_list(idata))//'chisqS.dat',form='unformatted')
      write(10) chisq_chain(:,:,idata)
      close(10)
    end do

    par(1,:,:,:) = chain(1:Npar,:,:)
    do idata = 1, Ndata
      open(10,file=trim(dirname)//trim(data_list(idata))//'par.dat',form='unformatted')
      write(10) par(:,:,:,idata)
      close(10)
    end do

    ! deallocate
    deallocate(par_chain,chisq_chain)
  end if

  deallocate(rnd)
  err = curandDestroyGenerator(gen)

contains

  double precision function rand_num()
    use curand
    implicit none
    if (rnd_i==0 .or. rnd_i==Nrnd) then
        rnd_i = 0
        err = curandGenerate(gen, rnd, Nrnd)
    end if
    rnd_i = rnd_i + 1
    rand_num = rnd(rnd_i)
  end function rand_num

  integer function rand_int_unif(rmin,rmax)
    implicit none
    integer, intent(in) :: rmin, rmax
    rand_int_unif = int(dble(rmax-rmin+1)*rand_num()) + rmin
  end function rand_int_unif


  double precision function rand_dble_unif(rmin,rmax)
    implicit none
    double precision, intent(in) :: rmin, rmax
    rand_dble_unif = (rmax-rmin)*rand_num() + rmin
  end function rand_dble_unif


  double precision function rand_dble_norm(mu,var)
    implicit none
    double precision, intent(in) :: mu, var
    rand_dble_norm = sqrt(-2.d0*log(rand_num()))*cos(2.d0*pi*rand_num()) * var + mu
  end function rand_dble_norm

end program main