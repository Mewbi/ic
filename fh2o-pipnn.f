!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,vpesmin=-176.0260035639990d0,
     % PI=3.141592653589793238d0,radian=PI/180.0d0,bohr=0.5291772d0
       integer,parameter::nbasis=18,nrbasissoc=147,npbasissoc=288
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nterm(1:nbasis),nindex(1:nbasis,1:100,1:6)
       integer nscale
       integer, allocatable::nodes(:)
       real*8 r1soccoef(1:nrbasissoc),r2soccoef(1:nrbasissoc)
       real*8 p2soccoef(1:npbasissoc)

       real*8, allocatable::weighta(:,:,:),biasa(:,:)
       real*8, allocatable::pdela(:),pavga(:)
       real*8, allocatable::weightb(:,:,:),biasb(:,:)
       real*8, allocatable::pdelb(:),pavgb(:)
       real*8, allocatable::weightc(:,:,:),biasc(:,:)
       real*8, allocatable::pdelc(:),pavgc(:)

       end module nnparam

       subroutine fh2oNN(ct,vpes,vpesa,vpesb,vpesc) !ct in HHFO order and angstrom
       use nnparam
       implicit none
       integer,parameter::ndim=6
       integer j,k
       real*8 rb(ndim),xbond(ndim),basis(0:nbasis-1),tmp1
       real*8 txinput(1:nbasis-1)
       real*8 ct(3,4),xvec(3,ndim),xct(3,4),vpes
       real*8 vpesa,vpesb,vpesc,r12,r13,r14,r23,r24,r34

       basis=0.d0
       xct=ct
       
       xvec(:,1)=xct(:,2)-xct(:,1) !  H2->H1
       xvec(:,2)=xct(:,3)-xct(:,1) !  F->H1
       xvec(:,3)=xct(:,4)-xct(:,1) !  O->H1
       xvec(:,4)=xct(:,3)-xct(:,2) !  H2->F
       xvec(:,5)=xct(:,4)-xct(:,2) !  H2->O
       xvec(:,6)=xct(:,4)-xct(:,3) !  F-O 
       rb(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
       rb(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
       rb(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3)))
       rb(4)=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
       rb(5)=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
       rb(6)=dsqrt(dot_product(xvec(:,6),xvec(:,6))) 

       r12=rb(1);r13=rb(2);r14=rb(3);r23=rb(4);r24=rb(5);r34=rb(6)

       if(minval(rb(:)).le.0.67d0.or.minval(rb(:)).gt.1.45d0)then
         vpes=5.0d0
         return
       endif

       if(r34.le.1.5d0.or.r12.le.0.8d0.or.minval(rb).eq.r12)then
         vpes=5.0d0
         return
       endif
      
      if(r14.lt.r24.and.r24.gt.1.9d0) then
        if(r14.ge.1.8d0.or.r23.ge.1.8d0)then   
         vpes=5.0d0
         return
        endif
      endif

      if(r24.lt.r14.and.r14.gt.1.9d0) then
        if(r24.ge.1.8d0.or.r13.ge.1.8d0)then  
         vpes=5.0d0
         return
        endif
      endif

      if(r23.lt.r13.and.r13.gt.1.9d0.and.r24.gt.1.7d0) then
        if(r14.ge.1.8d0.or.r23.ge.1.8d0)then   
         vpes=5.0d0
         return
        endif
      endif

      if(r13.lt.r23.and.r23.gt.1.9d0.and.r14.gt.1.7d0) then
        if(r13.ge.1.8d0.or.r24.ge.1.8d0)then  
         vpes=5.0d0
         return
        endif
      endif

      xbond(:)=dexp(-rb(:)/alpha)
       
!      do j=1,nbasis
!       tmp1=0.0d0
!       do k=1,nterm(j)
!        tmp1=tmp1+
!    $  (xbond(1)**nindex(j,k,1))*(xbond(2)**nindex(j,k,2))
!    & *(xbond(3)**nindex(j,k,3))*(xbond(4)**nindex(j,k,4))
!    & *(xbond(5)**nindex(j,k,5))*(xbond(6)**nindex(j,k,6))
!       enddo
!       basis(j)=basis(j)+tmp1
!      enddo

       call bemsav(xbond,basis)
      
       do j=1,nbasis-1
        txinput(j)=basis(j)
       enddo

       call getpota(txinput,vpesa)
       call getpotb(txinput,vpesb)
       call getpotc(txinput,vpesc)
       vpes=(vpesa+vpesb+vpesc)/3.0d0

       if(vpes.lt.-1.5d0)vpes=5.0d0
       vpes=min(vpes,5.0d0)
       
       return
        
       end subroutine fh2oNN

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling the getpot() subroutine
        subroutine pes_init
        use nnparam
        implicit none
        integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
        integer ibasis,npd,iterm,ib,nfile
        character f1*80,line

      open(4,file='biases.txt',status='old')
      rewind(4)
      read(4,'(a100)')line
109   read(4,'(i4,3i3,3x,6i2)',end=121)ibasis,npd,nterm(ibasis+1),
     $iterm,(nindex(ibasis+1,iterm,ib), ib=1,6)
      if(nindex(ibasis+1,iterm,1).eq.2)then
        goto 121
      else
        goto 109
      endif
121   continue

      read(4,'(a100)')line
      do i=1,nrbasissoc
       read(4,*)r1soccoef(i)
      enddo

      read(4,'(a100)')line
      do i=1,npbasissoc
       read(4,*)p2soccoef(i)
      enddo
      read(4,*)line

        nfile=7
        open(nfile,file='weights.txt')
        read(nfile,*)line
!       open(nfile,file='weights.txt-13')
!       open(4,file='biases.txt-13')
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(nodes(nlayer),pdela(nscale),pavga(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weighta(nodemax,nodemax,2:nlayer),
     %   biasa(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdela(i),i=1,nscale)
        read(nfile,*)(pavga(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weighta(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasa(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif

        read(4,*)line
        read(nfile,*)line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelb(nscale),pavgb(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightb(nodemax,nodemax,2:nlayer),
     %   biasb(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelb(i),i=1,nscale)
        read(nfile,*)(pavgb(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightb(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasb(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif

        read(4,*)line
        read(nfile,*)line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelc(nscale),pavgc(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightc(nodemax,nodemax,2:nlayer),
     % biasc(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelc(i),i=1,nscale)
        read(nfile,*)(pavgc(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightc(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasc(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        close(nfile)
        close(4)

        return

        end subroutine pes_init

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)
        return
        end subroutine getpota

        subroutine getpotb(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavgb(i))/pdelb(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightb(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightb(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelb(nscale)+pavgb(nscale)
        return
        end subroutine getpotb

        subroutine getpotc(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavgc(i))/pdelc(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightc(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weightc(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)
        return
        end subroutine getpotc

        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
c    ifunc=2, transfer function is log sigmoid function, 'logsig'
c    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if (ifunc.eq.1) then
        tranfun=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+exp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end

      function emsav(x,c) result(v)
      implicit none
      real*8,dimension(1:6)::x
      real*8,dimension(0:17)::c
      real*8::v
      ! ::::::::::::::::::::
      real*8,dimension(0:17)::p
      call bemsav(x,p)
      v = dot_product(p,c)
      return
      end function emsav
 
      subroutine bemsav(x,p)
      implicit none
      real*8,dimension(1:6),intent(in)::x
      real*8,dimension(0:17),intent(out)::p
      ! ::::::::::::::::::::
      real*8,dimension(0:10)::m
      call evmono(x,m)
      call evpoly(m,p)
      return
      end subroutine bemsav
 
      subroutine evmono(x,m)
      implicit none
      real*8,dimension(1:6),intent(in)::x
      real*8,dimension(0:10),intent(out)::m
 
      m(0)=1.D0
      m(1)=x(6)
      m(2)=x(5)
      m(3)=x(3)
      m(4)=x(4)
      m(5)=x(2)
      m(6)=x(1)
      m(7)=m(2)*m(3)
      m(8)=m(3)*m(4)
      m(9)=m(2)*m(5)
      m(10)=m(4)*m(5)
 
      return
      end subroutine evmono
 
      subroutine evpoly(m,p)
      implicit none
      real*8,dimension(0:10),intent(in)::m
      real*8,dimension(0:17),intent(out)::p
 
      p(0)=m(0)
      p(1)=m(1)
      p(2)=m(2)+m(3)
      p(3)=m(4)+m(5)
      p(4)=m(6)
      p(5)=p(1)*p(2)
      p(6)=m(7)
      p(7)=p(1)*p(3)
      p(8)=m(8)+m(9)
      p(9)=m(10)
      p(10)=p(2)*p(3)-p(8)
      p(11)=p(1)*p(4)
      p(12)=p(4)*p(2)
      p(13)=p(4)*p(3)
      p(14)=p(1)*p(1)
      p(15)=p(2)*p(2)-p(6)-p(6)
      p(16)=p(3)*p(3)-p(9)-p(9)
      p(17)=p(4)*p(4)
 
      return
      end subroutine evpoly
