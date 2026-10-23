c$Id:$
c-----[--.----+----.----+----.-----------------------------------------]
c     elmt16.f — Simplified copy of elmt15.f
c     3-D finite deformation Neo-Hookean solid (Updated Lagrangian)
c     NO damage, NO nonlocal / gradient-enhanced field
c
c     Control:   ndm = 3, ndf = 3, nen = 8
c     Material:  user 16
c                E, nu
c-----[--.----+----.----+----.-----------------------------------------]
      subroutine elmt16(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c     * * F E A P * * A Finite Element Analysis Program

      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eldatp.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'strnum.h'

      logical   errck
      integer   ndf,ndm,nst,isw, ix(*)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,nst),p(nst)

      real*8    td(5)
      real*8    e1,nu12,g12,blk

      save

c     Output element type
      if(isw.eq.0 .and. ior.lt.0) then
         write(*,*) '   Elmt  16: 3-D Finite Deformation Neo-Hookean'
         write(*,*) '             (no damage / no nonlocal field)'

c     Input material properties
      elseif(isw.eq.1) then
         istv  = 6
         pstyp = 3
         nh1   = 0
         write(iow,2000)
         if(ior.lt.0) write(*,2000)

         call dinput(td,2)

         e1   = td(1)
         nu12 = td(2)

         g12 = e1/(1.d0 + nu12)/2.d0
         blk = e1/(1.d0 - nu12*2.d0)/3.d0

         d(1)  = e1
         d(2)  = nu12
         d(21) = blk - 2.d0/3.d0*g12
         d(22) = g12

         if(ior.lt.0) then
            write(*,2010) e1,nu12,d(21),g12
         endif
         write(iow,2010) e1,nu12,d(21),g12

         errck = .false.
         if(ndf.lt.3 .or. nen.lt.8) then
            write(  *,3000) ndf,nen
            write(ilg,3000) ndf,nen
            write(iow,3000) ndf,nen
            errck = .true.
         endif
         if(errck) call plstop()

c     Residual / tangent / stress / projection
      elseif(isw.eq.3 .or. isw.eq.4 .or. isw.eq.6 .or. isw.eq.8) then
         call NeoHmodel(d,ul,xl,s,p,ndf,ndm,nst,isw)
      endif

 2000 format(
     &     /5x,'T h r e e     D i m e n s i o n a l       S o l i d'/
     &     'N e o - H o o k e a n   F i n i t e   D e f o r m a t i o n'
     &     /)

 2010 format(5x,'Model Parameter'//,
     &     10x,'Modulus E       ',1p,1e12.5/
     &     10x,'Poisson Ratio   ',0p,1f8.5/
     &     10x,'Lambda Parameter',1p,1e12.5/
     &     10x,'Shear Modulus   ',1p,1e12.5/1x)

 3000 format(' *ERROR* Problem control record incorrect:'/
     &     '         DOFs/Node (ndf) = ',i4,' - Should be 3 or more'/
     &     '         Nodes/Elm (nen) = ',i4,' - Should be 8')

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Updated Lagrangian Neo-Hookean driver (8-node hex, 2x2x2 Gauss)
c-----[--.----+----.----+----.-----------------------------------------]
      subroutine NeoHmodel(d,ul,xl,s,r,ndf,ndm,nst,isw)

      implicit none

      integer, parameter :: mgauss = 8
      integer, parameter :: nnode  = 8
      integer, parameter :: ntens  = 6

      real*8,  parameter :: zero = 0.d0

      include  'eldata.h'
      include  'cdata.h'
      include  'elplot.h'
      include  'iofile.h'
      include  'strnum.h'

      integer   ndm                 ! spatial dimension
      integer   ndf                 ! DOFs per node
      integer   nst                 ! element matrix size
      integer   isw                 ! FEAP solution switch
      integer   i,j,kk,mm,pp,qq     ! loop indices
      integer   intip               ! Gauss-point counter
      real*8    d(*)                ! material parameters
      real*8    ul(ndf,nen,*)       ! nodal solution (displacements)
      real*8    xl(ndm,*)           ! nodal reference coordinates
      real*8    s(nst,nst)          ! element tangent matrix
      real*8    r(nst)              ! element residual vector

      real*8    pgp(3,mgauss)       ! Gauss-point natural coords
      real*8    wgp(mgauss)         ! Gauss weights
      real*8    ksi(3)              ! current GP natural coords
      real*8    shapef(nnode)       ! shape functions N_I
      real*8    dshape(nnode,3)     ! dN/dxi (parent derivatives)
      real*8    shp3(nnode,mgauss)  ! N_I stored at each GP (plots)
      real*8    xsjdv(mgauss)       ! volume factor J*w at each GP

      real*8    xreff(3,nnode)      ! reference nodal coords X
      real*8    xcurr(3,nnode)      ! current nodal coords x = X+u
      real*8    xjacR(3,3)          ! Jacobian dX/dxi (reference)
      real*8    xjacinvR(3,3)       ! inverse of xjacR
      real*8    detxjacR            ! det(xjacR)
      real*8    xjacC(3,3)          ! Jacobian dx/dxi (current)
      real*8    xjacinvC(3,3)       ! inverse of xjacC
      real*8    detxjacC            ! det(xjacC)

      real*8    Fgrad(3,3)          ! deformation gradient F
      real*8    Fgradinv(3,3)       ! F^{-1}
      real*8    FgradT(3,3)         ! F^T
      real*8    detF                ! J = det(F)
      real*8    bTensor(3,3)        ! left Cauchy-Green b = F F^T
      real*8    Bmat(6,3*nnode)     ! spatial B-matrix (Voigt)

      real*8    stressE(3,3)        ! Cauchy stress sigma
      real*8    CmatE(3,3,3,3)      ! spatial elastic moduli c_ijkl
      real*8    CmatEVoigt(6,6)     ! moduli in Voigt notation
      real*8    sigVoigt(6)         ! Cauchy stress in Voigt
      real*8    sig8x8(6,mgauss)    ! stresses at all GPs (output)

      real*8    fact                ! integration weight det(xjacC)*w
      real*8    Fd(3*nnode)         ! assembled internal force
      real*8    felem(3*nnode)      ! GP contribution to force
      real*8    Kdd(3*nnode,3*nnode)     ! material stiffness at GP
      real*8    Kgeodd(3*nnode,3*nnode)  ! geometric stiffness at GP
      real*8    sumKdd(3*nnode,3*nnode)  ! accumulated K_uu

      integer   Nota(3,3)
      data Nota/1,4,6,
     *          4,2,5,
     *          6,5,3/

c     Initialize
      sumKdd = zero
      Fd     = zero
      xreff  = zero
      xcurr  = zero

      if(nel.ne.8) then
         write(*,*)    'ERROR in elmt16.f - only 8-node hex supported'
         write(iow,*)  'ERROR in elmt16.f - only 8-node hex supported'
         return
      endif

c     Integration points and weights
      call nhQuadrature3D8(pgp,wgp) ! 8-node hex, 2x2x2 Gauss

c     Reference / current geometry
      ! Reference geometry
      do i = 1,ndm
         do j = 1,nel
            xreff(i,j) = xl(i,j)
            xcurr(i,j) = xl(i,j) + ul(i,j,1)
         end do
      end do

c     -----------------------------------------------------------------
c     LOOP OVER INTEGRATION POINTS
c     -----------------------------------------------------------------
      do intip = 1,mgauss

         do i = 1,3
            ksi(i) = pgp(i,intip)
         end do

c        Evaluate shape functions
         call nhShape3D8(shapef,ksi,ndm)
         do i = 1,nel
            shp3(i,intip) = shapef(i)
         enddo
c        Evaluate shape function derivatives
         call nhDerivShape3D8(dshape,ksi,ndm)

c        Compute Jacobian matrices (reference and current)
         call nhJacobian3D(xjacR,detxjacR,xjacinvR,
     *        dshape,ndm,nel,xreff)
         call nhJacobian3D(xjacC,detxjacC,xjacinvC,
     *        dshape,ndm,nel,xcurr)

c        Compute deformation gradient tensor
         call nhDeformGrad3D(Fgrad,Fgradinv,FgradT,detF,
     *        xjacR,detxjacR,xjacinvR,xjacC,detxjacC,xjacinvC,ndm,nel)

c        Compute B matrix
         call nhBmat3D(Bmat,dshape,xjacC,detxjacC,xjacinvC,ndm,nel,ndm)

c        Compute stress tensor
         call nhStandardNeoh3(d,detF,Fgrad,stressE,CmatE,bTensor,ndm)

         fact = detxjacC*wgp(intip)
         xsjdv(intip) = fact

c        Tensor -> Voigt (11,22,33,12,23,13)
         do kk = 1,3
            do mm = 1,3
               do pp = 1,3
                  do qq = 1,3
                     CmatEVoigt(Nota(kk,mm),Nota(pp,qq)) =
     $                    CmatE(kk,mm,pp,qq)
                  end do
               end do
            end do
         end do

         sigVoigt(1) = stressE(1,1)
         sigVoigt(2) = stressE(2,2)
         sigVoigt(3) = stressE(3,3)
         sigVoigt(4) = stressE(1,2)
         sigVoigt(5) = stressE(2,3)
         sigVoigt(6) = stressE(1,3)

         do i = 1,6
            sig8x8(i,intip) = sigVoigt(i)
         enddo

         call nhGetKuu_Fu(Kdd,felem,nnode,ndm,ndm,fact,CmatEVoigt,
     *        Bmat,sigVoigt,ntens)
         call nhGetKuuGeom(Kgeodd,dshape,stressE,ndm,nnode,ntens,
     *        xjacinvC)

         do i = 1,ndm*nnode
            Fd(i) = Fd(i) - felem(i)
            do j = 1,ndm*nnode
               sumKdd(i,j) = sumKdd(i,j) + Kdd(i,j) + Kgeodd(i,j)*fact
            enddo
         enddo

      end do

c     Scatter to element arrays (ndf = 3, contiguous)
      do i = 1,ndm*nnode
         r(i) = Fd(i)
         do j = 1,ndm*nnode
            s(i,j) = sumKdd(i,j)
         enddo
      enddo

c     Stress projection for plots
      if(isw.eq.8) then
         call nhGaussValueProjection(sig8x8,shp3,xsjdv,
     *        mgauss,nel,6,r,s)
      endif

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compressible Neo-Hookean (spatial / Cauchy)
c     d: material parameters (input)
c     detf: determinant of deformation gradient tensor (input)
c     fgrad: deformation gradient tensor (input)
c     stressE: Kirchhoff stress tensor (output)
c     CmatE: constitutive matrix (output)
c     btensor: Left Cauchy-Green deformation tensor matrix (output)
c     ndim: dimension of the problem (input)
c     ____________________________________________________________
c-----[--.----+----.----+----.-----------------------------------------]
      subroutine nhStandardNeoh3(d,detf,fgrad,stressE,CmatE,btensor,ndim)

      implicit none

      real*8, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

      integer   i,j,k,l,ndim
      real*8    detf, d(*)
      real*8    btensor(3,3), CmatE(3,3,3,3), stressE(3,3)
      real*8    g12, lambda
      real*8    fgrad(3,3), fgradT(3,3)
      real*8    Id2(3,3), II4(3,3,3,3)

      do i = 1,3
         do j = 1,3
            Id2(i,j)     = zero
            btensor(i,j) = zero
         enddo
         Id2(i,i) = one
      enddo

      lambda = d(21)
      g12    = d(22)

      call nhTranspose3D(fgrad,3,3,fgradT)
      call nhMatrixMultiply(btensor,ndim,ndim,
     *     fgrad,ndim,ndim,fgradT,ndim,ndim)

c     Compute Kirchhoff stress tensor
      do i = 1,ndim
         do j = 1,ndim
            stressE(i,j) = g12*btensor(i,j)/detf +
     1           (lambda*dlog(detf)/detf - g12/detf)*Id2(i,j)
         enddo
      enddo

c     Compute constitutive matrix
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  CmatE(i,j,k,l) = zero
                  II4(i,j,k,l) = half*(Id2(i,k)*Id2(j,l) +
     $                 Id2(i,l)*Id2(j,k))
               enddo
            enddo
         enddo
      enddo

c     Compute constitutive matrix
      do i = 1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  CmatE(i,j,k,l) = lambda*Id2(i,j)*Id2(k,l)/detf +
     1                 two*(g12 - lambda*dlog(detf))*II4(i,j,k,l)/detf
               enddo
            enddo
         enddo
      enddo

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute material stiffness matrix and internal force vector
c     stiff: material stiffness matrix (output)
c     force: internal force vector (output)
c     nnode: number of nodes per element (input)
c     ndim: dimension of the problem (input)
c     ndof: number of degrees of freedom per node (input)
c     factor: integration weight (input)
c     Dmat: constitutive matrix (input)
c     Bmat: B matrix (input)
c     stress: stress tensor (input)
c     ntens: number of stress components (input)
c     ____________________________________________________________
      subroutine nhGetKuu_Fu(stiff,force,nnode,ndim,ndof,factor,
     *     Dmat,Bmat,stress,ntens)

      implicit none

      real*8, parameter :: zero = 0.d0

      integer nnode,ndim,ndof,ntens
      real*8  stiff(ndim*nnode,ndim*nnode)
      real*8  force(ndim*nnode)
      real*8  Bmat(ntens,ndim*nnode)
      real*8  Dmat(ntens,ntens)
      real*8  stress(ntens)
      real*8  factor

      real*8  db(ntens), sum, sumforce
      integer ndofnnode, ndimnnode
      integer i,j,k,l,m, istiff,jstiff, iforce,jforce

      ndofnnode = ndof*nnode
      ndimnnode = ndim*nnode

      do istiff = 1,ndofnnode
         force(istiff) = zero
         do jstiff = 1,ndofnnode
            stiff(istiff,jstiff) = zero
         end do
      end do

      do j = 1,ndimnnode
         do k = 1,ntens
            db(k) = zero
            do l = 1,ntens
               db(k) = db(k) + Dmat(k,l)*Bmat(l,j)*factor
            end do
         end do
         do i = 1,ndimnnode
            sum = zero
            do m = 1,ntens
               sum = sum + Bmat(m,i)*db(m)
            end do
            stiff(i,j) = stiff(i,j) + sum
         end do
      end do

      do jforce = 1,ndimnnode
         sumforce = zero
         do iforce = 1,ntens
            sumforce = sumforce + Bmat(iforce,jforce)*stress(iforce)*
     $           factor
         end do
         force(jforce) = sumforce
      end do

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute geometric stiffness matrix
c     Kgeodd: geometric stiffness matrix (output)
c     dshape: shape function derivatives (input)
c     stress: stress tensor (input)
c     ndim: dimension of the problem (input)
c     nnode: number of nodes per element (input)
c     ntens: number of stress components (input)
c     xjacinv: inverse of Jacobian matrix (input)
c     ____________________________________________________________
      subroutine nhGetKuuGeom(Kgeodd,dshape,stress,ndim,nnode,
     *     ntens,xjacinv)

      implicit none

      real*8, parameter :: zero = 0.d0

      integer ndim,nnode,ntens
      integer i,j,k,l,kk,ll
      real*8  dshape(nnode,ndim)
      real*8  dshapeT(ndim,nnode), dshapeC(ndim,nnode)
      real*8  Kgeodd(ndim*nnode,ndim*nnode)
      real*8  stress(ndim,ndim)
      real*8  xjacinv(ndim,ndim)
      real*8  scalar
      integer ndimnnode

      ndimnnode = ndim*nnode

      do i = 1,ndimnnode
         do j = 1,ndimnnode
            Kgeodd(i,j) = zero
         enddo
      enddo

      do i = 1,ndim
         do j = 1,nnode
            dshapeC(i,j) = zero
         enddo
      enddo

      call nhTranspose3D(dshape,ndim,nnode,dshapeT)

      do i = 1,nnode
         dshapeC(1,i) = dshapeT(1,i)*xjacinv(1,1) +
     *        dshapeT(2,i)*xjacinv(2,1) +
     *        dshapeT(3,i)*xjacinv(3,1)
         dshapeC(2,i) = dshapeT(1,i)*xjacinv(1,2) +
     *        dshapeT(2,i)*xjacinv(2,2) +
     *        dshapeT(3,i)*xjacinv(3,2)
         dshapeC(3,i) = dshapeT(1,i)*xjacinv(1,3) +
     *        dshapeT(2,i)*xjacinv(2,3) +
     *        dshapeT(3,i)*xjacinv(3,3)
      enddo

      kk = 1
      do k = 1,nnode
         ll = 1
         do l = 1,nnode
            scalar = zero
            do i = 1,ndim
               do j = 1,ndim
                  scalar = scalar + dshapeC(i,k)*stress(i,j)*dshapeC(j,l)
               enddo
            enddo
            Kgeodd(kk  ,ll  ) = scalar
            Kgeodd(kk+1,ll+1) = scalar
            Kgeodd(kk+2,ll+2) = scalar
            ll = ll + 3
         enddo
         kk = kk + 3
      enddo

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute Gauss value projection
c     sig: stress tensor (input)
c     shp3: shape functions (input)
c     xsj: Jacobian determinant (input)
c     lint: number of integration points (input)
c     nel: number of nodes per element (input)
c     nes: number of stress components (input)
c     r: internal force vector (output)
c     s: material stiffness matrix (output)
c     ____________________________________________________________  
      subroutine nhGaussValueProjection(sig,shp3,xsj,lint,nel,nes,r,s)

      implicit none

      include  'cdata.h'
      include  'strnum.h'

      integer   ii,jj,l,nel,lint,nes
      real*8    r(*),s(nen,*), sig(nes,lint), xj
      real*8    xsj(lint), shp3(nel,*)

      do ii = 1,nel
         r(ii) = 0.0d0
         do jj = 1,nes
            s(ii,jj) = 0.0d0
         end do
      end do

      do l = 1,lint
         do ii = 1,nel
            xj    = xsj(l)*shp3(ii,l)
            r(ii) = r(ii) + xj
            do jj = 1,nes
               s(ii,jj) = s(ii,jj) + sig(jj,l)*xj
            end do
         end do
      end do

      iste = nes

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Evaluate shape functions
c     N: shape functions (output)
c     xi: isoparametric coordinates (input)
c     ndim: dimension of the problem (input)
c     ____________________________________________________________
      subroutine nhShape3D8(N,xi,ndim)

      implicit none

      real*8, parameter :: one=1.d0, eight=8.d0
      integer ndim
      real*8  N(8), xi(ndim)

      N(1) = (one-xi(1))*(one-xi(2))*(one-xi(3))/eight
      N(2) = (one+xi(1))*(one-xi(2))*(one-xi(3))/eight
      N(3) = (one+xi(1))*(one+xi(2))*(one-xi(3))/eight
      N(4) = (one-xi(1))*(one+xi(2))*(one-xi(3))/eight
      N(5) = (one-xi(1))*(one-xi(2))*(one+xi(3))/eight
      N(6) = (one+xi(1))*(one-xi(2))*(one+xi(3))/eight
      N(7) = (one+xi(1))*(one+xi(2))*(one+xi(3))/eight
      N(8) = (one-xi(1))*(one+xi(2))*(one+xi(3))/eight

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Evaluate shape function derivatives
c     dNdxi: shape function derivatives (output)
c     xi: isoparametric coordinates (input)
c     ndim: dimension of the problem (input)
c     ____________________________________________________________
      subroutine nhDerivShape3D8(dNdxi,xi,ndim)

      implicit none

      real*8, parameter :: one=1.d0, eight=8.d0
      integer ndim
      real*8  dNdxi(8,ndim), xi(ndim)

      dNdxi(1,1) = -(one-xi(2))*(one-xi(3))/eight
      dNdxi(1,2) = -(one-xi(1))*(one-xi(3))/eight
      dNdxi(1,3) = -(one-xi(1))*(one-xi(2))/eight
      dNdxi(2,1) =  (one-xi(2))*(one-xi(3))/eight
      dNdxi(2,2) = -(one+xi(1))*(one-xi(3))/eight
      dNdxi(2,3) = -(one+xi(1))*(one-xi(2))/eight
      dNdxi(3,1) =  (one+xi(2))*(one-xi(3))/eight
      dNdxi(3,2) =  (one+xi(1))*(one-xi(3))/eight
      dNdxi(3,3) = -(one+xi(1))*(one+xi(2))/eight
      dNdxi(4,1) = -(one+xi(2))*(one-xi(3))/eight
      dNdxi(4,2) =  (one-xi(1))*(one-xi(3))/eight
      dNdxi(4,3) = -(one-xi(1))*(one+xi(2))/eight
      dNdxi(5,1) = -(one-xi(2))*(one+xi(3))/eight
      dNdxi(5,2) = -(one-xi(1))*(one+xi(3))/eight
      dNdxi(5,3) =  (one-xi(1))*(one-xi(2))/eight
      dNdxi(6,1) =  (one-xi(2))*(one+xi(3))/eight
      dNdxi(6,2) = -(one+xi(1))*(one+xi(3))/eight
      dNdxi(6,3) =  (one+xi(1))*(one-xi(2))/eight
      dNdxi(7,1) =  (one+xi(2))*(one+xi(3))/eight
      dNdxi(7,2) =  (one+xi(1))*(one+xi(3))/eight
      dNdxi(7,3) =  (one+xi(1))*(one+xi(2))/eight
      dNdxi(8,1) = -(one+xi(2))*(one+xi(3))/eight
      dNdxi(8,2) =  (one-xi(1))*(one+xi(3))/eight
      dNdxi(8,3) =  (one-xi(1))*(one+xi(2))/eight

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute Jacobian matrix (isoparametric mapping)
c     xjac: Jacobian matrix (output)
c     detxjac: determinant of Jacobian matrix (output)
c     xjacinv: inverse of Jacobian matrix (output)
c     dshape: shape function derivatives (input)
c     ndim: dimension of the problem (input)
c     nnode: number of nodes per element (input) 
c     coord: node coordinates (input)
c     ____________________________________________________________
      subroutine nhJacobian3D(xjac,detxjac,xjacinv,dshape,ndim,nnode,coord)

      implicit none

      include  'iofile.h'

      real*8, parameter :: zero=0.d0
      integer ndim,nnode, i,j,k
      real*8  xjac(ndim,ndim), xjacinv(ndim,ndim), detxjac
      real*8  coord(ndim,nnode), dshape(nnode,ndim)

      do i = 1,ndim
         do j = 1,ndim
            xjac(i,j)    = zero
            xjacinv(i,j) = zero
         enddo
      enddo

      do k = 1,nnode
         xjac(1,1) = xjac(1,1) + coord(1,k)*dshape(k,1)
         xjac(1,2) = xjac(1,2) + coord(1,k)*dshape(k,2)
         xjac(1,3) = xjac(1,3) + coord(1,k)*dshape(k,3)
         xjac(2,1) = xjac(2,1) + coord(2,k)*dshape(k,1)
         xjac(2,2) = xjac(2,2) + coord(2,k)*dshape(k,2)
         xjac(2,3) = xjac(2,3) + coord(2,k)*dshape(k,3)
         xjac(3,1) = xjac(3,1) + coord(3,k)*dshape(k,1)
         xjac(3,2) = xjac(3,2) + coord(3,k)*dshape(k,2)
         xjac(3,3) = xjac(3,3) + coord(3,k)*dshape(k,3)
      enddo

      call nhDeterminant3D(detxjac,xjac)

      if(detxjac.gt.zero) then
         call nhInverse3D(xjacinv,xjac,detxjac)
      else
         write(iow,*)'Warning: element has negative jacobian'
         write(*,*)  'Warning: element has negative jacobian'
      end if

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute determinant of a 3x3 matrix
c     det: determinant of the matrix (output)
c     A: matrix (input)
c     ____________________________________________________________
      subroutine nhDeterminant3D(det,A)

      implicit none

      real*8, parameter :: zero=0.d0
      real*8 A(3,3), det

      det = A(1,1)*A(2,2)*A(3,3) +
     *     A(1,3)*A(2,1)*A(3,2) +
     *     A(3,1)*A(1,2)*A(2,3) -
     *     A(3,1)*A(2,2)*A(1,3) -
     *     A(3,3)*A(1,2)*A(2,1) -
     *     A(1,1)*A(2,3)*A(3,2)

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute inverse of a 3x3 matrix
c     Ainv: inverse of the matrix (output)
c     A: matrix (input)
c     det: determinant of the matrix (input)
c     ____________________________________________________________
      subroutine nhInverse3D(Ainv,A,det)

      implicit none

      real*8 A(3,3), Ainv(3,3), det

      Ainv(1,1) =  (A(2,2)*A(3,3) - A(3,2)*A(2,3))/det
      Ainv(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))/det
      Ainv(1,3) =  (A(1,2)*A(2,3) - A(2,2)*A(1,3))/det
      Ainv(2,1) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))/det
      Ainv(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))/det
      Ainv(2,3) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))/det
      Ainv(3,1) =  (A(2,1)*A(3,2) - A(3,1)*A(2,2))/det
      Ainv(3,2) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))/det
      Ainv(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute deformation gradient tensor (reference to current configuration)
c     Fgrad: deformation gradient tensor (output)
c     Finv: inverse of deformation gradient tensor (output)
c     FgradT: transpose of deformation gradient tensor (output)
c     detF: determinant of deformation gradient tensor (output)
c     xjacR: Jacobian matrix of reference geometry (input)
c     detxjacR: determinant of Jacobian matrix of reference geometry (input)
c     xjacinvR: inverse of Jacobian matrix of reference geometry (input)
c     xjacC: Jacobian matrix of current geometry (input)
c     detxjacC: determinant of Jacobian matrix of current geometry (input)
c     xjacinvC: inverse of Jacobian matrix of current geometry (input)
c     ndim: dimension of the problem (input)
c     nnode: number of nodes per element (input)
c     ____________________________________________________________
      subroutine nhDeformGrad3D(Fgrad,Finv,FgradT,detF,xjacR,detxjacR,
     *     xjacinvR,xjacC,detxjacC,xjacinvC,ndim,nnode)

      implicit none

      real*8, parameter :: zero=0.d0
      integer ndim,nnode, i,j
      real*8  Fgrad(3,3), Finv(3,3), FgradT(3,3), detF
      real*8  xjacR(ndim,ndim), xjacinvR(ndim,ndim), detxjacR
      real*8  xjacC(ndim,ndim), xjacinvC(ndim,ndim), detxjacC

      do i = 1,ndim
         do j = 1,ndim
            Fgrad(i,j)  = zero
            Finv(i,j)   = zero
            FgradT(i,j) = zero
         enddo
      enddo

      ! Compute deformation gradient tensor as the product of the Jacobian matrix 
      ! of the current configuration and the inverse of the Jacobian matrix of the 
      ! reference configuration
      call nhMatrixMultiply(Fgrad,ndim,ndim,xjacC,ndim,ndim,
     *     xjacinvR,ndim,ndim)


      ! Compute determinant of deformation gradient tensor
      detF = detxjacC/detxjacR

      ! Compute inverse of deformation gradient tensor
      call nhInverse3D(Finv,Fgrad,detF)

      ! Compute transpose of deformation gradient tensor
      call nhTranspose3D(Fgrad,3,3,FgradT)

      end

c-----[--.----+----.----+----.-----------------------------------------]
      subroutine nhMatrixMultiply(A,a1,a2,B,b1,b2,C,c1,c2)

      implicit none

      integer a1,a2,b1,b2,c1,c2, i,j,k
      real*8  A(a1,a2), B(b1,b2), C(c1,c2)
      real*8, parameter :: zero = 0.d0
      real*8  sump

      do i = 1,a1
         do j = 1,a2
            A(i,j) = zero
         enddo
      enddo

      if(b2.ne.c1) then
         write(*,*) 'No matrix multiplication'
      end if

      do i = 1,b1
         do j = 1,c2
            sump = zero
            do k = 1,c1
               sump = sump + B(i,k)*C(k,j)
            end do
            A(i,j) = sump
         end do
      end do

      end

c-----[--.----+----.----+----.-----------------------------------------]
      subroutine nhTranspose3D(A,m,n,AT)

      implicit none

      real*8, parameter :: zero=0.d0
      integer n,m, i,j
      real*8  A(n,m), AT(m,n)

      do i = 1,m
         do j = 1,n
            AT(i,j) = zero
         end do
      end do

      do i = 1,n
         do j = 1,m
            AT(j,i) = A(i,j)
         end do
      end do

      end

c-----[--.----+----.----+----.-----------------------------------------]
c     Compute B matrix
c     Bmat: B matrix (output)
c     dshape: shape function derivatives (input)
c     xjac: Jacobian matrix (input)
c     detxjac: determinant of Jacobian matrix (input)
c     xjacinv: inverse of Jacobian matrix (input)
c     ndim: dimension of the problem (input)
c     nnode: number of nodes per element (input)
c     ndof: number of degrees of freedom per node (input)
c     ____________________________________________________________
      subroutine nhBmat3D(Bmat,dshape,xjac,detxjac,xjacinv,ndim,
     *     nnode,ndof)

      implicit none

      real*8, parameter :: zero=0.d0
      integer ndim,nnode,ndof, i,j,kint
      real*8  Bmat(6,nnode*ndof)
      real*8  dNidx, dNidy, dNidz
      real*8  xjac(ndim,ndim), xjacinv(ndim,ndim), detxjac
      real*8  dshape(nnode,ndim)

      do i = 1,6
         do j = 1,nnode*ndof
            Bmat(i,j) = zero
         enddo
      enddo

      kint = 1
      do i = 1,nnode
         dNidx = dshape(i,1)*xjacinv(1,1) +
     1        dshape(i,2)*xjacinv(2,1) +
     2        dshape(i,3)*xjacinv(3,1)
         dNidy = dshape(i,1)*xjacinv(1,2) +
     1        dshape(i,2)*xjacinv(2,2) +
     2        dshape(i,3)*xjacinv(3,2)
         dNidz = dshape(i,1)*xjacinv(1,3) +
     1        dshape(i,2)*xjacinv(2,3) +
     2        dshape(i,3)*xjacinv(3,3)

         Bmat(1,kint)   = dNidx
         Bmat(2,kint+1) = dNidy
         Bmat(3,kint+2) = dNidz
         Bmat(4,kint)   = dNidy
         Bmat(4,kint+1) = dNidx
         Bmat(5,kint+1) = dNidz
         Bmat(5,kint+2) = dNidy
         Bmat(6,kint)   = dNidz
         Bmat(6,kint+2) = dNidx

         kint = kint + 3
      end do

      end

c-----[--.----+----.----+----.-----------------------------------------]
      subroutine nhQuadrature3D8(pgp,wgp)

      implicit none

      real*8, parameter :: one = 1.d0, three = 3.d0
      integer i
      real*8  pgp(3,8), wgp(8), scalar

      do i = 1,8
         wgp(i) = one
      enddo

      scalar = one/dsqrt(three)

      pgp(1,1) = -scalar
      pgp(2,1) = -scalar
      pgp(3,1) = -scalar

      pgp(1,2) =  scalar
      pgp(2,2) = -scalar
      pgp(3,2) = -scalar

      pgp(1,3) =  scalar
      pgp(2,3) =  scalar
      pgp(3,3) = -scalar

      pgp(1,4) = -scalar
      pgp(2,4) =  scalar
      pgp(3,4) = -scalar

      pgp(1,5) = -scalar
      pgp(2,5) = -scalar
      pgp(3,5) =  scalar

      pgp(1,6) =  scalar
      pgp(2,6) = -scalar
      pgp(3,6) =  scalar

      pgp(1,7) =  scalar
      pgp(2,7) =  scalar
      pgp(3,7) =  scalar

      pgp(1,8) = -scalar
      pgp(2,8) =  scalar
      pgp(3,8) =  scalar

      end
