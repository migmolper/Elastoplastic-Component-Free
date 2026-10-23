c$Id:$
      subroutine elmt15(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c     * * F E A P * * A Finite Element Analysis Program

c.... Copyright (c) 1984-2008: Regents of the University of California
c     All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     for git purposes
c     Modification log                                Date (dd/mm/year)
c     Original version                                    01/11/2006
c     Implementation of elmt14                           09/09/2014 
c     GED model implemented                              03/13/2017 
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Three Dimensional Gradient Enhanced 
c     non-local damage model 
c     Inputs:
c     d(*)      - Element parameters
c     ul(ndf,*) - Current nodal solution parameters
c     xl(ndm,*) - Nodal coordinates
c     ix(*)     - Global nodal connections
c     ndf       - Degree of freedoms/node
c     ndm       - Mesh coordinate dimension
c     nst       - Element array dimension
c     isw       - Solution option switch

c     Outputs:
c     s(nst,*)  - Element array
c     p(ndf,*)  - Element vector
      
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'cdata.h'        !numnp,numel,nummat,nen,neq,ipr
      include  'eldata.h'       !n,ma,mct,iel,nel,pstyp,eltyp,elty2,elty3
      include  'eldatp.h' 
      include  'hdata.h'
      include  'comblk.h'
      include  'iofile.h'
      include  'strnum.h'

      LOGICAL   errck
      INTEGER   ndf,ndm,nst,isw, ix(*)
      REAL*8    d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,nst),p(nst) 

      REAL*8    td(15)          !input parameter vector
      REAL*8    e1,nu12,g12,blk !elastic modulus, poisson's ratio
c     User defined
      REAL*8    Kd, ndp               
      INTEGER   etype
      INTEGER   mgauss
      
      save
      
c     Extract type data                     _________
c     Gradient-Enhanced Damage Model flag, |_etype=1_|
c     Standard Damage Model flag, |_etype=2_| 
      etype = 2                       

c     Output element type
      if(isw.eq.0 .and. ior.lt.0) then  
         write(*,*) '   Elmt  14: 3-D Finite Deformation Element'
         write(*,*) ' Standard and Gradient Enhanced Damage Model'

c     Input material properties

      elseif(isw.eq.1) then
         istv = 8 
         pstyp = 3
         write(iow,2000)
         if(ior.lt.0) write(*,2000)
c     
         call dinput(td,5)
c     
         e1      = td(1)
         nu12    = td(2)
         Kd      = td(3)  
         ndp     = td(4)  
c     
         etype = td(5)                        

         g12 = e1/(1.d0 + nu12)/2.d0
         blk = e1/(1.d0 - nu12*2.d0)/3.d0
c     
         d(1)  = e1
         d(2)  = nu12
         d(21) = blk - 2.d0/3.d0*g12 ! lambda 
         d(22) = g12
         d(243)= Kd                        
         d(244)= ndp
         d(245)= etype
c     
         if(ior.lt.0) then
            write(*,2010) e1,nu12,d(21),g12,Kd,ndp
         endif
         
         write(iow,2010) e1,nu12,d(21),g12,Kd,ndp

         if(etype.eq.1) then             
            mgauss = 8          ! Gauss Quadrature 2x2x2 rule per element
            nh1 = mgauss
         elseif(etype.eq.2) then
            mgauss = 8
            nh1 = mgauss
         endif

c     Set plot sequence

c     Check for errors in problem size
         errck = .false.
         if(ndf.lt.3 .or. nen.lt.4) then
            write(  *,3000) ndf,nen
            write(ilg,3000) ndf,nen
            write(iow,3000) ndf,nen
            errck = .true.
         endif

         if(etype.eq.0 .or. etype.ge.3) then
            write(  *,4000)
            write(ilg,4000)
            write(iow,4000)
            errck = .true.
         endif

         if(errck) then
            call plstop()
         endif
         
c     Check element for errors

      elseif(isw.eq.2) then

c     Compute mass or geometric stiffness matrix

      elseif(isw.eq.5) then

c     Mass contribution

c     Compute damping matrix

      elseif(isw.eq.9) then

c     History manipulation

      elseif(isw.eq.12) then

         RETURN                 ! (R)-required step, so return is used. 

c     Body Force computation

      elseif(isw.eq.15 .or.isw.eq.23) then

c     Compute face normals

      elseif(isw.eq.24) then

c     Compute boundary nodes

      elseif(isw.eq.26) then

c     Compute residuals and tangents for parts

      elseif(isw.eq.3 .or.  isw.eq. 4  .or. isw.eq. 6 .or.
     &        (isw.ge.8 .and. isw.le.14) .or. isw.eq.16 ) then
         
         etype = d(245)
         if(etype.eq.1) then
c     Enhanced Strain Model
            call GEDmodel(d,ul,xl,s,p,ndf,ndm,nst,isw) !updated Lagrangian
            !write(7,*) "--------- GE Updated -------"
c     
         elseif(etype.eq.2) then
            call GEDmodelTL(d,ul,xl,s,p,ndf,ndm,nst,isw) !total Lagrangian
            !write(7,*) "--------- GE Updated TOTAL-------"         
         endif
      endif

c     Formats for input-output

 2000 format(
     &     /5x,'T h r e e     D i m e n s i o n a l       S o l i d'/
     &     'F i n i t e   D e f o r m a t i o n   E l e m e n t'/)

 2010 format(5x,'Model Parameter'//,
     &     10x,'Modulus E       ',1p,1e12.5/
     &     10x,'Poisson Ratio   ',0p,1f8.5/
     &     10x,'Lambda Parameter',1p,1e12.5/
     &     10x,'Shear Modulus   ',1p,1e12.5/
     &     10x,'Dam. Threshold  ',1p,1e12.5/
     &     10x,'relaxation no.  ',10x,1p,1e12.5/1x)

 3000 format(' *ERROR* Problem control record incorrect:'/
     &     '         DOFs/Node (ndf) = ',i4,' - Should be 3 or more'/
     &     '         Nodes/Elm (nen) = ',i4,' - Should be 4 or more')

 4000 format(' *ERROR* Element type does not exist.')

      end

c     $Id:$
      subroutine GEDmodel(d,ul,xl,s,r,ndf,ndm,nst,isw)

c     3-D finite deformation GRADIENT ENHANCED NON-LOCAL DAMAGE
c     MODEL 
c     
c     ul(ndf,*) - Current nodal solution parameters
c     xl(ndm,*) - Nodal coordinates
c     
c     ____________________________________________________________
c     Modification log                                Date (mm/year)
c     Implementation based on disp30.f CL element        09/2014    
c     
c     _____________________________________________________________
      implicit  none

      integer, parameter :: mgauss= 8 ! 8 -> Number of Gauss points 2x2x2
      INTEGER, parameter :: nhv = 1 ! number of history variables
      
      integer, parameter :: maxnel= 8 ! maximum number of nodes per element
      integer, parameter :: nnode = 8 ! # of nodes per element
      INTEGER, parameter :: nnodeGE = 8   
      integer, parameter :: nhtd  = mgauss
      INTEGER, parameter :: ntens = 6 ! size of stress tensor

      real*8,  parameter :: zero = 0.d0,quarter=0.25d0, half =0.5d0
      real*8,  parameter :: one = 1.d0, two = 2.d0, three=3.d0

      include  'eldata.h'       ! nel : # of nodes on element
      include  'cdata.h'        ! nen : max nodes per element
      include  'elplot.h'       ! tt(1000)
      include  'iofile.h'       ! iow
      include  'hdata.h'        ! nh1, nh2, nh3 
      include  'comblk.h'       ! hr

      INTEGER   i,j,k,l,m,ik,jk,kk,mm,pp,qq !INDICES
      integer   ndm             ! spatial dimension
      integer   ndf             ! max dof per node
      integer   isw  
      integer   nst             ! size of element matrix ?
      integer   imat            ! material type   

      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,nst),r(nst)   
!     --------------------------------------------------  
!     QUADRATURE - HISTORY
      integer  nn , ntm
      REAL*8   pgp(ndm,mgauss)  !position GP
      REAL*8   wgp(mgauss)      !weights GP
      INTEGER  intip            !counter for integration point
      REAL*8   ksi(ndm)
!     --------------------------------------------------  
      real*8   e1, nu12 ,lambda ,g12
!----------------------------------------------------
      REAL*8   Bmat(6,ndm*nnode)  
      REAL*8   fact
      REAL*8   CmatE(3,3,3,3)  
      REAL*8   stressE(3,3)
      REAL*8   sigDVoigt(6)
      REAL*8   sigD(3,3)  
      REAL*8   CmatD(3,3,3,3)  
      REAL*8   CmatDVoigt(6,6)  
!     --------------------------------------------------  
      REAL*8   shapef(nel)
      REAL*8   dshape(nel,ndm)
!----------------------------------------------------
      REAL*8   xreff(ndm,nel)   !REFERENCE COORDINATES
      REAL*8   xcurr(ndm,nel)   !CURRENT COORDINATES
!----------------------------------------------------
!Jacobian in the current and in the reference configuration
      REAL*8   xjacR(3,3)
      REAL*8   xjacinvR(3,3)
      REAL*8   detxjacR
      
      REAL*8   xjacC(3,3)
      REAL*8   xjacinvC(3,3)
      REAL*8   detxjacC
!--------------------------------------------------
! KINEMATICS
      REAL*8   Fgrad(3,3)       !Deformation gradient
      REAL*8   Fgradinv(3,3)    !Inverse deformation gradient
      REAL*8   FgradT(3,3)      !Deformation gradient transpose
      REAL*8   detF             !determinant deformation gradient 

      REAL*8   bTensor(3,3)
      REAL*8   EGL(3,3)         !Green-Lagrange strain
      REAL*8   EAlm(3,3)        !Euler Almansi strain
      REAL*8   Ctens(3,3)
!--------------------------------------------------
!--------------------------------------------------
      REAL*8   Fd(ndm*nel)
      REAL*8   felem(ndm*nel)
      
      REAL*8   Fg(nnodeGE)
      REAL*8   fgelem(nnodeGE,1)
      
      REAL*8   sumKdd(ndm*nel,ndm*nel)
      REAL*8   Kdd(ndm*nel,ndm*nel)
      REAL*8   Kgeodd(ndm*nel,ndm*nel)
      REAL*8   Kgd(nnodeGE,ndm*nel) 
      REAL*8   Kdg(ndm*nel,nnodeGE) !60x20 for Q2Q2
      REAL*8   sumKdg(ndm*nel,nnodeGE) !60x20 for Q2Q2
      REAL*8   sumKgg(nnodeGE,nnodeGE) !20x20 for Q2Q2
      REAL*8   Kgg(nnodeGE,nnodeGE) !20x20 for Q2Q2

!--------------------------------------------------
!--------------------------------------------------
!     GRADIENT ENHANCED MAGNITUDES
!--------------------------------------------------
!--------------------------------------------------
!     PROJECTION 
      REAL*8   shp3(nel,mgauss) ! (nnode,mgauss) shape functions
      REAL*8   xsjdv(mgauss)
      REAL*8   sigforp(1,mgauss)
      REAL*8   sig8x8(8,mgauss)    
!--------------------------------------------------
!  FOR output part of Local Newton
      REAL*8   yvector(3,1)
      REAL*8   yscalar
      REAL*8   dsigdphi(3,3)    ! for global K^{def,phi} and K^{phi,def}   
      REAL*8   dyscalardphi     ! for global K^{phi,phi}  nonlocal damage uncoupled part
      REAL*8   dydgradphi(3,3)  ! for global K^{phi,phi}  nonlocal damage uncoupled part 
!--------------------------------------------------
      real*8   fdfuncp 
      real*8   fdfunc
      REAL*8   Kd               ! Threshold parameter from input file
      REAL*8   ndp 
!--------------------------------------------------
!----------------C O U P L E D DEFINITIONS---------
!--------------------------------------------------
!Compute grad_x N^(phi)
      REAL*8   shapefGE(nnodeGE)
      REAL*8   dshapeGE(nnodeGE,ndm)
      REAL*8   GradshapePhi(nnodeGE,ndm)
      REAL*8   GradshapePhiT(ndm,nnodeGE)
      REAL*8   PhiInterp
      REAL*8   pVec(nnodeGE,1)
      REAL*8   Gradphi(ndm,1)
      REAL*8   GradphiT(1,ndm)

      INTEGER  A11,A12,B11,B12  
!     --------------------------------------------------

      integer Nota(3,3)         !notation for tensor transformation

      Data Nota/1,4,6,
     *     4,2,5,
     *     6,5,3/  
!     --------------------------------------------------
     
c     Get input parameters
      e1      =  d(1)   
      nu12    =  d(2)   
      lambda  =  d(21)  
      g12     =  d(22)   
      Kd      =  d(243)         !  growth parameter  
      ndp     =  d(244) 
      
! initialization aux matrices
      Kdd(1:ndm*nnode,1:ndm*nnode) = zero
      Kgeodd(1:ndm*nnode,1:ndm*nnode) = zero
      sumKdd(:,:) = zero

      Kgd(1:nnodeGE ,1:ndm*nnode) = zero

      Kdg(1:ndm*nnode,1:nnodeGE) = zero
      sumKdg(1:ndm*nnode,1:nnodeGE) = zero
      
      Kgg(1:nnodeGE,1:nnodeGE) = zero
      sumKgg(1:nnodeGE,1:nnodeGE) = zero
c     
      Fd(1:ndm*nnode) = zero
      felem(1:ndm*nnode) = zero
      Fg(1:nnodeGE)  = zero
      fgelem(1:nnodeGE,1)  = zero
c     
      xcurr(1:ndm,1:nel)=zero
      xreff(1:ndm,1:nel)=zero 
c     
c     Initilize history variables 
      if (isw.eq.14) then
         nn = 0
         do l=1,mgauss
            hr(nh1+nn) = Kd
!     hr(nh1+nn+1) = Kd              
            nn = nn+nhv
         end do                 !l
      endif

c     Get quadrature information - EVALUATE POSITIONS AND WEIGHTS

      if(nel.eq.8) then
         call quadrature3D8(pgp,wgp)
!call kgauss20Node(pgp,wgp,ndm,nnode,mgauss)
      else
         write(*,*)'ERROR in elmt14.f - wrong quadrature rule'
         write(iow,*) 'ERROR in elmt14.f - wrong quadrature rule'
      end if

c     Compute current geometry
      do i = 1,ndm
         do j = 1,nel
            xreff(i,j) = xl(i,j)
            xcurr(i,j) = xl(i,j) + ul(i,j,1)  
         end do                 ! j
      end do     
      
      do j=1,nnodeGE            ! pVec(nnodeGE)
         pVec(j,1) = ul(4,j,1)
      enddo
      
!call mprint(pVec,8,1,8,"PGI")
c     Compute shape functions

c     -----------------------------------------------------------------
c     LOOP OVER THE INTGRATION POINTS
c     -----------------------------------------------------------------
      
      nn = 0                    ! history variable update 
      do intip = 1,mgauss
         
!write(*,*) "INTEGRATION PART AT@" , intip

         do i=1,3
            ksi(i) = pgp(i,intip) !Position for each IP
         end do

!Evaluate shape functions 
         call getShape3D8(shapef,ksi,ndm)
         
         do i = 1,nel
            shp3(i,intip) = shapef(i)  
         enddo

!Evaluate shape functions derivatives
         call getDerivativeShape3D8(dshape,ksi,ndm)
         
!Compute Jacobian: REFERENCE CONFIGURATION
         call getJacobian3D(xjacR,detxjacR,xjacinvR,
     *        dshape,ndm,nel,xreff)

!Compute Jacobian: CURRENT CONFIGURATION

         call getJacobian3D(xjacC,detxjacC,xjacinvC,
     *        dshape,ndm,nel,xcurr)

!Compute deformation gradient F = je*Je^-1
         call getDeformGrad3D(Fgrad,Fgradinv,FgradT,detF,xjacR,detxjacR,
     *        xjacinvR,xjacC,detxjacC,xjacinvC,ndm,nel)
         
!call mprint(Fgrad,3,3,3,"F")
c     !Compute Deformation tensor
c     ! EGL --> Green-Langrage, C = FtF
c     ! EAlm --> Euler-Almansi, b = FFt
c     call kDefTensors(EGL,EAlm,Ctens,btens,Fgrad,
c     *  Fgradinv,FgradT,detF)    
c     !DEFORMATION GRADIENT DECOMPOSITION F=FeFg
c     call kDeformGradGrowth(Fgrowth,Felast,Fgrad,Fgradinv,
c     *  ndm,d,detFe,FelastT,Felastinv)
c     call kDefTensors(EGLE,EAlmE,CtensE,btensE,Felast,
c     *  Felastinv,FelastT,detFe)

!Compute B-operator
         call getBmat3D(Bmat,dshape,xjacC,detxjacC,xjacinvC,ndm,
     *        nel,ndm) 
!     ----------------------------------------------------
!     ----------- N O N L O C A L     P A R T  -----------
!     ----------------------------------------------------
         call getGradshapePhi(GradshapePhi,GradshapePhiT,dshape,
     *        xjacinvC,ndm,nnodeGE)

!     ----------------------------------------------------
!--------Interpolate----------------------------------    
c     \grad\phi
c     \phi
         call getGradPhiInterp(PhiInterp,Gradphi,GradphiT,GradshapePhi,
     1        pVec,shapef,dshape,nnodeGE,ndm)  
!---- -------------------------------------------------    
! MATERIAL MODEL 
         call StandardNeoh3(d,detF,Fgrad,stressE,CmatE,bTensor,ndm)
!     ----------------------------------------------------
         fact = detxjacC*wgp(intip)
         xsjdv(intip) = fact    ! for projection, volume element at every gauss point 

         ntm = nhv*mgauss       ! MAx number of history variables

c     damage evolution
         call LocalNewtonDamage(sigD,yvector,yscalar,CmatD,
     1        dsigdphi,dyscalardphi,dydgradphi,hr(nh1+nn),hr(nh2+nn),ntm
     $        ,d, Gradphi,PhiInterp,bTensor,detF,stressE,CmatE)

!call mprint(yvector,3,1,3,"yvector")
!call mprint(yscalar,1,1,1,"yscalar")
!call mprint(dsigdphi,3,3,3,"dsigdphi")

!     ---------------------------------------------
!     ------------ SWITCH TO VOIGT NOTATION
!     ---------------------------------------------
         do kk=1,3
            do mm=1,3
               do pp=1,3
                  do qq=1,3
                     CmatDVoigt(nota(kk,mm),nota(pp,qq))=CmatD(kk,mm,pp,
     $                    qq)
                  end do
               end do
            end do
         end do 


!Arrange stress in Voigt notation
         sigDVoigt(1) =sigD(1,1) 
         sigDVoigt(2) =sigD(2,2) 
         sigDVoigt(3) =sigD(3,3) 
         sigDVoigt(4) =sigD(1,2) 
         sigDVoigt(5) =sigD(2,3) 
         sigDVoigt(6) =sigD(1,3) 

         sig8x8(1,intip) = sigDVoigt(1)
         sig8x8(2,intip) = sigDVoigt(2)
         sig8x8(3,intip) = sigDVoigt(3)
         sig8x8(4,intip) = sigDVoigt(4)
         sig8x8(5,intip) = sigDVoigt(5)
         sig8x8(6,intip) = sigDVoigt(6)
         
!write(*,*) "history at current step: ", hr(nh2+nn)
         call damageFunction(fdfunc,fdfuncp,d,hr(nh2+nn))
         
!   write(*,*) "fdfunction: ",fdfunc,"Knp1 < Kn",
!1    hr(nh2+nn)<hr(nh1+nn)
         
         sig8x8(7,intip) = 1.d0 - fdfunc
!write(*,*) n,fdfunc
         
         nn = nn+nhv                     
!     ---------------------------------------------
!     ------------end / SWITCH TO VOIGT NOTATION
!     ---------------------------------------------

c     Compute MATerial STIFFNESS matrix
         call getKuu_Fu(Kdd,felem,nnode,ndm,ndm,fact,CmatDVoigt,
     *        Bmat,sigDVoigt,ntens) 
         
c     Compute Geometrical stiffness matrix
         call getKuuGeom(Kgeodd,dshape,sigD,ndm,nnode,ntens,
     *        xjacinvC)

!call mprint(Fd,24,1,24,"Fd")

!     1st part ------------------------
         do i=1,ndm*nnode
            Fd(i) = Fd(i) - felem(i)
            do j=1,ndm*nnode
               sumKdd(i,j) = sumKdd(i,j) + Kdd(i,j) + Kgeodd(i,j)*fact
            enddo
         enddo

!     ------------------------
!     ------------------------

!     2nd part ------------------------
!     ----- Coupling matrix Kdg between the displacements and the nonlocal part

         call getKuphi(Kdg,Bmat,dsigdphi,shapef,fact,
     *        ndm,nnode,nnodeGE)   

         do ik=1,ndm*nnode
            do jk=1,nnodeGE
               sumKdg(ik,jk) = sumKdg(ik,jk) + Kdg(ik,jk)
            enddo
         enddo 
!     ------------------------
!     ------------------------

!     3rd part ------------------------ 
!     ----- Coupling matrix Kdg between the displacements and the nonlocal part
! call kCoupMatPHIU(kgd,Bmat,dsigdphi,shapef,fact,
!*  ndm,nnode,nnode)
! 
!do ik=1,nnode
!  do jk=1,ndm*nnode
!   Kephid(ik,jk) = Kephid(ik,jk) + kphid(ik,jk)
!  enddo
!enddo 
         
!-------------------------
!-------------------------

!     4th part ------------------------ 
!-------------------------
         
         call getKphiphi(Kgg,GradshapePhi,GradshapePhiT,shapef,
     *        dyscalardphi,dydgradphi,fact,nnode,ndm) 
         
         do ik=1,nnode
            do jk=1,nnode
               sumKgg(ik,jk) = sumKgg(ik,jk) + Kgg(ik,jk)
            enddo
         enddo 
         
c     -------------------------    
c     Compute residual due to phi       
c     -------------------------    
         call getFphi(fgelem,btensor,yvector,yscalar,shapef,
     1        GradshapePhi,nnodeGE,ndm,fact)
         
         do i=1,nnodeGE
            Fg(i) = Fg(i) - fgelem(i,1)
         enddo

      end do                    ! intip

      do i =1,8
         tt(i)   =  sig8x8(7,i)
      enddo 

!     open (unit = 71, file = "Kdd.txt")
!     do, i=1,ndm*nnode
!     write(71,*) ( sumKdd(i,j), j=1,ndm*nnode )
!     enddo
!     close(71) 
!     
!call mprint(Fg,8,1,8,"Fgradient")
! write(*,*) 'Fg'
! write(*,'(1f20.10 )') (Fg(i), i =1,8)
! write(*,*) 'Fg'

      do i=1,nnode
         A11 = (ndm+1)*(i-1)+1
         A12 = ndm*(i-1)+1
         
! displacement
         r(A11  ) =  Fd(A12  ) 
         r(A11+1) =  Fd(A12+1)
         r(A11+2) =  Fd(A12+2)
         
! GE dof
         r(A11+3) =  Fg(i)
      enddo

      do i=1,nnode
         do j=1,nnode
            A11 = (ndm+1)*(i-1)+1
            A12 =  ndm*(i-1)+1
            B11 = (ndm+1)*(j-1)+1
            B12 =  ndm*(j-1)+1
            
!---------
! displacement
!---------
            s(A11,B11)     = sumKdd(A12,B12)
            s(A11,B11+1)   = sumKdd(A12,B12+1)
            s(A11,B11+2)   = sumKdd(A12,B12+2)
            
            s(A11+1,B11)   = sumKdd(A12+1,B12)
            s(A11+1,B11+1) = sumKdd(A12+1,B12+1)
            s(A11+1,B11+2) = sumKdd(A12+1,B12+2)

            s(A11+2,B11)   = sumKdd(A12+2,B12)
            s(A11+2,B11+1) = sumKdd(A12+2,B12+1)
            s(A11+2,B11+2) = sumKdd(A12+2,B12+2)
            
!---------
! Gradient enhanced 
!---------
            s(A11+3,B11+3) = sumKgg(i,j)
            
!---------
! displacement/GE Kedphi
!---------
            s(A11,B11+3)   = sumKdg(A12,j)
            s(A11+1,B11+3) = sumKdg(A12+1,j)
            s(A11+2,B11+3) = sumKdg(A12+2,j)
!---------
! GE/displacement Kephid
!---------
            s(A11+3,B11)   = sumKdg(B12,i)
            s(A11+3,B11+1) = sumKdg(B12+1,i)
            s(A11+3,B11+2) = sumKdg(B12+2,i)
            
         enddo
      enddo
      
!     do k =1,27
!     tt(27+k)   = one - sigforp(1,k)   ! extract element values.
!     enddo 

      
      if (isw.eq.8) then
c     write(*,*) 'sigforp'
c     write(*,'(7f20.10 /)') (sigforp(1,i), i =1,14)
c     call gaussValueProjection(sigforp, r, s, nel, xsjdv, shp3, nip)
         call gaussValueProjection(sig8x8, shp3, xsjdv, 
     1        mgauss, nel, 8, r, s)

      endif 
! if (isw.eq.8) then

!call gaussValueProjection(sigforp, r,s, nel,xsjdv,shp3,mgauss)

!endif 
      
c     
      END
!----------------------End Updated Lagrangian

!-----------------------Begin Total Lagrangian 
      subroutine GEDmodelTL(d,ul,xl,s,r,ndf,ndm,nst,isw)

c     3-D finite deformation GRADIENT ENHANCED NON-LOCAL DAMAGE
c     MODEL 
c     
c     ul(ndf,*) - Current nodal solution parameters
c     xl(ndm,*) - Nodal coordinates
c     
c     ____________________________________________________________
c     Modification log                                Date (mm/year)
c     Implementation based on disp30.f CL element        09/2014    
c     
c     _____________________________________________________________
      implicit  none

      integer, parameter :: mgauss= 8 ! 8 -> Number of Gauss points 2x2x2
      INTEGER, parameter :: nhv = 1 ! number of history variables
      
      integer, parameter :: maxnel= 8 ! maximum number of nodes per element
      integer, parameter :: nnode = 8 ! # of nodes per element
      INTEGER, parameter :: nnodeGE = 8   
      integer, parameter :: nhtd  = mgauss
      INTEGER, parameter :: ntens = 6 ! size of stress tensor

      real*8,  parameter :: zero = 0.d0,quarter=0.25d0, half =0.5d0
      real*8,  parameter :: one = 1.d0, two = 2.d0, three=3.d0

      include  'eldata.h'       ! nel : # of nodes on element
      include  'cdata.h'        ! nen : max nodes per element
      include  'elplot.h'       ! tt(1000)
      include  'iofile.h'       ! iow
      include  'hdata.h'        ! nh1, nh2, nh3 
      include  'comblk.h'       ! hr

      INTEGER   i,j,k,l,m,ik,jk,kk,mm,pp,qq !INDICES
      integer   ndm             ! spatial dimension
      integer   ndf             ! max dof per node
      integer   isw  
      integer   nst             ! size of element matrix ?
      integer   imat            ! material type   

      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,nst),r(nst)   
!     --------------------------------------------------  
!     QUADRATURE - HISTORY
      integer  nn , ntm
      REAL*8   pgp(ndm,mgauss)  !position GP
      REAL*8   wgp(mgauss)      !weights GP
      INTEGER  intip            !counter for integration point
      REAL*8   ksi(ndm)
!     --------------------------------------------------  
      real*8   e1, nu12 ,lambda ,g12
!----------------------------------------------------
      REAL*8   Bmat(6,ndm*nnode)  
      REAL*8   fact
      REAL*8   CmatE(3,3,3,3)  
      REAL*8   stressE(3,3)
      REAL*8   sigDVoigt(6)
      REAL*8   sigD(3,3)  
      REAL*8   CmatD(3,3,3,3)  
      REAL*8   CmatDVoigt(6,6)  

      REAL*8   BmatNL(6,ndm*nnode)

!     --------------------------------------------------  
      REAL*8   shapef(nel)
      REAL*8   dshape(nel,ndm)
!----------------------------------------------------
      REAL*8   xreff(ndm,nel)   !REFERENCE COORDINATES
      REAL*8   xcurr(ndm,nel)   !CURRENT COORDINATES
!----------------------------------------------------
!     Jacobian in the current and in the reference configuration
      REAL*8   xjacR(3,3)
      REAL*8   xjacinvR(3,3)
      REAL*8   detxjacR
      
      REAL*8   xjacC(3,3)
      REAL*8   xjacinvC(3,3)
      REAL*8   detxjacC
!--------------------------------------------------
!     KINEMATICS
      REAL*8   Fgrad(3,3)       !Deformation gradient
      REAL*8   Fgradinv(3,3)    !Inverse deformation gradient
      REAL*8   FgradT(3,3)      !Deformation gradient transpose
      REAL*8   detF             !determinant deformation gradient 

      REAL*8   CTensor(3,3)
!     REAL*8   bTensor(3,3)
      REAL*8   EGL(3,3)         !Green-Lagrange strain
      REAL*8   EAlm(3,3)        !Euler Almansi strain
      REAL*8   Ctens(3,3)
!--------------------------------------------------
!--------------------------------------------------
      REAL*8   Fd(ndm*nel)
      REAL*8   felem(ndm*nel)
      
      REAL*8   Fg(nnodeGE)
      REAL*8   fgelem(nnodeGE,1)
      
      REAL*8   sumKdd(ndm*nel,ndm*nel)
      REAL*8   Kdd(ndm*nel,ndm*nel)
      REAL*8   Kgeodd(ndm*nel,ndm*nel)
      REAL*8   Kgd(nnodeGE,ndm*nel) 
      REAL*8   Kdg(ndm*nel,nnodeGE)    !60x20 for Q2Q2
      REAL*8   sumKdg(ndm*nel,nnodeGE) !60x20 for Q2Q2
      REAL*8   sumKgg(nnodeGE,nnodeGE) !20x20 for Q2Q2
      REAL*8   Kgg(nnodeGE,nnodeGE) !20x20 for Q2Q2

!--------------------------------------------------
!--------------------------------------------------
!     GRADIENT ENHANCED MAGNITUDES
!--------------------------------------------------
!--------------------------------------------------
!     PROJECTION 
      REAL*8   shp3(nel,mgauss) ! (nnode,mgauss) shape functions
      REAL*8   xsjdv(mgauss)
      REAL*8   sigforp(1,mgauss)
      REAL*8   sig8x8(8,mgauss)    
!--------------------------------------------------
!     FOR output part of Local Newton
      REAL*8   yvector(3,1)
      REAL*8   yscalar
      REAL*8   dsigdphi(3,3)    ! for global K^{def,phi} and K^{phi,def}   
      REAL*8   dyscalardphi     ! for global K^{phi,phi}  nonlocal damage uncoupled part
      REAL*8   dydgradphi(3,3)  ! for global K^{phi,phi}  nonlocal damage uncoupled part 
!--------------------------------------------------
      real*8   fdfuncp 
      real*8   fdfunc
      REAL*8   Kd               ! Threshold parameter from input file
      REAL*8   ndp 
!--------------------------------------------------
!----------------C O U P L E D DEFINITIONS---------
!--------------------------------------------------
!     Compute grad_x N^(phi)
      REAL*8   shapefGE(nnodeGE)
      REAL*8   dshapeGE(nnodeGE,ndm)
      REAL*8   GradshapePhi(nnodeGE,ndm)
      REAL*8   GradshapePhiT(ndm,nnodeGE)
      REAL*8   PhiInterp
      REAL*8   pVec(nnodeGE,1)
      REAL*8   Gradphi(ndm,1)
      REAL*8   GradphiT(1,ndm)

      INTEGER  A11,A12,B11,B12  
!     --------------------------------------------------
!--------------EAS block added by Jose 07.04.2017
      REAL*8   shapef0(nel)     !shape function at the elemtn center 
      REAL*8   dshape0(nel,ndm) !derivative of shape function at the elemtn center 
      REAL*8   ksi0(3)          ! locations at the element center

!     Jacobian in the current and in the reference configuration
      REAL*8   xjacR0(3,3)
      REAL*8   xjacinvR0(3,3)
      REAL*8   detxjacR0

      REAL*8   kWloc(6,24)      !enhanced strain matrix interpolation: local setting
      REAL*8   kWmat(6,24)      !enhanced strain matrix interpolation: global setting


!     --------------------------------------------------
!     --------------------------------------------------
!     --------------------------------------------------
  
      integer Nota(3,3)         !notation for tensor transformation

      Data Nota/1,4,6,
     *          4,2,5,
     *          6,5,3/  
!     --------------------------------------------------
      
c     Get input parameters
      e1      =  d(1)   
      nu12    =  d(2)   
      lambda  =  d(21)  
      g12     =  d(22)   
      Kd      =  d(243)         !  growth parameter  
      ndp     =  d(244) 
      
!     initialization aux matrices
      Kdd(1:ndm*nnode,1:ndm*nnode) = zero
      Kgeodd(1:ndm*nnode,1:ndm*nnode) = zero
      sumKdd(:,:) = zero

      Kgd(1:nnodeGE ,1:ndm*nnode) = zero

      Kdg(1:ndm*nnode,1:nnodeGE) = zero
      sumKdg(1:ndm*nnode,1:nnodeGE) = zero
      
      Kgg(1:nnodeGE,1:nnodeGE) = zero
      sumKgg(1:nnodeGE,1:nnodeGE) = zero
c     
      Fd(1:ndm*nnode) = zero
      felem(1:ndm*nnode) = zero
      Fg(1:nnodeGE)  = zero
      fgelem(1:nnodeGE,1)  = zero
c     
      xcurr(1:ndm,1:nel)=zero
      xreff(1:ndm,1:nel)=zero 
c     
c     Initilize history variables 
      if (isw.eq.14) then
         nn = 0
         do l=1,mgauss
            hr(nh1+nn) = Kd
!     hr(nh1+nn+1) = Kd              
            nn = nn+nhv
         end do                 !l
      endif

c     Get quadrature information - EVALUATE POSITIONS AND WEIGHTS

      if(nel.eq.8) then
         call quadrature3D8(pgp,wgp)
!     call kgauss20Node(pgp,wgp,ndm,nnode,mgauss)
      else
         write(*,*)'ERROR in elmt14.f - wrong quadrature rule'
         write(iow,*) 'ERROR in elmt14.f - wrong quadrature rule'
      end if

c     Compute current geometry
      do i = 1,ndm
         do j = 1,nel
            xreff(i,j) = xl(i,j)
            xcurr(i,j) = xl(i,j) + ul(i,j,1)  
         end do                 ! j
      end do     
      
      do j=1,nnodeGE            ! pVec(nnodeGE)
         pVec(j,1) = ul(4,j,1)
      enddo
      
!     call mprint(pVec,8,1,8,"PGI")
c     Compute shape functions

!------------------------EAS block: added by Jose
      ksi0(1) = 0.d0
      ksi0(2) = 0.d0
      ksi0(3) = 0.d0

!     Evaluate shape functions at the element center    
      call getShape3D8(shapef0,ksi0,ndm)

!     Evaluate shape functions derivatives at the element center
      call getDerivativeShape3D8(dshape0,ksi,ndm)


!     Compute Jacobian: REFERENCE CONFIGURATION at the element center
      call getJacobian3DCenter(xjacR0,detxjacR0,
     *     dshape0,ndm,nel,xreff)

!---------------------------



c     -----------------------------------------------------------------
c     LOOP OVER THE INTGRATION POINTS
c     -----------------------------------------------------------------
      
      nn = 0                    ! history variable update 
      do intip = 1,mgauss
         
!     write(*,*) "INTEGRATION PART AT@" , intip

         do i=1,3
            ksi(i) = pgp(i,intip) !Position for each IP
         end do

!     Evaluate shape functions 
         call getShape3D8(shapef,ksi,ndm)
         
         do i = 1,nel
            shp3(i,intip) = shapef(i)  
         enddo

!     Evaluate shape functions derivatives
         call getDerivativeShape3D8(dshape,ksi,ndm)
         
!     Compute Jacobian: REFERENCE CONFIGURATION
         call getJacobian3D(xjacR,detxjacR,xjacinvR,
     *        dshape,ndm,nel,xreff)

!     Compute Jacobian: CURRENT CONFIGURATION

         call getJacobian3D(xjacC,detxjacC,xjacinvC,
     *        dshape,ndm,nel,xcurr)

!     Compute deformation gradient F = je*Je^-1
         call getDeformGrad3D(Fgrad,Fgradinv,FgradT,detF,xjacR,detxjacR,
     *        xjacinvR,xjacC,detxjacC,xjacinvC,ndm,nel)
         
!     call mprint(Fgrad,3,3,3,"F")
c     !Compute Deformation tensor
c     ! EGL --> Green-Langrage, C = FtF
c     ! EAlm --> Euler-Almansi, b = FFt
c     call kDefTensors(EGL,EAlm,Ctens,btens,Fgrad,
c     *  Fgradinv,FgradT,detF)    
c     !DEFORMATION GRADIENT DECOMPOSITION F=FeFg
c     call kDeformGradGrowth(Fgrowth,Felast,Fgrad,Fgradinv,
c     *  ndm,d,detFe,FelastT,Felastinv)
c     call kDefTensors(EGLE,EAlmE,CtensE,btensE,Felast,
c     *  Felastinv,FelastT,detFe
         
!     Compute B-operator: total Lagrangian uses the reference jacobian
         call getBmat3D(Bmat,dshape,xjacR,detxjacR,xjacinvR,ndm,
     *        nel,ndm) 


!     Compute the nonlinear B-operator = Bmat + BmatNl
         call getBmat3DNL(BmatNL,Bmat,Fgrad,xjacinvR)

      
!     -----------------------------------------------------------
!---------------EAS block
!----------------------------------------------------------
         
!     Interpolation of the EAS paramaters 
         call kEASMlocalGE3d(ksi(1),ksi(3),ksi(3),kWloc)

!     transformation from parametric to local setting
         call kEASTransfGE3D(kWmat,xjacR,xjacR0,kWloc,
     *        detxjacR0,detxjacR)

!     ----------------------------------------------------
!     ----------- N O N L O C A L     P A R T  -----------
!     ----------------------------------------------------
!-----------total Lagrangian uses the reference jacobian
         call getGradshapePhi(GradshapePhi,GradshapePhiT,dshape,
     *        xjacinvR,ndm,nnodeGE)

!     ----------------------------------------------------
!--------Interpolate----------------------------------    
c     \grad\phi
c     \phi
         call getGradPhiInterp(PhiInterp,Gradphi,GradphiT,GradshapePhi,
     1        pVec,shapef,dshape,nnodeGE,ndm)  
!---- -------------------------------------------------    
!     MATERIAL MODEL 
         call StandardNeoh3reference(d,detF,Fgrad,stressE,CmatE,
     $        CTensor,ndm)
!     ----------------------------------------------------
         fact = detxjacR*wgp(intip) ! changed to reference
         xsjdv(intip) = fact    ! for projection, volume element at every gauss point 
         
         ntm = nhv*mgauss       ! MAx number of history variables

c     damage evolution
         call LocalNewtonDamageReference(sigD,yvector,yscalar,CmatD,
     1        dsigdphi,dyscalardphi,dydgradphi,hr(nh1+nn),hr(nh2+nn),ntm
     $        ,d, Gradphi,PhiInterp,CTensor,detF,stressE,CmatE)
!     call mprint(yvector,3,1,3,"yvector")
!     call mprint(yscalar,1,1,1,"yscalar")
!     call mprint(dsigdphi,3,3,3,"dsigdphi")

!     ---------------------------------------------
!     ------------ SWITCH TO VOIGT NOTATION
!     ---------------------------------------------
         do kk=1,3
            do mm=1,3
               do pp=1,3
                  do qq=1,3
                 CmatDVoigt(nota(kk,mm),nota(pp,qq))=CmatD(kk,mm,pp,qq)
                  end do
               end do
            end do
         end do 


!     Arrange stress in Voigt notation
         sigDVoigt(1) =sigD(1,1) 
         sigDVoigt(2) =sigD(2,2) 
         sigDVoigt(3) =sigD(3,3) 
         sigDVoigt(4) =sigD(1,2) 
         sigDVoigt(5) =sigD(2,3) 
         sigDVoigt(6) =sigD(1,3) 

         sig8x8(1,intip) = sigDVoigt(1)
         sig8x8(2,intip) = sigDVoigt(2)
         sig8x8(3,intip) = sigDVoigt(3)
         sig8x8(4,intip) = sigDVoigt(4)
         sig8x8(5,intip) = sigDVoigt(5)
         sig8x8(6,intip) = sigDVoigt(6)
         
!     write(*,*) "history at current step: ", hr(nh2+nn)
         call damageFunction(fdfunc,fdfuncp,d,hr(nh2+nn))
         
!     write(*,*) "fdfunction: ",fdfunc,"Knp1 < Kn",
!     1    hr(nh2+nn)<hr(nh1+nn)
         
         sig8x8(7,intip) = 1.d0 - fdfunc
!     write(*,*) n,fdfunc
         
         nn = nn+nhv                     
!     ---------------------------------------------
!     ------------end / SWITCH TO VOIGT NOTATION
!     ---------------------------------------------

c     Compute MATerial STIFFNESS matrix: total Lagrangian uses the reference jacobian
         call getKuu_Fu(Kdd,felem,nnode,ndm,ndm,fact,CmatDVoigt,
     *        Bmat,sigDVoigt,ntens) 
         
c     Compute Geometrical stiffness matrix: total Lagrangian uses the reference jacobian
         call getKuuGeom(Kgeodd,dshape,sigD,ndm,nnode,ntens,
     *        xjacinvR)

!     call mprint(Fd,24,1,24,"Fd")

!     1st part ------------------------
         do i=1,ndm*nnode
            Fd(i) = Fd(i) - felem(i)
            do j=1,ndm*nnode
               sumKdd(i,j) = sumKdd(i,j) + Kdd(i,j) + Kgeodd(i,j)*fact
            enddo
         enddo

!     ------------------------
!     ------------------------

!     2nd part ------------------------
!     ----- Coupling matrix Kdg between the displacements and the nonlocal part

         call getKuphi(Kdg,Bmat,dsigdphi,shapef,fact,
     *        ndm,nnode,nnodeGE)   

         do ik=1,ndm*nnode
            do jk=1,nnodeGE
               sumKdg(ik,jk) = sumKdg(ik,jk) + Kdg(ik,jk)
            enddo
         enddo 
!     ------------------------
!     ------------------------

!     3rd part ------------------------ 
!     ----- Coupling matrix Kdg between the displacements and the nonlocal part
!     call kCoupMatPHIU(kgd,Bmat,dsigdphi,shapef,fact,
!     *  ndm,nnode,nnode)
!     
!     do ik=1,nnode
!     do jk=1,ndm*nnode
!     Kephid(ik,jk) = Kephid(ik,jk) + kphid(ik,jk)
!     enddo
!     enddo 
         
!-------------------------
!-------------------------

!     4th part ------------------------ 
!-------------------------
         
         call getKphiphi(Kgg,GradshapePhi,GradshapePhiT,shapef,
     *        dyscalardphi,dydgradphi,fact,nnode,ndm) 
         
         do ik=1,nnode
            do jk=1,nnode
               sumKgg(ik,jk) = sumKgg(ik,jk) + Kgg(ik,jk)
            enddo
         enddo 
         
c     -------------------------    
c     Compute residual due to phi       
c     -------------------------    
         call getFphiReference(fgelem,yvector,yscalar,shapef,
     1        GradshapePhi,nnodeGE,ndm,fact)
         
         do i=1,nnodeGE
            Fg(i) = Fg(i) - fgelem(i,1)
         enddo

      end do                    ! intip

      do i =1,8
         tt(i)   =  sig8x8(7,i)
      enddo 

!     open (unit = 71, file = "Kdd.txt")
!     do, i=1,ndm*nnode
!     write(71,*) ( sumKdd(i,j), j=1,ndm*nnode )
!     enddo
!     close(71) 
!     
!     call mprint(Fg,8,1,8,"Fgradient")
!     write(*,*) 'Fg'
!     write(*,'(1f20.10 )') (Fg(i), i =1,8)
!     write(*,*) 'Fg'

      do i=1,nnode
         A11 = (ndm+1)*(i-1)+1
         A12 = ndm*(i-1)+1
         
!     displacement
         r(A11  ) =  Fd(A12  ) 
         r(A11+1) =  Fd(A12+1)
         r(A11+2) =  Fd(A12+2)
         
!     GE dof
         r(A11+3) =  Fg(i)
      enddo

      do i=1,nnode
         do j=1,nnode
            A11 = (ndm+1)*(i-1)+1
            A12 =  ndm*(i-1)+1
            B11 = (ndm+1)*(j-1)+1
            B12 =  ndm*(j-1)+1
            
!---------
!     displacement
!---------
            s(A11,B11)     = sumKdd(A12,B12)
            s(A11,B11+1)   = sumKdd(A12,B12+1)
            s(A11,B11+2)   = sumKdd(A12,B12+2)
            
            s(A11+1,B11)   = sumKdd(A12+1,B12)
            s(A11+1,B11+1) = sumKdd(A12+1,B12+1)
            s(A11+1,B11+2) = sumKdd(A12+1,B12+2)

            s(A11+2,B11)   = sumKdd(A12+2,B12)
            s(A11+2,B11+1) = sumKdd(A12+2,B12+1)
            s(A11+2,B11+2) = sumKdd(A12+2,B12+2)
            
!---------
!     Gradient enhanced 
!---------
            s(A11+3,B11+3) = sumKgg(i,j)
            
!---------
!     displacement/GE Kedphi
!---------
            s(A11,B11+3)   = sumKdg(A12,j)
            s(A11+1,B11+3) = sumKdg(A12+1,j)
            s(A11+2,B11+3) = sumKdg(A12+2,j)
!---------
!     GE/displacement Kephid
!---------
            s(A11+3,B11)   = sumKdg(B12,i)
            s(A11+3,B11+1) = sumKdg(B12+1,i)
            s(A11+3,B11+2) = sumKdg(B12+2,i)
            
         enddo
      enddo
      
!     do k =1,27
!     tt(27+k)   = one - sigforp(1,k)   ! extract element values.
!     enddo 

      if (isw.eq.8) then
c     write(*,*) 'sigforp'
c     write(*,'(7f20.10 /)') (sigforp(1,i), i =1,14)
c     call gaussValueProjection(sigforp, r, s, nel, xsjdv, shp3, nip)
         call gaussValueProjection(sig8x8, shp3, xsjdv, 
     1        mgauss, nel, 8, r, s)

      endif 
!     if (isw.eq.8) then

!     call gaussValueProjection(sigforp, r,s, nel,xsjdv,shp3,mgauss)

!     endif
            
c     
      END

!--------------------End total Lagrangian

      
      subroutine StandardNeoh3(d,detf,fgrad,stressE,CmatE,btensor,ndim)
c     * * F E A P * * A Finite Element Analysis Program

c.... Copyright (c) 1984-2008: Regents of the University of California
c     All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Neo-Hookean model

c     Inputs:
c     d(100)     Material parameters
c     d(21)      Bulk  modulus
c     d(22)      Shear modulus
c     detf       Jacobian determinant at t_n+1
c     beVoigt(6)      Left Cauchy-Green tensor

c     Outputs:
c     sig(6)     CAUCHY stress tensor
c     aa(6,6)    CAUCHY (spatial) elastic moduli
c     estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0

      integer   i,jsw,j,k,l,ndim
      real*8    detf, press, ha,hp,hpp, estore
      real*8    xlamd,logj
      real*8    d(*)
      REAL*8    btensor(3,3)  
      REAL*8    CmatE(3,3,3,3)
      REAL*8    stressE(3,3)
      
      REAL*8    g12, lambda
      real*8    detfe
      real*8    fgrad(3,3)
      real*8    fgradT(3,3)

      REAL*8 Id2(3,3)
      REAL*8 II4(3,3,3,3)

!     Definition second order identity
      do i=1,3
         do j=1,3
            Id2(i,j) = zero
            btensor(i,j) = zero
         enddo
         Id2(i,i) = one     
      enddo

c     MATERIAL PARAMETERS

      lambda = d(21)
      g12    = d(22)

c     Compute Transpose of F
      call getTranspose3D(fgrad,3,3,fgradT)  
!     kCompTranspose20Node(A,m,n,AT) 
      
      call getMatrixMultiply(btensor,ndim,ndim,
     *     fgrad,ndim,ndim,fgradT,ndim,ndim)
      
      do i= 1,ndim
         do j= 1,ndim
            stressE(i,j) = g12*btensor(i,j) / detF + 
     1           (lambda * dlog(detF)/detF - g12 / detF)*Id2(i,j)
         enddo
      enddo

c     COMPUTE tangent

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  CmatE(i,j,k,l)=zero
                  II4(i,j,k,l)=(half)*(Id2(i,k)*Id2(j,l) +
     $                 Id2(i,l)*Id2(j,k))
               enddo
            enddo
         enddo
      enddo   

!     Compute material tangent
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  CmatE(i,j,k,l)= lambda*Id2(i,j)*Id2(k,l)/detF+
     1                 two*(g12 - lambda*dlog(detF))*II4(i,j,k,l)/detF
               enddo
            enddo
         enddo
      enddo  


      RETURN
      END

c     ---------------------------------------------------------         

      subroutine StandardNeoh3reference(d,detf,fgrad,stressE,CmatE,
     $     Ctensor, ndim)
c     * * F E A P * * A Finite Element Analysis Program

c.... Copyright (c) 1984-2008: Regents of the University of California
c     All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Finite deformation elasticity Neo-Hookean model

c     Inputs:
c     d(100)     Material parameters
c     d(21)      Bulk  modulus
c     d(22)      Shear modulus
c     detf       Jacobian determinant at t_n+1
c     beVoigt(6)      Left Cauchy-Green tensor

c     Outputs:
c     sig(6)     CAUCHY stress tensor
c     aa(6,6)    CAUCHY (spatial) elastic moduli
c     estore     Stored energy density
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0

      integer   i,jsw,j,k,l,ndim
      real*8    detf, press, ha,hp,hpp, estore
      real*8    xlamd,logj
      real*8    d(*)
      REAL*8    Ctensor(3,3)  
      REAL*8    CtensorInv(3,3)  
      REAL*8    CmatE(3,3,3,3)
      REAL*8    stressE(3,3)
      
      REAL*8    g12, lambda
      real*8    detC
      real*8    fgrad(3,3)
      real*8    fgradT(3,3)

      REAL*8 Id2(3,3)
      REAL*8 II4(3,3,3,3)

!     Definition second order identity
      do i=1,3
         do j=1,3
            Id2(i,j) = zero
            Ctensor(i,j) = zero
         enddo
         Id2(i,i) = one     
      enddo

c     MATERIAL PARAMETERS

      lambda = d(21)
      g12    = d(22)

c     Compute Transpose of F
      call getTranspose3D(fgrad,3,3,fgradT)  
!     kCompTranspose20Node(A,m,n,AT) 
      
      call getMatrixMultiply(Ctensor,ndim,ndim,
     *     fgradT,ndim,ndim,fgrad,ndim,ndim)

      detC = detF*detF
      call getInverse3D(CtensorInv,Ctensor,detC)

c$$$     call mprint(Ctensor,3,3,3,"C")
c$$$     call mprint(CtensorInv,3,3,3,"Cinv")

      do i= 1,ndim
         do j= 1,ndim
            stressE(i,j) = g12*id2(i,j) + 
     1           (lambda*dlog(detF) - g12)*CtensorInv(i,j)
         enddo
      enddo

c     COMPUTE tangent

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  CmatE(i,j,k,l)=zero
                  II4(i,j,k,l)=half*(CtensorInv(i,k)*CtensorInv(j,l) +
     $                 CtensorInv(i,l)*CtensorInv(j,k))
               enddo
            enddo
         enddo
      enddo   

!     Compute material tangent
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  CmatE(i,j,k,l)=lambda*CtensorInv(i,j)*CtensorInv(k,l)+
     1                 two*(g12 - lambda*dlog(detF))*II4(i,j,k,l)
               enddo
            enddo
         enddo
      enddo  


      RETURN
      END

c     ---------------------------------------------------------         
c     ---------------------------------------------------------  

      subroutine getKuu_Fu(stiff,force,nnode,ndim,ndof,factor,
     *     Dmat,Bmat,stress,ntens)
      
      IMPLICIT NONE
      
c     number definition
      REAL*8,PARAMETER :: zero = 0.d0
      REAL*8,PARAMETER :: one = 1.d0
      REAL*8,PARAMETER :: two = 2.d0
      REAL*8,PARAMETER :: three = 3.d0    
      REAL*8,PARAMETER :: half = 0.5d0 
      REAL*8,PARAMETER :: eight = 8.0d0
      
      INTEGER nnode,ndim,ndof,ntens
      
      REAL*8 stiff(ndim*nnode,ndim*nnode)
      REAL*8 force(ndim*nnode)
      REAL*8 Bmat(ntens,ndim*nnode)
      REAL*8 Dmat(ntens,ntens)
      REAL*8 stress(ntens)
      
      REAL*8 factor             !Jacobian
      
      REAL*8 db(ntens)
      REAL*8 sum
      integer ndofnnode,ndimnnode 

      INTEGER iforce,jforce
      REAL*8 sumforce
      
      INTEGER istiff,jstiff
      INTEGER i,j,k,l,m,n
      
      ndofnnode= ndof*nnode
      ndimnnode= ndim*nnode

!     Initialization
      do istiff=1,ndofnnode
         force(istiff)=zero
         do jstiff=1,ndofnnode
            stiff(istiff,jstiff)=zero
         end do
      end do

!     Computing the stiffness matrix in a condensed way

      do j=1,ndimnnode
         do k=1,ntens
            db(k) = zero        !initialization to zero
            do l=1,ntens
               db(k) = db(k) + Dmat(k,l)*Bmat(l,j)*factor
            end do
         end do
         do i=1,ndimnnode
            sum = zero
            do m=1,ntens
               sum = sum + bmat(m,i)*db(m)
            end do
            stiff(i,j) = stiff(i,j) + sum
         end do
      end do  
      
!     Computing the internal force vector
      
      do jforce =1,ndimnnode
         sumforce= zero
         do iforce = 1,ntens
            sumforce = sumforce + bmat(iforce,jforce)*stress(iforce)*
     $           factor
         end do
         force(jforce) = sumforce
      end do
      
      END

!------------------------------------------
!------------------------------------------------------

!---------------------------------------------------------
!---------------------------------------------------------


      subroutine getKuuGeom(Kgeodd,dshape,stress,ndim,nnode,
     *     ntens,xjacinv)
      
      Implicit none
      
      
c     number definition
      REAL*8,PARAMETER :: zero = 0.d0
      REAL*8,PARAMETER :: one = 1.d0
      REAL*8,PARAMETER :: two = 2.d0
      REAL*8,PARAMETER :: three = 3.d0    
      REAL*8,PARAMETER :: half = 0.5d0 
      REAL*8,PARAMETER :: eight = 8.0d0
      
      
      REAL*8 dshape(nnode,ndim)
      REAL*8 dshapeT(ndim,nnode)
      REAL*8 dshapeC(ndim,nnode)
      
      REAL*8 GRAD(ndim)
      REAL*8 Kgeodd(ndim*nnode,ndim*nnode)
      REAL*8 stress(ndim,ndim)  !Kirchoff stress in the current configuration
      
      REAL*8 xjacinv(ndim,ndim)
      integer ndofnnode,ndimnnode 
      
      REAL*8 scalar
      
      INTEGER ndim,nnode,ntens
      INTEGER i,j,k,l,kk,ll

      ndimnnode= ndim*nnode
      
!     initialization
      do i=1,ndimnnode
         do j=1,ndimnnode
            Kgeodd(i,j) = zero    
         enddo
      enddo    
      
!     initialization
      do i=1,ndim
         do j=1,nnode
            dshapeC(i,j) = zero
         enddo
      enddo    

!     Pay attention: dshape is stored in this element nnode*ndim (20*3), but 
!     this matrix needs to be trasnposed in order to be consistent with the Geometrical compputation
!     Transposing the shape derivative matrix --> shapeT
      call getTranspose3D(dshape,ndim,nnode,dshapeT) 
      
!     Compute derivatives wrt the current config. according to Wriggers definition 4.93 page 134
      do i=1,nnode
         dshapeC(1,i)=dshapeT(1,i)*xjacinv(1,1)+
     *        dshapeT(2,i)*xjacinv(2,1)+
     *        dshapeT(3,i)*xjacinv(3,1)
         
         dshapeC(2,i)=dshapeT(1,i)*xjacinv(1,2)+
     *        dshapeT(2,i)*xjacinv(2,2)+
     *        dshapeT(3,i)*xjacinv(3,2)
         
         dshapeC(3,i)=dshapeT(1,i)*xjacinv(1,3)+
     *        dshapeT(2,i)*xjacinv(2,3)+
     *        dshapeT(3,i)*xjacinv(3,3)
      enddo
      
      kk = 1   
      do k=1,nnode
         ll = 1   
         do l=1,nnode
            scalar = zero
            
            do i=1,ndim
               do j=1,ndim
                  scalar=scalar + dshapeC(i,k)*stress(i,j)*dshapeC(j,l)
               enddo
            enddo 
            
            Kgeodd(kk,ll) = scalar
            Kgeodd(kk+1,ll+1) = scalar
            Kgeodd(kk+2,ll+2) = scalar
            
            ll = ll + 3  
         enddo                  ! l
         
         kk = kk + 3     
      enddo                     ! k   
      
      end

c---------------------------------------------------------

      subroutine LocalNewtonDamage(sigd,yvector,yscalar,CmatD,
     1     dsigdphi,dyscalardphi,dydgradphi,hn,hn1,ntm,d,Gradphi,
     $     phinloc, bTensor,detF,stressE,CmatE)

      IMPLICIT NONE

      include 'iofile.h'
      include 'hdata.h'

      REAL*8,PARAMETER :: zero = 0.d0
      REAL*8,PARAMETER :: one = 1.d0
      REAL*8,PARAMETER :: two = 2.d0
      REAL*8,PARAMETER :: three = 3.d0    
      REAL*8,PARAMETER :: half = 0.5d0 
      REAL*8,PARAMETER :: eight = 8.0d0

!     History 
      real*8 hn(*), hn1(*)
      integer ntm               ! max number of history variables
!     input
      REAL*8 detF      
      REAL*8 stressE(3,3)       
      REAL*8 CmatE(3,3,3,3)       
      REAL*8 CmatD(3,3,3,3)       
      REAL*8 d(*)
      real*8 beVoigt(6)
      real*8 bTensor(3,3)
      real*8 detFe
      REAL*8 CmatVoigt(6,6)

!     LOCAL NEWTON ITERATION
      REAL*8 Knp1, Kn
      REAL*8 ndp                ! relaxation parameter of exponential damage function
      REAL*8 gammadp            ! non local switch 
      REAL*8 betadp             ! penalty parameter 
      REAL*8 trb
      REAL*8 phinloc
      real*8 fdfunc             ! value of - exponential- function
      real*8 fdfuncp            ! derivative of exponential-damage function
      real*8 q                  ! conjugate force to nonlocal damage

      REAL*8 fdfunctol
      REAL*8 chk                ! damage surface 
      integer damped            ! damped part 
      real*8 tol                ! local newton TOLERANCE
      real*8 ctol               ! local newton TOLERANCE
      REAL*8 gammak
      REAL*8 gamma1       
      REAL*8 residual           ! residula for local newton iteration
      REAL*8 dres               ! slope for local newton iteration

      REAL*8 dgamma             ! delta K = lambda * dt ----> increment 
      REAL*8 eta                ! damping constant
      REAL*8 psiE
!     material PARAMETERS
      real*8  e1, nu12 ,lambda ,g12,grthpar
      REAL*8  Kd                ! Threshold parameter from input file
      REAL*8 cdp                ! enhanced free energy part scalar parameter 
      REAL*8 thetadam           ! tangent moduli required parameter
      integer i,j,k,l,m         ! counter
      integer kk,mm,pp,qq 
!--------------------------------------------------
      integer Nota(3,3)         !notation for tensor transformation
      REAL*8 kc, kc2
      REAL*8 edamVoigt(6,6)
      REAL*8 edam(3,3,3,3)
      REAL*8 Id2(3,3)
      LOGICAL damageflg
!--------------------------------------------------
      REAL*8 Gradphi(3,1)
c     FOR output part
      REAL*8 sigd(3,3)  
      REAL*8 yvector(3,1)
      real*8 yscalar
!     O U T P U T 
      REAL*8 dsigdphi(3,3)      ! for global K^{def,phi} and K^{phi,def}   
      REAL*8 dyscalardphi       ! for global K^{phi,phi}  nonlocal damage uncoupled part
      REAL*8 dydgradphi(3,3)    ! for global K^{phi,phi}  nonlocal damage uncoupled part 
      
!--------------------------------------------------
      integer count

      Data Nota/1,4,6,
     *     4,2,5,
     *     6,5,3/  
      
c     !initialization
      count = 0

c     write(*,*) "phinloc" , phinloc
      
      do i=1,6
         do j=1,6
            edamVoigt(i,j) = zero
         enddo
      enddo

      do i=1,3
         do j=1,3
            sigD(i,j)= zero  
            Id2(i,j) = zero
         enddo
         Id2(i,i) = one
      enddo
c     ! end initilization 

c     ! get MAterial Parameters

      e1      =  d(1)   
      nu12    =  d(2)   
      lambda  =  d(21)  
      g12     =  d(22)   
      Kd      =  d(243)         !  growth parameter  
      ndp     =  d(244)

c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------        
c     (for each integration point)  Local Newton Iteration for local damage parameter update.
c     !!! WITHOUT DAMPING PART
c     write(*,*) "Kn1", Kn  
!     0.  given F, phi, Kd, gammad, nd, betad, 

      Kn = hn(1)                ! assign the 
c     Gradphi
c     phinloc
      
!     if (Kn.lt.Kd ) then
!     Kn = Kd
!     endif

      cdp = 10.0d0
      gammadp = one
      betadp = 1000.0d0
      
!     1.  calculate free energy psiE, use cauchy stress stressE, and elasticity tensor CmatVoigt
      trb = bTensor(1,1) + bTensor(2,2) + bTensor(3,3) 

      psiE = g12*half*(trb - three)-g12*dlog(detF) +
     &     lambda*half*dlog(detF)*dlog(detF)
      
!     3.COMPUTE q
      call damageFunction(fdfunc,fdfuncp,d,Kn)
      q = psiE + (betadp*gammadp)/(ndp*fdfunc)*(phinloc - gammadp*Kn)

!     check damage surface
      chk = q - Kn
      
!     write(*,'(3f20.10 )') chk,q,Kn

      fdfunctol = 0.05d0
      if ((chk.le.zero) .or. (fdfunc.lt.fdfunctol)) then
c     write(*,*) 'N O   D A M A G E!!!'
         damageflg = .false.
         Knp1 = Kn 
      else
         damageflg = .true.
!     write(*,*) ("DAMAGE !!!!! !")
c     damped = 0  gamma1  = zero  gammak = zero

         Knp1 = Kn              ! set Kn+1 =Kn
         tol = 1.0000E-7
         ctol = 50
         residual=1.0d0

         do while (abs(residual).gt.tol)
            count =count +1
c     (5a) calculate res
            
            if (count.gt.ctol) then
               write(*,*) count, "NO LOCAL CONVERGENCE!"
               write(*,*) residual
               stop
            endif
            
            call damageFunction(fdfunc,fdfuncp,d,Knp1)
            residual = (betadp*gammadp)/(ndp*fdfunc)
            residual = psiE + residual*(phinloc - gammadp*Knp1) - Knp1

c     calculate local tangent
            dres = (betadp*gammadp)/(ndp*fdfunc)
            dres = dres*(ndp*(phinloc - gammadp*Knp1)-gammadp) - one

!     write(*,*) "damage",count,residual
c     compute the increment 
            gamma1 = -residual/dres
c     !update K
            Knp1 = Knp1 + gamma1;
         end do                 ! do while loop 

      endif
      
      call damageFunction(fdfunc,fdfuncp,d,Knp1)
!     if (n.eq.2) then
!     write(*,*) n, fdfunc
!     endif

!     write(*,*) "Damage:" ,damageflg,"Kn:" ,Kn, "  Knp1",
!     1   Knp1,"Kd",Kd, "phinloc", phinloc   

      if (fdfunc.lt.0.01) then
         write(*,*) "_WARNING!!_",chk, Kn, Knp1, residual
      endif

      if (Kn.gt.Knp1) then
         write(*,*) "_________ HERE ____________"
         Knp1 = Kn 
         call damageFunction(fdfunc,fdfuncp,d,Knp1)
      endif 
!     Update history

      hn1(1) = Knp1   

c     step6 -> Compute flux/source terms
!     a stress ~ updated fdfunc  
!|||| 
!|||  FOR FORCE MATRIX
!|||| 
      do i =1,3
         yvector(i,1) = cdp/detF*Gradphi(i,1) 
         do j=1,3
            sigD(i,j)    = (fdfunc+0.0001)*stressE(i,j)
         enddo
      enddo
      yscalar = betadp*(phinloc - gammadp* Knp1)/detF 
      
!|||| 
!|||  FOR STIFFNESS MATRIX
!|||| 
c     step 7 -> Compute tangent moduli
      if (damageflg) then
!     a damage tangent
!     thetadam = gammadp*(one+ndp*Knp1)-ndp*phinloc
!     thetadam = betadp*gammadp*thetadam 
!     thetadam = - one - thetadam/(ndp*fdfunc)
!     thetadam = one/thetadam
         
         thetadam = -ndp*fdfunc - gammadp*gammadp*betadp + 
     1        gammadp*betadp*ndp*(phinloc - gammadp*Knp1)
         
         thetadam = ndp*fdfunc/thetadam

      else 
         thetadam = zero
      endif
c     write(iow,*) "thetadam",thetadam

      kc = ndp*fdfunc*thetadam*detF ! damage part tangent scaling multiplier

!     write(*,*) "kc",kc
      
      do i=1,3
         do j=1,3
            do k=1,3
               do m=1,3
                  edam(i,j,k,m) = kc*stressE(i,j)*stressE(k,m)          
               enddo
            enddo
         enddo
      enddo  


      do kk=1,3
         do mm=1,3
            do pp=1,3
               do qq=1,3
                  CmatD(kk,mm,pp,qq)=(fdfunc+0.0001)*CmatE(kk,mm,pp,qq)+ 
     &                 edam(kk,mm,pp,qq) 
               end do
            end do
         end do
      end do
!     b Second (middle) quantity needed for the couple part of global tangent 
      kc2 = betadp*gammadp*thetadam
      do i=1,3
         do j=1,3
            dsigdphi(i,j) = kc2*stressE(i,j) ! 2 dy/dg is same! leading to the symmetry of the tangent at the end. 
         enddo
      enddo
      
!     c Fourth (middle) quantity needed for the couple part of global tangent 
      dyscalardphi = betadp*gammadp*gammadp*thetadam  
      dyscalardphi = one + dyscalardphi/(ndp*fdfunc)   
      dyscalardphi = betadp*dyscalardphi/detF

c     write(*,*)  "chk" , chk,  "phiInterp",phinloc, "Kn+1", Knp1,
c     1  "thetadam", thetadam, "fdfunc", fdfunc, "ndp", ndp,
c     2  "detF", detF

!     d Fifth (middle) quantity needed for the couple part of global tangent 
      do i=1,3
         do j=1,3
            dydgradphi(i,j)= (cdp/detF)*bTensor(i,j) 
         enddo
      enddo
      
c     write(iow,*)  "tehtadam",thetadam, "dsigdphi"
c     write(iow,*)  "dsigdphi"
c     write(iow,'(3f20.10 /)')  ((dsigdphi(i,j),i=1,3),j=1,3)
c     write(iow,*)  "dydgradphi"
c     write(iow,'(3f20.10 /)')  ((dydgradphi(i,j),i=1,3),j=1,3)
c     write(iow,*)  "dyscalardphi",dyscalardphi*detf
c     c
c     
c     write(iow,*) "Knp1",Knp1, "phinloc", phinloc, "chk",chk

c     write(*,*) "hn(1 )", hn(1 )
c     write(*,*) "hn1(1)", hn1(1)
      

c     write(*,'(6f10.4 /)')  ((dsigdphi(i,j),i=1,3),j=1,3)

      return
      END
c     ----------------------------------------------------------------------------

      subroutine LocalNewtonDamageReference(sigd,Yvector,yscalar,CmatD,
     1     dsigdphi,dyscalardphi,dydgradphi,hn,hn1,ntm,d,Gradphi,
     $     phinloc, CTensor,detF,stressE,CmatE)

      IMPLICIT NONE

      include 'iofile.h'
      include 'hdata.h'

      REAL*8,PARAMETER :: zero = 0.d0
      REAL*8,PARAMETER :: one = 1.d0
      REAL*8,PARAMETER :: two = 2.d0
      REAL*8,PARAMETER :: three = 3.d0    
      REAL*8,PARAMETER :: half = 0.5d0 
      REAL*8,PARAMETER :: eight = 8.0d0

!     History 
      real*8 hn(*), hn1(*)
      integer ntm               ! max number of history variables
!     input
      REAL*8 detF      
      REAL*8 stressE(3,3)       
      REAL*8 CmatE(3,3,3,3)       
      REAL*8 CmatD(3,3,3,3)       
      REAL*8 d(*)
      real*8 beVoigt(6)
      real*8 CTensor(3,3)
      real*8 detFe
      REAL*8 CmatVoigt(6,6)

!     LOCAL NEWTON ITERATION
      REAL*8 Knp1, Kn
      REAL*8 ndp                ! relaxation parameter of exponential damage function
      REAL*8 gammadp            ! non local switch 
      REAL*8 betadp             ! penalty parameter 
      REAL*8 trb
      REAL*8 phinloc
      real*8 fdfunc             ! value of - exponential- function
      real*8 fdfuncp            ! derivative of exponential-damage function
      real*8 q                  ! conjugate force to nonlocal damage

      REAL*8 fdfunctol
      REAL*8 chk                ! damage surface 
      integer damped            ! damped part 
      real*8 tol                ! local newton TOLERANCE
      real*8 ctol               ! local newton TOLERANCE
      REAL*8 gammak
      REAL*8 gamma1       
      REAL*8 residual           ! residula for local newton iteration
      REAL*8 dres               ! slope for local newton iteration

      REAL*8 dgamma             ! delta K = lambda * dt ----> increment 
      REAL*8 eta                ! damping constant
      REAL*8 psiE
!     material PARAMETERS
      real*8  e1, nu12 ,lambda ,g12,grthpar
      REAL*8  Kd                ! Threshold parameter from input file
      REAL*8 cdp                ! enhanced free energy part scalar parameter 
      REAL*8 thetadam           ! tangent moduli required parameter
      integer i,j,k,l,m         ! counter
      integer kk,mm,pp,qq 
!--------------------------------------------------
      integer Nota(3,3)         !notation for tensor transformation
      REAL*8 kc, kc2
      REAL*8 edamVoigt(6,6)
      REAL*8 edam(3,3,3,3)
      REAL*8 Id2(3,3)
      LOGICAL damageflg
!--------------------------------------------------
      REAL*8 Gradphi(3,1)
c     FOR output part
      REAL*8 sigd(3,3)  
      REAL*8 yvector(3,1)
      real*8 yscalar
!     O U T P U T 
      REAL*8 dsigdphi(3,3)      ! for global K^{def,phi} and K^{phi,def}   
      REAL*8 dyscalardphi       ! for global K^{phi,phi}  nonlocal damage uncoupled part
      REAL*8 dydgradphi(3,3)    ! for global K^{phi,phi}  nonlocal damage uncoupled part 
      
!--------------------------------------------------
      integer count

      Data Nota/1,4,6,
     *     4,2,5,
     *     6,5,3/  
      
c     !initialization
      count = 0

c     write(*,*) "phinloc" , phinloc
      
      do i=1,6
         do j=1,6
            edamVoigt(i,j) = zero
         enddo
      enddo

      do i=1,3
         do j=1,3
            sigD(i,j)= zero  
            Id2(i,j) = zero
         enddo
         Id2(i,i) = one
      enddo
c     ! end initilization 

c     ! get MAterial Parameters

      e1      =  d(1)   
      nu12    =  d(2)   
      lambda  =  d(21)  
      g12     =  d(22)   
      Kd      =  d(243)         !  growth parameter  
      ndp     =  d(244)

c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------        
c     (for each integration point)  Local Newton Iteration for local damage parameter update.
c     !!! WITHOUT DAMPING PART
c     write(*,*) "Kn1", Kn  
!     0.  given F, phi, Kd, gammad, nd, betad, 

      Kn = hn(1)                ! assign the 
c     Gradphi
c     phinloc
      
!     if (Kn.lt.Kd ) then
!     Kn = Kd
!     endif

      cdp =  10.d0 !50.0d0
      gammadp = one
      betadp =  1000.0d0
      
!     1.  calculate free energy psiE, use cauchy stress stressE, and elasticity tensor CmatVoigt
      trb = CTensor(1,1) + CTensor(2,2) + CTensor(3,3) 

      psiE = g12*half*(trb - three)-g12*dlog(detF) +
     &     lambda*half*dlog(detF)*dlog(detF)
      
!     3.COMPUTE q
      call damageFunction(fdfunc,fdfuncp,d,Kn)
      q = psiE + (betadp*gammadp)/(ndp*fdfunc)*(phinloc - gammadp*Kn)

!     check damage surface
      chk = q - Kn
      
!     write(*,'(3f20.10 )') chk,q,Kn

      fdfunctol = 0.05d0
      if ((chk.le.zero) .or. (fdfunc.lt.fdfunctol)) then
c     write(*,*) 'N O   D A M A G E!!!'
         damageflg = .false.
         Knp1 = Kn 
      else
         damageflg = .true.
!     write(*,*) ("DAMAGE !!!!! !")
c     damped = 0  gamma1  = zero  gammak = zero

         Knp1 = Kn              ! set Kn+1 =Kn
         tol = 1.0000E-7
         ctol = 50
         residual=1.0d0

         do while (abs(residual).gt.tol)
            count =count +1
c     (5a) calculate res
            
            if (count.gt.ctol) then
               write(*,*) count, "NO LOCAL CONVERGENCE!"
               write(*,*) residual
               stop
            endif
            
            call damageFunction(fdfunc,fdfuncp,d,Knp1)
            residual = (betadp*gammadp)/(ndp*fdfunc)
            residual = psiE + residual*(phinloc - gammadp*Knp1) - Knp1

c     calculate local tangent
            dres = (betadp*gammadp)/(ndp*fdfunc)
            dres = dres*(ndp*(phinloc - gammadp*Knp1)-gammadp) - one

!     write(*,*) "damage",count,residual
c     compute the increment 
            gamma1 = -residual/dres
c     !update K
            Knp1 = Knp1 + gamma1;
         end do                 ! do while loop 

      endif
      
      call damageFunction(fdfunc,fdfuncp,d,Knp1)
!     if (n.eq.2) then
!     write(*,*) n, fdfunc
!     endif

!     write(*,*) "Damage:" ,damageflg,"Kn:" ,Kn, "  Knp1",
!     1   Knp1,"Kd",Kd, "phinloc", phinloc   

      if (fdfunc.lt.0.01) then
         write(*,*) "_WARNING!!_",chk, Kn, Knp1, residual
      endif

      if (Kn.gt.Knp1) then
         write(*,*) "_________ HERE ____________"
         Knp1 = Kn 
         call damageFunction(fdfunc,fdfuncp,d,Knp1)
      endif 
!     Update history

      hn1(1) = Knp1   

c     step6 -> Compute flux/source terms
!     a stress ~ updated fdfunc  
!|||| 
!|||  FOR FORCE MATRIX
!|||| 
      do i =1,3
         Yvector(i,1) = cdp*Gradphi(i,1) 
         do j=1,3
            sigD(i,j)    = (fdfunc+0.0001)*stressE(i,j)
         enddo
      enddo
      yscalar = betadp*(phinloc - gammadp* Knp1) 
      
!|||| 
!|||  FOR STIFFNESS MATRIX
!|||| 
c     step 7 -> Compute tangent moduli
      if (damageflg) then
!     
         thetadam = -ndp*fdfunc - gammadp*gammadp*betadp + 
     1        gammadp*betadp*ndp*(phinloc - gammadp*Knp1)
         
         thetadam = ndp*fdfunc/thetadam
!     
      else 
         thetadam = zero
      endif
c     write(iow,*) "thetadam",thetadam
!     
      kc = ndp*fdfunc*thetadam ! damage part tangent scaling multiplier
!     
!     write(*,*) "kc",kc
!     
      do i=1,3
         do j=1,3
            do k=1,3
               do m=1,3
                  edam(i,j,k,m) = kc*stressE(i,j)*stressE(k,m)          
               enddo
            enddo
         enddo
      enddo  


      do kk=1,3
         do mm=1,3
            do pp=1,3
               do qq=1,3
                  CmatD(kk,mm,pp,qq)=(fdfunc+0.0001)*CmatE(kk,mm,pp,qq)+ 
     &                 edam(kk,mm,pp,qq) 
               end do
            end do
         end do
      end do
!     b Second (middle) quantity needed for the couple part of global tangent 
      kc2 = betadp*gammadp*thetadam
      do i=1,3
         do j=1,3
            dsigdphi(i,j) = kc2*stressE(i,j) ! 2 dy/dg is same! leading to the symmetry of the tangent at the end. 
         enddo
      enddo
      
!     c Fourth (middle) quantity needed for the couple part of global tangent 
      dyscalardphi = betadp*gammadp*gammadp*thetadam  
      dyscalardphi = one + dyscalardphi/(ndp*fdfunc)   
      dyscalardphi = betadp*dyscalardphi

c     write(*,*)  "chk" , chk,  "phiInterp",phinloc, "Kn+1", Knp1,
c     1  "thetadam", thetadam, "fdfunc", fdfunc, "ndp", ndp,
c     2  "detF", detF

!     d Fifth (middle) quantity needed for the couple part of global tangent 
      do i=1,3
         do j=1,3
            dydgradphi(i,j)= (cdp)*id2(i,j) 
         enddo
      enddo
      
c     write(iow,*)  "tehtadam",thetadam, "dsigdphi"
c     write(iow,*)  "dsigdphi"
c     write(iow,'(3f20.10 /)')  ((dsigdphi(i,j),i=1,3),j=1,3)
c     write(iow,*)  "dydgradphi"
c     write(iow,'(3f20.10 /)')  ((dydgradphi(i,j),i=1,3),j=1,3)
c     write(iow,*)  "dyscalardphi",dyscalardphi*detf
c     c
c     
c     write(iow,*) "Knp1",Knp1, "phinloc", phinloc, "chk",chk

c     write(*,*) "hn(1 )", hn(1 )
c     write(*,*) "hn1(1)", hn1(1)
      

c     write(*,'(6f10.4 /)')  ((dsigdphi(i,j),i=1,3),j=1,3)

      return
      END
c     ----------------------------------------------------------------------------

c     ----------------------------------------------------------------------------
      
      subroutine getGradshapePhi(GradshapePhi,GradshapePhiT,dshapeGE,
     *     xjacinv,ndim,nnodeGE)
      
      IMPLICIT NONE
c     ------------------------------------------
c     ------------------------------------------
c     input - dshape (nnodeGE,ndim)
c     xjacinvC (ndim,ndim)
c     ndim 
c     nnodeGE
c     _ 
c     output- GradshapePhi  (nnodege,ndim) V_x PHI
c     GradshapePhiT  
c     ------------------------------------------
c     -----------------------------------------

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES 
      
!------O U T P U T  
      REAL*8 GradshapePhi(nnodeGE,ndim)
      REAL*8 GradshapePhiT(ndim,nnodeGE)
      REAL*8 dshapeGE(nnodeGE,ndim)
!------O U T P U T  
      
!------I N P U T  
      REAL*8 xjacinv(ndim,ndim)
      INTEGER ndim
      INTEGER nnodeGE           ! node number of GE nonlocal part
!------I N P U T  
      
      REAL*8 dNidx,dNidy,dNidz
      
!     initialization
      do i=1,nnodeGE
         do j=1,ndim
            GradshapePhi(i,j) = zero  
         enddo
      enddo 
      
      do  i=1,nnodeGE
         
         dNidx = dshapeGE(i,1)*xjacinv(1,1) + 
     *        dshapeGE(i,2)*xjacinv(2,1) + 
     *        dshapeGE(i,3)*xjacinv(3,1)
         
         dNidy = dshapeGE(i,1)*xjacinv(1,2) +
     *        dshapeGE(i,2)*xjacinv(2,2) + 
     *        dshapeGE(i,3)*xjacinv(3,2)
         
         dNidz = dshapeGE(i,1)*xjacinv(1,3) + 
     *        dshapeGE(i,2)*xjacinv(2,3) + 
     *        dshapeGE(i,3)*xjacinv(3,3) 
         
         
         GradshapePhi(i,1) =  dNidx
         GradshapePhi(i,2) =  dNidy
         GradshapePhi(i,3) =  dNidz
         
      enddo
!     Compute transpose
      call getTranspose3D(GradshapePhi,ndim,nnodeGE,GradshapePhiT)
      
      RETURN
      END


c     ----------------------------------------------------------------------------
!-----------------------------------------------------
      
      subroutine getGradPhiInterp(PhiInterp,Gradphi,GradphiT,
     1     GradshapePhi, pVec,shapef,dshape,nnodeGE,ndim)
      
      IMPLICIT NONE
      
c     ------------------------------------------
c     ------------------------------------------
c     input - GradshapePhi (nnodeGE,ndim)
c     pVec (nnodeGE,1)
c     shapef (nnodeGE)
c     dshape (nnodeGE,ndim)
c     nnodeGE 
c     ndim
c     
c     output- phiInterp
c     GradPhi  (ndim,1)
c     GradPhiT  (1,ndim)
c     
c     ------------------------------------------
c     -----------------------------------------

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      REAL*8 pVec(nnodeGE,1)
      REAL*8 Gradphi(ndim,1)
      REAL*8 GradphiT(1,ndim)
      REAL*8 PhiInterp
      
      REAL*8 dshape(nnodeGE,ndim)
      REAL*8 dshapeT(ndim,nnodeGE)
      REAL*8 shapef(nnodeGE)
      
      REAL*8 GradshapePhi(nnodeGE,ndim)
      REAL*8 GradshapePhiT(ndim,nnodeGE)

c     REAL*8 sum
      INTEGER nnodeGE,ndim
      
      do i=1,ndim
         Gradphi(i,1) = zero
         GradphiT(1,i) = zero
      enddo

      PhiInterp = zero     
      do i=1,nnodeGE
         PhiInterp  = PhiInterp + shapef(i)*pVec(i,1)
      enddo  
      
      call getTranspose3D(GradshapePhi,ndim,nnodeGE,GradshapePhiT)  
!     kCompTranspose20Node(A,m,n,AT)  

      call getMatrixMultiply(Gradphi,ndim,1,
     *     GradshapePhiT,ndim,nnodeGE,pVec,nnodeGE,1)

!     call mprint(GradshapePhiT,ndim,nnodeGE,ndim,"GradshapePhi")
!     write(*,*) "Gradphi"
!     write(*,'(1f20.10)') (Gradphi(i,1), i =1,3)

      
      do i=1,ndim
         GradphiT(1,i) = Gradphi(i,1)
      enddo
      
!     be careful about transpose!! 

      RETURN
      END
      
c-----------------------------------------------------
c-----------------------------------------------------

c-----------------------------------------------------
c-----------------------------------------------------
      
      subroutine getKuphi(Kdg,Bmat,dsigdphi,shapefGE,fact,
     *     ndim,nnode,nnodeGE)
      
      IMPLICIT NONE

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      REAL*8 Kdg(ndim*nnode,nnodeGE) !60x8
      REAL*8 Bmat(6,ndim*nnode)
      REAL*8 BmatT(nnode*ndim,6)
      
      REAL*8 dsigdphi(ndim,ndim)
      REAL*8 dsigdphiVoigt(6,1) 
      REAL*8 shapefGE(nnodeGE)
      REAL*8 shapefGEvect(1,nnodeGE)
      REAL*8 fact
      
      INTEGER ndim,nnode,nnodeGE
      
      do i=1,ndim*nnode
         do j=1,nnodeGE
            Kdg(i,j) =  zero   
         enddo
      enddo  
      
!     Compute BmatT
      call getTranspose3D(Bmat,nnode*ndim,6,BmatT)  
      
!     Arrange shapefGE in matrix form
      do i=1,nnodeGE
         shapefGEvect(1,i) = zero
         shapefGEvect(1,i) = shapefGE(i)
      enddo
      
!     Arrange stress in Voigt notation
      dsigdphiVoigt(1,1) =dsigdphi(1,1) 
      dsigdphiVoigt(2,1) =dsigdphi(2,2) 
      dsigdphiVoigt(3,1) =dsigdphi(3,3) 
      dsigdphiVoigt(4,1) =dsigdphi(1,2) 
      dsigdphiVoigt(5,1) =dsigdphi(2,3) 
      dsigdphiVoigt(6,1) =dsigdphi(1,3) 
      
c     This subroutine computing 
c     A = B*C*D  
      call  getMatrixTriply(Kdg,ndim*nnode,nnodeGE,
     *     BmatT,ndim*nnode,6,dsigdphiVoigt,6,1,shapefGEvect,1,nnodeGE)
      
      do i=1,ndim*nnode
         do j=1,nnodeGE
            Kdg(i,j) =   Kdg(i,j)*fact   
         enddo
      enddo  

      RETURN
      END
      


c---------------------------------
c---------------------------------
      subroutine  getMatrixTriply(Ares,a1,a2,B,b1,b2,
     &     C,c1,c2,Trd,d1,d2)  

c     This subroutine computing 
c     A = B*C*D
c     
      IMPLICIT NONE

      integer a1, a2
      integer b1, b2
      integer c1, c2
      integer d1, d2

      REAL*8 Ares(a1,a2)        !result matrix 
      REAL*8 B(b1,b2)           !first matrix
      REAL*8 C(c1,c2)           !second matrix        
      REAL*8 Trd(d1,d2)         !third matrix        

      REAL*8 Aux(c1,d2)         !auxiliar matrix for the computation

      REAL*8,parameter :: zero = 0.d0
      REAL*8,parameter :: one = 1.d0
      REAL*8,parameter :: two = 2.d0

      INTEGER i,j,k

      REAL*8 sump, sump2


!     initialization to zero of the auxiliar matrix
      do i=1,c1
         do j=1,d2
            Aux(i,j) = zero
         enddo
      enddo
      
      do i=1,a1
         do j=1,a2
            Ares(i,j) = zero
         enddo
      enddo

      if(b2-c1.ne.zero) then    !It can not be multiplied
         write(7,*)'No matrix multiplication'
      end if 
      
      if(c2-d1.ne.zero) then    !It can not be multiplied
         write(7,*)'No matrix multiplication'
      end if 
      
      
!     write(7,*) "Before first multiplication in triple"
      call getMatrixMultiply(Aux,c1,d2,C,c1,c2,Trd,d1,d2)
!     write(7,*) "Before second multiplication in triple"  
      call getMatrixMultiply(Ares,b1,d2,B,b1,c1,Aux,c1,d2)
!     write(7,*) "After second multiplication in triple"  



      return 
      end

!---------------------------------------------------------
!---------------------------------------------------------
!-----------------------------------------------------
      
      subroutine getKphiphi(kphiphi,GradshapePhi,GradshapePhiT,
     *     shapefGE,dyscalardphi,dydgradphi,fact,nnodeGE,ndim) 
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      
      REAL*8 kphiphi(nnodeGE,nnodeGE)
      REAL*8 kphiphiAux(nnodeGE,nnodeGE)
      REAL*8 GradshapePhi(nnodeGE,ndim)
      REAL*8 GradshapePhiT(ndim,nnodeGE)
      
      REAL*8 shapefGE(nnodeGE)
      REAL*8 shapefGEvect(1,nnodeGE)
      REAL*8 shapefGEvectT(nnodeGE,1)
      
      REAL*8 dyscalardphi
      REAL*8 dydgradphi(ndim,ndim)
      
      REAL*8 fact
      INTEGER ndim,nnodeGE
      
      
!     initialization
      do i=1,nnodeGE
         do j=1,nnodeGE
            kphiphi(i,j) =  zero  
            kphiphiAux(i,j) =  zero  
         enddo
      enddo  
      
      
!     Arrange shapefGE in matrix form
      do i=1,nnodeGE
         shapefGEvect(1,i) = zero
         shapefGEvect(1,i) = shapefGE(i)
         shapefGEvectT(i,1) = zero
         shapefGEvectT(i,1) = shapefGE(i)
      enddo
      
!     Compute first term N^phiT*N^phi
      call getMatrixMultiply(kphiphi,nnodeGE,nnodeGE,
     *     shapefGEvectT,nnodeGE,1,shapefGEvect,1,nnodeGE)
      
      
      do i=1,nnodeGE
         do j=1,nnodeGE
            kphiphi(i,j) =   kphiphi(i,j)*dyscalardphi*fact 
c     
c     write(*,*) "dyscalardphi",dyscalardphi
         enddo
      enddo 
      
      
!     Compute the second term:  grad N^phiT *dydgradphi*grad N^phi
!     This subroutine computing 
!     A = B*C*D   

      call  getMatrixTriply(kphiphiAux,nnodeGE,nnodeGE,
     *     GradshapePhi,nnodeGE,ndim,dydgradphi,ndim,ndim,
     *     GradshapePhiT,ndim,nnodeGE)
      
      
      do i=1,nnodeGE
         do j=1,nnodeGE
            kphiphiAux(i,j) =   kphiphiAux(i,j)*fact 
            
c     write(*,*) kphiphiAux(i,j) 

         enddo
      enddo 
      
c     write(*,*) kphiphiAux(i,j) 

      do i=1,nnodeGE
         do j=1,nnodeGE
            kphiphi(i,j) =   kphiphi(i,j) + kphiphiAux(i,j)  
         enddo
      enddo 
      
      
      RETURN
      END

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------
      
      subroutine getFphi(fphi,btensor,yvector,yscalar,shapef,
     *     dshape,nnodeGE,ndim,fact)

c     !  subroutine getFphi(fphi,yvector,dshape,nnodeGE,ndim,fact,yscalar)
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      INTEGER nnodeGE,ndim
      REAL*8 btensor(3,3)
      REAL*8 yscalar
      REAL*8 shapef(nnodeGE)
      REAL*8 fphi(nnodeGE,1)
      REAL*8 fphiExt(nnodeGE,1)
      REAL*8 yvector(ndim,1) 
c     REAL*8 yvectorAux(ndim,1) 
      REAL*8 dshape(nnodeGE,ndim)
      REAL*8 dshapeT(ndim,nnodeGE)
      REAL*8 BtensorDshape(ndim,nnodeGE)
      REAL*8 BtensorDshapeT(nnodeGE,ndim)
      REAL*8 fact
      
!initialization 
      do i=1,nnodeGE
         fphi(i,1) = zero
      enddo 

      BtensorDshape(1:ndim,1:nnodeGE) =zero
      BtensorDshapeT(1:nnodeGE,1:ndim) =zero
      dshapeT(1:ndim,1:nnodeGE) =zero

!     do i=1,ndim
!     yvectorAux(i,1) = yvector(i)
!     write(*,*) "yvector", yvector(i)
!     enddo

      call getTranspose3D(dshape,ndim,nnodeGE,dshapeT) 
      
      call getMatrixMultiply(BtensorDshape,ndim,nnodeGE,
     *     btensor,ndim,ndim,dshapeT,ndim,nnodeGE)

      call getTranspose3D(BtensorDshape,nnodeGE,ndim,BtensorDshapeT) 

      call getMatrixMultiply(fphi,nnodeGE,1,
     *     BtensorDshapeT,nnodeGE,ndim,yvector,ndim,1)
      

!call mprint(fphi,8,1,8,"Fphi") 
! write(*,*) "yvector"
! write(*,'(1f20.10 /)') (yvector(i,1), i =1,3)

      do i=1,nnodeGE
         fphiExt(i,1) =  yscalar*shapef(i)
      enddo 

      do i=1,nnodeGE
         fphi(i,1) = (fphi(i,1) + fphiExt(i,1))*fact
      enddo 

! do i=1,nnodeGE
!   fphi(i,1) = fphi(i,1)*fact
! enddo 
      
      RETURN 
      END

c     ----------------------------------------------------------------------------
!-----------------------------------------------------
      
      subroutine getFphiReference(fphi,yvector,yscalar,shapef,
     *     dshape,nnodeGE,ndim,fact)

c     !  subroutine getFphi(fphi,yvector,dshape,nnodeGE,ndim,fact,yscalar)
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      INTEGER nnodeGE,ndim
      REAL*8 yscalar
      REAL*8 shapef(nnodeGE)
      REAL*8 fphi(nnodeGE,1)
      REAL*8 fphiExt(nnodeGE,1)
      REAL*8 yvector(ndim,1) 
c     REAL*8 yvectorAux(ndim,1) 
      REAL*8 dshape(nnodeGE,ndim)
      REAL*8 dshapeT(ndim,nnodeGE)
      REAL*8 BtensorDshape(ndim,nnodeGE)
      REAL*8 BtensorDshapeT(nnodeGE,ndim)
      REAL*8 fact
      
!     initialization 
      do i=1,nnodeGE
         fphi(i,1) = zero
      enddo 

c$$$      BtensorDshape(1:ndim,1:nnodeGE) =zero
c$$$      BtensorDshapeT(1:nnodeGE,1:ndim) =zero
c$$$      dshapeT(1:ndim,1:nnodeGE) =zero

!     do i=1,ndim
!     yvectorAux(i,1) = yvector(i)
!     write(*,*) "yvector", yvector(i)
!     enddo

c$$$      call getTranspose3D(dshape,ndim,nnodeGE,dshapeT) 
c$$$      
c$$$      call getMatrixMultiply(BtensorDshape,ndim,nnodeGE,
c$$$     *     btensor,ndim,ndim,dshapeT,ndim,nnodeGE)
c$$$
c$$$      call getTranspose3D(BtensorDshape,nnodeGE,ndim,BtensorDshapeT) 
c$$$
c$$$      call getMatrixMultiply(fphi,nnodeGE,1,
c$$$     *     BtensorDshapeT,nnodeGE,ndim,yvector,ndim,1)

      call getMatrixMultiply(fphi,nnodeGE,1,
     *     dshape,nnodeGE,ndim,yvector,ndim,1)
      
      !call mprint(dshape,nnodeGE,ndim,nnodeGE,"dshape8,3")
!     call mprint(fphi,8,1,8,"Fphi") 
!     write(*,*) "yvector"
!     write(*,'(1f20.10 /)') (yvector(i,1), i =1,3)

      do i=1,nnodeGE
         fphiExt(i,1) =  yscalar*shapef(i)
      enddo 

      do i=1,nnodeGE
         fphi(i,1) = (fphi(i,1) + fphiExt(i,1))*fact
      enddo 

!     do i=1,nnodeGE
!     fphi(i,1) = fphi(i,1)*fact
!     enddo 
      
      RETURN 
      END

c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------


      subroutine  damageFunction(fdfunc,fdfuncp,d,K)
! DAMAGE FUNCTION

      IMPLICIT NONE

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0

! output 
      real*8 fdfuncp 
      real*8 fdfunc
! input
      REAL*8 d(*)
      real*8 K, Kd
      real*8 ndp
! material parameters 
      Kd = d(243)
      ndp = d(244)

      fdfunc = exp(ndp*(Kd-K))
      fdfuncp = - ndp *fdfunc
      
!write(*,*) "in fdamage"
!write(*,*) "ndp",ndp,"Kd",Kd,"K",K,"fdfunc",fdfunc
c     *       ndp  =  0.010d0
c     *       Kd   =  d(243) 
c     *       dinf = 0.98d0
c     *
c     *       fdfunc = (one -dinf*(one-exp(ndp*(Kd-K))))
c     *       fdfuncp = -dinf*ndp*exp(ndp*(Kd - K))
      end

c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------

!     (sig8x8, shp3, xsjdv, wgp, mgauss, nel, 8, r, s)
      subroutine gaussValueProjection(sig, shp3, xsj, lint, 
     1     nel, nes, r, s)

c     Outputs:
c     p(nen)       - Weights for 'lumped' projection
c     s(nen,*)     - Integral of variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
c     include  'qudshp.h'
      include  'strnum.h'       ! contains iste

      integer   ii, jj, l, nel, lint, nes
      real*8    r(*),s(nen,*), sig(nes,lint), xj
      real*8    xsj(lint), shp3(nel,*)

      save

c     Initialize arrays
      do ii = 1,nel
         r(ii)    = 0.0d0
!     s(ii,10) = 0.0d0
         do jj = 1,8            ! 6
            s(ii,jj) = 0.0d0
         end do                 ! jj
      end do                    ! ii

c     Compute projections: int ( sig * shp(i) * darea )
      do l = 1,lint
         do ii = 1,nel          ! nnode
            xj    = xsj(l)*shp3(ii,l)
            r(ii) = r(ii)   + xj

            s(ii,1) = s(ii,1) + sig(1,l)*xj
            s(ii,2) = s(ii,2) + sig(2,l)*xj
            s(ii,3) = s(ii,3) + sig(3,l)*xj
            s(ii,4) = s(ii,4) + sig(4,l)*xj
            s(ii,5) = s(ii,5) + sig(5,l)*xj
            s(ii,6) = s(ii,6) + sig(6,l)*xj
            s(ii,7) = s(ii,7) + sig(7,l)*xj
            s(ii,8) = s(ii,8) + sig(8,l)*xj

         end do                 ! ii
      end do                    ! l

c     write(*,*) 'subroutine inside'
c     write(*,'(7f10.4 /)') ((sig(jj,l) , jj = 1, 1),l = 1, lint)        !((s(ii,jj), ii = 1, nel),jj = 1, 1)

      iste = 8

      end

      subroutine getShape3D8(N,xi,ndim)

      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      REAL*8 N(8)
      REAL*8 xi(ndim)
      INTEGER ndim
      
! do i=1,8
! N(1:8) = zero
! end do
      
      N(1) = (one-xi(1))*(one-xi(2))*(one-xi(3))/eight;
      N(2) = (one+xi(1))*(one-xi(2))*(one-xi(3))/eight;
      N(3) = (one+xi(1))*(one+xi(2))*(one-xi(3))/eight;
      N(4) = (one-xi(1))*(one+xi(2))*(one-xi(3))/eight;
      N(5) = (one-xi(1))*(one-xi(2))*(one+xi(3))/eight;
      N(6) = (one+xi(1))*(one-xi(2))*(one+xi(3))/eight;
      N(7) = (one+xi(1))*(one+xi(2))*(one+xi(3))/eight;
      N(8) = (one-xi(1))*(one+xi(2))*(one+xi(3))/eight;
c     write(*,*) (N(i),i=1,20)

      
      end
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

!---------------------------------------------------------
      subroutine getDerivativeShape3D8(dNdxi,xi,ndim)

      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      REAL*8 dNdxi(8,ndim)
      REAL*8 xi(ndim)
      INTEGER ndim
      
!initialization 
! do i=1,20    ! #NODE
!     do j=1,ndim
!  dNdxi(i,j) = zero
!     enddo
! enddo         
      
      dNdxi(1,1) = -(one-xi(2))*(one-xi(3))/eight;
      dNdxi(1,2) = -(one-xi(1))*(one-xi(3))/eight;
      dNdxi(1,3) = -(one-xi(1))*(one-xi(2))/eight;
      dNdxi(2,1) =  (one-xi(2))*(one-xi(3))/eight;
      dNdxi(2,2) = -(one+xi(1))*(one-xi(3))/eight;
      dNdxi(2,3) = -(one+xi(1))*(one-xi(2))/eight;
      dNdxi(3,1) =  (one+xi(2))*(one-xi(3))/eight;
      dNdxi(3,2) =  (one+xi(1))*(one-xi(3))/eight;
      dNdxi(3,3) = -(one+xi(1))*(one+xi(2))/eight;
      dNdxi(4,1) = -(one+xi(2))*(one-xi(3))/eight;
      dNdxi(4,2) =  (one-xi(1))*(one-xi(3))/eight;
      dNdxi(4,3) = -(one-xi(1))*(one+xi(2))/eight;
      dNdxi(5,1) = -(one-xi(2))*(one+xi(3))/eight;
      dNdxi(5,2) = -(one-xi(1))*(one+xi(3))/eight;
      dNdxi(5,3) =  (one-xi(1))*(one-xi(2))/eight;
      dNdxi(6,1) =  (one-xi(2))*(one+xi(3))/eight;
      dNdxi(6,2) = -(one+xi(1))*(one+xi(3))/eight;
      dNdxi(6,3) =  (one+xi(1))*(one-xi(2))/eight;
      dNdxi(7,1) =  (one+xi(2))*(one+xi(3))/eight;
      dNdxi(7,2) =  (one+xi(1))*(one+xi(3))/eight;
      dNdxi(7,3) =  (one+xi(1))*(one+xi(2))/eight;
      dNdxi(8,1) = -(one+xi(2))*(one+xi(3))/eight;
      dNdxi(8,2) =  (one-xi(1))*(one+xi(3))/eight;
      dNdxi(8,3) =  (one-xi(1))*(one+xi(2))/eight;

      
      end
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine getJacobian3D(xjac,detxjac,xjacinv,
     *     dshape,ndim,nnode,coord)

      IMPLICIT NONE

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0

C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 
      include  'iofile.h'       ! iow

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      INTEGER ndim,nnode
      
      REAL*8 xjac(ndim,ndim)
      REAL*8 xjacinv(ndim,ndim)
      REAL*8 detxjac
      
      REAL*8 coord(ndim,nnode)
      REAL*8 dshape(nnode,ndim)
      
!initialization
      do i=1,ndim
         do j=1,ndim
            xjac(i,j) =zero 
            xjacinv(i,j) =zero 
         enddo
      enddo 
      
! detxjac = zero
      
!write(iow,*) "nnode in Jacob", nnode
!write(iow,*) "ndim in Jacob", ndim
      
!write(iow,*) "coord in Jacobian"
!write(iow,*) "", coord
      
      
      do k=1,nnode
!write(iow,*) "k", k
!write(iow,*) "", coord(1,k)
!write(iow,*) "", dshape(k,1)
         
         xjac(1,1) =  xjac(1,1) + coord(1,k)*dshape(k,1)
         xjac(1,2) =  xjac(1,2) + coord(1,k)*dshape(k,2) 
         xjac(1,3) =  xjac(1,3) + coord(1,k)*dshape(k,3) 
         
         xjac(2,1) =  xjac(2,1) + coord(2,k)*dshape(k,1)
         xjac(2,2) =  xjac(2,2) + coord(2,k)*dshape(k,2) 
         xjac(2,3) =  xjac(2,3) + coord(2,k)*dshape(k,3) 
         
         xjac(3,1) =  xjac(3,1) + coord(3,k)*dshape(k,1)
         xjac(3,2) =  xjac(3,2) + coord(3,k)*dshape(k,2) 
         xjac(3,3) =  xjac(3,3) + coord(3,k)*dshape(k,3)

      enddo
      
c     write(iow,*) "---  xjac in Jacob  ---"
c     write(iow,*) "", xjac
c     write(iow,*) "---    ---"
!COMPUTE DETERMINANT OF THE JACOBIAN
      call getDeterminant3D(detxjac,xjac)
      
      if(detxjac.gt.zero) then  !It can be formulated the inverse
!COMPUTE THE INVERSE 
         call getInverse3D(xjacinv,xjac,detxjac)
      else 
C     negative or zero jacobian
         write(iow,*)'Warning: element has negative jacobian'
         write(*,*)'Warning: element has negative jacobian'
      end if 
      

      END

      
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine getJacobian3DCenter(xjac,detxjac,
     *     dshape,ndim,nnode,coord)

      IMPLICIT NONE

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0

C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 
      include  'iofile.h'       ! iow

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      INTEGER ndim,nnode
      
      REAL*8 xjac(ndim,ndim)
!REAL*8 xjacinv(ndim,ndim)
      REAL*8 detxjac
      
      REAL*8 coord(ndim,nnode)
      REAL*8 dshape(nnode,ndim)
      
!initialization
      do i=1,ndim
         do j=1,ndim
            xjac(i,j) =zero 
!xjacinv(i,j) =zero 
         enddo
      enddo 
      
! detxjac = zero
      
!write(iow,*) "nnode in Jacob", nnode
!write(iow,*) "ndim in Jacob", ndim
      
!write(iow,*) "coord in Jacobian"
!write(iow,*) "", coord
      
      
      do k=1,nnode
!write(iow,*) "k", k
!write(iow,*) "", coord(1,k)
!write(iow,*) "", dshape(k,1)
         
         xjac(1,1) =  xjac(1,1) + coord(1,k)*dshape(k,1)
         xjac(1,2) =  xjac(1,2) + coord(1,k)*dshape(k,2) 
         xjac(1,3) =  xjac(1,3) + coord(1,k)*dshape(k,3) 
         
         xjac(2,1) =  xjac(2,1) + coord(2,k)*dshape(k,1)
         xjac(2,2) =  xjac(2,2) + coord(2,k)*dshape(k,2) 
         xjac(2,3) =  xjac(2,3) + coord(2,k)*dshape(k,3) 
         
         xjac(3,1) =  xjac(3,1) + coord(3,k)*dshape(k,1)
         xjac(3,2) =  xjac(3,2) + coord(3,k)*dshape(k,2) 
         xjac(3,3) =  xjac(3,3) + coord(3,k)*dshape(k,3)

      enddo
      
c     write(iow,*) "---  xjac in Jacob  ---"
c     write(iow,*) "", xjac
c     write(iow,*) "---    ---"
!COMPUTE DETERMINANT OF THE JACOBIAN
      call getDeterminant3D(detxjac,xjac)
      
!if(detxjac.gt.zero) then !It can be formulated the inverse
!COMPUTE THE INVERSE 
!  call getInverse3D(xjacinv,xjac,detxjac)
!else 
C     negative or zero jacobian
!     write(iow,*)'Warning: element has negative jacobian', detxjac
!     write(*,*)'Warning: element has negative jacobian', detxjac
!end if 
      

      END

      
!---------------------------------------------------------

!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getDeterminant3D(det,A)

!ROUTINE TO COMPUTE THE DETERMINANT 3x3 MATRIX
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      
      REAL*8 A(3,3)
      REAL*8 det

      
      det = zero
      
      det    = A(1,1)*A(2,2)*A(3,3) +
     *     A(1,3)*A(2,1)*A(3,2) +
     *     A(3,1)*A(1,2)*A(2,3) - 
     *     A(3,1)*A(2,2)*A(1,3) -
     *     A(3,3)*A(1,2)*A(2,1) -
     *     A(1,1)*A(2,3)*A(3,2)
      
c     write(*,*) 'Determinant jacobian', det

      END

!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getInverse3D(Ainv,A,det)
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0

      REAL*8 A(3,3),Ainv(3,3)
      REAL*8 det
      
! do i=1,3
!   do j=1,3
!     Ainv(i,j) = zero  
!   enddo
! enddo  
      
      Ainv(1,1) =  (A(2,2) *A(3,3) - A(3,2) *A(2,3))/det
      Ainv(1,2) = -(A(1,2) *A(3,3) - A(1,3) *A(3,2))/det
      Ainv(1,3) =  (A(1,2) *A(2,3) - A(2,2) *A(1,3))/det
      Ainv(2,1) = -(A(2,1) *A(3,3) - A(3,1) *A(2,3))/det
      Ainv(2,2) =  (A(1,1) *A(3,3) - A(1,3) *A(3,1))/det
      Ainv(2,3) = -(A(1,1) *A(2,3) - A(2,1) *A(1,3))/det
      Ainv(3,1) =  (A(2,1) *A(3,2) - A(3,1) *A(2,2))/det
      Ainv(3,2) = -(A(1,1) *A(3,2) - A(3,1) *A(1,2))/det
      Ainv(3,3) =  (A(1,1) *A(2,2) - A(1,2) *A(2,1))/det
      
      
      END
      
!---------------------------------------------------------
!---------------------------------------------------

      subroutine getDeformGrad3D(Fgrad,Finv,FgradT,detF,xjacR,detxjacR,
     *     xjacinvR,xjacC,detxjacC,xjacinvC,ndim,nnode)

      IMPLICIT NONE

      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four=4.d0
      REAL*8, parameter :: eight=8.d0

      include 'iofile.h'
c----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      REAL*8 Fgrad(3,3)         !Deformation gradient
      REAL*8 Finv(3,3)          !Inverse deformation gradient
      REAL*8 FgradT(3,3)        !Deformation gradient transpose
      REAL*8 detF               !determinant deformation gradient 
      
      REAL*8 xjacR(ndim,ndim)
      REAL*8 xjacinvR(ndim,ndim)
      REAL*8 detxjacR
      
      REAL*8 xjacC(ndim,ndim)
      REAL*8 xjacinvC(ndim,ndim)
      REAL*8 detxjacC
      
      INTEGER nnode,ndim
      
!initialization 
      do i=1,ndim
         do j=1,ndim
            Fgrad(i,j) = zero
            Finv(i,j) = zero
            FgradT(i,j) = zero
         enddo
      enddo  
      
! Fe = je*Je^-1 --> by definition according to Wriggers
      call getMatrixMultiply(Fgrad,ndim,ndim,xjacC,ndim,ndim,
     *     xjacinvR,ndim,ndim) 
      
!     write(99,*) "Fgrad"
!     do i=1,ndim
!     write(99,*) "",Fgrad(i,1),Fgrad(i,2),Fgrad(i,3)
!     enddo
!     write(99,*) "-------"
      
!detFe = det je/det Je
      detF = detxjacC/detxjacR 
      
!     write(iow,*) "detF", detF
      
!Compute the inverse of the deformation gradient
      call getInverse3D(Finv,Fgrad,detF)
      
c     write(iow,*) "Finv"
c     do i=1,ndim
c     write(iow,*) "",Finv(i,1),Finv(i,2),Finv(i,3)
c     enddo
c     write(iow,*) "-------"
      
!Compute transverse deformation gradient
      call getTranspose3D(Fgrad,3,3,FgradT)
      
      
      END

!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getMatrixMultiply(A,a1,a2,B,b1,b2,
     &     C,c1,c2)  

c     This subroutine stdcomputing 
c     A = B*C
c     
      IMPLICIT NONE

      include 'iofile.h'

      integer a1, a2
      integer b1, b2
      integer c1, c2

      REAL*8 A(a1,a2)           !result matrix 
      REAL*8 B(b1,b2)           !first matrix
      REAL*8 C(c1,c2)           !second matrix        

      REAL*8,parameter :: zero = 0.d0
      REAL*8,parameter :: one = 1.d0
      REAL*8,parameter :: two = 2.d0

      INTEGER i,j,k

      REAL*8 sump

!initialization to zero of the result matrix
      do i=1,a1
         do j=1,a2
            A(i,j) = zero    
         enddo
      enddo   

      if(b2-c1.ne.zero) then    !It can not be multiplied
         write (*,*)'No matrix multiplication'
      end if 

c     operation = factor1*factor2  !operation

      do i = 1,b1                
         do j = 1,c2
            sump = zero
            do k=1,c1
               sump=sump+B(i,k)*C(k,j)
            end do 
            A(i,j)=sump
         end  do 
      end do
      

      end

!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getTranspose3D(A,m,n,AT)  

      implicit none
      
      REAL*8 A(n,m)
      REAL*8 AT(m,n)
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      
      INTEGER n,m
      INTEGER i,j,k             !indices
      
      do i=1,m
         do j=1,n
            AT(i,j) = zero      !initialization
         end do
      end do
      
      do i=1,n
         do j=1,m
            AT(j,i) = A(i,j)
         end do
      end do
      
      END
!---------------------------------------------------


!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getBmat3D(Bmat,dshape,xjac,detxjac,xjacinv,ndim,
     *     nnode,ndof)
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      INTEGER ndim,nnode,ndof
      
      REAL*8  Bmat(6,nnode*ndof) 
      INTEGER kint
      REAL*8  dNidx, dNidy, dNidz
      
      REAL*8  xjac(ndim,ndim)
      REAL*8  xjacinv(ndim,ndim)
      REAL*8  detxjac
      
      REAL*8  dshape(nnode,ndim)
      
!initialization 
      do i=1,6
         do j=1,nnode*ndof
            Bmat(i,j) = zero
         enddo
      enddo
      
      
!Linear B-operator looks like that
! Nk,x    0       0|
!   0    Nk,y     0|
!   0     0    Nk,z| 

!Nk,z    0     Nk,x|
!   0    Nk,z  Nk,y|
!Nk,y    Nk,x     0|  
      

      kint = 1   

      do  i=1,nnode

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

c     Bmat(4,kint)   = dNidz
c     Bmat(4,kint+2) = dNidx
c     
c     Bmat(5,kint+1) = dNidz
c     Bmat(5,kint+2) = dNidy
c     
c     Bmat(6,kint)   = dNidy
c     Bmat(6,kint+1) = dNidx
c     
!update the counter
         kint = kint + 3

      end do
      
      
      END
      
!---------------------------------------------------------
!---------------------------------------------------------

      subroutine getBmat3DNL(BmatNL,Bop,Fgrad,xjacinv)


!----------Computation of nonlinear part of Bmat  
      
      IMPLICIT NONE
      
      REAL*8, parameter :: zero=0.d0
      REAL*8, parameter :: half=0.5d0
      REAL*8, parameter :: quarter=0.25d0
      REAL*8, parameter :: one=1.d0
      REAL*8, parameter :: two=2.d0
      REAL*8, parameter :: three=3.d0
      REAL*8, parameter :: four =4.d0
      
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C

      INTEGER i,j,k,l
      
      REAL*8 BmatNL(6,24)       !nonlinear B-operator
      REAl*8 Bop(6,24)          !Linear part B-operator


      REAL*8 Fgrad(3,3)         !Deformation gradient
      REAL*8 xjacinv(3,3)       !inverse of the jacobian
      

!     write(99,*) "DefGrad in BmatNL"
!     do i=1,3
!     write(99,*) "",DefGrad(i,1),DefGrad(i,2),DefGrad(i,3)
!     enddo
!     write(99,*) "-------"        

      
      do i=1,6
         do j=1,24
            BmatNL(i,j) = zero
         enddo
      enddo   
      
      k=1
      do i=1,8
! B(1,*) = F(1,1)*Nk,x
         BmatNL(1,k) = Fgrad(1,1)*Bop(1,k)
! B(1,*) = F(2,1)*Nk,x
         BmatNL(1,k+1) = Fgrad(2,1)*Bop(1,k)
! B(1,*) = F(3,1)*Nk,x
         BmatNL(1,k+2) = Fgrad(3,1)*Bop(1,k)

! B(2,*) = F(1,2)*Nk,y
         BmatNL(2,k) = Fgrad(1,2)*Bop(2,k+1)
! B(2,*) = F(2,2)*Nk,y
         BmatNL(2,k+1) = Fgrad(2,2)*Bop(2,k+1)
! B(2,*) = F(3,2)*Nk,y
         BmatNL(2,k+2) = Fgrad(3,2)*Bop(2,k+1)
         
! B(3,*) = F(1,3)*Nk,z
         BmatNL(3,k) = Fgrad(1,3)*Bop(3,k+2)
! B(3,*) = F(2,3)*Nk,z
         BmatNL(3,k+1) = Fgrad(2,3)*Bop(3,k+2)
! B(3,*) = F(3,3)*Nk,z
         BmatNL(3,k+2) = Fgrad(3,3)*Bop(3,k+2)

!12
! B(4,*) = F(1,1)*Nk,y + F(1,2)*Nk,x
         BmatNL(4,k) = Fgrad(1,1)*Bop(2,k+1) + Fgrad(1,2)*Bop(1,k)
! B(4,*) = F(2,1)*Nk,y + F(2,2)*Nk,x
         BmatNL(4,k+1) = Fgrad(2,1)*Bop(2,k+1) + Fgrad(2,2)*Bop(1,k)
! B(4,*) = F(3,1)*Nk,y + F(3,2)*Nk,x
         BmatNL(4,k+2) = Fgrad(3,1)*Bop(2,k+1) + Fgrad(3,2)*Bop(1,k)


!23
! B(5,*) = F(1,2)*Nk,z + F(1,3)*Nk,y
         BmatNL(5,k) = Fgrad(1,2)*Bop(3,k+2) + Fgrad(1,3)*Bop(2,k+1)
! B(5,*) = F(2,2)*Nk,z + F(2,3)*Nk,y
         BmatNL(5,k+1) = Fgrad(2,2)*Bop(3,k+2) + Fgrad(2,3)*Bop(2,k+1)
! B(5,*) = F(3,2)*Nk,z + F(3,3)*Nk,y
         BmatNL(5,k+2) = Fgrad(3,2)*Bop(3,k+2) + Fgrad(3,3)*Bop(2,k+1)
         
!13
! B(6,*) = F(1,1)*Nk,z + F(1,3)*Nk,x
         BmatNL(6,k) = Fgrad(1,1)*Bop(3,k+2) + Fgrad(1,3)*Bop(1,k)
! B(6,*) = F(2,1)*Nk,z + F(2,3)*Nk,x
         BmatNL(6,k+1) = Fgrad(2,1)*Bop(3,k+2) + Fgrad(2,3)*Bop(1,k)
!  B(6,*) = F(3,1)*Nk,z + F(3,3)*Nk,x
         BmatNL(6,k+2) = Fgrad(3,1)*Bop(3,k+2) + Fgrad(3,3)*Bop(1,k)

         
         k=k+3

      enddo 
      

      do i=1,6
         do j=1,24
            Bop(i,j) =  BmatNL(i,j) 
         enddo
      enddo

!     write(99,*) " BmatNL"
!     do i=1,3
!     write(99,*) "",BmatNL(i,1),BmatNL(i,2),BmatNL(i,3)
!     enddo
!     write(99,*) "-------"     

      RETURN
      END
      
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------


      
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine quadrature3D8(pgp,wgp)
C----------------------------------------------------C      
C     INDICES   & AUXILIARY VARIABLES           C 
C----------------------------------------------------C 
      IMPLICIT NONE

      INTEGER i,j,k,l,m,ik,jk,kk,mm !INDICES
      
      real*8, parameter :: zero = 0.d0,quarter=0.25d0, half =0.5d0
      real*8, parameter :: one = 1.d0, two = 2.d0, three=3.d0
      
      REAL*8 pgp(3,8)
      REAL*8 wgp(8) 
      REAL*8 scalar 

!do i=1,ndim
!  do j=1,8
!    pgp(j,i) = zero 
!  enddo
!enddo
      
      do i=1,8
!  wgp (i) = zero
         wgp (i) = one
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

      pgp(1,4) =  -scalar
      pgp(2,4) =   scalar
      pgp(3,4) = - scalar

      pgp(1,5) = -scalar
      pgp(2,5) = -scalar
      pgp(3,5) =  scalar

      pgp(1,6) =  scalar
      pgp(2,6) = -scalar
      pgp(3,6) =  scalar

      pgp(1,7) =  scalar
      pgp(2,7) =  scalar
      pgp(3,7) =  scalar

      pgp(1,8) =  -scalar
      pgp(2,8) =   scalar
      pgp(3,8) =   scalar

      END

c     ----------------------------------------------------------------------------


      subroutine kEASMlocalGE3d(xi1,xi2,xi3,kaa)
      
      IMPLICIT NONE 
      
      REAL*8 kaa(6,24) 
      INTEGER i,j
      
      INTEGER rr,rs,rt,ss,st,tt

      REAL*8 xi1,xi2,xi3
      
      call pzero(kaa,6*24)
      
      rr = 1
      ss = 2
      tt = 3
      rs = 4
      st = 5 
      rt = 3
      
!incompatible modes interpolation at local setting: see Klinkel 1998     
      
      kaa(rr,1) = xi1 
      kaa(ss,2) = xi2
      kaa(tt,3) = xi3
!----------------------
      kaa(rs,4) = xi1 
      kaa(rs,5) = xi2 

      kaa(st,6) = xi2 
      kaa(st,7) = xi3 
      
      kaa(rt,8) = xi1 
      kaa(rt,9) = xi3 

!----------------------
      kaa(rs,10) = xi1*xi3 
      kaa(rs,11) = xi2*xi3 

      kaa(st,12) = xi2*xi1 
      kaa(st,13) = xi3*xi1
      
      kaa(rt,14) = xi1*xi2 
      kaa(rt,15) = xi3*xi2

!-----------------------
      kaa(rr,16) = xi1*xi2
      kaa(rr,17) = xi1*xi3 

      kaa(ss,18) = xi2*xi1
      kaa(ss,19) = xi2*xi3

      kaa(tt,20) = xi3*xi1
      kaa(tt,21) = xi3*xi2
!-----------------------
      kaa(rs,22) = xi1 *xi2
      kaa(st,23) = xi2 *xi3
      kaa(rt,24) = xi1 *xi3

!-------------------------
!     kaa(rr,25) = xi1 *xi2*xi3
!     kaa(ss,26) = xi1 *xi2*xi3
!     kaa(tt,27) = xi1 *xi2*xi3
!     kaa(rs,28) = xi1 *xi2*xi3
!     kaa(st,29) = xi1 *xi2*xi3
!     kaa(rt,30) = xi1 *xi2*xi3

      end



      subroutine kEASTransfGE3D(Weas,xjacR,xjacR0,kaa,
     *     detr0,detr)
      
      IMPLICIT NONE
      
! number definition
      
      REAL*8 xjacR(3,3)
      REAL*8 xjacR0(3,3)
      
      REAL*8 T0mat(6,6)
      REAL*8 kaa(6,24)
      REAL*8 Weas(6,24)

      
      REAL*8 j11,j12,j13
      REAL*8 j21,j22,j23
      REAL*8 j31,j32,j33

      
      REAL*8 factor0,detr0,detr 
      
      INTEGER i,j,k,l
      
      REAL*8 TOL 
      
      TOL = 1.0d-15
      
      do i=1,6
         do j=1,6
            T0mat(i,j)=0.d0
         enddo
      enddo
      
      do i=1,6
         do j=1,24
            Weas(i,j) = 0.d0
         enddo
      enddo  

      factor0= detr0/detr       !relation between both determinants

      j11 = xjacR0(1,1)
      j12 = xjacR0(1,2)
      j13 = xjacR0(1,3)


      j21 = xjacR0(2,1)
      j22 = xjacR0(2,2)
      j23 = xjacR0(2,3)

      j31 = xjacR0(3,1)
      j32 = xjacR0(3,2)
      j33 = xjacR0(3,3)
      
      
      T0mat(1,1) = j11*j11
      T0mat(1,2) = j21*j21
      T0mat(1,3) = j31*j31
      T0mat(1,4) = 2.d0*j11*j21
      T0mat(1,5) = 2.d0*j21*j31
      T0mat(1,6) = 2.d0*j11*j31


      T0mat(2,1) = j12*j12
      T0mat(2,2) = j22*j22
      T0mat(2,3) = j32*j32
      T0mat(2,4) = 2.d0*j12*j22
      T0mat(2,5) = 2.d0*j22*j32
      T0mat(2,6) = 2.d0*j12*j32

      T0mat(3,1) = j13*j13
      T0mat(3,2) = j23*j23
      T0mat(3,3) = j33*j33
      T0mat(3,4) = 2.d0*j13*j23
      T0mat(3,5) = 2.d0*j23*j33
      T0mat(3,6) = 2.d0*j13*j33


      T0mat(4,1) = j11*j12
      T0mat(4,2) = j21*j22
      T0mat(4,3) = j31*j32
      T0mat(4,4) = j11*j22 + j21*j12
      T0mat(4,5) = j21*j32 + j31*j22
      T0mat(4,6) = j11*j32 + j31*j12

      T0mat(5,1) = j12*j13
      T0mat(5,2) = j23*j23
      T0mat(5,3) = j32*j33
      T0mat(5,4) = j12*j23 + j22*j13
      T0mat(5,5) = j22*j33 + j32*j23
      T0mat(5,6) = j12*j33 + j32*j13
      

      T0mat(6,1) = j11*j13
      T0mat(6,2) = j21*j23
      T0mat(6,3) = j31*j33
      T0mat(6,4) = j11*j23 + j21*j13
      T0mat(6,5) = j21*j33 + j31*j23
      T0mat(6,6) = j11*j33 + j31*j13



      do i=1,6
         do j=1,6
            
            if (abs(T0mat(i,j)).LT.TOL) then
               T0mat(i,j) =  0.d0
            endif  
            
         enddo
      enddo      
      
      
      
      call getMatrixMultiply(Weas,6,24,T0mat,6,6,
     *     kaa,6,24)  
      

      do i=1,6
         do j=1,24
            if (abs(Weas(i,j)).LE.TOL) then
               Weas(i,j) =  0.d0
            endif  
            
         enddo
      enddo  


      RETURN
      END   
!---------------------------------------------------------
