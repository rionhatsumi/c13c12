C
C      ESTE PROGRAMA CALCULA PENETRACAO DE BARREIRA COM NOSSO POTENCIAL GLOBAL
C      E com FOLDING PARA INTERACAO COULOMBIANA
C
      IMPLICIT REAL*8(A-H,O-Z)
      dimension vf(50000),vc(50000),Vef(50000,0:200),VREAL(50000),
     *IRL(0:200),RBL(0:200),VBL(0:200),HWL(0:200),TLWKB(0:200),
     *hwwkb(0:200),hwef(0:200),tlef(0:200),ff(0:1000),ffc(0:1000),
     *vblef(0:200),ro1(0:1000),roc1(0:1000),ro2(0:1000),roc2(0:1000),
     *TLW(0:200),TLHW(0:200)
      PI=4.d0*datan(1.d0)     
      a1=0.56
      ac1=0.53
      a2=0.56
      ac2=0.53
      dr=0.005
      dq=.02
      ircore=0
      nqmax=1000

      write(*,9518)
 9518 format('Put 0 for [mbarn] or 1 for [barn]')
      read(*,*) iopt1

      write(*,17)
 17   format('Put 1 to include "exp(gE)" in the S-factor equation;',1x,/
     *'Any other number makes exp(gE)=0')
      read(*,*) ioptg
      g=0.d0

      open(5,file='sfactor.dat',status='old')
      open(1,file='projectile.de1',status='unknown')
      open(2,file='target.de2',status='unknown')
      open(6,file='sfactor.out',status='unknown')
      open(4,file='fusion.out',status='unknown')
      open(3,file='verify.out',status='unknown')
c
c      Entrada das funcoes para barreiras efetivas
c
      READ(5,*) FMP,FMA,ZP,ZA
      FMI=FMP*FMA/(FMP+FMA)
      estvb=za*zp/(fma**0.333333+fmp**0.333333)
      estrb=1.4*(fma**0.333333+fmp**0.333333)
      read(5,*) iopde1,iopde2
      read(5,*) eemin,eemax
      npts=abs((eemax-eemin))/0.1
      write(3,*) npts
c
c     particulas identicas ou nao
c
      if(zp.eq.za.and.fmp.eq.fma) then
      iopart=1
      ldelta=2
      else
      iopart=0
      ldelta=1
      endif
      RC1=1.76*zp**.33333-.96
      RC2=1.76*za**.33333-.96
      RA1=1.31*fmp**.33333-.84
      RA2=1.31*fma**.33333-.84
c
c      Densities:
c
c      i) Projectile
c      
      if(iopde1.eq.0) then
c
c      density normalization
c
      rmax1=ra1+14.*a1
      nmax1=idint(rmax1/.02)
      if(nmax1.gt.1000) then
      write(3,900)
 900  format(1x,'Number of points for the density 1 is too large')
      stop
      endif
      soma=0.
      somac=0.
      do 1 i=0,nmax1
      r=dfloat(i)*.02+.01
      rho=1./(1.+dexp((r-ra1)/a1))
      rhoc=1./(1.+dexp((r-rc1)/ac1))
      soma=soma+r**2*rho
      somac=somac+r**2*rhoc
 1    continue
      soma=soma*4.*pi*.02
      rho01=fmp/soma
      somac=somac*4.*pi*.02
      rho0c1=zp/somac
c
c      density calculation
c
      do 2 i=0,nmax1
      r=dfloat(i)*.02+.01
      ro1(i)=rho01/(1.+dexp((r-ra1)/a1))
      roc1(i)=rho0c1/(1.+dexp((r-rc1)/ac1))
 2    continue
      endif
c
c      read external file
c
      if(iopde1.eq.1) then
      i=0
      somalu=0.
      somaluc=0.
 70   read(1,*,end=80) rnao,ro1(i),roc1(i)
      r=dfloat(i)*.02+.01      
      somalu=somalu+r**2*ro1(i)
      somaluc=somaluc+r**2*roc1(i)
      i=i+1
      go to 70
 80   nmax1=i-1
      if(nmax1.gt.1000) then
      write(3,900)
      stop
      endif
      somalu=somalu*4.*pi*0.02
      somaluc=somaluc*4.*pi*0.02
      if(dabs(somalu-fmp)/fmp.gt.0.01.or.dabs(somaluc-zp)/zp.gt..01) 
     *write(3,1050) somalu,somaluc
 1050 format(//,' warning: density 1 is not properly normalized, A= ',
     *f6.2,' Z= ',f6.2,//)
      endif
c
c      ii) Target
c
      if(iopde2.eq.0) then
      rmax2=ra2+14.*a2
      nmax2=idint(rmax2/.02)
      if(nmax2.gt.1000) then
      write(3,901)
 901  format(1x,'Number of points for the density 2 is too large')
      stop
      endif
      soma=0.
      somac=0.
      do 90 i=0,nmax2
      r=dfloat(i)*.02+.01
      rho=1./(1.+dexp((r-ra2)/a2))
      rhoc=1./(1.+dexp((r-rc2)/ac2))
      soma=soma+r**2*rho
      somac=somac+r**2*rhoc
 90   continue
      soma=soma*4.*pi*.02
      rho02=fma/soma
      somac=somac*4.*pi*.02
      rho0c2=za/somac
      do 103 i=0,nmax2
      r=dfloat(i)*.02+.01
      ro2(i)=rho02/(1.+dexp((r-ra2)/a2))
      roc2(i)=rho0c2/(1.+dexp((r-rc2)/ac2))
 103  continue
      endif

      if(iopde2.eq.1) then
      i=0
      somalu=0.
      somaluc=0.
 150  read(2,*,end=160) rnao,ro2(i),roc2(i)
      r=dfloat(i)*.02+.01      
      somalu=somalu+r**2*ro2(i)
      somaluc=somaluc+r**2*roc2(i)
      i=i+1
      go to 150
 160  nmax2=i-1
      if(nmax2.gt.1000) then
      write(3,901)
      stop
      endif
      somalu=somalu*4.*pi*0.02
      somaluc=somaluc*4.*pi*0.02
      if(dabs(somalu-fma)/fma.gt.0.01.or.dabs(somaluc-za)/za.gt..01) 
     *write(3,1061) somalu,somaluc
 1061 format(//,' warning: density 2 is not properly normalized, A= ',
     *f6.2,' Z= ',f6.2,//)
      endif
c
c      transformadas das densidades
c
      do 180 iq=0,1000
      q=dfloat(iq)*.02+.01
      soma1=0.
      somac1=0.
      do 190 i=0,nmax1
      r=dfloat(i)*.02+.01
      soma1=soma1+ro1(i)*r*dsin(q*r)
      somac1=somac1+roc1(i)*r*dsin(q*r)
 190  continue
      f1=4.*pi*soma1*.02/q
      fc1=4.*pi*somac1*.02/q
      soma2=0.
      somac2=0.
      do 220 i=0,nmax2
      r=dfloat(i)*.02+.01
      soma2=soma2+ro2(i)*r*dsin(q*r)
      somac2=somac2+roc2(i)*r*dsin(q*r)
 220  continue
      f2=4.*pi*soma2*.02/q
      fc2=4.*pi*somac2*.02/q
c
c      Multiplicacao das transformadas
c
      ff(iq)=f1*f2
      ffc(iq)=fc1*fc2
 180  continue
C
C      NAO ESTA PREPARADO PARA MAIS DE 50000 PASSOS DE INTEGRACAO
C
      rmax=ra1+ra2+35.
      retorno=1.43986*za*zp/eemin
      if(rmax.lt.retorno) rmax=retorno+20.
      NPR=IDINT(RMAX/DR)+2
      IF(NPR.GT.50000) then
      write(3,1000) nsist,fma,fmp,za,zp,npr
 1000 format(1x,i3,4(1x,f5.1),' NPR=',i6)
      npr=50000
      endif
c
c      calculo do potencial Coulombiano, real e imaginario
c
      DO 9 IR=1,NPR
      r=dfloat(ir)*dr
c
c     folding e Coulomb
c
      vf(ir)=0.
      vc(ir)=0.
      do 6 iq=0,nqmax
      q=dfloat(iq)*dq+dq/2.
      vf(ir)=vf(ir)+q*dsin(q*r)*ff(iq)
      vc(ir)=vc(ir)+dsin(q*r)*ffc(iq)/q
 6    continue
      vf(ir)=-456.*vf(ir)*dq/(2.*pi**2*r)
      vc(ir)=((1.43986*2.)/(pi*r))*vc(ir)*dq
 9    continue
C
C      Entrando com a energia
C
      ecm=eemin-0.1
      do 100 j=1,npts
      ecm=ecm+0.1
      if(ecm.le.estvb) then
      sigmae=0.01*(ecm/estvb)
      else
      sigmae=10.*pi*estrb**2*(1.-estvb/ecm)
      endif
      FK=DSQRT(ECM*FMI/20.913067)
c
c      estimativa de Lmax
c
      lmax=idint(dsqrt(sigmae*fk**2/(pi*10.)))+50
      if(iopart.eq.1) then
      lmax=(lmax/2)
      lmax=lmax*2
      endif
c
c      Nao esta preparado para mais de 200 ondas parciais
c
      if(lmax.gt.200) then
      write(3,2000) nsist,fma,fmp,za,zp,ecm,lmax
 2000 format(1x,i3,4(1x,f5.1),' Ecm=',f7.2,' Lmax=',i3)
      lmax=200
      endif
c
c     potencial local equivalente - por iteracoes
c
      do 91 ir=1,npr
      r=dfloat(ir)*dr
      vn1=vf(ir)
      ve2=2.*(ecm-vc(ir)-vf(ir))/(931.5*fmi)
      vn2=vf(ir)*dexp(-4.*ve2)
c
c     teste de convergencia
c
      itest=0
 7    itest=itest+1
      ve2=2.*(ecm-vc(ir)-vn2)/(931.5*fmi)
      vn1=vn2
      vn2=vf(ir)*dexp(-4.*ve2)
      test=2.*dabs(vn2-vn1)/(dabs(vn1)+dabs(vn2))
      if(test.gt..005) go to 7
c
c     potencial local equivalente
c
      vreal(ir)=vn2
C
C      CALCULO DE Vef
C
      RO=R*FK
      DO 8 L=0,LMAX,ldelta
      Vef(IR,L)=Ecm*DFLOAT(L*(L+1))/(RO**2)+VC(ir)+vreal(ir)
 8    CONTINUE
      if(j.eq.1) write(3,*) r,vc(ir),-vreal(ir),vef(ir,0)
 91   continue
c
c      calculo do raio da barreira para cada L
c
      lmaxsalva=lmax
      do 10 l=0,lmax,ldelta
      ir=npr
 20   ir=ir-1
      If(vef(ir,l).gt.vef(ir+1,l).and.ir.ne.5) go to 20
      irl(l)=ir
c
c      Se nao encontrou a barreira, adota outro lmax
c
      if(irl(l).eq.5.and.lmaxsalva.eq.lmax) then
      if(iopart.eq.0) lmaxsalva=l-1
      if(iopart.eq.1) lmaxsalva=l-2
      go to 93
      endif
c
c      Parametros das Barreiras
c
      RBL(l)=dfloat(irl(l))*dr
      VBL(l)=vef(irl(l),l)
c
c     segunda derivada de VEF
c
      dv=(-vef(ir-2,l)+16.*vef(ir-1,l)-30.*vbl(l)+16.*
     *vef(ir+1,l)-vef(ir+2,l))/(12.*dr**2)
      hwl(l)=dsqrt(dabs((197.3**2/(fmi*931.5))*dv))
 93   continue
 10   continue
      if(ecm.ge.(0.98*vbl(0)).and.ecm.le.(1.02*vbl(0))) write(3,*) ecm,
     *vbl(0),hwl(0)
      lmax=lmaxsalva
c
c      calculo dos demais HW
c
      do 30 l=0,lmax,ldelta
      if(ecm.ge.vbl(l)) then
      hwwkb(l)=hwl(l)
      go to 29
      endif
c
c      calculo dos pontos de retorno
c
      imax=npr
      do 25 ir=npr,irl(l),-1
      if((vef(ir,l)-ecm).lt.0.) imax=ir-1
 25   continue
 35   ir=ir-1
      if((vef(ir,l)-ecm).gt.0..and.ir.ge.2) go to 35 
      imin=ir+1
c
c     nuclear radius for S-factor
c
      ir=imin+10
 36   ir=ir-1
      if(vef(ir,0).ge.0.) go to 36
      if(ircore.le.0) ircore=ir+1

      if((imax-imin).lt.20) then
      hwwkb(l)=hwl(l)
      go to 29
      endif
      soma_jwkb = 0.
      do 45 ir=imin,imax
      soma_jwkb=soma_jwkb+dsqrt(0.191434*fmi*(vef(ir,l)-ecm))
 45   continue
      hwwkb(l)=2.*pi*(vbl(l)-ecm)/(soma_jwkb*dr)
 29   continue
c
c     Effective Parameters
c
      if(fmi.gt.8.) then
      vblef(l)=vbl(l)
      hwef(l)=hwwkb(l)*(1.+0.1*(fmi-8.))
      else
      vblef(l)=vbl(l)
      hwef(l)=hwwkb(l)
      endif
 30   continue
c
c      calculo dos TL
c
      do 40 l=0,lmax,ldelta
      ajuda=vblef(0)+dfloat(l*(l+1))*197.**2/(2.*fmi*931.*rbl(0)**2)
      x=2.*pi*(ajuda-ecm)/hwef(0)
c      if(x.ge.-30..and.x.le.30.) tlw(l)=1./(1.+dexp(x))
      if(x.ge.-30.) tlw(l)=1./(1.+dexp(x))
      if(x.lt.-30.) tlw(l)=0.
      x=2.*pi*(vbl(l)-ecm)/hwl(l)
      if(x.ge.-30..and.x.le.30.) tlhw(l)=1./(1.+dexp(x))
      if(x.lt.-30.) tlhw(l)=1.
      if(x.gt.30.) tlhw(l)=0.
      x=2.*pi*(vbl(l)-ecm)/hwwkb(l)
c      if(x.ge.-30..and.x.le.100.) tlwkb(l)=1./(1.+dexp(x))
      if(x.ge.-30.) tlwkb(l)=1./(1.+dexp(x))
      if(x.lt.-30.) tlwkb(l)=1.
c      if(x.gt.100.) tlwkb(l)=0.
      x=2.*pi*(vblef(l)-ecm)/hwef(l)
c      if(x.ge.-30..and.x.le.100.) tlef(l)=1./(1.+dexp(x))
      if(x.ge.-30.) tlef(l)=1./(1.+dexp(x))
      if(x.lt.-30.) tlef(l)=1.
c      if(x.gt.100.) tlef(l)=0.
 40   continue
c
c      calculo das secoes de choque
c
      Somaw=0.
      somahw=0.
      somawkb=0.
      somaef=0.
      do 50 l=0,lmax,ldelta
      somaw=somaw+(2.*dfloat(l)+1)*tlw(l)
      somahw=somahw+(2.*dfloat(l)+1)*tlhw(l)
      somawkb=somawkb+(2.*dfloat(l)+1)*tlwkb(l)
      somaef=somaef+(2.*dfloat(l)+1)*tlef(l)
 50   continue
      sw=(pi/fk**2)*somaw*10.
      shw=(pi/fk**2)*somahw*10.
      swkb=(pi/fk**2)*somawkb*10.
      sef=(pi/fk**2)*somaef*10.
      if(iopart.eq.1) then
      sw=sw*2.
      shw=shw*2.
      swkb=swkb*2.
      sef=sef*2.
      endif
      if(iopt1.eq.0) write(4,1234) ecm,swkb,sef
      if(iopt1.eq.1) write(4,1234) ecm,swkb/1000.,sef/1000.
c
c     S-Factor
c
      rcore=dfloat(ircore)*dr
      eta=(0.1574*zp*za)/(dsqrt(ecm/fmi))
c      sec=1.0495*dsqrt(zp*za*fmi*rcore)
      sec=0. 
      if(ioptg.eq.1) g=0.1215*dsqrt(fmi*rcore**3/(zp*za))
      eleva=2.*pi*eta-sec+g*ecm
      sfactor=swkb*ecm*dexp(eleva)
      sfactoref=sef*ecm*dexp(eleva)
c      SFACTOR=sfactor/(1.e16)
      if(iopt1.eq.0) write(6,1234) ecm,sfactor,sfactoref,swkb,sef
      if(iopt1.eq.1) write(6,1234) ecm,sfactor/1000.,sfactoref/1000.,
     & swkb,sef
 1234 format(1x,f7.3,4(1x,e14.7e3))
 100  continue
      write(3,*) rcore,vf(irl(0))
      close(5)
      close(1)
      close(2)
      close(6)
      close(4)
      close(3)
      stop
      end
