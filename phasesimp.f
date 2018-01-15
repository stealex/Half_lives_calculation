


       implicit double precision (a-h,o-z)
           dimension e(99),greal(99),gimag(99),freal(99),fimag(99)
           dimension greal1(99),gimag1(99),freal1(99),fimag1(99)
      DIMENSION X1x(200),B(200),C(200),D(200),Z(32),G(32)
      DIMENSION  A1A(200),DG(200)
      common/val/ca,ce(99),cgreal(99),cgimag(99),cfreal(99),cfimag(99),
     c cgreal1(99),cgimag1(99),cfreal1(99),cfimag1(99),nner
      common/en/wmax,amc2,raza
      external ai1
      amc2=0.51099891
cc hbar=6.58211928 10-16 eV s=6.58211928 10-22 MeV s=6.58211928 10-22/31622400 MeV Year
      hbar=6.58211928e-22/31622400
         vvv=3600*24*366

c daca Energia este in MeV si distanta in Fermi
c atunci p=radical(E^2-(mc^2)^)/c^2
c pR=pR/(2.9979245*6.58211928*10) in unitati hbar unde p este in MeV/c^2 si R in fm 
c deoarece c=2.9979246 10^23 fm/s
c          hbar=6.58211928 10^-22 MeV s



      qvalue=0.282+amc2
      qvalue=5.0+amc2
       print*,'valoarea lui E(beta)'
       read(*,*)veb
       qvalue=veb+amc2
c      a=20
       open(82,file='asiz.dat',status='old')
       read(82,*)a,z0
       close(82)
       raza=1.2*a**(1./3.)


      ei=0.0
       call gf(a,e,greal,gimag,freal,fimag,greal1,gimag1,freal1,fimag1)
       ca=a
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      numar de valori in energii ner (mai mic de 99)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ner=42
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       nner=ner
       do l=1,ner
       ce(l)=e(l)
       cgreal(l)=greal(l)
       cgimag(l)=gimag(l)
       cfreal(l)=freal(l)
       cfimag(l)=fimag(l)
       cgreal1(l)=greal1(l)
       cgimag1(l)=gimag1(l)
       cfreal1(l)=freal1(l)
       cfimag1(l)=fimag1(l)
       enddo
       do l=1,ner
      e1=ce(l)
      call fiuri(e1,f110)
      write(65,*)e1,f110,greal(l),gimag(l)
       enddo
c in aceasta versiune am zece energii n=11 intre 0.52 si 6
       alims=qvalue/amc2
       wmax=alims
       alimi=1.d0

       call dqg32(alimi,alims,ai1,ggg)
        print*,'ggg',ggg,'log ggg',dlog10(ggg)
       end



      double precision function ai1(w)
      implicit double precision (a-h,o-z)
      common/val/ca,ce(99),cgreal(99),cgimag(99),cfreal(99),cfimag(99),
     c cgreal1(99),cgimag1(99),cfreal1(99),cfimag1(99),nner
      common/en/wmax,amc2,raza
      common/atilda/atilda2
      p=dsqrt(w**2-1.d0)
      enn1=w*amc2
      call fiuri(enn1,f110)
       ai1=f110*p*w*(wmax-w)**2 !/2. /p**2
c      print*,'ai3,f110',ai3,f110
      return
      end



      double precision function ai2(w)
      implicit double precision (a-h,o-z)
      common/val/ca,ce(99),cgreal(99),cgimag(99),cfreal(99),cfimag(99),
     c cgreal1(99),cgimag1(99),cfreal1(99),cfimag1(99),nner
      common/en/wmax,amc2,raza
      common/atilda/atilda2
c pR=pR/(2.9979245*6.58211928*10) in unitati hbar unde p este in MeV/c^2 si R in fm 
      pen=dsqrt((w*amc2)**2-amc2**2)
      pR=pen*raza/(2.9979245*6.58211928*10)
      enn1=w*amc2
      call fiuri(enn1,f110)
      p=dsqrt(w**2-1.d0)
c nu calculez g-2 si f2
      alambda2=f110/2./p**2*(3./pR)**2
      almabda1=f110/2./p**2
      ai2=1./12.*p*w*(wmax-w)**2*(alambda1*(wmax-w)**2+alambda2*p**2)
c      print*,'ai3,f110',ai3,f110
      return
      end



      subroutine fiuri(e1,f110)
       implicit double precision (a-h,o-z)
           dimension e(50),greal(50),gimag(50),freal(50),fimag(50)
           dimension greal1(50),gimag1(50),freal1(50),fimag1(50)
      DIMENSION X1x(200),B(200),C(200),D(200),Z(32),G(32)
      DIMENSION  A1A(200),DG(200)
      common/val/ca,ce(99),cgreal(99),cgimag(99),cfreal(99),cfimag(99),
     c cgreal1(99),cgimag1(99),cfreal1(99),cfimag1(99),nner
      common/en/wmax,amc2,raza
       a=ca
       n=nner
       do l=1,n
       e(l)=ce(l)
       greal(l)=cgreal(l)
       gimag(l)=cgimag(l)
       freal(l)=cfreal(l)
       fimag(l)=cfimag(l)
       greal1(l)=cgreal1(l)
       gimag1(l)=cgimag1(l)
       freal1(l)=cfreal1(l)
       fimag1(l)=cfimag1(l)
       enddo


      do i3=1,n
      x1x(i3)=e(i3)
      enddo


      m=1
      j=1
      k=1

          z(1)=e1
          do i3=1,n
          a1a(i3)=greal(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          gr=g(1)
          do i3=1,n
          a1a(i3)=gimag(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          gi=g(1)
          do i3=1,n
          a1a(i3)=freal(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          fr=g(1)
          do i3=1,n
          a1a(i3)=fimag(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          fi=g(1)
          do i3=1,n
          a1a(i3)=greal1(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          grp=g(1)
          do i3=1,n
          a1a(i3)=gimag1(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          gip=g(1)
          do i3=1,n
          a1a(i3)=freal1(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          frp=g(1)
          do i3=1,n
          a1a(i3)=fimag1(i3)
          enddo
          call XSPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
          fip=g(1)


       gminus=gr**2+gi**2
       fplus=frp**2+fip**2
      f110=gminus+fplus



       return
       end





      subroutine gf(a,e,greal,gimag,freal,fimag,
     c                    greal1,gimag1,freal1,fimag1)

       implicit double precision (a-h,o-z)
       dimension xx(2),pp(2),qq(2)
           dimension e(99),greal(99),gimag(99),freal(99),fimag(99)
           dimension g2(99),f2(99)
           dimension greal1(99),gimag1(99),freal1(99),fimag1(99)
           dimension g21(99),f21(99)




             d=0
             ds=0



       r=1.2*a**(.33333333)
       open(1,file='minus1/e0.511/fort.1',status='old')
       open(11,file='minus1/e0.511/fort.11',status='old')
       e( 1)=0.51100
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  10
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 10    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 1)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 1)=dsin(-d-ds)*ppp*ano*ano1
       freal( 1)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 1)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.511/fort.1',status='old')
       open(11,file='plus1/e0.511/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  11
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 11    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 1)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 1)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 1)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 1)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5111/fort.1',status='old')
       open(11,file='minus1/e0.5111/fort.11',status='old')
       e( 2)=0.5111
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  20
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 20    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 2)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 2)=dsin(-d-ds)*ppp*ano*ano1
       freal( 2)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 2)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5111/fort.1',status='old')
       open(11,file='plus1/e0.5111/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  21
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 21    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 2)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 2)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 2)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 2)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5112/fort.1',status='old')
       open(11,file='minus1/e0.5112/fort.11',status='old')
       e( 3)=0.5112
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  30
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 30    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 3)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 3)=dsin(-d-ds)*ppp*ano*ano1
       freal( 3)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 3)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5112/fort.1',status='old')
       open(11,file='plus1/e0.5112/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  31
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 31    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 3)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 3)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 3)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 3)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5113/fort.1',status='old')
       open(11,file='minus1/e0.5113/fort.11',status='old')
       e( 4)=0.5113
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  40
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 40    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 4)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 4)=dsin(-d-ds)*ppp*ano*ano1
       freal( 4)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 4)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5113/fort.1',status='old')
       open(11,file='plus1/e0.5113/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  41
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 41    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 4)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 4)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 4)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 4)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5115/fort.1',status='old')
       open(11,file='minus1/e0.5115/fort.11',status='old')
       e( 5)=0.5115
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  50
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 50    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 5)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 5)=dsin(-d-ds)*ppp*ano*ano1
       freal( 5)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 5)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5115/fort.1',status='old')
       open(11,file='plus1/e0.5115/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  51
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 51    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 5)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 5)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 5)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 5)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5117/fort.1',status='old')
       open(11,file='minus1/e0.5117/fort.11',status='old')
       e( 6)=0.5117
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  60
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 60    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 6)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 6)=dsin(-d-ds)*ppp*ano*ano1
       freal( 6)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 6)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5117/fort.1',status='old')
       open(11,file='plus1/e0.5117/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  61
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 61    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 6)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 6)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 6)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 6)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.512/fort.1',status='old')
       open(11,file='minus1/e0.512/fort.11',status='old')
       e( 7)=0.512
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  70
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 70    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 7)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 7)=dsin(-d-ds)*ppp*ano*ano1
       freal( 7)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 7)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.512/fort.1',status='old')
       open(11,file='plus1/e0.512/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  71
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 71    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 7)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 7)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 7)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 7)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5125/fort.1',status='old')
       open(11,file='minus1/e0.5125/fort.11',status='old')
       e( 8)=0.5125
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  80
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 80    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 8)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 8)=dsin(-d-ds)*ppp*ano*ano1
       freal( 8)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 8)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5125/fort.1',status='old')
       open(11,file='plus1/e0.5125/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  81
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 81    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 8)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 8)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 8)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 8)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5130/fort.1',status='old')
       open(11,file='minus1/e0.5130/fort.11',status='old')
       e( 9)=0.5130
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  90
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 90    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 9)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 9)=dsin(-d-ds)*ppp*ano*ano1
       freal( 9)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 9)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5130/fort.1',status='old')
       open(11,file='plus1/e0.5130/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  91
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 91    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 9)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 9)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 9)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 9)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5135/fort.1',status='old')
       open(11,file='minus1/e0.5135/fort.11',status='old')
       e(10)=0.5135 
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  100
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 100    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 10)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 10)=dsin(-d-ds)*ppp*ano*ano1
       freal( 10)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 10)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5135/fort.1',status='old')
       open(11,file='plus1/e0.5135/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  101
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 101    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 10)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 10)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 10)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 10)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.5140/fort.1',status='old')
       open(11,file='minus1/e0.5140/fort.11',status='old')
       e(11)=0.5140
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  110
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 110    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 11)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 11)=dsin(-d-ds)*ppp*ano*ano1
       freal( 11)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 11)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.5140/fort.1',status='old')
       open(11,file='plus1/e0.5140/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  111
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 111    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 11)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 11)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 11)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 11)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.515/fort.1',status='old')
       open(11,file='minus1/e0.515/fort.11',status='old')
       e(12)=0.515
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  120
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 120    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 12)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 12)=dsin(-d-ds)*ppp*ano*ano1
       freal( 12)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 12)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.515/fort.1',status='old')
       open(11,file='plus1/e0.515/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  121
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 121    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 12)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 12)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 12)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 12)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.517/fort.1',status='old')
       open(11,file='minus1/e0.517/fort.11',status='old')
       e(13)=0.517
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  130
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 130    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 13)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 13)=dsin(-d-ds)*ppp*ano*ano1
       freal( 13)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 13)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.517/fort.1',status='old')
       open(11,file='plus1/e0.517/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  131
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 131    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 13)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 13)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 13)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 13)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.52/fort.1',status='old')
       open(11,file='minus1/e0.52/fort.11',status='old')
       e(14)=0.52
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  140
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 140    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 14)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 14)=dsin(-d-ds)*ppp*ano*ano1
       freal( 14)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 14)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.52/fort.1',status='old')
       open(11,file='plus1/e0.52/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  141
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 141    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 14)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 14)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 14)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 14)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.525/fort.1',status='old')
       open(11,file='minus1/e0.525/fort.11',status='old')
       e(15)=0.525
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  150
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 150    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 15)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 15)=dsin(-d-ds)*ppp*ano*ano1
       freal( 15)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 15)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.525/fort.1',status='old')
       open(11,file='plus1/e0.525/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  151
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 151    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 15)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 15)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 15)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 15)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.53/fort.1',status='old')
       open(11,file='minus1/e0.53/fort.11',status='old')
       e(16)=0.53
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  160
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 160    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 16)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 16)=dsin(-d-ds)*ppp*ano*ano1
       freal( 16)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 16)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.53/fort.1',status='old')
       open(11,file='plus1/e0.53/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  161
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 161    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 16)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 16)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 16)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 16)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.535/fort.1',status='old')
       open(11,file='minus1/e0.535/fort.11',status='old')
       e(17)=0.535
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  170
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 170    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 17)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 17)=dsin(-d-ds)*ppp*ano*ano1
       freal( 17)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 17)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.535/fort.1',status='old')
       open(11,file='plus1/e0.535/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  171
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 171    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 17)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 17)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 17)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 17)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.54/fort.1',status='old')
       open(11,file='minus1/e0.54/fort.11',status='old')
       e(18)=0.54
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  180
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 180    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 18)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 18)=dsin(-d-ds)*ppp*ano*ano1
       freal( 18)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 18)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.54/fort.1',status='old')
       open(11,file='plus1/e0.54/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  181
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 181    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 18)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 18)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 18)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 18)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.55/fort.1',status='old')
       open(11,file='minus1/e0.55/fort.11',status='old')
       e(19)=0.55
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  190
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 190    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 19)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 19)=dsin(-d-ds)*ppp*ano*ano1
       freal( 19)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 19)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.55/fort.1',status='old')
       open(11,file='plus1/e0.55/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  191
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 191    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 19)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 19)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 19)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 19)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.575/fort.1',status='old')
       open(11,file='minus1/e0.575/fort.11',status='old')
       e(20)=0.575
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  200
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 200    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 20)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 20)=dsin(-d-ds)*ppp*ano*ano1
       freal( 20)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 20)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.575/fort.1',status='old')
       open(11,file='plus1/e0.575/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  201
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 201    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 20)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 20)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 20)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 20)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.6/fort.1',status='old')
       open(11,file='minus1/e0.6/fort.11',status='old')
       e(21)=0.6
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  210
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 210    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 21)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 21)=dsin(-d-ds)*ppp*ano*ano1
       freal( 21)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 21)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.6/fort.1',status='old')
       open(11,file='plus1/e0.6/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  211
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 211    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 21)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 21)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 21)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 21)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.65/fort.1',status='old')
       open(11,file='minus1/e0.65/fort.11',status='old')
       e(22)=0.65
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  220
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 220    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 22)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 22)=dsin(-d-ds)*ppp*ano*ano1
       freal( 22)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 22)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.65/fort.1',status='old')
       open(11,file='plus1/e0.65/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  221
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 221    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 22)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 22)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 22)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 22)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.7/fort.1',status='old')
       open(11,file='minus1/e0.7/fort.11',status='old')
       e(23)=0.7
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  230
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 230    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 23)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 23)=dsin(-d-ds)*ppp*ano*ano1
       freal( 23)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 23)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.7/fort.1',status='old')
       open(11,file='plus1/e0.7/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  231
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 231    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 23)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 23)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 23)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 23)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.8/fort.1',status='old')
       open(11,file='minus1/e0.8/fort.11',status='old')
       e(24)=0.8
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  240
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 240    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 24)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 24)=dsin(-d-ds)*ppp*ano*ano1
       freal( 24)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 24)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.8/fort.1',status='old')
       open(11,file='plus1/e0.8/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  241
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 241    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 24)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 24)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 24)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 24)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e0.9/fort.1',status='old')
       open(11,file='minus1/e0.9/fort.11',status='old')
       e(25)=0.9
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  250
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 250    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 25)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 25)=dsin(-d-ds)*ppp*ano*ano1
       freal( 25)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 25)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e0.9/fort.1',status='old')
       open(11,file='plus1/e0.9/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  251
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 251    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 25)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 25)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 25)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 25)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e1/fort.1',status='old')
       open(11,file='minus1/e1/fort.11',status='old')
       e(26)=1.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  260
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 260    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 26)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 26)=dsin(-d-ds)*ppp*ano*ano1
       freal( 26)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 26)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e1/fort.1',status='old')
       open(11,file='plus1/e1/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  261
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 261    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 26)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 26)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 26)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 26)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e1.1/fort.1',status='old')
       open(11,file='minus1/e1.1/fort.11',status='old')
       e(27)=1.1
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  270
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 270    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 27)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 27)=dsin(-d-ds)*ppp*ano*ano1
       freal( 27)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 27)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e1.1/fort.1',status='old')
       open(11,file='plus1/e1.1/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  271
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 271    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 27)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 27)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 27)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 27)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e1.25/fort.1',status='old')
       open(11,file='minus1/e1.25/fort.11',status='old')
       e(28)=1.25
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  280
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 280    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 28)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 28)=dsin(-d-ds)*ppp*ano*ano1
       freal( 28)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 28)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e1.25/fort.1',status='old')
       open(11,file='plus1/e1.25/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  281
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 281    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 28)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 28)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 28)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 28)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e1.5/fort.1',status='old')
       open(11,file='minus1/e1.5/fort.11',status='old')
       e(29)=1.5
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  290
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 290    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 29)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 29)=dsin(-d-ds)*ppp*ano*ano1
       freal( 29)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 29)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e1.5/fort.1',status='old')
       open(11,file='plus1/e1.5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  291
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 291    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 29)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 29)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 29)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 29)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e1.75/fort.1',status='old')
       open(11,file='minus1/e1.75/fort.11',status='old')
       e(30)=1.75
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  300
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 300    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 30)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 30)=dsin(-d-ds)*ppp*ano*ano1
       freal( 30)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 30)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e1.75/fort.1',status='old')
       open(11,file='plus1/e1.75/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  301
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 301    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 30)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 30)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 30)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 30)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e2/fort.1',status='old')
       open(11,file='minus1/e2/fort.11',status='old')
       e(31)=2.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  310
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 310    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 31)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 31)=dsin(-d-ds)*ppp*ano*ano1
       freal( 31)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 31)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e2/fort.1',status='old')
       open(11,file='plus1/e2/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  311
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 311    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 31)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 31)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 31)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 31)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e2.1/fort.1',status='old')
       open(11,file='minus1/e2.1/fort.11',status='old')
       e(32)=2.1
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  320
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 320    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 32)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 32)=dsin(-d-ds)*ppp*ano*ano1
       freal( 32)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 32)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e2.1/fort.1',status='old')
       open(11,file='plus1/e2.1/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  321
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 321    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 32)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 32)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 32)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 32)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e2.2/fort.1',status='old')
       open(11,file='minus1/e2.2/fort.11',status='old')
       e(33)=2.2
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  330
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 330    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 33)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 33)=dsin(-d-ds)*ppp*ano*ano1
       freal( 33)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 33)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e2.2/fort.1',status='old')
       open(11,file='plus1/e2.2/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  331
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 331    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 33)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 33)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 33)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 33)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e2.5/fort.1',status='old')
       open(11,file='minus1/e2.5/fort.11',status='old')
       e(34)=2.5
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  340
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 340    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 34)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 34)=dsin(-d-ds)*ppp*ano*ano1
       freal( 34)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 34)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e2.5/fort.1',status='old')
       open(11,file='plus1/e2.5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  341
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 341    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 34)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 34)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 34)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 34)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e2.75/fort.1',status='old')
       open(11,file='minus1/e2.75/fort.11',status='old')
       e(35)=2.75
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  350
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 350    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 35)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 35)=dsin(-d-ds)*ppp*ano*ano1
       freal( 35)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 35)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e2.75/fort.1',status='old')
       open(11,file='plus1/e2.75/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  351
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 351    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 35)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 35)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 35)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 35)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e3/fort.1',status='old')
       open(11,file='minus1/e3/fort.11',status='old')
       e(36)=3.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  360
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 360    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 36)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 36)=dsin(-d-ds)*ppp*ano*ano1
       freal( 36)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 36)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e3/fort.1',status='old')
       open(11,file='plus1/e3/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  361
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 361    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 36)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 36)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 36)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 36)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e3.5/fort.1',status='old')
       open(11,file='minus1/e3.5/fort.11',status='old')
       e(37)=3.5
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  370
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 370    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 37)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 37)=dsin(-d-ds)*ppp*ano*ano1
       freal( 37)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 37)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e3.5/fort.1',status='old')
       open(11,file='plus1/e3.5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  371
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 371    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 37)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 37)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 37)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 37)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e4/fort.1',status='old')
       open(11,file='minus1/e4/fort.11',status='old')
       e(38)=4.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  380
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 380    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 38)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 38)=dsin(-d-ds)*ppp*ano*ano1
       freal( 38)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 38)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e4/fort.1',status='old')
       open(11,file='plus1/e4/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  381
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 381    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 38)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 38)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 38)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 38)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e4.5/fort.1',status='old')
       open(11,file='minus1/e4.5/fort.11',status='old')
       e(39)=4.5
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  390
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 390    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 39)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 39)=dsin(-d-ds)*ppp*ano*ano1
       freal( 39)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 39)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e4.5/fort.1',status='old')
       open(11,file='plus1/e4.5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  391
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 391    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 39)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 39)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 39)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 39)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e5/fort.1',status='old')
       open(11,file='minus1/e5/fort.11',status='old')
       e(40)=5.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  400
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 400    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 40)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 40)=dsin(-d-ds)*ppp*ano*ano1
       freal( 40)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 40)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e5/fort.1',status='old')
       open(11,file='plus1/e5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  401
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 401    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 40)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 40)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 40)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 40)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e5.5/fort.1',status='old')
       open(11,file='minus1/e5.5/fort.11',status='old')
       e(41)=5.5
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  410
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 410    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 41)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 41)=dsin(-d-ds)*ppp*ano*ano1
       freal( 41)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 41)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e5.5/fort.1',status='old')
       open(11,file='plus1/e5.5/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  411
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 411    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 41)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 41)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 41)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 41)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)

       open(1,file='minus1/e6/fort.1',status='old')
       open(11,file='minus1/e6/fort.11',status='old')
       e(42)=6.
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  420
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 420    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal( 42)=dcos(-d-ds)*ppp*ano*ano1
       gimag( 42)=dsin(-d-ds)*ppp*ano*ano1
       freal( 42)=dcos(-d-ds)*qqq*ano*ano2
       fimag( 42)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)
       open(1,file='plus1/e6/fort.1',status='old')
       open(11,file='plus1/e6/fort.11',status='old')
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(1)=x1x
       pp(1)=p*10**(2*(anorm-1))/x1x
       qq(1)=q*10**(2*(anorm-1))/x1x
       do i=1,10000
       read(1,*)x1x,p,q,anorm,aaas,pli
       xx(2)=x1x
       pp(2)=p*10**(2*(anorm-1))/x1x
       qq(2)=q*10**(2*(anorm-1))/x1x
       if(x1x.gt.r)then
       ppp=pp(1)+(pp(2)-pp(1))*(r-xx(1))/(xx(2)-xx(1))
       qqq=qq(1)+(qq(2)-qq(1))*(r-xx(1))/(xx(2)-xx(1))
       goto  421
         endif
       xx(1)=xx(2)
       pp(1)=pp(2)
       qq(1)=qq(2)
       enddo
 421    continue
       read(11,*)ano,ano1,ano2,ak,bohr,d,ds
       ano=ano/ak*bohr
       ano2=ano1
       greal1( 42)=dcos(-d-ds)*ppp*ano*ano1
       gimag1( 42)=dsin(-d-ds)*ppp*ano*ano1
       freal1( 42)=dcos(-d-ds)*qqq*ano*ano2
       fimag1( 42)=dsin(-d-ds)*qqq*ano*ano2
       close(1)
       close(11)


         return

        end





      subroutine mutim(ar,ai,br,bi,rr,ri)

      implicit double precision (a-h,o-z)
c multiply 2 complex number ar+iai and br+ibi
c and the result is in rr=ar br-ai bi and
c ri=ai br +ar bi
      rr=ar*br-ai*bi
      ri=ai*br+ar*bi
      return
      end



      SUBROUTINE SPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
C CHEAMA SUBROUTINA SPLAIS
C X,A,B,C,D,N,J,K AU ACEEASI SEMNIFICATIE CA LA SPLAIS
C Z-VECTOR DE INTRARE DE DIMENSIUNE M CONTININD NODURILE Zi IN
C   CARE VREM SA CALCULAM VALORILE INTERPOLATE Y(Zi) 1=<i=<M
C G-VECTOR DE IESIRE DIMENSIUNE M CONTININD VALORILE INTERPOLATE
C   Y(Zi) 1=<i=<M
C M-INTREG DE INTRARE CONTININD NUMARUL DE PUNCTE Zi
C DG derivatele lui G
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     common/vdsp/der
      DIMENSION X(200),A(200),B(200),C(200),D(200),Z(32),G(32)
      DIMENSION  A1A(200),DG(200),X1X(200)
      NA=N-1
       do ihj=1,n
       a(ihj)=a1a(ihj)
       x(ihj)=x1x(ihj)
       b(ihj)=0.d0
       c(ihj)=0.d0
       d(ihj)=0.d0
       enddo
      CALL SPLAIS(X,A,B,C,D,N,J,K)
       do ihj=1,n
       a1a(ihj)=a(ihj)
       x1x(ihj)=x(ihj)
       enddo
      DO 1 I=1,M
      W=Z(I)
      L=1
      DO 2 IX=2,NA
      IF(W.GT.X(IX))L=IX
2     CONTINUE
      P=W-X(L)
      Q=P*P
      R=Q*P
      E=C(L)*Q+D(L)*R
      G(I)=A(L)+B(L)*P+E
      der=b(l)+2*c(l)*p+3*d(l)*q
      DG(I)=DER
c     print*,'in  splaks  L,P=W-X(L),G(I),I,A(L),X(L),B(L),C(L),D'
c     print*,L,P,G(I),I,A(L),X(L),B(L),C(L),D(L)
1     CONTINUE
      RETURN
      END

      SUBROUTINE SPLAIS(X,A,B,C,D,N,J,K)
C X-VECTOR DE INTRARE DIMENSIUNE N CONTININD NODURILE RETELEI Xi
C   1.LE.i.LE.N
C A-VECTOR DE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   VALORILE IN NODURILE DATE Yi 1.LE.i.LE.N IAR LA IESIRE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C B-VECTORDE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   CAZUL UNU: B(1)=Y'(X1), B(N)=Y'(Xn)
C   CAZUL DOI: B(1)=Y''(X1),B(N)=Y''(Xn)
C   IN REST NU SE CERE INITIALIZARE. LA IESIRE CONTINE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C C,D-VECTORI DE IESIRE DIMENSIUNE N FIECARE, CONTIN COEFICENTII
C   POLINOAMELOR DE INTERPOLARE ASTFEL INCIT:
C      Y(X)=Pi(X)=A(I)+B(I)*(X-Xi)+C(I)*(X-Xi)**2+D(I)*(X-Xi)**3
C      daca x inclus in [Xi,Xi+1] 1.LE.i.LE.N-1
C N-INTREG DE INTRARE CARE SPECIFICA NR. DE NODURI DATE >=5
C J,K-INTREGI DE INTRARE AVIND VALORILE 1,2,3 DUPA CUM URMEAZA
C     J=1-NU DISPUNEM DE Y'(X1) SAU Y''(X1)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(X1)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(X1)
C     K=1-NU DISPUNEM DE Y'(Xn) SAU Y''(Xn)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(Xn)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(Xn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200),A(200),B(200),C(200),D(200)
C     DATA T/3.D0/,U/1.5D0/,V/0.25D0/,Z/0.5D0/
      T=3.D0
      U=1.5D0
      V=.25D0
      Z=.5D0
      M=N-1
      L=M-1
      D(1)=X(2)-X(1)
      IL=1
      DO 1 I=2,M
      IU=I+1
      D(I)=X(IU)-X(I)
      E=D(I)/D(IL)
      F=(A(IU)-A(I))/E
      G=(A(I)-A(IL))*E
      B(I)=T*(F+G)
1     IL=I
      GOTO(2,3,4),J
2     CONTINUE
      R=D(1)+D(2)
      E=T*D(1)+D(2)+D(2)
      F=D(2)/D(1)
      G=(A(3)-A(2))/F
      P=(A(2)-A(1))*F
      E=E*P+G*D(I)
      B(1)=E/R
      B(2)=B(2)-B(1)
      GO TO 5
3     CONTINUE
      R=D(1)+D(2)
      R=R+R
      B(2)=B(2)-D(2)*B(1)
      GO TO 5
4     CONTINUE
      R=D(1)+D(1)+U*D(2)
      E=(A(2)-A(1))/D(1)
      F=V*D(1)*B(1)
      Q=U*E-F
      B(2)=B(2)-D(2)*Q
5     CONTINUE
      GO TO(6,7,8),K
6     CONTINUE
      S=D(L)+D(M)
      E=T*D(M)+D(L)+D(L)
      F=D(L)/D(M)
      G=(A(M)-A(L))/F
      P=(A(N)-A(M))*F
      E=E*P+G*D(M)
      B(N)=E/S
      B(M)=B(M)-B(N)
      GO TO 9
7     CONTINUE
      S=D(L)+D(M)
      S=S+S
      B(M)=B(M)-D(L)*B(N)
      GO TO 9
8     CONTINUE
      S=D(M)+D(M)+U*D(L)
      E=(A(N)-A(M))/D(M)
      F=V*D(M)*B(N)
      W=U*E+F
      B(M)=B(M)-D(L)*W
9     CONTINUE
      C(2)=R
      DO 10 I=3,L
      IL=I-1
      E=D(IL)+D(I)
      F=D(I)/C(IL)
      G=F*D(IL-1)
      C(I)=E+E-G
      B(I)=B(I)-F*B(IL)
10    CONTINUE
      F=D(M)/C(L)
      C(M)=S-F*D(L-1)
      B(M)=B(M)-F*B(L)
      I=L
      B(M)=B(M)/C(M)
11    CONTINUE
      E=B(I)-D(I-1)*B(I+1)
      B(I)=E/C(I)
      I=I-1
      IF(I.GE.2) GO TO 11
      GO TO (12,13,14),J
12    CONTINUE
      E=B(1)-R*B(2)
      B(1)=E/D(2)
      GO TO 13
14    CONTINUE
      B(1)=Q-Z*B(2)
13    CONTINUE
      GO TO(15,16,17),K
15    CONTINUE
      E=B(N)-S*B(M)
      B(N)=E/D(L)
      GO TO 16
17    CONTINUE
      B(N)=W-Z*B(M)
16    CONTINUE
      DO 18 I=1,M
      IU=I+1
      H=D(I)
      E=H*H
      G=(A(IU)-A(I))/E
      P=B(IU)+B(I)
      Q=G/H
      F=(P+B(I))/H
      C(I)=T*G-F
      F=P/E
      D(I)=F-Q-Q
18    CONTINUE
c      print*,'in splais x a b c d'
c     do iuy=1,21
c      print*,x(iuy),a(iuy),b(iuy),c(iuy),d(iuy)
c     enddo
      RETURN
      END


      SUBROUTINE DQG32(XL,XU,FCT,Y)
C
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49863193092474078D0*B
      Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
      C=.49280575577263417D0*B
      Y=Y+.8137197365452835D-2*(FCT(A+C)+FCT(A-C))
      C=.48238112779375322D0*B
      Y=Y+.12696032654631030D-1*(FCT(A+C)+FCT(A-C))
      C=.46745303796886984D0*B
      Y=Y+.17136931456510717D-1*(FCT(A+C)+FCT(A-C))
      C=.44816057788302606D0*B
      Y=Y+.21417949011113340D-1*(FCT(A+C)+FCT(A-C))
      C=.42468380686628499D0*B
      Y=Y+.25499029631188088D-1*(FCT(A+C)+FCT(A-C))
      C=.39724189798397120D0*B
      Y=Y+.29342046739267774D-1*(FCT(A+C)+FCT(A-C))
      C=.36609105937014484D0*B
      Y=Y+.32911111388180923D-1*(FCT(A+C)+FCT(A-C))
      C=.33152213346510760D0*B
      Y=Y+.36172897054424253D-1*(FCT(A+C)+FCT(A-C))
      C=.29385787862038116D0*B
      Y=Y+.39096947893535153D-1*(FCT(A+C)+FCT(A-C))
      C=.25344995446611470D0*B
      Y=Y+.41655962113473378D-1*(FCT(A+C)+FCT(A-C))
      C=.21067563806531767D0*B
      Y=Y+.43826046502201906D-1*(FCT(A+C)+FCT(A-C))
      C=.16593430114106382D0*B
      Y=Y+.45586939347881942D-1*(FCT(A+C)+FCT(A-C))
      C=.11964368112606854D0*B
      Y=Y+.46922199540402283D-1*(FCT(A+C)+FCT(A-C))
      C=.7223598079139825D-1*B
      Y=Y+.47819360039637430D-1*(FCT(A+C)+FCT(A-C))
      C=.24153832843869158D-1*B
      Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END


      SUBROUTINE DQG321(XL,XU,FCT,Y)
C
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49863193092474078D0*B
      Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
      C=.49280575577263417D0*B
      Y=Y+.8137197365452835D-2*(FCT(A+C)+FCT(A-C))
      C=.48238112779375322D0*B
      Y=Y+.12696032654631030D-1*(FCT(A+C)+FCT(A-C))
      C=.46745303796886984D0*B
      Y=Y+.17136931456510717D-1*(FCT(A+C)+FCT(A-C))
      C=.44816057788302606D0*B
      Y=Y+.21417949011113340D-1*(FCT(A+C)+FCT(A-C))
      C=.42468380686628499D0*B
      Y=Y+.25499029631188088D-1*(FCT(A+C)+FCT(A-C))
      C=.39724189798397120D0*B
      Y=Y+.29342046739267774D-1*(FCT(A+C)+FCT(A-C))
      C=.36609105937014484D0*B
      Y=Y+.32911111388180923D-1*(FCT(A+C)+FCT(A-C))
      C=.33152213346510760D0*B
      Y=Y+.36172897054424253D-1*(FCT(A+C)+FCT(A-C))
      C=.29385787862038116D0*B
      Y=Y+.39096947893535153D-1*(FCT(A+C)+FCT(A-C))
      C=.25344995446611470D0*B
      Y=Y+.41655962113473378D-1*(FCT(A+C)+FCT(A-C))
      C=.21067563806531767D0*B
      Y=Y+.43826046502201906D-1*(FCT(A+C)+FCT(A-C))
      C=.16593430114106382D0*B
      Y=Y+.45586939347881942D-1*(FCT(A+C)+FCT(A-C))
      C=.11964368112606854D0*B
      Y=Y+.46922199540402283D-1*(FCT(A+C)+FCT(A-C))
      C=.7223598079139825D-1*B
      Y=Y+.47819360039637430D-1*(FCT(A+C)+FCT(A-C))
      C=.24153832843869158D-1*B
      Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END


      SUBROUTINE DQG322(XL,XU,FCT,Y)
C
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49863193092474078D0*B
      Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
      C=.49280575577263417D0*B
      Y=Y+.8137197365452835D-2*(FCT(A+C)+FCT(A-C))
      C=.48238112779375322D0*B
      Y=Y+.12696032654631030D-1*(FCT(A+C)+FCT(A-C))
      C=.46745303796886984D0*B
      Y=Y+.17136931456510717D-1*(FCT(A+C)+FCT(A-C))
      C=.44816057788302606D0*B
      Y=Y+.21417949011113340D-1*(FCT(A+C)+FCT(A-C))
      C=.42468380686628499D0*B
      Y=Y+.25499029631188088D-1*(FCT(A+C)+FCT(A-C))
      C=.39724189798397120D0*B
      Y=Y+.29342046739267774D-1*(FCT(A+C)+FCT(A-C))
      C=.36609105937014484D0*B
      Y=Y+.32911111388180923D-1*(FCT(A+C)+FCT(A-C))
      C=.33152213346510760D0*B
      Y=Y+.36172897054424253D-1*(FCT(A+C)+FCT(A-C))
      C=.29385787862038116D0*B
      Y=Y+.39096947893535153D-1*(FCT(A+C)+FCT(A-C))
      C=.25344995446611470D0*B
      Y=Y+.41655962113473378D-1*(FCT(A+C)+FCT(A-C))
      C=.21067563806531767D0*B
      Y=Y+.43826046502201906D-1*(FCT(A+C)+FCT(A-C))
      C=.16593430114106382D0*B
      Y=Y+.45586939347881942D-1*(FCT(A+C)+FCT(A-C))
      C=.11964368112606854D0*B
      Y=Y+.46922199540402283D-1*(FCT(A+C)+FCT(A-C))
      C=.7223598079139825D-1*B
      Y=Y+.47819360039637430D-1*(FCT(A+C)+FCT(A-C))
      C=.24153832843869158D-1*B
      Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END


      SUBROUTINE xSPLAKS(X1X,A1A,B,C,D,Nt,J,K,Z,G,nd,DG)
      implicit double precision (a-h,o-z)
      DIMENSION B(200),C(200),D(200),Z(32),G(32)
      DIMENSION  A1A(200),DG(200),X1X(200)

c     subroutine lplazvd(x1x,a1a,nt,z,g,nd,dg)
c linear interpolation
c nt mumar de noduri
c x1x valoarea nodurilor
c a1a valoarea in noduri
c z valorile in care se interpoleaza
c g valorile interpolate
c derg derivata in punctele interpolate
c numar de puncte interpolate
c     dimension x1x(200),a1a(200),
c    c   z(32),g(32),dg(200)
      do l=1,nd
      v=z(l)
      if(v.le.x1x(1))then
      g(l)=a1a(1)+(a1a(2)-a1a(1))*(x1x(1)-v)/(x1x(2)-x1x(1))
      dg(l)=(a1a(2)-a1a(1))/(x1x(2)-x1x(1))
      else
      if(v.ge.x1x(nt))then
      g(l)=a1a(nt)+(a1a(nt)-a1a(nt-1))*(v-x1x(nt))/(x1x(nt)-x1x(nt-1))
      dg(l)=(a1a(nt)-a1a(nt-1))/(x1x(nt)-x1x(nt-1))
      else
      ipos=0
      do i=2,nt
      if(v.le.x1x(i).and.ipos.eq.0)then
      ipos=1
      nrm=i
      endif
      enddo
      g(l)=a1a(nrm-1)+(a1a(nrm)-a1a(nrm-1))*(v-x1x(nrm-1))/(x1x(nrm)-
     c                  x1x(nrm-1))
      dg(l)=(a1a(nrm)-a1a(nrm-1))/(x1x(nrm)-x1x(nrm-1))
      endif
      endif
      enddo
      return
      end






