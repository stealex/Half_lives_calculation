cc  Calculez functia csi care este modul patrat a fiecarei functii de unda
      subroutine suprafata(a0,s0,s,raps)
c calculez suprafata nucleului
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,deltax
      common/r116/R0
       common/suprafete/s1,s2,s3
      pi=3.141592645d0
      s0=4.*pi*(R0*a0**(1./3.))**2
      
c partea mediana doar in cazul in care ro3\ne 0
      S3=0
      if(r3.gt.300.d0)then
      radic=dsqrt(1.+((y201-y101)/(x201-x101))**2)
      s3=2*pi*(x201-x101)*y101*radic-
     c   2*pi*(y201-y101)*x101*radic+
     c   pi*(y201-y101)/(x201-x101)*(x201**2-x101**2)*radic
      goto 10
      endif
      if(ro3.eq.0.or.r3.lt.0.003d0)goto 10
c am doua cazuri x1=u1 x2=u2 si x1\ne u1 x2\ne u2
      print*,'sp',sp
       dif3=ro3-r3
      if(sp.eq.-1.or.(sp.eq.1..and.dif3.gt.0.d0))then
      S3=2*pi*r3*(x1-x2)*sp
      S3=S3+2*pi*ro3*r3*(dasin((x2-c3)/r3)-dasin((x1-c3)/r3))
      else
      S3=2*pi*r3*(u2-x2)*sp
      S3=S3+2*pi*ro3*r3*(dasin((x2-c3)/r3)-dasin((u2-c3)/r3))
      S3=S3+2*pi*r3*(x1-u1)*sp
      S3=S3+2*pi*ro3*r3*(dasin((u1-c3)/r3)-dasin((x1-c3)/r3))
      endif
10    continue


c suprafata fragment 1 partea stanga
c  am trei cazuri sfere, prolate, oblate
      diferenta=1.d0-b1p/a1p
      if(dabs(diferenta).lt.0.001d0)then !cazul sferic
      S1=2.*pi*a1p*(x1-c1+a1p)
      goto 20
      else
      if(diferenta.gt.0.d0)then ! cazul prolate
      rap=1.d0-b1p**2/a1p**2
      raport=1.d0/(1.d0-b1p**2/a1p**2)
      radical=dsqrt(1.d0/(1.d0-b1p**2/a1p**2))
      S1=2*pi*a1p*b1p*dsqrt(rap)*(
     c  (x1-c1)/a1p*dsqrt(raport-((x1-c1)/a1p)**2)/2.d0+
     c  raport/2.d0*dasin((x1-c1)/a1p/radical)+
     c  dsqrt(raport-1)/2.d0+
     c  raport/2.d0*dasin(1/radical))
      else
      rap=b1p**2/a1p**2-1
      raport=1.d0/(b1p**2/a1p**2-1)
      radical=dsqrt(1.d0/(b1p**2/a1p**2-1))
      S1=2*pi*a1p*b1p*dsqrt(rap)*(
     c   (x1-c1)/a1p*dsqrt(raport+((x1-c1)/a1p)**2)/2.d0+
     c   raport/2.d0*dlog((x1-c1)/a1p+dsqrt(((x1-c1)/a1p)**2+raport))+
     c   dsqrt(raport+1)/2.d0-
     c   raport/2.d0*dlog(-1+dsqrt(1+raport)))
      endif
      endif
20    continue

c suprafata fragment 2 partea dreapta
c se face la fel ca fragment 1 facandu-se urmatoarele substitutii
      xx1=-x2
      cc1=-c2
      bb1p=b2p
      aa1p=a2p
c  am trei cazuri sfere, prolate, oblate
      diferenta=1.d0-bb1p/aa1p
      if(dabs(diferenta).lt.0.001d0)then !cazul sferic
      S2=2.*pi*aa1p*(xx1-cc1+aa1p)
      goto 30
      else
      if(diferenta.gt.0.d0)then ! cazul prolate
      rap=1.d0-bb1p**2/aa1p**2
      raport=1.d0/(1.d0-bb1p**2/aa1p**2)
      radical=dsqrt(1.d0/(1.d0-bb1p**2/aa1p**2))
      S2=2*pi*aa1p*bb1p*dsqrt(rap)*(
     c  (xx1-cc1)/aa1p*dsqrt(raport-((xx1-cc1)/aa1p)**2)/2.d0+
     c  raport/2.d0*dasin((xx1-cc1)/aa1p/radical)+
     c  dsqrt(raport-1)/2.d0+
     c  raport/2.d0*dasin(1/radical))
      else
      rap=bb1p**2/aa1p**2-1
      raport=1.d0/(bb1p**2/aa1p**2-1)
      radical=dsqrt(1.d0/(bb1p**2/aa1p**2-1))
      S2=2*pi*aa1p*bb1p*dsqrt(rap)*(
     c   (xx1-cc1)/aa1p*dsqrt(raport+((xx1-cc1)/aa1p)**2)/2.d0+
     c   raport/2.d0*dlog((xx1-cc1)/aa1p+
     c                  dsqrt(((xx1-cc1)/aa1p)**2+raport))+
     c   dsqrt(raport+1)/2.d0-
     c   raport/2.d0*dlog(-1+dsqrt(1+raport)))
      endif
      endif
30    continue


      s=s1+s2+s3
      raps=s/s0


      return
      end


c             subroutine main
      implicit double precision (a-h,o-z)
      dimension enou(2925),een(2925)
      dimension tray(50),trax(50),trayy(50)
      dimension amc2(150),zmc2(150)

      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,deltax
      common/dltw/dltw
      common/r116/R0
      common/ier3parelp/ier3
      common/gradlib/ga0,geps1,geps2,ga1pa2,gcr3,gdelt
      common/lasciziune/delatsciz
      common/dvtr/dereps1,dereps2,dercr3,dera
c     common/rapginte/rapgint 
       common/suprafete/s1,s2,s3
          common/rapsup/rapsu
        common/efep/efep,dufp,dpfp,duep,dpep,dufn,dpfn,duen,dpen
        common/initfin/a2final,z2final,initfin
c initfin=1 atunci cand calculez valorile indicelui care-mi determina prima sau a doua
c groapa, nrfinale de masa si de sarcina a fragmentului usor a2final,z2final

      dimension ccc(7)
      DIMENSION ENn(2925),eimpn(10)
      dimension infinn(2925),nrspinn(2925)
      DIMENSION ENp(2925),eimp(10)
      dimension infinp(2925),nrspinp(2925)

      dimension x1x(200),a1a(200),b1b(200),c1c(200),d1d(200)
      dimension z1z(32),g1g(32),dg(200)
      dimension crti(28),eps1i(28),eps2i(28),a1pa2i(28)
      dimension ddd(28),ddda(28)
           
      dimension rnou(50),crnou(50)
      dimension na2(200),nz2(200),naa2(200),eps(200)
      open(22,file='asiz.dat',status='old')
      read(22,*)numarmasa,numarsarcina
      close(22)
c     open(79,file='partitii.dat',status='old')
        n56=72
        n112=144
c      do k=1,n56
c      read(79,*)na2(k),na1,nz2(k),nz1
c      enddo
c      close(79)

                 initfin=0

               nvmax=-1000
               nvmin=1000
c      do k=1,n112
c      read(83,*)aaa2,zzz2,eps(k)
c      naa2(k)=aaa2
c      if(nvmax.lt.naa2(k))nvmax=naa2(k)
c      if(nvmin.gt.naa2(k))nvmin=naa2(k)
c      enddo


       
             eps1fi=0.6
             eps2fi=0.

           sciziunea0=1.16*((248.)**(1./3.)/(1-eps1fi**2)**(1./3.)
     c      +48.**(1./3.)/(1-eps2fi**2)**(1./3.))
           sciziunea=sciziunea0
        print*,sciziunea
           raportul0=((248.)/48.)**(1./3.)
            raportul=raportul0
           touching0=sciziunea0

 
       trax(1)=0.9
       trax(2)=2.
       trax(3)=3.04
       trax(4)=4.811
       trax(5)=6.24 
       trax(6)=7.15 
       trax(7)=10.227
       trax(8)=18.   
       trax(9)=100.
       trayy(1)=0.02
       trayy(2)=0.055
       trayy(3)=0.15 
       trayy(4)=0.533+0.05 
       trayy(5)=0.92   
       trayy(6)=0.973   
       trayy(7)=0.973  
       trayy(8)=0.973 
       trayy(9)=0.973  


c            do l=1,9 
c            tray(l)=296*(1-tray(l))/2.
c            if(l.gt.5)tray(l)=4.d0
c            a0=296
c      tray(l)=(((a0-tray(l))/tray(l))*(1-eps2f**2)/(1-eps1f**2))**(1./3.)
c            enddo


      sciz0=sciziunea0

           slin0=raportul0



       rnou(1)=0.1/sciziunea
       rnou(2)=5.52/sciziunea
       rnou(3)=7.68/sciziunea
       rnou(4)=11.13/sciziunea
       rnou(5)=12.2/sciziunea
       rnou(6)=sciziunea/sciziunea 
       rnou(7)=10000.     
       print*,'sciziune',sciziunea,'raportul',raportul
       crnou(1)=-0.07
       crnou(2)=-0.07
       crnou(3)=0.18 
       crnou(4)=5!0.400
       crnou(5)=10!0.500
       crnou(6)=500   
       crnou(7)=5000  


!       do ia1pa2=1,36
!        do ia2=4,148,2
        do idelta=1,1




       delta=.01

   




ccccccccccccc              ma intereseaza valoarea lui ia2


               a0=numarmasa 
               z0=numarsarcina
  


            a1pa2=1.


      a2e=(numarmasa-2)/2
      z2e=(numarsarcina-2)/2
      ia2=a2e

                 n10nou=7!16


           aaa2=ia2

     
             eps1fi=0.
             eps2fi=0.

           zzz2=ia2/2
 

       scizmom=1.16*((a0-aaa2)**(1./3.)/(1-eps1fi**2)**(1./3.)
     c   +aaa2**(1./3.)/(1-eps2fi**2)**(1./3.))

                a2final=aaa2
                z2final=zzz2    


     

           cr3=-0.05


            dereps1=0
            dereps2=0


        ! calculez un eps1final care sa fie functie de masa
        ! fiicei (este zero pentru plumb) si este 0.7 pentru 232
          
          
             
               eps1f=eps1fi
               eps2f=eps2fi
               dereps1=0
               dereps2=0
              

          touching=scizmom                                   
            eps1=eps1f
            eps2=eps2f


             eps0=0





        
!c cu formula urmatoarea ma asigur ca voi avea acelasi raport al
!c  volumelor ca in cazul in care nu am introdus eccentricitati
!        a1pa2=a1pa2*((1.-eps2**2)/(1.-eps1**2))**(1./3.)
!           aaa2m=ia2
!      a1pa2=(((a0-aaa2m)/aaa2m)*(1-eps2f**2)/(1-eps1f**2))**(1./3.)

     
c excentricitate zero=sfera
c parametrizez eps1,eps2,a1pa2 (excentricitati si raport intre semi-axele a1 si a2)
      r0=1.16 ! valoarea lui r0 pentru calculul macroscopic

c      call vefinel(a0,z0,a2e,z2e,eps0,eps1f,eps2f,cr3,delta,
c     c eps1,eps2,a1pa2,a2m,z2m)

      call vefinelx(a0,z0,a1pa2,deltai,eps0,eps1f,eps2f,
     c a2m,z2m)


          amij=a0/2
          aaa2lim1=aaa2-1
          aaa2lim2=aaa2+1
          if(aaa2lim1.lt.amij.and.aaa2lim2.gt.amij)goto 921
       z2m=zzz2+(z0/2.-zzz2)*(a2m-aaa2)/(a0/2.-aaa2)
921     continue
c      if(a2m.gt.48.d0)then
c      z2m=20+(z2m-20)*(a2m-48)/(a0/2.-48)
c      else
c      z2m=20+(z2m-20)*(a2m-48)/(a0/2.-48)
c      z2m=2+(z2m-2)*(a2m-4)/(48-2)
c      endif

         naaaa2=dint(a2m+0.1)
         nzzzz2=dint(z2m+0.1)
      a2m=naaaa2
      z2m=nzzzz2




       a2e=a2m
       z2e=z2m
!       deltai=dabs(deltai)+0.003
      delta0=delta
!         if(delta.le.deltai)delta=deltai
c         deltap=delta*1.16/1.27
          deltap=delta
        if(deltap.le.0.01)deltap=0.01
         cr3p=cr3*1.27/1.16
      call parelp2(a0,eps1,eps2,a1pa2,cr3p,deltap,a1,b1,a2,b2)
      cc3=cr3
      call psciziune(a1,b1,a2,b2,cc3,deltasciz,ier1,ier2,ier3,ier4)
      delatsciz=deltasciz
      delatsciz=scizmom
c          print*,'deltasciz',deltasciz
       write(88,*)delta,a1,b1,a2,b2,cc3,a1pa2,deltasciz,dercr3,dera
        print*,'delta,cr3',delta,cr3
      call suprafata(a0,sup0,sup,raps)
          rapsu=raps
c           rapsu=1.d0 !elimin dependenta de suprafata a lui G
c          goto 111
c       write(88,*)delta,a1,b1,a2,b2,cc3,deltasciz
      if(ier1.ne.0.or.ier2.ne.0.or.ier3.ne.0.or.ier4.ne.0)then
       print*,'a1,b1,a2,b2,cc3,deltasciz'
       print*,a1,b1,a2,b2,cc3,deltasciz
       print*,'ier1,...ier4',ier1,ier2,ier3,ier4
       print*,'idelta',idelta
c      stop
      endif
        nump=0
c       if(delta.lt.deltai)nump=1
        print*,'ier !!!!!!!!!!!!!!',ier11,ier22,ier3
      if(ier11.ne.0.or.ier22.ne.0.or.ier3.ne.0.or.nump.eq.1)then
       poten=1
       itest=i
      efeps1p=1
      efeps2p=1
      efa1pa2p=1
      efcr3p=1
      efdeltp=1
      efeps1n=1
      efeps2n=1
      efa1pa2n=1
      efcr3n=1
      efdeltn=1

      efeps1eps2p=1
      efeps1a1pa2p=1
      efeps1cr3p=1
      efeps1deltp=1
      efeps2a1pa2p=1
      efeps2cr3p=1
      efeps2deltp=1
      efa1pa2cr3p=1
      efa1pa2deltp=1
      efcr3deltp=1
      efeps1eps2n=1
      efeps1a1pa2n=1
      efeps1cr3n=1
      efeps1deltn=1
      efeps2a1pa2n=1
      efeps2cr3n=1
      efeps2deltn=1
      efa1pa2cr3n=1
      efa1pa2deltn=1
      efcr3deltn=1
      goto 8745
      endif
      delt=delta0
      CALL MPE(A0,Z0,A2m,Z2m,EC0,EC,EY0,ENY,EV,ED0,ED,ew,V)
               vldm=v
           call grafic(a0,z0,a2e,z2e,eps1,eps2,a1pa2,cr3,delta)  
      qvalue=91.19161      
      vqvalue=-qvalue-6.56+(Z0-Z2m)*Z2m/delta*1.4399764+0.5
      r0=1.27 ! valoarea lui r0 pentru calculul microscopic
c      call vefinel(a0,z0,a2e,z2e,eps0,eps1f,eps2f,cr3,delta,
c     c eps1,eps2,a1pa2,a2mm,z2mm)
      call vefinelx(a0,z0,a1pa2,deltai,eps0,eps1f,eps2f,
     c a2mm,z2mm)
      delta0=delta
!         deltai=dabs(deltai)+0.003
!         if(delta.le.deltai)delta=deltai

      ga0=a0
      geps1=eps1
      geps2=eps2
      ga1pa2=a1pa2
      gcr3=cr3
      gdelt=delta
             print*,'intru nrasimtp'
c                   goto 99
      NMAX=12
      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
        print*,'ier !!!!!!!!!!!!!!',ier11,ier22,ier3


      nrimpa=0
      ga0=a0
      geps1=eps1
      geps2=eps2
      ga1pa2=a1pa2
      gcr3=cr3
      gdelt=delta
          dv=0
          if(delta.gt.30.d0)goto 85
       call micros(a0,z0,nmax,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
85     continue
        poten=vldm+dv
c       if(delta.gt.touching)then
c       pppi=poten
c       poten=vqvalue
c       toucj3=touching+3
c       if(delta.lt.toucj3)then
c       poten=pppi+(vqvalue-pppi)*(delta-touching)/(toucj3-touching)
c       endif
c       endif

        amaseftot=(efeps1)*dereps1**2+
     c (efeps2)*dereps2**2+(efa1pa2)*dera**2+
     c (efcr3)*dercr3**2+efdelt+
     c 2*(efeps1eps2)*dereps1*dereps2+
     c 2*(efeps1a1pa2)*dereps1*dera+
     c 2*(efeps1cr3)*dereps1*dercr3+
     c 2*(efeps1delt)*dereps1+
     c 2*(efeps2a1pa2)*dereps2*dera+
     c 2*(efeps2cr3)*dereps2*dercr3+
     c 2*(efeps2delt)*dereps2+
     c 2*(efa1pa2cr3)*dera*dercr3+
     c 2*(efa1pa2delt)*dera+
     c 2*(efcr3delt)*dercr3
c99             continue
8745          continue
      amasr=efdelt
      amasc=(efcr3)*dercr3**2
      amasa1pa2=(efa1pa2)*dera**2
         if(initfin.eq.1)goto 435
      write(78,*)delta,amasr,amasc,amasa1pa2,efcr3,efa1pa2,dercr3,dera
        write(24,*)efep,dufp,dpfp,duep,dpep,dufn,dpfn,duen,dpen

      write(60,*)delta,A0,Z0,A2m,Z2m,EC0,EC,EY0,ENY,EV,ED0,ED,ew,Vldm
      write(59,*)a2m,z2m
      write(57,*)k1,iasim,idelta,aaa2,zzz2,cr3,delta,poten,amaseftot
      write(18,*)delta,cr3,eps1,eps2,a1pa2,dupro,dppro,du,dp,vldm,dv
      write(51,*)dercr3,dera,dereps1,dereps2
      write(41,*)efeps1,efeps2,efa1pa2,efcr3,efdelt,nimi
      write(42,*)efeps1eps2,
     c  efeps1a1pa2,efeps1cr3,efeps1delt,efeps2a1pa2,nimi
      write(43,*)efeps2cr3,
     c  efeps2delt,efa1pa2cr3,efa1pa2delt,efcr3delt,nimi
435      continue
  55             continue
111      continue
               enddo
c               enddo
c              enddo

      end
      


       subroutine micros(a0,z0,nmaxim,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
! intrari: numere parinte, numar parturi majore<=14, numere emis,
! distanta intre centre, gatuire, deformari, raport a1/a2
! iesire energii,spini,nrfinale pentru protoni (dimensiune 2925)
!  energii, spini, nrfinale pentru neutroni (dimensiune 2925)
! (energii ordonate descrescator)
! numar de nivele la iesire, corectii paturi si imprechere protoni,
! corectii paturi si imperechere neutroni, corectii totale,
! energie Fermi si gap protoni, energie Fermi si ngap neutroni,
! 15 mase efective
      implicit double precision (a-h,o-z)
      dimension enp(2925),nrspinp(2925),infinp(2925)
      dimension enn(2925),nrspinn(2925),infinn(2925)
      dimension enpm(2925),nrspinpm(2925),infinpm(2925)
      dimension ennm(2925),nrspinnm(2925),infinnm(2925)
!!!!! prin urmatorul common se trimit valorile initiale
      common/gradlib/ga0,geps1,geps2,ga1pa2,gcr3,gdelt
      common/enfermin/enfermip,enfermin


            aaa0=a0
            zzz0=z0
            aa2e=a2e
            zz2e=z2e
            eps11=eps1
            eps22=eps2
            a1pa22=a1pa2
            deltaa=delta
            cr33=cr3
            ga0=aaa0
            geps1=eps11
            geps2=eps22
            ga1pa2=a1pa22
            gcr3=cr33
            gdelt=deltaa
           xrapx=a1pa2
           if(xrapx.ge.1.0455d0.or.xrapx.le.0.9565d0)then
      CALL microscopic(a0,z0,nmaxim,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
           do inn=1,nrniv
           write(61,*)inn,delta,enp(inn),nrspinp(inn),infinp(inn)
           write(62,*)inn,delta,enn(inn),nrspinn(inn),infinn(inn)
           enddo
c          write(63,*)delta,enfermip
c          write(64,*)delta,enfermin
       return
           endif
          if(xrapx.ge.0.999d0.and.xrapx.le.1.001d0)then
              ga1pa2=1.
             a2ef=a0/2.
             z2ef=z0/2.
              a1pa2sim=1.
      CALL microscopic(a0,z0,nmaxim,a2ef,z2ef,delta,cr3,
     c eps1,eps2,a1pa2sim,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
           do inn=1,nrniv
           write(61,*)inn,delta,enp(inn),nrspinp(inn),infinp(inn)
           write(62,*)inn,delta,enn(inn),nrspinn(inn),infinn(inn)
           enddo
c          write(63,*)delta,enfermip
c          write(64,*)delta,enfermin
           return

           endif
             if(xrapx.lt.1.d0)then
             a1pa2max=1.
             a1pa2min=0.9565
            ga0=aaa0
            geps1=eps11
            geps2=eps22
             ga1pa2=a1pa2max
            gcr3=cr33
            gdelt=deltaa
             ga1pa2=a1pa2max
      CALL microscopic(a0,z0,nmaxim,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2max,
     c enpm,nrspinpm,infinpm,ennm,nrspinnm,infinnm,nrniv,
     c duprom,dpprom,dum,dpm,dvm,sollampm,soldelpm,sollamm,soldelm,
     c efeps1m,efeps2m,efa1pa2m,efcr3m,
     c efdeltm,efeps1eps2m,efeps1a1pa2m,efeps1cr3m,efeps1deltm,
     c efeps2a1pa2m,efeps2cr3m,efeps2deltm,efa1pa2cr3m,
     c efa1pa2deltm,efcr3deltm)
      print*,'delta,cr3,eps1,eps2,a1pa2max'
      print*,delta,cr3,eps1,eps2,a1pa2max
      print*,'duprom,dpprom,dum,dpm,dvm'
      print*,duprom,dpprom,dum,dpm,dvm
        enfermipmax=enfermip
        enferminmax=enfermin
            ga0=aaa0
            geps1=eps11
            geps2=eps22
             ga1pa2=a1pa2min
            gcr3=cr33
            gdelt=deltaa
             ga1pa2=a1pa2min
            delta=deltaa
            cr3=cr33
            eps1=eps11
            eps2=eps22
      CALL microscopic(aaa0,zzz0,nmaxim,aa2e,zz2e,delta,cr3,
     c eps1,eps2,a1pa2min,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
      print*,'delta,cr3,eps1,eps2,a1pa2'
      print*,delta,cr3,eps1,eps2,a1pa2
      print*,'dupro,dppro,du,dp,dv'
      print*,dupro,dppro,du,dp,dv
          do i=1,nrniv
       enp(i)=enp(i)+(enpm(i)-enp(i))*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       enn(i)=enn(i)+(ennm(i)-enn(i))*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
          enddo
       dupro=dupro+(duprom-dupro)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dppro=dppro+(dpprom-dppro)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       du=du+(dum-du)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dp=dp+(dpm-dp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dv=dv+(dvm-dv)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       sollamp=sollamp+(sollampm-sollamp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       soldelp=soldelp+(soldelpm-soldelp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       sollam=sollam+(sollamm-sollam)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       soldel=soldel+(soldelm-soldel)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps1=efeps1+(efeps1m-efeps1)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps2=efeps2+(efeps2m-efeps2)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efa1pa2=efa1pa2+(efa1pa2m-efa1pa2)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efcr3=efcr3+(efcr3m-efcr3)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efdelt=efdelt+(efdeltm-efdelt)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps1eps2=efeps1eps2+(efeps1eps2m-efeps1eps2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1a1pa2=efeps1a1pa2+(efeps1a1pa2m-efeps1a1pa2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1cr3=efeps1cr3+(efeps1cr3m-efeps1cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1delt=efeps1delt+(efeps1deltm-efeps1delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2a1pa2=efeps2a1pa2+(efeps2a1pa2m-efeps2a1pa2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2cr3=efeps2cr3+(efeps2cr3m-efeps2cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2delt=efeps2delt+(efeps2deltm-efeps2delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efa1pa2cr3=efa1pa2cr3+(efa1pa2cr3m-efa1pa2cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efa1pa2delt=efa1pa2delt+(efa1pa2deltm-efa1pa2delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efcr3delt=efcr3delt+(efcr3deltm-efcr3delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
           do inn=1,nrniv
           write(61,*)inn,delta,enp(inn),nrspinp(inn),infinp(inn)
           write(62,*)inn,delta,enn(inn),nrspinn(inn),infinn(inn)
           enddo
        enfermip=enfermip+(enfermipmax-enfermip)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
        enfermin=enfermin+(enferminmax-enfermin)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
c          write(63,*)delta,enfermip
c          write(64,*)delta,enfermin
       return
       else

             a1pa2max=1.0455
             a1pa2min=1.
      print*,'fac a1pa2max'
            ga0=aaa0
            geps1=eps11
            geps2=eps22
             ga1pa2=a1pa2max
            gcr3=cr33
            gdelt=deltaa
      CALL microscopic(a0,z0,nmaxim,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2max,
     c enpm,nrspinpm,infinpm,ennm,nrspinnm,infinnm,nrniv,
     c duprom,dpprom,dum,dpm,dvm,sollampm,soldelpm,sollamm,soldelm,
     c efeps1m,efeps2m,efa1pa2m,efcr3m,
     c efdeltm,efeps1eps2m,efeps1a1pa2m,efeps1cr3m,efeps1deltm,
     c efeps2a1pa2m,efeps2cr3m,efeps2deltm,efa1pa2cr3m,
     c efa1pa2deltm,efcr3deltm)
      print*,'delta,cr3,eps1,eps2,a1pa2max'
      print*,delta,cr3,eps1,eps2,a1pa2max
      print*,'duprom,dpprom,dum,dpm,dvm'
      print*,duprom,dpprom,dum,dpm,dvm
        enfermipmax=enfermip
        enferminmax=enfermin
            ga0=aaa0
            geps1=eps11
            geps2=eps22
             ga1pa2=a1pa2min
            gcr3=cr33
            gdelt=deltaa
            delta=deltaa
            cr3=cr33
            eps1=eps11
            eps2=eps22
      CALL microscopic(aaa0,zzz0,nmaxim,aa2e,zz2e,delta,cr3,
     c eps1,eps2,a1pa2min,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
      print*,'delta,cr3,eps1,eps2,a1pa2'
      print*,delta,cr3,eps1,eps2,a1pa2
      print*,'dupro,dppro,du,dp,dv'
      print*,dupro,dppro,du,dp,dv
c     write(93,*)delta,cr3,eps1,eps2,a1pa2max,duprom,dpprom,dum,dpm,dvm
c     write(93,*)delta,cr3,eps1,eps2,a1pa2min,dupro,dppro,du,dp,dv
          do i=1,nrniv
       enp(i)=enp(i)+(enpm(i)-enp(i))*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       enn(i)=enn(i)+(ennm(i)-enn(i))*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
          enddo
       dupro=dupro+(duprom-dupro)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dppro=dppro+(dpprom-dppro)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       du=du+(dum-du)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dp=dp+(dpm-dp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       dv=dv+(dvm-dv)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       sollamp=sollamp+(sollampm-sollamp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       soldelp=soldelp+(soldelpm-soldelp)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       sollam=sollam+(sollamm-sollam)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       soldel=soldel+(soldelm-soldel)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps1=efeps1+(efeps1m-efeps1)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps2=efeps2+(efeps2m-efeps2)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efa1pa2=efa1pa2+(efa1pa2m-efa1pa2)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efcr3=efcr3+(efcr3m-efcr3)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efdelt=efdelt+(efdeltm-efdelt)*(xrapx-a1pa2min)/
     c        (a1pa2max-a1pa2min)
       efeps1eps2=efeps1eps2+(efeps1eps2m-efeps1eps2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1a1pa2=efeps1a1pa2+(efeps1a1pa2m-efeps1a1pa2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1cr3=efeps1cr3+(efeps1cr3m-efeps1cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps1delt=efeps1delt+(efeps1deltm-efeps1delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2a1pa2=efeps2a1pa2+(efeps2a1pa2m-efeps2a1pa2)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2cr3=efeps2cr3+(efeps2cr3m-efeps2cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efeps2delt=efeps2delt+(efeps2deltm-efeps2delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efa1pa2cr3=efa1pa2cr3+(efa1pa2cr3m-efa1pa2cr3)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efa1pa2delt=efa1pa2delt+(efa1pa2deltm-efa1pa2delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
       efcr3delt=efcr3delt+(efcr3deltm-efcr3delt)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)

           do inn=1,nrniv
           write(61,*)inn,delta,enp(inn),nrspinp(inn),infinp(inn)
           write(62,*)inn,delta,enn(inn),nrspinn(inn),infinn(inn)
           enddo
        enfermip=enfermip+(enfermipmax-enfermip)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
        enfermin=enfermin+(enferminmax-enfermin)*
     c        (xrapx-a1pa2min)/(a1pa2max-a1pa2min)
c          write(63,*)delta,enfermip
c          write(64,*)delta,enfermin
        return
        endif
       end


       subroutine microscopic(a0,z0,nmaxim,a2e,z2e,delta,cr3,
     c eps1,eps2,a1pa2,
     c enp,nrspinp,infinp,enn,nrspinn,infinn,nrniv,
     c dupro,dppro,du,dp,dv,sollamp,soldelp,sollam,soldel,
     c efeps1,efeps2,efa1pa2,efcr3,
     c efdelt,efeps1eps2,efeps1a1pa2,efeps1cr3,efeps1delt,
     c efeps2a1pa2,efeps2cr3,efeps2delt,efa1pa2cr3,
     c efa1pa2delt,efcr3delt)
! intrari: numere parinte, numar parturi majore<=14, numere emis,
! distanta intre centre, gatuire, deformari, raport a1/a2
! iesire energii,spini,nrfinale pentru protoni (dimensiune 2925)
!  energii, spini, nrfinale pentru neutroni (dimensiune 2925)
! (energii ordonate descrescator)
! numar de nivele la iesire, corectii paturi si imprechere protoni,
! corectii paturi si imperechere neutroni, corectii totale,
! energie Fermi si gap protoni, energie Fermi si ngap neutroni,
! 15 mase efective
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,deltax
      common/dltw/dltw
      common/r116/R0
      common/ier3parelp/ier3
      common/gradlib/ga0,geps1,geps2,ga1pa2,gcr3,gdelt
      common/nrizospin/nizospin
      common/enfermin/enfermip,enfermin
      common/sumg/sumeps1,sumeps2,suma1pa2,sumcr3,sumdelt,sumeps1eps2,
     c sumeps1a1pa2,sumeps1cr3,sumeps1delt,sumeps2a1pa2,sumeps2cr3,
     c sumeps2delt,suma1pa2cr3,suma1pa2delt,sumcr3delt
        common/initfin/a2final,z2final,initfin
       common/enrde/enrgd(2925),nspind(2925),infnpd(2925)
! numere asimptotice
      common/nrasimptotice/nrasimptp(2925),nrasimptn(2925)
      common/lasciziune/delatsciz
       common/nrnivc/gpair,nrnivcpair,iierpair
       common/nrnivcmass/gmass,nrnivcmass,iiermass
      common/hw0ptg/hw00
        common/partnrs/eonri(2925),eimpari(10),nrnivei,nrimparei,nrocupi
        common/partnr2/eonr2(2925),eimpar2(10),nrnive2,nrimpare2,nrocup2
        common/efep/efep,dufp,dpfp,duep,dpep,dufn,dpfn,duen,dpen
        common/denstt/sollam0,soldel0,sollam10,soldl10,rrsc,nmaxy
      common/vectoripro/VECPRO(25,325,325),LELMAX(25)
          common/distanta/distanta
cccc Common cu valori densitate uniparticula in noduri Gauss Legendre
      common/densitup/dens(2,100,100),n1,n2,n3
c  n1=1 sau 2 pentru protoni sau neutroni
c n1 si n2 array pentru z si rho
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      dimension ccc(7)
      DIMENSION EN(2925),eimp(10)
      dimension infin(2925),nrspin(2925)
      dimension enp(2925),nrspinp(2925),infinp(2925)
      dimension enn(2925),nrspinn(2925),infinn(2925)
!     dimension enimpar(2925)
      dimension x1x(200),a1a(200),b1b(200),c1c(200),d1d(200)
      dimension z1z(32),g1g(32),dg(200)
      dimension crti(18),eps1i(18),eps2i(18),a1pa2i(18)
      dimension enimp(2925)
      dimension enin(2925)
      dimension gelm(300,300),gelm1(300,300),gelm2(300,300),nordin(2925)
      dimension nordinp(2925),nordinn(2925)
      dimension vecintermp(25,325,325),vecintermn(25,325,325)
            dimension jfinp(2925),jfinn(2925) !nr cuantice gropi finale ordonate     
            dimension ig1(2925),ig2(2925)
           nimi=0
      R0=1.27d0
      NMAX=nmaxim


c excentricitate zero=sfera
c parametrizez eps1,eps2,a1pa2 (excentricitati si raport intre semi-axele a1 si a2)

      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
        print*,'ier !!!!!!!!!!!!!!',ier11,ier22,ier3
          print*,'a0',a0,'dupa parelp2'

      if(ier11.ne.0.or.ier22.ne.0.or.ier3.ne.0)then
       poten=1000
       itest=i
      efeps1p=1
      efeps2p=1
      efa1pa2p=1
      efcr3p=1
      efdeltp=1
      efeps1n=1
      efeps2n=1
      efa1pa2n=1
      efcr3n=1
      efdeltn=1

      efeps1eps2p=1
      efeps1a1pa2p=1
      efeps1cr3p=1
      efeps1deltp=1
      efeps2a1pa2p=1
      efeps2cr3p=1
      efeps2deltp=1
      efa1pa2cr3p=1
      efa1pa2deltp=1
      efcr3deltp=1
      efeps1eps2n=1
      efeps1a1pa2n=1
      efeps1cr3n=1
      efeps1deltn=1
      efeps2a1pa2n=1
      efeps2cr3n=1
      efeps2deltn=1
      efa1pa2cr3n=1
      efa1pa2deltn=1
      efcr3deltn=1
      goto 8745
      endif
      delt=delta0
      IZ=0
        nizospin=iz
      nmax1=nmax-1
c  calculez valorile polinoamelor Laguerre
      call clag(nmax1)
       print*,'a1,a2,a1pa2',a1,a2,a1pa2
      CALL ENERG(A0,Z0,A2E,Z2E,NMAX,IZ,NRNIV,EN,NFIT,
     C                  IAS,infin,nrspin,nnlel)
          print*,'a0',a0,'dupa energ'
      print*,'iesit din energ'
           do iop=1,nrniv
           enp(iop)=en(iop)
           nrspinp(iop)=nrspin(iop)
           infinp(iop)=infin(iop)
c          print*,'enp(iop)',enp(iop),iop
           enrgd(iop)=en(iop)
           nspind(iop)=nrspin(iop)
           infnpd(iop)=infin(iop)
           nordinp(iop)=iop
           enddo
            call donare(nrniv,en,infin,nrspin,nordinp)
            do iop=1,nrniv
            jfinp(iop)=infin(iop)
            enddo
        
      print*,'iesit din donare'
            nnfer=nrniv+1-(z0+1)/2
            enfermip=en(nnfer)
          if(en(nrniv).lt.-300.d0)then                   
           poten=1000.d0
            goto 9092
          endif
      idp=iz
      gamma=1.4
      hw0=41*a0**(-.33333333333d0)
         nrimpa=0
            print*,'a0',a0,'intru emedienou'
      call emedienou(idp,a0,z0,nrniv,gamma,en,hw0,
     c       nrimpa,eimp,
     c       u,se,du,alambd)
            print*,'idp,z0,nrniv,gamma,hw0,nrimpa'
            print*,idp,z0,nrniv,gamma,hw0,nrimpa
      call pairnou(idp,a0,z0,en,
     c   gamma,dp,sollam,soldel)
          hw00=hw0

c avem grija de valorile ce provin din common/partnrs si vor ajunge
c prin common/partnr2 in calculul maselor efective
          do kk=1,nrniv
          eonr2(kk)=eonri(kk)
          enddo
          do kk=1,10
          eimpar2(kk)=eimpari(kk)
          enddo
          nrnive2=nrnivei
          nrimpare2=nrimparei
          nrocup2=nrocupi

            gmass=gpair
            nrnivcmass=nrnivcpair
            iiermass=iierpair
      du=hw0*du
      dp=hw0*dp
              dcc=du+dp
         print*,'a0,z0,du0,dp0',a0,z0,du,dp
      sollam=sollam*hw0
      soldel=soldel*hw0
      sollam0=sollam
      soldel0=soldel
      sollam0p=sollam
      soldel0p=soldel
      sollamp=sollam
      soldelp=soldel

            dupro=du
            dppro=dp
      dv=du+dp


           if(initfin.eq.1)then
           do iop=1,nrniv
           nrasimptp(iop)=infinp(iop)
           enddo
           else
           do iop=1,nrniv
           infinp(iop)=nrasimptp(iop)
           enddo
           endif

c corectii de paturi si imprechere pe fragment 1
           do iop=1,nrniv
           enin(iop)=enp(iop)
           if(infinp(iop).eq.2)enin(iop)=1000
           enddo
            call donare(nrniv,enin,infin,nrspin,nordin)
      a01=a0-a2final+0.1
      z01=z0-z2final+0.1
      call emedienou(idp,a01,z01,nrniv,gamma,enin,hw0,
     c       nrimpa,eimp,
     c       u1,se1,du1,alambd1)
      call pairnou(idp,a01,z01,enin,
     c   gamma,dp1,sollam1,soldel1)
               dufp=hw0*du1
               dpfp=hw0*dp1
            sollam1=sollam1*hw0
            soldel1=soldel1*hw0
            sollam10=sollam1
            soldel10=soldel1
            sollam10p=sollam1
            soldel10p=soldel1
c corectii de paturi si imprechere pe fragment 2
           do iop=1,nrniv
           enin(iop)=enp(iop)
           if(infinp(iop).eq.1)enin(iop)=1000
           enddo
            call donare(nrniv,enin,infin,nrspin,nordin)
      a01=a2final+0.1
      z01=z2final+0.1
      call emedienou(idp,a01,z01,nrniv,gamma,enin,hw0,
     c       nrimpa,eimp,
     c       u2,se2,du2,alambd2)
      call pairnou(idp,a01,z01,enin,
     c   gamma,dp2,sollam2,soldel2)
               duep=hw0*du2
               dpep=hw0*dp2




c           do iomg=1,nmax
c           nelmx=(nmax+1-iomg)*(nmax+2-iomg)/2
c           do ll=1,nelmx

c     do iz=1,16 
c          xcx=15.+0.5*delta
c     z=-xcx+(iz-1)*2.-3
c     do iro=1,16  
c     ro=-15+(iro-1)*2. 

c     vall=csi(z,ro,nmax,iomg,ll)


!c            print*,'intru in densitati'
c     call densitateuni(z,ro,sollam,soldel,sollam1,soldel1,rrsc,
c    c                  nmax,densit)
c     write(97,*)z,ro,vall    
c     enddo
c     enddo
c         enddo
c         enddo




cccc Fac modificar astfel incat sa calculez elementele de matrice gelm(300,300)
c pentru interactia de pairing numai pentru primele 300 de nivele in ordine crescatoare
      rrsc=delatsciz
            nmaxy=nmax
               n1=1
      print*,'intru intepridid'

           do i=1,100
           distanta=0.2*(i-1)

               call intetripid(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt)
          print*,'ies intepridid'
          print*,zt,'zt',zp1,zp2,'zp1 zp2'
               call intetripid2(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt2)
          print*,zt2,'zt2',zp1,zp2,'zp1 zp2'
          write(97,*)distanta,zt,zt2
           enddo


                stop
              do iu1=1,25
              do iu2=1,325
              do iu3=1,325
              vecintermp(iu1,iu2,iu3)=vecpro(iu1,iu2,iu3)
              enddo
              enddo
              enddo









!!!!!!!!!! Daca nmas=0 calculez cu model cranking
!!!!!!!!!! Daca nmas nu e zero calcules cu model GOA
!!!!!!!!!! Cand utilizez GOA sa nu uit sa pun common/dvtr/
!!!!!!!!! in programul principal cu derivatele coordonatelor
!!!!!!!!! generalizate
!!!!!!1!!      common/dvtr/dereps1,dereps2,dercr3,dera
       nmas=1 
          if(nmas.eq.0)then
      call AMEF(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1p,efeps2p,
     c                    efa1pa2p,efcr3p,efdeltp,efeps1eps2p,
     c  efeps1a1pa2p,efeps1cr3p,efeps1deltp,efeps2a1pa2p,efeps2cr3p,
     c  efeps2deltp,efa1pa2cr3p,efa1pa2deltp,efcr3deltp)
          else
          if(nmas.eq.1)then
      call GOA(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1p,efeps2p,
     c                    efa1pa2p,efcr3p,efdeltp,efeps1eps2p,
     c  efeps1a1pa2p,efeps1cr3p,efeps1deltp,efeps2a1pa2p,efeps2cr3p,
     c  efeps2deltp,efa1pa2cr3p,efa1pa2deltp,efcr3deltp)
      sumeps1p=sumeps1
      sumeps2p=sumeps2
      suma1pa2p=suma1pa2
      sumcr3p=sumcr3
      sumdeltp=sumdelt
      sumeps1eps2p=sumeps1eps2
      sumeps1a1pa2p=sumeps1a1pa2
      sumeps1cr3p=sumeps1cr3
      sumeps1deltp=sumeps1delt
      sumeps2a1pa2=sumeps2a1pa2
      sumeps2cr3p=sumeps2cr3
      sumeps2deltp=sumeps2delt
      suma1pa2cr3p=suma1pa2cr3
      suma1pa2deltp=suma1pa2delt
      sumcr3deltp=sumcr3delt
          else
      call JPG(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1p,efeps2p,
     c                    efa1pa2p,efcr3p,efdeltp,efeps1eps2p,
     c  efeps1a1pa2p,efeps1cr3p,efeps1deltp,efeps2a1pa2p,efeps2cr3p,
     c  efeps2deltp,efa1pa2cr3p,efa1pa2deltp,efcr3deltp)
          endif
          endif

cccccccc      NMAX=14
      IZ=1
        nizospin=iz
      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      CALL ENERG(A0,Z0,A2E,Z2E,NMAX,IZ,NRNIV,EN,NFIT,
     C                  IAS,infin,nrspin,nnlel)
           do iop=1,nrniv
           enn(iop)=en(iop)
           nrspinn(iop)=nrspin(iop)
           infinn(iop)=infin(iop)
           enrgd(iop)=en(iop)
           nspind(iop)=nrspin(iop)
           infnpd(iop)=infin(iop)
           nordinn(iop)=iop
           enddo
      call donare(nrniv,en,infin,nrspin,nordinn)
            
            do iop=1,nrniv
            jfinn(iop)=infin(iop)
            enddo
            nnfer=nrniv+1-(a0+1-z0)/2
            enfermin=en(nnfer)
c              do ki=1,nrniv
c              iki=nrniv+1-ki
c              if(iki.ge.nnfer)then
c              enimp(ki+1)=en(ki)
c              else
c              enimp(ki)=en(ki)
c              endif
c              enddo
c              en(1)=1000.
c              do ki=2,nrniv
c              en(ki)=enimp(ki)
c              enddo
      idp=iz
      gamma=1.4
      hw0=41*a0**(-.33333333333d0)
      call emedienou(idp,a0,z0,nrniv,gamma,en,hw0,
     c       nrimpa,eimp,
     c       u,se,du,alambd)
      call pairnou(idp,a0,z0,en,
     c   gamma,dp,sollam,soldel)
          hw00=hw0

c avem grija de valorile ce provin din common/partnrs si vor ajunge
c prin common/partnr2 in calculul maselor efective
          do kk=1,nrniv
          eonr2(kk)=eonri(kk)
          enddo
          do kk=1,10
          eimpar2(kk)=eimpari(kk)
          enddo
          nrnive2=nrnivei
          nrimpare2=nrimparei
          nrocup2=nrocupi

            gmass=gpair
            nrnivcmass=nrnivcpair
            iiermass=iierpair
      du=hw0*du
      dp=hw0*dp
      sollam=sollam*hw0
      soldel=soldel*hw0
      sollam0=sollam
      soldel0=soldel
      sollam0n=sollam
      soldel0n=soldel
         print*,'a0,z0,du0,dp0',a0,z0,du,dp

       print*,'lambda delta in pair',sollam,soldel
       print*,'corectii de paturi in MeV du si dp',du,dp
c          stop
      dv=dv+du+dp


           if(initfin.eq.1)then
           do iop=1,nrniv
           nrasimptn(iop)=infinn(iop)
           enddo
           else
           do iop=1,nrniv
           infinn(iop)=nrasimptn(iop)
           enddo
           endif


c corectii de paturi si imprechere pe fragment 1
           do iop=1,nrniv
           enin(iop)=enn(iop)
           if(infinn(iop).eq.2)enin(iop)=1000
           enddo
            call donare(nrniv,enin,infin,nrspin,nordin)
      a01=a0-a2final+0.1
      z01=z0-z2final+0.1
      call emedienou(idp,a01,z01,nrniv,gamma,enin,hw0,
     c       nrimpa,eimp,
     c       u1,se1,du1,alambd1)
      call pairnou(idp,a01,z01,enin,
     c   gamma,dp1,sollam1,soldel1)
               dufn=hw0*du1
               dpfn=hw0*dp1
            sollam1=sollam1*hw0
            soldel1=soldel1*hw0
            sollam10=sollam1
            soldel10=soldel1
            sollam10n=sollam1
            soldel10n=soldel1
c corectii de paturi si imprechere pe fragment 2
           do iop=1,nrniv
           enin(iop)=enn(iop)
           if(infinn(iop).eq.1)enin(iop)=1000
           enddo
            call donare(nrniv,enin,infin,nrspin,nordin)
      a01=a2final+0.1
      z01=z2final+0.1
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!emedie'
      call emedienou(idp,a01,z01,nrniv,gamma,enin,hw0,
     c       nrimpa,eimp,
     c       u2,se2,du2,alambd2)
       print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!pairnou'
      call pairnou(idp,a01,z01,enin,
     c   gamma,dp2,sollam2,soldel2)
               duen=hw0*du2
               dpen=hw0*dp2

        amin=delatsciz-1.0000000000000005
        amax=delatsciz-1.000000005


        efep=dufp+dpfp+duep+dpep+dufn+dpfn+duen+dpen
        if(a2final.lt.18.d0)efep=dufp+dpfp+dufn+dpfn
             if(delta.ge.delatsciz)dv=efep
             if(delta.gt.amin.and.delta.lt.amax)then
c                   ddvv=dv
c       ain=dcos(3.14d0*(delta-amin)/(amax-amin))
c       ain=1-((ain+1.d0)/2.)
c       dv=ddvv+(efep-ddvv)*ain
                     ddvv=dv
             dv=ddvv+(efep-ddvv)*(delta-amin)/(amax-amin)
             endif


      do iz=1,121
      do iro=1,101
           xcx=15.+0.5*delta
      z=-xcx+(iz-1)*0.3-3
      ro=-15+(iro-1)*0.3
c     call densitateuni(z,ro,sollam,soldel,sollam1,soldel1,rrsc,
c    c                  nmax,densit)
c     write(98,*)z,ro,densit
      enddo
      enddo

c           do iomg1=1,nmax
c           nelmx1=(nmax+1-iomg1)*(nmax+2-iomg1)/2
c           do ll1=1,nelmx1
c           do iomg2=1,nmax
c           nelmx2=(nmax+1-iomg2)*(nmax+2-iomg2)/2
c           do ll2=1,nelmx2
c           print*,'l1 l2',ll1,ll2
c     call intetrip(sollam,soldel,sollam1,soldel1,
c    c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt)

c     write(97,*)iomg1,iomg2,ll1,ll2,zp1,zp2,zt
c     enddo
c     enddo
c         enddo
c         enddo

      rrsc=delatsciz
              if(initfin.eq.1)goto 8435
            nmaxy=nmax
            n100=100
               n1=2
               call intetripid(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt)
          print*,'ies intepridid'
          print*,zt,'zt',zp1,zp2,'zp1 zp2'
                stop

              do iu1=1,25
              do iu2=1,325
              do iu3=1,325
              vecintermn(iu1,iu2,iu3)=vecpro(iu1,iu2,iu3)
              enddo
              enddo
              enddo




            sollam0=sollam0p
            soldel0=soldel0p
            sollam10=sollam10p
            soldel10=soldel10p

              do iu1=1,25
              do iu2=1,325
              do iu3=1,325
              vecpro(iu1,iu2,iu3)=vecintermp(iu1,iu2,iu3)
              enddo
              enddo
              enddo

c                   lord=0
c                   do iopp=1,n100
c                   iop=nrniv+1-iop
c                   if(jfinp(iop).eq.2)then
c                   lord=lord+1
c                   nordin(lord)=nordinp(lord)
c                   do ip=iop-lord+1,nrniv-lord
c                   nordinp(iop-lord+1)=nordinp(iop+1-lord+1)
c                   enddo
c                   endif
c                   enddo
c                   do iop=1,lord
c                   nordinp(nrniv-lord+iop)=nordin(iop) 
c                   enddo

                  do idist=0,0
                  suma=0
                  distanta=idist*0.2      
                     nll1=0
                   do nr1=1,n100
                   jjj1=nrniv+1-nr1
                   nrr1=nordinp(jjj1)
                   ig1(nr1)=jfinp(jjj1)
                     nll2=0


                     nl1=0
            do iomg1=1,nmax
            nelmx1=(nmax+1-iomg1)*(nmax+2-iomg1)/2
            do ll1=1,nelmx1
                     nl1=nl1+1




      call intetrip(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt)
          gelm(nr1,nr2)=zt
          gelm1(nr1,nr2)=zp1
          gelm2(nr1,nr2)=zp2
      print*,'nr1,nr2,zt,zp1,zp2',nr1,nr2,zt,zp1,zp2
      print*,'iomg1,iomg2,ll1,ll2',iomg1,iomg2,ll1,ll2
               suma=suma+zt
      enddo
      enddo
            
            enddo
           write(97,*)distanta,suma
            enddo



              do nr1=1,n100
              do nr2=1,n100
      if(nr2.lt.nr1)then
c     write(97,*)nr2,nr1,ig2(nr2),ig1(nr1),gelm1(nr2,nr1),gelm2(nr2,nr1)
c    c ,gelm(nr2,nr1)
      else
c     write(97,*)nr1,nr2,ig1(nr1),ig2(nr2),gelm1(nr1,nr2),gelm2(nr1,nr2)
c    c ,gelm(nr1,nr2)
      endif
              enddo
              enddo



            sollam0=sollam0n
            soldel0=soldel0n
            sollam10=sollam10n
            soldel10=soldel10n

              do iu1=1,25
              do iu2=1,325
              do iu3=1,325
              vecpro(iu1,iu2,iu3)=vecintermn(iu1,iu2,iu3)
              enddo
              enddo
              enddo


c                   lord=0
c                   do iop=1,nrniv
c                   if(jfinn(iop).eq.1)then
c                   lord=lord+1
c                   nordin(lord)=nordinn(lord)
c                   do ip=iop-lord+1,nrniv-lord
c                   nordinn(iop-lord+1)=nordinn(iop+1-lord+1)
c                   enddo
c                   endif
c                   enddo
c                   do iop=1,lord
c                   nordinn(nrniv-lord+iop)=nordin(iop) 
c                   enddo



                     nll1=0
                   do nr1=1,n100
                   jjj1=nrniv+1-nr1
                   nrr1=nordinn(jjj1)
                   ig1(nr1)=jfinn(jjj1)
                     nll2=0
                   do nr2=1,n100
           if(nr2.lt.nr1)goto 6677
                   jjj2=nrniv+1-nr2
                   nrr2=nordinn(jjj2)
                   ig2(nr2)=jfinn(jjj2)


                     nl1=0
            do iomg1=1,nmax
            nelmx1=(nmax+1-iomg1)*(nmax+2-iomg1)/2
            do ll1=1,nelmx1
                     nl1=nl1+1

            if(nl1.eq.nrr1)then


                     nl2=0
            do iomg2=1,nmax
            nelmx2=(nmax+1-iomg2)*(nmax+2-iomg2)/2
            do ll2=1,nelmx2
                     nl2=nl2+1

            if(nl2.eq.nrr2)then
             n1=2
c     call intetrip(sollam,soldel,sollam1,soldel1,
c    c    sollam2,soldel2,nmax,iomg1,iomg2,ll1,ll2,zp1,zp2,zt)
          gelm(nr1,nr2)=zt
          gelm1(nr1,nr2)=zp1
          gelm2(nr1,nr2)=zp2
                     endif

      enddo
      enddo
             endif
          enddo
          enddo

6677   continue
            enddo
            enddo



              do nr1=1,n100
              do nr2=1,n100
      if(nr2.lt.nr1)then
c     write(97,*)nr2,nr1,ig2(nr2),ig1(nr1),gelm1(nr2,nr1),gelm2(nr2,nr1)
c    c ,gelm(nr2,nr1)
      else
c     write(97,*)nr1,nr2,ig1(nr1),ig2(nr2),gelm1(nr1,nr2),gelm2(nr1,nr2)
c    c ,gelm(nr1,nr2)
      endif
              enddo
              enddo

8435    continue


          if(nmas.eq.0)then
      call AMEF(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1n,efeps2n,
     c                    efa1pa2n,efcr3n,efdeltn,efeps1eps2n,
     c  efeps1a1pa2n,efeps1cr3n,efeps1deltn,efeps2a1pa2n,efeps2cr3n,
     c  efeps2deltn,efa1pa2cr3n,efa1pa2deltn,efcr3deltn)
      else
         if(nmas.eq.1)then
      call GOA(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1n,efeps2n,
     c                    efa1pa2n,efcr3n,efdeltn,efeps1eps2n,
     c  efeps1a1pa2n,efeps1cr3n,efeps1deltn,efeps2a1pa2n,efeps2cr3n,
     c  efeps2deltn,efa1pa2cr3n,efa1pa2deltn,efcr3deltn)
      else
      call JPG(sollam,soldel,iz,a0,z0,nmax,nrniv,efeps1n,efeps2n,
     c                    efa1pa2n,efcr3n,efdeltn,efeps1eps2n,
     c  efeps1a1pa2n,efeps1cr3n,efeps1deltn,efeps2a1pa2n,efeps2cr3n,
     c  efeps2deltn,efa1pa2cr3n,efa1pa2deltn,efcr3deltn)
      endif
      endif
9092  continue
8745  continue
          if(nmas.eq.0)then
      efeps1=efeps1p+efeps1n
      efeps2=efeps2p+efeps2n
      efa1pa2=efa1pa2p+efa1pa2n
      efcr3=efcr3p+efcr3n
      efdelt=efdeltp+efdeltn
      efeps1eps2=efeps1eps2p+efeps1eps2n
      efeps1a1pa2=efeps1a1pa2p+efeps1a1pa2n
      efeps1cr3=efeps1cr3p+efeps1cr3n
      efeps1delt=efeps1deltp+efeps2deltn
      efeps2a1pa2=efeps2a1pa2p+efeps2a1pa2n
      efeps2cr3=efeps2cr3p+efeps2cr3n
      efeps2delt=efeps2deltp+efeps2deltn
      efa1pa2cr3=efa1pa2cr3p+efa1pa2cr3n
      efa1pa2delt=efa1pa2deltp+efa1pa2deltn
      efcr3delt=efcr3deltp+efcr3deltn
         else
         if(nmas.eq.1)then
      efeps1=(sumeps1p+sumeps1)**2*efeps1p*efeps1n/
     c (sumeps1p**2*efeps1p+sumeps1**2*efeps1n)
      efeps2=(sumeps2p+sumeps2)**2*efeps2p*efeps2n/
     c (sumeps2p**2*efeps2p+sumeps2**2*efeps2n)
      efa1pa2=(suma1pa2+suma1pa2)**2*efa1pa2p*efa1pa2n/
     c (suma1pa2p**2*efa1pa2p+suma1pa2**2*efa1pa2n)
      efcr3=(sumcr3p+sumcr3)**2*efcr3p*efcr3n/
     c (sumcr3p**2*efcr3p+sumcr3**2*efcr3n)
      efdelt=(sumdeltp+sumdelt)**2*efdeltp*efdeltn/
     c (sumdeltp**2*efdeltp+sumdelt**2*efdeltn)
      efeps1eps2=(sumeps1eps2p+sumeps1eps2)**2*efeps1eps2p*efeps1eps2n/
     c (sumeps1eps2p**2*efeps1eps2p+sumeps1eps2**2*efeps1eps2n)
      efeps1a1pa2=(sumeps1a1pa2p+sumeps1a1pa2)**2*
     c efeps1a1pa2p*efeps1a1pa2n/
     c (sumeps1a1pa2p**2*efeps1a1pa2p+sumeps1a1pa2*efeps1a1pa2n)
      efeps1cr3=(sumeps1cr3p+sumeps1cr3)**2*efeps1cr3p*efeps1cr3n/
     c (sumeps1cr3p**2*efeps1cr3p+sumeps1cr3**2*efeps1cr3n)
      efeps1delt=(sumeps1deltp+sumeps1delt)**2*efeps1deltp*efeps2deltn/
     c (sumeps1deltp**2*efeps1deltp+sumeps1delt**2*efeps1deltn)
      efeps2a1pa2=(sumeps2a1pa2p+sumeps2a1pa2)**2*
     c efeps2a1pa2p*efeps2a1pa2n/
     c (sumeps2a1pa2p**2*efeps2a1pa2p+sumeps2a1pa2**2*efeps2a1pa2n)
      efeps2cr3=(sumeps2cr3p+sumeps2cr3)**2*efeps2cr3p*efeps2cr3n/
     c (sumeps2cr3**2*efeps2cr3p+sumeps2cr3**2*efeps2cr3n)
      efeps2delt=(sumeps2deltp+sumeps2delt)**2*efeps2deltp*efeps2deltn/
     c (sumeps2deltp**2*efeps2deltp+sumeps2delt**2*efeps2deltn)
      efa1pa2cr3=(suma1pa2cr3p+suma1pa2cr3)**2*efa1pa2cr3p*efa1pa2cr3n/
     c (suma1pa2cr3p**2*efa1pa2cr3p+suma1pa2cr3**2*efa1pa2cr3n)
      efa1pa2delt=(suma1pa2deltp+suma1pa2delt)**2*
     c efa1pa2deltp*efa1pa2deltn/
     c (suma1pa2deltp**2*efa1pa2deltp+suma1pa2delt**2*efa1pa2deltn)
      efcr3delt=(sumcr3deltp+sumcr3delt)**2*efcr3deltp*efcr3deltn/
     c (sumcr3deltp**2*efcr3deltp+sumcr3delt*efcr3deltn)
      efeps1=0
      efeps2=0
      efa1pa2=0
      efcr3=0
      efeps1eps2=0
      efeps1a1pa2=0
      efeps1cr3=0
      efeps1delt=0
      efeps2a1pa2=0
      efeps2cr3=0
      efeps2delt=0
      efa1pa2cr3=0
      efa1pa2delt=0
      efcr3delt=0
         else
      efeps1=efeps1p+efeps1n
      efeps2=efeps2p+efeps2n
      efa1pa2=efa1pa2p+efa1pa2n
      efcr3=efcr3p+efcr3n
      efdelt=efdeltp+efdeltn
      efeps1eps2=efeps1eps2p+efeps1eps2n
      efeps1a1pa2=efeps1a1pa2p+efeps1a1pa2n
      efeps1cr3=efeps1cr3p+efeps1cr3n
      efeps1delt=efeps1deltp+efeps2deltn
      efeps2a1pa2=efeps2a1pa2p+efeps2a1pa2n
      efeps2cr3=efeps2cr3p+efeps2cr3n
      efeps2delt=efeps2deltp+efeps2deltn
      efa1pa2cr3=efa1pa2cr3p+efa1pa2cr3n
      efa1pa2delt=efa1pa2deltp+efa1pa2deltn
      efcr3delt=efcr3deltp+efcr3deltn
         endif
         endif
       return
      end















        subroutine schema
c am luat programul coul.f si incerc sa fac un program
c pentru schema de nivele neutronice


c  am inlocuit deja ws

      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,delta
      common/dltw/dltw
      common/r116/R0
      dimension ccc(7)
      DIMENSION EN(2925)
      dimension infin(2925),nrspin(2925)
      R0=1.16d0
      ccc(1)=-0.07
      ccc(2)=-0.05
      ccc(3)=-0.02
      ccc(4)=-0.01
      ccc(5)=-0.007
      ccc(6)=-0.004
      ccc(7)=-0.002
      mmmm=0


      nmax=24
      nmax1=nmax-1
c  calculez valorile polinoamelor Laguerre
      call clag(nmax1)

     

      do j= 1,70!70 ! corespunde delta
      mmmm=mmmm+1

      delta=.01d0 +(j-1)*0.3


      if(delta.lt.7.d0)then
      cr3=-0.04
      else
      if(delta.gt.15.d0)then
      cr3=0.1
      else
      cr3=-0.04+(0.1+0.04)*(delta-7.d0)/(15.d0-7.d0)
      endif
      endif

                    cr3=1.0

      dcr3=dabs(cr3)
      a0=235
      z0=92
      a2e=100
      z2e=40
      eps0=0
      eps1f=0.6
      eps2f=0.5
           c  a0=26
           c  z0=13
           c  a2e=12
           c  z2e=6

             eps0=0
             eps1f=0
             eps2f=0
c      eps1=0.5
c      eps2=0.5
c      a1pa2=1.5
                 eps1f=0
                 eps2f=0
c excentricitate zero=sfera
c parametrizez eps1,eps2,a1pa2 (excentricitati si raport intre semi-axele a1 si a2)
      r0=1.27 ! valoarea lui r0 pentru calculul microscopic
      delta0=delta
      call vefinel(a0,z0,a2e,z2e,eps0,eps1f,eps2f,cr3,delta,
     c eps1,eps2,a1pa2,a2m,z2m)
      delta0=delta
      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      delt=delta0
      write(88,*)a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,delta
      r0=1.27 ! valoarea lui r0 pentru calculul microscopic

c      NMAX=4!14!7!5
      IZ=1
      CALL ENERG(A0,Z0,A2E,Z2E,NMAX,IZ,NRNIV,EN,NFIT,
     C                  IAS,infin,nrspin,nnlel)
        print*,'iesit din energ'
      idp=iz
      gamma=1.4
           do ii=1,nnlel
           write(51,*)ii,delta,cr3,en(ii),nrspin(ii)
           enddo


      enddo
      end








           subroutine grafic(a0,z0,a2e,z2e,eps1,eps2,a1pa2,cr3,delta)  
c am luat programul coul.f si incerc sa fac un program
c care sa dea valori pentru potentialul wood-saxon (fac niste grafice)
!!      double precision function ws(z,rho)
!!c ne da valorile potentialului ws functie de z si ro (coordonatele cilindrice)
      implicit double precision (a-h,o-z)
      common/r116/R0
      common/adancimi/vpo1,vpo2
      common/adifz/adiz
      common/paramel/a1v,b1v,a2v,b2v,c1,c2,c3,x1p,
     c               x2p,r3,u1,u2,ro3,s,delt
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2 ! valori calculate in NENGLS
      COMMON/ALF12/ALF1,ALF2
      common/oadancime/v000
      common/CKST/CKST
      common/lasciziune/delatsciz
      r0=1.27 ! valoarea lui r0 pentru calculul microscopic
c      call vefinel(a0,z0,a2e,z2e,eps0,eps1f,eps2f,cr3,delta,
c     c eps1,eps2,a1pa2,a2m,z2m)
      delta0=delta
c     call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      print*,'a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2'
      print*,a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2
          i=1
      call parw_s(i,a0,z0,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      adiz=adz
      v00=vn0
      v000=vn0
          a1=a0-a2e
          z1=z0-z2e
      call parw_s(i,a1,z1,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      v01=vn0
          w2=a2e
          z2=z2e
      call parw_s(i,w2,z2,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      v02=vn0


      DPAR=delatsciz
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      IF(DELTA.GT.DPAR)THEN
      VP1=V01
      VP2=V02
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      VP1=V00
      VP2=V00
      ELSE
      VP1=V00+(DELTA-DPAR07)/(DPAR-DPAR07)*(V01-V00)
      VP2=V00+(DELTA-DPAR07)/(DPAR-DPAR07)*(V02-V00)
      ENDIF
      ENDIF
      vpo1=vp1
      vpo2=vp2


      CTw=c2/dabs(c1)
      cx1=-DELT/(1.D0+CTw)
      cx2=CT*dabs(cx1)
C DISTANTA DINTRE PLANUL DE INTERSECTIE A PARAMETRIZARII SEMISIMETRICE  SI
C PLANUL DE INTERSECTIE A PARAMTERIZARII ASIMETRICE
      Z0P=cx1-c1
      CKST=41.d0
      r0=1.16
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      W2RO=W0*R0/b2
      ALF1ro=DSQRT(MPH*W1ro)*1.D-15
      ALF2ro=DSQRT(MPH*W2ro)*1.D-15
      HWRO1=CKST*R0/b1
      HWZ1=b1/a1*HWRO1
      HWRO2=CKST*R0/b2
      HWZ2=b2/a2*HWRO2

      w1z1=w0*r0/a1
      w2z2=w0*r0/a2
      ALF1=DSQRT(MPH*W1z1)*1.D-15
      ALF2=DSQRT(MPH*W2z2)*1.D-15

          difr=50./50.
c     do l1=1,51
c     do l2=1,51
c     z=-25+difr*(l1-1)
c     ro=25-difr*(l2-1)
      do iz=1,121
           xcx=15.+0.5*delta
      z=-xcx+(iz-1)*0.3-3
      do iro=1,101
      ro=-15+(iro-1)*0.3


      woods=ws(z,ro)
c     difere= dif_ws_os(z,ro)
c     call dervVx(z,ro,derddz,derdro)
c     am=dsqrt(derddz**2+derdro**2)
      vdens=rhodens(z,ro)

      write(93,*)z,woods,vdens
      enddo
      enddo

                                return

                                end








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    MASE EFECTIVE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine AMEF(efermi,delaaa,iz,a0,z0,nmax,nrniv,efeps1,efeps2,
     c                    efa1pa2,efcr3,efdelt,efeps1eps2,
     c  efeps1a1pa2,efeps1cr3,efeps1delt,efeps2a1pa2,efeps2cr3,
     c  efeps2delt,efa1pa2cr3,efa1pa2delt,efcr3delt)
c iz=0 protoni : iz=1 neutroni
c deltag=gapu-ul, efermi este energia Fermi
      implicit double precision (a-h,o-z)
c calculez masele efective
c aici am elementele de matrice obtinute din nenglscumef
c primul indice iomg, al doilea este lel si al treilea este nel
c emasef sunt energii functie de iomg si lel
      dimension en(2925),infin(2925),nrspin(2925),nordin(2925)
      dimension u(25,325),v(25,325)
      common/nivselect2/eselect(2925),nnnnn
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/ginter/ginte ! constanta de interactie
      common/encuazi/eq(25,325)
      common/conmef/peps1,peps2,pa1pa2,pcr3,pdelt,peps1eps2,
     c peps1a1pa2,peps1cr3,peps1delt,peps2a1pa2,peps2cr3,
     c peps2delt,pa1pa2cr3,pa1pa2delt,pcr3delt
c urmatoarele doua common-uri provin din function bcs1nou
      common/nrnivcmass/g,nrnivc,iier
        common/partnr2/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
c dimesiunea maxima pentru nel si lel este unde iomg incepe de la 1
C PANA AICI NEL=LEL=(Nmax+1-IOMG)*(Nmax+2-IOMG)/2

c creez vectorii cu energii, numerele lel si iomg care vor fi
c pusi in ordine
      i=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do l=1,nelmax
      i=i+1
      en(i)=emasef(iomg,l)
      infin(i)=l
      nrspin(i)=iomg
      enddo
      enddo
      call donare(nrniv,en,infin,nrspin,nordin)

c      nrocup=(a0-z0)/2
c      if(iz.eq.0)nrocup=z0/2
c      ginte=27
c      if(iz.eq.0)ginte=35
c      numniv=30
c      numniv2=-numniv+1
      nnnnn=0
      do i=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+i
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
c      nnnnn=2*numniv
c      call bcssolu(a0,delaaa,efermi)

c          print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
          print*,delaaa,efermi
      energiamax=eselect(1)
      energiamin=eselect(nnnnn)
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j=1,nelmax
      eq(iomg,j)=dsqrt((emasef(iomg,j)-efermi)**2+delaaa**2)
      if(emasef(iomg,j).gt.energiamax)then
      u(iomg,j)=1.d0
      v(iomg,j)=0.d0
      else
      if(emasef(iomg,j).lt.energiamin)then
      u(iomg,j)=0.d0
      v(iomg,j)=1.d0
      else
      u(iomg,j)=dsqrt(0.5d0*(1.d0+(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      v(iomg,j)=dsqrt(0.5d0*(1.d0-(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      endif
      endif
c      print*,'iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)'
c      print*,iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)
      enddo
      enddo
                  

      sumeps1=0
      sumeps2=0
      suma1pa2=0
      sumcr3=0
      sumdelt=0
      sumeps1eps2=0
      sumeps1a1pa2=0
      sumeps1cr3=0
      sumeps1delt=0
      sumeps2a1pa2=0
      sumeps2cr3=0
      sumeps2delt=0
      suma1pa2cr3=0
      suma1pa2delt=0
      sumcr3delt=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      do j2=1,nelmax
!!!!!!!!!!!!!!!!  TERMENII DIAGONALI
      sumeps1=sumeps1+elmateps1(iomg,j1,j2)*elmateps1(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps2=sumeps2+elmateps2(iomg,j1,j2)*elmateps2(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      suma1pa2=suma1pa2+elmata1pa2(iomg,j1,j2)*elmata1pa2(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumcr3=sumcr3+elmatcr3(iomg,j1,j2)*elmatcr3(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumdelt=sumdelt+elmatdelt(iomg,j1,j2)*elmatdelt(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
c      print*,'iomg,j1,j2,sumdelt,elmatdelt(iomg,j1,j2),eq(iomg,j1)'
c      print*,'eq(iomg,j2),u(iomg,j1),u(iomg,j2),v(iomg,j1),v(iomg,j2)'
c      print*,iomg,j1,j2,sumdelt,elmatdelt(iomg,j1,j2),eq(iomg,j1)
c      print*,eq(iomg,j2),u(iomg,j1),u(iomg,j2),v(iomg,j1),v(iomg,j2)
!!!!!!!!!!!!!!!!!! TERMENII NEDIAGONALI
      sumeps1eps2=sumeps1eps2+elmateps1(iomg,j1,j2)*
     c  elmateps2(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps1a1pa2=sumeps1a1pa2+elmateps1(iomg,j1,j2)*
     c  elmata1pa2(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps1cr3=sumeps1cr3+elmateps1(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps1delt=sumeps1delt+elmateps1(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps2a1pa2=sumeps2a1pa2+elmateps2(iomg,j1,j2)*
     c  elmata1pa2(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps2cr3=sumeps2cr3+elmateps2(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumeps2delt=sumeps2delt+elmateps2(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      suma1pa2cr3=suma1pa2cr3+elmata1pa2(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      suma1pa2delt=suma1pa2delt+elmata1pa2(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      sumcr3delt=sumcr3delt+elmatcr3(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)/
     c  (eq(iomg,j1)+eq(iomg,j2))**3*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
      enddo
      enddo
      enddo      
          print*,delaaa,efermi
       call pereche(delaaa,efermi,nmax) 
      efeps1=2*sumeps1  +peps1
      efeps2=2*sumeps2  +peps2
      efa1pa2=2*suma1pa2  +pa1pa2
      efcr3=2*sumcr3  +pcr3
      efdelt=2*sumdelt  +pdelt
      efeps1eps2=2*sumeps1eps2  +peps1eps2
      efeps1a1pa2=2*sumeps1a1pa2  +peps1a1pa2
      efeps1cr3=2*sumeps1cr3  +peps1cr3
      efeps1delt=2*sumeps1delt  +peps1delt
      efeps2a1pa2=2*sumeps2a1pa2  +peps2a1pa2
      efeps2cr3=2*sumeps2cr3  +peps2cr3
      efeps2delt=2*sumeps2delt  +peps2delt
      efa1pa2cr3=2*suma1pa2cr3 +pa1pa2cr3
      efa1pa2delt=2*suma1pa2delt  +pa1pa2delt
      efcr3delt=2*sumcr3delt  +pcr3delt


      return
      end


      subroutine pereche(delta,efermi,nmax)
c calculez corectia P la masele efective
      implicit double precision (a-h,o-z)
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/encuazi/eq(25,325)
c contributii mase efective
      common/conmef/peps1,peps2,pa1pa2,pcr3,pdelt,peps1eps2,
     c peps1a1pa2,peps1cr3,peps1delt,peps2a1pa2,peps2cr3,
     c peps2delt,pa1pa2cr3,pa1pa2delt,pcr3delt
c 
             gap=delta
      cceps1=0
      cdeps1=0
      cceps2=0
      cdeps2=0
      cca1pa2=0
      cda1pa2=0
      cccr3=0
      cdcr3=0
      ccdelt=0
      cddelt=0
      ca=0
      cb=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      ca=ca+2*delta/eq(iomg,j1)**3
      cb=cb+2*(emasef(iomg,j1)-efermi)/eq(iomg,j1)**3
      cceps1=cceps1+delta*elmateps1(iomg,j1,j1)/eq(iomg,j1)**3
      cdeps1=cdeps1+(emasef(iomg,j1)-efermi)*elmateps1(iomg,j1,j1)/
     c                 eq(iomg,j1)**3

      cceps2=cceps2+delta*elmateps2(iomg,j1,j1)/eq(iomg,j1)**3
      cdeps2=cdeps2+(emasef(iomg,j1)-efermi)*elmateps2(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      cca1pa2=cca1pa2+delta*elmata1pa2(iomg,j1,j1)/eq(iomg,j1)**3
      cda1pa2=cda1pa2+(emasef(iomg,j1)-efermi)*elmata1pa2(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      cccr3=cccr3+delta*elmatcr3(iomg,j1,j1)/eq(iomg,j1)**3
      cdcr3=cdcr3+(emasef(iomg,j1)-efermi)*elmatcr3(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      ccdelt=ccdelt+delta*elmatdelt(iomg,j1,j1)/eq(iomg,j1)**3
      cddelt=cddelt+(emasef(iomg,j1)-efermi)*elmatdelt(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      enddo
      enddo
         patrat=ca**2+cb**2
       dferdeps1=(ca*cceps1+cb*cdeps1)/patrat
       dgapdeps1=(cb*cceps1-ca*cdeps1)/patrat
       dferdeps2=(ca*cceps2+cb*cdeps2)/patrat
       dgapdeps2=(cb*cceps2-ca*cdeps2)/patrat
       dferda1pa2=(ca*cca1pa2+cb*cda1pa2)/patrat
       dgapda1pa2=(cb*cca1pa2-ca*cda1pa2)/patrat
       dferdcr3=(ca*cccr3+cb*cdcr3)/patrat
       dgapdcr3=(cb*cccr3-ca*cdcr3)/patrat
       dferddelt=(ca*ccdelt+cb*cddelt)/patrat
       dgapddelt=(cb*ccdelt-cb*cddelt)/patrat



      peps1=0
      peps2=0
      pa1pa2=0
      pcr3=0
      pdelt=0
      peps1eps2=0
      peps1a1pa2=0
      peps1cr3=0
      peps1delt=0
      peps2a1pa2=0
      peps2cr3=0
      peps2delt=0
      pa1pa2cr3=0
      pa1pa2delt=0
      pcr3delt=0
      

      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      peps11=delta**2*dferdeps1*dferdeps1+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdeps1+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdeps1+
     +  dferdeps1*dgapdeps1)-
     -  gap**2*(dferdeps1*elmateps1(iomg,j1,j1)+
     +  dferdeps1*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmateps1(iomg,j1,j1)+dgapdeps1*elmateps1(iomg,j1,j1))
      peps1=peps1+peps11/eq(iomg,j1)**5/4.
      peps21=delta**2*dferdeps2*dferdeps2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapdeps2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapdeps2+
     +  dferdeps2*dgapdeps2)-
     -  gap**2*(dferdeps2*elmateps2(iomg,j1,j1)+
     +  dferdeps2*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmateps2(iomg,j1,j1)+dgapdeps2*elmateps2(iomg,j1,j1))
      peps2=peps2+peps21/eq(iomg,j1)**5/4.
      pa1pa21=delta**2*dferda1pa2*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapda1pa2+
     +  dferda1pa2*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmata1pa2(iomg,j1,j1))
      pa1pa2=pa1pa2+pa1pa21/eq(iomg,j1)**5/4.
      pcr31=delta**2*dferdcr3*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdcr3*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdcr3*dgapdcr3+
     +  dferdcr3*dgapdcr3)-
     -  gap**2*(dferdcr3*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmatcr3(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdcr3*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmatcr3(iomg,j1,j1))
      pcr3=pcr3+pcr31/eq(iomg,j1)**5/4.
      pdelt1=delta**2*dferddelt*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapddelt*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferddelt*dgapddelt+
     +  dferddelt*dgapddelt)-
     -  gap**2*(dferddelt*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmatdelt(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapddelt*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmatdelt(iomg,j1,j1))
      pdelt=pdelt+pdelt1/eq(iomg,j1)**5/4.
      peps1eps21=delta**2*dferdeps1*dferdeps2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdeps2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdeps2+
     +  dferdeps2*dgapdeps1)-
     -  gap**2*(dferdeps1*elmateps2(iomg,j1,j1)+
     +  dferdeps2*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmateps2(iomg,j1,j1)+dgapdeps2*elmateps1(iomg,j1,j1))
      peps1eps2=peps1eps2+peps1eps21/eq(iomg,j1)**5/4.
      peps1a1pa21=delta**2*dferdeps1*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapda1pa2+
     +  dferda1pa2*dgapdeps1)-
     -  gap**2*(dferdeps1*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmateps1(iomg,j1,j1))
      peps1a1pa2=peps1a1pa2+peps1a1pa21/eq(iomg,j1)**5/4.
      peps1cr31=delta**2*dferdeps1*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdcr3+
     +  dferdcr3*dgapdeps1)-
     -  gap**2*(dferdeps1*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmateps1(iomg,j1,j1))
      peps1cr3=peps1cr3+peps1cr31/eq(iomg,j1)**5/4.
      peps1delt1=delta**2*dferdeps1*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapddelt+
     +  dferddelt*dgapdeps1)-
     -  gap**2*(dferdeps1*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmateps1(iomg,j1,j1))
      peps1delt=peps1delt+peps1delt1/eq(iomg,j1)**5/4.
      peps2a1pa21=delta**2*dferdeps2*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapda1pa2+
     +  dferda1pa2*dgapdeps2)-
     -  gap**2*(dferdeps2*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmateps2(iomg,j1,j1))
      peps2a1pa2=peps2a1pa2+peps2a1pa21/eq(iomg,j1)**5/4.
      peps2cr31=delta**2*dferdeps2*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapdcr3+
     +  dferdcr3*dgapdeps2)-
     -  gap**2*(dferdeps2*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmateps2(iomg,j1,j1))
      peps2cr3=peps2cr3+peps2cr31/eq(iomg,j1)**5/4.
      peps2delt1=delta**2*dferdeps2*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapddelt+
     +  dferddelt*dgapdeps2)-
     -  gap**2*(dferdeps2*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmateps2(iomg,j1,j1))
      peps2delt=peps2delt+peps2delt1/eq(iomg,j1)**5/4.
      pa1pa2cr31=delta**2*dferda1pa2*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapdcr3+
     +  dferdcr3*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmata1pa2(iomg,j1,j1))
      pa1pa2cr3=pa1pa2cr3+pa1pa2cr31/eq(iomg,j1)**5/4.
      pa1pa2delt1=delta**2*dferda1pa2*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapddelt+
     +  dferddelt*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmata1pa2(iomg,j1,j1))
      pa1pa2delt=pa1pa2delt+pa1pa2delt1/eq(iomg,j1)**5/4.
      pcr3delt1=delta**2*dferdcr3*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdcr3*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdcr3*dgapddelt+
     +  dferddelt*dgapdcr3)-
     -  gap**2*(dferdcr3*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmatcr3(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdcr3*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmatcr3(iomg,j1,j1))
      pcr3delt=pcr3delt+pcr3delt1/eq(iomg,j1)**5/4.

      enddo
      enddo

      return
      end      



      subroutine GOA(efermi,delaaa,iz,a0,z0,nmax,nrniv,efeps1,efeps2,
     c                    efa1pa2,efcr3,efdelt,efeps1eps2,
     c  efeps1a1pa2,efeps1cr3,efeps1delt,efeps2a1pa2,efeps2cr3,
     c  efeps2delt,efa1pa2cr3,efa1pa2delt,efcr3delt)
c iz=0 protoni : iz=1 neutroni
c deltag=gapu-ul, efermi este energia Fermi
      implicit double precision (a-h,o-z)
c calculez masele efective
c aici am elementele de matrice obtinute din nenglscumef
c primul indice iomg, al doilea este lel si al treilea este nel
c emasef sunt energii functie de iomg si lel
      dimension en(2925),infin(2925),nrspin(2925),nordin(2925)
      dimension u(25,325),v(25,325)
      dimension pjkdelt(25,325,325),pjkcr3(25,325,325),
     c pjkeps1(25,325,325),pjkeps2(25,325,325),pjka1pa2(25,325,325)
      dimension pjkdelta(25,325,325)
      common/nivselect2/eselect(2925),nnnnn
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/ginter/ginte ! constanta de interactie
      common/encuazi/eq(25,325)
      common/conmef/peps1,peps2,pa1pa2,pcr3,pdelt,peps1eps2,
     c peps1a1pa2,peps1cr3,peps1delt,peps2a1pa2,peps2cr3,
     c peps2delt,pa1pa2cr3,pa1pa2delt,pcr3delt
c urmatoarele doua common-uri provin din function bcs1nou
      common/nrnivcmass/g,nrnivc,iier
        common/partnr2/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/derparam/dferdeps1,dgapdeps1,dferdeps2,dgapdeps2,
     c dferda1pa2,dgapda1pa2,dferdcr3,dgapdcr3,dferddelt,dgapddelt
      common/sumg/sumeps1,sumeps2,suma1pa2,sumcr3,sumdelt,sumeps1eps2,
     c sumeps1a1pa2,sumeps1cr3,sumeps1delt,sumeps2a1pa2,sumeps2cr3,
     c sumeps2delt,suma1pa2cr3,suma1pa2delt,sumcr3delt

      common/dvtr/dereps1,dereps2,dercr3,dera
c dimesiunea maxima pentru nel si lel este unde iomg incepe de la 1
C PANA AICI NEL=LEL=(Nmax+1-IOMG)*(Nmax+2-IOMG)/2

c creez vectorii cu energii, numerele lel si iomg care vor fi
c pusi in ordine
      i=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do l=1,nelmax
      i=i+1
      en(i)=emasef(iomg,l)
      infin(i)=l
      nrspin(i)=iomg
      enddo
      enddo
      call donare(nrniv,en,infin,nrspin,nordin)

c      nrocup=(a0-z0)/2
c      if(iz.eq.0)nrocup=z0/2
c      ginte=27
c      if(iz.eq.0)ginte=35
c      numniv=30
c      numniv2=-numniv+1
      nnnnn=0
      do i=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+i
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
c      nnnnn=2*numniv
c      call bcssolu(a0,delaaa,efermi)

c          print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
          print*,delaaa,efermi
      energiamax=eselect(1)
      energiamin=eselect(nnnnn)
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j=1,nelmax
      eq(iomg,j)=dsqrt((emasef(iomg,j)-efermi)**2+delaaa**2)
      if(emasef(iomg,j).gt.energiamax)then
      u(iomg,j)=1.d0
      v(iomg,j)=0.d0
      else
      if(emasef(iomg,j).lt.energiamin)then
      u(iomg,j)=0.d0
      v(iomg,j)=1.d0
      else
      u(iomg,j)=dsqrt(0.5d0*(1.d0+(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      v(iomg,j)=dsqrt(0.5d0*(1.d0-(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      endif
      endif
c      print*,'iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)'
c      print*,iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)
      enddo
      enddo
                  
      call pereche2(delaaa,efermi,nmax) 
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      do j2=1,nelmax
      pjkdelt(iomg,j1,j2)=-elmatdelt(iomg,j1,j2)*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))/
     c  (eq(iomg,j1)+eq(iomg,j2))
      pjkcr3(iomg,j1,j2)=-elmatcr3(iomg,j1,j2)*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))/
     c  (eq(iomg,j1)+eq(iomg,j2))
      pjkeps1(iomg,j1,j2)=-elmateps1(iomg,j1,j2)*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))/
     c  (eq(iomg,j1)+eq(iomg,j2))
      pjkeps2(iomg,j1,j2)=-elmateps2(iomg,j1,j2)*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))/
     c  (eq(iomg,j1)+eq(iomg,j2))
      pjka1pa2(iomg,j1,j2)=-elmata1pa2(iomg,j1,j2)*
     c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))/
     c  (eq(iomg,j1)+eq(iomg,j2))
      if(j1.eq.j2)then
      pjkdelt(iomg,j1,j1)=pjkdelt(iomg,j1,j1)+
     c  0.5/eq(iomg,j1)**2*(delaaa*dferddelt+
     c  (emasef(iomg,j1)-efermi)*dgapddelt)
      pjkcr3(iomg,j1,j1)=pjkcr3(iomg,j1,j1)+
     c  0.5/eq(iomg,j1)**2*(delaaa*dferdcr3+
     c  (emasef(iomg,j1)-efermi)*dgapdcr3)
      pjkeps1(iomg,j1,j1)=pjkeps1(iomg,j1,j1)+
     c  0.5/eq(iomg,j1)**2*(delaaa*dferdeps1+
     c  (emasef(iomg,j1)-efermi)*dgapdeps1)
      pjkeps2(iomg,j1,j1)=pjkeps2(iomg,j1,j1)+
     c  0.5/eq(iomg,j1)**2*(delaaa*dferdeps2+
     c  (emasef(iomg,j1)-efermi)*dgapdeps2)
      pjka1pa2(iomg,j1,j1)=pjka1pa2(iomg,j1,j1)+
     c  0.5/eq(iomg,j1)**2*(delaaa*dferda1pa2+
     c  (emasef(iomg,j1)-efermi)*dgapda1pa2)
      endif
      pjkdelta(iomg,j1,j2)=pjkdelt(iomg,j1,j2)+
     c pjkcr3(iomg,j1,j2)*dercr3+pjkeps1(iomg,j1,j2)*dereps1+
     c pjkeps2(iomg,j1,j2)*dereps2+pjka1pa2(iomg,j1,j2)*dera
      enddo
      enddo
      enddo

      sumeps1=0
      sumeps1p=0
      sumeps2=0
      sumeps2p=0
      suma1pa2=0
      suma1pa2p=0
      sumcr3=0
      sumcr3p=0
      sumdelt=0
      sumdeltp=0
      sumeps1eps2=0
      sumeps1eps2p=0
      sumeps1a1pa2=0
      sumeps1a1pa2p=0
      sumeps1cr3=0
      sumeps1cr3p=0
      sumeps1delt=0
      sumeps1deltp=0
      sumeps2a1pa2=0
      sumeps2a1pa2p=0
      sumeps2cr3=0
      sumeps2cr3p=0
      sumeps2delt=0
      sumeps2deltp=0
      suma1pa2cr3=0
      suma1pa2cr3p=0
      suma1pa2delt=0
      suma1pa2deltp=0
      sumcr3delt=0
      sumcr3deltp=0

      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      do j2=1,nelmax
c        if(v(iomg,j1).eq.1.d0.or.v(iomg,j1).eq.0.d0.or.
c     c  v(iomg,j2).eq.1.d0.or.v(iomg,j2).eq.0.d0)goto 5151
!!!!!!!!!!!!!!!!  TERMENII DIAGONALI
      sumeps1=sumeps1+pjkeps1(iomg,j1,j2)*pjkeps1(iomg,j2,j1)
      sumeps1p=sumeps1p+pjkeps1(iomg,j1,j2)*pjkeps1(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps2=sumeps2+pjkeps2(iomg,j1,j2)*pjkeps1(iomg,j2,j1)
      sumeps2p=sumeps2p+pjkeps2(iomg,j1,j2)*pjkeps2(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      suma1pa2=suma1pa2+pjka1pa2(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)
      suma1pa2p=suma1pa2p+pjka1pa2(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumcr3=sumcr3+pjkcr3(iomg,j1,j2)*pjkcr3(iomg,j2,j1)
      sumcr3p=sumcr3p+pjkcr3(iomg,j1,j2)*pjkcr3(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumdelt=sumdelt+pjkdelta(iomg,j1,j2)*pjkdelta(iomg,j2,j1)
      sumdeltp=sumdeltp+pjkdelta(iomg,j1,j2)*pjkdelta(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))


!!!!!!!!!!!!!!!!!! TERMENII NEDIAGONALI
      sumeps1eps2=sumeps1eps2+pjkeps1(iomg,j1,j2)*pjkeps2(iomg,j2,j1)
      sumeps1eps2p=sumeps1eps2p+pjkeps1(iomg,j1,j2)*pjkeps2(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps1a1pa2=sumeps1a1pa2+pjkeps1(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)
      sumeps1a1pa2p=sumeps1a1pa2p+
     c pjkeps1(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps1cr3=sumeps1cr3+pjkeps1(iomg,j1,j2)*pjkcr3(iomg,j2,j1)
      sumeps1cr3p=sumeps1cr3p+pjkeps1(iomg,j1,j2)*pjkcr3(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps1delt=sumeps1delt+pjkeps1(iomg,j1,j2)*pjkdelt(iomg,j2,j1)
      sumeps1deltp=sumeps1deltp+pjkeps1(iomg,j1,j2)*pjkdelt(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps2a1pa2=sumeps2a1pa2+pjkeps2(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)
      sumeps2a1pa2p=sumeps2a1pa2p+
     c pjkeps2(iomg,j1,j2)*pjka1pa2(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps2cr3=sumeps2cr3+pjkeps2(iomg,j1,j2)*pjkcr3(iomg,j2,j1)
      sumeps2cr3p=sumeps2cr3p+pjkeps2(iomg,j1,j2)*pjkcr3(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumeps2delt=sumeps2delt+pjkeps2(iomg,j1,j2)*pjkdelt(iomg,j2,j1)
      sumeps2deltp=sumeps2deltp+pjkeps2(iomg,j1,j2)*pjkdelt(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      suma1pa2cr3=suma1pa2cr3+pjka1pa2(iomg,j1,j2)*pjkcr3(iomg,j2,j1)
      suma1pa2cr3p=suma1pa2cr3p+pjka1pa2(iomg,j1,j2)*pjkcr3(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      suma1pa2delt=suma1pa2delt+pjka1pa2(iomg,j1,j2)*pjkdelt(iomg,j2,j1)
      suma1pa2deltp=suma1pa2deltp+
     c pjka1pa2(iomg,j1,j2)*pjkdelt(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
      sumcr3delt=sumcr3delt+pjkcr3(iomg,j1,j2)*pjkdelt(iomg,j2,j1)
      sumcr3deltp=sumcr3deltp+
     c pjkcr3(iomg,j1,j2)*pjkdelt(iomg,j2,j1)*
     c (eq(iomg,j1)+eq(iomg,j2))
c5151  continue
      enddo
      enddo
      enddo      
          print*,delaaa,efermi
c      call pereche(delaaa,efermi,nmax) 
      efeps1=2*sumeps1**2/sumeps1p ! +peps1
      efeps2=2*sumeps2**2/sumeps2p ! +peps2
      efa1pa2=2*suma1pa2**2/suma1pa2p ! +pa1pa2
      efcr3=2*sumcr3**2/sumcr3p ! +pcr3
      efdelt=2*sumdelt**2/sumdeltp ! +pdelt
      efeps1eps2=2*sumeps1eps2**2/sumeps1eps2p ! +peps1eps2
      efeps1a1pa2=2*sumeps1a1pa2**2/sumeps1a1pa2p ! +peps1a1pa2
      efeps1cr3=2*sumeps1cr3**2/sumeps1cr3p ! +peps1cr3
      efeps1delt=2*sumeps1delt**2/sumeps1deltp ! +peps1delt
      efeps2a1pa2=2*sumeps2a1pa2**2/sumeps2a1pa2p ! +peps2a1pa2
      efeps2cr3=2*sumeps2cr3**2/sumeps2cr3p ! +peps2cr3
      efeps2delt=2*sumeps2delt**2/sumeps2deltp ! +peps2delt
      efa1pa2cr3=2*suma1pa2cr3**2/suma1pa2cr3p !+pa1pa2cr3
      efa1pa2delt=2*suma1pa2delt**2/suma1pa2deltp ! +pa1pa2delt
      efcr3delt=2*sumcr3delt**2/sumcr3deltp  !+pcr3delt


      return
      end


      subroutine pereche2(delta,efermi,nmax)
c calculez corectia P la masele efective
      implicit double precision (a-h,o-z)
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/encuazi/eq(25,325)
c contributii mase efective
      common/conmef/peps1,peps2,pa1pa2,pcr3,pdelt,peps1eps2,
     c peps1a1pa2,peps1cr3,peps1delt,peps2a1pa2,peps2cr3,
     c peps2delt,pa1pa2cr3,pa1pa2delt,pcr3delt
      common/derparam/dferdeps1,dgapdeps1,dferdeps2,dgapdeps2,
     c dferda1pa2,dgapda1pa2,dferdcr3,dgapdcr3,dferddelt,dgapddelt
c 
             gap=delta
      cceps1=0
      cdeps1=0
      cceps2=0
      cdeps2=0
      cca1pa2=0
      cda1pa2=0
      cccr3=0
      cdcr3=0
      ccdelt=0
      cddelt=0
      ca=0
      cb=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      ca=ca+2*delta/eq(iomg,j1)**3
      cb=cb+2*(emasef(iomg,j1)-efermi)/eq(iomg,j1)**3
      cceps1=cceps1+delta*elmateps1(iomg,j1,j1)/eq(iomg,j1)**3
      cdeps1=cdeps1+(emasef(iomg,j1)-efermi)*elmateps1(iomg,j1,j1)/
     c                 eq(iomg,j1)**3

      cceps2=cceps2+delta*elmateps2(iomg,j1,j1)/eq(iomg,j1)**3
      cdeps2=cdeps2+(emasef(iomg,j1)-efermi)*elmateps2(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      cca1pa2=cca1pa2+delta*elmata1pa2(iomg,j1,j1)/eq(iomg,j1)**3
      cda1pa2=cda1pa2+(emasef(iomg,j1)-efermi)*elmata1pa2(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      cccr3=cccr3+delta*elmatcr3(iomg,j1,j1)/eq(iomg,j1)**3
      cdcr3=cdcr3+(emasef(iomg,j1)-efermi)*elmatcr3(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      ccdelt=ccdelt+delta*elmatdelt(iomg,j1,j1)/eq(iomg,j1)**3
      cddelt=cddelt+(emasef(iomg,j1)-efermi)*elmatdelt(iomg,j1,j1)/
     c                 eq(iomg,j1)**3
      enddo
      enddo
         patrat=ca**2+cb**2
       dferdeps1=(ca*cceps1+cb*cdeps1)/patrat
       dgapdeps1=(cb*cceps1-ca*cdeps1)/patrat
       dferdeps2=(ca*cceps2+cb*cdeps2)/patrat
       dgapdeps2=(cb*cceps2-ca*cdeps2)/patrat
       dferda1pa2=(ca*cca1pa2+cb*cda1pa2)/patrat
       dgapda1pa2=(cb*cca1pa2-ca*cda1pa2)/patrat
       dferdcr3=(ca*cccr3+cb*cdcr3)/patrat
       dgapdcr3=(cb*cccr3-ca*cdcr3)/patrat
       dferddelt=(ca*ccdelt+cb*cddelt)/patrat
       dgapddelt=(cb*ccdelt-cb*cddelt)/patrat


               return
      peps1=0
      peps2=0
      pa1pa2=0
      pcr3=0
      pdelt=0
      peps1eps2=0
      peps1a1pa2=0
      peps1cr3=0
      peps1delt=0
      peps2a1pa2=0
      peps2cr3=0
      peps2delt=0
      pa1pa2cr3=0
      pa1pa2delt=0
      pcr3delt=0
      

      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      peps11=delta**2*dferdeps1*dferdeps1+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdeps1+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdeps1+
     +  dferdeps1*dgapdeps1)-
     -  gap**2*(dferdeps1*elmateps1(iomg,j1,j1)+
     +  dferdeps1*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmateps1(iomg,j1,j1)+dgapdeps1*elmateps1(iomg,j1,j1))
      peps1=peps1+peps11/eq(iomg,j1)**5/4.
      peps21=delta**2*dferdeps2*dferdeps2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapdeps2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapdeps2+
     +  dferdeps2*dgapdeps2)-
     -  gap**2*(dferdeps2*elmateps2(iomg,j1,j1)+
     +  dferdeps2*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmateps2(iomg,j1,j1)+dgapdeps2*elmateps2(iomg,j1,j1))
      peps2=peps2+peps21/eq(iomg,j1)**5/4.
      pa1pa21=delta**2*dferda1pa2*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapda1pa2+
     +  dferda1pa2*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmata1pa2(iomg,j1,j1))
      pa1pa2=pa1pa2+pa1pa21/eq(iomg,j1)**5/4.
      pcr31=delta**2*dferdcr3*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdcr3*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdcr3*dgapdcr3+
     +  dferdcr3*dgapdcr3)-
     -  gap**2*(dferdcr3*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmatcr3(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdcr3*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmatcr3(iomg,j1,j1))
      pcr3=pcr3+pcr31/eq(iomg,j1)**5/4.
      pdelt1=delta**2*dferddelt*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapddelt*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferddelt*dgapddelt+
     +  dferddelt*dgapddelt)-
     -  gap**2*(dferddelt*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmatdelt(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapddelt*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmatdelt(iomg,j1,j1))
      pdelt=pdelt+pdelt1/eq(iomg,j1)**5/4.
      peps1eps21=delta**2*dferdeps1*dferdeps2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdeps2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdeps2+
     +  dferdeps2*dgapdeps1)-
     -  gap**2*(dferdeps1*elmateps2(iomg,j1,j1)+
     +  dferdeps2*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmateps2(iomg,j1,j1)+dgapdeps2*elmateps1(iomg,j1,j1))
      peps1eps2=peps1eps2+peps1eps21/eq(iomg,j1)**5/4.
      peps1a1pa21=delta**2*dferdeps1*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapda1pa2+
     +  dferda1pa2*dgapdeps1)-
     -  gap**2*(dferdeps1*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmateps1(iomg,j1,j1))
      peps1a1pa2=peps1a1pa2+peps1a1pa21/eq(iomg,j1)**5/4.
      peps1cr31=delta**2*dferdeps1*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapdcr3+
     +  dferdcr3*dgapdeps1)-
     -  gap**2*(dferdeps1*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmateps1(iomg,j1,j1))
      peps1cr3=peps1cr3+peps1cr31/eq(iomg,j1)**5/4.
      peps1delt1=delta**2*dferdeps1*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps1*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps1*dgapddelt+
     +  dferddelt*dgapdeps1)-
     -  gap**2*(dferdeps1*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmateps1(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps1*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmateps1(iomg,j1,j1))
      peps1delt=peps1delt+peps1delt1/eq(iomg,j1)**5/4.
      peps2a1pa21=delta**2*dferdeps2*dferda1pa2+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapda1pa2+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapda1pa2+
     +  dferda1pa2*dgapdeps2)-
     -  gap**2*(dferdeps2*elmata1pa2(iomg,j1,j1)+
     +  dferda1pa2*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmata1pa2(iomg,j1,j1)+dgapda1pa2*elmateps2(iomg,j1,j1))
      peps2a1pa2=peps2a1pa2+peps2a1pa21/eq(iomg,j1)**5/4.
      peps2cr31=delta**2*dferdeps2*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapdcr3+
     +  dferdcr3*dgapdeps2)-
     -  gap**2*(dferdeps2*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmateps2(iomg,j1,j1))
      peps2cr3=peps2cr3+peps2cr31/eq(iomg,j1)**5/4.
      peps2delt1=delta**2*dferdeps2*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdeps2*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdeps2*dgapddelt+
     +  dferddelt*dgapdeps2)-
     -  gap**2*(dferdeps2*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmateps2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdeps2*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmateps2(iomg,j1,j1))
      peps2delt=peps2delt+peps2delt1/eq(iomg,j1)**5/4.
      pa1pa2cr31=delta**2*dferda1pa2*dferdcr3+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapdcr3+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapdcr3+
     +  dferdcr3*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmatcr3(iomg,j1,j1)+
     +  dferdcr3*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmatcr3(iomg,j1,j1)+dgapdcr3*elmata1pa2(iomg,j1,j1))
      pa1pa2cr3=pa1pa2cr3+pa1pa2cr31/eq(iomg,j1)**5/4.
      pa1pa2delt1=delta**2*dferda1pa2*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapda1pa2*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferda1pa2*dgapddelt+
     +  dferddelt*dgapda1pa2)-
     -  gap**2*(dferda1pa2*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmata1pa2(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapda1pa2*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmata1pa2(iomg,j1,j1))
      pa1pa2delt=pa1pa2delt+pa1pa2delt1/eq(iomg,j1)**5/4.
      pcr3delt1=delta**2*dferdcr3*dferddelt+
     +  (emasef(iomg,j1)-efermi)**2*dgapdcr3*dgapddelt+
     +  delta*(emasef(iomg,j1)-efermi)*(dferdcr3*dgapddelt+
     +  dferddelt*dgapdcr3)-
     -  gap**2*(dferdcr3*elmatdelt(iomg,j1,j1)+
     +  dferddelt*elmatcr3(iomg,j1,j1))-
     -  delta*(emasef(iomg,j1)-efermi)*(dgapdcr3*
     *  elmatdelt(iomg,j1,j1)+dgapddelt*elmatcr3(iomg,j1,j1))
      pcr3delt=pcr3delt+pcr3delt1/eq(iomg,j1)**5/4.

      enddo
      enddo

      return
      end      



      subroutine JPG(efermi,delaaa,iz,a0,z0,nmax,nrniv,efeps1,efeps2,
     c                    efa1pa2,efcr3,efdelt,efeps1eps2,
     c  efeps1a1pa2,efeps1cr3,efeps1delt,efeps2a1pa2,efeps2cr3,
     c  efeps2delt,efa1pa2cr3,efa1pa2delt,efcr3delt)
c iz=0 protoni : iz=1 neutroni
c deltag=gapu-ul, efermi este energia Fermi
      implicit double precision (a-h,o-z)
c calculez masele efective cu metoda din JPG
c aici am elementele de matrice obtinute din nenglscumef
c primul indice iomg, al doilea este lel si al treilea este nel
c emasef sunt energii functie de iomg si lel
      dimension en(2925),infin(2925),nrspin(2925),nordin(2925)
      dimension u(25,325),v(25,325)
      dimension tk(25,325),Emm(25,325),Ejl(25,325,325)
      dimension deltam(25,325),efermim(25,325)
      dimension rho(25,325),akap(25,325)
      dimension u2(25,325),v2(25,325),eq2(25,325),eq22(25,325)
      dimension tk2(25,325),Ediftk(25,325,325)
      dimension rho2(25,325),akap2(25,325)
      dimension u22(25,325),v22(25,325),rho22(25,325)
      dimension akap22(25,325),tk22(25,325),edifmm(25,325)
      dimension ssumtk2(25,325,325)
      common/nivselect2/eselect(2925),nnnnn
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/ginter/ginte ! constanta de interactie
      common/ginter2/ginte2 ! constanta de interactie
      common/encuazi/eq(25,325)
      common/conmef/peps1,peps2,pa1pa2,pcr3,pdelt,peps1eps2,
     c peps1a1pa2,peps1cr3,peps1delt,peps2a1pa2,peps2cr3,
     c peps2delt,pa1pa2cr3,pa1pa2delt,pcr3delt
      common/derparam/dferdeps1,dgapdeps1,dferdeps2,dgapdeps2,
     c dferda1pa2,dgapda1pa2,dferdcr3,dgapdcr3,dferddelt,dgapddelt
c urmatoarele doua common-uri provin din function bcs1nou
      common/nrnivcmass/gnorm,nrnivc,iier
      common/hw0ptg/hw00
        common/partnr2/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
c prin common voi transmite energia de vibratie de ground
      common/envib0/energiav0
      common/bcsefermi/efermi0
      common/informare/informare
c dimesiunea maxima pentru nel si lel este unde iomg incepe de la 1
C PANA AICI NEL=LEL=(Nmax+1-IOMG)*(Nmax+2-IOMG)/2

c creez vectorii cu energii, numerele lel si iomg care vor fi
c pusi in ordine

             g=gnorm*hw00
c            print*,'g,gnorm,hw00',g,gnorm,hw00
             efermix=0
             efermi6x=0
             efermi66x=0
      i=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do l=1,nelmax
      i=i+1
      en(i)=emasef(iomg,l)
      infin(i)=l
      nrspin(i)=iomg
      enddo
      enddo
      call donare(nrniv,en,infin,nrspin,nordin)

      nrocup=(a0-z0)/2
      if(iz.eq.0)nrocup=z0/2
      ginte=g

      nnnnn=0
      nrnive=nrniv
c      print*,'nrnive,nrniv,nrnivc,nrocup'
c      print*,nrnive,nrniv,nrnivc,nrocup
      do i=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+i
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
c     
c     print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
c         print*,delaaa,efermi
c     print*,'nrnivc',nrnivc,' nnnnn',nnnnn
      call bcssolu(a0,delaaa,efermi)
           efermi0=efermi
c         print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
c         print*,delaaa,efermi
c       write(95,*)delaaa,efermi,nnnnn,' pentru e0'
                   
      energiamax=eselect(1)
      energiamin=eselect(nnnnn)
      energia0=-10000
c        write(23,*)delaaa,efermi


! scad numarul de nivele din interactie de perechi cu 2

      nnnnn=0
      nrnive=nrniv
c      print*,'nrnive,nrniv,nrnivc,nrocup'
c      print*,nrnive,nrniv,nrnivc,nrocup
       nrnivcm1=nrnivc-1
      do i=-nrnivcm1,nrnivcm1-1
      n1=nrnive+1-nrocup+i
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
c     
c     print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
c         print*,delaaa,efermi
c     print*,'nrnivc',nrnivc,' nnnnn',nnnnn
      call bcssolu(a0,delaaa,efermi)
           efermi0=efermi
c         print*,'!!!!!!!!!!!!!! !!!!!!!!!! delaaa,efermi'
c         print*,delaaa,efermi
c       write(95,*)delaaa,efermi,nnnnn,' pentru e0'
                   
      energiamax0=eselect(1)
      energiamin0=eselect(nnnnn)
      energia0=-10000
c        write(23,*)delaaa,efermi




          efermix=efermi


       E0=0.
       E00=0
       sumrho=0.
       delta0=0
      sumtk=0
      an0=0
      an02=0
         afermi=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j=1,nelmax
      eq(iomg,j)=dsqrt((emasef(iomg,j)-efermi)**2+delaaa**2)
c       write(23,*)'1',emasef(iomg,j),eq(iomg,j),iomg,j
      if(emasef(iomg,j).gt.energiamax0)then
      u(iomg,j)=1.d0
      v(iomg,j)=0.d0
      else
      if(emasef(iomg,j).lt.energiamin0)then
      u(iomg,j)=0.d0
      v(iomg,j)=1.d0
      else
      u(iomg,j)=dsqrt(0.5d0*(1.d0+(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      v(iomg,j)=dsqrt(0.5d0*(1.d0-(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      endif
      endif
         rho(iomg,j)=v(iomg,j)**2
         akap(iomg,j)=v(iomg,j)*u(iomg,j)
c      print*,'iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)'
c      print*,iomg,j,u(iomg,j),v(iomg,j),eq(iomg,j),emasef(iomg,j)
      tk(iomg,j)=0
       an0=an0+2*rho(iomg,j)

       if(emasef(iomg,j).ge.energiamin0.
     c   and.emasef(iomg,j).le.energiamax0)then
      tk(iomg,j)=2*rho(iomg,j)*(emasef(iomg,j)-efermix)-
     c 2*Ginte2*rho(iomg,j)**2+
     c akap(iomg,j)*delaaa*(rho(iomg,j)**2/akap(iomg,j)**2-1.)
       else
      tk(iomg,j)=0.         !2*rho(iomg,j)*(emasef(iomg,j)-efermix)
       endif

       if(emasef(iomg,j).le.energiamax0.and.
     c    emasef(iomg,j).ge.energia0)then     
       E0=E0+2*rho(iomg,j)*(emasef(iomg,j))
       E00=E00+2*rho(iomg,j)*(emasef(iomg,j))-tk(iomg,j)
       if(emasef(iomg,j).ge.energiamin0)then
c       E00=E00+2*rho(iomg,j)*(emasef(iomg,j)-efermix)-tk(iomg,j)
       delta0=delta0+akap(iomg,j)
       sumrho=sumrho+rho(iomg,j)**2
       an02=an02+2*rho(iomg,j)
       afermi=afermi+efermi*2*rho(iomg,j)
       endif
      sumtk=sumtk+tk(iomg,j)
       endif
      enddo
      enddo
           delta0=delta0*ginte2
c       print*,'E0,E00,delaaa,ginte2,sumrho'   
c       print*,E0,E00,delaaa,ginte2,sumrho       
      E0=E0-delta0**2/ginte2-ginte2*sumrho               
      E00=E00-delta0**2/ginte2-ginte2*sumrho

       etest=e00+sumtk
c     print*,'e0,e00,sumtk,etest',e0,e00,sumtk,etest
c       print*,'delaaa,delta0',delaaa,delta0,an0,an02


        
  
! stari de senioritate 2 pe acelasi nivel
      do iomg1=1,nmax
       nelmax1=(nmax+1-iomg1)*(nmax+2-iomg1)/2
      do j1=1,nelmax1
       Emm(iomg1,j1)=0
       deltamm=0
       sumrhomm=0 
       Edifmm(iomg1,j1)=0 
      sumtk22=0
c stari senioritate 2 pe acelasi nivel
       if(emasef(iomg1,j1).le.energiamax.and.
     c    emasef(iomg1,j1).ge.energiamin)then     

      i5=0
      do iomg5=1,nmax
      nelmax5=(nmax+1-iomg5)*(nmax+2-iomg5)/2
      do j5=1,nelmax5
         if(iomg5.eq.iomg1.and.j5.eq.j1)goto 9177
      i5=i5+1
      en(i5)=emasef(iomg5,j5)
      infin(i5)=j5
      nrspin(i5)=iomg5
9177       continue
      enddo
      enddo
      nrnivm1=nrniv-1
      call donare(nrnivm1,en,infin,nrspin,nordin)
      nnnnn=0
      do i9=-nrnivc+2,nrnivc-1!1 ! plus 1 apare 
                              ! fiindca sunt 1 nivele mai putin
      n1=(nrnive-1)+1-nrocup+i9
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
      call bcssolu(a0,delaaa6,efermi66)
c          print*,'!!!!!! 0 delaaa6,efermi66',delaaa6,efermi66,nnnnn
c          print*,'1111111 iomg1,iomg1,j1,j1',iomg1,iomg1,j1,j1
      efermi66x=efermi66   

      energiamax11=eselect(1)
      energiamin11=eselect(nnnnn)


      do iomg6=1,nmax
      nelmax6=(nmax+1-iomg6)*(nmax+2-iomg6)/2
      do j6=1,nelmax6
      eq22(iomg6,j6)=dsqrt((emasef(iomg6,j6)-efermi66)**2+delaaa6**2)
      if(emasef(iomg6,j6).gt.energiamax11)then
      u22(iomg6,j6)=1.d0
      v22(iomg6,j6)=0.d0
      else
      if(emasef(iomg6,j6).lt.energiamin11)then
      u22(iomg6,j6)=0.d0
      v22(iomg6,j6)=1.d0
      else
      u22(iomg6,j6)=dsqrt(0.5d0*(1.d0+(emasef(iomg6,j6)-efermi66)/
     c    eq22(iomg6,j6)))
      v22(iomg6,j6)=dsqrt(0.5d0*(1.d0-(emasef(iomg6,j6)-efermi66)/
     c    eq22(iomg6,j6)))
      endif
      endif
         rho22(iomg6,j6)=v22(iomg6,j6)**2
         akap22(iomg6,j6)=v22(iomg6,j6)*u22(iomg6,j6)
      tk22(iomg6,j6)=0
       if(emasef(iomg6,j6).le.energiamax11)then 
       if(emasef(iomg6,j6).ge.energiamin11)then          
      tk22(iomg6,j6)=2*rho22(iomg6,j6)*(emasef(iomg6,j6)-efermi66x)-
     c 2*Ginte2*rho22(iomg6,j6)**2+
     c akap22(iomg6,j6)*delaaa6*
     c (rho22(iomg6,j6)**2/akap22(iomg6,j6)**2-1.)
       else
      tk22(iomg6,j6)=2*rho22(iomg6,j6)*(emasef(iomg6,j6)-efermi66x)
       endif
       endif

       if(emasef(iomg6,j6).le.energiamax11)then 
       if(j1.eq.j6.and.iomg1.eq.iomg6)goto 22
       if(emasef(iomg6,j6).ge.energiamin11)then
       edifmm(iomg1,j1)=Emm(iomg1,j1)+
     c      2*rho22(iomg6,j6)*(emasef(iomg6,j6))
     c      -tk22(iomg6,j6)
       Emm(iomg1,j1)=Emm(iomg1,j1)+
     c      2*rho22(iomg6,j6)*(emasef(iomg6,j6))
       deltamm=deltamm+akap22(iomg6,j6)
       sumrhomm=sumrhomm+rho22(iomg6,j6)**2
       else
       edifmm(iomg1,j1)=Emm(iomg1,j1)+
     c      2*rho22(iomg6,j6)*(emasef(iomg6,j6))
     c      -tk22(iomg6,j6)
       Emm(iomg1,j1)=Emm(iomg1,j1)+
     c      2*rho22(iomg6,j6)*(emasef(iomg6,j6))
      endif
      sumtk22=sumtk22+tk22(iomg6,j6)
22    continue
       endif
     
      enddo
      enddo

c pana aici stari senioritate 2 pe acelasi nivel
        deltamm=ginte2*deltamm
        deltam(iomg1,j1)=deltamm
        efermim(iomg1,j1)=efermi66
        Emm(iomg1,j1)=Emm(iomg1,j1)+
     c      2*(emasef(iomg1,j1))-
     c      deltamm**2/ginte2
     c      -ginte2*sumrhomm

       edifmm(iomg1,j1)=edifmm(iomg1,j1)+
     c      2*(emasef(iomg1,j1))-
     c      deltamm**2/ginte2
     c      -ginte2*sumrhomm

        endif
        enddo
        enddo


              enmax=energiamax+5



      do iomg1=1,nmax
       nelmax1=(nmax+1-iomg1)*(nmax+2-iomg1)/2
      do j1=1,nelmax1                    
       do iomg2=1,nmax
       nelmax2=(nmax+1-iomg2)*(nmax+2-iomg2)/2
      do j2=1,nelmax2
c  aici trebuie sa gasesc u si v pentru fiecare stare de
c senioritate 2
        if(j2.ne.j1.and.iomg1.eq.iomg2)then 

      Ejl(iomg1,j1,j2)=0
      Ediftk(iomg1,j1,j2)=1000000000.
      deltajl=0
      sumrhojl=0
      sumtk2=0

       an=0
       an2=0
       afermi2=0

      if(emasef(iomg1,j1).gt.energiamax.and.
     c   emasef(iomg2,j2).gt.energiamax)goto 52
      if(emasef(iomg1,j1).lt.energiamin.and.
     c   emasef(iomg2,j2).lt.energiamin)goto 52
      if(emasef(iomg1,j1).gt.enmax)goto 52
      if(emasef(iomg2,j2).gt.enmax)goto 52
      

       if(emasef(iomg1,j1).le.energiamax)then
       if(emasef(iomg2,j2).le.energiamax)then 
        if(j2.ne.j1.or.iomg1.ne.iomg2)then 

      Ediftk(iomg1,j1,j2)=0
      i5=0
      do iomg5=1,nmax
      nelmax5=(nmax+1-iomg5)*(nmax+2-iomg5)/2
      do j5=1,nelmax5
         if((iomg5.eq.iomg1.and.j5.eq.j1).or.
     c       (iomg5.eq.iomg2.and.j5.eq.j2))goto 177
      i5=i5+1
      en(i5)=emasef(iomg5,j5)
      infin(i5)=j5
      nrspin(i5)=iomg5
177       continue
      enddo
      enddo
      nrnivm2=nrniv-2
      call donare(nrnivm2,en,infin,nrspin,nordin)
      nnnnn=0
      do i9=-nrnivc+2,nrnivc-1 ! plus 2 apare 
                              ! fiindca sunt 2 nivele mai putin
      n1=(nrnive-2)+1-nrocup+i9
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
      call bcssolu(a0,delaaa6,efermi6)

       efermi6x=efermi6
         

      energiamax1=eselect(1)
      energiamin1=eselect(nnnnn)

  
      sumtk2=0
      do iomg6=1,nmax
      nelmax6=(nmax+1-iomg6)*(nmax+2-iomg6)/2
      do j6=1,nelmax6
      eq2(iomg6,j6)=dsqrt((emasef(iomg6,j6)-efermi6)**2+delaaa6**2)
      if(emasef(iomg6,j6).gt.energiamax1)then
      u2(iomg6,j6)=1.d0
      v2(iomg6,j6)=0.d0
      else
      if(emasef(iomg6,j6).lt.energiamin1)then
      u2(iomg6,j6)=0.d0
      v2(iomg6,j6)=1.d0
      else
      u2(iomg6,j6)=dsqrt(0.5d0*(1.d0+(emasef(iomg6,j6)-efermi6)/
     c    eq2(iomg6,j6)))
      v2(iomg6,j6)=dsqrt(0.5d0*(1.d0-(emasef(iomg6,j6)-efermi6)/
     c    eq2(iomg6,j6)))
      endif
      endif
         rho2(iomg6,j6)=v2(iomg6,j6)**2
         akap2(iomg6,j6)=v2(iomg6,j6)*u2(iomg6,j6)
      tk2(iomg6,j6)=0


      
       if(emasef(iomg6,j6).le.energiamax1)then

       if(emasef(iomg6,j6).ge.energiamin1)then
      tk2(iomg6,j6)=2*rho2(iomg6,j6)*(emasef(iomg6,j6)-efermi6x)-
     c 2*Ginte2*rho2(iomg6,j6)**2+
     c akap2(iomg6,j6)*delaaa6*
     c (rho2(iomg6,j6)**2/akap2(iomg6,j6)**2-1.)
      else
      tk2(iomg6,j6)=0.    !!!!!!!!2*rho2(iomg6,j6)*(emasef(iomg6,j6)-efermi6x)
      endif

         if((iomg6.eq.iomg1.and.j6.eq.j1).or.
     c       (iomg6.eq.iomg2.and.j6.eq.j2))goto 1771
      if(emasef(iomg6,j6).ge.energiamin1)then
      Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))
      Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))-tk2(iomg6,j6)
         deltajl=deltajl+akap2(iomg6,j6)
         sumrhojl=sumrhojl+rho2(iomg6,j6)**2
      sumtk2=sumtk2+tk2(iomg6,j6)      
        an2=an2+2*rho2(iomg6,j6)
        afermi2=afermi2+efermi6*2*rho2(iomg6,j6)
      else
      Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))
      Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))-tk2(iomg6,j6)
      endif
        an=an+2*rho2(iomg6,j6)     
1771  continue

c      if(iomg1.eq.1.and.iomg2.eq.1.and.j1.eq.78.and.j2.eq.79)then
c      print*,'sumtk2,ediftk,iomg6,j6',
c     c sumtk2,Ediftk(iomg1,iomg2,j1,j2),iomg6,j6
c      endif

         endif
   
      enddo
      enddo
           sss=sumtk-sumtk2
c     print*,'an,an2,nfer2',an,an2,afermi2,' an0,an02,afermi',
c    c    an0,an02,afermi,' sumtk-sumtk2',sss
              ssumtk2(iomg1,j1,j2)=sumtk2
           deltajl=ginte2*deltajl

c      print*,'Ejl,Ediftk,deltajl,ginte2,sumrhojl'

c      print*,Ejl(iomg1,j1,j2),Ediftk(iomg1,j1,j2),deltajl,
c     c ginte2,sumrhojl

        Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)
     c   -deltajl**2/ginte2-
     c   ginte2*sumrhojl+emasef(iomg1,j1)+
     c   emasef(iomg2,j2)             
      
        Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)
     c   -deltajl**2/ginte2-
     c   ginte2*sumrhojl+(emasef(iomg1,j1))+
     c  (emasef(iomg2,j2))
c      print*,'Ejl,Ediftk,sumtk2',Ejl(iomg1,iomg2,j1,j2),
c     c  Ediftk(iomg1,iomg2,j1,j2),sumtk2

c       if(iomg1.eq.iomg2)then
c               cncn=Ediftk(iomg1,iomg2,j1,j2)+sumtk2
c               cnc=Ejl(iomg1,iomg2,j1,j2)-sumtk2+sumtk
cc      print*,'Ejl,ediftk,sumtk2,sumtk',Ejl(iomg1,j1,j2),
c     c Ediftk(iomg1,j1,j2),sumtk2,sumtk,iomg1,iomg2,j1,j2
c      print*,'Ediftk+sumtk2',cncn
c      print*,'Ejl(iomg1,j1,j2)-sum2tk+sumtk',cnc
cc      stop
c        endif         
           goto 52               
         endif
         endif
         endif




       an=0
       an2=0
       afermi2=0

        if(j2.ne.j1.and.iomg1.eq.iomg2)then 

      Ediftk(iomg1,j1,j2)=0
      i5=0
      do iomg5=1,nmax
      nelmax5=(nmax+1-iomg5)*(nmax+2-iomg5)/2
      do j5=1,nelmax5
         if((iomg5.eq.iomg1.and.j5.eq.j1).or.
     c       (iomg5.eq.iomg2.and.j5.eq.j2))goto 1177
      i5=i5+1
      en(i5)=emasef(iomg5,j5)
      infin(i5)=j5
      nrspin(i5)=iomg5
1177       continue
      enddo
      enddo
      nrnivm2=nrniv-2
      call donare(nrnivm2,en,infin,nrspin,nordin)
      nnnnn=0
      do i9=-nrnivc+2,nrnivc-1 ! plus 1 si -2 apar 
                              ! fiindca este 1 nivel mai putin
      n1=(nrnive-2)+1-nrocup+i9
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
      call bcssolu(a0,delaaa6,efermi6)

         efermi6x=efermi6         

      energiamax1=eselect(1)
      energiamin1=eselect(nnnnn)

  
      sumtk2=0
      do iomg6=1,nmax
      nelmax6=(nmax+1-iomg6)*(nmax+2-iomg6)/2
      do j6=1,nelmax6
      eq2(iomg6,j6)=dsqrt((emasef(iomg6,j6)-efermi6)**2+delaaa6**2)
      if(emasef(iomg6,j6).gt.energiamax1)then
      u2(iomg6,j6)=1.d0
      v2(iomg6,j6)=0.d0
      else
      if(emasef(iomg6,j6).lt.energiamin1)then
      u2(iomg6,j6)=0.d0
      v2(iomg6,j6)=1.d0
      else
      u2(iomg6,j6)=dsqrt(0.5d0*(1.d0+(emasef(iomg6,j6)-efermi6)/
     c    eq2(iomg6,j6)))
      v2(iomg6,j6)=dsqrt(0.5d0*(1.d0-(emasef(iomg6,j6)-efermi6)/
     c    eq2(iomg6,j6)))
      endif
      endif
         rho2(iomg6,j6)=v2(iomg6,j6)**2
         akap2(iomg6,j6)=v2(iomg6,j6)*u2(iomg6,j6)
      tk2(iomg6,j6)=0

      
       if(emasef(iomg6,j6).le.energiamax1)then

       if(emasef(iomg6,j6).ge.energiamin1)then
      tk2(iomg6,j6)=2*rho2(iomg6,j6)*(emasef(iomg6,j6)-efermi6x)-
     c 2*Ginte2*rho2(iomg6,j6)**2+
     c akap2(iomg6,j6)*delaaa6*
     c (rho2(iomg6,j6)**2/akap2(iomg6,j6)**2-1.)
      else
      tk2(iomg6,j6)=0.     !!!!!!!!!!2*rho2(iomg6,j6)*(emasef(iomg6,j6)-efermi6x)
      endif

         if((iomg6.eq.iomg1.and.j6.eq.j1).or.
     c       (iomg6.eq.iomg2.and.j6.eq.j2))goto 1671
      if(emasef(iomg6,j6).ge.energiamin1)then
      Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))
      Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))-tk2(iomg6,j6)
         deltajl=deltajl+akap2(iomg6,j6)
         sumrhojl=sumrhojl+rho2(iomg6,j6)**2
        
        an2=an2+2*rho2(iomg6,j6)

        afermi2=afermi2+efermi6*2*rho2(iomg6,j6)
      sumtk2=sumtk2+tk2(iomg6,j6)  
      else
      Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))
      Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)+
     c        2*rho2(iomg6,j6)*(emasef(iomg6,j6))-tk2(iomg6,j6)
      endif 
        an=an+2*rho2(iomg6,j6)   
1671  continue



         endif
   
      enddo
      enddo
              ssumtk2(iomg1,j1,j2)=sumtk2
            deltajl=ginte2*deltajl
                 sss=sumtk-sumtk2
c     print*,'an,an2,nfer',an,an2,afermi2, 
c    c    ' an0,an02,nfer',an0,an02,afermi,' sumtk-sumtk2',sss
c     print*,'emasef(iomgo,jo),emasef(iomgo,joo),energiamax,energiamin'
c     print*,emasef(iomg1,j1),emasef(iomg2,j2),energiamax,energiamin

c      if(iomg1.eq.iomg2.and.iomg1.eq.3.and.j2.eq.61.and.j1.eq.71)then
c     print*,'Ejl,Ediftk,deltajl,ginte2,sumrhojl'
c     print*,Ejl(iomg1,j1,j2),Ediftk(iomg1,j1,j2),deltajl,
c    c ginte2,sumrhojl
c     print*,'emasef(iomg1,j1),emasef(iomg2,j2)',emasef(iomg1,j1),
c    c emasef(iomg2,j2),iomg1,j1,iomg2,j2
c     endif

        Ejl(iomg1,j1,j2)=Ejl(iomg1,j1,j2)
     c   -deltajl**2/ginte2-
     c   ginte2*sumrhojl+emasef(iomg1,j1)+
     c   emasef(iomg2,j2)             
      
        Ediftk(iomg1,j1,j2)=Ediftk(iomg1,j1,j2)
     c   -deltajl**2/ginte2-
     c   ginte2*sumrhojl+(emasef(iomg1,j1))+
     c  (emasef(iomg2,j2))
c      print*,'Ejl,Ediftk,sumtk2',Ejl(iomg1,iomg2,j1,j2),
c     c  Ediftk(iomg1,j1,j2),sumtk2

c      if(iomg1.eq.iomg2.and.iomg1.eq.3.and.j2.eq.61.and.j1.eq.71)then
c      print*,'e0',e0
c              cncn=Ediftk(iomg1,j1,j2)+sumtk2
c              cnc=Ejl(iomg1,j1,j2)-sumtk2+sumtk
c      print*,'Ejl,ediftk,sumtk2,sumtk',Ejl(iomg1,j1,j2),
c    c Ediftk(iomg1,j1,j2),sumtk2,sumtk,iomg1,iomg2,j1,j2
c     print*,'Ediftk+sumtk2',cncn
c     print*,'Ejl(iomg1,j1,j2)-sum2tk+sumtk',cnc
!      stop
c       endif         
               
         endif
         

52     continue
          endif
       enddo
       enddo
       enddo
       enddo





              do iomgo=1,nmax
         nelmaxo=(nmax+1-iomgo)*(nmax+2-iomgo)/2
         do jo=1,nelmaxo
         do joo=1,nelmaxo

         dif=Ejl(iomgo,jo,joo)-E0
         dif2=ediftk(iomgo,jo,joo)-E00
         if(jo.eq.joo)then
c         ddif=emm(iomgo,jo)-E0
c         ddif2=edifmm(iomgo,jo)-E0
       if(emasef(iomgo,jo).le.energiamax.and.
     c    emasef(iomgo,jo).ge.energiamin)then  
c         print*,'ddif,ddif2',ddif,ddif2,iomgo,jo,joo
       endif
         if(Ddif.le.0.d0.or.Ddif2.le.0.d0)then
c         print*,'ddif,ddif2',ddif,ddif2,iomgo,jo,joo
         endif
         endif
         if(jo.eq.joo)goto 4454
       if(emasef(iomgo,jo).gt.enmax.or.  
     c   emasef(iomgo,joo).Gt.enmax)goto 4454
      if(emasef(iomgo,jo).gt.energiamax.and.
     c emasef(iomgo,joo).gt.energiamax)goto 4454      

      if(emasef(iomgo,jo).lt.energiamin.and.
     c emasef(iomgo,joo).lt.energiamin)goto 4454      
c        print*,'dif,dif2i,e1,e2,emax,emin',dif,dif2,iomgo,jo,joo,
c    c emasef(iomgo,jo),emasef(iomgo,joo),energiamax,energiamin
c         print*,Ejl(iomgo,iomgo,jo,joo),E00
4454   continue

      enddo
      enddo
      enddo
                             

      sumeps1=0
      sumeps2=0
      suma1pa2=0
      sumcr3=0
      sumdelt=0
      sumeps1eps2=0
      sumeps1a1pa2=0
      sumeps1cr3=0
      sumeps1delt=0
      sumeps2a1pa2=0
      sumeps2cr3=0
      sumeps2delt=0
      suma1pa2cr3=0
      suma1pa2delt=0
      sumcr3delt=0
      call pereche2(delaaa,efermi,nmax) 
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j1=1,nelmax
      do j2=1,nelmax
c      if(emasef(iomg,j1).gt.energiamax.
c     c    or.emasef(iomg,j1).lt.energiamin)goto 90
c      if(emasef(iomg,j2).gt.energiamax.
c     c    or.emasef(iomg,j2).lt.energiamin)goto 90
      if(emasef(iomg,j1).gt.enmax)goto 90
      if(emasef(iomg,j2).gt.enmax)goto 90
            if(j1.eq.j2)goto 89
!!!!!!!!!!!!!!!!  TERMENII DIAGONALI
c           if(Ejl(iomg,iomg,j1,j2).eq.0)goto 89
c      rapeneg=(Ejl(iomg,iomg,j1,j2)-E0)/
c     c  (Ejl(iomg,iomg,j1,j2)+Tk(iomg,j1)+Tk(iomg,j2)-E0)**2/
c     c  (emasef(iomg,j1)-emasef(iomg,j2))**2

c ar trebui sa incerc

      rapeneg=(Ejl(iomg,j1,j2)-E0)/
     c  (Ediftk(iomg,j1,j2)-E00)**2/
!     c  (Ejl(iomg,j1,j2)-E0-ssumtk2(iomg,j1,j2)+sumtk)**2/          !(Ediftk(iomg,j1,j2)-E00)**2/
     c  (emasef(iomg,j1)-emasef(iomg,j2))**2


c      rap1=akap(iomg,j1)*dsqrt(rho(iomg,j1))*dabs(akap(iomg,j2))/
c    c    dabs(akap(iomg,j1))/dsqrt(rho(iomg,j2))
c      rap2=akap(iomg,j2)*dsqrt(rho(iomg,j2))*dabs(akap(iomg,j1))/
c    c    dabs(akap(iomg,j2))/dsqrt(rho(iomg,j1))
c      rapeneg=(rap1-rap2)**2*rapeneg

      rrr=(u(iomg,j1)*v(iomg,j2)-u(iomg,j2)*v(iomg,j1))**2
      rapeneg=rapeneg*rrr
      sumeps1=sumeps1+elmateps1(iomg,j1,j2)*elmateps1(iomg,j2,j1)*
     c        rapeneg
      sumeps2=sumeps2+elmateps2(iomg,j1,j2)*elmateps2(iomg,j2,j1)*
     c        rapeneg
      suma1pa2=suma1pa2+elmata1pa2(iomg,j1,j2)*elmata1pa2(iomg,j2,j1)*
     c        rapeneg
      sumcr3=sumcr3+elmatcr3(iomg,j1,j2)*elmatcr3(iomg,j2,j1)*
     c        rapeneg
      sumdelt=sumdelt+elmatdelt(iomg,j1,j2)*elmatdelt(iomg,j2,j1)*
     c        rapeneg
c     sumdelt=sumdelt+elmatdelt(iomg,j1,j2)*elmatdelt(iomg,j2,j1)/
c    c  (eq(iomg,j1)+eq(iomg,j2))**3*
c    c  (u(iomg,j1)*v(iomg,j2)+u(iomg,j2)*v(iomg,j1))**2
            edif1=Ejl(iomg,j1,j2)-E0
            edif2=Ediftk(iomg,j1,j2)-E00
            edif3=emasef(iomg,j1)-emasef(iomg,j2)
            ssss=sumtk-ssumtk2(iomg,j1,j2)
c      write(23,*)sumdelt,iomg,j1,j2,energiamax,energiamin,
c    c   emasef(iomg,j1),emasef(iomg,j2),edif1,edif2,edif3
c      write(23,*)'difsumtk',ssss
!!!!!!!!!!!!!!!!!! TERMENII NEDIAGONALI
      sumeps1eps2=sumeps1eps2+elmateps1(iomg,j1,j2)*
     c  elmateps2(iomg,j2,j1)*
     c        rapeneg
      sumeps1a1pa2=sumeps1a1pa2+elmateps1(iomg,j1,j2)*
     c  elmata1pa2(iomg,j2,j1)*
     c        rapeneg
      sumeps1cr3=sumeps1cr3+elmateps1(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)*
     c        rapeneg
      sumeps1delt=sumeps1delt+elmateps1(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)*
     c        rapeneg
      sumeps2a1pa2=sumeps2a1pa2+elmateps2(iomg,j1,j2)*
     c  elmata1pa2(iomg,j2,j1)*
     c        rapeneg
      sumeps2cr3=sumeps2cr3+elmateps2(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)*
     c        rapeneg
      sumeps2delt=sumeps2delt+elmateps2(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)*
     c        rapeneg
      suma1pa2cr3=suma1pa2cr3+elmata1pa2(iomg,j1,j2)*
     c  elmatcr3(iomg,j2,j1)*
     c        rapeneg
      suma1pa2delt=suma1pa2delt+elmata1pa2(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)*
     c        rapeneg
      sumcr3delt=sumcr3delt+elmatcr3(iomg,j1,j2)*
     c  elmatdelt(iomg,j2,j1)*
     c        rapeneg
          goto 90
89        continue
cccccccc  fac o aproximatie pentru termenii diagonali
c         if(u(iomg,j1).gt.0.d0.and.v(iomg,j1).gt.0.d0)then
       goto 110
      rapeneg=(Emm(iomg,j1)-E0)/
     c  (edifmm(iomg,j1)-E00)**2



c      print*,rapeneg,Emm(iomg,j1),E0,Tk(iomg,j1),iomg,j1
      aeps1=0.5/eq(iomg,j1)**2*
     c (delaaa*(elmateps1(iomg,j1,j1)-dferdeps1)-
     c (emasef(iomg,j1)-efermi)*dgapdeps1)
      aeps2=0.5/eq(iomg,j1)**2*
     c (delaaa*(elmateps2(iomg,j1,j1)-dferdeps2)-
     c (emasef(iomg,j1)-efermi)*dgapdeps2)
      aa1pa2=0.5/eq(iomg,j1)**2*
     c (delaaa*(elmata1pa2(iomg,j1,j1)-dferda1pa2)-
     c (emasef(iomg,j1)-efermi)*dgapda1pa2)
      acr3=0.5/eq(iomg,j1)**2*
     c (delaaa*(elmatcr3(iomg,j1,j1)-dferdcr3)-
     c (emasef(iomg,j1)-efermi)*dgapdcr3)
      adelt=0.5/eq(iomg,j1)**2*
     c (delaaa*(elmatdelt(iomg,j1,j1)-dferddelt)-
     c (emasef(iomg,j1)-efermi)*dgapddelt)
c      print*,'adelt,elmatdelt(iomg,j1,j1),eq(iomg,j1),deltam(iomg,j1)'
c      print*,adelt,elmatdelt(iomg,j1,j1),eq(iomg,j1),deltam(iomg,j1)
c      print*,'j1,j2',j1,j2
c     write(91,*)adelt,elmatdelt(iomg,j1,j1),eq(iomg,j1),deltam(iomg,j1)
c    c,dferddelt,emasef(iomg,j1),dgapddelt,iomg,j1,j2
      sumeps1=sumeps1+aeps1*aeps1*
     c        rapeneg
      sumeps2=sumeps2+aeps2*aeps2*
     c        rapeneg
      suma1pa2=suma1pa2+aa1pa2*aa1pa2*
     c        rapeneg
      sumcr3=sumcr3+acr3*acr3*
     c        rapeneg
      sumdelt=sumdelt+adelt*adelt*
     c        rapeneg
c      print*,'sumdelt',sumdelt
      write(91,*),sumdelt,rapeneg,Emm(iomg,j1),E0,Tk(iomg,j1),iomg,j1
!!!!!!!!!!!!!!!!!! TERMENII NEDIAGONALI
      sumeps1eps2=sumeps1eps2+aeps1*
     c  aeps2*
     c        rapeneg
      sumeps1a1pa2=sumeps1a1pa2+aeps1*
     c  aa1pa2*
     c        rapeneg
      sumeps1cr3=sumeps1cr3+aeps1*
     c  acr3*
     c        rapeneg
      sumeps1delt=sumeps1delt+aeps1*
     c  adelt*
     c        rapeneg
      sumeps2a1pa2=sumeps2a1pa2+aeps2*
     c  aa1pa2*
     c        rapeneg
      sumeps2cr3=sumeps2cr3+aeps2*
     c  acr3*
     c        rapeneg
      sumeps2delt=sumeps2delt+aeps2*
     c  adelt*
     c        rapeneg
      suma1pa2cr3=suma1pa2cr3+aa1pa2*
     c  acr3*
     c        rapeneg
      suma1pa2delt=suma1pa2delt+aa1pa2*
     c  adelt*
     c        rapeneg
      sumcr3delt=sumcr3delt+acr3*
     c  adelt*
     c        rapeneg
c        endif
110      continue
90       continue
      enddo
      enddo
      enddo      
c          print*,delaaa,efermi
      efeps1=2*sumeps1  !+peps1
      efeps2=2*sumeps2  !+peps2
      efa1pa2=2*suma1pa2  !+pa1pa2
      efcr3=2*sumcr3  !+pcr3
      efdelt=2*sumdelt  !+pdelt
      efeps1eps2=2*sumeps1eps2  !+peps1eps2
      efeps1a1pa2=2*sumeps1a1pa2  !+peps1a1pa2
      efeps1cr3=2*sumeps1cr3  !+peps1cr3
      efeps1delt=2*sumeps1delt  !+peps1delt
      efeps2a1pa2=2*sumeps2a1pa2  !+peps2a1pa2
      efeps2cr3=2*sumeps2cr3  !+peps2cr3
      efeps2delt=2*sumeps2delt  !+peps2delt
      efa1pa2cr3=2*suma1pa2cr3 !+pa1pa2cr3
      efa1pa2delt=2*suma1pa2delt  !+pa1pa2delt
      efcr3delt=2*sumcr3delt  !+pcr3delt
           print*,'efdelt',efdelt

      return
      end

      subroutine bcssolu(a0,delaaa,efermi)
      implicit double precision (a-h,o-z)
      dimension en(2925)
      dimension ea(2925)
      common/nivselect2/eselect(2925),nnnnn


c sa determin energiile in punctul x
      do i=1,nnnnn
      ea(i)=eselect(i)
      enddo
       

c trebuie sa punem toate nivelele in ordine
 
      jh=nnnnn
      call ordoneza(ea,jh)
c am ordonat toate nivelele

              do ip=1,jh
              en(ip)=ea(ip)
              enddo

               nnn=jh
      a=a0
      call bcseq(a,en,nnn,delaaa,efermi)

      return
      end




      Subroutine bcseq(a,en,n,delta,fermi)
      implicit double precision (a-h,o-z)
      real eps
      dimension en(2925)
c en nivele dublu ocupate in ordine crescatoare
c n numar de nivele deasup[ra si sub energia Fermi dublu ocupate
      
      common/bcsequ/a0,e(2925),efm,delt,ier2,nn
      common/bcsefermi/efermi0
      common/informare/informare
      common/ginter/ginte ! constanta de interactie
      common/ginter2/ginte2 ! constanta de interactie
      external bcs1,bcsdelta0
          informare=0
      ginte2=ginte
      a0=a
      nn=n
      do i=1,n
      e(i)=en(i)
      enddo
           inf=0      
c cautam valorile lui efermi pentru a rezolva sistemul
1000  continue
           inf=inf+1
      n1=n/2
      n2=n1+1
      ef=(e(n1)+e(n2))/2.d0
      eff=ef
      do j=1,50
      ef1=ef
      ef2=ef1+0.1d0
      ef=ef2
      eps=1.e-5
c       print*,'1'
      call drtmi(x,v,bcs1,ef1,ef2,eps,100,ier)
cc       print*,'+ier,ier2',ier,ier2
      if(ier.eq.0.and.ier2.eq.0)then
      fermi=x
      delta=delt
      return
      endif
      enddo
      ef=eff
      do j=1,50
      ef1=ef
      ef2=ef1-0.1d0
      ef=ef2
      eps=1.e-5
      call drtmi(x,v,bcs1,ef2,ef1,eps,100,ier)
cc       print*,'-ier,ier2',ier,ier2
      if(ier.eq.0.and.ier2.eq.0)then
      fermi=x
      delta=delt
      return
      endif
      enddo

      if(ier.eq.2.or.ier2.eq.2)then
       ginte2=ginte2*1.05
      if(inf.eq.20)goto 3000
      goto 1000
      endif
3000  continue


         informare=1
        print*,'bcseq,ier,ier2',ier,ier2,' inf',inf
          fermi=eff
          delta=0.001d0
cc caut o energie fermi pentru un delta=0.001
c           ef1=eff-1
c           ef2=eff+1
c      call drtmi(x,v,bcsdelta0,ef1,ef2,eps,100,ier)
c          fermi=x
c      fermi=efermi0
c      if(fermi.lt.eff)fermi=eff
      return
      end   

      double precision function bcsdelta0(ef)
      implicit double precision (a-h,o-z)
      common/bcsequ/a0,e(2925),efc,delti,ier2,nn     
      delta=0.001d0
      dlt=delta
      sum=0
      do i=1,nn
      enr=e(i)-ef
      rad=dsqrt(dlt**2+enr**2)
      rap=enr/rad
      v2=0.5d0*(1.d0-rap)
      sum=sum+v2
      enddo
      bcsdelta0=2*sum-nn
      return
      end
      

      double precision function bcs1(ef1)
      
      implicit double precision (a-h,o-z)
      real eps
      common/bcsequ/a0,e(2925),ef,delti,ier2,nn
      common/v2/v2(2925)
      external bcs2

c       print*,'intru bcs1'

      eps=1.e-5
      ef=ef1
      delt=0.00001d0
      do j=1,50
      delt1=delt
      delt2=delt1+0.1
      delt=delt2
      call ddrtmi(x,v,bcs2,delt1,delt2,eps,100,ier)
      ier2=ier
      if(ier.eq.0)goto 10
      enddo
 10   continue
      sum=0
      delti=x
      do i=1,nn
      sum=sum+v2(i)
      enddo
      bcs1=2*sum-nn
cc      print*,'sum,nn,delti,ef1',sum,nn,delti,ef1
      return
      end


      double precision function bcs2(dlt)
      implicit double precision (a-h,o-z)
      common/bcsequ/a0,e(2925),ef,delt,ier2,nn
      common/v2/v2(2925)
      common/ginter2/ginte ! constanta de interactie

c       print*,'intru bcs2'


      delt=dlt
      sum=0.d0
      do j=1,nn
      enr=e(j)-ef
c      print*,'enr',enr
      rad=dsqrt(dlt**2+enr**2)
      rap=enr/rad
      sum=sum+1.d0/rad
      v2(j)=0.5d0*(1.d0-rap)
      enddo
c     g=ginte/a0
                    g=ginte
c       print*,'sum',sum
      bcs2=sum-2.d0/g
      g2=2.d0/g
cc      print*,'sum,g2,delt,ef',sum,g2,delt,ef
      return
      end


      subroutine ordoneza(ea,nr)
      implicit double precision (a-h,o-z)
      dimension ea(2925),e(2925)


      do i=1,nr
      e(i)=ea(i)
      enddo

      do i=1,nr-1
      do j=i+1,nr
      if(e(i).gt.e(j))then
      ee=e(i)
      e(i)=e(j)
      e(j)=ee
      endif
      enddo
      enddo

      do i=1,nr
      ea(i)=e(i)
      enddo
      return
      end



      subroutine dpot5(xll,xuu)
! calculez derivatele potential fata de mase efective
      implicit double precision (a-h,o-z)
      dimension ro(32),z(32),zz1(16),zz2(16)
      common/gradlib/ga0,geps1,geps2,ga1pa2,gcr3,gdelt
      common/dder11eps1/der11eps1(32,32)
      common/dder12eps1/der12eps1(32,32)
      common/dder21eps1/der21eps1(32,32)
      common/dder22eps1/der22eps1(32,32)
      common/dder11eps2/der11eps2(32,32)
      common/dder12eps2/der12eps2(32,32)
      common/dder21eps2/der21eps2(32,32)
      common/dder22eps2/der22eps2(32,32)
      common/dder11a1pa2/der11a1pa2(32,32)
      common/dder12a1pa2/der12a1pa2(32,32)
      common/dder21a1pa2/der21a1pa2(32,32)
      common/dder22a1pa2/der22a1pa2(32,32)
      common/dder11cr3/der11cr3(32,32)
      common/dder12cr3/der12cr3(32,32)
      common/dder21cr3/der21cr3(32,32)
      common/dder22cr3/der22cr3(32,32)
      common/dder11delt/der11delt(32,32)
      common/dder12delt/der12delt(32,32)
      common/dder21delt/der21delt(32,32)
      common/dder22delt/der22delt(32,32)

      common/rezmasef/dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta

      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/nrizospin/nizospin
c     ckst=41.d0

      a0=ga0
      eps1=geps1
      eps2=geps2
      a1pa2=ga1pa2
      cr3=gcr3
      delta=gdelt

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      WRO=W0*R0/b1
      alfro=DSQRT(MPH*Wro)*1.D-15
       if(nizospin.eq.0)then
      call deri_gl(a0,eps1,eps2,a1pa2,delta,cr3)
       endif


      ro(1)=.11175139809793770D3
      ro(2)=.9882954286828397D2    
      ro(3)=.8873534041789240D2
      ro(4)=.8018744697791352D2
      ro(5)=.7268762809066271D2
      ro(6)=.65975377287935053D2
      ro(7)=.59892509162134018D2
      ro(8)=.54333721333396907D2
      ro(9)=.49224394987308639D2
      ro(10)=.44509207995754938D2
      ro(11)=.40145719771539442D2
      ro(12)=.36100494805751974D2
      ro(13)=.32346629153964737D2
      ro(14)=.28862101816323475D2
      ro(15)=.25628636022459248D2
      ro(16)=.22630889013196774D2
      ro(17)=.19855860940336055D2
      ro(18)=.17292454336715315D2
      ro(19)=.14931139755522557D2
      ro(20)=.12763697986742725D2
      ro(21)=.10783018632539972D2
      ro(22)=.8982940924212596D1
      ro(23)=.7358126733186241D1
      ro(24)=.59039585041742439D1
      ro(25)=.46164567697497674D1
      ro(26)=.34922132730219945D1
      ro(27)=.25283367064257949D1
      ro(28)=.17224087764446454D1
      ro(29)=.10724487538178176D1
      ro(30)=.57688462930188643D0
      ro(31)=.23452610951961854D0
      ro(32)=.44489365833267018D-1
              
      z(1)=.10526123167960546D2
      z(2)=.9895287586829539D1
      z(3)=.9373159549646721D1
      z(4)=.8907249099964770D1
      z(5)=.8477529083379863D1
      z(6)=.8073687285010225D1
      z(7)=.7689540164040497D1
      z(8)=.7321013032780949D1
      z(9)=.69652411205511075D1
      z(10)=.66201122626360274D1
      z(11)=.62840112287748282D1
      z(12)=.59556663267994860D1
      z(13)=.56340521643499721D1
      z(14)=.53183252246332709D1
      z(15)=.50077796021987682D1
      z(16)=.47018156474074998D1
      z(17)=.43999171682281376D1
      z(18)=.41016344745666567D1
      z(19)=.38065715139453605D1
      z(20)=.35143759357409062D1
      z(21)=.32247312919920357D1
      z(22)=.29373508230046218D1
      z(23)=.26519724354306350D1
      z(24)=.23683545886324014D1
      z(25)=.20862728798817620D1
      z(26)=.18055171714655449D1
      z(27)=.15258891402098637D1
      z(28)=.12472001569431179D1
      z(29)=.9692694230711780D0
      z(30)=.69192230581004458D0
      z(31)=.41498882412107868D0
      z(32)=.13830224498700972D0

        xu=0.d0
        xl=xll
	Az1=.5D0*(XU+XL)
	Bz1=XU-XL
	zz1(1)=.49863193092474078D0*Bz1
	zz1(2)=.49280575577263417D0*Bz1
c	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	zz1(3)=.48238112779375322D0*Bz1
	zz1(4)=.46745303796886984D0*Bz1
	zz1(5)=.44816057788302606D0*Bz1
	zz1(6)=.42468380686628499D0*Bz1
	zz1(7)=.39724189798397120D0*Bz1
	zz1(8)=.36609105937014484D0*Bz1
	zz1(9)=.33152213346510760D0*Bz1
	zz1(10)=.29385787862038116D0*Bz1
	zz1(11)=.25344995446611470D0*Bz1
	zz1(12)=.21067563806531767D0*Bz1
	zz1(13)=.16593430114106382D0*Bz1
	zz1(14)=.11964368112606854D0*Bz1
	zz1(15)=.7223598079139825D-1*Bz1
	zz1(16)=.24153832843869158D-1*Bz1
	
        xu=xuu
        xl=0.d0
	Az2=.5D0*(XU+XL)
	Bz2=XU-XL
	zz2(1)=.49863193092474078D0*Bz2
	zz2(2)=.49280575577263417D0*Bz2
c	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	zz2(3)=.48238112779375322D0*Bz2
	zz2(4)=.46745303796886984D0*Bz2
	zz2(5)=.44816057788302606D0*Bz2
	zz2(6)=.42468380686628499D0*Bz2
	zz2(7)=.39724189798397120D0*Bz2
	zz2(8)=.36609105937014484D0*Bz2
	zz2(9)=.33152213346510760D0*Bz2
	zz2(10)=.29385787862038116D0*Bz2
	zz2(11)=.25344995446611470D0*Bz2
	zz2(12)=.21067563806531767D0*Bz2
	zz2(13)=.16593430114106382D0*Bz2
	zz2(14)=.11964368112606854D0*Bz2
	zz2(15)=.7223598079139825D-1*Bz2
	zz2(16)=.24153832843869158D-1*Bz2
	

      do i=1,32
      rhoo=dsqrt(ro(i))/alfro
      do j=1,32
      z1p=(-z(j)-zet1)/alf1
      z2p=(z(j)+zet2)/alf2


      call dervate(z1p,rhoo)
      der11eps1(j,i)=dvdeps1
      der11eps2(j,i)=dvdeps2
      der11a1pa2(j,i)=dvda1pa2
      der11cr3(j,i)=dvdcr3
      der11delt(j,i)=dvddelta
c        print*,'1'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta
      call dervate(z2p,rhoo)
      der12eps1(j,i)=dvdeps1
      der12eps2(j,i)=dvdeps2
      der12a1pa2(j,i)=dvda1pa2
      der12cr3(j,i)=dvdcr3
      der12delt(j,i)=dvddelta

c        print*,'2'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta

      enddo
      enddo

      do i=1,32
      rhoo=dsqrt(ro(i))/alfro
      do j=1,16
      zz1p=(az1+zz1(j))/alf1
        jj=2*(j-1)+1
      call dervate(zz1p,rhoo)
      der21eps1(jj,i)=dvdeps1
      der21eps2(jj,i)=dvdeps2
      der21a1pa2(jj,i)=dvda1pa2
      der21cr3(jj,i)=dvdcr3
      der21delt(jj,i)=dvddelta

 
c        print*,'3'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta
      zz1p=(az1-zz1(j))/alf1
        jj=jj+1
      call dervate(zz1p,rhoo)
      der21eps1(jj,i)=dvdeps1
      der21eps2(jj,i)=dvdeps2
      der21a1pa2(jj,i)=dvda1pa2
      der21cr3(jj,i)=dvdcr3
      der21delt(jj,i)=dvddelta


c        print*,'4'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta
      zz2p=(az2+zz2(j))/alf2
        jj=2*(j-1)+1
      call dervate(zz2p,rhoo)
      der22eps1(jj,i)=dvdeps1
      der22eps2(jj,i)=dvdeps2
      der22a1pa2(jj,i)=dvda1pa2
      der22cr3(jj,i)=dvdcr3
      der22delt(jj,i)=dvddelta


c        print*,'5'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta
      zz2p=(az2-zz2(j))/alf2
        jj=jj+1
      call dervate(zz2p,rhoo)
      der22eps1(jj,i)=dvdeps1
      der22eps2(jj,i)=dvdeps2
      der22a1pa2(jj,i)=dvda1pa2
      der22cr3(jj,i)=dvdcr3
      der22delt(jj,i)=dvddelta


c        print*,'6'
c      print*,dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta

      enddo
      enddo




      return
      end

      subroutine dervate(zi,rho)
      !! calculez dervivatele potentialului
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,s,delt
      common/prtdta/z0,ro0 
      common/dltw/deltax
      common/dereps1/a1eps1,a2eps1,b1eps1,b2eps1,c1eps1,c2eps1,c3eps1,
     c               ro3eps1,x101eps1,x201eps1,y101eps1,y201eps1
      common/dereps2/a1eps2,a2eps2,b1eps2,b2eps2,c1eps2,c2eps2,c3eps2,
     c               ro3eps2,x101eps2,x201eps2,y101eps2,y201eps2
      common/dera1pa2/a1a1pa2,a2a1pa2,b1a1pa2,b2a1pa2,c1a1pa2,c2a1pa2,
     c       c3a1pa2,ro3a1pa2,x101a1pa2,x201a1pa2,y101a1pa2,y201a1pa2
      common/dercr3/a1cr3,a2cr3,b1cr3,b2cr3,c1cr3,c2cr3,c3cr3,ro3cr3,
     c               x101cr3,x201cr3,y101cr3,y201cr3
      common/derdelta/a1delta,a2delta,b1delta,b2delta,c1delta,c2delta,
     c       c3delta,ro3delta,x101delta,x201delta,y101delta,y201delta
      common/adancimi/vp1,vp2
      common/oadancime/v000
      common/adifz/adiz
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
!
      common/rezmasef/dvdeps1,dvdeps2,dvda1pa2,dvdcr3,dvddelta
      external delta1,delta2

        
        ro0=dabs(rho)
        ro=ro0
        z=zi-z0p
        z0=z
      V0=vp1
      ca=adiz
      


ccc      call deri_gl(a0,eps1,eps2,a1pa2,delta,cr3)



c deosebim mai multe cazuri
c      ro=dabs(rho)
c      z0=z
c      ro0=ro
      a1pa2pr3=a1+a2+2*r3-1.d-6
      if(a1pa2pr3.lt.delt)then ! cazul elipselor complet separate
! in acest caz utilizam planul de separare z=0
         if(z.lt.0.d0)then
! elipsa din partea stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
                  if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1

                        if(z.gt.c1)then
                        delta11=z-(c1+a1)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1delta-a1delta)
         return
                        else
                        delta11=c1-a1-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1delta-a1delta)
         return
                        endif

                  else      
                  if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1  
        delta11=ro0-b1
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1delta)
        return
                  else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
      x00=x0
      if(z.gt.c1)then
      diff=c1-x0
      x00=c1+diff
      endif
                              !!!!!!!!! eps1 !!!!!!!!
       dxdeps1=-2*b1/a1*(b1eps1/a1-b1*a1eps1/a1**2)+c1eps1/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps1=dxdeps1+(b1eps1/a1**2-2.d0*b1/a1**3*a1eps1)/rad
       dxdeps1=dxdeps1+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps1/a1-(x00-c1)/a1**2*a1eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps1=b1eps1*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps1-c1eps1)/a1-(x00-c1)/a1**2*a1eps1)

                              !!!!!!!!! eps2 !!!!!!!!
       dxdeps2=-2*b1/a1*(b1eps2/a1-b1*a1eps2/a1**2)+c1eps2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps2=dxdeps2+(b1eps2/a1**2-2.d0*b1/a1**3*a1eps2)/rad
       dxdeps2=dxdeps2+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps2/a1-(x00-c1)/a1**2*a1eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps2=b1eps2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps2-c1eps2)/a1-(x00-c1)/a1**2*a1eps2)
                              !!!!!!!!! a1pa2 !!!!!!!!
       dxda1pa2=-2*b1/a1*(b1a1pa2/a1-b1*a1a1pa2/a1**2)+c1a1pa2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxda1pa2=dxda1pa2+(b1a1pa2/a1**2-2.d0*b1/a1**3*a1a1pa2)/rad
       dxda1pa2=dxda1pa2+b1/a1**2*(x00-c1)/a1*
     !         (-c1a1pa2/a1-(x00-c1)/a1**2*a1a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyda1pa2=b1a1pa2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxda1pa2-c1a1pa2)/a1-(x00-c1)/a1**2*a1a1pa2)
                              !!!!!!!!! cr3 !!!!!!!!
       dxdcr3=-2*b1/a1*(b1cr3/a1-b1*a1cr3/a1**2)+c1cr3/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdcr3=dxdcr3+(b1cr3/a1**2-2.d0*b1/a1**3*a1cr3)/rad
       dxdcr3=dxdcr3+b1/a1**2*(x00-c1)/a1*
     !         (-c1cr3/a1-(x00-c1)/a1**2*a1cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydcr3=b1cr3*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdcr3-c1cr3)/a1-(x00-c1)/a1**2*a1cr3)
                              !!!!!!!!! delta !!!!!!!!
       dxddelta=-2*b1/a1*(b1delta/a1-b1*a1delta/a1**2)+c1delta/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxddelta=dxddelta+(b1delta/a1**2-2.d0*b1/a1**3*a1delta)/rad
       dxddelta=dxddelta+b1/a1**2*(x00-c1)/a1*
     !         (-c1delta/a1-(x00-c1)/a1**2*a1delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyddelta=b1delta*rad-b1/rad*(x00-c1)/a1*
     !         ((dxddelta-c1delta)/a1-(x00-c1)/a1**2*a1delta)
        endif
        endif
! pana aici elipsa stanga
          else
! elipsa din partea dreapta (2)
        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1


                        if(z.gt.c2)then
                        delta11=z-(c2+a2)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2delta-a2delta)
         return
                        else
                        delta11=c2-a2-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2delta-a2delta)
         return
                       endif
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
  
        delta11=ro0-b2
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2delta)
        return
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
      x00=x0
      if(z.gt.c2)then
      diff=c2-x0
      x00=c2+diff
      endif
                !!!!!!!!!!!!!! eps1 !!!!!!!!!!!!!!!!
       dxdeps1=-2*b2/a2*(b2eps1/a2-b2*a2eps1/a2**2)+c2eps1/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps1=dxdeps1+(b2eps1/a2**2-2.d0*b2/a2**3*a2eps1)/rad
       dxdeps1=dxdeps1+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps1/a2-(x00-c2)/a2**2*a2eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps1=b2eps1*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps1-c2eps1)/a2-(x00-c2)/a2**2*a2eps1)
    

                !!!!!!!!!!!!!! eps2 !!!!!!!!!!!!!!!!
       dxdeps2=-2*b2/a2*(b2eps2/a2-b2*a2eps2/a2**2)+c2eps2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps2=dxdeps2+(b2eps2/a2**2-2.d0*b2/a2**3*a2eps2)/rad
       dxdeps2=dxdeps2+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps2/a2-(x00-c2)/a2**2*a2eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps2=b2eps2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps2-c2eps2)/a2-(x00-c2)/a2**2*a2eps2)

                !!!!!!!!!!!!!! a1pa2 !!!!!!!!!!!!!!!!
       dxda1pa2=-2*b2/a2*(b2a1pa2/a2-b2*a2a1pa2/a2**2)+c2a1pa2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxda1pa2=dxda1pa2+(b2a1pa2/a2**2-2.d0*b2/a2**3*a2a1pa2)/rad
       dxda1pa2=dxda1pa2+b2/a2**2*(x00-c2)/a2*
     !         (-c2a1pa2/a2-(x00-c2)/a2**2*a2a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyda1pa2=b2a1pa2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxda1pa2-c2a1pa2)/a2-(x00-c2)/a2**2*a2a1pa2)

                !!!!!!!!!!!!!! cr3 !!!!!!!!!!!!!!!!
       dxdcr3=-2*b2/a2*(b2cr3/a2-b2*a2cr3/a2**2)+c2cr3/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdcr3=dxdcr3+(b2cr3/a2**2-2.d0*b2/a2**3*a2cr3)/rad
       dxdcr3=dxdcr3+b2/a2**2*(x00-c2)/a2*
     !         (-c2cr3/a2-(x00-c2)/a2**2*a2cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydcr3=b2cr3*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdcr3-c2cr3)/a2-(x00-c2)/a2**2*a2cr3)

                !!!!!!!!!!!!!! delta !!!!!!!!!!!!!!!!
       dxddelta=-2*b2/a2*(b2delta/a2-b2*a2delta/a2**2)+c2delta/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxddelta=dxddelta+(b2delta/a2**2-2.d0*b2/a2**3*a2delta)/rad
       dxddelta=dxddelta+b2/a2**2*(x00-c2)/a2*
     !         (-c2delta/a2-(x00-c2)/a2**2*a2delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyddelta=b2delta*rad-b2/rad*(x00-c2)/a2*
     !         ((dxddelta-c2delta)/a2-(x00-c2)/a2**2*a2delta)
          endif
          endif
! pana aici elipsa dreapta
 
          endif
        delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps1)+2.*(ro0-y0)*(-dydeps1))

      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps2)+2.*(ro0-y0)*(-dydeps2))

      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxda1pa2)+2.*(ro0-y0)*(-dyda1pa2))

      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdcr3)+2.*(ro0-y0)*(-dydcr3))

      dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxddelta)+2.*(ro0-y0)*(-dyddelta))
      return



! Pana aici elipse complet separate
      endif

      if(r3.ge.400.9999d0)then ! cazul r3 infinit
! in acest caz trebuie calculate liniile de separare intre cele trei zone
!      alpha=atan((y101-201)/(x201-x101))
!      alpha=pi/2-gamma ! gamma da panta planelor de separate tan(gamma)=cotg(alpha)
! tan(gamma)=(x201-x101)/(y101-y201)
      !calculez primul segment de separare 
        if(y101.eq.y201)then
        if(z.le.x101)goto 4441 ! atunci am prima elipsa
        if(z.ge.x201)goto 4443 ! atunci am a doua elipsa
        goto 4442 ! cazul intermediar
        endif
       tang=(x201-x101)/(y101-y201)
      x0prim=x101-1.d0/tang*y101 ! calculat originea
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
          if(rodif.le.ro)goto 4441 ! atunci am prima elipsa 
      x0prim=x201-1.d0/tang*y201 ! calculat originea pentru a doua separare
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
           if(rodif.ge.ro)goto 4443 ! atunci am a doua elipsa

           goto 4442

4441      continue
! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1


                        if(z.gt.c1)then
                        delta11=z-(c1+a1)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1delta-a1delta)
         return
                        else
                        delta11=c1-a1-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1delta-a1delta)
         return
                        endif
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
  
        delta11=ro0-b1
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1delta)
        return
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)

      x00=x0
      if(z.gt.c1)then
      diff=c1-x0
      x00=c1+diff
      endif
                              !!!!!!!!! eps1 !!!!!!!!
       dxdeps1=-2*b1/a1*(b1eps1/a1-b1*a1eps1/a1**2)+c1eps1/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps1=dxdeps1+(b1eps1/a1**2-2.d0*b1/a1**3*a1eps1)/rad
       dxdeps1=dxdeps1+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps1/a1-(x00-c1)/a1**2*a1eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps1=b1eps1*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps1-c1eps1)/a1-(x00-c1)/a1**2*a1eps1)

                              !!!!!!!!! eps2 !!!!!!!!
       dxdeps2=-2*b1/a1*(b1eps2/a1-b1*a1eps2/a1**2)+c1eps2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps2=dxdeps2+(b1eps2/a1**2-2.d0*b1/a1**3*a1eps2)/rad
       dxdeps2=dxdeps2+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps2/a1-(x00-c1)/a1**2*a1eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps2=b1eps2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps2-c1eps2)/a1-(x00-c1)/a1**2*a1eps2)
                              !!!!!!!!! a1pa2 !!!!!!!!
       dxda1pa2=-2*b1/a1*(b1a1pa2/a1-b1*a1a1pa2/a1**2)+c1a1pa2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxda1pa2=dxda1pa2+(b1a1pa2/a1**2-2.d0*b1/a1**3*a1a1pa2)/rad
       dxda1pa2=dxda1pa2+b1/a1**2*(x00-c1)/a1*
     !         (-c1a1pa2/a1-(x00-c1)/a1**2*a1a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyda1pa2=b1a1pa2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxda1pa2-c1a1pa2)/a1-(x00-c1)/a1**2*a1a1pa2)
                              !!!!!!!!! cr3 !!!!!!!!
       dxdcr3=-2*b1/a1*(b1cr3/a1-b1*a1cr3/a1**2)+c1cr3/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdcr3=dxdcr3+(b1cr3/a1**2-2.d0*b1/a1**3*a1cr3)/rad
       dxdcr3=dxdcr3+b1/a1**2*(x00-c1)/a1*
     !         (-c1cr3/a1-(x00-c1)/a1**2*a1cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydcr3=b1cr3*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdcr3-c1cr3)/a1-(x00-c1)/a1**2*a1cr3)
                              !!!!!!!!! delta !!!!!!!!
       dxddelta=-2*b1/a1*(b1delta/a1-b1*a1delta/a1**2)+c1delta/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxddelta=dxddelta+(b1delta/a1**2-2.d0*b1/a1**3*a1delta)/rad
       dxddelta=dxddelta+b1/a1**2*(x00-c1)/a1*
     !         (-c1delta/a1-(x00-c1)/a1**2*a1delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyddelta=b1delta*rad-b1/rad*(x00-c1)/a1*
     !         ((dxddelta-c1delta)/a1-(x00-c1)/a1**2*a1delta)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         
         ws=-V0/(1+dexp(delta/ca))

      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps1)+2.*(ro0-y0)*(-dydeps1))

      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps2)+2.*(ro0-y0)*(-dydeps2))

      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxda1pa2)+2.*(ro0-y0)*(-dyda1pa2))

      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdcr3)+2.*(ro0-y0)*(-dydcr3))

       dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxddelta)+2.*(ro0-y0)*(-dyddelta))               
! pana aici elipsa din stanga 1

           goto 4445
4443       continue
! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1

                        if(z.gt.c2)then
                        delta11=z-(c2+a2)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2delta-a2delta)
         return
                        else
                        delta11=c2-a2-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2delta-a2delta)
         return
                       endif
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1

  
        delta11=ro0-b2
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2delta)
        return
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)

      x00=x0
      if(z.gt.c2)then
      diff=c2-x0
      x00=c2+diff
      endif
                !!!!!!!!!!!!!! eps1 !!!!!!!!!!!!!!!!
       dxdeps1=-2*b2/a2*(b2eps1/a2-b2*a2eps1/a2**2)+c2eps1/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps1=dxdeps1+(b2eps1/a2**2-2.d0*b2/a2**3*a2eps1)/rad
       dxdeps1=dxdeps1+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps1/a2-(x00-c2)/a2**2*a2eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps1=b2eps1*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps1-c2eps1)/a2-(x00-c2)/a2**2*a2eps1)
    

                !!!!!!!!!!!!!! eps2 !!!!!!!!!!!!!!!!
       dxdeps2=-2*b2/a2*(b2eps2/a2-b2*a2eps2/a2**2)+c2eps2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps2=dxdeps2+(b2eps2/a2**2-2.d0*b2/a2**3*a2eps2)/rad
       dxdeps2=dxdeps2+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps2/a2-(x00-c2)/a2**2*a2eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps2=b2eps2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps2-c2eps2)/a2-(x00-c2)/a2**2*a2eps2)

                !!!!!!!!!!!!!! a1pa2 !!!!!!!!!!!!!!!!
       dxda1pa2=-2*b2/a2*(b2a1pa2/a2-b2*a2a1pa2/a2**2)+c2a1pa2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxda1pa2=dxda1pa2+(b2a1pa2/a2**2-2.d0*b2/a2**3*a2a1pa2)/rad
       dxda1pa2=dxda1pa2+b2/a2**2*(x00-c2)/a2*
     !         (-c2a1pa2/a2-(x00-c2)/a2**2*a2a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyda1pa2=b2a1pa2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxda1pa2-c2a1pa2)/a2-(x00-c2)/a2**2*a2a1pa2)

                !!!!!!!!!!!!!! cr3 !!!!!!!!!!!!!!!!
       dxdcr3=-2*b2/a2*(b2cr3/a2-b2*a2cr3/a2**2)+c2cr3/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdcr3=dxdcr3+(b2cr3/a2**2-2.d0*b2/a2**3*a2cr3)/rad
       dxdcr3=dxdcr3+b2/a2**2*(x00-c2)/a2*
     !         (-c2cr3/a2-(x00-c2)/a2**2*a2cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydcr3=b2cr3*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdcr3-c2cr3)/a2-(x00-c2)/a2**2*a2cr3)

                !!!!!!!!!!!!!! delta !!!!!!!!!!!!!!!!
       dxddelta=-2*b2/a2*(b2delta/a2-b2*a2delta/a2**2)+c2delta/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxddelta=dxddelta+(b2delta/a2**2-2.d0*b2/a2**3*a2delta)/rad
       dxddelta=dxddelta+b2/a2**2*(x00-c2)/a2*
     !         (-c2delta/a2-(x00-c2)/a2**2*a2delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyddelta=b2delta*rad-b2/rad*(x00-c2)/a2*
     !         ((dxddelta-c2delta)/a2-(x00-c2)/a2**2*a2delta)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
               
         ws=-V0/(1+dexp(delta/ca))

      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps1)+2.*(ro0-y0)*(-dydeps1))

      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps2)+2.*(ro0-y0)*(-dydeps2))

      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxda1pa2)+2.*(ro0-y0)*(-dyda1pa2))

      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdcr3)+2.*(ro0-y0)*(-dydcr3))

      dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxddelta)+2.*(ro0-y0)*(-dyddelta))
! pana aici elipsa din dreapta

           goto 4445
4442  continue ! a ramas partea intermediara
      z0=z
      ro0=ro
      if(y101.eq.y201)then
      delta=ro0-y101

         v000=V0
         ws=-V0/(1+dexp(delta/ca))
                
       !!! R3 este infinit cu simetrie la reflectie

         dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-y101eps1)
c       print*,'dvdeps1,y101eps1',dvdeps1,y101eps1
         dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-y101eps2)
         dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-y101a1pa2)
         dvdcr3=0.d0
         dvdelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-y101delta)
      else
      alpha=datan((y101-y201)/(x201-x101))
      roc=y101+(z-x101)*(y201-y101)/(x201-x101)!calculez punctul 
                                                !corespunzator pe suprafata 
      delta=(ro0-roc)*dcos(alpha)

         v000=V0
         ws=-V0/(1+dexp(delta/ca))
       rap=(y201-y101)/(x201-x101)
       radrap=dsqrt(1.d0+rap**2)
      coseps1=rap/radrap**3*((y201eps1-y101eps1)/(x201-x101)-
     !                       rap/(x201-x101)*(x201eps1-x101eps1))
      roceps1=y101eps1-x101eps1*rap+(z-x101)*((y201eps1-y101eps1)
     !         /(x201-x101)-rap/(x201-x101)*(x201eps1-x101eps1))
      coseps2=rap/radrap**3*((y201eps2-y101eps2)/(x201-x101)-
     !                       rap/(x201-x101)*(x201eps2-x101eps2))
      roceps2=y101eps2-x101eps2*rap+(z-x101)*((y201eps2-y101eps2)
     !         /(x201-x101)-rap/(x201-x101)*(x201eps2-x101eps2))
      cosa1pa2=rap/radrap**3*((y201a1pa2-y101a1pa2)/(x201-x101)-
     !                       rap/(x201-x101)*(x201a1pa2-x101a1pa2))
      roca1pa2=y101a1pa2-x101a1pa2*rap+(z-x101)*((y201a1pa2-y101a1pa2)
     !         /(x201-x101)-rap/(x201-x101)*(x201a1pa2-x101a1pa2))
      cosdelta=rap/radrap**3*((y201delta-y101delta)/(x201-x101)-
     !                       rap/(x201-x101)*(x201delta-x101delta))
      rocdelta=y101delta-x101delta*rap+(z-x101)*((y201delta-y101delta)
     !         /(x201-x101)-rap/(x201-x101)*(x201delta-x101delta))

         dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-roceps1*dcos(alpha)+(ro0-roc)*coseps1)
         dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-roceps2*dcos(alpha)+(ro0-roc)*coseps2)
         dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-roca1pa2*dcos(alpha)+(ro0-roc)*cosa1pa2)
           dvdcr3=0.d0
         dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !           (-rocdelta*dcos(alpha)+(ro0-roc)*cosdelta)


      endif
4445  continue
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
c              write(33,*)'crau'
      return
! pana aici cazul r3 infinit
      endif


      if(s.gt.0.d0)then
      if(dabs(ro).gt.ro3)then
      xx02=c3
      goto 6688
      endif
! calculez y11 si y22 punctele de pe suprafata elipsei la trecerea
! dintre regiuni
! folosesc rapoarte pentru determinarea valorilor x functie de ro
! pe suprafetele de separare
! daca z se afla intre cele doua limite atunci am cazul regiunii mediane
                       ! s pozitiv
       y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
       xx01=(x1-c3)*(ro3-dabs(ro))/(ro3-y11)+c3

       y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)
       xx02=(x2-c3)*(ro3-dabs(ro))/(ro3-y22)+c3
           if(r3.lt.1.d-3)then
           xx01=c3
           xx02=c3
           endif
     
         if(z.ge.xx01.and.z.le.xx02)goto 6677 ! regiunea mediana

             else
                        ! s negativ
      y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
      y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)

      xx01=(x1-c3)*(dabs(ro)+s*ro3)/(s*ro3+y11)+c3
      xx02=(x2-c3)*(dabs(ro)+s*ro3)/(s*ro3+y22)+c3
         if(z.le.xx02.and.z.ge.xx01)goto 6677 ! regiunea mediana
             endif
     
      goto 6688 ! mergi intr-una din elipse


6677        continue ! regiunea mediana
      V0=(vp1+vp2)/2.d0
      ! ramane cazul intermediar pentru r3
      ! ro3 are semn negativ pentru s=-1 si pozitiv pentru s=1
                delta=-s*(dsqrt((c3-z)**2+(dabs(ro)-ro3)**2)-r3)
c                if(s.gt.0.d0)delta=-delta
      delta11=-s*dsqrt((c3-z)**2+(dabs(ro)-ro3)**2)

      dvdeps1=-((z-c3)*c3eps1+(ro0-ro3)*ro3eps1)/delta11
      dvdeps2=-((z-c3)*c3eps2+(ro0-ro3)*ro3eps2)/delta11
      dvda1pa2=-((z-c3)*c3a1pa2+(ro0-ro3)*ro3a1pa2)/delta11
      dvdcr3=-((z-c3)*c3cr3+(ro0-ro3)*ro3cr3)/delta11-r3**2
      dvddelta=-((z-c3)*c3delta+(ro0-ro3)*ro3delta)/delta11

c        print*,'s,delta,r3,dvdcr3',s,delta,r3,dvdcr3
         
         ws=-V0/(1+dexp(delta/ca))

      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*dvdeps1
      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*dvdeps2
      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*dvda1pa2
      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*dvdcr3
      dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*dvddelta
   
c        print*,'dvdcr3',dvdcr3   

                goto 5555
6688  continue
c      if(z.gt.z3.and.dabs(ro).gt.ro3)goto 5553
      if(z.ge.xx02)goto 5553
!5551  continue      

! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1


                        if(z.gt.c1)then
                        delta11=z-(c1+a1)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c1delta-a1delta)
         return
                        else
                        delta11=c1-a1-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps1-a1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1eps2-a1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1a1pa2-a1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1cr3-a1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c1delta-a1delta)
         return
                        endif
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
  
        delta11=ro0-b1
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b1delta)
        return
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)


      x00=x0
      if(z.gt.c1)then
      diff=c1-x0
      x00=c1+diff
      endif
                              !!!!!!!!! eps1 !!!!!!!!
       dxdeps1=-2*b1/a1*(b1eps1/a1-b1*a1eps1/a1**2)+c1eps1/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps1=dxdeps1+(b1eps1/a1**2-2.d0*b1/a1**3*a1eps1)/rad
       dxdeps1=dxdeps1+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps1/a1-(x00-c1)/a1**2*a1eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps1=b1eps1*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps1-c1eps1)/a1-(x00-c1)/a1**2*a1eps1)

                              !!!!!!!!! eps2 !!!!!!!!
       dxdeps2=-2*b1/a1*(b1eps2/a1-b1*a1eps2/a1**2)+c1eps2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdeps2=dxdeps2+(b1eps2/a1**2-2.d0*b1/a1**3*a1eps2)/rad
       dxdeps2=dxdeps2+b1/a1**2*(x00-c1)/a1*
     !         (-c1eps2/a1-(x00-c1)/a1**2*a1eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydeps2=b1eps2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdeps2-c1eps2)/a1-(x00-c1)/a1**2*a1eps2)
                              !!!!!!!!! a1pa2 !!!!!!!!
       dxda1pa2=-2*b1/a1*(b1a1pa2/a1-b1*a1a1pa2/a1**2)+c1a1pa2/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxda1pa2=dxda1pa2+(b1a1pa2/a1**2-2.d0*b1/a1**3*a1a1pa2)/rad
       dxda1pa2=dxda1pa2+b1/a1**2*(x00-c1)/a1*
     !         (-c1a1pa2/a1-(x00-c1)/a1**2*a1a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyda1pa2=b1a1pa2*rad-b1/rad*(x00-c1)/a1*
     !         ((dxda1pa2-c1a1pa2)/a1-(x00-c1)/a1**2*a1a1pa2)
                              !!!!!!!!! cr3 !!!!!!!!
       dxdcr3=-2*b1/a1*(b1cr3/a1-b1*a1cr3/a1**2)+c1cr3/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxdcr3=dxdcr3+(b1cr3/a1**2-2.d0*b1/a1**3*a1cr3)/rad
       dxdcr3=dxdcr3+b1/a1**2*(x00-c1)/a1*
     !         (-c1cr3/a1-(x00-c1)/a1**2*a1cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dydcr3=b1cr3*rad-b1/rad*(x00-c1)/a1*
     !         ((dxdcr3-c1cr3)/a1-(x00-c1)/a1**2*a1cr3)
                              !!!!!!!!! delta !!!!!!!!
       dxddelta=-2*b1/a1*(b1delta/a1-b1*a1delta/a1**2)+c1delta/(x00-c1)-
     !         (z-c1)/(x00-c1)**2*c1delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c1)/a1)**2)
       dxddelta=dxddelta+(b1delta/a1**2-2.d0*b1/a1**3*a1delta)/rad
       dxddelta=dxddelta+b1/a1**2*(x00-c1)/a1*
     !         (-c1delta/a1-(x00-c1)/a1**2*a1delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c1)/(x00-c1)**2-
     !         b1/a1**2*(x00-c1)/a1**2/rad**3)
       dyddelta=b1delta*rad-b1/rad*(x00-c1)/a1*
     !         ((dxddelta-c1delta)/a1-(x00-c1)/a1**2*a1delta)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         
         ws=-V0/(1+dexp(delta/ca))

      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps1)+2.*(ro0-y0)*(-dydeps1))

      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps2)+2.*(ro0-y0)*(-dydeps2))

      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxda1pa2)+2.*(ro0-y0)*(-dyda1pa2))

      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdcr3)+2.*(ro0-y0)*(-dydcr3))

      dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxddelta)+2.*(ro0-y0)*(-dyddelta))
! pana aici elipsa din stanga 1

      goto 5555

5553  continue      

! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1

                        if(z.gt.c2)then
                        delta11=z-(c2+a2)
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-c2delta-a2delta)
         return
                        else
                        delta11=c2-a2-z
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps1-a2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2eps2-a2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2a1pa2-a2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2cr3-a2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (c2delta-a2delta)
         return
                       endif
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0.and.z0.gt.x0)semn=1

  
        delta11=ro0-b2
         v000=V0
         ws=-V0/(1+dexp(delta11/ca))
         dvdeps1=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps1)
         dvdeps2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2eps2)
         dvda1pa2=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2a1pa2)
         dvdcr3=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2cr3)
         dvddelta=-ws/(1+dexp(delta11/ca))/ca*dexp(delta11/ca)*
     !           (-b2delta)
        return
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)


      x00=x0
      if(z.gt.c2)then
      diff=c2-x0
      x00=c2+diff
      endif
                !!!!!!!!!!!!!! eps1 !!!!!!!!!!!!!!!!
       dxdeps1=-2*b2/a2*(b2eps1/a2-b2*a2eps1/a2**2)+c2eps1/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps1
       dxdeps1=1/ro0*dxdeps1
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps1=dxdeps1+(b2eps1/a2**2-2.d0*b2/a2**3*a2eps1)/rad
       dxdeps1=dxdeps1+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps1/a2-(x00-c2)/a2**2*a2eps1)/rad**3
       dxdeps1=dxdeps1/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps1=b2eps1*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps1-c2eps1)/a2-(x00-c2)/a2**2*a2eps1)
    

                !!!!!!!!!!!!!! eps2 !!!!!!!!!!!!!!!!
       dxdeps2=-2*b2/a2*(b2eps2/a2-b2*a2eps2/a2**2)+c2eps2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2eps2
       dxdeps2=1/ro0*dxdeps2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdeps2=dxdeps2+(b2eps2/a2**2-2.d0*b2/a2**3*a2eps2)/rad
       dxdeps2=dxdeps2+b2/a2**2*(x00-c2)/a2*
     !         (-c2eps2/a2-(x00-c2)/a2**2*a2eps2)/rad**3
       dxdeps2=dxdeps2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydeps2=b2eps2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdeps2-c2eps2)/a2-(x00-c2)/a2**2*a2eps2)

                !!!!!!!!!!!!!! a1pa2 !!!!!!!!!!!!!!!!
       dxda1pa2=-2*b2/a2*(b2a1pa2/a2-b2*a2a1pa2/a2**2)+c2a1pa2/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2a1pa2
       dxda1pa2=1/ro0*dxda1pa2
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxda1pa2=dxda1pa2+(b2a1pa2/a2**2-2.d0*b2/a2**3*a2a1pa2)/rad
       dxda1pa2=dxda1pa2+b2/a2**2*(x00-c2)/a2*
     !         (-c2a1pa2/a2-(x00-c2)/a2**2*a2a1pa2)/rad**3
       dxda1pa2=dxda1pa2/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyda1pa2=b2a1pa2*rad-b2/rad*(x00-c2)/a2*
     !         ((dxda1pa2-c2a1pa2)/a2-(x00-c2)/a2**2*a2a1pa2)

                !!!!!!!!!!!!!! cr3 !!!!!!!!!!!!!!!!
       dxdcr3=-2*b2/a2*(b2cr3/a2-b2*a2cr3/a2**2)+c2cr3/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2cr3
       dxdcr3=1/ro0*dxdcr3
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxdcr3=dxdcr3+(b2cr3/a2**2-2.d0*b2/a2**3*a2cr3)/rad
       dxdcr3=dxdcr3+b2/a2**2*(x00-c2)/a2*
     !         (-c2cr3/a2-(x00-c2)/a2**2*a2cr3)/rad**3
       dxdcr3=dxdcr3/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dydcr3=b2cr3*rad-b2/rad*(x00-c2)/a2*
     !         ((dxdcr3-c2cr3)/a2-(x00-c2)/a2**2*a2cr3)

                !!!!!!!!!!!!!! delta !!!!!!!!!!!!!!!!
       dxddelta=-2*b2/a2*(b2delta/a2-b2*a2delta/a2**2)+c2delta/(x00-c2)-
     !         (z-c2)/(x00-c2)**2*c2delta
       dxddelta=1/ro0*dxddelta
       rad=dsqrt(1.d0-((x0-c2)/a2)**2)
       dxddelta=dxddelta+(b2delta/a2**2-2.d0*b2/a2**3*a2delta)/rad
       dxddelta=dxddelta+b2/a2**2*(x00-c2)/a2*
     !         (-c2delta/a2-(x00-c2)/a2**2*a2delta)/rad**3
       dxddelta=dxddelta/(-1/ro0*(z-c2)/(x00-c2)**2-
     !         b2/a2**2*(x00-c2)/a2**2/rad**3)
       dyddelta=b2delta*rad-b2/rad*(x00-c2)/a2*
     !         ((dxddelta-c2delta)/a2-(x00-c2)/a2**2*a2delta)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
               
         ws=-V0/(1+dexp(delta/ca))

      dvdeps1=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps1)+2.*(ro0-y0)*(-dydeps1))

      dvdeps2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdeps2)+2.*(ro0-y0)*(-dydeps2))

      dvda1pa2=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxda1pa2)+2.*(ro0-y0)*(-dyda1pa2))

      dvdcr3=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxdcr3)+2.*(ro0-y0)*(-dydcr3))

      dvddelta=-ws/(1+dexp(delta/ca))/ca*dexp(delta/ca)*
     !   1./2./delta*(2.*(z-x00)*(-dxddelta)+2.*(ro0-y0)*(-dyddelta))
      goto 5555
5555  continue
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
      return
      end
            



      subroutine deri_gl(a0,eps1,eps2,a1pa2,delta,cr3)
      !! calculez dervivatele parametrilor a1,a2,b1,b2,c1,c2,c3,ro3,r3,c3 
      !! fata de grade
      !! de libertate si le trimit prin common mai departe
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1x,b1x,a2x,b2x,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,deltax
      common/dereps1/a1eps1,a2eps1,b1eps1,b2eps1,c1eps1,c2eps1,c3eps1,
     c               ro3eps1,x101eps1,x201eps1,y101eps1,y201eps1
      common/dereps2/a1eps2,a2eps2,b1eps2,b2eps2,c1eps2,c2eps2,c3eps2,
     c               ro3eps2,x101eps2,x201eps2,y101eps2,y201eps2
      common/dera1pa2/a1a1pa2,a2a1pa2,b1a1pa2,b2a1pa2,c1a1pa2,c2a1pa2,
     c       c3a1pa2,ro3a1pa2,x101a1pa2,x201a1pa2,y101a1pa2,y201a1pa2
      common/dercr3/a1cr3,a2cr3,b1cr3,b2cr3,c1cr3,c2cr3,c3cr3,ro3cr3,
     c               x101cr3,x201cr3,y101cr3,y201cr3
      common/derdelta/a1delta,a2delta,b1delta,b2delta,c1delta,c2delta,
     c       c3delta,ro3delta,x101delta,x201delta,y101delta,y201delta





      ! fac variatie dupa eps1
       h=1.d-3
         if(dabs(eps1).le.2.d0*h+1.d-4)then
         semn=1.
         if(eps1.lt.0.d0)semn=-1.
      eps1p1=eps1+semn*h+semn*1.d-4
      call parelp2(a0,eps1p1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      eps10=eps1+semn*1.d-4
      call parelp2(a0,eps10,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a10=a1
      a20=a2
      b10=b1
      b20=b2
      c10=c1-cm
      c20=c2-cm
      c30=c3-cm
      x10=x1-cm
      x20=x2-cm
      r30=r3
      u10=u1-cm
      u20=u2-cm
      ro30=ro3
      x1010=x101-cm
      x2010=x201-cm
      y1010=y101
      y2010=y201

       a1eps1=(a1p1-a10)/(semn*h)
       a2eps1=(a2p1-a20)/(semn*h)
       b1eps1=(b1p1-b10)/(semn*h)
       b2eps1=(b2p1-b20)/(semn*h)
       c1eps1=(c1p1-c10)/(semn*h)
       c2eps1=(c2p1-c20)/(semn*h)
       c3eps1=(c3p1-c30)/(semn*h)
       ro3eps1=(ro3p1-ro30)/(semn*h)
       x101eps1=(x101p1-x1010)/(semn*h)
       x201eps1=(x201p1-x2010)/(semn*h)
       y101eps1=(y101p1-y1010)/(semn*h)
       y201eps1=(y201p1-y2010)/(semn*h)
         else
      eps1m2=eps1-2*h
      call parelp2(a0,eps1m2,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      eps1m1=eps1-1*h
      call parelp2(a0,eps1m1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      eps1p1=eps1+1*h
      call parelp2(a0,eps1p1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      eps1p2=eps1+2*h
      call parelp2(a0,eps1p2,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1eps1=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2eps1=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1eps1=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2eps1=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1eps1=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2eps1=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3eps1=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3eps1=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101eps1=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201eps1=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101eps1=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201eps1=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h
         endif


      ! fac variatie dupa eps2
      h=1.d-3
         if(dabs(eps2).le.2.d0*h)then
         semn=1.
         if(eps2.lt.0.d0)semn=-1
      eps2p1=eps2+semn*h+semn*5.d-4
      call parelp2(a0,eps1,eps2p1,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      eps20=eps2+semn*5.d-4
      call parelp2(a0,eps1,eps20,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a10=a1
      a20=a2
      b10=b1
      b20=b2
      c10=c1-cm
      c20=c2-cm
      c30=c3-cm
      x10=x1-cm
      x20=x2-cm
      r30=r3
      u10=u1-cm
      u20=u2-cm
      ro30=ro3
      x1010=x101-cm
      x2010=x201-cm
      y1010=y101
      y2010=y201

       a1eps2=(a1p1-a10)/(semn*h)
       a2eps2=(a2p1-a20)/(semn*h)
       b1eps2=(b1p1-b10)/(semn*h)
       b2eps2=(b2p1-b20)/(semn*h)
       c1eps2=(c1p1-c10)/(semn*h)
       c2eps2=(c2p1-c20)/(semn*h)
       c3eps2=(c3p1-c30)/(semn*h)
       ro3eps2=(ro3p1-ro30)/(semn*h)
       x101eps2=(x101p1-x1010)/(semn*h)
       x201eps2=(x201p1-x2010)/(semn*h)
       y101eps2=(y101p1-y1010)/(semn*h)
       y201eps2=(y201p1-y2010)/(semn*h)
         else
      eps2m2=eps2-2*h
      call parelp2(a0,eps1,eps2m2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      eps2m1=eps2-1*h
      call parelp2(a0,eps1,eps2m1,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      eps2p1=eps2+1*h
      call parelp2(a0,eps1,eps2p1,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      eps2p2=eps2+2*h
      call parelp2(a0,eps1,eps2p2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1eps2=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2eps2=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1eps2=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2eps2=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1eps2=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2eps2=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3eps2=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3eps2=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101eps2=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201eps2=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101eps2=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201eps2=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h
       endif

      ! fac variatie dupa a1pa2 (fii atent a1pa2>1)
      h=1.d-3
      if(a1pa2.lt.1.d0+2*h)then
      a1pa2p1=a1pa2+1*h
      call parelp2(a0,eps1,eps2,a1pa2p1,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
       a1a1pa2=(a1p1-a1)/h
       a2a1pa2=(a2p1-a2)/h
       b1a1pa2=(b1p1-b1)/h
       b2a1pa2=(b2p1-b2)/h
       c1a1pa2=(c1p1-(c1-cm))/h
       c2a1pa2=(c2p1-(c2-cm))/h
       c3a1pa2=(c3p1-(c3-cm))/h
       ro3a1pa2=(ro3p1-ro3)/h
       x101a1pa2=(x101p1-(x101-cm))/h
       x201a1pa2=(x201p1-(x201-cm))/h
      y101a1pa2=(y101p1-y101)/h
      y201a1pa2=(y201p1-y201)/h
      else
      a1pa2m2=a1pa2-2*h
      call parelp2(a0,eps1,eps2,a1pa2m2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      a1pa2m1=a1pa2-1*h
      call parelp2(a0,eps1,eps2,a1pa2m1,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      a1pa2p1=a1pa2+1*h
      call parelp2(a0,eps1,eps2,a1pa2p1,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      a1pa2p2=a1pa2+2*h
      call parelp2(a0,eps1,eps2,a1pa2p2,cr3,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1a1pa2=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2a1pa2=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1a1pa2=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2a1pa2=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1a1pa2=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2a1pa2=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3a1pa2=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3a1pa2=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101a1pa2=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201a1pa2=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101a1pa2=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201a1pa2=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h
      endif

      ! fac variatie dupa cr3 (fii atent la intervale)
      if(dabs(cr3).lt.0.0025)then
                 h=1.d-4
                 cmr3=0.0025
              cmr3p=cmr3
      cr3m2=cmr3-2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      cr3m1=cmr3-1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      cr3p1=cmr3+1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      cr3p2=cmr3+2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1cr3p=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2cr3p=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1cr3p=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2cr3p=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1cr3p=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2cr3p=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3cr3p=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3cr3p=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101cr3p=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201cr3p=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101cr3p=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201cr3p=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h

      cmr3=-0.0025
              cmr3m=cmr3
      cr3m2=cmr3-2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      cr3m1=cmr3-1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      cr3p1=cmr3+1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      cr3p2=cmr3+2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1cr3m=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2cr3m=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1cr3m=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2cr3m=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1cr3m=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2cr3m=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3cr3m=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3cr3m=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101cr3m=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201cr3m=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101cr3m=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201cr3m=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h
         a1cr3=a1cr3m+(a1cr3p-a1cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         a2cr3=a2cr3m+(a2cr3p-a2cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         b1cr3=b1cr3m+(b1cr3p-b1cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         b2cr3=b2cr3m+(b2cr3p-b2cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         c1cr3=c1cr3m+(c1cr3p-c1cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         c2cr3=c2cr3m+(c2cr3p-c2cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         c3cr3=c3cr3m+(c3cr3p-c3cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         ro3cr3=ro3cr3m+(ro3cr3p-ro3cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         x101cr3=x101cr3m+(x101cr3p-x101cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         x201cr3=x201cr3m+(x201cr3p-x201cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         y101cr3=y101cr3m+(y101cr3p-y101cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
         y201cr3=y201cr3m+(y201cr3p-y201cr3m)*(cr3-cmr3m)/(cmr3p-cmr3m)
      else
c deoarece pentru valori dabs(cr3) mici, am mase efective ciudate,
c am sa fixez prin cmr3 o valoare minima a lui dabs(cr3) la 7.d-3
      h=1.d-3
                       cmr3=cr3
       if(dabs(cr3).le.1.d-2)then
       h=1.d-4+(1.d-3-1.d-4)*(dabs(cr3)-0.5d-2)/(1.d-2-0.5d-2)
       if(dabs(cr3).lt.0.5d-2)h=1.d-4
       endif
      cr3m2=cmr3-2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      cr3m1=cmr3-1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3m1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3m1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      cr3p1=cmr3+1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p1,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p1,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      cr3p2=cmr3+2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3p2,delta,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,delta,cr3p2,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1cr3=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2cr3=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1cr3=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2cr3=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1cr3=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2cr3=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3cr3=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3cr3=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101cr3=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201cr3=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101cr3=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201cr3=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h
      endif

      ! fac variatie dupa delta
      h=1.d-3
      deltam2=delta-2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3,deltam2,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,deltam2,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m2=a1
      a2m2=a2
      b1m2=b1
      b2m2=b2
      c1m2=c1-cm
      c2m2=c2-cm
      c3m2=c3-cm
      x1m2=x1-cm
      x2m2=x2-cm
      r3m2=r3
      u1m2=u1-cm
      u2m2=u2-cm
      ro3m2=ro3
      x101m2=x101-cm
      x201m2=x201-cm
      y101m2=y101
      y201m2=y201
      deltam1=delta-1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3,deltam1,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,deltam1,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1m1=a1
      a2m1=a2
      b1m1=b1
      b2m1=b2
      c1m1=c1-cm
      c2m1=c2-cm
      c3m1=c3-cm
      x1m1=x1-cm
      x2m1=x2-cm
      r3m1=r3
      u1m1=u1-cm
      u2m1=u2-cm
      ro3m1=ro3
      x101m1=x101-cm
      x201m1=x201-cm
      y101m1=y101
      y201m1=y201
      deltap1=delta+1*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3,deltap1,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,deltap1,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p1=a1
      a2p1=a2
      b1p1=b1
      b2p1=b2
      c1p1=c1-cm
      c2p1=c2-cm
      c3p1=c3-cm
      x1p1=x1-cm
      x2p1=x2-cm
      r3p1=r3
      u1p1=u1-cm
      u2p1=u2-cm
      ro3p1=ro3
      x101p1=x101-cm
      x201p1=x201-cm
      y101p1=y101
      y201p1=y201
      deltap2=delta+2*h
      call parelp2(a0,eps1,eps2,a1pa2,cr3,deltap2,a1,b1,a2,b2)
      call centrumasa(cm,a1,b1,a2,b2,deltap2,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
      a1p2=a1
      a2p2=a2
      b1p2=b1
      b2p2=b2
      c1p2=c1-cm
      c2p2=c2-cm
      c3p2=c3-cm
      x1p2=x1-cm
      x2p2=x2-cm
      r3p2=r3
      u1p2=u1-cm
      u2p2=u2-cm
      ro3p2=ro3
      x101p2=x101-cm
      x201p2=x201-cm
      y101p2=y101
      y201p2=y201
       a1delta=(a1m2-8*a1m1+8*a1p1-a1p2)/12./h
       a2delta=(a2m2-8*a2m1+8*a2p1-a2p2)/12./h
       b1delta=(b1m2-8*b1m1+8*b1p1-b1p2)/12./h
       b2delta=(b2m2-8*b2m1+8*b2p1-b2p2)/12./h
       c1delta=(c1m2-8*c1m1+8*c1p1-c1p2)/12./h
       c2delta=(c2m2-8*c2m1+8*c2p1-c2p2)/12./h
       c3delta=(c3m2-8*c3m1+8*c3p1-c3p2)/12./h
       ro3delta=(ro3m2-8*ro3m1+8*ro3p1-ro3p2)/12./h
       x101delta=(x101m2-8*x101m1+8*x101p1-x101p2)/12./h
       x201delta=(x201m2-8*x201m1+8*x201p1-x201p2)/12./h
       y101delta=(y101m2-8*y101m1+8*y101p1-y101p2)/12./h
       y201delta=(y201m2-8*y201m1+8*y201p1-y201p2)/12./h



    !! refac valorile initiale
      call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
      a10=a1
      a20=a2
      b10=b1
      b20=b2
      c10=c1
      c20=c2
      c30=c3
      x10=x1
      x20=x2
      r30=r3
      u10=u1
      u20=u2
      ro30=ro3
      return
      end   




      subroutine centrumasa(cm,a1,b1,a2,b2,delta,cr3,c1,c2,
     c       x1,x2,c3,ro3,r3,u1,u2)
c calculez volumul elipsei
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/volume/vol1,vol2,vol31,vol32,vol3
      s=1
      if(cr3.lt.0)s=-1
      pi=3.141592645d0
c volumul intre capatul din stanga pana la x1
      v1=pi*b1**2*(x1-c1+a1-(x1-c1)**3/(3.d0*a1**2)-a1/3.d0)
      rapz=(x1-c1)/a1
!      ai1=pi*a1*b1**2*(a1/2.d0*(rapz**2+1)-a1/4.d0*(rapz**4+1)+
!     c   c1*(rapz+1)-c1/3.d0*(rapz**3+1))-c1*v1  !fata de c1 
      ai1=pi*b1**2*(0.5*((x1-c1)**2-a1**2)
     c            -0.25/a1**2*((x1-c1)**4-a1**4))

c volumul intre capatul din dreapta pana la x2
      v2=pi*b2**2*(-(x2-c2)+a2+(x2-c2)**3/(3.d0*a2**2)-a2/3.d0)
      rapz=(x2-c2)/a2
!      ai2=pi*a2*b2**2*(a2/2.d0*(1-rapz**2)-a2/4.d0*(1-rapz**4)+
!     c   c2*(1-rapz)-c2/3.d0*(1-rapz**3))-c2*v2  !fata de c2
      ai2=pi*b2**2*(0.5*(a2**2-(x2-c2)**2)
     c            -0.25/a2**2*(a2**4-(x2-c2)**4))
      vol1=v1
      vol2=v2
      if(r3.gt.500.d0)then !cazul r3 infinit
      R3=10000.D0
      v3=pi/3.d0*(x201-x101)*(y101**2+y201**2+y101*y201)
      v=v1+v2+v3
      v11=v1+v3
      v22=v2
      v00=v
      vol3=v3
      ray=(y201-y101)/y101
      ai3=pi*y101**2*(x201-x101)**2*(0.5d0+2.d0/3.d0*ray+0.25d0*ray**2) ! fata de x101
      cm=(c1*vol1+ai1+x101*vol3+ai3+c2*vol2+ai2)/(vol1+vol3+vol2)
      return
      endif
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
           a1pa2pr3=2.d0*r3+a1+a2-1.d-6
           if(a1pa2pr3.le.delta)then ! cele doua elipse sunt complet separate
           v1=pi*b1**2*(2.d0*a1-2.d0/3.d0*a1)
           v2=pi*b2**2*(2.d0*a2-2.d0/3.d0*a2)
           v=v1+v2
           vol1=v1
           vol2=v2
           v11=v1
           v22=v2
           v00=v
           cm=(c1*vol1+c2*vol2)/(vol1+vol2)
           return
           endif
       if(ro3.lt.r3.and.s.gt.0.d0.and.c3.gt.c1.and.c3.lt.c2)then !u1<>u2<>c3
      ARGi=(u1 -c3)/R3
      ARG2=(x1-c3)/R3
      V31=(RO3*RO3+R3*R3)*(U1-X1)
     &     +((X1-C3)**3 )/3.D0
     &     +RO3*((X1-C3)*DSQRT(DABS(R3*R3-(X1-C3)**2))
     &    +R3*R3*DASIN(  ARG2)   )
     &      -((U1-C3)**3  )/3.D0
     &   -RO3*((U1 -C3)*DSQRT(DABS(R3*R3-(U1 -C3)**2))
     &    +R3*R3*DASIN(ARGi))
      ARGU=(X2-C3)/R3
      ARGi=(u2-c3)/R3
      V2S=(RO3*RO3+R3*R3)*(X2-U2)
     &      -((X2-C3)**3  )/3.D0
     &   -RO3*((X2-C3)*DSQRT(DABS(R3*R3-(X2-C3)**2))
     &    +R3*R3*DASIN(ARGU))
     &      +((U2 -C3)**3  )/3.D0
     &   +RO3*((U2 -C3)*DSQRT(DABS(R3*R3-(U2 -C3)**2))
     &    +R3*R3*DASIN(ARGi))
        v32=v2s
      aI31=pi*((RO3*RO3+R3*R3 )/2.D0*((U1-c3)**2 -(x1-c3)**2) -
     &   ((U1-c3)**4 -(x1-c3)**4 )/4.D0+
     &  2.D0/3.D0*RO3*R3**3 *((1.D0-(U1-c3)**2 /R3/R3 )**1.5
     &  -(1.D0-(x1-c3)**2 /R3/R3 )**1.5))
      aI32=pi*((RO3*RO3+R3*R3)/2.D0*((x2-c3)**2 -(U2-c3)**2 )-
     &  ((x2-c3)**4 -(U2-c3)**4 )/4.D0+
     &  2.D0/3.D0*RO3*R3**3 *((1.D0-(x2-c3)**2 /R3/R3 )**1.5
     &  -(1.D0-(U2-c3)**2 /R3/R3 )**1.5))
c       print*,'3 v31 v32 ',v31,v32,' u1,u2,c3',u1,u2,c3

       v=v1+v2+pi*(v31+v32)
       v11=v1+pi*v31
       v22=v2+pi*v32
       v00=v
       vol31=pi*v31
       vol32=pi*v32
       cm=(c1*vol1+ai1+c3*(vol31+vol32)+ai31+ai32+c2*vol2+ai2)/
     c    (vol1+vol2+vol31+vol32)
       return
       endif
       vvv1=(r3**2-(x2-c3)**2)
       if(vvv1.lt.1.d-50)vvv1=1.d-50
       vvv2=(r3**2-(x1-c3)**2)
       if(vvv2.lt.1.d-50)vvv2=1.d-50
       vvv4=(x2-c3)/r3
       if(vvv4.gt.1.d0)  vvv4=1.d0
       if(vvv4.lt.-1.d0)  vvv4=-1.d0
       vvv5=(x1-c3)/r3
       if(vvv5.gt.1.d0)  vvv5=1.d0
       if(vvv5.lt.-1.d0)  vvv5=-1.d0
       v3=(ro3**2+r3**2)*(x2-x1)-1.d0/3.d0*((x2-c3)**3-(x1-c3)**3)-
     - s*ro3*((x2-c3)*dsqrt(vvv1)+r3**2*dasin(vvv4)-
     -        (x1-c3)*dsqrt(vvv2)-r3**2*dasin(vvv5))
c    - s*ro3*((x2-c3)*dsqrt(r3**2-(x2-c3)**2)+r3**2*dasin((x2-c3)/r3)-
c    -        (x1-c3)*dsqrt(r3**2-(x1-c3)**2)-r3**2*dasin((x1-c3)/r3))
       v3=pi*v3
       vol3=v3
c       print*,'4 v3',v3
       v=v1+v2+v3
       v00=v
       if(c3.gt.c1.and.c2.gt.c3.and.s.gt.0.d0)then  ! pentru calculul volumelor celor doua elipse
c deosebim 2 cazuri c1<c3<c2  si c2<c3
      v31=(ro3**2+r3**2)*(c3-x1)-1.d0/3.d0*(-(x1-c3)**3)-
     - s*ro3*(-(x1-c3)*dsqrt(r3**2-(x1-c3)**2)-r3**2*dasin((x1-c3)/r3))
       v31=pi*v31
       v32=(ro3**2+r3**2)*(x2-c3)-1.d0/3.d0*((x2-c3)**3)-
     - s*ro3*((x2-c3)*dsqrt(r3**2-(x2-c3)**2)+r3**2*dasin((x2-c3)/r3))
       v32=pi*v32
       vol31=v31
       vol32=v32
      aI31=pi*((RO3*RO3+R3*R3 )/2.D0*((U1-c3)**2 -(x1-c3)**2) -
     &   ((U1-c3)**4 -(x1-c3)**4 )/4.D0+
     &  2.D0/3.D0*RO3*R3**3 *((1.D0-(U1-c3)**2 /R3/R3 )**1.5
     &  -(1.D0-(x1-c3)**2 /R3/R3 )**1.5))
      aI32=pi*((RO3*RO3+R3*R3)/2.D0*((x2-c3)**2 -(U2-c3)**2 )-
     &  ((x2-c3)**4 -(U2-c3)**4 )/4.D0+
     &  2.D0/3.D0*RO3*R3**3 *((1.D0-(x2-c3)**2 /R3/R3 )**1.5
     &  -(1.D0-(U2-c3)**2 /R3/R3 )**1.5))
       cm=(c1*vol1+ai1+c3*(vol31+vol32)+ai31+ai32+c2*vol2+ai2)/
     c    (vol1+vol2+vol31+vol32)
       v11=v1+v31
       v22=v2+v32
       else
       v11=v1+v3
       v22=v2
      ai3=pi*((RO3*RO3 +R3*R3 )/2.D0*((x2-c3)**2 -(x1-c3)**2 )
     &  -1.D0/4.D0*((x2-c3)**4 -(x1-c3)**4 )+S*2.D0/3.D0*
     &   RO3*R3**3*((1.D0-(x2-c3)**2 /R3/R3 )**1.5-
     &   (1.D0-(x1-c3)**2 /R3/R3 )**1.5))
       cm=(c1*vol1+ai1+c3*vol3+ai3+c2*vol2+ai2)/(vol1+vol3+vol2)
       endif
       return
      end


      subroutine ameps1i(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ameps1_0i,ameps1_1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ameps1_1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
      numdqh64=0
           CALL IGAUSS2(ameps1_1i,xl,xu,Y2)
      call dqh64(ameps1_0i,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function ameps1_0i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps64/ameps064(64,0:22,0:22,0:22,0:22)
      external ameps1_11,ameps1_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (ameps1_11,Y)
      ameps1_0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (ameps1_12,Y)
      ameps1_0i=Y*c2*c2p/alf2
      endif
      ameps064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function ameps1_1i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps32/ameps132(64,0:22,0:22,0:22,0:22)
      external ameps1_21,ameps1_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (ameps1_21,Y)
      ameps1_1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (ameps1_22,Y)
c integrala se face pe (alfa*ro)**2
      ameps1_1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      ameps132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end




      subroutine ameps1(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2

      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ameps1_0,ameps1_1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p

      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ameps1_1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(ameps1_1,xl,xu,Y2)
      numdqh64=0
      call dqh64(ameps1_0,pert)
            
      pert=pert+y1+y2

      return
      end

      double  precision function ameps1_0(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps64/ameps064(64,0:22,0:22,0:22,0:22)
      external ameps1_11,ameps1_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1

c           CALL DQL32 (ameps1_11,Y)
      y=ameps064(numdqh64,nro1x,nro2x,m1x,m2x)
      ameps1_0=Y*C1*C1P/alf1
      else

      zp=(z+zet2)/alf2

c           CALL DQL32 (ameps1_12,Y)
      y=ameps064(numdqh64,nro1x,nro2x,m1x,m2x)
      ameps1_0=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function ameps1_1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps32/ameps132(64,0:22,0:22,0:22,0:22)
      external ameps1_21,ameps1_22
      numigauss32=numigauss32+1

           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p

c           CALL DQL32 (ameps1_21,Y)
      y=ameps132(numigauss32,nro1x,nro2x,m1x,m2x)
      ameps1_1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))

      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p

c           CALL DQL32 (ameps1_22,Y)
c integrala se face pe (alfa*ro)**2
      y=ameps132(numigauss32,nro1x,nro2x,m1x,m2x)
      ameps1_1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))

      endif
      return
      end



      double precision function ameps1_11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder11eps1/der11eps1(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

       f=der11eps1(nrnodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps1_11=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ameps1_12(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder12eps1/der12eps1(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro

       f=der12eps1(nrnodzz,nrnod)

      pertur22=f         !*pla1*pla2

      ameps1_12=pertur22*(rho)**(iabs(mp1))

      return
      end


      double precision function ameps1_21(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder21eps1/der21eps1(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der21eps1(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps1_21=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ameps1_22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1

      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder22eps1/der22eps1(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der22eps1(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps1_22=pertur11*(rho)**(iabs(mp1))

      return
      end

      subroutine ama1pa2i(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ama1pa2_0i,ama1pa2_1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ama1pa2_1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(ama1pa2_1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(ama1pa2_0i,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function ama1pa2_0i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ama1p064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ama1p64/ama1p064(64,0:22,0:22,0:22,0:22)
      external ama1pa2_11,ama1pa2_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (ama1pa2_11,Y)
      ama1pa2_0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (ama1pa2_12,Y)
      ama1pa2_0i=Y*c2*c2p/alf2
      endif
      ama1p064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function ama1pa2_1i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ama1p132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ama1p32/ama1p132(64,0:22,0:22,0:22,0:22)
      external ama1pa2_21,ama1pa2_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (ama1pa2_21,Y)
      ama1pa2_1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (ama1pa2_22,Y)
c integrala se face pe (alfa*ro)**2
      ama1pa2_1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      ama1p132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end



      subroutine ama1pa2(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ama1pa2_0,ama1pa2_1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p

      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ama1pa2_1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(ama1pa2_1,xl,xu,Y2)
      numdqh64=0
      call dqh64(ama1pa2_0,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function ama1pa2_0(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ama1p064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ama1p64/ama1p064(64,0:22,0:22,0:22,0:22)
      external ama1pa2_11,ama1pa2_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c           CALL DQL32 (ama1pa2_11,Y)
      y=ama1p064(numdqh64,nro1x,nro2x,m1x,m2x)
      ama1pa2_0=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
c           CALL DQL32 (ama1pa2_12,Y)
      y=ama1p064(numdqh64,nro1x,nro2x,m1x,m2x)
      ama1pa2_0=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function ama1pa2_1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ama1p132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ama1p32/ama1p132(64,0:22,0:22,0:22,0:22)
      external ama1pa2_21,ama1pa2_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c           CALL DQL32 (ama1pa2_21,Y)
      y=ama1p132(numigauss32,nro1x,nro2x,m1x,m2x)
      ama1pa2_1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c           CALL DQL32 (ama1pa2_22,Y)
c integrala se face pe (alfa*ro)**2
      y=ama1p132(numigauss32,nro1x,nro2x,m1x,m2x)
      ama1pa2_1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end



      double precision function ama1pa2_11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder11a1pa2/der11a1pa2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

       f=der11a1pa2(nrnodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ama1pa2_11=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ama1pa2_12(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder12a1pa2/der12a1pa2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro

       f=der12a1pa2(nrnodzz,nrnod)

      pertur22=f         !*pla1*pla2

      ama1pa2_12=pertur22*(rho)**(iabs(mp1))

      return
      end


      double precision function ama1pa2_21(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder21a1pa2/der21a1pa2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der21a1pa2(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ama1pa2_21=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ama1pa2_22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder22a1pa2/der22a1pa2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der22a1pa2(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ama1pa2_22=pertur11*(rho)**(iabs(mp1))

      return
      end

      
      subroutine ameps2i(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ameps2_0i,ameps2_1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ameps2_1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(ameps2_1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(ameps2_0i,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function ameps2_0i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps264/ameps064(64,0:22,0:22,0:22,0:22)
      external ameps2_11,ameps2_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (ameps2_11,Y)
      ameps2_0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (ameps2_12,Y)
      ameps2_0i=Y*c2*c2p/alf2
      endif
      ameps064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function ameps2_1i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps232/ameps132(64,0:22,0:22,0:22,0:22)
      external ameps2_21,ameps2_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (ameps2_21,Y)
      ameps2_1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (ameps2_22,Y)
c integrala se face pe (alfa*ro)**2
      ameps2_1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      ameps132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end



      subroutine ameps2(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2

      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external ameps2_0,ameps2_1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p

      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(ameps2_1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(ameps2_1,xl,xu,Y2)
      numdqh64=0
      call dqh64(ameps2_0,pert)
            
      pert=pert+y1+y2

      return
      end

      double  precision function ameps2_0(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps264/ameps064(64,0:22,0:22,0:22,0:22)
      external ameps2_11,ameps2_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1

c           CALL DQL32 (ameps2_11,Y)
      y=ameps064(numdqh64,nro1x,nro2x,m1x,m2x)
      ameps2_0=Y*C1*C1P/alf1
      else

      zp=(z+zet2)/alf2

c           CALL DQL32 (ameps2_12,Y)
      y=ameps064(numdqh64,nro1x,nro2x,m1x,m2x)
      ameps2_0=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function ameps2_1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real ameps132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/ameps232/ameps132(64,0:22,0:22,0:22,0:22)
      external ameps2_21,ameps2_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p

c           CALL DQL32 (ameps2_21,Y)
      y=ameps132(numigauss32,nro1x,nro2x,m1x,m2x)
      ameps2_1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))

      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p

c           CALL DQL32 (ameps2_22,Y)
c integrala se face pe (alfa*ro)**2
      y=ameps132(numigauss32,nro1x,nro2x,m1x,m2x)
      ameps2_1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))

      endif
      return
      end



      double precision function ameps2_11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder11eps2/der11eps2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

       f=der11eps2(nrnodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps2_11=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ameps2_12(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder12eps2/der12eps2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro

       f=der12eps2(nrnodzz,nrnod)

      pertur22=f         !*pla1*pla2

      ameps2_12=pertur22*(rho)**(iabs(mp1))

      return
      end


      double precision function ameps2_21(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder21eps2/der21eps2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der21eps2(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps2_21=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function ameps2_22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder22eps2/der22eps2(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der22eps2(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      ameps2_22=pertur11*(rho)**(iabs(mp1))

      return
      end


      subroutine amcr3i(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external amcr3_0i,amcr3_1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(amcr3_1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(amcr3_1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(amcr3_0i,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function amcr3_0i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amcr3064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amcr364/amcr3064(64,0:22,0:22,0:22,0:22)
      external amcr3_11,amcr3_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (amcr3_11,Y)
      amcr3_0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (amcr3_12,Y)
      amcr3_0i=Y*c2*c2p/alf2
      endif
      amcr3064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function amcr3_1i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amcr3132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amcr332/amcr3132(64,0:22,0:22,0:22,0:22)
      external amcr3_21,amcr3_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (amcr3_21,Y)
      amcr3_1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (amcr3_22,Y)
c integrala se face pe (alfa*ro)**2
      amcr3_1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      amcr3132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end



      subroutine amcr3(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external amcr3_0,amcr3_1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(amcr3_1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(amcr3_1,xl,xu,Y2)
      numdqh64=0
      call dqh64(amcr3_0,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function amcr3_0(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amcr3064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amcr364/amcr3064(64,0:22,0:22,0:22,0:22)
      external amcr3_11,amcr3_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c           CALL DQL32 (amcr3_11,Y)
      y=amcr3064(numdqh64,nro1x,nro2x,m1x,m2x)
      amcr3_0=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
c           CALL DQL32 (amcr3_12,Y)
      y=amcr3064(numdqh64,nro1x,nro2x,m1x,m2x)
      amcr3_0=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function amcr3_1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amcr3132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amcr332/amcr3132(64,0:22,0:22,0:22,0:22)
      external amcr3_21,amcr3_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c           CALL DQL32 (amcr3_21,Y)
      y=amcr3132(numigauss32,nro1x,nro2x,m1x,m2x)
      amcr3_1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c           CALL DQL32 (amcr3_22,Y)
c integrala se face pe (alfa*ro)**2
      y=amcr3132(numigauss32,nro1x,nro2x,m1x,m2x)
      amcr3_1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end



      double precision function amcr3_11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder11cr3/der11cr3(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

       f=der11cr3(nrnodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amcr3_11=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function amcr3_12(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder12cr3/der12cr3(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro

       f=der12cr3(nrnodzz,nrnod)

      pertur22=f         !*pla1*pla2

      amcr3_12=pertur22*(rho)**(iabs(mp1))

      return
      end


      double precision function amcr3_21(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder21cr3/der21cr3(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der21cr3(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amcr3_21=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function amcr3_22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder22cr3/der22cr3(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der22cr3(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amcr3_22=pertur11*(rho)**(iabs(mp1))

      return
      end



      subroutine amdelti(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external amdelt_0i,amdelt_1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(amdelt_1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(amdelt_1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(amdelt_0i,pert)
            
      pert=pert+y1+y2
      return
      end

      double  precision function amdelt_0i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amdel064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amdel64/amdel064(64,0:22,0:22,0:22,0:22)
      external amdelt_11,amdelt_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (amdelt_11,Y)
      amdelt_0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (amdelt_12,Y)
      amdelt_0i=Y*c2*c2p/alf2
      endif
      amdel064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function amdelt_1i(z)

c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amdel132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amdel32/amdel132(64,0:22,0:22,0:22,0:22)
      external amdelt_21,amdelt_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (amdelt_21,Y)
      amdelt_1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (amdelt_22,Y)
c integrala se face pe (alfa*ro)**2
      amdelt_1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      amdel132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end

  

      subroutine amdelt(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external amdelt_0,amdelt_1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(amdelt_1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(amdelt_1,xl,xu,Y2)
      numdqh64=0
      call dqh64(amdelt_0,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function amdelt_0(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amdel064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amdel64/amdel064(64,0:22,0:22,0:22,0:22)
      external amdelt_11,amdelt_12
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c           CALL DQL32 (amdelt_11,Y)
      y=amdel064(numdqh64,nro1x,nro2x,m1x,m2x)
      amdelt_0=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
c           CALL DQL32 (amdelt_12,Y)
      y=amdel064(numdqh64,nro1x,nro2x,m1x,m2x)
      amdelt_0=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function amdelt_1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real amdel132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/amdel32/amdel132(64,0:22,0:22,0:22,0:22)
      external amdelt_21,amdelt_22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c           CALL DQL32 (amdelt_21,Y)
      y=amdel132(numigauss32,nro1x,nro2x,m1x,m2x)
      amdelt_1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c           CALL DQL32 (amdelt_22,Y)
c integrala se face pe (alfa*ro)**2
      y=amdel132(numigauss32,nro1x,nro2x,m1x,m2x)
      amdelt_1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end



      double precision function amdelt_11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder11delt/der11delt(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

       f=der11delt(nrnodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amdelt_11=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function amdelt_12(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/dder12delt/der12delt(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro

       f=der12delt(nrnodzz,nrnod)

      pertur22=f         !*pla1*pla2

      amdelt_12=pertur22*(rho)**(iabs(mp1))

      return
      end


      double precision function amdelt_21(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder21delt/der21delt(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der21delt(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amdelt_21=pertur11*(rho)**(iabs(mp1))

      return
      end

      double precision function amdelt_22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/dder22delt/der22delt(32,32)

      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro

      f=der22delt(nodzz,nrnod)

      pertur11=f                     !*pla1*pla2

      amdelt_22=pertur11*(rho)**(iabs(mp1))

      return
      end

    
!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!! adaug aici un potential coulombian pe care il diagonalizez
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! calculez potentialul Coulomb !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! intra subrutinele qgauss si delke !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine coulini(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/zilim12/z1lim,z2lim
      common/numintegr/numigauss32,numdqh64
      external couli1i,coulo0i  
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(couli1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(couli1i,xl,xu,Y2)
      pert=y1+y2
      numdqh64=0
           CALL DQH64(coulo0i,pp)
      pert=pert+pp
      return
      end

      double  precision function coulo0i(z)
      implicit double precision (a-h,o-z)
      real coulo064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/coulo64/coulo064(64,0:22,0:22,0:22,0:22)
      external coul01,coul02
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (coul01,Y)
      coulo0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (coul02,Y)
      coulo0i=Y*c2*c2p/alf2
      endif
      coulo064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function couli1i(z)
      implicit double precision (a-h,o-z)
      real coulo132
      double precision lzsz11,lzsz22
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/coulo32/coulo132(64,0:22,0:22,0:22,0:22)
      external coul11,coul22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (coul11,Y)
      couli1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (coul22,Y)
c integrala se face pe (alfa*ro)**2
      couli1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      coulo132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end


      subroutine coulin(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
      double precision lzsz1
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/zilim12/z1lim,z2lim
      common/numintegr/numigauss32,numdqh64
      external couli1,coulo0  
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
c      z1lim=cc1-a1x-4.d0
c      z2lim=cc2+a2x+4.d0
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(couli1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(couli1,xl,xu,Y2)
      pert=y1+y2
      numdqh64=0
           CALL DQH64(coulo0,pp)
      pert=pert+pp
      return
      end

      double  precision function coulo0(z)
      implicit double precision (a-h,o-z)
      double precision lzsz11,lzsz22
      real coulo064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/coulo64/coulo064(64,0:22,0:22,0:22,0:22)
      external coul01,coul02
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c           CALL DQL32 (coul01,Y)
      y=coulo064(numdqh64,nro1x,nro2x,m1x,m2x)
      coulo0=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
c           CALL DQL32 (coul02,Y)
      y=coulo064(numdqh64,nro1x,nro2x,m1x,m2x)
      coulo0=Y*c2*c2p/alf2
      endif
      return
      end

      double precision function coul01(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delta
      common/CKST/CKST
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/volumel/v00,v11,v22
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vcoul1ext/coulo1ext(32,32) 
      common/lasciziune/delatsciz  
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
ccc      call vc(z,ro,f)
           f=coulo1ext(nrnodzz,nrnod)
! trebuie densitate *e  ro*e=Z*e**2/V*r**2 (r**2 vine din integrala)
! e**2=1.4399764 MeV fm
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      RO0=zZ0/aA0
      RO1F=(zZ0-Zz2)/(aA0-aA2)
      RO2F=zZ2/aA2
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z1M=zZ0-RO2*V22/V00*aA0
      A1M=aA0-V22/V00*aA0
      RO1=Z1M/A1M
      ENDIF
      ENDIF
! ro*e**2=Z/A*e**2/(4*pi/3*r0**3)
        cost=ro1*1.4399764/(4.d0*3.141592645d0/3.d0*r0**3)
      f=cost*f
       coul01=f*(rho)**(iabs(mp1))
      return
      end

      double precision function coul02(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delta
      common/CKST/CKST
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/volumel/v00,v11,v22
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vcoul2ext/coulo2ext(32,32)  
      common/lasciziune/delatsciz  
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
ccc      call vc(z,ro,f)
           f=coulo2ext(nrnodzz,nrnod)
! trebuie densitate *e  ro*e=Z*e**2/V*r**2 (r**2 vine din integrala)
! e**2=1.4399764 MeV fm
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      RO0=zZ0/aA0
      RO1F=(zZ0-Zz2)/(aA0-aA2)
      RO2F=zZ2/aA2
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z1M=zZ0-RO2*V22/V00*aA0
      A1M=aA0-V22/V00*aA0
      RO1=Z1M/A1M
      ENDIF
      ENDIF
! ro*e**2=Z/A*e**2/(4*pi/3*r0**3)
        cost=ro2*1.4399764/(4.d0*3.141592645d0/3.d0*r0**3)
      f=cost*f
       coul02=f*(rho)**(iabs(mp1))
      return
      end




      double  precision function couli1(z)
      implicit double precision (a-h,o-z)
      double precision lzsz11,lzsz22
      real coulo132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/coulo32/coulo132(64,0:22,0:22,0:22,0:22)
      external coul11,coul22
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
c           CALL DQL32 (coul11,Y)
      y=coulo132(numigauss32,nro1x,nro2x,m1x,m2x)
      couli1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
c           CALL DQL32 (coul22,Y)
c integrala se face pe (alfa*ro)**2
      y=coulo132(numigauss32,nro1x,nro2x,m1x,m2x)
      couli1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end

      double precision function coul11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delta
      common/CKST/CKST
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/volumel/v00,v11,v22
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vcoul1/coulo1(32,32)  
      common/lasciziune/delatsciz    

c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
ccc      call vc(z,ro,f)
           f=coulo1(nodzz,nrnod)
! trebuie densitate *e  ro*e=Z*e**2/V*r**2 (r**2 vine din integrala) 
! e**2=1.4399764 MeV fm
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2
      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      RO0=zZ0/aA0
      RO1F=(zZ0-Zz2)/(aA0-aA2)
      RO2F=zZ2/aA2
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z1M=zZ0-RO2*V22/V00*aA0
      A1M=aA0-V22/V00*aA0
      RO1=Z1M/A1M
      ENDIF
      ENDIF
! ro*e**2=Z/A*e**2/(4*pi/3*r0**3)
        cost=ro1*1.4399764/(4.d0*3.141592645d0/3.d0*r0**3)
      f=cost*f 
c      n=n1
c      k=mp1
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur11=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
       coul11=pertur11*(rho)**(iabs(mp1))
c      coul11=pertur11*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))/ro !*dexp(-rho)
      return
      end

      double precision function coul22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delta
      common/CKST/CKST
      common/micoul/aa0,zz0,aa2,zz2,volu0
      common/volumel/v00,v11,v22
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vcoul2/coulo2(32,32) 
      common/lasciziune/delatsciz      

c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro
ccccc      call vc(z,ro,f)
           f=coulo2(nodzz,nrnod)
! trebuie densitate *e  ro*e=Z*e**2/V*r**2 (r**2 vine din integrala) 
! e**2=1.4399764 MeV fm
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz  
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      RO0=zZ0/aA0
      RO1F=(zZ0-Zz2)/(aA0-aA2)
      RO2F=zZ2/aA2
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z1M=zZ0-RO2*V22/V00*aA0
      A1M=aA0-V22/V00*aA0
      RO1=Z1M/A1M
      ENDIF
      ENDIF
! ro*e**2=Z/A*e**2/(4*pi/3*r0**3)
        cost=ro1*1.4399764/(4.d0*3.141592645d0/3.d0*r0**3)
      f=cost*f 
c      n=n1
c      k=mp1
c      rho2=rho**2
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur22=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
       coul22=pertur22*(rho)**(iabs(mp1))
c      coul22=pertur22*2.d0*alf2ro**2*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c dsqrt((rho)**(iabs(mp1)+iabs(mp1p)))*dsqrt(rho)/ro !*dexp(-rho)
      return
      end

      subroutine vc(zi,ro,fs)
C calculez potentialul Coulomb fata de z si rho
      implicit double precision (a-h,o-z)
      external vcinteg
      common/zzroro/zz,rro
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/cx12/cx1,cx2
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
!
      z=zi-z0p
      raza0=raza00
      zz=z
      rro=ro 
        if(u1.lt.u2)then
      ali=cc1-a1
      als=u1
      call qgauss(vcinteg,ali,als,nn,y)
      ali=u2
      als=cc2+a2
      call qgauss(vcinteg,ali,als,nn,y1)
      y=y+y1
        else
      ali=cc1-a1
      als=cc2+a2
      call qgauss(vcinteg,ali,als,nn,y)
         endif
      fs=y
      return
      end


      double precision function vcinteg(zprim)
c integrant pentru potentialul Coulomb uniparticula
      implicit double precision (a-h,o-z)
      common/zzroro/zz,rro
      z=zz
      ro=rro
! subrutina care calculeaza roprim si derivata lui roprim**2 fata de zprim
      call PARELC(Zprim,dROpri,ROprim)
      a=dsqrt((zprim-z)**2+(roprim+ro)**2)
      b=dsqrt((zprim-z)**2+(roprim-ro)**2)
      ak2=(a**2-b**2)/a**2
            CALL DELKE(ak2,DK,DE)
      fab=dk/a
      eab=de*a
c     vcinteg=(roprim**2-ro**2-(zprim-z)**2-(zprim-z)*(-2*zprim))*fab !pentru sfera
      vcinteg=(roprim**2-ro**2-(zprim-z)**2-(zprim-z)*(dropri))*fab
     +        +eab
      return
      end

      SUBROUTINE PARELC(Z,RO2DZ,ROFZ)
C ACEASTA SUBROUTINA NE DA DEPENDENTA VARIABILEI
C RO LA PATRAT ( dRO**2Z) FATA DE VARIABILA Z
C SI A VARIABILEI RO FATA DE Z FARA NORMALIZARE PENTRU POTENTIALUL COULOMB
C DIN HAMILTONIANUL UNIPARTICULA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/limitesol/x101,x201,y101,y201

C ACEASTA FUNCTIE ARE DIFERITE FORME IN FUNCTIE DE INTERVALUL
C PE CARE ESTE DEFINITA. AVEM TREI INTERVALE [-1,ZC1],
C [ZC1,ZC2], [ZC2,1]. IN INTERVALUL [ZC1,ZC2] AVEM
C SITUATIE SPECIALA CIND R3=FORV(4) (NENORMALIZAT) TINDE LA
C INFINIT
      IF (Z.LE.X1) GO TO 200
      IF (Z.LT.X2) GO TO 100
C FIE INTERVALUL [X2N,1]
      ZZ2=(Z-C2)/A2**2
      RO2DZ=B2**2*(-ZZ2-ZZ2)
      ROFZ=B2*DSQRT(1.D0-((Z-C2)/A2)**2)
      RETURN
C FIE INTERVALUL [-1,X1N]
200   ZZ1=(Z-C1)/A1**2
      RO2DZ=B1**2*(-ZZ1-ZZ1)
      ROFZ=B1*DSQRT(1.D0-((Z-C1)/A1)**2)
      RETURN
C FIE INTERVALUL [X1N,X2N], AVEM SI POSIBILITATEA CA R3.GT.1.D6
100   CONTINUE
      DABR3=DABS(R3)
      IF( DABR3 .GT. 500.D0) GO TO 300
      ZZ3=Z-C3
      EXI=S*DSQRT(R3*R3-ZZ3*ZZ3)
      RO2DZ=(-2.D0+2.D0*RO3/EXI)*ZZ3
      ROFZ=RO3-EXI
      RETURN
C CAZUL R3 FOARTE MARE
300   CONTINUE
      ROFZ=Y101+(Y201-Y101)*(Z-X101)/(X201-X101)
      RO2DZ=2.D0*ROFZ*(Y201-Y101)/(X201-X101)
      RETURN
      END



      double precision function ceifk(ak)
c compute the eliptic integral of first kind up to pi/2
      implicit double precision (a-h,o-z)
      pi=3.141592645d0/2.d0
      sum=1
      dfm1=-1.d0 
      df=0.d0
      x=1.d0
      s=1.d0
1     continue
      dfm1=dfm1+2.d0
      df=df+2.d0
      s=s*(dfm1/df)**2*ak**2
      sum=sum+s
      prec=dabs(sum*1.d-12)
      ds=dabs(s)
      if(ds.gt.prec)goto 1
      ceifk=pi*sum
      return
      end

      double precision function ceisk(ak)
c compute the eliptic integral of second kind up to pi/2
      implicit double precision (a-h,o-z)
      pi=3.141592645d0/2.d0
      sum=1
      dfm1=1.d0
      df=2.d0
      x=1.d0
      s=1.d0
      n=1
1     continue
      n=n+1
      an=2*n-1.d0
      dfm1=dfm1+2.d0
      df=df+2.d0
      s=s*(dfm1/df)**2*ak**2/an
      sum=sum+s
      prec=dabs(sum*1.d-12)
      ds=dabs(s)
      if(ds.gt.prec)goto 1
      ceisk=pi*(1.d0-sum)
      return
      end

      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANA AICI COULOMB

      subroutine parw_s(i,a0,z0,alp,aln,rsp,rsn,r0p,r0n,ak,vp,vn,a)
      implicit double precision (a-h,o-z)
c ne da parametrii wood saxon functie de optiunea i
c i=101 ROST, 104 Universal, in rest Blomqvist si Wahlborn
       if(i.eq.101)goto 101
       if(i.eq.104)goto 104
       alp=32.d0 !(lambda protoni) - constanta spin orbita
       aln=32.d0 !(lambda neutroni)
       rsp=1.27d0 ! r0 pentru spin orbita protoni
       rsn=1.27d0 ! r0 pentru spin orbita neutroni
       r0p=1.27d0 ! r0 pentru protoni
       r0n=1.27d0 ! r0 pentru neutroni
       ak=0.67 ! parametrizeaza adancimea V=V0*(1\pm ak(N-Z)/(N+Z)) plus=protoni minus neutroni
       v0=51.d0
       vp=v0*(1.d0+ak*(a0-z0-z0)/a0)
       vn=v0*(1.d0-ak*(a0-z0-z0)/a0)
       a=0.67 ! constanta difuseness parameter of the central potential
       return
101    continue
       alp=17.8d0
       aln=31.5d0
       rsp=0.932d0
       rsn=1.280d0
       r0p=1.275d0
       r0n=1.347d0
       ak=0.86 
       v0=49.60
       vp=v0*(1.d0+ak*(a0-z0-z0)/a0)
       vn=v0*(1.d0-ak*(a0-z0-z0)/a0)
       a=0.70
       return
104    continue
       alp=36.0d0
       aln=35.0d0
       rsp=1.200d0
       rsn=1.310d0
       r0p=1.275d0
       r0n=1.347d0
       ak=0.86
       v0=49.60
       vp=v0*(1.d0+ak*(a0-z0-z0)/a0)
       vn=v0*(1.d0-ak*(a0-z0-z0)/a0)
       a=0.70
       return
       end



      SUBROUTINE ENERG(A0,Z0,A2E,Z2E,NMAX,IZ,NRNIV,EN,NFIT,
     C                  IAS,infin,nrspin,nnlel)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PARAMETRII DE INTRARE SINT A0-NR.MASA PARINTE, Z0-NR.ATOMIC PARINTE
C Z2E-NR.ATOMIC FRAGMENT EMIS, A2E-NR.MASA FRAGMENT EMIS,
C NMAX<25 Numar de nivele majore de la 1 la 25
C IZ OPTIUNE 0 pentru protoni, alt numar intreg in rest
C NRNIV ESTE NUMARUL DE ORBITALI PENTRU NMAX (FIECARE VA FI OCUPAT CU
C DOI NUCLEONI
C PARAMETRII DE IESIRE EN(2925) NIVELE -ENERGETICE (EN(2925=NRNIV MAXIM))
c daca nfit = 5 atunci se face fitarea initiala
c infin(2925) indexare nivele cu 1 sau 2 (prima sau a 2 a groapa)
c NMAX=25
      DIMENSION EN(2925),CO1D(0:24),CO2D(0:24),X1D(0:24),
     C          X2D(0:24)
      dimension infin(2925),nrspin(2925),nordin(2925)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/REOIER/IER
      common/infin/infinu(2925)
      common/r116/R0
      common/CKST/CKST
      common/adifz/adiz
      common/micoul/aa0,zz0,aa2,zz2,volu0


      aa0=a0
      zz0=z0
      aa2=a2e
      zz2=z2e
      nrniv=nmax*(nmax+1)*(nmax+2)/6 ! NUMAR DE NIVELE CALCULATE
      PI=3.14159265358979324D0
c     R0=1.16D0
      CKST=41.D0*1.2
      if(iz.eq.1)CKST=41.d0
      i=1
      call parw_s(i,a0,z0,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      volu0=r0n0*a0**(0.333333)*4.*PI/3.d0
      if(iz.eq.0)volu0=r0p0*a0**(0.333333)*4.*PI/3.d0
!calculez constantele spin orbita si l2
c in general variatie cu masa nu are decat adancimea potentialului
      adiz=adz
          if(iz.eq.0)then
      alp0=alp0
      raso0=rsp0 
      ra0=r0p0
      v00=vp0
          else
      alp0=aln0
      raso0=rsn0 
      ra0=r0n0
      v00=vn0
          endif
      ai=a2e
      zi=z2e
      call parw_s(i,ai,zi,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      adiz=adz
          if(iz.eq.0)then
      alp2=alp0
      raso2=rsp0 
      ra2=r0p0
      v02=vp0
          else
      alp2=aln0
      raso2=rsn0 
      ra2=r0n0
      v02=vn0
          endif
      ai=a0-a2e
      zi=z0-z2e
      call parw_s(i,ai,zi,alp0,aln0,rsp0,rsn0,r0p0,r0n0,ak0,vp0,vn0,adz)
      adiz=adz
      z1e=z0-z2e
      a1e=a0-a2e
          if(iz.eq.0)then
      alp1=alp0
      raso1=rsp0 
      ra1=r0p0
      v01=vp0
           else
      alp1=aln0
      raso1=rsn0 
      ra1=r0n0
      v01=vn0
           endif
      npr=0
c trebuie sa calculez r0f astfel incat r1f=a1 si r2f=a2 pentru delt
      r1f=a1
      r2f=a2 
c         xrapx=a2/a1
c          if(xrapx.lt.1.04d0.and.xrapx.gt.0.96d0)then
c     CALL NENGLScumlimita(A0,a2e,R1F,R2F,V00,V01,V02,alp0,raso0,ra0,
c    C     NMAX,IZ,EN,X1D,X2D,CO1D,CO2D,npr,
c    c     nfit,IAS,nrspin,nnlel)
c          endif

      CALL NENGLScumef(A0,a2e,R1F,R2F,V00,V01,V02,alp0,raso0,ra0,
     C     NMAX,IZ,EN,X1D,X2D,CO1D,CO2D,npr,
     c     nfit,IAS,nrspin,nnlel,infin)
!      do i9=1,nnlel!nrniv
!      infin(i9)=infinu(i9)
!     print*,'en(i9)',en(i9),i9
!      enddo
!      call donare(nrniv,en,infin,nrspin,nordin)
C EN SUNT ENERGIILE
      RETURN
      END


      subroutine laguerre(x,n,ik,pla)
c calculeaza pla=L_{n}^{k}(x) (polinom laguerre asociat)
c ordin maxim 100
      implicit double precision (a-h,o-z)
      k=iabs(ik)
      npk=n+k
      call vect(npk,anpk) ! inmultesc la sfarsit cu (n+k)!
      call vect(n,an)
      call vect(k,ak)
      am=1.d0
      isemn=1
      sum=1.d0/(an*ak)
      xx=1.d0
      do m=1,n
      isemn=-isemn
      am=am*m
      kk=k+m
      ak=ak*kk
      nn=n-m+1
      an=an/nn
      xx=x*xx
      sum=sum+isemn*xx/(an*ak*am)
      enddo
      pla=sum*anpk
      return
      end

      double precision function dif_ws_os(z,rho)
c calculeaza diferenta dintre wood-saxon si oscilator
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2 ! valori calculate in NENGLS
      COMMON/ALF12/ALF1,ALF2
      common/r116/R0    
      common/oadancime/v000
      common/CKST/CKST
      common/adancimea/vadanc
      common/cx12/cx1,cx2
! caut ws
 
      woods=ws(z,rho)
c        woods=0
c     CKST=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      W2RO=W0*R0/b2
      ALF1ro=DSQRT(MPH*W1ro)*1.D-15
      ALF2ro=DSQRT(MPH*W2ro)*1.D-15
      if(z.le.0.d0)then
      Eosc=0.5d0*(hwz1*alf1**2*(z+dabs(cx1))**2+hwro1*alf1ro**2*rho**2)
      else
      Eosc=0.5d0*(hwz2*alf2**2*(z-dabs(cx2))**2+hwro1*alf1ro**2*rho**2)
      endif

c      print*,'eosc',eosc,'hwz1,hwz2,hwro1',hwz1,hwz2,hwro1
c      print*,'alf1,alf2,alf1ro', alf1,alf2,alf1ro      
      
      dif_ws_os=woods-eosc+vadanc
    
c      print*,'dif,woods,eosc,v000',dif_ws_os,woods,eosc,vadanc                  
      return
      end




      double precision function rhodens(zi,rho)
c ne da valorile potentialului de imperechere functie de z si ro (coordonatele cilindrice)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/prtdta/z0,ro0 
      common/limitesol/x101,x201,y101,y201
      common/dltw/delta
      common/oadancime/v000
      common/adifz/adiz
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
!
      external delta1,delta2

c  rhoc densitatea la saturatie, vrho este valoarea cu care se imnulteste,
c gamma este exponentul
      rhoc=0.16
      vrho=-999
      gamma=1
      vp1=rhoc
      vp2=rhoc


        z=zi-z0p
      V0=vp1
      ca=adiz
c     print*,'vp1,vp2,adiz',vp1,vp2,adiz
c deosebim mai multe cazuri
      ro=dabs(rho)
      z0=z
      ro0=ro
      a1pa2pr3=a1+a2+2*r3-1.d-6
      if(a1pa2pr3.lt.delt)then ! cazul elipselor complet separate
! in acest caz utilizam planul de separare z=0
         if(z.lt.0.d0)then
! elipsa din partea stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
! pana aici elipsa stanga
          else
! elipsa din partea dreapta (2)
        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
          endif
          endif
! pana aici elipsa dreapta
 
          endif
        delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         v000=V0
         ws=rhoc/(1+dexp(delta/ca))
        rhodens=vrho*(1.d0-(ws/rhoc)**gamma)
      return
! Pana aici elipse complet separate
      endif
      if(r3.ge.500.d0)then ! cazul r3 infinit
! in acest caz trebuie calculate liniile de separare intre cele trei zone
!      alpha=atan((y101-201)/(x201-x101))
!      alpha=pi/2-gamma ! gamma da panta planelor de separate tan(gamma)=cotg(alpha)
! tan(gamma)=(x201-x101)/(y101-y201)
      !calculez primul segment de separare 
        if(y101.eq.y201)then
        if(z.le.x101)goto 4441 ! atunci am prima elipsa
        if(z.ge.x201)goto 4443 ! atunci am a doua elipsa
        goto 4442 ! cazul intermediar
        endif
       tang=(x201-x101)/(y101-y201)
      x0prim=x101-1.d0/tang*y101 ! calculat originea
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
          if(rodif.le.ro)goto 4441 ! atunci am prima elipsa 
      x0prim=x201-1.d0/tang*y201 ! calculat originea pentru a doua separare
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
           if(rodif.ge.ro)goto 4443 ! atunci am a doua elipsa

           goto 4442

4441      continue
! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
! pana aici elipsa din stanga 1

           goto 4445
4443       continue
! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         
! pana aici elipsa din dreapta

           goto 4445
4442  continue ! a ramas partea intermediara
      z0=z
      ro0=ro
      if(y101.eq.y201)then
      delta=ro0-y101
      else
      alpha=datan((y101-y201)/(x201-x101))
      roc=y101+(z-x101)*(y201-y101)/(x201-x101)!calculez punctul corespunzator pe suprafata 
      delta=(ro0-roc)*dcos(alpha)
      endif
4445  continue
         v000=V0
         ws=rhoc/(1+dexp(delta/ca))
        rhodens=vrho*(1.d0-(ws/rhoc)**gamma)
c              write(33,*)'crau'
      return
! pana aici cazul r3 infinit
      endif


      if(s.gt.0.d0)then
      if(dabs(ro).gt.ro3)then
      xx02=c3
      goto 6688
      endif
! calculez y11 si y22 punctele de pe suprafata elipsei la trecerea
! dintre regiuni
! folosesc rapoarte pentru determinarea valorilor x functie de ro
! pe suprafetele de separare
! daca z se afla intre cele doua limite atunci am cazul regiunii mediane
                       ! s pozitiv
       y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
       xx01=(x1-c3)*(ro3-dabs(ro))/(ro3-y11)+c3

       y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)
       xx02=(x2-c3)*(ro3-dabs(ro))/(ro3-y22)+c3
           if(r3.lt.1.d-3)then
           xx01=c3
           xx02=c3
           endif
     
         if(z.ge.xx01.and.z.le.xx02)goto 6677 ! regiunea mediana

             else
                        ! s negativ
      y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
      y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)

      xx01=(x1-c3)*(dabs(ro)+s*ro3)/(s*ro3+y11)+c3
      xx02=(x2-c3)*(dabs(ro)+s*ro3)/(s*ro3+y22)+c3
         if(z.le.xx02.and.z.ge.xx01)goto 6677 ! regiunea mediana
             endif
     
      goto 6688 ! mergi intr-una din elipse


6677        continue ! regiunea mediana
      V0=(vp1+vp2)/2.d0
      ! ramane cazul intermediar pentru r3
      ! ro3 are semn negativ pentru s=-1 si pozitiv pentru s=1
                delta=dsqrt((c3-z)**2+(dabs(ro)-ro3)**2)-r3
                if(s.gt.0.d0)delta=-delta
                goto 5555
6688  continue
c      if(z.gt.z3.and.dabs(ro).gt.ro3)goto 5553
      if(z.ge.xx02)goto 5553
!5551  continue      

! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
! pana aici elipsa din stanga 1

      goto 5555

5553  continue      

! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         

      goto 5555
5555  continue
         v000=V0
         ws=rhoc/(1+dexp(delta/ca))
        rhodens=vrho*(1.d0-(ws/rhoc)**gamma)
      return
      end
            

      double precision function ws(zi,rho)
c ne da valorile potentialului ws functie de z si ro (coordonatele cilindrice)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/prtdta/z0,ro0 
      common/limitesol/x101,x201,y101,y201
      common/dltw/delta
      common/adancimi/vp1,vp2
      common/oadancime/v000
      common/adifz/adiz
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
!
      external delta1,delta2
        z=zi-z0p
      V0=vp1
      ca=adiz
c     print*,'vp1,vp2,adiz',vp1,vp2,adiz
c deosebim mai multe cazuri
      ro=dabs(rho)
      z0=z
      ro0=ro
      a1pa2pr3=a1+a2+2*r3-1.d-6
      if(a1pa2pr3.lt.delt)then ! cazul elipselor complet separate
! in acest caz utilizam planul de separare z=0
         if(z.lt.0.d0)then
! elipsa din partea stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
! pana aici elipsa stanga
          else
! elipsa din partea dreapta (2)
        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
          endif
          endif
! pana aici elipsa dreapta
 
          endif
        delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
      return
! Pana aici elipse complet separate
      endif
      if(r3.ge.500.d0)then ! cazul r3 infinit
! in acest caz trebuie calculate liniile de separare intre cele trei zone
!      alpha=atan((y101-201)/(x201-x101))
!      alpha=pi/2-gamma ! gamma da panta planelor de separate tan(gamma)=cotg(alpha)
! tan(gamma)=(x201-x101)/(y101-y201)
      !calculez primul segment de separare 
        if(y101.eq.y201)then
        if(z.le.x101)goto 4441 ! atunci am prima elipsa
        if(z.ge.x201)goto 4443 ! atunci am a doua elipsa
        goto 4442 ! cazul intermediar
        endif
       tang=(x201-x101)/(y101-y201)
      x0prim=x101-1.d0/tang*y101 ! calculat originea
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
          if(rodif.le.ro)goto 4441 ! atunci am prima elipsa 
      x0prim=x201-1.d0/tang*y201 ! calculat originea pentru a doua separare
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
           if(rodif.ge.ro)goto 4443 ! atunci am a doua elipsa

           goto 4442

4441      continue
! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
! pana aici elipsa din stanga 1

           goto 4445
4443       continue
! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-5,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         
! pana aici elipsa din dreapta

           goto 4445
4442  continue ! a ramas partea intermediara
      z0=z
      ro0=ro
      if(y101.eq.y201)then
      delta=ro0-y101
      else
      alpha=datan((y101-y201)/(x201-x101))
      roc=y101+(z-x101)*(y201-y101)/(x201-x101)!calculez punctul corespunzator pe suprafata 
      delta=(ro0-roc)*dcos(alpha)
      endif
4445  continue
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
c              write(33,*)'crau'
      return
! pana aici cazul r3 infinit
      endif


      if(s.gt.0.d0)then
      if(dabs(ro).gt.ro3)then
      xx02=c3
      goto 6688
      endif
! calculez y11 si y22 punctele de pe suprafata elipsei la trecerea
! dintre regiuni
! folosesc rapoarte pentru determinarea valorilor x functie de ro
! pe suprafetele de separare
! daca z se afla intre cele doua limite atunci am cazul regiunii mediane
                       ! s pozitiv
       y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
       xx01=(x1-c3)*(ro3-dabs(ro))/(ro3-y11)+c3

       y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)
       xx02=(x2-c3)*(ro3-dabs(ro))/(ro3-y22)+c3
           if(r3.lt.1.d-3)then
           xx01=c3
           xx02=c3
           endif
     
         if(z.ge.xx01.and.z.le.xx02)goto 6677 ! regiunea mediana

             else
                        ! s negativ
      y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
      y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)

      xx01=(x1-c3)*(dabs(ro)+s*ro3)/(s*ro3+y11)+c3
      xx02=(x2-c3)*(dabs(ro)+s*ro3)/(s*ro3+y22)+c3
         if(z.le.xx02.and.z.ge.xx01)goto 6677 ! regiunea mediana
             endif
     
      goto 6688 ! mergi intr-una din elipse


6677        continue ! regiunea mediana
      V0=(vp1+vp2)/2.d0
      ! ramane cazul intermediar pentru r3
      ! ro3 are semn negativ pentru s=-1 si pozitiv pentru s=1
                delta=dsqrt((c3-z)**2+(dabs(ro)-ro3)**2)-r3
                if(s.gt.0.d0)delta=-delta
                goto 5555
6688  continue
c      if(z.gt.z3.and.dabs(ro).gt.ro3)goto 5553
      if(z.ge.xx02)goto 5553
!5551  continue      

! elipsa din stanga (1)
      V0=vp1
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c1)then
      dif=z0-c1
      z0=c1-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c1-a1
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c1-z0.lt.0.0001d0)then
      x0=c1
      y0=b1
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c1-a1+0.0001
      xri=c1-0.0001
      CALL DRTMI(X,F,delta1,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b1**2+(c1-z0)**2/a1**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b1*dsqrt(1.d0-(x0-c1)**2/a1**2)
        endif
        endif
          
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
! pana aici elipsa din stanga 1

      goto 5555

5553  continue      

! elipsa din dreapta (2)


        v0=vp2
         z0=z
         ro0=dabs(ro)
      if(z0.gt.c2)then
      dif=z0-c2
      z0=c2-dif
      endif
      if(ro0.lt.0.0001)then
      y0=0
      x0=c2-a2
      semn=-1
      if(z0.lt.x0)semn=1
      else      
      if(c2-z0.lt.0.0001d0)then
      x0=c2
      y0=b2
      semn=-1
      if(ro0.gt.y0)semn=1
      else
      xli=c2-a2+0.0001
      xri=c2-0.0001
      CALL DRTMI(X,F,delta2,XLI,XRI,1.e-6,200,IER)
      semn=1.d0
      ssss=ro0**2/b2**2+(c2-z0)**2/a2**2
      if(ssss.lt.1.d0)semn=-1
      x0=x
      y0=b2*dsqrt(1.d0-(x0-c2)**2/a2**2)
        endif
        endif
                delta=semn*dsqrt((z0-x0)**2+(ro0-y0)**2)
         

      goto 5555
5555  continue
         v000=V0
         ws=-V0/(1+dexp(delta/ca))
      return
      end
            


      double precision function delta1(x) ! determina x
      implicit double precision (a-h,o-z) ! pe suprafata
      common/prtdta/z0,ro0 
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt

      sx=1.d0-(b1/a1)**2-(z0-c1)/(x-c1)
      sx=sx/ro0
      ss=b1/a1**2/dsqrt(1.d0-((x-c1)/a1)**2)
      delta1=sx+ss
      return
      end
            
      double precision function delta2(x) ! determina x 
      implicit double precision (a-h,o-z) ! pe suprafata
      common/prtdta/z0,ro0 
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      sx=1.d0-(b2/a2)**2-(z0-c2)/(x-c2)
      sx=sx/ro0
      ss=b2/a2**2/dsqrt(1.d0-((x-c2)/a2)**2)
      delta2=sx+ss
      return
      end

      subroutine donare(nrnive,en,infin,nrspin,nordin)
      implicit double precision (a-h)
      dimension en(2925),e(2925),infin(2925),iiff(2925),nspin(2925),
     c          nrspin(2925),nordin(2925),norm(2925)     
      do i=1,nrnive
      e(i)=en(i)
      iiff(i)=infin(i)
      nspin(i)=nrspin(i)
      norm(i)=nordin(i)
      enddo
      do i=1,nrnive
      ii=nrnive+1-i
      emin=e(ii)
      iemin=iiff(ii)
      iespin=nspin(ii)
      no=norm(ii)
      do j=1,nrnive-i
      if(emin.gt.e(j))then
      iein=iiff(j)
      ispin=nspin(j)
      nom=norm(j)
      iiff(j)=iemin
      nspin(j)=iespin
      norm(j)=no
      iemin=iein
      iespin=ispin
      no=nom
      ein=e(j)
      e(j)=emin
      emin=ein
      endif
      enddo
      iiff(ii)=iemin
      nspin(ii)=iespin
      norm(ii)=no
      e(ii)=emin
      enddo
      do ij=1,nrnive
      en(ij)=e(ij)
      infin(ij)=iiff(ij)
      nrspin(ij)=nspin(ij)
      nordin(ij)=norm(ij)
      enddo
      return
      end

      subroutine ordonar(nrnive,en)
      implicit double precision (a-h)
      dimension en(2925),e(2925)
      do i=1,nrnive
      e(i)=en(i)
      enddo
      do i=1,nrnive
      ii=nrnive+1-i
      emin=e(ii)
      do j=1,nrnive-i
      if(emin.gt.e(j))then
      ein=e(j)
      e(j)=emin
      emin=ein
      endif
      enddo
      e(ii)=emin
      enddo
      do ij=1,nrnive
      en(ij)=e(ij)
      enddo
      return
      end

      subroutine orniv(niv,e1,e2,indfin)
c creem vector infin cu valoare 1 daca nivelul este localizat
c in 1 groapa si 2 pt a 2 a groapa la fiecare omega
      implicit double precision (a-h,o-z)
      dimension e1(2925),e2(2925),indfin(2925)
      call ordonar(niv,e1)
      call ordonar(niv,e2)
      do i=1,niv
      niv1=niv+1-i
      ee1=e1(niv1)
      do j=1,niv
      niv2=niv+1-j
      ee2=e2(niv2)
      if(ee1.gt.ee2)then
      indfin(niv1)=2
      e2(niv2)=10000000.d0
      do jj=2,niv1
      e1(jj-1)=e1(jj)
      enddo
      e1(niv1)=ee2
      goto 1
      endif
      enddo
      indfin(niv1)=1
1     continue
      enddo
      return
      end

      subroutine finale(a,a2,ckst,ak1,ak2,ae1,ae2,nmax,infin,hw1,hw2)
      implicit double precision (a-h,o-z)
      dimension e1(2925),e2(2925),indfin(2925),infin(2925)
c vector infin cu indice 1 sau 2 pentru localizari in gropi
c acesti indici sunt pentru primul subspatiu omega, apoi pentru urmatorul 
c si asa mai departe, deci avem ordine in fiecare dintre subspatii
      hw1=ckst*(a-a2)**(-1.d0/3.d0)
      hw2=ckst*a2**(-1.d0/3.d0)
      nivvv=0
c pentru un subspatiu omega initializem toate valorile=10000
      do i=1,325
      e1(i)=10000
      e2(i)=10000
      enddo
c luam un subspatiu omega
      do nom=1,nmax
      omega=nom-.5d0
      niv=0
c luam numarul cuantic nro
      do nnn=1,nmax
      n=nnn-1
      np2=n/2
      no2=np2+np2
      lin=0
      if(no2.ne.n)lin=1
c luam numarul cuantic l
      do l=lin,n,2
      aj1=l-.5d0
      aj2=l+.5d0
      aj11=aj1+.01d0
      if(aj11.gt.omega)then
      NIV=NIV+1
      e1(niv)=hw1*((n+1.5d0)-ak1*(aj1*(aj1+1.d0)-l*(l+1.d0)-.75d0)-
     -        ak1*ae1*(l*(l+1)-.5d0*n*(n+3.d0)))
      e2(niv)=hw2*((n+1.5d0)-ak2*(aj1*(aj1+1.d0)-l*(l+1.d0)-.75d0)-
     -        ak2*ae2*(l*(l+1)-.5d0*n*(n+3.d0)))
      endif
      aj22=aj2+.01d0
      if(aj22.gt.omega)then
      niv=niv+1
      e1(niv)=hw1*((n+1.5d0)-ak1*(aj2*(aj2+1.d0)-l*(l+1.d0)-.75d0)-
     -        ak1*ae1*(l*(l+1)-.5d0*n*(n+3.d0)))
      e2(niv)=hw2*((n+1.5d0)-ak2*(aj2*(aj2+1.d0)-l*(l+1.d0)-.75d0)-
     -        ak2*ae2*(l*(l+1)-.5d0*n*(n+3.d0)))
      endif
      enddo
      enddo
      call orniv(niv,e1,e2,indfin)
      do jk=1,niv
      nivk=jk+nivvv
      infin(nivk)=indfin(jk)
      enddo
      nivvv=nivvv+niv
      enddo
      return
      end

ccc ATENTIE MODIFICARI FACUTE VALABILE NUMAI PENTRU ALFA
      subroutine densitateuni(z,ro,efermi0,delaaa0,efermi1,delaaa1,rrsc,
     c                  nmax,densit)
      implicit double precision (a-h,o-z)
      dimension en(2925),infin(2925),nrspin(2925),nordin(2925)
      dimension u(25,325),v(25,325),eq(25,325)
      dimension indice(25,325)
      common/nivselect2/eselect(2925),nnnnn

      common/nrcuantice/numcz1(25,325),numcm1(25,325),numcro1(25,325),
     c          numcsp1(25,325),
     c          numcz2(25,325),numcm2(25,325),numcro2(25,325),
     c          numcsp2(25,325)
      common/vectoripro/VECPRO(25,325,325),LELMAX(25)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      common/numecuantce/X1D(0:24),X2D(0:24),CO1D(0:24),CO2D(0:24)
      common/ckst/ckst
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)
      common/nrnivcmass/g,nrnivc,iier
        common/partnr2/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
       common/enrde/enrgd(2925),nspind(2925),infnpd(2925)

cccc  determin u si v din BCS la fel ca la mase efective
c creez vectorii cu energii, numerele lel si iomg care vor fi
c pusi in ordine
          nrniv=nrnive
      i=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do l=1,nelmax
      i=i+1
      en(i)=emasef(iomg,l)
      infin(i)=infnpd(i)
      nrspin(i)=iomg
      indice(iomg,l)=infnpd(i)
      enddo
      enddo
      efermi=efermi0
      delaaa=delaaa0
      rrsc=100000.
c     if(delt.gt.rrsc)then
c     efermi=efermi1
c     delaaa=delaaa1
c     endif


      call donare(nrniv,en,infin,nrspin,nordin)
      nnnnn=0
      do i=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+i
      nnnnn=nnnnn+1
      eselect(nnnnn)=en(n1)
      enddo
      energiamax=eselect(1)
      energiamin=eselect(nnnnn)
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do j=1,nelmax
      eq(iomg,j)=dsqrt((emasef(iomg,j)-efermi)**2+delaaa**2)
      if(emasef(iomg,j).gt.energiamax)then
      u(iomg,j)=1.d0
      v(iomg,j)=0.d0
      else
      if(emasef(iomg,j).lt.energiamin)then
      u(iomg,j)=0.d0
      v(iomg,j)=1.d0
      else
      u(iomg,j)=dsqrt(0.5d0*(1.d0+(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
      v(iomg,j)=dsqrt(0.5d0*(1.d0-(emasef(iomg,j)-efermi)/
     c    eq(iomg,j)))
c         if(indice(iomg,j).eq.2)then
c         u(iomg,j)=0
c         v(iomg,j)=1
c         endif
            
      endif
      endif
c           print*,'v(iomg,j)',v(iomg,j),iomg,j
      enddo
      enddo


      densit=0
      do iomg=1,nmax
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do ll=1,nelmax
              densit1=0
              densit2=0
      do l=1,nelmax
c calculeaza pla=L_{n}^{k}(x) (polinom laguerre asociat)
cccc constanta de normare laguerre
      fact1=factorial(numcro1(iomg,l))
      fact2=factorial(numcro1(iomg,l)+numcm1(iomg,l))
      cl=dsqrt(2.d0*fact1/fact2)
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1R0=W0*R0/b1
      alf1ro=DSQRT(MPH*W1r0)*1.D-15
      rho=(alf1ro*ro)**2
      n=numcro1(iomg,l)
      k=iabs(numcm1(iomg,l))
      call laguerre(rho,n,k,pla1)
c        print*,'rho,n,k,pla1',rho,n,k,pla1
      rnm=cl*alf1ro*dexp(-rho/2.)*(rho)**(k/2.)*pla1
        v1=dexp(-rho/2.)
        v2=(rho)**(k/2.) 
c      print*,'cl,alf1ro,rho,v1,v2',cl,alf1ro,rho,v1,v2
c         print*,'rnm,iomg,l,ll',rnm,iomg,l,ll
cccc impart la radical din pi pentru funtia unghiulara
      rnm=rnm/dsqrt(2.d0*3.14)
cccc calculez partea dupa z
          if(z.lt.0)then
          zz=alf1*z+dabs(zet1)
          anr1=x1d(numcz1(iomg,l))
          call dherm(anr1,zz,h1,d1)
c              print*,'anr1,zz,h1,d1',anr1,zz,h1,d1
          fz=co1d(numcz1(iomg,l))*dexp(-zz**2/2.)*h1
          else
          zz=-(z*alf2-zet2)
          anr2=x2d(numcz1(iomg,l))
          call dherm(anr2,zz,h2,d2)
c              print*,'anr2,zz,h2,d2',anr1,zz,h1,d1
          fz=co2d(numcz1(iomg,l))*dexp(-zz**2/2.)*h2
          endif
       rnmf=v(iomg,ll)*rnm*fz*vecpro(iomg,l,ll)
c         print*,'v(iomg,ll)*rnm*fz*vecpro(iomg,l,ll)',v(iomg,ll),rnm,fz
c    c ,vecpro(iomg,l,ll)
c         print*,'rnmf,iomg,l,ll',rnmf,iomg,l,ll
                 if(numcsp1(iomg,l).eq.-1)then
                 densit1=densit1+rnmf
                 else
                 densit2=densit2+rnmf
                 endif
c         print*,'densit',densit
      enddo
       densit=densit+densit1**2+densit2**2
c             stop
      enddo
      enddo

      return
      end


      double precision function factorial(numar)
      implicit double precision (a-h,o-z)
      if(numar.eq.0)then
      factorial=1.d0
      return
      endif
      factorial=1.d0
      do i=1,numar
      factorial=factorial*i
      enddo
      return
      end
      
      


      subroutine stangadreapta(aniu1,aniu1p,aniu2,aniu2p,
     c const1,const1p,const2,const2p,stang,dreapt)
c calculam separat pe partea dreapta si partea stanga
c c^2*integrala functia de unda la patrat
c ca sa aflam unde este localizata functia de unda
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      ZETA1=ZET1
      ZETA2=ZET2
      call jcptn(aj1,aniu1,aniu1p,zeta1)
      stang=const1*const1p*aj1/alf1
      call jcptn(aj2,aniu2,aniu2p,zeta2)
      dreapt=const2*const2p*aj2/alf2
c      print*,aniu1,aniu2,const1,const2,stang,dreapt,aj1,aj2
      return
      end

      SUBROUTINE NENGLScumef(A0I,A2E,R1F,R2F,V00,V01,V02,
     C alp0,raso0,ra0,
     C NMAX,IZ,EN,X1D,X2D,CO1D,CO2D,npr,nfit,IAS,nrspin,nnlel,nfinal)

C IZ=0 protoni, IZ=1 neutroni


C       NMAX NR. NIVELE MAXIM CONSIDERATE < 14
C CELELALTE VALORI DE INTRARE SINT INTRODUSE PRIN COMMON/FORV/FORV(18)
C OUTPUT EN NIVELE ENERGETICE OBTINUTE DIMENSIUNE EN(168)
C        X1D, X2D NR. NIU PE AXA Z PENTRU 1 SI 2 DIMENSIUNE 7
C        CO1D, CO2D CONST. NORMARE FUNCTII DE UNDA DIMENSIUNE 7
C      PARAMETER (N1=14, N2=N1-1, N3=N1*(N1+1)*(N1+2)/6, N4=N6+(N6*N6-N6)/2,
C                 N6=N1*(N+1)/2,  N5=2*N1-1, N7=N6*N6)
C                 N1=14 in acest program

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      REAL ECC,R
      dimension dreapta(25,325,325),stanga(25,325,325)
      dimension ecc(10000000),r(10000000)
      common/nrcuantice/numcz1(25,325),numcm1(25,325),numcro1(25,325),
     c          numcsp1(25,325),
     c          numcz2(25,325),numcm2(25,325),numcro2(25,325),
     c          numcsp2(25,325)
      dimension aemeps1(25,325,325),aemeps2(25,325,325),
     c         aema1pa2(25,325,325),aemcr3(25,325,325),
     c         aemdelt(25,325,325)
      dimension nparitate(0:24)      
      dimension cnorm1(325),reanc1(325),cnorm2(325),reanc2(325)

       dimension X1D(0:24),X2D(0:24),EN(2925),ENI(25,25),
     C          CO1D(0:24),CO2D(0:24),EI(25,25),!ECC(36100),
     c          nrspin(2925),nfinal(2925)
     C         ,NR1(0:49),NR2(0:49)
      common/numecuantce/X1Dc(0:24),X2Dc(0:24),CO1Dc(0:24),CO2Dc(0:24)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/vectoripro/VECPRO(25,325,325),LELMAX(25)

      common/limitesol/x101,x201,y101,y201

      common/elmdef/elmateps1(25,325,325),elmateps2(25,325,325),
     c         elmata1pa2(25,325,325),elmatcr3(25,325,325),
     c         elmatdelt(25,325,325),emasef(25,325)

      common/elmdeflim/elmateps1lim(25,325,325),elmateps2lim(25,325,325)
     c         ,elmata1pa2lim(25,325,325),elmatcr3lim(25,325,325),
     c         elmatdeltlim(25,325,325),emaseflim(25,325)
      common/elim/energlim(25,325,325)

      common/variat/variat
c      common/eig/r(5565)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/W12/W1,W2
      COMMON/ZET1S2/Z1,Z2,ZN1,ZN2
      common/nrfin/nrfin1,nrfin2
      common/ptrads/detai,xdr1i,xdr2i,ndri
c urmatorul commmon ddimmensionat pt nmmax=25
      common/nrint/k1,k2
      common/r116/R0
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/optest/optest
      common/adancimi/vpo1,vpo2
      common/CKST/CKST
      common/nznnll/nznn,nzll
      common/zilim12/z1lim,z2lim
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
      common/adancimea/vadanc
      common/cx12/cx1,cx2
      common/lasciziune/delatsciz
      vadanc=v00

! pentru inceput redefinesc coordonatele functie de originea potentialului cu 
! doua centre
      CT=a2/dabs(a1)
      cx1=-DELT/(1.D0+CT)
      cx2=CT*dabs(cx1)
C DISTANTA DINTRE PLANUL DE INTERSECTIE A PARAMETRIZARII SEMISIMETRICE  SI
C PLANUL DE INTERSECTIE A PARAMTERIZARII ASIMETRICE
      Z0P=cx1-c1
c           print*,'z0p',z0p



      ndri=0
      npr=0
      ALP=ALP0
      nrfin1=0
      nrfin2=0
c calculez HW pentru elipsoizi, trebuie sa am semi-axele initiale, la suprafata trebuie
c sa am mereu CKST=41 MeV
c     CKST=41.D0
      HWRO1=CKST*R0/b1
      HWZ1=b1/a1*HWRO1
      HWRO2=CKST*R0/b2
      HWZ2=b2/a2*HWRO2

      DELTA=DELT  
      R1FI=A1 
      R2FI=A2 
      NMAX0=NMAX
      CALL ORDO(R1FI,R2FI,NMAX0,NR1,NR2) !da nr1(0:27) si nr2(0:27)
C CALCULAM NR DE STARI CU ENERGII DIFERITE NMAXX PENTRU NR NIVELE NMAX 
      NMAXX=NMAX*(NMAX+1)*(NMAX+2)/6

      R0C=R0*A0I**(1.D0/3.D0)
      Z1=DABS(Cx1)
      Z2=Cx2
      CALL P2C(R0C,A1,A2,Z1,Z2)
                  ! calculeaza zet1,zet2,alf1,alf2,w1,w2
C parametrizez adancimile potentialelor, constantele
C spin orbita nu trebuie parametrizate deoarece sunt aceleasi
C pentru orice A si Z
c      DPAR=R0*((A0I-A2E)**(1./3.)+A2E**(1./3.))+2.d0*R3
      DPAR=delatsciz
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      IF(DELTA.GT.DPAR)THEN
      VP1=V01
      VP2=V02
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      VP1=V00
      VP2=V00
      ELSE
      VP1=V00+(DELTA-DPAR07)/(DPAR-DPAR07)*(V01-V00)
      VP2=V00+(DELTA-DPAR07)/(DPAR-DPAR07)*(V02-V00)
      ENDIF
      ENDIF
      vpo1=vp1
      vpo2=vp2

      DO 2 I=1,NMAX
      L=I-1
      NRFIN1=NR1(L)
      NRFIN2=NR2(L)
      xdr1i=x1d(l)
      xdr2i=x2d(l)
      R11=A1
      R22=A2
      CE1=DABS(Cx1)
      CE2=Cx2
      CALL NRADAC(L,R11,R22,DELTA,CE1,CE2,X1,X2)
      CALL KNORM(X1,X2,CO1,CO2)
      nparitate(l)=1
      xrapx=R22/R11
      if(xrapx.lt.1.04d0.and.xrapx.gt.0.96d0)then
      CALL NRADAC1(L,R11,DELTA,X1)
      x2=x1
      CALL KNORM1(X1,X2,CO1,CO2)
      nparitate(l)=(-1)**(l)
      endif
c     print*,'x1x2co1co2',x1,x2,co1,co2
      X1D(L)=X1
      X2D(L)=X2
      CO1D(L)=CO1
      CO2D(L)=CO2
2     CONTINUE
            do lc=0,nmax-1
            x1dc(lc)=x1d(lc)
            x2dc(lc)=x2d(lc)
            co1dc(lc)=co1d(lc)
            co2dc(lc)=co2d(lc)
            enddo
c     write(61,*)x1d(0),x1d(1),x1d(2),x1d(3),x1d(4),x1d(5),x1d(6),
c    c  x1d(7),x1d(8),x1d(9),x1d(10),x1d(11),x1d(12),x1d(13)
c     write(62,*)x2d(0),x2d(1),x2d(2),x2d(3),x2d(4),x2d(5),x2d(6),
c    c  x2d(7),x2d(8),x2d(9),x2d(10),x2d(11),x2d(12),x2d(13)
c     write(63,*)co1d(0),co1d(1),co1d(2),co1d(3),
c    c   co1d(4),co1d(5),co1d(6),
c    c  co1d(7),co1d(8),co1d(9),co1d(10),co1d(11),co1d(12),co1d(13)
c     write(64,*)co2d(0),co2d(1),co2d(2),co2d(3),
c    c   co2d(4),co2d(5),co2d(6),
c    c  co2d(7),co2d(8),co2d(9),co2d(10),co2d(11),co2d(12),co2d(13)

      z1lim=cx1
      z2lim=cx2
                 nmaxm1=nmax-1
      xl=alf1*z1lim
      xu=0.d0
      xl1=xl
        call herval1(nmaxm1,x1d,xl,xu) 
      xu=alf2*z2lim
      xl=0.d0
        call herval2(nmaxm1,x2d,xl,xu) 
      call fherm(nmaxm1,x1d,x2d)    
      call vpot(xl1,xu)
      call dpot5(xl1,xu)

       do mmm1=0,nmaxm1
       mxn1=(nmaxm1-mmm1)/2
       do mmm2=0,nmaxm1
       mxn2=(nmaxm1-mmm2)/2
       do nnn1=0,mxn1
       do nnn2=0,mxn2
       call pertwsi(nnn1,mmm1,nnn2,mmm2)
       call lzszi(nnn1,mmm1,nnn2,mmm2)
       if(iz.eq.0)call coulini(nnn1,mmm1,nnn2,mmm2)
       call ameps1i(nnn1,mmm1,nnn2,mmm2)
       call ameps2i(nnn1,mmm1,nnn2,mmm2)
       call ama1pa2i(nnn1,mmm1,nnn2,mmm2)
       call amcr3i(nnn1,mmm1,nnn2,mmm2)
       call amdelti(nnn1,mmm1,nnn2,mmm2)
       enddo
       enddo
       enddo
       enddo


      NNEL=0
      LLEL=0
C CALCULAM DIFERITE VALORI PENTRU OMEGA 0.5,1.5,2.5,3.5,...(NMAX-1)/2
      DO 9876 IOMG=1,NMAX!2!NMAX
      OMG=IOMG-.5d0
c          omg=-omg


      LEL=0
      DO 7 J=1,NMAX
      DO 4 JJ=1,J
C INDEXAREA NUMARULUI CUANTIC NZ1 CU VALOAREA PENTRU SFERA
      NN=JJ-1
C NUMARUL N AL NIVELULUI
      N=J-1
C NUMARUL  M=NFI  AL NIVELULUI
      DO 4 MM1=-N,N
      M1=-MM1
C NUMARUL RO AL NIVELULUI
      NRO1=N-NN-IABS(M1)
      NRO1=NRO1/2
      IF(NRO1.LT.0) GO TO 444
      IF(NRO1+NRO1+IABS(M1)+NN.NE.N) GO TO 444
C INDEXARE DUPA SPINUL 1 NS1
      DO 4444 NS1S=-1,1,2
      SPIN1=NS1S/2.D0
      omeg1=dabs(m1+spin1)
      OMEG11=M1+SPIN1
      IF(OMEG11.NE.OMG)GOTO 4444
C INDEXAREA MATRICII DUPA PRIMUL INDICE
      LEL=LEL+1

       numcz1(iomg,lel)=nn
       numcm1(iomg,lel)=m1
       numcro1(iomg,lel)=nro1
       numcsp1(iomg,lel)=ns1s     




C
C SE INDEXEAZA LA FEL SI NUMERELE CUANTICE PENTRU AL 2-LEA INDICE
C AL MATRICII
      NEL=0
      DO 3 K=1,NMAX
      DO 3 KK=1,K
      LL=KK-1
      L=K-1
      DO 3 MM2=-L,L
      M2=-MM2
      NRO2=L-LL-IABS(M2)
      NRO2=NRO2/2
      IF(NRO2.LT.0)GOTO333
      IF(NRO2+NRO2+IABS(M2)+LL.NE.L) GO TO 333
      DO 3333 NS2S=-1,1,2
      SPIN2=NS2S/2.D0
C INDEXAREA MATRICII DUPA AL DOILEA INDICE
      omeg2=dabs(m2+spin2)
      OMEG22=M2+SPIN2
      IF(OMEG22.NE.OMG)GOTO 3333
      NEL=NEL+1

       numcz2(iomg,nel)=ll
       numcm2(iomg,nel)=m2
       numcro2(iomg,nel)=nro2
       numcsp2(iomg,nel)=ns2s     


C
           nzll=ll
           nznn=nn
      k1=nn
      k2=ll
      prs=1.d0
c matricea fiind simetrica, calculam numai partea triunghiulara
      IF(NEL.LT.LEL) GO TO 3333


       
      aemeps1(iomg,nel,lel)=0.d0
      aemeps1(iomg,lel,nel)=0.d0      
      aemeps2(iomg,nel,lel)=0.d0
      aemeps2(iomg,lel,nel)=0.d0       
      aema1pa2(iomg,nel,lel)=0.d0
      aema1pa2(iomg,lel,nel)=0.d0       
      aemcr3(iomg,nel,lel)=0.d0
      aemcr3(iomg,lel,nel)=0.d0       
      aemdelt(iomg,nel,lel)=0.d0
      aemdelt(iomg,lel,nel)=0.d0


      stanga(iomg,nel,lel)=0.d0
      stanga(iomg,lel,nel)=0.d0
      dreapta(iomg,nel,lel)=0.d0
      dreapta(iomg,lel,nel)=0.d0



c     print*,'nel,lel,omg',nel,lel,omg

      enlel=0.D0
      enlsim=0.d0
      enel=0.d0
      enell=0.d0
c calculez partea simetrica
      if(spin1.eq.spin2)then
         if(m1.eq.m2.and.nro1.eq.nro2.and.nn.eq.ll)then

      enlel=hwz1*(x1d(nn)+0.5)+
     + hwro1*(2*nro1+iabs(m1)+1.d0)-vadanc
c                print*,'!!!!!!!!!!!!!!! vadanc',vadanc

      enlsim=enlel 
         endif
C CALCULAM PERTURBATIILE
c perturbatia datorita potentialului wood saxon este pentru m1=m2
         if(m1.eq.m2)then
      aniu1=x1d(nn)
      aniu2=x2d(nn)
      const1=co1d(nn)
      const2=co2d(nn)
      aniu1p=x1d(ll)
      aniu2p=x2d(ll)
      const1p=co1d(ll)
      const2p=co2d(ll)
! calculez perturbatia wood-saxon
           nzll=ll
           nznn=nn

      call pertws(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,pert)

c      print*,'am iesit cu pert(n1,n2,m1,m2)',pert,nro1,nro2,m1,m2
c      write(25,*)nn,ll,nro1,nro2,m1,m2,pert
c       print*,'pert',pert
c       stop
         enlel=enlel+pert
         enlsim=enlsim+pert
         pertlz=0
         call lzsz(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,pertlz)
!c     alp !constanta spin orbita si semn minus in fata
       pertlz=-alp*spin1*m1*pertlz*(197.327)**2/2.d0/(931.49432d0)**2 ! hbar*c=197.327 MeV fm
c      write(25,*)nn,ll,nro1,nro2,m1,m2,pertlz
         enlel=enlel+pertlz
         enlsim=enlsim+pertlz
c          print*,'pertlz',pertlz
           percou=0
              if(iz.eq.0)then
              call coulin(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,percou)
c          if(dabs(percou).ge.1.d0)print*,percou,nn,ll,nro1,nro2,m1,m2
              enlel=enlel+percou
              enlsim=enlsim+percou
              endif               
       

!!!!!!!!!! aici trebuie introduse elementele de matrice ale maselor efective
      call ameps1(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,emeps1)
      call ameps2(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,emeps2)
      call ama1pa2(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,ema1pa2)
      call amcr3(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,emcr3)
      call amdelt(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,emdelt)
c      print*,'emeps1,emeps2,ema1pa2,emcr3,emdelt'
c      print*,emeps1,emeps2,ema1pa2,emcr3,emdelt

       
      aemeps1(iomg,nel,lel)=emeps1
      aemeps1(iomg,lel,nel)=emeps1       
      aemeps2(iomg,nel,lel)=emeps2
      aemeps2(iomg,lel,nel)=emeps2       
      aema1pa2(iomg,nel,lel)=ema1pa2
      aema1pa2(iomg,lel,nel)=ema1pa2       
      aemcr3(iomg,nel,lel)=emcr3
      aemcr3(iomg,lel,nel)=emcr3       
      aemdelt(iomg,nel,lel)=emdelt
      aemdelt(iomg,lel,nel)=emdelt







        endif


      endif

      spin1p=spin1-1.d0
      spin2p=spin1+1.d0
      m1plus=m1+1
      m1min=m1-1

      if(spin1p.eq.spin2.and.m1plus.eq.m2)then
      aniu1=x1d(nn)
      aniu2=x2d(nn)
      const1=co1d(nn)
      const2=co2d(nn)
      aniu1p=x1d(ll)
      aniu2p=x2d(ll)
      const1p=co1d(ll)
      const2p=co2d(ll)
c     print*,'intru in lplus cu spin1p,spin2,m1plus,m2'
c     print*,spin1p,spin2,m1plus,m2
      pertlp=0
      call lplussmin(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,pertlp)
c      print*,'am iesit cu pertlp(n1,n2,m1,m2)',pertlp,nro1,nro2,m1,m2
!c     alp !constanta spin orbita
c      write(25,*)nn,ll,nro1,nro2,m1,m2,pertlp
      pertlp=-alp*0.5*pertlp*(197.327)**2/2.d0/(931.49432d0)**2 ! hbar*c=197.327 MeV fm
         enlel=enlel+pertlp
         enlsim=enlsim+pertlp





      endif


      if(spin2p.eq.spin2.and.m1min.eq.m2)then
      aniu1=x1d(nn)
      aniu2=x2d(nn)
      const1=co1d(nn)
      const2=co2d(nn)
      aniu1p=x1d(ll)
      aniu2p=x2d(ll)
      const1p=co1d(ll)
      const2p=co2d(ll)
!c     print*,'intru in lmin cu spin1p,spin2,m1plus,m2'
!c     print*,spin1p,spin2,m1plus,m2
      pertlm=0
      call lminplus(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro2,m2,pertlm)
c      print*,'am iesit cu pertlm(n1,n2,m1,m2)',pertlm,nro1,nro2,m1,m2
!c     alp !constanta spin orbita
c      write(25,*)nn,ll,nro1,nro2,m1,m2,pertlm
      pertlm=-alp*0.5*pertlm*(197.327)**2/2.d0/(931.49432d0)**2 ! hbar*c=197.327 MeV fm
         enlel=enlel+pertlm
         enlsim=enlsim+pertlm



      endif

      ENLEL=(ENLSIM+ENLEL)/2.D0


c calculez patratul functiei de unda pe partea 
c stanga si partea dreapta ca sa vad unde este localizata
      stang=0
      dreapt=0
      call stangadreapta(aniu1,aniu1p,aniu2,aniu2p,
     c const1,const1p,const2,const2p,stang,dreapt)
      stanga(iomg,nel,lel)=stang
      stanga(iomg,lel,nel)=stang
      dreapta(iomg,nel,lel)=dreapt
      dreapta(iomg,lel,nel)=dreapt
         sumsd=stang+dreapt
c       print*,'sumsd iomg nel,lel',sumsd,iomg,nel,lel,stang,dreapt



c     if(xrapx.lt.1.04d0.and.xrapx.gt.0.96d0)then
c       if(xrapx.ge.1.d0)enlel=enlel+
c    c               (energlim(iomg,nel,lel)-enlel)*(xrapx-1.)/(1.04-1.)
c       if(xrapx.lt.1.d0)enlel=enlel+
c    c               (energlim(iomg,nel,lel)-enlel)*(xrapx-1.)/(0.96-1.)
c     endif
        
        
      KEI=LEL+(NEL*NEL-NEL)/2
      R(KEI)=ENLEL
3333  continue
333   continue
3     CONTINUE
4444  continue
444   continue
4     CONTINUE
7     continue
       DO INDX=1,LEL
       DO INDY=1,LEL
       VECPRO(IOMG,INDX,INDY)=0.
       ENDDO
       ENDDO
       LELLE=LEL**2
       DO IW=1,LELLE
       ECC(IW)=0.
       ENDDO
       nnlel=lel
      LELLE=LEL**2
      DO IW=1,LELLE
      ECC(IW)=0.
      ENDDO
      CALL EIGEN(R,ECC,LEL,0)
      DO 130 IOI=1,LEL
      LOI=IOI+(IOI*IOI-IOI)/2
      NRSPIN(IOI+LLEL)=IOMG
      emasef(iomg,ioi)=r(loi)
c     if(xrapx.lt.1.04d0.and.xrapx.gt.0.96d0)then
c       if(xrapx.ge.1.d0)emasef(iomg,ioi)=emasef(iomg,ioi)+
c    c     (emaseflim(iomg,ioi)-emasef(iomg,ioi))*(xrapx-1.)/(1.04-1.)
c       if(xrapx.lt.1.d0)emasef(iomg,ioi)=emasef(iomg,ioi)+
c    c     (emaseflim(iomg,ioi)-emasef(iomg,ioi))*(xrapx-1.)/(0.96-1.)
c     endif
130   EN(IOI+LLEL)=R(LOI)
      IW=0

      DO INDY=1,LEL
      DO INDX=1,LEL
      IW=IW+1
      VECPRO(IOMG,INDX,INDY)=ECC(IW)
c      vecpro(iomg,indy,indx)=ecc(iw)
      ENDDO
      ENDDO
      LELMAX(IOMG)=LEL
C TEST PENTRU VERIFICAREA VALORILOR PROPRII si pentru determinarea
c localizarii functiei de unda
      DO INDX=1,LEL
      SUMAVE=0.
       sdre=0
       ssta=0
            do indx1=1,lel
            do indy=1,lel
       sdre=sdre+dreapta(iomg,indy,indx1)*VECPRO(IOMG,INDy,INDx)*
     c                                   vecpro(iomg,indx1,indx)
       ssta=ssta+stanga(iomg,indy,indx1)*VECPRO(IOMG,INDy,INDx)*
     c                                   vecpro(iomg,indx1,indx)
            enddo
            enddo
      DO INDY=1,LEL
       SUMAVE=SUMAVE+VECPRO(IOMG,INDY,INDX)**2
c      print*,'dr st',dreapta(iomg,indy,indx),stanga(iomg,indy,indx)
      ENDDO
        nfinal(INDX+LLEL)=1
        if(sdre.gt.ssta)nfinal(indx+llel)=2
             nnfinal=nfinal(indx+llel)
C VALOAREA LUI SUMAVE TREBUIE SA FIE 1
c      PRINT*,'SUMAVE IOMG', SUMAVE,IOMG
c         sumave2=sdre+ssta
c       print*,'nnfinal,indx,sumave2',nnfinal,indx,sumave2    
      ENDDO
                               
C PANA AICI NEL=LEL=(Nmax+1-IOMG)*(Nmax+2-IOMG)/2
      NNEL=NNEL+NEL
         nnlel=nnel
      LLEL=LLEL+LEL

!! aici facem elementele de matrice intre doua functii lel,nel
      do i1=1,nel
      do i2=1,lel
c            if(lel.gt.nel)goto 7373
          sumeps1=0
          sumeps2=0
          suma1pa2=0
          sumcr3=0
          sumdelt=0
      do j1=1,nel
      do j2=1,lel
       sumeps1=sumeps1+vecpro(iomg,j1,i1)*vecpro(iomg,j2,i2)*
     c aemeps1(iomg,j1,j2)
       sumeps2=sumeps2+vecpro(iomg,j1,i1)*vecpro(iomg,j2,i2)*
     c aemeps2(iomg,j1,j2)
       suma1pa2=suma1pa2+vecpro(iomg,j1,i1)*vecpro(iomg,j2,i2)*
     c aema1pa2(iomg,j1,j2)
       sumcr3=sumcr3+vecpro(iomg,j1,i1)*vecpro(iomg,j2,i2)*
     c aemcr3(iomg,j1,j2)
       sumdelt=sumdelt+vecpro(iomg,j1,i1)*vecpro(iomg,j2,i2)*
     c aemdelt(iomg,j1,j2)
      enddo
      enddo
       elmateps1(iomg,i1,i2)=sumeps1
       elmateps2(iomg,i1,i2)=sumeps2
       elmata1pa2(iomg,i1,i2)=suma1pa2
       elmatcr3(iomg,i1,i2)=sumcr3
       elmatdelt(iomg,i1,i2)=sumdelt
      enddo
      enddo
             
9876  CONTINUE
                                          
      ndri=npr
      RETURN
      END






      double precision function delta1x(d) ! determina distanta delta
      implicit double precision (a-h,o-z) ! normala pe suprafata pana la un punct z,ro
      common/prtdta/z0,ro0 
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      delta1x=ro0**2*(a1+d)**2+(z0-c1)**2*(b1+d)**2-(a1+d)**2*(b1+d)**2
      return
      end
            
      double precision function delta2x(d) ! determina distanta delta
      implicit double precision (a-h,o-z) ! normala pe suprafata pana la un punct z,ro
      common/prtdta/z0,ro0 
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      delta2x=ro0**2*(a2+d)**2+(z0-c2)**2*(b2+d)**2-(a2+d)**2*(b2+d)**2
      return
      end


       subroutine dervV(z,rho,derddz,derdro)
c calculez derivatele potentialului WS fata de z si de ro
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/prtdta/z0,ro0
      common/limitesol/x101,x201,y101,y201
      common/dltw/delta
      common/adancimi/vp1,vp2
      common/adifz/adiz
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
      external delta1x,delta2x
      z=z+z0p
      if(z.lt.0.d0)then
      V0=vp1
      ca=adiz
      delta=dsqrt((z-c1)**2+rho**2)
      derddz=(z-c1)/delta
      derdro=rho/delta
      else
      v0=vp2
      ca=adiz
      delta=dsqrt((z-c2)**2+rho**2)
      derddz=(z-c2)/delta
      derdro=rho/delta
      endif
      derddz=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derddz
      derdro=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derdro
      return
      end
     
      subroutine dervVx(zi,rho,derddz,derdro)
c calculez derivatele potentialului WS fata de z si de ro
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/prtdta/z0,ro0
      common/limitesol/x101,x201,y101,y201
      common/dltw/delta
      common/adancimi/vp1,vp2
      common/adifz/adiz
      common/zop/z0p !valoarea care se adauga la z pentru potential ws (z0p este negativ)
      external delta1x,delta2x
      z=zi-z0p
      V0=vp1
      ca=adiz
c deosebim mai multe cazuri
      ro=dabs(rho)
      z0=z
      ro0=ro
      a1pa2pr3=a1+a2+2*r3-1.d-6
      if(a1pa2pr3.lt.delt)then ! cazul elipselor complet separate
! in acest caz utilizam planul de separare z=0
         if(z.lt.0.d0)then
      V0=vp1
         z0=z
         ro0=ro
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0 !daca este 0 se cauta si solutii in interior elipsa
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c1)
      if(xxx.gt.a1)nmaxio=1
      if(xxx.lt.a1)then
      yyy=b1*dsqrt(1.d0-xxx**2/a1**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta1x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 7277
      xri1=1.d0*i
      d2xr1=delta1x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 9161
      endif
      xli1=xri1
      d2xi1=d2xr1
7277  continue
       if(nmaxio.eq.0)then ! distanta delta nu poate fi mai mare decat b1
      xri2=-1.d0*i
      alimin=b1+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b1
      endif
      d2xr2=delta1x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 9161
      endif
      xli2=xri2
      d2xi2=d2xr2
       endif
      enddo
9161    continue
      CALL DRTMI(X,F,delta1x,XLI,XRI,1.e-4,100,IER)
      if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER1=',IER
      delta=x
! daca a1 sau b1 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a1.gt.dddd1.and.a1.le.dddd.and.a1.eq.b1)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a1.gt.dddd1.and.a1.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b1+delta)**2/(ro**2/(b1+delta)**3)
      derdro=0.d0
      else
      if(b1.gt.dddd1.and.b1.le.dddd)then
      derdro=1.d0
!      derddz=(z-c1)/(a1+delta)**2/((z-c1)**2/(a1+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c1)/(a1+delta)**2
      derddz=derddz/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      derdro=ro/(b1+delta)**2
      derdro=derdro/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      endif
      endif
      endif
         else
      V0=vp2
         z0=z
         ro0=ro
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c2)
      if(xxx.gt.a2)nmaxio=1
      if(xxx.lt.a2)then
      yyy=b2*dsqrt(1.d0-xxx**2/a2**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta2x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 1178
      xri1=1.d0*i
      d2xr1=delta2x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 6263
      endif
      xli1=xri1
      d2xi1=d2xr1
1178  continue
       if(nmaxio.eq.0)then
      xri2=-1.d0*i
      alimin=b2+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b2
      endif
      d2xr2=delta2x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 6263
      endif
      xli2=xri2
      d2xi2=d2xr2
        endif
      enddo
6263  continue
      CALL DRTMI(X,F,delta2x,XLI,XRI,1.e-4,100,IER)
      if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER2=',IER
         delta=x
! daca a2 sau b2 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a2.gt.dddd1.and.a2.le.dddd.and.a2.eq.b2)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a2.gt.dddd1.and.a2.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b2+delta)**2/(ro**2/(b2+delta)**3)
      derdro=0.d0
      else
      if(b2.gt.dddd1.and.b2.le.dddd)then
      derdro=1.d0
!      derddz=(z-c2)/(a2+delta)**2/((z-c2)**2/(a2+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c2)/(a2+delta)**2
      derddz=derddz/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      derdro=ro/(b2+delta)**2
      derdro=derdro/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      endif
      endif
      endif
         endif
      derddz=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derddz
      derdro=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derdro
      return
      endif
            !!!!!!!!!!!!!!  cazul r3 infinit
      if(r3.ge.500.d0)then ! cazul r3 infinit
! in acest caz trebuie calculate liniile de separare intre cele trei zone
!      alpha=atan((y101-201)/(x201-x101))
!      alpha=pi/2-gamma ! gamma da panta planelor de separate tan(gamma)=cotg(alpha)
! tan(gamma)=(x201-x101)/(y101-y201)
      !calculez primul segment de separare
        if(y101.eq.y201)then
        if(z.le.x101)goto 4441 ! atunci am prima elipsa
        if(z.ge.x201)goto 4443 ! atunci am a doua elipsa
        goto 4442 ! cazul intermediar
        endif
       tang=(x201-x101)/(y101-y201)
      x0prim=x101-1.d0/tang*y101 ! calculat originea
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
          if(rodif.le.ro)goto 4441 ! atunci am prima elipsa
      x0prim=x201-1.d0/tang*y201 ! calculat originea pentru a doua separare
      zdif=(z-x0prim)
      rodif=zdif*tang !calculez panta de separare
           if(rodif.ge.ro)goto 4443 ! atunci am a doua elipsa

           goto 4442

4441      continue
          z0=z
          ro0=ro
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0 !daca este 0 se cauta si solutii in interior elipsa
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c1)
      if(xxx.gt.a1)nmaxio=1
      if(xxx.lt.a1)then
      yyy=b1*dsqrt(1.d0-xxx**2/a1**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta1x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 7117
      xri1=1.d0*i
      d2xr1=delta1x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 6111
      endif
      xli1=xri1
      d2xi1=d2xr1
7117  continue
       if(nmaxio.eq.0)then ! distanta delta nu poate fi mai mare decat b1
      xri2=-1.d0*i
      alimin=b1+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b1
      endif
      d2xr2=delta1x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 6111
      endif
      xli2=xri2
      d2xi2=d2xr2 
       endif
      enddo
6111    continue
      CALL DRTMI(X,F,delta1x,XLI,XRI,1.e-4,100,IER)
           if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER3=',IER
           delta=x
! daca a1 sau b1 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a1.gt.dddd1.and.a1.le.dddd.and.a1.eq.b1)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a1.gt.dddd1.and.a1.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b1+delta)**2/(ro**2/(b1+delta)**3)
      derdro=0.d0
      else
      if(b1.gt.dddd1.and.b1.le.dddd)then
      derdro=1.d0
!      derddz=(z-c1)/(a1+delta)**2/((z-c1)**2/(a1+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c1)/(a1+delta)**2
      derddz=derddz/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      derdro=ro/(b1+delta)**2
      derdro=derdro/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      endif
      endif
      endif
           goto 4445
4443       continue
           z0=z
           ro0=ro
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c2)
      if(xxx.gt.a2)nmaxio=1
      if(xxx.lt.a2)then
      yyy=b2*dsqrt(1.d0-xxx**2/a2**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta2x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 7118
      xri1=1.d0*i
      d2xr1=delta2x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 6222
      endif
      xli1=xri1
      d2xi1=d2xr1
7118  continue
       if(nmaxio.eq.0)then
      xri2=-1.d0*i
      alimin=b2+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b2
      endif
      d2xr2=delta2x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 6222
      endif
      xli2=xri2
      d2xi2=d2xr2
        endif
      enddo
6222  continue
      CALL DRTMI(X,F,delta2x,XLI,XRI,1.e-4,100,IER)
           if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER4=',IER
           delta=x
! daca a2 sau b2 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a2.gt.dddd1.and.a2.le.dddd.and.a2.eq.b2)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a2.gt.dddd1.and.a2.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b2+delta)**2/(ro**2/(b2+delta)**3)
      derdro=0.d0
      else
      if(b2.gt.dddd1.and.b2.le.dddd)then
      derdro=1.d0
!      derddz=(z-c2)/(a2+delta)**2/((z-c2)**2/(a2+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c2)/(a2+delta)**2
      derddz=derddz/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      derdro=ro/(b2+delta)**2
      derdro=derdro/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      endif
      endif
      endif
           goto 4445
4442  continue ! a ramas partea intermediara
      z0=z
      ro0=ro
      if(y101.eq.y201)then
      delta=ro0-y101
      derdro=1.
      derddz=0.d0
      else
      alpha=datan((y101-y201)/(x201-x101))
      roc=y101+(z-x101)*(y201-y101)/(x201-x101)!calculez punctul corespunzator pe suprafata
      delta=(ro0-roc)*dcos(alpha)
      derddz=-(y201-y101)/(x201-x101)*dcos(alpha)
      derdro=dcos(alpha)
      endif
4445  continue
      derddz=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derddz
      derdro=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derdro
      return
      endif
                !!!!!!!!!! ultimul caz cu r3 diferit de infinit
      if(s.gt.0.d0)then
      if(dabs(ro).gt.ro3)then
      xx02=c3
      goto 6688
      endif
! calculez y11 si y22 punctele de pe suprafata la trecerea
! dintre regiuni
! folosesc rapoarte pentru determinarea valorilor x functie de ro
! pe suprafetele de separare
! daca z se afla intre cele doua limite atunci am cazul regiunii mediane
                       ! s pozitiv
       y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
       xx01=(x1-c3)*(ro3-dabs(ro))/(ro3-y11)+c3

       y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)
       xx02=(x2-c3)*(ro3-dabs(ro))/(ro3-y22)+c3
           if(r3.lt.1.d-3)then
           xx01=c3
           xx02=c3
           endif

     
         if(z.ge.xx01.and.z.le.xx02)goto 6677 ! regiunea mediana

             else
                        ! s negativ
      y11=b1*dsqrt(1.d0-((x1-c1)/a1)**2)
      y22=b2*dsqrt(1.d0-((x2-c2)/a2)**2)

      xx01=(x1-c3)*(dabs(ro)+s*ro3)/(s*ro3+y11)+c3
      xx02=(x2-c3)*(dabs(ro)+s*ro3)/(s*ro3+y22)+c3
         if(z.le.xx02.and.z.ge.xx01)goto 6677 ! regiunea mediana
             endif
     
      goto 6688 ! mergi intr-una din elipse

6677  continue ! regiunea mediana

      ! ramane cazul intermediar pentru r3
       V0=(vp1+vp2)/2.d0
                delta=dsqrt((c3-z)**2+(ro-ro3)**2)-r3
                if(s.gt.0.d0)delta=-delta
      !derddz=-(z-c3)/delta
      !derdro=-(ro-ro3)/delta
        if(s.gt.0.d0)then ! am impresia ca ro3 are semn
        derddz=-(z-c3)/dsqrt((c3-z)**2+(ro-ro3)**2)
        derdro=-(ro-ro3)/dsqrt((c3-z)**2+(ro-ro3)**2)
        else
        derddz=(z-c3)/dsqrt((c3-z)**2+(ro-ro3)**2)
        derdro=(ro-ro3)/dsqrt((c3-z)**2+(ro-ro3)**2)
        endif
                goto 5555
6688  continue
      if(z.gt.xx02)goto 5553
5551  continue
      V0=vp1
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0 !daca este 0 se cauta si solutii in interior elipsa
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c1)
      if(xxx.gt.a1)nmaxio=1
      if(xxx.lt.a1)then
      yyy=b1*dsqrt(1.d0-xxx**2/a1**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta1x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 7177
      xri1=1.d0*i
      d2xr1=delta1x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 6161
      endif
      xli1=xri1
      d2xi1=d2xr1
7177  continue
       if(nmaxio.eq.0)then ! distanta delta nu poate fi mai mare decat b1
      xri2=-1.d0*i
      alimin=b1+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b1
      endif
      d2xr2=delta1x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 6161
      endif
      xli2=xri2
      d2xi2=d2xr2
       endif
      enddo
6161    continue
      CALL DRTMI(X,F,delta1x,XLI,XRI,1.e-4,100,IER)
      if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER5=',IER
      if(IER.eq.2)stop
      delta=x
! daca a1 sau b1 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a1.gt.dddd1.and.a1.le.dddd.and.a1.eq.b1)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a1.gt.dddd1.and.a1.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b1+delta)**2/(ro**2/(b1+delta)**3)
      derdro=0.d0
      else
      if(b1.gt.dddd1.and.b1.le.dddd)then
      derdro=1.d0
!      derddz=(z-c1)/(a1+delta)**2/((z-c1)**2/(a1+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c1)/(a1+delta)**2
      derddz=derddz/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      derdro=ro/(b1+delta)**2
      derdro=derdro/((z-c1)**2/(a1+delta)**3+ro**2/(b1+delta)**3)
      endif
      endif
      endif
      goto 5555

5553  continue
      V0=vp2
         ! incep sa caut limitele pentru rezolvare
      nmaxio=0
! fac un test sa vad daca punctul este in afara elipsei
      xxx=dabs(z-c2)
      if(xxx.gt.a2)nmaxio=1
      if(xxx.lt.a2)then
      yyy=b2*dsqrt(1.d0-xxx**2/a2**2)
      if(yyy.lt.ro)nmaxio=1
      endif
! fac o bucla cu toate valorile -50 la 50
         xli1=0.d0
         xli2=0.d0
         d2xi1=delta2x(xli1)
         d2xi2=d2xi1
      do i=1,50
      if(nmaxio.eq.0)goto 7178
      xri1=1.d0*i
      d2xr1=delta2x(xri1)
      p1=d2xi1*d2xr1
      if(p1.le.0.d0)then
      xli=xli1
      xri=xri1
      goto 6262
      endif
      xli1=xri1
      d2xi1=d2xr1
7178  continue
       if(nmaxio.eq.0)then
      xri2=-1.d0*i
      alimin=b2+xri2
      if(alimin.lt.0.d0)then
      nmaxio=1
      xri2=-b2
      endif
      d2xr2=delta2x(xri2)
      p2=d2xi2*d2xr2
      if(p2.le.0.d0)then
      xli=xri2
      xri=xli2
      goto 6262
      endif
      xli2=xri2
      d2xi2=d2xr2
        endif
      enddo
6262  continue
      CALL DRTMI(X,F,delta2x,XLI,XRI,1.e-4,100,IER)
      if(IER.ne.0)print*,'!!!!!!!!!!!!!! IER6=',IER
      if(IER.eq.2)stop
      delta=x
! daca a2 sau b2 egal cu -delta
      dddd=-delta
      dddd1=-delta-1.d-4
      if(a2.gt.dddd1.and.a2.le.dddd.and.a2.eq.b2)then !caz special sfera
      derddz=1.d0
      derdro=1.d0
      else
      if(a2.gt.dddd1.and.a2.le.dddd)then
      derddz=1.d0
!      derdro=ro/(b2+delta)**2/(ro**2/(b2+delta)**3)
      derdro=0.d0
      else
      if(b2.gt.dddd1.and.b2.le.dddd)then
      derdro=1.d0
!      derddz=(z-c2)/(a2+delta)**2/((z-c2)**2/(a2+delta)**3)
      derdrz=0.d0
      else
      derddz=(z-c2)/(a2+delta)**2
      derddz=derddz/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      derdro=ro/(b2+delta)**2
      derdro=derdro/((z-c2)**2/(a2+delta)**3+ro**2/(b2+delta)**3)
      endif
      endif
      endif
      goto 5555
5555  continue
      derddz=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derddz
      derdro=+V0/ca*dexp(delta/ca)/(1+dexp(delta/ca))**2*derdro
      return
      end

      subroutine lminplus(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
      double precision lminus,lminus0
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2

      common/zilim12/z1lim,z2lim
      external lminus,lminus0
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
c      z1lim=cc1-a1x-4.d0
c      z2lim=cc2+a2x+4.d0
      xl=alf1*z1lim
      xu=0.d0
           CALL IGAUSS1(lminus,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(lminus,xl,xu,Y2)
           call dqh64(lminus0,p)
      pert=y1+y2+p
      return
      end

      double  precision function lminus0(z)
      implicit double precision (a-h,o-z)
      double precision lminu1,lminu2
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      common/zetal/zeta
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/fdez/valfz
      external lminu1,lminu2 
           nro1x=n1
           nro2x=n1p
           m1x=mp1
           m2x=mp1p     
c     print*,' '
c      zeta=z    
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
      zeta=zp*alf1
c      zz=z
c      a1=an1
c      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
      zz=zeta+dabs(zet1)
      a1=an1
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an1-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(-alf1*zz-2.d0*an1*alf1*raport)
            CALL DQL32 (lminu1,Y)
c      lminus=y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      lminus0=y*C1*C1P/alf1
      else
c      zz=-z
      zp=(z+zet2)/alf2
      zeta=zp*alf2
c      a1=an2
c      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
      zz=-(zeta-zet2)
      a1=an2
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an2-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(alf2*zz+2.d0*an2*alf2*raport)
            CALL DQL32 (lminu2,Y)
c integrala se face pe (alfa*ro)**2
c      lminus=y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      lminus0=y*c2*c2p/alf2
      endif
      return
      end

      double  precision function lminus(z)
      implicit double precision (a-h,o-z)
      double precision lminu11,lminu22
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      common/zetal/zeta
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/fdez/valfz
      external lminu11,lminu22 
           nro1x=n1
           nro2x=n1p
           m1x=mp1
           m2x=mp1p     
c     print*,' '
      zeta=z    
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
      zz=zeta+dabs(zet1)
      a1=an1
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an1-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(-alf1*zz-2.d0*an1*alf1*raport)
            CALL DQL32 (lminu11,Y)
c      lminus=y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      lminus=y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
      zz=-(zeta-zet2)
      a1=an2
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an2-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(alf2*zz+2.d0*an2*alf2*raport)
            CALL DQL32 (lminu22,Y)
c integrala se face pe (alfa*ro)**2
c      lminus=y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      lminus=y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end

      double precision function lminu1(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz1(nrnodzz,nrnod)
       derro=derdro1(nrnodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=-derz/ro*(mp1+iabs(mp1)+2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=zeta+dabs(zet1z)
c      a1=an1
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an1-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(-alf1z*zz-2.d0*an1*alf1z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lminu1=pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
! 

      double precision function lminu2(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz2(nrnodzz,nrnod)
       derro=derdro2(nrnodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=-derz/ro*(mp1+iabs(mp1)+2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=-(zeta-zet2z)
c      a1=an2
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an2-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(alf2z*zz+2.d0*an2*alf2z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lminu2=pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end

      double precision function lminu11(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
c      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz11(nodzz,nrnod)
       derro=derdro11(nodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=-derz/ro*(mp1+iabs(mp1)+2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=zeta+dabs(zet1z)
c      a1=an1
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an1-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(-alf1z*zz-2.d0*an1*alf1z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lminu11=pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
! 

      double precision function lminu22(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
c      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz22(nodzz,nrnod)
       derro=derdro22(nodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=-derz/ro*(mp1+iabs(mp1)+2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=-(zeta-zet2z)
c      a1=an2
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an2-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(alf2z*zz+2.d0*an2*alf2z*raport)
      pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lminu22=pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
! 
! 

      subroutine lplussmin(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
      double precision lplus,lplus0
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2

      common/zilim12/z1lim,z2lim
      external lplus,lplus0
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
c      z1lim=cc1-a1x-4.d0
c      z2lim=cc2+a2x+4.d0
      xl=alf1*z1lim
      xu=0.d0
           CALL IGAUSS1(lplus,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(lplus,xl,xu,Y2)
           call dqh64(lplus0,p)
      pert=y1+y2+p
      return
      end

      double  precision function lplus0(z)
      implicit double precision (a-h,o-z)
      double precision lplus1,lplus2
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      common/zetal/zeta
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/fdez/valfz
      external lplus1,lplus2
           nro1x=n1
           nro2x=n1p
           m1x=mp1
           m2x=mp1p
c     print*,' '
c      zeta=z
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
      zeta=zp*alf1
c      zz=z
c      a1=an1
c      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
      zz=zeta+dabs(zet1)
      a1=an1
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an1-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(-alf1*zz-2.d0*an1*alf1*raport)
            CALL DQL32 (lplus1,Y)
c      lminus=y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      lplus0=y*C1*C1P/alf1
      else
c      zz=-z
      zp=(z+zet2)/alf2
      zeta=zp*alf2
c      a1=an2
c      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
      zz=-(zeta-zet2)
      a1=an2
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an2-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(alf2*zz+2.d0*an2*alf2*raport)
            CALL DQL32 (lplus2,Y)
c integrala se face pe (alfa*ro)**2
c      lminus=y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      lplus0=y*c2*c2p/alf2
      endif
      return
      end


      double  precision function lplus(z)
      implicit double precision (a-h,o-z)
      double precision lplus11,lplus22
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      common/zetal/zeta
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/fdez/valfz
      external lplus11,lplus22
           nro1x=n1
           nro2x=n1p
           m1x=mp1
           m2x=mp1p
c     print*,' '
      zeta=z
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
      zz=zeta+dabs(zet1)
      a1=an1
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an1-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(-alf1*zz-2.d0*an1*alf1*raport)
            CALL DQL32 (lplus11,Y)
c      lminus=y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      lplus=y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
      zz=-(zeta-zet2)
      a1=an2
      CALL DHERM(A1,ZZ,H11,D1)
      a1=an2-1.d0
      CALL DHERM(A1,ZZ,H12,D1)
      h112=h11
      if(dabs(h112).lt.1.d-12)h112=1.d-12
      raport=h12/h112
      valfz=(alf2*zz+2.d0*an2*alf2*raport)
            CALL DQL32 (lplus22,Y)
c integrala se face pe (alfa*ro)**2
c      lminus=y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      lplus=y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end

      double precision function lplus1(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz1(nrnodzz,nrnod)
       derro=derdro1(nrnodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=derz/ro*(mp1-iabs(mp1)-2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=zeta+dabs(zet1z)
c      a1=an1
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an1-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(-alf1z*zz-2.d0*an1*alf1z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lplus1=-pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
!

      double precision function lplus2(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz2(nrnodzz,nrnod)
       derro=derdro2(nrnodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=derz/ro*(mp1-iabs(mp1)-2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=-(zeta-zet2z)
c      a1=an2
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an2-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(alf2z*zz+2.d0*an2*alf2z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lplus2=-pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
!

      double precision function lplus11(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
c      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz11(nodzz,nrnod)
       derro=derdro11(nodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=derz/ro*(mp1-iabs(mp1)-2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=zeta+dabs(zet1z)
c      a1=an1
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an1-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(-alf1z*zz-2.d0*an1*alf1z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)

c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lplus11=-pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
!

      double precision function lplus22(rho)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      COMMON/ALF12/ALF1z,ALF2z
      COMMON/ZET12/ZET1z,ZET2z
      common/zetal/zeta
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
c      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
      common/fdez/valfz
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,derz,derro)
       derz=derddz22(nodzz,nrnod)
       derro=derdro22(nodzz,nrnod)
      n=n1
      k=iabs(mp1)
      pl1=plag(n,k,nrnod)
c      call laguerre(rho,n,k,pla1)
      n=n1-1
      pl2=0
      if(n.ge.0)      pl2=plag(n,k,nrnod)
c      if(n.lt.k.or.n.lt.0.d0)then
c      plam=0.d0
c      else
c      call laguerre(rho,n,k,plam)
c      endif
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
c      roal=rho
c      pertur11=-derz*alf1ro*(-roal+(2.d0*n1+1+mp1)/roal-
c     -  2.d0*(n1+iabs(mp1))/roal*plam/pla1)
       ccoef=0
       if(n1.gt.0)ccoef=dsqrt(1.d0*n1/(n1+iabs(mp1)))
       pertur11=derz/ro*(mp1-iabs(mp1)-2*n1)+
     c         derz*alf1ro**2*ro+
     c         derz*2*(n1+iabs(mp1))/ro*pl2/pl1*ccoef
! calculam si functiile Hermite
c      zz=-(zeta-zet2z)
c      a1=an2
c      CALL DHERM(A1,ZZ,H11,D1)
c      a1=an2-1.d0
c      CALL DHERM(A1,ZZ,H12,D1)
c      raport=h12/h11
c      pertur11=pertur11+derro*(alf2z*zz+2.d0*an2*alf2z*raport)
       pertur11=pertur11+derro*valfz
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lminu1=pertur11*2.d0*alf1ro**2*
c     c dsqrt(aan1*aan1p/(amm1*amm2))*pla1*pla2*
c     c (rho)**(0.5d0*(iabs(mp1)+iabs(mp1p)))*dsqrt(rho) !*dexp(-rho)
      lplus22=-pertur11*
     c (dsqrt(rho))**(iabs(mp1)+iabs(mp1p)) !*dexp(-rho)
      return
      end
!

! 

      subroutine lzszi(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
      double precision lzsz1i,lzsz0i
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      common/zilim12/z1lim,z2lim
      common/numintegr/numigauss32,numdqh64
      external lzsz1i,lzsz0i  
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(lzsz1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(lzsz1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(lzsz0i,p)
      pert=p+y1+y2
      return
      end



      double  precision function lzsz0i(z)
      implicit double precision (a-h,o-z)
      double precision lzsz11,lzsz22
      real alzsz064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/lzsz64/alzsz064(64,0:22,0:22,0:22,0:22)
      external lzsz11,lzsz22
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (lzsz11,Y)
      lzsz0i=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (lzsz22,Y)
c integrala se face pe (alfa*ro)**2
      lzsz0i=Y*c2*c2p/alf2
      endif
      alzsz064(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end


      double  precision function lzsz1i(z)
      implicit double precision (a-h,o-z)
      double precision lzsz110,lzsz220
      real alzsz132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/lzsz32/alzsz132(64,0:22,0:22,0:22,0:22)
      external lzsz110,lzsz220
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (lzsz110,Y)
      lzsz1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (lzsz220,Y)
c integrala se face pe (alfa*ro)**2
      lzsz1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      alzsz132(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end

      subroutine lzsz(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
      double precision lzsz1,lzsz0
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      common/zilim12/z1lim,z2lim
      common/numintegr/numigauss32,numdqh64
      external lzsz1,lzsz0  
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
c      z1lim=cc1-a1x-4.d0
c      z2lim=cc2+a2x+4.d0
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(lzsz1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(lzsz1,xl,xu,Y2)
      numdqh64=0
      call dqh64(lzsz0,p)
      pert=p+y1+y2
      return
      end



      double  precision function lzsz0(z)
      implicit double precision (a-h,o-z)
      double precision lzsz11,lzsz22
      real alzsz064
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/lzsz64/alzsz064(64,0:22,0:22,0:22,0:22)
      external lzsz11,lzsz22
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c      zz=z
c      a1=an1
c      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (lzsz11,Y)
c      lzsz1=Y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      y=alzsz064(numdqh64,nro1x,nro2x,m1x,m2x)
      lzsz0=Y*C1*C1P/alf1
      else
c      zz=-z
      zp=(z+zet2)/alf2
c      a1=an2
c      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (lzsz22,Y)
c integrala se face pe (alfa*ro)**2
c      lzsz1=Y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      y=alzsz064(numdqh64,nro1x,nro2x,m1x,m2x)
      lzsz0=Y*c2*c2p/alf2
      endif
      return
      end


      double  precision function lzsz1(z)
      implicit double precision (a-h,o-z)
      double precision lzsz110,lzsz220
      real alzsz132
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/lzsz32/alzsz132(64,0:22,0:22,0:22,0:22)
      external lzsz110,lzsz220
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)     
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (lzsz110,Y)
c      lzsz1=Y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      y=alzsz132(numigauss32,nro1x,nro2x,m1x,m2x)
      lzsz1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (lzsz220,Y)
c integrala se face pe (alfa*ro)**2
c      lzsz1=Y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      y=alzsz132(numigauss32,nro1x,nro2x,m1x,m2x)
      lzsz1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      return
      end

      double precision function lzsz11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,der1,f)
       f=derdro1(nrnodzz,nrnod)
c      n=n1
c      k=mp1
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur11=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lzsz11=pertur11*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      lzsz11=pertur11*
     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      return
      end

      double precision function lzsz22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro
c      call dervVx(z,ro,der1,f)
       f=derdro2(nrnodzz,nrnod)
c      n=n1
c      k=mp1
c      rho2=rho**2
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur22=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lzsz22=pertur22*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      lzsz22=pertur22*
     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      return
      end


      double precision function lzsz110(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
c      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,der1,f)
       f=derdro11(nodzz,nrnod)
c      n=n1
c      k=mp1
c      rho2=rho**2
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur22=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lzsz22=pertur22*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      lzsz110=pertur22*
     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      return
      end

      double precision function lzsz220(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
c      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c      call dervVx(z,ro,der1,f)
       f=derdro22(nodzz,nrnod)
c      n=n1
c      k=mp1
c      rho2=rho**2
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=mp1p
c      call laguerre(rho,n,k,pla2)
      pertur22=f
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
c      lzsz22=pertur22*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      lzsz220=pertur22*
     c (rho)**(iabs(mp1))/ro !*dexp(-rho)
      return
      end


      subroutine pertwsi(nro1,m1,
     c  nro1p,m1p)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2
      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external perturbi,perturb1i
      n1=nro1
      mp1=m1
      n1p=nro1p
      mp1p=m1p
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(perturb1i,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(perturb1i,xl,xu,Y2)
      numdqh64=0
      call dqh64(perturbi,pert)
      pert=pert+y1+y2
      return
      end

      double  precision function perturbi(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real perturb64
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/pertur64/perturb64(64,0:22,0:22,0:22,0:22)
      external pertur11,pertur22
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
            CALL DQL32 (pertur11,Y)
      perturbi=Y*C1*C1P/alf1
      else
      zp=(z+zet2)/alf2
            CALL DQL32 (pertur22,Y)
c integrala se face pe (alfa*ro)**2
      perturbi=Y*c2*c2p/alf2
      endif
      perturb64(numdqh64,nro1x,nro2x,m1x,m2x)=y
      return
      end

      double  precision function perturb1i(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real perturb32
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/pertur32/perturb32(64,0:22,0:22,0:22,0:22)
      external pertur00,pertur02
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
            CALL DQL32 (pertur00,Y)
      perturb1i=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
            CALL DQL32 (pertur02,Y)
c integrala se face pe (alfa*ro)**2
      perturb1i=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      endif
      perturb32(numigauss32,nro1x,nro2x,m1x,m2x)=y
      return
      end



      subroutine pertws(aniu1,aniu2,const1,const2,nro1,m1,
     c  aniu1p,aniu2p,const1p,const2p,nro1p,m1p,pert)
      implicit double precision (a-h,o-z)
c calculez perturbatiile wood-saxon
c trebuie sa fac o integrare a functiei perturb1
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/hwuri/HWRO1,HWZ1,HWRO2,HWZ2
      common/adancimi/vpo1,vpo2

      common/zilim12/z1lim,z2lim
      common/oadancime/v000
      common/adancimea/vadanc
      common/numintegr/numigauss32,numdqh64
      external perturb,perturb1
      an1=aniu1
      an2=aniu2
      c1=const1
      c2=const2
      n1=nro1
      mp1=m1
      an1p=aniu1p
      an2p=aniu2p
      c1p=const1p
      c2p=const2p
      n1p=nro1p
      mp1p=m1p
c      z1lim=cc1-a1x-4.d0
c      z2lim=cc2+a2x+4.d0

c call jcptn(rez,x1n,x2n,zet) I2N1N2(rez,x1n,x2n,zet)
      xl=alf1*z1lim
      xu=0.d0
      numigauss32=0
           CALL IGAUSS1(perturb1,xl,xu,Y1)
      xu=alf2*z2lim
      xl=0.d0
           CALL IGAUSS2(perturb1,xl,xu,Y2)
c      pert=y1+y2!+v1+v2


      numdqh64=0
      call dqh64(perturb,pert)
            
      pert=pert+y1+y2

c      ai1=0
c      if(n1.eq.n1p.or.n1p.eq.(n1+1).or.n1p.eq.(n1-1))then
c      zzz=(z1lim-cc1)*alf1
c      call i2n1n2(va1,an1p,an1,zzz)
c      call jcptn(va,an1p,an1,zzz)
c      if(n1.eq.n1p)then
c      ai1=c1*c1p/alf1*hwz1/2.d0*va1
c      ai1=ai1-c1*c1p/alf1*vadanc*va
c      ai1=ai1+c1*c1p/alf1*hwro1/2.d0*va*(1+2*n1+iabs(mp1))
c      endif
c      if(n1+1.eq.n1p)then
c      ai1=-c1*c1p/alf1*hwro1/2.d0*va*dsqrt((n1+1.d0)*(n1+iabs(mp1)+1))
c      endif
c      if(n1-1.eq.n1p)then
c      ai1=-c1*c1p/alf1*hwro1/2.d0*va*dsqrt(n1*(n1+iabs(mp1))*1.d0)
c      endif
c      endif
c      v1=-ai1

c           print*,'va1',va1,an1p,an1
c      ai2=0
c      if(n1.eq.n1p.or.n1p.eq.(n1+1).or.n1p.eq.(n1-1))then
c      zzz=-(z2lim-cc2)*alf2
c      call i2n1n2(va1,an2p,an2,zzz)
c      call jcptn(va,an2p,an2,zzz)
c      if(n1.eq.n1p)then
c      ai2=c2*c2p/alf2*hwz2/2.d0*va1
c      ai2=ai2-c2*c2p/alf2*vadanc*va
c      ai2=ai2+c2*c2p/alf2*hwro1/2.d0*va*(1+2*n1+iabs(mp1))
c      endif
c      if(n1+1.eq.n1p)then
c      ai2=-c2*c2p/alf2*hwro1/2.d0*va*dsqrt((n1+1.d0)*(n1+iabs(mp1)+1))
c      endif
c      if(n1-1.eq.n1p)then
c      ai2=-c2*c2p/alf2*hwro1/2.d0*va*dsqrt(n1*(n1+iabs(mp1))*1.d0)
c      endif
c      endif
c      v2=-ai2
c            pert=pert+v1+v2
c          print*,'va1',va1,an2p,an2
c      stop
      return
      end

      double  precision function perturb(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real perturb64
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/pertur64/perturb64(64,0:22,0:22,0:22,0:22)
      external pertur11,pertur22
         numdqh64=numdqh64+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=(z-zet1)/alf1
c      zz=z
c      a1=an1
c      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (pertur11,Y)
c      print*,'!!!!!!!!!! y',y
c      perturb=Y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
c      perturb=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      y=perturb64(numdqh64,nro1x,nro2x,m1x,m2x)
      perturb=Y*C1*C1P/alf1
      else
c      zz=-z
      zp=(z+zet2)/alf2
c      a1=an2
c      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (pertur22,Y)
c integrala se face pe (alfa*ro)**2
c      perturb=Y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
c      perturb=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      y=perturb64(numdqh64,nro1x,nro2x,m1x,m2x)
      perturb=Y*c2*c2p/alf2
      endif
      return
      end

      double  precision function perturb1(z)
c aici z este notatia pentru zeta=alfa z
      implicit double precision (a-h,o-z)
      real perturb32
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/paramel/a1x,b1x,a2x,b2x,cc1,cc2,c3x,x1x,x2x,
     c            r3x,u1x,u2x,ro3x,sx,deltx
      common/zpert/zp
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET12/ZET1,ZET2
      common/optest/optest
      common/nrom/nro1x,nro2x,m1x,m2x
      common/numintegr/numigauss32,numdqh64
      common/pertur32/perturb32(64,0:22,0:22,0:22,0:22)
      external pertur00,pertur02
      numigauss32=numigauss32+1
           nro1x=n1
           nro2x=n1p
           m1x=iabs(mp1)
           m2x=iabs(mp1p)
      if(z.lt.0.d0)then
      zp=z/alf1
      zz=z+dabs(zet1)
      a1=an1
      a2=an1p
c      CALL DHERM(A1,ZZ,H1,D1)
c      CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (pertur00,Y)
c      print*,'!!!!!!!!!! y',y
c      perturb=Y*C1*C1P/alf1*H1*H2*dexp(-z**2-zet1**2+2*dabs(z*zet1))
      y=perturb32(numigauss32,nro1x,nro2x,m1x,m2x)
      perturb1=Y*C1*C1P/alf1*dexp(-z**2-zet1**2+2*dabs(z*zet1))
c      perturb=Y*C1*C1P/alf1*dexp(-zet1**2+2*dabs(z*zet1))*h1*h2
      else
      zz=-(z-zet2)
      zp=z/alf2
      a1=an2
      a2=an2p
c       CALL DHERM(A1,ZZ,H1,D1)
c       CALL DHERM(A2,ZZ,H2,D2)
c            CALL DQL32 (pertur02,Y)
c integrala se face pe (alfa*ro)**2
c      perturb=Y*c2*c2p/alf2*h1*h2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
      y=perturb32(numigauss32,nro1x,nro2x,m1x,m2x)
      perturb1=Y*c2*c2p/alf2*dexp(-z**2-zet2**2+2*dabs(z*zet2))
c      perturb1=Y*c2*c2p/alf2*dexp(-zet2**2+2*dabs(z*zet2))*h1*h2
      endif
      return
      end



      double precision function pertur11(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vpoten1/val1(32,32) !z si rho (z negativ)
      common/vpoten2/val2(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c       f=dif_ws_os(z,ro)
       f=val1(nrnodzz,nrnod)
c      n=n1
c      k=iabs(mp1)
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=iabs(mp1p)
c      call laguerre(rho,n,k,pla2)
      pertur11=f                     !*pla1*pla2
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
      pertur11=pertur11*(rho)**(iabs(mp1))
c      pertur11=pertur11*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1)) !*dexp(-rho)
      return
      end

      double precision function pertur22(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
      common/lagnod/nrnod,nrnodzz,nodzz
      common/vpoten1/val1(32,32) !z si rho (z negativ)
      common/vpoten2/val2(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W2RO=W0*R0/b1
      alf2ro=DSQRT(MPH*W2ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf2ro
c      f=dif_ws_os(z,ro)
       f=val2(nrnodzz,nrnod)
c      n=n1
c      k=iabs(mp1)
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=iabs(mp1p)
c      call laguerre(rho,n,k,pla2)
      pertur22=f         !*pla1*pla2
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
      pertur22=pertur22*(rho)**(iabs(mp1))
c      pertur22=pertur22*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1)) !*dexp(-rho)
      return
      end


      double precision function pertur00(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/vpoten11/vall1(32,32) !z si rho (z negativ)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c       f=dif_ws_os(z,ro)
      f=vall1(nodzz,nrnod)
c      n=n1
c      k=iabs(mp1)
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=iabs(mp1p)
c      call laguerre(rho,n,k,pla2)
      pertur11=f                     !*pla1*pla2
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
      pertur00=pertur11*(rho)**(iabs(mp1))
c      pertur11=pertur11*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1)) !*dexp(-rho)
      return
      end

      double precision function pertur02(rho)
      implicit double precision (a-h,o-z)
      common/ppertur/an1,an2,c1,c2,n1,mp1
      common/pperturp/an1p,an2p,c1p,c2p,n1p,mp1p
      common/zpert/zp
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/lagnod/nrnod,nrnodzz,nodzz
      common/CKST/CKST
      common/vpoten22/vall2(32,32) ! z si rho (z pozitiv)
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1RO=W0*R0/b1
      alf1ro=DSQRT(MPH*W1ro)*1.D-15
      z=zp
      ro=dsqrt(rho)/alf1ro
c       f=dif_ws_os(z,ro)
      f=vall2(nodzz,nrnod)
c      n=n1
c      k=iabs(mp1)
c      call laguerre(rho,n,k,pla1)
c      n=n1p
c      k=iabs(mp1p)
c      call laguerre(rho,n,k,pla2)
      pertur11=f                     !*pla1*pla2
c      call vect(n1,aan1)
c      call vect(n1p,aan1p)
c      mm=n1+iabs(mp1)
c      call vect(mm,amm1)
c      mm=n1p+iabs(mp1p)
c      call vect(mm,amm2)
      pertur02=pertur11*(rho)**(iabs(mp1))
c      pertur11=pertur11*dsqrt(aan1*aan1p/(amm1*amm2))*
c     c (rho)**(iabs(mp1)) !*dexp(-rho)
      return
      end

C
C     ..................................................................
C
C        SUBROUTINE DQH64
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
C                               -INFINITY TO +INFINITY).
C
C        USAGE
C           CALL DQH64 (FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 64-POINT GAUSSIAN-HERMITE
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 127.
C           FOR REFERENCE, SEE
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
C           TR00.1100 (MARCH 1964), PP.213-214.
C
C     ..................................................................
C
      SUBROUTINE DQH64(FCT,Y)
C
C
      DOUBLE PRECISION X,Y,Z,FCT
        double precision herm1,herm2

      common/fherdqh/herm1(0:24,32),herm2(0:24,32)
        common/nznnll/nznn,nzll
C
      common/lagnod/nrnod,nrnodzz,nodzz
           nrnodzz=1
      X=.10526123167960546D2
      Z=-X
      Y=.55357065358569428D-48*(FCT(X)*herm2(nznn,1)*herm2(nzll,1)+
     c                          FCT(Z)*herm1(nznn,1)*herm1(nzll,1))
           nrnodzz=2
      X=.9895287586829539D1
      Z=-X
      Y=Y+.16797479901081592D-42*(FCT(X)*herm2(nznn,2)*herm2(nzll,2)+
     c                          FCT(Z)*herm1(nznn,2)*herm1(nzll,2))
           nrnodzz=3
      X=.9373159549646721D1
      Z=-X
      Y=Y+.34211380112557405D-38*(FCT(X)*herm2(nznn,3)*herm2(nzll,3)+
     c                          FCT(Z)*herm1(nznn,3)*herm1(nzll,3))
           nrnodzz=4
      X=.8907249099964770D1
      Z=-X
      Y=Y+.15573906246297638D-34*(FCT(X)*herm2(nznn,4)*herm2(nzll,4)+
     c                          FCT(Z)*herm1(nznn,4)*herm1(nzll,4))
           nrnodzz=5
      X=.8477529083379863D1
      Z=-X
      Y=Y+.25496608991129993D-31*(FCT(X)*herm2(nznn,5)*herm2(nzll,5)+
     c                          FCT(Z)*herm1(nznn,5)*herm1(nzll,5))
           nrnodzz=6
      X=.8073687285010225D1
      Z=-X
      Y=Y+.19291035954649669D-28*(FCT(X)*herm2(nznn,6)*herm2(nzll,6)+
     c                          FCT(Z)*herm1(nznn,6)*herm1(nzll,6))
           nrnodzz=7
      X=.7689540164040497D1
      Z=-X
      Y=Y+.7861797788925910D-26*(FCT(X)*herm2(nznn,7)*herm2(nzll,7)+
     c                          FCT(Z)*herm1(nznn,7)*herm1(nzll,7))
           nrnodzz=8
      X=.7321013032780949D1
      Z=-X
      Y=Y+.19117068833006428D-23*(FCT(X)*herm2(nznn,8)*herm2(nzll,8)+
     c                          FCT(Z)*herm1(nznn,8)*herm1(nzll,8))
           nrnodzz=9
      X=.69652411205511075D1
      Z=-X
      Y=Y+.29828627842798512D-21*(FCT(X)*herm2(nznn,9)*herm2(nzll,9)+
     c                          FCT(Z)*herm1(nznn,9)*herm1(nzll,9))
           nrnodzz=10
      X=.66201122626360274D1
      Z=-X
      Y=Y+.31522545665037814D-19*(FCT(X)*herm2(nznn,10)*herm2(nzll,10)+
     c                          FCT(Z)*herm1(nznn,10)*herm1(nzll,10))
           nrnodzz=11
      X=.62840112287748282D1
      Z=-X
      Y=Y+.23518847106758191D-17*(FCT(X)*herm2(nznn,11)*herm2(nzll,11)+
     c                          FCT(Z)*herm1(nznn,11)*herm1(nzll,11))
           nrnodzz=12
      X=.59556663267994860D1
      Z=-X
      Y=Y+.12800933913224380D-15*(FCT(X)*herm2(nznn,12)*herm2(nzll,12)+
     c                          FCT(Z)*herm1(nznn,12)*herm1(nzll,12))
           nrnodzz=13
      X=.56340521643499721D1
      Z=-X
      Y=Y+.52186237265908475D-14*(FCT(X)*herm2(nznn,13)*herm2(nzll,13)+
     c                          FCT(Z)*herm1(nznn,13)*herm1(nzll,13))
           nrnodzz=14
      X=.53183252246332709D1
      Z=-X
      Y=Y+.16283407307097204D-12*(FCT(X)*herm2(nznn,14)*herm2(nzll,14)+
     c                          FCT(Z)*herm1(nznn,14)*herm1(nzll,14))
           nrnodzz=15
      X=.50077796021987682D1
      Z=-X
      Y=Y+.39591777669477239D-11*(FCT(X)*herm2(nznn,15)*herm2(nzll,15)+
     c                          FCT(Z)*herm1(nznn,15)*herm1(nzll,15))
           nrnodzz=16
      X=.47018156474074998D1
      Z=-X
      Y=Y+.7615217250145451D-10*(FCT(X)*herm2(nznn,16)*herm2(nzll,16)+
     c                          FCT(Z)*herm1(nznn,16)*herm1(nzll,16))
           nrnodzz=17
      X=.43999171682281376D1
      Z=-X
      Y=Y+.11736167423215493D-8*(FCT(X)*herm2(nznn,17)*herm2(nzll,17)+
     c                          FCT(Z)*herm1(nznn,17)*herm1(nzll,17))
           nrnodzz=18
      X=.41016344745666567D1
      Z=-X
      Y=Y+.14651253164761094D-7*(FCT(X)*herm2(nznn,18)*herm2(nzll,18)+
     c                          FCT(Z)*herm1(nznn,18)*herm1(nzll,18))
           nrnodzz=19
      X=.38065715139453605D1
      Z=-X
      Y=Y+.14955329367272471D-6*(FCT(X)*herm2(nznn,19)*herm2(nzll,19)+
     c                          FCT(Z)*herm1(nznn,19)*herm1(nzll,19))
           nrnodzz=20
      X=.35143759357409062D1
      Z=-X
      Y=Y+.12583402510311846D-5*(FCT(X)*herm2(nznn,20)*herm2(nzll,20)+
     c                          FCT(Z)*herm1(nznn,20)*herm1(nzll,20))
           nrnodzz=21
      X=.32247312919920357D1
      Z=-X
      Y=Y+.8788499230850359D-5*(FCT(X)*herm2(nznn,21)*herm2(nzll,21)+
     c                          FCT(Z)*herm1(nznn,21)*herm1(nzll,21))
           nrnodzz=22
      X=.29373508230046218D1
      Z=-X
      Y=Y+.51259291357862747D-4*(FCT(X)*herm2(nznn,22)*herm2(nzll,22)+
     c                          FCT(Z)*herm1(nznn,22)*herm1(nzll,22))
           nrnodzz=23
      X=.26519724354306350D1
      Z=-X
      Y=Y+.25098369851306249D-3*(FCT(X)*herm2(nznn,23)*herm2(nzll,23)+
     c                          FCT(Z)*herm1(nznn,23)*herm1(nzll,23))
           nrnodzz=24
      X=.23683545886324014D1
      Z=-X
      Y=Y+.10363290995075777D-2*(FCT(X)*herm2(nznn,24)*herm2(nzll,24)+
     c                          FCT(Z)*herm1(nznn,24)*herm1(nzll,24))
           nrnodzz=25
      X=.20862728798817620D1
      Z=-X
      Y=Y+.36225869785344588D-2*(FCT(X)*herm2(nznn,25)*herm2(nzll,25)+
     c                          FCT(Z)*herm1(nznn,25)*herm1(nzll,25))
           nrnodzz=26
      X=.18055171714655449D1
      Z=-X
      Y=Y+.10756040509879137D-1*(FCT(X)*herm2(nznn,26)*herm2(nzll,26)+
     c                          FCT(Z)*herm1(nznn,26)*herm1(nzll,26))
           nrnodzz=27
      X=.15258891402098637D1
      Z=-X
      Y=Y+.27203128953688918D-1*(FCT(X)*herm2(nznn,27)*herm2(nzll,27)+
     c                          FCT(Z)*herm1(nznn,27)*herm1(nzll,27))
           nrnodzz=28
      X=.12472001569431179D1
      Z=-X
      Y=Y+.58739981964099435D-1*(FCT(X)*herm2(nznn,28)*herm2(nzll,28)+
     c                          FCT(Z)*herm1(nznn,28)*herm1(nzll,28))
           nrnodzz=29
      X=.9692694230711780D0
      Z=-X
      Y=Y+.10849834930618684D0*(FCT(X)*herm2(nznn,29)*herm2(nzll,29)+
     c                          FCT(Z)*herm1(nznn,29)*herm1(nzll,29))
           nrnodzz=30
      X=.69192230581004458D0
      Z=-X
      Y=Y+.17168584234908370D0*(FCT(X)*herm2(nznn,30)*herm2(nzll,30)+
     c                          FCT(Z)*herm1(nznn,30)*herm1(nzll,30))
           nrnodzz=31
      X=.41498882412107868D0
      Z=-X
      Y=Y+.23299478606267805D0*(FCT(X)*herm2(nznn,31)*herm2(nzll,31)+
     c                          FCT(Z)*herm1(nznn,31)*herm1(nzll,31))
           nrnodzz=32
      X=.13830224498700972D0
      Z=-X
      Y=Y+.27137742494130398D0*(FCT(X)*herm2(nznn,32)*herm2(nzll,32)+
     c                          FCT(Z)*herm1(nznn,32)*herm1(nzll,32))
      RETURN
      END


      subroutine fherm(nm,x1d,x2d)
      implicit double precision (a-h,o-z)
        dimension x1d(0:24),x2d(0:24)
      common/fherdqh/herm1(0:24,32),herm2(0:24,32)
        COMMON/ZET12/ZET1,ZET2

 
             do k=0,nm
        a1=x1d(k)
        a2=x2d(k)
               
      X=.10526123167960546D2
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,1)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,1)=h1

      X=.9895287586829539D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,2)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,2)=h1

      X=.9373159549646721D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,3)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,3)=h1

      X=.8907249099964770D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,4)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,4)=h1


      X=.8477529083379863D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,5)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,5)=h1

      X=.8073687285010225D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,6)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,6)=h1

      X=.7689540164040497D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,7)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,7)=h1

      X=.7321013032780949D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,8)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,8)=h1

      X=.69652411205511075D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,9)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,9)=h1

      X=.66201122626360274D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,10)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,10)=h1

      X=.62840112287748282D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,11)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,11)=h1

      X=.59556663267994860D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,12)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,12)=h1

      X=.56340521643499721D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,13)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,13)=h1

      X=.53183252246332709D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,14)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,14)=h1

      X=.50077796021987682D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,15)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,15)=h1

      X=.47018156474074998D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,16)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,16)=h1

      X=.43999171682281376D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,17)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,17)=h1

      X=.41016344745666567D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,18)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,18)=h1

      X=.38065715139453605D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,19)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,19)=h1

      X=.35143759357409062D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,20)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,20)=h1

      X=.32247312919920357D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,21)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,21)=h1

      X=.29373508230046218D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,22)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,22)=h1

      X=.26519724354306350D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,23)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,23)=h1

      X=.23683545886324014D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,24)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,24)=h1

      X=.20862728798817620D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,25)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,25)=h1

      X=.18055171714655449D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,26)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,26)=h1

      X=.15258891402098637D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,27)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,27)=h1

      X=.12472001569431179D1
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,28)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,28)=h1

      X=.9692694230711780D0
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,29)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,29)=h1

      X=.69192230581004458D0
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,30)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,30)=h1

      X=.41498882412107868D0
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,31)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,31)=h1

      X=.13830224498700972D0
      zz=-x
      a2=x2d(k) 
      call dherm(a2,zz,h1,d1)
      herm2(k,32)=h1
      a1=x1d(k)
      call dherm(a1,zz,h1,d1)
      herm1(k,32)=h1
      enddo
      RETURN
      END


      subroutine vpot(xll,xuu)
c voi calcula valorile potentialului in cele 64 puncte z
c si cele 32 puncte ro pentru dqh64 si dql32

       implicit double precision (a-h,o-z)
      dimension ro(32),z(32),zz1(16),zz2(16)
      common/vpoten1/val1(32,32) !z si rho (z negativ)
      common/vpoten2/val2(32,32) ! z si rho (z pozitiv)
      common/vpoten11/vall1(32,32) !z si rho (z negativ)
      common/vpoten22/vall2(32,32) ! z si rho (z pozitiv)

      common/vcoul1/coulo1(32,32)
      common/vcoul2/coulo2(32,32)       
      common/vcoul1ext/coulo1ext(32,32)
      common/vcoul2ext/coulo2ext(32,32)



      common/vderiv1/derddz1(32,32),derdro1(32,32) ! z si rho (z negativ)
      common/vderiv2/derddz2(32,32),derdro2(32,32) ! z si rho (z pozitiv)
      common/vderiv11/derddz11(32,32),derdro11(32,32) ! z si rho (z negativ)
      common/vderiv22/derddz22(32,32),derdro22(32,32) ! z si rho (z pozitiv)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      common/r116/R0
      common/paramel/a1,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/CKST/CKST
c     ckst=41.d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      WRO=W0*R0/b1
      alfro=DSQRT(MPH*Wro)*1.D-15
      ro(1)=.11175139809793770D3
      ro(2)=.9882954286828397D2    
      ro(3)=.8873534041789240D2
      ro(4)=.8018744697791352D2
      ro(5)=.7268762809066271D2
      ro(6)=.65975377287935053D2
      ro(7)=.59892509162134018D2
      ro(8)=.54333721333396907D2
      ro(9)=.49224394987308639D2
      ro(10)=.44509207995754938D2
      ro(11)=.40145719771539442D2
      ro(12)=.36100494805751974D2
      ro(13)=.32346629153964737D2
      ro(14)=.28862101816323475D2
      ro(15)=.25628636022459248D2
      ro(16)=.22630889013196774D2
      ro(17)=.19855860940336055D2
      ro(18)=.17292454336715315D2
      ro(19)=.14931139755522557D2
      ro(20)=.12763697986742725D2
      ro(21)=.10783018632539972D2
      ro(22)=.8982940924212596D1
      ro(23)=.7358126733186241D1
      ro(24)=.59039585041742439D1
      ro(25)=.46164567697497674D1
      ro(26)=.34922132730219945D1
      ro(27)=.25283367064257949D1
      ro(28)=.17224087764446454D1
      ro(29)=.10724487538178176D1
      ro(30)=.57688462930188643D0
      ro(31)=.23452610951961854D0
      ro(32)=.44489365833267018D-1
              
      z(1)=.10526123167960546D2
      z(2)=.9895287586829539D1
      z(3)=.9373159549646721D1
      z(4)=.8907249099964770D1
      z(5)=.8477529083379863D1
      z(6)=.8073687285010225D1
      z(7)=.7689540164040497D1
      z(8)=.7321013032780949D1
      z(9)=.69652411205511075D1
      z(10)=.66201122626360274D1
      z(11)=.62840112287748282D1
      z(12)=.59556663267994860D1
      z(13)=.56340521643499721D1
      z(14)=.53183252246332709D1
      z(15)=.50077796021987682D1
      z(16)=.47018156474074998D1
      z(17)=.43999171682281376D1
      z(18)=.41016344745666567D1
      z(19)=.38065715139453605D1
      z(20)=.35143759357409062D1
      z(21)=.32247312919920357D1
      z(22)=.29373508230046218D1
      z(23)=.26519724354306350D1
      z(24)=.23683545886324014D1
      z(25)=.20862728798817620D1
      z(26)=.18055171714655449D1
      z(27)=.15258891402098637D1
      z(28)=.12472001569431179D1
      z(29)=.9692694230711780D0
      z(30)=.69192230581004458D0
      z(31)=.41498882412107868D0
      z(32)=.13830224498700972D0

        xu=0.d0
        xl=xll
	Az1=.5D0*(XU+XL)
	Bz1=XU-XL
	zz1(1)=.49863193092474078D0*Bz1
	zz1(2)=.49280575577263417D0*Bz1
c	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	zz1(3)=.48238112779375322D0*Bz1
	zz1(4)=.46745303796886984D0*Bz1
	zz1(5)=.44816057788302606D0*Bz1
	zz1(6)=.42468380686628499D0*Bz1
	zz1(7)=.39724189798397120D0*Bz1
	zz1(8)=.36609105937014484D0*Bz1
	zz1(9)=.33152213346510760D0*Bz1
	zz1(10)=.29385787862038116D0*Bz1
	zz1(11)=.25344995446611470D0*Bz1
	zz1(12)=.21067563806531767D0*Bz1
	zz1(13)=.16593430114106382D0*Bz1
	zz1(14)=.11964368112606854D0*Bz1
	zz1(15)=.7223598079139825D-1*Bz1
	zz1(16)=.24153832843869158D-1*Bz1
	
        xu=xuu
        xl=0.d0
	Az2=.5D0*(XU+XL)
	Bz2=XU-XL
	zz2(1)=.49863193092474078D0*Bz2
	zz2(2)=.49280575577263417D0*Bz2
c	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	zz2(3)=.48238112779375322D0*Bz2
	zz2(4)=.46745303796886984D0*Bz2
	zz2(5)=.44816057788302606D0*Bz2
	zz2(6)=.42468380686628499D0*Bz2
	zz2(7)=.39724189798397120D0*Bz2
	zz2(8)=.36609105937014484D0*Bz2
	zz2(9)=.33152213346510760D0*Bz2
	zz2(10)=.29385787862038116D0*Bz2
	zz2(11)=.25344995446611470D0*Bz2
	zz2(12)=.21067563806531767D0*Bz2
	zz2(13)=.16593430114106382D0*Bz2
	zz2(14)=.11964368112606854D0*Bz2
	zz2(15)=.7223598079139825D-1*Bz2
	zz2(16)=.24153832843869158D-1*Bz2
	

      do i=1,32
      rhoo=dsqrt(ro(i))/alfro
      do j=1,32
      z1p=(-z(j)-zet1)/alf1
      z2p=(z(j)+zet2)/alf2

       val1(j,i)=dif_ws_os(z1p,rhoo)
       val2(j,i)=dif_ws_os(z2p,rhoo)
      call dervVx(z1p,rhoo,derddz1(j,i),derdro1(j,i))
      call dervVx(z2p,rhoo,derddz2(j,i),derdro2(j,i))
        call vc(z1p,rhoo,coulo1ext(j,i))
        call vc(z2p,rhoo,coulo2ext(j,i))

      

      enddo
      enddo

      do i=1,32
      rhoo=dsqrt(ro(i))/alfro
      do j=1,16
      zz1p=(az1+zz1(j))/alf1
        jj=2*(j-1)+1
        vall1(jj,i)=dif_ws_os(zz1p,rhoo)
        call vc(zz1p,rhoo,vvv)
        coulo1(jj,i)=vvv
      call dervVx(zz1p,rhoo,derddz11(jj,i),derdro11(jj,i))
      zz1p=(az1-zz1(j))/alf1
        jj=jj+1
        vall1(jj,i)=dif_ws_os(zz1p,rhoo)
        call vc(zz1p,rhoo,vvv)
        coulo1(jj,i)=vvv
      call dervVx(zz1p,rhoo,derddz11(jj,i),derdro11(jj,i))
      zz2p=(az2+zz2(j))/alf2
        jj=2*(j-1)+1
        vall2(jj,i)=dif_ws_os(zz2p,rhoo)
        call vc(zz2p,rhoo,vvv)
        coulo2(jj,i)=vvv
      call dervVx(zz2p,rhoo,derddz22(jj,i),derdro22(jj,i))
      zz2p=(az2-zz2(j))/alf2
        jj=jj+1
        vall2(jj,i)=dif_ws_os(zz2p,rhoo)
        call vc(zz2p,rhoo,vvv)
        coulo2(jj,i)=vvv
      call dervVx(zz2p,rhoo,derddz22(jj,i),derdro22(jj,i))

      enddo
      enddo

       return
       end



C
C     ..................................................................
C
C        SUBROUTINE DQL32
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X
C                               FROM 0 TO INFINITY).
C
C        USAGE
C           CALL DQL32 (FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSSIAN-LAGUERRE
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 63.
C           FOR REFERENCE, SEE
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
C           TR00.1100 (MARCH 1964), PP.24-25.
C
C     ..................................................................
C
      SUBROUTINE DQL32(FCT,Y)
C
C
      DOUBLE PRECISION X,Y,FCT,PLAG
C

      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
      common/nrom/nro1x,nro2x,m1x,m2x

      common/lagnod/nrnod,nrnodzz,nodzz
      n1=nro1x
      n2=nro2x
      k1=iabs(m1x)
      k2=iabs(m2x)
      X=.11175139809793770D3
      nrnod=1
      Y=.45105361938989742D-47*FCT(X)*plag(n1,k1,1)*plag(n2,k2,1)
      X=.9882954286828397D2
      nrnod=2
      Y=Y+.13386169421062563D-41*FCT(X)*plag(n1,k1,2)*plag(n2,k2,2)
      X=.8873534041789240D2
      nrnod=3
      Y=Y+.26715112192401370D-37*FCT(X)*plag(n1,k1,3)*plag(n2,k2,3)
      X=.8018744697791352D2
      nrnod=4
      Y=Y+.11922487600982224D-33*FCT(X)*plag(n1,k1,4)*plag(n2,k2,4)
      X=.7268762809066271D2
      nrnod=5
      Y=Y+.19133754944542243D-30*FCT(X)*plag(n1,k1,5)*plag(n2,k2,5)
      X=.65975377287935053D2
      nrnod=6
      Y=Y+.14185605454630369D-27*FCT(X)*plag(n1,k1,6)*plag(n2,k2,6)
      X=.59892509162134018D2
      nrnod=7
      Y=Y+.56612941303973594D-25*FCT(X)*plag(n1,k1,7)*plag(n2,k2,7)
      X=.54333721333396907D2
      nrnod=8
      Y=Y+.13469825866373952D-22*FCT(X)*plag(n1,k1,8)*plag(n2,k2,8)
      X=.49224394987308639D2
      nrnod=9
      Y=Y+.20544296737880454D-20*FCT(X)*plag(n1,k1,9)*plag(n2,k2,9)
      X=.44509207995754938D2
      nrnod=10
      Y=Y+.21197922901636186D-18*FCT(X)*plag(n1,k1,10)*plag(n2,k2,10)
      X=.40145719771539442D2
      nrnod=11
      Y=Y+.15421338333938234D-16*FCT(X)*plag(n1,k1,11)*plag(n2,k2,11)
      X=.36100494805751974D2
      nrnod=12
      Y=Y+.8171823443420719D-15*FCT(X)*plag(n1,k1,12)*plag(n2,k2,12)
      X=.32346629153964737D2
      nrnod=13
      Y=Y+.32378016577292665D-13*FCT(X)*plag(n1,k1,13)*plag(n2,k2,13)
      X=.28862101816323475D2
      nrnod=14
      Y=Y+.9799379288727094D-12*FCT(X)*plag(n1,k1,14)*plag(n2,k2,14)
      X=.25628636022459248D2
      nrnod=15
      Y=Y+.23058994918913361D-10*FCT(X)*plag(n1,k1,15)*plag(n2,k2,15)
      X=.22630889013196774D2
      nrnod=16
      Y=Y+.42813829710409289D-9*FCT(X)*plag(n1,k1,16)*plag(n2,k2,16)
      X=.19855860940336055D2
      nrnod=17
      Y=Y+.63506022266258067D-8*FCT(X)*plag(n1,k1,17)*plag(n2,k2,17)
      X=.17292454336715315D2
      nrnod=18
      Y=Y+.7604567879120781D-7*FCT(X)*plag(n1,k1,18)*plag(n2,k2,18)
      X=.14931139755522557D2
      nrnod=19
      Y=Y+.7416404578667552D-6*FCT(X)*plag(n1,k1,19)*plag(n2,k2,19)
      X=.12763697986742725D2
      nrnod=20
      Y=Y+.59345416128686329D-5*FCT(X)*plag(n1,k1,20)*plag(n2,k2,20)
      X=.10783018632539972D2
      nrnod=21
      Y=Y+.39203419679879472D-4*FCT(X)*plag(n1,k1,21)*plag(n2,k2,21)
      X=.8982940924212596D1
      nrnod=22
      Y=Y+.21486491880136419D-3*FCT(X)*plag(n1,k1,22)*plag(n2,k2,22)
      X=.7358126733186241D1
      nrnod=23
      Y=Y+.9808033066149551D-3*FCT(X)*plag(n1,k1,23)*plag(n2,k2,23)
      X=.59039585041742439D1
      nrnod=24
      Y=Y+.37388162946115248D-2*FCT(X)*plag(n1,k1,24)*plag(n2,k2,24)
      X=.46164567697497674D1
      nrnod=25
      Y=Y+.11918214834838557D-1*FCT(X)*plag(n1,k1,25)*plag(n2,k2,25)
      X=.34922132730219945D1
      nrnod=26
      Y=Y+.31760912509175070D-1*FCT(X)*plag(n1,k1,26)*plag(n2,k2,26)
      X=.25283367064257949D1
      nrnod=27
      Y=Y+.70578623865717442D-1*FCT(X)*plag(n1,k1,27)*plag(n2,k2,27)
      X=.17224087764446454D1
      nrnod=28
      Y=Y+.12998378628607176D0*FCT(X)*plag(n1,k1,28)*plag(n2,k2,28)
      X=.10724487538178176D1
      nrnod=29
      Y=Y+.19590333597288104D0*FCT(X)*plag(n1,k1,29)*plag(n2,k2,29)
      X=.57688462930188643D0
      nrnod=30
      Y=Y+.23521322966984801D0*FCT(X)*plag(n1,k1,30)*plag(n2,k2,30)
      X=.23452610951961854D0
      nrnod=31
      Y=Y+.21044310793881323D0*FCT(X)*plag(n1,k1,31)*plag(n2,k2,31)
      X=.44489365833267018D-1
      nrnod=32
      Y=Y+.10921834195238497D0*FCT(X)*plag(n1,k1,32)*plag(n2,k2,32)
      RETURN
      END

      subroutine clag(nm)
      implicit double precision (a-h,o-z)
      dimension plag(0:24,0:24,32) !n,m,nrnod
      common/plag/plag
           !!!!!!! valorile clag se calculeaza o singura data in program principal
      do n=0,nm
      do k=0,nm
      call vect(n,aan1)
      mm=n+k
      call vect(mm,amm1)
      rad=dsqrt(aan1/amm1)
      X=.11175139809793770D3
             call laguerre(x,n,k,pla2)
      plag(n,k,1)=pla2*rad
      X=.9882954286828397D2
            call laguerre(x,n,k,pla2)
      plag(n,k,2)=pla2*rad    
      X=.8873534041789240D2
            call laguerre(x,n,k,pla2)      
      plag(n,k,3)=pla2*rad
      X=.8018744697791352D2
            call laguerre(x,n,k,pla2)
      plag(n,k,4)=pla2*rad
      X=.7268762809066271D2
            call laguerre(x,n,k,pla2)
      plag(n,k,5)=pla2*rad
      X=.65975377287935053D2
            call laguerre(x,n,k,pla2)
      plag(n,k,6)=pla2*rad
      X=.59892509162134018D2
            call laguerre(x,n,k,pla2)
      plag(n,k,7)=pla2*rad
      X=.54333721333396907D2
           call laguerre(x,n,k,pla2)
      plag(n,k,8)=pla2*rad
      X=.49224394987308639D2
            call laguerre(x,n,k,pla2)
      plag(n,k,9)=pla2*rad
      X=.44509207995754938D2
            call laguerre(x,n,k,pla2)
      plag(n,k,10)=pla2*rad
      X=.40145719771539442D2
           call laguerre(x,n,k,pla2)
      plag(n,k,11)=pla2*rad
      X=.36100494805751974D2
            call laguerre(x,n,k,pla2)
      plag(n,k,12)=pla2*rad
      X=.32346629153964737D2
            call laguerre(x,n,k,pla2)
      plag(n,k,13)=pla2*rad
      X=.28862101816323475D2
            call laguerre(x,n,k,pla2)
      plag(n,k,14)=pla2*rad
      X=.25628636022459248D2
            call laguerre(x,n,k,pla2)
      plag(n,k,15)=pla2*rad
      X=.22630889013196774D2
            call laguerre(x,n,k,pla2)
      plag(n,k,16)=pla2*rad
      X=.19855860940336055D2
            call laguerre(x,n,k,pla2)
      plag(n,k,17)=pla2*rad
      X=.17292454336715315D2
            call laguerre(x,n,k,pla2)
      plag(n,k,18)=pla2*rad
      X=.14931139755522557D2
            call laguerre(x,n,k,pla2)
      plag(n,k,19)=pla2*rad
      X=.12763697986742725D2
            call laguerre(x,n,k,pla2)
      plag(n,k,20)=pla2*rad
      X=.10783018632539972D2
            call laguerre(x,n,k,pla2)
      plag(n,k,21)=pla2*rad
      X=.8982940924212596D1
            call laguerre(x,n,k,pla2)
      plag(n,k,22)=pla2*rad
      X=.7358126733186241D1
            call laguerre(x,n,k,pla2)
      plag(n,k,23)=pla2*rad
      X=.59039585041742439D1
           call laguerre(x,n,k,pla2)
      plag(n,k,24)=pla2*rad
      X=.46164567697497674D1
            call laguerre(x,n,k,pla2)
      plag(n,k,25)=pla2*rad
      X=.34922132730219945D1
            call laguerre(x,n,k,pla2)
      plag(n,k,26)=pla2*rad
      X=.25283367064257949D1
            call laguerre(x,n,k,pla2)
      plag(n,k,27)=pla2*rad
      X=.17224087764446454D1
            call laguerre(x,n,k,pla2)
      plag(n,k,28)=pla2*rad
      X=.10724487538178176D1
            call laguerre(x,n,k,pla2)
      plag(n,k,29)=pla2*rad
      X=.57688462930188643D0
            call laguerre(x,n,k,pla2)
      plag(n,k,30)=pla2*rad
      X=.23452610951961854D0
            call laguerre(x,n,k,pla2)
      plag(n,k,31)=pla2*rad
      X=.44489365833267018D-1
      call laguerre(x,n,k,pla2)
      plag(n,k,32)=pla2*rad
      enddo
      enddo      
      RETURN
      END


      SUBROUTINE NRADAC(L,R1,R2,DELTA,C1,C2,X1,X2)
C ACEASTA SUBRUTINA OBTINE RADACINA NIU PENTRU VALOAREA L
C LA DOUA SFERE INTERSECTATE DACA SE CUNOSC R1,R2,DELTA
C SOLUTIILE SINT X1 SI X2
c aici R1 si R2 sunt semi-axele A1 si A2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real cucu1,cucu2
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      COMMON/ZET12/ZET1,ZET2
c     common/nrfin/nrfin1,nrfin2
      common/ptrads/detai,xdr1i,xdr2i,ndri ! vechile radacini ! ndri=0 se incepe calculul
                                           ! de la parametrizarea sfera, daca nu de la ultimele valori
      common/CKST/CKST
      EXTERNAL FNIU1,FNIU2
      ndri=0
c     CKST=41.D0
      n=l
      VINT1=.007D0
      VINT2=.001D0
      XST=L
c     print*,'XST,L',XST,L
      N=L
      AL=L
      LP2=L/2
      EAL=AL/2.D0-LP2
c     print*,'intru in reparm cu r1,r2,delta ',r1,r2,delta
      CALL REPARM(R1,R2,DELTA,VC,V2C,R22C,R0C)
c     print*,'ies din reparm cu VC,V2C,R22C,R0C',VC,V2C,R22C,R0C
      ALFA0=1.071186277015679D0/DSQRT(R0C) 
      DIF=1.D-1*ALFA0
      CT=R2/R1
      Z1=DELTA/(1.D0+CT)
      Z2=CT*Z1
c     print*,'1 z1,z2,ct',z1,z2,ct
c Z1 si Z2 sunt centrele celor doua elipse
      Z1=DABS(C1)
      Z2=DABS(C2)
      CT=Z2/Z1
c     print*,'2 z1,z2,ct',z1,z2,ct
      

      Z1C1=0.D0
      Z2C2=0.D0
c     print*,'ndri',ndri
         if(ndri.ne.0)then
      deti=detai
c     print*,'intru in REPARE(r0c,vc,v2c,r22c,deti,r1in,r2in'
c     print*,r0c,vc,v2c,r22c,deti,r1in,r2in
      CALL REPARE(r0c,vc,v2c,r22c,deti,r1in,r2in)
c     print*,'ies REPARE ',r0c,vc,v2c,r22c,deti,r1in,r2in
      cti=r2in/r1in
      z1in=deti/(1.d0+cti)
      llin=z1in/dif
          endif
          IF(EAL.LT..25D0) THEN ! L ESTE PAR
      if(ndri.ne.0)xst=xdr1i
      LL=Z1/DIF
          IF(LL.EQ.0) GO TO 2
          DO 1 J=1,LL
      Z1C1=Z1C1+DIF
      Z2C2=CT*Z1C1
      D=(Z1C1+Z2C2)
          if(ndri.ne.0.and.j.le.llin) go to 1
c     print*,'intru in REPARE(r0c,vc,v2c,r22c,d,r1c,r2c'
c     print*,r0c,vc,v2c,r22c,d,r1c,r2c
      CALL REPARE (R0C,VC,V2C,R22C,D,R1C,R2C)
c     print*,'ies REPARE ',r0c,vc,v2c,r22c,d,r1c,r2c
      CT1=R2C/R1C
      Z1C=D/(1.D0+CT1)
      Z2C=CT1*Z1C
      CALL P2C(R0C,R1C,R2C,Z1C,Z2C)
c     print*,'ies din p2c cu r0c,r1c,r2c,z1c,z2c'
c     print*,r0c,r1c,r2c,z1c,z2c
      XST1=XST
      XST2=XST
      FNST11=FNIU1(XST1)
c     print*,'FNST11(XST1)',FST11,XST1 
      FNST22=FNST11
       DO 101 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU1(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1 
c     print*,'FNST1,F1,XST1',FNST1,F1,XST1
C     F1=FNIU1(XST11)*FNIU1(XST1)
C     PRINT*,'F1',F1,'  XST1',XST1
          IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
c     print*,'als,ali',als,ali
      GO TO 3
          ENDIF
      FNST2=FNIU1(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
C     F1=FNIU1(XST22)*FNIU1(XST2)
C     PRINT*,'F1',F1,'    XST2',XST2
         IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
c     print*,'als,ali',als,ali
      GO TO 3
          ENDIF
101   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL, L=',L,' LL=',J
         RETURN
3     CONTINUE
      cucu1=1.e-5
      ncu=200
c       print*,'intru in drtmi ali,als',ali,als
      CALL DRTMI(X,FF,FNIU1,ALI,ALS,cucu1,ncu,IER)
c       print*,'x,ff,ier',x,ff,ier,'  PUN STOP'
      IF(IER.EQ.2)PRINT*,'IER1=',IER,'IN NRADAC L D',L,DELTA
      XST=X
1     CONTINUE
2     CONTINUE
      CALL P2C(R0C,R1,R2,Z1,Z2) 
c     print*,'ies din p2c cu r1c,r2,z1c,z2c'
c     print*,r0c,r1,r2,z1,z2
      XST1=XST
      XST2=XST
      FNST11=FNIU1(XST1)
c     print*,'FNST11(XST1)',FST11,XST1 
      FNST22=FNST11 
         DO 102 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU1(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
c     print*,'FNST1,F1,XST1',FNST1,F1,XST1
C     F1=FNIU1(XST11)*FNIU1(XST1)
C     PRINT*,'F1',F1
           IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 4
          ENDIF
      FNST2=FNIU1(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
c     print*,'FNST2,F1,XST2',FNST2,F1,XST2
C     F1=FNIU1(XST22)*FNIU1(XST2)
C     PRINT*,'F1',F1
          IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 4
          ENDIF
102   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL, L=',L
      RETURN
4     CONTINUE
      cucu2=1.e-7
      ncu=200
      CALL DRTMI(X,FF,FNIU1,ALI,ALS,cucu2,ncu,IER)
           IF(IER.EQ.2)PRINT*,'IER2=',IER,'IN NRADAC L DELTA',L,DELTA
      X1=X
      X2=(W1*(X1+.5D0))/W2-.5D0
c     print*,'x1,x2,w1,w2',x1,x2,w1,w2
c     stop
      RETURN
           ELSE
           if(ndri.ne.0)xst=xdr2i
      LL=Z2/DIF
      IF(LL.EQ.0) GO TO 5
      DO 6 J=1,LL
      Z2C2=Z2C2+DIF
      Z1C1=Z2C2/CT
      D=Z2C2+Z1C1
          if(ndri.ne.0.and.j.le.llin) goto6
      CALL REPARE(R0C,VC,V2C,R22C,D,R1C,R2C)
      CT1=R2C/R1C
      Z1C=D/(1.D0+CT1)
      Z2C=CT1*Z1C
      CALL P2C(R0C,R1C,R2C,Z1C,Z2C)
      XST1=XST
      XST2=XST
      FNST11=FNIU2(XST1)
      FNST22=FNST11
      DO 103 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2 
      FNST1=FNIU2(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
C     F1=FNIU2(XST11)*FNIU2(XST1)
C     PRINT*,'F1',F1,'   XXST1',XST1
            IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 7
           ENDIF
      FNST2=FNIU2(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
C     F1=FNIU2(XST22)*FNIU2(XST2)
C     PRINT*,'F1',F1,'   XST2',XST2
          IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 7
          ENDIF
103   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL, L=',L,'LL=',J
      RETURN
7     CONTINUE
      cucu1=1.e-5
      ncu=200
      CALL DRTMI(X,FF,FNIU2,ALI,ALS,cucu1,ncu,IER)
        IF(IER.EQ.2)PRINT*,'IER3=',IER,'IN NRADAC L DELTA',L,DELTA
      XST=X
6     CONTINUE
5     CONTINUE
      CALL P2C(R0C,R1,R2,Z1,Z2)
      XST1=XST
      XST2=XST
      FNST11=FNIU2(XST1)
      FNST22=FNST11
      DO 104 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU2(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
C     F1=FNIU2(XST11)*FNIU2(XST1)
C     PRINT*,'F1',F1
      IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 8
           ENDIF
      FNST2=FNIU2(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
C     F1=FNIU2(XST22)*FNIU2(XST2)
C     PRINT*,'F1',F1
          IF(F1.LT.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 8
          ENDIF
104   CONTINUE
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL, L=',L
      RETURN
8     CONTINUE
      cucu2=1.e-7
      ncu=200
      CALL DRTMI(X,FF,FNIU2,ALI,ALS,cucu2,ncu,IER)
      IF(IER.EQ.2)PRINT*,'IER4=',IER,'IN NRADAC L DELTA',L,DELTA
      X2=X
      X1=(W2*(X2+.5D0))/W1-.5D0
      RETURN
          ENDIF
      END


      SUBROUTINE REPARM(R1,R2,D,V,V2F,R2F,R0)
C ACEASTA SUBRUTINA PORNESTE DE LA PARAMETRII R1,R2 DELTA
C SI REPARAMETRIZEAZA R0,VOLUMUL TOTAL,R2F,V2F, IN
C PARAMETRIZAREA V2=CONSTANT (DOUA SFERE) PENTRU OBTINEREA
C SOLUTIILOR NIU IN POTENTIALUL CU DOUA CENTRE. S-A
C NORMALIZAT CU PI/3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(R1+R2.LE.D)GO TO 12
      R12=R1*R1
      R22=R2*R2
      D2=D*D
      D22=D+D
      C1=(R12+D2-R22)/D22 
      C2=(R22+D2-R12)/D22
      P1=R1+C1
      P11=P1*P1
      R11=R1+R1-C1
      P2=R2+C2
      P22=P2*P2
      R222=R2+R2-C2
      V2F=P22*R222
      V=P11*R11+V2F
      TH=1.D0/3.D0
      R0=(V/4.D0)**TH
      R2F=(V2F/4.D0)**TH
      RETURN
12    CONTINUE
      R222=R2*R2*R2 
      R222=R222+R222
      V2F=R222+R222
      R111=R1*R1*R1
      R111=R111+R111
      V1F=R111+R111
      V=V1F+V2F
      R0=(V/4.D0)**(.33333333333333333D0)
      R2F=R2
      RETURN
      END


      SUBROUTINE REPARE(R0,V,V2F,R2F,DELTA,R1,R2)
C ACEASTA SUBRUTINA PORNESTE DE LA PARAMETRII R0,V,V2F,R2C
C CALCULATI DE REPARM SI CALCULEAZA PE R1 SI R2 
C INTERMEDIARI PENTRU SOLUTIILE NIU (NORMALIZAT CU PI/3) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real cucu2
      COMMON/VREPC/R,V0,V2,D,R2D,R2A
      COMMON/RVREPC/IER2
      EXTERNAL VREP 
      R1FI=((V-V2F)/4.D0)**(1.D0/3.D0)
      IF(R1FI+R2F.LE.DELTA)THEN
      R1=R1FI 
      R2=R2F
      RETURN
      ENDIF
      R=R0
      V0=V
      V2=V2F
      D=DELTA
      R2D=R2F
      ALIMS=R0
      ALIMI=R1FI
      AI=(ALIMS+ALIMI)/2.D0
      AI1=(ALIMI+AI)/2.D0
      AI2=(ALIMS+AI)/2.D0
      NR1=1
      BLIMI=ALIMI
      BLIMS=AI1
55    CONTINUE
      cucu2=1.e-7
      ncu=200
      CALL DRTMI(X,FF,VREP,BLIMI,BLIMS,cucu2,ncu,IER)
      IF(IER.EQ.2.AND.NR1.EQ.1)THEN
      BLIMI=AI1
      BLIMS=AI
      NR1=2
      GO TO 55
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.2)THEN
      BLIMI=AI
      BLIMS=AI2
      NR1=3
      GO TO 55
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.3) THEN
      BLIMI=AI2
      BLIMS=ALIMS
      NR1=4
      GO TO 55
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.4) THEN
      BLIMI=ALIMS-.01D0
      BLIMS=ALIMS+.1D0
      NR1=5
      GO TO 55
      ENDIF
      ier22=ier2
      IF(IER2.NE.0)then
      PRINT*,'IER2',IER2,'IN VREP'
      endif
      IF(IER.NE.0)then
      PRINT*,'IER=',IER,'IN REPARE'
      endif
      IF(IER.NE.0.OR.IER22.NE.0)then
      PRINT*,'R0,V,V2F,DELTA,R,V0,V2,D,R2D,R2A,R2F'
      PRINT*,R0,V,V2F,DELTA,R,V0,V2,D,R2D,R2A,R2F
      endif
      R1=X
      R2=R2A
      RETURN
      END


      DOUBLE PRECISION FUNCTION VREP(X)
C FUNCTIE CARE CONDUCE LA CONSERVAREA VOLUMULUI TOTAL  
C NORMALIZATA CU PI/3  PARAMETRIZARE 2 SFERE CU V2F=CONSTANT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real cucu3
      COMMON/VREPC/R0,V0,V2F,DELTA,R2C,R2 
C     COMMON/IERVREPC/IER
      COMMON/RVREPC/IER
      COMMON/VREP2C/R1,R12,R22,D2,D22
      EXTERNAL VREP2
      R1=X
      ALIMS=R0
      ALIMI=R2C
      R12=R1*R1
      D2=DELTA*DELTA
      D22=DELTA+DELTA
      AI=(ALIMI+ALIMS)/2.D0
      AI1=(AI+ALIMI)/2.D0
      AI2=(AI+ALIMS)/2.D0
      NR1=1
      BLIMI=ALIMI
      BLIMS=AI1
77    CONTINUE
      cucu3=1.e-8
      ncu=200
      CALL DDRTMI(XX,FF,VREP2,BLIMI,BLIMS,cucu3,ncu,IER)
      IF(IER.EQ.2.AND.NR1.EQ.1) THEN
      NR1=2
      BLIMI=AI1
      BLIMS=AI
      GO TO 77
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.2) THEN
      NR1=3
      BLIMI=AI
      BLIMS=AI2
      GO TO 77
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.3) THEN
      NR1=4
      BLIMI=AI2
      BLIMS=ALIMS
      GO TO 77
      ENDIF
      IF(IER.EQ.2.AND.NR1.EQ.4) THEN
      NR1=5
      BLIMI=ALIMS-.01D0
      BLIMS=ALIMS+.1D0
      GO TO 77
      ENDIF
      R2=XX
      C1=(R12+D2-R22)/D22
      IF(C1.GE.R1) C1=R1
      IF(C1.LE.-R1) C1=-R1
      P1=R1+C1
      P11=P1*P1
      R11=R1+R1-C1
      VREP=V0-V2F-P11*R11 
      RETURN
      END


      DOUBLE PRECISION FUNCTION VREP2(X)
C ACEASTA FUNCTIE CONDUCE LA CONSERVAREA VOLUMULUI 2
C INPARAMETRIZAREA A DOUA SFERE INTERSECTATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VREPC/R0,V0,V2F,DELTA,R2C,R2
      COMMON/VREP2C/R1,R12,R22,D2,D22
      R22=X*X
      C2=(R22+D2-R12)/D22
      IF(C2.GE.X) C2=.9999999999D0*X
      IF(C2.LE.-X) C2=-.99999999999D0*X
      P2=X+C2
      P22=P2*P2
      R222=X+X-C2
      VREP2=V2F-P22*R222
      RETURN
      END


      DOUBLE PRECISION FUNCTION FNIU1(X)
C FUNCTIE CARE TREBUIE REZOLVATA PENTRU A OBTINE VALORILE 
C PROPRII        CAZ PAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/XN2/XN2
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      XN1=X 
      XN2=(W1*(XN1+.5D0))/W2-.5D0
      RZ =ZET1/ZET2
      IF( DABS(DABS(RZ)-1.D0).LT.2.D-7) THEN
      CALL DHERM(X,ZET1,H,D)
      CALL DHERM(X-1.D0,ZET1,H1,D1)
      FNIU1=ZET1*H+(X+X)*H1
      RETURN
      ELSE
      CALL DHERM(X,ZET1,H11,D11)
      CALL DHERM(X-1.D0,ZET1,H12,D12)
      CALL DHERM(XN2,ZET2,H21,D21)
      CALL DHERM(XN2-1.D0,ZET2,H22,D22)
      XN22=XN2+XN2
      XN11=XN1+XN1
      FNIU1=ALF2*H11*(ZET2*H21+XN22*H22)+ALF1*(ZET1*H11+XN11*H12)
     C *H21
      RETURN
      ENDIF
      END


      DOUBLE PRECISION FUNCTION FNIU2(X)
C FUNCTIE CARE TREBUIE REZOLVATA PENTRU A OBTINE VALORILE 
C PROPRII        CAZ IMPAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/XN1/XN1
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      XN2=X 
      XN1=(W2*(XN2+.5D0))/W1-.5D0
      RZ =ZET1/ZET2
      IF( DABS(DABS(RZ)-1.D0).LT.2.D-7) THEN
      CALL DHERM(X,ZET2,H1,D1)
      FNIU2=H1
      RETURN
      ELSE
      CALL DHERM(X,ZET2,H21,D21)
      CALL DHERM(X-1.D0,ZET2,H22,D22)
      CALL DHERM(XN1,ZET1,H11,D11)
      CALL DHERM(XN1-1.D0,ZET1,H12,D12)
      XN22=XN2+XN2
      XN11=XN1+XN1
      FNIU2=ALF2*H11*(ZET2*H21+XN22*H22)+ALF1*(ZET1*H11+XN11*H12)
     C *H21
      RETURN
      ENDIF
      END


      SUBROUTINE KNORM(ANIU1,ANIU2,C1,C2)
C     SE DETERMINA CONSTANTELE DE NORMARE C1 SI C2 PENTRU CAZ ASIMETRIC
C     INPUT VALORILE NIU1 SI NIU2 ( ANIU )
C     OUTPUT VALORILE C1 SI C2
C     PRIN COMMON SE OBTIN ZET1,ZET2,ALFA1,ALFA2 SI VALORILE FUNCTIILOR
C     HERMITE H1 SI H2 (PE RIND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/nrad/nnn
      common/nrfin/nrfin1,nrfin2
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/HPRN/H
c     print*,'zet1,zet2,nnn',zet1,zet2,nnn
      ZETA1=ZET1
      ZETA2=ZET2
c modificare 1 ->  0
          if((zeta2.gt.14.2d0).or.(zeta1.gt.3.2d0.and.nnn.eq.0))then
      call  dherm(aniu1,zet1,hr1,dr1)
      call  dherm(aniu2,zet2,hr2,dr2)
      expon1=dexp(-.5d0*zet1*zet1)
      expon2=dexp(-.5d0*zet2*zet2)
        if(nrfin1.ge.0)then
      nr=nrfin1
      call vect(nr,ve)
      c1=dsqrt(alf1/(2**nrfin1*ve*dsqrt(3.14159265358979324d0)))
      nrn=nnn
      nnn2=nrn/2
      nnnp=nnn2+nnn2
      ss=1.d0
      if(nnnp.ne.nrn)ss=-1.d0
c      c1=ss*dabs(c1)
      c2=c1*expon1/expon2*hr1/hr2
      c2=0.d0
        else
      nr=nrfin2
      call vect(nr,ve)
      c2=dsqrt(alf2/(2**nrfin2*ve*dsqrt(3.14159265358979324d0)))
      c1=c2*expon2/expon1*hr2/hr1
      c1=0.d0
        endif
           else
      call jcptn(aj1,aniu1,aniu1,zeta1)
      call dherm(aniu1,zeta1,h,dnr)
      H1=H
      call jcptn(aj2,aniu2,aniu2,zeta2)
      call dherm(aniu2,zeta2,h,dnr)
      H2=H
      E1=DEXP(-ZETA1*ZETA1*.5D0)*H1
      E2=DEXP(-ZETA2*ZETA2*.5D0)*H2
      E21=E2/E1
      E212=E21*E21
          if(nrfin2.lt.0)then
      nrn=nnn
      nnn2=nrn/2
      nnnp=nnn2+nnn2
      ss=1.d0
         if(nnnp.ne.nrn) ss=-1.d0
      c1=1.d0/dsqrt(dabs(aj1/alf1+aj2/alf2/e212))
c      c1=ss*dabs(c1)
      c2=c1/e21
           ss=1.
           if(c2.lt.0.d0)then
           ss=-1.
           c1=ss*c1
           c2=ss*c2
           endif
      else
      nrn=nnn
      nnn2=nrn/2
      nnnp=nnn2+nnn2
      ss=1.d0
         if(nnnp.ne.nrn) ss=-1.d0
      C2=1.D0/DSQRT(dabs(AJ2/ALF2+AJ1/ALF1*E212))
      C1=C2*E21
c      c1=ss*dabs(c2*e21)
           ss=1.
           if(c2.lt.0.d0)then
           ss=-1.
           c1=ss*c1
           c2=ss*c2
           endif
          endif
          endif
      RETURN
      END
 
      subroutine vect(lc,alc)
      implicit double precision( a)
c calculeaza alc=lc!
      if(lc.eq.0)then
      alc=1.d0
      else
      alc=1.d0
      do 1 i=1,lc
      alc=alc*i
1     continue
      endif
      return
      end





      SUBROUTINE I2N1N2(AI2,ANIU1,ANIU2,ZETA)
C FUNCTIA I(NIU1,2,NIU2,ZETA)=AI2 OUTPUT
C ANIU1,ANIU2,ZETA INPUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VPIURI/HV1,HV2,HV3,HV4,EXPO
      common/nrint/k1,k2
      CALL JCPTN(AJ1,ANIU1,ANIU2,ZETA)
      H1=HV1
      H2=HV2
      H3=HV3
      H4=HV4
      AJ=.5D0*(AJ1+EXPO*(H1*H3*(-ZETA)+
     +    ANIU1*H2*H3+ANIU2*H1*H4))
      if(k1.eq.0.or.k1.eq.1)then
      aj2=0.
      else
      CALL JCPTN(AJ2,ANIU1-2.D0,ANIU2,ZETA)
      endif
      ANIU21=ANIU2-1.D0
      ANIU11=ANIU1-1.D0
      if(k1.eq.0.or.k2.eq.0.or.k1.eq.1.or.k2.eq.1)then
      aj3=0.
      else
      CALL JCPTN(AJ3,ANIU11,ANIU21,ZETA)
      endif
      if(k2.eq.0.or.k2.eq.1)then
      aj5=0.
      else
      CALL JCPTN(AJ5,ANIU1,ANIU2-2.D0,ZETA)
      endif
      AJ4=AJ3
      AI2=AJ+ANIU1*(ANIU11*AJ2+ANIU2*AJ3)
     +     +ANIU2*(ANIU1*AJ4+ANIU21*AJ5)
      RETURN
      END


      SUBROUTINE I1N1N2(AI1,ANIU1,ANIU2,ZETA)
C FUNCTIA I(NIU1,1,NIU2,ZETA)=AI1 OUTPUT
C ANIU1, ANIU2, ZETA INPUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/nrint/k1,k2
      CALL JCPTN(AJ1,ANIU1,ANIU2+1,ZETA)
      if(k2.eq.0)then
      aj2=0
      else
      CALL JCPTN(AJ2,ANIU1,ANIU2-1.D0,ZETA)
      endif
      CALL DHERM(ANIU1,ZETA,H1,D1)
      CALL DHERM(ANIU2,ZETA,H2,D2)
      AI1=.5D0*AJ1+ANIU2*AJ2
      RETURN
      END

      SUBROUTINE P2C(R0I,R1,R2,Z1,Z2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/W12/W1,W2
      COMMON/ALF12/ALF1,ALF2
      common/r116/R0
      COMMON/CKST/CKST
c R1 si R2 sunt semi-axele A1 si A2
c     CKST=41.d0
c     R0=1.16d0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1=W0*R0/R1
      W2=W0*R0/R2
      ALF1=DSQRT(MPH*W1)*1.D-15
      ALF2=DSQRT(MPH*W2)*1.D-15
      ZET1=ALF1*Z1
      ZET2=ALF2*Z2
      RETURN
      END

      subroutine hermasimp(aniu,zeta,plus)
c reprezentarea asimptotice pentru functia Hermite zeta tinde la infinit
      implicit double precision (a-h,o-z)
c pentru zeta tinde la plus infinit
      plus=0.d0
      minus=1
      v=-aniu
      akf=1.d0
      zek=1.d0
      t=1.d0
      gam=1.d0
      k=0
1     continue
      plus=t+plus
      k=k+1
      ak2k=k+k-aniu-1
      ak2km=ak2k-1
      gam=ak2k*ak2km
      t=-t*gam/k/(zeta+zeta)**2
      alms=1.d-13*plus
      if(dabs(t).gt.dabs(alms))goto 1
      plus=plus*(2.d0*zeta)**(aniu)
      return
      end
     
      subroutine gama(z,r)
c gamma for x<0
      implicit double precision (a-h,o-z)
      n=0
      if(z.eq.0.d0)then
      r=1.d0
      return
      endif
      n=0
      if(z.lt.0.d0)then
        zz=z
1     continue
      if(zz.lt.0.d0)then
        n=n+1
        zz=zz+1.d0
        goto 1
      endif
      call gamma(zz,gam,ier)
      do i=1,n
        zz=zz-1.d0
        gam=gam/zz
      enddo
      r=gam
      return
      endif
      call gamma(z,r,ier)
      return
      end

 
      SUBROUTINE DHERM(A,Z,R1,R2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   COMPUTE THE HERMITE FUNCTION 'R1' AND ITS DERIVATIVE 'R2'
C
c compute the negative of z
      DIMENSION YHER(50)
      QPI=0.177245385090552D1
      TWO=0.69314718055994D0
      A1=-A*.5D0
      A2=A1+.5D0
      X=Z*Z
      AN=A+1.D-7
      NA=AN
      EA=AN-NA
! introduc si dezvoltarea asimptotica
      zpoz=-z ! deoarece se calculeaza pentru negativ lui z
      alsi=6.535d0!+dabs(a)
      alss=alsi+1.d0
       if(zpoz.ge.alsi)then
       call hermasimp(a,zpoz,plus)
c     print*,'plus 1 ',plus
       if(zpoz.ge.alss)then
       r1=plus
       return
       endif 
       endif

      IF(EA.LT.2.D-7.AND.A.GT.-1.D-7)GO TO 10
c     CALL F(A1,0.5D0,X,F1)
c     CALL F(A2,.15D1,X,F2)
      CALL FDUBLU(A1,0.5D0,X,F1,A2,.15D1,X,F2)
      CALL GAMAI(A1,G1)
      CALL GAMAI(A2,G2)
      CALL DHIPER(A1,.5D0,X,D1)
      CALL DHIPER(A2,.15D1,X,D2)
      H1=DERG(A1)
      H2=DERG(A2)
      P1=F1*G2+(Z+Z)*F2*G1
      P2=.2D1**A*QPI
      R1=P1*P2
       if(zpoz.ge.alsi.and.zpoz.le.alss)then
c     print*,'r1,plus,zpoz',r1,plus,zpoz
       r1=r1+(plus-r1)*(zpoz-alsi)
c     print*,'r1',r1
       endif 
      R2=P2*(TWO*P1-D1*G2*.5D0+F1*H2*.5D0+Z*(F2*H1-D2*G1))
      RETURN
10    CONTINUE
      EA2=AN/2.D0-NA/2
      IF(EA2.LT..25D0) THEN
      aniu=na
      eps=a-aniu
      es=1.d-7
      if(eps.lt.0.d0)es=-1.d-7
      a1=-aniu*0.5d0
      a2=a1+0.5d0
      CALL F(A1,.5D0,X,F1)
      CALL GAMAI(A2,G2)
      CALL DHIPER(A1,.5D0,X,D1)
      H2=DERG(A2)
      P1=F1*G2
      P2=.2D1**nA*QPI
      R1=P1*P2
      R2=P2*(TWO*P1-D1*G2*.5D0+F1*H2*.5D0)
      zm=-z
      call dhep(yher,zm,25)
      r1=yher(na+1)
      ann=aniu+es
      a1=-ann*0.5d0
      a2=a1+0.5d0
c     CALL F(A1,0.5D0,X,F1)
c     CALL F(A2,.15D1,X,F2)
      CALL FDUBLU(A1,0.5D0,X,F1,A2,.15D1,X,F2)
      CALL GAMAI(A1,G1)
      CALL GAMAI(A2,G2)
      CALL DHIPER(A1,.5D0,X,D1)
      CALL DHIPER(A2,.15D1,X,D2)
      H1=DERG(A1)
      H2=DERG(A2)
      P1=F1*G2+(Z+Z)*F2*G1
      P2=.2D1**Ann*QPI
      R11=P1*P2
      R22=P2*(TWO*P1-D1*G2*.5D0+F1*H2*.5D0+Z*(F2*H1-D2*G1))
      eps=dabs(eps)
      r1=r1+(r11-r1)*eps/1.d-7
      r2=r2+(r22-r2)*eps/1.d-7
      RETURN
      ELSE
      aniu=na
      eps=a-aniu
      es=1.d-7
      if(eps.lt.0.d0)es=-1.d-7
      a1=-aniu*0.5d0
      a2=a1+0.5d0
      CALL F(A2,.15D1,X,F2)
      CALL GAMAI(A1,G1)
      CALL DHIPER(A2,.15D1,X,D2)
      H1=DERG(A1)
      P1=(Z+Z)*F2*G1
      P2=.2D1**nA*QPI
      R1=P1*P2
      R2=P2*(TWO*P1+Z*(F2*H1-D2*G1))
      zm=-z
      call dhep(yher,zm,25)
      r1=yher(na+1)
      ann=aniu+es
      a1=-ann*0.5d0
      a2=a1+0.5d0
c     CALL F(A1,0.5D0,X,F1)
c     CALL F(A2,.15D1,X,F2)
      CALL FDUBLU(A1,0.5D0,X,F1,A2,.15D1,X,F2)
      CALL GAMAI(A1,G1)
      CALL GAMAI(A2,G2)
      CALL DHIPER(A1,.5D0,X,D1)
      CALL DHIPER(A2,.15D1,X,D2)
      H1=DERG(A1)
      H2=DERG(A2)
      P1=F1*G2+(Z+Z)*F2*G1
      P2=.2D1**Ann*QPI
      R11=P1*P2
      R22=P2*(TWO*P1-D1*G2*.5D0+F1*H2*.5D0+Z*(F2*H1-D2*G1))
      eps=dabs(eps)
      r1=r1+(r11-r1)*eps/1.d-7
      r2=r2+(r22-r2)*eps/1.d-7
      RETURN
      ENDIF
      END

      SUBROUTINE GAMAI(X,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C   COMPUTE THE RECIPROCAL 'R' OF THE GAMMA FUNCTION FOR THE ARGUMENT X
      IF(X.LE.0.D0) GO TO 1
      CALL GAMMA(X,DG1,IER)
      R=.1D1/DG1
      GO TO 2
1     K=DABS(X)
      P=X
      IF(K.EQ.0) GO TO 4
      DO 3 I=1,K
3     P=P*(X+I)
4     CALL GAMMA(X+(K+1),DG2,IER)
      R=P/DG2
2     RETURN
      END

      SUBROUTINE F(A,B,X,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C   COMPUTE THE HYPERGEOMETRIC FUNCTION
      EPS=1.D-13
      R=.1D1
      T=A*X/B
      AK=0.D0
1     R=R+T
      AK=AK+.1D1
      T=T*(A+AK)*X/((B+AK)*(AK+.1D1))
      IF(DABS(T).GT.(EPS*DABS(R))) GO TO 1
      RETURN
      END

      SUBROUTINE FDUBLU(A,B,X,R,AA1,BB1,XX1,RR1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C   COMPUTE THE HYPERGEOMETRIC FUNCTION
      EPS=1.D-13
      R=.1D1
      T=A*X/B
      AK=0.D0
      RR1=.1D1
      TT1=AA1*XX1/BB1
      AK1=0.D0
1     R=R+T
      RR1=RR1+TT1
      AK=AK+.1D1
      AK1=AK1+.1D1
      T=T*(A+AK)*X/((B+AK)*(AK+.1D1))
      TT1=TT1*(AA1+AK1)*XX1/((BB1+AK1)*(AK1+.1D1))
      IF(DABS(T).GT.(EPS*DABS(R))) GO TO 1
      RETURN
      END

      SUBROUTINE DHIPER(A,B,X,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C COMPUTE THE DERIVATIVE R=DF/DA OF THE HYPERGEOMETRIC FUNCTION
C F(A,B,X)
      EPS=.1D-12
      S=X/B
      F=A*S
      R=S
      AK=.1D1
1     C1=A+AK
      C2=(B+AK)*(AK+.1D1)
      S=X*(C1*S+F)/C2
      R=R+S
      IF(DABS(S).LE.(EPS*DABS(R))) GO TO 2
      F=C1*X*F/C2
      AK=AK+.1D1
      GO TO 1
2     RETURN
      END

      DOUBLE PRECISION FUNCTION DERG(X)
C COMPUTE THE RATIO PSI/GAMMA FUNCTION FOR ARGUMENT X
C BY SUBROUTINE DHERO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(X.LT.0.D0) GO TO 1
      P=X
      S=.1D1
      Y=X+.1D1
      GO TO 2
1     K=DABS(X)
      K=K+2
      Y=X+K
      P=.1D1
      DO 3 I=1,K
3     P=P*(X+(I-1))
      S=0.D0
      DO 4 I=1,K
      R=.1D1
      DO 5 J=1,K
      IF(J.NE.I) GO TO 6
      Z=.1D1
      GO TO 7
6     Z=X+(J-1)
7     R=R*Z
5     CONTINUE
4     S=S+R
2     CALL GAMMA(Y,DG1,IER)
      DERG=(P*PSIDR(Y)-S)/DG1
      RETURN
      END


      DOUBLE PRECISION FUNCTION PSIDR(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  COMPUTE THE PSI FUNCTION
      IF(X.GT.0.2D2) GO TO 5
      I=X
      Y=X-I
      R=0.D0
      IF(X.GE.0.3D1) GO TO 1
      K=3-I
      DO 3 J=1,K
3     R=R-.1D1/(X+(J-1))
      GO TO 2
1     IF(X.LT.0.4D1) GO TO 2
      K=I-1
      DO 4 J=3,K
4     R=R+.1D1/(Y+J)
2     PSIDR=.9227843350984667D0+Y*(.3949340668482262D0+Y*(-.77056903
     *15959405D-1+Y*(.1982323371113089D-1+Y*(-.5677755143195022D-2
     1+Y*(.1718061982054804D-2+Y*(-.5367773607419286D-3+Y*(.171106
     20679388364D-3+Y*(-.5526724702502909D-4+Y*(.180107006063411D-4
     3+Y*(-.5902470118113279D-5+Y*(.1936396575047679D-5+Y*(-.628631
     40427852799D-6+Y*(.19611156217856D-6+Y*(-.5539299852288D-7+Y*
     5(.128312147968D-7+Y*(-.208305913856D-8+Y*.17179869184D-9))))
     6))))))))))))
      PSIDR=PSIDR+R
      GO TO 6
5     Z1=.1D1/X
      Z2=Z1*Z1
      PSIDR=DLOG(X)-.5D0*Z1-.83333333333333333D-1*Z2*(.1D1+Z2*(-.1D0
     *+Z2*(.4761904761904762D-1+Z2*(-.5D-1+Z2*.9090909090909091D-1
     1))))
6     RETURN
      END


C        SUBROUTINE GAMMA
C
C        PURPOSE
C           COMPUTES THE DOUBLE PRECISION OF THE
C           GAMMA FUNCTION OF A GIVEN DOUBLE PRECISION ARGUMENT.
C
C        USAGE
C           CALL GAMMA(XX,DLNG,IER)
C
C        DESCRIPTION OF PARAMETERS
C           XX   - THE DOUBLE PRECISION ARGUMENT FOR THE GAMMA
C                  FUNCTION.
C           DLNG - THE RESULTANT DOUBLE PRECISION GAMMA FUNCTION
C                  VALUE.
C           IER  - RESULTANT ERROR CODE WHERE
C                  IER= 0----NO ERROR.
C                  IER=-1----XX IS WITHIN 10**(-9) OF BEING ZERO OR XX
C                            IS NEGATIVE.  LOG(GAMMA) IS SET TO -1.OD38.
C                  IER=+1----XX IS GREATER THAN 10**35. LOG(GAMMA) IS SET TO
C                            +1.OD38.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           THE EULER-MCLAURIN EXPANSION TO THE SEVENTH DERIVATIVE TERM
C           IS USED, AS GIVEN BY M. ABRAMOWITZ AND I.A. STEGUN,
C           'HANDBOOK OF MATHEMATICAL FUNCTIONS', U. S. DEPARTMENT OF
C           COMMERCE, NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS
C           SERIES, 1966, EQUATION 6.1.41.
C
C     ..................................................................
C
      SUBROUTINE GAMMA(XX,DLNG,IER)
      DOUBLE PRECISION XX,ZZ,TERM,RZ2,DLNG
      IER=0
      ZZ=XX
      IF(XX-1.D10) 2,2,1
    1 IF(XX-1.D35) 8,9,9
C
C        SEE IF XX IS NEAR ZERO OR NEGATIVE
C
    2 IF(XX-1.D-9) 3,3,4
    3 IER=-1
      DLNG=-1.D38
      GO TO 10
C
C        XX GREATER THAN ZERO AND LESS THAN OR EQUAL TO 1.D+10
C
    4 TERM=1.D0
    5 IF(ZZ-18.D0) 6,6,7
    6 TERM=TERM*ZZ
      ZZ=ZZ+1.D0
      GO TO 5
    7 RZ2=1.D0/ZZ**2
      DLNG =(ZZ-0.5D0)*DLOG(ZZ)-ZZ +0.9189385332046727 -DLOG(TERM)+
     1(1.D0/ZZ)*(.8333333333333333D-1 -(RZ2*(.2777777777777777D-2 +(RZ2*
     2(.7936507936507936D-3 -(RZ2*(.5952380952380952D-3)))))))
      DLNG =DEXP(DLNG)
      GO TO 10
C
C        XX GREATER THAN 1.D+10 AND LESS THAN 1.D+35
C
    8 DLNG=ZZ*(DLOG(ZZ)-1.D0)
      DLNG=DEXP(DLNG)
      GO TO 10
C
C        XX GREATER THAN OR EQUAL TO 1.D+35
C
    9 IER=+1
      DLNG=1.D38
   10 RETURN
      END
C
C     ..................................................................
C
C        SUBROUTINE EIGEN
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C        matrix
C        USAGE
C           CALL EIGEN(A,R,N,MV)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN DESCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C                   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C                   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C                       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C                       SEQUENCE)
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN 'MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS', EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
      SUBROUTINE sEIGEN(A,R,N,MV)
      DIMENSION A(1),R(1)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
      DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
     1                 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
    5 RANGE=1.0D-12
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0
      IF(I-J) 20,15,20
   15 R(IJ)=1.0
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM) 165,165,40
   40 ANORM=1.414*dSQRT(ANORM)
      ANRMX=ANORM*RANGE/FLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/FLOAT(N)
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
   62 IF(dABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5*(A(LL)-A(MM))
   68 Y=-A(LM)/dSQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 SINX=Y/dSQRT(2.0*(1.0+(dSQRT(1.0-Y*Y))))
      SINX2=SINX*SINX
   78 COSX=dSQRT(1.0-SINX2)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(LL)-A(MM)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
      RETURN
      END
C
C     ..................................................................
C
C        SUBROUTINE DQG12
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)
C
C        USAGE
C           CALL DQG12 (XL,XU,FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 12-POINT GAUSS QUADRATURE
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 23
C           EXACTLY. FOR REFERENCE, SEE
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.
C
C     ..................................................................
C
      SUBROUTINE DQG12(XL,XU,FCT,Y)
      EXTERNAL FCT
C
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49078031712335963D0*B
      Y=.23587668193255914D-1*(FCT(A+C)+FCT(A-C))
      C=.45205862818523743D0*B
      Y=Y+.53469662997659215D-1*(FCT(A+C)+FCT(A-C))
      C=.38495133709715234D0*B
      Y=Y+.8003916427167311D-1*(FCT(A+C)+FCT(A-C))
      C=.29365897714330872D0*B
      Y=Y+.10158371336153296D0*(FCT(A+C)+FCT(A-C))
      C=.18391574949909010D0*B
      Y=Y+.11674626826917740D0*(FCT(A+C)+FCT(A-C))
      C=.62616704255734458D-1*B
      Y=B*(Y+.12457352290670139D0*(FCT(A+C)+FCT(A-C)))
      RETURN
      END

      SUBROUTINE DDRTMI(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C SUBRUTINA DE REZOLVARE X=SOLUTIA,FCT=FUNCTIA,XLI SI XRI=LIMITELE
C INTERVALULUI PE CARE SE CAUTA SOLUTII
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
      EXTERNAL FCT
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
1     FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
2     FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCES FOR FUNCTIONS VALUES.
C
3     I=0
      TOLF=100.D0*EPS
C
C     START ITERATION LOOP
C
4     I=I+1
C
C     START BISECTION LOOP
C
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
5     IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
C
6     TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
7     TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
8     IF(I-IEND)17,17,9
9     XR=X
      FR=F
C
C     TEST ON SATISFATORY ACCURACY IN BISECTION LOOP
C
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
10    TOL=TOL*A
11    IF(DABS(XR-XL)-TOL)12,12,13
12    IF(DABS(FR-FL)-TOLF)14,14,13
13    CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS.ERROR RETURN.
C
      IER=1
14    IF(DABS(FR)-DABS(FL))16,16,15
15    X=XL
      F=FL
16    RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
C
17    A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFATORY ACCURACY ON ITERATION LOOP
18    TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
19    TOL=TOL*A
20    IF(DABS(DX)-TOL)21,21,22
21    IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
C
22    IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
23    XR=X
      FR=F
      GO TO 4
24    XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C
C     END OF ITERATION LOOP
C     
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
C
25    IER=2
      RETURN
      END

      SUBROUTINE DDDRTM(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C SUBRUTINA DE REZOLVARE X=SOLUTIA,FCT=FUNCTIA,XLI SI XRI=LIMITELE
C INTERVALULUI PE CARE SE CAUTA SOLUTII
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
      EXTERNAL FCT
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
1     FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
2     FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCES FOR FUNCTIONS VALUES.
C
3     I=0
      TOLF=100.D0*EPS
C
C     START ITERATION LOOP
C
4     I=I+1
C
C     START BISECTION LOOP
C
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
5     IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
C
6     TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
7     TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
8     IF(I-IEND)17,17,9
9     XR=X
      FR=F
C
C     TEST ON SATISFATORY ACCURACY IN BISECTION LOOP
C
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
10    TOL=TOL*A
11    IF(DABS(XR-XL)-TOL)12,12,13
12    IF(DABS(FR-FL)-TOLF)14,14,13
13    CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS.ERROR RETURN.
C
      IER=1
14    IF(DABS(FR)-DABS(FL))16,16,15
15    X=XL
      F=FL
16    RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
C
17    A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFATORY ACCURACY ON ITERATION LOOP
18    TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
19    TOL=TOL*A
20    IF(DABS(DX)-TOL)21,21,22
21    IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
C
22    IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
23    XR=X
      FR=F
      GO TO 4
24    XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C
C     END OF ITERATION LOOP
C     
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
C
25    IER=2
      RETURN
      END

      DOUBLE PRECISION FUNCTION FIN11(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALF12/ALF1,ALF2 
      COMMON/W12/W1,W2
      COMMON/FORV/FORV(18)
      COMMON/ZET1S2/Z1,Z2,A1,A2
      RO3=FORV(13)
      R3=FORV(4) 
      R1=FORV(5)
      ZC1=FORV(8)-Z1
      Z3=FORV(12)-Z1
      ZC2=FORV(9)-Z1
      SR3=FORV(3)
      ALFP=ALF1*ALF1
      A=Z-ALF1*(Z1+Z3) 
      AA=DSQRT(ALFP*R3*R3-A*A)
      BB=ALF1*RO3-SR3*AA
      CALL DHERM(A1,Z,H1,D1)
      CALL DHERM(A2,Z,H2,D2)
      FIN11=-((BB*sr3*A/AA)+z)*DEXP(-Z*Z)*H1*H2
      RETURN
      END





      DOUBLE PRECISION FUNCTION FIN1(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALF12/ALF1,ALF2 
      COMMON/W12/W1,W2
      COMMON/FORV/FORV(18)
      COMMON/ZET1S2/Z1,Z2,A1,A2
      RO3=FORV(13)
      R3=FORV(4) 
      R1=FORV(5)
      ZC1=FORV(8)-Z1
      Z3=FORV(12)-Z1
      ZC2=FORV(9)-Z1
      SR3=FORV(3)
      ALFP=ALF1*ALF1
      A=Z-ALF1*(Z1+Z3) 
      AA=DSQRT(ALFP*R3*R3-A*A)
      BB=ALF1*RO3-SR3*AA
      BBP=BB*BB
      CC=Z-ALF1*(Z1+ZC1)
      CCP=CC
      ZZZ=ZC2-ZC1
      W2PW1=W2/W1
      EE=1.D0+(W2PW1*W2PW1-1.D0)*CCP/(ALF1*ZZZ)
      CALL DHERM(A1,Z,H1,D1)
      CALL DHERM(A2,Z,H2,D2)
      FIN1=EE*BBP*DEXP(-Z*Z)*H1*H2
      RETURN
      END



      DOUBLE PRECISION FUNCTION FIN2(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALF12/ALF1,ALF2
      COMMON/W12/W1,W2
      COMMON/FORV/FORV(18)
      COMMON/ZET1S2/Z1,Z2,A1,A2
      RO3=FORV(13)
      R3=FORV(4)
      R1=FORV(5)
      ZC1=FORV(8)-Z1
      ZC2=FORV(9)-Z1
      Z3=FORV(12)-Z1
      SR3=FORV(3)
      ALFP=ALF2*ALF2
      A=Z+ALF2*(Z2-Z3)
      AA=DSQRT(ALFP*R3*R3-A*A)
      BB=ALF2*RO3-SR3*AA
      BBP=BB*BB
      CC=Z+ALF2*(Z2-ZC1)
      CCP=CC
      ZZZ=ZC2-ZC1
C      print*,'w2,w1,ccp,alf2,zzz'
C      print*,w2,w1,ccp,alf2,zzz
      w2pw1=w2/w1
      EE=1.D0+(W2PW1*W2PW1-1.D0)*CCP/(ALF2*ZZZ)
      ZS=-Z
      CALL DHERM(A1,ZS,H1,D1)
      CALL DHERM(A2,ZS,H2,D2)
      FIN2=EE*BBP*DEXP(-Z*Z)*H1*H2
      RETURN
      END

      DOUBLE PRECISION FUNCTION FIN3(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALF12/ALF1,ALF2
      COMMON/W12/W1,W2
      COMMON/FORV/FORV(18)
      COMMON/ZET1S2/Z1,Z2,A1,A2
      RO3=FORV(13)
      R3=FORV(4)
      R1=FORV(5)
      ZC1=FORV(8)-Z1
      ZC2=FORV(9)-Z1
      Z3=FORV(12)-Z1
      SR3=FORV(3)
      ALFP=ALF2*ALF2
      A=Z+ALF2*(Z2-Z3)
      AA=DSQRT(ALFP*R3*R3-A*A)
      BB=ALF2*RO3-SR3*AA
      BBP=BB*BB
      CC=Z+ALF2*(Z2-ZC2)
      CCP=-CC
      ZZZ=ZC2-ZC1
      W1PW2=W1/W2
      EE=1.D0+(W1PW2*W1PW2-1.D0)*CCP/(ALF2*ZZZ)
      ZS=-Z
      CALL DHERM(A1,ZS,H1,D1)
      CALL DHERM(A2,ZS,H2,D2)
      FIN3=EE*BBP*DEXP(-Z*Z)*H1*H2
      RETURN
      END

      DOUBLE PRECISION FUNCTION FIN4(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INFINT/AAA,BBB,ZZCC11,ZZCC22
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET1S2/Z1,Z2,A1,A2
      COMMON/W12/W1,W2
      D1=Z+ALF2*(Z2-ZZCC11)
      D1=D1/ALF2
      D=AAA*D1+BBB
      W2W1=W2/W1
      WP=W2W1*W2W1
      D3=1.D0+(WP-1.D0)*D1/(ZZCC22-ZZCC11)
      D2=D*D
      D2=D3*D2
      ZS=-Z
      CALL DHERM(A1,ZS,H1,D1)
      CALL DHERM(A2,ZS,H2,DF)
      FIN4=D2*H1*H2*DEXP(-Z*Z)
      RETURN
      END
     

      DOUBLE PRECISION FUNCTION FIN5(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INFINT/AAA,BBB,ZZCC11,ZZCC22
      COMMON/ALF12/ALF1,ALF2
      COMMON/W12/W1,W2
      COMMON/ZET1S2/Z1,Z2,A1,A2
      D1=Z+ALF2*(Z2-ZZCC11)
      D1=D1/ALF2
      D=AAA*D1+BBB
      W1W2=W1/W2
      WP=W1W2*W1W2
      D4=-Z+ALF2*(ZZCC22-Z2)
      D4=D4/ALF2
      D3=1.D0+(WP-1.D0)*D4/(ZZCC22-ZZCC11)
      D2=D*D
      D2=D3*D2
      ZS=-Z
      CALL DHERM(A1,ZS,H1,D1)
      CALL DHERM(A2,ZS,H2,DF)
      FIN5=D2*H1*H2*DEXP(-Z*Z)
      RETURN
      END
     


      DOUBLE PRECISION FUNCTION FIN6(Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INFINT/AAA,BBB,ZZCC11,ZZCC22
      COMMON/ALF12/ALF1,ALF2
      COMMON/ZET1S2/Z1,Z2,A1,A2
      COMMON/W12/W1,W2
      D1=Z-ALF1*(Z1+ZZCC11)
      D1=D1/ALF1
      D=AAA*D1+BBB
      W2W1=W2/W1
      WP=W2W1*W2W1
      D3=1.D0+(WP-1.D0)*D1/(ZZCC22-ZZCC11)
      D2=D*D
      D2=D3*D2
      ZS=Z
      CALL DHERM(A1,ZS,H1,D1)
      CALL DHERM(A2,ZS,H2,DF)
      FIN6=D2*H1*H2*DEXP(-Z*Z)
      RETURN
      END
     
	SUBROUTINE IGAUSS(FCT,XL,XU,Y)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
	DOUBLE PRECISION XL,XU,Y,A,B,C,FCT

        EXTERNAL FCT
	A=.5D0*(XU+XL)
	B=XU-XL
	C=.49863193092474078D0*B
	Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
	C=.49280575577263417D0*B
	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	C=.48238112779375322D0*B
	Y=.12696032654631030D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.46745303796886984D0*B
	Y=.17136931456510717D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.44816057788302606D0*B
	Y=.21417949011113340D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.42468380686628499D0*B
	Y=.25499029631188088D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.39724189798397120D0*B
	Y=.29342046739267774D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.36609105937014484D0*B
	Y=.32911111388180923D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.33152213346510760D0*B
	Y=.36172897054424253D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.29385787862038116D0*B
	Y=.39096947893535153D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.25344995446611470D0*B
	Y=.41655962113473378D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.21067563806531767D0*B
	Y=.43826046502201906D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.16593430114106382D0*B
	Y=.45586939347881942D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.11964368112606854D0*B
	Y=.46922199540402283D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.7223598079139825D-1*B
	Y=.47819360039637430D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.24153832843869158D-1*B
	Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
	RETURN
	END

	SUBROUTINE IGAUSS1(FCT,XL,XU,Y)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
	DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
        double precision hvpl1,hvpl2,hvmi1,hvmi2
      
        dimension hvpl1(0:24,16),hvmi1(0:24,16)
        dimension hvpl2(0:24,16),hvmi2(0:24,16)
        common/hval1/hvpl1,hvmi1
        common/hval2/hvpl2,hvmi2
        common/nznnll/nznn,nzll
      common/lagnod/nrnod,nrnodzz,nodzz
        EXTERNAL FCT
	A=.5D0*(XU+XL)
	B=XU-XL
        n=nznn
        l=nzll

	C=.49863193092474078D0*B
        nodzz=1
	Y=.35093050047350483D-2*FCT(A+C)*hvpl1(n,1)*hvpl1(l,1)
        nodzz=2
	Y=.35093050047350483D-2*FCT(A-C)*hvmi1(n,1)*hvmi1(l,1)+Y

	C=.49280575577263417D0*B
        nodzz=3
	Y=.8137197365452835D-2*FCT(A+C)*hvpl1(n,2)*hvpl1(l,2)+Y
        nodzz=4
	Y=.8137197365452835D-2*FCT(A-C)*hvmi1(n,2)*hvmi1(l,2)+Y

	C=.48238112779375322D0*B
        nodzz=5
	Y=.12696032654631030D-1*FCT(A+C)*hvpl1(n,3)*hvpl1(l,3)+Y
        nodzz=6
	Y=.12696032654631030D-1*FCT(A-C)*hvmi1(n,3)*hvmi1(l,3)+Y
        
	C=.46745303796886984D0*B
        nodzz=7
	Y=.17136931456510717D-1*FCT(A+C)*hvpl1(n,4)*hvpl1(l,4)+Y
        nodzz=8
	Y=.17136931456510717D-1*FCT(A-C)*hvmi1(n,4)*hvmi1(l,4)+Y
        
	C=.44816057788302606D0*B
        nodzz=9
	Y=.21417949011113340D-1*FCT(A+C)*hvpl1(n,5)*hvpl1(l,5)+Y
        nodzz=10
	Y=.21417949011113340D-1*FCT(A-C)*hvmi1(n,5)*hvmi1(l,5)+Y
      
	C=.42468380686628499D0*B
        nodzz=11
	Y=.25499029631188088D-1*FCT(A+C)*hvpl1(n,6)*hvpl1(l,6)+Y
        nodzz=12
	Y=.25499029631188088D-1*FCT(A-C)*hvmi1(n,6)*hvmi1(l,6)+Y
        
	C=.39724189798397120D0*B
        nodzz=13
	Y=.29342046739267774D-1*FCT(A+C)*hvpl1(n,7)*hvpl1(l,7)+Y
        nodzz=14
	Y=.29342046739267774D-1*FCT(A-C)*hvmi1(n,7)*hvmi1(l,7)+Y
        
	C=.36609105937014484D0*B
        nodzz=15
	Y=.32911111388180923D-1*FCT(A+C)*hvpl1(n,8)*hvpl1(l,8)+Y
        nodzz=16
	Y=.32911111388180923D-1*FCT(A-C)*hvmi1(n,8)*hvmi1(l,8)+Y
       
	C=.33152213346510760D0*B
        nodzz=17
	Y=.36172897054424253D-1*FCT(A+C)*hvpl1(n,9)*hvpl1(l,9)+Y
        nodzz=18
	Y=.36172897054424253D-1*FCT(A-C)*hvmi1(n,9)*hvmi1(l,9)+Y
        
	C=.29385787862038116D0*B
        nodzz=19
	Y=.39096947893535153D-1*FCT(A+C)*hvpl1(n,10)*hvpl1(l,10)+Y
        nodzz=20
	Y=.39096947893535153D-1*FCT(A-C)*hvmi1(n,10)*hvmi1(l,10)+Y
        
	C=.25344995446611470D0*B
        nodzz=21
	Y=.41655962113473378D-1*FCT(A+C)*hvpl1(n,11)*hvpl1(l,11)+Y
        nodzz=22
	Y=.41655962113473378D-1*FCT(A-C)*hvmi1(n,11)*hvmi1(l,11)+Y
        
	C=.21067563806531767D0*B
        nodzz=23
	Y=.43826046502201906D-1*FCT(A+C)*hvpl1(n,12)*hvpl1(l,12)+Y
        nodzz=24
	Y=.43826046502201906D-1*FCT(A-C)*hvmi1(n,12)*hvmi1(l,12)+Y
       
	C=.16593430114106382D0*B
        nodzz=25
	Y=.45586939347881942D-1*FCT(A+C)*hvpl1(n,13)*hvpl1(l,13)+Y
        nodzz=26
	Y=.45586939347881942D-1*FCT(A-C)*hvmi1(n,13)*hvmi1(l,13)+Y
        
	C=.11964368112606854D0*B
        nodzz=27
	Y=.46922199540402283D-1*FCT(A+C)*hvpl1(n,14)*hvpl1(l,14)+Y
        nodzz=28
	Y=.46922199540402283D-1*FCT(A-C)*hvmi1(n,14)*hvmi1(l,14)+Y
        
	C=.7223598079139825D-1*B
        nodzz=29
	Y=.47819360039637430D-1*FCT(A+C)*hvpl1(n,15)*hvpl1(l,15)+Y
        nodzz=30
	Y=.47819360039637430D-1*FCT(A-C)*hvmi1(n,15)*hvmi1(l,15)+Y
       
	C=.24153832843869158D-1*B
        nodzz=31
	Y=Y+.48270044257363900D-1*FCT(A+C)*hvpl1(n,16)*hvpl1(l,16)
        nodzz=32
	Y=B*(Y+.48270044257363900D-1*FCT(A-C)*hvmi1(n,16)*hvmi1(l,16))
	RETURN
	END

	SUBROUTINE IGAUSS2(FCT,XL,XU,Y)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
	DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
        double precision hvpl1,hvpl2,hvmi1,hvmi2
      
        dimension hvpl1(0:24,16),hvmi1(0:24,16)
        dimension hvpl2(0:24,16),hvmi2(0:24,16)
        common/hval1/hvpl1,hvmi1
        common/hval2/hvpl2,hvmi2
        common/nznnll/nznn,nzll
      common/lagnod/nrnod,nrnodzz,nodzz
        
        EXTERNAL FCT
        n=nznn
        l=nzll
	A=.5D0*(XU+XL)
	B=XU-XL
      

	C=.49863193092474078D0*B
        nodzz=1
	Y=.35093050047350483D-2*FCT(A+C)*hvpl2(n,1)*hvpl2(l,1)
        nodzz=2
	Y=.35093050047350483D-2*FCT(A-C)*hvmi2(n,1)*hvmi2(l,1)+Y

	C=.49280575577263417D0*B
        nodzz=3
	Y=.8137197365452835D-2*FCT(A+C)*hvpl2(n,2)*hvpl2(l,2)+Y
        nodzz=4
	Y=.8137197365452835D-2*FCT(A-C)*hvmi2(n,2)*hvmi2(l,2)+Y

	C=.48238112779375322D0*B
        nodzz=5
	Y=.12696032654631030D-1*FCT(A+C)*hvpl2(n,3)*hvpl2(l,3)+Y
        nodzz=6
	Y=.12696032654631030D-1*FCT(A-C)*hvmi2(n,3)*hvmi2(l,3)+Y
        
	C=.46745303796886984D0*B
        nodzz=7
	Y=.17136931456510717D-1*FCT(A+C)*hvpl2(n,4)*hvpl2(l,4)+Y
        nodzz=8
	Y=.17136931456510717D-1*FCT(A-C)*hvmi2(n,4)*hvmi2(l,4)+Y
        
	C=.44816057788302606D0*B
        nodzz=9
	Y=.21417949011113340D-1*FCT(A+C)*hvpl2(n,5)*hvpl2(l,5)+Y
        nodzz=10
	Y=.21417949011113340D-1*FCT(A-C)*hvmi2(n,5)*hvmi2(l,5)+Y
      
	C=.42468380686628499D0*B
        nodzz=11
	Y=.25499029631188088D-1*FCT(A+C)*hvpl2(n,6)*hvpl2(l,6)+Y
        nodzz=12
	Y=.25499029631188088D-1*FCT(A-C)*hvmi2(n,6)*hvmi2(l,6)+Y
        
	C=.39724189798397120D0*B
        nodzz=13
	Y=.29342046739267774D-1*FCT(A+C)*hvpl2(n,7)*hvpl2(l,7)+Y
        nodzz=14
	Y=.29342046739267774D-1*FCT(A-C)*hvmi2(n,7)*hvmi2(l,7)+Y
        
	C=.36609105937014484D0*B
        nodzz=15
	Y=.32911111388180923D-1*FCT(A+C)*hvpl2(n,8)*hvpl2(l,8)+Y
        nodzz=16
	Y=.32911111388180923D-1*FCT(A-C)*hvmi2(n,8)*hvmi2(l,8)+Y
       
	C=.33152213346510760D0*B
        nodzz=17
	Y=.36172897054424253D-1*FCT(A+C)*hvpl2(n,9)*hvpl2(l,9)+Y
        nodzz=18
	Y=.36172897054424253D-1*FCT(A-C)*hvmi2(n,9)*hvmi2(l,9)+Y
        
	C=.29385787862038116D0*B
        nodzz=19
	Y=.39096947893535153D-1*FCT(A+C)*hvpl2(n,10)*hvpl2(l,10)+Y
        nodzz=20
	Y=.39096947893535153D-1*FCT(A-C)*hvmi2(n,10)*hvmi2(l,10)+Y
        
	C=.25344995446611470D0*B
        nodzz=21
	Y=.41655962113473378D-1*FCT(A+C)*hvpl2(n,11)*hvpl2(l,11)+Y
        nodzz=22
	Y=.41655962113473378D-1*FCT(A-C)*hvmi2(n,11)*hvmi2(l,11)+Y
        
	C=.21067563806531767D0*B
        nodzz=23
	Y=.43826046502201906D-1*FCT(A+C)*hvpl2(n,12)*hvpl2(l,12)+Y
        nodzz=24
	Y=.43826046502201906D-1*FCT(A-C)*hvmi2(n,12)*hvmi2(l,12)+Y
       
	C=.16593430114106382D0*B
        nodzz=25
	Y=.45586939347881942D-1*FCT(A+C)*hvpl2(n,13)*hvpl2(l,13)+Y
        nodzz=26
	Y=.45586939347881942D-1*FCT(A-C)*hvmi2(n,13)*hvmi2(l,13)+Y
        
	C=.11964368112606854D0*B
        nodzz=27
	Y=.46922199540402283D-1*FCT(A+C)*hvpl2(n,14)*hvpl2(l,14)+Y
        nodzz=28
	Y=.46922199540402283D-1*FCT(A-C)*hvmi2(n,14)*hvmi2(l,14)+Y
        
	C=.7223598079139825D-1*B
        nodzz=29
	Y=.47819360039637430D-1*FCT(A+C)*hvpl2(n,15)*hvpl2(l,15)+Y
        nodzz=30
	Y=.47819360039637430D-1*FCT(A-C)*hvmi2(n,15)*hvmi2(l,15)+Y
       
	C=.24153832843869158D-1*B
        nodzz=31
	Y=Y+.48270044257363900D-1*FCT(A+C)*hvpl2(n,16)*hvpl2(l,16)
        nodzz=32
	Y=B*(Y+.48270044257363900D-1*FCT(A-C)*hvmi2(n,16)*hvmi2(l,16))
	RETURN
	END


        subroutine herval1(nm,x1d,xl,xu)
c calculez valori functii hermite pe cuadraturi
        implicit double precision (a-h,o-z)
        dimension hvpl1(0:24,16),hvmi1(0:24,16)
        dimension x1d(0:24)
        common/hval1/hvpl1,hvmi1
        COMMON/ZET12/ZET1,ZET2
	A=.5D0*(XU+XL)
	B=XU-XL
             do k=0,nm
        a1=x1d(k)
               
	C=.49863193092474078D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,1)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,1)=h1

          
	C=.49280575577263417D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,2)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,2)=h1

          
	C=.48238112779375322D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,3)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,3)=h1

	C=.46745303796886984D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,4)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,4)=h1

	C=.44816057788302606D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,5)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,5)=h1

	C=.42468380686628499D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,6)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,6)=h1

	C=.39724189798397120D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,7)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,7)=h1

	C=.36609105937014484D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,8)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,8)=h1

	C=.33152213346510760D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,9)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,9)=h1

	C=.29385787862038116D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,10)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,10)=h1
	
	C=.25344995446611470D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,11)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,11)=h1
	
	C=.21067563806531767D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,12)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,12)=h1
	
	C=.16593430114106382D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,13)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,13)=h1
	
	C=.11964368112606854D0*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,14)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,14)=h1
	
	C=.7223598079139825D-1*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,15)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,15)=h1
	
	C=.24153832843869158D-1*B
        z=(a+c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvpl1(k,16)=h1
        z=(a-c)
        zz=z+dabs(zet1)
        call dherm(a1,zz,h1,d1)
        hvmi1(k,16)=h1
	enddo
        return
        end

        subroutine herval2(nm,x2d,xl,xu)
c calculez valori functii hermite pentru cuadraturi
        implicit double precision (a-h,o-z)
        dimension hvpl2(0:24,16),hvmi2(0:24,16)
        dimension x2d(0:24)
        common/hval2/hvpl2,hvmi2
        COMMON/ZET12/ZET1,ZET2
	A=.5D0*(XU+XL)
	B=XU-XL
             do k=0,nm
        a1=x2d(k)         
	C=.49863193092474078D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,1)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,1)=h1

	C=.49280575577263417D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,2)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,2)=h1

	C=.48238112779375322D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,3)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,3)=h1

	C=.46745303796886984D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,4)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,4)=h1

	C=.44816057788302606D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,5)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,5)=h1

	C=.42468380686628499D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,6)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,6)=h1

	C=.39724189798397120D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,7)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,7)=h1

	C=.36609105937014484D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,8)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,8)=h1

	C=.33152213346510760D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,9)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,9)=h1

	C=.29385787862038116D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,10)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,10)=h1
	
	C=.25344995446611470D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,11)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,11)=h1
	
	C=.21067563806531767D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,12)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,12)=h1
	
	C=.16593430114106382D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,13)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,13)=h1
	
	C=.11964368112606854D0*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,14)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,14)=h1
	
	C=.7223598079139825D-1*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,15)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,15)=h1
	
	C=.24153832843869158D-1*B
        z=(a+c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvpl2(k,16)=h1
        z=(a-c)
        zz=-(z-zet2)
        call dherm(a1,zz,h1,d1)
        hvmi2(k,16)=h1
	enddo
        return
        end

       subroutine ordo(r1,r2,nmax,n,m)
c n  la iesire este un numar care este negativ daca nivelul este localizat
c in a doua groapa dupa sciziune iar daca este localizat in prima groapa are 
c valoarea unui intreg care ne da pe nz
c m la iesire este negativ daca nivelul este localizat in a doua groapa
c PENTRU DETERMINARE AM NEVOIE DE VALORILE semiaxelor A1 si A2 finale
c dar le voi folosi pe cele momentane (din timpul deformarii) fiindca aceste
c valori ma intereseaza mai mult la sciziune unde am aceste valori bune
       implicit double precision (a-h,o-z)
       dimension x(0:49),n(0:49),m(0:49)
       w1=1.d0/r1
       w2=1.d0/r2
       do 6345 i=0,nmax-1
       x(i)=i
       x(i+nmax)=w2/w1*(i+.5)-.5
6345   continue
       nnn=nmax+nmax-1
       do 6346 i=0,nnn
       xx=x(i)
       do 6347 j=i+1,nnn
       if(xx.gt.x(j)) then
       aa=xx
       xx=x(j)
       x(j)=aa
       endif
6347   continue
       x(i)=xx
       n(i)=xx
6346   continue


            kkk=0
            do kk=0,nnn
            difa=dabs(x(kk)-dint(x(kk)))
            if(difa.lt.1.d-3)then
            n(kk)=x(kk)
            m(kk)=-1
            else
            n(kk)=-1
            m(kk)=kkk
            kkk=kkk+1
            endif
            enddo


       return
       end







C
C     ..................................................................
C
C        SUBROUTINE EIGEN
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C           MATRIX
C
C        USAGE
C           CALL EIGEN(A,R,N,MV)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN DESCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C                   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C                   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C                       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C                       SEQUENCE)
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN 'MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS', EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
      SUBROUTINE EIGEN(A,R,N,MV)
      DIMENSION A(10000000),R(10000000)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
      DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
     1                 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
    5 RANGE=1.0E-12
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0
      IF(I-J) 20,15,20
   15 R(IJ)=1.0
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM) 165,165,40
   40 ANORM=1.414*dSQRT(ANORM)
      ANRMX=ANORM*RANGE/FLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/FLOAT(N)
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
   62 IF( dABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5*(A(LL)-A(MM))
   68 Y=-A(LM)/ dSQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 SINX=Y/ dSQRT(2.0*(1.0+( dSQRT(1.0-Y*Y))))
      SINX2=SINX*SINX
   78 COSX= dSQRT(1.0-SINX2)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(LL)-A(MM)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
      RETURN
      END


C        SUBROUTINE EIGEN
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C        matrix
C        USAGE
C           CALL EIGEN(A,R,N,MV)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN DESCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C                   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C                   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C                       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C                       SEQUENCE)
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN 'MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS', EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
      SUBROUTINE DEIGEN(R,N,MV)
      DIMENSION R(36100)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
C    1                 COSX2,SINCS,RANGE,ewqter
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
      common/eig/a(5565)
    5 RANGE=1.0D-6
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0
      IF(I-J) 20,15,20
   15 R(IJ)=1.0
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM) 165,165,40
   40 ANORM=1.414*SQRT(ANORM)
      ANRMX=ANORM*RANGE/FLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/FLOAT(N)
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
   62 IF(ABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5*(A(LL)-A(MM))
   68 Y=-A(LM)/SQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 SINX=Y/SQRT(2.0*(1.0+(SQRT(1.0-Y*Y))))
      SINX2=SINX*SINX
   78 COSX=SQRT(1.0-SINX2)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(LL)-A(MM)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
      RETURN
      END


      SUBROUTINE JCPTNN(AJ,ANIU1,ANIU2,ZETA)
C FUNCTIA J(NIU1,NIU2,ZETA)
C SE DETERMINA J=AJ OUTPUT
C ANIU1,ANIU2,ZETA INPUT
      IMPLICIT double precision  (A-H,O-Z)
      double precision aj,aniu1,aniu2,zeta,h1,h2,h3,h4,
     c expo,zet1,zet2,aniu,
     c anniu,d1,d2,d3,d4,ga,a1,a2,aa11,aa22,g11,g22,psidr
      COMMON/VPIURI/H1,H2,H3,H4,EXPO
      common/zet12/zet1,zet2
      common/nrfin/nrfin1,nrfin2
      IF(ANIU1+1.D-7.GT.ANIU2.AND.ANIU1-1.D-7.LT.ANIU2)THEN
      ANIU=.5D0*(ANIU1+ANIU2)
      naniu=aniu
      anniu=naniu
      if(anniu.eq.aniu) ANiu=aniu+1.d-7
      CALL DHERM(ANIU,ZETA,H1,D1)
      CALL DHERM(ANIU-1.D0,ZETA,H2,D2)
      A1=-ANIU/2.D0
      A2=A1+.5D0
      semn=1.d0
      if(a1.lt.0.d0)semn=-1.d0
      if(dabs(a1).lt.1.d-8) a1=1.d-8*semn
      semn=1.d0
      if(a2.lt.0.d0)semn=-1.d0
      if(dabs(a2).lt.1.d-8) a2=1.d-8*semn
      na1=a1
      aa11=na1
      if(a1.eq.aa11) a1=a1+1.d-8
      na2=a2
      aa22=na2
      if(a2.eq.aa22) a2=a2+1.d-8
      CALL GAMAI(A1,G11)
      CALL GAMAI(A2,G22)
      G1=1.D0/G11
      G2=1.D0/G22
      P1=PSIDR(A1)
      P2=PSIDR(A2)
      Z2=-ZETA-ZETA
      E1=G1*P1
      E2=G2*P2
      EX1=E1
      EX2=E2
      ANF=0.D0
      ZLBNF=1.D0
1     CONTINUE
      ANF=ANF+1.D0
      ZLBNF=ZLBNF*Z2/ANF
      E1=E2
      FT1=-E1*ZLBNF
      EX1=EX1+FT1
      P1=P1+1.D0/A1
      G1=G1*A1
      A1=A1+1.D0
      E2=P1*G1
      EX2=EX2-E2*ZLBNF
      ANF=ANF+1.D0
      ZLBNF=ZLBNF*Z2/ANF
      E1=E2
      FT2=E1*ZLBNF
      EX1=EX1+FT2
      P2=P2+1.D0/A2
      G2=A2*G2
      A2=A2+1.D0
      E2=P2*G2
      EX2=EX2+E2*ZLBNF
      FT=FT1+FT2
      YER=G2*ZLBNF
      IF(ANF.LT.140)THEN
      IF(ANF.LT.125.D0)THEN
      IF(FT.GT.1.D-8.OR.ANF.LT.15.D0) GO TO 1
      ELSE
      AANNFF=1.D-8/DEXP(-ZETA*ZETA)*DEXP(-6.D0)
      IF(FT.GT.AANNFF) GO TO 1
      ENDIF
      ENDIF
      CALL GAMAI(-ANIU,GA)
      EXPO=DEXP(-ZETA*ZETA)
      AJ=EXPO*(ANIU*H2*EX1+H1*EX2)*GA/4.D0
      HHH=H1
      H3=HHH
      HHH=H2
      H4=HHH
      RETURN
      ENDIF
      CALL DHERM(ANIU1,ZETA,H1,D1)
      CALL DHERM(ANIU1-1.D0,ZETA,H2,D2)
      CALL DHERM(ANIU2,ZETA,H3,D3)
      CALL DHERM(ANIU2-1.D0,ZETA,H4,D4)
      EXPO=DEXP(-ZETA*ZETA)
      if(aniu1.eq.0.d0)then
      aj=expo*h4
      return
      endif
      if(aniu2.eq.0.d0)then
      aj=expo*h2
      return
      endif
      AJ=EXPO*(ANIU1*H2*H3-ANIU2*H1*H4)/(ANIU1-ANIU2)
      RETURN
      END


      DOUBLE PRECISION FUNCTION ERFU(X)
c functia eroare 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PP(2)
                  if(x.lt.-5.d0)then
                  erfu=-1.d0
                  return
                  endif
      IF(X.GT.5.D0) GOTO 100
      XX=X*X
      SS=X
      S=X
      N=0
      ZZZ=DABS(X-XX*X/3.D0)*1.D-10
1     CONTINUE
      N=N+1
      S=-S*(XX/N)
      SSS=S/(N+N+1)
      SS=SS+SSS
      PP(1)=PP(2)
      PP(2)=SSS
      PPP=DABS(PP(1)+PP(2))
      IF(PPP.GT.ZZZ.OR.N.LT.10) GO TO 1
      ERFU=(SS+SS)/DSQRT(3.14159265358979324D0)
      RETURN
100   CONTINUE
      P=.3275911D0
      A1=.254829592D0
      A2=-.284496736D0
      A3=1.421413741D0
      A4=-1.453152027D0
      A5=1.061405429D0
      T=1.D0/(1+P*X)
      TT=T*T
      T3=TT*T
      T4=TT*TT
      T5=TT*T3
      ERFU=1.D0-(A1*T+A2*TT+A3*T3+A4*T4+A5*T5)*DEXP(-X*X)
      RETURN
      END

      subroutine norma(zeta,n,m,r)
c aceasta  subrutina calculeazaintegrala Hn*Hm*exp(-z2)  de la
c  -zeta  la infinit pt. n si m >= 0
      implicit double precision  (a-h,o-z)
      dimension  xh(50)
      s=0.d0
      if(n.eq.m)  then

      if(n.eq.0)  then
      erfi=erfu(zeta)
      r=dsqrt(3.1415926535897324d0)/2.d0*(1.d0+erfi)
      return
      else
      erfi=erfu(zeta)
      call vect(m,rv)
      s=2**(m-1)*rv*(1.d0+erfi)*dsqrt(3.14159265358979324d0)
      call dhep(xh,-zeta,m+1)
      exp=dexp(-zeta*zeta)
      s1=exp*xh(m)*xh(m+1)
      aam=1.d0
      do 9173 i=1,m-1
      am=(m-i+1)
      aam=2.d0*am*aam
      s1=s1+aam*(exp*xh(m-i)*xh(m-i+1))
9173  continue
      r=s+s1
      return
      endif

      else

      if(n.eq.0.or.m.eq.0) then
      nm=m
      if(m.eq.0)nm=n
      call dhep(xh,-zeta,nm)
      exp=dexp(-zeta*zeta)
      r=(exp*xh(nm))
      return
      else

      if(m.gt.n)then
      nm=m
      mn=n
      else
      nm=n
      mn=m
      endif

      mmm=nm-mn
      call dhep(xh,-zeta,nm+1)
      exp=dexp(-zeta*zeta)
      aam=1.d0
      s1=0.d0
      do 6139 i=2,mn
      am=(nm-i+2)
      aam=2.d0*am*aam
      s1=s1+aam*(xh(mn-i+1)*xh(nm-i+2))
6139  continue
      aam=2.d0*(nm-mn+1)*aam
      s1=(s1+aam*xh(nm-mn)+xh(mn)*xh(nm+1))*exp
      r=s+s1
      return
      endif

      endif

      end

      subroutine jcptn(aj,aniu1,aniu2,zeta)
      implicit double precision (a-h,o-z)
      common/vpiuri/h1,h2,h3,h4,expo
c      print*,'aniu1,aniu2',aniu1,aniu2
C MODIFICARE
c      print*,'aniu1,aniu2,zeta',aniu1,aniu2,zeta
      if(dabs(aniu1).lt.1.d-8.and.dabs(aniu2).lt.1.d-8.and.zeta.gt.
     c  7.3d0)then
      call norma(zeta,0,0,rt)
      aj=rt
      return
      endif
C PINA AICI
      if(aniu1+1.d-7.gt.aniu2.and.aniu1-1.d-7.lt.aniu2)then
      aniu=(aniu1+aniu2)/2.d0
        if(aniu.gt.-2.d-7) then
        aniup=aniu+2.d-7
        naniu=aniup
        anan=naniu
        eps=aniu-anan
         if(dabs(eps).le.2.d-7) then
         call norma(zeta,naniu,naniu,r1)
         semn=1.d0
         if(eps.lt.0.d0)semn=-1.d0
         a1=anan+2.d-7*semn
         call jcptnn(r2,a1,a1,zeta)
         r=r1+(r2-r1)*dabs(eps)/2.d-7
         an1=aniu-1.d-7
         an2=aniu+1.d-7
         call jcptnn(r3,an1,an2,zeta)
         epsi=dabs(aniu-aniu1)
         aj=r+(r3-r)*epsi/1.d-7
         else
         call jcptnn(r1,aniu,aniu,zeta)
         an1=aniu-1.d-7
         an2=aniu+1.d-7
         call jcptnn(r2,an1,an2,zeta)
         epsi=dabs(aniu-aniu1)
         aj=r1+(r2-r1)*epsi/1.d-7
         endif
        else
        call jcptnn(r1,aniu,aniu,zeta)
        an1=aniu-1.d-7
        an2=aniu+1.d-7
        call jcptnn(r2,an1,an2,zeta)
        epsi=dabs(aniu-aniu1)
        aj=r1+(r2-r1)*epsi/1.d-7
        endif
      call dherm(aniu1,zeta,h1,d1)
      call dherm(aniu1-1.d0,zeta,h2,d2)
      call dherm(aniu2,zeta,h3,d3)
      call dherm(aniu2-1.d0,zeta,h4,d4)
      expo=dexp(-zeta*zeta)
      return
      else
      call jcptnn(aj,aniu1,aniu2,zeta)
      return
      endif
      end

      subroutine normr(zeta,an,m,r)
c aceasta  subrutina calculeazaintegrala Han*Hm*exp(-z2)  de la
c  -zeta  la infinit, an real, m intreg
      implicit double precision (a-h,o-z)
      dimension xh(50),yh(50)
      if(m.eq.0) then
      call dherm(an-1.d0,zeta,r1,r2)
      exp=dexp(-zeta*zeta)
      r=r1*exp
      return
      else
      n=an
      if(m.lt.n)then
      mn=n-m
      call vect(m,rv)
      call dhep(xh,zeta,m+1)
      call dhep(yh,0.d0,m+1)
      exp=dexp(-zeta*zeta)
      and=an-m-1
      call dherm(and,zeta,xhnm,r2e)
      s=2**(m)*rv*(-exp*xhnm)
      call dherm(an,zeta,xhp,r2e)
      s1=-exp*xh(mn)*xhp
      aam=1.d0
      do 5287 i=1,mn-2
      am=(mn-i+1)
      aam=2.d0*am*aam
      ani=an-i
      call dherm(ani,zeta,xhnm,r2e)
      s1=s1+aam*(-exp*xh(mn-i)*xhnm)
5287  continue
      r=s+s1
      return
      endif
      endif
      end










      SUBROUTINE NRADAC1(L,R1,DELTA,X1)
c simetric

C ACEASTA SUBRUTINA OBTINE RADACINA NIU PENTRU VALOAREA L
C LA DOUA SFERE INTERSECTATE DACA SE CUNOSC R1,R2,DELTA
C SOLUTIILE SINT X1 SI X2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real cucu1,cucu2
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      COMMON/ZET12/ZET1,ZET2
      common/nrfin/nrfin1,nrfin2
      common/ptrads/detai,xdr1i,xdr2i,ndri
      EXTERNAL FNIU11,FNIU22
      common/CKST/CKST
c     CKST=41.D0
      n=l
      VINT1=.007D0
      VINT2=.001D0
      XST=L
      N=L
      AL=L
      LP2=L/2
      EAL=AL/2.D0-LP2
      
      CALL REPARM1(R1,DELTA,VC,V2C,R22C,R0C)
c     ALFA0=1.071186277015679D0/DSQRT(R0C)*dsqrt(ckst/41.d0) 
      ALFA0=1.071186277015679D0/DSQRT(R0C) 
      DIF=1.D-1*ALFA0
      r2=r1
      CT=R2/R1
      Z1=DELTA/(1.D0+CT)
      Z2=CT*Z1
      Z1C1=0.D0
      Z2C2=0.D0
      if(ndri.ne.0)then
      deti=detai
      call repare1(r0c,vc,deti,r1in,r2in)
      cti=r2in/r1in
      z1in=deti/(1.d0+cti)
      llin=z1in/dif
      endif
      IF(EAL.LT..25D0) THEN
      if(ndri.ne.0)xst=xdr1i
      LL=Z1/DIF
      IF(LL.EQ.0) GO TO 2
      DO 1 J=1,LL
      Z1C1=Z1C1+DIF
      Z2C2=CT*Z1C1
      D=(Z1C1+Z2C2)
      if(ndri.ne.0.and.j.le.llin) go to 1
      CALL REPARE1 (R0C,VC,D,R1C,R2C)
      CT1=R2C/R1C
      Z1C=D/(1.D0+CT1)
      Z2C=CT1*Z1C
      CALL P2C1(R0C,R1C,Z1C)
      XST1=XST
      XST2=XST
      FNST11=FNIU11(XST1)
      FNST22=FNST11
      DO 101 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU11(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1 
      IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 3
      ENDIF
      FNST2=FNIU11(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
      IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 3
      ENDIF
101   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL, L=',L,' LL=',J
      RETURN
3     CONTINUE
      cucu1=1.e-5
      ncu=200
      CALL DRTMI(X,FF,FNIU11,ALI,ALS,cucu1,ncu,IER)
      IF(IER.EQ.2)PRINT*,'IER1=',IER,'IN NRADAC1 L D',L,DELTA
      XST=X
1     CONTINUE
2     CONTINUE
      CALL P2C1(R0C,R1,Z1) 
      XST1=XST
      XST2=XST
      FNST11=FNIU11(XST1)
      FNST22=FNST11 
      DO 102 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU11(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
      IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 4
      ENDIF
      FNST2=FNIU11(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
      IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 4
      ENDIF
102   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL nradac1, L=',L
      RETURN
4     CONTINUE
      cucu2=1.e-7
      ncu=200
      CALL DRTMI(X,FF,FNIU11,ALI,ALS,cucu2,ncu,IER)
      IF(IER.EQ.2)PRINT*,'IER2=',IER,'IN NRADAC1 L DELTA',L,DELTA
      X1=X
      X2=x
      RETURN
      ELSE
      if(ndri.ne.0)xst=xdr2i
      LL=Z2/DIF
      IF(LL.EQ.0) GO TO 5
      DO 6 J=1,LL
      Z2C2=Z2C2+DIF
      Z1C1=Z2C2/CT
      D=Z2C2+Z1C1
      if(ndri.ne.0.and.j.le.llin) goto6
      CALL REPARE1(R0C,VC,D,R1C,R2C)
      CT1=R2C/R1C
      Z1C=D/(1.D0+CT1)
      Z2C=CT1*Z1C
      CALL P2C1(R0C,R1C,Z1C)
      XST1=XST
      XST2=XST
      FNST11=FNIU22(XST1)
      FNST22=FNST11
      DO 103 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2 
      FNST1=FNIU22(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
      IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 7
      ENDIF
      FNST2=FNIU22(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
      IF(F1.LE.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 7
      ENDIF
103   CONTINUE 
      PRINT*,' NU AM GASIT SOLUTIE IN INTERVAL (nradac1), L=',L,'LL=',J
      RETURN
7     CONTINUE
      cucu1=1.e-5
      ncu=200
      CALL DRTMI(X,FF,FNIU22,ALI,ALS,cucu1,ncu,IER)
      IF(IER.EQ.2)PRINT*,'IER3=',IER,'IN NRADAC1 L DELTA',L,DELTA
      XST=X
6     CONTINUE
5     CONTINUE
      CALL P2C1(R0C,R1,Z1)
      XST1=XST
      XST2=XST
      FNST11=FNIU22(XST1)
      FNST22=FNST11
      DO 104 I=1,300
      XST11=XST1
      XST22=XST2
      XST1=XST1-VINT1
      XST2=XST2+VINT2
      FNST1=FNIU22(XST1)
      F1=FNST11*FNST1
      FNST11=FNST1
      IF(F1.LE.0.D0) THEN
      ALS=XST11
      ALI=XST1
      GO TO 8
      ENDIF
      FNST2=FNIU22(XST2)
      F1=FNST22*FNST2
      FNST22=FNST2
      IF(F1.LT.0.D0) THEN
      ALS=XST2
      ALI=XST22
      GO TO 8
      ENDIF
104   CONTINUE
      PRINT*,' NO SOLUTION FOUND IN RANGE, L=',L
      RETURN
8     CONTINUE
      CUCU2=1.e-7
      NCU=200
      CALL DRTMI(X,FF,FNIU22,ALI,ALS,CUCU2,NCU,IER)
      IF(IER.EQ.2)PRINT*,'IER4=',IER,'IN NRADAC1 L DELTA',L,DELTA
      X2=X
      X1=x
      RETURN
      ENDIF
      END


      SUBROUTINE REPARM1(R1,D,V,V2F,R2F,R0)
c simetric
C ACEASTA SUBRUTINA PORNESTE DE LA PARAMETRII R1,R2 DELTA
C SI REPARAMETRIZEAZA R0,VOLUMUL TOTAL,R2F,V2F, IN
C PARAMETRIZAREA V2=CONSTANT (DOUA SFERE) PENTRU OBTINEREA
C SOLUTIILOR NIU IN POTENTIALUL CU DOUA CENTRE. S-A
C NORMALIZAT CU PI/3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     print*,'in reparm1 r1,d,v,v2f,r2f,r0'
c     print*,r1,d,v,v2f,r2f,r0
      r2=r1
      IF(R1+R2.LE.D)GO TO 12
      R12=R1*R1
      D2=D*D
      D22=D+D
      C1=d/2.d0
      P1=R1+C1
      P11=P1*P1
      R11=R1+R1-C1
      V=2.d0*P11*R11
      TH=1.D0/3.D0
      R0=(V/4.D0)**TH
      R2F=(V/8.D0)**TH
      RETURN
12    CONTINUE
      R111=R1*R1*R1
      R111=R111+R111
      V1F=R111+R111
      V=V1F+V1F
      R0=(V/4.D0)**(.33333333333333333D0)
      R2F=R1
      RETURN
      END


      DOUBLE PRECISION FUNCTION FNIU11(X)
C SIMETRIC
C FUNCTIE CARE TREBUIE REZOLVATA PENTRU A OBTINE VALORILE 
C PROPRII        CAZ PAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/XN2/XN2
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      CALL DHERM(X,ZET1,H,D)
      CALL DHERM(X-1.D0,ZET1,H1,D1)
      FNIU11=ZET1*H+(X+X)*H1
      RETURN
      END


      DOUBLE PRECISION FUNCTION FNIU22(X)
C SIMETRIC
C FUNCTIE CARE TREBUIE REZOLVATA PENTRU A OBTINE VALORILE 
C PROPRII        CAZ IMPAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/XN1/XN1
      COMMON/NRAD/N
      COMMON/W12/W1,W2
      CALL DHERM(X,ZET2,H1,D1)
      FNIU22=H1
      RETURN
      END

      SUBROUTINE P2C1(R0,R1,Z1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ZET12/ZET1,ZET2
      COMMON/W12/W1,W2
      COMMON/ALF12/ALF1,ALF2
      common/r116/R_0
      common/CKST/CKST
c     CKST=41.D0
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst *R_0/R0/H
      W1=W0*R0/R1
      W2=w1
      ALF1=DSQRT(MPH*W1)*1.D-15
      ALF2=alf1
      ZET1=ALF1*Z1
      ZET2=zet1
      RETURN
      END


      SUBROUTINE REPARE1(R0,V,DELTA,R1,R2)
c simetric
C ACEASTA SUBRUTINA PORNESTE DE LA PARAMETRII R0,V,
C CALCULATI DE REPARM SI CALCULEAZA PE R1 SI R2 
C INTERMEDIARI PENTRU SOLUTIILE NIU (NORMALIZAT CU PI/3) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL CUCU2
      COMMON/VREPC/R,V0,V2,D,R2D,R2A
      COMMON/RVREPC/IER2
      EXTERNAL VREP1 
c     print*,'in repare1 r0,v,delta,r1,r2,r,v0,v2,d,r2d,r2a'
c     print*,r0,v,delta,r1,r2,r,v0,v2,d,r2d,r2a
      R1FI=(V/8.D0)**(1.D0/3.D0)
      R2F=R1FI
      IF(R1FI+R1fi.LE.DELTA)THEN
      R1=R1FI 
      R2=r1
      RETURN
      ENDIF
      R=R0
      V0=V
      V2=V0/2.d0
      D=DELTA
      R2D=R2F
      ALIMS=R0
      ALIMI=R1FI
      CUCU2=1.E-7
      NCU=200
      CALL DRTMI(X,FF,VREP1,ALIMI,ALIMS,CUCU2,NCU,IER)
      IF(IER2.NE.0)PRINT*,'IER2',IER2,'IN VREP1'
      IF(IER.NE.0)PRINT*,'IER=',IER,'IN REPARE1'
      IF(IER.NE.0.OR.IER2.NE.0) THEN
      PRINT*,'R0,V,DELTA,R,V0,V2,D,R2D,R2A,R2F'
      PRINT*,R0,V,DELTA,R,V0,V2,D,R2D,R2A,R2F
      ENDIF
      R1=X
      R2=R1
      RETURN
      END


      DOUBLE PRECISION FUNCTION VREP1(X)
c simetric
C FUNCTIE CARE CONDUCE LA CONSERVAREA VOLUMULUI TOTAL  
C NORMALIZATA CU PI/3  PARAMETRIZARE 2 SFERE CU V2F=CONSTANT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VREPC/R0,V0,V2F,DELTA,R2C,R2 
      R1=X
      VREP1=V0/2.d0-(r1+delta/2.)**2*(2*r1-delta/2.)
      RETURN
      END


      SUBROUTINE KNORM1(ANIU1,ANIU2,C1,C2)
C     SE DETERMINA CONSTANTELE DE NORMARE C1 SI C2 PENTRU CAZ ASIMETRIC
C     INPUT VALORILE NIU1 SI NIU2 ( ANIU )
C     OUTPUT VALORILE C1 SI C2
C     PRIN COMMON SE OBTIN ZET1,ZET2,ALFA1,ALFA2 SI VALORILE FUNCTIILOR
C     HERMITE H1 SI H2 (PE RIND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/nrad/nnn
      common/nrfin/nrfin1,nrfin2
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      COMMON/HPRN/H
      ZETA1=ZET1
      ZETA2=ZET2
      zeta1=zet1
      call jcptn(aj1,aniu1,aniu1,zeta1)
      cn=dsqrt(dabs(alf1/aj1/2.))
      c2=cn
      c1=(-1)**(nnn)*cn
      RETURN
      END
 
      SUBROUTINE I3N1N2(AI3,ANIU1,ANIU2,ZETA)
C FUNCTIA I(NIU1,3,NIU2,ZETA)=AI3 OUTPUT
C ANIU1,ANIU2,ZETA INPUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL I2N1N2(AI21,ANIU1-1.D0,ANIU2,ZETA)
      CALL I2N1N2(AI22,ANIU1,ANIU2-1.D0,ZETA)
      CALL DHERM(ANIU1,ZETA,H1,D1)
      CALL DHERM(ANIU2,ZETA,H2,D2)
      CALL I1N1N2(AI1,ANIU1,ANIU2,ZETA)
      ZET2=ZETA*ZETA
      AI3=.5D0*ZET2*DEXP(-ZET2)*H1*H2+
     C    ANIU1*AI21+ANIU2*AI22+AI1
      RETURN
      END

      SUBROUTINE COSO(A,AKN,AKP,ANN,ANP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C CALCULEZ CONSTANTE SPIN ORBITA SI L**2
      IF(A.LE.25.D0)THEN
      AKP=.08D0
      AKN=.08D0
      ANP=0.D0
      ANN=0.D0
      RETURN
      ENDIF
      IF(A.LE.103.D0)THEN
      AKP=.08D0+(.0686D0-.08D0)*(A-25.D0)/(103.D0-25.D0)
      AKN=.08D0+(.0638D0-.08D0)*(A-25.D0)/(103.D0-25.D0)
      ANP=.560*(A-25.D0)/(103.D0-25.D0)
      ANN=.497*(A-25.D0)/(103.D0-25.D0)
      RETURN
      ENDIF
      IF(A.LE.110.D0)THEN
      AKP=.0686D0+(.070D0-.0686D0)*(A-103.D0)/(110.D0-103.D0)
      AKN=.0638D0+(.066D0-.0638D0)*(A-103.D0)/(110.D0-103.D0)
      ANP=.560D0+(.4D0-.560D0)*(A-103.D0)/(110.D0-103.D0)
      ANN=.497D0+(.35D0-.497D0)*(A-103.D0)/(110.D0-103.D0)
      RETURN
      ENDIF
      IF(A.LE.122.D0)THEN
      AKP=.070D0+(.0671D0-.070D0)*(A-110.D0)/(122.D0-110.D0)
      AKN=.066D0+(.0638D0-.066D0)*(A-110.D0)/(122.D0-110.D0)
      ANP=.40D0+(.572D0-.40D0)*(A-110.D0)/(122.D0-110.D0)
      ANN=.35D0+(.493D0-.35D0)*(A-110.D0)/(122.D0-110.D0)
      RETURN
      ENDIF
      IF(A.LE.140.D0)THEN
      AKP=.0671D0+(.0657D0-.0671D0)*(A-122.D0)/(140.D0-122.D0)
      AKN=.0638D0+(.0637D0-.0638D0)*(A-122.D0)/(140.D0-122.D0)
      ANP=.572D0+(.584D0-.572D0)*(A-122.D0)/(140.D0-122.D0)
      ANN=.493D0+(.451D0-.493D0)*(A-122.D0)/(140.D0-122.D0)
      RETURN
      ENDIF
      IF(A.LE.165.D0)THEN
      AKP=.0657D0+(.0637D0-.0657D0)*(A-140.D0)/(165.D0-140.D0)
      AKN=.0637D0+(.0637D0-.0637D0)*(A-140.D0)/(165.D0-140.D0)
      ANP=.584D0+(.6D0-.584D0)*(A-140.D0)/(165.D0-140.D0)
      ANN=.451D0+(.42D0-.451D0)*(A-140.D0)/(165.D0-140.D0)
      RETURN
      ENDIF
      IF(A.LE.242.D0)THEN
      AKP=.0637D0+(.0577D0-.0637D0)*(A-165.D0)/(242.D0-165.D0)
      AKN=.0637D0+(.0635D0-.0637D0)*(A-165.D0)/(242.D0-165.D0)
      ANP=.6D0+(.65D0-.6D0)*(A-165.D0)/(242.D0-165.D0)
      ANN=.42D0+(.325D0-.42D0)*(A-165.D0)/(242.D0-165.D0)
      RETURN
      ENDIF
      AKP=.0577D0
      AKN=.0635D0
      ANP=.65
      ANN=.325
      RETURN
      END


C     ..................................................................
C
C        SUBROUTINE DHEP
C
C        PURPOSE
C           COMPUTE THE VALUES OF THE HERMITE POLYNOMIALS H(N,X)
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.
C
C        USAGE
C           CALL DHEP(Y,X,N)
C
C        DESCRIPTION OF PARAMETERS
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUES
C                   OF HERMITE POLYNOMIALS OF ORDER 0 UP TO N
C                   FOR GIVEN ARGUMENT X.
C                   DOUBLE PRECISION VECTOR.
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER
C           X     - ARGUMENT OF HERMITE POLYNOMIAL
C                   DOUBLE PRECISION VARIABLE.
C           N     - ORDER OF HERMITE POLYNOMIAL
C
C        REMARKS
C           N LESS THAN 0 IS TREATED AS IF N WERE 0
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR
C           HERMITE POLYNOMIALS H(N,X)
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X))
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,
C           THE SECOND IS THE ARGUMENT.
C           STARTING VALUES ARE H(0,X)=1, H(1,X)=2*X.
C
C     ..................................................................
C
      SUBROUTINE DHEP(Y,X,N)
C
      DIMENSION Y(50)
      DOUBLE PRECISION Y,X,F
C
C        TEST OF ORDER
      Y(1)=1.D0
      IF(N)1,1,2
    1 RETURN
C
    2 Y(2)=X+X
      IF(N-1)1,1,3
C
    3 DO 4 I=2,N
      F=X*Y(I)-DFLOAT(I-1)*Y(I-1)
    4 Y(I+1)=F+F
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! POTENTIAL MACROSCOPIC !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MPE(A0,Z0,A2,Z2,EC0,EC,EY0,EN,EV,ED0,ED,ew,V)
      implicit double precision (a-h,o-z)
      common/lasciziune/delatsciz
      CALL GRULED(32)
      CALL COULOM(A0,Z0,A2,Z2,EC0,EC)
      print*,'a0,z0,a2,z2,ec0,ec'
      print*,a0,z0,a2,z2,ec0,ec

      CALL ENVOL(A0,Z0,A2,Z2,EV)
      print*,'EV',EV
      CALL ENYPE(A0,z0,A2,Z2,EY0,EN)
      print*,'EY0,EN',EY0,EN
      CALL DIFUZC(A0,Z0,A2,Z2,ED0,ED)
      Print*,'ED0,ED',ED0,ED
      call wignerplusa0(a0,z0,a2,z2,ew)
      print*,'ew',ew
      V=EC-EC0+EN-EY0+EV+ED-ED0+ew
      print*,'V',V
c     stop
      return
      end

      SUBROUTINE COULOM(A0,Z0,A2E,Z2E,EC0,EC)
C ACEASTA SUBRUTINA OBTINE ENERGIA COULOMBIANA A NUCLEULUI 
C NOTATA EC. DATELE DE INTRARE SINT FORV(18)=DIMENSIUNILE
C CARE PARAMETRIZEAZA FORMA, A0=NUMARUL DE MASA AL 
C FRAGMENTULUI PARINTE, Z0=NUMARUL ATOMIC AL PARINTELUI
C Z2=NUMARUL ATOMIC AL FRAGMENTULUI EMIS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/volumel/v00,v11,v22
      common/r116/R0
      COMMON/NGAUSS/NY
      common/lasciziune/delatsciz
      EXTERNAL FFCOUL
      PI=3.141592645d0
      Z1=Z0-Z2
      SCON=1.4399764D0 ! e^2 in MeV fm
      DELTA=DELT
      print*,'a1,delta,a2',a1,delta,a2
      ALUNG=(A1+DELTA+A2)/2.D0 ! SEMILUNGIMEA NUCLEU
      DIFER=C1-A1+ALUNG ! POZITIA MIJLOCULUI
      print*,'alung,difer',alung,difer
      CALL RENORM(ALUNG,DIFER)
      ! test
      akl=c1n-a1n
      bkl=c2n+a2n
      print*,'lim,c1n,x1n,x2n,c2n,lim',akl,c1n,x1n,x2n,c2n,bkl

      RR=R0*A0**(1.D0/3.D0) ! RAZA INITIALA
C DETERMINAM CONSTANTA CU CARE SE INMULTESC INTEGRALELE
      ACON=.625D0*(ALUNG/RR)**5/PI
      EC0=.6D0*SCON*Z0*Z0/RR ! 	ENERGIA COULOMB A SFEREI
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz!-2 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
c     DPAR07=delatsciz-4
      RO0=Z0/A0
      RO1F=(Z0-Z2E)/(A0-A2E)
      RO2F=Z2E/A2E
      print*,'delta,dpar,dpar07'
      print*,delta,dpar,dpar07
      print*,'ro0.ro1f,ro2f',ro0,ro1f,ro2f
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z1M=Z0-RO2*V22/V00*A0
      A1M=A0-V22/V00*A0
      print*,'z1m,a1m',z1m,a1m
      RO1=Z1M/A1M
      ENDIF
      ENDIF
      print*,'ro1,ro2,ro0',ro1,ro2,ro0
      CO1R=(RO1/RO0)**2
      CO2R=(RO2/RO0)**2
      CO3R=RO1*RO2/RO0**2
      print*,'co1r,co2r,co3r',co1r,co2r,co3r
C SE DEOSEBESC MAI MULTE CAZURI

      A1PA2PR3=A1+A2+2.D0*R3-1.D-6
c     a1pa2pr3=delatsciz
      IF(A1PA2PR3.LT.DELTA)THEN ! ELIPSE COMPLET SEPARATE
      CALL WGAUSS(FFCOUL,-1.D0,U1N,R1)
      CALL SGAUSS(FFCOUL,-1.D0,U1N,U2N,1.D0,RR3)
      CALL WGAUSS(FFCOUL,U2N,1.D0,R2)
      R1=ACON*R1*CO1R
      R2=ACON*R2*CO2R
      RR3=2.D0*ACON*RR3*CO3R
      EC=EC0*(R1+R2+RR3)
      print*,'caz ELIPSE COMPLET SEPARATE'
      RETURN
      ENDIF
      IF(S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3.AND.R3.LT.500.d0)THEN !C1<C3<C2
      IF(RO3.LT.R3)THEN ! U1 NU E EGAL CU U2
      CALL INTEGQ(FFCOUL,-1.D0,X1N,U1N,U2N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FFCOUL,-1.D0,X1N,U1N,NY,R1)
      CALL SINTQ(FFCOUL,U2N,X2N,1.D0,NY,R2)
      R1=ACON*R1*CO1R
      R2=ACON*R2*CO2R
      RR3=2.D0*ACON*RR3*CO3R
      print*,'r1,r2,rr3',r1,r2,rr3
      EC=EC0*(R1+R2+RR3)
      print*,'x1n,u1n,c3n,u2n,x2n'
      print*,x1n,u1n,c3n,u2n,x2n
      print*,'caz U1 NU E EGAL CU U2'
      RETURN
      ELSE ! C1<C3<C2 DAR CELE DOUA FRAGMENTE COMUNICA
      CALL INTEGQ(FFCOUL,-1.D0,X1N,C3N,C3N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FFCOUL,-1.D0,X1N,C3N,NY,R1)
      CALL SINTQ(FFCOUL,C3N,X2N,1.D0,NY,R2)
      R1=ACON*R1*CO1R
      R2=ACON*R2*CO2R
      RR3=2.D0*ACON*RR3*CO3R
      EC=EC0*(R1+R2+RR3)
      print*,'caz S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3'
      RETURN
      ENDIF
      ELSE
      ! VOLUMUL DIN MIJLOC APARTINE PARTII STANGI
      APINT=(1.d0+X2N)/2.D0
      CALL SINTQ(FFCOUL,-1.D0,X1N,X2N,NY,R1)
      CALL INTEGQ(FFCOUL,-1.D0,X1N,X2N,X2N,APINT,1.D0,NY,RR3)
      CALL WGAUSS(FFCOUL,X2N,1.D0,R2)
      print*,'acon,r1,r2,r3',acon,r1,r2,r33
      R1=ACON*R1*CO1R
      R2=ACON*R2*CO2R
      RR3=2.D0*ACON*RR3*CO3R
      print*,'x1n,x2n',x1n,x2n
      print*,'r1,r2,rr3',r1,r2,rr3
      EC=EC0*(R1+R2+RR3)
      print*,'caz VOLUMUL DIN MIJLOC APARTINE PARTII STANGI'
      RETURN
      ENDIF
      END


      DOUBLE PRECISION FUNCTION FFCOUL(Z1,Z2)
C REPREZINTA INTEGRANTUL NECESAR CALCULULUI ENERGIEI POTENTIALE
C A UNUI NUCLEU DEFORMAT. FORMULA ESTE DATA DE DAVIES SI
C SIERK, 1975, J. COMPUT. PHYS., 18, 311
C O VARIABILA ESTE DATA PRIN COMMON/FFCOUL/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
C CALCULAM VALOAREA FUNCTIEI IN PUNCTUL IN CARE APARE 
C NEDETERMINAREA Z1=Z2. PENTRU A CISTIGA TIMP INITIALIZAM
C RZ1=ROFZ(Z1), RZ2=ROFZ(Z2), R2DZ1=RO2DZ(Z1), R2DZ2=RO2DZ(Z2)
      Z1Z=Z1
      CALL VARPAR(Z1Z,R2DZ1,RZ1)
c     print*,'z1z,r2dz1,rz1',z1z,r2dz1,rz1
C INITIALIZAM RZ1RZ1, Z1Z2, RZ2RZ2 PENTRU ECONOMIE TIMP RULARE
      RZ1RZ1=RZ1*RZ1
      IF(Z1.NE.Z2) GO TO 1
      FFCOUL=RZ1RZ1*RZ1/0.75D0
      RETURN
1     Z2Z=Z2
      CALL VARPAR(Z2Z,R2DZ2,RZ2)
      RZ2RZ2=RZ2*RZ2
      Z1Z2=Z1-Z2
      Z1Z2LP=Z1Z2*Z1Z2
      RZ1Z2R=RZ1+RZ2
      A=RZ1Z2R*RZ1Z2R+Z1Z2LP
      AK2M=4.D0*RZ1*RZ2/A
      if(ak2m.gt.0.999999999999d0)ak2m=0.999999999999d0
      CALL DELKE(AK2M,AK,D)
      D=(AK-D)/AK2M
C     FORMULA DAVIES-SIERK ESTE
      FFCOUL=(RZ1*RZ2*(AK-2.D0*D)/3.D0*
     &     (2.D0*(RZ1RZ1+RZ2RZ2)-Z1Z2LP+
     &      1.5D0*Z1Z2*(R2DZ2-R2DZ1))+
     &      AK*(RZ1RZ1*RZ2RZ2/3.D0+
     &      (RZ1RZ1-.5D0*Z1Z2*R2DZ1)*
     &      (RZ2RZ2+.5D0*Z1Z2*R2DZ2)))/DSQRT(A)
      RETURN
      END


      
      SUBROUTINE VARPAR(Z,RO2DZ,ROFZ)
C ACEASTA SUBROUTINA NE DA DEPENDENTA VARIABILEI
C RO LA PATRAT ( dRO**2Z) FATA DE VARIABILA Z
C SI A VARIABILEI RO FATA DE Z ( DIMENSIUNILE
C SINT NORMALIZATE )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
c     common/limitesol/x101,x201,y101,y201
C ACEASTA FUNCTIE ARE DIFERITE FORME IN FUNCTIE DE INTERVALUL
C PE CARE ESTE DEFINITA. AVEM TREI INTERVALE [-1,ZC1], 
C [ZC1,ZC2], [ZC2,1]. IN INTERVALUL [ZC1,ZC2] AVEM
C SITUATIE SPECIALA CIND R3=FORV(4) (NENORMALIZAT) TINDE LA
C INFINIT
      IF (Z.LE.X1N) GO TO 200
      IF (Z.LT.X2N) GO TO 100
C FIE INTERVALUL [X2N,1]
      ZZ2=(Z-C2N)/A2N**2
      RO2DZ=B2N**2*(-ZZ2-ZZ2)
      ROFZ=B2N*DSQRT(1.D0-((Z-C2N)/A2N)**2)
      RETURN
C FIE INTERVALUL [-1,X1N]
200   ZZ1=(Z-C1N)/A1N**2
      RO2DZ=B1N**2*(-ZZ1-ZZ1)
      ROFZ=B1N*DSQRT(1.D0-((Z-C1N)/A1N)**2)
      RETURN
C FIE INTERVALUL [X1N,X2N], AVEM SI POSIBILITATEA CA R3.GT.1.D6
100   CONTINUE
      DABR3=DABS(R3) 
      IF( DABR3 .GT. 500.D0) GO TO 300
      ZZ3=Z-C3N
      EXI=S*DSQRT(R3N*R3N-ZZ3*ZZ3)
      RO2DZ=(-2.D0+2.D0*RHO3N/EXI)*ZZ3
      ROFZ=RHO3N-EXI
      RETURN
C CAZUL R3 FOARTE MARE
300   CONTINUE
      ROFZ=Y101N+(Y201N-Y101N)*(Z-X101N)/(X201N-X101N)
      RO2DZ=2.D0*ROFZ*(Y201N-Y101N)/(X201N-X101N)
      RETURN
      END

      SUBROUTINE DELKE(X,DK,DE)
C COMPLETE ELLIPTIC INTEGRALS DK AND DE OF ARGUMENT X=K**2
C CHEBYSHEV APPROXIMATION. SEE J.CODY, MATH.OF.COMPUTATION
C 19,(1965),105
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AV(8),BV(8),CV(8),DV(8)
      COMMON/REOIER/IERM
      DATA AV/.39684709020989782D-2,.10795990490591634D-1,
     &        .10589953620989358D-1,.75193867218083810D-2,
     &        .89266462945564662D-2,.14942029142282078D-1,
     &        .30885173001899709D-1,.96573590301742528D-1/
      DATA BV/.17216147097986521D-2,.92811603829686042D-2,
     &        .20690240005100840D-1,.29503729348688713D-1,
     &        .37335546682286029D-1,.48827155048118010D-1,
     &        .70312495459546608D-1,.12499999999764066D0/
      DATA CV/.43025377747931166D-2,.11785841008733935D-1,
     &        .11841925995501249D-1,.90355277375408818D-2,
     &        .11716766944657723D-1,.21836131405486897D-1,
     &        .56805223329308289D-1,.44314718058336814D0/
      DATA DV/.18645379184063363D-2,.10087958494375100D-1,
     &        .22660309891604122D-1,.32811069172721062D-1,
     &        .42672510126591752D-1,.58592707184265274D-1,
     &        .93749995116367067D-1,.24999999999746142D0/
      E=1.D0-X
      IF(E)1,2,3
1     PRINT *,'ARGUMENT IN DELKE OUT OF RANGE -1,+1'
      RETURN
2     DK=1.D35
      DE=1.D0
      PRINT *,'ARGUMENT IN DELKE = 1'
      IERM=2
      RETURN
3     EL=DLOG(1.D0/E)
      AA=E*.30072519903686484D-3
      BB=E*.66631752464607315D-4
      CC=E*.32519201550639042D-3
      DD=E*.72031696345715460D-4
      DO 5 I=1,8
      AA=E*(AA+AV(I))
      BB=E*(BB+BV(I))
      CC=E*(CC+CV(I))
5     DD=E*(DD+DV(I))
      DE=1.D0+CC+EL*DD
      DK=1.3862943611198906D0+AA+EL*(.5D0+BB)
      RETURN
      END

	SUBROUTINE WGAUSS(FCT,XL,XU,YY)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
C INTEGREAZA FUNCTII SIMETRICE IN DOUA ARGUMENTE INTEGRATE
C PE ACELASI INTERVAL
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION C(16),Y(16)
        EXTERNAL FCT
	A=.5D0*(XU+XL)
	B=XU-XL
	C(1)=.49863193092474078D0*B
	Y(1)=.35093050047350483D-2
	C(2)=.49280575577263417D0*B
	Y(2)=.8137197365452835D-2
	C(3)=.48238112779375322D0*B
	Y(3)=.12696032654631030D-1
	C(4)=.46745303796886984D0*B
	Y(4)=.17136931456510717D-1
	C(5)=.44816057788302606D0*B
	Y(5)=.21417949011113340D-1
	C(6)=.42468380686628499D0*B
	Y(6)=.25499029631188088D-1
	C(7)=.39724189798397120D0*B
	Y(7)=.29342046739267774D-1
	C(8)=.36609105937014484D0*B
	Y(8)=.32911111388180923D-1
	C(9)=.33152213346510760D0*B
	Y(9)=.36172897054424253D-1
	C(10)=.29385787862038116D0*B
	Y(10)=.39096947893535153D-1
	C(11)=.25344995446611470D0*B
	Y(11)=.41655962113473378D-1
	C(12)=.21067563806531767D0*B
	Y(12)=.43826046502201906D-1
	C(13)=.16593430114106382D0*B
	Y(13)=.45586939347881942D-1
	C(14)=.11964368112606854D0*B
	Y(14)=.46922199540402283D-1
	C(15)=.7223598079139825D-1*B
	Y(15)=.47819360039637430D-1
	C(16)=.24153832843869158D-1*B
        Y(16)=.482700442573639D-1
      YY=0.D0
      DO 1 I=1,15
      CI=C(I)
      AP=A+CI
      AM=A-CI
      YI=Y(I)
      YY=YI*YI*(FCT(AP,AP)+FCT(AM,AM)
     &             +FCT(AP,AM)+FCT(AM,AP))+YY
      DO 2 J=I+1,16
      CJ=C(J)
      BP=A+CJ
      BM=A-CJ
      YY=YY+2.D0*YI*Y(J)*(FCT(AP,BP)+FCT(AM,BM)
     &                     +FCT(AP,BM)+FCT(AM,BP))
2     CONTINUE
1     CONTINUE
      YY=YY+Y(16)*Y(16)*(FCT(A+C(16),A+C(16))+FCT(A-C(16),A-C(16))
     &                  +FCT(A+C(16),A-C(16))+FCT(A-C(16),A+C(16)))
      YY=B*B*YY
	RETURN
	END

      SUBROUTINE GRULED(N)
C ZERO AND WEIGHT OF LEGENDRE POLYNOMIALS FOR GAUSSIAN
C QUADRATURE IN (-1,1).AFTER P.J.DAVIS, P.RABINOWITZ,
C METHODS OF NUMERICAL INTEGRATION, ACADEMIC PRESS, 1975.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(32),W(32)
      COMMON/GRUCOM/X,W
      COMMON/NGAUSS/NNN
      NNN=N
      IF(N.LE.32) GO TO 2
      N=32
      PRINT *,'ATTEMPT TO SURPASS N=32 IN GRULED'
2     N1=N+1
      M=N1/2
      E1=N*N1
      DO 1 I=1,M
      T=(4*I-1)*3.141592653589793D0/(4*N+2)
      X0=(1.D0-(1.D0-1.D0/N)/(8.D0*N*N))*DCOS(T)
      PKM1=1.D0
      PK=X0
      DO 3 K=2,N
      T1=X0*PK
      PKP1=T1-PKM1-(T1-PKM1)/K+T1
      PKM1=PK
3     PK=PKP1
      DEN=1.D0-X0*X0
      D1=N*(PKM1-X0*PK)
      DPN=D1/DEN
      D2PN=(2.D0*X0*DPN-E1*PK)/DEN
      D3PN=(4.D0*X0*D2PN+(2.D0-E1)*DPN)/DEN
      D4PN=(6.D0*X0*D3PN+(6.D0-E1)*D2PN)/DEN
      U=PK/DPN
      V=D2PN/DPN
      HU=V+U*(V*V-U*D3PN/(3.D0*DPN))
      H=-U*(1.D0+.5D0*U*HU)
      PU=D2PN+H/3.D0*(D3PN+.25D0*H*D4PN)
      P=PK+H*(DPN+.5D0*H*PU)
      DP=DPN+H*(D2PN+.5D0*H*(D3PN+H*D4PN/3.D0))
      H=H-P/DP
      X(I)=X0+H
      FU=D2PN+.25D0*H*(D3PN+.2D0*H*D4PN)
      FX=D1-H*E1*(PK+.5D0*H*(DPN+H/3.D0*FU))
      W(I)=2.D0*(1.D0-X(I)*X(I))/(FX*FX)
      IC=N1-I
      X(IC)=-X(I)
      W(IC)=W(I)
1     CONTINUE
      IF(M+M.GT.N)X(M)=.0D0
      RETURN
      END

	SUBROUTINE SGAUSS(FCT,XL,XU,XL1,XL2,YY)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
C INTEGREAZA FUNCTII SIMETRICE IN DOUA ARGUMENTE INTEGRATE
C PE INTERVALE DIFERITE
	DOUBLE PRECISION XL,XU,Y,A,B,C,FCT,YY,XL1,XL2,A1,B1,C1
        DIMENSION C(16),Y(16),C1(16)
        EXTERNAL FCT
	A=.5D0*(XU+XL)
	B=XU-XL
        A1=.5D0*(XL2+XL1)
      	B1=-XL1+XL2
	C(1)=.49863193092474078D0*B
        C1(1)=.49863193092474078D0*B1
	Y(1)=.35093050047350483D-2
	C(2)=.49280575577263417D0*B
      	C1(2)=.49280575577263417D0*B1
	Y(2)=.8137197365452835D-2
	C(3)=.48238112779375322D0*B
      	C1(3)=.48238112779375322D0*B1
	Y(3)=.12696032654631030D-1
	C(4)=.46745303796886984D0*B
      	C1(4)=.46745303796886984D0*B1
	Y(4)=.17136931456510717D-1
	C(5)=.44816057788302606D0*B
      	C1(5)=.44816057788302606D0*B1
	Y(5)=.21417949011113340D-1
	C(6)=.42468380686628499D0*B
      	C1(6)=.42468380686628499D0*B1
	Y(6)=.25499029631188088D-1
	C(7)=.39724189798397120D0*B
      	C1(7)=.39724189798397120D0*B1
	Y(7)=.29342046739267774D-1
	C(8)=.36609105937014484D0*B
      	C1(8)=.36609105937014484D0*B1
	Y(8)=.32911111388180923D-1
	C(9)=.33152213346510760D0*B
      	C1(9)=.33152213346510760D0*B1
	Y(9)=.36172897054424253D-1
	C(10)=.29385787862038116D0*B
      	C1(10)=.29385787862038116D0*B1
	Y(10)=.39096947893535153D-1
	C(11)=.25344995446611470D0*B
      	C1(11)=.25344995446611470D0*B1
	Y(11)=.41655962113473378D-1
	C(12)=.21067563806531767D0*B
      	C1(12)=.21067563806531767D0*B1
	Y(12)=.43826046502201906D-1
	C(13)=.16593430114106382D0*B
      	C1(13)=.16593430114106382D0*B1
	Y(13)=.45586939347881942D-1
	C(14)=.11964368112606854D0*B
      	C1(14)=.11964368112606854D0*B1
	Y(14)=.46922199540402283D-1
	C(15)=.7223598079139825D-1*B
      	C1(15)=.7223598079139825D-1*B1
	Y(15)=.47819360039637430D-1
	C(16)=.24153832843869158D-1*B
      	C1(16)=.24153832843869158D-1*B1
        Y(16)=.482700442573639D-1
      YY=0.D0
      DO 1 I=1,16
      DO 2 J=1,16
      YY=YY+Y(I)*Y(J)*(FCT(A+C(I),A1+C1(J))+FCT(A-C(I),A1-C1(J))
     &                     +FCT(A+C(I),A1-C1(J))+FCT(A-C(I),A1+C1(J)))
2     CONTINUE
1     CONTINUE
      YY=B*B1*YY
	RETURN
	END

      SUBROUTINE INTEGQ(FCT,AL0,AL1,AL2,AL3,AL4,AL5,N,R3)
C INTEGREAZA FUNCTII SIMETRICE PRIN GAUSS LEGENDRE IN N PUNCTE
C PE INTERVALELE AL1-AL0,...,AL4-AL3 (N-INTREG)
C REZ1,REZ2,REZ3 SINT REZULTATELE
C SE CHEAMA INITIAL GRULED(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(32),W(32),XX(128),WW(128),C(128)
      COMMON/GRUCOM/X,W
      EXTERNAL FCT
      K=0
      DO 2 JJ=1,4
      DO 1 I=1,N
      K=(JJ-1)*N+I
      WW(K)=W(I)
      XX(K)=X(I)
1     CONTINUE
2     CONTINUE
      AL10=.5D0*(AL1-AL0)
      AL21=.5D0*(AL2-AL1)
      AL43=.5D0*(AL4-AL3)
      AL54=.5D0*(AL5-AL4)
      BL21=.5D0*(AL2+AL1)
      BL10=.5D0*(AL1+AL0)
      BL54=.5D0*(AL5+AL4)
      BL43=.5D0*(AL4+AL3)
      N2=2*N
      N3=3*N
      DO 3 L=1,N
      C(L)=XX(L)*AL10
      C(N+L)=XX(N+L)*AL21
      C(N2+L)=XX(N2+L)*AL43
      C(N3+L)=XX(N3+L)*AL54
3     CONTINUE
      E1=0.D0
      E2=0.D0
      E3=0.D0
      E4=0.D0
      DO 4 L1=1,N
      DO 5 L2=1,N
      WP=W(L1)*W(L2)
      E1=E1+WP*FCT(BL10+C(L1),BL43+C(N2+L2))
      E2=E2+WP*FCT(BL21+C(N+L1),BL43+C(N2+L2))
      E3=E3+WP*FCT(BL10+C(L1),BL54+C(N3+L2))
      E4=E4+WP*FCT(BL21+C(N+L1),BL54+C(N3+L2))
5     CONTINUE
4     CONTINUE
      E1=AL10*AL43*E1
      E2=AL21*AL43*E2
      E3=AL10*AL54*E3
      E4=AL21*AL54*E4
      R3=E1+E2+E3+E4      
      RETURN 
      END
            
      SUBROUTINE RENORM(ALUNG,DIFER)
C RENORMEAZA LUNGIMEA NUCLEULUI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/limitesol/x101,x201,y101,y201
      a1n=a1/alung
      b1n=b1/alung
      a2n=a2/alung
      b2n=b2/alung
      c1n=(c1-difer)/alung
      c2n=(c2-difer)/alung
      c3n=(c3-difer)/alung
      x1n=(x1-difer)/alung
      x2n=(x2-difer)/alung
      r3n=r3/alung
      u1n=(u1-difer)/alung
      u2n=(u2-difer)/alung
      rho3n=ro3/alung
      sn=s
      x101n=(x101-difer)/alung
      x201n=(x201-difer)/alung
      y101n=y101/alung
      y201n=y201/alung
      RETURN
      END
      
      SUBROUTINE SINTQ(FCT,AL0,AL1,AL2,N,R)
C INTEGREAZA O FUNCTIE SIMETRICA IN ARGUMENTE PE INTERVALUL AL0-AL2
C PE DOUA INTERVALE AL0-AL1 AL1-AL2 PRIN GAUSS LEGENDRE IN N 
C PUNCTE PE FIECARE SUBINTERVAL . SE CHEAMA CALL GRULED(N) INITIAL.
C REZULTATUL ESTE R. N=INTREG, IN REST TOT E IN D.P.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(32),W(32),XX(64),C(64)
      COMMON/GRUCOM/X,W
      EXTERNAL FCT
      K=0
      DO 2 JJ=1,2
      DO 1 I=1,N
      K=(JJ-1)*N+I
      XX(K)=X(I)
1     CONTINUE
2     CONTINUE
      AL10=.5D0*(AL1-AL0)
      AL21=.5D0*(AL2-AL1)
      BL10=.5D0*(AL1+AL0)
      BL21=.5D0*(AL1+AL2)
      DO 3 L=1,N
      C(L)=XX(L)*AL10
      C(L+N)=XX(L+N)*AL21
3     CONTINUE
      G1=0.D0
      G2=0.D0
      DO 4 K=1,N
      WP=W(K)*W(K)
      G1=G1+WP*FCT(BL10+C(K),BL10+C(K))
      G2=G2+WP*FCT(BL21+C(N+K),BL21+C(N+K))
      IF(K.EQ.N) GO TO 5
      DO 6 M=K+1,N
      WP=W(K)*W(M)
      G1=G1+2.D0*WP*FCT(BL10+C(K),BL10+C(M))
      G2=G2+2.D0*WP*FCT(BL21+C(N+K),BL21+C(M+N))
6     CONTINUE
5     CONTINUE
4     CONTINUE
      G1=AL10*AL10*G1
      G2=AL21*AL21*G2
      G3=0.D0
      DO 7 M1=1,N
      DO 8 M2=1,N
      G3=G3+W(M1)*W(M2)*FCT(BL10+C(M1),BL21+C(M2+N))
8     CONTINUE
7     CONTINUE
      G3=AL21*AL10*G3
      R=G1+G2+2.D0*G3
      RETURN
      END



      SUBROUTINE ENVOL(A0,Z0,A2E,Z2E,EV)
C ACEASTA SUBRUTINA CALCULEAZA ENERGIA NUCLEARA DE VOLUM EV
C PARAMETRII DE INTRARE SINT A: NR. DE MASA A PARINTELUI,
C A2:NU. DE MASA AL FRAGMENTULUI EMIS , Z: NR. ATOMIC AL
C PARINTELUI, Z2: NR. ATOMIC AL FRAGMENTTULUI EMIS
      IMPLICIT DOUBLE PRECISION (A-Z)
      common/paramel/a1,b1,a2,b2,c11,c22,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/volumel/v00,v11,v22
      common/lasciziune/delatsciz
      DELTA=DELT 
      C1=15.9937D0
      C2=1.927D0
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2 (LA FEL CA PENTRU
C ENERGIA COULOMB) PENTRU A AFLA IZOSPINII AI (SAU denumite si excesele de neutroni)
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz!-2 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
c     DPAR07=delatsciz-4
      RO0=Z0/A0
      RO1F=(Z0-Z2E)/(A0-A2E)
      RO2F=Z2E/A2E
      AI0=(A0-2*Z0)/A0
      AI2F=(A2E-2*Z2E)/A2E
      AI1F=((A0-A2E)-2*(Z0-Z2E))/(A0-A2E)
c      print*,'delta,dpar,dpar07'
c      print*,delta,dpar,dpar07
c      print*,'ro0.ro1f,ro2f',ro0,ro1f,ro2f
      A2M=V22/V00*A0
      A1M=A0-V22/V00*A0
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      AI2=AI2F
      AI1=AI1F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      AI2=AI0
      AI1=AI0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      A2M=V22/V00*A0
      Z2M=RO2*A2M
      Z1M=Z0-RO2*V22/V00*A0
      A1M=A0-V22/V00*A0
      print*,'z1m,a1m',z1m,a1m
      RO1=Z1M/A1M
      AI2=(A2M-2*Z2M)/A2M
      AI1=(A1M-2*Z1M)/A1M
      ENDIF
      ENDIF
c      print*,'in envol c1,c2,pa1,ai1,ai2'
c      print*,c1,c2,pa1,ai1,ai2
c      print*,'pa2,pz2',pa2,pz2
      EV=C1*C2*(A1M*AI1*AI1+A2M*AI2*AI2-A0*AI0*AI0)
      RETURN
      END


      subroutine wignerplusa0(a0,z0,a2e,z2e,ew)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c11,c22,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/volumel/v00,v11,v22
      common/lasciziune/delatsciz
      W=30
      ca0=2.615
      DELTA=DELT 
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2 (LA FEL CA PENTRU
C ENERGIA COULOMB) PENTRU A AFLA IZOSPINII AI (SAU denumite si excesele de neutroni)
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz+2 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      DPAR07=delatsciz-2
      RO0=Z0/A0
      RO1F=(Z0-Z2E)/(A0-A2E)
      RO2F=Z2E/A2E
      AI0=(A0-2*Z0)/A0
      AI2F=(A2E-2*Z2E)/A2E
      AI1F=((A0-A2E)-2*(Z0-Z2E))/(A0-A2E)
c      print*,'delta,dpar,dpar07'
c      print*,delta,dpar,dpar07
c      print*,'ro0.ro1f,ro2f',ro0,ro1f,ro2f
      A2M=V22/V00*A0
      A1M=A0-V22/V00*A0
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      AI2=AI2F
      AI1=AI1F
           ea0=ca0
           ew0=W*(AI1+AI2-AI0)
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      AI2=AI0
      AI1=AI0
            ea0=0
            ew0=0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      A2M=V22/V00*A0
      Z2M=RO2*A2M
      Z1M=Z0-RO2*V22/V00*A0
      A1M=A0-V22/V00*A0
      print*,'z1m,a1m',z1m,a1m
      RO1=Z1M/A1M
      AI2=(A2M-2*Z2M)/A2M
      AI1=(A1M-2*Z1M)/A1M
            ea0=ca0*(delta-dpar07)/(dpar-dpar07)
            ew0=W*(AI1+AI2-AI0)*(delta-dpar07)/(dpar-dpar07)
      ENDIF
      ENDIF
c      print*,'in envol c1,c2,pa1,ai1,ai2'
c      print*,c1,c2,pa1,ai1,ai2
c      print*,'pa2,pz2',pa2,pz2
      EW=ew0+ea0
      print*,'AI1+AI2-AI0,ea0',ai1,ai2,ai0,ea0
      return
      end

      SUBROUTINE ENYPE(A0,z0,A2E,Z2E,EY0,EN)
C ACEASTA SUBRUTINA NE DA ENERGIA NUCLEARA FOLOSIND UN POTENTIAL
C YUKAWA PLUS EXPONENTIALA. A ESTE NUMARUL DE MASA AL PARINTELUI,
C A2 ESTE NUMARUL DE MASA AL FRAGMENTULUI EMIS, Z ESTE NUMARUL
C ATOMIC AL PARINTELUI, Z2  AL FRAGMENTULUI EMIS
C EN ( OUTPUT ) ESTE ENERGIA NUCLEARA
c aici o modificare, folosesc diferite constante de
c difuzivitate pe suprafata ccare se transmit prin
c common/cdfn/adn0,adn1,adn2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/volumel/v00,v11,v22
      COMMON/NGAUSS/NY
      COMMON/ALUNGI/ALUNG,RR,ADN
      common/r116/R0
      common/lasciziune/delatsciz
      EXTERNAL FEYPE
      PI=3.141592645d0
      RR=R0*A0**(1.D0/3.D0) ! RAZA INITIALA
      DELTA=DELT
      print*,'a1,delta,a2',a1,delta,a2
      ALUNG=(A1+DELTA+A2)/2.D0 ! SEMILUNGIMEA NUCLEU
      DIFER=C1-A1+ALUNG ! POZITIA MIJLOCULUI
      print*,'alung,difer',alung,difer
      CALL RENORM(ALUNG,DIFER)
C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2 (LA FEL CA PENTRU
C ENERGIA COULOMB) PENTRU A AFLA IZOSPINII AI (SAU denumite si excesele de neutroni)
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz!-2 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
c     DPAR07=delatsciz-4
      RO0=Z0/A0
      RO1F=(Z0-Z2E)/(A0-A2E)
      RO2F=Z2E/A2E
      AI0=(A0-2*Z0)/A0
      AI2F=(A2E-2*Z2E)/A2E
      AI1F=((A0-A2E)-2*(Z0-Z2E))/(A0-A2E)
c      print*,'delta,dpar,dpar07'
c      print*,delta,dpar,dpar07
c      print*,'ro0.ro1f,ro2f',ro0,ro1f,ro2f
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      AI2=AI2F
      AI1=AI1F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      AI2=AI0
      AI1=AI0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      A2M=V22/V00*A0
      Z2M=RO2*A2M
      Z1M=Z0-RO2*V22/V00*A0
      A1M=A0-V22/V00*A0
      print*,'z1m,a1m',z1m,a1m
      RO1=Z1M/A1M
      AI2=(A2M-2*Z2M)/A2M
      AI1=(A1M-2*Z1M)/A1M
      ENDIF
      ENDIF
C ADN ESTE CONSTANTA DOMENIULUI FINIT AL FORTELOR NUCLEARE
c     R0=1.16d0
      ADN=.68D0
      AS=21.13D0
      HI=2.3D0
      A21=AS*(1.D0-HI*AI1*AI1)
c     print*,'a21',a21
      A22=AS*(1.D0-HI*AI2*AI2)
      A20=AS*(1.D0-HI*AI0*AI0)
      ADNAD0=ADN*ADN


C CALCULAM ENERGIA NUCLEARA YUKAWA PLUS EXPONENTIALA PENTRU PARINTE SFERIC
c     RR=1.16d0*A**(1.d0/3.d0)
      EY0=A20*A0**.6666666666666666D0*(1.D0-3.D0*ADNAD0/RR/RR+
     +   (1.D0+RR/ADN)*(2.D0+3.D0*ADN/RR*(1.D0+ADN/RR))*
     *   DEXP(-2.D0*RR/ADN))

C CALCULAM O CONSTANTA CU CARE SE INMULTESC INTEGRALELE : CONST=
C =-d**4*(ro/2/a**2)*A20*RR*A0/EY0
      CONST=-(ALUNG/RR)**4*R0/2.D0/ADNAD0*A20*RR*A0/EY0
C SE DEOSEBESC MAI MULTE CAZURI

      A1PA2PR3=A1+A2+2.D0*R3-1.D-6
      IF(A1PA2PR3.LT.DELTA)THEN ! ELIPSE COMPLET SEPARATE
      CALL WGAUSS(FEYPE,-1.D0,U1N,R1)
      CALL SGAUSS(FEYPE,-1.D0,U1N,U2N,1.D0,RR3)
      CALL WGAUSS(FEYPE,U2N,1.D0,R2)
      R1=CONST*R1*A21/A20
      R2=CONST*R2*A22/A20
      RR3=2.D0*CONST*RR3*DSQRT(A21*A22)/A20
      EN=EY0*(R1+R2+RR3)
      print*,'caz ELIPSE COMPLET SEPARATE'
      RETURN
      ENDIF
      IF(S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3.AND.R3.LT.500.d0)THEN !C1<C3<C2
      IF(RO3.LT.R3)THEN ! U1 NU E EGAL CU U2
      CALL INTEGQ(FEYPE,-1.D0,X1N,U1N,U2N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FEYPE,-1.D0,X1N,U1N,NY,R1)
      CALL SINTQ(FEYPE,U2N,X2N,1.D0,NY,R2)
      R1=CONST*R1*A21/A20
      R2=CONST*R2*A22/A20
      RR3=2.D0*CONST*RR3*DSQRT(A21*A22)/A20
      EN=EY0*(R1+R2+RR3)
      print*,'r1,r2,rr3',r1,r2,rr3
      print*,'x1n,u1n,c3n,u2n,x2n'
      print*,x1n,u1n,c3n,u2n,x2n
      print*,'caz U1 NU E EGAL CU U2'
      RETURN
      ELSE ! C1<C3<C2 DAR CELE DOUA FRAGMENTE COMUNICA
      CALL INTEGQ(FEYPE,-1.D0,X1N,C3N,C3N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FEYPE,-1.D0,X1N,C3N,NY,R1)
      CALL SINTQ(FEYPE,C3N,X2N,1.D0,NY,R2)
      R1=CONST*R1*A21/A20
      R2=CONST*R2*A22/A20
      RR3=2.D0*CONST*RR3*DSQRT(A21*A22)/A20
      EN=EY0*(R1+R2+RR3)
      print*,'caz S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3'
      RETURN
      ENDIF
      ELSE
      ! VOLUMUL DIN MIJLOC APARTINE PARTII STANGI
      APINT=(1.d0+X2N)/2.D0
      CALL SINTQ(FEYPE,-1.D0,X1N,X2N,NY,R1)
      CALL INTEGQ(FEYPE,-1.D0,X1N,X2N,X2N,APINT,1.D0,NY,RR3)
      CALL WGAUSS(FEYPE,X2N,1.D0,R2)
      print*,'acon,r1,r2,r3',acon,r1,r2,r3
      R1=CONST*R1*A21/A20
      R2=CONST*R2*A22/A20
      RR3=2.D0*CONST*RR3*DSQRT(A21*A22)/A20
      EN=EY0*(R1+R2+RR3)
      print*,'x1n,x2n',x1n,x2n
      print*,'r1,r2,rr3',r1,r2,rr3
      print*,'caz VOLUMUL DIN MIJLOC APARTINE PARTII STANGI'
      RETURN
      ENDIF
      END



      DOUBLE PRECISION FUNCTION FEYPE(Z1,Z2)
C AC. F. ESTE INT. EN. YUK. PLUS EXP. DUPA CE S-A INTEGRAT DUPA W
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/Z1Z2/X1,X2
      EXTERNAL FEY
      X1=Z1
      X2=Z2
      CALL QGAUSS(FEY,0.D0,1.D0,NN,R)
      FEYPE=R
      RETURN
      END

      DOUBLE PRECISION FUNCTION FEY(W)
C ACEASTA FUNCTIE REPREZINTA INTEGRANTUL ENERGIEI NUCLEARE
C DE DEFORMARE IN APROX YUKAWA PLUS EXPONENTIALA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ALUNGI/AL,R0,ADN
      COMMON/Z1Z2/Z1,Z2
      PI=3.141592645d0
      ALL=AL/R0
C ALL=SEMILUNGIME IMPARTITA LA RAZA ZER0, R0=RAZA ZERO
      PI2=PI+PI
      Z1Z=Z1
      Z2Z=Z2
      CALL VARPAR(Z1Z,RDZ1,RO1)
      CALL VARPAR(Z2Z,RDZ2,RO2)
C CALCULAM EXPRESII INTERMEDIARE CARE SE REPETA PENTRU A REDUCE
C TIMPUL DE CALCUL
      RO1RO1=RO1*RO1
      RO1RO2=RO1*RO2
      RO2RO2=RO2*RO2
      Z1Z2=Z1-Z2
      FI=PI2*W
      COSFI=DCOS(FI)
      RR12FI=RO1RO2*COSFI
      RO1FI=RO1RO1-RR12FI
      RO2FI=RO2RO2-RR12FI
      F1=RO1FI-.5D0*Z1Z2*RDZ1
      F2=RO2FI+.5D0*Z1Z2*RDZ2
      P=RO1FI+RO2FI+Z1Z2*Z1Z2
      RP=DSQRT(P)
      aaa=ADN
c      S=A/AL
      s=aaa/al
      RPS=RP/S
      SP=S/P
      Q=SP*SP*((RPS*(RPS+2.D0)+2.D0)*DEXP(-RPS)-2.D0)
      FEY=F1*F2*Q
      RETURN
      END


	SUBROUTINE QGAUSS(FCT,XL,XU,NN,Y)
C SUBRUTINA DE INTEGRARE PRIN GAUSS IN 32 PUNCTE
	DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
        EXTERNAL FCT
        nn=32
	A=.5D0*(XU+XL)
	B=XU-XL
	C=.49863193092474078D0*B
	Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
	C=.49280575577263417D0*B
	Y=.8137197365452835D-2*(FCT(A+C)+FCT(A-C))+Y
	C=.48238112779375322D0*B
	Y=.12696032654631030D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.46745303796886984D0*B
	Y=.17136931456510717D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.44816057788302606D0*B
	Y=.21417949011113340D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.42468380686628499D0*B
	Y=.25499029631188088D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.39724189798397120D0*B
	Y=.29342046739267774D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.36609105937014484D0*B
	Y=.32911111388180923D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.33152213346510760D0*B
	Y=.36172897054424253D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.29385787862038116D0*B
	Y=.39096947893535153D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.25344995446611470D0*B
	Y=.41655962113473378D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.21067563806531767D0*B
	Y=.43826046502201906D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.16593430114106382D0*B
	Y=.45586939347881942D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.11964368112606854D0*B
	Y=.46922199540402283D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.7223598079139825D-1*B
	Y=.47819360039637430D-1*(FCT(A+C)+FCT(A-C))+Y
	C=.24153832843869158D-1*B
	Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
	RETURN
	END



      SUBROUTINE DIFUZC(A0,Z0,A2E,Z2E,EDIF0,EDIF)
C ACEASTA SUBRUTINA NE DA ENERGIA DIFUZIVITATE COULOMBIANA
C A ESTE NUMARUL DE MASA AL PARINTELUI,
C A2 ESTE NUMARUL DE MASA AL FRAGMENTULUI EMIS, Z ESTE NUMARUL
C ATOMIC AL PARINTELUI, Z2  AL FRAGMENTULUI EMIS
C EN ( OUTPUT ) ESTE ENERGIA NUCLEARA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      common/parameln/a1n,b1n,a2n,b2n,c1n,c2n,c3n,x1n,x2n,r3n,u1n,u2n,
     c               rho3n,sn,x101n,x201n,y101n,y201n
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1,x2,r3,u1,u2,ro3,s,delt
      common/volumel/v00,v11,v22
      COMMON/NGAUSS/NY
      COMMON/ALUNGI/ALUNG,RR,ADN
      common/r116/R0
      common/lasciziune/delatsciz
      

      EXTERNAL FDIFUZ
      PI=3.141592645d0
      RR=R0*A0**(1.D0/3.D0) ! RAZA INITIALA
      DELTA=DELT
      print*,'a1,delta,a2',a1,delta,a2
      ALUNG=(A1+DELTA+A2)/2.D0 ! SEMILUNGIMEA NUCLEU
      DIFER=C1-A1+ALUNG ! POZITIA MIJLOCULUI
      print*,'alung,difer',alung,difer
      CALL RENORM(ALUNG,DIFER)
C ADN ESTE CONSTANTA DOMENIULUI FINIT AL FORTELOR NUCLEARE
      ADN=.99D0/DSQRT(2.D0)
      SCON=1.4399764D0
c     R0=1.16d0

      AL2=ADN*ADN
      ADNADN=AL2
      PIPI=PI*PI
      NG=NY
      RAZA0=RR
      R0R0=R0*R0
      AR0=ADN/RAZA0
      AR02=AR0*AR0
      AR03=AR02*AR0
C CALCULAM EN. DIF. COUL. PT. PARINTE SFERIC
c     RAZA0=1.16d0*A**(1.d0/3.d0)
      RAZA0=RR
      EDIF0=-3.D0*Z0*Z0*SCON*AR02/RAZA0*
     &      (1.D0-15.D0/8.D0*AR0+21.D0/8.D0*AR03-.75D0*DEXP(-2.D0/AR0)*
     &      (1.D0+4.5D0*AR0+7.D0*AR02+3.5D0*AR03))


C PARAMETRIZEZ DENSITATILE RO0, RO1 SI RO2 (LA FEL CA PENTRU
C ENERGIA COULOMB) PENTRU A AFLA IZOSPINII AI (SAU denumite si excesele de neutroni)
c      DPAR=A1+A2+2.D0*R3 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      DPAR=delatsciz-2 ! DUPA ACEASTA VALOARE AM DENSITATILE FINALE
      if(s.lt.0.d0)dpar=1000.
      DPAR07=0.7D0*DPAR
      DPAR07=delatsciz-4
c      RO0=Z0/A0
c      RO1F=(Z0-Z2E)/(A0-A2E)
c      RO2F=Z2E/A2E
c spre deosebire de celelalte subrutine, aici densitatile de sarcina
c se calculeaza ca raport Z/volum ( nu raport Z/A)
      ro0=z0/v00
      ro2f=z2e/(4.d0/3.d0*pi*(R0*a2e**.333333d0)**3)
      ro1f=(z0-z2e)/(4.d0/3.d0*pi*(R0*(a0-a2e)**.333333d0)**3)
c      print*,'delta,dpar,dpar07'
c      print*,delta,dpar,dpar07
c      print*,'ro0.ro1f,ro2f',ro0,ro1f,ro2f
      IF(DELTA.GT.DPAR)THEN
      RO1=RO1F
      RO2=RO2F
      ELSE
      IF(DELTA.LT.DPAR07)THEN
      RO1=RO0
      RO2=RO0
      ELSE
      RO2=RO0+(RO2F-RO0)*(DELTA-DPAR07)/(DPAR-DPAR07)
      Z2M=RO2*V22
      Z1M=Z0-Z2m
      RO1=Z1M/V11
      ENDIF
      ENDIF
C CALCULAM O CONSTANTA CU CARE SE INMULTESC INTEGRALELE 
      ACON=-4.D0*PIPI*AL2*R0R0*SCON*R0
      CONST=(ALUNG/RR)**3*A0*ACON
C SE DEOSEBESC MAI MULTE CAZURI

      A1PA2PR3=A1+A2+2.D0*R3-1.D-6
      IF(A1PA2PR3.LT.DELTA)THEN ! ELIPSE COMPLET SEPARATE
      CALL WGAUSS(FDIFUZ,-1.D0,U1N,R1)
      CALL SGAUSS(FDIFUZ,-1.D0,U1N,U2N,1.D0,RR3)
      CALL WGAUSS(FDIFUZ,U2N,1.D0,R2)
      R1=CONST*R1*RO1*RO1
      R2=CONST*R2*RO2*RO2
      RR3=2.D0*CONST*RR3*RO1*RO2
      EDIF=R1+R2+RR3
      print*,'caz ELIPSE COMPLET SEPARATE'
      RETURN
      ENDIF
      IF(S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3.AND.R3.LT.500.d0)THEN !C1<C3<C2
      IF(RO3.LT.R3)THEN ! U1 NU E EGAL CU U2
      CALL INTEGQ(FDIFUZ,-1.D0,X1N,U1N,U2N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FDIFUZ,-1.D0,X1N,U1N,NY,R1)
      CALL SINTQ(FDIFUZ,U2N,X2N,1.D0,NY,R2)
      R1=CONST*R1*RO1*RO1
      R2=CONST*R2*RO2*RO2
      RR3=2.D0*CONST*RR3*RO1*RO2
      EDIF=R1+R2+RR3
      print*,'Difuziv coul, r1,r2,rr3',r1,r2,rr3
      print*,'difuz coul-x1n,u1n,c3n,u2n,x2n'
      print*,x1n,u1n,c3n,u2n,x2n
      print*,'difuz coul -caz U1 NU E EGAL CU U2'
      RETURN
      ELSE ! C1<C3<C2 DAR CELE DOUA FRAGMENTE COMUNICA
      CALL INTEGQ(FDIFUZ,-1.D0,X1N,C3N,C3N,X2N,1.D0,NY,RR3)
      CALL SINTQ(FDIFUZ,-1.D0,X1N,C3N,NY,R1)
      CALL SINTQ(FDIFUZ,C3N,X2N,1.D0,NY,R2)
      R1=CONST*R1*RO1*RO1
      R2=CONST*R2*RO2*RO2
      RR3=2.D0*CONST*RR3*RO1*RO2
      EDIF=R1+R2+RR3
      print*,'dif coul-caz S.GT.0.D0.AND.C1.LT.C3.AND.C2.GT.C3'
      RETURN
      ENDIF
      ELSE
      ! VOLUMUL DIN MIJLOC APARTINE PARTII STANGI
      APINT=(1.d0+X2N)/2.D0
      CALL SINTQ(FDIFUZ,-1.D0,X1N,X2N,NY,R1)
      print*,'x1n,x2n,x101n,x201n,y101n,y201n'
      print*,x1n,x2n,x101n,x201n,y101n,y201n
      CALL INTEGQ(FDIFUZ,-1.D0,X1N,X2N,X2N,APINT,1.D0,NY,RR3)
      CALL WGAUSS(FDIFUZ,X2N,1.D0,R2)
      print*,'dif coul-acon,r1,r2,r3',acon,r1,r2,rr3
      R1=CONST*R1*RO1*RO1
      R2=CONST*R2*RO2*RO2
      RR3=2.D0*CONST*RR3*RO1*RO2
      EDIF=R1+R2+RR3
      print*,'dif coul=x1n,x2n',x1n,x2n
      print*,'dif coul-r1,r2,rr3',r1,r2,rr3
      print*,'dif coul-caz VOLUMUL DIN MIJLOC APARTINE PARTII STANGI'
      RETURN
      ENDIF
      END





      DOUBLE PRECISION FUNCTION FDIFUZ(Z1,Z2)
C AC. F. ESTE INT. EN. DIF. COUL. DUPA CE S-A INTEGRAT DUPA W
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/Z1Z2/X1,X2
      EXTERNAL FDIFZ
      X1=Z1
      X2=Z2
      CALL QGAUSS(FDIFZ,0.D0,1.D0,NN,R)
      FDIFUZ=R
      RETURN
      END

      DOUBLE PRECISION FUNCTION FDIFZ(W)
C ACEASTA FUNCTIE REPREZINTA INTEGRANTUL ENERGIEI 
C DE DIFUZIVITATE COULOMBIANA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/Z1Z2/Z1,Z2
      COMMON/ALUNGI/AL,R0,A
      ALL=AL/R0
C ALL=SEMILUNGIME IMPARTITA LA RAZA ZER0, R0=RAZA ZERO
      PI=3.14159265358979324D0
      Z1Z=Z1
      Z2Z=Z2
      CALL VARPAR(Z1Z,RDZ1,RO1)
      CALL VARPAR(Z2Z,RDZ2,RO2)
C CALCULAM EXPRESII INTERMEDIARE CARE SE REPETA PENTRU A REDUCE
C TIMPUL DE CALCUL
      RO1RO1=RO1*RO1
      RO1RO2=RO1*RO2
      RO2RO2=RO2*RO2
      Z1Z2=Z1-Z2
      FI=2.D0*PI*W
      COSFI=DCOS(FI)
      F1=RO1RO1-RO1RO2*COSFI-.5D0*Z1Z2*RDZ1
      F2=RO2RO2-RO1RO2*COSFI+.5D0*Z1Z2*RDZ2
      P=RO1RO1+RO2RO2-2.D0*RO1RO2*COSFI+Z1Z2*Z1Z2
      RP=DSQRT(P)
      S=A/R0/ALL
      Q=((.25D0/S*P+1.5D0*RP+2.5*S)*DEXP(-RP/S)-2.5*S+RP)/P/P         
      FDIFZ=F1*F2*Q
      RETURN
      END






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PARTEA DE PARAMETRIZARE CU ELIPSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine zzz
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,delt
      common/r116/R0
      dimension ccc(7)
      R0=1.16
      ccc(1)=-0.07
      ccc(2)=-0.05
      ccc(3)=-0.02
      ccc(4)=-0.01
      ccc(5)=-0.007
      ccc(6)=-0.004
      ccc(7)=-0.002
      mmmm=0
      do i=1,30
      do j=01,50
      mmmm=mmmm+1
c     if(mmmm.eq.  7)then
      delta=4.3d0 +(j-1)*0.3
      cr=-1.4+(i-1)*0.2
      if(i.le.7)then
      cr3=ccc(i)
      else
      cr3=(i-8)*1.d-1
      endif
      dcr3=dabs(cr3)
      a0=235
      z0=92
      a2e=100
      z2e=40
      eps0=0.7
      eps1f=0.6
      eps2f=0.5
c      eps1=0.5
c      eps2=0.5
c      a1pa2=1.5
c excentricitate zero=sfera
c parametrizez eps1,eps2,a1pa2 (excentricitati si raport intre semi-axele a1 si a2)
      call vefinel(a0,z0,a2e,z2e,eps0,eps1f,eps2f,cr3,delta,
     c eps1,eps2,a1pa2,a2m,z2m)

      delta0=delta
       call parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
c     call paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
      delt=delta0
      CALL MPE(A0,Z0,A2E,Z2E,EC0,EC,EY0,EN,EV,ED0,ED,ew,V)
      cc1=0
      cc2=0
      write(11,*)c1,cc1,a1,b1,ier11,ier22
      write(22,*)c2,cc2,a2,b2
      write(33,*)c3,ro3,r3
      x10=x101
      x20=x201
      y1=y101
      y2=y201
      write(55,*)v,ier11,ier22
      write(77,*)x10,x20,y1,y2
      write(88,*)cr3,ier11,ier22
      write(99,*)cr3,c1,x1,x10,c2,x2,x20
      write(66,*)delta0,ec0,ec,ey0,en,ev,ed0,ed,v
c     if(mmmm.eq.2)stop
c     endif
      enddo
      enddo
      end


       subroutine testparamael
c fac un test in care rezolv inainte de toate cazul r3 infinit
c ca sa gasesc limita intervalului in care se poate gasi 
c x1 si x2 pentru orice r3


      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/paramel/a1p,b1p,a2p,b2p,c1,c2,c3,x1,x2,r3,u1,u2,
     c               ro3,sp,delt
      common/r116/R0
      dimension ccc(7)
      R0=1.16
      ccc(1)=-0.07
      ccc(2)=-0.05
      ccc(3)=-0.02
      ccc(4)=-0.01
      ccc(5)=-0.007
      ccc(6)=-0.004
      ccc(7)=-0.002
      do j=1 ,200
      do i=1,100
      do k=1,10
      a1=7.
      b1=4.5
      b1=4.1
      a2=R0*235.d0**(1.d0/3.d0)/1.5d0
      b2=2.5
      b2=b1-(k-1)*0.1 
      delta=(a1-a2)+1.d-3 +(j-1)*0.1
      delta=13.3d0 +(j-1)*0.1
      delta=4.+(j-1)*0.1
      if(i.le.7)then
      cr3=ccc(i)
      else
      cr3=2.d-3+(i-8)*9.d-3
      endif
c      cr3=1+0.1*i
      dcr3=dabs(cr3)
c      if(dcr3.lt.2.d-3)goto2
c      call paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
c      call volels(a1,b1,a2,b2,delta,cr3,v)
      a0=235
       call parelp(a0,b1,a2,b2,cr3,delta,a1)
c     call paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
      cc1=0
      cc2=0
      write(11,*)c1,cc1,a1,b1,ier11,ier22
      write(22,*)c2,cc2,a2,b2
      write(33,*)c3,ro3,r3
      x10=x101
      x20=x201
      y1=y101
      y2=y201
      write(55,*)v,ier11,ier22
      write(77,*)x10,x20,y1,y2
      write(88,*)cr3,ier11,ier22
      write(99,*)cr3,c1,x1,x10,c2,x2,x20
c2     continue
      enddo
      enddo
      enddo
      end


      subroutine parelp2(a0,eps1,eps2,a1pa2,cr3,delta,a1,b1,a2,b2)
c rezolva ecuatia pentru volumul elipsei cand se dau la intrare
c a0 numar de masa parinte
c eps1-eccentricity (excentricitatea) elipsei 1 eps=radical(1-b^2/a^2)
c             excentricitatea este zero pentru cerc
c eps2-excentricitatea elipsei 2
c a1pa2 raport intre a1 si a2
c cr3 curbura gatuirii
c delta distanta dintre centre
c se mai da prin common R_0
c output parametrii elipsei a1,b1,a2,b2
      implicit double precision (a-h,o-z)
      common/a0elp/a00,b11,a22,b22,delt,cr33
      common/ier3parelp/ier3
      common/epsurisir/ep1,ep2,ra12
      common/ierelips/ier11,ier22
      common/limitesol/x101,x201,y101,y201
      common/paramel/a1p,b1p,a2p,b2p,c1p,c2p,c3p,x1p,x2p,r3p,u1p,u2p,
     c               ro3p,sp,deltap
      common/volumel/v00,v11,v22
      common/r116/R_0
      external velidf
      a00=a0
      ep1=eps1
      ep2=eps2
      ra12=a1pa2
      delt=delta
      cr33=cr3
      r0=R_0*a0**(1.d0/3.d0)
      asup=r0/(1.-eps1**2)**(1.d0/3.d0)+0.1
      ainf=0.2*asup
      ddd=(asup-ainf)/100.d0
c daca incep cu limita inferioara pentru gasirea solutiei
c      xx1=ainf
c      vv1=velidf(xx1)
c daca incep cu limita superioara pentru gasirea solutiei
      xx1=asup
      vv1=velidf(xx1)
      do i=1,100
c inferior
c      xx2=xx1+ddd
c      vv2=velidf(xx2)
c superior
      xx2=xx1-ddd
      delt=delta
      vv2=velidf(xx2)
      pr=vv2*vv1
c     print*,'r0,i,xx1,xx2,vv1,vv2 ',r0,i,xx1,xx2,vv1,vv2
c     print*,'ier11,ier22',ier11,ier22
      if(pr.le.0.d0.and.ier11.eq.0.and.ier22.eq.0)then
c inferior
c      ainf=xx1
c      asup=xx2
c superior
      ainf=xx2
      asup=xx1
      goto 1
      endif
      xx1=xx2
      vv1=vv2
      enddo
1     continue
      delt=delta
      call drtmi3(x,f,velidf,ainf,asup,1.e-5,100,ier)
      a1=x
      b1=b11
      a2=a22
      b2=b22
      ier3=ier
      print*,'iesit din drtmi3 a1,b1,a2,b2',a1,b1,a2,b2
      dltp=delta
      call paramelips(a1,b1,a2,b2,dltp,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
      a1p=a1
      b1p=b1
      a2p=a2
      b2p=b2
      c1p=c1
      c2p=c2
      c3p=c3
      x1p=x1
      x2p=x2
      r3p=r3
      u1p=u1
      u2p=u2
      ro3p=ro3
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
      sp=s
      deltap=delta
      return
      end

      double precision function velidf(a1)
c diferenta intre volumul v0 dat de nr de masa a0 si
c voulumul elipselor pentru input excentricitati
      implicit double precision (a-h,o-z)
      common/a0elp/a0,b11,a22,b22,delt,cr33
      common/epsurisir/eps1,eps2,ra1a2
      common/r116/R_0
      r0=R_0*a0**(1.d0/3.d0)
      pi=3.141592645d0
      v0=4.d0/3.d0*r0**3*pi
        if(eps1.gt.0.d0)then
      b11=a1*dsqrt(1.d0-eps1**2)
        else
      b11=a1/dsqrt(1.d0-eps1**2)
        endif
      a22=a1/ra1a2
        if(eps2.gt.0.d0)then
      b22=a22*dsqrt(1.d0-eps2**2)
        else
      b22=a22/dsqrt(1.d0-eps2**2)
        endif
      delta=delt
      cr3=cr33
      b1=b11
      a2=a22
      b2=b22
c     print*,'intru in volels cu a1,b1,a2,b2,delta,cr3(,v)'
c     print*,a1,b1,a2,b2,delta,cr3
      call volels(a1,b1,a2,b2,delta,cr3,v)
c     print*,'IES DIN VOLELS',A1,B1,A2,B2,DELTA,CR3,V
      velidf=v-v0
c     print*,'velidf,v,v0',velidf,v,v0
      b11=b1
      b22=b2
      a22=a2
      return
      end



      subroutine vefinelx(a0,z0,a1pa2i,deltai,eps0,eps1f,eps2f,
     c a2m,z2m)
c calculez eps1 si eps2 si a1pa2 (o parametrizare pe care
c o aleg--- liniara) 
      implicit double precision (a-h,o-z)
      common/r116/R_0




      r0=R_0*a0**(1.d0/3.d0)
      a00=r0/(1-eps1f**2)**(1.d0/3.d0)
      a2=a00/a1pa2i
      deltai=0.05+(a00-a2)
      a2f=r0/(a1pa2i**3*(1.-eps1f**2)+(1.-eps2f**2))**(1./3.)
      a1f=a1pa2i*a2f

      b1f=a1f*dsqrt(1.d0-eps1f**2)
      
      b2f=a2f*dsqrt(1.d0-eps2f**2)

      a2m=a2f*b2f**2/r_0**3
      z2m=z0/a0*a2m

c      z2m=2.+(z2m-2)*(a2m-2)/(a0-2)
       return
       end



      subroutine vefinel(a0,z0,a2,z2,eps0,eps1f,eps2f,cr3,delta,
     c eps1,eps2,a1pa2,a2m,z2m)
c calculez eps1 si eps2 si a1pa2 (o parametrizare pe care
c o aleg--- liniara) 
      implicit double precision (a-h,o-z)
      common/r116/R_0
      r0=R_0*a0**(1.d0/3.d0)
      r2f=R_0*a2**(1.d0/3.d0)
      a1=(a0-a2)
      r1f=R_0*a1**(1.d0/3.d0)
c valorile semi-axelor a finale pentru excentricitatile finale
      a2f=r2f/(1-eps2f**2)**(1.d0/3.d0)
      a1f=r1f/(1-eps1f**2)**(1.d0/3.d0)
c raportul final intre semiaxe
      a1pa2f=a1f/a2f
c facem o interpolare liniara
      deltaim=6
      if(delta.le.deltaim-1)then
       eps1=eps0
       eps2=eps0
       a1pa2=1.d0
       a2m=a0/2.
       z2m=z0/2.
      return
      endif
       if(delta.ge.deltaim)then
       eps1=eps1f
       eps2=eps2f
       a1pa2=a1pa2f
       a2m=a2
       z2m=z2
       return
       endif
       eps1=eps0+(eps1f-eps0)*(delta-deltaim+1)/(1)
       eps2=eps0+(eps2f-eps0)*(delta-deltaim+1)/(1)
       a1pa2=1.d0+(a1pa2f-1.d0)*(delta-deltaim+1)/(1)
       a2m=a0/2.+(a2-a0/2.)*(delta-deltaim+1)/(1)
       z2m=z0/2.+(z2-z0/2.)*(delta-deltaim+1)/(1)
       return
       end


      subroutine parelp(a0,b1,a2,b2,cr3,delta,a1)
c rezolva ecuatia pentru volumul elipsei pentru a0=nr masa parinte
c b1,a2,b2,cr3,delta = parametrii input
c a1 = parametru output
      implicit double precision (a-h,o-z)
      common/a0elp/a00,b11,a22,b22,delt,cr33
      common/ier3parelp/ier3
      common/limitesol/x101,x201,y101,y201
      common/paramel/a1p,b1p,a2p,b2p,c1p,c2p,c3p,x1p,x2p,r3p,u1p,u2p,
     c               ro3p,sp,deltap
      common/volumel/v00,v11,v22
      common/r116/R_0
      external veldif
      a00=a0
      b11=b1
      a22=a2
      b22=b2
      delt=delta
      cr33=cr3
      r0=R_0*a0**(1.d0/3.d0)
      asup=2*r0
      ainf=0.2*r0
      ddd=(asup-ainf)/20.d0
      
      xx1=ainf
      vv1=veldif(xx1)
      do i=1,20
      xx2=xx1+ddd
      vv2=veldif(xx2)
      pr=vv2*vv1
      print*,'xx1,xx2,vv1,vv2',xx1,xx2,vv1,vv2
      if(pr.le.0.d0)then
      ainf=xx1
      asup=xx2
      goto 1
      endif
      xx1=xx2
      vv1=vv2
      enddo
1     continue
      call drtmi3(x,f,veldif,ainf,asup,1.e-3,100,ier)
      a1=x
      ier3=ier
      call paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
      a1p=a1
      b1p=b1
      a2p=a2
      b2p=b2
      c1p=c1
      c2p=c2
      c3p=c3
      x1p=x1
      x2p=x2
      r3p=r3
      u1p=u1
      u2p=u2
      ro3p=ro3
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
      sp=s
      deltap=delta
      return
      end

      double precision function veldif(a1)
c diferenta intre volumul v0 dat de nr de masa a0 si
c voulumul elipselor
      implicit double precision (a-h,o-z)
      common/a0elp/a0,b11,a22,b22,delt,cr33
      common/r116/R_0
      r0=R_0*a0**(1.d0/3.d0)
      pi=3.141592645d0
      v0=4.d0/3.d0*r0**3*pi
      b1=b11
      a2=a22
      b2=b22
      delta=delt
      cr3=cr33
      call volels(a1,b1,a2,b2,delta,cr3,v)
      veldif=v-v0
      return
      end

      subroutine volels(a1,b1,a2,b2,deltai,cr3,v)
c calculez volumul elipsei
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/volumel/v00,v11,v22
      delta=deltai
      call paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2)
c     print*,'paramelips a1,b1,a2,b2,delta'
c     print*,' cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2'
c     print*,a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,r3,u1,u2
      pi=3.141592645d0
c volumul intre capatul din stanga pana la x1
      v1=pi*b1**2*(x1-c1+a1-(x1-c1)**3/(3.d0*a1**2)-a1/3.d0)
c volumul intre capatul din dreapta pana la x2
      v2=pi*b2**2*(-(x2-c2)+a2+(x2-c2)**3/(3.d0*a2**2)-a2/3.d0)
c     print*,'v1,v2',v1,v2
      dcr3=dabs(cr3)
      if(dcr3.lt.2.d-3)then !cazul r3 infinit
      R3=10000.D0
      v3=pi/3.d0*(x201-x101)*(y101**2+y201**2+y101*y201)
c     print*,'v3 1',v3
      v=v1+v2+v3
      v11=v1+v3
      v22=v2
      v00=v
      return
      endif
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
           a1pa2pr3=2.d0*r3+a1+a2-1.d-6
           if(a1pa2pr3.le.delta)then ! cele doua elipse sunt complet separate
           v1=pi*b1**2*(2.d0*a1-2.d0/3.d0*a1)
           v2=pi*b2**2*(2.d0*a2-2.d0/3.d0*a2)
           v=v1+v2
           v11=v1
           v22=v2
           v00=v
c          print*,'2 v=v1+v2',v
           return
           endif
       if(ro3.lt.r3.and.s.gt.0.d0.and.c3.gt.c1.and.c3.lt.c2)then !u1<>u2<>c3
c      v31=(ro3**2+r3**2)*(u1-x1)-1.d0/3.d0*((u1-c3)**3-(x1-c3)**3)-
c    - ro3*r3**2*(dasin((u1-c3)/r3)-dasin((x1-c3)/r3))+
c    + ro3*r3**2*(u1-c3)/r3*dsqrt(1.d0-((u1-c3)/r3)**2)-
c    - ro3*r3**2*(x1-c3)/r3*dsqrt(1.d0-((x1-c3)/r3)**2)
c      v32=(ro3**2+r3**2)*(x2-u2)-1.d0/3.d0*((x2-c3)**3-(u2-c3)**3)-
c    - ro3*r3**2*(dasin((x2-c3)/r3)-dasin((u2-c3)/r3))-
c    - ro3*r3**2*(x2-c3)/r3*dsqrt(1.d0-((x2-c3)/r3)**2)+
c    + ro3*r3**2*(u2-c3)/r3*dsqrt(1.d0-((u2-c3)/r3)**2)
c       print*,'3 v31 v32 ',v31,v32,' u1,u2,c3',u1,u2,c3

      ARGi=(u1 -c3)/R3
      ARG2=(x1-c3)/R3
      V31=(RO3*RO3+R3*R3)*(U1-X1)
     &     +((X1-C3)**3 )/3.D0
     &     +RO3*((X1-C3)*DSQRT(DABS(R3*R3-(X1-C3)**2))
     &    +R3*R3*DASIN(  ARG2)   )
     &      -((U1-C3)**3  )/3.D0
     &   -RO3*((U1 -C3)*DSQRT(DABS(R3*R3-(U1 -C3)**2))
     &    +R3*R3*DASIN(ARGi))
      ARGU=(X2-C3)/R3
      ARGi=(u2-c3)/R3
      V2S=(RO3*RO3+R3*R3)*(X2-U2)
     &      -((X2-C3)**3  )/3.D0
     &   -RO3*((X2-C3)*DSQRT(DABS(R3*R3-(X2-C3)**2))
     &    +R3*R3*DASIN(ARGU))
     &      +((U2 -C3)**3  )/3.D0
     &   +RO3*((U2 -C3)*DSQRT(DABS(R3*R3-(U2 -C3)**2))
     &    +R3*R3*DASIN(ARGi))
        v32=v2s
c       print*,'3 v31 v32 ',v31,v32,' u1,u2,c3',u1,u2,c3

       v=v1+v2+pi*(v31+v32)
       v11=v1+pi*v31
       v22=v2+pi*v32
       v00=v
       return
       endif
       vvv1=(r3**2-(x2-c3)**2)
       if(vvv1.lt.1.d-50)vvv1=1.d-50
       vvv2=(r3**2-(x1-c3)**2)
       if(vvv2.lt.1.d-50)vvv2=1.d-50
       vvv4=(x2-c3)/r3
       if(vvv4.gt.1.d0)  vvv4=1.d0
       if(vvv4.lt.-1.d0)  vvv4=-1.d0
       vvv5=(x1-c3)/r3 
       if(vvv5.gt.1.d0)  vvv5=1.d0
       if(vvv5.lt.-1.d0)  vvv5=-1.d0
       v3=(ro3**2+r3**2)*(x2-x1)-1.d0/3.d0*((x2-c3)**3-(x1-c3)**3)-
     - s*ro3*((x2-c3)*dsqrt(vvv1)+r3**2*dasin(vvv4)-
     -        (x1-c3)*dsqrt(vvv2)-r3**2*dasin(vvv5))
c    - s*ro3*((x2-c3)*dsqrt(r3**2-(x2-c3)**2)+r3**2*dasin((x2-c3)/r3)-
c    -        (x1-c3)*dsqrt(r3**2-(x1-c3)**2)-r3**2*dasin((x1-c3)/r3))
       v3=pi*v3
c       print*,'4 v3',v3
       v=v1+v2+v3
       v00=v
       if(c3.gt.c1.and.c2.gt.c3)then  ! pentru calculul volumelor celor doua elipse
c deosebim 2 cazuri c1<c3<c2  si c2<c3
       v31=(ro3**2+r3**2)*(c3-x1)-1.d0/3.d0*(-(x1-c3)**3)-
     - s*ro3*(-(x1-c3)*dsqrt(r3**2-(x1-c3)**2)-r3**2*dasin((x1-c3)/r3))
       v31=pi*v31
       v32=(ro3**2+r3**2)*(x2-c3)-1.d0/3.d0*((x2-c3)**3)-
     - s*ro3*((x2-c3)*dsqrt(r3**2-(x2-c3)**2)+r3**2*dasin((x2-c3)/r3))
       v32=pi*v32
       v11=v1+v31
       v22=v2+v32
       else
       v11=v1+v3
       v22=v2
       endif
       return
      end

      subroutine paramelips(a1,b1,a2,b2,delta,cr3,c1,c2,x1,x2,c3,ro3,
     c                     r3,u1,u2)
c SE PARE CA MERGE BINE PANA LA R3=1.d6 cu s pozitiv
c                               R3=1.d3 cu s negativ    
c R3 MAXIM DEDUS DIN TESTE ESTE 500 adica dabs(cr3)> 0.002
c dupa care trebuie facuta o interpolare cu r3 infinit
c SE PARE CA NU MERGE DELOC SUB DELTA (a1-a2)+1.d-3
c merge daca b1 este mai mare cu 1 la suta decat b2 sau b1=b2
      implicit double precision (a-h,o-z)
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
       
      call elipsa(a1,b1,a2,b2,delta,c1,c2)
c as putea sa rezolv ecuatia volumului pentru doua cazuri extreme si sa   
c fitez pentru obtinerea parametrilor
      dcr3=dabs(cr3)
      alicr3=2.d-3 ! limita cr3 care mai functioneaza bine
      if(dcr3.lt.alicr3)goto 2 !cazul r3 infinit
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
      r3=1.d0/cr3
      r3=dabs(r3)
             a1pa2pr3=2.d0*r3+a1+a2-5.d-6
           if(a1pa2pr3.le.delta)then
      ier11=0
      ier22=0
c cele doua elipse sunt complet separate
      u1=c1+a1
      u2=c2-a2
      x1=u1
      x2=u2
      ro3=0.d0
      dr=(delta-a1-a2)/2.d0
      c3=u1+dr
      return
           else
c     print*,'intru in solelipd pentru a gasi limitele in care 
c    c se gaseste solutia'
c     PRINT*,'intru in solelipd',a1,b1,a2,b2,c1,c2
      call solelipd(a1,b1,a2,b2,c1,c2,delta,x10,x20,y1,y2)
      x101=x10
      x201=x20
      y101=y1
      y201=y2
c     print*,'intru in solelips x1x2y1y2',x10,x20,y1,y2
      call solelips(a1,b1,a2,b2,c1,c2,delta,s,r3,x1,x2,c3,ro3,u1,u2)
      return
           endif
2     continue
      if(dcr3.lt.1.d-6)cr3=1.d-6
      s=1.d0
      if(cr3.lt.0.d0)s=-1.d0
      call solelipd(a1,b1,a2,b2,c1,c2,delta,x10,x20,y1,y2)
      x101=x10
      x201=x20
      y101=y1
      y201=y2
      alir3=1.d0/alicr3 
      R3=10000.D0
c     r3=alir3+(1.d6-alir3)*(cr3-s*alicr3)/(s*1.d-6-s*alicr3)
c     call solelips(a1,b1,a2,b2,c1,c2,delta,s,alir3,x1,x2,c3,ro3,u1,u2)
c     x1=x1-(x1-x10)*(cr3-s*alicr3)/(s*1.d-6-s*alicr3)
c     x2=x2-(x2-x20)*(cr3-s*alicr3)/(s*1.d-6-s*alicr3)
      ro3=10000.d0
      x1=x101
      x2=x201
      u1=0.d0
      u2=0.d0
      return
      end

      subroutine elipsa(a1,b1,a2,b2,delta,c1,c2)
c calculez c1 si c2 -centrii elipseo
c cei sase  parametrii independenti sunt:
c semiaxele celor doua elipse, raza gatuirii cu semn si distanta dintre centre
      implicit double precision (a-h,o-z)
c conditia ca delta>a1+a2
      a1pa2=a1+a2
      if(delta.gt.a1pa2)then
      c1=-a1*delta/(a1pa2)
      c2=delta+c1
      else
c conditia ca b1=b2 (cazul simetric)
      if(b1.eq.b2)then
      c1=-a1*delta/(a1+a2)
      c2=delta+c1
      else
c conditia dabs(a1-a2)>delta
      if(dabs(a1-a2).gt.delta)then
      print*,'conditia dabs(a1-a2) nesatisfacuta'
c      if(a1.gt.a2)delta=a1-a2+1.d-10
c      if(a2.gt.a1)delta=a2-a1+1.d-10 
      delta=dabs(a1-a2)+1.d-10
      endif
c calculez c1 si c2 centrele elipselor
      b2a2=b2/a2
      b2a22=b2a2*b2a2
      b1a1=b1/a1
      b1a12=b1a1*b1a1
      del2=delta*delta
      rad=del2*b2a22*b2a22-(b2a22-b1a12)*(b1**2-b2**2+b2a22*del2)
c      print*,'rad',rad
      rad=dsqrt(rad)
c      print*,'rad ',rad,' b2a22 b1a12 ',b2a22,b1a12
      b2a2i=b2a2-1.d-6*b2a2
      b2a2s=b2a2+1.d-6*b2a2
      if(b1a1.ge.b2a2s.or.b1a1.le.b2a2i)then
      c1=(-delta*b2a22+rad)/(b2a22-b1a12)
      else
      c1=0.5d0*((b2**2-b1**2)*a1**2/delta/b1**2-delta)
      endif
      c2=delta+c1
c      print*,'c1,c2,delta',c1,c2,delta
c      test=-delta*b2a22
c      test2=(b2a22-b1a12)
      endif
      endif
      return
      end

c PARTEA R3 INFINIT
      subroutine solelipd(a1,b1,a2,b2,c1,c2,delta,x1,x2,y1,y2)
      implicit double precision (a-h,o-z)
      common/fdex1x2i/ac1,bc1,ac2,bc2,cc1,cc2,sc,rc3,xc1,xc2,cc3,deltc
      external fdex1d
c gaseste punctele x1,x2 pentru r3=infinit
      if(b1.eq.b2)then
c caz in care b1=b2
      x1=-a1*delta/(a1+a2)
      x2=delta+x1
      y1=b1
      y2=b2
      else
      ac1=a1
      bc1=b1
      ac2=a2
      bc2=b2
      cc1=c1
      cc2=c2
      deltc=delta
c se cauta solutia in intervalul
      alimi=c1+1.d-7 
      alims=c1+a1-1.d-16
      if(alims.gt.0.d0)alims=0.d0

      itest=0

      goto 7777
6666  continue

        itest=2

        if(b2.gt.b1)then ! se schimba intervalul)
        alimi=c1-a1+1.d-7
        alims=c1-1.d-7
        endif
7777  continue

      call drtmi(x,f,fdex1d,alimi,alims,1.e-10,100,ier)

        if(ier.ne.0.and.b2.gt.b1.and.itest.eq.0)goto 6666


c     print*,'ier drtmi pentru r3 infinit',ier
      x1=x
      x2=xc2
      y1=b1*dsqrt(1.d0-(x1-c1)**2/a1**2)
      y2=b2*dsqrt(1.d0-(x2-c2)**2/a2**2)
      endif
      return
      end

      double precision function fdex1d(x1)
      implicit double precision (a-h,o-z)
      common/fdex1x2i/a1,b1,a2,b2,c1,c2,sc,rc3,xc1,x2,cc3,delta
      x1mc1=x1-c1
      ra1=1.d0-((x1mc1)/a1)**2
      x1mc12=(x1mc1)**2
      x2=x1mc1*ra1+(b2*a1/(a2*b1))**2*c2*ra1+x1/a1**2*x1mc12
      x2=x2/((b2*a1/(a2*b1))**2*ra1+x1mc12/a1**2)
      x1mc1pa1=((x1-c1)/a1)**2
      a1pa2=a1+a2
      if(delta.ge.a1pa2)then
c     print*,'delta,a1pa2',delta,a1pa2
      if(x1mc1pa1.ge.0.9999999d0)x1mc1pa1=0.9999999d0
      endif
      f=       b1*(x1-c1)/a1**2/dsqrt(1.d0-x1mc1pa1)
      x2mc2pa2=((x2-c2)/a2)**2
c     a1pa2=a1+a2
      if(delta.ge.a1pa2)then
c     print*,'delta,a1pa2',delta,a1pa2
      if(x2mc2pa2.ge.0.9999999d0)x2mc2pa2=0.9999999d0
      endif
      if(x2mc2pa2.ge.0.9999999d0)x2mc2pa2=0.9999999d0
      fdex1d=f-b2*(x2-c2)/a2**2/dsqrt(1.d0-x2mc2pa2)
c     print*,'x2mc2pa2 mai mic ca 1',x2mc2pa2
c      print*,'f',f
c      print*,'x1,x2,fdex1d',x1,x2,fdex1d
      return
      end
      
C PARTEA IN CARE R3 ESTE DIFERIT DE INFINIT
      subroutine solelips(a1,b1,a2,b2,c1,c2,delta,s,r3,x1,x2,c3,ro3,
     c u1,u2)
      implicit double precision (a-h,o-z)
      common/fdex1x2/ac1,bc1,ac2,bc2,cc1,cc2,sc,rc3,xc1,xc2,cc3,dt,ier2
      common/limitesol/x101,x201,y101,y201
      common/ierelips/ier11,ier22
      common/nopindex/noind,noind1,noind2
      external fdex1
c gaseste punctele x1,x2,c3,ro3 
c
c b1>=b2 a1>=a2
      noind2=0
      ac1=a1
      bc1=b1
      ac2=a2
      bc2=b2
      cc1=c1
      cc2=c2
      sc=s
      rc3=r3
      dt=delta 
      Ni=100
c solutia x1 se afla in diferite intervale
      if(s.gt.0.d0)then
c impart intervalul in Ni parti:
!                         S POZITIV 
      alimii=c1+1.d-16
      alimss=c1+a1-1.d-16
      alimii=x101-1.d-10
c     alimii=x101
      if(alimss.gt.0.d0)alimss=0.d0
c pentru solutie intervalul nu poate fi mai mare decat r3
c fata de origine cu exceptia cazului special 
c CAZUL SPECIAL IN CARE CEI DOI SFEROIZI SE SEPARA S=1
c ATUNCI LIMITELE SOLUTIILOR NU POT EXISTA DECAT IN INTERVALUL
c [C1,C1+a1] si [c2,c2-a2]
      a1pa2=a1+a2
      if(delta.lt.a1pa2)then
      alimin=-r3
      else
c cazul in care cercul 3 se gaseste intre cele doua elipse si
c cele doua elipse sunt separate. Scad din 2*r3 diferenta (a1-a2)
c si maresc intervalul in care caut aceasta solutie cu (2*r3-dabs(a1-a2))
c      alimin=-r3-dabs(2.d0*r3-dabs(a1-a2))
      alimin=-r3-dabs(a1-a2)/2.d0
      endif
      if(alimin.gt.alimii)alimii=alimin
      dddd=(alimss-alimii)/Ni
      do ii=1,Ni                 !+3
      alimi=alimii+(ii-1)*dddd
      alims=alimii+ii*dddd
                            !alimi=alimii+(ii-3)*dddd
                            !alims=alimii+(ii-2)*dddd
      if(alims.gt.0.d0)alims=0.d0
      noind=0
c          a1pla2=a1+a2
c          if(delta.ge.a1pla2)then
c          c1pla1=c1+a1-1.d-16
c          if(alims.gt.c1pla1)alims=c1pla1
c          c1pla1m=c1pla1-1.d-1
c          if(alimi.gt.c1pla1m)alimi=c1pla1m
c          endif
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      noind=1
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      if(noind1.eq.0)then
      noind=2
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      noind=3
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
            endif      !de la noind1.eq.0
      enddo

! MAIN INCERC
      alimii2=c1-a1
      dddd=(alimss-alimii2)/Ni
            print*,'alimii2,dddd,alimss',alimii2,dddd,alimss,c1,a1
      do ii=1,Ni                 !+3
      alimi=alimii2+(ii-1)*dddd
      alims=alimii2+ii*dddd
                            !alimi=alimii+(ii-3)*dddd
                            !alims=alimii+(ii-2)*dddd
      if(alims.gt.0.d0)alims=0.d0
      noind=0
c       print*,'alimi,alims',alimi,alims
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'1 x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
c      if(ier.eq.0.and.ier2.eq.0)stop
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      noind=1
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'2 x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      if(noind1.eq.0)then
      noind=2
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'3 x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      noind=3
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
c     print*,'4 x1,x2,fdex1 ',x,x2c,f,' ier,ier2 ',ier,ier2
c     print*,' a fost alimi,alims ',alimi,alims
      if(ier.eq.0.and.ier2.eq.0)goto 1312
            endif      !de la noind1.eq.0
      enddo

      noind2=1
c cazul in care cercul 3 se gaseste intre cele doua elipse si
c cele doua elipse sunt separate. Mai fac o incercare fiindca
c s-a dovedit ca n-a functionat bine
      alimi=c1+a1-2.d-3
      alims=c1+a1-1.d-8
      if(alims.gt.0.d0)then
      alims=0.d0
      alimi=-2.d-3
      endif
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
      if(ier.eq.0.and.ier2.eq.0)goto 1312
      alimi=c1+a1-1.d-2
      alims=c1+a1-1.d-8
      if(alims.gt.0.d0)then
      alims=-2.d-3
      alimi=-1.d-1
      endif
      call drtmi(x,f,fdex1,alimi,alims,1.e-6,100,ier)
1312  continue
      else
c impart intervalul in Ni parti: 
!                        S NEGATIV
! caut solutia intre alimii=c1-a1+1.d-16 si alimss=c1 (primul interval) si
!                    alimii=c1 si alimss=x101+1.d-5  (al doilea interval)
      if(s.lt.0.d0.and.r3.gt.100.d0)then
      Ni=400
      endif
          nprimacaut=0
      alimii=c1-a1-1.d-16
      alimss=c1+a1-1.d-16
      alimss=x101
      if(alimss.gt.0.d0)alimss=0.d0
      alimss=c1+1.d-16
2131  continue
      dddd=(alimss-alimii)/Ni
      do ii=1,Ni+1
      alimi=alimii+(ii-1)*dddd
      alims=alimii+ii*dddd
      noind=0
      if(s.lt.0.d0.and.r3.gt.100.d0)then
      call drtmi(x,f,fdex1,alimi,alims,1.e-2,100,ier)
      else
      call drtmi(x,f,fdex1,alimi,alims,1.e-5,100,ier)
      endif
      if(ier.eq.0.and.ier2.eq.0)goto 1212
c!!!!!!!! ca limita superioara ar trebui sa pun punctul de intersectie al elipselor
      enddo
        if(nprimacaut.eq.0)then
         alimss=c1+a1
         nprimacaut=1
         goto 2131
         endif
c fac urmatorul algoritm, iau pe x101, il impart la zece si caut solutia x1
c in cincizeci de intervale [i*x101/10,(i-1)*x101/10] (i=1,Ni). La fel procedez si pentru
c solutia x2
      dddd=(c1-a1-x101)/(Ni-1) ! ar trebui sa fie negativ
      do ii=1,Ni
      alims=x101+(ii-2)*dddd
      alimi=x101+(ii-1)*dddd
c     print*,'INTRU drtmi cu alimi,alims ',alimi,alims,' x101,x201 ',
c    cx101,x201
      noind=1
      if(s.lt.0.d0.and.r3.gt.100.d0)then
      call drtmi(x,f,fdex1,alimi,alims,1.e-2,100,ier)
      else
      call drtmi(x,f,fdex1,alimi,alims,1.e-5,100,ier)
      endif
c     print*,'IES din drtmi cu '
c     print*,'x1,x2 ',x,xc2,' f,alimi,alims ',f,alimi,alims,' c1,x101 ',
c    c c1,x101,' ier,ier2 ',ier,ier2
      if(ier.eq.0.and.ier2.eq.0)goto 1212
      enddo
1212  continue
      endif
      ier11=ier
      ier22=ier2
c      if(ier.eq.0.and.ier2.eq.0)write(74,*)r3,ier,ier2
      x1=x
      x2=xc2
      c3=x1+s*r3*b1/a1**2*(x1-c1)/dsqrt(1.d0+
     + (b1**2/a1**4-1.d0/a1**2)*(x1-c1)**2)
      ro3=b1*dsqrt(1.d0-(x1-c1)**2/a1**2)+
     + s*dsqrt(r3**2-(x1-c3)**2)
       gh3=1.d0-(x1-c1)**2/a1**2
       gh4=r3**2-(x1-c3)**2
           if(gh3.le.0.d0)gh3=0.d0
           if(gh4.le.0.d0)gh4=0.d0
        ro3=b1*dsqrt(gh3)+s*dsqrt(gh4)


      u1=c3
      u2=c3
      if(ro3.lt.r3.and.s.gt.0.d0.and.c3.gt.c1.and.c3.lt.c2)then !u1<>u2<>c3
c calculez u1 si u2
      dif=dsqrt(r3*r3-ro3*ro3)
      u1=c3-dif
      u2=c3+dif
      endif
      return
      end

      double precision function fdex1(x1)
      implicit double precision (a-h,o-z)
      common/fdex1x2/a1,b1,a2,b2,c1,c2,s,r3,xx1,x2,c3,delta,ier2
      common/limitesol/x101,x201,y101,y201
      common/nopindex/noind,noind1,noind2
      external fdex2
      xx1=x1      
      Ni=100
c solutia x2 se afla in diferite intervale      
      if(s.lt.0.d0)then
!                           S NEGATIV
        if(noind.eq.0)then
c impart intervalul in Ni parti: 
      alimss=c2+a2+1.d-16
      alimii=c2
      alimii=x201
      if(alimii.lt.0.d0)alimii=0.d0
      dddd=(alimss-alimii)/Ni
      do ii=1,Ni+2
      alimi=alimii+(ii-2)*dddd
      alims=alimii+(ii-1)*dddd
c     print*,'   INTRU drtmi2 cu alimi,alims ',alimi,alims
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8,100,ier)
c     print*,' IES drtmi2, x2,fdex2,ier',x,f,ier
      ier2=ier
      if(ier.eq.0)goto 1312
      enddo
         else
      dddd=(a2+c2-x201)/(Ni-2)
      do ii=1,Ni
      alimi=a2+c2-(ii-1)*dddd
      alims=a2+c2-(ii-2)*dddd
c     print*,'   INTRU drtmi2 cu alimi,alims ',alimi,alims
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8,100,ier)
c     print*,' IES drtmi2, x2,fdex2,ier',x,f,ier
      ier2=ier
      if(ier.eq.0)goto 1312
      enddo
        endif
1312  continue
      else
!!!! ca limita inferioara ar trebui punctul intersectie al elipselor
c impart intervalul in Ni parti: 
!                         S POZITIV
                if(noind2.eq.1)then
c cazul in care cercul 3 se gaseste intre cele doua elipse si
c cele doua elipse sunt separate. Mai fac o incercare fiindca
c s-a dovedit ca n-a functionat bine
      alimi=c2-a2+1.d-6
      alims=c2-a2+5.d-2
      c2ma2=c2-a2
      if(c2ma2.lt.0.d0)then
      alimi=1.d-20
      alims=5.d-2
      endif
      ali=alimi
      dddd=(alims-alimi)/20.d0
      do ia=1,100
      alimi=ali+(ia-1)*dddd
      alims=ali+ia*dddd
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
      do ia=1,20
      alimi=ali+20*(ia-1)*dddd
      alims=ali+20*ia*dddd
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
      goto 1212
                else
      alimii=c2-a2+1.d-16
      if(alimii.lt.0.d0)alimii=0.d0
      alimss=c2+a2-1.d-16
      alms=c2+a2
      alimss=x201+1.d-10
c     alimss=x201
      a1pa2=a1+a2
      if(delta.lt.a1+a2)then
      alimsu=2*r3
      if(alimsu.gt.alms)alimsu=alms
      else
c      alimsu=r3+dabs(2.d0*r3-dabs(a1-a2))
      alimsu=2*r3+dabs(a1-a2)/2.d0
      if(alimsu.gt.alms)alimsu=alms
      endif
c            a1pla2=a1+a2
c            c2ma2=c2-a2+1.d-16
c            c2ma2s=c2ma2+1.d-1
c CAZUL SPECIAL IN CARE CEI DOI SFEROIZI SE SEPARA S=1
c ATUNCI LIMITELE SOLUTIILOR NU POT EXISTA DECAT IN INTERVALUL
c [C1,C1+a1] si [c2,c2-a2]
c     a1pa2=a1+a2
c     if(delta.lt.a1pa2)then
      if(alimsu.lt.alimss)alimss=alimsu 
c     endif
c AM PATRU INTERVALE IN CARE LUCREZ X2<C2 X2>C2 X2<C3 X2>C3
c si trebuie cautata solutia in aceste intervale, fiecare
c interval este indexat cu un numar intreg noind
c Deoarece pe c3 nu-l cunosc utilizez doar doua intervale
c [alimii,c2] si [c2,alimss] si ambele intervale vor fi baleiate
c in ordine crescatoare si descrescatoare
c mai fac un indice care este 1 daca c2 nu este in intervalul [alimii,alimss]
      dddd=(alimss-alimii)/Ni

         if(c2.le.alimss.and.c2.ge.alimii)then
         noind1=0
         dddd=(c2-alimii)/Ni
         dddd1=(alimss-c2)/Ni
         else 
         noind1=1
         dddd=(alimss-alimii)/Ni
         endif
                if(noind1.eq.0)then  
                if(noind.eq.0)then
! FAC BUCLA CRESCATOR IN PRIMUL INTERVAL                 
      do ii=1,Ni                  !+2
      alimi=alimii+(ii-1)*dddd               !(ii-2)*dddd
      alims=alimii+ii*dddd                   !(ii-1)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
           if(itest.eq.71.or.itest.eq.72.or.itest.eq.73)then
       print*,'pt fdex2 alimii,ddd,alimi,alims',alimi,alims
       endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
               else
               if(noind.eq.1)then
! FAC BUCLA INVERS IN PRIMUL INTERVAL

      do ii=1,Ni                  !+2
      alims=c2-(ii-1)*dddd                !(ii-2)*dddd
      alimi=c2-ii*dddd                    !(ii-1)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
              else
              if(noind.eq.2)then
! FAC BUCLA CRESCATOR IN AL DOILEA INTERVAL
         dddd=(alimss-c2)/Ni
      do ii=1,Ni
      alimi=c2+(ii-1)*dddd
      alims=c2+(ii)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
               else
! FAC BUCLA DESCRESCATOR IN AL DOILEA INTERVAL
         dddd=(alimss-c2)/Ni
      do ii=1,Ni
      alimi=alimss-(ii)*dddd
      alims=alimss-(ii-1)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
               endif
               endif
               endif
               endif ! de la noind1=0
                           if(noind1.eq.1)then  
                if(noind.eq.0)then
! FAC BUCLA CRESCATOR IN PRIMUL INTERVAL
                 dddd=(alimss-alimii)/Ni                 
      do ii=1,Ni                  !+2
      alimi=alimii+(ii-1)*dddd               !(ii-2)*dddd
      alims=alimii+ii*dddd                   !(ii-1)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
           if(itest.eq.71.or.itest.eq.72.or.itest.eq.73)then
       print*,'pt fdex2 alimii,ddd,alimi,alims',alimi,alims
       endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo
               else
! FAC BUCLA INVERS IN PRIMUL INTERVAL

                 dddd=(alimss-alimii)/Ni
      do ii=1,Ni                  !+2
      alims=alimss-(ii-1)*dddd                !(ii-2)*dddd
      alimi=alimss-ii*dddd                    !(ii-1)*dddd
c            if(delta.gt.a1pla2)then
c            if(alims.lt.c2ma2s)alims=c2ma2s
c            if(alimi.lt.c2ma2)alimi=c2ma2
c            endif
      call drtmi2(x,f,fdex2,alimi,alims,1.e-8 ,100,ier)
c        print*,'drtmi2 alimi,alims ',alimi,alims,' ier2 ',ier
      ier2=ier
      if(ier.eq.0)goto 1212
      enddo  
              endif            
                      endif ! de la noind1=1
      endif
1212  continue
      endif
c     print*,'ier din drtmi2 s=1',ier
      x2=x
      ra1=1.d0+
     + (b1**2/a1**4-1.d0/a1**2)*(x1-c1)**2
      ra2=1.d0+
     + (b2**2/a2**4-1.d0/a2**2)*(x2-c2)**2
      if(ra1.lt.1.d-50)ra1=1.d-50
      if(ra2.lt.1.d-50)ra2=1.d-50
      x1mc3=-s*r3*b1/a1**2*( x1-c1)/dsqrt(ra1)
      x2mc3=-s*r3*b2/a2**2*(x2-c2)/dsqrt(ra2)
c      x1mc3=-s*r3*b1/a1**2*( x1-c1)/dsqrt(1.d0+
c     + (b1**2/a1**4-1.d0/a1**2)*(x1-c1)**2)
c      x2mc3=-s*r3*b2/a2**2*(x2-c2)/dsqrt(1.d0+
c     + (b2**2/a2**4-1.d0/a2**2)*(x2-c2)**2)
             
             r1=1.d0-(x1-c1)**2/a1**2
             r13=(r3**2-x1mc3**2)
             r2=(1.d0-(x2-c2)**2/a2**2)
             r23=(r3**2-x2mc3**2)
       if(r1.lt.1.d-50)r1=1.d-50
       if(r13.lt.1.d-50)r13=1.d-50
       if(r2.lt.1.d-50)r2=1.d-50
       if(r23.lt.1.d-50)r23=1.d-50
      fdex1=b1*dsqrt(r1)+
     + s*dsqrt(r13)-
     -      b2*dsqrt(r2)-
     - s*dsqrt(r23) 
c      fdex1=b1*dsqrt(1.d0-(x1-c1)**2/a1**2)+
c     + s*dsqrt(r3**2-x1mc3**2)-
c     -      b2*dsqrt(1.d0-(x2-c2)**2/a2**2)-
c     - s*dsqrt(r3**2-x2mc3**2)
c     print*,'x1,x2,fdex1,fdex2',x1,x2,fdex1,f,' ier2',ier 
c     print*,'x1mc3,x2mc3',x1mc3,x2mc3
c     print*,'r1,r13,r2,r23',r1,r13,r2,r23
      return
      end

      double precision function fdex2(x2)
      implicit double precision (a-h,o-z)
      common/fdex1x2/a1,b1,a2,b2,c1,c2,s,r3,x1,x22,c3,delta,ier2
      r1=1.d0+(b2**2/a2**4-1.d0/a2**2)*
     * (x2-c2)**2
c     if(r1.lt.1.d-50)r1=1.d-50
      r2=1.d0+(b1**2/a1**4-1.d0/a1**2)*
     * (x1-c1)**2
c     if(r2.lt.1.d-50)r2=1.d-50
c     fdex2=x2-x1-s*r3*(
c    ( -b2/a2**2*(x2-c2)/dsqrt(r1)+
c    + b1/a1**2*(x1-c1)/dsqrt(r2))
       sr1=1.d0
       if(r1.lt.0.d0)sr1=-1.d0
       sr2=1.d0
       if(r2.lt.0.d0)sr2=-1.d0
c       r1=dabs(r1)
c       r2=dabs(r2)
      fdex2=sr1*sr2*dsqrt(r1*r2)*(x2-x1)-s*r3*(
     ( -b2/a2**2*(x2-c2)*dsqrt(r2)*sr2+
     + b1/a1**2*(x1-c1)*dsqrt(r1)*sr1)
c      fdex2=x2-x1-s*r3*(
c     ( -b2/a2**2*(x2-c2)/dsqrt(1.d0+(b2**2/a2**4-1.d0/a2**2)*
c     * (x2-c2)**2)+
c     +  b1/a1**2*(x1-c1)/dsqrt(1.d0+(b1**2/a1**4-1.d0/a1**2)*
c     * (x1-c1)**2))
      return
      end

C
C     ..................................................................
C
C        SUBROUTINE DRTMI
C
C        PURPOSE
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM FCT(X)=0
C           BY MEANS OF MUELLER-S ITERATION METHOD.
C
C        USAGE
C           CALL DRTMI (X,F,FCT,XLI,XRI,EPS,IEND,IER)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.
C
C        DESCRIPTION OF PARAMETERS
C           X      - DOUBLE PRECISION RESULTANT ROOT OF EQUATION
C                    FCT(X)=0.
C           F      - DOUBLE PRECISION RESULTANT FUNCTION VALUE
C                    AT ROOT X.
C           FCT    - NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           XLI    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE
C                    INITIAL LEFT BOUND OF THE ROOT X.
C           XRI    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE
C                    INITIAL RIGHT BOUND OF THE ROOT X.
C           EPS    - SINGLE PRECISION INPUT VALUE WHICH SPECIFIES THE
C                    UPPER BOUND OF THE ERROR OF RESULT X.
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS
C                     IER=0 - NO ERROR,
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS
C                             FOLLOWED BY IEND SUCCESSIVE STEPS OF
C                             BISECTION,
C                     IER=2 - BASIC ASSUMPTION FCT(XLI)*FCT(XRI) LESS
C                             THAN OR EQUAL TO ZERO IS NOT SATISFIED.
C
C        REMARKS
C           THE PROCEDURE ASSUMES THAT FUNCTION VALUES AT INITIAL
C           BOUNDS XLI AND XRI HAVE NOT THE SAME SIGN. IF THIS BASIC
C           ASSUMPTION IS NOT SATISFIED BY INPUT VALUES XLI AND XRI, THE
C           PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           SOLUTION OF EQUATION FCT(X)=0 IS DONE BY MEANS OF MUELLER-S
C           ITERATION METHOD OF SUCCESSIVE BISECTIONS AND INVERSE
C           PARABOLIC INTERPOLATION, WHICH STARTS AT THE INITIAL BOUNDS
C           XLI AND XRI. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF
C           FCT(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP
C           REQUIRES TWO EVALUATIONS OF FCT(X). FOR TEST ON SATISFACTORY
C           ACCURACY SEE FORMULAE (3,4) OF MATHEMATICAL DESCRIPTION.
C           FOR REFERENCE, SEE G. K. KRISTIANSEN, ZERO OF ARBITRARY
C           FUNCTION, BIT, VOL. 3 (1963), PP.205-206.
C
C     ..................................................................
C
      SUBROUTINE DRTMI(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C
C
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
C
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
    1 FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
    2 FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCE FOR FUNCTION VALUES.
    3 I=0
      TOLF=100.*EPS
C
C
C     START ITERATION LOOP
    4 I=I+1
C
C     START BISECTION LOOP
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
    5 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
    6 TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
    7 TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
    8 IF(I-IEND)17,17,9
    9 XR=X
      FR=F
C
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
   10 TOL=TOL*A
   11 IF(DABS(XR-XL)-TOL)12,12,13
   12 IF(DABS(FR-FL)-TOLF)14,14,13
   13 CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.
      IER=1
   14 IF(DABS(FR)-DABS(FL))16,16,15
   15 X=XL
      F=FL
   16 RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
   17 A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP
   18 TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
   19 TOL=TOL*A
   20 IF(DABS(DX)-TOL)21,21,22
   21 IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
   22 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
   23 XR=X
      FR=F
      GO TO 4
   24 XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C     END OF ITERATION LOOP
C
C
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
   25 IER=2
      RETURN
      END
      

      SUBROUTINE DRTMI2(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C
C
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
C
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
    1 FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
    2 FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCE FOR FUNCTION VALUES.
    3 I=0
      TOLF=100.*EPS
C
C
C     START ITERATION LOOP
    4 I=I+1
C
C     START BISECTION LOOP
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
    5 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
    6 TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
    7 TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
    8 IF(I-IEND)17,17,9
    9 XR=X
      FR=F
C
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
   10 TOL=TOL*A
   11 IF(DABS(XR-XL)-TOL)12,12,13
   12 IF(DABS(FR-FL)-TOLF)14,14,13
   13 CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.
      IER=1
   14 IF(DABS(FR)-DABS(FL))16,16,15
   15 X=XL
      F=FL
   16 RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
   17 A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP
   18 TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
   19 TOL=TOL*A
   20 IF(DABS(DX)-TOL)21,21,22
   21 IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
   22 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
   23 XR=X
      FR=F
      GO TO 4
   24 XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C     END OF ITERATION LOOP
C
C
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
   25 IER=2
      RETURN
      END
      

      SUBROUTINE DRTMI3(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C
C
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
C
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
    1 FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
    2 FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCE FOR FUNCTION VALUES.
    3 I=0
      TOLF=100.*EPS
C
C
C     START ITERATION LOOP
    4 I=I+1
C
C     START BISECTION LOOP
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
    5 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
    6 TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
    7 TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
    8 IF(I-IEND)17,17,9
    9 XR=X
      FR=F
C
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
   10 TOL=TOL*A
   11 IF(DABS(XR-XL)-TOL)12,12,13
   12 IF(DABS(FR-FL)-TOLF)14,14,13
   13 CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.
      IER=1
   14 IF(DABS(FR)-DABS(FL))16,16,15
   15 X=XL
      F=FL
   16 RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
   17 A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP
   18 TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
   19 TOL=TOL*A
   20 IF(DABS(DX)-TOL)21,21,22
   21 IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
   22 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
   23 XR=X
      FR=F
      GO TO 4
   24 XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C     END OF ITERATION LOOP
C
C
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
   25 IER=2
      RETURN
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! TERMINAT PARTE PARAMETRIZARE CU ELIPSE !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!EFECTE DE PATURI SI DE IMPERECHERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pairnou(idpsn,a0,z0,ess,
     c   gammas,dp,sollam,soldel)
c calculeaza efectul de pairing
c idpsn=1 pt neutroni si 0 pentru protoni
      implicit double precision (a-h,o-z)
       dimension ess(2925)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/g,nrnivc,iier
       common/deltaaa/del,fermi
      common/enmoment/emom(21,325),eim(2,21)
c     common/rapginte/rapgint 
          common/rapsup/rapsu

       external bcs1nou,partnrnou

          nrnivc=28 !parametru de taiere
      if(idpsn.eq.1)then
      nrocup=((a0-z0)-nrimpare)/2
      else
      nrocup=(z0-nrimpare)/2
      endif
          if(nrnivc.gt.nrocup)nrnivc=nrocup
      
      do i=1,nrnive
      eonr(i)=ess(i)/hw0
      enddo

      gamma=gammas
      nrpart=nrocup
      cdet=12.d0
c      print*,'hw0',hw0
            
      do i=1,nrnive-1
      j=nrnive-i+1 
      jk1=j
      xli=eonr(jk1)
      jk2=j-1
      xri=eonr(jk2)
      gsd=partnrnou(xli)*partnrnou(xri)
      if(gsd.le.0.d0)goto 142
      enddo
142    continue
      call drtmi(xx,ff,partnrnou,xli,xri,1.e-4,80,ier)
      if(ier.ne.0)then
      print*,'par ier nu e zero cand calculam lambda tilda'
      print*,'ier=',ier,'  xx=',xx,' ff=',ff,'  xli=',xli,'  xri=',xri
      endif
       alambd=xx
c    
c    
      gtilam=gtildnou(alambd)
      print*,'gtilam,alambd',gtilam,alambd

          nrnicc=hw0*gtilam/2
         ncnv=nrnivc
          if(ncnv.lt.nrnicc)ncnv=nrnicc
          if(ncnv.gt.nrocup)ncnv=nrocup
         print*,'nrnivc,nrnicc,nrocup',nrnivc,nrnicc,nrocup

      omega2=nrnivc/gtilam

      deltil=cdet/dsqrt(a0)/hw0

      g=1.d0/gtilam/dlog(2.*omega2/deltil)
         print*,'gvechi',g
           g=g*rapsu       
C FOLOSESC GVECHI
        goto 91

      nrniv=nrnive
      do j=1,nrnive!-nrnivc,nrnivc-1
      emom(1,j)=ess(j)
      enddo
      call ginterac(a0,0,nrniv,nrocup,nrnivc,hw0,gval)
      g=gval
         g=g/hw0
      print*,'g din renormalizare',g,'omega2,hw0',omega2,hw0

91      continue
      
                 noption=0 
 6000     continue
      d=.001d0
      dmin=d
      v1=bcs1nou(dmin)
         ii1=iier
1000  continue
      dmax=dmin+.001d0
      v2=bcs1nou(dmax)
         ii2=iier
      prod=v1*v2
c     print*,dmin,dmax,v1,v2,prod
      if(prod.le.0.d0.and.ii1.eq.0.and.ii2.eq.0)goto 2000
      v1=v2
      ii1=ii2
      dmin=dmax
      if(dmin.gt.1.2d0.and.noption.eq.0)goto 2001
      if(dmin.gt.1.2d0.and.noption.eq.2)goto 2000
!      if(dmin.gt.0.8d0.and.noption.eq.0)goto 5000
!      if(dmin.gt.0.8d0.and.noption.eq.1)goto 8000
      goto 1000
2001    continue
!5000  continue
!c     g=1.d0/gtilam/dlog(2.*omega2/(20/dsqrt(a0)/hw0))
!       nrnivc=ncnv
!      noption=1
!      goto 6000
!8000  continue
       g=g*1.05 !!!!!!!!!!! FAC O VALOAREA MAI MICA
       noption=2
       goto 6000



2000  continue

      call drtmi(xx,ff,bcs1nou,dmin,dmax,1.e-4,50,ier)
      
      if(iier.ne.0.or.ier.ne.0) then 
      print*,'nu se rezolva bine bcs2 sau bcs1'
      print*,'ier in rezolvarea bcs1',ier
      print*,'ier in rezolvarea bcs2',iier
c      soldel=0.01
c      fermi=0.5*(eonr(nrnive-nrocup)+eonr(nrnive+1-nrocup))
      xx=0.001
      v2=bcs1nou(xx)
      endif

      soldel=xx
      sollam=fermi
       print*,'soldel sollam',soldel,sollam
        sssss=soldel*hw0
        ddddd=sollam*hw0
        print*,'delta fermi',sssss,ddddd,'   hw0',hw0

c se face sumarea pentru P
      sum1=0.d0


      do j=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+j
     
      ene=eonr(n1)-sollam
      eps=dsqrt(ene*ene+soldel*soldel)
      vk2=0.5*(1.-ene/eps)
      if(n1.ge.nrnive+1-nrocup)then
      v1=2*eonr(n1)
      else
      v1=0
      endif
      sum1=sum1+2*vk2*eonr(n1)-v1
c      eps2=ene*ene+.5d0*soldel*soldel
c      sum1=sum1+dabs(ene)-eps2/eps
      
      enddo
      P=sum1-soldel*soldel/g
               ffsf=soldel*soldel/g
               print*,'delta*delta/g',ffsf
c      print*,'p=sum1',p


      pc=-1.d0/2.d0*deltil*deltil*gtilam
      print*,'deltil,gtilam',deltil,gtilam
      dp=p-pc
       print*,'dp p pc',dp,p,pc
      cim=0
      do i=1,nrimpare
      cim=cim+dsqrt((eimpar(i)-sollam)**2+soldel**2)
      enddo
      dp=dp+cim
      return
      end



      double precision function bcs1nou(d)
      implicit double precision (a-h,o-z)


        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/g,nrnivc,iier
       common/deltaaa/del,fermi
      external bcs2nou

      del=d
      xmed=eonr(nrnive+1-nrocup)

      super=1000.d0
      
      sodelt=d
     
     

      xmin1=xmed
      xmax1=xmed
      
      
1000  continue
      xmin2=xmin1-.001d0
      vmin2=bcs2nou(xmin2)
      xmax2=xmax1+.001d0
      vmax2=bcs2nou(xmax2)
      prod=vmin2*vmax2
      if(prod.le.0.d0)goto 2000
      xmin1=xmin2
      xmax1=xmax2
      if(dabs(xmax1-xmin1).gt.5.d0)goto 2000
      goto 1000
2000  continue

      call ddrtmi(xx,ff,bcs2nou,xmin2,xmax2,1.e-4,80,ier)
      al=xx
c daca iier transmis prin common in bcs1 este zero atunci solutia nu e buna
      iier=ier
 

      fermi=al
      sum=0.d0
      do j=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+j
 
      ene1=eonr(n1)-al
   
      eps1=1.d0/dsqrt(ene1*ene1+d*d)
  
      sum=sum+eps1
      enddo
      ginv=2.d0/g
      zzz=2.d0/g-sum
c      print*,'bcs1',zzz,'lambda',all,'delta',d
      bcs1nou=2.d0/g-sum
      return
      end
      
      double precision function bcs2nou(al)
      implicit double precision (a-h,o-z)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/gf,nrnivc,iier
       common/deltaaa/del,fmmm
      d=del
      sum=0.d0
      do j=-nrnivc,nrnivc-1
      n1=nrnive+1-nrocup+j

      ene1=eonr(n1)-al
      
      eps1=1.d0/dsqrt(ene1*ene1+d*d)
      sum=sum+ene1*eps1

      enddo
      bcs2nou=sum
      return
      end
      

      double precision function gtildanou(epsilo)
c calculeaza densitatea de nivele la energia epsilo
      implicit double precision (a-h,o-z)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/gf,nrnivc,iier
       common/deltaaa/del,fmmm
      g=0
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnive+1-nrnivmaxim,nrnive!1,nrnivmaxim
      gamms=gamma
      xi=(epsilo-eonr(j))/gamms
c     if(dabs(xi).gt.3.d0)goto 1
      yi=xi*xi
      g=g+dexp(-yi)*(2.1875d0+yi*(yi*(1.75d0-yi/6.d0)-4.375d0))/gamms
c1     continue
      enddo
      gtildanou=g/(1.77245385)
      return
      end



      double precision function partnrnou(alambd)
c conservarea numarului de particule pentru obtinerea energiei Fermi
c mediata
      implicit double precision (a-h,o-z)
      dimension hp(50)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/gff,nrnivc,iier
c deoarece nivelele sunt in ordine descrescatoare, sumarea se va face
c de la ultimul nivel pana ce se bucleaza toate particulele (nrpart)
      pi=3.14159265358979324d0
      spi=dsqrt(pi)
      anr=0
      m=3
      mm=m+m
      nrprtm=0
c     print*,'in function partn'
c     print*,'nrnis,nrniv,nrpart',nrnis,nrniv,nrpart
c     print*,'penultim nivel ocupat,ultim nivel si primul liber'
  
        anr=0
        nrprtm=0
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnive+1-nrnivmaxim,nrnive!1,nrnivmaxim
      i=nrnive+1-j
      gammj=gamma
            
      xif=(alambd-eonr(j))/gammj
          
      if(i.le.nrocup)nrprtm=nrprtm+2
c     if(xif.lt.-6.d0)goto 11
      sss=1.d0+erfu(xif)
      a2k=1.d0
      call dhep(hp,xif,mm)
      suma2k=0.d0
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      anr=anr+sss+suma2k
c11    continue
      enddo

      partnrnou=nrprtm-anr

      return
      end
      

      double precision function partnnou(alambd)
c conservarea numarului de particule pentru obtinerea energiei Fermi
c mediata
      implicit double precision (a-h,o-z)
      dimension hp(50)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
c deoarece nivelele sunt in ordine descrescatoare, sumarea se va face
c de la ultimul nivel pana ce se bucleaza toate particulele (nrpart)
      pi=3.14159265358979324d0
      spi=dsqrt(pi)
      anr=0
      m=3
      mm=m+m
      nrprtm=0
c     print*,'in function partn'
c     print*,'nrnis,nrniv,nrpart',nrnis,nrniv,nrpart
c     print*,'penultim nivel ocupat,ultim nivel si primul liber'
                
        anr=0
        nrprtm=0
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnive+1-nrnivmaxim,nrnive!1,nrnivmaxim
      i=nrnive+1-j
      gammj=gamma
            
      xif=(alambd-eonr(j))/gammj
          
      if(i.le.nrocup)nrprtm=nrprtm+2
c     if(xif.lt.-6.d0)goto 11
      sss=1.d0+erfu(xif)
      a2k=1.d0
      call dhep(hp,xif,mm)
      suma2k=0.d0
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      anr=anr+sss+suma2k
c11    continue
      enddo

      partnnou=nrprtm-anr


c          PARTEA IMPARA
            if(nrimpare.eq.1)then
        anr=0
      do j=1,nrimpare
      gammj=gamma
      xif=(alambd-eimpar(j))/gammj
      if(xif.lt.-6.d0)goto 111
      sss=1.d0+erfu(xif)
      a2k=1.d0
      call dhep(hp,xif,mm)
      suma2k=0.d0
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      anr=anr+sss+suma2k
111    continue
      enddo

      partnnou=partnnou+nrimpare-anr
         endif

      return
      end
     



      
      subroutine ginterac(a0,m,nrniv,nrocup,nrnivc,hw,gval)
c daca m este 1 atunci avem sistem cu o pereche sparta
c nrnivc este numar de nivele deasupra si sub nivel fermi
      implicit double precision (a-h,o-z)
      dimension ess(2925)
c prin urmatorul commom se transmit nivele selctionate
c pentru imperechere (in jurul fermi) (pana la o suta) pentru diferite
c configuratii de excitare (pana la 21 configuratii) precum si cele 2
c nivele cu excitare pentru fiecare configuratie
      common/enmoment/emom(21,325),eim(2,21)
      common/conspartnr/eonr(2925),nrnive,nrocupa
      hw0=hw
      nrnive=nrniv
      nrocupa=nrocup
      if(m.eq.1)nrnive=nrniv-2

      do i=1,nrniv
      ess(i)=emom(1,i)       
      enddo
c     nrocup=46
      nrnive=nrniv
      nrperechi=nrnivc*2
      call  ordoneza(ess,nrniv)
      call phenomen(a0,ess,nrperechi,gval,hw0)
      return
      end



      subroutine phenomen(a0,ess,nrperechi,ggg,hw)
c calculeaza g interactie pentru pairing
      implicit double precision (a-h,o-z)
      dimension ess(2925)
      common/conspartnr/eonr(2925),nrnive,nrocup
c     common/gamnasi/gamma,hw0
      common/gamnas/gamma,hw0
      external conspart         
c     hw0=hw      
c      hw0=41.*a0**(-.33333333) 
      do i=1,nrnive    
      eonr(i)=ess(i)/hw0
      enddo
c     gammas=1.4
c     gamma=gammas
      nrpart=nrocup
      cdet=12.d0
            
      do i=1,nrnive-1
      j=i 
      jk1=j
      xli=eonr(jk1)
      jk2=j+1
      xri=eonr(jk2)
      gsd=conspart(xli)*conspart(xri)
      if(gsd.le.0.d0)goto 142
      enddo
142    continue
      call drtmi(xx,ff,conspart,xli,xri,1.e-4,80,ier)
      if(ier.ne.0)then
      endif
       alambd=xx


c    
c    
      gapm=cdet/dsqrt(a0)
      gtilam=gtildanou(alambd)/hw0
      c=nrperechi/2./gtilam/gapm
      c=gtilam*dlog(dsqrt(c**2+1)+c)
      ggg=1.12/c
      return
      end


      double precision function gtildnou(epsilo)
c calculeaza densitatea de nivele la energia epsilo
      implicit double precision (a-h,o-z)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw0
       common/nrnivc/gf,nrnivc,iier
       common/deltaaa/del,fmmm
      g=0
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnive+1-nrnivmaxim,nrnive!1,nrnivmaxim
      gamms=gamma
      xi=(epsilo-eonr(j))/gamms
c     if(dabs(xi).gt.3.d0)goto 1
      yi=xi*xi
      g=g+dexp(-yi)*(2.1875d0+yi*(yi*(1.75d0-yi/6.d0)-4.375d0))/gamms
c1     continue
      enddo
      gtildnou=g/(1.77245385)
      return
      end



!      double precision function gtildnou(epsilo)
!c calculeaza densitatea de nivele la energia epsilo
!      implicit double precision (a-h,o-z)
!      common/conspartnr/eonr(2925),nrnive,nrocup
!      common/gamnasi/gamma,hw0
!      g=0
!      do i=1,nrnive
!      gamms=gamma
!      xi=(epsilo-eonr(i))/gamms
!      if(dabs(xi).gt.3.d0)goto 1
!      yi=xi*xi
!      g=g+dexp(-yi)*(2.1875d0+yi*(yi*(1.75d0-yi/6.d0)-4.375d0))/gamms
!1     continue
!      enddo
!      gtildnou=g/(1.77245385)
!      return
!      end



      double precision function conspart(alambd)
c conservarea numarului de particule pentru obtinerea energiei Fermi
c mediata
      implicit double precision (a-h,o-z)
      dimension hp(50)
      common/conspartnr/eonr(2925),nrnive,nrocup
c     common/gamnasi/gamma,hw0
      common/gamnas/gamma,hw0
c deoarece nivelele sunt in ordine descrescatoare, sumarea se va face
c de la ultimul nivel pana ce se bucleaza toate particulele (nrpart)
      pi=3.14159265358979324d0
      spi=dsqrt(pi)
      anr=0
      m=3
      mm=m+m
      nrprtm=0
        anr=0
        nrprtm=0
c          print*,'nrnive,nrocup,gamma,hw0 in conspart'
c          print*,nrnive,nrocup,gamma,hw0 
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnive+1-nrnivmaxim,nrnive!1,nrnivmaxim
c     do j=1,nrnive
      i=nrnive+1-j
      gammj=gamma
      xif=(alambd-eonr(i))/gammj
      if(i.le.nrocup)nrprtm=nrprtm+2
c     if(xif.lt.-6.d0)goto 11
      sss=1.d0+erfu(xif)
      a2k=1.d0
      call dhep(hp,xif,mm)
      suma2k=0.d0
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      anr=anr+sss+suma2k
c11    continue
      enddo
      conspart=nrprtm-anr
      return
      end

      subroutine lplazv(x1x,a1a,nt,z,g,nd)
      implicit double precision (a-h,o-z)
c linear interpolation
c nt mumar de noduri
c x1x valoarea nodurilor
c a1a valoarea in noduri
c z valorile in care se interpoleaza
c g valorile interpolate
c derg derivata in punctele interpolate
c numar de puncte interpolate

      dimension x1x(1000),a1a(1000),
     c   z(1),g(1)
      do l=1,nd
      v=z(l)
      if(v.le.x1x(1))then
      g(l)=a1a(1)+(a1a(2)-a1a(1))*(x1x(1)-v)/(x1x(2)-x1x(1))
      else
      if(v.ge.x1x(nt))then
      g(l)=a1a(nt)+(a1a(nt)-a1a(nt-1))*(v-x1x(nt))/(x1x(nt)-x1x(nt-1))
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
      endif
      endif
      enddo
      return
      end     
 
      
      subroutine emedienou(idpsn,a0,z0,nrnivs,gammas,es,hw0,
     c       nrimpa,eimp,
     c       u,se,du,alambd)
c calculeaza energia medie 
c idpsn=1 pt neutroni si 0 pentru protoni
c a0
c nrnivs numar de nivele total < 560
c gammas constanta gamma
c es(2925) nivele in MeV
c hw0 distanta intre paturi pentru parinte
c nrimpa nr de nucleoni impari<10
c eimp(10) energiile nucleoni impari
c alambd energia Fermi mediata
c se suma nivele uniparticula
c du efecte de paturi
c ut set dut se refera la corectii fara excitari
      implicit double precision (a-h,o-z)
      dimension es(2925),hp(50)
      dimension eimp(10),eo(2925)
        common/partnrs/eonr(2925),eimpar(10),nrnive,nrimpare,nrocup
      common/gamnas/gamma,hw01
      external partnnou
                 hw0=41.d0*(a0)**(-1.d0/3.d0)
                 hw01=hw0
                 print*,'!!!!!!!!!!!! hw0',hw0
         do i=1,nrimpa
         eimpar(i)=eimp(i)/hw0
         write(92,*)'eimp nrimpa',eimp(i),nrimpa
         enddo
         nrnive=nrnivs
         nrimpare=nrimpa     
         gamma=gammas


c      if(idpsn.eq.1)then
c      nrocup=((a0-z0)-nrimpa)/2
c      else
c      nrocup=(z0-nrimpa)/2
c      endif

      if(idpsn.eq.1)then
      nrocup=(a0-z0)/2
      else
      nrocup=z0/2
      endif

      do i=1,nrnivs
      eonr(i)=es(i)/hw0
      enddo

                  !      print*,'AAAAAAAAAici'

            nvalin=nrocup/2
      do i=nvalin,nrnivs-1
      jk1=nrnivs-i+1 
      xli=eonr(jk1)
      jk2=jk1-1
      xri=eonr(jk2)
      gsd=partnnou(xli)*partnnou(xri)
      if(gsd.le.0.d0)goto 14
      enddo
                      


14    continue

                   !     print*,'AAAAAAAAAici'


c      print*,'in emedie xli,xri',xli,xri
      call drtmi(xx,ff,partnnou,xli,xri,1.e-4,80,ier)
c     print*,'solutia drtmi la emedie xx,ff,ier',xx,ff,ier
      if(ier.ne.0)then
      print*,'ier nu e zero cand calculam lambda tilda'
      print*,'ier=',ier,'  xx=',xx,' ff=',ff,'  xli=',xli,'  xri=',xri
      endif



      alambd=xx

c deoarece nivelele sunt in ordine descrescatoare, sumarea se va face
c de la ultimul nivel pana ce se bucleaza toate particulele (nrpart)
      pi=3.14159265358979324d0
      spi=dsqrt(pi)
      eps1=0.d0
      m=3
      se=0
      mm=m+m
           nrnivmaxim=2*nrocup+nrocup/2
      do j=nrnivs+1-nrnivmaxim,nrnivs!1,nrnivmaxim
      i=nrnivs+1-j
      gammj=gamma
      xif=(alambd-eonr(j))/gammj
      if(i.le.nrocup)se=se+2.d0*eonr(j)
c      if(xif.lt.-6.d0)goto 11
      sss=1.d0+erfu(xif)
      a2k=1.d0
      suma2k=0.d0
      call dhep(hp,xif,mm)
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      eps1=eps1+(sss+suma2k)*eonr(j)
c11    continue
      enddo
c  
      eps2=0.d0
      do j=nrnivs+1-nrnivmaxim,nrnivs!1,nrnivmaxim
      i=nrnivs+1-j
      gammj=gamma
      xif=(alambd-eonr(j))/gammj
c      if(xif.lt.-6.d0) goto 12
      sss=-dexp(-xif*xif)/2.d0
      a2k=1.d0
      suma2k=0.d0
      call dhep(hp,xif,mm)
      do k=1,m
      n2km1=k+k
      n2km2=n2km1-1
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*(xif*hp(n2km1)+hp(n2km2))
      enddo
      eps2=eps2+(sss+suma2k)*gammj
c12    continue
      enddo
      eps2=2.d0/spi*eps2
c    
      u=eps1+eps2


      du=se-u

c             print*,'eps1,eps2,se',eps1,eps2,se
      if(nrimpare.eq.1)then
          print*,'calc nrinpare'
                 se=0
                 do j=1,nrimpare
                 se=se+eimpar(j)
                 write(92,*)'eimpar',eimpar(j),j
                 enddo

      eps1=0
      mm=m+m

      do j=1,nrimpare
      gammj=gamma
      xif=(alambd-eimpar(j))/gammj
      if(xif.lt.-6.d0)goto 111
      sss=1.d0+erfu(xif)
      a2k=1.d0
      suma2k=0.d0
      call dhep(hp,xif,mm)
      do k=1,m
      n2km1=2*k
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*hp(n2km1)
      enddo
      suma2k=2.d0/spi*suma2k
      eps1=eps1+0.5*(sss+suma2k)*eimpar(j)
111    continue
      enddo
c  
      eps2=0.d0
      do j=1,nrimpare
      gammj=gamma
      xif=(alambd-eimpar(j))/gammj
      if(xif.lt.-6.d0) goto 121
      sss=-dexp(-xif*xif)/2.d0
      a2k=1.d0
      suma2k=0.d0
      call dhep(hp,xif,mm)
      do k=1,m
      n2km1=k+k
      n2km2=n2km1-1
      a2k=-a2k/(4.d0*k)
      suma2k=suma2k-a2k*dexp(-xif*xif)*(xif*hp(n2km1)+hp(n2km2))
      enddo
      eps2=eps2+(sss+suma2k)*gammj
121    continue
      enddo
      eps2=1.d0/spi*eps2

           du=du+se-2*(eps1+eps2)
       write(92,*)se,eps1,eps2


      endif


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


      subroutine lplazvd(x1x,a1a,nt,z,g,nd,dg)
c linear interpolation
c nt mumar de noduri
c x1x valoarea nodurilor
c a1a valoarea in noduri
c z valorile in care se interpoleaza
c g valorile interpolate
c derg derivata in punctele interpolate
c numar de puncte interpolate
      implicit double precision (a-h,o-z)
      dimension x1x(200),a1a(200),
     c   z(32),g(32),dg(200)
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



      subroutine descrestere(enrg,numar,nelmax,nuint,einte)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension enrg(2925),numar(2925),eint(2925),nuint(2925)
      dimension einte(2925)
      nrniv=2*nelmax
c le pun in ordine
       do j=1,nrniv
       aa=enrg(j)
       na=numar(j)
       do i=j+1,nrniv
       if(aa.lt.enrg(i))then
       bb=enrg(i)
       enrg(i)=aa
       aa=bb
       nb=numar(i)
       numar(i)=na
       na=nb
       endif
       enddo
       enrg(j)=aa
       numar(j)=na
       enddo
       do j=1,nelmax
       nuint(j)=numar(j+nelmax)
       einte(j)=enrg(j+nelmax)
       enddo
      return
      end




      subroutine corectii2(a0,z0,gms,en,nrspin,iz,nf,nrniv,du,dp,dv)
      implicit double precision (a-h,o-z)
      dimension en(2925),enp(2925),nrspin(2925),infin(2925),
     c   nrspinp(2925),nordin(2925)
      dimension infinp(2925),eimp(10)
! numere asimptotice
      common/paramel/a1px,b1px,a2px,b2px,c1x,c2x,c3x,x1x,x2x,r3x,u1x,u2x
     c               ,ro3x,spx,deltax
      common/nrasimptotice/nrasimptp(2925),nrasimptn(2925)
           do iop=1,nrniv
       if(iz.eq.0)then
           infinp(iop)=nrasimptp(iop)
       else
           infinp(iop)=nrasimptn(iop)
       endif

!!!! modific numarul corespunzator fragment alfa astfel incat sa tina
!cont de intersectii de nivele interzise 
      if(iz.eq.0.and.deltax.lt.14.9d0)then
      infinp(88)=2 !numere identificate cu ajutorul graficului
      infinp(96)=1
      endif 



           enp(iop)=1000.
           enddo
           io=0
           do iop=1,nrniv 
           if(infinp(iop).eq.nf)then
           io=io+1
           enp(io)=en(iop)
           nrspinp(io)=nrspin(iop)
           endif
           enddo
            call donare(nrniv,enp,infinp,nrspinp,nordin)
      print*,'iesit din donare'
            if(iz.eq.0)then
            nnfer=nrniv+1-(z0+1)/2
            else
            nnfer=nrniv+1-(a0-z0)/2
            endif
            enfermip=en(nnfer)
      idp=iz
      gamma=gms
      hw0=41*a0**(-.33333333333d0)
         nrimpa=0
            print*,'a0',a0,'intru emedienou'
      call emedienou(idp,a0,z0,nrniv,gamma,enp,hw0,
     c       nrimpa,eimp,
     c       u,se,du,alambd)
            print*,'idp,z0,nrniv,gamma,hw0,nrimpa'
            print*,idp,z0,nrniv,gamma,hw0,nrimpa
      call pairnou(idp,a0,z0,enp,
     c   gamma,dp,sollam,soldel)
      du=hw0*du
      dp=hw0*dp
      sollam=sollam*hw0
      soldel=soldel*hw0
      sollamp=sollam
      soldelp=soldel
       print*,'lambda delta in pair',sollamp,soldelp
       print*,'corectii de paturi in MeV du si dp',du,dp
            dupro=du
            dppro=dp
      dv=du+dp
          return
          end


      subroutine corectii3(a0,z0,gms,en,nrspin,iz,nf,nrniv,du,dp,dv)
c fac corectii pe nivelele asimptotice (spre deosebire de corectii2
c unde se fac corectii pe nivelele extrase din schema de nivele
      implicit double precision (a-h,o-z)
      dimension en(2925),enp(2925),nrspin(2925),infin(2925),
     c   nrspinp(2925),nordin(2925)
      dimension infinp(2925),eimp(10)
! numere asimptotice
      common/paramel/a1px,b1px,a2px,b2px,c1x,c2x,c3x,x1x,x2x,r3x,u1x,u2x
     c               ,ro3x,spx,deltax
      common/nrasimptotice/nrasimptp(2925),nrasimptn(2925)
      common/energiiasim/epas1(2925),epas2(2925),enas1(2925),enas2(2925)
c nf = nr fraGMENT 
           io=0
           do iop=1,nrniv 
           io=io+1
           IF(NF.EQ.1.AND.IZ.EQ.0)THEn
           enp(io)=EPAS1(iop)
c             print*,'epas1(iop)',epas1(iop),iop
           endif
           IF(NF.EQ.1.AND.IZ.EQ.1)THEn
           enp(io)=EnAS1(iop)
c             print*,'enas1(iop)',enas1(iop),iop
           endif
           nrspinp(io)=1
           IF(NF.EQ.2.AND.IZ.EQ.0)THEn
           enp(io)=EPAS2(iop)
           endif
           IF(NF.EQ.2.AND.IZ.EQ.1)THEn
           enp(io)=EnAS2(iop)
           endif
           enddo
            call donare(nrniv,enp,infinp,nrspinp,nordin)
      print*,'iesit din donare'
            if(iz.eq.0)then
            nnfer=nrniv+1-(z0+1)/2
            else
            nnfer=nrniv+1-(a0-z0)/2
            endif
            enfermip=en(nnfer)
c      print*,'iz,nf,nnfer,a0,z0',iz,nf,nnfer,a0,z0
c              stop
      idp=iz
      gamma=gms
      hw0=41*a0**(-.33333333333d0)
         nrimpa=0
            print*,'a0',a0,'intru emedienou'
      call emedienou(idp,a0,z0,nrniv,gamma,enp,hw0,
     c       nrimpa,eimp,
     c       u,se,du,alambd)
            print*,'idp,z0,nrniv,gamma,hw0,nrimpa'
            print*,idp,z0,nrniv,gamma,hw0,nrimpa
      call pairnou(idp,a0,z0,enp,
     c   gamma,dp,sollam,soldel)
      du=hw0*du
      dp=hw0*dp
      sollam=sollam*hw0
      soldel=soldel*hw0
      sollamp=sollam
      soldelp=soldel
       print*,'lambda delta in pair',sollamp,soldelp
       print*,'corectii de paturi in MeV du si dp',du,dp
            dupro=du
            dppro=dp
      dv=du+dp
          return
          end



      subroutine psciziune(a1,b1,a2,b2,cr3,dltsciz,ier1,ier2,ier3,ier4)
c calculam elongatia de sciziune
      implicit double precision (a-h,o-z)
      common/ptsc/a,b,r33,dd,iier
      common/lasciziune/delatsciz
      external sciziune
       r3=50
       ier1=0
       ier2=0
       ier3=0
       ier4=0
      dltsciz=100.
      deltatsciz=dltsciz
      if(cr3.le.0.02)return
      if(cr3.gt.0.02)r3=1.d0/cr3
c     print*,'calculez sciziunea cu a1,b1,a2,b2,r3'
c     print*,a1,b1,a2,b2,r3
c calculez pentru fragment 1
      a=a1
      b=b1
      r33=r3
       ef2=r3+a1
       ef1=a1-1
        if(ef1.lt.0.d0)ef1=0.d0
          dife=(ef2-ef1)/50.  
  ! caut delta
        ef11=ef1
        vvv1=sciziune(ef11)
        valmin=dabs(vvv1)
        efmin=ef1
       do ii=1,50
        ef11=ef1+dife*(ii-1)
        ef22=ef11+dife
        vvv2=sciziune(ef22)
        if(valmin.gt.dabs(vvv2))then
        valmin=dabs(vvv2)
        efmin=vvv2
        endif
        prod=vvv1*vvv2
c       print*,'ef11,ef22,vvv1,vvv2',ef11,ef22,vvv1,vvv2
        vvv1=vvv2
        if(prod.le.0.d0)goto 3
        enddo
        x1=efmin
        x1=a1+.5*r3
              goto 4
3     continue
      call drtmi(x1,v,sciziune,ef11,ef22,1.e-3,100,ier)
c          print*,'ier drtmi (1)',ier
       ier1=ier
       ier2=iier
4      continue
c calculez pentru fragment 2
      a=a2
      b=b2
      r33=r3
       ef2=r3+a2
       ef1=a2-1
        if(ef2.lt.0.d0)ef2=0.d0
          dife=(ef2-ef1)/50.
  ! caut delta
        ef11=ef1
        vvv1=sciziune(ef11)
        valmin=dabs(vvv1)
        efmin=ef1
       do ii=1,50
        ef11=ef1+dife*(ii-1)
        ef22=ef11+dife
        vvv2=sciziune(ef22)
        if(valmin.gt.dabs(vvv2))then
        valmin=dabs(vvv2)
        efmin=vvv2
        endif
        prod=vvv1*vvv2
c       print*,'ef11,ef22,vvv1,vvv2',ef11,ef22,vvv1,vvv2
        vvv1=vvv2
        if(prod.le.0.d0)goto 33
        enddo
        x2=efmin
        x2=a2+.5*r3
              goto 44
33     continue
      call drtmi(x2,v,sciziune,ef1,ef2,1.e-3,100,ier)
      ier3=ier
      ier4=iier
44    continue
c          print*,'ier drtmi (2)',ier
      dltsciz=x1+x2
      delatsciz=dltsciz
      return
      end
      double precision function sciziune(d)
      implicit double precision (a-h,o-z)
      common/ptsc/a,b,r33,dd,iier
      external sciziune2
       dd=d
       ef1=0.0001
       ef2=d
c rezolvam derivatele
      call ddrtmi(x,v,sciziune2,ef1,ef2,1.e-3,100,ier)
c          print*,'ier2 ddrtmi',ier
           iier=ier
      zz=x
      val1=zz/a
      if(val1.gt.0.9999d0)val1=0.9999d0
      ro1=b*dsqrt(1.d0-val1**2)
      val2=z-dd
      r33min=r33-1.d-4
      if(dabs(val2).gt.r33min)val2=r33min    
      ro2=r33-dsqrt(r33**2-val2**2)
      sciziune=ro1-ro2
c         print*,'d,sciziune,ro1,ro2',d,sciziune,ro1,ro2
      return
      end
      double precision function sciziune2(z)
      implicit double precision (a-h,o-z)
      common/ptsc/a,b,r33,dd,iier
      val1=z/a
c      print*,'z,a,val1',z,a,val1
      if(val1.gt.0.9999d0)val1=0.9999
      if(val1.lt.-0.9999d0)val1=-0.9999
c      print*,'val1',val1
      rad1=(dsqrt(1.d0-val1**2))**(3./2.)
c      print*,'rad1',rad1
      der1=-b*val1/rad1/a**2
      val2=z-dd
      r33min=r33-1.d-7
      if(val2.gt.r33min)val2=r33min    
      if(val2.lt.-r33min)val2=-r33min    
      rad2=(dsqrt(r33**2-val2**2))**(3./2.)
      der2=val2/rad2
      sciziune2=der1-der2
c      print*,'sciziune2,der1,der2',sciziune2,der1,der2
c      print*,'val1,val2,rad1,rad2',val1,va2,rad1,rad2
      return
      end 
 



      double precision function csi(z,ro,nmax,iomg,ll)
      implicit double precision (a-h,o-z)
ccc Calculez modul patrat din functia de unda cu numar de proiectie a spin
ccc iomg si numar nivel in  spatiul proiectiei de spin ll

   
      common/nrcuantice/numcz1(25,325),numcm1(25,325),numcro1(25,325),
     c          numcsp1(25,325),
     c          numcz2(25,325),numcm2(25,325),numcro2(25,325),
     c          numcsp2(25,325)
      common/vectoripro/VECPRO(25,325,325),LELMAX(25)
      COMMON/ZET12/ZET1,ZET2
      COMMON/ALF12/ALF1,ALF2
      common/numecuantce/X1D(0:24),X2D(0:24),CO1D(0:24),CO2D(0:24)
      common/ckst/ckst
      common/r116/R0
      common/paramel/a1x,b1,a2,b2,cc1,cc2,c3,x1,x2,r3,u1,u2,ro3,s,delt
    

          csi=0.d0
          densit1=0
          densit2=0
      nelmax=(nmax+1-iomg)*(nmax+2-iomg)/2
      do l=1,nelmax
c calculeaza pla=L_{n}^{k}(x) (polinom laguerre asociat)
cccc constanta de normare laguerre
      fact1=factorial(numcro1(iomg,l))
      fact2=factorial(numcro1(iomg,l)+numcm1(iomg,l))
      cl=dsqrt(2.d0*fact1/fact2)
      H=6.583021D-22
      MPH=1.58823D7
      W0=ckst/H
      W1R0=W0*R0/b1
      alf1ro=DSQRT(MPH*W1r0)*1.D-15
      rho=(alf1ro*ro)**2
      n=numcro1(iomg,l)
      k=iabs(numcm1(iomg,l))
      call laguerre(rho,n,k,pla1)
c        print*,'rho,n,k,pla1',rho,n,k,pla1
      rnm=cl*alf1ro*dexp(-rho/2.)*(rho)**(k/2.)*pla1
        v1=dexp(-rho/2.)
        v2=(rho)**(k/2.)
c      print*,'cl,alf1ro,rho,v1,v2',cl,alf1ro,rho,v1,v2
c         print*,'rnm,iomg,l,ll',rnm,iomg,l,ll
cccc impart la radical din pi pentru funtia unghiulara
      rnm=rnm/dsqrt(2.d0*3.14)
cccc calculez partea dupa z
          if(z.lt.0)then
          zz=alf1*z+dabs(zet1)
          anr1=x1d(numcz1(iomg,l))
          call dherm(anr1,zz,h1,d1)
c              print*,'anr1,zz,h1,d1',anr1,zz,h1,d1
          fz=co1d(numcz1(iomg,l))*dexp(-zz**2/2.)*h1
          else
          zz=-(z*alf2-zet2)
          anr2=x2d(numcz1(iomg,l))
          call dherm(anr2,zz,h2,d2)
c              print*,'anr2,zz,h2,d2',anr1,zz,h1,d1
          fz=co2d(numcz1(iomg,l))*dexp(-zz**2/2.)*h2
          endif
       rnmf=rnm*fz*vecpro(iomg,l,ll)
c         print*,'v(iomg,ll)*rnm*fz*vecpro(iomg,l,ll)',v(iomg,ll),rnm,fz
c    c ,vecpro(iomg,l,ll)
c         print*,'rnmf,iomg,l,ll',rnmf,iomg,l,ll
                 if(numcsp1(iomg,l).eq.-1)then
                 densit1=densit1+rnmf
                 else
                 densit2=densit2+rnmf
                 endif
c         print*,'densit',densit
      enddo
       csi=densit1**2+densit2**2

      return
      end


      subroutine intetrip(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1n,iomg2n,ll1n,ll2n,zp1,zp2,zt)
      implicit double precision (a-h,o-z)
c integrez densitatea electrica pe volum


      common/numciogll/iomg1,iomg2,ll1,ll2
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/asiz/s1,s2
      common/densitup/dens(2,100,100),n1,n2,n3
      common/distanta/distanta
      external fuctic3

       iomg1=iomg1n
       iomg2=iomg2n
       ll1=ll1n
       ll2=ll2n
! pentru inceput redefinesc coordonatele functie de originea potentialului cu 
! doua centre
      CT=a2/dabs(a1)
      cx1=-DELT/(1.D0+CT)
      cx2=CT*dabs(cx1)
C DISTANTA DINTRE PLANUL DE INTERSECTIE A PARAMETRIZARII SEMISIMETRICE  SI
C PLANUL DE INTERSECTIE A PARAMTERIZARII ASIMETRICE
      Z0P=cx1-c1

      c1i=c1+z0p
      c2i=c2+z0p
      c3i=c3+z0p
      x1pi=x1p+z0p
      x2pi=x2p+z0p
      x3pi=0.
      u1i=u1+z0p
      u2i=u2+z0p
      sollami=sollam
      soldeli=soldel
      sollami1=sollam1
      soldeli1=soldel1
      sollami2=sollam2
      soldeli2=soldel2
      nmaxi=nmax
      n2=0 
      if(c3i.gt.c2i.or.c3i.lt.c1i)then
      alim1=c1i-a1-2
      alim2=x2pi
      alim3=c2i+a2+2
c             print*,'0'
      call dqg32(alim1,alim2,fuctic3,zp1)
      call dqg32(alim2,alim3,fuctic3,zp2)
      zt=zp1+zp2 
      else
      if(u1i.ge.u2i)then
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=c2i+a2+2
c          print*,'1'
      call dqg32(alim1,alim2,fuctic3,zp1)
      call dqg32(alim2,alim3,fuctic3,zp2)
      zt=zp1+zp2 
      else
c           print*,'2'
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=x3pi
      alim4=c2i+a2+2
      call dqg32(alim1,alim2,fuctic3,zp1)
      call dqg32(alim3,alim4,fuctic3,zp2)
      zt=zp1+zp2 
      endif
      endif
        pi2=dasin(1.d0)
c        print*,'pi2',pi2
        zt=4*pi2*zt*2
        zp1=4*pi2*zp1*2
        zp2=4*pi2*zp2*2
c         write(96,*)s1,s2,' z1 z2'
      return
      end


      double precision function fuctic3(z)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxa
      common/densitup/dens(2,100,100),n1,n2,n3
      external fuctic23
      n2=n2+1
      n3=0
        zzzz=z
      if(z.le.x1pi)then
      rad=(1.-((z-c1i)/a1)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b1*dsqrt(rad)                
      else
      if(z.ge.x2pi)then
      rad=(1.-((z-c2i)/a2)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b2*dsqrt(rad)                   
      else
      rad=r3**2-(z-c3i)**2
      if(rad.lt.1.e-8)rad=1.e-8
      rho=ro3-s*dsqrt(rad)
      if(rho.le.0.d0)rho=0.
      endif
      endif
       alim=rho+3.
       call dqg322(0.d0,alim,fuctic23,y)
       fuctic3=y
c       print*,'cccccccccccccccccccccccccccccccc'
c       print*,fuctic3
      return
      end


      double precision function fuctic23(rh)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/numciogll/iomg1,iomg2,ll1,ll2
      common/asiz/s1,s2
        common/denstt/sollam,soldel,sollam1,soldl1,rrsc,nmaxy
      common/densitup/dens(2,100,100),n1,n2,n3
      common/distanta/distanta
      n3=n3+1
      z=zzzz
       djos=dsqrt((distanta-z)**2+rh**2)
      ro=rh
      sollam=sollami
      soldel=soldeli
      sollam1=sollami1
      soldel1=soldeli1
      sollam2=sollami2
      soldel2=soldeli2
      nmax=nmaxi
      v1=csi(z,ro,nmax,iomg1,ll1)
      v2=csi(z,ro,nmax,iomg2,ll2)
c     pim=rhodens(z,ro)
        sollamm=sollam
        soldelm=soldel
        sollam1m=sollam1
        soldel1m=soldel1
        rrscm=rrsc
        nmaxm=nmaxy
c     call densitateuni(z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,
c    c                  nmaxm,densit)
      densit=dens(1,n2,n3)+dens(2,n2,n3)
      vrho=-999.d0
      gam=1.d0
      rhoc=0.16d0
      pim=vrho*(1.d0-(densit/rhoc)**gam)
c       call densitateuni(z,ro,sollam,soldel,sollam1,soldel1,
c     c     sollam2,soldel2,nmax,densit,suma1,suma2)
c            print*,'z,ro,densit',z,ro,densit
c      s1=suma1
c      s2=suma2
c        fuctic23=ro*densit
c       fuctic23=-v1*v2*pim
        fuctic23=dens(1,n2,n3)!*1.44/djos
c       print*,'v1*v2*pim',v1,v2,pim,n1,n2,fuctic23
c       print*,'iomg1,ll1,iomg2,ll2,nmax',iomg1,ll1,iomg2,ll2,nmax
       return
       end
    

 


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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C calculez valorile densitatii uniparticula in nodurile 
c Gauss Legendre

      subroutine intetripid(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1n,iomg2n,ll1n,ll2n,zp1,zp2,zt)
      implicit double precision (a-h,o-z)
c integrez densitatea electrica pe voluma


      common/densitup/dens(2,100,100),n1,n2,n3
c  n1=1 sau 2 pentru protoni sau neutroni
c n1 si n2 array pentru z si rho

      common/numciogll/iomg1,iomg2,ll1,ll2
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/asiz/s1,s2
      external fuctic3d


         print*,'n1,n2,n3 1',n1,n2,n3

       iomg1=iomg1n
       iomg2=iomg2n
       ll1=ll1n
       ll2=ll2n
! pentru inceput redefinesc coordonatele functie de originea potentialului cu 
! doua centre
      CT=a2/dabs(a1)
      cx1=-DELT/(1.D0+CT)
      cx2=CT*dabs(cx1)
C DISTANTA DINTRE PLANUL DE INTERSECTIE A PARAMETRIZARII SEMISIMETRICE  SI
C PLANUL DE INTERSECTIE A PARAMTERIZARII ASIMETRICE
      Z0P=cx1-c1

      c1i=c1+z0p
      c2i=c2+z0p
      c3i=c3+z0p
      x1pi=x1p+z0p
      x2pi=x2p+z0p
      x3pi=0.
      u1i=u1+z0p
      u2i=u2+z0p
      sollami=sollam
      soldeli=soldel
      sollami1=sollam1
      soldeli1=soldel1
      sollami2=sollam2
      soldeli2=soldel2
      nmaxi=nmax
                  n2=0
      if(c3i.gt.c2i.or.c3i.lt.c1i)then
      alim1=c1i-a1-2
      alim2=x2pi
      alim3=c2i+a2+2
c             print*,'0'
      call dqg32(alim1,alim2,fuctic3d,zp1)
      call dqg32(alim2,alim3,fuctic3d,zp2)
      zt=zp1+zp2
      else
      if(u1i.ge.u2i)then
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=c2i+a2+2
c          print*,'1'
      call dqg32(alim1,alim2,fuctic3d,zp1)
      call dqg32(alim2,alim3,fuctic3d,zp2)
      zt=zp1+zp2
      else
c           print*,'2'
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=x3pi
      alim4=c2i+a2+2
      call dqg32(alim1,alim2,fuctic3d,zp1)
      call dqg32(alim3,alim4,fuctic3d,zp2)
      zt=zp1+zp2
      endif
      endif
        pi2=dasin(1.d0)
c        print*,'pi2',pi2
        zt=4*pi2*zt
        zp1=4*pi2*zp1
        zp2=4*pi2*zp2
c         write(96,*)s1,s2,' z1 z2'
      return
      end


      double precision function fuctic3d(z)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxa
      common/densitup/dens(2,100,100),n1,n2,n3
      external fuctic23d
       n2=n2+1
       n3=0
         print*,'n1,n2,n3 1',n1,n2,n3
        zzzz=z
      if(z.le.x1pi)then
      rad=(1.-((z-c1i)/a1)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b1*dsqrt(rad)
      else
      if(z.ge.x2pi)then
      rad=(1.-((z-c2i)/a2)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b2*dsqrt(rad)
      else
      rad=r3**2-(z-c3i)**2
      if(rad.lt.1.e-8)rad=1.e-8
      rho=ro3-s*dsqrt(rad)
      if(rho.le.0.d0)rho=0.
      endif
      endif
       alim=rho+3.
       call dqg322(0.d0,alim,fuctic23d,y)
       fuctic3=y
      return
      end


      double precision function fuctic23d(rh)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/numciogll/iomg1,iomg2,ll1,ll2
      common/asiz/s1,s2
        common/denstt/sollam,soldel,sollam1,soldl1,rrsc,nmaxy
      common/densitup/dens(2,100,100),n1,n2,n3
      n3=n3+1
      z=zzzz
      ro=rh
      sollam=sollami
      soldel=soldeli
      sollam1=sollami1
      soldel1=soldeli1
      sollam2=sollami2
      soldel2=soldeli2

      nmax=nmaxi
c     pim=rhodens(z,ro)
        sollamm=sollam
        soldelm=soldel
        sollam1m=sollam1
        soldel1m=soldel1
        rrscm=rrsc
        nmaxm=nmaxy
c     print*,'z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit'
c     print*,z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit
      call densitateuni(z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,
     c                  nmaxm,densit)
c     print*,z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit
       dens(n1,n2,n3)=2*densit*ro
c              stop

        fuctic23d=2*densit*ro
ccccc       (n1,n2,n3)',dens(n1,n2,n3),n1,n2,n3
      vrho=-999.d0
      gam=1.d0
      rhoc=0.16d0
      pim=vrho*(1.d0-(densit/rhoc)**gam)
   
       return
       end


      subroutine intetripid2(sollam,soldel,sollam1,soldel1,
     c    sollam2,soldel2,nmax,iomg1n,iomg2n,ll1n,ll2n,zp1,zp2,zt)
      implicit double precision (a-h,o-z)
c integrez densitatea electrica pe voluma


      common/densitup/dens(2,100,100),n1,n2,n3
c  n1=1 sau 2 pentru protoni sau neutroni
c n1 si n2 array pentru z si rho

      common/numciogll/iomg1,iomg2,ll1,ll2
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/asiz/s1,s2
      external fuctic3d2


         print*,'n1,n2,n3 1',n1,n2,n3

       iomg1=iomg1n
       iomg2=iomg2n
       ll1=ll1n
       ll2=ll2n
! pentru inceput redefinesc coordonatele functie de originea potentialului cu 
! doua centre
      CT=a2/dabs(a1)
      cx1=-DELT/(1.D0+CT)
      cx2=CT*dabs(cx1)
C DISTANTA DINTRE PLANUL DE INTERSECTIE A PARAMETRIZARII SEMISIMETRICE  SI
C PLANUL DE INTERSECTIE A PARAMTERIZARII ASIMETRICE
      Z0P=cx1-c1

      c1i=c1+z0p
      c2i=c2+z0p
      c3i=c3+z0p
      x1pi=x1p+z0p
      x2pi=x2p+z0p
      x3pi=0.
      u1i=u1+z0p
      u2i=u2+z0p
      sollami=sollam
      soldeli=soldel
      sollami1=sollam1
      soldeli1=soldel1
      sollami2=sollam2
      soldeli2=soldel2
      nmaxi=nmax
                  n2=0
      if(c3i.gt.c2i.or.c3i.lt.c1i)then
      alim1=c1i-a1-2
      alim2=x2pi
      alim3=c2i+a2+2
c             print*,'0'
      call dqg32(alim1,alim2,fuctic3d2,zp1)
      call dqg32(alim2,alim3,fuctic3d2,zp2)
      zt=zp1+zp2
      else
      if(u1i.ge.u2i)then
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=c2i+a2+2
c          print*,'1'
      call dqg32(alim1,alim2,fuctic3d2,zp1)
      call dqg32(alim2,alim3,fuctic3d2,zp2)
      zt=zp1+zp2
      else
c           print*,'2'
      alim1=c1i-a1-2
      alim2=x3pi
      alim3=x3pi
      alim4=c2i+a2+2
      call dqg32(alim1,alim2,fuctic3d2,zp1)
      call dqg32(alim3,alim4,fuctic3d2,zp2)
      zt=zp1+zp2
      endif
      endif
        pi2=dasin(1.d0)
c        print*,'pi2',pi2
        zt=4*pi2*zt
        zp1=4*pi2*zp1
        zp2=4*pi2*zp2
c         write(96,*)s1,s2,' z1 z2'
      return
      end


      double precision function fuctic3d2(z)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxa
      common/densitup/dens(2,100,100),n1,n2,n3
      external fuctic23d2
       n2=n2+1
       n3=0
         print*,'n1,n2,n3 1',n1,n2,n3
        zzzz=z
      if(z.le.x1pi)then
      rad=(1.-((z-c1i)/a1)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b1*dsqrt(rad)
      else
      if(z.ge.x2pi)then
      rad=(1.-((z-c2i)/a2)**2)
      if(rad.lt.1.e-8)rad=1.e-8
      rho=b2*dsqrt(rad)
      else
      rad=r3**2-(z-c3i)**2
      if(rad.lt.1.e-8)rad=1.e-8
      rho=ro3-s*dsqrt(rad)
      if(rho.le.0.d0)rho=0.
      endif
      endif
       alim=rho+3.
       call dqg322(0.d0,alim,fuctic23d2,y)
       fuctic3=y
      return
      end


      double precision function fuctic23d2(rh)
      implicit double precision (a-h,o-z)
      common/paramel/a1,b1,a2,b2,c1,c2,c3,x1p,x2p,r3,u1,u2,ro3,s,delt
      common/integrc1/c1i,c2i,c3i,x1pi,x2pi,u1i,u2i,zzzz
      common/integrc2/sollami,soldeli,nmaxi
      common/integrc3/sollami1,soldeli1,sollami2,soldeli2
      common/numciogll/iomg1,iomg2,ll1,ll2
      common/asiz/s1,s2
        common/denstt/sollam,soldel,sollam1,soldl1,rrsc,nmaxy
      common/densitup/dens(2,100,100),n1,n2,n3
      common/distanta/distanta
      n3=n3+1
      z=zzzz
      ro=rh
      sollam=sollami
      soldel=soldeli
      sollam1=sollami1
      soldel1=soldeli1
      sollam2=sollami2
      soldel2=soldeli2

      nmax=nmaxi
c     pim=rhodens(z,ro)
        sollamm=sollam
        soldelm=soldel
        sollam1m=sollam1
        soldel1m=soldel1
        rrscm=rrsc
        nmaxm=nmaxy
c     print*,'z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit'
c     print*,z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit
      call densitateuni(z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,
     c                  nmaxm,densit)
c     print*,z,ro,sollamm,soldelm,sollam1m,soldel1m,rrscm,nmaxm,densit
       dens(n1,n2,n3)=2*densit*ro
c              stop

        fuctic23d2=2*densit*ro/dsqrt((distanta-z)**2+rh**2)*1.44
ccccc       (n1,n2,n3)',dens(n1,n2,n3),n1,n2,n3
      vrho=-999.d0
      gam=1.d0
      rhoc=0.16d0
      pim=vrho*(1.d0-(densit/rhoc)**gam)
   
       return
       end



