    cc Program written by: Lindsay King
      
      PROGRAM nfwsim
      implicit real*8(a-h,o-z)
      parameter (N=600, cv=3.d5, ho=70.d0,ngriddy=50)
      parameter (gg=6.67e-11, om=0.3d0, ol=0.7d0)
      parameter	(nmax=100000)
      complex*16 ep,gamman
      dimension p(2),xi(2,2),xx(100,100),x(nmax),y(nmax)
      complex*16 epl,epld
      integer nn(2)
      common/pass1/epld(nmax),phild(nmax),thld(nmax),nk
      common/pass2/sigcrit,rhoc,dol,rin,rout
c      common/pass3/sig
      external dangg


      open(unit=33,file="parameters3.dat",form="formatted")

cc ellipticity dispersion

        sig=0.2d0


cc density of galaxies per arcminsq

        endens=300.



c inner and outer radii -- rin and rout in Mpc

         rout=2.d0

         rin=0.05d0


c         rin=0.12d0
c         rout=0.8d0


cc seeds for random number generation - might not all be used
	idum1=-10
	idum2=-30
	idum3=-40
	idum4=-50
	idum5=-33
	idum=-66
        idummm=-66
	

        xxmin=10000000.

	pi=4.*atan(1.)



        zl=0.2d0

        zs=1.0d0


cc       angular diam distances and sig crit


 	call qromb(dangg,0.d0,zl,dd)
 
  	dd=(dd*cv/ho)/(1.d0+zl)

	dol=dd

         call qromb(dangg,0.d0,zs,ds)

    	 ds=(ds*cv/ho)/(1.d0+zs)

         ddn=dd*(1.d0+zl)
	
 	 dsn=ds*(1.d0+zs)

	 dds=(1.d0/(1.d0+zs))*(dsn-ddn)

c        	write(*,*)dd,ds,dds

          dang=ds/(dd*dds)

          sigcrit=(cv**2/(4.d0*pi*GG))*dang



ccc*** function w, lensing efficiency for Miyoung

ccc         if (zs.gt.zl) then 

ccc          weff=dds/ds

ccc          else

ccc            weff=0.

ccc            endif

c          write(*,*)sigcrit, 'sigcrit'

cc rhocrit

         rhoc=(3.d0/(8.d0*pi*GG))*ho*ho*((om*((1.d0+zl)**3))+ol)


cc area of square containing circle radius rout

         rea=4*((rout**2))

cc physical scale of 1 arcmin anglular scale at lens redshift

         physlenlens=dol*(pi/(60.*180.))

c         write(*,*)physlenlens*0.6, physlenlens*4.0

c physical scale at source redshift corresponding to 1 arcmin

c         physlen=ds*(pi/(60.*180.))


cc galaxies per arcminsq to per Mpc^2

         engal=(endens*rea)/(physlenlens**2.)

         ndensity=int(engal/rea)

         engal=int(engal)

c         write(*,*)engal

cc catalogue generation -- 1 catalogue

          do ik=1,100


cc zero out arrays

            do iik=1,nmax
            thld(iik)=0.
            epld(iik)=complex(0.,0.)
            phild(iik)=0.
            x(iik)=0.
            y(iik)=0.
            enddo

	 nk=0

         nmag=0

c Poisson noise on galaxy count

         gal=poidev(engal,idummm)

c         write(*,*)'poisson gal', gal

c make sure it is an integer! ;-)
c Will: Change this to the size of the 
         ngal=int(gal)

c now loop over these galaxies, generating x, y coords
c coords go from -rout to +rout
c th is distance from centre, ph is position angle of
c galaxy location
c Will: ran2 is built in thing in FORTRAN to make random numbers
	do in=1,ngal

         xxx=(-1.*rout)+(ran2(idum)*2*rout)

         yyy=(-1.*rout)+(ran2(idum)*2*rout)

         th=((xxx**2)+(yyy**2))**0.5

         ph = atan2(yyy,xxx)


c generate components of ellipticity
c Will: gasdev is a gaussian distrbitution
  	ep1=gasdev(idum3)*sig/(2.d0**0.5)
	ep2=gasdev(idum2)*sig/(2.d0**0.5)

	ep=complex(ep1,ep2)


c call the nfw routine to get reduced shear, magnification at
c galaxy location -- returns g1,g2, gamma, appa, emu

	call enfw(th,ph,sigcrit,rhoc,g1,g2,gamma,appa,emu)
	do ii=1,1000
        read (99,*) x,y,e1,e2
	
	gamman=complex(g1,g2)


c put limit on reduced shear and location of galaxies to be 
c analysed


c must be outside inner radius
            if (th.ge.rin) then

c must be inside outer radius

               if (th.le.rout) then

c check there is no crazy high magnification
 
            if (emu.lt.3.) then

c check there is no high reduced shear

  	       if (abs(gamman).le.0.5d0) then


c now comes the bit to reject a fraction of galaxies depending on their magnifications - more 
c highly magnified, more rejected
c since space also gets magnified (n=no * mu**(beta-1))

c so generate random number to attach to galaxy         

	     rannum=ran2(idum5)

c work out mu**(beta-1)

   	    emupow=emu**(-0.5)


c rannum is between 0 and 1, so if emupow is bigger there is less chance that random number is bigger than it!
c see e.g. numerical recipes for this so-called rejection technique

           if (rannum.le.emupow) then

              nk=nk+1

c              write(*,*)th,emupow,nk

c then work out the lensed ellipticity using unlensed and reduced shear, and put corresponding phi and r in arrays

               epld(nk)=(ep+gamman)/(1.d0+(conjg(gamman)*ep))

               phild(nk)=ph

               thld(nk)=th

c              write(*,*)nk, ep, gamman


           endif
           endif
           endif
           endif
           endif

         	enddo



cccc end of generation of lensed cats.

c        write(*,*)'out', nk


c now either use a minimization technique or a gridded likelihood

	ftol=1.d-8
	

ccc
cccc uncomment next lines for gridded stuff
cc        cmin=5.5
cc        cmax=6.5
cc        r200min=1.65
cc        r200max=1.85
cc        deec=(cmax-cmin)/(ngriddy-1.)
cc        deer=(r200max-r200min)/(ngriddy-1.)
cc        do i=1,ngriddy
cc           p(2)=cmin + (i*deec)
cc        do jj=1,ngriddy
cc           p(1)= r200min + (jj*deer)


cccc parameters for powell minimization

c	p(2)=4.
c	p(1)=1.
c	xi(1,1)=0.5
c	xi(2,2)=0.5
c	xi(1,2)=0.000001
c	xi(2,1)=0.000001


	p(2)=4.
	p(1)=1.
	xi(1,1)=0.5
	xi(2,2)=0.5
	xi(1,2)=0.1
	xi(2,1)=0.1

cc        write(*,*) 'Using Max Like'

cc this goes and calculates the function func, which powell is going to minimize
cc comment this line if you don't want to run powell!

	call powell(p,xi,2,2,ftol,iter,fret)


	egaparsec=3.08568025d22
	fac=1.d6/gg
	disfac=1.d6*egaparsec
	cosfac=(om*((1.d0+zl)**3))+ol

c just getting mass if needed

	eM=rhoc*(800./3.)*pi*(p(1)**3)*egaparsec*1.d6/1.9889d45

	write(*,*)ik, p(1),p(2)

	write(33,*)p(1),p(2)

c        write(*,*)'here'


cc uncomment for gridded stuff

c         xxval=func(p)

c         xxmin=min(xxval,xxmin)

c         xx(i,jj)=xxval

c         if (xxval.eq.xxmin) then
            
c            p1min=p(1)
c            p2min=p(2)

c            endif

c  	write(99,*)p(2),p(1),xxmin


ccc uncomment these two enddos to do likelihood on a grid

c	enddo
c        enddo

c        write(*,*)p1min,p2min

c        do iii=1,ngriddy
c           do jjj=1,ngriddy
c              write(99,*)jjj,iii,xx(iii,jjj), 2.*(xx(iii,jjj)-xxmin)
c              enddo

              enddo




      STOP
      END


cccc this subroutine is calculating the nfw pararameters for the lens at the 
cc location of the galaxy, like shear and so on


      subroutine enfw(th,ph,sigcrit,rhoc,g1,g2,gamma,appa,emu)
      implicit real*8(a-h,o-z)
      parameter (PI=3.1415926536d0,ho=70.d0,gg=6.67e-11,cv=3.e5)
      parameter (om=0.3d0,ol=0.7d0,n=600)      
c      external dangg, func

ccc can change to get different cluster
cc nfw concentration and r_200 in Mpc

  	c=4.d0

	r_200=1.5d0

        deltac=(200.d0/3.d0)*(c**3)/(log(1.+c)-(c/(1.d0+c)))
	  
        r=th

        x=r*c/r_200

        xsq=x*x

        xsqm1=xsq-1.d0

        rxsqm1=xsqm1**0.5

        pmxsqm1=1.d0-xsq

        rpmxsqm1=pmxsqm1**0.5

       fac1=2.d0*(r_200/c)*deltac*rhoc/sigcrit

       faci=((1.d0-x)/(1.d0+x))**0.5

c      fac2=atanh(((1.d0-x)/(1.d0+x))**0.5)

      fac2=0.5d0*log((-1.d0-faci)/(faci-1.d0))  

      fac3=atan(((x-1.d0)/(1.d0+x))**0.5)
c	These conditional statements find the values for SigmaNFW
      if (x.lt.1.) then
         appa = (fac1/xsqm1)
     +        *(1.d0-(2.d0*fac2/rpmxsqm1))
      endif
      if (x.eq.1.) then
         appa=fac1/3.d0
      endif
      if (x.gt.1.) then
         appa = (fac1/xsqm1)
     +        *(1.d0-(2.d0*fac3/rxsqm1))
      endif

c	These find the Mean of Sigma NFW
      if (x.lt.1.) then
         appam = (2.d0*fac1/xsq)
     +        *((2.d0*fac2/rpmxsqm1)+log(x/2.))
      endif

      if (x.eq.1.) then
         appam =2.d0*fac1*(1.d0+log(0.5d0))
      endif

      if (x.gt.1.) then
         appam = (2.d0*fac1/xsq)
     +        *((2.d0*fac3/(rxsqm1))+log(x/2.))
      endif
c	gamma = (MeanSigmaNFW - SigmaNFW ) / SigmaCrit
	gamma=appam-appa

	emu=abs(1.d0/((1.d0-appa)**2-gamma**2))

c        write(98,*)r,emu,appa,gamma

 	g=gamma/(1.-appa)

	g1=-g*cos(2.*ph)
	g2=-g*sin(2.*ph)

	end


       function func(p)
       implicit real*8(a-h,o-z)
       parameter (PI=3.1415926536)
       parameter (ho=70.d0)
       parameter (GG=6.67e-11)
       parameter (cv=3.e5,n=2048)
       parameter (nmax=100000)
       dimension p(2)
       complex*16 epl,ep,epld,g
       common/pass1/epld(nmax),phild(nmax),thld(nmax),nk
       common/pass2/sigcrit,rhoc,dol,rin,rout
c       common/pass3/sig

       sig=0.2d0
    

       r_200=p(1)
       c=p(2)

       uss = 1./(sig**2)
 
       summep=0.d0
       summ = 0.d0
       suma = 0.d0
       sumt=  0.d0
       sum=   0.d0

       deltac=(200.d0/3.d0)*(c**3)/(log(1.+c)-(c/(1.d0+c)))

c       write(*,*)p(1),p(2)

        do ii=1,nk

             th=thld(ii)
             epl=epld(ii)
             phi=phild(ii)
  
           r=th
             x=r*c/r_200
             xsq=x*x
      xsqm1=xsq-1.
      rxsqm1=xsqm1**0.5
      pmxsqm1=1.-xsq
      rpmxsqm1=pmxsqm1**0.5
      fac1=2.d0*(r_200/c)*deltac*rhoc/sigcrit
      faci=((1.d0-x)/(1.d0+x))**0.5
      fac2=0.5d0*log((-1.-faci)/(faci-1.))  
      fac3=atan(((x-1.)/(1.+x))**0.5)
      if (x.lt.1.) then
         appa = (fac1/xsqm1)
     +        *(1.d0-(2.d0*fac2/rpmxsqm1))
      endif
      if (x.eq.1.) then
         appa=fac1/3.d0
      endif
      if (x.gt.1.) then
         appa = (fac1/xsqm1)
     +        *(1.d0-(2.d0*fac3/rxsqm1))
      endif


      if (x.lt.1.) then
         appam = (2.d0*fac1/xsq)
     +        *((2.d0*fac2/rpmxsqm1)+log(x/2.))
      endif
      if (x.eq.1.) then
         appam=2.d0*fac1*(1.d0+log(0.5))
      endif
      if (x.gt.1.) then
         appam = (2.d0*fac1/xsq)
     +        *((2.d0*fac3/(rxsqm1))+log(x/2.))
      endif


      glx=8.*fac2/(xsq*rpmxsqm1)
     +     +4.d0*log(x/2.)/xsq
     +     -2.d0/xsqm1
     +     +4.d0*fac2/(xsqm1*rpmxsqm1)

      ggx=8.*fac3/(xsq*rxsqm1)
     +     +4.d0*log(x/2.)/xsq
     +     -2.d0/xsqm1
     +     +4.d0*fac3/(xsqm1**1.5)


      if (x.lt.1.) then
         gamma=(fac1/2.d0)*glx
      endif
      if (x.eq.1.) then
         gamma=(fac1/2.d0)*(4.d0*log(0.5)+10.d0/3.d0)
      endif
      if (x.gt.1.) then
         gamma=(fac1/2.d0)*ggx
      endif

      deta=((1.d0-appa)**2)-(gamma**2)

      emu=abs(1.d0/deta)


c      if (((appa+gamma)-1.)**2.le.0.1) then 
c         write(85,*)c,r_200,emu,th
c        endif
c        if (((appa-gamma)-1.)**2.le.0.1) then
c        write(86,*)c,r_200,emu,th

c        endif



c        if (ii.le.int(nk/4)) then

        gr=gamma/(1.d0-appa)

        g1=-gr*cos(2.d0*phi)

        g2=-gr*sin(2.d0*phi)

         g=complex(g1,g2)        

         uss=(1.d0/sig)**2

         if (abs(g).le.1.0) then
         Ya = (((abs(g)**2) - 1.)**2)/(abs((epl*conjg(g))-1.)**4)
         ep=(epl-g)/(1.-(conjg(g)*epl))
         elseif (abs(g).gt.1.0) then
         Ya = ((abs(g)**2 - 1.)**2)/(abs(epl-g)**4)
         ep=(1.-(g*conjg(epl)))/(conjg(epl)-conjg(g))
         endif

          pref= ya/(pi*(sig**2.)*(1.-exp(-uss)))
          Pepl = pref*exp(-(abs(ep)/sig)**2)
          summep = summep + log(Pepl)


       enddo

        func=-summep

       return
       end



      function dangg(zz)
      implicit real*8(a-h,o-z)

    	om=0.3d0
    	ol=0.7d0
        dangg=1.d0/((om*((1.d0+zz)**3)) + ol)**0.5

      return
      end

        


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END


      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL*8 data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END


      FUNCTION poidev(xm,idum)
      INTEGER idum
      REAL*8 poidev,xm,PI
      PARAMETER (PI=3.141592654)
CU    USES gammln,ran1
      REAL*8 alxm,em,g,oldm,sq,t,y,gammln,ran1
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
 2             em=em+1.
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
 1             y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END


      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
 11         continue
      gammln=tmp+log(stp*ser/x)
      return
      END

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
      INTEGER iter,n,np,NMAX,ITMAX
      REAL*8 fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=20,ITMAX=200)
CU    USES func,linmin
      INTEGER i,ibig,j
      REAL*8 del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END


      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      REAL*8 fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.e-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL*8 brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      write(*,*) 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END


      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL*8 f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      REAL*8 pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL*8 dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
