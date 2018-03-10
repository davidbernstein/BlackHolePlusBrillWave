c Copyright 1992-2018, David Bernstein
c
c This file is part of BlackHolePlusBrillWave.
c   
c BlackHolePlusBrillWave is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.
c    
c BlackHolePlusBrillWave is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c    
c You should have received a copy of the GNU General Public License
c along with BlackHolePlusBrillWave.  If not, see <http://www.gnu.org/licenses/>.

      program onebh10
      implicit double precision (a-h,o-z)
c
c        An axisymmetric black hole evolution code.
c        This code uses the multigrid solver usmg2.
c        Created 9 March, 1992.
c
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c     waytogo='batch'
      waytogo='interact'
c     waytogo='ivpcrank'
c
      if (waytogo.eq.'batch') then
      call batchrun
      elseif (waytogo.eq.'interact') then
      call interact
      elseif (waytogo.eq.'ivpcrank') then
      call ivpcrank
      endif
c
c
      stop
      end
c
c
c
      subroutine ivpcrank
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      character*8 qrun,waytogo,dervans,timesche,timescal
      common /charblck/ qrun,waytogo,dervans,timesche,timescal
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
      real amp(40),range(40),width(40)
c
      da1=0.
      dr1=0.
      dw1=1.
      da2=0.
      dr2=0.
      dw2=1.
      call paramtrs
c
c        open file for output from initial value problem solver
      open(67,file='ivp.dat',status='unknown')
c
c        do experiment
      amplow=-1.
      amphigh=3.
      damp=0.25
      kamp=nint((amphigh-amplow)/damp)+1
      do 5 k=1,kamp
      amp(k)=amplow+(k-1.)*damp
    5 continue
c
      ranglow=0.
      ranghigh=2.5
      drang=0.25
      krang=nint((ranghigh-ranglow)/drang)+1
      do 6 k=1,krang
      range(k)=ranglow+(k-1.)*drang
    6 continue
c
c        write header for ivp.dat
c     write(67,*) kamp
c     write(67,*) krang
c     do 7 k=1,kamp
c   7 write(67,*) amp(k)
c     do 8 k=1,krang
c   8 write(67,*) range(k)
c
      dw1=1.
c
      do 10 ki=1,kamp
      da1=amp(ki)
      write(67,*)
c
      do 10 kj=1,krang
      dr1=range(kj)
c
      call ivp
c
      call admmass(ans1) 
      ans1=ans1/scaleprm
c
c     call ivptest1(ans2)
c
      write(67,*) da1,dr1,ans1
c
c     print*,da1,dr1,dw1,ans1,ans2
   10 continue
c
  100 format(5e15.6)
c
      close(67)
c
c
      return
      end
c
c
c
      subroutine batchrun
      include 'param.h'
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        set distortion coefficients
      da1=-1.
      dr1=2.
      dw1=1.
c
      da2=0.0
      dr2=0.
      dw2=1.
c
c        set spatial differencing option (s,f,s24)
      dervans='f'
c        set temporal differencing option (b,m)
      timesche='b'
c
      qrun='y'
c
c        total time of the run is in units of the ADM mass
      time=70.
c
      call go
c
c
      return
      end
c
c
c
      subroutine interact
      include 'param.h'
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        set distortion constants da and dr for IVP
      write(5,*)
      write(5,*)'distortion parameters (magnitude,range,width)?'
      write(5,*)'first set:'
      read(5,*) da1,dr1,dw1
      write(5,*)
      write(5,*)'second set?'
      read(5,100) qrun
      if (qrun.eq.'y') then
      read(5,*) da2,dr2,dw2
      else
      da2=0.
      dr2=0.
      dw2=1.
      endif
c
      write(5,*)
      write(5,*)'spatial differencing scheme (s,f,s24)?'
      read(5,100) dervans
c
      write(5,*)
      write(5,*)'temporal differencing scheme (b,m)?'
      read(5,100) timesche
c
      write(5,*)
      write(5,*)'compute complete run (y) or only initial data (n)?'
      read(5,100) qrun
  100 format(8a)
c
      call go
c
c
      return
      end
c
c
c
      subroutine go
      include 'param.h'
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
      common /equblock/ nequ,nfequ,fiequ
      common /horiznbk/ nah,nfah,fiah
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /rpatblck/ nradpat,nfirp,firp
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        set space and time steps, open data files, etc. 
      call paramtrs
c
c        initialize variables and compute some constants
      call initial
c
c        output initial data
      call allball2
c
      call brillmas(ambrill)
      ratio=ambrill/adm-1.
      write(5,*)'(Brill mass-ADM)/ADM is   ',ratio
c
      call hwkgmass(amhwk)
      ratio=amhwk/adm-1.
      write(5,*)'(Hawking mass-ADM)/ADM is  ',ratio
c
      call ivptest1(ans)
      write(5,*)
      write(5,*)'hmass/ADM is ',ans/adm
      write(5,*)dlog10(ans)
c
      call ivptest2
c
c        find initial apparent horizon
      call firstah
      call circumf(ans)
      write(5,*)
      write(5,*)'polar to equatorial circumf of apparent horizon ',ans
      call ahmass(ans)
      write(5,*)'mass of apparent horizon (/Madm) ',ans/adm
      write(5,*)'max radiation efficiency ',1.-ans/adm
      write(5,*)
      call ahtest1(ans)
      write(5,*)'max ah residual test1 ',ans,dlog10(ans)
      write(5,*)
c
c     call wojtkiew
c
c        event horizon estimate
      call eventest(ans)
      write(5,*)'location of mass M surface at eta= ',ans
      write(5,*)
c
      if (nwaves.gt.1) call zerilli
      if (nequ.gt.1)  call eqout
      if (nframes.gt.1) call output
      if (nah.gt.1) call ahout
      if (nmwaves.gt.1) call zermovie
      if (nradpat.gt.1) call rpout
c
c
      if (qrun.eq.'n') go to 99
c
c        Go dogs go!
      do 10 jtime=0,ntime-1
c
      if (timesche.eq.'m') then
      call mcormack
      elseif (timesche.eq.'b') then
      call brlvskya
      endif
c
      if (jtime.eq.nfequ.and.nequ.gt.1) call eqout
c
      if (jtime.eq.nf.and.nframes.gt.1) call output
c
      if (jtime.eq.nfah.and.nah.gt.1) call ahout
c
      if (jtime.eq.nfwave.and.nwaves.gt.1) call zerilli
c
      if (jtime.eq.nfmwave.and.nmwaves.gt.1) call zermovie
c
      if(jtime.eq.nfirp.and.nradpat.gt.1) call rpout
c
c     if (mod(jtime,nint(ntime/real(nmovie))).eq.0) call movieout
c
   10 continue
c
c        output final data
      if (nwaves.gt.1) call zerilli
      if (nequ.gt.1) call eqout
      if (nframes.gt.1) call output
      if (nah.gt.1) call ahout
      if (nmwaves.gt.1) call zermovie
      if (nradpat.gt.1) call rpout
c
      call ahtest1(ans)
      write(5,*)'max ah residual test1 ',ans,dlog10(ans)
      write(5,*)
c
      call alldone
c
   99 continue
c
c
      return
      end
c
c
c
      subroutine brlvskya
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /fdotblck/ adot(id,jd),bdot(id,jd),
     &                  cdot(id,jd),ddot(id,jd)
      common /hdotblck/ hadot(id,jd),hbdot(id,jd),
     &                  hcdot(id,jd),hddot(id,jd)
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real  ta(id,jd), tb(id,jd), tc(id,jd), td(id,jd),
     &     tha(id,jd),thb(id,jd),thc(id,jd),thd(id,jd)
c
c
c        store variables temporarily
      do 10 j=1,nthtamax
      do 10 i=1,netamax
c
      ta(i,j)=a(i,j)
      tb(i,j)=b(i,j)
      tc(i,j)=c(i,j)
      td(i,j)=d(i,j)
c
      tha(i,j)=ha(i,j)
      thb(i,j)=hb(i,j)
      thc(i,j)=hc(i,j)
      thd(i,j)=hd(i,j)
   10 continue
c
      call allball1
c
      do 20 j=1,nthtamax
      do 20 i=1,netamax
c     
      a(i,j)=a(i,j)+adot(i,j)*dt
      b(i,j)=b(i,j)+bdot(i,j)*dt
      c(i,j)=c(i,j)+cdot(i,j)*dt
      d(i,j)=d(i,j)+ddot(i,j)*dt
c
      ha(i,j)=ha(i,j)+hadot(i,j)*dt
      hb(i,j)=hb(i,j)+hbdot(i,j)*dt
      hc(i,j)=hc(i,j)+hcdot(i,j)*dt
      hd(i,j)=hd(i,j)+hddot(i,j)*dt
   20 continue
c
      call allball2
c
      do 30 j=1,nthtamax
      do 30 i=1,netamax
c     
      a(i,j)=ta(i,j)+adot(i,j)*dt
      b(i,j)=tb(i,j)+bdot(i,j)*dt
      c(i,j)=tc(i,j)+cdot(i,j)*dt
      d(i,j)=td(i,j)+ddot(i,j)*dt
c
      ha(i,j)=tha(i,j)+hadot(i,j)*dt
      hb(i,j)=thb(i,j)+hbdot(i,j)*dt
      hc(i,j)=thc(i,j)+hcdot(i,j)*dt
      hd(i,j)=thd(i,j)+hddot(i,j)*dt
   30 continue
c
c        save wave zone boundary values for next time step
      do 40 j=1,nthtamax
      abnd(j)=a(netamax,j)
      bbnd(j)=b(netamax,j)
      cbnd(j)=c(netamax,j)
      dbnd(j)=d(netamax,j)
c
      habnd(j)=ha(netamax,j)
      hbbnd(j)=hb(netamax,j)
      hcbnd(j)=hc(netamax,j)
      hdbnd(j)=hd(netamax,j)
   40 continue
c
c
      return
      end
c
c
c
      subroutine mcormack
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /fdotblck/ adot(id,jd),bdot(id,jd),
     &                  cdot(id,jd),ddot(id,jd)
      common /hdotblck/ hadot(id,jd),hbdot(id,jd),
     &                  hcdot(id,jd),hddot(id,jd)
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real  ta(id,jd), tb(id,jd), tc(id,jd), td(id,jd),
     &     tha(id,jd),thb(id,jd),thc(id,jd),thd(id,jd)
c
c
c        store variables temporarily
      do 10 j=1,nthtamax
      do 10 i=1,netamax
c
      ta(i,j)=a(i,j)
      tb(i,j)=b(i,j)
      tc(i,j)=c(i,j)
      td(i,j)=d(i,j)
c
      tha(i,j)=ha(i,j)
      thb(i,j)=hb(i,j)
      thc(i,j)=hc(i,j)
      thd(i,j)=hd(i,j)
   10 continue
c
      call allball1
c
      do 20 j=1,nthtamax
      do 20 i=1,netamax
c     
      a(i,j)=a(i,j)+adot(i,j)*dt
      b(i,j)=b(i,j)+bdot(i,j)*dt
      c(i,j)=c(i,j)+cdot(i,j)*dt
      d(i,j)=d(i,j)+ddot(i,j)*dt
c
      ha(i,j)=ha(i,j)+hadot(i,j)*dt
      hb(i,j)=hb(i,j)+hbdot(i,j)*dt
      hc(i,j)=hc(i,j)+hcdot(i,j)*dt
      hd(i,j)=hd(i,j)+hddot(i,j)*dt
   20 continue
c
      call allball2
c
      do 30 j=1,nthtamax
      do 30 i=1,netamax
c     
      a(i,j)=0.5*(a(i,j)+ta(i,j)+adot(i,j)*dt)
      b(i,j)=0.5*(b(i,j)+tb(i,j)+bdot(i,j)*dt)
      c(i,j)=0.5*(c(i,j)+tc(i,j)+cdot(i,j)*dt)
      d(i,j)=0.5*(d(i,j)+td(i,j)+ddot(i,j)*dt)
c
      ha(i,j)=0.5*(ha(i,j)+tha(i,j)+hadot(i,j)*dt)
      hb(i,j)=0.5*(hb(i,j)+thb(i,j)+hbdot(i,j)*dt)
      hc(i,j)=0.5*(hc(i,j)+thc(i,j)+hcdot(i,j)*dt)
      hd(i,j)=0.5*(hd(i,j)+thd(i,j)+hddot(i,j)*dt)
   30 continue
c
c        save wave zone boundary values for next time step
      do 40 j=1,nthtamax
      abnd(j)=a(netamax,j)
      bbnd(j)=b(netamax,j)
      cbnd(j)=c(netamax,j)
      dbnd(j)=d(netamax,j)
c
      habnd(j)=ha(netamax,j)
      hbbnd(j)=hb(netamax,j)
      hcbnd(j)=hc(netamax,j)
      hdbnd(j)=hd(netamax,j)
   40 continue
c
c
      return
      end
c
c
c
      subroutine allball1
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        construct the inverse of e
      do 10 j=1,nthtamax
      do 10 i=1,netamax
   10 em1(i,j)=1./(a(i,j)*b(i,j)-c(i,j)*c(i,j))
c
      if (dervans.eq.'s') then
      call diff2
      elseif (dervans.eq.'f') then
      call diff4
      elseif (dervans.eq.'s24') then
      call diff24
      endif
c
c     call ricci
      call riccinoc
c
      call maxlapse
c     call alglapse
c
      call shiftc0
c
      call dot
c     call macsmdot
c
      return
      end
c
c
c
      subroutine allball2
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        construct the inverse of e
      do 10 j=1,nthtamax
      do 10 i=1,netamax
   10 em1(i,j)=1./(a(i,j)*b(i,j)-c(i,j)*c(i,j))
c
      if (dervans.eq.'s') then
      call diff2
      elseif (dervans.eq.'f') then
      call diff4
      elseif (dervans.eq.'s24') then
      call diff24
      endif
c
c     call ricci
      call riccinoc
c
      call dot
c     call macsmdot
c
      return
      end
c
c
c
      subroutine diff2
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c        compute metric derivatives
      call scnddrvs(a,a1n,a1t,a2n,a2t,a2m, 1.,abnd,
     &              netamax,nthtamax)
      call scnddrvs(b,b1n,b1t,b2n,b2t,b2m, 1.,bbnd,
     &              netamax,nthtamax)
      call scnddrvs(c,c1n,c1t,c2n,c2t,c2m,-1.,cbnd,
     &              netamax,nthtamax)
      call scnddrvs(d,d1n,d1t,d2n,d2t,d2m, 1.,dbnd,
     &              netamax,nthtamax)
c
c        and extrinsic curvature first derivatives
      call frstdrv2(ha,ha1n,ha1t, 1.,habnd,netamax,nthtamax)
      call frstdrv2(hb,hb1n,hb1t, 1.,hbbnd,netamax,nthtamax)
      call frstdrv2(hc,hc1n,hc1t,-1.,hcbnd,netamax,nthtamax)
      call frstdrv2(hd,hd1n,hd1t, 1.,hdbnd,netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine diff4
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c        compute metric derivatives
      call frthdrvs(a,a1n,a1t,a2n,a2t,a2m, 1.,abnd,
     &              netamax,nthtamax)
      call frthdrvs(b,b1n,b1t,b2n,b2t,b2m, 1.,bbnd,
     &              netamax,nthtamax)
      call frthdrvs(c,c1n,c1t,c2n,c2t,c2m,-1.,cbnd,
     &              netamax,nthtamax)
      call frthdrvs(d,d1n,d1t,d2n,d2t,d2m, 1.,dbnd,
     &              netamax,nthtamax)
c
c        and extrinsic curvature first derivatives
      call frstdrv4(ha,ha1n,ha1t, 1.,habnd,netamax,nthtamax)
      call frstdrv4(hb,hb1n,hb1t, 1.,hbbnd,netamax,nthtamax)
      call frstdrv4(hc,hc1n,hc1t,-1.,hcbnd,netamax,nthtamax)
      call frstdrv4(hd,hd1n,hd1t, 1.,hdbnd,netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine diff24
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c        compute metric derivatives
      call scheme24(a,a1n,a1t,a2n,a2t,a2m, 1.,abnd,
     &              netamax,nthtamax)
      call scheme24(b,b1n,b1t,b2n,b2t,b2m, 1.,bbnd,
     &              netamax,nthtamax)
      call scheme24(c,c1n,c1t,c2n,c2t,c2m,-1.,cbnd,
     &              netamax,nthtamax)
      call scheme24(d,d1n,d1t,d2n,d2t,d2m, 1.,dbnd,
     &              netamax,nthtamax)
c
c        and extrinsic curvature first derivatives
      call frstdv24(ha,ha1n,ha1t, 1.,habnd,netamax,nthtamax)
      call frstdv24(hb,hb1n,hb1t, 1.,hbbnd,netamax,nthtamax)
      call frstdv24(hc,hc1n,hc1t,-1.,hcbnd,netamax,nthtamax)
      call frstdv24(hd,hd1n,hd1t, 1.,hdbnd,netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine riccinoc
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /chrisblk/ g111(id,jd),g112(id,jd),g122(id,jd),
     &                  g133(id,jd),g211(id,jd),g212(id,jd),
     &                  g222(id,jd),g233(id,jd),g313(id,jd),
     &                  g323(id,jd)
      common /ricciblk/ ricci11(id,jd),ricci12(id,jd),
     &                  ricci22(id,jd),ricci33(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      real r1(14),r2(14),r3(16),r4(17)
      real r1212(id,jd),r1313(id,jd),r2323(id,jd),r1323(id,jd)
c
c
c        construct Christoffel symbols and Ricci tensor
c        components off axis
c        (stored in ricci33 is actually ricci33/sin(theta)**2)
c        (and g133 and g233 are similarly divided by sin(theta)**2)
c
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      g111(i,j)=0.5*(a1n(i,j)/a(i,j)+phi1n(i,j))
      g112(i,j)=0.5*(a1t(i,j)/a(i,j)+phi1t(i,j))
      g122(i,j)=-0.5*(b1n(i,j)+b(i,j)*phi1n(i,j))/a(i,j)
      g211(i,j)=-0.5*(a1t(i,j)+a(i,j)*phi1t(i,j))/b(i,j)
      g212(i,j)=0.5*(b1n(i,j)/b(i,j)+phi1n(i,j))
      g222(i,j)=0.5*(b1t(i,j)/b(i,j)+phi1t(i,j))
c
      g313(i,j)=0.5*(d1n(i,j)/d(i,j)+phi1n(i,j))
      g323(i,j)=0.5*(d1t(i,j)/d(i,j)+phi1t(i,j))+ct(j)
c
      g133(i,j)=-0.5*(d1n(i,j)+d(i,j)*phi1n(i,j))/a(i,j)
      g233(i,j)=-0.5*(d1t(i,j)+d(i,j)*phi1t(i,j))/b(i,j)
     &          -d(i,j)*ct(j)/b(i,j)
   10 continue
c
      do 11 j=2,nthtamax
      do 11 i=1,netamax
      r1(1)=-0.5*a(i,j)*phi2t(i,j)
      r1(2)=0.25*a(i,j)*b1t(i,j)*phi1t(i,j)/b(i,j)
      r1(3)=-0.25*a1t(i,j)*phi1t(i,j)
      r1(4)=-0.5*b(i,j)*phi2n(i,j)
      r1(5)=0.25*a1n(i,j)*b(i,j)*phi1n(i,j)/a(i,j)
      r1(6)=-0.25*b1n(i,j)*phi1n(i,j)
      r1(7)=0.25*a1t(i,j)*b1t(i,j)/b(i,j)
      r1(8)=0.25*b1n(i,j)*b1n(i,j)/b(i,j)
      r1(9)=0.25*a1n(i,j)*b1n(i,j)/a(i,j)
      r1(10)=0.25*a1t(i,j)*a1t(i,j)/a(i,j)
      r1(11)=-0.5*b2n(i,j)
      r1(12)=-0.5*a2t(i,j)
c
      r1212(i,j)=r1( 1)+r1( 2)+r1( 3)+r1( 4)+r1( 5)+r1( 6)+r1( 7)+r1( 8)
     &          +r1( 9)+r1(10)+r1(11)+r1(12)
   11 continue
c
      do 12 j=2,nthtamax
      do 12 i=1,netamax
      r2(1)=-0.5*a(i,j)*d(i,j)*phi1t(i,j)*ct(j)/b(i,j)
      r2(2)=-0.5*a1t(i,j)*d(i,j)*ct(j)/b(i,j)
      r2(3)=-0.25*a(i,j)*d(i,j)*phi1t(i,j)*phi1t(i,j)/b(i,j)
      r2(4)=-0.25*a(i,j)*d1t(i,j)*phi1t(i,j)/b(i,j)
      r2(5)=-0.25*a1t(i,j)*d(i,j)*phi1t(i,j)/b(i,j)
      r2(6)=-0.5*d(i,j)*phi2n(i,j)
      r2(7)=-0.25*d1n(i,j)*phi1n(i,j)
      r2(8)=0.25*a1n(i,j)*d(i,j)*phi1n(i,j)/a(i,j)
      r2(9)=-0.25*a1t(i,j)*d1t(i,j)/b(i,j)
      r2(10)=0.25*a1n(i,j)*d1n(i,j)/a(i,j)
      r2(11)=-0.5*d2n(i,j)
      r2(12)=0.25*d1n(i,j)*d1n(i,j)/d(i,j)
c
      r1313(i,j)=r2( 1)+r2( 2)+r2( 3)+r2( 4)+r2( 5)+r2( 6)+r2( 7)+r2( 8)
     &     +r2( 9)+r2(10)+r2(11)+r2(12)
   12 continue
c
      do 13 j=2,nthtamax
      do 13 i=1,netamax
      r3(1)=-0.5*d(i,j)*phi1t(i,j)*ct(j)
      r3(2)=0.5*b1t(i,j)*d(i,j)*ct(j)/b(i,j)
      r3(3)=-d1t(i,j)*ct(j)
      r3(4)=-0.5*d(i,j)*phi2t(i,j)
      r3(5)=-0.25*d1t(i,j)*phi1t(i,j)
      r3(6)=0.25*b1t(i,j)*d(i,j)*phi1t(i,j)/b(i,j)
      r3(7)=-0.25*b(i,j)*d(i,j)*phi1n(i,j)*phi1n(i,j)/a(i,j)
      r3(8)=-0.25*b(i,j)*d1n(i,j)*phi1n(i,j)/a(i,j)
      r3(9)=-0.25*b1n(i,j)*d(i,j)*phi1n(i,j)/a(i,j)
      r3(10)=0.25*b1t(i,j)*d1t(i,j)/b(i,j)
      r3(11)=-0.25*b1n(i,j)*d1n(i,j)/a(i,j)
      r3(12)=-0.5*d2t(i,j)
      r3(13)=0.25*d1t(i,j)*d1t(i,j)/d(i,j)
      r3(14)=d(i,j)
c
      r2323(i,j)=r3( 1)+r3( 2)+r3( 3)+r3( 4)+r3( 5)+r3( 6)+r3( 7)+r3( 8)
     &     +r3( 9)+r3(10)+r3(11)+r3(12)+r3(13)+r3(14)
   13 continue
c
      do 14 j=2,nthtamax
      do 14 i=1,netamax
      r4(1)=0.5*b1n(i,j)*d(i,j)*ct(j)/b(i,j)
      r4(2)=-0.5*d1n(i,j)*ct(j)
      r4(3)=0.25*d(i,j)*phi1n(i,j)*phi1t(i,j)
      r4(4)=0.25*b1n(i,j)*d(i,j)*phi1t(i,j)/b(i,j)
      r4(5)=-0.5*d(i,j)*phi2m(i,j)
      r4(6)=0.25*a1t(i,j)*d(i,j)*phi1n(i,j)/a(i,j)
      r4(7)=0.25*b1n(i,j)*d1t(i,j)/b(i,j)
      r4(8)=0.25*a1t(i,j)*d1n(i,j)/a(i,j)
      r4(9)=0.25*d1n(i,j)*d1t(i,j)/d(i,j)
      r4(10)=-0.5*d2m(i,j)
c
      r1323(i,j)=r4( 1)+r4( 2)+r4( 3)+r4( 4)+r4( 5)+r4( 6)+r4( 7)+r4( 8)
     &          +r4( 9)+r4(10)
   14 continue
c
      do 15 j=2,nthtamax
      do 15 i=1,netamax
      ricci11(i,j)=r1212(i,j)/b(i,j)+r1313(i,j)/d(i,j)
      ricci12(i,j)=r1323(i,j)/d(i,j)
      ricci22(i,j)=r1212(i,j)/a(i,j)+r2323(i,j)/d(i,j)
      ricci33(i,j)=r1313(i,j)/a(i,j)+r2323(i,j)/b(i,j)
   15 continue
c
c        and on axis
      do 20 i=1,netamax
      r1(1)=-0.5*a(i,1)*phi2t(i,1)
      r1(2)=-0.5*b(i,1)*phi2n(i,1)
      r1(3)=0.25*a1n(i,1)*b(i,1)*phi1n(i,1)/a(i,1)
      r1(4)=-0.25*b1n(i,1)*phi1n(i,1)
      r1(5)=0.25*b1n(i,1)*b1n(i,1)/b(i,1)
      r1(6)=0.25*a1n(i,1)*b1n(i,1)/a(i,1)
      r1(7)=-0.5*b2n(i,1)
      r1(8)=-0.5*a2t(i,1)
c
      r1212(i,1)=r1(1)+r1(2)+r1(3)+r1(4)+r1(5)+r1(6)+r1(7)+r1(8)
   20 continue
c
      do 21 i=1,netamax
      r2(1)=-0.5*a(i,1)*phi2t(i,1)
      r2(2)=-0.5*a2t(i,1)
      r2(3)=-0.5*d(i,1)*phi2n(i,1)
      r2(4)=-0.25*d1n(i,1)*phi1n(i,1)
      r2(5)=0.25*a1n(i,1)*d(i,1)*phi1n(i,1)/a(i,1)
      r2(6)=0.25*a1n(i,1)*d1n(i,1)/a(i,1)
      r2(7)=-0.5*d2n(i,1)
      r2(8)=0.25*d1n(i,1)*d1n(i,1)/d(i,1)
c
      r1313(i,1)=r2(1)+r2(2)+r2(3)+r2(4)+r2(5)+r2(6)+r2(7)+r2(8)
   21 continue
c
      do 22 i=1,netamax
      r3(1)=-0.5*d(i,1)*phi2t(i,1)
      r3(2)=0.5*b2t(i,1)
      r3(3)=-1.5*d2t(i,1)
      r3(4)=-0.5*d(i,1)*phi2t(i,1)
      r3(5)=-0.25*b(i,1)*d(i,1)*phi1n(i,1)*phi1n(i,1)/a(i,1)
      r3(6)=-0.25*b(i,1)*d1n(i,1)*phi1n(i,1)/a(i,1)
      r3(7)=-0.25*b1n(i,1)*d(i,1)*phi1n(i,1)/a(i,1)
      r3(8)=-0.25*b1n(i,1)*d1n(i,1)/a(i,1)
      r3(9)=d(i,1)
c
      r2323(i,1)=r3(1)+r3(2)+r3(3)+r3(4)+r3(5)+r3(6)+r3(7)+r3(8)+r3(9)
   22 continue
c
      do 23 i=1,netamax
      ricci11(i,1)=r1212(i,1)/b(i,1)+r1313(i,1)/d(i,1)
      ricci22(i,1)=r1212(i,1)/a(i,1)+r2323(i,1)/d(i,1)
      ricci33(i,1)=r1313(i,1)/a(i,1)+r2323(i,1)/b(i,1)
      ricci12(i,1)=0.
   23 continue
c
c
      return
      end
c
c
c
      subroutine dot
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /eblock  / em1(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /chrisblk/ g111(id,jd),g112(id,jd),g122(id,jd),
     &                  g133(id,jd),g211(id,jd),g212(id,jd),
     &                  g222(id,jd),g233(id,jd),g313(id,jd),
     &                  g323(id,jd)
      common /ricciblk/ ricci11(id,jd),ricci12(id,jd),
     &                  ricci22(id,jd),ricci33(id,jd)
      common /hdotblck/ hadot(id,jd),hbdot(id,jd),
     &                  hcdot(id,jd),hddot(id,jd)
      common /fdotblck/ adot(id,jd),bdot(id,jd),
     &                  cdot(id,jd),ddot(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /shiftblk/ beta1(id,jd),beta2(id,jd),omega(id,jd)
      common /shftdblk/ beta1n(id,jd),beta1t(id,jd),
     &                  beta2n(id,jd),beta2t(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
c
c
c        construct time derivatives off axis
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      trk=em1(i,j)*(b(i,j)*ha(i,j)+a(i,j)*hb(i,j)-2.*c(i,j)*hc(i,j))
     &   +hd(i,j)/d(i,j)
c
      r1=alpha(i,j)*ricci11(i,j)-alpha2n(i,j)+
     &   g111(i,j)*alpha1n(i,j)+g211(i,j)*alpha1t(i,j)
      r2=2.*em1(i,j)*(ha(i,j)*(b(i,j)*ha(i,j)-c(i,j)*hc(i,j))+
     &                hc(i,j)*(a(i,j)*hc(i,j)-c(i,j)*ha(i,j)))
      hageomdt=r1*psim4(i,j)+alpha(i,j)*(ha(i,j)*trk-r2)         
c
      r3=alpha(i,j)*ricci22(i,j)-alpha2t(i,j)+
     &   g122(i,j)*alpha1n(i,j)+g222(i,j)*alpha1t(i,j)
      r4=2.*em1(i,j)*(hc(i,j)*(b(i,j)*hc(i,j)-c(i,j)*hb(i,j))+
     &                hb(i,j)*(a(i,j)*hb(i,j)-c(i,j)*hc(i,j)))
      hbgeomdt=r3*psim4(i,j)+alpha(i,j)*(hb(i,j)*trk-r4)         
c
      r5=alpha(i,j)*ricci12(i,j)-alpha2m(i,j)+
     &   g112(i,j)*alpha1n(i,j)+g212(i,j)*alpha1t(i,j)
      r6=2.*em1(i,j)*(hc(i,j)*(b(i,j)*ha(i,j)-c(i,j)*hc(i,j))+
     &                hb(i,j)*(a(i,j)*hc(i,j)-c(i,j)*ha(i,j)))
      hcgeomdt=r5*psim4(i,j)+alpha(i,j)*(hc(i,j)*trk-r6)         
c
      r7=alpha(i,j)*ricci33(i,j)+
     &   g133(i,j)*alpha1n(i,j)+g233(i,j)*alpha1t(i,j)
      r8=2.*hd(i,j)*hd(i,j)/d(i,j)
      hdgeomdt=r7*psim4(i,j)+alpha(i,j)*(hd(i,j)*trk-r8)
c
      t1=beta1(i,j)*(ha1n(i,j)+ha(i,j)*phi1n(i,j))
     &  +beta2(i,j)*(ha1t(i,j)+ha(i,j)*phi1t(i,j))
      t2=2.*(ha(i,j)*beta1n(i,j)+hc(i,j)*beta2n(i,j))
      hashftdt=t1+t2
c
      t3=beta1(i,j)*(hb1n(i,j)+hb(i,j)*phi1n(i,j))
     &  +beta2(i,j)*(hb1t(i,j)+hb(i,j)*phi1t(i,j))
      t4=2.*(hc(i,j)*beta1t(i,j)+hb(i,j)*beta2t(i,j))
      hbshftdt=t3+t4
c
      t5=beta1(i,j)*(hc1n(i,j)+hc(i,j)*phi1n(i,j))
     &  +beta2(i,j)*(hc1t(i,j)+hc(i,j)*phi1t(i,j))
      t6=hc(i,j)*(beta1n(i,j)+beta2t(i,j))+ha(i,j)*beta1t(i,j)
     &                                    +hb(i,j)*beta2n(i,j)
      hcshftdt=t5+t6
c
      hdshftdt=beta1(i,j)*(hd1n(i,j)+hd(i,j)*phi1n(i,j))
     &        +beta2(i,j)*(hd1t(i,j)+hd(i,j)*(2.*ct(j)+phi1t(i,j)))
c
      hadot(i,j)=hageomdt+hashftdt
      hbdot(i,j)=hbgeomdt+hbshftdt
      hcdot(i,j)=hcgeomdt+hcshftdt
      hddot(i,j)=hdgeomdt+hdshftdt
c
      ageomdt=-2.*alpha(i,j)*ha(i,j)
      bgeomdt=-2.*alpha(i,j)*hb(i,j)
      cgeomdt=-2.*alpha(i,j)*hc(i,j)
      dgeomdt=-2.*alpha(i,j)*hd(i,j)
c
      ashftdt=2.*(a(i,j)*beta1n(i,j)+c(i,j)*beta2n(i,j))
     &       +beta1(i,j)*(a1n(i,j)+a(i,j)*phi1n(i,j))
     &       +beta2(i,j)*(a1t(i,j)+a(i,j)*phi1t(i,j))
c
      bshftdt=2.*(c(i,j)*beta1t(i,j)+b(i,j)*beta2t(i,j))
     &       +beta1(i,j)*(b1n(i,j)+b(i,j)*phi1n(i,j))
     &       +beta2(i,j)*(b1t(i,j)+b(i,j)*phi1t(i,j))
c
      cshftdt=c(i,j)*(beta1n(i,j)+beta2t(i,j))
     &       +b(i,j)*beta2n(i,j)+a(i,j)*beta1t(i,j)
     &       +beta1(i,j)*(c1n(i,j)+c(i,j)*phi1n(i,j))
     &       +beta2(i,j)*(c1t(i,j)+c(i,j)*phi1t(i,j))
c
      dshftdt=beta1(i,j)*(d1n(i,j)+d(i,j)*phi1n(i,j))
     &       +beta2(i,j)*(d1t(i,j)+d(i,j)*(2.*ct(j)+phi1t(i,j)))
c
      adot(i,j)=ageomdt+ashftdt
      bdot(i,j)=bgeomdt+bshftdt
c     cdot(i,j)=cgeomdt+cshftdt
      ddot(i,j)=dgeomdt+dshftdt
   10 continue
c
c        and on axis
      do 20 i=1,netamax
      trk=ha(i,1)/a(i,1)+2.*hb(i,1)/b(i,1)
c
      g111axis=0.5*(a1n(i,1)/a(i,1)+phi1n(i,1))
      g122axis=-0.5*(b1n(i,1)+b(i,1)*phi1n(i,1))/a(i,1)
     &        +c1t(i,1)/a(i,1)
c
      hageomdt=psim4(i,1)*(alpha(i,1)*ricci11(i,1)-alpha2n(i,1)
     &                    +g111axis*alpha1n(i,1))
     &        +alpha(i,1)*(ha(i,1)*trk-2.*ha(i,1)*ha(i,1)/a(i,1))
c
      hbgeomdt=psim4(i,1)*(alpha(i,1)*ricci22(i,1)-alpha2t(i,1)
     &                    +g122axis*alpha1n(i,1))
     &        +alpha(i,1)*(hb(i,1)*trk-2.*hb(i,1)*hb(i,1)/b(i,1))
c
      hashftdt=beta1(i,1)*(ha1n(i,1)+ha(i,1)*phi1n(i,1))
     &        +2.*ha(i,1)*beta1n(i,1)
c
      hbshftdt=beta1(i,1)*(hb1n(i,1)+hb(i,1)*phi1n(i,1))
     &        +2.*hb(i,1)*beta2t(i,1)
c
      hadot(i,1)=hageomdt+hashftdt
      hbdot(i,1)=hbgeomdt+hbshftdt
      hddot(i,1)=hbdot(i,1)
      hcdot(i,1)=0.
c
c
      ageomdt=-2.*alpha(i,1)*ha(i,1)
      bgeomdt=-2.*alpha(i,1)*hb(i,1)
c
      ashftdt=beta1(i,1)*(a1n(i,1)+a(i,1)*phi1n(i,1))
     &       +2.*a(i,1)*beta1n(i,1)
c
      bshftdt=beta1(i,1)*(b1n(i,1)+b(i,1)*phi1n(i,1))
     &       +2.*b(i,1)*beta2t(i,1)
c
      adot(i,1)=ageomdt+ashftdt
      bdot(i,1)=bgeomdt+bshftdt
      cdot(i,1)=0.
      ddot(i,1)=bdot(i,1)
   20 continue
c
c
      return
      end
c
c
c
      subroutine macsmdot
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /hdotblck/ hadot(id,jd),hbdot(id,jd),
     &                  hcdot(id,jd),hddot(id,jd)
      common /fdotblck/ adot(id,jd),bdot(id,jd),
     &                  cdot(id,jd),ddot(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /shiftblk/ beta1(id,jd),beta2(id,jd),omega(id,jd)
      common /shftdblk/ beta1n(id,jd),beta1t(id,jd),
     &                  beta2n(id,jd),beta2t(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      real r(38),s(38),p(38),t(38)
c
c
c        construct time derivatives off axis
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      r(1) =-0.5*a(i,j)*alpha(i,j)*phi1t(i,j)*ct(j)*psim4(i,j)/b(i,j)
      r(2) =-0.5*a1t(i,j)*alpha(i,j)*psim4(i,j)*ct(j)/b(i,j)
      r(3) =-0.5*a(i,j)*alpha(i,j)*psim4(i,j)*phi2t(i,j)/b(i,j)
      r(4) =-0.25*a(i,j)*alpha(i,j)*psim4(i,j)*(phi1t(i,j)**2.)/b(i,j)
      r(5) =-0.25*a(i,j)*alpha(i,j)*d1t(i,j)*psim4(i,j)*phi1t(i,j)
     &     /(b(i,j)*d(i,j))
      r(6) =0.25*a(i,j)*alpha(i,j)*b1t(i,j)*psim4(i,j)*phi1t(i,j)
     &     /(b(i,j)*b(i,j))
      r(7) =-0.5*a(i,j)*alpha1t(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      r(8) =-0.5*a1t(i,j)*alpha(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      r(9) =beta2(i,j)*ha(i,j)*phi1t(i,j)
      r(10)=-alpha(i,j)*psim4(i,j)*phi2n(i,j)
      r(11)=-0.25*alpha(i,j)*d1n(i,j)*psim4(i,j)*phi1n(i,j)/d(i,j)
      r(12)=-0.25*alpha(i,j)*b1n(i,j)*psim4(i,j)*phi1n(i,j)/b(i,j)
      r(13)=0.5*alpha1n(i,j)*psim4(i,j)*phi1n(i,j)
      r(14)=0.5*a1n(i,j)*alpha(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      r(15)=beta1(i,j)*ha(i,j)*phi1n(i,j)
      r(16)=-0.25*a1t(i,j)*alpha(i,j)*d1t(i,j)*psim4(i,j)
     &     /(b(i,j)*d(i,j))
      r(17)=-0.5*alpha(i,j)*d2n(i,j)*psim4(i,j)/d(i,j)
      r(18)=0.25*alpha(i,j)*(d1n(i,j)**2.)*psim4(i,j)/(d(i,j)*d(i,j))
      r(19)=0.25*a1n(i,j)*alpha(i,j)*d1n(i,j)*psim4(i,j)/(a(i,j)*d(i,j))
      r(20)=0.25*a1t(i,j)*alpha(i,j)*b1t(i,j)*psim4(i,j)/(b(i,j)*b(i,j))
      r(21)=-0.5*alpha(i,j)*b2n(i,j)*psim4(i,j)/b(i,j)
      r(22)=0.25*alpha(i,j)*(b1n(i,j)**2.)*psim4(i,j)/(b(i,j)*b(i,j))
      r(23)=0.25*a1n(i,j)*alpha(i,j)*b1n(i,j)*psim4(i,j)/(a(i,j)*b(i,j))
      r(24)=-0.5*a1t(i,j)*alpha1t(i,j)*psim4(i,j)/b(i,j)
      r(25)=-0.5*a2t(i,j)*alpha(i,j)*psim4(i,j)/b(i,j)
      r(26)=0.25*(a1t(i,j)**2.)*alpha(i,j)*psim4(i,j)/(a(i,j)*b(i,j))
      r(27)=-alpha2n(i,j)*psim4(i,j)
      r(28)=0.5*a1n(i,j)*alpha1n(i,j)*psim4(i,j)/a(i,j)
      r(29)=alpha(i,j)*ha(i,j)*hd(i,j)/d(i,j)
      r(30)=-2.*alpha(i,j)*hc(i,j)*hc(i,j)/b(i,j)
      r(31)=2.*beta2n(i,j)*hc(i,j)
      r(32)=alpha(i,j)*ha(i,j)*hb(i,j)/b(i,j)
      r(33)=beta2(i,j)*ha1t(i,j)
      r(34)=beta1(i,j)*ha1n(i,j)
      r(35)=-alpha(i,j)*ha(i,j)*ha(i,j)/a(i,j)
      r(36)=2.*beta1n(i,j)*ha(i,j)
c
      hadot(i,j)= r(1)+ r(2)+ r(3)+ r(4)+ r(5)+ r(6)+ r(7)+ r(8)+ r(9)
     &    +r(10)+r(11)+r(12)+r(13)+r(14)+r(15)+r(16)+r(17)+r(18)+r(19)
     &    +r(20)+r(21)+r(22)+r(23)+r(24)+r(25)+r(26)+r(27)+r(28)+r(29)
     &    +r(30)+r(31)+r(32)+r(33)+r(34)+r(35)+r(36)
   10 continue
c
      do 11 j=2,nthtamax
      do 11 i=1,netamax
      t(1) =-0.5*alpha(i,j)*psim4(i,j)*phi1t(i,j)*ct(j)
      t(2) =-alpha(i,j)*d1t(i,j)*psim4(i,j)*ct(j)/d(i,j)
      t(3) =0.5*alpha(i,j)*b1t(i,j)*psim4(i,j)*ct(j)/b(i,j)
      t(4) =-alpha(i,j)*psim4(i,j)*phi2t(i,j)
      t(5) =-0.25*alpha(i,j)*d1t(i,j)*psim4(i,j)*phi1t(i,j)/d(i,j)
      t(6) =0.5*alpha(i,j)*b1t(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      t(7) =0.5*alpha1t(i,j)*psim4(i,j)*phi1t(i,j)
      t(8) =-0.25*a1t(i,j)*alpha(i,j)*psim4(i,j)*phi1t(i,j)/a(i,j)
      t(9) =beta2(i,j)*hb(i,j)*phi1t(i,j)
      t(10)=-0.5*alpha(i,j)*b(i,j)*psim4(i,j)*phi2n(i,j)/a(i,j)
      t(11)=-0.25*alpha(i,j)*b(i,j)*psim4(i,j)*(phi1n(i,j)**2.)/a(i,j)
      t(12)=-0.25*alpha(i,j)*b(i,j)*d1n(i,j)*psim4(i,j)*phi1n(i,j)
     &     /(a(i,j)*d(i,j))
      t(13)=-0.5*alpha(i,j)*b1n(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      t(14)=-0.5*alpha1n(i,j)*b(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      t(15)=0.25*a1n(i,j)*alpha(i,j)*b(i,j)*psim4(i,j)*phi1n(i,j)
     &     /(a(i,j)*a(i,j))
      t(16)=beta1(i,j)*hb(i,j)*phi1n(i,j)
      t(17)=-0.5*alpha(i,j)*d2t(i,j)*psim4(i,j)/d(i,j)
      t(18)=0.25*alpha(i,j)*d1t(i,j)*d1t(i,j)*psim4(i,j)/(d(i,j)*d(i,j))
      t(19)=0.25*alpha(i,j)*b1t(i,j)*d1t(i,j)*psim4(i,j)/(b(i,j)*d(i,j))
      t(20)=-0.25*alpha(i,j)*b1n(i,j)*d1n(i,j)*psim4(i,j)
     &     /(a(i,j)*d(i,j))
      t(21)=0.5*alpha1t(i,j)*b1t(i,j)*psim4(i,j)/b(i,j)
      t(22)=0.25*a1t(i,j)*alpha(i,j)*b1t(i,j)*psim4(i,j)/(a(i,j)*b(i,j))
      t(23)=-0.5*alpha(i,j)*b2n(i,j)*psim4(i,j)/a(i,j)
      t(24)=0.25*alpha(i,j)*(b1n(i,j)**2.)*psim4(i,j)/(a(i,j)*b(i,j))
      t(25)=-0.5*alpha1n(i,j)*b1n(i,j)*psim4(i,j)/a(i,j)
      t(26)=0.25*a1n(i,j)*alpha(i,j)*b1n(i,j)*psim4(i,j)/(a(i,j)*a(i,j))
      t(27)=-alpha2t(i,j)*psim4(i,j)
      t(28)=-0.5*a2t(i,j)*alpha(i,j)*psim4(i,j)/a(i,j)
      t(29)=0.25*(a1t(i,j)**2.)*alpha(i,j)*psim4(i,j)/(a(i,j)*a(i,j))
      t(30)=alpha(i,j)*psim4(i,j)
      t(31)=alpha(i,j)*hb(i,j)*hd(i,j)/d(i,j)
      t(32)=-2.*alpha(i,j)*hc(i,j)*hc(i,j)/a(i,j)
      t(33)=2.*beta1t(i,j)*hc(i,j)
      t(34)=beta2(i,j)*hb1t(i,j)
      t(35)=beta1(i,j)*hb1n(i,j)
      t(36)=-alpha(i,j)*hb(i,j)*hb(i,j)/b(i,j)
      t(37)=alpha(i,j)*ha(i,j)*hb(i,j)/a(i,j)
      t(38)=2.*beta2t(i,j)*hb(i,j)
c
      hbdot(i,j)= t(1)+ t(2)+ t(3)+ t(4)+ t(5)+ t(6)+ t(7)+ t(8)+ t(9)
     &    +t(10)+t(11)+t(12)+t(13)+t(14)+t(15)+t(16)+t(17)+t(18)+t(19)
     &    +t(20)+t(21)+t(22)+t(23)+t(24)+t(25)+t(26)+t(27)+t(28)+t(29)
     &    +t(30)+t(31)+t(32)+t(33)+t(34)+t(35)+t(36)+t(37)+t(38)
   11 continue
c
      do 12 j=2,nthtamax
      do 12 i=1,netamax
      s(1) =-0.5*alpha(i,j)*d1n(i,j)*psim4(i,j)*ct(j)/d(i,j)
      s(2) =0.5*alpha(i,j)*b1n(i,j)*psim4(i,j)*ct(j)/b(i,j)
      s(3) =0.25*alpha(i,j)*psim4(i,j)*phi1n(i,j)*phi1t(i,j)
      s(4) =0.25*alpha(i,j)*b1n(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      s(5) =0.5*alpha1n(i,j)*psim4(i,j)*phi1t(i,j)
      s(6) =beta2(i,j)*hc(i,j)*phi1t(i,j)
      s(7) =-0.5*alpha(i,j)*psim4(i,j)*phi2m(i,j)
      s(8) =0.5*alpha1t(i,j)*psim4(i,j)*phi1n(i,j)
      s(9) =0.25*a1t(i,j)*alpha(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      s(10)=beta1(i,j)*hc(i,j)*phi1n(i,j)
      s(11)=0.25*alpha(i,j)*d1n(i,j)*d1t(i,j)*psim4(i,j)/(d(i,j)*d(i,j))
      s(12)=0.25*alpha(i,j)*b1n(i,j)*d1t(i,j)*psim4(i,j)/(b(i,j)*d(i,j))
      s(13)=-0.5*alpha(i,j)*d2m(i,j)*psim4(i,j)/d(i,j)
      s(14)=0.25*a1t(i,j)*alpha(i,j)*d1n(i,j)*psim4(i,j)/(a(i,j)*d(i,j))
      s(15)=0.5*alpha1t(i,j)*b1n(i,j)*psim4(i,j)/b(i,j)
      s(16)=-alpha2m(i,j)*psim4(i,j)
      s(17)=0.5*a1t(i,j)*alpha1n(i,j)*psim4(i,j)/a(i,j)
      s(18)=alpha(i,j)*hc(i,j)*hd(i,j)/d(i,j)
      s(19)=beta2(i,j)*hc1t(i,j)
      s(20)=beta1(i,j)*hc1n(i,j)
      s(21)=-alpha(i,j)*hb(i,j)*hc(i,j)/b(i,j)
      s(22)=-alpha(i,j)*ha(i,j)*hc(i,j)/a(i,j)
      s(23)=beta2t(i,j)*hc(i,j)
      s(24)=beta1n(i,j)*hc(i,j)
      s(25)=beta2n(i,j)*hb(i,j)
      s(26)=beta1t(i,j)*ha(i,j)
c
      hcdot(i,j)= s(1)+ s(2)+ s(3)+ s(4)+ s(5)+ s(6)+ s(7)+ s(8)+ s(9)
     &    +s(10)+s(11)+s(12)+s(13)+s(14)+s(15)+s(16)+s(17)+s(18)+s(19)
     &    +s(20)+s(21)+s(22)+s(23)+s(24)+s(25)+s(26)
   12 continue
c
      do 13 j=2,nthtamax
      do 13 i=1,netamax
      p(1) =-alpha(i,j)*d(i,j)*psim4(i,j)*phi1t(i,j)*ct(j)/b(i,j)
      p(2) =-alpha(i,j)*d1t(i,j)*psim4(i,j)*ct(j)/b(i,j)
      p(3) =0.5*alpha(i,j)*b1t(i,j)*d(i,j)*psim4(i,j)*ct(j)
     &   /(b(i,j)*b(i,j))
      p(4) =-alpha1t(i,j)*d(i,j)*psim4(i,j)*ct(j)/b(i,j)
      p(5) =-0.5*a1t(i,j)*alpha(i,j)*d(i,j)*psim4(i,j)*ct(j)
     &   /(a(i,j)*b(i,j))
      p(6) =2.*beta2(i,j)*hd(i,j)*ct(j)
      p(7) =-0.5*alpha(i,j)*d(i,j)*psim4(i,j)*phi2t(i,j)/b(i,j)
      p(8) =-0.25*alpha(i,j)*d(i,j)*psim4(i,j)*(phi1t(i,j)**2.)/b(i,j)
      p(9) =-0.5*alpha(i,j)*d1t(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      p(10)=0.25*alpha(i,j)*b1t(i,j)*d(i,j)*psim4(i,j)*phi1t(i,j)
     &   /(b(i,j)*b(i,j))
      p(11)=-0.5*alpha1t(i,j)*d(i,j)*psim4(i,j)*phi1t(i,j)/b(i,j)
      p(12)=-0.25*a1t(i,j)*alpha(i,j)*d(i,j)*psim4(i,j)*phi1t(i,j)
     &   /(a(i,j)*b(i,j))
      p(13)=beta2(i,j)*hd(i,j)*phi1t(i,j)
      p(14)=-0.5*alpha(i,j)*d(i,j)*psim4(i,j)*phi2n(i,j)/a(i,j)
      p(15)=-0.25*alpha(i,j)*d(i,j)*psim4(i,j)*(phi1n(i,j)**2.)/a(i,j)
      p(16)=-0.5*alpha(i,j)*d1n(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      p(17)=-0.25*alpha(i,j)*b1n(i,j)*d(i,j)*psim4(i,j)*phi1n(i,j)
     &   /(a(i,j)*b(i,j))
      p(18)=-0.5*alpha1n(i,j)*d(i,j)*psim4(i,j)*phi1n(i,j)/a(i,j)
      p(19)=0.25*a1n(i,j)*alpha(i,j)*d(i,j)*psim4(i,j)*phi1n(i,j)
     &   /(a(i,j)*a(i,j))
      p(20)=beta1(i,j)*hd(i,j)*phi1n(i,j)
      p(21)=-0.5*alpha(i,j)*d2t(i,j)*psim4(i,j)/b(i,j)
      p(22)=0.25*alpha(i,j)*(d1t(i,j)**2.)*psim4(i,j)/(b(i,j)*d(i,j))
      p(23)=0.25*alpha(i,j)*b1t(i,j)*d1t(i,j)*psim4(i,j)/(b(i,j)**2.)
      p(24)=-0.5*alpha1t(i,j)*d1t(i,j)*psim4(i,j)/b(i,j)
      p(25)=-0.25*a1t(i,j)*alpha(i,j)*d1t(i,j)*psim4(i,j)
     &     /(a(i,j)*b(i,j))
      p(26)=-0.5*alpha(i,j)*d2n(i,j)*psim4(i,j)/a(i,j)
      p(27)=0.25*alpha(i,j)*(d1n(i,j)**2.)*psim4(i,j)/(a(i,j)*d(i,j))
      p(28)=-0.25*alpha(i,j)*b1n(i,j)*d1n(i,j)*psim4(i,j)
     &     /(a(i,j)*b(i,j))
      p(29)=-0.5*alpha1n(i,j)*d1n(i,j)*psim4(i,j)/a(i,j)
      p(30)=0.25*a1n(i,j)*alpha(i,j)*d1n(i,j)*psim4(i,j)/(a(i,j)**2.)
      p(31)=alpha(i,j)*d(i,j)*psim4(i,j)/b(i,j)
      p(32)=beta2(i,j)*hd1t(i,j)
      p(33)=beta1(i,j)*hd1n(i,j)
      p(34)=-alpha(i,j)*hd(i,j)*hd(i,j)/d(i,j)
      p(35)=alpha(i,j)*hb(i,j)*hd(i,j)/b(i,j)
      p(36)=alpha(i,j)*ha(i,j)*hd(i,j)/a(i,j)
c
      hddot(i,j)= p(1)+ p(2)+ p(3)+ p(4)+ p(5)+ p(6)+ p(7)+ p(8)+ p(9)
     &    +p(10)+p(11)+p(12)+p(13)+p(14)+p(15)+p(16)+p(17)+p(18)+p(19)
     &    +p(20)+p(21)+p(22)+p(23)+p(24)+p(25)+p(26)+p(27)+p(28)+p(29)
     &    +p(30)+p(31)+p(32)+p(33)+p(34)+p(35)+p(36)
   13 continue
c
      do 14 j=2,nthtamax
      do 14 i=1,netamax
      q1=a(i,j)*beta2(i,j)*phi1t(i,j)
      q2=a(i,j)*beta1(i,j)*phi1n(i,j)
      q3=-2.*alpha(i,j)*ha(i,j)
      q4=a1t(i,j)*beta2(i,j)
      q5=2.*a(i,j)*beta1n(i,j)
      q6=a1n(i,j)*beta1(i,j)
c
      adot(i,j)=q1+q2+q3+q4+q5+q6
c
      q1=b(i,j)*beta2(i,j)*phi1t(i,j)
      q2=b(i,j)*beta1(i,j)*phi1n(i,j)
      q3=-2.*alpha(i,j)*hb(i,j)
      q4=2.*b(i,j)*beta2t(i,j)
      q5=b1t(i,j)*beta2(i,j)
      q6=b1n(i,j)*beta1(i,j)
c
      bdot(i,j)=q1+q2+q3+q4+q5+q6
c
      cdot(i,j)=0.
c
      q1=2.*beta2(i,j)*d(i,j)*ct(j)
      q2=beta2(i,j)*d(i,j)*phi1t(i,j)
      q3=beta1(i,j)*d(i,j)*phi1n(i,j)
      q4=-2.*alpha(i,j)*hd(i,j)
      q5=beta2(i,j)*d1t(i,j)
      q6=beta1(i,j)*d1n(i,j)
c
      ddot(i,j)=q1+q2+q3+q4+q5+q6
   14 continue
c
c        and on axis
      do 20 i=1,netamax
      r(1) =-0.5*a(i,1)*alpha(i,1)*phi2t(i,1)*psim4(i,1)/b(i,1)
      r(2) =-0.5*a2t(i,1)*alpha(i,1)*psim4(i,1)/b(i,1)
      r(3) =-0.5*a(i,1)*alpha(i,1)*psim4(i,1)*phi2t(i,1)/b(i,1)
      r(4) =-0.25*a(i,1)*alpha(i,1)*psim4(i,1)*(phi1t(i,1)**2.)/b(i,1)
      r(5) =-0.25*a(i,1)*alpha(i,1)*d1t(i,1)*psim4(i,1)*phi1t(i,1)
     &     /(b(i,1)*d(i,1))
      r(6) =0.25*a(i,1)*alpha(i,1)*b1t(i,1)*psim4(i,1)*phi1t(i,1)
     &     /(b(i,1)*b(i,1))
      r(7) =-0.5*a(i,1)*alpha1t(i,1)*psim4(i,1)*phi1t(i,1)/b(i,1)
      r(8) =-0.5*a1t(i,1)*alpha(i,1)*psim4(i,1)*phi1t(i,1)/b(i,1)
      r(9) =beta2(i,1)*ha(i,1)*phi1t(i,1)
      r(10)=-alpha(i,1)*psim4(i,1)*phi2n(i,1)
      r(11)=-0.25*alpha(i,1)*d1n(i,1)*psim4(i,1)*phi1n(i,1)/d(i,1)
      r(12)=-0.25*alpha(i,1)*b1n(i,1)*psim4(i,1)*phi1n(i,1)/b(i,1)
      r(13)=0.5*alpha1n(i,1)*psim4(i,1)*phi1n(i,1)
      r(14)=0.5*a1n(i,1)*alpha(i,1)*psim4(i,1)*phi1n(i,1)/a(i,1)
      r(15)=beta1(i,1)*ha(i,1)*phi1n(i,1)
      r(16)=-0.25*a1t(i,1)*alpha(i,1)*d1t(i,1)*psim4(i,1)
     &     /(b(i,1)*d(i,1))
      r(17)=-0.5*alpha(i,1)*d2n(i,1)*psim4(i,1)/d(i,1)
      r(18)=0.25*alpha(i,1)*(d1n(i,1)**2.)*psim4(i,1)/(d(i,1)*d(i,1))
      r(19)=0.25*a1n(i,1)*alpha(i,1)*d1n(i,1)*psim4(i,1)/(a(i,1)*d(i,1))
      r(20)=0.25*a1t(i,1)*alpha(i,1)*b1t(i,1)*psim4(i,1)/(b(i,1)*b(i,1))
      r(21)=-0.5*alpha(i,1)*b2n(i,1)*psim4(i,1)/b(i,1)
      r(22)=0.25*alpha(i,1)*(b1n(i,1)**2.)*psim4(i,1)/(b(i,1)*b(i,1))
      r(23)=0.25*a1n(i,1)*alpha(i,1)*b1n(i,1)*psim4(i,1)/(a(i,1)*b(i,1))
      r(24)=-0.5*a1t(i,1)*alpha1t(i,1)*psim4(i,1)/b(i,1)
      r(25)=-0.5*a2t(i,1)*alpha(i,1)*psim4(i,1)/b(i,1)
      r(26)=0.25*(a1t(i,1)**2.)*alpha(i,1)*psim4(i,1)/(a(i,1)*b(i,1))
      r(27)=-alpha2n(i,1)*psim4(i,1)
      r(28)=0.5*a1n(i,1)*alpha1n(i,1)*psim4(i,1)/a(i,1)
      r(29)=alpha(i,1)*ha(i,1)*hd(i,1)/d(i,1)
      r(30)=-2.*alpha(i,1)*hc(i,1)*hc(i,1)/b(i,1)
      r(31)=2.*beta2n(i,1)*hc(i,1)
      r(32)=alpha(i,1)*ha(i,1)*hb(i,1)/b(i,1)
      r(33)=beta2(i,1)*ha1t(i,1)
      r(34)=beta1(i,1)*ha1n(i,1)
      r(35)=-alpha(i,1)*ha(i,1)*ha(i,1)/a(i,1)
      r(36)=2.*beta1n(i,1)*ha(i,1)
c
      hadot(i,1)= r(1)+ r(2)+ r(3)+ r(4)+ r(5)+ r(6)+ r(7)+ r(8)+ r(9)
     &    +r(10)+r(11)+r(12)+r(13)+r(14)+r(15)+r(16)+r(17)+r(18)+r(19)
     &    +r(20)+r(21)+r(22)+r(23)+r(24)+r(25)+r(26)+r(27)+r(28)+r(29)
     &    +r(30)+r(31)+r(32)+r(33)+r(34)+r(35)+r(36)
   20 continue
c
      do 21 i=1,netamax
      t(1) =-0.5*alpha(i,1)*psim4(i,1)*phi2t(i,1)
      t(2) =-alpha(i,1)*d2t(i,1)*psim4(i,1)/d(i,1)
      t(3) =0.5*alpha(i,1)*b2t(i,1)*psim4(i,1)/b(i,1)
      t(4) =-alpha(i,1)*psim4(i,1)*phi2t(i,1)
      t(5) =-0.25*alpha(i,1)*d1t(i,1)*psim4(i,1)*phi1t(i,1)/d(i,1)
      t(6) =0.5*alpha(i,1)*b1t(i,1)*psim4(i,1)*phi1t(i,1)/b(i,1)
      t(7) =0.5*alpha1t(i,1)*psim4(i,1)*phi1t(i,1)
      t(8) =-0.25*a1t(i,1)*alpha(i,1)*psim4(i,1)*phi1t(i,1)/a(i,1)
      t(9) =beta2(i,1)*hb(i,1)*phi1t(i,1)
      t(10)=-0.5*alpha(i,1)*b(i,1)*psim4(i,1)*phi2n(i,1)/a(i,1)
      t(11)=-0.25*alpha(i,1)*b(i,1)*psim4(i,1)*(phi1n(i,1)**2.)/a(i,1)
      t(12)=-0.25*alpha(i,1)*b(i,1)*d1n(i,1)*psim4(i,1)*phi1n(i,1)
     &     /(a(i,1)*d(i,1))
      t(13)=-0.5*alpha(i,1)*b1n(i,1)*psim4(i,1)*phi1n(i,1)/a(i,1)
      t(14)=-0.5*alpha1n(i,1)*b(i,1)*psim4(i,1)*phi1n(i,1)/a(i,1)
      t(15)=0.25*a1n(i,1)*alpha(i,1)*b(i,1)*psim4(i,1)*phi1n(i,1)
     &     /(a(i,1)*a(i,1))
      t(16)=beta1(i,1)*hb(i,1)*phi1n(i,1)
      t(17)=-0.5*alpha(i,1)*d2t(i,1)*psim4(i,1)/d(i,1)
      t(18)=0.25*alpha(i,1)*d1t(i,1)*d1t(i,1)*psim4(i,1)/(d(i,1)*d(i,1))
      t(19)=0.25*alpha(i,1)*b1t(i,1)*d1t(i,1)*psim4(i,1)/(b(i,1)*d(i,1))
      t(20)=-0.25*alpha(i,1)*b1n(i,1)*d1n(i,1)*psim4(i,1)
     &     /(a(i,1)*d(i,1))
      t(21)=0.5*alpha1t(i,1)*b1t(i,1)*psim4(i,1)/b(i,1)
      t(22)=0.25*a1t(i,1)*alpha(i,1)*b1t(i,1)*psim4(i,1)/(a(i,1)*b(i,1))
      t(23)=-0.5*alpha(i,1)*b2n(i,1)*psim4(i,1)/a(i,1)
      t(24)=0.25*alpha(i,1)*(b1n(i,1)**2.)*psim4(i,1)/(a(i,1)*b(i,1))
      t(25)=-0.5*alpha1n(i,1)*b1n(i,1)*psim4(i,1)/a(i,1)
      t(26)=0.25*a1n(i,1)*alpha(i,1)*b1n(i,1)*psim4(i,1)/(a(i,1)*a(i,1))
      t(27)=-alpha2t(i,1)*psim4(i,1)
      t(28)=-0.5*a2t(i,1)*alpha(i,1)*psim4(i,1)/a(i,1)
      t(29)=0.25*(a1t(i,1)**2.)*alpha(i,1)*psim4(i,1)/(a(i,1)*a(i,1))
      t(30)=alpha(i,1)*psim4(i,1)
      t(31)=alpha(i,1)*hb(i,1)*hd(i,1)/d(i,1)
      t(32)=-2.*alpha(i,1)*hc(i,1)*hc(i,1)/a(i,1)
      t(33)=2.*beta1t(i,1)*hc(i,1)
      t(34)=beta2(i,1)*hb1t(i,1)
      t(35)=beta1(i,1)*hb1n(i,1)
      t(36)=-alpha(i,1)*hb(i,1)*hb(i,1)/b(i,1)
      t(37)=alpha(i,1)*ha(i,1)*hb(i,1)/a(i,1)
      t(38)=2.*beta2t(i,1)*hb(i,1)
c
      hbdot(i,1)= t(1)+ t(2)+ t(3)+ t(4)+ t(5)+ t(6)+ t(7)+ t(8)+ t(9)
     &    +t(10)+t(11)+t(12)+t(13)+t(14)+t(15)+t(16)+t(17)+t(18)+t(19)
     &    +t(20)+t(21)+t(22)+t(23)+t(24)+t(25)+t(26)+t(27)+t(28)+t(29)
     &    +t(30)+t(31)+t(32)+t(33)+t(34)+t(35)+t(36)+t(37)+t(38)
   21 continue
c
      do 22 i=1,netamax
      hcdot(i,1)=0.
c
      hddot(i,1)=hbdot(i,1)
c
      q1=a(i,1)*beta2(i,1)*phi1t(i,1)
      q2=a(i,1)*beta1(i,1)*phi1n(i,1)
      q3=-2.*alpha(i,1)*ha(i,1)
      q4=a1t(i,1)*beta2(i,1)
      q5=2.*a(i,1)*beta1n(i,1)
      q6=a1n(i,1)*beta1(i,1)
c
      adot(i,1)=q1+q2+q3+q4+q5+q6
c
      q1=b(i,1)*beta2(i,1)*phi1t(i,1)
      q2=b(i,1)*beta1(i,1)*phi1n(i,1)
      q3=-2.*alpha(i,1)*hb(i,1)
      q4=2.*b(i,1)*beta2t(i,1)
      q5=b1t(i,1)*beta2(i,1)
      q6=b1n(i,1)*beta1(i,1)
c
      bdot(i,1)=q1+q2+q3+q4+q5+q6
c
      cdot(i,1)=0.
c
      ddot(i,1)=bdot(i,1)
   22 continue
c
c
      return
      end
c
c
c
      subroutine maxlapse
      include 'param.h'
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c        construct coefficient functions used in the maximal slicing
c        equation
      call setlpsf
c
      call setmc5(1.0,alphabnd)
c
      tol=1.0e-9
c
      call multigrid(netamax,1,netamax,nthtamax,1,nthtamax,
     &               cc,cn,cs,cw,ce,cnw,cne,csw,cse,alpha,rhs,tol,
     &               0,iflag)
c
      if (iflag.eq.0) then
      write(5,*)'multigrid did not converge for maximal slicing'
      stop
      endif
c
c        compute lapse derivatives
      if (dervans.eq.'s') then
      call scnddrvs(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      elseif (dervans.eq.'f') then
      call frthdrvs(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      elseif (dervans.eq.'s24') then
      call scheme24(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      endif
c
c
      return
      end
c
c
c
      subroutine setlpsf
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /chrisblk/ g111(id,jd),g112(id,jd),g122(id,jd),
     &                  g133(id,jd),g211(id,jd),g212(id,jd),
     &                  g222(id,jd),g233(id,jd),g313(id,jd),
     &                  g323(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        off axis coefficients
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      f11(i,j)=b(i,j)*em1(i,j)
      f22(i,j)=xi*xi*a(i,j)*em1(i,j)
      f12(i,j)=-0.5*xi*c(i,j)*em1(i,j)
c
      f1(i,j)=-0.5*deta*(em1(i,j)*(g111(i,j)*b(i,j)+g122(i,j)*a(i,j)
     &                            -2.*g112(i,j)*c(i,j))
     &                  +g133(i,j)/d(i,j))
c
      f2(i,j)=-0.5*deta*xi*(em1(i,j)*(g211(i,j)*b(i,j)+g222(i,j)*a(i,j)
     &                               -2.*g212(i,j)*c(i,j))
     &                     +g233(i,j)/d(i,j))
c
      r1=b(i,j)*b(i,j)*ha(i,j)*ha(i,j)
      r2=a(i,j)*a(i,j)*hb(i,j)*hb(i,j)
      r3=2.*(a(i,j)*b(i,j)+c(i,j)*c(i,j))*hc(i,j)*hc(i,j)
      r4=-4.*b(i,j)*c(i,j)*ha(i,j)*hc(i,j)
      r5=-4.*a(i,j)*c(i,j)*hb(i,j)*hc(i,j)
      r6=2.*c(i,j)*c(i,j)*ha(i,j)*hb(i,j)
      r7=hd(i,j)*hd(i,j)/(d(i,j)*d(i,j))
c
      sum=r1+r2+r3+r4+r5+r6
c
      f(i,j)=-deta*deta*(em1(i,j)*em1(i,j)*sum+r7)/psim4(i,j)
c
      g(i,j)=0.
   10 continue
c
c        on axis coefficients
      do 20 i=1,netamax
      f11(i,1)=1./a(i,1)
      f22(i,1)=2.*xi*xi/b(i,1)
      f12(i,1)=0.
c
      g111axis=0.5*(a1n(i,1)/a(i,1)+phi1n(i,1))
      g122axis=0.5*(2.*c1t(i,1)-b1n(i,1)-b(i,1)*phi1n(i,1))/a(i,1)
c
      f1(i,1)=-0.5*deta*(g111axis/a(i,1)+2.*g122axis/b(i,1))
      f2(i,1)=0.
c
      r1=(ha(i,1)/a(i,1))**2+2.*(hb(i,1)/b(i,1))**2
      f(i,1)=-deta*deta*r1/psim4(i,1)
c
      g(i,1)=0.
   20 continue
c
c
      return
      end
c
c
c
      subroutine alglapse
      include 'param.h'
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
      call computeq
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      detginit=exp(4.*f(i,j))
      r=a(i,j)*b(i,j)*d(i,j)/detginit
c
      alpha(i,j)=1.+r-1./r
   10 continue
c
c        compute lapse derivatives
      if (dervans.eq.'s') then
      call scnddrvs(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      elseif (dervans.eq.'f') then
      call frthdrvs(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      elseif (dervans.eq.'s24') then
      call scheme24(alpha,alpha1n,alpha1t,alpha2n,alpha2t,alpha2m,
     &              1.,alphabnd,netamax,nthtamax)
      endif
c
c
      return
      end
c
c
c
      subroutine shiftc0
      include 'param.h'
      common /shiftblk/ beta1(id,jd),beta2(id,jd),omega(id,jd)
      common /shftdblk/ beta1n(id,jd),beta1t(id,jd),
     &                  beta2n(id,jd),beta2t(id,jd)
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        this subroutine computes the up components of the
c        shift which diagonalizes the 3-metric.
c
      call setshftf
c
      call setmc5(-1.0,omegabnd)
c
      tol=1.0e-10
c
      call multigrid(netamax,1,netamax,nthtamax,1,nthtamax,
     &               cc,cn,cs,cw,ce,cnw,cne,csw,cse,omega,rhs,tol,
     &               0,iflag)
c
      if (iflag.eq.0) then
      write(5,*)'multigrid did not converge for shift'
      stop
      endif
c
c        compute derivatives
      if (dervans.eq.'s') then
      call scnddrvs(omega,beta2,beta1,beta2n,beta1t,beta1n,
     &              -1.,omegabnd,netamax,nthtamax)
      elseif (dervans.eq.'f') then
      call frthdrvs(omega,beta2,beta1,beta2n,beta1t,beta1n,
     &              -1.,omegabnd,netamax,nthtamax)
      elseif (dervans.eq.'s24') then
      call scheme24(omega,beta2,beta1,beta2n,beta1t,beta1n,
     &              -1.,omegabnd,netamax,nthtamax)
      endif
c
      do 50 j=1,nthtamax
      do 50 i=1,netamax
   50 beta2t(i,j)=beta1n(i,j)
c
c
      return
      end
c
c
c
      subroutine setshftf
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /lapseblk/ alpha(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      f11(i,j)=b(i,j)
      f22(i,j)=xi*xi*a(i,j)
      f12(i,j)=0.5*xi*c(i,j)
c
      f1(i,j)=0.5*deta*(c1t(i,j)*c(i,j)*phi1t(i,j))
      f2(i,j)=0.5*deta*xi*(c1n(i,j)+c(i,j)*phi1n(i,j))
c
      f(i,j)=0.
c
      g(i,j)=-2.*alpha(i,j)*hc(i,j)*deta*deta
   10 continue
c
c
      return
      end
c
c
c
      subroutine setmc5(sym,bnd)
      include 'param.h'
      dimension bnd(*)
c
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
c        set rhs
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      rhs(i,j)=-g(i,j)
   10 continue
      do 11 j=1,nthtamax
      rhs(netamax,j)=-bnd(j)*(f11(netamax,j)+f1(netamax,j))-g(netamax,j)
   11 continue
c
c        zero all components
      do 20 j=1,nthtamax
      do 20 i=1,netamax
      cc(i,j)=0.
      ce(i,j)=0.
      cs(i,j)=0.
      cw(i,j)=0.
      cn(i,j)=0.
      cnw(i,j)=0.
      cne(i,j)=0.
      csw(i,j)=0.
      cse(i,j)=0.
   20 continue
c
c        all elements of cc
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      cc(i,j)=-2.*(f11(i,j)+f22(i,j))+f(i,j)
   30 continue
c
c        axis
      ce(1,1)=f11(1,1)+f1(1,1)+sym*(f11(1,1)-f1(1,1))
      cn(1,1)=f22(1,1)+f2(1,1)+sym*(f22(1,1)-f2(1,1))
c
      do 40 i=2,netamax-1
      ce(i,1)=f11(i,1)+f1(i,1)
      cw(i,1)=f11(i,1)-f1(i,1)
      cn(i,1)=f22(i,1)+f2(i,1)+sym*(f22(i,1)-f2(i,1))
   40 continue
c
      cw(netamax,1)=f11(netamax,1)-f1(netamax,1)
      cn(netamax,1)=f22(netamax,1)+f2(netamax,1)
     &             +sym*(f22(netamax,1)-f2(netamax,1))
c
c        body
      do 60 j=2,nthtamax-1
      cs(1,j)=f22(1,j)-f2(1,j)
      cn(1,j)=f22(1,j)+f2(1,j)
      ce(1,j)=f11(1,j)+f1(1,j)+sym*(f11(1,j)-f1(1,j))
c
      do 50 i=2,netamax-1
      cn(i,j)=f22(i,j)+f2(i,j)
      cs(i,j)=f22(i,j)-f2(i,j)
      ce(i,j)=f11(i,j)+f1(i,j)
      cw(i,j)=f11(i,j)-f1(i,j)
   50 continue
c
      cn(netamax,j)=f22(netamax,j)+f2(netamax,j)
      cs(netamax,j)=f22(netamax,j)-f2(netamax,j)
      cw(netamax,j)=f11(netamax,j)-f1(netamax,j)
   60 continue
c
      cs(1,nthtamax)=f22(1,nthtamax)-f2(1,nthtamax)
     &              +sym*(f22(1,nthtamax)+f2(1,nthtamax))
      ce(1,nthtamax)=f11(1,nthtamax)+f1(1,nthtamax)
     &              +sym*(f11(1,nthtamax)-f1(1,nthtamax))
c
      do 70 i=2,netamax-1
      cs(i,nthtamax)=f22(i,nthtamax)-f2(i,nthtamax)
     &              +sym*(f22(i,nthtamax)+f2(i,nthtamax))
      cw(i,nthtamax)=f11(i,nthtamax)-f1(i,nthtamax)
      ce(i,nthtamax)=f11(i,nthtamax)+f1(i,nthtamax)
   70 continue
c
      cs(netamax,nthtamax)=f22(netamax,nthtamax)-f2(netamax,nthtamax)
     &                +sym*(f22(netamax,nthtamax)+f2(netamax,nthtamax))
      cw(netamax,nthtamax)=f11(netamax,nthtamax)-f1(netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine setmcivr
      include 'param.h'
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        this subroutine implements Robin boundary conditions
c        for the conformal factor
c
c        set rhs 
      bnd=2.*deta*sqrthfm*exp(0.5*eta(netamax))
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      rhs(i,j)=0.
   10 continue
      do 11 j=1,nthtamax
      rhs(netamax,j)=-bnd*f11(netamax,j)
   11 continue
c
c        zero all components
      do 20 j=1,nthtamax
      do 20 i=1,netamax
      cc(i,j)=0.
      ce(i,j)=0.
      cs(i,j)=0.
      cw(i,j)=0.
      cn(i,j)=0.
      cnw(i,j)=0.
      cne(i,j)=0.
      csw(i,j)=0.
      cse(i,j)=0.
   20 continue
c
c        all elements of cc
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      cc(i,j)=-2.*(f11(i,j)+f22(i,j))+f(i,j)
   30 continue
      do 31 j=1,nthtamax
      cc(netamax,j)=cc(netamax,j)-deta*f11(netamax,j)
   31 continue
c
c        axis
      ce(1,1)=2.*f11(1,1)
      cn(1,1)=2.*f22(1,1)
c
      do 40 i=2,netamax-1
      ce(i,1)=f11(i,1)
      cw(i,1)=f11(i,1)
      cn(i,1)=2.*f22(i,1)
   40 continue
c
      cw(netamax,1)=2.*f11(netamax,1)
      cn(netamax,1)=2.*f22(netamax,1)
c
c        body
      do 60 j=2,nthtamax-1
      cs(1,j)=f22(1,j)-f2(1,j)
      ce(1,j)=2.*f11(1,j)
      cn(1,j)=f22(1,j)+f2(1,j)
c
      do 50 i=2,netamax-1
      cn(i,j)=f22(i,j)+f2(i,j)
      cs(i,j)=f22(i,j)-f2(i,j)
      ce(i,j)=f11(i,j)
      cw(i,j)=f11(i,j)
   50 continue
c
      cs(netamax,j)=f22(netamax,j)-f2(netamax,j)
      cn(netamax,j)=f22(netamax,j)+f2(netamax,j)
      cw(netamax,j)=2.*f11(netamax,j)
   60 continue
c
      cs(1,nthtamax)=2.*f22(1,nthtamax)
      ce(1,nthtamax)=2.*f11(1,nthtamax)
c
      do 70 i=2,netamax-1
      cs(i,nthtamax)=2.*f22(i,nthtamax)
      ce(i,nthtamax)=f11(i,nthtamax)
      cw(i,nthtamax)=f11(i,nthtamax)
   70 continue
c
      cs(netamax,nthtamax)=2.*f22(netamax,nthtamax)
      cw(netamax,nthtamax)=2.*f11(netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine setmciv3
      include 'param.h'
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        this subroutine implements Robin boundary conditions
c        for the conformal factor
c
c        set rhs
      fac1=1./(1.+exp(-eta(netamax)))
      fac2=2.*deta*fac1
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      rhs(i,j)=-g(i,j)
   10 continue
      do 11 j=1,nthtamax
      rhs(netamax,j)=-g(netamax,j)-fac2*(f11(netamax,j)+f1(netamax,j))
   11 continue
c
c        zero all components
      do 20 j=1,nthtamax
      do 20 i=1,netamax
      cc(i,j)=0.
      ce(i,j)=0.
      cs(i,j)=0.
      cw(i,j)=0.
      cn(i,j)=0.
      cnw(i,j)=0.
      cne(i,j)=0.
      csw(i,j)=0.
      cse(i,j)=0.
   20 continue
c
c        all elements of cc
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      cc(i,j)=-2.*(f11(i,j)+f22(i,j))+f(i,j)
   30 continue
      do 35 j=1,nthtamax
      cc(netamax,j)=cc(netamax,j)-fac2*(f11(netamax,j)+f1(netamax,j))
   35 continue
c
c        axis
      ce(1,1)=2.*f11(1,1)
      cn(1,1)=2.*f22(1,1)
c
      do 40 i=2,netamax-1
      ce(i,1)=f11(i,1)+f1(i,1)
      cw(i,1)=f11(i,1)-f1(i,1)
      cn(i,1)=2.*f22(i,1)
   40 continue
c
      cw(netamax,1)=2.*f11(netamax,1)
      cn(netamax,1)=2.*f22(netamax,1)
c
c        body
      do 60 j=2,nthtamax-1
      cn(1,j)=f22(1,j)+f2(1,j)
      cs(1,j)=f22(1,j)-f2(1,j)
      ce(1,j)=2.*f11(1,j)
c
      do 50 i=2,netamax-1
      cn(i,j)=f22(i,j)+f2(i,j)
      cs(i,j)=f22(i,j)-f2(i,j)
      ce(i,j)=f11(i,j)+f1(i,j)
      cw(i,j)=f11(i,j)-f1(i,j)
   50 continue
c
      cn(netamax,j)=f22(netamax,j)+f2(netamax,j)
      cs(netamax,j)=f22(netamax,j)-f2(netamax,j)
      cw(netamax,j)=2.*f11(netamax,j)
   60 continue
c
      cs(1,nthtamax)=2.*f22(1,nthtamax)
      ce(1,nthtamax)=2.*f11(1,nthtamax)
c
      do 70 i=2,netamax-1
      cs(i,nthtamax)=2.*f22(i,nthtamax)
      ce(i,nthtamax)=f11(i,nthtamax)+f1(i,nthtamax)
      cw(i,nthtamax)=f11(i,nthtamax)-f1(i,nthtamax)
   70 continue
c
      cs(netamax,nthtamax)=2.*f22(netamax,nthtamax)
      cw(netamax,nthtamax)=2.*f11(netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine scnddrvs(ux,u1n,u1t,u2n,u2t,u2m,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt), u1n(ne,nt),u1t(ne,nt),
     &          u2n(ne,nt),u2t(ne,nt),u2m(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      call scnddeta(u1n,u2n,netamax,nthtamax)
      call scndthta(u1t,u2t,netamax,nthtamax)
      call scndmixd(u2m    ,netamax,nthtamax)
c
      return
      end
c
c
c
      subroutine scheme24(ux,u1n,u1t,u2n,u2t,u2m,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt), u1n(ne,nt),u1t(ne,nt),
     &          u2n(ne,nt),u2t(ne,nt),u2m(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      call scnddeta(u1n,u2n,netamax,nthtamax)
      call frththta(u1t,u2t,netamax,nthtamax)
      call mixedd24(u2m    ,netamax,nthtamax)
c
      return
      end
c
c
c
      subroutine frthdrvs(ux,u1n,u1t,u2n,u2t,u2m,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt), u1n(ne,nt),u1t(ne,nt),
     &          u2n(ne,nt),u2t(ne,nt),u2m(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      call frthdeta(u1n,u2n,netamax,nthtamax)
      call frththta(u1t,u2t,netamax,nthtamax)
      call mixedd44(u2m    ,netamax,nthtamax)
c
      return
      end
c
c
c
      subroutine scnddeta(u1n,u2n,ne,nt)
      include 'param.h'
      dimension u1n(ne,nt),u2n(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      ip1=i+1
      im1=i-1
c
      u1n(i,j)=(ub(ip1,j)-ub(im1,j))*dn1o2
      u2n(i,j)=(ub(ip1,j)-2.*ub(i,j)+ub(im1,j))*dn2o2
   10 continue
c
c
      return
      end
c
c
c
      subroutine scndthta(u1t,u2t,ne,nt)
      include 'param.h'
      dimension u1t(ne,nt),u2t(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      do 10 j=1,nthtamax
      jp1=j+1
      jm1=j-1
c
      do 10 i=1,netamax
      u1t(i,j)=(ub(i,jp1)-ub(i,jm1))*dt1o2
      u2t(i,j)=(ub(i,jp1)-2.*ub(i,j)+ub(i,jm1))*dt2o2
   10 continue
c
c
      return
      end
c
c
c
      subroutine scndmixd(u2m,ne,nt)
      include 'param.h'
      dimension u2m(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      do 10 j=1,nthtamax
      jp1=j+1
      jm1=j-1
c
      do 10 i=1,netamax
      ip1=i+1
      im1=i-1
c
      u2m(i,j)=(ub(ip1,jp1)-ub(ip1,jm1)
     &         -ub(im1,jp1)+ub(im1,jm1))*dm2o2
   10 continue
c
c
      return
      end
c
c
c
      subroutine mixedd24(u2m,ne,nt)
      include 'param.h'
      dimension u2m(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      dm=dn1o2*dt1o4
c
      do 10 j=1,nthtamax
      jp1=j+1
      jp2=j+2
      jm1=j-1
      jm2=j-2
c
      do 10 i=1,netamax
      ip1=i+1
      im1=i-1
c
      sp1=ub(ip1,jm2)-ub(ip1,jp2)+8.*(ub(ip1,jp1)-ub(ip1,jm1))
      sm1=ub(im1,jm2)-ub(im1,jp2)+8.*(ub(im1,jp1)-ub(im1,jm1))
      u2m(i,j)=(sp1-sm1)*dm
   10 continue
c
c
      return
      end
c
c
c
      subroutine mixedd44(u2m,ne,nt)
      include 'param.h'
      dimension u2m(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      dm=dn1o2*dt1o4
c
      do 10 j=1,nthtamax
      jp1=j+1
      jp2=j+2
      jm1=j-1
      jm2=j-2
c
      do 10 i=1,netamax
      ip1=i+1
      ip2=i+2
      im1=i-1
      im2=i-2
c
      sm2=ub(im2,jm2)-ub(im2,jp2)+8.*(ub(im2,jp1)-ub(im2,jm1))
      sm1=ub(im1,jm2)-ub(im1,jp2)+8.*(ub(im1,jp1)-ub(im1,jm1))
      sp1=ub(ip1,jm2)-ub(ip1,jp2)+8.*(ub(ip1,jp1)-ub(ip1,jm1))
      sp2=ub(ip2,jm2)-ub(ip2,jp2)+8.*(ub(ip2,jp1)-ub(ip2,jm1))
      u2m(i,j)=(sm2-sp2+8.*(sp1-sm1))*dm2o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine frththta(u1t,u2t,ne,nt)
      include 'param.h'
      dimension u1t(ne,nt),u2t(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      do 10 j=1,nthtamax
      jp1=j+1
      jp2=j+2
      jm1=j-1
      jm2=j-2
c
      do 10 i=1,netamax
      u1t(i,j)=(ub(i,jm2)-ub(i,jp2)+8.*(ub(i,jp1)-ub(i,jm1)))*dt1o4
      u2t(i,j)=(-ub(i,jm2)-ub(i,jp2)+16.*(ub(i,jp1)+ub(i,jm1))
     &          -30.*ub(i,j))*dt2o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine frthdeta(u1n,u2n,ne,nt)
      include 'param.h'
      dimension u1n(ne,nt),u2n(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      im2=i-2
      im1=i-1
      ip1=i+1
      ip2=i+2
c
      u1n(i,j)=(ub(im2,j)-ub(ip2,j)+8.*(ub(ip1,j)-ub(im1,j)))*dn1o4
      u2n(i,j)=(-ub(im2,j)-ub(ip2,j)+16.*(ub(ip1,j)+ub(im1,j))
     &          -30.*ub(i,j))*dn2o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine frthdvbk(ux,u1n,u1t,u2n,u2t,u2m,sym,ne,nt)
      include 'param.h'
      dimension  ux(ne,nt),u1n(ne,nt),u1t(ne,nt),
     &          u2n(ne,nt),u2t(ne,nt),u2m(ne,nt)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
      real bnd(110)
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      jm2=j-2
      jm1=j-1
      jp1=j+1
      jp2=j+2
c
      do 10 i=1,netamax-2
      im2=i-2
      im1=i-1
      ip1=i+1
      ip2=i+2
c
      u1n(i,j)=(ub(im2,j)-ub(ip2,j)+8.*(ub(ip1,j)-ub(im1,j)))*dn1o4
      u1t(i,j)=(ub(i,jm2)-ub(i,jp2)+8.*(ub(i,jp1)-ub(i,jm1)))*dt1o4
c
      u2n(i,j)=(-ub(im2,j)-ub(ip2,j)+16.*(ub(ip1,j)+ub(im1,j))
     &          -30.*ub(i,j))*dn2o4
      u2t(i,j)=(-ub(i,jm2)-ub(i,jp2)+16.*(ub(i,jp1)+ub(i,jm1))
     &          -30.*ub(i,j))*dt2o4
c
      sm2=ub(im2,jm2)-ub(im2,jp2)+8.*(ub(im2,jp1)-ub(im2,jm1))
      sm1=ub(im1,jm2)-ub(im1,jp2)+8.*(ub(im1,jp1)-ub(im1,jm1))
      sp1=ub(ip1,jm2)-ub(ip1,jp2)+8.*(ub(ip1,jp1)-ub(ip1,jm1))
      sp2=ub(ip2,jm2)-ub(ip2,jp2)+8.*(ub(ip2,jp1)-ub(ip2,jm1))
      u2m(i,j)=(sm2-sp2+8.*(sp1-sm1))*dm2o4
   10 continue
c
      do 20 j=1,nthtamax
      jm2=j-2
      jm1=j-1
      jp1=j+1
      jp2=j+2
c
      i=netamax-1
      im1=i-1
      im2=i-2
      im3=i-3
      ip1=i+1
c
      u1t(i,j)=(ub(i,jm2)-ub(i,jp2)+8.*(ub(i,jp1)-ub(i,jm1)))*dt1o4
      u2t(i,j)=(-ub(i,jm2)-ub(i,jp2)+16.*(ub(i,jp1)+ub(i,jm1))
     &          -30.*ub(i,j))*dt2o4
c
      u1n(i,j)=(-ub(im3,j)+6.*ub(im2,j)-18.*ub(im1,j)
     &          +10.*ub(i,j)+3.*ub(ip1,j))*dn1o4
c
      u2n(i,j)=(-ub(im3,j)+4.*ub(im2,j)+6.*ub(im1,j)-20.*ub(i,j)
     &          +11.*ub(ip1,j))*dn2o4
c
      sm3=ub(im3,jm2)-8.*ub(im3,jm1)+8.*ub(im3,jp1)-ub(im3,jp2)
      sm2=ub(im2,jm2)-8.*ub(im2,jm1)+8.*ub(im2,jp1)-ub(im2,jp2)
      sm1=ub(im1,jm2)-8.*ub(im1,jm1)+8.*ub(im1,jp1)-ub(im1,jp2)
      s  =ub(i  ,jm2)-8.*ub(i  ,jm1)+8.*ub(i  ,jp1)-ub(i  ,jp2)
      sp1=ub(ip1,jm2)-8.*ub(ip1,jm1)+8.*ub(ip1,jp1)-ub(ip1,jp2)
      u2m(i,j)=(-sm3+6.*sm2-18.*sm1+10.*s+3.*sp1)*dm2o4
c
      i=netamax
      im1=i-1
      im2=i-2
      im3=i-3
      im4=i-4
c
      u1t(i,j)=(ub(i,jm2)-ub(i,jp2)+8.*(ub(i,jp1)-ub(i,jm1)))*dt1o4
      u2t(i,j)=(-ub(i,jm2)-ub(i,jp2)+16.*(ub(i,jp1)+ub(i,jm1))
     &          -30.*ub(i,j))*dt2o4
c
      u1n(i,j)=(3.*ub(im4,j)-16.*ub(im3,j)+36.*ub(im2,j)
     &         -48*ub(im1,j)+25.*ub(i,j))*dn1o4
      u2n(i,j)=(11.*ub(im4,j)-56.*ub(im3,j)+114.*ub(im2,j)
     &          -104.*ub(im1,j)+35.*ub(i,j))*dn2o4
c
      sm4=ub(im4,jm2)-8.*ub(im4,jm1)+8.*ub(im4,jp1)-ub(im4,jp2)
      sm3=ub(im3,jm2)-8.*ub(im3,jm1)+8.*ub(im3,jp1)-ub(im3,jp2)
      sm2=ub(im2,jm2)-8.*ub(im2,jm1)+8.*ub(im2,jp1)-ub(im2,jp2)
      sm1=ub(im1,jm2)-8.*ub(im1,jm1)+8.*ub(im1,jp1)-ub(im1,jp2)
      s  =ub(i  ,jm2)-8.*ub(i  ,jm1)+8.*ub(i  ,jp1)-ub(i  ,jp2)
      u2m(i,j)=(3.*sm4-16.*sm3+36.*sm2-48.*sm1+25.*s)*dm2o4
   20 continue
c
c
      return
      end
c
c
c
      subroutine frstdrv2(ux,u1n,u1t,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt),u1n(ne,nt),u1t(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      u1n(i,j)=(ub(i+1,j)-ub(i-1,j))*dn1o2
      u1t(i,j)=(ub(i,j+1)-ub(i,j-1))*dt1o2
   10 continue
c
c
      return
      end
c
c
c
      subroutine frstdv24(ux,u1n,u1t,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt),u1n(ne,nt),u1t(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      u1n(i,j)=(ub(i+1,j)-ub(i-1,j))*dn1o2
      u1t(i,j)=(ub(i,j-2)-ub(i,j+2)+8.*(ub(i,j+1)-ub(i,j-1)))*dt1o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine frstdrv4(ux,u1n,u1t,sym,bnd,ne,nt)
      include 'param.h'
      dimension ux(ne,nt),u1n(ne,nt),u1t(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
c
c
      call boundry2(ux,bnd,sym,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      u1n(i,j)=(ub(i-2,j)-ub(i+2,j)+8.*(ub(i+1,j)-ub(i-1,j)))*dn1o4
      u1t(i,j)=(ub(i,j-2)-ub(i,j+2)+8.*(ub(i,j+1)-ub(i,j-1)))*dt1o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine boundry2(ux,bnd,sym,ne,nt)
      include 'param.h'
      dimension ux(ne,nt),bnd(*)
c
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
      do 5 j=1,nthtamax
      do 5 i=1,netamax
    5 ub(i,j)=ux(i,j)
c
      ntm2=nthtamax-2
      ntp2=nthtamax+2
c
      do 10 i=1,netamax
      ub(i,0)=sym*ux(i,2)
      ub(i,-1)=sym*ux(i,3)
      ub(i,ntp1)=sym*ux(i,ntm1)
      ub(i,ntp2)=sym*ux(i,ntm2)
   10 continue
c
      do 20 j=1,nthtamax
      ub(0,j)=sym*ux(2,j)
      ub(-1,j)=sym*ux(3,j)
   20 continue
c
      ub(0,0)=ux(2,2)
      ub(-1,0)=ux(3,2)
      ub(0,-1)=ux(2,3)
      ub(-1,-1)=ux(3,3)
c
      ub(0,ntp1)=ux(2,ntm1)
      ub(0,ntp2)=ux(2,ntm2)
      ub(-1,ntp1)=ux(3,ntm1)
      ub(-1,ntp2)=ux(3,ntm2)
c
      nep2=netamax+2
      do 30 j=1,nthtamax
      ub(nep1,j)=bnd(j)
      ub(nep2,j)=bnd(j)
   30 continue
c
      ub(nep1, 0)=sym*bnd(2)
      ub(nep1,-1)=sym*bnd(3)
      ub(nep2, 0)=sym*bnd(2)
      ub(nep2,-1)=sym*bnd(3)
c
      ntp2=nthtamax+2
      ntm2=nthtamax-2
      ub(nep1,ntp1)=sym*bnd(ntm1)
      ub(nep1,ntp2)=sym*bnd(ntm2)
      ub(nep2,ntp1)=sym*bnd(ntm1)
      ub(nep2,ntp2)=sym*bnd(ntm2)
c
c
      return
      end
c
c
c
      subroutine hwkgmass(hmass)
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /spinblck/ rho(id,jd),amu(id,jd),sigma(id,jd),
     &                  alambda(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /func1blk/ h(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
      character*8 answer
c
      call spincoef
c
      ns=netamax
c
c        compute surface area of coordinate sphere
      do 10 j=1,nthtamax
   10 h(j)=sqrt(b(ns,j)*d(ns,j))*sn(j)/psim4(ns,j)
c
      call simpson(h,dtheta,nthtamax,ans)
      area=ans*8.*hfpi
c
c        compute surface integral over mu*rho
      do 20 j=1,nthtamax
      ds=sqrt(b(ns,j)*d(ns,j))*sn(j)/psim4(ns,j)
      h(j)=rho(ns,j)*amu(ns,j)*ds
   20 continue
c
      call simpson(h,dtheta,nthtamax,ans)
      rhomuint=ans*8.*hfpi
c
      hmass=((8.*hfpi)**(-1.5))*sqrt(area)*(4.*hfpi-rhomuint)
c
      go to 99
      write(5,*)
      write(5,*) 'areas?'
      read(5,101) answer
  101 format(8a)
      if (answer.eq.'y') then
c        compute surface area of coordinate sphere
      do 89 i=1,netamax/2
      do 88 j=1,nthtamax
   88 h(j)=sqrt(b(i,j)*d(i,j))*sn(j)/psim4(i,j)
      call simpson(h,dtheta,nthtamax,ans)
      area=ans*8.*hfpi
      write(5,*)eta(i),area
   89 continue
      endif
c
c
   99 return
      end
c
c
c
      subroutine eventest(ans)
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /func1blk/ h(500)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real am(400)
c
c
      imin=1
c
c        compute surface area of coordinate sphere
      do 20 i=1,netamax
      do 10 j=1,nthtamax
   10 h(j)=sqrt(b(i,j)*d(i,j))*sn(j)/psim4(i,j)
      call simpson(h,dtheta,nthtamax,ans1)
      area=ans1*8.*hfpi
      am(i)=sqrt(area/(32.*hfpi))/adm
      if (am(i).lt.1.0) imin=i
   20 continue
c
      ax=(am(imin+1)-1.)/(am(imin+1)-am(imin))
      bx=(1.-am(imin))/(am(imin+1)-am(imin))
      ans=ax*eta(imin)+bx*eta(imin+1)
c
c
      return
      end
c
c
c
      subroutine admmass(ans)
      include 'param.h'
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /func1blk/ h(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        this subroutine computes the ADM mass of the initial
c        slice.
c
      ns=netamax
c
      do 10 j=1,nthtamax
      h(j)=0.5*psi(ns,j)*(0.5*phi1n(ns,j)-1.)*sn(j)
   10 continue
c
      call simpson(h,dtheta,nthtamax,ans)
c
      adm=-2.*sqrthfm*exp(0.5*eta(ns))*ans
      ans=adm
c
c
      return
      end
c
c
c
      subroutine brillmas(brillmsc)
      include 'param.h'
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /func1blk/ h(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      real angleint(500)
c
c
c        this subroutine computes the brill mass of the initial
c        slice.
c
      ns=netamax
c
c        compute the volume integral
      do 20 i=1,ns
      do 10 j=1,nthtamax
   10 h(j)=((phi1n(i,j)-2.0)**2.+phi1t(i,j)**2.)*sn(j)
c
      call simpson(h,dtheta,nthtamax,ans)
      angleint(i)=2.*ans
   20 continue
c
      do 30 i=1,ns
      e=eta(i)
      h(i)=exp(e)*angleint(i)
   30 continue
c
      call simpson(h,deta,ns,volint)
c
c        compute q term
      call computeq
      do 50 j=1,nthtamax
   50 h(j)=f(1,j)*sn(j)
c
      call simpson(h,dtheta,nthtamax,ans)
      qterm=2.*ans
c
      brillms=scaleprm*(0.5+qterm/8.+volint/32.)
c
c        compute brill mass with 1/R correction
      bigr=0.5*scaleprm*exp(eta(netamax))
      brillmsc=brillms/(1.-0.5*brillms/bigr)
c
c
      return
      end
c
c
c
      subroutine simpson(u,dx,nxmax,ans)
      dimension u(*)
c
c
      ans=0.
      do 10 k=5,nxmax-4
   10 ans=ans+u(k)
c
      ans=ans+(u(1)+u(nxmax  ))*(17./48.)
     &       +(u(2)+u(nxmax-1))*(59./48.)
     &       +(u(3)+u(nxmax-2))*(43./48.)
     &       +(u(4)+u(nxmax-3))*(49./48.)
c
      ans=ans*dx
c
c
      return
      end
c
c
c
      subroutine abrahams
      include 'param.h'
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
c
      real radius(3)
c
c
      iw=nint(jtime/fiwave)
      nfwave=nint(fiwave*(1.+iw))
c
      thetime=dt*jtime/adm
c
      radius(1)=15.*scaleprm
      radius(2)=30.*scaleprm
      radius(3)=60.*scaleprm
c
      do 200 iradius=1,3
      r=radius(iradius)
      etad=dlog(2.*r/scaleprm)
      i=nint(etad/deta)+1
c
      call amplitud(a022,a202,a222,a242,a022r,a202r,a222r,a242r,i)
c
      ifile=50+iradius
      write(ifile,*) thetime
      write(ifile,*) a022
      write(ifile,*) a202
      write(ifile,*) a222
      write(ifile,*) a242
      write(ifile,*) a022r
      write(ifile,*) a202r
      write(ifile,*) a222r
      write(ifile,*) a242r
  200 continue
c
c
      return
      end
c
c
c
      subroutine zerilli
      include 'param.h'
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
c The following common blocks and arrays added by Ed Seidel 12/7/90
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
c
      real r2int(110),a2int(110),h2int2(110),gint2(110),kint2(110),
     1 h2int4(110),gint4(110),kint4(110),gint2r(110),kint2r(110),
     2 gint4r(110),kint4r(110),h1int2(110),h1int4(110)
      real radius(4)
      real k2,k4,k2r,k4r,g2r,g4r
c
c  This subroutine computes the l=2 and l=4 Zerilli functions, as
c  well as Andrew's flat space gauge invariant function.
c  l=2 has been tested against Andrew's expressions, but not l=4.
c  It is written to be completely general in that it should
c  work with any three-metric (doesn't assume it to be diagonal)
c ***WARNING:  IT HAS NEVER BEEN TESTED WITH NON-DIAGONAL 3-METRIC***
c  It has been completely rewritten by Ed Seidel in December, 1990
c  from the earlier version from David Bernstein.
c
c
      iw=nint(jtime/fiwave)
      nfwave=nint(fiwave*(1.+iw))
c
      thetime=dt*jtime/adm
c
      pi = 2.*acos(0.)
c
      radius(1)=15.*scaleprm
      radius(2)=30.*scaleprm
      radius(3)=45.*scaleprm
      radius(4)=60.*scaleprm
c
      do 10 iradius=1,4
      r=radius(iradius)
      etad=dlog(2.*r/scaleprm)
      i=int(etad/deta)+1
c
      e=eta(i)
      ri=0.5*scaleprm*exp(e)
c test against other measures of radius
c this was the isotropic radius
      dedr = exp(e)/(exp(2.*e)-1.)
      do 11 j=1,nthtamax
c Note: I did not include the factor sin^2 in gphph to avoid
c having the code blow up on me at theta = 0
      grr = a(i,j)/psim4(i,j)
      gthth = b(i,j)/psim4(i,j)
      grth = c(i,j)/psim4(i,j)
      gphph = d(i,j)/psim4(i,j)
      gththr = dedr/psim4(i,j)*(b1n(i,j)+b(i,j)*phi1n(i,j))
      gphphr = dedr/psim4(i,j)*(d1n(i,j)+d(i,j)*phi1n(i,j))
c
c compute (l = 0) spherical background metric quantities
c
      a2int(j) = sn(j)*grr
      r2int(j) = sn(j)*gthth
c
c compute the l=2 expression
c
      h1int2(j) = acs(j)*sn(j)*sn(j)*grth
      h2int2(j) = sn(j)*grr*(3.*acs(j)**2-1.)
      gint2(j) =  sn(j)**3*(gthth-gphph)
      kint2(j) = sn(j)*(gphph*(4.-9.*sn(j)**2) + gthth*
     1            (4.-3.*sn(j)**2))
      gint2r(j) =  sn(j)**3*(gththr-gphphr)
      kint2r(j) = sn(j)*(gphphr*(4.-9.*sn(j)**2) + gththr*
     1            (4.-3.*sn(j)**2))
c
c compute the l=4 expressions
c
      h1int4(j) = (3.-7.*acs(j)**2)*sn(j)*acs(j)*grth
      h2int4(j) = sn(j)*grr*3./(16.*sqrt(2.*hfpi))*
     1  (3.-30.*acs(j)**2+35.*acs(j)**4)
      gint4(j) = (7.*acs(j)**2-1.)*(gthth-gphph)*sn(j)**3
      kint4(j) = sn(j)*(gthth*(35.*sn(j)**4-60.*sn(j)**2+24.) +
     1  gphph*(175.*sn(j)**4-180.*sn(j)**2+24.))
      gint4r(j) = (7.*acs(j)**2-1.)*(gththr-gphphr)*sn(j)**3
      kint4r(j) = sn(j)*(gththr*(35.*sn(j)**4-60.*sn(j)**2+24.) +
     1  gphphr*(175.*sn(j)**4-180.*sn(j)**2+24.))
11    continue
c
c compute the projected spherical areal radius
c and grr metric component if desired
c
      call simpson(r2int,dtheta,nthtamax,ans)
      r = sqrt(ans)
c      call simpson(a2int,dtheta,nthtamax,ans)
c      ax = sqrt(ans)
c  Testing of small amplitude Brill wave (amp=0.1) shows
c  Schwarzschild is a good approximation for ax, so use it.
c
c here are some background quantities
      s = 1.-2.*scaleprm/r
      ax = sqrt(r**2/(1.-2.*scaleprm/r))
c
c Compute l = 2 expressions
c
      call simpson(h1int2,dtheta,nthtamax,ans)
      h12 = -sqrt(10.*hfpi)*ans
      call simpson(h2int2,dtheta,nthtamax,ans)
      h22 = sqrt(10.*hfpi)/(ax**2)*ans
      call simpson(gint2,dtheta,nthtamax,ans)
      g2 = sqrt(10.*hfpi)/(4.*r**2)*ans
      call simpson(kint2,dtheta,nthtamax,ans)
      k2 = sqrt(10.*hfpi)/(4.*r**2)*ans
      call simpson(gint2r,dtheta,nthtamax,ans)
      g2r = sqrt(10.*hfpi)/(4.*r**2)*ans-2./r*g2
      call simpson(kint2r,dtheta,nthtamax,ans)
      k2r = sqrt(10.*hfpi)/(4.*r**2)*ans-2./r*k2
c
c Now put the pieces together for the l = 2 Zerilli function
c
      xk12 = k2 + s*r*g2r - 2.*s*h12/r
      xk22 = h22/(2.*s)-.5/s*(r*k2r+k2*(r-3.*scaleprm)/(r-2.*scaleprm))
      q12 = 4.*r*s*s*xk22 + 6.*r*xk12
      xlam2 = 1. + 3.*scaleprm/(2.*r)
      psi2 = q12/xlam2
c
c Compute l = 4 expressions
c
      call simpson(h1int4,dtheta,nthtamax,ans)
      h14 = 2.5*sqrt(2.*hfpi)*ans
      call simpson(h2int4,dtheta,nthtamax,ans)
      h24 = 8.*hfpi/ax**2*ans
      call simpson(gint4,dtheta,nthtamax,ans)
      g4 = sqrt(2.*hfpi)/(8.*r**2)*ans
      call simpson(kint4,dtheta,nthtamax,ans)
      k4 = sqrt(2.*hfpi)/(8.*r**2)*ans
      call simpson(gint4r,dtheta,nthtamax,ans)
      g4r = sqrt(2.*hfpi)/(8.*r**2)*ans-2./r*g4
      call simpson(kint4r,dtheta,nthtamax,ans)
      k4r = sqrt(2.*hfpi)/(8.*r**2)*ans-2./r*k4
c
c Now put the pieces together for the l = 4 Zerilli function
c ***CURRENTLY HAVE NOT ATTEMPTED TO NORMALIZE IT FOR COMPARISON
C TO THE L = 2 ZERILLI FUNCTION***
c
      xk14 = k4 + s*r*g4r - 2.*s*h14/r
      xk24 = h24/(2.*s)-.5/s*(r*k4r+k4*(r-3.*scaleprm)/(r-2.*scaleprm))
      q14 = 4.*r*s*s*xk24 + 20.*r*xk14
      xlam4 = 1. + scaleprm/(3.*r)
      psi4 = q14/xlam4
c
c        normalize for energy radiated
      psi2=psi2/sqrt(12.)
      psi4=psi4/sqrt(180.)
c
c compute Andrew's l = 2 flat space extraction
c
      call amplitud(a022,a202,a222,a242,a022r,a202r,a222r,a242r,i)
      gia = sqrt(8.*hfpi/5.)*r*(6.*a022-9.*a222+5.*a242
     1      -2.*a022r+a202r+a222r+a242r)
      rstar = r + 2.*scaleprm*log(r/(2.*scaleprm)-1.)
c write to zerilli detectors at fixed eta
      ifile=50+iradius
      write(ifile,*) thetime,r/adm,rstar,psi2,psi4
   10 continue
c
c
      return
      end
c
c
c
      subroutine zermovie
      include 'param.h'
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
c The following common blocks and arrays added by Ed Seidel 12/7/90
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
c
      real r2int(110),a2int(110),h2int2(110),gint2(110),kint2(110),
     1 h2int4(110),gint4(110),kint4(110),gint2r(110),kint2r(110),
     2 gint4r(110),kint4r(110),h1int2(110),h1int4(110)
      real k2,k4,k2r,k4r,g2r,g4r
c
c  This subroutine computes the l=2 and l=4 Zerilli functions, as
c  well as Andrew's flat space gauge invariant function.
c  l=2 has been tested against Andrew's expressions, but not l=4.
c  It is written to be completely general in that it should
c  work with any three-metric (doesn't assume it to be diagonal)
c ***WARNING:  IT HAS NEVER BEEN TESTED WITH NON-DIAGONAL 3-METRIC***
c  It has been completely rewritten by Ed Seidel in December, 1990
c  from the earlier version from David Bernstein.
c
      iw=nint(jtime/fimwave)
      nfmwave=nint(fimwave*(1.+iw))
c
      thetime=dt*jtime/adm
      write(61,*) thetime
      write(62,*) thetime
      pi = 2.*acos(0.)
c
c        points outside horizon
      do 5 i=1,netamax
      do 4 j=1,nthtamax
      gthth = b(i,j)/psim4(i,j)
      r2int(j) = sn(j)*gthth
    4 continue
      call simpson(r2int,dtheta,nthtamax,ans)
      r = sqrt(abs(ans))
      if(r.gt.2.*scaleprm) go to 6
    5 continue
    6 imin=i+1
      write(61,*) netamax-imin+1
      write(62,*) netamax-imin+1
c
      do 10 i=imin,netamax
c
      e=eta(i)
      ri=0.5*scaleprm*exp(e)
c test against other measures of radius
c this was the isotropic radius
      dedr = exp(e)/(exp(2.*e)-1.)
      do 11 j=1,nthtamax
c Note: I did not include the factor sin^2 in gphph to avoid
c having the code blow up on me at theta = 0
      grr = a(i,j)/psim4(i,j)
      gthth = b(i,j)/psim4(i,j)
      grth = c(i,j)/psim4(i,j)
      gphph = d(i,j)/psim4(i,j)
      gththr = dedr/psim4(i,j)*(b1n(i,j)+b(i,j)*phi1n(i,j))
      gphphr = dedr/psim4(i,j)*(d1n(i,j)+d(i,j)*phi1n(i,j))
c
c compute (l = 0) spherical background metric quantities
c
      a2int(j) = sn(j)*grr
      r2int(j) = sn(j)*gthth
c
c compute the l=2 expression
c
      h1int2(j) = acs(j)*sn(j)*sn(j)*grth
      h2int2(j) = sn(j)*grr*(3.*acs(j)**2-1.)
      gint2(j) =  sn(j)**3*(gthth-gphph)
      kint2(j) = sn(j)*(gphph*(4.-9.*sn(j)**2) + gthth*
     1            (4.-3.*sn(j)**2))
      gint2r(j) =  sn(j)**3*(gththr-gphphr)
      kint2r(j) = sn(j)*(gphphr*(4.-9.*sn(j)**2) + gththr*
     1            (4.-3.*sn(j)**2))
c
c compute the l=4 expressions
c
      h1int4(j) = (3.-7.*acs(j)**2)*sn(j)*acs(j)*grth
      h2int4(j) = sn(j)*grr*3./(16.*sqrt(2.*hfpi))*
     1  (3.-30.*acs(j)**2+35.*acs(j)**4)
      gint4(j) = (7.*acs(j)**2-1.)*(gthth-gphph)*sn(j)**3
      kint4(j) = sn(j)*(gthth*(35.*sn(j)**4-60.*sn(j)**2+24.) +
     1  gphph*(175.*sn(j)**4-180.*sn(j)**2+24.))
      gint4r(j) = (7.*acs(j)**2-1.)*(gththr-gphphr)*sn(j)**3
      kint4r(j) = sn(j)*(gththr*(35.*sn(j)**4-60.*sn(j)**2+24.) +
     1  gphphr*(175.*sn(j)**4-180.*sn(j)**2+24.))
11    continue
c
c compute the projected spherical areal radius
c and grr metric component if desired
c
      call simpson(r2int,dtheta,nthtamax,ans)
      r = sqrt(abs(ans))
c      call simpson(a2int,dtheta,nthtamax,ans)
c      ax = sqrt(abs(ans))
c  Testing of small amplitude Brill wave (amp=0.1) shows
c  Schwarzschild is a good approximation for ax, so use it.
c
c here are some background quantities
      s = 1.-2.*scaleprm/r
      ax = sqrt(abs(r**2/(1.-2.*scaleprm/r)))
c
c Compute l = 2 expressions
c
      call simpson(h1int2,dtheta,nthtamax,ans)
      h12 = -sqrt(10.*hfpi)*ans
      call simpson(h2int2,dtheta,nthtamax,ans)
      h22 = sqrt(10.*hfpi)/(ax**2)*ans
      call simpson(gint2,dtheta,nthtamax,ans)
      g2 = sqrt(10.*hfpi)/(4.*r**2)*ans
      call simpson(kint2,dtheta,nthtamax,ans)
      k2 = sqrt(10.*hfpi)/(4.*r**2)*ans
      call simpson(gint2r,dtheta,nthtamax,ans)
      g2r = sqrt(10.*hfpi)/(4.*r**2)*ans-2./r*g2
      call simpson(kint2r,dtheta,nthtamax,ans)
      k2r = sqrt(10.*hfpi)/(4.*r**2)*ans-2./r*k2
c
c Now put the pieces together for the l = 2 Zerilli function
c
      xk12 = k2 + s*r*g2r - 2.*s*h12/r
      xk22 = h22/(2.*s)-.5/s*(r*k2r+k2*(r-3.*scaleprm)/(r-2.*scaleprm))
      q12 = 4.*r*s*s*xk22 + 6.*r*xk12
      xlam2 = 1. + 3.*scaleprm/(2.*r)
      psi2 = q12/xlam2
c
c Compute l = 4 expressions
c
      call simpson(h1int4,dtheta,nthtamax,ans)
      h14 = 2.5*sqrt(2.*hfpi)*ans
      call simpson(h2int4,dtheta,nthtamax,ans)
      h24 = 8.*hfpi/ax**2*ans
      call simpson(gint4,dtheta,nthtamax,ans)
      g4 = sqrt(2.*hfpi)/(8.*r**2)*ans
      call simpson(kint4,dtheta,nthtamax,ans)
      k4 = sqrt(2.*hfpi)/(8.*r**2)*ans
      call simpson(gint4r,dtheta,nthtamax,ans)
      g4r = sqrt(2.*hfpi)/(8.*r**2)*ans-2./r*g4
      call simpson(kint4r,dtheta,nthtamax,ans)
      k4r = sqrt(2.*hfpi)/(8.*r**2)*ans-2./r*k4
c
c Now put the pieces together for the l = 4 Zerilli function
c ***CURRENTLY HAVE NOT ATTEMPTED TO NORMALIZE IT FOR COMPARISON
C TO THE L = 2 ZERILLI FUNCTION***
c
      xk14 = k4 + s*r*g4r - 2.*s*h14/r
      xk24 = h24/(2.*s)-.5/s*(r*k4r+k4*(r-3.*scaleprm)/(r-2.*scaleprm))
      q14 = 4.*r*s*s*xk24 + 20.*r*xk14
      xlam4 = 1. + scaleprm/(3.*r)
      psi4 = q14/xlam4
c
c        normalize for energy radiated
      psi2=psi2/sqrt(12.)
      psi4=psi4/sqrt(180.)
c
c compute Andrew's l = 2 flat space extraction
c
      call amplitud(a022,a202,a222,a242,a022r,a202r,a222r,a242r,i)
      gia = sqrt(8.*hfpi/5.)*r*(6.*a022-9.*a222+5.*a242
     1      -2.*a022r+a202r+a222r+a242r)
      rstar = r + 2.*scaleprm*dlog(abs(r/(2.*scaleprm)-1.))
      write(61,*)r/adm,rstar,psi2,psi4,gia,eta(i)
      write(62,*)r/adm,h22,g2,k2,h24,g4,k4
   10 continue
c
c
      return
      end
c
c
c
      subroutine amplitud(a022 ,a202 ,a222 ,a242 ,
     &                    a022r,a202r,a222r,a242r,i)
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /func1blk/ h(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
c
      real  q1(500), q2(500), q3(500), q4(500),
     &      q5(500), q6(500), q7(500), q8(500)
      real radius(3)
c
c
      e=eta(i)
      r=0.5*scaleprm*exp(e)
c
      do 10 j=1,nthtamax
      fac=1./(psim4(i,j)*r*r)
c
      cs2=acs(j)*acs(j)
      trig1=3.*cs2-1.
      trig2=2.-3.*cs2
      trig3=3.-7.*cs2
      trig4=1.-5.*cs2
c
c        radial derivatives multiplied by r
      ap1n=a1n(i,j)+a(i,j)*phi1n(i,j)-2.*a(i,j)
      bp1n=b1n(i,j)+b(i,j)*phi1n(i,j)-2.*b(i,j)
      dp1n=d1n(i,j)+d(i,j)*phi1n(i,j)-2.*d(i,j)
c
      q1(j)=fac*(a(i,j)+b(i,j)+d(i,j))*trig1*sn(j)
      q2(j)=fac*(a(i,j)*trig1+b(i,j)*trig2-d(i,j))*sn(j)
      q3(j)=fac*(a(i,j)*trig1-b(i,j)+d(i,j)*trig2)*sn(j)
      q4(j)=fac*(4.*trig1*a(i,j)+trig3*b(i,j)+trig4*d(i,j))*sn(j)
c
      q5(j)=fac*(ap1n+bp1n+dp1n)*trig1*sn(j)
      q6(j)=fac*(ap1n*trig1+bp1n*trig2-dp1n)*sn(j)
      q7(j)=fac*(ap1n*trig1-bp1n+dp1n*trig2)*sn(j)
      q8(j)=fac*(4.*trig1*ap1n+trig3*bp1n+trig4*dp1n)*sn(j)
   10 continue
c
      call simpson(q1,dtheta,nthtamax,ans)
      a022=ans*5./6.
c
      call simpson(q2,dtheta,nthtamax,ans)
      a202=ans/3.
c
      call simpson(q3,dtheta,nthtamax,ans)
      a222=ans*10./21.
c
      call simpson(q4,dtheta,nthtamax,ans)
      a242=ans*3./14.
c
      call simpson(q5,dtheta,nthtamax,ans)
      a022r=ans*5./6.
c
      call simpson(q6,dtheta,nthtamax,ans)
      a202r=ans/3.
c
      call simpson(q7,dtheta,nthtamax,ans)
      a222r=ans*10./21.
c
      call simpson(q8,dtheta,nthtamax,ans)
      a242r=ans*3./14.
c
c
      return
      end
c
c
c
      subroutine rpout
      include 'param.h'
      common /rpatblck/ nradpat,nfirp,firp
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /elmagblk/ e11(id,jd),e12(id,jd),e22(id,jd),
     &                  e33(id,jd),b13(id,jd),b23(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
      iw=nint(jtime/firp)
      nfirp=nint(firp*(1.+iw))
c
      call elmagnoc
      call belrobvc(f1,f2,f,netamax,nthtamax)
c
      r=15.*scaleprm
      etad=dlog(2.*r/scaleprm)
      i1=int(etad/deta)+1
c
      r=30.*scaleprm
      etad=dlog(2.*r/scaleprm)
      i2=int(etad/deta)+1
c
      r=45.*scaleprm
      etad=dlog(2.*r/scaleprm)
      i3=int(etad/deta)+1
c
      r=60.*scaleprm
      etad=dlog(2.*r/scaleprm)
      i4=int(etad/deta)+1
c
      do 10 j=1,nthtamax
      write(71,101) f1(i1,j),f1(i2,j),f1(i3,j),f1(i4,j)
   10 continue
c
  101 format(4e16.5)
c
c
      return
      end
c
c
c
      subroutine spincoef
      include 'param.h'
      common /spinblck/ rho(id,jd),amu(id,jd),sigma(id,jd),
     &                  alambda(id,jd)
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
      sqrt8m1=0.5/sqrt2
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      p2=psi(i,j)*psi(i,j)
c
      r1=0.5*phi1n(i,j)/(sqrt(2.*a(i,j))*p2)
      r2=d1n(i,j)/(4.*sqrt(2.*a(i,j))*d(i,j)*p2)
      r3=b1n(i,j)/(4.*sqrt(2.*a(i,j))*b(i,j)*p2)
      r4=hd(i,j)/(2.*sqrt(2.)*d(i,j))
      r5=hb(i,j)/(2.*sqrt(2.)*b(i,j))
c
      rho(i,j)=-r1-r2-r3+r4+r5
c
      amu(i,j)=-r1-r2-r3-r4-r5
c
      sigma(i,j)=r2-r3-r4+r5
c
      alambda(i,j)=r2-r3+r4-r5
   10 continue
c
c
      return
      end
c
c
c
      subroutine elmagnoc
      include 'param.h'
      common /elmagblk/ e11(id,jd),e12(id,jd),e22(id,jd),
     &                  e33(id,jd),b13(id,jd),b23(id,jd)
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /chrisblk/ g111(id,jd),g112(id,jd),g122(id,jd),
     &                  g133(id,jd),g211(id,jd),g212(id,jd),
     &                  g222(id,jd),g233(id,jd),g313(id,jd),
     &                  g323(id,jd)
      common /ricciblk/ ricci11(id,jd),ricci12(id,jd),
     &                  ricci22(id,jd),ricci33(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        construct electric and magnetic components of the
c        Riemann tensor
c        (a factor of psi^4 is taken out of all the electric
c         components and sin(theta)^2 is taken out of e33.
c         a factor of sin(theta)*psi^2 is taken out of all the
c         magnetic components)
c
c        the 3-metric is assumed to be diagonal in this subroutine
c
c        off axis
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      e11(i,j)=-ricci11(i,j)*psim4(i,j)
     &         -ha(i,j)*(hb(i,j)/b(i,j)+hd(i,j)/d(i,j))
     &         +hc(i,j)*hc(i,j)/b(i,j)
c
      e12(i,j)=-ricci12(i,j)*psim4(i,j)-hc(i,j)*hd(i,j)/d(i,j)
c
      e22(i,j)=-ricci22(i,j)*psim4(i,j)
     &         -hb(i,j)*(ha(i,j)/a(i,j)+hd(i,j)/d(i,j))
     &         +hc(i,j)*hc(i,j)/a(i,j)
c
      e33(i,j)=-ricci33(i,j)*psim4(i,j)
     &         -hd(i,j)*(ha(i,j)/a(i,j)+hb(i,j)/b(i,j))
c
      fac=1./sqrt(a(i,j)*b(i,j)*d(i,j))
c
      x1=d(i,j)*(hc1n(i,j)+hc(i,j)*phi1n(i,j)
     &          -ha1t(i,j)-ha(i,j)*phi1t(i,j))
      x2=a(i,j)*(hd1t(i,j)+hd(i,j)*phi1t(i,j)+2.*hd(i,j)*ct(j))
c
      x3=d(i,j)*(-ha(i,j)*g112(i,j)+hb(i,j)*g211(i,j)
     &           +hc(i,j)*(g111(i,j)-g212(i,j)))
      x4=a(i,j)*(-hc(i,j)*g133(i,j)-hb(i,j)*g233(i,j)
     &           +hd(i,j)*g323(i,j))
c
      b13(i,j)=-0.5*fac*(x1+x2-x3-x4)
c
      r1=d(i,j)*(-hc1t(i,j)-hc(i,j)*phi1t(i,j)
     &           +hb1n(i,j)+hb(i,j)*phi1n(i,j))
      r2=-b(i,j)*(hd1n(i,j)+hd(i,j)*phi1n(i,j))
c
      r3=d(i,j)*(-ha(i,j)*g122(i,j)+hb(i,j)*g212(i,j)
     &           +hc(i,j)*(g112(i,j)-g222(i,j)))
      r4=b(i,j)*( ha(i,j)*g133(i,j)+hc(i,j)*g233(i,j)
     &           -hd(i,j)*g313(i,j))
c
      b23(i,j)=-0.5*fac*(r1+r2-r3-r4)
   10 continue
c
c        and on axis
      do 20 i=1,netamax
c
      e11(i,1)=-ricci11(i,1)*psim4(i,1)-2.*ha(i,1)*hb(i,1)/b(i,1)
c
      e12(i,1)=-ricci12(i,1)*psim4(i,1)
c
      e22(i,1)=-ricci22(i,1)*psim4(i,1)
     &         -hb(i,1)*(ha(i,1)/a(i,1)+hb(i,1)/b(i,1))
c
      e33(i,1)=e22(i,1)
c
      b13(i,1)=0.
      b23(i,1)=0.
   20 continue
c
c
      return
      end
c
c
c
      subroutine paramtrs
      include 'param.h'
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
      common /rpatblck/ nradpat,nfirp,firp
      common /equblock/ nequ,nfequ,fiequ
      common /horiznbk/ nah,nfah,fiah
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        set the scale parameter m
      scaleprm=2.0
c
c        set the maximum value of eta and the number of radial
c        grid points and calculate deta
      etamax=6.
      netamax=id
      deta=etamax/(netamax-1)
c
c        the number of angular grid points is set to make deta and
c        dtheta as close as possible
      hfpi=dacos(0.)
      nthtamax=iidnnt(hfpi/deta)+1
      if (nthtamax.ne.jd) then
      write(5,*) 'no match between jd and nthtamax calculation'
      stop
      endif
      dtheta=hfpi/(nthtamax-1)
      xi=deta/dtheta
c
      do 10 i=1,netamax
      eta(i)=(i-1)*deta
   10 continue
c
      theta(1)=0.
      sn(1)=0.
      acs(1)=1.
      do 20 j=2,nthtamax
      t=(j-1)*dtheta
      theta(j)=t
      sn(j)=sin(t)
      acs(j)=cos(t)
      ct(j)=acs(j)/sn(j)
   20 continue
c
c        compute total number of grid points
      ngridpts=netamax*nthtamax
c
      nem1=netamax-1
      nep1=netamax+1
      ntm1=nthtamax-1
      ntp1=nthtamax+1
c
c        a number which occurs often is the square root of M/2
      sqrthfm=sqrt(0.5*scaleprm)
c
c        compute differences used in the finite difference schemes
      dn1o2=0.5/deta
      dn2o2=1./(deta*deta)
      dt1o2=0.5/dtheta
      dt2o2=1./(dtheta*dtheta)
      dm2o2=0.25/(dtheta*deta)
c
      dn1o4=1./(12.*deta)
      dn2o4=1./(12.*deta*deta)
      dt1o4=1./(12.*dtheta)
      dt2o4=1./(12.*dtheta*dtheta)
      dm2o4=1./(144.*deta*dtheta)
c
c        compute or read in conformal factor
      call ivp
      call admmass(adm)
      write(5,*)'ADM mass (/m) at grid boundary is  ',adm/scaleprm
c
      if (qrun.eq.'y') then
      if (waytogo.eq.'interact') then
      write(5,*)
      write(5,*)'time of run (/ADM)'
      read(5,*) time
      endif
      else
      time=1000.
      endif
c
c        set time step and calculate the total number of steps
      dt=1.0*deta
      ntime=nint(time*adm/dt)
c
c        set number of frames of various types of output
      nframes=11
      if (nframes.gt.1) then
      fi=real(ntime)/real(nframes-1)
      nf=nint(fi)
      endif
c
      nmovie=2
c
c        wave amplitude output
      nwaves=300
      if (ntime.lt.nwaves) nwaves=ntime-2
      if (nwaves.gt.1) then
      fiwave=real(ntime)/real(nwaves-1)
      nfwave=nint(fiwave)
      endif
c
c        and for movie
      nmwaves=50
      if (nmwaves.gt.1) then
      fimwave=real(ntime)/real(nmwaves-1)
      nfmwave=nint(fimwave)
      endif
c
c        equator data output
      nequ=2
      if (ntime.lt.nequ) nequ=ntime-2
      if (nequ.gt.1) then
      fiequ=real(ntime)/real(nequ-1)
      nfequ=nint(fiequ)
      endif
c
c        apparent horizon data output
      nah=300
      if (ntime.lt.nah) nah=ntime-2
      if (nah.gt.1) then
      fiah=real(ntime)/real(nah-1)
      nfah=nint(fiah)
      endif
c
c        radiation 2-sphere pattern output
      nradpat=300
      if (ntime.lt.nradpat) nradpat=ntime-2
      if (nradpat.gt.1) then
      firp=real(ntime)/real(nradpat-1)
      nfirp=nint(firp)
      endif
c
c        oh, yeah
      sqrt2=sqrt(2.)
c
      if (waytogo.eq.'ivpcrank') go to 99
c        write vital run info to terminal
      write(5,*)
      write(5,*)'           *** Run Info ***'
      write(5,*)
      write(5,*) 'grid size (eta X theta)  ',netamax,nthtamax
      write(5,*) 'maximum eta value   ',etamax
      write(5,*) 'scale parameter m   ',scaleprm
      write(5,*) 'time of run (/ADM)  ',time
      write(5,*)
      write(5,*) 'number of time steps  ',ntime
      write(5,*) 'number of data outputs  ',nframes
      write(5,*) 'number of apparent horizon outputs  ',nah
      write(5,*) 'number of wave amp outputs  ',nwaves
      write(5,*) 'number of equator outputs  ',nequ
      write(5,*)
      write(5,*) 'delta eta     ',deta
      write(5,*) 'delta theta   ',dtheta
      write(5,*) 'delta t       ',dt
      write(5,*)
      write(5,*) 'distortion coefficients:  ','da1=',da1,'   dr1=',dr1,
     &           '    dw1=',dw1
      write(5,*) '                          ','da2=',da2,'   dr2=',dr2,
     &           '    dw2=',dw2
      write(5,*)
      write(5,*)
c
      if (dervans.eq.'s') then
      write(5,*)'spatial differencing is second order'
      elseif (dervans.eq.'f') then
      write(5,*)'spatial differencing is fourth order'
      elseif (dervans.eq.'s24') then
      write(5,*)'spatial differencing is scheme24'
      endif
c
      if (timesche.eq.'m') then
      write(5,*)'time differencing is McCormack'
      elseif(timesche.eq.'b') then
      write(5,*)'time differencing is Brailovskaya'
      endif
      write(5,*)
c
c     write(5,*)'*** Code is evolving metric C with no shift ***'
c     write(5,*)
c
c
   99 return
      end
c
c
c
      subroutine initial
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /shiftblk/ beta1(id,jd),beta2(id,jd),omega(id,jd)
      common /geobndbk/  abnd(500), bbnd(500), cbnd(500), dbnd(500),
     &                  habnd(500),hbbnd(500),hcbnd(500),hdbnd(500)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /waveblck/ nwaves,nfwave,fiwave,nmwaves,
     &                  nfmwave,fimwave
      common /equblock/ nequ,nfequ,fiequ
      common /horiznbk/ nah,nfah,fiah
      common /rpatblck/ nradpat,nfirp,firp
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
c        open files for output
      open(1,file='a.dat',status='unknown')
      open(2,file='b.dat',status='unknown')
      open(3,file='c.dat',status='unknown')
      open(4,file='d.dat',status='unknown')
      open(6,file='ha.dat',status='unknown')
      open(7,file='hb.dat',status='unknown')
      open(8,file='hc.dat',status='unknown')
      open(9,file='hd.dat',status='unknown')
      open(10,file='lpse.dat',status='unknown')
      open(11,file='omga.dat',status='unknown')
      open(12,file='stf1.dat',status='unknown')
      open(21,file='equa.dat',status='unknown')
      open(30,file='movi.dat',status='unknown')
      open(31,file='ah.dat',status='unknown')
      open(32,file='ahd.dat',status='unknown')
      open(33,file='pe.dat',status='unknown')
      open(34,file='mah.dat',status='unknown')
      open(35,file='gcah.dat',status='unknown')
      open(36,file='emtr.dat',status='unknown')
c
c
c        initialize metric and extrinsic curvature
      call computeq
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      q=f(i,j)
      e2q=exp(2.*q)
c
      a(i,j)=e2q
      b(i,j)=e2q
      c(i,j)=0.0
      d(i,j)=1.0
c
      ha(i,j)=0.0
      hb(i,j)=0.0
      hc(i,j)=0.0
      hd(i,j)=0.0
   30 continue
c
c        set initial outer boundary values
      schwzlps=tanh(0.5*(eta(netamax)+deta))
      do 40 j=1,nthtamax
      abnd(j)=1.0
      bbnd(j)=1.0
      cbnd(j)=0.0
      dbnd(j)=1.0
c
      habnd(j)=0.0
      hbbnd(j)=0.0
      hcbnd(j)=0.0
      hdbnd(j)=0.0
c
      alphabnd(j)=schwzlps
      omegabnd(j)=1.0e-20
   40 continue
c
c        set initial lapse, shift potential, and shift vector
      bnd=tanh(0.5*(eta(netamax)+deta))
      do 50 j=1,nthtamax
      do 50 i=1,netamax
      alpha(i,j)=bnd
      omega(i,j)=1.0e-20
   50 continue
c
c        compute initial lapse 
      call allball2
      call maxlapse
c
c        compute initial shift
      call shiftc0
c
c        open files for Andrew's wave detector
      if (.true.) then
      open(51,file='r15.dat',status='unknown')
      open(52,file='r30.dat',status='unknown')
      open(53,file='r45.dat',status='unknown')
      open(54,file='r60.dat',status='unknown')
c
c        header for file is radius in units of M
      r1=15.*scaleprm
      etad=dlog(2.*r1/scaleprm)
      idx=int(etad/deta)+1
      e=eta(idx)
      r1=0.5*scaleprm*exp(e)
c
      r2=30.*scaleprm
      etad=dlog(2.*r2/scaleprm)
      idx=int(etad/deta)+1
      e=eta(idx)
      r2=0.5*scaleprm*exp(e)
c
      r3=45.*scaleprm
      etad=dlog(2.*r3/scaleprm)
      idx=int(etad/deta)+1
      e=eta(idx)
      r3=0.5*scaleprm*exp(e)
c
      r4=60.*scaleprm
      etad=dlog(2.*r4/scaleprm)
      idx=int(etad/deta)+1
      e=eta(idx)
      r4=0.5*scaleprm*exp(e)
c
      write(51,*) nwaves
      write(52,*) nwaves
      write(53,*) nwaves
      write(54,*) nwaves
c
      write(51,*) r1
      write(52,*) r2
      write(53,*) r3
      write(54,*) r4
      endif
c
c        open file for zerilli function
      open(61,file='zerl.dat',status='unknown')
      open(62,file='zer2.dat',status='unknown')
      write(61,*) nmwaves
      write(62,*) nmwaves
      write(61,*) netamax
      write(62,*) netamax
      do 35 i=1,netamax
      write(61,*) eta(i)
      write(62,*) eta(i)
   35 continue
c
c        prepare other data files
      do 82 ifile=1,12
      if (ifile.ne.5) then
      if (qrun.eq.'y') then
      write(ifile,*) nframes
      else
      write(ifile,*) 1
      endif
      write(ifile,*) netamax
      write(ifile,*) nthtamax
      do 80 i=1,netamax
   80 write(ifile,*) eta(i)
      do 81 j=1,nthtamax
   81 write(ifile,*) theta(j)
      endif
   82 continue
c
c        prepare output for movie data
      write(30,*) nmovie
      write(30,*) netamax
      write(30,*) nthtamax
      do 83 i=1,netamax
   83 write(30,*) eta(i)
      do 84 j=1,nthtamax
   84 write(30,*) theta(j)
c
c        prepare output for equ data
      write(21,*) nequ
      write(21,*) netamax
      do 85 i=1,netamax
   85 write(21,*) eta(i)
c
c        prepare output for apparent horizon
      write(31,*) nah
      write(33,*) nah
      write(34,*) nah
      write(35,*) nah
      write(36,*) nah
c
      write(31,*) nthtamax
      write(35,*) nthtamax
      write(36,*) nthtamax
c
      if (qrun.eq.'y') then
      write(32,*) nframes
      else
      write(32,*) 1
      endif
      write(32,*) nthtamax
c
      do 86 j=1,nthtamax
      write(31,*) theta(j)
      write(32,*) theta(j)
      write(35,*) theta(j)
      write(36,*) theta(j)
   86 continue
c
c        prepare files for radiation pattern data
      open(71,file='rpat.dat',status='unknown')
      write(71,*) nradpat
      write(71,*) nthtamax
      do 87 j=1,nthtamax
   87 write(71,*) theta(j)
c
c
      return
      end
c
c
c
      subroutine computeq
      include 'param.h'
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /distblck/ da1,dr1,dw1,da2,dr2,dw2
c
c
c        compute function q and put it in array f
      do 10 j=1,nthtamax
c        specify angular form of q
c
c        sin(theta)^2 dependence
      qth=sn(j)**2.
      q1th=2.*acs(j)*sn(j)
      q2th=2.*(acs(j)*acs(j)-sn(j)*sn(j))
c
c        sin(theta)^4 dependence
c     qth=sn(j)**4.
c     q1th=4.*acs(j)*sn(j)**3.
c     q2th=12.*(sn(j)*acs(j))**2.-4.*sn(j)**4.
c
c        (sin(theta)*cos(theta))**2 dependence
c     qth=(sn(j)*acs(j))**2.
c     q1th=2.*sn(j)*acs(j)*(acs(j)**2.-sn(j)**2.)
c     q2th=2.*(acs(j)**4.-6.*(sn(j)*acs(j))**2.+sn(j)**4.)
c
      do 10 i=1,netamax
      e=eta(i)
c
      em1=(e-dr1)/dw1
      ep1=(e+dr1)/dw1
c
      expp1=exp(-ep1*ep1)
      expm1=exp(-em1*em1)
c
      em2=(e-dr2)/dw2
      ep2=(e+dr2)/dw2
c
      expp2=exp(-ep2*ep2)
      expm2=exp(-em2*em2)
c
c        put q here
      q1=da1*qth*(expp1+expm1)
      q2=da2*qth*(expp2+expm2)
c
c        eta derivative
      q1n1=da1*qth*(-2.*em1*expm1-2.*ep1*expp1)/dw1
      q1n2=da2*qth*(-2.*em2*expm2-2.*ep2*expp2)/dw2
c
c        theta derivative
      q1t1=da1*q1th*(expp1+expm1)
      q1t2=da2*q1th*(expp2+expm2)
c
c        eta second derivative
      q2n1=da1*qth*(expp1*(-2.+4.*ep1*ep1)
     &             +expm1*(-2.+4.*em1*em1))/(dw1*dw1)
      q2n2=da2*qth*(expp2*(-2.+4.*ep2*ep2)
     &             +expm2*(-2.+4.*em2*em2))/(dw2*dw2)
c
c        theta second derivative
      q2t1=da1*q2th*(expp1+expm1)
      q2t2=da2*q2th*(expp2+expm2)
c
c        mixed second derivative
      q2m1=da1*q1th*(-2.*em1*expm1-2.*ep1*expp1)/dw1
      q2m2=da2*q1th*(-2.*em2*expm2-2.*ep2*expp2)/dw2
c
      f(i,j)=q1+q2
      f1(i,j)=q1n1+q1n2
      f2(i,j)=q1t1+q1t2
      f11(i,j)=q2n1+q2n2
      f22(i,j)=q2t1+q2t2
      f12(i,j)=q2m1+q2m2
   10 continue
c
c     call frthdvbk(f,f1,f2,f11,f22,f12,1.,netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine swzlapse
      include 'param.h'
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      e=eta(i)
      th=tanh(0.5*e)
c
      alpha(i,j)=th
      alpha1n(i,j)=0.5*(1.-th*th)
      alpha2n(i,j)=-0.5*th*(1.-th*th)
c
      alpha1t(i,j)=0.0
      alpha2t(i,j)=0.0
      alpha2m(i,j)=0.0
   10 continue
c
      return
      end
c
c
c
      subroutine ivp
      include 'param.h'
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c
      call ivpaxis1(ir1)
      if (ir1.eq.1) call ivptest1(ans1)
c
      call ivpaxis3(ir3)
      if (ir3.eq.1) call ivptest1(ans3)
      write(5,*)ir1,ans1,ir3,ans3
c
      if (ir1.eq.1.or.ir3.eq.1) ir=1
c
      if (ans3.le.ans1) then
      call ivpaxis3(ix)
      method=3
      else
      call ivpaxis1(ix)
      method=1
      endif     
c
      if (waytogo.ne.'ivpcrank') then
      if (ir.eq.1) then
      write(5,*)'ivp method used:  ',method
      go to 99 
      else
      write(5,*)'ivp solver failure: code aborted (bub)'
      call alldone
      stop
      endif
      endif
c
c
   99 return
      end
c
c
c
      subroutine ivpaxis1(iflag)
      include 'param.h'
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
c        construct functions f
      call setivpfr
c
c        construct coefficent matrices
      call setmcivr
c
c        compute initial guess for psi
      do 20 j=1,nthtamax
      do 20 i=1,netamax
      e=eta(i)
      psi(i,j)=2.*sqrthfm*cosh(0.5*e)
   20 continue
c
      tol=1.0e-10
c
      call multigrid(netamax,1,netamax,nthtamax,1,nthtamax,
     &               cc,cn,cs,cw,ce,cnw,cne,csw,cse,psi,rhs,tol,
     &               1,iflag)
c
      if (iflag.eq.0) go to 99
c
c        compute derivatives of the conformal factor
c        (phi is temporarily stored in f)
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      cf=psi(i,j)
      psim4(i,j)=cf**(-4.)
      f(i,j)=4.*dlog(cf)
   30 continue
c
      call frthdvbk(f,phi1n,phi1t,phi2n,phi2t,phi2m,1.,
     &              netamax,nthtamax)
c
c        output conformal factor
      open(97,file='conf.dat',status='unknown')
      write(97,*) 1
      write(97,*) netamax
      write(97,*) nthtamax
      do 50 i=1,netamax
   50 write(97,*) eta(i)
      do 51 j=1,nthtamax
   51 write(97,*) theta(j)
      do 53 j=1,nthtamax
      do 53 i=1,netamax
      write(97,*) psi(i,j)
   53 continue
c
      close(97)
c
c
   99 return
      end
c
c
c
      subroutine ivpaxis3(iflag)
      include 'param.h'
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /mcoefblk/ cc(id,jd),cn(id,jd),cs(id,jd),cw(id,jd),
     &               ce(id,jd),cnw(id,jd),cne(id,jd),csw(id,jd),
     &                  cse(id,jd),rhs(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /sclbndbk/ alphabnd(500),omegabnd(500),phibnd(500)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
c        construct functions f
      call setivpf3
c
c        construct coefficient matrices
      call setmciv3
c
c        compute initial guess for psi (f1 used temporarily)
      do 20 j=1,nthtamax
      do 20 i=1,netamax
      f1(i,j)=1.
   20 continue
c
      tol=1.0e-10
c
      call multigrid(netamax,1,netamax,nthtamax,1,nthtamax,
     &               cc,cn,cs,cw,ce,cnw,cne,csw,cse,f1,rhs,tol,
     &               1,iflag)
c
      if (iflag.eq.0) go to 99
c
c        compute derivatives of the conformal factor
c        (phi is temporarily stored in f)
      do 30 j=1,nthtamax
      do 30 i=1,netamax
      e=eta(i)
      sscf=2.*sqrthfm*cosh(0.5*e)
c
      cf=sscf*f1(i,j)
      psi(i,j)=cf
      psim4(i,j)=cf**(-4.)
      f(i,j)=4.*dlog(cf)
   30 continue
c
      e=eta(netamax)
      fac=sqrthfm*exp(0.5*e)
      do 40 j=1,nthtamax
      psibnd=psi(netamax-1,j)+2.*deta*(fac-0.5*psi(netamax,j))
      phibnd(j)=4.*dlog(psibnd)
   40 continue
c
c     call scnddrvs(f,phi1n,phi1t,phi2n,phi2t,phi2m,1.,phibnd,
c    &              netamax,nthtamax)
      call frthdvbk(f,phi1n,phi1t,phi2n,phi2t,phi2m,1.,
     &              netamax,nthtamax)
c
c        output conformal factor
      open(97,file='conf.dat',status='unknown')
      write(97,*) 1
      write(97,*) netamax
      write(97,*) nthtamax
      do 50 i=1,netamax
   50 write(97,*) eta(i)
      do 51 j=1,nthtamax
   51 write(97,*) theta(j)
      do 53 j=1,nthtamax
      do 53 i=1,netamax
      write(97,*) psi(i,j)
   53 continue
c
      close(97)
c
c
   99 return
      end
c
c
c
      subroutine setivpfr
      include 'param.h'
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
      call computeq
c
      do 20 j=1,nthtamax
      do 20 i=1,netamax
   20 f(i,j)=0.25*deta*deta*(f11(i,j)+f22(i,j)-1.0)
c
      do 30 j=2,nthtamax
      do 30 i=1,netamax
      f11(i,j)=1.
      f22(i,j)=xi*xi
      f12(i,j)=0.
c
      f1(i,j)=0.
      f2(i,j)=0.5*xi*deta*ct(j)
c
      g(i,j)=0.
   30 continue
c
      do 40 i=1,netamax
      f11(i,1)=1.
      f22(i,1)=2.*xi*xi
      f12(i,1)=0.
c
      f1(i,1)=0.
      f2(i,1)=0.
c
      g(i,1)=0.
   40 continue
c
      return
      end
c
c
c
      subroutine setivpf3
      include 'param.h'
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
      call computeq
c
      do 20 j=1,nthtamax
      do 20 i=1,netamax
   20 f(i,j)=0.25*deta*deta*(f11(i,j)+f22(i,j))
c
      do 30 j=2,nthtamax
      do 30 i=1,netamax
      e=eta(i)
c
      f11(i,j)=1.
      f22(i,j)=xi*xi
      f12(i,j)=0.
c
      f1(i,j)=0.5*deta*tanh(0.5*e)
      f2(i,j)=0.5*xi*deta*ct(j)
c
      g(i,j)=0.
   30 continue
c
      do 40 i=1,netamax
      e=eta(i)
c
      f11(i,1)=1.
      f22(i,1)=2.*xi*xi
      f12(i,1)=0.
c
      f1(i,1)=0.5*deta*tanh(0.5*e)
      f2(i,1)=0.
c
      g(i,1)=0.
   40 continue
c
      return
      end
c
c
c
      subroutine output
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /ricciblk/ ricci11(id,jd),ricci12(id,jd),
     &                  ricci22(id,jd),ricci33(id,jd)
      common /lapseblk/ alpha(id,jd)
      common /lpsedblk/ alpha1n(id,jd),alpha1t(id,jd),
     &                  alpha2n(id,jd),alpha2t(id,jd),
     &                  alpha2m(id,jd)
      common /shiftblk/ beta1(id,jd),beta2(id,jd),omega(id,jd)
      common /shftdblk/ beta1n(id,jd),beta1t(id,jd),
     &                  beta2n(id,jd),beta2t(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /fdotblck/ adot(id,jd),bdot(id,jd),
     &                  cdot(id,jd),ddot(id,jd)
      common /hdotblck/ hadot(id,jd),hbdot(id,jd),
     &                  hcdot(id,jd),hddot(id,jd)
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /spinblck/ rho(id,jd),amu(id,jd),sigma(id,jd),
     &                  alambda(id,jd)
      common /elmagblk/ e11(id,jd),e12(id,jd),e22(id,jd),
     &                  e33(id,jd),b13(id,jd),b23(id,jd)
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
      io=nint(jtime/fi)
      nf=nint(fi*(1.+io))
c
      thetime=dt*jtime/adm
      write (5,101) thetime,time
  101 format(f7.3,' ADM out of ',f7.3,' ADM')
c
c     call spincoef
c     call elmagnoc
c
c     call tracek(f2,netamax,nthtamax)
c     call hamilton(f1,netamax,nthtamax)
c     call belrobvc(f11,f22,f12,netamax,nthtamax)
c     call momentum(f11,f22,netamax,nthtamax)
      call computeq
c
c        find apparent horizon for output
      call cook
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      write(1,100) a(i,j)
      write(2,100) b(i,j)
c     write(3,100) c(i,j)
      write(4,100) d(i,j)
      write(6,100) ha(i,j)
      write(7,100) hb(i,j)
      write(8,100) hc(i,j)
      write(9,100) hd(i,j)
      write(10,100) alpha(i,j)
      write(11,100) omega(i,j)
c
c        construct incidental quantities
      rhat=-2.*exp(-2.*f(i,j)-2.*eta(i))
     &    *(f11(i,j)+f22(i,j))
      write(12,100) rhat
   10 continue
c
      do 20 j=1,nthtamax
      write(32,100) h(j)
   20 continue
  100 format(e22.10e4)
c
c
      return
      end
c
c
c
      subroutine eqout
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
      common /equblock/ nequ,nfequ,fiequ
c
c
      iw=nint(jtime/fiequ)
      nfequ=nint(fiequ*(1.+iw))
c
      thetime=dt*jtime/adm
c
      write(21,*) thetime
c
      do 10 i=1,netamax
      write(21,101) a(i,nthtamax),b(i,nthtamax),
     &              c(i,nthtamax),d(i,nthtamax)
      write(21,101) ha(i,nthtamax),hb(i,nthtamax),
     &              hc(i,nthtamax),hd(i,nthtamax)
c
      write(21,101) a(i,ntm1),b(i,ntm1),
     &              c(i,ntm1),d(i,ntm1)
      write(21,101) ha(i,ntm1),hb(i,ntm1),
     &              hc(i,ntm1),hd(i,ntm1)
   10 continue
c
  101 format(4e22.14)
c
c
      return
      end
c
c
c
      subroutine movieout
      include 'param.h'
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /spinblck/ rho(id,jd),amu(id,jd),sigma(id,jd),
     &                  alambda(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
      call elmagnoc
      call belrobvc(f,f22,f12,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      write(30,100) f(i,j)
   10 continue
c
  100 format(e22.10e4)
c
c
      return
      end
c
c
c
      subroutine ahout
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /emtrpblk/ gcah(110),curvet(110),arcleng(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /horiznbk/ nah,nfah,fiah
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
      iw=nint(jtime/fiah)
      nfah=nint(fiah*(1.+iw))
c
      thetime=jtime*dt/adm
c
      call cook
c
c        output coordinate location of apparent horizon
      write(31,*) thetime
      do 10 j=1,nthtamax
   10 write(31,*) h(j)
c
c        output other data
      call circumf(ans)
      write(33,*) thetime,ans
c
      call ahmass(ans)
      write(34,*) thetime,ans/adm
c
      call embdtrp
      write(35,*) thetime
      write(36,*) thetime
      do 20 j=1,nthtamax
      write(35,*) gcah(j)
      write(36,*) arcleng(j),curvet(j)
   20 continue
c
c
      return
      end
c
c
c
      subroutine tracek(u,ne,nt)
      include 'param.h'
      dimension u(ne,nt)
c
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
c        this subroutine constructs the trace of the extrinsic 
c        curvature and puts the result in array u.
c
      do 10 j=1,nthtamax
      do 10 i=1,netamax
      u(i,j)=em1(i,j)*(b(i,j)*ha(i,j)+a(i,j)*hb(i,j)-2.*c(i,j)*hc(i,j))
     &      +hd(i,j)/d(i,j)
   10 continue
c
c
      return
      end
c
c
c
      subroutine hamilton(u,ne,nt)
      include 'param.h'
      dimension u(ne,nt)
c
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /eblock  / em1(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /ricciblk/ ricci11(id,jd),ricci12(id,jd),
     &                  ricci22(id,jd),ricci33(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        this subroutine constructs the Hamiltonian constraint and
c        returns it in array u.
c
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      r1=em1(i,j)*(b(i,j)*ricci11(i,j)+a(i,j)*ricci22(i,j)
     &            -2.*c(i,j)*ricci12(i,j))
     &  +ricci33(i,j)/d(i,j)
c
      scalcurv=psim4(i,j)*r1
c
      trk=em1(i,j)*(b(i,j)*ha(i,j)+a(i,j)*hb(i,j)-2.*c(i,j)*hc(i,j))
     &   +hd(i,j)/d(i,j)
c
      r1=b(i,j)*b(i,j)*ha(i,j)*ha(i,j)
      r2=a(i,j)*a(i,j)*hb(i,j)*hb(i,j)
      r3=2.*(a(i,j)*b(i,j)+c(i,j)*c(i,j))*hc(i,j)*hc(i,j)
      r4=-4.*b(i,j)*c(i,j)*ha(i,j)*hc(i,j)
      r5=-4.*a(i,j)*c(i,j)*hb(i,j)*hc(i,j)
      r6=2.*c(i,j)*c(i,j)*ha(i,j)*hb(i,j)
      r7=hd(i,j)*hd(i,j)/(d(i,j)*d(i,j))
      sum=r1+r2+r3+r4+r5+r6
      squarek=em1(i,j)*em1(i,j)*sum+r7
c
      u(i,j)=scalcurv-squarek+trk*trk
   10 continue
c
c
c        and on axis
      do 20 i=1,netamax
      scalcurv=psim4(i,1)*(ricci11(i,1)/a(i,1)+ricci22(i,1)/b(i,1)
     &                    +ricci33(i,1)/d(i,1))
c
      trk=ha(i,1)/a(i,1)+2.*hb(i,1)/b(i,1)
c
      squarek=(ha(i,1)/a(i,1))**2+2.*(hb(i,1)/b(i,1))**2
c
      u(i,1)=scalcurv-squarek+trk*trk
   20 continue
c
c
      return
      end
c
c
c
      subroutine momentum(u1,u2,ne,nt)
      include 'param.h'
      dimension u1(ne,nt),u2(ne,nt)
c
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /chrisblk/ g111(id,jd),g112(id,jd),g122(id,jd),
     &                  g133(id,jd),g211(id,jd),g212(id,jd),
     &                  g222(id,jd),g233(id,jd),g313(id,jd),
     &                  g323(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c        this subroutine constructs the down components of the
c        momentum constraint and returns the eta component in the
c        array u1 and the theta component in array u2.
c        (the 3-metric is assumed to be diagonal in this subroutine)
c
c
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      r1= hc1t(i,j)+hc(i,j)*phi1t(i,j)
      r2=-hb1n(i,j)-hb(i,j)*phi1n(i,j)
      r3=-g122(i,j)*ha(i,j)-g222(i,j)*hc(i,j)
     &   +g112(i,j)*hc(i,j)+g212(i,j)*hb(i,j)
c
      r4=-hd1n(i,j)-hd(i,j)*phi1n(i,j)
      r5=-g133(i,j)*ha(i,j)-g233(i,j)*hc(i,j)+g313(i,j)*hd(i,j)
c
      u1(i,j)=(r1+r2+r3)/b(i,j)+(r4+r5)/d(i,j)
c
      r1= hc1n(i,j)+hc(i,j)*phi1n(i,j)
      r2=-ha1t(i,j)-ha(i,j)*phi1t(i,j)
      r3=-g211(i,j)*hb(i,j)-g111(i,j)*hc(i,j)
     &   +g112(i,j)*ha(i,j)+g212(i,j)*hc(i,j)
c
      r4=-hd1t(i,j)-hd(i,j)*phi1t(i,j)-2.*ct(j)*hd(i,j)
      r5=-g133(i,j)*hc(i,j)-g233(i,j)*hb(i,j)+g323(i,j)*hd(i,j)
c
      u2(i,j)=(r1+r2+r3)/a(i,j)+(r4+r5)/d(i,j)
   10 continue
c
c
c        and on axis
      do 20 i=1,netamax
      g122axis=-0.5*(b1n(i,1)+b(i,1)*phi1n(i,1))/a(i,1)
      g212axis= 0.5*(b1n(i,1)/b(i,1)+phi1n(i,1))
c
      r1=hc1t(i,1)-hb1n(i,1)-hb(i,1)*phi1n(i,1)
      r2=-g122axis*ha(i,1)+g212axis*hb(i,1)
      u1(i,1)=2.*(r1+r2)/b(i,1)
      u2(i,1)=0.
   20 continue
c
c
      return
      end
c
c
c
      subroutine belrobvc(u1,u2,u3,ne,nt)
      include 'param.h'
      dimension u1(ne,nt),u2(ne,nt),u3(ne,nt)
c
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /elmagblk/ e11(id,jd),e12(id,jd),e22(id,jd),
     &                  e33(id,jd),b13(id,jd),b23(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
c
c        this subroutine constructs the components of the 
c        Bel-Robinson vector.  The up eta component is put in u1,
c        the up theta component is put in u2, and the norm is in u3.
c        (the 3-metric is assumed to be diagonal)
c        Each is appropriate units of the ADM mass
c
      fac1=adm**5.
      fac2=adm**4.
c
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      fac=1./sqrt(a(i,j)*b(i,j)*d(i,j))
c
      r1=-b13(i,j)*e12(i,j)/a(i,j)
      r2=-b23(i,j)*e22(i,j)/b(i,j)
      r3= b23(i,j)*e33(i,j)/d(i,j)
c
      bre=fac*(r1+r2+r3)*psim4(i,j)
      u1(i,j)=fac1*bre
c
      r1= b13(i,j)*e11(i,j)/a(i,j)
      r2= b23(i,j)*e12(i,j)/b(i,j)
      r3=-b13(i,j)*e33(i,j)/d(i,j)
c
      brt=fac*(r1+r2+r3)*psim4(i,j)
      u2(i,j)=fac1*brt
c
      r1=sqrt(bre*bre*a(i,j)+brt*brt*b(i,j))
      u3(i,j)=fac2*psi(i,j)*psi(i,j)*r1
   10 continue
c
c        and on axis
      do 20 i=1,netamax
      u1(i,1)=0.
      u2(i,1)=0.
      u3(i,1)=0.
   20 continue
c
c
      return
      end
c
c
c
      subroutine ivptest1(hmass)
      include 'param.h'
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real h1(500),h2(500),hres(id,jd)
      real psi1n(id,jd),psi1t(id,jd),psi2n(id,jd),
     &     psi2t(id,jd),psi2m(id,jd)
c
c
c        this subroutine tests the solution of the Hamiltonian
c        constraint by fourth order differencing and performing
c        a volume integral over the initial slice
c
      call frthdvbk(psi,psi1n,psi1t,psi2n,psi2t,psi2m,1.,
     &              netamax,nthtamax)
c
      call computeq
c
      do 10 j=2,nthtamax
      do 10 i=1,netamax
      r1=psi2n(i,j)+psi2t(i,j)+psi1t(i,j)*ct(j)
     &  +0.25*psi(i,j)*(f11(i,j)+f22(i,j)-1.)
      r2=8.*psi(i,j)*sn(j)
      hres(i,j)=abs(r1*r2)
   10 continue
c
      do 11 i=1,netamax
   11 hres(i,1)=0.
c
      do 30 i=1,netamax
      do 25 j=1,nthtamax
   25 h1(j)=hres(i,j)
      call simpson(h1,dtheta,nthtamax,ans)
      h2(i)=ans
   30 continue
c
      call simpson(h2,deta,netamax,ans)
c
      hmass=2.*ans
c
c
      return
      end
c
c
c
      subroutine ivptest2
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /func2blk/ f(id,jd),f1(id,jd),f2(id,jd),
     &                  f11(id,jd),f22(id,jd),f12(id,jd),
     &                  g(id,jd)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real psi1n(id,jd),psi1t(id,jd),psi2n(id,jd),
     &     psi2t(id,jd),psi2m(id,jd)
c
c
c        this subroutine tests the solution of the Hamiltonian
c        constraint by the use of alternative differencing
c        (from that used to solve the problem).
c
      call frthdvbk(psi,psi1n,psi1t,psi2n,psi2t,psi2m,1.,
     &              netamax,nthtamax)
c
      call computeq
c
      hmax=-1.
      sum=0.
c
      j=1
      do 10 i=1,netamax
      res=abs(psi2n(i,j)+2.0*psi2t(i,j)
     &   +0.25*psi(i,j)*(f11(i,j)+f22(i,j)-1.))
      res=8.*res/(a(i,j)*psi(i,j)**5)
      hmax=dmax1(hmax,res)
      sum=sum+res
   10 continue
c
c
      do 11 j=2,nthtamax
      do 11 i=1,netamax
      res=abs(psi2n(i,j)+psi2t(i,j)+psi1t(i,j)*ct(j)
     &   +0.25*psi(i,j)*(f11(i,j)+f22(i,j)-1.))
      res=8.*res/(a(i,j)*psi(i,j)**5)
      hmax=dmax1(hmax,res)
      sum=sum+res
   11 continue
c
      write(5,*)
      write(5,*)'log maximum residual of finite differenced equation'
      write(5,*)dlog10(hmax/(scaleprm**2.))
c
      sum=sum/(float(netamax*nthtamax))
      write(5,*)'log average residual of finite differenced equation'
      write(5,*)dlog10(sum/(scaleprm**2.))
c
c
      return
      end
c
c
c
      subroutine firstah
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
      character*8 qrun,waytogo,dervans,timesche
      common /charblck/ qrun,waytogo,dervans,timesche
c
c        give initial guess for first apparent horizon
      if (waytogo.eq.'interact') then
      write(5,*)
      write(5,*)'guess for initial apparent horizon?'
      read(5,*) ah
      else
      ah=0.3
      endif
c
      do 10 j=1,nthtamax
   10 h(j)=ah
c
      call cook
c
      return
      end
c
c
c
      subroutine cook
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /diagblck/ dd(110),du(110),dl(110)
      common /residblk/ res(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real rhs(110)
c
c        this subroutine computes the location of the apparent 
c        horizon on a general slice.  It uses a method due to G. Cook
c        (see his PhD thesis for details).
c
c
c        iterate Newton-Raphson procedure
      do 90 inr=1,50
c
c        compute diagonals in linearized coefficient matrix
      call lincoeff
c
c        get residuals of current estimate
      call diffh2
      call residual
c
c        solve tridiagonal system for adjustment vector
      do 10 j=1,nthtamax
   10 rhs(j)=-res(j)
      call tridiag(dd,du,dl,rhs,dh)
c
      call newh
c
c        check to see if equation is satisfied
      call check(iflag)
      if (iflag.eq.1) then
c     write(5,*)'equation satisfied on iteration ',inr
      go to 99
      endif
   90 continue
c
      write(5,*)'Newton-Raphson in trouble'
c
   99 continue
c
c
      return
      end
c
c
c
      subroutine newh
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
c        construct new apparent horizon
c
      isum=0
c
      do 10 j=1,nthtamax
      h(j)=h(j)+dh(j)
c
      i=int(sign(1.,h(j)))
      isum=i+isum
   10 continue
c
      if (isum.eq.-nthtamax) then
      do 20 j=1,nthtamax
   20 h(j)=-h(j)
      endif
c
c
      return
      end
c
c
c
      subroutine lincoeff
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /diagblck/ dd(110),du(110),dl(110)
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c        compute second order derivatives
      call diffh2
c
c        get interpolated values of all variables along
c        current guess of apparent horizon
c     call interpl1
      call interpl3
c
c        main diagonal
      do 10 j=2,nthtamax
      h2=h1t(j)*h1t(j)
      h3=h2*h1t(j)
c
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
      trk1n=ha1ni(j)/ai(j)-hai(j)*a1ni(j)/(ai(j)*ai(j))
     &     +hb1ni(j)/bi(j)-hbi(j)*b1ni(j)/(bi(j)*bi(j))
     &     +hd1ni(j)/di(j)-hdi(j)*d1ni(j)/(di(j)*di(j))
c
      f=ai(j)*ai(j)*h2/bi(j)+ai(j)
      fh=sqrt(f)
      f1n=2.*ai(j)*a1ni(j)*h2/bi(j)-b1ni(j)*h2*(ai(j)/bi(j))**2.
     &   +a1ni(j)
c
      s1 =-2./(dtheta*dtheta)
c
      r1 =a1ni(j)*ct(j)/bi(j)
      r2 =-ai(j)*b1ni(j)*ct(j)/(bi(j)*bi(j))
      r3 =0.5*a1ni(j)*d1ti(j)/(bi(j)*di(j))
      r4 =0.5*ai(j)*d2mi(j)/(bi(j)*di(j))
      r5 =-0.5*ai(j)*d1ti(j)*b1ni(j)/(bi(j)*bi(j)*di(j))
      r6 =-0.5*ai(j)*d1ti(j)*d1ni(j)/(bi(j)*di(j)*di(j))
      r7 =0.5*a2mi(j)/bi(j)
      r8 =-0.5*a1ti(j)*b1ni(j)/(bi(j)*bi(j))
      r9 =phi2mi(j)*ai(j)/bi(j)
      r10=phi1ti(j)*a1ni(j)/bi(j)
      r11=-phi1ti(j)*ai(j)*b1ni(j)/(bi(j)*bi(j))
      s2=h3*(r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11)
c
      r1=-0.5*d2ni(j)/di(j)
      r2=0.5*d1ni(j)*d1ni(j)/(di(j)*di(j))
      r3=-b2ni(j)/bi(j)
      r4=b1ni(j)*b1ni(j)/(bi(j)*bi(j))
      r5=0.5*a2ni(j)/ai(j)
      r6=-0.5*a1ni(j)*a1ni(j)/(ai(j)*ai(j))
      r7=-phi2ni(j)
      s3=h2*(r1+r2+r3+r4+r5+r6+r7)
c
      r1=0.5*d2mi(j)/di(j)
      r2=-0.5*d1ti(j)*d1ni(j)/(di(j)*di(j))
      r3=-0.5*b2mi(j)/bi(j)
      r4=0.5*b1ti(j)*b1ni(j)/(bi(j)*bi(j))
      r5=a2mi(j)/ai(j)
      r6=-a1ti(j)*a1ni(j)/(ai(j)*ai(j))
      r7=phi2mi(j)
      s4=h1t(j)*(r1+r2+r3+r4+r5+r6+r7)
c
      r1=-0.5*b1ni(j)*d1ni(j)/(ai(j)*di(j))
      r2=-0.5*bi(j)*d2ni(j)/(ai(j)*di(j))
      r3=0.5*bi(j)*d1ni(j)*a1ni(j)/(ai(j)*ai(j)*di(j))
      r4=0.5*bi(j)*d1ni(j)*d1ni(j)/(ai(j)*di(j)*di(j))
      r5=-0.5*b2ni(j)/ai(j)
      r6=0.5*b1ni(j)*a1ni(j)/(ai(j)*ai(j))
      r7=-phi1ni(j)*b1ni(j)/ai(j)
      r8=-phi2ni(j)*bi(j)/ai(j)
      r9=phi1ni(j)*bi(j)*a1ni(j)/(ai(j)*ai(j))
      s5=r1+r2+r3+r4+r5+r6+r7+r8+r9
c
      r1=-bi(j)*hai(j)/(ai(j)*ai(j))
      r2=-h2*hbi(j)/bi(j)
      r3=2.*hci(j)*h1t(j)/ai(j)
      s6=0.5*psii(j)*psii(j)*phi1ni(j)*fh*(r1+r2+r3)
c
      s7=0.5*psii(j)*psii(j)*(f1n/fh)*(r1+r2+r3)
c
      r1=-b1ni(j)*hai(j)/(ai(j)*ai(j))
      r2=-bi(j)*ha1ni(j)/(ai(j)*ai(j))
      r3=2.*bi(j)*hai(j)*a1ni(j)/(ai(j)*ai(j)*ai(j))
      r4=-h2*hb1ni(j)/bi(j)
      r5=h2*hbi(j)*b1ni(j)/(bi(j)*bi(j))
      r6=2.*hc1ni(j)*h1t(j)/ai(j)
      r7=-2.*hci(j)*a1ni(j)*h1t(j)/(ai(j)*ai(j))
      s8=psii(j)*psii(j)*fh*(r1+r2+r3+r4+r5+r6+r7)
c
      r1=b1ni(j)*psii(j)*psii(j)/(ai(j)*ai(j))
      r2=0.5*psii(j)*psii(j)*phi1ni(j)*bi(j)/(ai(j)*ai(j))
      r3=-2.*psii(j)*psii(j)*bi(j)*a1ni(j)/(ai(j)**3.)
      s9=(fh**3.)*trk*(r1+r2+r3)
c
      s10=1.5*bi(j)*psii(j)*psii(j)*f1n*fh*trk/(ai(j)*ai(j))
c
      s11=bi(j)*psii(j)*psii(j)*(fh**3.)*trk1n/(ai(j)*ai(j))
c
      dd(j)=s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11
   10 continue
c
      j=1
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
      trk1n=ha1ni(j)/ai(j)-hai(j)*a1ni(j)/(ai(j)*ai(j))
     &     +hb1ni(j)/bi(j)-hbi(j)*b1ni(j)/(bi(j)*bi(j))
     &     +hd1ni(j)/di(j)-hdi(j)*d1ni(j)/(di(j)*di(j))
c
      r1=-4./(dtheta*dtheta)
      r2=-b2ni(j)/ai(j)
      r3=a1ni(j)*b1ni(j)/(ai(j)*ai(j))
      r4=-b1ni(j)*phi1ni(j)/ai(j)
      r5=-bi(j)*phi2ni(j)/ai(j)
      r6=a1ni(j)*bi(j)*phi1ni(j)/(ai(j)*ai(j))
      r7=-0.5*psii(j)*psii(j)*phi1ni(j)*bi(j)*hai(j)/(ai(j)**1.5)
      r8=-psii(j)*psii(j)*b1ni(j)*hai(j)/(ai(j)**1.5)
      r9=-psii(j)*psii(j)*bi(j)*ha1ni(j)/(ai(j)**1.5)
      r10=1.5*psii(j)*psii(j)*bi(j)*hai(j)*a1ni(j)/(ai(j)**2.5)
      r11=b1ni(j)*psii(j)*psii(j)*trk/sqrt(ai(j))
      r12=0.5*psii(j)*psii(j)*phi1ni(j)*bi(j)*trk/sqrt(ai(j))
      r13=psii(j)*psii(j)*bi(j)*trk1n/sqrt(ai(j))
      r14=-0.5*a1ni(j)*psii(j)*psii(j)*bi(j)*trk/(ai(j)**1.5)
c
      dd(j)=r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14
c
c        upper and lower diagonalz
      dtm1=1./dtheta
      do 20 j=2,nthtamax-1
      h2=h1t(j)*h1t(j)
      h3=h2*h1t(j)
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
c
      f=ai(j)*ai(j)*h2/bi(j)+ai(j)
      fh=sqrt(f)
c
      s1=1./(dtheta*dtheta)
c
      r1=ai(j)*ct(j)/bi(j)
      r2=0.5*ai(j)*d1ti(j)/(bi(j)*di(j))
      r3=0.5*a1ti(j)/bi(j)
      r4=phi1ti(j)*ai(j)/bi(j)
      s2=1.5*h2*dtm1*(r1+r2+r3+r4)
c
      r1=-0.5*d1ni(j)/di(j)
      r2=-b1ni(j)/bi(j)
      r3=0.5*a1ni(j)/ai(j)
      r4=-phi1ni(j)
      s3=h1t(j)*dtm1*(r1+r2+r3+r4)
c
      r1=ct(j)
      r2=0.5*d1ti(j)/di(j)
      r3=-0.5*b1ti(j)/bi(j)
      r4=a1ti(j)/ai(j)
      r5=phi1ti(j)
      s4=0.5*dtm1*(r1+r2+r3+r4+r5)
c
      r1=-bi(j)*hai(j)/(ai(j)*ai(j))
      r2=-h2*hbi(j)/bi(j)
      r3=2.*hci(j)*h1t(j)/ai(j)
      s5=(0.5*ai(j)*ai(j)*psii(j)*psii(j)*h1t(j)*dtm1/(bi(j)*fh))
     &  *(r1+r2+r3)
c
      r1=-hbi(j)*h1t(j)/bi(j)
      r2=hci(j)/ai(j)
      s6=psii(j)*psii(j)*fh*dtm1*(r1+r2)
c
      s7=1.5*psii(j)*psii(j)*fh*h1t(j)*dtm1*trk
c
      du(j)=s1+s2+s3+s4+s5+s6+s7
      dl(j)=s1-s2-s3-s4-s5-s6-s7
   20 continue
      du(1)=4./(dtheta*dtheta)
      dl(nthtamax)=2./(dtheta*dtheta)
c
c
      return
      end
c
c
c
      subroutine tridiag(dg,dgu,dgl,rhs,ans)
      dimension dg(*),dgu(*),dgl(*),rhs(*),ans(*)
c
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
       real gam(110)
c
c
c        this subroutine inverts a tridiagonal matric with
c        diagonal dg, upper diagonal dgu, lower diagonal dgl,
c        and right hand side vector rhs.  the result is returned
c        in array ans.
c
c
      bet=dg(1)
      ans(1)=rhs(1)/bet
c
      do 10 j=2,nthtamax
      gam(j)=dgu(j-1)/bet
      bet=dg(j)-dgl(j)*gam(j)
      if (bet.eq.0) write(5,*)'bet equals zero, bub',j
      ans(j)=(rhs(j)-dgl(j)*ans(j-1))/bet
   10 continue
c
      do 20 j=nthtamax-1,1,-1
      ans(j)=ans(j)-gam(j+1)*ans(j+1)
   20 continue
c
c
      return
      end
c
c
c
      subroutine diffh4
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
      real hx(-1:110)
c
c
      do 5 j=1,nthtamax
    5 hx(j)=h(j)
c
      hx(0)=h(2)
      hx(-1)=h(3)
      ntp2=nthtamax+2
      ntm2=nthtamax-2
      hx(ntp1)=h(ntm1)
      hx(ntp2)=h(ntm2)
c
      do 10 j=1,nthtamax
      jm1=j-1
      jp1=j+1
      jm2=j-2
      jp2=j+2
c
      h1t(j)=(hx(jm2)-hx(jp2)+8.*(hx(jp1)-hx(jm1)))*dt1o4
      h2t(j)=(-hx(jm2)-hx(jp2)+16.*(hx(jp1)+hx(jm1))
     &        -30.*hx(j))*dt2o4
   10 continue
c
c
      return
      end
c
c
c
      subroutine diffh2
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
      do 5 j=2,nthtamax-1
      h1t(j)=(h(j+1)-h(j-1))*dt1o2
      h2t(j)=(h(j+1)-2.*h(j)+h(j-1))*dt2o2
    5 continue
      h1t(1)=0.
      h1t(nthtamax)=0.
      h2t(1)=2.*(h(2)-h(1))*dt2o2
      h2t(nthtamax)=2.*(h(nthtamax-1)-h(nthtamax))*dt2o2
c
c
      return
      end
c
c
c
      subroutine residual
      include 'param.h'
      common /residblk/ res(110)
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
c
c     call interpl1
      call interpl3
c
      do 10 j=2,nthtamax-1
      h2=h1t(j)*h1t(j)
      h3=h1t(j)*h2
      f=ai(j)*ai(j)*h2/bi(j)+ai(j)
      fh=sqrt(f)
      fh3=fh*fh*fh
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
c
      s0=h2t(j)
c
      r1=ai(j)*ct(j)/bi(j)
      r2=0.5*ai(j)*d1ti(j)/(bi(j)*di(j))
      r3=0.5*a1ti(j)/bi(j)
      r4=phi1ti(j)*ai(j)/bi(j)
      s1=h3*(r1+r2+r3+r4)
c
      r1=-0.5*d1ni(j)/di(j)
      r2=-b1ni(j)/bi(j)
      r3=0.5*a1ni(j)/ai(j)
      r4=-phi1ni(j)
      s2=h2*(r1+r2+r3+r4)
c
      r1=ct(j)
      r2=0.5*d1ti(j)/di(j)
      r3=-0.5*b1ti(j)/bi(j)
      r4=a1ti(j)/ai(j)
      r5=phi1ti(j)
      s3=h1t(j)*(r1+r2+r3+r4+r5)
c
      r1=-0.5*bi(j)*d1ni(j)/(ai(j)*di(j))
      r2=-0.5*b1ni(j)/ai(j)
      r3=-bi(j)*phi1ni(j)/ai(j)
      s4=r1+r2+r3
c
      r1=-bi(j)*hai(j)/(ai(j)*ai(j))
      r2=-h2*hbi(j)/bi(j)
      r3=2.*hci(j)*h1t(j)/ai(j)
      s5=psii(j)*psii(j)*fh*(r1+r2+r3)
c
      s6=bi(j)*psii(j)*psii(j)*fh3*trk/(ai(j)*ai(j))
c
      res(j)=s0+s1+s2+s3+s4+s5+s6
   10 continue
c
      j=1
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
      r1=2.*h2t(j)
      r2=-b1ni(j)/ai(j)
      r3=-bi(j)*phi1ni(j)/ai(j)
      r4=-psii(j)*psii(j)*bi(j)*hai(j)/(ai(j)**1.5)
      r5=psii(j)*psii(j)*bi(j)*trk/sqrt(ai(j))
      res(j)=r1+r2+r3+r4+r5
c
      j=nthtamax
      trk=hai(j)/ai(j)+hbi(j)/bi(j)+hdi(j)/di(j)
      r1=h2t(j)
      r2=-0.5*bi(j)*d1ni(j)/(ai(j)*di(j))
      r3=-0.5*b1ni(j)/ai(j)
      r4=-bi(j)*phi1ni(j)/ai(j)
      r5=-psii(j)*psii(j)*bi(j)*hai(j)/(ai(j)**1.5)
      r6=psii(j)*psii(j)*bi(j)*trk/sqrt(ai(j))
      res(j)=r1+r2+r3+r4+r5+r6
c
c
      return
      end
c
c
c
      subroutine check(iflag)
      include 'param.h'
      common /residblk/ res(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
c        get residual
      call diffh2
      call residual
c
      rmax=-1.
c
      do 20 j=1,nthtamax
      r=abs(res(j))
      rmax=dmax1(r,rmax)
   20 continue
c
       epsilon=1.0e-8
c
       if (rmax.lt.epsilon) then
       iflag=1.
       else
       iflag=0.
       endif
c      write(5,*)'rmax  ',rmax
c
c
      return
      end
c
c
c
      subroutine ahtest1(rmax)
      include 'param.h'
      common /residblk/ res(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
c        get residual
      call diffh4
      call residual
c
      rmax=-1.
c
      do 20 j=1,nthtamax
      r=abs(res(j))
      rmax=dmax1(r,rmax)
   20 continue
c
c
      return
      end
c
c
c
      subroutine interpl1
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
      do 10 j=1,nthtamax
      sgn=sign(1.,h(j))
      hh=abs(h(j))
      i=int(hh/deta)+1
c
      if (i.gt.netamax-1)  i=netamax-1
c
      p=(hh-eta(i))/deta
      r=1.-p
c
      ip1=i+1
c
      ai(j)=r*a(i,j)+p*a(ip1,j)
      bi(j)=r*b(i,j)+p*b(ip1,j)
      di(j)=r*d(i,j)+p*d(ip1,j)
      psii(j)=r*psi(i,j)+p*psi(ip1,j)
c
      hai(j)=r*ha(i,j)+p*ha(ip1,j)
      hbi(j)=r*hb(i,j)+p*hb(ip1,j)
      hci(j)=sgn*(r*hc(i,j)+p*hc(ip1,j))
      hdi(j)=r*hd(i,j)+p*hd(ip1,j)
c
      a1ni(j)=sgn*(r*a1n(i,j)+p*a1n(ip1,j))
      a1ti(j)=r*a1t(i,j)+p*a1t(ip1,j)
      a2ni(j)=r*a2n(i,j)+p*a2n(ip1,j)
      a2ti(j)=r*a2t(i,j)+p*a2t(ip1,j)
      a2mi(j)=sgn*(r*a2m(i,j)+p*a2m(ip1,j))
c
      b1ni(j)=sgn*(r*b1n(i,j)+p*b1n(ip1,j))
      b1ti(j)=r*b1t(i,j)+p*b1t(ip1,j)
      b2ni(j)=r*b2n(i,j)+p*b2n(ip1,j)
      b2ti(j)=r*b2t(i,j)+p*b2t(ip1,j)
      b2mi(j)=sgn*(r*b2m(i,j)+p*b2m(ip1,j))
c
      d1ni(j)=sgn*(r*d1n(i,j)+p*d1n(ip1,j))
      d1ti(j)=r*d1t(i,j)+p*d1t(ip1,j)
      d2ni(j)=r*d2n(i,j)+p*d2n(ip1,j)
      d2ti(j)=r*d2t(i,j)+p*d2t(ip1,j)
      d2mi(j)=sgn*(r*d2m(i,j)+p*d2m(ip1,j))
c
      phi1ni(j)=sgn*(r*phi1n(i,j)+p*phi1n(ip1,j))
      phi1ti(j)=r*phi1t(i,j)+p*phi1t(ip1,j)
      phi2ni(j)=r*phi2n(i,j)+p*phi2n(ip1,j)
      phi2ti(j)=r*phi2t(i,j)+p*phi2t(ip1,j)
      phi2mi(j)=sgn*(r*phi2m(i,j)+p*phi2m(ip1,j))
c
      ha1ni(j)=sgn*(r*ha1n(i,j)+p*ha1n(ip1,j))
      ha1ti(j)=r*ha1t(i,j)+p*ha1t(ip1,j)
c
      hb1ni(j)=sgn*(r*hb1n(i,j)+p*hb1n(ip1,j))
      hb1ti(j)=r*hb1t(i,j)+p*hb1t(ip1,j)
c
      hc1ni(j)=r*hc1n(i,j)+p*hc1n(ip1,j)
      hc1ti(j)=sgn*(r*hc1t(i,j)+p*hc1t(ip1,j))
c
      hd1ni(j)=sgn*(r*hd1n(i,j)+p*hd1n(ip1,j))
      hd1ti(j)=r*hd1t(i,j)+p*hd1t(ip1,j)
   10 continue
c
c
      return
      end
c
c
c
      subroutine interpl3
      include 'param.h'
      common /metrcblk/ a(id,jd),b(id,jd),c(id,jd),d(id,jd)
      common /excrvblk/ ha(id,jd),hb(id,jd),hc(id,jd),hd(id,jd)
      common /cnfmlblk/ psim4(id,jd),psi(id,jd),phi1n(id,jd),
     &                  phi2n(id,jd),phi1t(id,jd),
     &                  phi2t(id,jd),phi2m(id,jd)
      common /gdervblk/ a1n(id,jd),b1n(id,jd),c1n(id,jd),
     &                  d1n(id,jd),a1t(id,jd),b1t(id,jd), 
     &                  c1t(id,jd),d1t(id,jd),a2n(id,jd), 
     &                  b2n(id,jd),c2n(id,jd),d2n(id,jd), 
     &                  a2t(id,jd),b2t(id,jd),c2t(id,jd), 
     &                  d2t(id,jd),a2m(id,jd),b2m(id,jd), 
     &                  c2m(id,jd),d2m(id,jd) 
      common /hdervblk/ ha1n(id,jd),ha1t(id,jd),hb1n(id,jd),
     &                  hb1t(id,jd),hc1n(id,jd),hc1t(id,jd),
     &                  hd1n(id,jd),hd1t(id,jd)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
c
      call spline3( a, ai,1., 1.,netamax,nthtamax)
      call spline3( b, bi,1., 1.,netamax,nthtamax)
      call spline3( d, di,1., 1.,netamax,nthtamax)
      call spline3(ha,hai,0., 1.,netamax,nthtamax)
      call spline3(hb,hbi,0., 1.,netamax,nthtamax)
      call spline3(hd,hdi,0., 1.,netamax,nthtamax)
      call spline3(hc,hci,0.,-1.,netamax,nthtamax)
c
      call spline3( a1n, a1ni,0.,-1.,netamax,nthtamax)
      call spline3( b1n, b1ni,0.,-1.,netamax,nthtamax)
      call spline3( d1n, d1ni,0.,-1.,netamax,nthtamax)
      call spline3(ha1n,ha1ni,0.,-1.,netamax,nthtamax)
      call spline3(hb1n,hb1ni,0.,-1.,netamax,nthtamax)
      call spline3(hd1n,hd1ni,0.,-1.,netamax,nthtamax)
      call spline3(hc1n,hc1ni,0., 1.,netamax,nthtamax)
c
      call spline3( a1t, a1ti,0., 1.,netamax,nthtamax)
      call spline3( b1t, b1ti,0., 1.,netamax,nthtamax)
      call spline3( d1t, d1ti,0., 1.,netamax,nthtamax)
      call spline3(ha1t,ha1ti,0., 1.,netamax,nthtamax)
      call spline3(hb1t,hb1ti,0., 1.,netamax,nthtamax)
      call spline3(hd1t,hd1ti,0., 1.,netamax,nthtamax)
      call spline3(hc1t,hc1ti,0.,-1.,netamax,nthtamax)
c
      call spline3(a2n,a2ni,0., 1.,netamax,nthtamax)
      call spline3(b2n,b2ni,0., 1.,netamax,nthtamax)
      call spline3(d2n,d2ni,0., 1.,netamax,nthtamax)
c
      call spline3(a2t,a2ti,0., 1.,netamax,nthtamax)
      call spline3(b2t,b2ti,0., 1.,netamax,nthtamax)
      call spline3(d2t,d2ti,0., 1.,netamax,nthtamax)
c
      call spline3(a2m,a2mi,0.,-1.,netamax,nthtamax)
      call spline3(b2m,b2mi,0.,-1.,netamax,nthtamax)
      call spline3(d2m,d2mi,0.,-1.,netamax,nthtamax)
c
      p1=psi(netamax,1)
      p2=phi1n(netamax,1)
      p3=phi1t(netamax,1)
      p4=phi2n(netamax,1)
      p5=phi2t(netamax,1)
      p6=phi2m(netamax,1)
      call spline3(psi  ,psii  ,p1, 1.,netamax,nthtamax)
      call spline3(phi1n,phi1ni,p2,-1.,netamax,nthtamax)
      call spline3(phi1t,phi1ti,p3, 1.,netamax,nthtamax)
      call spline3(phi2n,phi2ni,p4, 1.,netamax,nthtamax)
      call spline3(phi2t,phi2ti,p5, 1.,netamax,nthtamax)
      call spline3(phi2m,phi2mi,p6,-1.,netamax,nthtamax)
c
c
      return
      end
c
c
c
      subroutine spline3(u,ui,bnd,sym,ne,nt)
      include 'param.h'
      dimension u(ne,nt),ui(*)
c
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /dervablk/ ub(-1:idp2,-1:jdp2)
      common /diffblck/ dn1o2,dn2o2,dt1o2,dt2o2,dm2o2,
     &                  dn1o4,dn2o4,dt1o4,dt2o4,dm2o4
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
      real sbnd(110)
c
c
      do 1 j=1,nthtamax
    1 sbnd(j)=bnd
c
      call boundry2(u,sbnd,sym,netamax,nthtamax)
c
      do 10 j=1,nthtamax
      i=int(h(j)/deta)+1
      ip1=i+1
      im1=i-1
      ip2=i+2
c
c        form radial derivatives at i and i+1
      ub1ni=(ub(ip1,j)-ub(im1,j))*dn1o2
      ub1nip1=(ub(ip2,j)-ub(i,j))*dn1o2
c
c        form cubic spline coefficients at i 
      sp=abs((ub(ip1,j)-ub(i,j))/deta)
      sm=abs((ub(i,j)-ub(im1,j))/deta)
      rmin=3.*dmin1(sp,sm)
      fmin=dmin1(0.,ub1ni)
      fmax=dmax1(0.,ub1ni)
      filtmin=dmin1(fmax,rmin)
      filtmax=dmax1(fmin,-rmin)
      sigma=sign(1.,ub1ni)
      if (sigma.lt.0.0) then
      c2i=filtmax
      else
      c2i=filtmin
      endif
c     c2i=cvmgm(filtmax,filtmin,sigma)
c
      sp=abs((ub(ip2,j)-ub(ip1,j))/deta)
      sm=abs((ub(ip1,j)-ub(i,j))/deta)
      rmin=3.*dmin1(sp,sm)
      fmin=dmin1(0.,ub1nip1)
      fmax=dmax1(0.,ub1nip1)
      filtmin=dmin1(fmax,rmin)
      filtmax=dmax1(fmin,-rmin)
      sigma=sign(1.,ub1nip1)
      if (sigma.lt.0.0) then
      c2ip1=filtmax
      else
      c2ip1=filtmin
      endif
c     c2ip1=cvmgm(filtmax,filtmin,sigma)
c
      sp=(ub(ip1,j)-ub(i,j))/deta
      c3i=(3.*sp-c2ip1-2.*c2i)/deta
      c4i=(-2.*sp+c2ip1+c2i)/(deta*deta)
c
      d=h(j)-eta(i)
      d2=d*d
      d3=d2*d
c
      ui(j)=ub(i,j)+d*c2i+d2*c3i+d3*c4i
   10 continue
c
c
      return
      end
c
c
c
      subroutine wojtkiew
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
c
      call area(areahor)
      write(5,*)
      write(5,*)'horizon at axis,equator ',h(1),h(nthtamax)
      write(5,*)'area of apparent horizon ',areahor
c
      do 10 j=1,nthtamax
   10 h(j)=0.
      call area(areausm)
      write(5,*)'area of unstable surface ',areausm
      write(5,*)'difference and log of difference of areas'
      write(5,*) areausm-areahor,dlog10(abs(areausm-areahor))
c
      return
      end
c
c
c
      subroutine ahmass(ans)
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      call area(ar)
      ans=0.25*sqrt(ar/(2.*hfpi))
c
      return
      end
c
c
c
      subroutine area(ans)
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      real arg(110)
c
c
      call diffh4
c     call interpl1
      call interpl3
c
      do 10 j=1,nthtamax
      r1=(psii(j)**4.)*sqrt(bi(j)*di(j))     
      r2=sqrt(1.+(ai(j)/bi(j))*h1t(j)*h1t(j))
c
      arg(j)=r1*r2*sn(j)
   10 continue
c
      call simpson(arg,dtheta,nthtamax,ans)
      ans=4.*hfpi*2.*ans
c
c
      return
      end
c
c
c
      subroutine circumf(ans)
      include 'param.h'
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
c
      real arg(110)
c
c
      call diffh4
c     call interpl1
      call interpl3
c
c        equatorial circumference of trapped surface
      circe=4.*hfpi*psii(nthtamax)*psii(nthtamax)*sqrt(di(nthtamax))
c
c        polar circumference of trapped surface
      do 10 j=1,nthtamax
      arg(j)=psii(j)*psii(j)*sqrt(ai(j)*h1t(j)*h1t(j)+bi(j))
   10 continue
c
      call simpson(arg,dtheta,nthtamax,circ)
      circp=4.*circ
c
      ans=circp/circe
c
c
      return
      end
c
c
c
      subroutine embdtrp
      include 'param.h'
      common /emtrpblk/ gcah(110),curvet(110),arcleng(110)
      common /horizblk/ h(110),h1t(110),h2t(110),dh(110)
      common /intplblk/ ai(110),bi(110),di(110),hai(110),hbi(110),
     &                  hci(110),hdi(110),psii(110),a1ni(110),
     &                  b1ni(110),d1ni(110),a1ti(110),b1ti(110),
     &                  d1ti(110),a2ni(110),b2ni(110),d2ni(110),
     &                  a2ti(110),b2ti(110),d2ti(110),a2mi(110),
     &                  b2mi(110),d2mi(110),
     &                  phi1ni(110),phi1ti(110),phi2ni(110),
     &                  phi2ti(110),phi2mi(110),ha1ni(110),ha1ti(110),
     &                  hb1ni(110),hb1ti(110),hc1ni(110),hc1ti(110),
     &                  hd1ni(110),hd1ti(110)
      common /gridblck/ eta(500),theta(500),deta,dtheta,xi,
     &                  netamax,nthtamax,ngridpts,nem1,ntm1,nep1,ntp1
      common /constblk/ sn(500),acs(500),ct(500),sqrthfm,hfpi,sqrt2
      common /timeblck/ time,dt,scaleprm,adm,ntime,jtime,nf,fi,
     &                  nframes,nmovie
c
      real r1212(110),y(110),y1d(110)
c
c
c        r1212 has a factor of sin(theta)^2*psi^4 taken out
c
      call interpl3
      call diffh4
c
c        on axis
      j=1
c
      r1 =-2.*bi(j)*di(j)*di(j)*phi2ti(j)
      r2 =-2.*bi(j)*di(j)*di(j)*h2t(j)*phi1ni(j)
      r3 =4.*ai(j)*di(j)*di(j)*(h2t(j)**2.)
      r4 =-4.*bi(j)*di(j)*d1ni(j)*h2t(j)
      r5 =2.*b1ni(j)*di(j)*di(j)*h2t(j)
      r6 =-4.*bi(j)*di(j)*d2ti(j)
      r7 =2.*b2ti(j)*di(j)*di(j)
      r8 =-2.*bi(j)*di(j)*di(j)*phi2ti(j)
      r9 =-2.*bi(j)*di(j)*di(j)*h2t(j)*phi1ni(j)
      r10=-2.*bi(j)*di(j)*d1ni(j)*h2t(j)
      r11=-2.*bi(j)*di(j)*d2ti(j)
      r12=4.*bi(j)*di(j)*di(j)
c
      fac=1./(4.*di(j)*bi(j))
c
      r1212(j)=fac*(r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12)
c
c        off axis
      do 10 j=2,nthtamax
      r1 =-2.*ai(j)*di(j)*di(j)*(h1t(j)**2.)*phi1ti(j)*ct(j)
      r2 =-2.*bi(j)*di(j)*di(j)*phi1ti(j)*ct(j)
      r3 =-2.*ai(j)*di(j)*di(j)*(h1t(j)**3.)*phi1ni(j)*ct(j)
      r4 =-2.*bi(j)*di(j)*di(j)*h1t(j)*phi1ni(j)*ct(j)
      r5 =4.*ai(j)*di(j)*di(j)*h1t(j)*h2t(j)*ct(j)
      r6 =-4.*ai(j)*di(j)*d1ni(j)*(h1t(j)**3.)*ct(j)
      r7 =2.*a1ni(j)*di(j)*di(j)*(h1t(j)**3.)*ct(j)
      r8 =-4.*ai(j)*di(j)*d1ti(j)*(h1t(j)**2.)*ct(j)
      r9 =2.*a1ti(j)*di(j)*di(j)*(h1t(j)**2.)*ct(j)
c
      r10=-4.*bi(j)*di(j)*d1ni(j)*h1t(j)*ct(j)
      r11=2.*b1ni(j)*di(j)*di(j)*h1t(j)*ct(j)
      r12=-4.*bi(j)*di(j)*d1ti(j)*ct(j)
      r13=2.*b1ti(j)*di(j)*di(j)*ct(j)
      r14=-2.*ai(j)*di(j)*di(j)*(h1t(j)**2.)*phi2ti(j)
      r15=-2.*bi(j)*di(j)*di(j)*phi2ti(j)
      r16=2.*ai(j)*di(j)*di(j)*h1t(j)*h2t(j)*phi1ti(j)
      r17=-ai(j)*di(j)*d1ni(j)*(h1t(j)**3.)*phi1ti(j)
      r18=a1ni(j)*di(j)*di(j)*(h1t(j)**3.)*phi1ti(j)
      r19=-ai(j)*di(j)*d1ti(j)*(h1t(j)**2.)*phi1ti(j)
c
      r20=a1ti(j)*di(j)*di(j)*(h1t(j)**2.)*phi1ti(j)
      r21=-bi(j)*di(j)*d1ni(j)*h1t(j)*phi1ti(j)
      r22=b1ni(j)*di(j)*di(j)*h1t(j)*phi1ti(j)
      r23=-bi(j)*di(j)*d1ti(j)*phi1ti(j)
      r24=b1ti(j)*di(j)*di(j)*phi1ti(j)
      r25=-2.*ai(j)*di(j)*di(j)*(h1t(j)**4.)*phi2ni(j)
      r26=-2.*bi(j)*di(j)*di(j)*(h1t(j)**2.)*phi2ni(j)
      r27=-4.*ai(j)*di(j)*di(j)*(h1t(j)**3.)*phi2mi(j)
      r28=-4.*bi(j)*di(j)*di(j)*h1t(j)*phi2mi(j)
      r29=-2.*bi(j)*di(j)*di(j)*h2t(j)*phi1ni(j)
c
      r30=-ai(j)*di(j)*d1ni(j)*(h1t(j)**4.)*phi1ni(j)
      r31=a1ni(j)*di(j)*di(j)*(h1t(j)**4.)*phi1ni(j)
      r32=-ai(j)*di(j)*d1ti(j)*(h1t(j)**3.)*phi1ni(j)
      r33=a1ti(j)*di(j)*di(j)*(h1t(j)**3.)*phi1ni(j)
      r34=-bi(j)*di(j)*d1ni(j)*(h1t(j)**2.)*phi1ni(j)
      r35=b1ni(j)*di(j)*di(j)*(h1t(j)**2.)*phi1ni(j)
      r36=-bi(j)*di(j)*d1ti(j)*h1t(j)*phi1ni(j)
      r37=b1ti(j)*di(j)*di(j)*h1t(j)*phi1ni(j)
      r38=2.*ai(j)*di(j)*d1ti(j)*h1t(j)*h2t(j)
      r39=-2.*bi(j)*di(j)*d1ni(j)*h2t(j)
c
      r40=-2.*ai(j)*di(j)*d2ni(j)*(h1t(j)**4.)
      r41=ai(j)*d1ni(j)*d1ni(j)*(h1t(j)**4.)
      r42=a1ni(j)*di(j)*d1ni(j)*(h1t(j)**4.)
      r43=2.*ai(j)*d1ni(j)*d1ti(j)*(h1t(j)**3.)
      r44=a1ni(j)*di(j)*d1ti(j)*(h1t(j)**3.)
      r45=-4.*ai(j)*di(j)*d2mi(j)*(h1t(j)**3.)
      r46=a1ti(j)*di(j)*d1ni(j)*(h1t(j)**3.)
      r47=-2.*ai(j)*di(j)*d2ti(j)*(h1t(j)**2.)
      r48=ai(j)*d1ti(j)*d1ti(j)*(h1t(j)**2.)
      r49=a1ti(j)*di(j)*d1ti(j)*(h1t(j)**2.)
c
      r50=-2.*bi(j)*di(j)*d2ni(j)*(h1t(j)**2.)
      r51=bi(j)*d1ni(j)*d1ni(j)*(h1t(j)**2.)
      r52=b1ni(j)*di(j)*d1ni(j)*(h1t(j)**2.)
      r53=4.*ai(j)*di(j)*di(j)*(h1t(j)**2.)
      r54=2.*bi(j)*d1ni(j)*d1ti(j)*h1t(j)
      r55=b1ni(j)*di(j)*d1ti(j)*h1t(j)
      r56=-4.*bi(j)*di(j)*d2mi(j)*h1t(j)
      r57=b1ti(j)*di(j)*d1ni(j)*h1t(j)
      r58=-2.*bi(j)*di(j)*d2ti(j)
      r59=bi(j)*d1ti(j)*d1ti(j)
c
      r60=b1ti(j)*di(j)*d1ti(j)
      r61=4.*bi(j)*di(j)*di(j)
c
      fac=1./(4.*di(j)*(ai(j)*h1t(j)*h1t(j)+bi(j)))
c
      r1212(j)=fac*(      r1 +r2 +r3 +r4 +r5 +r6 +r7 +r8 +r9 
     &              +r10+r11+r12+r13+r14+r15+r16+r17+r18+r19 
     &              +r20+r21+r22+r23+r24+r25+r26+r27+r28+r29 
     &              +r30+r31+r32+r33+r34+r35+r36+r37+r38+r39 
     &              +r40+r41+r42+r43+r44+r45+r46+r47+r48+r49 
     &              +r50+r51+r52+r53+r54+r55+r56+r57+r58+r59 
     &              +r60+r61)
   10 continue
c
      do 20 j=1,nthtamax
      p4=psii(j)**4.
      fac=1./((ai(j)*h1t(j)*h1t(j)+bi(j))*di(j)*p4)
      gcah(j)=fac*r1212(j)
   20 continue
c
c        compute derivative of angular arclength
      do 30 j=1,nthtamax
   30 arcleng(j)=psii(j)*psii(j)*sqrt(ai(j)*h1t(j)*h1t(j)+bi(j))
c
c        construct normal curvature of arclength
c        on axis
      inoexist=0
c
      gausscrv=r1212(1)/(bi(1)*di(1))
      if (gausscrv.lt.0.0) then
      inoexist=1
      gausscrv=1.0e-6
      endif
      curvet(1)=sqrt(gausscrv)/(psii(1)*psii(1))
c
      do 40 j=2,nthtamax
      pxt=phi1ti(j)+h1t(j)*phi1ni(j)
      dxt=d1ti(j)+h1t(j)*d1ni(j)
      geodcurv=(0.5*pxt+0.5*dxt/di(j)+ct(j))
     &        /sqrt(ai(j)*h1t(j)*h1t(j)+bi(j))
c
      radical=1./(di(j)*sn(j)*sn(j))-geodcurv*geodcurv
      if (radical.lt.0.0) then
      inoexist=1
      radical=1.0e-6
      endif
      curvnorm=sqrt(radical)/(psii(j)*psii(j))
c
      gausscrv=r1212(j)*(psii(j)**(-4.))
     &        /((ai(j)*h1t(j)*h1t(j)+bi(j))*di(j))
c
      curvet(j)=gausscrv/curvnorm
   40 continue
c
      if (inoexist.eq.1) then
      write(5,*)'AH embedding does not exist at time ',dt*jtime/scaleprm
      endif
c
c
      return
      end
c
c
c
      subroutine alldone
c
c
      do 10 i=1,12
      if (i.ne.5) close(i)
   10 continue
c
      close(21)
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)
      close(51)
      close(52)
      close(53)
      close(54)
      close(61)
      close(62)
c
      write(5,*)
      write(5,*)'you are done, dudley-do-right'
      write(5,*)
c
      return
      end
c
c
c
      subroutine multigrid (idim,ilower,iupper,jdim,jlower,jupper,
     &             cc,cn,cs,cw,ce,cnw,cne,csw,cse,u,rhs,tol,
     &             ifmg,iflag)
************************************************************************
*
*
*     This routine is a wrapper for the multigrid solver.  The input
*  arrays are the finite difference coefficients on the 2D grid with
*  the boundary conditions absorbed into them.
*
****  Parameters
*
*  INTEGERS:
*     idim,jdim :
*          These define sizes of the coefficient arrays passed to this
*          routine.
*
*     ilower,iupper,
*     jlower,jupper :
*          These are the indices of the computational grid that
*          correspond to a fictional upper and lower index limits for
*          points on which computations are actually done.
*
*     ifmg :
*          ifmg = 0 - The full multigrid algorithm is not used to
*                     obtain a good initial guess on the fine grid.
*                     (use this if you can provide a good initial guess)
*          ifmg = 1 - The full multigrid algorithm is used to obtain a
*                     good initial guess on the fine grid.
*
*     type :
*          Defines the type of equation for the two black hole code
*              type = 1 -> lapse equation -> tol = 1.e-9
*              type = 2 -> shift equation -> tol = 1.e-10
*     Note: I have changed this so that tol is a passable parameter
*     D.B. 12 March 1992
*
*  REALS:
*     cc(idim,jdim),cn(idim,jdim),
*     cs(idim,jdim),cw(idim,jdim),ce(idim,jdim),
*     cnw(idim,jdim),cne(idim,jdim),
*     csw(idim,jdim),cse(idim,jdim) :
*          These are the finite difference coefficient arrays for a
*          nine-point stencil on the two dimensional grid as follows:
*
*
*               cnw    cn    cne
*            ^
*            |  cw     cc    ce
* increasing |
*  j values  |  csw    cs    cse
*   (theta)  |
*             --------> increasing i values (eta)
*
*  Ie. the coresspondence is : nw : i-1,j+1
*                              n  : i,j+1
*                              ne : i+1,j+1
*                              w  : i-1,j
*                              c  : i,j
*                              e  : i+1,j
*                              sw : i-1,j-1
*                              s  : i,j-1
*                              se : i+1,j-1
*
*     u :
*          Input: this contains the initial guess to the solution of the
*                 equation
*          Output: This contains the final approximation to the solution
*                 determined by the multigrid solver.
*     
*     rhs :
*          This array contains the values of the right hand side of the
*          equation at every point on the rwo dimensional grid.
*
************************************************************************
*
*

      integer idim,ilower,iupper,jdim,jlower,jupper,ifmg,type
      real cc(idim,jdim),cn(idim,jdim),cs(idim,jdim),cw(idim,jdim),
     &     ce(idim,jdim),cnw(idim,jdim),cne(idim,jdim),csw(idim,jdim),
     &     cse(idim,jdim),u(idim,jdim),rhs(idim,jdim)

*
************************************************************************
*
*  Variable definitions:
*
*  Integers:
*     ifd59 :
*          ifd59 = 5 - means a 5-point finite difference stencil
*                      (ac,an,as,aw,ae) is defined on the finest grid.
*          ifd59 = 9 - means a 9-point finite difference stencil
*                      (ac,an,as,aw,ae,anw,ane,asw,ase) is defined on
*                      the finest grid by the user.
*          (NOTE: This routine is coded to use ifd59=9 ONLY!)
*
*     ifmg :
*          ifmg = 0 - The full multigrid algorithm is not used to
*                     obtain a good initial guess on the fine grid.
*                     (use this if you can provide a good initial guess)
*          ifmg = 1 - The full multigrid algorithm is used to obtain a
*                     good initial guess on the fine grid.
*
*     ncyc :
*          The maximum number of multigrid "v"-cycles to be used. If
*          the maximum norm of the residual is not less than tol at the
*          end of ncyc cycles, the algorithm is terminated.
*          (NOTE: ncyc <= 40 )
*          
*     id5 :
*          Dimension of the arrays ac,aw,as,ae,an,q and f. id5 is the
*          total number of grid points on the finest grid and all
*          coarser grids.
*
*     id9 :
*          Dimension of the arrays asw,ase,ane,anw. If ifd59=5 then
*          id9=idi.  If ifd59=9 then id9=id5.
*
*     idi :
*          Dimension of the work arrays pu and pd. idi is the total
*          number of grid points on all of the coarser grids.
*
*     idg :
*          Dimension of the work array gam. It is set to the value im,
*          the number of grid points in the i-direction on the finest
*          grid.
*     np(20,8) :
*          Input: When the iskip=1,-1 or -2 option is used, np2 is
*                 assumed to contain the grid information for umgs2.
*          Output: When the iskip=0 option is used, the grid
*                 information for umgs2 is returned in np2.
*          (NOTE: This is only useful for multiple instance problems)
*
*     iskip :
*          iskip = 0 - The coarse grid information, coarse grid
*                      operators and interpolation coefficients are
*                      calculated by umgs2.  This information is stored
*                      in the arrays ac, aw, as, asw, ase, pu, pd, np2
*                      and the variable.
*          iskip = 1 - The calculation of the coarse grid information,
*                      coarse grid operators and interpolation
*                      coefficients is skipped.  This option would be
*                      used when umgs2 has been called with iskip=0 and
*                      is being called again to solve a system of
*                      equations with the same matrix. This would be
*                      the case in, say, parabolic problems with time
*                      independent coefficients.
*          iskip =-1 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and the given matrix
*                      equation is solved.  This option would be used
*                      when umgs2 has been called with iskip=0 and is
*                      being called again to solve a system of equations
*                      with a different matrix of the same dimensions.
*                      This would be the case for, say, parabolic
*                      problems with time dependent coefficients.
*          iskip =-2 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and returned.
*                      No matrix solve.
*
*     ipc :
*          ipc = 0 or 1.
*          ipc is a multigrid parameter which determines the type of
*          interpolation to be used.  Usually ipc=1 is best.  However,
*          if the boundary condition equations have been absorbed into
*          the interior equations then ipc=0 can be used which results
*          in a slightly more efficient algorithm.
*
*     nman :
*          nman = 0 usually.
*          nman =1 signals that the fine grid equations are singular for
*          the case when homogeneous Neumann boundary conditions are
*          applied along the entire boundary.  In this case, the
*          difference equations are singular and the condition that the
*          integral of q over the domain be zero is added to the set of
*          difference equations.  This condition is satisfied by adding
*          the appropriate constant vector to q on the fine grid.  It is
*          assumed, in this case, that a well-defined problem has been
*          given to mgss2, i.e. the integral of f over the domain is
*          zero.
*
*     im :
*          The number of grid points in the x-direction (including two
*          ficticious points)
*     jm :
*          The number of grid points in the y-direction (including two
*          ficticious points)
*
*     linp :
*          This is a dummy argument left over from the authors
*          development of the code
*             Use:  common /io/ linp,lout
*
*     lout :
*          lout = unit number of output file into which the maximum norm
*          of the residual after each multigrid v-cycle" is printed.
*             Use:  common /io/ linp,lout
*
*     iscale :
*          Flag to indicate whether problem can be diagonnaly scaled to
*          speed convergence of the multigrid solver.
*
*     REALS:
*
*     ac(id5),an(id5),as(id5),aw(id5),ae(id5),
*     anw(id9),ane(id9),asw(id9),ase(id9) :
*          Input: ac, an, as, aw, ae, anw, ane, asw and ase contain the
*                 stencil coefficients for the difference operator on
*                 the finest grid. When the iskip=1 option is used,
*                 these arrays also are assumed to contain the coarse
*                 grid difference stencil coeficients.
*          Output: when the iskip=0 option is used, the coarse grid
*                 stencil coeficients are returned in ac, an, as, aw,
*                 ae, anw, ane, asw and ase.
*
*     ru(idi),rd(idi),rc(idi) :
*          Real work arrays.
*
*     pu(idi),pd(idi),pc(idi) :
*          Real work arrays.
*          Input: when the iskip=1 option is used, these arrays are
*                 assumed to contain the interpolation coefficients used
*                 in the semi-coarsening multigrid algorithm.
*          Output: when the iskip=0 option is used, the interpolation
*                 coeficients are returned in pu and pd.
*
*     f(id5) :
*          f contains the right hand side vector of the matrix
*          equation to be solved by umgs2.
*
*     q(id5) :
*          If ifmg=0, q contains the initial guess on the fine grid.
*          If ifmg=1, the initial guess on the fine grid is determined
*                     by the full multigrid process and the value of
*                     q on input to umgs2 not used.
*
*     tol :
*          tol > 0.  The maximum norm of the residual is calculated at
*                    the end of each multigrid cycle. The algorithm is
*                    terminated when this maximum becomes less than tol
*                    or when the maximum number of iterations (see ncyc)
*                    is exceeded.  It is up to the user to provide a
*                    meaningfull tolerance criteria for the particular
*                    problem being solved.
*          tol = 0.  Perform ncyc multigrid cycles.  Calculate and print
*                    the maximum norm of the residual after each cycle.
*          tol =-1.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the final residual is calculated and returned in
*                    the variable rmax in the calling list of umgs2.
*          tol =-2.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the residual is not calculated.
*
*     rmax :
*          If tol.ge.-1., the final residual norm is returned in rmax.
*
************************************************************************
*
*

      integer ifd59,ncyc
      parameter (ifd59=9,ncyc=40)

      integer id5,id9,idi,idg
c  This is for a 103x28 grid
c     parameter(id5=6732,id9=6732,idi=3774,idg=102)
c  This is for a 203x56 grid
c     parameter(id5=24846,id9=24846,idi=13534,idg=202)
c  This is for a 403x118 grid
      parameter(id5=100098,id9=100098,idi=52662,idg=402)

      integer np2(20,8)
      integer iskip,ipc,nman
      parameter (iskip=0,ipc=1,nman=0)
      integer irc,irurd,im,jm
      integer linp,lout
      common /io/ linp,lout

      real ac(id5),an(id5),as(id5),aw(id5),ae(id5),
     &     anw(id9),ane(id9),asw(id9),ase(id9),
     &     q(id5),f(id5),gam(idg)

      real ru(idi),rd(idi),rc(idi),pu(idi),pd(idi),pc(idi)

      real tol,rmax

      integer iscale

*
* Set some parameters for multigrid solver
*
      lout=20
*     rewind(unit=lout)
      irc=0
      irurd=0
      im=iupper-ilower+3
      jm=jupper-jlower+3
*     if (type.eq.1) tol=1.e-9
*     if (type.eq.2) tol=1.e-10

*
*  Set up coefficients into vectors with correct indexing
*
      do 110 j=jlower,jupper
      do 100 i=ilower,iupper
      n=(j)*im + i+1-(ilower-1)
      ac(n)=cc(i,j)
      an(n)=cn(i,j)
      as(n)=cs(i,j)
      aw(n)=cw(i,j)
      ae(n)=ce(i,j)
      anw(n)=cnw(i,j)
      ane(n)=cne(i,j)
      asw(n)=csw(i,j)
      ase(n)=cse(i,j)
      q(n)=u(i,j)
      f(n)=rhs(i,j)
 100  continue
 110  continue

*
*  Determine whether we can diagonal scale the problem to speed
*  convergence. Can only be done if there are no zeros on the main
*  diagonal (ie. central difference coefficient).
*
      iscale=1
      do 200 j=jlower,jupper
      do 205 i=ilower,iupper
      n=(j)*im + i+1-(ilower-1)
      if (ac(n) .eq. 0.) then
        iscale=0
      endif
 205  continue
 200  continue

*
*  Do the diagonal scaling if we can.
*
      if (iscale.eq.1) then
      do 210 j=jlower,jupper
      do 215 i=ilower,iupper
        n=(j)*im + i+1-(ilower-1)
        f(n)=f(n)/ac(n)
        ae(n)=ae(n)/ac(n)
        aw(n)=aw(n)/ac(n)
        as(n)=as(n)/ac(n)
        an(n)=an(n)/ac(n)
        ase(n)=ase(n)/ac(n)
        asw(n)=asw(n)/ac(n)
        ane(n)=ane(n)/ac(n)
        anw(n)=anw(n)/ac(n)
        ac(n)=1.
 215  continue
 210  continue
      endif

c 
c   Now call the multigrid routine
      
*     write(5,*)'equation is type : ',type
      call umgs2 (
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,pu,pd,pc,ru,rd,rc,gam,np2,
     + ifd59,ifmg,ncyc,tol,nman,im,jm,id5,id9,idi,m,iskip,rmax,
     + ipc,irc,irurd)

      iflag=1
      if (rmax.gt.tol) then
*       write(5,*)'Did not converge type=',type
*       write(5,*)' maximum residual    = ',rmax
*       write(5,*)' tolerance           = ',tol
        iflag=0
*     stop
      endif

*     write(5,*)' maximum residual    = ',rmax
*     write(5,*)

*
*  Convert the solution back to the 2D array form
*
      do 510 j=jlower,jupper
      do 500 i=ilower,iupper
      n=(j)*im + i+1-(ilower-1)
      u(i,j)=q(n)
 500  continue
 510  continue

      return
      end
*
c----------------------------------------------------------------------
      subroutine ugrdfn(m,ifd59,is5,is9,isi,np2,imx,jmx)
c----------------------------------------------------------------------
c  Given imx, jmx and ifd59 (See comments in mgss2), ugrdfn calculates
c  the number of grids that will be needed.  Pointers into the arrays
c  ac, aw, as, asw, ase, q, f, pu, pd, pc, ru, rd and rc and the size
c  of each grid is calculated and stored in the array np2.  The
c  subroutine ukey is called to retrieve the grid information.
c .....................................................................
      parameter(n5=1,n9=2,ni=3,jm=4,i9=5,j9=6,ifd=7,jred=8)
      dimension np2(20,8)
      common /cs/ icorstr,iprint
      iq5=1
      iq9=1
      iqi=1
      m=1
      np2(m,1)=jmx
      np2(m,2)=3
   10 if(np2(m,1).le.3) go to 20
      m=m+1
      np2(m,1)=np2(m-1,1)/2+1
      if(np2(m-1,2).eq.2.and.mod(np2(m-1,1),2).eq.1)
     + np2(m,1)=np2(m,1)+1
      np2(m,2)=2
      go to 10
   20 do 30 k=1,m
      np2(m-k+1,jm)=np2(k,1)
   30 np2(m-k+1,jred)=np2(k,2)
      do 40 k=m,1,-1
      ktot=imx*np2(k,jm)
      np2(k,n5)=iq5
      iq5=iq5+ktot
      np2(k,n9)=iq9
      if(k.lt.m.or.ifd59.eq.9) iq9=iq9+ktot
      np2(k,ni)=iqi
   40 if(k.lt.m) iqi=iqi+ktot
      do 50 k=1,m
      np2(k,i9)=imx
      np2(k,j9)=np2(k,jm)
   50 np2(k,ifd)=9
         if(ifd59.eq.5) then
      np2(m,i9)=1
      np2(m,j9)=1
      np2(m,ifd)=5
         endif
      is5=iq5-1
      is9=iq9-1
      isi=iqi-1
      return
      end
c----------------------------------------------------------------------
      subroutine uintad(q,qc,pu,pd,im,jm,jmc,iadd,jred,ipc)
c----------------------------------------------------------------------
c  iadd=1:
c  Interpolates and adds the coarse grid (kf-1) correction, qc, to the
c  fine grid (kf) approximation, q, at the black y-lines.
c  iadd=0:
c  In the full multigrid algorithm, the solution to the coarse grid
c  (kf-1) difference equation is interpolated to the fine grid (kf)
c  to be used as the initial guess vector for kf=2,3,...,m.
c  Interpolation is at black y-lines only.
c .....................................................................
      dimension q(im,jm),qc(im,jmc),pu(im,jmc),pd(im,jmc)
      im1=im-1
      jm1=jm-1
      jblack=5-jred
c                                    add correction to next finer grid
 1000    if(iadd.eq.1) then
      jc=3-jred
      do 10 j=jblack,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 q(i,j)=q(i,j)+pd(i,jc)*qc(i,jc)+pu(i,jc)*qc(i,jc+1)
c
c                       interpolate solution to next finer grid in fmg
 1001    else
      jc=3-jred
      do 40 j=jblack,jm1,2
      jc=jc+1
      do 40 i=2,im1
   40 q(i,j)=pd(i,jc)*qc(i,jc)+pu(i,jc)*qc(i,jc+1)
 1002    endif
      return
      end
c----------------------------------------------------------------------
      subroutine ukey(k,np2,nn5,nn9,nni,jjm,ii9,jj9,iifd,jjred)
c----------------------------------------------------------------------
c  Returns the grid pointers and dimension variables for grid k.  The
c  information is stored in the array np2.
c......................................................................
      parameter(n5=1,n9=2,ni=3,jm=4,i9=5,j9=6,ifd=7,jred=8)
      dimension np2(20,8)
      nn5=np2(k,n5)
      nn9=np2(k,n9)
      nni=np2(k,ni)
      jjm=np2(k,jm)
      ii9=np2(k,i9)
      jj9=np2(k,j9)
      iifd=np2(k,ifd)
      jjred=np2(k,jred)
      return
      end
c----------------------------------------------------------------------
      subroutine umgs2(
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,pu,pd,pc,ru,rd,rc,gam,np2,
     + ifd59,ifmg,ncyc,tol,nman,im,jm,id5,id9,idi,m,iskip,rmax,
     + ipc,irc,irurd)
c----------------------------------------------------------------------
c** SUBROUTINE UMGS2
c
c** COPYRIGHT: Ecodynamics Research Associates, Inc.
c
c** Date written: June, 1990
c** Author:  Steve Schaffer
c            Mathematics Department
c            New Mexico Tech
c            Socorro, NM  87801
c            505-835-5811
c
c** DESCRIPTION:
c     umgs2 is a black box symmetric matrix solver.  It is written
c     in unsymmetric storage mode and can be used to solve mildly
c     nonsymmetric problems.  The user provides a matrix and right hand
c     side vector corresponding to a 5 or 9 point finite difference/
c     finite volume discretization of a symmetric second order PDE.
c     umgs2 will construct a sequence of coarse grids and coarse
c     grid operators and then solve the matrix equation using a
c     y-direction semi-coarsening multigrid algorithm and return
c     the solution vector.  If a sequence of matrix problems are
c     to be solved using the same matrix, computational time can
c     be saved by skipping the construction of the coarse grid
c     information in subsequent calls to umgs2.
c
c     The matrix on the finest grid is stored in the arrays ac,aw,as,
c     an,asw,ase,ane and anw.  The difference stencil at the point
c     (i,j) given by
c
c         nw  n  ne        anw(i,j)   an(i,j)   ane(i,j)
c          w  c  e    =    aw(i,j)    ac(i,j)   aw(i,j)
c         sw  s  se        asw(i,j)   as(i,j)   ase(i,j)
c
c     If the difference stencil on the fine grid is a 5 point stencil
c     then the arrays asw,ase,ane,anw are not used and the
c     stencil is given by
c
c             n                     an(i,j)
c          w  c  e    =    aw(i,j)  ac(i,j)   ae(i,j)
c             s                     as(i,j)
c
c     However, asw,ase,ane,anw still need to be dimensioned (by id9)
c     in the calling program as they are used in the coarse grid
c     calculations.
c
c** STORAGE
c     It is assumed that a set of ficticious points have been defined
c     along the entire boundary.  These points have nothing to do with
c     the solution and are used for programming convenience and
c     vectorization purposes.  Storage is allocated for the stencil
c     elements ac,aw,as,asw,ase, the solution vector, q, and the
c     right hand side vector, f, at these ficticious points.  The
c     stencils at these ficticious points and all stencil connections
c     to them are set to zero in the subroutine useta which is called
c     by umgs2.  The computational grid is depicted by
c
c         x   x   x   x              x   x   x   x
c
c         x   *   *   *              *   *   *   x
c
c         x   *   *   *              *   *   *   x
c    .
c    .
c    .
c         x   *   *   *              *   *   *   x
c
c         x   *   *   *              *   *   *   x
c
c         x   x   x   x              x   x   x   x
c
c     where x depicts the ficticious points and * depicts the interior
c     points.  The total storage requirements for the fine grid problem
c     is then 5*im*jm for 5 point stencils and 7*im*jm for 9 point
c     stencils.  The total storage requirements for the multigrid
c     solution is approximately 2 to 3 times that of the storage
c     requirements of the fine grid problem. (See DIMENSION PARAMETERS).
c     Note:  The first im*jm elements of the arrays ac,aw,as,[asw,ase],
c     q and f correspond to the finest grid.
c
c** DIMENSION PARAMETERS
c     The arrays ac,aw,ae,asw,ase,q,f,pu,pd and pc are dimensioned as on
c     dimensional arrays in the calling program.  They are dimensioned
c     as two dimensional arrays in the working subroutines.  The one
c     dimensional storage of the arrays, say q, follows: n=(j-1)*jm+i,
c     where n is the element location in the one dimensional storage of
c     q corresponding to the (i,j)th element of the two dimensional
c     storage of q and jm is the number of grid points in the j
c     direction (including the two ficticious points).
c
c     The dimension parameters are id5, id9, idi and idg.  They can be
c     determined by running the companion program MSS2DIM.F.
c      id5 - Integer variable.
c            Dimension of the arrays ac,aw,as,ae,an,q and f in the
c            calling program.  id5 is the total number of grid points
c            on the finest grid and all coarser grids.
c      id9 - Integer variable.
c            Dimension of the arrays asw,ase,ane,anw in the calling
c            program.  If ifd59=5 then id9=idi.  If ifd59=9 then
c            id9=id5.
c      idi - Integer variable.
c            Dimension of the work arrays pu and pd in the calling
c            program.  idi is the total number of grid points on all
c            of the coarser grids.
c      idg - Integer variable.
c            Dimension of the work array gam in the calling program.
c            It is set to the value im, the number of grid points
c            in the i-direction on the finest grid.
c
c** INPUT
c            (Note: all variable types are set implicitly)
c ac,aw,as
c    ae,an - Real arrays.  Dimensioned (id5) in calling program.
c            See comments in DESCRIPTION and DIMENSION PARAMETERS.
c  asw,ase
c  ane,anw - Real arrays.  Dimensioned (id9) in calling program.
c            See comments in DESCRIPTION and DIMENSION PARAMETERS.
c        f - Real array.  Dimensioned (id5) in calling program.
c            f contains the right hand side vector of the matrix
c            equation to be solved by umgs2.
c        q - Real array.  Dimensioned (id5) in calling program.
c            If ifmg=0, q contains the initial guess on the fine
c            grid.  If ifmg=1, the initial guess on the fine grid
c            is determined by the full multigrid process and the
c            value of q on input to umgs2 not used.
c   ifd59 -  Integer variable.
c            =5 - means a 5-point finite difference stencil (ac,aw and
c                 as) is defined on the finest grid by the user.
c            =9 - means a 9-point finite difference stencil (ac,aw,as,
c                 asw, ase) is defined on the finest grid by the user.
c    ifmg - Integer variable.
c           =0 - The full multigrid algorithm is not used to obtain a
c                good initial guess on the fine grid.
c           =1 - The full multigrid algorithm is used to obtain a good
c                initial guess on the fine grid.
c    ncyc - Integer variable.
c           The maximum number of multigrid "v"-cycles to be used.
c           If the maximum norm of the residual is not less than tol
c           at the end of ncyc cycles, the algorithm is terminated.
c     tol - Real variable.
c           >0 - The maximum norm of the residual is calculated at the
c           end of each multigrid cycle.  The algorithm is terminated
c           when this maximum becomes less than tol or when the maximum
c           number of iterations (see ncyc) is exceeded.  It is up to
c           the user to provide a meaningfull tolerance criteria for
c           the particular problem being solved.
c           =0 - Perform ncyc multigrid cycles.  Calculate and print
c           the maximum norm of the residual after each cycle.
c           =-1. - Perform ncyc multigrid cycles.  The maximum norm of
c           the final residual is calculated and returned in the
c           variable rmax in the calling list of umgs2.
c           =-2. - Perform ncyc multigrid cycles.  The maximum norm of
c           the residual is not calculated.
c   iskip - Integer variable.
c           =0 - The coarse grid information, coarse grid operators
c                and interpolation coefficients are calculated by
c                umgs2.  This information is stored in the arrays
c                ac, aw, as, asw, ase, pu, pd, np2 and the variable m
c                and returned to the calling program.
c           =1 - The calculation of the coarse grid information, coarse
c                grid operators and interpolation coefficients is
c                skipped.  This option would be used when umgs2 has
c                been called with iskip=0 and is being called again
c                to solve a system of equations with the same matrix.
c                This would be the case in, say, parabolic problems
c                with time independent coefficients.
c           =-1 -The set up of pointers (ugrdfn) is skipped.  Coarse gri
c                operators and interpolation coefficients are calculated
c                and the given matrix equation is solved.  This option
c                would be used when umgs2 has been called with iskip=0
c                and is being called again to solve a system of
c                equations with a different matrix of the same
c                dimensions.  This would be the case for, say,
c                parabolic problems with time dependent coefficients.
c           =-2 -The set up of pointers (ugrdfn) is skipped.  Coarse gri
c                operators and interpolation coefficients are calculated
c                and returned to the calling program.  No matrix solve.
c     ipc - Integer variable.
c           =0 or 1.
c           ipc is a multigrid parameter which determines the type of
c           interpolation to be used.  Usually ipc=1 is best.  However,
c           the boundary contition equations have been absorbed into the
c           interior equations then ipc=0 can be used which results in a
c           slightly more efficient algorithm.
c    nman - Integer variable.
c           =0 usually.
c           =1 signals that the fine grid equations are singular for
c           the case when homogeneous Neumann boundary conditions are
c           applied along the entire boundary.  In this case, the
c           difference equations are singular and the condition that
c           the integral of q over the domain be zero is added to the
c           set of difference equations.  This condition is satisfied
c           by adding the appropriate constant vector to q on the fine
c           grid.  It is assumed, in this case, that a well-defined
c           problem has been given to mgss2, i.e. the integral of f
c           over the domain is zero.
c      im - Integer variable.
c           The number of grid points in the x-direction (including two
c           ficticious points)
c      jm - Integer variable.
c           The number of grid points in the y-direction (including two
c           ficticious points)
c    lout - Integer variable.
c           = unit number of output file into which the maximum norm
c             of the residual after each multigrid v-cycle" is printed.
c             Use:  common /iout/ lout
c
c** INPUT/OUTPUT
c       q - Real array.  Dimensioned (id5)
c           On input, if ifmg=0, q contains the initial guess on the
c           finest grid for umgs2.  On output, q contains the final
c           solution on the finest grid.
c  ac-anw - Real arrays.  See DIMENSION.
c           On input, ac, aw, as, [asw and ase] contain the stencil
c           coefficients for the difference operator on the finest
c           grid.  When the iskip=1 option is used, these arrays
c           also are assumed to contain the coarse grid difference
c           stencil coeficients.
c           On output, when the iskip=0 option is used, the coarse
c           grid stencil coeficients are returned in ac - ase.
c
c ru,rd,rc - Real work arrays.  Dimensioned (idi)
c
c pu,pd,pc - Real work arrays.  Dimensioned (idi).
c           On input, when the iskip=1 option is used, these arrays
c           are assumed to contain the interpolation coefficients
c           used in the semi-coarsening multigrid algorithm.
c           On output, when the iskip=0 option is used, the
c           interpolation coeficients are returned in pu and pd.
c     np2 - Integer work array.  Dimensioned np2(20,8).
c           On input, when the iskip=1,-1 or -2 option is used, np2 is
c           assumed to contain the grid information for umgs2.
c           On output, when the iskip=0 option is used, the grid
c           information for umgs2 is returned in np2.
c** OUTPUT
c    rmax - If tol.ge.-1., the final residual norm is returned in rmax.
c
c** SUBROUTINES CALLED BY UMGS2
c
c    - ugrdfn, ukey, uintad, urelax, urscal, ursrhs, useta
c
c** END OF DESCRIPTION OF UMGS2
c .....................................................................
      dimension ac(id5),aw(id5),as(id5),ae(id5),an(id5),asw(id9),
     + ase(id9),ane(id9),anw(id9),q(id5),f(id5)
      dimension pu(idi),pd(idi),pc(idi),gam(im),np2(20,8)
      dimension ru(idi),rd(idi),rc(idi)
      dimension resid(0:40),confac(0:40)
      common /io/ linp,lout
c
c-time           tsu0=second()
         if(iskip.eq.0) then
      call ugrdfn(m,ifd59,is5,is9,isi,np2,im,jm)
         iquit=0
           if(m.gt.20) then
      iquit=1
      write(lout,*) ' m=',m,' > 20 - np2 is dimensioned np2(m=20,8)'
           endif
           if(is5.gt.id5) then
      iquit=1
      write(lout,*) ' id5=',id5,' too small.  Should be set to',is5
           endif
           if(is9.gt.id9) then
      iquit=1
      write(lout,*) ' id9=',id9,' too small.  Should be set to',is9
           endif
           if(isi.gt.idi) then
      iquit=1
      write(lout,*) ' idi=',idi,' too small.  Should be set to',isi
            endif
            if(is5.lt.2*im*jm) then
      iquit=1
      write(lout,*) ' id5.lt.2*im*jm can cause problems in useta'
      write(lout,*) ' this can be remedied by setting id5 larger'
            endif
      if(iquit.eq.1) return
         endif
         if(iskip.le.0) then
c     ----------  interpolation and coarse grid operators -----------
      do 5 k=m-1,1,-1
      call ukey(k+1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
      if(k.eq.m-1) n5cqf=n5c
    5 call useta(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),ac(n5c),aw(n5c),as(n5c),ae(n5c),an(n5c),asw(n9c),
     + ase(n9c),ane(n9c),anw(n9c),pu(nic),pd(nic),pc(nic),ru(nic),
     + rd(nic),rc(nic),q(n5cqf),f(n5cqf),gam,
     + im,jm,jmc,ifd,i9,j9,nman,k+1,m,jr,ipc,irc,irurd)
         endif
      if(iskip.eq.-2) return
c
         if(ifmg.ge.1) then
      do 6 k=m-1,1,-1
      call ukey(k+1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
    6 call ursrhs(f(n5),f(n5c),pu(nic),pd(nic),pc(nic),ru(nic),
     + rd(nic),rc(nic),im,jm,jmc,m,k+1,jr,irc)
         endif
c-time            tsu1=second()
c-time            write(lout,*) ' time for setup =',tsu1-tsu0
      l=1
      if(ifmg.eq.0) l=m
      k=l
      mcyc=0
      rmaxo=1.
c     ----------   begin multigrid cycling  ----------------------------
c
      if(l.eq.1) go to 20
   10 call ukey(k,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
      call ukey(k-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,tol,rmax,ipc,irc)
      if(k.eq.m.and.rmax.lt.tol) go to 60
         if(k.eq.m.and.tol.ge.-.5) then
      if(rmaxo.ne.0.) rate=rmax/rmaxo
      rmaxo=rmax
      if(mcyc.eq.0) rmax0=rmax
      resid(mcyc)=rmax
      confac(mcyc)=rate
         endif
      if(tol.eq.-.5) write(lout,*) ' down ',k,rmax
      k=k-1
      if(k.gt.1) go to 10
c     ---------  solve coarsest grid  ----------------------------------
c
   20 call ukey(1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
         if(l.eq.1) go to 40
c     ----------  interpolate correction to next finer grid  -----------
c
   30 k=k+1
      call ukey(k,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
      call uintad(
     + q(n5),q(n5c),pu(nic),pd(nic),im,jm,jmc,1,jr,ipc)
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
         if(tol.eq.-.5) then
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,tol,rmax,ipc,irc)
      write(lout,*) ' up   ',k,rmax
         endif
      if(k.lt.l) go to 30
      if(l.eq.m) go to 50
c     ----------  interpolate solution to new finest grid l+1 in fmg  --
c
   40 l=l+1
      k=l
      call ukey(l,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(l-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
      call uintad(
     + q(n5),q(n5c),pu(nic),pd(nic),im,jm,jmc,0,jr,0)
      go to 10
c
   50 if(nman.eq.1) call uneuman(q(n5),im,jm)
      mcyc=mcyc+1
c     ----------  Cycle ncyc times on grid m  --------------------------
      if(mcyc.lt.ncyc) go to 10
c-time            tmg1=second()
c-time            write(lout,*) ' time in ',ncyc,' cycles =',tmg1-tsu1
c
c     ----------  print out final residual and work units  -------------
        if(tol.ge.-1.) then
      call ukey(m,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(m-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,1.,rmax,ipc,irc)
      resid(mcyc)=rmax
      confac(mcyc)=rmax/rmaxo
      nb=0
      ne=min0(6,mcyc)
 2029 write(lout,2033) (mc,mc=nb,ne)
      write(lout,2032) (resid(mc),mc=nb,ne)
      write(lout,2031) (confac(mc),mc=nb,ne)
      nb=ne+1
      ne=ne+min0(6,mcyc-ne)
      if(nb.le.ne) go to 2029
      fconfac=(rmax/rmax0)**(1./float(mcyc))
      write(lout,2034) fconfac
 2034 format(30x,6(1h*)/,' average convergence factor =',f7.3,/,
     + 30x,6(1h*))
 2033 format(7(4x,i2,4x))
 2031 format(7(1x,f9.3))
 2032 format(7(1x,e9.3))
         endif
      return
   60 write(lout,1003) mcyc,resid(mcyc),tol
      return
 1003 format(' cyc=',i2,' max(res)=',1pe8.2/
     + ' tolerance condition tol=',1pe8.2,' satisfied')
      end
c----------------------------------------------------------------------
      subroutine uneuman(q,im,jm)
c----------------------------------------------------------------------
c  For problems with homogeneous Neumann boundary contitions, the
c  condition that the integral of q over the domain be zero is added
c  to the set of difference equations in order to obtain a unique
c  solution.
c......................................................................
      dimension q(im,jm)
      im1=im-1
      jm1=jm-1
      con=0.
      do 10 j=2,jm1
      do 10 i=2,im1
   10 con=con+q(i,j)
      con=con/((im-2)*(jm-2))
      do 20 j=2,jm1
      do 20 i=2,im1
   20 q(i,j)=q(i,j)-con
      return
      end
c**********************************************************************
      subroutine uoutpt(q,im,jm)
c**********************************************************************
c  Sample output subroutine.  Prints out the values of q at the
c  interior points of the finest grid.
c**********************************************************************
      common /io/ linp,lout
      dimension q(im,jm)
      im1=im-1
      jm1=jm-1
      ie=1
 20   ib=ie+1
      ie=ib+min0(5,im1-ib)
      do 10 j=jm1,2,-1
   10 write(lout,100) j,(q(i,j),i=ib,ie)
      if(ie.lt.im1) go to 20
  100 format(1x,i2,1x,6(1x,f10.4))
      return
      end
c**********************************************************************
      subroutine uputf(ac,aw,as,ae,an,f,nx,ny,
     + lo,nxd,nyd,i32su)
c**********************************************************************
      dimension ac(lo:nxd,lo:nyd),aw(lo:nxd,lo:nyd),
     + as(lo:nxd,lo:nyd),ae(lo:nxd,lo:nyd),
     + an(lo:nxd,lo:nyd),f(lo:nxd,lo:nyd)
      dimension a(20),b(20),ab(20),il(20),ir(20),jb(20),jt(20)
      dimension ibc(4)
      common /io/ linp,lout
c
c dcell is the value assigned to the diagonal element of a dead
c cell, i.e. a cell that has 0 conections to all its neighbors.
c icendif determines the differencing scheme for the first order terms
c icendif=0 - central differencing, =1 - forward differencing.
c
      dcell = 1.
      do 321 j=lo,nyd
      do 321 i=lo,nxd
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
 321  an(i,j)=0.
c
      nx1=nx+1
      ny1=ny+1
      read(linp,*) icendif
      read(linp,*) ibc(1),ibc(2),ibc(3),ibc(4)
      read(linp,*) nreg
      hx=1./nx
      hy=1./ny
      hx2=hx*hx
      hy2=hy*hy
      hxy2=hx2*hy2
      dcell=hxy2*dcell
      write(lout,1011) ibc(1),ibc(2),ibc(3),ibc(4),icendif,dcell,hx,hy
 1011 format(' ibc_l ibc_b ibc_r ibc_t   icendiff   dcell   hx      hy'/
     + 4x,i1,5x,i1,5x,i1,5x,i1,7x,i1,6x,e8.2,2x,f6.4,2x,f6.4/)
c
      do 10 irg=1,nreg
      read(linp,*) il(irg),ir(irg),jb(irg),jt(irg)
      read(linp,*) xk,yk,sreg,freg
      read(linp,*) a(irg),b(irg),ab(irg)
      write(lout,1000) il(irg),ir(irg),jb(irg),jt(irg),xk,yk,sreg,freg
 1000 format(1x,i3,',',i3,' X ',i3,',',i3,2x,1pe8.2,1x,1pe8.2,2x,
     + 1pe8.1,1x,1pe8.1)
      write(lout,1001) a(irg),b(irg),ab(irg)
 1001 format(17x,' a=',1pe10.3,'  b=',1pe10.3,' ab=',1pe10.3)
      xk=xk*hy2
      yk=yk*hx2
      sreg=sreg*hxy2
      freg=freg*hxy2
      a(irg)=a(irg)*hx2*hy
      b(irg)=b(irg)*hx*hy2
      ab(irg)=ab(irg)*hx*hy
         if(icendif.eq.0) then
      a(irg)=a(irg)/2.
      b(irg)=b(irg)/2.
      ab(irg)=ab(irg)/4.
         endif
      if(il(irg).eq.1) il(irg)=0
      if(ir(irg).eq.nx) ir(irg)=nx1
      if(jb(irg).eq.1) jb(irg)=0
      if(jt(irg).eq.ny) jt(irg)=ny1
      do 20 i=il(irg),ir(irg)
      do 20 j=jb(irg),jt(irg)
      aw(i,j)=xk
      as(i,j)=yk
      ac(i,j)=sreg
   20 f(i,j)=freg
   10 continue
      write(lout,*) ' - - - - - - - - - - - - - - - - - - - - - - - - -'
c                             defining coeficients by harmonic averaging
      do 30 i=1,nx
      asio=as(i,0)
      do 30 j=1,ny1
      aa=as(i,j)*asio
         if(aa.gt.0.) then
      t=2.*aa/(as(i,j)+asio)
      asio=as(i,j)
      as(i,j)=t
         else
      asio=as(i,j)
      as(i,j)=0.
         endif
   30 continue
      do 40 j=1,ny
      awoj=aw(0,j)
      do 40 i=1,nx1
      aa=aw(i,j)*awoj
         if(aa.gt.0.) then
      t=2.*aa/(aw(i,j)+awoj)
      awoj=aw(i,j)
      aw(i,j)=t
         else
      awoj=aw(i,j)
      aw(i,j)=0.
         endif
   40 continue
      do 45 i=0,nx
      do 45 j=0,ny
      ae(i,j)=aw(i+1,j)
   45 an(i,j)=as(i,j+1)
      do 50 i=1,nx
      do 50 j=1,ny
      ac(i,j)=ac(i,j)-aw(i,j)-as(i,j)-ae(i,j)-
     + an(i,j)
   50 if(ac(i,j).eq.0.) ac(i,j)=dcell
c                                        adding on the unsymmetric terms
      do 51 irg=1,nreg
c                                           icendif=0 ==> central diff'g
         if(icendif.eq.0) then
      do 52 i=il(irg),ir(irg)
      do 52 j=jb(irg),jt(irg)
      aw(i,j)=aw(i,j)-a(irg)
      ae(i,j)=ae(i,j)+a(irg)
      an(i,j)=an(i,j)+b(irg)
   52 as(i,j)=as(i,j)-b(irg)
c                                        icendif=1 ==> upstream diff's
         elseif(icendif.eq.1) then
      do 54 i=il(irg),ir(irg)
      do 54 j=jb(irg),jt(irg)
      ac(i,j)=ac(i,j)-a(irg)
      ae(i,j)=ae(i,j)+a(irg)
      an(i,j)=an(i,j)+b(irg)
   54 ac(i,j)=ac(i,j)-b(irg)
         endif
   51 continue
c                           set boundary conditions for 5 point operator
      do 60 j=1,ny
c                                                          left boundary
      ae(0,j)=aw(1,j)
         if(ibc(1).eq.1) then
      ac(0,j)=aw(1,j)
      f(0,j)=2.*aw(1,j)*0.0
         elseif(ibc(1).eq.2) then
      ac(0,j)=-aw(1,j)
      f(0,j)=hx*aw(1,j)*0.0
         endif
         if(ac(0,j).eq.0.) then
      ae(0,j)=0.
      ac(0,j)=dcell
      f(0,j)=0.
         endif
c                                                         right boundary
      aw(nx1,j)=ae(nx,j)
         if(ibc(3).eq.1) then
      ac(nx1,j)=ae(nx,j)
      f(nx1,j)=2.*ae(nx,j)*0.
         elseif(ibc(3).eq.2) then
      ac(nx1,j)=-ae(nx,j)
      f(nx1,j)=hx*ae(nx,j)*0.
         endif
         if(ac(nx1,j).eq.0.) then
      aw(nx1,j)=0.
      ac(nx1,j)=dcell
      f(nx1,j)=0.
         endif
   60 continue
c
      do 80 i=1,nx
c                                                         lower boundary
      an(i,0)=as(i,1)
         if(ibc(2).eq.1) then
      ac(i,0)=as(i,1)
      f(i,0)=2.*as(i,1)*0.
         elseif(ibc(2).eq.2) then
      ac(i,0)=-as(i,1)
      f(i,0)=hy*as(i,1)*0.
         endif
         if(ac(i,0).eq.0.) then
      an(i,0)=0.
      ac(i,0)=dcell
      f(i,0)=0.
         endif
c                                                         upper boundary
      as(i,ny1)=an(i,ny)
         if(ibc(4).eq.1) then
      ac(i,ny1)=an(i,ny)
      f(i,ny1)=2.*an(i,ny)*0.
         elseif(ibc(4).eq.2) then
      ac(i,ny1)=-an(i,ny)
      f(i,ny1)=2.*an(i,ny)*0.
         endif
         if(ac(i,ny1).eq.0.) then
      as(i,ny1)=0.
      ac(i,ny1)=dcell
      f(i,ny1)=0.
         endif
   80 continue
c                     connections between "ghost" boundary points zeroed
      do 83 j=1,ny1
      as(0,j)=0.
      an(0,j-1)=0.
      as(nx1,j)=0.
   83 an(nx1,j-1)=0.
      do 86 i=1,nx1
      aw(i,0)=0.
      ae(i-1,0)=0.
      aw(i,ny1)=0.
   86 ae(i-1,ny1)=0.
c                                        corner stencils and rhs defined
         if(i32su.eq.32) then
      do 90 j=0,ny1,ny1
      do 90 i=0,nx1,nx1
      ac(i,j)=dcell
      aw(i,j)=0.
      ae(i,j)=0.
      as(i,j)=0.
      an(i,j)=0.
   90 f(i,j)=0.
         endif
c                                i32su=22 - boundary conditions absorbed
         if(i32su.eq.22) then
      do 100 j=1,ny
      awac=aw(1,j)/ac(0,j)
      ac(1,j)=ac(1,j)-awac*ae(0,j)
      aw(1,j)=0.
      f(1,j)=f(1,j)-awac*f(0,j)
      ac(0,j)=0.
      ae(0,j)=0.
      f(0,j)=0
      awac=aw(nx1,j)/ac(nx1,j)
      ac(nx,j)=ac(nx,j)-awac*ae(nx1,j)
      ae(nx,j)=0.
      f(nx,j)=f(nx,j)-awac*f(nx1,j)
      ac(nx1,j)=0.
      aw(nx1,j)=0.
  100 f(nx1,j)=0.
c
      do 110 i=1,nx
      asac=as(i,1)/ac(i,0)
      ac(i,1)=ac(i,1)-asac*an(i,0)
      as(i,1)=0.
      f(i,1)=f(i,1)-asac*f(i,0)
      ac(i,0)=0.
      an(i,0)=0.
      f(i,0)=0.
      anac=an(i,ny)/ac(i,ny1)
      ac(i,ny)=ac(i,ny)-anac*as(i,ny1)
      an(i,ny)=0.
      f(i,ny)=f(i,ny)-anac*f(i,ny1)
      ac(i,ny1)=0.
      as(i,ny1)=0.
  110 f(i,ny1)=0.
         endif
      return
      end
c----------------------------------------------------------------------
      subroutine urelax(ac,aw,as,ae,an,asw,ase,ane,anw,f,q,gam,
     + im,jm,i9,j9,ifd,nman,k,m,jred,ipc,iprcud)
c----------------------------------------------------------------------
c  Performs red/black x-line relaxation.  The Thomas algorithm is used
c  to solve the tridiagonal matrices.
c** INPUT -
c       ac-anw=  finite difference operator coeficients
c            q=  initial approximation
c            f=  right hand side vector
c        im,jm=  the number of grid points in the x,y-directions
c        i9,j9=  the i,j-dimensions of the arrays asw,ase
c          ifd=  5 or 9 - the size of the stencil
c         nman-  =0 usually.
c                =1 signals that the fine grid equations are singular
c                for the case when Neumann boundary conditions are
c                applied along the entire boundary.  In this case, the
c                equations on the coarsest grid (consisting of a single
c                line of unknowns) is a singular tridiagonal system
c                and the Thomas algorithm is modified on this grid to
c                obtain a solution with an arbitrary constant vector
c                component.  This constant vector is removed on the
c                finest grid by the call to subroutine uneuman.
c** OUTPUT -
c            q=  final approximation after a red/black relaxation sweep
c .....................................................................
      dimension ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9)
      dimension f(im,jm),q(im,jm),gam(im)
      jm1=jm-1
      im1=im-1
      im2=im-2
      jblack=5-jred
c                                              usual red/black relaxatio
      nrel=2
      jrb=jred
c                                             ipc ..brbr relaxation swee
 1000    if(iprcud.eq.1) then
      nrel=ipc
      if(mod(ipc,2).eq.0) jrb=jblack
c                                     1 black relax for calc'g pu,pd,ru,
 1001    elseif(iprcud.eq.2) then
      nrel=1
      jrb=jblack
 1002    endif
c
c
      do 109 nrr=1,nrel
 5000    if(jrb.eq.jblack) then
c                                                             black rela
 6000    if(jblack.le.jm1) then
 1400    if(iprcud.ne.2) then
c
      do 110 j=jblack,jm1,2
      do 110 i=2,im1
  110 q(i,j)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)
 7000    if(ifd.eq.9) then
      do 120 j=jblack,jm1,2
      do 120 i=2,im1
  120 q(i,j)=q(i,j)-asw(i,j)*q(i-1,j-1)-ase(i,j)*q(i+1,j-1)-
     + anw(i,j)*q(i-1,j+1)-ane(i,j)*q(i+1,j+1)
 7001    endif
 1401    endif
c                                                black tridiagonal solve
c**
c**  Moved calculation of loop 129 from loop 130 for vectorization
c**  on vector machines (ie. Cray)
c**  By: John Towns 2/6/92
c**
      do 129 j=jblack,jm1,2        
      if (abs(ac(2,j)).lt.1.0e-65) then
        write(5,*)' central coef(2,',j+1,') = ',ac(2,j),q(2,j)
        stop
      endif
  129 q(2,j)=q(2,j)/ac(2,j)

c**
c**  Changed bet=(quantity) to bet=1./(quantity) to trade two divisions
c**  for one division and two multiplies (more efficient on all
c**  machines)
c**  By: John Towns 2/6/92
c**
      do 130 j=jblack,jm1,2
      bet=1./ac(2,j)
      do 140 i=3,im1
      gam(i)=ae(i-1,j)*bet
      bet=1./(ac(i,j)-aw(i,j)*gam(i))
  140 q(i,j)=(q(i,j)-aw(i,j)*q(i-1,j))*bet
      do 150 i=im2,2,-1
  150 q(i,j)=q(i,j)-gam(i+1)*q(i+1,j)
  130 continue
 6001    endif
c                                                             red relax
 5001    else
c
      do 210 j=jred,jm1,2
      do 210 i=2,im1
  210 q(i,j)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)
 1100    if(ifd.eq.9) then
      do 220 j=jred,jm1,2
      do 220 i=2,im1
  220 q(i,j)=q(i,j)-asw(i,j)*q(i-1,j-1)-ase(i,j)*q(i+1,j-1)-
     + anw(i,j)*q(i-1,j+1)-ane(i,j)*q(i+1,j+1)
 1101    endif
c                                                      tridiagonal solve
c                          nman=1 ==> avoid singularity on coarsest grid
      imm=im1
          if(nman.eq.1.and.k.eq.1) then
      imm=im-2
      q(im1,2)=0.
      gam(im1)=0.
          endif
c
c**
c**  Moved calculation of loop 229 from loop 230 for vectorization
c**  on vector machines (ie. Cray)
c**  By: John Towns 2/6/92
c**
      do 229 j=jred,jm1,2
      if (abs(ac(2,j)).lt.1.e-65) then
        write(5,*)' central coef(2,',j+1,') = ',ac(2,j),q(2,j)
        stop
      endif
  229 q(2,j)=q(2,j)/ac(2,j)

c**
c**  Changed bet=(quantity) to bet=1./(quantity) to trade two divisions
c**  for one division and two multiplies (more efficient on all
c**  machines)
c**  By: John Towns 2/6/92
c**
      do 230 j=jred,jm1,2
      bet=1./ac(2,j)
      do 240 i=3,imm
      gam(i)=ae(i-1,j)*bet
      bet=1./(ac(i,j)-aw(i,j)*gam(i))
  240 q(i,j)=(q(i,j)-aw(i,j)*q(i-1,j))*bet
      do 250 i=im2,2,-1
  250 q(i,j)=q(i,j)-gam(i+1)*q(i+1,j)
  230 continue
 5002    endif
      jrb=5-jrb
  109 continue
      return
      end
c----------------------------------------------------------------------
      subroutine urscal(
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,fc,qc,rc,
     + im,jm,jmc,ifd,i9,j9,kf,m,jred,tol,rmax,ipc,irc)
c----------------------------------------------------------------------
c  Defines the grid kf-1 right hand side, fc, as the restriction of the
c  grid kf residual.  The restriction operator is the transpose of the
c  interpolation operator.  Note:  The grid kf residual is zero at the
c  black lines (j-direction) as a result of red/black relaxation.
c  Thus, the restriction is simple injection.  The initial guess, qc,
c  for the coarse grid correction equation is set to zero.  The
c  maximum norm of the residual is calculated and returned in rmax.
c......................................................................
      dimension ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9)
      dimension f(im,jm),q(im,jm),fc(im,jmc),qc(im,jmc)
      dimension rc(im,jmc)
      rmax=0.
      im1=im-1
      jm1=jm-1
      jmc1=jmc-1
      jc=1
      do 10 j=jred,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 fc(i,jc)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)-
     + aw(i,j)*q(i-1,j)-ae(i,j)*q(i+1,j)-ac(i,j)*q(i,j)
 1000    if(ifd.eq.9) then
      jc=1
      do 20 j=jred,jm1,2
      jc=jc+1
      do 20 i=2,im1
   20 fc(i,jc)=fc(i,jc)-asw(i,j)*q(i-1,j-1)-ane(i,j)*q(i+1,j+1)-
     + ase(i,j)*q(i+1,j-1)-anw(i,j)*q(i-1,j+1)
 1001    endif
c                                           zero out qc as initial guess
      do 25 jc=1,jmc
      do 25 i=1,im
   25 qc(i,jc)=0.
c                                        if kf=m calculate residual norm
 2000    if((kf.eq.m.and.tol.ge.0.).or.tol.eq.-.5) then
      do 30 jc=2,jmc1
      do 30 i=2,im1
      resmax=abs(fc(i,jc))
   30 if(resmax.gt.rmax) rmax=resmax
 2001    endif
c                                                 weight rhs if irc.ge.1
 3000    if(irc.eq.1.and.ipc.ge.1) then
      do 40 jc=2,jmc1
      do 40 i=2,im1
   40 fc(i,jc)=rc(i,jc)*fc(i,jc)
 3001    endif
c
      return
      end
c----------------------------------------------------------------------
      subroutine ursrhs(f,fc,ru,rd,rc,im,jm,jmc,m,kf,jred,irc)
c----------------------------------------------------------------------
c  Restricts the right hand side vector on grid kf onto grid kf-1 when
c  the full multigrid (ifmg>0) option is used.  The restriction operator
c  is NOT necessarily the transpose of the interpolation operator.
c......................................................................
      dimension f(im,jm),fc(im,jmc),ru(im,jmc),rd(im,jmc),rc(im,jmc)
      jm1=jm-1
      im1=im-1
      jc=1
 1000    if(irc.eq.0) then
      do 10 j=jred,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 fc(i,jc)=ru(i,jc-1)*f(i,j-1)+rd(i,jc)*f(i,j+1)+f(i,j)
 1001    else
      do 20 j=jred,jm1,2
      jc=jc+1
      do 20 i=2,im1
   20 fc(i,jc)=ru(i,jc-1)*f(i,j-1)+rd(i,jc)*f(i,j+1)+
     + rc(i,jc)*f(i,j)
 1002    endif
      return
      end
c----------------------------------------------------------------------
      subroutine useta(
     + ac,aw,as,ae,an,asw,ase,ane,anw,acc,awc,asc,aec,
     + anc,aswc,asec,anec,anwc,pu,pd,pc,ru,rd,rc,qw,fw,gam,
     + im,jm,jmc,ifd,i9,j9,nman,kf,m,jred,ipc,irc,irurd)
c----------------------------------------------------------------------
c     Calculates the interpolation coefficients from grid kf-1 to
c     grid kf and the coarse grid operator on grid kf-1.
c** INPUT -
c    ac - anw = fine grid (kf) array stencil coeficients
c            m=  total number of grids
c           kf=  grid number of the fine grid
c          ifd=  the size of the fine grid stencil (= 5 or 9)
c        i9,j9=  the i,j-dimensions of the arrays asw,ase
c        qw,fw=  coarse grid portions of q and f used for work arrays he
c   (See comments in MGSS2 for details)
c** OUTPUT -
c  acc - anwc = coarse grid (kf-1) array stencil coeficients
c        pu,pd=  arrays of interpolation coefficients from grid kf-1
c                to grid kf
c .....................................................................
      dimension ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9),
     + ru(im,jmc),rd(im,jmc),rc(im,jmc),
     + pu(im,jmc),pd(im,jmc),pc(im,jmc),gam(im)
      dimension acc(im,jmc),awc(im,jmc),asc(im,jmc),aec(im,jmc),
     + anc(im,jmc),aswc(im,jmc),asec(im,jmc),anec(im,jmc),anwc(im,jmc)
      dimension qw(im,jm),fw(im,jm)
      common /io/ linp,lout
      common /prsol/ iprsol
c
      pcscale=.001
c
      im1=im-1
      jm1=jm-1
      jmc1=jmc-1
      jblack=5-jred
c                           zeroing out connections to fictitious points
      do 1 j=1,jm
      do 2 i=1,im,im1
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
    2 an(i,j)=0.
      aw(2,j)=0.
    1 ae(im1,j)=0.
      do 3 i=1,im
      do 4 j=1,jm,jm1
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
    4 an(i,j)=0.
      as(i,2)=0.
    3 an(i,jm1)=0.
 1000    if(ifd.eq.9) then
      do 5 j=1,jm
      do 6 i=1,im,im1
      asw(i,j)=0.
      ase(i,j)=0.
      ane(i,j)=0.
    6 anw(i,j)=0.
      asw(2,j)=0.
      anw(2,j)=0.
      ase(im1,j)=0.
    5 ane(im1,j)=0.
      do 7 i=1,im
      do 8 j=1,jm,jm1
      asw(i,j)=0.
      ase(i,j)=0.
      ane(i,j)=0.
    8 anw(i,j)=0.
      ase(i,2)=0.
      asw(i,2)=0.
      ane(i,jm1)=0.
    7 anw(i,jm1)=0.
 1001 endif
c
      do 9 jc=1,jmc
      do 9 i=1,im
      pc(i,jc)=0.
      pu(i,jc)=0.
      pd(i,jc)=0.
      rc(i,jc)=0.
      ru(i,jc)=0.
    9 rd(i,jc)=0.
c
c                               calculation of interpolation coeficients
c

c                                                              define pc
 2000    if(ipc.ge.1) then
      do 20 j=2,jm1
      do 20 i=2,im1
      fw(i,j)=0.
   20 qw(i,j)=1.
c
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,1)
c                                                              scale pc
      pcmax=0.
      jc=1
      do 40 j=jred,jm1,2
      jc=jc+1
      do 40 i=2,im1
      pc(i,jc)=qw(i,j)
   40 pcmax=dmax1(pcmax,abs(qw(i,j)))
      do 50 jc=2,jmc1
      do 50 i=2,im1
      if(pc(i,jc).eq.0.) pc(i,jc)=pcscale
   50 pc(i,jc)=pc(i,jc)/pcmax
c
 2001    else
      do 55 jc=2,jmc1
      do 55 i=2,im1
   55 pc(i,jc)=1.
 2002    endif
c
c                                                              define pu
      jc=3-jred
      do 60 j=jblack,jm1,2
      jc=jc+1
 4000    if(ipc.eq.0) then
      do 70 i=2,im1
   70 qw(i,j)=-an(i,j)
 5000    if(ifd.eq.9) then
      do 80 i=2,im1
   80 qw(i,j)=qw(i,j)-ane(i,j)-anw(i,j)
 5001    endif
 4001    else
      do 90 i=2,im1
   90 qw(i,j)=-an(i,j)*pc(i,jc+1)
 6000    if(ifd.eq.9) then
      do 100 i=2,im1
  100 qw(i,j)=qw(i,j)-ane(i,j)*pc(i+1,jc+1)-anw(i,j)*pc(i-1,jc+1)
 6001    endif
 4002    endif
   60 continue
c                                                          solve for pu
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c

      jc=3-jred
      do 102 j=jblack,jm1,2
      jc=jc+1
 3020    if(j.lt.jm1) then
      do 103 i=2,im1
  103 pu(i,jc)=qw(i,j)
 3021    endif
  102 continue
c
c                                                              define pd
      jc=3-jred
      do 106 j=jblack,jm1,2
      jc=jc+1
 8000    if(ipc.eq.0) then
      do 130 i=2,im1
  130 qw(i,j)=-as(i,j)
 9000    if(ifd.eq.9) then
      do 140 i=2,im1
  140 qw(i,j)=qw(i,j)-ase(i,j)-asw(i,j)
 9001    endif
c
 8001    else
c
      do 150 i=2,im1
  150 qw(i,j)=-as(i,j)*pc(i,jc)
 1100    if(ifd.eq.9) then
      do 160 i=2,im1
  160 qw(i,j)=qw(i,j)-ase(i,j)*pc(i+1,jc)-asw(i,j)*pc(i-1,jc)
 1101    endif
 8002    endif
  106 continue
c                                                          solve for pd
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c
      jc=3-jred
      do 105 j=jblack,jm1,2
      jc=jc+1
 7010    if(j.gt.2) then
      do 104 i=2,im1
  104 pd(i,jc)=qw(i,j)
 7011    endif
  105 continue
c
c                                            define restriction operator
c
c                                                              define rc
 1200    if(irc.eq.1) then
      do 500 jc=2,jmc1
      do 500 i=2,im1
  500 rc(i,jc)=pc(i,jc)
         else
      do 502 jc=2,jmc1
      do 502 i=2,im1
  502 rc(i,jc)=1.
 1201    endif
c
c                                           compute qw = -Cb(inv) * eb*
 1300    if(irurd.ge.1) then
      jc=3-jred
 3300    if(irurd.eq.1) then
      do 560 j=jblack,jm1,2
      jc=jc+1
      do 560 i=2,im1
  560 qw(i,j)=1.
 3301    elseif(irurd.eq.2) then
      do 561 j=jblack,jm1,2
      jc=jc+1
      do 561 i=2,im1
  561 qw(i,j)=(pd(i,jc)*pc(i,jc)+pu(i,jc)*pc(i,jc+1))
 3302    endif
c
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c
      jc=3-jred
      do 566 j=jblack,jm1,2
      jc=jc+1
c                                           compute ru = "-b(j+1)" * qw
 1400    if(j.lt.jm1) then
      do 570 i=2,im1
  570 ru(i,jc)=-as(i,j+1)*qw(i,j)
 1500    if(ifd.eq.9) then
      do 580 i=2,im1
  580 ru(i,jc)=ru(i,jc)-ase(i,j+1)*qw(i+1,j)-asw(i,j+1)*qw(i-1,j)
 1501    endif
 1401    endif
c                             compute rd = "-a(j-1)" * "c(j)(inv)" * qw
 1600    if(j.gt.2) then
      do 650 i=2,im1
  650 rd(i,jc)=-an(i,j-1)*qw(i,j)
 1700    if(ifd.eq.9) then
      do 660 i=2,im1
  660 rd(i,jc)=rd(i,jc)-ane(i,j-1)*qw(i+1,j)-anw(i,j-1)*qw(i-1,j)
 1701    endif
 1601    endif
  566 continue
c
 1301    else
c                                              else set ru=pu and rd=pd
      jc=3-jred
      do 670 j=jblack,jm1,2
      jc=jc+1
      do 670 i=2,im1
      ru(i,jc)=pu(i,jc)
  670 rd(i,jc)=pd(i,jc)
 1303    endif
c
c                                   calculating the coarse grid operator
c
 1800    if(ipc+irc+irurd.eq.0) then
      j=jred-2
      do 200 jc=2,jmc1
      j=j+2
      do 200 i=2,im1
      acc(i,jc)=ac(i,j)+an(i,j-1)*pu(i,jc-1)+as(i,j+1)*pd(i,jc)+
     + pu(i,jc-1)*(as(i,j)+ac(i,j-1)*pu(i,jc-1))+
     + pd(i,jc)*(an(i,j)+ac(i,j+1)*pd(i,jc))
      awc(i,jc)=aw(i,j)+pd(i-1,jc)*aw(i,j+1)*pd(i,jc)+pu(i-1,jc-1)*
     + aw(i,j-1)*pu(i,jc-1)
      asc(i,jc)=as(i,j)*pd(i,jc-1)+pu(i,jc-1)*(as(i,j-1)+
     + ac(i,j-1)*pd(i,jc-1))
      aec(i,jc)=ae(i,j)+pd(i+1,jc)*ae(i,j+1)*pd(i,jc)+pu(i+1,jc-1)*
     + ae(i,j-1)*pu(i,jc-1)
      anc(i,jc)=an(i,j)*pu(i,jc)+pd(i,jc)*(an(i,j+1)+
     + ac(i,j+1)*pu(i,jc))
      aswc(i,jc)=pd(i-1,jc-1)*aw(i,j-1)*pu(i,jc-1)
      asec(i,jc)=pd(i+1,jc-1)*ae(i,j-1)*pu(i,jc-1)
      anec(i,jc)=pu(i+1,jc)*ae(i,j+1)*pd(i,jc)
  200 anwc(i,jc)=pu(i-1,jc)*aw(i,j+1)*pd(i,jc)
 1900    if(ifd.eq.9) then
      j=jred-2
      do 210 jc=2,jmc1
      j=j+2
      do 210 i=2,im1
      awc(i,jc)=awc(i,jc)+asw(i,j+1)*pd(i,jc)+anw(i,j-1)*pu(i,jc-1)+
     + pd(i-1,jc)*anw(i,j)+pu(i-1,jc-1)*asw(i,j)
      aec(i,jc)=aec(i,jc)+ase(i,j+1)*pd(i,jc)+ane(i,j-1)*pu(i,jc-1)+
     + pd(i+1,jc)*ane(i,j)+pu(i+1,jc-1)*ase(i,j)
      aswc(i,jc)=aswc(i,jc)+asw(i,j-1)*pu(i,jc-1)+pd(i-1,jc-1)*asw(i,j)
      asec(i,jc)=asec(i,jc)+ase(i,j-1)*pu(i,jc-1)+pd(i+1,jc-1)*ase(i,j)
      anec(i,jc)=anec(i,jc)+ane(i,j+1)*pd(i,jc)+pu(i+1,jc)*ane(i,j)
  210 anwc(i,jc)=anwc(i,jc)+anw(i,j+1)*pd(i,jc)+pu(i-1,jc)*anw(i,j)
 1901    endif
c
 1801    else
c
      j=jred-2
      do 300 jc=2,jmc1
      j=j+2
      do 300 i=2,im1
      acc(i,jc)=rc(i,jc)*(ac(i,j)*pc(i,jc)+
     +                    as(i,j)*pu(i,jc-1)+
     +                    an(i,j)*pd(i,jc))+
     +        ru(i,jc-1)*(ac(i,j-1)*pu(i,jc-1)+
     +                    an(i,j-1)*pc(i,jc))+
     +          rd(i,jc)*(ac(i,j+1)*pd(i,jc)+
     +                    as(i,j+1)*pc(i,jc))
      awc(i,jc)=rc(i,jc)*aw(i,j)*pc(i-1,jc)+
     +          rd(i,jc)*aw(i,j+1)*pd(i-1,jc)+
     +        ru(i,jc-1)*aw(i,j-1)*pu(i-1,jc-1)
      asc(i,jc)=rc(i,jc)*as(i,j)*pd(i,jc-1)+
     +       ru(i,jc-1)*(ac(i,j-1)*pd(i,jc-1)+
     +                   as(i,j-1)*pc(i,jc-1))
      aec(i,jc)=rc(i,jc)*ae(i,j)*pc(i+1,jc)+
     +          rd(i,jc)*ae(i,j+1)*pd(i+1,jc)+
     +        ru(i,jc-1)*ae(i,j-1)*pu(i+1,jc-1)
      anc(i,jc)=rc(i,jc)*an(i,j)*pu(i,jc)+
     +         rd(i,jc)*(ac(i,j+1)*pu(i,jc)+
     +                   an(i,j+1)*pc(i,jc+1))
      aswc(i,jc)=ru(i,jc-1)*aw(i,j-1)*pd(i-1,jc-1)
      asec(i,jc)=ru(i,jc-1)*ae(i,j-1)*pd(i+1,jc-1)
      anec(i,jc)=rd(i,jc)*ae(i,j+1)*pu(i+1,jc)
  300 anwc(i,jc)=rd(i,jc)*aw(i,j+1)*pu(i-1,jc)
 2100    if(ifd.eq.9) then
      j=jred-2
      do 310 jc=2,jmc1
      j=j+2
      do 310 i=2,im1
      awc(i,jc)=awc(i,jc)+(rd(i,jc)*asw(i,j+1)+
     +                   ru(i,jc-1)*anw(i,j-1))*pc(i-1,jc)+
     +                    rc(i,jc)*(anw(i,j)*pd(i-1,jc)+
     +                              asw(i,j)*pu(i-1,jc-1))
      aec(i,jc)=aec(i,jc)+(rd(i,jc)*ase(i,j+1)+
     +                   ru(i,jc-1)*ane(i,j-1))*pc(i+1,jc)+
     +                    rc(i,jc)*(ane(i,j)*pd(i+1,jc)+
     +                              ase(i,j)*pu(i+1,jc-1))
      aswc(i,jc)=aswc(i,jc)+ru(i,jc-1)*asw(i,j-1)*pc(i-1,jc-1)+
     +                        rc(i,jc)*asw(i,j)*pd(i-1,jc-1)
      asec(i,jc)=asec(i,jc)+ru(i,jc-1)*ase(i,j-1)*pc(i+1,jc-1)+
     +                        rc(i,jc)*ase(i,j)*pd(i+1,jc-1)
      anec(i,jc)=anec(i,jc)+rd(i,jc)*ane(i,j+1)*pc(i+1,jc+1)+
     +                      rc(i,jc)*ane(i,j)*pu(i+1,jc)
  310 anwc(i,jc)=anwc(i,jc)+rd(i,jc)*anw(i,j+1)*pc(i-1,jc+1)+
     +                      rc(i,jc)*anw(i,j)*pu(i-1,jc)
 2101    endif
 1802    endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         if(iprsol.eq.2.and.kf.eq.m.and.ifd.eq.5) then
      do 111 j=2,jm1
      do 111 i=2,im1
               write(lout,1010) im,jm,an(i,j)
               write(lout,1012) i,j,aw(i,j),ac(i,j),ae(i,j)
 111           write(lout,1011) as(i,j)
 1010 format(2(1x,i2),14x,f12.5)
 1011 format(20x,f12.5)
 1012 format(2(1x,i2),3(1x,f12.5))
          endif
         if(iprsol.eq.2.and.kf.eq.m.and.ifd.eq.9) then
      do 115 j=2,jm1
      do 115 i=2,im1
               write(lout,1017) im,jm,anw(i,j),an(i,j),ane(i,j)
               write(lout,1017) i,j,aw(i,j),ac(i,j),ae(i,j)
 115           write(lout,1016) asw(i,j),as(i,j),ase(i,j)
 1016  format(6x,3(1x,f12.5))
 1017  format(2(1x,i2),3(1x,f12.5))
          endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(iprsol.eq.2.or.iprsol.eq.3) then
      write(lout,*) ' coarse grid kf-1 ==',kf-1
      do 211 jc=2,jmc1
      do 211 i=2,im1
      write(lout,2015) im,jmc,anwc(i,jc),anc(i,jc),anec(i,jc),pd(i,jc),
     + rd(i,jc)
      write(lout,2013) i,jc,awc(i,jc),acc(i,jc),aec(i,jc),pc(i,jc),
     + rc(i,jc)
      write(lout,2014) aswc(i,jc),asc(i,jc),asec(i,jc),pu(i,jc-1),
     + ru(i,jc-1)
 211  write(lout,*) '                          m, kf-1=',m,kf-1
 2014 format(9x,3(1x,f11.5),3x,2(f11.5,1x))
 2013 format(3x,2(1x,i2),3(1x,f11.5),3x,2(f11.5,1x))
 2015 format(3x,2(1x,i2),3(1x,f11.5),3x,2(f11.5,1x))
          endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      return
      end
