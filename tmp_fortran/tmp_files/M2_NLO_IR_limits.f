      function M2_S(i,xs,xp,wgt,Wsoft,xj,nit,extra,ierr)
c     single-soft limit S_(i) * Wsoft
c     it returns 0 if i is not a gluon
      implicit none
c      include 'setup.inc'
      include 'nexternal.inc'
      include 'damping_factors.inc'
      integer i,l,m,lb,mb,ierr,nit
      double precision M2_S,pref,M2tmp,wgt,wgtpl,Wsoft,xj,xjCS
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,ccBLO(nexternal-1,nexternal-1),extra
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sil,sim,slm,y,z
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_S=0d0
      M2tmp=0d0
      ierr=0
c
c     soft limit only for gluons
      if(.not.isgNLO(i))return
c
c     overall kernel prefix
c     TODO: read pi*alphas
      pref=-8d0*pi*alphas
c
c     eikonal double sum
c     TODO: multiply mult_factors = 2d0
c     TODO: m and l have to run on a list
c           of colored partons of the process
      do m=1,nexternal
         if(m.eq.i)cycle
         do l=1,m-1
            if(l.eq.i)cycle
c
c     determine indices in the n-1 body kinematics
c     TODO: write include file for imap labels
            lb=imap(l,i,l,0,0,npartNLO)
            mb=imap(m,i,l,0,0,npartNLO)
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the single-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,l,m,xp,xpb,nexternal,xjCS)
            if(xjCS.eq.0d0)goto 999
            call invariants_from_p(xpb,nexternal-1,xsb,ierr)
            if(ierr.eq.1)goto 999
c     TODO: look at phase_space_CS_inv(), invariants_from_p()
c
c     possible cuts
c     TODO: look at docut()
            if(docut(xpb,nexternal-1))cycle
c
c     invariant quantities
            sil=xs(i,l)
            sim=xs(i,m)
            slm=xs(l,m)
            y=sil/(sil+sim+slm)
            z=sim/(sim+slm)
c     TODO: add other variables for damping
c
c     safety check
            if(sil*sim.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_S',sil,sim
               goto 999
            endif
c
c     call colour-connected Born
c     TODO: look at cc_Born_LO()
            call cc_Born_LO(xsb,ccBLO,ierr)
            if(ierr.eq.1)goto 999
c
c     eikonal
            if(m.ge.2.and.l.ge.2)then
               M2tmp=ccBLO(lb,mb)*2d0*slm/(sil*sim)
     &            *((1d0-y)*(1d0-z))**alpha
            elseif(m.ge.2.and.l.le.2)then
c              M2tmp=... ()()**alpha
            elseif(m.le.2.and.l.le.2)then
c              M2tmp=... ()()**alpha
            endif
            M2_S=M2_S+pref*M2tmp*Wsoft*extra
c
c     plot
            wgtpl=-pref*M2tmp*Wsoft*extra*xj*wgt/nit
            if(doplot)call histo_fill(xpb,xsb,npartLO,wgtpl)
c
         enddo
      enddo
c
c     sanity check
      if(abs(M2_S).ge.huge(1d0).or.isnan(M2_S))then
         write(77,*)'Exception caught in M2_S',M2_S
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_H_C_FgFg(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib) - S_(ib)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'damping_factors.inc'
      integer ia,ib,ir,ierr,nit
      double precision M2_H_C,pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,BLOkp,pkt(nexternal-1)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),ktkt
      double precision sab,sar,sbr,x,y,xinit
      double precision wa,wb,wr
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_H_C=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
c     TODO: check only first time
      if(ia.eq.ib.or.ia.eq.ir.or.ib.eq.ir)then
         write(*,*)'Wrong input indices in M2_H_C'
         write(*,*)ia,ib,ir
         stop
      endif
c
c     possible cuts
      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     coefficients of kt
c     kt = wa pa + wb pb + wr pr 
      wa=x
      wb=-(1d0-x)
      wr=(1d0-2d0*x)*sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_H_C',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
c     TODO: look at Born_LO()
      call Born_LO(xsb,BLO,ierr)
      if(ierr.eq.1)goto 999
c
c     TODO: this was get_eps()
      call get_kt(ia,ib,ir,xp,xpb,nexternal,wa,wb,wr,pkt,ktkt)
      call Born_LO_kp(xsb,pkt,ktkt,BLOkp,ierr)
      if(ierr.eq.1)goto 999
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         M2tmp=CA*2d0*(2d0/sab*BLOkp+
     &                x/(1d0-x)*(1d0-x**alpha)*BLO+
     &                (1d0-x)/x*(1d0-(1d0-x)**alpha)*BLO)
     &           *(1d0-y)**beta_FF
      else
         M2tmp=CA*2d0*(2d0/sab*BLOkp+
     &                x/(1d0-x)*(1d0-x**alpha)*BLO+
     &                (1d0-x)/x*(1d0-(1d0-x)**alpha)*BLO)
     &           *xinit**beta_FI
      endif
c
      M2_H_C=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C*xj*wgt/nit
      if(doplot)call histo_fill(xpb,xsb,npartLO,wgtpl)
c
c     sanity check
      if(abs(M2_H_C).ge.huge(1d0).or.isnan(M2_H_C))then
         write(77,*)'Exception caught in M2_H_C',M2_H_C
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
  
                  
      function M2_H_C_FgFq(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib) - S_(ia)C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'damping_factors.inc'
      integer ia,ib,ir,ierr,nit
      double precision M2_H_C,pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1)
      double precision sab,sar,sbr,x,y,xinit
      double precision wa,wb,wr
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_H_C=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
c     TODO: check only first time
      if(ia.eq.ib.or.ia.eq.ir.or.ib.eq.ir)then
         write(*,*)'Wrong input indices in M2_H_C'
         write(*,*)ia,ib,ir
         stop
      endif
c
c     possible cuts
c     TODO 
      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_H_C',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
c     TODO: look at Born_LO()
      call Born_LO(xsb,BLO,ierr)
      if(ierr.eq.1)goto 999
c
c     TODO: split the various flavours (without if)
c     q -> gq case
      if(ir.ge.2)then
         M2tmp=BLO*CF*((1d0-x)+
     &                 2d0*x/(1d0-x)*(1d0-x**alpha))
     &            *(1d0-y)**beta_FF
      else
         M2tmp=BLO*CF*((1d0-x)+
     &                 2d0*x/(1d0-x)*(1d0-x**alpha))
     &            *xinit**beta_FI
c
      M2_H_C=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C*xj*wgt/nit
c     TODO: look at histo_fill()
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C).ge.huge(1d0).or.isnan(M2_H_C))then
         write(77,*)'Exception caught in M2_H_C',M2_H_C
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_H_C_FqFqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,ierr)
c     hard-collinear limit C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      implicit none
      include 'nexternal.inc'
      include 'damping_factors.inc'
      integer ia,ib,ir,ierr,nit
      double precision M2_H_C,pref,M2tmp,wgt,wgtpl,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision BLO,BLOkp,pkt(nexternal-1)
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),ktkt
      double precision sab,sar,sbr,x,y,xinit
      double precision wa,wb,wr
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_H_C=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
c     TODO: check only first time
      if(ia.eq.ib.or.ia.eq.ir.or.ib.eq.ir)then
         write(*,*)'Wrong input indices in M2_H_C'
         write(*,*)ia,ib,ir
         stop
      endif
c
c     possible cuts
      if(docut(xpb,nexternal-1))return
c
c     overall kernel prefix
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      x=sar/(sar+sbr)
      y=sab/(sab+sar+sbr)
      xinit = 1d0 - sab/(sar+sbr)
c
c     coefficients of kt
c     kt = wa pa + wb pb + wr pr 
      wa=x
      wb=-(1d0-x)
      wr=(1d0-2d0*x)*sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_H_C',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
c     TODO: look at Born_LO()
      call Born_LO(xsb,BLO,ierr)
      if(ierr.eq.1)goto 999
c
c     g -> qq case
      call get_kt(ia,ib,ir,xp,xpb,nexternal,wa,wb,wr,pkt,ktkt)
      call Born_LO_kp(xsb,pkt,ktkt,BLOkp,ierr)
      if(ierr.eq.1)goto 999
      if(ir.ge.2)then
         M2tmp=TR*(BLO-4d0/sab*BLOkp)*(1d0-y)**beta_FF
      else
         M2tmp=TR*(BLO-4d0/sab*BLOkp)*xinit**beta_FI
      endif
c
      M2_H_C=M2tmp*pref/sab*extra
c
c     plot
      wgtpl=-M2_H_C*xj*wgt/nit
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_H_C).ge.huge(1d0).or.isnan(M2_H_C))then
         write(77,*)'Exception caught in M2_H_C',M2_H_C
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end