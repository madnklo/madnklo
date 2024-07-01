c     BEGIN FUNCTIONS FOR THE PURE DOUBLE-UNRESOLVED COUNTERTERM
      function M2_SS(i,j,xs,xp,iA2,wgt,Wsoft,xj,ierr)
c     returns S_(i,j) * Wsoft
c     where Wsoft is the SS limit of the W function
c     it returns 0 if (i,j) is not (q,q) nor (g,g)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer i,j,l,m,jb,lb,mb,lbb,mbb,ierr,iA2
      double precision M2_SS,pref,M2tmp,wgt,wgtpl,Wsoft,xj,xjCS1,xjCS2
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xsbb(-2:maxdim,-2:maxdim)
      double precision ccBLO(maxdim,maxdim)
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision xpbb(0:3,-2:npartLO)
      double precision sij,sil,sim,sjl,sjm,slm
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_SS=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
      if(i.eq.j)then
         write(*,*)'Wrong input indices in M2_SS'
         write(*,*)i,j
         stop
      endif
c
c     double-soft limit only for gluon or quark pairs
      if((.not.isggNNLO(i,j)).and.(.not.is_NNLO_SS_qqb_pair(i,j)))return
c
c     overall kernel prefix
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do l=1,npartNNLO
         if(l.eq.i.or.l.eq.j)cycle
         do m=1,npartNNLO
            if(m.eq.i.or.m.eq.j.or.m.eq.l)cycle
c
c     in the definition of barred Sij the case l=m is excluded
c     and recovered through explicit terms I_ll and I_mm
c
c     determine indices in the n-1 and in the n-body kinematics
            jb=imap(j,i,j,0,0,npartNNLO)
            lb=imap(l,i,j,0,0,npartNNLO)
            mb=imap(m,i,j,0,0,npartNNLO)
            lbb=imap(lb,jb,lb,0,0,npartNLO)
            mbb=imap(mb,jb,lb,0,0,npartNLO)
c
c     check that id(l) = id(lb) = id(lbb), that id(m) = id(mb) = id(mbb) 
c     and that id(jb) = 21 (the remapped j is a gluon both for (i,j) =
c     (q,q), and = (g,g))
            if( idNLO(lb).ne.idNNLO(l).or.idNLO(mb).ne.idNNLO(m)
     &      .or.idLO(lbb).ne.idNLO(lb).or.idLO(mbb).ne.idNLO(mbb)
     &      .or.idNLO(jb).ne.21 )then
               write(*,*)'Inconsistent imap 1 in M2_SS'
               write(*,*)i,j,l,m,jb,lb,mb,lbb,mbb
               write(*,*)idLO(lbb),idNLO(lb),idNNLO(l),idLO(mbb),
     &                   idNLO(mb),idNNLO(m),idNLO(jb),idNNLO(j)
               stop
            endif
c
c     phase-space mapping according to l and m, at fixed radiation
c     phase-space point: the singular kernel is in the same point
c     as the double-real, ensuring numerical stability, while the
c     underlying Born configuration is remapped
            call phase_space_CS_inv(i,j,l,xp,xpb,npartNNLO,xjCS1)
            if(xjCS1.eq.0d0)goto 999
            call invariants_from_p(xpb,npartNLO,xsb,ierr)
            if(ierr.eq.1)goto 999
            call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,npartNLO,xjCS2)
            if(xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,npartLO,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
      if(docut(xpbb,npartLO))return
c
c     invariant quantities
            sij=xs(i,j)
            sil=xs(i,l)
            sim=xs(i,m)
            sjl=xs(j,l)
            sjm=xs(j,m)
            slm=xs(l,m)
c
c     safety check
            if(sil+sjl.le.0d0.or.sim+sjm.le.0d0.or.sij.le.0d0)then
               write(77,*)'Inaccuracy 1 in M2_SS',sil+sjl,sim+sjm,sij
               goto 999
            endif
c
c     call colour-connected Born
            call cc_Born_LO(xsbb,ccBLO,ierr)
            if(ierr.eq.1)goto 999
c
c     Catani-Grazzini double-soft kernel
            if(is_NNLO_SS_qqb_pair(i,j))then
c
c     given I_lm = 2*TR*(sil*sjm+sim*sjl-sij*slm)/(sij**2*(sil+sjl)*(sim+sjm))
c     SS = (I_lm - I_ll/2 - I_mm/2)*ccBLO
               M2tmp=2d0*TR*((sim*sjl-sil*sjm)**2-sij*slm*(sil+sjl)*(sim+sjm))/
     &         (sij**2*(sil+sjl)**2*(sim+sjm)**2)*ccBLO(lbb,mbb)*Wsoft
            else
               write(*,*)'M2_SS not yet implememted for gluons'
               stop
            endif
            M2_SS=M2_SS+pref*M2tmp
c
c     plot
            wgtpl=-pref*M2tmp*xj*wgt/nitRR
            if(doplot)call histo_fill(xpbb,xsbb,npartLO,wgtpl)
c
         enddo
      enddo
c
c     sanity check
      if(abs(M2_SS).ge.huge(1d0).or.isnan(M2_SS))then
         write(77,*)'Exception caught in M2_SS',M2_SS
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_HH_CC(ia,ib,ic,ir,xs,xp,xsb,xpb,xsbb,xpbb,wgt,Wcoll,Wsoftcoll,xj,ierr)
c     returns C_(ia,ib,ic) * Wcoll - S_(ia,ib) C_(ia,ib,ic) * Wsoftcoll
c     or C_(ia,ib,ic) * Wcoll if the sector is not double-soft singular
c     Wcoll and Wsoftcoll are the CC and SS_CC limits of the W function
c     ic and ir define mapping labels
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer ia,ib,ic,ir,ierr
      double precision M2_HH_CC,pref,M2tmp,wgt,wgtpl,Wcoll,Wsoftcoll,xj
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xsbb(-2:maxdim,-2:maxdim)
      double precision BLO
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision xpbb(0:3,-2:npartLO)
      double precision sab,sac,sbc,sar,sbr,scr,sabc
      double precision tabc,tacb,tcba,za,zb,zc
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_HH_CC=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
      if(ia.eq.ib.or.ia.eq.ic.or.ia.eq.ir.or.
     &   ib.eq.ic.or.ib.eq.ir.or.ic.eq.ir)then
         write(*,*)'Wrong input indices in M2_HH_CC'
         write(*,*)ia,ib,ic,ir
         stop
      endif
c
c     possible cuts
      if(docut(xpbb,npartLO))return
c
c     overall kernel prefix
      pref=32d0*pi**2*alphas**2
c
c     invariant quantities
      sab=xs(ia,ib)
      sac=xs(ia,ic)
      sbc=xs(ib,ic)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      scr=xs(ic,ir)
      sabc=sab+sac+sbc
      za=sar/(sar+sbr+scr)
      zb=sbr/(sar+sbr+scr)
      zc=scr/(sar+sbr+scr)
c
c     safety check
      if(sab.le.0d0.or.sabc.le.0d0.or.sar+sbr.le.0d0.or.sar+scr.le.0d0
     &.or.sbr+scr.le.0d0.or.sar+sbr+scr.le.0d0.or.za.le.0d0.or.za.gt.1d0
     &.or.zb.le.0d0.or.zb.gt.1d0.or.zc.le.0d0.or.zc.gt.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HH_CC',
     &   sab,sabc,sar+sbr,sar+scr,sbr+scr,sar+sbr+scr,za,zb,zc
         goto 999
      endif
c
c     call Born
      call Born_LO(xsbb,BLO,ierr)
      if(ierr.eq.1)goto 999
c
c     q -> q q' q' case
      if(isqorqbNNLO(ia).and.isqorqbNNLO(ib).and.isqorqbNNLO(ic).and.
     &idNNLO(ia).eq.-idNNLO(ib).and.abs(idNNLO(ia)).ne.abs(idNNLO(ic)))then
         tabc=2d0*(za*sbc-zb*sac)/(za+zb)+(za-zb)/(za+zb)*sab
         M2tmp=CF*TR/(sab*sabc)*(-tabc**2/(sab*sabc)+(4d0*zc+(za-zb)**2)
     &   /(za+zb)+(za+zb-sab/sabc))*BLO*Wcoll
c     subtract double-soft singularity
         if(is_NNLO_SS_singular_pair(ia,ib))
     &   M2tmp=M2tmp-4d0*CF*TR*(-(sar*sbc-sbr*sac)**2+sab*scr*(sac+sbc)*
     &   (sar+sbr))/(sab**2*(sac+sbc)**2*(sar+sbr)**2)*BLO*Wsoftcoll
      endif
c
      if(isqorqbNNLO(ia).and.isqorqbNNLO(ib).and.isqorqbNNLO(ic).and.
     &idNNLO(ia).eq.-idNNLO(ic).and.abs(idNNLO(ia)).ne.abs(idNNLO(ib)))then
         tacb=2d0*(za*sbc-zc*sab)/(za+zc)+(za-zc)/(za+zc)*sac
         M2tmp=CF*TR/(sac*sabc)*(-tacb**2/(sac*sabc)+(4d0*zb+(za-zc)**2)
     &   /(za+zc)+(za+zc-sac/sabc))*BLO*Wcoll
c     subtract double-soft singularity
         if(is_NNLO_SS_singular_pair(ia,ic))
     &   M2tmp=M2tmp-4d0*CF*TR*(-(sar*sbc-scr*sab)**2+sac*sbr*(sab+sbc)*
     &   (sar+scr))/(sac**2*(sab+sbc)**2*(sar+scr)**2)*BLO*Wsoftcoll
      endif
c
      if(isqorqbNNLO(ia).and.isqorqbNNLO(ib).and.isqorqbNNLO(ic).and.
     &idNNLO(ib).eq.-idNNLO(ic).and.abs(idNNLO(ia)).ne.abs(idNNLO(ib)))then
         tcba=2d0*(zc*sab-zb*sac)/(zc+zb)+(zc-zb)/(zc+zb)*sbc
         M2tmp=CF*TR/(sbc*sabc)*(-tcba**2/(sbc*sabc)+(4d0*za+(zc-zb)**2)
     &   /(zc+zb)+(zc+zb-sbc/sabc))*BLO*Wcoll
c     subtract double-soft singularity
         if(is_NNLO_SS_singular_pair(ic,ib))
     &   M2tmp=M2tmp-4d0*CF*TR*(-(scr*sab-sbr*sac)**2+sbc*sar*(sac+sab)*
     &   (scr+sbr))/(sbc**2*(sac+sab)**2*(scr+sbr)**2)*BLO*Wsoftcoll
      endif
c
      M2_HH_CC=pref*M2tmp
c
c     plot
      wgtpl=-pref*M2tmp*xj*wgt/nitRR
      if(doplot)call histo_fill(xpbb,xsbb,npartLO,wgtpl)
c
c     sanity check
      if(abs(M2_HH_CC).ge.huge(1d0).or.isnan(M2_HH_CC))then
         write(77,*)'Exception caught in M2_HH_CC',M2_HH_CC
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
c     END FUNCTIONS FOR THE PURE DOUBLE-UNRESOLVED COUNTERTERM


c     BEGIN FUNCTIONS FOR THE MIXED DOUBLE-UNRESOLVED COUNTERTERM
      function M2_SS_C(i,j,r,xs,xp,xsb,xpb,wgt,Wsoftcoll,xjj,ierr)
c     returns S_(i,j) * C_(i,j) * Wsoftcoll
c     where Wsoftcoll is the SS_C limit of the W function
c     it returns 0 if (i,j) is not (q,q) nor (g,g)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer i,j,r,l,m,jb,lb,mb,lbb,mbb,ierr,mu
      double precision M2_SS_C,pref,M2tmp,wgt,wgtpl,Wsoftcoll,xjj,xjCS1,xjCS2
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xsbb(-2:maxdim,-2:maxdim)
      double precision ccBLO(maxdim,maxdim)
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision xpbb(0:3,-2:npartLO)
      double precision sij,sir,sjr,sbjl,sbjm,sblm
      double precision xi,xj,Qij,ktsq,ktkbl,ktkbm
      double precision wa,wb,wr
      double precision kt(0:3),kbl(0:3),kbm(0:3),dot
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_SS_C=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
      if(i.eq.j.or.i.eq.r.or.j.eq.r)then
         write(*,*)'Wrong input indices in M2_SS_C'
         write(*,*)i,j,r
         stop
      endif
c
c     double-soft single-collinear limit only for gluon or quark pairs
      if((.not.isggNNLO(i,j)).and.(.not.is_NNLO_SS_qqb_pair(i,j)))return
c
c     overall kernel prefix
      pref=32d0*pi**2*alphas**2
c
c     eikonal double sum
      do l=1,npartNNLO
         if(l.eq.i.or.l.eq.j)cycle
         do m=1,npartNNLO
            if(m.eq.i.or.m.eq.j.or.m.eq.l)cycle
c     in the definition of barred Sij Cij the case l=m
c     is excluded and recovered through explicit terms
c
c     determine indices in the n-1 and in the n-body kinematics
            jb=imap(j,i,j,0,0,npartNNLO)
            lb=imap(l,i,j,0,0,npartNNLO)
            mb=imap(m,i,j,0,0,npartNNLO)
            lbb=imap(lb,jb,lb,0,0,npartNLO)
            mbb=imap(mb,jb,lb,0,0,npartNLO)
c
c     check that id(l) = id(lb) = id(lbb), that id(m) = id(mb) = id(mbb) 
c     and that id(jb) = 21 (the remapped j is a gluon both for (i,j) =
c     (q,q), and = (g,g))
            if( idNLO(lb).ne.idNNLO(l).or.idNLO(mb).ne.idNNLO(m)
     &      .or.idLO(lbb).ne.idNLO(lb).or.idLO(mbb).ne.idNLO(mbb)
     &      .or.idNLO(jb).ne.21 )then
               write(*,*)'Inconsistent imap 1 in M2_SS_C'
               write(*,*)i,j,l,m,jb,lb,mb,lbb,mbb
               write(*,*)idLO(lbb),idNLO(lb),idNNLO(l),idLO(mbb),
     &                   idNLO(mb),idNNLO(m),idNLO(jb),idNNLO(j)
               stop
            endif
c
c     phase-space formo n+2 to n+1 is (ijr) (or (jir)), hence it is the one that was
c     used to generate xpb and xsb: retrieve them instead of recalculating. From n+1
c     to n the mapping depends on the eikonal indices
            call phase_space_CS_inv(jb,lb,mb,xpb,xpbb,npartNLO,xjCS2)
            if(xjCS2.eq.0d0)goto 999
            call invariants_from_p(xpbb,npartLO,xsbb,ierr)
            if(ierr.eq.1)goto 999
c
c     possible cuts
            if(docut(xpbb,npartLO))return
c
c     invariant quantities and transverse momentum
            sij=xs(i,j)
            sir=xs(i,r)
            sjr=xs(j,r)
            sbjl=xsb(jb,lb)
            sbjm=xsb(jb,mb)
            sblm=xsb(lb,mb)
            xi=sir/(sir+sjr)
            xj=sjr/(sir+sjr)
c
c     coefficients of kperp
c     these are written *exactly* in the same way as done in the
c     M2_H_C_NNLO to avoid unnecessary accuracy losses when kt(mu)
c     is small, and guarantee a good numerical cancellation with
c     the single collinear limit
            wa=xj
            wb=-(1d0-xj)
            wr=(1d0-2d0*xj)*sij/(sir+sjr)
            do mu=0,3
               kt(mu)=wa*xp(mu,i)+wb*xp(mu,j)+wr*xp(mu,r)
               kbl(mu)=xpb(mu,lb)
               kbm(mu)=xpb(mu,mb)
            enddo
            ktsq=-xi*(1d0-xi)*sij
            ktkbl=dot(kt,kbl)
            ktkbm=dot(kt,kbm)
c
c     safety check
            if(xi.le.0d0.or.xi.ge.1d0.or.sij.le.0d0.or.sbjl.le.0d0.or.
     &         sbjm.le.0d0.or.ktsq.ge.0d0)then
               write(77,*)'Inaccuracy 1 in M2_SS_C',xi,sij,sbjl,sbjm,ktsq
               goto 999
            endif
c
c     call colour-connected Born
            call cc_Born_LO(xsbb,ccBLO,ierr)
            if(ierr.eq.1)goto 999
c
c     q q case
            if(is_NNLO_SS_qqb_pair(i,j))then
               Qij=TR*2d0*xi*(1d0-xi)
               M2tmp=-1d0/sij*(TR*2d0*sblm/(sbjl*sbjm)+Qij/ktsq*
     &         (2d0*ktkbl/sbjl-2d0*ktkbm/sbjm)**2)*ccBLO(lbb,mbb)*Wsoftcoll
            else
               write(*,*)'M2_SS_C not yet implememted for gluons'
               stop
            endif
            M2_SS_C=M2_SS_C+pref*M2tmp
c
c     plot
            wgtpl=pref*M2tmp*xjj*wgt/nitRR
            if(doplot)call histo_fill(xpbb,xsbb,npartLO,wgtpl)
c
         enddo
      enddo
c
c     sanity check
      if(abs(M2_SS_C).ge.huge(1d0).or.isnan(M2_SS_C))then
         write(77,*)'Exception caught in M2_SS_C',M2_SS_C
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end


      function M2_HH_CC_C(ia,ib,ic,ir,xs,xp,xsb,xpb,xsbb,xpbb,wgt,Wcoll,Wsoftcoll,xjac,itopo,ierr)
c     returns C_(ia,ib) [ C_(ia,ib,ic) * Wcoll - S_(ia,ib) C_(ia,ib,ic) * Wsoftcoll ]
c     or C_(ia,ib) C_(ia,ib,ic) * Wcoll if the sector is not double-soft singular
c     Wcoll and Wsoftcoll are the C_CC and C_SS_CC limits of the W function
c     ic and ir define mapping labels
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer ia,ib,ic,ir,iab,ibb,icb,irb,ierr,mu,itopo
      double precision M2_HH_CC_C,pref,M2tmp,wgt,wgtpl,Wcoll,Wsoftcoll,xjac
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xsbb(-2:maxdim,-2:maxdim)
      double precision BLO
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision xpbb(0:3,-2:npartLO)
      double precision sij,sik,sir,sjk,sjr,skr,sbjk
      double precision xi,xj,xjp,Pij,Qij,Pbjk,Qbjk,ktsq,ktpsq,ktktp
      double precision wa,wb,wr
      double precision kt(0:3),ktp(0:3),dot
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_HH_CC_C=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
      if(ia.eq.ib.or.ia.eq.ic.or.ia.eq.ir.or.
     &   ib.eq.ic.or.ib.eq.ir.or.ic.eq.ir)then
         write(*,*)'Wrong input indices in M2_HH_CC_C'
         write(*,*)ia,ib,ic,ir
         stop
      endif
c
c     possible cuts
      if(docut(xpbb,npartLO))return
c
c     overall kernel prefix
      pref=32d0*pi**2*alphas**2
c
c     remapped labels. Recall that this function is called with
c     arguments isec = ia, jsec = ib, lsec = ic, iref = ir for itopo = 1.
c     Note that it is not necessary to call remapped invariants, since
c     these are already the xsb passed as argument of this function.
      iab=imap(ia,ia,ib,0,0,npartNNLO)
      ibb=imap(ib,ia,ib,0,0,npartNNLO)
      icb=imap(ic,ia,ib,0,0,npartNNLO)
      irb=imap(ir,ia,ib,0,0,npartNNLO)
c
c     check remapped identities
      if(idNNLO(ic).ne.idNLO(icb).or.idNNLO(ir).ne.idNLO(irb).or.
     &   (isgNNLO(ia).and.idNNLO(ib).ne.idNLO(ibb)).or.
     &   (isgNNLO(ib).and.idNNLO(ia).ne.idNLO(iab)).or.
     &   (isqorqbNNLO(ia).and.isqorqbNNLO(ib).and..not.isgNLO(ibb)))then
         write(*,*)'Inconsistent imap 1 in M2_HH_CC_C'
         write(*,*)ia,ib,ic,ir,iab,ibb,icb,irb
         write(*,*)idNNLO(ia),idNNLO(ib),idNNLO(ic),idNNLO(ir)
         write(*,*)idNLO(iab),idNLO(ibb),idNLO(icb),idNLO(irb)
         stop
      endif
c
c     invariant quantities
      sij=xs(ia,ib)
      sik=xs(ia,ic)
      sir=xs(ia,ir)
      sjk=xs(ib,ic)
      sjr=xs(ib,ir)
      skr=xs(ic,ir)
      xi=sir/(sir+sjr)
      xj=sjr/(sir+sjr)
c      xjp=xsb(ibb,irb)/(xsb(ibb,irb)+xsb(icb,irb))
      xjp=(sir+sjr)/(sir+sjr+skr)
c
c     coefficients of kperp
c     these are written *exactly* in the same way as done in the
c     M2_H_C_NNLO to avoid unnecessary accuracy losses when kt(mu)
c     is small, and guarantee a good numerical cancellation with
c     the single collinear limit
      wa=xj
      wb=-(1d0-xj)
      wr=(1d0-2d0*xj)*sij/(sir+sjr)
      do mu=0,3
         kt(mu)=wa*xp(mu,ia)+wb*xp(mu,ib)+wr*xp(mu,ir)
      enddo
c
      sbjk=xsb(ibb,icb)
      ktsq=-xi*(1d0-xi)*sij
      ktpsq=-xjp*(1d0-xjp)*sbjk
c     in ktp the components along bar(kj) and bar(kr) give 0
c     when contracted with kt
      ktktp=xjp*dot(kt,xpb(0,icb))
c
c     safety check
      if(xi.le.0d0.or.xi.ge.1d0.or.xjp.le.0d0.or.xjp.ge.1d0.or.
     &   ktsq.ge.0d0.or.ktpsq.ge.0d0)then
         write(77,*)'Inaccuracy 1 in M2_HH_CC_C',xi,xjp,ktsq,ktpsq
         goto 999
      endif
c
c     call Born
      call Born_LO(xsbb,BLO,ierr)
      if(ierr.eq.1)goto 999
c
c     q -> q q' q' case
      if(isqorqbNNLO(ia).and.isqorqbNNLO(ib).and.isqorqbNNLO(ic).and.
     &idNNLO(ia).eq.-idNNLO(ib).and.abs(idNNLO(ia)).ne.abs(idNNLO(ic)))then
         Qij=TR*2d0*xi*(1d0-xi)
         Pbjk=CF*(1d0+(1d0-xjp)**2)/xjp
         M2tmp=2d0/sij/sbjk*(TR*Pbjk-CF*Qij*xjp-CF*Qij*(1d0-xjp)/xjp
     &        *(2d0*ktktp)**2/ktsq/ktpsq)*BLO*Wcoll
c
c     subtract double-soft singularity: this function is called with ia =
c     isec (first sector index) and ib = jsec (second sector index)
c     Here topology is checked. Check that this is fine in all cases
         if(is_NNLO_SS_qqb_pair(ia,ib).and.itopo.eq.1)then
            M2tmp=M2tmp-4d0/sij/sbjk*CF*(1d0-xjp)/xjp*
     &      (TR-Qij/2d0*(2d0*ktktp)**2/ktsq/ktpsq)*BLO*Wsoftcoll
         endif
      endif
c
      M2_HH_CC_C=pref*M2tmp
c
c     plot
      wgtpl=pref*M2tmp*xjac*wgt/nitRR
      if(doplot)call histo_fill(xpbb,xsbb,npartLO,wgtpl)
c
c     sanity check
      if(abs(M2_HH_CC_C).ge.huge(1d0).or.isnan(M2_HH_CC_C))then
         write(77,*)'Exception caught in M2_HH_CC_C',M2_HH_CC_C
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
c     END FUNCTIONS FOR THE MIXED DOUBLE-UNRESOLVED COUNTERTERM


c     BEGIN FUNCTIONS FOR THE SINGLE-UNRESOLVED COUNTERTERM
      function M2_H_C_NNLO(ia,ib,ir,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,Wcoll,ierr)
c     returns C_(ia,ib) * Wcoll - S_(ia) C_(ia,ib)
c     where Wcoll is the C limit of the W function
c     (the Wsoftcoll weight is = 1)
      implicit none
      include 'dims.inc'
      include 'setup.inc'
      integer ia,ib,ir,ierr,nit
      integer iU1,iS1,iB1,iU2,iS2,iB2
      double precision M2_H_C_NNLO,pref,M2tmp,wgt,wgtpl,Wcoll,xj
      double precision xs(-2:maxdim,-2:maxdim),xsb(-2:maxdim,-2:maxdim)
      double precision xsbb(-2:maxdim,-2:maxdim)
      double precision RNLO,RNLOkp,peps(-2:npartNLO),epseps,z
      double precision xp(0:3,-2:npartNNLO),xpb(0:3,-2:npartNLO)
      double precision xpbb(0:3,-2:npartLO)
      double precision sab,sar,sbr
      double precision wa,wb,wr
      logical doplot
      common/cdoplot/doplot
c
c     initialise
      M2_H_C_NNLO=0d0
      M2tmp=0d0
      ierr=0
c
c     input check
      if(ia.eq.ib.or.ia.eq.ir.or.ib.eq.ir)then
         write(*,*)'Wrong input indices in M2_H_C_NNLO'
         write(*,*)ia,ib,ir
         stop
      endif
c
c     possible cuts
      if(docut(xpb,npartNLO))return
c
c     overall kernel prefix
      pref=8d0*pi*alphas
c
c     invariant quantities
      sab=xs(ia,ib)
      sar=xs(ia,ir)
      sbr=xs(ib,ir)
      z=sbr/(sar+sbr)
c
c     coefficients of kperp
      wa=z
      wb=-(1d0-z)
      wr=(1d0-2d0*z)*sab/(sar+sbr)
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.z.le.0d0.or.z.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_H_C_NLO',sab,sar+sbr,z
         goto 999
      endif
c
c     call real
      call real_NLO(xsb,RNLO,ierr)
      if(ierr.eq.1)goto 999
c
c     q -> gq case
      if(isgNNLO(ia).and.isqorqbNNLO(ib))then
         write(*,*)'implement in H_C_NNLO'
         stop
c
c     q -> qg case
      elseif(isgNNLO(ib).and.isqorqbNNLO(ia))then
         write(*,*)'implement in H_C_NNLO'
         stop
c
c     g -> gg case
      elseif(isggNNLO(ia,ib))then
         write(*,*)'implement in H_C_NNLO'
         stop
c
c     g -> qq case
      elseif(isqqbNNLO(ia,ib))then
         call get_eps(ia,ib,ir,xp,xpb,npartNNLO,wa,wb,wr,peps,epseps,ierr)
         if(ierr.eq.1)goto 999
         call real_NLO_kp(xsb,peps,epseps,RNLOkp,ierr)
         if(ierr.eq.1)goto 999
         M2tmp=TR*(RNLO-4d0/sab*RNLOkp)/sab*Wcoll
      endif
c
      M2_H_C_NNLO=pref*M2tmp
c
c     plot
      wgtpl=-pref*M2tmp*xj*wgt/nitRR
      if(doplot)call histo_fill(xpb,xsb,npartNLO,wgtpl)
c
c     sanity check
      if(abs(M2_H_C_NNLO).ge.huge(1d0).or.isnan(M2_H_C_NNLO))then
         write(77,*)'Exception caught in M2_H_C_NNLO',M2_H_C_NNLO
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end
c     END FUNCTIONS FOR THE SINGLE-UNRESOLVED COUNTERTERM
