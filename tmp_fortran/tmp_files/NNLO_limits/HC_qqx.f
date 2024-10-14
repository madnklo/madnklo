

      double precision function M2_HC_qqx(ia,ib,ir,xs,xp,xsb,xpb,wgt,xj,nit,extra,wgt_chan,ierr)
c     hard-collinear limit C_(ia,ib)
c     this is meant to represent the full hard-collinear
c     for sectors (ia,ib)+(ib,ia)
      USE SECTORS2_MODULE
      USE SECTORS4_MODULE
      implicit none
      include 'nexternal.inc'
      INCLUDE 'coupl.inc'
      include 'math.inc'
      include 'damping_factors.inc'
      include 'nsqso_born.inc'
      include 'leg_PDGs.inc'
      INCLUDE 'input.inc'
      INCLUDE 'run.inc'      
      integer ia,ib,ir,ic,id,ierr,nit,parent,sec_index(2)
      double precision pref,M2tmp,wgt,wgtpl,wgt_chan,xj,extra
      double precision xs(nexternal,nexternal),xsb(nexternal-1,nexternal-1)
      double precision RNLO,KKRNLO
      double precision xp(0:3,nexternal),xpb(0:3,nexternal-1),kt(0:3)
      double precision sab,sar,sbr,x,y,xinit,damp
      double precision wa,wb,wr
      double precision ANS(0:NSQSO_BORN)
      integer mapped_labels(nexternal),mapped_flavours(nexternal)
      integer, parameter :: hel = - 1
c     set logical doplot
      logical doplot
      common/cdoplot/doplot
      double precision sCM
      common/cscm/sCM
      logical docut
      integer %(proc_prefix_rr)s_fl_factor
      common/%(proc_prefix_rr)s_flavour_factor/%(proc_prefix_rr)s_fl_factor
      double precision alphas,alpha_qcd
      double precision %(proc_prefix_HC_qqx)s_get_kkblo
      integer %(proc_prefix_rr)s_den
      common/%(proc_prefix_rr)s_iden/%(proc_prefix_rr)s_den
      integer %(proc_prefix_HC_qqx)s_den
      common/%(proc_prefix_HC_qqx)s_iden/%(proc_prefix_HC_qqx)s_den
      INTEGER ISEC,JSEC,KSEC,LSEC
      COMMON/CSECINDICES/ISEC,JSEC,KSEC,LSEC
      INTEGER BORN_LEG_PDGS(NEXTERNAL-1)
      INTEGER UNDERLYING_LEG_PDGS(NEXTERNAL-1)
      integer mapped_sec(2,nexternal)
      integer i,j,k
c
c     initialise
      M2_HC_qqx=0d0
      M2tmp=0d0
      ierr=0
      damp=0d0
      ic=0
      id=0
c
c     initial checks and label assignment
      if(lsec.eq.0)then
         if((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)) then
            ic = ksec
         elseif((ia.eq.isec.and.ib.eq.ksec).or.(ia.eq.ksec.and.ib.eq.isec)) then
            ic = jsec
         elseif((ia.eq.jsec.and.ib.eq.ksec).or.(ia.eq.ksec.and.ib.eq.jsec)) then
            ic = isec
         else
            write(*,*)'Wrong indices 1 in M2_HC_qqx'
            write(*,*)ia,ib,isec,jsec,ksec
            stop
         endif
      else
         if((ia.eq.isec.and.ib.eq.jsec).or.(ia.eq.jsec.and.ib.eq.isec)) then
            ic = ksec
            id = lsec
         elseif((ia.eq.ksec.and.ib.eq.lsec).or.(ia.eq.lsec.and.ib.eq.ksec)) then
            ic = isec
            id = jsec
         else
            write(*,*)'Wrong indices 2 in M2_HC_qqx'
            write(*,*)ia,ib,isec,jsec,ksec,lsec
            stop
         endif
      endif
c
c     possible cuts
c      call GET_BORN_PDGS(ISEC,JSEC,NEXTERNAL-1,BORN_LEG_PDGS)
      call GET_UNDERLYING_PDGS(ISEC,JSEC,KSEC,LSEC,NEXTERNAL-1,UNDERLYING_LEG_PDGS)

      IF(DOCUT(XPB,NEXTERNAL-1,UNDERLYING_LEG_PDGS,0))RETURN
c
c     overall kernel prefix
      alphas=alpha_QCD(asmz,nloop,scale)
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
      wa= 1d0 - x
      wb= - x
      wr= - (1d0-2d0*x)*sab/(sar+sbr)
      kt(:) = wa*xp(:,ia) + wb*xp(:,ib) + wr*xp(:,ir) 
c
c     safety check
      if(sab.le.0d0.or.sar+sbr.le.0d0.or.x.le.0d0.or.x.ge.1d0)then
         write(77,*)'Inaccuracy 1 in M2_HC_qqx',sab,sar+sbr,x
         goto 999
      endif
c
c     call Born
      call %(proc_prefix_HC_qqx)s_ME_ACCESSOR_HOOK(xpb,hel,alphas,ANS)
      RNLO = ANS(0)
c
c     assign mapped labels and flavours
      call get_collinear_mapped_labels(ia,ib,nexternal,leg_PDGs,mapped_labels,mapped_flavours)
      parent = mapped_labels(ib)
      if(mapped_flavours(ib).ne.21)then
         write(*,*) 'Wrong parent particle label!', ib, mapped_flavours(ib)
         stop
      endif
c
c     call remapped sector function
      CALL GET_SIG2(XSB,1D0,NEXTERNAL-1)
      if(lsec.eq.0)then
         sec_index(1) = parent
         sec_index(2) = mapped_labels(ic)
      else
         sec_index(1) = mapped_labels(ic)
         sec_index(2) = mapped_labels(id)
      endif

c     Fill  mapped_sec_list(2,nexternal) with the pairs of all the
c     final state particles after mapping n+2 --> n+1

      k = 0
      do i=3,nexternal-1
         do j=i+1,nexternal
            k=k+1
            mapped_sec(1,k) = mapped_labels(i)
            mapped_sec(2,k) = mapped_labels(j)
         enddo
      enddo

      
      CALL GET_ZHC_NNLO(sec_index(1),sec_index(2),mapped_sec)
c
      KKRNLO = %(proc_prefix_HC_qqx)s_GET_KKBLO(parent,xpb,kt)
c     TODO: improve ktmuktnuBmunu / kt^2
      M2tmp=TR*(RNLO-4d0/sab*KKRNLO)*Z_HC_NNLO
c     Including correct multiplicity factor
      M2tmp = M2tmp*dble(%(proc_prefix_HC_qqx)s_den)/dble(%(proc_prefix_rr)s_den)
c     account for different damping factors according to
c     recoiler position (ir)
      if(ir.ge.2)then
         damp=(1d0-y)**beta_FF
      else
         damp=xinit**beta_FI
      endif
      M2tmp=M2tmp*damp
      M2_HC_qqx=M2tmp*pref/sab*xj*extra
c     apply flavour factor
      M2_HC_qqx=M2_HC_qqx*%(proc_prefix_rr)s_fl_factor
c
c     plot
      wgtpl=-M2_HC_qqx*wgt/nit*wgt_chan
      if(doplot)call histo_fill(xpb,xsb,nexternal-1,wgtpl)
c
c     sanity check
      if(abs(M2_HC_qqx).ge.huge(1d0).or.isnan(M2_HC_qqx))then
         write(77,*)'Exception caught in M2_HC_qqx',M2_HC_qqx
         goto 999
      endif
c
      return
 999  ierr=1
      return
      end

