      double precision function pdf_alphas(q)
c     MadNkLO, hadron beams: alpha_s(q) of the LHAPDF set (ren_scale.f)
      implicit none
      double precision q
      double precision alphasPDF
      external alphasPDF
      pdf_alphas=alphasPDF(q)
      return
      end
