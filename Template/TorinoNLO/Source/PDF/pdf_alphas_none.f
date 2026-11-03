      double precision function pdf_alphas(q)
c     MadNkLO, lepton beams (library compiled without LHAPDF): alpha_s
c     is computed from the param_card value, never from a PDF set
      implicit none
      double precision q
      write(*,*)'pdf_alphas: no LHAPDF in this build (lepton beams); ',
     &     'with hadron beams, set lhapdf in input/madnklo_config.txt ',
     &     'and rerun bin/run_all.py'
      stop 1
      pdf_alphas=0d0
      return
      end
