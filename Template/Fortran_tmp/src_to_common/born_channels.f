c=======================================================================
c     Integration channels of the n-body (Born-like) integrands: the
c     phase-space configurations of configs.f (MAPCONFIG(0) of them),
c     with the multi-channel weights of the diagrams that have one.
c     Used by the Born (LO_B.f, driver_n.f) and by the virtual
c     (NLO_V_*.f, driver_v.f). Diagrams without a configuration (e.g.
c     with a four-gluon vertex only, as in g g > g g) have no channel of
c     their own: looping over all diagrams instead (N_MAX_CG) gives
c     wrong channel weights and an out-of-range configuration index.
c=======================================================================
      integer function n_born_channels()
c     Number of integration channels: the phase-space configurations of
c     configs.f (MAPCONFIG(0)). Diagrams without a configuration (e.g.
c     with a four-gluon vertex only) have no channel of their own.
      implicit none
      include 'maxconfigs.inc'
      integer mapconfig(0:lmaxconfigs),this_config
      common/to_mconfigs/mapconfig,this_config
      call configs_born
      n_born_channels=mapconfig(0)
      return
      end


      double precision function born_chan_weight(ich)
c     Multi-channel weight of channel ich (configuration ich of
c     configs.f, diagram MAPCONFIG(ich)), from the |AMP|^2 of the last
c     matrix-element call (AMP2, normalised over all the diagrams), so
c     that the weights of the channels add up to 1. When every diagram
c     has its own configuration in the same order, this is AMP2(ich).
      implicit none
      include 'maxconfigs.inc'
      include 'ngraphs.inc'
      integer ich,k
      integer mapconfig(0:lmaxconfigs),this_config
      common/to_mconfigs/mapconfig,this_config
      integer NGRAPHS2
      double precision amp2(N_MAX_CG)
      COMMON/TO_AMP2/AMP2,NGRAPHS2
      logical identity
      double precision sum
c
      identity=(mapconfig(0).eq.N_MAX_CG)
      do k=1,mapconfig(0)
         if(mapconfig(k).ne.k)identity=.false.
      enddo
      if(identity)then
         born_chan_weight=amp2(ich)
      else
         sum=0d0
         do k=1,mapconfig(0)
            sum=sum+amp2(mapconfig(k))
         enddo
         born_chan_weight=0d0
         if(sum.gt.0d0)born_chan_weight=amp2(mapconfig(ich))/sum
      endif
      return
      end
