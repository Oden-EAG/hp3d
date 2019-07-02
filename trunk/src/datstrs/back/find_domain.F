c----------------------------------------------------------------------
c
c   routine name       - find_domain
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c   purpose            - routine finds domain for an element
c
c   arguments :
c     in:
c             Mdle     - an element middle node, identified with
c                        the element
c     out:
c             Ndomain  - GMP domain number 
c
c   required  routines - 
c
c----------------------------------------------------------------------
c 
      subroutine find_domain(Mdle, Ndom)
c
      use GMP
      use data_structure3D
c
      iprint=0
c
      nfath = NODES(Mdle)%father
      do while(nfath.gt.0)
        nfath = NODES(nfath)%father
      enddo
      nel = abs(nfath)
      Ndom = Ndomain_gmp(ELEMS(nel)%GMPblock)
c
c
      end

      

