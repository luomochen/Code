      program likheng
ctc 
c   Space group finder/checker for crystalline structures
c   Christmas 1990
c
c  convention  {R|t}r = Rr + t
ctc
      implicit double precision (a-h,o-z)
      parameter(nal=200,na=200)
      dimension a(3,3),coor(na,3),id(na),x(3),org(3)
      dimension rot(3,3),bvec(3,3)
      dimension  dex(nal),irank(nal),index(nal)
      dimension idxx(nal)
      dimension  cx(nal,3),cxx(nal,3),cyy(nal,3),dx(nal,3)
      dimension  gro(48,3,3), dxt(48,3)
      common/rl/aij(3,3)
      common /traf/ ro(48,3,3), tnp(48,3),ntrans
      character*5 grponc
      data small/1.d-10/, zero/0.d0/
c234567
c     open(unit=5,file='groupin',status='old',form='formatted')
c     open(unit=6,file='groupout',status='unknown',form='formatted')
c
c    read in lattice vectors and atomic positions
c
      read(5,*) iguess,ipr
      read(5,*) (org(j),j=1,3)
      read(5,*) ((a(i,j),j=1,3),i=1,3)
czw   write(6,*) ' origin '
czw   write(6,611) (org(j),j=1,3)
czw   write(6,*)  ' a vectors '
czw   write(6,610) ((a(i,j),j=1,3),i=1,3)
      call rla(a)
 612  format(3f10.8,i5)
czw   write(6,*)  ' atomic positions (in lattice coordinates) '
      read(5,*) natom
      if(natom.gt.na) stop
      do 12 i=1,natom
      read(5,*) (coor(i,j),j=1,3), id(i)
      do 122 j=1,3
      coor(i,j) = coor(i,j) - org(j)
 122  continue
      do 123 j=1,3
      if( abs(coor(i,j)) .gt. 1.d0 ) then
czw      write(6,*) ' please put atoms inside unit cell ' 
         stop
      endif
      if( coor(i,j).lt. 0.d0) coor(i,j) = coor(i,j) + 1.d0
 123  continue
czw   write(6,612)  (coor(i,j), j=1,3), id(i)
 12    continue
c
c read in or generate operations, "rotation" part only
c
      if(iguess.eq.1) then
       call autogen
       else
       call trafo
      endif
c
c    generate a cluster
c
      nn12 = 100*100
      nn1  = 100
      iall = natom
      if(natom.gt.nal) stop
      do 22  j=1,natom
      do 24 k=1,3
 24   cx(j,k) = coor(j,k)
      dex(j) = cx(j,1)* nn12 + cx(j,2)*nn1 + cx(j,3)
 22   continue
 20   continue
c
c    order the cluster in some way
c
      call tri(dex,index  ,natom,irank)
c    
czw   write(6,*) ' cluster after ordering '
czw   write(6,*) ' no of atoms in cluster = ', natom
      do 40 i=1,natom
      ii = index(i)
      do 44 k=1,3
      cxx(i,k) = cx(ii,k)
      idxx(i) = id(ii)
 44   continue
czw   write(6,460)  (cxx(i,k),k=1,3),idxx(i),dex(ii)
 40   continue
 460  format(3f10.5,i5,f20.5)
c
c   loop over rotation
c
      ngroup = 0
      nonsym = 0
      do 80 ng=1,ntrans

czw   write(6,*)  ' '
czw   write(6,*) ' operation no ', ng
czw   write(6,610) ((ro(ng,i,k),k=1,3),i=1,3)  

      if(iguess.ne.1) then
       call checkgr(ng,iflag)
czw   if(iflag.gt.0) 
czw  $ write(6,*) ' not group operation for lattice, iflag=', iflag
      if(iflag.gt.0) goto 80
      endif
  
      do 82 i=1,natom
c
c operate on the atoms in the cluster
c      
      do 822 j=1,3
      cx(i,j) = zero
      do 822 k=1,3
      cx(i,j) = cxx(i,k)*ro(ng,k,j) + cx(i,j)
 822  continue
c
c translate back to first quadrant
c
      do 824 k=1,3
      if( cx(i,k) .lt. 0.d0 ) cx(i,k) = cx(i,k) + 1.d0
 824  continue
 82   continue
c
c find all possible non-primitive translations 
c
      nany = 1
      do 86 i=1,natom
      do 862 k=1,3
c//sign
      tt = cxx(i,k) - cx(nany,k)
      dx(i,k) = round(tt,1)
czw   if( abs(dx(i,k)). gt. 0.90) write(6,*) tt
 862  continue
 86   continue
c
c see if any of the non-primitive translations is good
c
      do 88 j=1,natom
czw   if(ipr.eq.1) write(6,8661) (dx(j,k),k=1,3)
      do 882 i=1,natom
      do 884 k=1,3
c//sign
      tt = cx(i,k) + dx(j,k)
      cyy (i,k) = round(tt,2)
 884  continue
      dex(i) = cyy(i,1)* nn12 + cyy(i,2)*nn1 + cyy(i,3)
 882  continue
c
c order them
c
      call tri(dex,index  ,natom,irank)
c
c see if the transformed cluster (rotation + translation) is equivalent
c to the original cluster
c
      nogood = 0
c- check species first
      do 83 i=1,natom
      ii = index(i)
      if(idxx(i).ne.idxx(ii)) then
        nogood = 1
        goto 833
      endif
 83   continue
c- then check position
      do 84 i=1,natom
      ii = index(i)
cc      if(idxx(i).ne.idxx(ii)) nogood = nogood + 1
      do 842 k=1,3
      del = cxx(i,k) - cyy(ii,k)
      del = abs(del)
      if( del.lt. 1.d-5 .or.  abs(del-1.d0).lt. 1.d-5) goto 842 
      nogood = nogood + 1
 842   continue
czw   if(ipr.eq.1) 
czw  $ write(6,8662) (cyy(ii,k),k=1,3),idxx(ii),(cxx(i,k),k=1,3),idxx(i)
 84   continue
 8661 format('  dx ', 3f10.5)
 8662 format(3f10.5,i5,3f10.5,i5)
 833  continue
      if(nogood.eq.0) indtr = j
      if(nogood.eq.0) goto 888
 88   continue

 888  continue
      if(nogood.eq.0) then
      ngroup = ngroup + 1
czw   write(6,6101)  (dx(indtr,k), k =1,3)
      do 900 i=1,3
      dxt(ngroup,i)=dx(indtr,i)
      do 900 j=1,3
      gro(ngroup,i,j)=ro(ng,i,j)
 900  continue
      sum = abs(dx(indtr,1)) + abs(dx(indtr,2)) + abs(dx(indtr,3))
      if(sum.gt.1.d-5) nonsym=nonsym + 1
      else
czw   write(6,*)  ' not a group operation for the basis '
      endif
 80   continue

 610  format(3f10.5)
 6101 format( '      [',3f10.5,'   ]')
 611  format(8f10.5)
 
czw   write(6,*) ' '
czw   if(iguess.eq.1) then
czw   write(6,*) ' no of operations allowed in autogen ', ntrans
czw   else
czw   write(6,*) ' no of operations checked from input ', ntrans
czw   endif
czw   write(6,*) ' no of group operations found ', ngroup
czw   write(6,*) ' no of symmorphic group operations ', ngroup-nonsym
c
c write in format suitable for TRAFO input
c
czw   write(6,*) ' '
czw   write(6,*) ' for TRAFO input '
      grponc='group'
      write(6,6010)   ngroup
      write(6,6011) (((gro(i,j,k),k=1,3),j=1,3),i=1,ngroup)  
czw   write(6,6012) ( ( dxt(i,l),l=1,3),i=1,ngroup )     
 6010 format (  i3 )
 6011 format ( 18f4.0 ) 
 6012 format  (8f10.5  )  
czw  the following three lines added by C. Z. Wang. 
      read (5,*)
      read (5,*) ((bvec(i,j),j=1,3),i=1,3)
      write(6,777) ((bvec(i,j),j=1,3),i=1,3)
 777  format (3f8.4) 
      stop
      end

      double precision function round(tt,mm)
      implicit double precision (a-h,o-z)
      t = tt 
      n = t 
      if( t .lt. 0.d0 ) n = n - 1
      round = t - n
      if(mm.eq.2.and. abs(round-1.d0) .lt. 1.d-6) round = 0.d0
      return
      end

       subroutine checkgr(ng,iflag)
       implicit double precision (a-h,o-z)
       common /traf/ ro(48,3,3), tnp(48,3),ntrans
       dimension rot(3,3)
       data small / 1.d-8/
c
c check if the input operation (1) constitute a lattice transformation
c                              (2) norm-conserving (proper or improper)
c
c step (1) : 
       iflag = 0
       do 135 i=1,3
       do 135 j=1,3
       rot(i,j)=ro(ng,i,j)
       kr = rot(i,j)
       if( abs( kr - rot(i,j) ) .gt. small) iflag = 1  
  135  continue
       if(iflag.eq.1) return
c step (2)
       call unita(rot,iu)
       if(iu.ne.1) iflag = 2
      return
      end
c
      subroutine trafo
      implicit double precision (a-h,o-z)
c
c     subroutine reads and write s transformation matrices t
c     and non-primitive translations tnp
c
      common /traf/ t(48,3,3), tnp(48,3),ntrans
c
 10   format (  a5, i3 )
 11   format ( 18f4.0 )
   12 format(8f10.5  )
 20   format (30x,38htransformation matrices of point group,a5,1x,22hin
     &lattice coordinates/ 29x,68(1h-) //// )
   23 format(3(7x,3f4.1//) )
   21 format(5x,i5,3f10.4/)
 24   format ( 1h1 )
 25   format(35x,49hnon-primitive translations in lattice coordinates/34
     &x,51(1h-)///)
 26   format (2x,12i10//)
 27   format ( 5x,12f10.4 )
c
      read(5,10) grponc, ntrans
      read(5,11) (((t(i,j,k),k=1,3),j=1,3),i=1,ntrans)
c      read(5,12) ( ( tnp(i,l),l=1,3),i=1,ntrans )
      iprin = 0
c      write(6,20) grponc
czw   write(6,*) ' no of operation read in from trafo =', ntrans
      if(iprin.eq.0) return
czw   write(6,25)
      do 89 mm=1,ntrans
czw   write(6,21) mm,(tnp(mm,ll),ll=1,3)
czw   write(6,23) ((t(mm,jj,kk),kk=1,3),jj=1,3)
   89 continue
      return
      end
c
c
      double precision function prosca ( a,b,gij )
      implicit double precision (a-h,o-z)
c
c     subroutine calculates scalar product in arbitrary metrik
c
      dimension a(3), gij(3,3), b(3)
      prosca = 0.d0
      do 1   j = 1,3
      do 1   i = 1,3
      prosca = prosca + a(i) * gij(i,j) * b(j)
 1    continue
      return
      end

      subroutine rla(a)
      implicit double precision (a-h,o-z)
c form the metrix matrix
      dimension a(3,3)
      common/rl/aij(3,3)
      do 10 i=1,3
      do 10 j=1,3
      aij(i,j) = 0.d0
      do 10 k=1,3
      aij(i,j) = aij(i,j) + a(i,k)*a(j,k)
 10   continue
      return
      end

      subroutine autogen
      implicit double precision (a-h,o-z)
      common/rl/aij(3,3)
      common /traf/ ro(48,3,3), tnp(48,3),ntrans
      dimension x(3),nc(3),rot(3,3),xx(3,27,3)
      data small/1.d-5/
c
c try to guess the group operations, by finding those operations that
c are "norming preserving", it will then be tested in the main routine
c to select a subset that is the symmetry group
c
c step 1 : select lattice vector transformation that preserves length
c          for the 3 lattice vectors, this is a necessary condition 
c         
c  n=1 should be  good enough if a 'sensible' unit cell is chosen
c  if funny unit cells are chosen, larger n's are needed
c
      n = 1
      do 5   nx = 1,3

      do 8  j=1,3
 8    x(j)= 0.d0
      x(nx)= 1.d0
      ddv = prosca(x,x,aij)

      ic = 0
      do 10  n1 = -n, n 
      do 10  n2 = -n, n 
      do 10  n3 = -n, n 
      x(1) = n1
      x(2) = n2
      x(3) = n3
      dd = prosca(x,x,aij)
      if( abs(dd-ddv).gt. small) goto 10
      ic = ic + 1
      do 102 j=1,3
 102  xx(nx,ic,j)= x(j)
 10   continue
      nc(nx) = ic
  5   continue
c
czw   write(6,*) ' subroutine autogen: '
czw   write(6,*) ' no of norm-conserving transformation for lat vec'
czw   write(6,*) (nc(i),i=1,3)
czw   write(6,*) ' all together ', nc(1)*nc(2)*nc(3)
c
c step 2 : construct transformation matrices from lattice vector 
c transformations found in step 1, and select those that are 'norm-conserving'
c
      ntrans = 0
      do 40 ix=1,nc(1)
      do 40 iy=1,nc(2)
      do 40 iz=1,nc(3)

      do 42 j=1,3
      rot(1,j) = xx(1,ix,j)
      rot(2,j) = xx(2,iy,j)
      rot(3,j) = xx(3,iz,j)
 42   continue
       
      call unita(rot,iu)
      if(iu.ne.1) goto 40
      ntrans = ntrans + 1
      if(ntrans.gt.48) then
czw    write(6,*) ' ntrans too large '
       stop
      endif
     
      do 44 i=1,3
      do 44 j=1,3
      ro(ntrans,i,j) = rot(i,j)
 44   continue

czw   write(6,*) ' ntrans=', ntrans
czw   write(6,610) ((rot(i,k),k=1,3),i=1,3)  
 610  format(3f10.5)

 40   continue

      return
      end
c
      subroutine unita(r,iu)
      implicit double precision (a-h,o-z)
c
c check if operation 'r' is norm-conserving
c  r A A+ r+ = A A+ , where A+ = A transpose
c where AA+ = aij = metric matrix of lattice vectors
c
      common/rl/aij(3,3)
      dimension r(3,3),ar(3,3)
      do 10 i=1,3
      do 10 j=1,3
      ar(i,j) = 0.d0
      do 20 k1=1,3
      do 20 k2=1,3
      ar(i,j) = ar(i,j) + r(i,k1) * r(j,k2) * aij(k1,k2)
 20   continue
 10   continue
      iu = 1
      do 30 i=1,3
      do 30 j=1,3
      if( abs( aij(i,j)-ar(i,j) ) .gt. 1.d-8) iu = 0
 30   continue
      return
      end
     
c - sorter
C
      subroutine tri(array,iindex,n,irank)
      implicit double precision (a-h,o-z)

c*******************************************************************************
c     Sorting of an array with the HEAPSORT algorithm (cf. W. H. Press et al., 
c     "Numerical Recipes", Cambridge(1986), p. 229 - 235.
c*******************************************************************************

      dimension array(n),iindex(n),irank(n)

c *** "indexing"
      call sindex(n,array,iindex)

c *** "ranking"
      call srank(n,iindex,irank)

      return
      end

*###############################################################################

      subroutine sindex(n,array,iindex)

      implicit double precision (a-h,o-z)
c*******************************************************************************
c     indexes an array ARRAY of length N, i. e. outputs the array IINDEX such 
c     that ARRAY(IINDEX(J)) is in ascending order for J=1,2,...,N. The input 
c     quantities N and ARRAY are not changed.
c*******************************************************************************

      dimension array(n),iindex(n)

c *** initialize the index array with consecutive integers
      do 1, j=1,n
        iindex(j) = j
    1 continue
      
ctc
      if(n.eq.1) return

c *** heapsort algorithm
      l = n/2+1
      ir = n
   10 continue
        if (l.gt.1) then
          l = l-1
          indext = iindex(l)
          q = array(indext)
        else
          indext = iindex(ir)
          q = array(indext)
          iindex(ir) = iindex(1)
          ir = ir-1
          if (ir.eq.1) then
            iindex(1) = indext
            return
          end if
        end if
        i = l
        j = l+l
   20   if (j.le.ir) then
            if (j.lt.ir) then
              if (array(iindex(j)).lt.array(iindex(j+1))) j = j+1
            end if
            if (q.lt.array(iindex(j))) then
            iindex(i) = iindex(j)
            i = j
            j = j+j
          else
            j = ir+1
          end if
          go to 20
        end if
        iindex(i) = indext
      go to 10

      end

*###############################################################################

      subroutine srank(n,iindex,irank)
      implicit double precision (a-h,o-z)

c*******************************************************************************
c     Given IINDEX of length N as output from the routine INDEX, this routine
c     returns an array IRANK, the corresponding table of ranks
c*******************************************************************************

      dimension iindex(n), irank(n)

      do 1 j=1,n
        irank(iindex(j)) = j
    1 continue
      return

      end
