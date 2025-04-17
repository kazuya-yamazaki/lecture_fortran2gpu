!C
!C***
!C*** INPUT_GRID
!C***
!C
      subroutine INPUT_GRID (OLDtoNEWpe, ITERkk)
      use hpcmw_all

      integer(kind=kint ), parameter  ::  l_err = 6

      integer, parameter   ::  ndepth = 1
      integer(kind=kint ), parameter  ::  neibpetot_max = 26

      integer(kind=kint )        ::  nx_all , ny_all , nz_all
      integer(kind=kint ), save  ::  npp, ndx, ndy, ndz           
      integer(kind=kint ), save  ::  ITERmax

      integer(kind=kint )  ::  nxi, nyi, nzi  
      integer(kind=kint )  ::  nx , ny , nz 

      integer(kind=kint )  ::  ipe    , jpe    , kpe    , pe_id
      integer(kind=kint )  ::  inp    , jnp    , knp
      integer(kind=kint )  ::  inp_st , jnp_st , knp_st ,               &
     &                         inp_end, jnp_end, knp_end

      integer(kind=kint )  ::  ioff   , joff   , koff
      integer(kind=kint )  ::  i      , j      , k  
      integer(kind=kint )  ::  i_st   , j_st   , k_st   ,               &
     &                         i_end  , j_end  , k_end
      integer(kind=kint )  ::  is     , js     , ks     ,               &
     &                         ie     , je     , ke   

      integer(kind=kint )  ::  nodtot,   intnodtot, elmtot, enod
      integer(kind=kint )  ::  item_pos, item_tot

      character(len=12 )   ::  WORKFIL

      integer(kind=kint ), dimension(:,:,:), allocatable  :: node_id_lc
      integer(kind=kint ), dimension(:,:,:), allocatable  :: node_id_gl
      integer(kind=kint ), dimension(:), allocatable  :: Pindex

      integer(kind=kint )  ::  node_id, element_id, element_id_gl

      integer(kind=kint), dimension(:), allocatable :: nWORK

      integer(kind=kint )  ::  elmgrptot, nodgrptot, sufgrptot

      real   (kind=kreal)  ::  forc
      integer(kind=kint )  ::  i1 , i2, i3, i4, i5, i6, i7, i8

      real   (kind=kreal)                :: eps         = 1d-7
      integer(kind=kint )                :: max_itr   
      real   (kind=kreal)                :: large_value = 1.0e30_kreal
      integer(kind=kint )                :: nparm       =   3
      real   (kind=kreal)                :: p_corner, p_line, p_mid

      integer(kind=kint ), parameter            ::  nindex = 6

      integer :: CONFIGopt
      integer, dimension(:,:), allocatable :: NEIBregion
      integer, dimension(:  ), allocatable :: NEIBregionNUM
      integer, dimension(:,:), allocatable :: IWKG

      integer, dimension(PETOT) :: OLDtoNEWpe

      write(WORKFIL,'(a,i6.6)') 'wk.p_', my_rank


      if (ITERkk.eq.1) then
      if (my_rank.eq.0) then
        open (11, file='mesh.inp',status='unknown', form='formatted')
        read (11,*)  npx, npy, npz
        read (11,*)  ndx, ndy, ndz
        read (11,*)  PEsmpTOT
        read (11,*)  ITERmax
        read (11,*)  NCOLORk
        write (*,'(3i10)') npx, npy, npz
        write (*,'(3i10)') ndx, ndy, ndz
        write (*,'( i10)') PEsmpTOT
        write (*,'( i10)') NCOLORk
        close (11)
      endif

      call mpi_bcast (npx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (npy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (npz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (ndx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (ndy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (ndz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (PEsmpTOT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (ITERmax , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (NCOLORk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      endif      

      if (allocated(IWKG))       deallocate (IWKG)
      if (allocated(NEIBregion)) deallocate (NEIBregion)

      if (allocated(node_id_lc)) deallocate (node_id_lc) 
      if (allocated(node_id_gl)) deallocate (node_id_gl) 
      if (allocated(Pindex))     deallocate (Pindex)

      NEIBPE= 0
      NOD_STACK_IMPORT= 0
      NOD_STACK_EXPORT= 0
      if (.not. allocated(NOD_IMPORT)) then
        allocate(NOD_IMPORT(300000*4))
        allocate(NOD_EXPORT(300000*4))
        allocate(NOD_EXPORT_NEW(300000*4))
      endif
      NOD_IMPORT      = 0
      NOD_EXPORT      = 0

      if (allocated (ICELNOD)) deallocate (ICELNOD)
      if (allocated (BC_STACK)) deallocate (BC_STACK)
      if (allocated (BC_NOD  )) deallocate (BC_NOD)

      nxi= npx/ndx
      nyi= npy/ndy
      nzi= npz/ndz
 
      if (mod(ITERkk,3).eq.1) CONFIGopt= 1
      if (mod(ITERkk,3).eq.2) CONFIGopt= 2
      if (mod(ITERkk,3).eq.0) CONFIGopt= 3

      CONFIGopt=1

      nx_all= ndx * nxi
      ny_all= ndy * nyi
      nz_all= ndz * nzi

      hpcmwIarray(1)= ITERmax

!C
!C +---------------+
!C | PE reordering |
!C +---------------+
!C===
      OLDtoNEWpe= 0

      if (CONFIGopt.eq.1) then
!C
!C== NORMAL
      if (my_rank.eq.0) write (*,'(a)') '### NORMAL'
      do ip= 1, PEtot
        OLDtoNEWpe(ip)= ip
      enddo
!C==
      endif

      if (CONFIGopt.eq.2) then
!C
!C== RED-BLACK
      if (my_rank.eq.0) write (*,'(a)') '### RED-and-BLACK'
      icou = 0
      icouG= 0
      do k= 1, ndz
        do j= 1, ndy
          do i= 1, ndx
            icou= icou + 1
            imod= mod(i+j+k,2)
            if (imod.eq.1) then
              icouG= icouG + 1
              OLDtoNEWpe(icou)= icouG
            endif
          enddo
        enddo
      enddo

      icou= 0
      do k= 1, ndz
        do j= 1, ndy
          do i= 1, ndx
            icou= icou + 1
            imod= mod(i+j+k,2)
            if (imod.eq.0) then
              icouG= icouG + 1
              OLDtoNEWpe(icou)= icouG
            endif
          enddo
        enddo
      enddo
!C==
      endif

      if (CONFIGopt.eq.3) then
!C
!C== HYPER-PLANE
      if (my_rank.eq.0) write (*,'(a)') '### HYPER-PLANE'

      if (allocated(NEIBregion)) deallocate (NEIBregion)
      if (allocated(IWKG))       deallocate (IWKG)
      allocate (NEIBregion(PETOT,6))
      allocate (IWKG(PETOT,6))

      NEIBregion= 0
      IWKG      = 0

      icou = 0
      do k= 1, ndz
        do j= 1, ndy
          do i= 1, ndx
            icou= icou + 1
           
            NEIBregion(icou,1)= icou - 1
            NEIBregion(icou,2)= icou + 1
            NEIBregion(icou,3)= icou - ndx
            NEIBregion(icou,4)= icou + ndx
            NEIBregion(icou,5)= icou - ndx*ndy
            NEIBregion(icou,6)= icou + ndx*ndy

            if (i.eq.1  ) NEIBregion(icou,1)= 0
            if (i.eq.ndx) NEIBregion(icou,2)= 0
            if (j.eq.1  ) NEIBregion(icou,3)= 0
            if (j.eq.ndy) NEIBregion(icou,4)= 0
            if (k.eq.1  ) NEIBregion(icou,5)= 0
            if (k.eq.ndz) NEIBregion(icou,6)= 0
          enddo
        enddo
      enddo

      icouG= 1
      IWKG(1,1)= 1
      IWKG(1,2)= 1
      do icol= 2, PETOT
        do ip= 1, PETOT
          IWKG(ip,3)= 0
        enddo
        do ip= 1, PETOT
          if (IWKG(ip,1).eq.0) then
            do k= 1, 6
              ic= NEIBregion(ip,k)
              if (ic.ne.0) then
                if (IWKG(ic,1).eq.icol) goto 100
              endif
            enddo

            do k= 1, 6
              ic= NEIBregion(ip,k)
              if (ic.ne.0) then
                if (IWKG(ic,1).eq.icol-1) then
                  icouG= icouG + 1
                  IWKG(ip,1)= icol
                  IWKG(ip,2)= icouG
                  if (icouG.eq.PEtot) goto 200
                  do kk= 1, 6
                    ick= NEIBregion(ip,k)                   
                    if (ick.ne.0) IWKG(ick,3)= 1
                  enddo
                  goto 100
                endif
              endif
            enddo
          endif
 100            continue
        enddo
      enddo
 200    continue
      do ip= 1, PEtot
        OLDtoNEWpe(ip)= IWKG(ip,2)
      enddo
      deallocate (NEIBregion, IWKG)
!C==
      endif

!
! ***** set internal node count
!

      nxi = nx_all / ndx
      nyi = ny_all / ndy
      nzi = nz_all / ndz
!
! ***** allocate nodal id table
      if (allocated(node_id_lc)) deallocate (node_id_lc)
      if (allocated(node_id_gl)) deallocate (node_id_gl)
      allocate ( node_id_lc(nxi+2*ndepth,nyi+2*ndepth,nzi+2*ndepth) )
      allocate ( node_id_gl(nxi+2*ndepth,nyi+2*ndepth,nzi+2*ndepth) )
!
! **********   domain loop for each pe   **********
!
      pe_id = 1

!      write (*,*) 'INITIAL ID ?'
!      read  (*,*) INIT_ID 
      INIT_ID= 0
      INIT_ID= INIT_ID + 1


      nnp= ndz * ndy * ndx
      if (allocated(Pindex)) deallocate (Pindex)
      allocate (Pindex(0:nnp))
      Pindex= 0

      nnpp= nnp / PETOT
      nr  = nnp - nnpp * PETOT
      do i= 1, PETOT
        Pindex(i)= nnpp
        if (i.le.nr) then
          Pindex(i)= nnpp + 1
        endif
      enddo

      do i= 1, PETOT
        Pindex(i)= Pindex(i-1) + Pindex(i)
      enddo

      do kpe=1,ndz
        do jpe=1,ndy
          do ipe=1,ndx
!
! ***** open output file
!
             if (OLDtoNEWpe(pe_id).ge.Pindex(my_rank)+1 .and.           &
     &           OLDtoNEWpe(pe_id).le.Pindex(my_rank+1)) then

!
! ***** set and write basic local model parameters
!                                       .. pe nod per 1 line
                             nx = nxi
            if (ipe /=   1)  nx = nx + ndepth
            if (ipe /= ndx)  nx = nx + ndepth

                             ny = nyi
            if (jpe /=   1)  ny = ny + ndepth
            if (jpe /= ndy)  ny = ny + ndepth

                             nz = nzi
            if (kpe /=   1)  nz = nz + ndepth
            if (kpe /= ndz)  nz = nz + ndepth

!                                       .. total node
            nodtot    = nx  * ny  * nz
!                                       .. internal nodes
            intnodtot = nxi * nyi * nzi
!                                       .. search neighbor pe
                            inp_st  = -1
            if (ipe ==   1) inp_st  =  0
                            inp_end =  1
            if (ipe == ndx) inp_end =  0

                            jnp_st  = -1
            if (jpe ==   1) jnp_st  =  0
                            jnp_end =  1
            if (jpe == ndy) jnp_end =  0

                            knp_st  = -1
            if (kpe ==   1) knp_st  =  0
                            knp_end =  1
            if (kpe == ndz) knp_end =  0
!                                       .. set neighbor pe
            neibpetot = 0

            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

                  neibpetot = neibpetot  + 1
                  neibpe     (neibpetot) =  pe_id    +                  &
     &                                  inp + jnp*ndx + knp*ndx*ndy
                enddo
              enddo
            enddo

            do neib= 1, NEIBPETOT
              ip= NEIBPE(neib)
              NEIBPE(neib)= OLDtoNEWpe(ip) - 1
            enddo

            NP= nodtot
            N = intnodtot

!
! ***** set coordinate off set (starting corner for pe node)
!

!
! ***** set nodal position off set (i,j,k starting position -1)
!
                        ioff = (ipe-1)*nxi
            if (ipe/=1) ioff = ioff - ndepth
                        joff = (jpe-1)*nyi
            if (jpe/=1) joff = joff - ndepth
                        koff = (kpe-1)*nzi
            if (kpe/=1) koff = koff - ndepth
!
! ***** set and write coordinate for internal nodes
!
            node_id = 0
                            i_st  =     1         + ndepth
            if (ipe ==   1) i_st  =     1
                            i_end =  i_st + nxi-1

                            j_st  =     1         + ndepth
            if (jpe ==   1) j_st  =     1
                            j_end =  j_st + nyi-1

                            k_st  =     1         + ndepth
            if (kpe ==   1) k_st  =     1
                            k_end =  k_st + nzi-1

            do k=k_st,k_end
              do j=j_st,j_end
                do i=i_st,i_end

                  node_id = node_id + 1

                  node_id_lc(i,j,k) =  node_id
                  node_id_gl(i,j,k) = (ioff+i  ) +                      &
     &                                (joff+j-1)*nx_all +               &
     &                                (koff+k-1)*nx_all*ny_all 
                enddo
              enddo
            enddo
!
! ***** set and write coordinate for sleeve area nodes
!
!!          write(l_out,'(a)') ' '

            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=      1
                  if ( inp == -1 )  ie=      1
                  if ( jnp == -1 )  js=      1
                  if ( jnp == -1 )  je=      1
                  if ( knp == -1 )  ks=      1
                  if ( knp == -1 )  ke=      1

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end+1
                  if ( inp ==  1 )  ie=i_end+1
                  if ( jnp ==  1 )  js=j_end+1
                  if ( jnp ==  1 )  je=j_end+1
                  if ( knp ==  1 )  ks=k_end+1
                  if ( knp ==  1 )  ke=k_end+1

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie

                      node_id = node_id + 1

                      node_id_lc(i,j,k) =  node_id
                      node_id_gl(i,j,k) = (ioff+i  ) +                  &
     &                                    (joff+j-1)*nx_all +           &
     &                                    (koff+k-1)*nx_all*ny_all 
                    enddo
                   enddo
                  enddo

                enddo
              enddo
            enddo
!
! ..... write 2.2 element (connection)
!

            elmtot = (nx-1)*(ny-1)*(nz-1)
            enod   =  8
          
            ICELTOT= elmtot
            if (allocated(ICELNOD)) deallocate (ICELNOD)            
            allocate (ICELNOD(ICELTOT,8))
            element_id = 0

            do k=1,nz-1
              do j=1,ny-1
                do i=1,nx-1

                  element_id    =  element_id + 1
                  element_id_gl =  (ioff+i  ) +                         &
     &                             (joff+j-1)*(nx_all-1) +              &
     &                             (koff+k-1)*(nx_all-1)*(ny_all-1)

                  i1 = node_id_lc( i  , j  , k   )
                  i2 = node_id_lc( i+1, j  , k   )
                  i3 = node_id_lc( i  , j+1, k   )
                  i4 = node_id_lc( i+1, j+1, k   )
                  i5 = node_id_lc( i  , j  , k+1 )
                  i6 = node_id_lc( i+1, j  , k+1 )
                  i7 = node_id_lc( i  , j+1, k+1 )
                  i8 = node_id_lc( i+1, j+1, k+1 )

                  ICELNOD(element_id,1)= i1
                  ICELNOD(element_id,2)= i2
                  ICELNOD(element_id,3)= i4
                  ICELNOD(element_id,4)= i3
                  ICELNOD(element_id,5)= i5
                  ICELNOD(element_id,6)= i6
                  ICELNOD(element_id,7)= i8
                  ICELNOD(element_id,8)= i7

                enddo
              enddo
            enddo
!
! ..... write 3.import / export information
!
!
! ***** set and write import nodes
!                                     .... count nodes 

            NOD_STACK_IMPORT= 0
            NOD_STACK_EXPORT= 0

            neibpetot = 0
            node_id   = 0

            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=      1
                  if ( inp == -1 )  ie=      1
                  if ( jnp == -1 )  js=      1
                  if ( jnp == -1 )  je=      1
                  if ( knp == -1 )  ks=      1
                  if ( knp == -1 )  ke=      1

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end+1
                  if ( inp ==  1 )  ie=i_end+1
                  if ( jnp ==  1 )  js=j_end+1
                  if ( jnp ==  1 )  je=j_end+1
                  if ( knp ==  1 )  ks=k_end+1
                  if ( knp ==  1 )  ke=k_end+1

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      node_id = node_id + 1
                    enddo
                   enddo
                  enddo

                  neibpetot = neibpetot  + 1
                  NOD_STACK_IMPORT(neibpetot) = node_id

                enddo
              enddo
            enddo

!                                     .... write nodes 
            neibpetot = 0
            icou= 0
            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=      1
                  if ( inp == -1 )  ie=      1
                  if ( jnp == -1 )  js=      1
                  if ( jnp == -1 )  je=      1
                  if ( knp == -1 )  ks=      1
                  if ( knp == -1 )  ke=      1

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end+1
                  if ( inp ==  1 )  ie=i_end+1
                  if ( jnp ==  1 )  js=j_end+1
                  if ( jnp ==  1 )  je=j_end+1
                  if ( knp ==  1 )  ks=k_end+1
                  if ( knp ==  1 )  ke=k_end+1

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  neibpetot = neibpetot  + 1

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      icou= icou + 1
                    enddo
                   enddo
                  enddo

                enddo
              enddo
            enddo


           if (allocated(nWORK)) deallocate (nWORK)
           allocate (nWORK(icou))

            neibpetot= 0
            icou     = 0
            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=      1
                  if ( inp == -1 )  ie=      1
                  if ( jnp == -1 )  js=      1
                  if ( jnp == -1 )  je=      1
                  if ( knp == -1 )  ks=      1
                  if ( knp == -1 )  ke=      1

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end+1
                  if ( inp ==  1 )  ie=i_end+1
                  if ( jnp ==  1 )  js=j_end+1
                  if ( jnp ==  1 )  je=j_end+1
                  if ( knp ==  1 )  ks=k_end+1
                  if ( knp ==  1 )  ke=k_end+1

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  neibpetot = neibpetot  + 1

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      icou= icou + 1
                      nWORK(icou)= node_id_lc(i,j,k)
                    enddo
                   enddo
                  enddo

                enddo
              enddo
            enddo

            if (neibpetot.ne.0) then
              nn= NOD_STACK_IMPORT(neibpetot)               
              do kk= 1, nn
                NOD_IMPORT(kk)= nWORK(kk)
              enddo
              deallocate (nWORK) 
            endif
!
! ***** set and write export nodes
!                                     .... count nodes 

            neibpetot = 0
            node_id   = 0

            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=i_st
                  if ( inp == -1 )  ie=i_st
                  if ( jnp == -1 )  js=j_st
                  if ( jnp == -1 )  je=j_st
                  if ( knp == -1 )  ks=k_st
                  if ( knp == -1 )  ke=k_st

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end
                  if ( inp ==  1 )  ie=i_end
                  if ( jnp ==  1 )  js=j_end
                  if ( jnp ==  1 )  je=j_end
                  if ( knp ==  1 )  ks=k_end
                  if ( knp ==  1 )  ke=k_end

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      node_id = node_id + 1
                    enddo
                   enddo
                  enddo

                  neibpetot = neibpetot  + 1
                  NOD_STACK_EXPORT(neibpetot) = node_id

                enddo
              enddo
            enddo

!                                     .... write nodes 
            neibpetot = 0
            icou      = 0
            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=i_st
                  if ( inp == -1 )  ie=i_st
                  if ( jnp == -1 )  js=j_st
                  if ( jnp == -1 )  je=j_st
                  if ( knp == -1 )  ks=k_st
                  if ( knp == -1 )  ke=k_st

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end
                  if ( inp ==  1 )  ie=i_end
                  if ( jnp ==  1 )  js=j_end
                  if ( jnp ==  1 )  je=j_end
                  if ( knp ==  1 )  ks=k_end
                  if ( knp ==  1 )  ke=k_end

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  neibpetot = neibpetot  + 1

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      icou= icou + 1
                    enddo
                   enddo
                  enddo

                enddo
              enddo
            enddo

            neibpetot = 0
            if (allocated(nWORK)) deallocate (nWORK)
            allocate (nWORK(icou))
            icou      = 0
            do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
                do inp=inp_st,inp_end

                  if ((inp==0).and.(jnp==0).and.(knp==0)) cycle

!                                       .. start side

                  if ( inp == -1 )  is=i_st
                  if ( inp == -1 )  ie=i_st
                  if ( jnp == -1 )  js=j_st
                  if ( jnp == -1 )  je=j_st
                  if ( knp == -1 )  ks=k_st
                  if ( knp == -1 )  ke=k_st

!                                       .. finish side

                  if ( inp ==  1 )  is=i_end
                  if ( inp ==  1 )  ie=i_end
                  if ( jnp ==  1 )  js=j_end
                  if ( jnp ==  1 )  je=j_end
                  if ( knp ==  1 )  ks=k_end
                  if ( knp ==  1 )  ke=k_end

!                                       .. line pattern

                  if ( inp ==  0 )  is=i_st
                  if ( inp ==  0 )  ie=i_end
                  if ( jnp ==  0 )  js=j_st
                  if ( jnp ==  0 )  je=j_end
                  if ( knp ==  0 )  ks=k_st
                  if ( knp ==  0 )  ke=k_end

                  neibpetot = neibpetot  + 1

                  do k=ks,ke
                   do j=js,je
                    do i=is,ie
                      icou= icou + 1
                      nWORK(icou)= node_id_lc(i,j,k)
                    enddo
                   enddo
                  enddo

                enddo
              enddo
            enddo

            nn= NOD_STACK_EXPORT(neibpetot)
            if (neibpetot.ne.0) then
              do kk= 1, nn
                NOD_EXPORT(kk)= nWORK(kk)
              enddo
              deallocate (nWORK) 
            endif
!
! ..... write 4.group information
!

! 
! ***** write boundary condition (x,y,z=0 plane sym., x-force)
!
            elmgrptot = 0
            nodgrptot = 0
            sufgrptot = 0
!                                       ... node    group

!                                        .. count node group and stack

            nodgrptot = 4
            if (allocated(BC_STACK)) deallocate (BC_STACK)            
            allocate (BC_STACK(0:nodgrptot))
            BC_STACK= 0

!                                        .. count nodal item length
!                                                 .. XMIN
            item_tot = 0
            item_pos = 0

            item_pos= 1
            BC_STACK(item_pos)= BC_STACK(item_pos-1)
            if (ipe == 1) then 
              item_tot = item_tot +  ny*nz
              BC_STACK(item_pos) = item_tot
            endif
!                                                 .. YMIN
            item_pos= 2
            BC_STACK(item_pos)= BC_STACK(item_pos-1)
            if (jpe == 1) then 
              item_tot = item_tot +  nx*nz
              BC_STACK(item_pos) = item_tot
            endif
!                                                 .. ZMIN
            item_pos= 3
            BC_STACK(item_pos)= BC_STACK(item_pos-1)
            if (kpe == 1) then 
              item_tot = item_tot +  nx*ny
              BC_STACK(item_pos) = item_tot
            endif

!                                                 .. ZMAX
            item_pos= 4
            BC_STACK(item_pos)= BC_STACK(item_pos-1)
            if (kpe == ndz) then 
              item_tot = item_tot +  nx*ny
              BC_STACK(item_pos) = item_tot
            endif

            item_pos = 0

            nn= BC_STACK(4)
            if (allocated(BC_NOD)) deallocate (BC_NOD)            
            allocate (BC_NOD(nn))

!                                                 .. XMIN
            if (ipe == 1) then 
              icc= 0
              do k= 1, nz
                do j= 1, ny
                  icc= icc + 1
                  in = BC_STACK(0) + icc
                  BC_NOD(in)= node_id_lc(1,j,k)
                enddo
              enddo
            endif
!                                                 .. YMIN
            if (jpe == 1) then 
              icc= 0
              do k= 1, nz
                do i= 1, nx
                  icc= icc + 1
                  in = BC_STACK(1) + icc
                  BC_NOD(in)= node_id_lc(i,1,k)
                enddo
              enddo
            endif
!                                                 .. ZMIN
            if (kpe == 1) then 
              icc= 0
              do j= 1, ny
                do i= 1, nx
                  icc= icc + 1
                  in = BC_STACK(2) + icc
                  BC_NOD(in)= node_id_lc(i,j,1)
                enddo
              enddo
            endif
!                                                 .. ZMAX
            if (kpe == ndz) then 
              icc= 0
              do j= 1, ny
                do i= 1, nx
                  icc= icc + 1
                  in = BC_STACK(3) + icc
                  BC_NOD(in)= node_id_lc(i,j,nz)
                enddo
              enddo
            endif

            endif

            pe_id = pe_id + 1


          enddo
        enddo
      enddo

      return
      end
