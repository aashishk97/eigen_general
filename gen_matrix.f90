program matrix
      implicit none
      complex*16, allocatable, dimension(:,:) :: A, C, U, U_inv, D, E, G, Q
      integer :: i,j,m,n
      real*8 :: r1,r2
      real*8, allocatable :: eigv(:)

      write(*,*)'Enter the size you want for square matrix A :'
      read(*,*)m
      n=m
      allocate (A(m,n))
      allocate (C(m,n))
      allocate (U(m,n))
      allocate (U_inv(n,m))
      allocate (eigv(m))

!###################...Define the Matrix A...############################          
      write(*,*)''                                                             
      write(*,*)'For making the matrix A.......'                                
      write(*,*)''                                                              
      A = dcmplx(0.d0,0.d0)
      call random_seed()
      do j=1,n
         do i=j,m
            !write(*,*)'Enter the complex elements of matrix A one by one',i,j 
            !read(*,*) r1,r2                                                     
            ! r1=rand()
            ! r2=rand()
            call random_number(r1)
            call random_number(r2)
            if(i.eq.j) then
                A(i,j)=dcmplx(r1,0.d0)                                                
            else
                A(i,j)=dcmplx(r1,r2)                                                
                A(j,i)=dconjg(A(i,j))
            end if
         end do                                                                 
      end do
      write(*,*)''
      write(*,*)'matrix A :'
      do i=1,m
           write(*,*)A(i,:)
           C(i,:)= A(i,:)
      end do

      open(10,file='original_matrix',status='unknown')
      open(11,file='unitary_matrix',status='unknown')
      open(12,file='unitary_inv_matrix',status='unknown')

      do i=1,m
           write(10,*) C(i,:)
      enddo
      close(10)
      write(*,*)''
!###################...Diagonalization...################################

      call diagonalize(n, A, eigv)

      write(*,*)'The eigenvalues computed are...'
      write(*,*) eigv(1:m) 
      write(*,*)'matrix U is...'
      do i=1,m
      write(*,*)A(i,:)  !Here the original matrix A will be over written into Unitary matrix when JOBZ='V'.
      U(i,:)=A(i,:)
      enddo

      do i=1,m
      write(11,*) U(i,:)
      enddo
      write(*,*)''
!###################...Matrix Inversion...###############################
      call cmatinv(n, A)
      write(*,*)'matrix U_inv :'
      do i=1,m
         write(*,*) A(i,:)
         U_inv(i,:)= A(i,:)
      enddo
      write(*,*)''

      do i=1,m
         write(12,*) U_inv(i,:)
      enddo
      write(*,*)''
!#######...Matrix multiplication between U, C & U_inv...#################
      E=MATMUL(C, U)
      D=MATMUL(U_inv, E)
      write(*,*)'diagonal matrix D :'
      write(*,*)''
      do i=1,m
         write(*,*) D(i,:)
      end do
      write(*,*)''
      G=MATMUL(U, U_inv)
      do i=1,m
      write(*,*)G(i,:) 
      end do
      write(*,*)''
      close(11)
      close(12)
end program matrix
!########################################################################
subroutine diagonalize(N, A, evals)

character*1 :: JOBZ, UPLO
integer :: N, LDA, LWORK, INFO
real*8 :: evals(N), RWORK(3*N-2)
complex*16 :: A(N,N), WORK(2*N-1)

JOBZ='V'
UPLO='L'
LDA=N
LWORK=2*N-1


call ZHEEV(JOBZ, UPLO, N, A, LDA, evals, WORK, LWORK, RWORK, INFO)

write(6,*) 'INFO=',INFO

return
end subroutine diagonalize
!########################################################################
subroutine cmatinv(n,a)
  implicit none
  integer :: n,INFO,NRHS,LDA,LDB,i,j
  complex*16 :: a(n,n),x(n,n),b(n,n)
  integer :: IPIV(n)
  complex*16,parameter :: zero=dcmplx(0.d0,0.d0),one=cmplx(1.d0,0.d0)

  NRHS=N

  LDA=N

  B=A

  LDB=N

  do i=1,n
    do j=1,n
       A(i,j)=zero
    end do
    A(i,i)=one
  end do

  call ZGESV( N, NRHS, B, LDA, IPIV, A, LDB, INFO )

return
end subroutine cmatinv
!########################################################################
