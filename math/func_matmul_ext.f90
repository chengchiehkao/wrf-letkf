
function matmul_ext(A,B,m,k,n) result(C)

implicit none
integer m,n,k
real(kind=8) :: A(m,k),B(k,n),C(m,n)
!================================================
call dgemm('N','N',m,n,k,1.d0,A,m,B,k,0.d0,C,m)
!================================================
return
end function matmul_ext

