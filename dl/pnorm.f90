program norm
!hmb data
implicit none

integer,parameter :: nexpt=318
integer :: i,j,k(nexpt,402), ru(nexpt,18),n
real :: ln(85),d(nexpt,402),x,y,m(nexpt,nexpt)

open(unit=85,file='allh.ph',status="old")
 do i=1,402
  read (85,*) k(:,i)
!  endif
 enddo
close(85)

do i=1,nexpt
 ru(i,1)=sum(k(i,5:19))  !cyanobacteria
 ru(i,2)=sum(k(i,20:69))+sum(k(i,72:79)) !firmicutes
 ru(i,3)=sum(k(i,82:88)) !chloroflexi
 ru(i,4)=sum(k(i,90:92)) !synergistetes 
 ru(i,5)=sum(k(i,94:95)) ! Deinococcus
 ru(i,6)=sum(k(i,97:99)) ! thermotogae
 ru(i,7)=sum(k(i,100:149)) !actinobacteria 
 ru(i,8)=sum(k(i,151:248)) !b/g proteobacteria
 ru(i,9)=sum(k(i,255:256))  !acidobacteria
 ru(i,10)=sum(k(i,258:261))  !aquificae 
 ru(i,11)=sum(k(i,265:291)) !delta / epsilon proteobacteria
 ru(i,12)=sum(k(i,293:350)) !alpha proteobacteria
 ru(i,13)=sum(k(i,352:355)) !spirochaetes
 ru(i,14)=sum(k(i,358:361)) !planctomycetes
 ru(i,15)=sum(k(i,362:368)) !chlamydiae+verrucomicrobia
 ru(i,16)=sum(k(i,370:402)) !bacteroidetes
 ru(i,17)=k(i,2)+k(i,3)+k(i,4)+k(i,70)+k(i,71)+k(i,80)+k(i,81)+k(i,89)+k(i,93)+k(i,96)+k(i,150)+sum(k(i,249:254))+k(i,257)+ &
  k(i,292)+k(i,351)+sum(k(i,356:357))+k(i,369) !multiphyla
 ru(i,18)=k(i,1) !root
enddo

!do i=1,18
! print *,ru(:,i)  !phylogenetic roll-ups
!enddo

do i=1,nexpt
 x=sum(k(i,2:402))+1
 d(i,:)=k(i,:)/x
! print *,d(i,:)
enddo

!do i=2,402
! print *,d(:,i)*1000
!enddo

do i=1,nexpt
x=0.
 do j=2,402
  y=d(i,j)
  x=x+y*y
 enddo
 d(i,:)=d(i,:)/sqrt(x)
enddo

m(:,:)=0
do i=1,nexpt
 do j=1,nexpt
  do n=2,402
   m(i,j)=d(i,n)*d(j,n)+m(i,j)
  enddo
 enddo
enddo

do i=1,nexpt
 print *,m(i,:)  !dot products
enddo

end program norm
