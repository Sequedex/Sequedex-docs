program norm
!hmb data
implicit none

integer,parameter :: nexpt=547
integer :: i,j,n
real :: k(nexpt,964), ru(nexpt,30)
real :: ln(85),d(nexpt,964),x,y,m(nexpt,nexpt)

open(unit=85,file='what',status="old")
 do i=1,964
  read (85,*) k(:,i)
!  endif
 enddo
close(85)

do i=1,nexpt
 ru(i,1)=sum(k(i,1:54))  !amino acids and derivatives
 ru(i,2)=sum(k(i,55:153)) !carbohydrates
 ru(i,3)=sum(k(i,154:162)) !cell division and cell cycle
 ru(i,4)=sum(k(i,163:213)) !cell wall and capsule 
 ru(i,5)=sum(k(i,214:357)) ! clustering-based
 ru(i,6)=sum(k(i,358:397)) !cofactors, vitamins, prosthetic groups, pigments 
 ru(i,7)=sum(k(i,398:431)) !dna metabolism 
 ru(i,8)=sum(k(i,432:445)) !dormancy and sporulation
 ru(i,9)=sum(k(i,446:464))  !fatty acids, lipids, and isoprenoids
 ru(i,10)=sum(k(i,465:488))  !Iron acquisition and metabolism 
 ru(i,11)=sum(k(i,489:534)) !membrane transport
 ru(i,12)=sum(k(i,535:564)) !metablism of aeromatic compounds
 ru(i,13)=sum(k(i,565:579)) !miscellaneous
 ru(i,14)=sum(k(i,580:587)) !motility and chemotaxis
 ru(i,15)=sum(k(i,588:597)) !nitrogen metabolism
 ru(i,16)=sum(k(i,598:615)) !nucleosides and nucleotides
 ru(i,17)=sum(k(i,616:637)) !phages, prophages, transposable elements
 ru(i,18)=sum(k(i,638:643)) !phosphorus metabolism
 ru(i,19)=sum(k(i,644:652)) !photosynthesis
 ru(i,20)=sum(k(i,653:656)) !potassium metabolism
 ru(i,21)=sum(k(i,657:722)) !protein metabolism
 ru(i,22)=sum(k(i,723:753)) !RNA metabolism
 ru(i,23)=sum(k(i,754:791)) !regulation and cell signaling
 ru(i,24)=sum(k(i,792:827)) !respiration
 ru(i,25)=sum(k(i,828:843)) !secondary metabolism
 ru(i,26)=sum(k(i,844:881)) !stress response
 ru(i,27)=sum(k(i,882:894)) !sulfur metabolism
 ru(i,28)=sum(k(i,895:962)) !virulence
 ru(i,29)=k(i,963) !ribosome
 ru(i,30)=k(i,964) !no match
enddo

open(unit=85,file='froll',status="replace")
do i=1,30
  write(85,*) ru(:,i)  ! roll-ups
enddo
close(85)

do i=1,nexpt
 x=sum(k(i,1:963))+1
 d(i,:)=k(i,:)/x
! print *,d(i,:)
enddo

open(unit=85,file='fnorm',status="replace")
do i=1,963
  write(85,*) d(:,i)*1000
enddo
close(85)

do i=1,nexpt
x=0.
 do j=1,963
  y=d(i,j)
  x=x+y*y
 enddo
 d(i,:)=d(i,:)/sqrt(x)
enddo

m(:,:)=0
do i=1,nexpt
 do j=1,nexpt
  do n=1,963
   m(i,j)=d(i,n)*d(j,n)+m(i,j)
  enddo
 enddo
enddo

open(unit=85,file='fdot',status="replace")
do i=1,nexpt
  write(85,*) m(i,:)  !dot products
enddo
close(85)

end program norm
