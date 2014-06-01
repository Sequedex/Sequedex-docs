program norm
!hmb data
implicit none

integer,parameter :: nexpt=547
integer :: i,j,k(nexpt,2549), ru(nexpt,40),n
real :: ln(85),d(nexpt,2549),x,y,m(nexpt,nexpt)

open(unit=85,file='who',status="old")
 do i=1,2549
  read (85,*) k(:,i)
!  endif
 enddo
close(85)

do i=1,nexpt
ru(i,1)=k(i,7) !Elusimicrobia
ru(i,2)=sum(k(i,9:25)) !Fusobacteriales
ru(i,3)=sum(k(i,26:75)) !Spirochaetes
ru(i,4)=sum(k(i,78:89)) !Chlamydiae
ru(i,5)=sum(k(i,91:99)) !Verrucomicrobia
ru(i,6)=sum(k(i,100:111)) !Planctomycetales
ru(i,7)=sum(k(i,116:125)) !Chlorobi
ru(i,8)=sum(k(i,128:300)) !Bacteroidetes
ru(i,9)=sum(k(i,306:341)) !Epsilon_proteobacteria
ru(i,10)=sum(k(i,343:353)) !Aquificae
ru(i,11)=sum(k(i,356:363)) !Acidobacteria
ru(i,12)=sum(k(i,365:367)) !Deferribacteres
ru(i,13)=sum(k(i,369:427)) !Delta_proteobacteria
ru(i,14)=sum(k(i,429:691)) !Alpha_proteobacteria
ru(i,15)=sum(k(i,693:1083)) !Gamma_proteobacteria
ru(i,16)=sum(k(i,1084:1219)) !Beta_proteobacteria
ru(i,17)=sum(k(i,1222:1511)) !Actinobacteria
ru(i,18)=sum(k(i,1518:1521)) !Nitrospirae
ru(i,19)=k(i,1523) !Thermodesulfobacteria
ru(i,20)=sum(k(i,1525:1541)) !Deinococcus-Thermus
ru(i,21)=sum(k(i,1548:1562)) !Thermotogae
ru(i,22)=sum(k(i,1563:1572)) !Synergistetes
ru(i,23)=sum(k(i,1573:1590)) !Chloroflexi
ru(i,24)=sum(k(i,1592:1820)) !Clostridia
ru(i,25)=sum(k(i,1822:2065)) !Bacilli
ru(i,26)=sum(k(i,2067:2081)) !Erysipelotrichi
ru(i,27)=sum(k(i,2082:2120)) !Mollicutes
ru(i,28)=sum(k(i,2121:2213)) !Cyanobacteria
ru(i,29)=sum(k(i,2215:2219)) !Archaea
ru(i,30)=sum(k(i,2220:2253)) !Crenarchaeota
ru(i,31)=sum(k(i,2254:2259)) !Thaumarchaeota
ru(i,32)=sum(k(i,2260:2401)) !Euryarchaeota
ru(i,33)=k(i,2402) !Eukaryotes
ru(i,34)=sum(k(i,2403:2448)) !Protozoa
ru(i,35)=sum(k(i,2449:2482)) !Metazoa
ru(i,36)=sum(k(i,2483:2549)) !Fungi

ru(i,37)=sum(k(i,3:6))+k(i,8)+sum(k(i,76:77))+k(i,90)+sum(k(i,112:115))+sum(k(i,126:127))+sum(k(i,301:302))+ &
  k(i,342)+sum(k(i,354:355))+k(i,364)+k(i,368)+k(i,428)+k(i,692) !Gram-negative
ru(i,38)=sum(k(i,1220:1221))+sum(k(i,1512:1517))+k(i,1522)+k(i,1524)+sum(k(i,1542:1547))+k(i,1591)+k(i,1821)+k(i,2066)  !Gram-positive
ru(i,39)=k(i,2) !Bacteria
ru(i,40)=k(i,1)+k(i,2214) !Root
enddo

open(unit=85,file='proll',status="replace")
do i=1,40
 write(85,*) ru(:,i)  !phylogenetic roll-ups
enddo
close(85)

do i=1,nexpt
 x=sum(k(i,3:2549))+1
 d(i,:)=k(i,:)/x
! print *,d(i,:)
enddo

open(unit=85,file='pnorm',status="replace")
do i=3,2549
 write(85,*) d(:,i)*1000
enddo
close(85)

do i=1,nexpt
x=0.
 do j=3,2549
  y=d(i,j)
  x=x+y*y
 enddo
 d(i,:)=d(i,:)/sqrt(x)
enddo

m(:,:)=0
do i=1,nexpt
 do j=1,nexpt
  do n=3,2549
   m(i,j)=d(i,n)*d(j,n)+m(i,j)
  enddo
 enddo
enddo

open(unit=85,file='pdot',status="replace")
do i=1,nexpt
 write(85,*) m(i,:)  !dot products
enddo
close(85)

end program norm
