CalculiX installation on gianthk Debian - 09.01.2021
install dir tree:
/usr/local/CalculiX
/usr/local/SPOOLES.2.2
/usr/local/ARPACK

- follow step by step instructions at:
http://www.dhondt.de/ccx_2.17.README.INSTALL


- SPOOLES.2.2. :
	http://www.netlib.org/linalg/spooles/spooles.2.2.html
	http://www.netlib.org/linalg/spooles/Install.ps.gz
	http://www.netlib.org/linalg/spooles/spooles.2.2.tgz

	see modified Makefile!

- ARPACK
	follow README

	See modified Makefile!
	See modified second.f in UTIL!
	https://ubuntuforums.org/showthread.php?t=1714367

	https://www.linuxquestions.org/questions/linux-software-2/install-arpack-fortran-libraries-440794/#post4855296

- CalculiX
	- use the Makefile_MT ONLY: do not compile ccx but only ccx_MT

______________________________________________________________________________________________

- other usefull links (documentation, bug-fixes...)
	- https://www.libremechanics.com/?q=node/9
	- https://github.com/calculix/ccx_fff#downloads
	- https://precice.discourse.group/t/calculix-parallelization/167/5
	- 

