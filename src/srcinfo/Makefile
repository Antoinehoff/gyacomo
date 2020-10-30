#Default and consistent setings for finding libraries and modules
include ../../local/dirs.inc
include ../../local/make.inc

#VERSION=`svn info $(SRCDIR) | grep Revision | awk '{print $$2}'`
VERSION=`git describe --always --dirty`
BRANCH=`git rev-parse --abbrev-ref HEAD`
AUTHOR=`whoami`
DATE=`date`
HOST=`hostname`

all: srcinfo

srcinfo:
	@echo 'character(len=40) VERSION' | cat > srcinfo.h
	@echo 'character(len=40) BRANCH' | cat >> srcinfo.h
	@echo 'character(len=20) AUTHOR'| cat >> srcinfo.h
	@echo 'character(len=40) EXECDATE'| cat >> srcinfo.h
	@echo 'character(len=40) HOST'| cat >> srcinfo.h
	@echo parameter \(VERSION=\'${VERSION}\'\) >> srcinfo.h
	@echo parameter \(BRANCH=\'${BRANCH}\'\) >> srcinfo.h
	@echo parameter \(AUTHOR=\'${AUTHOR}\'\) >> srcinfo.h
	@echo parameter \(EXECDATE=\'${DATE}\'\) >> srcinfo.h
	@echo parameter \(HOST =\'${HOST}\'\) >> srcinfo.h
	@cp srcinfo.h $(SRCDIR)
clean:
	@rm -f srcinfo.h
