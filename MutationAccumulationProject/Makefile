# readindata.h
# Authors: Trevor Graham, Ruchira S. Datta
# Copyright (c) 2012, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# o Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# o Redistributions in binary form mus reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# o Neither the name of the University of California, San Francisco nor the
# names of its contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PUPROSE ARE
# DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMTIED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# macros
CC = mpicc
#CFLAGS = -Wall -pipe
CFLAGS = -O3 -Wall -pipe
#CFLAGS = -g -Wall -pipe
INCLUDE = -I/netapp/home/ruchira/include -I/netapp/home/ruchira/include/marray -I/netapp/home/ruchira/sprng2.0/include
#LFLAGS = -L/netapp/home/ruchira/sprng2.0/lib -L/netapp/home/ruchira/debug/lib 
LFLAGS = -L/netapp/home/ruchira/sprng2.0/lib -L/netapp/home/ruchira/lib 
LIBS = -lsprng -lgmp -lgsl -lgslcblas -lm -lmarray
#CFLAGS = -g

objects = bw9.o readindata.o

exe = /netapp/home/ruchira/bin/bw9

#${exe} :  ${objects}
#	export LD_RUN_PATH="/netapp/home/ruchira/debug/lib:$LD_RUN_PATH"; ${CC} ${CFLAGS} ${INCLUDE} ${LFLAGS} -o ${exe} ${objects} ${LIBS}
${exe} :  ${objects}
	${CC} ${CFLAGS} ${INCLUDE} ${LFLAGS} -o ${exe} ${objects} ${LIBS}
readindata.o : readindata.c
	${CC}  ${CFLAGS} ${INCLUDE} -c readindata.c
bw9.o : bw9.c
	${CC}  ${CFLAGS} ${INCLUDE} -c bw9.c

.PHONY : clean
clean :
	-rm $(exe) $(objects)
