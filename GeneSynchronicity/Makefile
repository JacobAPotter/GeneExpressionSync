# Makefile for running a program that uses Computer vision Techniques to
# score the syncronicity of gene expression data
# Jacob Potter
# 3/20/2019


CC	= g++
CFLAGS 	= -g
TARGET	= SyncScore
OBJS	= SyncScore.o

$(TARGET):	$(OBJS)
		$(CC) -o $(TARGET) $(OBJS)

readExpr.o:	SyncScore.cpp readExpr.cpp readExpr.h Operations.cpp Operations.h timer.h 
		$(CC) $(CFLAGS) -c SyncScore.cpp

clean:
		/bin/rm -f *.o $(TARGET)

run1:
	./SyncScore ExpressionData/M13_24genes.txt 24 67
run2:
	./SyncScore ExpressionData/test5genes10indiv.txt 5 10
run3:
	./SyncScore RandomData/MyersRandom_364indiv_100genes.txt 100 364
