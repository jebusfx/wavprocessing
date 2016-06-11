from copy import deepcopy
from matplotlib import pyplot
import numpy as np

class Matrix:
	def __init__(self,mat,zeropos):
		self.mat = mat
		self.zeropos = zeropos		
		self.length = len(mat)

	def getPeriod(self):
		for i in xrange(1,self.length):
			if self.mat[i] == self.mat[0]:
				self.period = i
				break
		self.mat = self.mat[:self.period]
		self.length = len(self.mat)

	def pad(self,pad,length=0):
		for i in xrange(abs(length-pad)):
			self.mat.append(0)

	@staticmethod
	def ceil(x,y):
		return (x+y-1)/y

	@staticmethod
	def shiftVector(vector,i):
		return vector[len(vector)-i:] + vector[:len(vector)-i]

class Convolution:
	CIRCULAR, FINITE, PERIODIC = xrange(3)	
	def __init__(self,matrixA,matrixB,convtype):
		self.matrixA = matrixA
		self.matrixB = matrixB
		self.matrixAOriginal = deepcopy(matrixA)
		self.matrixBOriginal = deepcopy(matrixB)		
		self.isPeriodic = (True if convtype==Convolution.PERIODIC else False)
		if convtype == Convolution.CIRCULAR:
			self.padding = self.matrixA.length - self.matrixB.length
			self.padFunction = self.padMatrixCircular
			self.greaterLength = (matrixA.length if (matrixA.length >= matrixB.length) else matrixB.length)			
		elif convtype == Convolution.FINITE:
			self.padding = self.matrixA.length + self.matrixB.length - 1
			self.padFunction = self.padMatrixFinite
			self.greaterLength = self.padding	
		elif convtype == Convolution.PERIODIC:
			self.periodicResult = []
			self.matrixA.getPeriod()
			self.padding = self.matrixA.length + self.matrixB.length - 1
			self.padFunction = self.padMatrixFinite
			self.greaterLength = self.padding			

	def calculatePeriodicRes(self):
		newpad = self.matrixA.length * Matrix.ceil(self.matrixRes.length,self.matrixA.length)
		auxMatrixRes = deepcopy(self.matrixRes)
		auxMatrixRes.pad(auxMatrixRes.length,newpad)
		auxMatrixRes.length = len(auxMatrixRes.mat)
		self.periodicResult = [0 for i in xrange(self.matrixA.period)]
		for i in xrange(auxMatrixRes.length):
			modulo = i%self.matrixA.period
			self.periodicResult[modulo] += auxMatrixRes.mat[i]	
		self.perResult = Matrix(self.periodicResult[:self.matrixA.period],0)	
		print self.periodicResult[:self.matrixA.period]
		self._plotPer()

	def padMatrixCircular(self):	
		if self.padding > 0:			
			self.matrixB.pad(self.padding)
		elif self.padding < 0:
			self.matrixA.pad(abs(self.padding))

	def padMatrixFinite(self):
		self.matrixA.pad(self.padding,self.matrixA.length)
		self.matrixB.pad(self.padding,self.matrixB.length)	

	def concatenateMatrix(self):				
		temp = []
		for i in xrange(self.greaterLength):								
			aux = (Matrix.shiftVector(self.matrixA.mat,i) if i else self.matrixA.mat)
			temp.append(aux)
		self.matrixA.mat = temp					

	def convolute(self):		
		self.padFunction()
		self.concatenateMatrix()		
		self.matrixA.mat = np.matrix(self.matrixA.mat).transpose()
		self.matrixB.mat = np.matrix(self.matrixB.mat).transpose()		
		self.matrixRes = Matrix(np.dot(self.matrixA.mat,self.matrixB.mat).transpose().tolist()[0],3)
		if self.isPeriodic:
			self.calculatePeriodicRes()
		return np.matrix(self.matrixRes.mat).transpose()

	def _plotPeriodic(self):
		xPer = np.arange(self.perResult.zeropos*-1,self.perResult.length-self.perResult.zeropos,1.0)
		pyplot.plot([xPer,xPer],[[0 for i in xrange(self.perResult.length)],self.perResult.mat],'k-',lw=2)	
		pyplot.grid(True)
		pyplot.axis([xPer[len(xPer)-1]*-3,xPer[len(xPer)-1]*3,max(self.perResult.mat)*-2,max(self.perResult.mat)*2])			
		pyplot.show()

	def _plot(self):
		xA = np.arange(self.matrixAOriginal.zeropos*-1,self.matrixAOriginal.length-self.matrixAOriginal.zeropos,1.0)
		xB = np.arange(self.matrixBOriginal.zeropos*-1,self.matrixBOriginal.length-self.matrixBOriginal.zeropos,1.0)
		xRes = np.arange(self.matrixRes.zeropos*-1,self.matrixRes.length-self.matrixRes.zeropos,1.0)
		pyplot.figure(1)	
		plotX = pyplot.subplot(311)
		plotX.plot([xA,xA],[[0 for i in xrange(self.matrixAOriginal.length)],self.matrixAOriginal.mat],'k-',lw=2)	
		plotX.set_title('X(n)')
		plotX.grid(True)
		plotX.axis([xA[len(xA)-1]*-3,xA[len(xA)-1]*3,max(self.matrixAOriginal.mat)*-2,max(self.matrixAOriginal.mat)*2])		
		plotY = pyplot.subplot(312)
		plotY.plot([xB,xB],[[0 for i in xrange(self.matrixBOriginal.length)],self.matrixBOriginal.mat],'k-',lw=2)	
		plotY.set_title('Y(n)')
		plotY.grid(True)
		plotY.axis([xB[len(xB)-1]*-3,xB[len(xB)-1]*3,max(self.matrixAOriginal.mat)*-2,max(self.matrixAOriginal.mat)*2])		
		plotRes = pyplot.subplot(313)
		plotRes.plot([xRes,xRes],[[0 for i in xrange(self.matrixRes.length)],self.matrixRes.mat],'k-',lw=2)	
		plotRes.set_title('Resultado')
		plotRes.grid(True)
		plotRes.axis([xRes[len(xRes)-1]*-3,xRes[len(xRes)-1]*3,max(self.matrixRes.mat)*-2,max(self.matrixRes.mat)*2])
		pyplot.show()

if __name__ == "__main__":	
	while(True):
		opc = input('1. Convolucion Circular\n2. Convolucion Finita\n3. Convolucion Periodica\n4. Salir\n')
		if opc == 1:
			x = map(float, raw_input('Teclea una secuencia para x: ').split())
			originx = input('Origen: ')
			y = map(float, raw_input('Teclea una secuencia para y: ').split())
			originy = input('Origen: ')
			conv = Convolution(Matrix(x,originx),Matrix(y,originy),Convolution.CIRCULAR)
		elif opc == 2:
			x = map(float, raw_input('Teclea una secuencia para x: ').split())
			originx = input('Origen: ')
			y = map(float, raw_input('Teclea una secuencia para y: ').split())
			originy = input('Origen: ')
			conv = Convolution(Matrix(y,originx),Matrix(x,originy),Convolution.FINITE)
		elif opc == 3:
			x = map(float, raw_input('Teclea una secuencia periodica para x: ').split())
			originx = input('Origen: ')
			y = map(float, raw_input('Teclea una secuencia para y: ').split())
			originy = input('Origen: ')
			conv = Convolution(Matrix(x,originx),Matrix(y,originy),Convolution.PERIODIC)
		else:
			exit(0)
		print conv.convolute()
		conv._plot()
	