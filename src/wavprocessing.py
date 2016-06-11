from datetime import datetime
from copy import deepcopy
from convolution import Convolution, Matrix
from cmath import exp,pi
from matplotlib import pyplot
from os import path
from tempfile import mkdtemp
from Tkinter import *
from tkFileDialog import askopenfilename
from sys import getsizeof
from scipy.io import wavfile
from uuid import uuid1
import numpy as np

class WavProcessing:
	OPTIONS = {}
	OPTIONS['defaultextension'] = '.wav'
	OPTIONS['filetypes'] = [('WAVE files','.wav')]

	def __init__(self):
		global OPTIONS
		root = Tk()
		root.withdraw()
		self.filename = askopenfilename(**self.OPTIONS)         		
		self.sampleFrequency, self.audioData = wavfile.read(self.filename)		
		self.audioDataOriginal = deepcopy(self.audioData[:])

	def amplify(self,factor):
		self.audioData *= factor

	def diminish(self,factor):
		self.audioData /= factor

	def moveleft(self,factor):		
		zeros = np.zeros((factor,2),dtype=np.int16)		
		self.audioData = np.concatenate((self.audioData[factor:],zeros[:]))		

	def moveright(self,factor):
		zeros = np.zeros((factor,2),dtype=np.int16)						
		self.audioData = np.concatenate((zeros[:],self.audioData[:self.audioData.shape[0]-factor]))

	def reflex(self):
		self.audioData *= -1

	def diezmado(self,factor):				
		limit = self.audioData.shape[0]		
		newAudioData = np.zeros((limit,2),dtype=np.int16)		
		for i in xrange(limit):
			diezm = factor*i
			if diezm > -1 and diezm < limit:
				newAudioData[i] = self.audioData[diezm]
		self.audioData = newAudioData[:]

	def interpoladocero(self):
		temp = np.zeros((self.audioData.shape[0]*2,self.audioData.shape[1]),dtype=np.int16)
		for i in xrange(self.audioData.shape[0]):
			temp[i] = self.audioData[i]
			temp[i+1] = 0
		self.audioData = deepcopy(temp)

	def interpoladoescalon(self):
		temp = np.zeros((self.audioData.shape[0]*2,self.audioData.shape[1]),dtype=np.int16)
		for i in xrange(self.audioData.shape[0]):
			temp[i] = self.audioData[i]
			temp[i+1] = self.audioData[i]
		self.audioData = deepcopy(temp)

 	def interpoladolineal(self,factor):
		temp = np.zeros((self.audioData.shape[0]*factor,self.audioData.shape[1]),dtype=np.int16)		
		for i in xrange(self.audioData.shape[0]-1):
			temp[i] = self.audioData[i]
			prom = abs(self.audioData[i+1] - self.audioData[i])/factor			
			for j in xrange(1,factor+1):
				temp[i + j] = self.audioData[i] + prom
		self.audioData = deepcopy(temp)

	def convoluteWavs(self,wav):
		BLOCKSIZE = 1024				
		def preprocess():
			if self.audioData.shape[0] < wav.audioData.shape[0]:
				audioSize = self.audioData.shape[0]				
				if audioSize%BLOCKSIZE != 0:
					lastblock = self.audioData[(self.audioData.shape[0]-(self.audioData.shape[0]%BLOCKSIZE)):]
					temp = self.audioData[:(self.audioData.shape[0]-(self.audioData.shape[0]%BLOCKSIZE))]
					self.audioData = np.concatenate((temp,WavProcessing.pad(len(lastblock),lastblock,BLOCKSIZE)))
					resultConv = np.zeros((self.audioData.shape[0],2),dtype=np.int16)
			elif self.audioData.shape[0] >= wav.audioData.shape[0]:
				audioSize = wav.audioData.shape[0]
				if audioSize%BLOCKSIZE != 0:
					lastblock = wav.audioData[(wav.audioData.shape[0]-(wav.audioData.shape[0]%BLOCKSIZE)):]
					temp = wav.audioData[:(wav.audioData.shape[0]-(wav.audioData.shape[0]%BLOCKSIZE))]
					wav.audioData = np.concatenate((temp,WavProcessing.pad(len(lastblock),lastblock,BLOCKSIZE)))					
					resultConv = np.zeros((wav.audioData.shape[0],2),dtype=np.int16)
			return resultConv,Matrix.ceil(audioSize,BLOCKSIZE)
		
		resultConv, nblocks = preprocess()			
		for i in xrange(nblocks):
			start = BLOCKSIZE*i
			end = BLOCKSIZE*(i+1)
			conv = Convolution(Matrix(self.audioData[start:end,0].tolist(),0),Matrix(wav.audioData[start:end,0].tolist(),0),Convolution.CIRCULAR)
			resultConv[start:end] = conv.convolute()
		self.audioData = resultConv[:]
		conv._plot()
                               	
	@staticmethod
	def exponential(k,n):
			return exp(k*(-2j*pi/n))

	@staticmethod
	def dft(x):			
		x = np.asarray(x,dtype=complex)
		n = x.shape[0]		
		wMatrix = np.ones((n,n),dtype=complex)		
		precomputeW = [WavProcessing.exponential(i,n) for i in xrange(n)]
		for i in xrange(1,n):
			for j in xrange(1,n):
				k = (i*j)%n
				wMatrix[i][j] = precomputeW[k]		
		return np.dot(wMatrix,x)

	@staticmethod
	def fft(x):				
		x = np.asarray(x,dtype=complex)
		audioSize = x.shape[0]
		if audioSize <= 64:
			return WavProcessing.dft(x)
		else:
			even, odd = WavProcessing.evenoddview(x)
			even = WavProcessing.fft(even)
			odd = WavProcessing.fft(odd)
			ponder = [WavProcessing.exponential(i,audioSize) for i in xrange(audioSize)]			
			return np.concatenate([even + ponder[:audioSize/2] * odd, even + ponder[audioSize/2:] * odd])

	def preprocess_fft(self):		
		if self.audioData.shape[0]%(2**(self.audioData.shape[0].bit_length())) != 0: 
			nextPowerOfTwo = self.audioData.shape[0].bit_length()+1
			self.audioData = WavProcessing.pad(self.audioData.shape[0],self.audioData,2**nextPowerOfTwo)
		numberofbytes = (self.audioData.shape[0]-1).bit_length()
		temp = np.zeros((self.audioData.shape),dtype=complex)
		for i in xrange(self.audioData.shape[0]):
			temp[i] = self.audioData[WavProcessing.reversebits(i,numberofbytes)]
		self.audioData = temp[:]

	def process_fft(self):
		BLOCKSIZE = 2**16
		w.preprocess_fft()
		nblocks = self.audioData.shape[0]/BLOCKSIZE
		for i in xrange(nblocks):
			start = BLOCKSIZE*i
			end = BLOCKSIZE*(i+1)
			w.audioData[start:end,0] = WavProcessing.fft(w.audioData[start:end,0])
			w.audioData[start:end,1] = WavProcessing.fft(w.audioData[start:end,1])

	def plot_fft(self):		
		auxplot = np.fft.fftfreq(self.audioData.shape[0])
		pyplot.plot(auxplot, self.audioData.real,auxplot, self.audioData.imag)
		pyplot.show()
		
	def plotAudio(self):
		time = np.arange(0, self.audioData.shape[0], 1)
		time = time / self.sampleFrequency
		timeOriginal = np.arange(0, self.audioDataOriginal.shape[0], 1)
		timeOriginal = timeOriginal / self.sampleFrequency
		
		pyplot.figure(1)	
		plotF = pyplot.subplot(211)
		plotF.plot(time,self.audioData[:,0],color='k')
		plotF.set_ylabel('Amplitud')
		plotF.set_xlabel('Tiempo')		
		plotOriginal = pyplot.subplot(212)
		plotOriginal.plot(timeOriginal,self.audioDataOriginal[:,0],color='k')		
		plotOriginal.set_ylabel('Amplitud')
		plotOriginal.set_xlabel('Tiempo')		
		pyplot.show()

	def writeWav(self):		
		wavfile.write('{}{}.wav'.format(path.split(self.filename)[-1].split('.')[0],str(uuid1())),self.sampleFrequency,self.audioData)

	@staticmethod
	def evenoddview(x):
		numberofbytes = (x.shape[0]-1).bit_length()
		temp = np.zeros((x.shape),dtype=complex)
		evens,odds = [],[]
		for i in xrange(x.shape[0]):	
			no = WavProcessing.reversebits(i,numberofbytes)
			if no%2 == 0:
				evens.append(no)
			elif no%2 != 0:
				odds.append(no)
		evens.sort()
		odds.sort()
		midpoint = x.shape[0]/2
		for i in xrange(len(evens)):
			temp[i] = x[evens[i]]
		for i in xrange(len(odds)):
			temp[midpoint + i] = x[odds[i]]

		return temp[:x.shape[0]/2], temp[x.shape[0]/2:]
		#return x[::2],x[1::2]

	@staticmethod
	def reversebits(x,numberofbytes):						
		r = bin(x)[2:].zfill(numberofbytes)	
		return int(r[::-1], 2)

	@staticmethod
	def pad(shape,vector,blocksize):					
		x = Matrix.ceil(shape,blocksize)			
		topad = abs(shape - blocksize*x)
		temp = np.zeros((shape+topad,2),dtype=np.int16)			
		temp[:shape] = vector[:]
		return temp

def test_fft():
	x = np.random.random(2**16)
	print np.allclose(WavProcessing.fft(x),np.fft.fft(x))

if __name__ == "__main__":	
	w = WavProcessing()
	opc = input('1. Amplificar\n2. Atenuar\n3. Desplazar a la izquierda\n4. Desplazar a la derecha\n5. Reflejar\n6. Diezmado\n7. Interpolado cero\n8. Interpolado escalon\n9. Interpolado lineal\n10. Convolucion\n11. FFT\n')
	if opc == 1:
		factor = input('Teclea el factor de amplificacion: ')
		w.amplify(factor)
	elif opc == 2:
		factor = input('Teclea el factor de atenuacion: ')
		w.diminish(factor)
	elif opc == 3:
		factor = input('Teclea el factor de desplazamiento: ')
		w.moveleft(factor)
	elif opc == 4:
		factor = input('Teclea el factor de desplazamiento: ')
		w.moveright(factor)
	elif opc == 5:
		w.reflex()
	elif opc == 6:
		factor = input('Teclea el factor de diezmado: ')
		w.diezmado(factor)
	elif opc == 7:		
		w.interpoladocero()
	elif opc == 8:		
		w.interpoladoescalon()
	elif opc == 9:
		factor = input('Teclea el factor de interpolado: ')
		w.interpoladolineal(factor)
	elif opc == 10:
		newW = WavProcessing()
		w.convoluteWavs(newW)
	elif opc == 11:
		w.process_fft()
		w.plot_fft()
		exit(0)
	else:
		exit(0)	
	w.writeWav()
	w.plotAudio()
	exit(0)
