

def fate_single_mut():

	import numpy as np
	import matplotlib.pyplot as plt
	import random

	N=1e5
	WT = []
	Mut = []
	T = []
	s = 0.1
	mu = 0
	b = 1.0 + s
	d = 1.0

	WT.append(N/2)
	Mut.append(1)
	T.append(0)

	for i in range(100000):

		#Ask if the mutant pop. went extint
		if(Mut[-1]==0):
			mu=1e-5
		else:
			mu=0

		sum_rate = b*Mut[-1] + d*Mut[-1] + mu*WT[-1]
		#print(sum_rate)

		#Produce random numbers
		r1 = random.random()
		r2 = random.random()

		#sample time of the next event
		tau = -(1/sum_rate)*np.log(1-r1)
		T.append(T[-1]+tau)

		#what happened?

		if(0 < r2 <= ((b*Mut[-1])/sum_rate) ):
			Mut.append(Mut[-1]+1)
			#print("birth event")
		elif( ((b*Mut[-1])/sum_rate) < r2 <= ((b*Mut[-1] + d*Mut[-1])/sum_rate)  ):
			Mut.append(Mut[-1]-1)
			#print("death event")
		else:
			Mut.append(1)
			#print("mutation event")

	plt.figure(figsize=(12,8))
	plt.plot(T,Mut)
	plt.plot(T, np.exp(s*np.array(T)))
	plt.hlines(1/s,0,max(T))
	plt.vlines((0.577216/s),0.1,1/s)
	plt.yscale('log')
	plt.xlabel('t', fontsize = 14)
	plt.ylabel(r'$\log{N_{\mathrm{mut}}}$', fontsize = 14)
	plt.title('Fate of a single mutant', fontsize = 15)
	plt.yticks([1/s], [r"$\frac{1}{s}$"], fontsize=16)
	plt.xticks([(0.577216/s)], [r"$<\tau_{\mathrm{est}}>$"], fontsize=14)
	plt.show()
