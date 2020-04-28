import sys
import numpy as np
import matplotlib.pyplot as plt
import random

def fate_single_mut(s, n_pop, short_t = False, long_t = True):

	m = 500
	if short_t == True:
		factor = 0.05
	else:
		factor = 4

	T_total = (1/s)*factor
	T_avg = np.linspace(0, T_total, m)
	Mut_avg = np.zeros(shape = (m))
	n_pop = n_pop
	N=1e5
	fig, ax = plt.subplots(figsize=(12,8))
	fig.suptitle(r'Fate of a single mutant ; $\left< n|\mathrm{not\;extinct}\right>$', fontsize = 15)
	for n in range(n_pop): 

		#WT_temp = np.zeros(shape = (1, 1000))
		WT = np.array([])
		Mut = np.array([])
		T = np.array([])
		s = s
		mu = 0
		b = 1.0 + s
		d = 1.0

		WT = np.append(WT,N/2)
		Mut = np.append(Mut,1)
		T = np.append(T,0)

		i = 0
		while T[i] < T_total*1.1:

			#Ask if the mutant pop. went extint
			if(Mut[-1]==0):
				mu=1e-5

				WT = np.array([])
				Mut = np.array([])
				T = np.array([])

				WT = np.append(WT,N/2)
				Mut = np.append(Mut,1)
				T = np.append(T,0)

				i = 0
			else:
				mu=0

			sum_rate = b*Mut[-1] + d*Mut[-1] + mu*WT[-1]
			#print(sum_rate)

			#Produce random numbers
			r1 = random.random()
			r2 = random.random()

			#sample time of the next event
			tau = -(1/sum_rate)*np.log(1-r1)
			T = np.append(T, T[-1]+tau)

			#what happened?

			if(0 < r2 <= ((b*Mut[-1])/sum_rate) ):
				Mut = np.append(Mut, Mut[-1]+1)
				#print("birth event")
			elif( ((b*Mut[-1])/sum_rate) < r2 <= ((b*Mut[-1] + d*Mut[-1])/sum_rate)  ):
				Mut = np.append(Mut, Mut[-1]-1)
				#print("death event")
			else:
				Mut = np.append(Mut, 1)
				#print("mutation event")

			i+=1

		k = 0
		j = 0
		while (T[j]<=T_total):
			while (T[j]>=T_avg[k]):
				Mut_avg[k] += Mut[j]
				k+=1
				if k >= m:
					break
			j+=1

		ax.plot(T,Mut, 'g', alpha = 0.05)
	
	if short_t == True:
		ax.plot(T_avg[:-1], Mut_avg[:-1]/n_pop, 'g-', ms = 30, label = r'$\left<n\right>$')
		ax.plot(T_avg, 1+T_avg, 'b--', label = r'$1+t$')
		ax.plot(T_avg, np.exp(s*T_avg) + (np.exp(s*T_avg)-1)/(s), 'b-', label = r'$e^{st} + \frac{e^{st}-1}{s}$')
		ax.hlines(1/s,0,max(T))
		ax.vlines((0.577216/s),0.5,1/s)
		#ax.set_yscale('log')
		ax.set_xlim(0,(1/s)*0.02)
		ax.set_ylim(0.8, (1/s)*0.04)
		ax.set_xlabel('t', fontsize = 14)
		ax.set_ylabel(r'$N_{\mathrm{mut}}$', fontsize = 14)
		ax.set_title(r'$t\ll\frac{1}{s}$', fontsize=16)
		ax.set_yticks([1/s])
		ax.set_yticklabels([r"$\frac{1}{s}$"])
		ax.set_xticks([(0.577216/s)])
		ax.set_xticklabels([r"$<\tau_{\mathrm{est}}>$"])
		ax.tick_params(labelsize = 16)
		ax.legend(fontsize = 16)
	else:
		ax.plot(T_avg[:-1], Mut_avg[:-1]/n_pop, 'g-', label = r'$\left<n\right>$')
		ax.plot(T_avg, np.exp(s*T_avg)/s, 'b--', label = r'$\frac{1}{s}e^{st}$')
		ax.plot(T_avg, np.exp(s*T_avg) + (np.exp(s*T_avg)-1)/(s), 'b-', label = r'$e^{st} + \frac{e^{st}-1}{s}$')
		ax.hlines(1/s,0,max(T))
		ax.vlines((0.577216/s),0.5,1/s)
		ax.set_yscale('log')
		ax.set_xlabel('t', fontsize = 14)
		ax.set_ylabel(r'$\log{N_{\mathrm{mut}}}$', fontsize = 14)
		ax.set_title(r'$t\gg\frac{1}{s}$', fontsize=16)
		ax.set_yticks([1/s])
		ax.set_yticklabels([r"$\frac{1}{s}$"])
		ax.set_xticks([(0.577216/s)])
		ax.set_xticklabels([r"$<\tau_{\mathrm{est}}>$"])
		ax.tick_params(labelsize = 16)
		ax.legend(fontsize = 16)

	#plt.show()
