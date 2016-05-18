#Clase metodo diferencias finitas
# Para ejecutar: mpiexec -n 4 python test.py

from mpi4py import MPI
import time
import numpy
import itertools

class Diferencias_finitas:
	num_proc = MPI.COMM_WORLD.Get_rank()
	num_proc_total = MPI.COMM_WORLD.Get_size()
	def __init__(self,*args):
		if len(args) == 0:
			self.tolerancia = 10**(-6)
			self.iteraciones = 10000
			self.a = 0
			self.b = 1
			self.fa = 0
			self.fb = 1
			self.n = 10000
		else:
			self.tolerancia = args[0]
			self.iteraciones = args[1]
			self.a = args[2]
			self.b = args[3]
			self.fa = args[4]
			self.fb = args[5]
			self.n = args[6]
	def show(self):
		print "Ajustes:\n -Tolerancia = %f" %(self.tolerancia)
		print " -Iteraciones = %d" %(self.iteraciones)
		print " -a = %f" %(self.a)
		print " -b = %f" %(self.b)
		print " -fa = %f" %(self.fa)
		print " -fb = %f" %(self.fb)
		print " -n = %d" %(self.n)
		print " -numero de procesador = %d/%d" %(self.num_proc, self.num_proc_total)
	def f(self,x):
		return 12*x*x
	def metodo(self):
		# En primer lugar reparto los nodos entre los procesos
		num_puntos = self.n/self.num_proc_total
		if self.num_proc == 0:
			print "Cada nodo se ocupa de %d puntos." %(num_puntos)
		# Calculo que intervalo trabaja cada nodo
		punto_inicio = (self.num_proc)*num_puntos+1
		punto_final = (self.num_proc+1)*num_puntos
		print "El procesador %d trabaja en el intervalo (%d,%d)" \
			%(self.num_proc,punto_inicio,punto_final)
		time.sleep(0.1)
		# Inicializo los vectores que contienen los valores de las funciones f y u
		h = float(1.0)/(self.n-1)
		vf = []
		vu = []
		for i in range(punto_inicio, punto_final+1):
			x = self.a+(i-1)*h
			vu.append(self.fa+(i-1)*(float(self.fb-self.fa))*h)
			vf.append(self.f(x))
		# Inicializo un vector para contener las diferencias
		dif = []
		for i in range(0,num_puntos):
			dif.append(1)
		dif[0]=0
		dif[num_puntos-1]=0
		time.sleep(0.1)
		# Coloco las condiciones de frontera
		if(self.num_proc == 0):
			vu[0] = self.fa
		if(self.num_proc == self.num_proc_total-1):
			vu[num_puntos-1] = self.fb
		time.sleep(0.3)		
		# Entro en las iteraciones
		for i in range(1,self.iteraciones+1):
			# Si no soy el primer procesador, necesito pasar informacion al nodo de la izquierda
			if self.num_proc > 0:
				MPI.COMM_WORLD.send(vu[0],dest=self.num_proc-1)
			# Si no soy el ultimo procesador, necesito pasar informacion a la derecha
			if self.num_proc < self.num_proc_total-1:
				MPI.COMM_WORLD.send(vu[num_puntos-1],dest=self.num_proc+1)
			# Si no soy el primer procesador, recibo informacion de la izquierda	
			if self.num_proc > 0:
				aux_izq = MPI.COMM_WORLD.recv(source=self.num_proc-1)
			# Si no soy el ultimo procesador, recibo informacion de la derecha
			if self.num_proc < self.num_proc_total-1:
				aux_dcha = MPI.COMM_WORLD.recv(source=self.num_proc+1)
			# Compruebo si tengo que hacer Jacobi o he alcanzado la tolerancia
			if max(dif)>self.tolerancia:
				#Calculo una iteracion de Jacobi
				# Si soy el primer procesador ya tengo un extremo calculado.
				if self.num_proc == 0:	
					for j in range(1,num_puntos-1):
						vu_antiguo = vu[j]
						vu[j] = 0.5*(vu[j-1]+vu[j+1]+h*h*vf[j])
						dif[j] = abs(vu_antiguo-vu[j])
					# Ahora calculo el extremo derecho
					vu_antiguo = vu[num_puntos-1]
					vu[num_puntos-1] = 0.5*(vu[num_puntos-2]+aux_dcha+h*h*vf[num_puntos-1])
					dif[num_puntos-1] = abs(vu_antiguo-vu[num_puntos-1])
				# Si soy el ultimo procesador ya tengo un extremo calculado
				if self.num_proc == self.num_proc_total-1:
					for j in range(1,num_puntos-1):
						vu_antiguo = vu[j]
						vu[j] = 1/2.0*(vu[j-1]+vu[j+1]+h*h*vf[j])
						dif[j] = abs(vu_antiguo-vu[j])
					# Ahora calculo el extremo izquierdo
					vu_antiguo = vu[0]
					vu[0] = 1/2.0*(aux_izq+vu[1]+h*h*vf[0])
					dif[0] = abs(vu_antiguo-vu[0])
				# Si soy un procesador de enmedio
				if self.num_proc > 0 and  self.num_proc < self.num_proc_total-1:
					for j in range(1,num_puntos-1):
						vu_antiguo = vu[j]
						vu[j] = 1/2.0*(vu[j-1]+vu[j+1]+h*h*vf[j])
						dif[j] = abs(vu_antiguo-vu[j])
					# Ahora calculo el extremo derecho
					vu_antiguo = vu[num_puntos-1]
					vu[num_puntos-1] = 1/2.0*(vu[num_puntos-2]+aux_dcha+h*h*vf[num_puntos-1])
					dif[num_puntos-1] = abs(vu_antiguo-vu[num_puntos-1])
					# Ahora calculo el extremo izquierdo
					vu_antiguo = vu[0]
					vu[0] = 1/2.0*(aux_izq+vu[0]+h*h*vf[0])
					dif[0] = abs(vu_antiguo-vu[0])
		time.sleep(0.3)
		# Ya se ha aplicado Jacobi. Ahora unimos todos los vectores
		#Creo en el procesador 0 un vector para almacenar el vector u
		
		u = MPI.COMM_WORLD.gather(vu,root=0)
		if(self.num_proc == 0):
			return list(itertools.chain.from_iterable(u))
			
def Modulo_Cargado():
	print "El modulo se ha cargado correctamente"
