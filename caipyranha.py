# -*- coding: utf-8 -*-

#### Trabajo por hacer
#
# Percentil 99/95 para eliminar singletons [ver numpy array]
#	mayus = open(url2+"/alignment_mayus.fas" ,"r")
#	alignment = AlignIO.read(mayus, "fasta")
#	for a in alignment:
#		print a.seq
#
# Abrir alineamientos con clustalx
# Conversor MEGA
# rawinput de la url
# alinear processed and processed_contig
#

import sys, getopt, os

from commands import *
import commands

from Bio import AlignIO, SeqIO, pairwise2
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq


###### Lectura de argumentos y mensaje de ayuda

def usage():
	print """
NOMBRE
\t Caipyranha - Procesador y analizador de secuencias\n
MODO DE EMPLEO
\t python caipyranha.py [-e] [-m] [-a] [-f] [-b]\n
\t Los nombres de los ficheros y carpetas deben ir sin espacios.\n
\t Se debe respetar el orden marcado arriba excepto para la opción -b, pudiéndose saltar alguna
\t opción excepto -e.\n
\t Debe indicar la pareja de primers utilizada para amplificar el producto.\n
\t Debe indicar el intervalo de la etiqueta que quiere utilizar para hacer comparaciones.\n
OPCIONES
\t -e \t Extrae el inserto del plásmido si la puntuación alineamiento local de uno de los primers
\t\t con la secuencia es mayor del 90% a la máxima. Utiliza la función "pairwise2" de Biopython.
\t -m \t Si dos etiquetas son iguales en el intervalo seleccionado, se unen las secuencia extraídas
\t\t por las regiones solapantes. Utiliza "merger" de la suite EMBOSS.
\t -a \t Alineamiento de las secuencias extraídas y/o unidas mediante MAFFT.
\t -f \t Obtener árboles filogenéticos mediante los métodos de Neighbor-Joining y Máxima Verosimilitud
\t\t mediante PHYLIPnew
\t -b \t Realizar un BLAST a través de la web y devuelve el número de resultados indicado.
"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "abc:ehfm", ["ali", "blast", "coli", "extract", "help", "filo", "merger"])
except getopt.GetoptError, err:
	# print help information and exit:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

###### Invocar al terminal

def run_command(cmd):
	getstatusoutput(cmd)

##### Función merger

def merger():

	salida = open(url2 + "/processed_contig.fas","w") #crea archivo en la carpeta indicada

	record = list(SeqIO.parse(open(url2 + "/processed.fas"), format="fasta"))

	for a in range(0,len(record)):
		if a < len(record)-1:
			b = a+1
			seq1 = record[a]
			seq2 = record[b]

			if seq1.id[0:5] == seq2.id[0:5]:

				open_seq1 = open(url2 + "/seq1.fas","w")
				open_seq1.write("%s" % (">"+ seq1.id +"\n"+ seq1.seq))
				open_seq1.close()
	
				open_seq2 = open(url2 + "/seq2.fas","w")
				open_seq2.write("%s" % (">"+ seq2.id +"\n"+ seq2.seq))
				open_seq2.close()

	#			comando_merger = str("merger -asequence " + url2 + "/seq1.fas -bsequence " + url2 + "/seq2.fas -outfile " + url2 + "/contig.merger -outseq " + url2 + "/contig.fas")
				comando_merger = str("megamerger -asequence " + url2 + "/seq1.fas -bsequence " + url2 + "/seq2.fas -wordsize 20 -outfile " + url2 + "/contig.merger -outseq " + url2 + "/contig.fas")
				run_command(comando_merger)

				overlap = open(url2 + "/contig.fas", "r")
				salida.write("%s" % (str.upper(overlap.read()) + "\n\n"))

				overlap.close()
				os.remove(url2 + "/seq1.fas")
				os.remove(url2 + "/seq2.fas")
				os.remove(url2 + "/contig.merger")
				os.remove(url2 + "/contig.fas")

#			elif:
#				seq1.

	salida.close()

##### Función BLAST

def blastear():

	print "\n\nComenzando búsqueda en BLAST... por favor, espere"

	try:
		record = list(SeqIO.parse(open(url2 + "/processed_contig.fas"), format="fasta"))
		fasta_string = open(url2 + "/processed_contig.fas").read()
		result_handle = NCBIWWW.qblast("blastn", "nr", fasta_string)
	except:
		record = list(SeqIO.parse(open(url2 + "/processed.fas"), format="fasta"))
		fasta_string = open(url2 + "/processed.fas").read()
		result_handle = NCBIWWW.qblast("blastn", "nr", fasta_string)

	save_file = open(url2 + "/my_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

	result_handle = open(url2 + "/my_blast.xml")
	blast_records = NCBIXML.parse(result_handle) #NCBIXML.parse para múltiples secuencias

	blast_save = open(url2 + "/blast", "w")

	for blast_record in blast_records:
		for alignment in blast_record.alignments[0:1]:
			for hsp in alignment.hsps:
				blast_save.write("****Alignment****\n")
				blast_save.write("sequence: " + alignment.title + "\n")
				blast_save.write("length: " + str(alignment.length) + "\n")
				blast_save.write("e value: " + str(hsp.expect) + "\n")
				blast_save.write(hsp.query[0:75] + "...\n")
				blast_save.write(hsp.match[0:75] + "...\n")
				blast_save.write(hsp.sbjct[0:75] + "...\n")
				blast_save.write("\n")

	print "BLAST realizado con éxito"

#	os.remove(url2 + "/my_blast.xml")

###### Función alinear

def alinear():

	try:
		open(url2 + "/processed_contig.fas", "r")
		aliin = "processed_contig.fas"
		aliout = "alignment_contig.fas"

	except:
		open(url2 + "/processed.fas", "r")
		aliin = "processed.fas"
		aliout = "alignment.fas"

	comando_mafft = str("mafft " + url2 + "/" + aliin + " > " + url2 + "/" + aliout)
	print "\n\nRealizando alineamiento con MAFFT para " + aliin
	run_command(comando_mafft)

	mayus_open = open(url2 + "/" + aliout, "r")
	mayus_leer = mayus_open.read()

	str_mayus = str.upper(mayus_leer)

	f_salida = url2 + "/" + aliout
	mayus = open(f_salida ,"w")
	mayus.write(str_mayus)
	mayus_open.close()
	mayus.close()

#	aliout = list(aliout)
#	aliout = aliout[:len(aliout)-4]
#	aliout = ''.join(str(n) for n in aliout)

#	mega = open(url2 + "/" + aliout + ".meg", "w")
#	mega_open = open(f_salida, "r")
#	mega_leer = mega_open.read()
#	str_mega = str(mega_leer).replace(">", "#")
#	mega.write("#Mega\n")
#	mega.write("Title: " + aliin + "\n\n")
#	mega.write(str_mega)

#	mega_open.close()
#	mega.close()

##### Función filogenia

def filogenia():
	try:
		mayus = open(url2 + "/alignment_contig.fas" ,"r")
		alignment = AlignIO.read(mayus, "fasta")
		phylip = url2 + "/alignment_contig.phy"
		phy = open(phylip ,"w")
	except:
		mayus = open(url2 + "/alignment.fas" ,"r")
		alignment = AlignIO.read(mayus, "fasta")
		phylip = url2 + "/alignment.phy"
		phy = open(phylip ,"w")

#	nexus = url2 + "/prueba.nex"
#	nexus = open(phylip ,"w")

#	for i, record in enumerate(alignment):
#		record.id = "seq%i" % i

	AlignIO.write([alignment], phy, "phylip") # convertir a formato phylip# Blastear cada secuencia y dar la salida	
	#AlignIO.write([alignment], nexus, "nexus")

#	a[0:0]= [" ", " "]
#	a[5:5]= [" ", " "]

	phy.close()
	mayus.close()
#	nexus.close()

	comando_fdnadist = "fdnadist -sequence " + url2 + "/alignment_contig.phy -method f -outfile " + url2 + "/alignment.fdnadist" # matriz distancias
	comando_fneighbor = "fneighbor -datafile " + url2 + "/alignment.fdnadist -outfile " + url2 + "/alignment.fneighbor -outtreefile " + url2 +  "/neighbor.treefile" # Neighbor-joinig
	comando_fdrawtreeN = "fdrawtree -intreefile " + url2 + "/neighbor.treefile -plotfile " + url2 + "/neighbour.ps" # árbol NJ

	comando_fdnaml = "fdnaml -sequence " + url2 + "/alignment_contig.phy -outfile " + url2 + "/alignment.fdnaml -intreefile " + url2 + "/neighbor.treefile -outtreefile " + url2 + "/max_like.treefile" # maximum likelihood
	comando_fdrawtreeM = "fdrawtree -intreefile " + url2 + "/max_like.treefile -plotfile " + url2 + "/max_like.ps" # árbol ML

#	comando_fdnapars = "fdnapars -sequence " + url2 + "/alignment_contig.phy -outfile " + url2 + "/alignment.fdnapars -outtreefile " + url2 + "/parsimony.treefile" # pasimony
#	comando_fdrawtreeP = "fdrawtree -intreefile " + url2 + "/parsimony.treefile -plotfile " + url2 + "/parsimony.ps"

	run_command(comando_fdnadist)
	run_command(comando_fneighbor)
	run_command(comando_fdrawtreeN)
	run_command(comando_fdnaml)
	run_command(comando_fdrawtreeM)
#	run_command(comando_fdnapars)
#	run_command(comando_fdrawtreeP)

#####################################################
######creamos archivo con listado de imagenes########
#####################################################

# Indicar la ruta en que se encuentran la secuencia a procesar

#url = "/home/fruano/science/programs/pylab/caipyranha/1"
url = "/home/fruano/pylab/caipyranha/2"

home = str(commands.getoutput("echo $HOME"))

url_temp_folders = str(home + "/.temp_folders")		#poner oculto

url_temp_files = str(home + "/.temp_files")				#poner oculto

comando1 = str("ls " + url + " > " + home + "/.temp_folders")		#oculto contiene lista de folders

run_command(comando1)

################################################################
######creamos lista que contiene direcciones de carpetas########
################################################################

lista_folders = []

file_read_folders = open(url_temp_folders)

for linea in file_read_folders:
				
		lista_folders.append(linea)


###Con este bucle eliminamos el caracter de retorno de carro \n 

i = 0

for folder in lista_folders:

		ruta = folder[0:-1]
		lista_folders[i] = str(ruta)
		i += 1

################################################################
######definimos las extensiones con que vamos a trabajar #######
################################################################

extension = ".txt"

#####################################################
######Aplicamos el algoritmo_paco a cada una ########
#####################################################

for o, a in opts:
	if o in ("-h", "--help"):
		usage()
		sys.exit()

print "\nCaipyranha"
print "\nThis is the Funny Program for extracting inserts of vectors... and more\n"

##Definimos el algoritmo que se aplica a cada fichero

primer_forward = "GCTTTTGTACACACCGCCCGTCGC"
primer_forward_revcom = Seq(primer_forward)
primer_forward_revcom = primer_forward_revcom.reverse_complement()
primer_forward_revcom = "".join(primer_forward_revcom)

primer_reverse = "ATATGCTTAAATTCAGCGGG"
primer_reverse_revcom = Seq(primer_reverse)
primer_reverse_revcom = primer_reverse_revcom.reverse_complement()
primer_reverse_revcom = "".join(primer_reverse_revcom)


###### Función que separa el inserto del vector

def algoritmo(sequencia):

	print sequencia # para conocer el progreso

	items = []
	seq = ''
	name = ''
	index = 0

	direccion = str(url2 +"/" + sequencia)

	secu = open(direccion , "r").readlines()

	for line in secu:
		if line.startswith(">"):
			if index >= 1:
					items.append(seq)
					seq = ''
			esp = line.find(".")
			index += 1
			name = ">" + line[15:esp] + "\n"
		else:
			seq += line[:-1]

	localforward = pairwise2.align.localms(seq, primer_forward, 2, -1, -1, -.1)
	bestforward = localforward[0]
	localreverse = pairwise2.align.localms(seq, primer_reverse, 2, -1, -1, -.1)
	bestreverse = localreverse[0]

	if bestforward[2] > len(primer_forward)*2*0.95:
#		print bestforward[3]
		b = 650
		a = bestforward[3]
		its_seq = seq[a:b]
	elif bestreverse[2] > len(primer_reverse)*2*0.95:
#		print bestreverse[3]
		b = 650
		a = bestreverse[3]
		its_seq = seq[a:b]
		its_seq = Seq(its_seq)
		its_seq = its_seq.reverse_complement()
	else:
		its_seq = "N"
#				error += name
		print ("Error en %s") % name
		errores.write ("%s" % (url2 + "/" + sequencia + "\n"))

	bestforward[3] == -1
	bestreverse[3] == -1

#	if its_seq != "N":
	salida.write("%s" % (name + its_seq + "\n\n"))



#############################################################
######Creamos la lista de archivos y aplicamos algoritmo#####
#############################################################



for folder in lista_folders:

	url2 = url + "/" + folder

	print "Leyendo carpeta \"%s\" \n" % folder

	comando2 = str("ls " + url2 + " | grep '" + extension + "'" + " > " + home + "/.temp_files")	#oculto contiene lista de files

	run_command(comando2)

	lista_files = []

	file_read_files = open(url_temp_files)

#creamos la lista_files que contiene direcciones de imagenes

	for linea in file_read_files:

		linea = str(linea[0:-1])		#eliminamos el caracter de retorno de carro \n		
		lista_files.append(linea)		

	###Comprobacion de carpeta. Si esta vacia se sale del programa

	if len(lista_files) == 0:
		break

###Creamos archivo que contiene los calculos 

	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit()
		elif o in ("-e", "--extract"):
			salida = open(url2 + "/processed.fas","w") #crea archivo en la carpeta indicada
			errores = open(url2 + "/errores", "w")
			errores.write ("No se encontraron los primers seleccionados en estos ficheros:\n\n")

			for sequencia in lista_files:
				algoritmo(str(sequencia))
			salida.close()
			errores.close()
			print "\n\nInsertos extraídos del plásmido"
		elif o in ("-m", "--merger"):
			merger()
			print "\n\nSecuencias solapadas"
		elif o in ("-b", "--blast"):
			blastear()
		elif o in ("-a", "--ali"):
			alinear()
		elif o in ("-f", "--filo"):
			filogenia()
		else:
			assert False, "error"

	print "\n\nTodo correcto, a la siguiente carpeta..."
	print "@"*60
	print "\n"

#####################################################
######Borramos archivos temporales y ceramos ########
#####################################################

file_read_folders.close()
file_read_files.close()

comando_borrar_folders = "rm " + url_temp_folders
comando_borrar_files = "rm " + url_temp_files  

run_command(comando_borrar_folders)
run_command(comando_borrar_files)

print "\nCaipyranha ha finalizado la tarea. Gracias por utilizarlo. (o'.'o)\n"

