#***MAKE FILE***#

#VARIABILI
inputMain=main_iplib.c
outputMain=main

#Se l'immagine è nello stesso percorso del file eseguibile, 
#basta cancellare il valore da assegnare alla variabile (es: dirIn= )
dirIn=img/original/
dirOut=img/results/
img1 = $(dirIn)flower.bmp

#Immagini per il bland
img2 = $(dirIn)mongolfiere.bmp
img3 = $(dirIn)flower2.bmp

concat_images = 1  #Flag per concatenare le immagini

#Decommentare se si vuole eseguire un check della memoria con Valgrind
MemCheck = #valgrind -v --leak-check=full



#TEST ALL-IN-ONE
runMainTest: compile corrupt gray brighten blend sharp edge avg emboss gauss



#COMPILAZIONE MAIN
compile: bmp.o ip_lib.o bmp.h ip_lib.h
	gcc $(inputMain) bmp.o ip_lib.o -o $(outputMain) -Wall -lm -g -O1 --ansi --pedantic -g3 -O3 -std=gnu89 -Wextra

#ESEGUZIONE MAIN
corrupt:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) corrupt $(dirOut)corrupt.bmp $(concat_images) 50
gray:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) gray $(dirOut)gray.bmp $(concat_images)
brighten:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) brighten $(dirOut)brighten.bmp $(concat_images) 50
blend:
	$(MemCheck) ./$(outputMain) $(img2) $(img3) blend $(dirOut)blend.bmp $(concat_images) 0 0.5
sharp:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) sharp $(dirOut)sharp.bmp $(concat_images)
edge:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) edge $(dirOut)edge.bmp $(concat_images)
emboss:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) emboss $(dirOut)emboss.bmp $(concat_images)
avg:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) avg $(dirOut)avg.bmp $(concat_images) 5
gauss:
	$(MemCheck) ./$(outputMain) $(img1) $(img2) gauss $(dirOut)gauss.bmp $(concat_images) 3 5


#DIPENDENZE
ip_lib.o: ip_lib.c ip_lib.h bmp.h
	gcc ip_lib.c -o ip_lib.o -c -Wall -lm -g -O1 --ansi --pedantic -g3 -O3 -std=gnu89 -Wextra

bmp.o: bmp.c bmp.h
	gcc bmp.c -o bmp.o -Wall -c 


#CLEAN
clean:
	@rm -f bmp.o ip_lib.o $(outputMain)
	@rm -r $(dirOut)
	@mkdir $(dirOut)
