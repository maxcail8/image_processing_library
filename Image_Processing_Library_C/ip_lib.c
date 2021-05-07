
/*
 ***Created by Lohackers***
    @Massimo Cailotto: 880763
    @Giovanni Costa: 880892
    @Matteo Minardi: 880895
    @Andrea Munarin: 879607
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"


/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

/*Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
  Inoltre crea un vettore di stats per contenere le statische sui singoli canali*/
ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned  int k, float v){
    
    ip_mat *nuova = (ip_mat *)malloc(sizeof(ip_mat));
    stats *st = (stats *)malloc(sizeof(stats)*k);
    
    /*Generazione righe*/
    float ***newData=(float***)malloc(h*sizeof(float**)); /*Puntatore ad un vettore di puntatori, ogni elemento punta ad un vettore di puntatori*/
    
    unsigned int r, c, p; /*profondità, riga, colonna*/

    if(nuova==NULL || st==NULL || newData==NULL){
        printf("Creazione fallita!\n");
        exit(1);
    }
    
    /*Inizializzazione dimensioni*/
    nuova->h = h;
    nuova->w = w;
    nuova->k = k;
    
    nuova->stat= st;
    
    /*Inizializzazione Data*/
    nuova->data = newData;
    for (r=0; r<h; r++){
        /*Generazione colonne*/
        newData[r]=(float **) malloc (w*sizeof(float*)); /*Puntatore ad un vettore di puntatori, ogni elemento punta a float*/
        if(newData[r]==NULL){
            printf("Creazione fallita!\n");
            exit(1);
        }
        
        for (c=0; c<w; c++){
            /*Generazione profondità*/
            newData[r][c]=(float *) malloc (k*sizeof(float)); /*Puntarore a float*/
            if(newData[r][c]==NULL){
                printf("Creazione fallita!\n");
                exit(1);
            }
            
            for(p=0; p<k; p++){
                set_val(nuova, r, c, p, v);
            }
        }
    }

    /*Inizializzazione del vettore di stats*/
    compute_stats(nuova);

    return nuova;    
}



/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){

    if(a!=NULL){

        unsigned int r, c;
        /*Libera il vettore di stats*/
        free(a->stat);
        
        for(r=0; r<a->h; r++){
            for(c=0; c<a->w; c++){
                /*Libera la memoria allocata al puntatore a float*/
                free(a->data[r][c]);
            }
            /*Libera la memoria allocata al puntatore "colonna"*/
            free(a->data[r]);
        }
        /*Libera la memoria allocata al puntatore della matrice 3D "data"*/
        free(a->data);
        
        /*Libera la struttura ip_mat allocata*/
        free(a);
    }

}



/* Restituisce il valore in posizione i, j, k */
float get_val(ip_mat * a, unsigned int i, unsigned int j, unsigned int k){

    if(a!=NULL && i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!\n");
        exit(1);
    }
}



/* Setta il valore in posizione i, j, k a v*/
void set_val(ip_mat * a, unsigned int i, unsigned int j, unsigned int k, float v){
    
    if(a!=NULL && i<a->h && j<a->w && k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!\n");
        exit(1);
    }
}



/* Calcola il valore minimo, il massimo e la media per ogni canale
   e li salva dentro la struttura ip_mat stats */
void compute_stats(ip_mat * t){
    
    if(t!=NULL){
        unsigned int r, c, p;
        float min, max, somma, mean;
    
        /*Pattern del minimo/massimo*/
        min=get_val(t, 0, 0, 0);
        max=min;
        somma=0;
        
        for(p=0; p<t->k; p++){
            for(r=0; r<t->h; r++){
                for(c=0; c<t->w; c++){
                    float tmp = get_val(t, r, c, p);
                    if(min>tmp)
                        min=tmp;
                    if(max<tmp)
                        max=tmp;
                    somma += tmp;
                }
            }

            /*Assegnazione degli stati*/
            t->stat[p].min=min;
            t->stat[p].max=max;
            mean = somma/(t->w*t->h); /*Calcolo della media*/
            t->stat[p].mean=mean;

            /*Inizializzazione variabili per il canale successivo (ad eccezione dell'ultimo canale/dimensione)*/
            if(p<t->k-1){
                min=get_val(t, 0, 0, p+1);
                max=min;
                somma=0;
            }
        }
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}



/* Inizializza una ip_mat con dimensioni w h e k.
   Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var){

    if(t!=NULL){
        unsigned int r, c, p;

        for (r=0; r<t->h; r++){
            for (c=0; c<t->w; c++){
                for(p=0; p<t->k; p++){
                    float n = get_normal_random(mean, var);
                    set_val(t, r, c, p, n);
                }
            }
        }

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

    /*compute_stats(t);*/
}



/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in){
    
    if(in!=NULL){
        unsigned int r, c, p;
        ip_mat * nuova =  ip_mat_create(in->h, in->w, in->k, 0.0);
        
        for(r=0; r<in->h; r++)
            for(c=0; c<in->w; c++)
                for(p=0; p<in->k; p++)
                    set_val(nuova, r, c, p, get_val(in, r, c, p));
        
        compute_stats(nuova);

        return nuova;

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}



/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    
    /*Controllo correttezza parametri dati*/
    if (t!=NULL && row_start<row_end && col_start<col_end && row_end<=t->h && col_end<=t->w){
        
        ip_mat * nuova;
        unsigned int r, c, p;
        nuova =  ip_mat_create(row_end-row_start, col_end-col_start, t->k, 0.0);
        
        for(r=0; r<row_end-row_start; r++)
            for(c=0; c<col_end-col_start; c++)
                for(p=0; p<t->k; p++)
                    set_val(nuova, r, c, p, get_val(t, r+row_start, c+col_start, p));

        return nuova;
    
    }else{
        printf("\nDimensioni date non valide\n");
        exit(1);
    }
    
}



/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
*/
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){

    if(a!=NULL && b!=NULL){
        ip_mat * nuova=NULL;
        unsigned int r, c, p;
        
        switch(dimensione){
            case 0: /*Concatenazine matrice rispetto alle righe*/
                if(a->w==b->w && a->k==b->k){
                    unsigned int tmp = a->h;
                    nuova = ip_mat_create((a->h)+(b->h), a->w, a->k, 0.0);
                    for(r=0; r<nuova->h; r++){
                        for(c=0; c<nuova->w; c++){
                            for(p=0; p<nuova->k; p++){
                                if(r<tmp) /*Scorrimento prima matrice*/
                                    set_val(nuova, r, c, p, get_val(a, r, c, p));
                                else /*Scorrimento seconda matrice*/
                                    set_val(nuova, r, c, p, get_val(b, (r-tmp), c, p));
                            }
                        }
                    }
                }
                else{
                    printf("Dimensioni non valide!\n");
                    exit(1);
                }

            break;
            
            
            case 1: /*Concatenzazione matrice rispetto alle colonne*/
                if(a->h==b->h && a->k==b->k){
                    unsigned int tmp = a->w;
                    nuova = ip_mat_create(a->h, (a->w)+(b->w), a->k, 0.0);
                    for(r=0; r<nuova->h; r++){
                        for(c=0; c<nuova->w; c++){
                            for(p=0; p<nuova->k; p++){
                                if(c<tmp) /*Scorrimento prima matrice*/
                                set_val(nuova, r, c, p, get_val(a, r, c, p));
                                else /*Scorrimento secondo matrice*/
                                    set_val(nuova, r, c, p, get_val(b, r, (c-tmp), p));
                            }
                        }
                    }
                }
                else{
                    printf("Dimensioni non valide!\n");
                    exit(1);
                }

            break;
            
            
            case 2: /*Concatenazione matrice rispetto ai canali*/
                if(a->h==b->h && a->w==b->w){
                    unsigned int tmp = a->k;
                    nuova = ip_mat_create(a->h, a->w, (a->k)+(b->k), 0.0);
                    for(r=0; r<nuova->h; r++){
                        for(c=0; c<nuova->w; c++){
                            for(p=0; p<nuova->k; p++){
                                if(p<tmp) /*Scorrimento prima matrice*/
                                    set_val(nuova, r, c, p, get_val(a, r, c, p));
                                else /*Scorrimento seconda matrice*/
                                    set_val(nuova, r, c, p, get_val(b, r, c, (p-tmp)));
                            }
                        }
                    }
                }
                else{
                    printf("Dimensioni non valide!\n");
                    exit(1);
                }

            break;

            default:
                printf("\nValore %d non disponibile\n", dimensione);
                exit(1);
            break;
        }


        return nuova;

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}



/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    
    /*Verifica se possibile effettuare l'operazione*/
    if(a!=NULL && b!=NULL && a->h==b->h && a->w==b->w && a->k==b->k){
        ip_mat * nuova=NULL;
        unsigned int r, c, p;
        nuova= ip_mat_create(a->h, a->w, a->k, 0.0);
        for(r=0; r<nuova->h; r++)
            for(c=0; c<nuova->w; c++)
                for(p=0; p<nuova->k; p++)
                    set_val(nuova, r, c, p, (get_val(a, r, c, p)+get_val(b, r, c, p)));
        
        return nuova;
    
    }else{
        printf("\nLe dimensioni delle matrici sono diverse! Failed\n");
        exit(1);
    }
    
}


/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    
    /*Verifica se possibile effettuare l'operazione*/
    if(a!=NULL && b!=NULL && a->h==b->h && a->w==b->w && a->k==b->k){
        ip_mat * nuova=NULL;
        unsigned int r, c, p;
        nuova= ip_mat_create(a->h, a->w, a->k, 0.0);
        for(r=0; r<nuova->h; r++)
            for(c=0; c<nuova->w; c++)
                for(p=0; p<nuova->k; p++)
                    set_val(nuova, r, c, p, (get_val(a, r, c, p)-get_val(b, r, c, p)));
        
        return nuova;
        
    }else{
        printf("\nLe dimensioni delle matrici sono diverse! Failed\n");
        exit(1);
    }
    
}


/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){

    if(a!=NULL){
        ip_mat * nuova=NULL;
        unsigned int r, col, p;
        
        nuova= ip_mat_create(a->h, a->w, a->k, 0.0);
        for(r=0; r<nuova->h; r++)
            for(col=0; col<nuova->w; col++)
                for(p=0; p<nuova->k; p++)
                    set_val(nuova, r, col, p, (get_val(a, r, col, p)*c));

        /*compute_stats(nuova);*/
        return nuova;
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}


/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){

    if(a!=NULL){
        ip_mat * nuova=NULL;
        unsigned int r, col, p;
        
        nuova= ip_mat_create(a->h, a->w, a->k, 0.0);
        for(r=0; r<nuova->h; r++)
            for(col=0; col<nuova->w; col++)
                for(p=0; p<nuova->k; p++)
                    set_val(nuova, r, col, p, (get_val(a, r, col, p)+c));

        return nuova;
    
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
    
}


/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    
    /*Verifica se possibile effettuare l'operazione*/
    if(a!=NULL && b!= NULL && a->h==b->h && a->w==b->w && a->k==b->k){
        ip_mat * nuova = NULL;
        ip_mat * tmp = NULL;
        /*Somma i valori delle due matrici cella per cella e li divide per 2*/
        tmp=ip_mat_sum(a, b);
        nuova=ip_mat_mul_scalar(tmp, 0.5);
        ip_mat_free(tmp);

        return nuova;
        
    }else{
        printf("\nLe dimensioni delle matrici sono diverse! Failed\n");
        exit(1);
    }

}



/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    
    if(in!=NULL){
        unsigned int r, c, p;
        float somma, media;
        ip_mat * nuova= ip_mat_create(in->h, in->w, in->k, 0.0);
        
        /*Pattern della media*/
        somma = 0.0;
        for (r=0; r<in->h; r++){
            for(c=0; c<in->w; c++){
                for(p=0; p<in->k; p++){
                    somma+=get_val(in, r, c, p);
                }
                media=somma/in->k;
                for(p=0; p<in->k; p++){
                    set_val(nuova, r, c, p, media);
                }
                somma=0.0;
            }
        }

        return nuova;
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}



/* Effettua la fusione (combinazione convessa) di due immagini
   FORMULA:
    Blend = alpha * A + (1-alpha)* B;
    alpha è in range [0,1]
*/
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){

    /*Verifica se possibile effettuare l'operazione*/
    if (a!=NULL && b!=NULL && a->h == b->h && a->w == b->w && a->k == b->k){
        ip_mat * nuova = NULL;
        /*Controllo correttezza parametri*/
        if(alpha >= 0 && alpha <=1){
            ip_mat *aux1 = ip_mat_mul_scalar(a, alpha); /*alpha * A */
            ip_mat *aux2 = ip_mat_mul_scalar(b, (1-alpha)); /*(1-alpha)* B */
            nuova = ip_mat_sum(aux1, aux2);
        
            ip_mat_free(aux1);
            ip_mat_free(aux2);
        }

        return nuova;
        
    }else{
        printf("\nLe dimensioni delle immagini diverse! Failed\n");
        exit(1);
    }
        
}



/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    
    if(a!=NULL){
        ip_mat * nuova = ip_mat_add_scalar(a, bright);

        return nuova;
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}



/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 *
 * out = a + gauss_noise*amount
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){

    if(a!=NULL){
        ip_mat * gauss = ip_mat_create(a->h, a->w, a->k, 0.0);
        ip_mat * nuova;
        
        /*Inizializzazione matrice gaussiana (gauss_noise)*/
        ip_mat_init_random(gauss, 0, amount/2.0); /* gauss_noise/2-->copre 95,45% della distribuzione normale(centro della gaussiana in 0)*/
        nuova = ip_mat_sum(a, gauss); /*a + gauss_noise*amount*/
        
        ip_mat_free(gauss);

        return nuova;
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}



/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){

    if(a!=NULL && f!=NULL){
        unsigned int i, j, r, c, p;
        unsigned int dim_filtro = f->h; /*Dimensione del filtro*/
        float somma=0;
        ip_mat *a_new; /*Matrice con padding*/
        ip_mat *nuova; /*Matrice finale*/
        
        if(f->h != f->w){ /*Controllo che il filtro sia quadrato*/
            printf("Il filtro deve essere quadrato! Errore\n");
            exit(1);
        }

        /*Calcolo padding da aggiungere tramite formula, per una dimensione corretta dell'ip_mat finale*/
        a_new = ip_mat_padding(a, (dim_filtro-1)/2, (dim_filtro-1)/2);
    
        nuova = ip_mat_create(a->h, a->w, a->k, 0.0);
        
        for(p=0; p<a_new->k; p++){ /*Scorrimento canale per canale*/

            /*Scorrimento per ottenere tutte le sottomatrici*/
            for(i=0; i<=a_new->h-dim_filtro; i++){
                for(j=0; j<=a_new->w-dim_filtro; j++){
                    
                    /*Sottomatrice da moltiplicare*/
                    ip_mat *sub_mat = ip_mat_subset(a_new, i, dim_filtro+i, j, dim_filtro+j);
                    
                    /*Moltiplicazione sottomatrice per filtro*/
                    for(c=0; c<sub_mat->w; c++)
                        for(r=0; r<sub_mat->h; r++)
                            if(f->k>1) /*Se il filtro ha più canali, esegue l'operazione con il canale specifico*/
                                somma += (get_val(sub_mat, r, c, p)*get_val(f, r, c, p));
                            else
                                somma += (get_val(sub_mat, r, c, p)*get_val(f, r, c, 0));                   
                    
                    /*Risultato calcolo filtro*sottomatrice e inserimento valore nella matrice finale*/
                    set_val(nuova, i, j, p, somma);
                    somma=0.0;
                    ip_mat_free(sub_mat);
                }
            }
        }
        
        ip_mat_free(a_new);
        
        return nuova;
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}


/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
*/
ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){

    if(a!=NULL){
        unsigned int r, c, p;
        ip_mat *nuova = ip_mat_create(a->h+2*pad_h, a->w+2*pad_w, a->k, 0.0);
        
        for(r=pad_h; r<nuova->h - pad_h; r++)
            for(c=pad_w; c<nuova->w - pad_w; c++)
                for(p=0; p<nuova->k; p++)
                    set_val(nuova, r, c, p, get_val(a, r-pad_h, c-pad_w, p));

        return nuova;
        
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }
}


/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter(){
    ip_mat *nuova = ip_mat_create(3, 3, 3, 0.0);
    unsigned int p;

    for(p=0; p<nuova->k; p++){
        set_val(nuova, 0, 1, p, -1.0);
        set_val(nuova, 1, 0, p, -1.0);
        set_val(nuova, 1, 1, p, 5.0);
        set_val(nuova, 1, 2, p, -1.0);
        set_val(nuova, 2, 1, p, -1.0);
    }

    return nuova;
}


/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    unsigned int p;
    ip_mat *nuova = ip_mat_create(3, 3, 3, -1.0);

    for(p=0; p<nuova->k; p++){
        set_val(nuova, 1, 1, p, 8.0);
    }

    return nuova;
}


/* Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter(){
    unsigned int p;
    ip_mat *nuova = ip_mat_create(3, 3, 3, 1.0);

    for(p=0; p<nuova->k; p++){
        set_val(nuova, 0, 0, p, -2.0);
        set_val(nuova, 0, 1, p, -1.0);
        set_val(nuova, 0, 2, p, 0.0);
        set_val(nuova, 1, 0, p, -1.0);
        set_val(nuova, 2, 0, p, 0.0);
        set_val(nuova, 2, 2, p, 2.0);
    }

    return nuova;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
    float c = (1.0)/(w*h); /*FORMULA*/
    ip_mat *nuova = ip_mat_create(h, w, k, c);

    return nuova;
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma){
    int cx, cy;
    unsigned int r, c, p;
    float somma=0.0;
    ip_mat *nuova=ip_mat_create(h, w, k, 0.0);

    /*Coordinate del centro del kernel (filtro)*/
    cx=w/2;
    cy=h/2;

    /*Somma di tutti i valori dei pixel dell'immagine a cui applichiamo la funzione G(x, y), 
    per la normalizzazione della matrice finale*/
    for(r=0; r<h; r++){
        for(c=0; c<w; c++){

                /*Calcolo della distanza dal centro*/
                int x = r-cx; 
                int y = c-cy;
                
                /*FUNZIONE G(x, y) data*/
                float esp = exp(-((pow(x, 2.0)+pow(y, 2.0))/(2.0*pow(sigma, 2.0))));
                float fattore = (1.0/(2.0*PI*pow(sigma, 2.0)));
                
                somma += (fattore*esp);
        }
    }

    /*Calcolo dei valori dei pixel dell'immagine a cui applichiamo la funzione G(x, y), 
    poi la dividiamo per la somma calcolata prima per normalizzarli*/
    for(r=0; r<h; r++){
        for(c=0; c<w; c++){
            for(p=0; p<k; p++){

                /*Calcolo della distanza dal centro*/
                int x = r-cx;
                int y = c-cy;

                /*FUNZIONE G(x, y) data*/
                float esp = exp(-((pow(x, 2.0)+pow(y, 2.0))/(2.0*pow(sigma, 2.0))));
                float fattore = (1.0/(2.0*PI*pow(sigma, 2.0)));
                float f = (fattore*esp);

                /*Inserimento valori normalizzati*/
                set_val(nuova, r, c, p, (f/somma));
            }
        }
    }
    

    return nuova;
}



/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula (valore-min)/(max-min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat * t, float new_max){
    
    if(t!=NULL){
        unsigned int p, r, c;
        
        compute_stats(t);
        for(p=0; p<t->k; p++){
            float min, max;
            min=t->stat[p].min;
            max=t->stat[p].max;
            for(r=0; r<t->h; r++)
                for(c=0; c<t->w; c++){
                    float val = get_val(t, r, c, p);
                    float scal = (val-min)/(max-min);
                    set_val(t, r, c, p, scal*new_max);               
                }
        }

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}



/* Nell'operazione di clamping i valori <low si convertono in low e i valori>high in high.*/
void clamp(ip_mat * t, float low, float high){
    
    if(t!=NULL){
        unsigned int p, r, c;
        
        for(r=0; r<t->h; r++)
            for(c=0; c<t->w; c++)
                for(p=0; p<t->k; p++){

                    float val=get_val(t, r, c, p);
                    
                    if (val > high)
                        set_val(t, r, c, p, high);
                    else if (val < low)
                        set_val(t, r, c, p, low);
                }
                
    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}


void ip_mat_show(ip_mat * t){

    if (t != NULL) {
        unsigned int i,l,j;
        printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
        for (l = 0; l < t->k; l++) {
            printf("Slice %d\n", l);
            for(i=0;i<t->h;i++) {
                for (j = 0; j < t->w; j++) {
                    printf("%f ", get_val(t,i,j,l));
                }
                printf("\n");
            }
            printf("\n");
        }

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}


void ip_mat_show_stats(ip_mat * t){

    if (t != NULL) {
        unsigned int k;

        compute_stats(t);

        for(k=0; k<t->k; k++){
            printf("Channel %d:\n", k);
            printf("\t Min: %f\n", t->stat[k].min);
            printf("\t Max: %f\n", t->stat[k].max);
            printf("\t Mean: %f\n", t->stat[k].mean);
        }

    }else{
        printf("L'elemento passato è NULL\n");
        exit(1);
    }

}


ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;
    unsigned char R,G,B;
    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}


Bitmap * ip_mat_to_bitmap(ip_mat * t){
    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}


float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}