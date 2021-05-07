#include <stdio.h>
#include <string.h>
#include "ip_lib.h"
#include "bmp.h"

int main(){

    /*Create*/
    /*ip_mat *rnd = ip_mat_create(5, 5, 3, 0.0);
    ip_mat *mat = ip_mat_create(5, 5, 3, 1.0);*/

    /*Inizializzazione Random*/
    /*ip_mat_init_random(rnd, 100, 1);
    ip_mat_show(rnd);
    printf("------------------------\n");*/
    /*Stats*/
    /*compute_stats(rnd);
    ip_mat_show_stats(rnd);
    printf("------------------------\n");*/

    /*Riduciamo matrice mat-> h=2, w=5, k=3*/
    /*mat=ip_mat_subset(mat, 0, 2, 0, 5);
    ip_mat_show(mat);
    printf("------------------------\n");*/

    /*Concateniamo matrice mat con rnd*/
    /*ip_mat *conc = ip_mat_concat(rnd, mat, 0);
    ip_mat_show(conc);
    printf("------------------------\n");*/

    /*Liberiamo la memoria*/
    /*ip_mat_free(rnd);
    ip_mat_free(mat);
    ip_mat_free(conc);*/

    /*PARTE 1*/
    /*ip_mat *a = NULL;
    ip_mat *b = NULL;
    ip_mat *c = NULL;
    printf("\n\nPARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT\n\n");
    a = ip_mat_create(5, 5, 3, 0.0);
    b = ip_mat_create(5, 5, 3, 0.0);
    
    ip_mat_init_random(a, 20, 3);
    printf("A:\n");
    ip_mat_show(a);
    printf("------------------------\n");
    
    ip_mat_init_random(b, 60, 3);
    printf("B:\n");
    ip_mat_show(b);
    printf("------------------------\n");

    c = ip_mat_sum(a, b);
    printf("Sum:\n");
    ip_mat_show(c);
    printf("------------------------\n");
    ip_mat_free(c);

    c = ip_mat_sub(a, b);
    printf("Sub:\n");
    ip_mat_show(c);
    printf("------------------------\n");
    ip_mat_free(c);

    c = ip_mat_mul_scalar(a, 7.0);
    printf("Mul_scal:\n");
    ip_mat_show(c);
    printf("------------------------\n");
    ip_mat_free(c);

    c = ip_mat_add_scalar(a, 7.0);
    printf("Add_scal:\n");
    ip_mat_show(c);
    printf("------------------------\n");
    ip_mat_free(c);

    c = ip_mat_mean(a,b);
    printf("Mat_mean:\n");
    ip_mat_show(c);
    printf("------------------------\n");
    
    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(c);*/
    
    /*PARTE 2*/
    /*ip_mat *a = NULL;
    ip_mat *b = NULL;
    ip_mat *c = NULL;
    
    Bitmap *flowerImg;
    Bitmap * mongolfImg;
    Bitmap *ris;

    flowerImg = bm_load("img/original/flower.bmp");
    mongolfImg = bm_load("img/original/mongolfiere.bmp");
    
    
    printf("PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI\n\n");
    
    a = bitmap_to_ip_mat(flowerImg);
    b = bitmap_to_ip_mat(mongolfImg);*/

    /*GRAYSCALE*/
    /*c=ip_mat_to_gray_scale(a);
    ris = ip_mat_to_bitmap(c);
    bm_save(ris, "img/gray.bmp");
    bm_free(ris);
    ip_mat_free(c);*/
    
    /*BLAND*/
    /*c = ip_mat_blend(a, b, 0.5);
    ris = ip_mat_to_bitmap(c);
    bm_save(ris, "img/blend.bmp");
    bm_free(ris);
    ip_mat_free(c);*/

    /*BRIGHTEN*/
    /*c = ip_mat_brighten(a, 70.0);
    ris = ip_mat_to_bitmap(c);
    bm_save(ris, "img/brighten.bmp");
    bm_free(ris);
    ip_mat_free(c);*/
    
    /*CORRUPT*/
    /*c = ip_mat_corrupt(a, 100.0);
    ris = ip_mat_to_bitmap(c);
    bm_save(ris, "img/corrupt.bmp");
    bm_free(ris);
    
    bm_free(flowerImg);
    bm_free(mongolfImg);
    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Corrupt created...\n\n");*/
    
    /*PARTE 3*/
    /*ip_mat *a = NULL;
    ip_mat *b = NULL;
    ip_mat *c=NULL;
    
    Bitmap *flowerImg;
    Bitmap *ris;

    flowerImg = bm_load("img/original/flower.bmp");
    
    a=bitmap_to_ip_mat(flowerImg);*/
    
    /*SHARPEN*/
    /*c= create_sharpen_filter();
    printf("MATRICE FILTRO:\n");
    ip_mat_show(c);

    b = ip_mat_convolve(a, c);
    
    ris = ip_mat_to_bitmap(b);
    bm_save(ris, "img/results/sharpen.bmp");
    bm_free(ris);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Sharpen created...\n\n");*/

    /*EDGE*/
    /*c= create_edge_filter();
    printf("MATRICE FILTRO:\n");
    ip_mat_show(c);

    b = ip_mat_convolve(a, c);

    ris = ip_mat_to_bitmap(b);
    bm_save(ris, "img/results/edge.bmp");
    bm_free(ris);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Edge created...\n\n");*/

    /*EMBOSS*/
    /*c= create_emboss_filter();
    printf("MATRICE FILTRO:\n");
    ip_mat_show(c);

    b = ip_mat_convolve(a, c);
    
    ris = ip_mat_to_bitmap(b);
    bm_save(ris, "img/results/emboss.bmp");
    bm_free(ris);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Emboss created...\n\n");*/

    /*AVERAGE*/
    /*c= create_average_filter(a->w, a->h, a->k);
    printf("MATRICE FILTRO:\n");
    ip_mat_show(c);

    b = ip_mat_convolve(a, c);
    
    ris = ip_mat_to_bitmap(b);
    bm_save(ris, "img/results/average.bmp");
    bm_free(ris);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Average created...\n\n");*/    

    /*GAUSSIAN*/
    /*c= create_gaussian_filter(13, 13, 3, 5);
    printf("MATRICE FILTRO:\n");
    ip_mat_show(c);

    b = ip_mat_convolve(a, c);
    clamp(b, 0, 255);
    
    ris = ip_mat_to_bitmap(b);
    bm_save(ris, "img/gaussian.bmp");
    bm_free(ris);
    ip_mat_free(b);
    ip_mat_free(c);
    printf("Gaussian created...\n\n");*/

    
    /*TODO: DUBBI:
        ip_mat_init_random --> formula
        corrupt->mean e var come sono legate a amount?
        differenza tra rescale e clamp, quando è preferibile usare una piuttosta che un l'atra?
        OK: sulle immagine con convolve viene un bordo, questo è dato da cosa?
    */


    return 0;
}
