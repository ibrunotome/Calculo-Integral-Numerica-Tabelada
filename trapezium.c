/*
    Nome: trapezium.c
    Copyright (c) 2014 IFMG. All rights reserved.

    =================== Membros ===================

    Bruno Tomé - 0011254 - ibrunotome@gmail.com
    Matheus Calixto - 0011233 - calixtinn@gmail.com

    =========== Instruções de compilação ==========

    Abra o Terminal e digite:

    cd <DIRETÓRIO>
    gcc trapezium.c -otrapezium.bin
    ./trapezium.bin <ARQUIVO ENTRADA.txt> <ARQUIVO SAIDA.r>
	Rscript <ARQUIVO SAIDA.r>

    ========== Ambiente de Desenvolvimento ========

    Bruno Tomé

    Sistema Operacional: OS X 10.9.4
    Hardware: Core i5 3ª Geração 2.5 GHz| 16 GB RAM 1600Mhz
    Desenvolvido no editor de textos Sublime Text 3
    Compilado no Terminal do OS X

    GCC

    Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
    Target: x86_64-apple-darwin13.3.0
    Thread model: posix

    Matheus Calixto

    Sistema Operacional: Linux Ubuntu 13.10 
    Hardware: Core i5 3ª Geração 1.8 ~ 2.5 GHz| 8GB RAM 1600 Mhz
    Desenvolvido no NetBeans 8.0
    Compilado no Terminal do Ubuntu.

    GCC (Ubuntu 4.8.2-19ubuntu1) 4.8.2

    ============ Objetivo do programa ============

    Criação de um software de interpolação e integração numérica
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int number_trapezium = 0, number_interpol = 0, cont = 0, contadorx = 0, contadorp = 0, grau = 0;
float limdown = 0, limtop = 0, valuesofx[256], valuesofy[256], valuesofp[256], auxinterp[256], pinterpolados[256], interpolador[256], valor_integral;
char valuestemp[256], figura[256];

/*
    Requisito 6 - Avaliaçãode pn(x) usando Horner
    Implementação Método Horner
*/

float metodohorner(float c[],float a){
    float y;
    int i = 1, n;
    n = grau;
    y = c[0];
    
    do{
        y =  y * a + c[i];
        i++;
    }while(i<(n+1));
    
    return y; 
}

// Requisito 7 - Cálculo da Área do Trapézio
// Função areatrapezio calcula a área do trapézio a partir da
// fórmula = Base maior + base menor * altura / 2

float areatrapezio(float B, float b, float h){
    return (((B+b)*h)/2);
}

// Fim areatrapezio

// Requisito 8 - Cálculo da Integral definida no intervalo [a,b] utilizando a regra do trapézio.

float calculointegral(float a, float b, int m) { // a = Limite inferior // b = limite superior // m = quantidade de trapézios.

    float somatorio, altura, valoresx[m + 1], valoresy[m + 1], aux;
    int i = 0, k = 0;

    altura = (b - a) / m; //Calcula o tamanho de cada subintervalo (altura dos trapézios))

    valoresx[0] = a;

    for (i = 1; i < m + 1; i++) {
        valoresx[i] = (valoresx[i - 1] + altura); //Adiciona num vetor todos os valores de x que estão no subintervalo

    }

    for (i = 0; i < m + 1; i++) {

        valoresy[i] = metodohorner(interpolador, valoresx[i]); // Salva num vetor Y os valores dos X's do subintervalo.
    }

    for (i = 0; i < m; i++) {
        somatorio = somatorio + areatrapezio(valoresy[i], valoresy[i + 1], altura);
    }

    return (somatorio);
}

/* Função para calcular a potência de um número */

float potencia(float a, int b){
    float pot = 1;
    int i = 0;
    
    for (i = 0; i < b; i++){
        pot *= a;
    }
    
    return pot;
}

/* Requisito - 04  Determinação do Polinômio Interpolador.
    -Criação da Matriz de Vandermond.
    -Resolução do sistema por Eliminação de Gauss
*/

void gausselimination(float a[], float b[]) {

    int l, c, indice, i, k, j, aux = -1, n;
    float matrixvanderm[contadorx][contadorx], vetorcoef[contadorx], x[contadorx], m, soma;

    // - Tópico 1 - Criação da Matriz de Vandermond.

    for (i = 0; i < contadorx; i++) { // Preenche a variável Vetor Coef com os valores do vetor Y
        vetorcoef[i] = b[i];
    }

    for (i = 0; i < contadorx; i++) { // Preenche os alementos da matriz[i][0] com 1.
        matrixvanderm[i][0] = 1;
    }

    for (l = 0; l < contadorx; l++) {

        indice = 1;
        aux++;

        for (c = 1; c < contadorx; c++) {
            matrixvanderm[l][c] = potencia(a[aux], indice);
            indice++;
        }
    }

    // - Tópico 2 - Requisito 04 - Resolução do Sistema por Eliminação de Gauss.

    n = contadorx; // n é a ordem da matriz.


    for (k = 0; k < (n - 1); k++) {
        //Metodo da eliminacao de Gauss (Triangularizacao)

        for (i = k + 1; i < n; i++) {
            m = matrixvanderm[i][k] / matrixvanderm[k][k]; //calcula o multiplicador

            for (j = 1; j < n; j++) //comeca do elemento apos o valor a ser zerado
            {
                matrixvanderm[i][j] = matrixvanderm[i][j] - (m * matrixvanderm[k][j]);
            }
            vetorcoef[i] = vetorcoef[i] - (m * vetorcoef[k]);
        }
    }

    //Requisito 05 - Calcula o vetor solução x do sistema linear Ax=b usando Substituição Retroativa

    auxinterp[n] = vetorcoef[n] / matrixvanderm[n][n];

    for (i = n - 1; i >= 0; i--) {
        soma = 0;
        for (j = (i + 1); j < n; j++) {
            soma = soma + (matrixvanderm[i][j] * auxinterp[j]);
        }
        auxinterp[i] = (vetorcoef[i] - soma) / matrixvanderm[i][i];
    }

} // Fim gausselimination

// Lendo arquivo de entrada

void lendoarquivo(char **argv) {
    
    int i = 0, k = 0, u = 0, tamanho = 0, tamanhox = 0,tamanhoy = 0, tamanhop = 0, tamanhot = 0;
    char *veriflinha, linha[256], pvalues[256], inferior[256], superior[256], numtrapezios[256], npinterpol[256], tvalues[256], xvalues[256], yvalues[256];
    
    FILE *arquivoin;
    
    arquivoin = fopen(argv[1],"r");
    
    while (!feof(arquivoin)){
        veriflinha = fgets (linha, 256, arquivoin);
        if (veriflinha){   
            if (linha[0] == 'a'){
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    inferior[k]=linha[i];
                    k++;
                } // end for verificação 'a'              
            } // end if linha == 'a'
            
            else if (linha[0] == 'b'){
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++){
                    superior[k]=linha[i];
                    k++;
                } // end for verificação 'b'             
            } // end if linha == 'b'
            
            else if (linha[0] == 't'){
                tamanhot = strlen(linha);
                k = 0;
                for(i=2;i<(tamanhot-1);i++) {
                    numtrapezios[k]=linha[i];
                    k++;
                } // end for verificação 't'      
            } // end if linha == 't'
            
            else if (linha[0] == 'i'){
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++){
                    npinterpol[k]=linha[i];
                    k++;
                } // end for verificação 'i'
            } // end if linha == 'i'

            else if (linha[0] == 'x'){
                tamanhox = strlen(linha);
                k = 0;
                for(i=2;i<(tamanhox-1);i++){
                    xvalues[k] = linha[i];
                    k++;
                } // end for verificação 'x'
            } // end if verificação 'x'

            else if (linha[0] == 'y'){
                tamanhoy= strlen(linha);
                k = 0;
                for(i=2;i<(tamanhoy-1);i++){
                    yvalues[k] = linha[i];
                    k++;
                } // end for verificação 'y'
            } // end if verificação 'y'

            else if (linha[0] == 'p'){
                tamanhop = strlen(linha);
                k = 0;
                for(i=2;i<(tamanhop-1);i++){
                    pvalues[k] = linha[i];
                    k++;
                } // end for verificação 'p'
            } // end if verificação 'p'
        } // end veriflinha
    } // end while
    
    fclose(arquivoin); // Fechando o arquivo de Leitura.
    
    limdown = atof(inferior); // Salva o Limite Inferior em uma variável do tipo float
    limtop = atof(superior); // Salva o Limite Superior em uma variável do tipo float
    number_interpol = atoi(npinterpol); // Salva o Número de interpolações em uma variável do tipo int
    number_trapezium = atoi(numtrapezios); // Salva o Número de trapézios em uma variável do tipo int
    xvalues[tamanhox-3] = xvalues[tamanhox-3] + ' '; // Adiciona um espaço no final da string, para controle.
    yvalues[tamanhoy-3] = yvalues[tamanhoy-3] + ' '; // Adiciona um espaço no final da string, para controle.
    pvalues[tamanhop-3] = pvalues[tamanhop-3] + ' '; // Adiciona um espaço no final da string, para controle.
    i = 0;
    u = 0;
    cont = 0;
    contadorx = 0;

  // Salvar os valores de X, Y e P em vetores.

    // Salvando valores de X. 

    while (xvalues[i] != '\0') {

        if (xvalues[i] != ' ') {
            valuestemp[u] = xvalues[i];
            i++;
            u++;
        }
        else {
            valuesofx[cont] = atof(valuestemp);
            contadorx++;
            i++;
            u = 0;
            cont++;
            strcpy(valuestemp, "");
        }
    }

    // Salvando valores de Y.

    i = 0;
    u = 0;
    cont = 0;

    while (yvalues[i] != '\0') {

        if (yvalues[i] != ' ') {
            valuestemp[u] = yvalues[i];
            i++;
            u++;
        }

        else {
            valuesofy[cont] = atof(valuestemp);
            i++;
            u = 0;
            cont++;
            strcpy(valuestemp, "");
        }
    }

    // Salvando valores de P.

    i = 0;
    u = 0;
    cont = 0;
    contadorp = 0;
    strcpy(valuestemp, "");

    while (pvalues[i] != '\0') {

        if (pvalues[i] != ' ') {
            valuestemp[u] = pvalues[i];
            i++;
            u++;
        }

        else {
            valuesofp[cont] = atof(valuestemp);
            contadorp++;
            i++;
            u = 0;
            cont++;
            strcpy(valuestemp, "");
        }
    }

    grau = contadorx - 1; // Grau do Polinômio interpolador

} // Terminou função lendoarquivo

/*
    Criar e escrever no arquivo de saída, a função irá receber o argumento 2
    como parâmetro e fará a escrita no arquivo a partir desta e de outras
    funções chamadas no corpo da função.
*/

void escrevendoarquivo(char **argv){
    FILE *rscript;
    
    int i = 0, j = 2;

    strcpy(figura, argv[2]); // copiando o nome do arquivo de saída pra usar no nome da png

    figura[(strlen(figura)-1)] = '\0'; // tirando o "r" do final do parâmetro 2

    rscript=fopen(argv[2],"w+");
    
    // Impressão no arquivo

    fprintf(rscript, "######################################################################\n");
    fprintf(rscript, "# Script automatico gerado por 'trapezium', software de interpolacao #\n");
    fprintf(rscript, "# e integracao numerica                                              #\n");
    fprintf(rscript, "######################################################################\n");

    
    fprintf(rscript, "\n# Nome da figura");
    fprintf(rscript, "\nnome <- '%spng'",figura);


    fprintf(rscript, "\n\n# Dados tabelados");
    fprintf(rscript, "\nx.tab <- c(");
    fprintf(rscript, "%f",valuesofx[0]);
        for (i=1;i<contadorx;i++){
            fprintf(rscript, ", %f", valuesofx[i]);
        }
    fprintf(rscript, ");");

    fprintf(rscript, "\ny.tab <- c(");
    fprintf(rscript, "%f", valuesofy[0]);
        for (i=1;i<contadorx;i++){
            fprintf(rscript, ", %f", valuesofy[i]);
        }
    fprintf(rscript, ");");


    fprintf(rscript, "\n\n# Pontos interpolados, calculados pelo 'trapezium'");
    fprintf(rscript, "\nx.int <- c(");
    fprintf(rscript, "%f", valuesofp[0]);
        for (i=1;i<contadorp;i++){
            fprintf(rscript, ", %f", valuesofp[i]);
        }
    fprintf(rscript, ");");

    fprintf(rscript, "\ny.int <- c(");
    fprintf(rscript, "%f", pinterpolados[0]);
        for (i=1;i<contadorp;i++){
            fprintf(rscript, ", %f", pinterpolados[i]);
        }
    fprintf(rscript, ");");


    fprintf(rscript, "\n\n# Coeficientes do polinomio interpolador");
    fprintf(rscript, "\ncoef <- c(");
    fprintf(rscript, "%f", interpolador[contadorx-1]);
        for (i=contadorx-2;i>=0;i--){
            fprintf(rscript, ", %f", interpolador[i]);
        }
    fprintf(rscript, ");");


    fprintf(rscript, "\n\n# Numero de pontos da tabela");
    fprintf(rscript, "\nn.tab <- %d;", contadorx);


    fprintf(rscript, "\n\n# Numero de pontos a interpolar");
    fprintf(rscript, "\nn.int <- %d;", number_interpol);


    fprintf(rscript, "\n\n# Numero de trapezios");
    fprintf(rscript, "\nn.tpz <- %d;", number_trapezium);


    fprintf(rscript, "\n\n# Titulo");
    fprintf(rscript, "\ntitulo <- 'P(x) = %f", interpolador[contadorx-1]);
    fprintf(rscript, " + %f * x\n", interpolador[contadorx-2]);
        for (i=contadorx-3;i>1;i--){
            fprintf(rscript," + %f * x^%d",interpolador[i],j);
            j++;
        }
    fprintf(rscript, "\n + %f * x^%d",interpolador[0],j+1);
    fprintf(rscript, "'");
    
    /*
        Parte Estática do script R
        copiado do enunciado
    */

    fprintf(rscript, "\n\n# Calcula o valor interpolado para o pto x");
    fprintf(rscript, "\npolinomio <- function(x, coef, n){");
    fprintf(rscript, "\n  resultado <- 0;");
    fprintf(rscript, "\n  for(i in 1:n){");
    fprintf(rscript, "\n      resultado <- resultado + coef[i]*(x^(i-1));");
    fprintf(rscript, "\n  }");
    fprintf(rscript, "\n  return(resultado);");
    fprintf(rscript, "\n}");

    fprintf(rscript, "\n\n#");
    fprintf(rscript, "\n# Aqui comecam os comandos para plotar os resultados");
    fprintf(rscript, "\n#");

    fprintf(rscript, "\n\n# Cria o arquivo .png");
    fprintf(rscript, "\npng(nome);");

    fprintf(rscript, "\n\n# Gerando figura com 100 pontos");
    fprintf(rscript, "\ngap <- (max(x.tab) - min(x.tab)) / 100;");
    fprintf(rscript, "\nx <- seq(min(x.tab), max(x.tab), gap);");
    fprintf(rscript, "\ny <- polinomio(x, coef, n.tab);");
    fprintf(rscript, "\nplot(x,y,type='l', main=titulo);");

    fprintf(rscript, "\n\n# Plota os trapezios");
    fprintf(rscript, "\nh <- (max(x.tab) - min(x.tab)) / n.tpz;");
    fprintf(rscript, "\nxp <- seq(min(x.tab), max(x.tab), h);");
    fprintf(rscript, "\nyp <- polinomio(xp, coef, n.tab);");
    fprintf(rscript, "\nfor(i in 1:(n.tpz)){");
    fprintf(rscript, "\n  polygon(c(xp[i], xp[i], xp[i+1], xp[i+1], xp[i]), c(0, yp[i], yp[i+1], 0, 0), col='yellow', border='black', lty=2, lwd=1.3);");
    fprintf(rscript, "\n}");

    fprintf(rscript, "\n\n# Pontos da tabela");
    fprintf(rscript, "\nfor(i in 1:n.tab){");
    fprintf(rscript, "\n  points(x.tab[i], y.tab[i], col='red', pch=19);");
    fprintf(rscript, "\n}");

    fprintf(rscript, "\n\n# Pontos interpolados");
    fprintf(rscript, "\nfor(i in 1:n.tab){");
    fprintf(rscript, "\n  points(x.int[i], y.int[i], col='blue', pch=19);");
    fprintf(rscript, "\n}");

    fprintf(rscript, "\n\n# Encerra a criacao do arquivo .png");
    fprintf(rscript, "\ndev.off()");

    fclose(rscript);
} // Terminou função escrevendoarquivo

int main(int argc, char **argv){

    /*
        Requisito 1 - Entrada de Dados: Verificando se o número de argumentos passado está correto
        Espera-se que sejam passados 3 argumentos: Nome do Programa, Arquivo Input e Arquivo Output
    */

    int k = 0, i = 0, j = 2;
    if(argc != 3){
        printf("\n Número errado de argumentos!\n");
        printf("\n ./trapezium.bin <inputFile> <outputFile>\n\n");
    }

    else{
        lendoarquivo(argv); // leitura do Arquivo de entrada.
        gausselimination(valuesofx, valuesofy); // Chama a função de Eliminação de Gauss para calcular o polinômio interpolador.
        
        k = contadorx - 1;

        for (i = 0; i < contadorx; i++) {   

        interpolador[i] = auxinterp[k];       // Inverte o vetor do polinomio interpolador, para ser jogado no método de Horner corretamente.
        k--;
        
        }
        
        for (i = 0; i < contadorp; i++) {
        
            pinterpolados[i] = metodohorner(interpolador, valuesofp[i]); // Calcula-se os pontos a interpolar utilizando Horner (Requisito 6).
        }
        
        valor_integral = calculointegral(limdown, limtop, number_trapezium);
        
        printf("Trapezium: Interpolador/Integrador Numerico \n\n");
        
        printf("Polinomio Interpolador: \n\n");

        printf("P(x) = %f", interpolador[contadorx-1]);
        printf(" + %f * x", interpolador[contadorx-2]);
            for (i=contadorx-3;i>=0;i--){
                printf(" + %f * x^%d",interpolador[i],j);
                j++;
            }
        
        printf("\n\nValores interpolados : \n\n");

        for(i=0;i<contadorp;i++){
            printf("P(%.2f) = %f \n",valuesofp[i],pinterpolados[i]);
        }
        
        printf("\n\nIntegral em [%.2f, %.2f] = %f ",limdown,limtop,valor_integral);

        escrevendoarquivo(argv);
            
        return 0;
    }
}
