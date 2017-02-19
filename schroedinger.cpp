/*

    - Ao alterar a constante de Planck de 7.61... (para distâncias em ângstrons) para 0.0761... (em nanometros),
    a função "schroedinger" passou a fornecer valores crescentes! Na primeira situação ela ficava sempre em torno
    da unidade. O que está acontecendo?

    - Na região classicamente proibida, y, y' e y'' [=F(x,y)] crescem exponencialmente, de modo que é importante escolher
    um fator multiplicativo (associado a 2m/hbar^2) pequeno.

 */

#include "RungeKuttaCashKarp.h"
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <fstream>

#include "TimeoutException.cpp"
#include "RungeKuttaException.cpp"

//#define SELF_ENERGY

using std::cout;


double energia = -0.0009111255;


double **psi;
//double **potential;

ofstream fileSchroedinger( "schroedinger.log" );

//------------------------------------------------------------------------------
// Equação de Schrödinger: y'' = 2m/hbar^2 * (V-E) * y.

// obs.: the value of <hbar2> was obtained by this way:
// hbar2 = (0.6582e-15 eV s)*(1.055e-34 J s)
// Now, exchange "J" for "kg m^2/s^2" and convert "kg" to "9.109e-31 kg" and "m" to "Å".
//------------------------------------------------------------------------------
double schroedinger(double x0, double *y0){

    const double mEff  = 0.067,
                 hbar2 = 7.6199682 * mEff; // É só dividir por 100 que pára de funcionar!

    //double V = 1/( 1 + exp((x0-40)/10) ) + 1/( 1 + exp((60-x0)/10) );
    double V = ( x0 > -8 && x0 < 8 ? -0.001 : 0 );
    fileSchroedinger << x0 << "\t" << y0[0] << "\t" << (2 * mEff / hbar2)*(V-energia)*y0[0] << endl;

    return( /*(2 * mEff / hbar2)*/(V-energia)*y0[0] );
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main( int argc, char *argv[] ){

    int nPtos = 10000, ordem = 2;

    double xI  = -20, // micro-metros
           xF  = +20, // micro-metros
           y0[2] = { 0, 1e-5 };

    double dx = (xF-xI)/(nPtos-1), next_x;
                                   
    //------------------------------------
    // Abre o arquivo para escrita.
    //------------------------------------
    FILE *arquivo;

    if( (arquivo = fopen( argv[1], "w+" )) == NULL ){
        printf( "I could not open the output file %s\n", argv[1] );
        return( -1 );
    }

    //------------------------------------
    // Monta os vetores *xE e *xD.
    //------------------------------------
    int ptosComuns = 3;
    if( nPtos < ptosComuns ){
        if( nPtos % 2 != 0 ) ptosComuns = nPtos;
        else ptosComuns = nPtos-1;
    }

    int meio = (nPtos-1) >> 1;

    ptosComuns = ptosComuns >> 1;
    if( ptosComuns == 0 ) ptosComuns = 1;

    const int nPtos_E =     1 + meio + ptosComuns,
              nPtos_D = nPtos - meio + ptosComuns;

    double *xE = new double[nPtos_E],
           *xD = new double[nPtos_D];

    for( register int i = 0; i <= nPtos-1; i++ ){

        next_x = xI + i*dx;

	    if( i <= meio+ptosComuns ) xE[i] = next_x;
	    if( i >= meio-ptosComuns ) xD[nPtos-1-i] = next_x;

    }

    //------------------------------------
    // Resolve a eq. de Schrödinger.
    //------------------------------------
    RungeKuttaCashKarp *EDO = new RungeKuttaCashKarp( ordem, schroedinger, true );
    EDO->setError( 1.0E-2 );
    EDO->setTimeout(10.0);
    /*
     * Define e aloca o estado inicial dos vetores-solução.
     */
    double **yE = new double * [ordem],
           **yD = new double * [ordem];

    for( register int i = 0; i <= ordem-1; i++ ){
        yE[i] = new double[nPtos];
        yD[i] = new double[nPtos];
    }

    #ifdef SELF_ENERGY
    const double E_i = 1e-3, E_f = 0.4;
    const int N = 100;

    for( register int iEnergy = 0; iEnergy <= N-1; iEnergy++ ){

        energia = E_i + iEnergy * ( E_f - E_i ) / ( N - 1 );
        cout << energia << "\n";
    #endif
        for( register int i = 0; i <= ordem-1; i++ )
            for( register int j = 0; j <= nPtos-1; j++ ){
                yE[i][j] = 0;
            }

        try{
            EDO->solve( nPtos_E, xE, yE, y0, "esquerda.log" );
            EDO->solve( nPtos_D, xD, yD, y0, "direita.log" );

            #ifdef SELF_ENERGY
            double avaliacao = fabs(yE[1][meio] / yE[0][meio] - yD[1][meio] / yD[0][meio] );
            cout << "avaliacao = " << avaliacao;


            /*
             * Constrói a função de onda.

            for( register int j = 0; j <= nPtos-1; j++ ){
                if( j < meio )
                    psi[j] = yE[];
                else
                    psi[j] = yD[];
            } */

            /*
             * Normaliza a função de onda.
             */
            //normalize( psi );

            fprintf( arquivo, "%f\t%f\n", energia, avaliacao );
            #endif
        }
        catch( std::bad_alloc ){
            cout << "not enough memory to run.\n";
            exit( -1 );
        }
        catch( TimeoutException e ){
            cout << e.getMessage() << "\n";
        }
        catch( RungeKuttaException e ){
            cout << e.getMessage() << "\n";
        }

    #ifdef SELF_ENERGY
    }
    #endif

    delete EDO;

    /*
     * Imprime a solução no arquivo.
     */
    #ifndef SELF_ENERGY
    for( register int i = 0; i <= nPtos-1; i++ ){

        if( i <= meio+ptosComuns )
            fprintf( arquivo, "%f\t%f\t", xE[i], yE[0][i] );

   	    if( i >= meio-ptosComuns ){
            if( i > meio+ptosComuns ){
                fprintf( arquivo, "%f\t\t", xD[nPtos-1-i] );
            }
            fprintf( arquivo, "%f", yD[0][nPtos-1-i] );
        }
        fprintf( arquivo, "\n" );
    }
    #endif

    fclose( arquivo );
    fileSchroedinger.close();

    return( 0 );
}


