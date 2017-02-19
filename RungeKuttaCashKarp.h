//------------------------------------------------------------------------------
//   Arquivo: RungeKuttaCashKarp.h (versão 1.5)
//   Autor:   Ivan Ramos Pagnossin (irpagnossin@usp.br)
//   Data:    dezembro de 2006.
//
//   A classe RungeKuttaCashKarp resolve uma equação diferencial ordinária (EDO)
// de ordem qualquer através do consagrado método de Runge, Kutta e CashKarp  de
// 5ª ordem utilizando passos  adaptativos para atender ao erro solicitado. Veja
// o arquivo exemplo.c, fornecido com este, para outros detalhes.
//------------------------------------------------------------------------------
#ifndef RUNGEKUTTACashKarp_H
#define RUNGEKUTTACashKarp_H

class RungeKuttaCashKarp{
public:
     RungeKuttaCashKarp( const unsigned int, double (*)(double, double *) );
     RungeKuttaCashKarp( const unsigned int, double (*)(double, double *), const bool );
    ~RungeKuttaCashKarp( void );

    void setTimeout( double );
    void setError( double );

    double getTimeout( void ) const;
    double getError( void ) const;

    void solve( unsigned int, const double *, double **, double *, char * );

private:
    void checkOverflow( double );

    double error, eps, timeout, (*F)(double, double *);
    const int order;
    const bool isLinear;
    bool hasTimeout;
};

#endif
//------------------------------- FIM DA DECLARAÇÃO DA CLASSE RUNGEKUTTACashKarp

