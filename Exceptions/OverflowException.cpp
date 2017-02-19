//------------------------------------------------------------------------------
// Arquivo: OverflowException.cpp (versão 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.12.08
//------------------------------------------------------------------------------
#ifndef OVERFLOWEXCEPTION_CPP
#define OVERFLOWEXCEPTION_CPP

#include "RungeKuttaException.cpp"

class OverflowException : public RungeKuttaException{
public:
    OverflowException() : RungeKuttaException( "OVERFLOW" ) {}
};

#endif
//--------------- FIM DA DECLARAÇÃO E IMPLEMENTAÇÃO DA CLASSE OVERFLOWTEXCEPTION
 