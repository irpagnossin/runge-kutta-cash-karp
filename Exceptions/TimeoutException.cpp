//------------------------------------------------------------------------------
// Arquivo: TimeoutException.cpp (vers�o 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.12.08
//------------------------------------------------------------------------------
#ifndef TIMEOUTEXCEPTION_CPP
#define TIMEOUTEXCEPTION_CPP

#include "RungeKuttaException.cpp"

class TimeoutException : public RungeKuttaException{
public:
    TimeoutException() : RungeKuttaException( "TIMEOUT.\n" ) {}
};

#endif
//----------------- FIM DA DECLARA��O E IMPLEMENTA��O DA CLASSE TIMEOUTEXCEPTION
 