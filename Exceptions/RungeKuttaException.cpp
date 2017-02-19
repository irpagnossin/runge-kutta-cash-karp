//------------------------------------------------------------------------------
// Arquivo: RungeKuttaException.cpp (versão 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.12.08
//------------------------------------------------------------------------------
#ifndef RUNGEKUTTAEXCEPTION_CPP
#define RUNGEKUTTAEXCEPTION_CPP

class RungeKuttaException{
public:
    RungeKuttaException( char * message ) : message( message ){}
    const char * getMessage( void ) const{ return( this->message ); }
private:
    const char * message;
};

#endif
//-------------- FIM DA DECLARAÇÃO E IMPLEMENTAÇÃO DA CLASSE RUNGEKUTTAEXCEPTION
