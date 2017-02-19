package cepa.edu.math.ode;

/**
 *
 * @author ivan.pagnossin
 */
public abstract class ODE {

    private int order = 0;
    
    public ODE( int order ){
        setOrder(order);
    }
    public void setOrder( int order ){
        if ( order >= 0 ) this.order = order;
    }
    public int getOrder(){
        return order;
    }
    
    /*
     * Funcional que define a equação diferencial ordinária. Por exemplo,
     * y''' = F(x,y,y',y'')
     * 
     * x é a variável independente e y[j] a j-ésima derivada de y.
     * O objetivo da EDO é determinar y[n], sendo n a ordem da EDO.
     */
    public abstract double F( double x, double[] y );
}
