
package MixtureModel;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.text.FieldPosition;
import java.io.PrintWriter;

public class Vector {

    public int n_dim;
    
    public double[] V;
    
    public Vector(int n_dim) {
        this.n_dim = n_dim;
        V = new double[n_dim];
    }                

    public Vector(double[] V, int n_dim) {
        this.n_dim = n_dim;
        this.V = V;
    }                

    public Vector(int n_dim, int s) {
        int i;

        this.n_dim = n_dim;
        V = new double[n_dim];
        for (i = 0;i < n_dim;i++)
            V[i] = s;
    }                
    
    void zero() {
        int i;
        
        for (i = 0;i < n_dim;i++)
            V[i] = 0;
    }

    Vector dup() {
        Vector nv = new Vector(n_dim);
        int i;
        
        for (i = 0;i < n_dim;i++)
            nv.V[i] = V[i];

        return nv;
    }

    void sub(Vector v) {
        int i;

        for (i = 0;i < n_dim;i++)
            V[i] -= v.V[i];
    }

    public void print(int w, int d) {
          print(new PrintWriter(System.out,true),w,d); 
    }

    public void print(PrintWriter output, int w, int d) {
        DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        print(output,format,w+2);
    }
    
    public void print (PrintWriter output, NumberFormat format, int width) {
        output.println();  // start on new line.
        for (int i = 0; i < n_dim; i++) {
            String s = format.format(V[i]); // format the number
            int padding = Math.max(1,width-s.length()); // At _least_ 1 space
            for (int k = 0; k < padding; k++)
                output.print(' ');
            output.print(s);
            output.println();
        }
        output.println();   // end with blank line.
    }

}