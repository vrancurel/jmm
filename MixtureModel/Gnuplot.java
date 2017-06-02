
package MixtureModel;

import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.OutputStreamWriter;

public class Gnuplot {

    Process process;
    BufferedWriter output;
    //PrintWriter output;

    public Gnuplot()
        throws java.io.IOException {
        System.out.println("launching gnuplot");
        process = Runtime.getRuntime().exec("gnuplot");
        output = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
        //output = new PrintWriter(System.out,true);
    }

    public void draw(String str) 
        throws java.io.IOException {
        output.write(str);
    }

    public void flush()
        throws java.io.IOException {
        output.flush();
    }

}