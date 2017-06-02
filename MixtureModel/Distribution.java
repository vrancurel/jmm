
package MixtureModel;

import java.io.PrintWriter;

public abstract class Distribution {

    abstract void DensityInit();
    abstract double DensityCompute(Range ranges[], double sample[]);

    abstract void UpdateParameters(MixtureModel mix, int kernel);

    abstract void RandParameters(Range ranges[], int k);
    
    //print parameters
    abstract void print(PrintWriter writer);

    abstract boolean isDrawable();
    abstract void draw(Gnuplot gnuplot) throws java.io.IOException;
}