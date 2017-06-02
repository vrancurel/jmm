
package MixtureModel;

import java.io.PrintWriter;

public class UniformDistribution extends Distribution {

    int n_dim;

    public UniformDistribution(int n_dim) {
        this.n_dim = n_dim;
    }

    void DensityInit() {
    }

    double DensityCompute(Range ranges[], double sample[]) {
        double prod, result;
        int i;

        prod = 1;
        for (i = 0;i < n_dim;i++)
            prod *= ranges[i].end - ranges[i].start;

        result = 1/prod;

        return result;
    }

    void UpdateParameters(MixtureModel mix, int kernel) {
    }

    void RandParameters(Range ranges[], int k) {
    }

    public void print(PrintWriter output) {
    }

    boolean isDrawable() {
        return false;
    }

    public void draw(Gnuplot gnuplot) {
    }
    
}
