
package MixtureModel;

import java.io.PrintWriter;

public class GaussianDistribution extends Distribution {

    int n_dim;

    Matrix cov;
    Vector mu;

    public GaussianDistribution(int n_dim) {
        this.n_dim = n_dim;
        mu = new Vector(n_dim);
        cov = new Matrix(n_dim, n_dim);
    }

    //temporary variables for density calculus
    double det;
    Matrix invcov;

    void DensityInit() {
        det = cov.det();
        
        //System.out.print("det = " + det + "\n");
        
        cov.print(3, 2);
        invcov = cov.inverse();

        //System.out.print("invcov ");
        //invcov.print(5, 10);
    }

    double DensityCompute(Range ranges[], double sample[]) {
        Vector csample, r1;
        double r2, result;
        int j, r;
        double k;

        k = n_dim;

        //centered sample
        csample = new Vector(sample, n_dim);
        csample.sub(mu); 

        //csample.print(3, 2);

        r1 = new Vector(n_dim, 0);

        // compute r1 = csample^T * invcov
        for (j = 0;j < csample.n_dim;j++)
            for (r = 0;r < csample.n_dim;r++)
                r1.V[j] += csample.V[r]*invcov.A[r][j];

        //r1.print(3, 2);

        //compute r2 = r1 * csample
        r2 = 0;
        for (r = 0;r < r1.n_dim;r++)
            r2 += r1.V[r]*csample.V[r];

        result = 1.0/(Math.pow(2.0*Math.PI,k/2.0)*Math.sqrt(det))*Math.exp((-1.0/2.0)*r2);
        
        //System.out.println(result);

        return result;
    }

    void UpdateParameters(MixtureModel mix, int kernel) {
        Vector tmpmu;
        Matrix tmpcov;
        int i, j, k;
        double tmp;
        
        tmpmu = new Vector(n_dim, 0);
        tmpcov = new Matrix(n_dim, n_dim, 0);

        for (i = 0;i < mix.n_samples;i++)
            for (j = 0;j < n_dim;j++)
                tmpmu.V[j] += mix.kernels[kernel].probs[i] * mix.data[j][i];

        tmp = mix.n_samples * mix.kernels[kernel].weight;

        for (i = 0;i < n_dim;i++)
            mu.V[i] = tmpmu.V[i]/tmp;

        //estimation of the covariance matrix
        for (i = 0;i < n_dim;i++)
            for (j = 0;j < n_dim;j++) {
                tmpcov.A[i][j] = 0;
                for (k = 0;k < mix.n_samples;k++) {
                    tmpcov.A[i][j] += 
                        mix.kernels[kernel].probs[k] *
                        (mix.data[i][k]-mu.V[i])*
                        (mix.data[j][k]-mu.V[j]);
                }
                cov.A[i][j] = tmpcov.A[i][j]/tmp;

                //variance limiting
                if (i == j && 
                    cov.A[i][j] < 0.000001)
                    cov.A[i][j] = 0.000001;
            }
    }

    void RandParameters(Range ranges[], int k) {
        
        cov.zero();
        
        for (int i = 0;i < n_dim;i++) {
            double length;

            length = Math.abs(ranges[i].end - ranges[i].start);
            mu.V[i] = ranges[i].start + length * Math.random();
            //mu.V[i] = ranges[i].start + length * 0.1 * k;
            cov.A[i][i] = (length/4 + 1.0) * (length/4 + 1.0);
        }
    }

    public void print(PrintWriter output) {

        output.println("mu:");
        mu.print(output, 3, 2);
        output.println("cov:");
        cov.print(output, 3, 2);
    }

    boolean isDrawable() {
        return true;
    }

    public void draw(Gnuplot gnuplot) 
        throws java.io.IOException {
        if (1 == n_dim) {
            //FIXME
        } else if (2 == n_dim) {
            String str = new String();

            String nstr = str.format("%f+cos(t)*%f,%f+sin(t)*%f", 
                                     mu.V[0], Math.sqrt(cov.A[0][0]),
                                     mu.V[1], Math.sqrt(cov.A[1][1]));
            gnuplot.draw(nstr);

        } else if (3 == n_dim) {
            String str = new String();

            //use Mercator parametrization
            //FIXME represent correlation matrix
            String nstr = str.format("%f+%f*sech(v)*cos(u),%f+%f*sech(v)*sin(u),%f+%f*tanh(v)", 
                                     mu.V[0], Math.sqrt(cov.A[0][0]),
                                     mu.V[1], Math.sqrt(cov.A[1][1]),
                                     mu.V[2], Math.sqrt(cov.A[2][2]));
            gnuplot.draw(nstr);

        } else {
            System.out.println("not drawable\n");
        }
    }
    
}
