
package MixtureModel;

import java.io.PrintWriter;

public class MixtureModel {

    class Kernel {
        Distribution distrib;
        double weight;
        int n_samples;
        double probs[]; // n_samples

        void reset_probs() {
            for (int i = 0;i < n_samples;i++)
                probs[i] = 0;
        }

        Kernel(Distribution distrib, int n_samples) {
            this.distrib = distrib;
            this.n_samples = n_samples;
            probs = new double[n_samples];
            reset_probs();
        }
    }
                       
    int n_dim;
    Range ranges[];
    int n_samples;
    double data[][];
    int n_kernels;
    Kernel kernels[];
    double init_weights[];

    double probs[];

    void reset_probs() {
        for (int i = 0;i < n_samples;i++)
            probs[i] = 0;
    }

    double init_sum;
    
    public MixtureModel(int n_dim, 
                        Range ranges[],
                        int n_samples,
                        double data[][],
                        int n_distribs,
                        Distribution distribs[],
                        double init_weights[]) {
        
        this.n_dim = n_dim;
        this.ranges = ranges;
        this.n_samples = n_samples;
        this.data = data;
        n_kernels = n_distribs;
        this.init_weights = init_weights;

        kernels = new Kernel[n_kernels];

        init_sum = 0.0;
        for (int i = 0;i < n_distribs;i++) {
            kernels[i] = new Kernel(distribs[i], n_samples);
            init_sum += init_weights[i];
        }
        
        for (int i = 0;i < n_distribs;i++) {
            kernels[i].weight = init_weights[i] / init_sum;
        }
        
        probs = new double[n_samples];
    }

    /** 
     * random parameters
     * 
     */
    public void rand() {
        for (int i = 0;i < n_kernels;i++) {
            kernels[i].distrib.RandParameters(ranges, i);
        }
    }

    public void reset() {
        reset_probs();
        for (int i = 0;i < n_kernels;i++) {
            kernels[i].reset_probs();
            kernels[i].weight = init_weights[i] / init_sum;
        }
        
    }

    public void prob() {
        reset_probs();

        for (int i = 0;i < n_kernels;i++) {

            kernels[i].reset_probs();
            kernels[i].distrib.DensityInit();

            for (int j = 0;j < n_samples;j++) {
                
                double sample[] = new double[n_dim];
                for (int k = 0;k < n_dim;k++) {
                    sample[k] = data[k][j];
                    //System.out.println(sample[k]);
                }

                double result = kernels[i].distrib.DensityCompute(ranges, sample);
                
                //System.out.println(result);

                kernels[i].probs[j]= kernels[i].weight * result;

                //compute total probabilities for each sample
                probs[j] += kernels[i].probs[j];
            }
        }
    }

    public void iter() {

        prob();

        //adjust probabilities and recompute weight
        double sum = 0;
        for (int i = 0; i < n_kernels;i++) {

            kernels[i].weight = 0;
            for (int j = 0;j < n_samples;j++) {
                kernels[i].probs[j] /= probs[j];
                kernels[i].weight += kernels[i].probs[j];
            }

            kernels[i].weight /= n_samples;
            
            sum += init_weights[i];
        }

        //update parameters
        for (int i = 0;i < n_kernels;i++) {

            kernels[i].distrib.UpdateParameters(this, i);

            kernels[i].weight = 
                0.9 * kernels[i].weight +
                0.1 * init_weights[i] / sum;
        }
    }

    public double likelihood() {
        double tmp = 0;
        for (int i = 0;i < n_samples;i++)
            tmp += Math.log(probs[i]);
        
        return tmp / n_samples;
    }

    public void print(PrintWriter output) {
        int i;
        
        for (i = 0;i < n_kernels;i++) {
            output.println("kernel " + i + " weight " + kernels[i].weight);
            kernels[i].distrib.print(output);
        }
        
        output.println("likelihood " + likelihood());
    }

    public void redraw(Gnuplot gnuplot, String datfile) 
        throws java.io.IOException {
        int i;
        boolean first = true;
        
        if (!(n_dim == 2 ||
              n_dim == 3)) {
            throw new IllegalArgumentException("Unable to draw in such dimensions");
        }

        if (n_dim == 2)
            gnuplot.draw("plot \"" + datfile + "\""); 
        else if (n_dim == 3)
            gnuplot.draw("splot \"" + datfile + "\""); 

        for (i = 0;i < n_kernels;i++) {
            
            if (kernels[i].distrib.isDrawable()) {
                
                if (first)
                    gnuplot.draw(",");
                else
                    first = false;
                
                kernels[i].distrib.draw(gnuplot);
            }
        }
        gnuplot.draw("\n");
        gnuplot.flush();
    }
}