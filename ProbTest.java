
import MixtureModel.*;
import java.io.*;
import java.util.Scanner;
import java.util.Locale;

/*
 * javac ProbTest.java
 * java -cp . ProbTest
 */

public class ProbTest {

    static void test_mat_inv()
    {
        double [][] A1 = 
            {{-2, 2, -3},
             {-1, 1, 3},
             {2, 0, 1}};
        Matrix mat = new Matrix(A1, 3, 3);

        System.out.println("det = " + mat.det());

        mat.print(3, 2);

        Matrix inv = mat.inverse();

        inv.print(3, 2);
    }

    static void test_chol_decomp()
    {
        final int N_DIM = 3;
        PrintWriter output = new PrintWriter(System.out,true);

        double [][] A1 = 
            {{1,    0.6,  0.30},
             {0.6,  1,    0.5},
             {0.30, 0.5,  1}};
        Matrix mat = new Matrix(A1, N_DIM, N_DIM);
        
        mat.print(output, 3, 2);
        Matrix mat2 = mat.chol().getL2();
        mat2.print(output, 3, 2);
    }

    static void test_mvmix(boolean reload)
        throws java.io.IOException {
        MixtureModel mix;
        Gnuplot gnuplot;
        Distribution distribs[];
        double init_weights[];
        double V[];
        double A[][];
        Scanner input = new Scanner(System.in);
        File tmpfile;
        PrintWriter output = new PrintWriter(System.out,true);
        PrintWriter foutput;
        FileReader fin;
        Scanner finput;
        Range ranges[];
        final int N_SAMPLES = 300;
        final int LENGTH = 500;
        final int N_DIM = 2;
        final int N_KERNELS = 4;
        String DATFILE = new String("/tmp/mvmix.dat");

        distribs = new Distribution[N_KERNELS];
        distribs[0] = new UniformDistribution(N_DIM);
        distribs[1] = new GaussianDistribution(N_DIM);
        distribs[2] = new GaussianDistribution(N_DIM);
        distribs[3] = new GaussianDistribution(N_DIM);

        init_weights = new double[N_KERNELS];
        init_weights[0] = 0.1;
        init_weights[1] = 1;
        init_weights[2] = 1;
        init_weights[3] = 1;


        gnuplot = new Gnuplot();
        gnuplot.draw("set parametric\n");

        if (!reload) {
            finput = null;
            tmpfile = new File(DATFILE);
            tmpfile.createNewFile();
            foutput = new PrintWriter(new FileWriter(tmpfile));
        } else {
            foutput = null;
            fin = new FileReader(DATFILE);
            finput = new Scanner(fin);
            finput.useLocale(Locale.US);
        }
        
        ranges = new Range[N_DIM];
        ranges[0] = new Range(0, LENGTH);
        ranges[1] = new Range(0, LENGTH);

        double [][] data = new double[N_DIM][N_SAMPLES];

        Utils utils = new Utils();

        for (int i = 0;i < N_SAMPLES/3;i++) {

            if (reload) {
                if (!finput.hasNext())
                    break ;
            }

            Vector x1;
            if (reload) {
                x1 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x1.V[0] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x1.V[2] = Double.parseDouble(finput.next());
            } else {
                x1 = utils.bvnsamp(100, 20, 100, 20, 0.34);
            }
            data[0][3*i] = x1.V[0];
            data[1][3*i] = x1.V[1];

            Vector x2;
            if (reload) {
                x2 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x2.V[0] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x2.V[2] = Double.parseDouble(finput.next());
            } else {
                x2 = utils.bvnsamp(200, 50, 300, 10, -0.21);
            }
            data[0][3*i+1] = x2.V[0];
            data[1][3*i+1] = x2.V[1];

            Vector x3;
            if (reload) {
                x3 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x3.V[0] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x3.V[2] = Double.parseDouble(finput.next());
            } else {
                x3 = utils.bvnsamp(400, 15, 400, 15, 0.15);
            }
            data[0][3*i+2] = x3.V[0];
            data[1][3*i+2] = x3.V[1];

            if (!reload) {
                foutput.println(data[0][3*i] + " " + data[1][3*i]);
                foutput.println(data[0][3*i+1] + " " + data[1][3*i+1]);
                foutput.println(data[0][3*i+2] + " " + data[1][3*i+2]);
            }
        }

        if (!reload)
            foutput.close();

        //gnuplot.draw("plot \"" + DATFILE + "\"\n");

        mix = new MixtureModel(N_DIM, ranges, N_SAMPLES, data, N_KERNELS, distribs, init_weights);

        mix.rand();
        mix.print(output);
        
        mix.redraw(gnuplot, DATFILE);
        input.nextLine();        

        for (int i = 0;i < 100;i++) {
            mix.iter();
            mix.redraw(gnuplot, DATFILE);
            mix.print(output);
            input.nextLine();
        }

        try {
            Thread.sleep(50000);
        } catch (java.lang.InterruptedException e) {
        }
    }
 
    static void test_mvmix2(boolean reload) 
        throws java.io.IOException {
        MixtureModel mix;
        Gnuplot gnuplot;
        Distribution distribs[];
        double init_weights[];
        Vector mu1, sigma1, mu2, sigma2, mu3, sigma3, x;
        Matrix rho1, rho2, rho3;
        double V[];
        double A[][];
        Scanner input = new Scanner(System.in);
        File tmpfile;
        PrintWriter output = new PrintWriter(System.out,true);
        PrintWriter foutput;
        FileReader fin;
        Scanner finput;
        Range ranges[];
        final int N_SAMPLES = 900;
        final int LENGTH = 500;
        final int N_DIM = 3;
        final int N_KERNELS = 4;
        String DATFILE = new String("/tmp/mvmix1bis.dat");

        distribs = new Distribution[N_KERNELS];
        distribs[0] = new UniformDistribution(N_DIM);
        distribs[1] = new GaussianDistribution(N_DIM);
        distribs[2] = new GaussianDistribution(N_DIM);
        distribs[3] = new GaussianDistribution(N_DIM);

        init_weights = new double[N_KERNELS];
        init_weights[0] = 0.1;
        init_weights[1] = 1;
        init_weights[2] = 1;
        init_weights[3] = 1;

        double [] V1a = {100, 90, 120};
        mu1 = new Vector(V1a, N_DIM);
        double [] V1b = {30, 20, 70};
        sigma1 = new Vector(V1b, N_DIM);
        double [][] A1 = 
            {{1,    0.26, 0.31},
             {0.26, 1,    0.5},
             {0.31, 0.5,  1}};
        rho1 = new Matrix(A1, N_DIM, N_DIM);

        //mu1.print(output, 3, 2);
        //sigma1.print(output, 3, 2);
        //rho1.print(output, 3, 2);

        double [] V2a = {300, 50, 320};
        mu2 = new Vector(V2a, N_DIM);
        double [] V2b = {10, 10, 25};
        sigma2 = new Vector(V2b, N_DIM);
        double [][] A2 = 
            {{1,    0.25, 0.19},
             {0.25, 1,    0.42},
             {0.19, 0.42, 1}};
        rho2 = new Matrix(A2, N_DIM, N_DIM);

        double [] V3a = {320, 450, 350};
        mu3 = new Vector(V3a, N_DIM);
        double [] V3b = {70, 20, 55};
        sigma3 = new Vector(V3b, N_DIM);
        double [][] A3 = 
            {{1,   0,   0},
             {0,   1,   0}, 
             {0,   0,   1}};
        rho3 = new Matrix(A3, N_DIM, N_DIM);

        gnuplot = new Gnuplot();
        gnuplot.draw("set parametric\n");
        gnuplot.draw("set isosamples 30\n");
        gnuplot.draw("sech(x) = 1/cosh(x)\n");

        if (!reload) {
            finput = null;
            tmpfile = new File(DATFILE);
            tmpfile.createNewFile();
            foutput = new PrintWriter(new FileWriter(tmpfile));
        } else {
            foutput = null;
            fin = new FileReader(DATFILE);
            finput = new Scanner(fin);
            finput.useLocale(Locale.US);
        }
        
        ranges = new Range[N_DIM];
        ranges[0] = new Range(0, LENGTH);
        ranges[1] = new Range(0, LENGTH);
        ranges[2] = new Range(0, LENGTH);

        double [][] data = new double[N_DIM][N_SAMPLES];

        Utils utils = new Utils();

        for (int i = 0;i < N_SAMPLES/3;i++) {

            if (reload) {
                if (!finput.hasNext())
                    break ;
            }

            Vector x1;
            if (reload) {
                x1 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x1.V[0] = Double.parseDouble(finput.next());
                x1.V[1] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x1.V[2] = Double.parseDouble(finput.next());
            } else {
                x1 = utils.mvnsamp(N_DIM, mu1, sigma1, rho1);
            }
            data[0][3*i] = x1.V[0];
            data[1][3*i] = x1.V[1];
            data[2][3*i] = x1.V[2];

            Vector x2;
            if (reload) {
                x2 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x2.V[0] = Double.parseDouble(finput.next());
                x2.V[1] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x2.V[2] = Double.parseDouble(finput.next());
            } else {
                x2 = utils.mvnsamp(N_DIM, mu2, sigma2, rho2);
            }
            data[0][3*i+1] = x2.V[0];
            data[1][3*i+1] = x2.V[1];
            data[2][3*i+1] = x2.V[2];

            Vector x3;
            if (reload) {
                x3 = new Vector(N_DIM);
                finput.useDelimiter(" ");
                x3.V[0] = Double.parseDouble(finput.next());
                x3.V[1] = Double.parseDouble(finput.next());
                finput.useDelimiter("\n");
                x3.V[2] = Double.parseDouble(finput.next());
            } else {
                x3 = utils.mvnsamp(N_DIM, mu3, sigma3, rho3);
            }
            data[0][3*i+2] = x3.V[0];
            data[1][3*i+2] = x3.V[1];
            data[2][3*i+2] = x3.V[2];

            if (!reload) {
                foutput.println(data[0][3*i] + " " + data[1][3*i] + " " + data[2][3*i]);
                foutput.println(data[0][3*i+1] + " " + data[1][3*i+1] + " " + data[2][3*i+1]);
                foutput.println(data[0][3*i+2] + " " + data[1][3*i+2] + " " + data[2][3*i+2]);
            }
        }

        if (!reload)
            foutput.close();

        //gnuplot.draw("splot \"" + DATFILE + "\"\n");

        mix = new MixtureModel(N_DIM, ranges, N_SAMPLES, data, N_KERNELS, distribs, init_weights);

        mix.rand();
        mix.print(output);
        
        mix.redraw(gnuplot, DATFILE);
        input.nextLine();        

        for (int i = 0;i < 100;i++) {
            mix.iter();
            mix.redraw(gnuplot, DATFILE);
            mix.print(output);
            input.nextLine();
        }

        try {
            Thread.sleep(50000);
        } catch (java.lang.InterruptedException e) {
        }
    }
    
    public static void main(String[] args) throws IOException {

        if (args.length != 1)
            throw new IllegalArgumentException("Usage: ProbTest nb");

        int test = Integer.parseInt(args[0]);
        
        if (test == 8)
            test_chol_decomp();
        else if (test == 16)
            test_mat_inv();
        else if (test == 20)
            test_mvmix(false);
        else if (test == 21)
            test_mvmix2(false);
        else if (test == 210)
            test_mvmix2(true);
        else
            throw new IllegalArgumentException("No such test nb");
        return  ;
    }
}
 

