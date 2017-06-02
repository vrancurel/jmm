package MixtureModel;

public class Utils {

   /** 
    * sqrt(a^2 + b^2) without under/overflow. 
    */
    public static double hypot(double a, double b) {
        double r;
        if (Math.abs(a) > Math.abs(b)) {
            r = b/a;
            r = Math.abs(a)*Math.sqrt(1+r*r);
        } else if (b != 0) {
            r = a/b;
            r = Math.abs(b)*Math.sqrt(1+r*r);
        } else {
            r = 0.0;
        }
        return r;
    }

    /** 
     * Check correlation matrix
     * 
     * @param rho 
     * 
     * @return true if matrix is OK
     */
    boolean rho_check(Matrix rho) {
        int i, j;

        if (rho.m != rho.n)
            return false;
        
        for (i = 0; i < rho.m;i++)
            if (rho.A[i][i] != 1)
                return false;

        for (i = 0;i < rho.m;i++)
            for (j = 0;j < i;j++)
                if (rho.A[i][j] != rho.A[j][i])
                    return false;
        
        return true;
    }
    
    /** 
     * Compute the covariance matrix
     * 
     * @param n_dim 
     * @param sigma deviation vector
     * @param rho correlation matrix (not checked)
     * 
     * @return the covariance matrix
     */
    Matrix comp_cov(int n_dim, Vector sigma, Matrix rho) {
        Matrix cov = new Matrix(n_dim, n_dim);
        
        int i, j;

        for (i = 0;i < n_dim;i++)
            for (j = 0;j < n_dim;j++)
                cov.A[i][j] = rho.A[i][j] * sigma.V[i] * sigma.V[j];
        
        return cov;
    }

    
    /** 
     * Peter J. Acklam, compute the inverse normal cdf
     * 
     * @param p probability
     * 
     * @return z, x random and P(x<=z) = p
     */
    double invnorm(double p)
    {
        /* Coefficients in rational approximations. */
        double [] a =
            {
                -3.969683028665376e+01,
                2.209460984245205e+02,
                -2.759285104469687e+02,
                1.383577518672690e+02,
                -3.066479806614716e+01,
                2.506628277459239e+00
            };
        
        double [] b =
            {
                -5.447609879822406e+01,
                1.615858368580409e+02,
                -1.556989798598866e+02,
                6.680131188771972e+01,
                -1.328068155288572e+01
            };
        
        double [] c =
            {
                -7.784894002430293e-03,
                -3.223964580411365e-01,
                -2.400758277161838e+00,
                -2.549732539343734e+00,
                4.374664141464968e+00,
                2.938163982698783e+00
            };
        
        double [] d =
            {
                7.784695709041462e-03,
                3.224671290700398e-01,
                2.445134137142996e+00,
                3.754408661907416e+00
            };
        
        final double LOW = 0.02425;
        final double HIGH = 0.97575;
        
        double q, r;
        
        if (p < 0 || p > 1)
            throw new IllegalArgumentException("Domain error");
        else if (p == 0)
            throw new IllegalArgumentException("Range error");
        else if (p == 1)
            throw new IllegalArgumentException("Range error");
        else if (p < LOW) {
            /* Rational approximation for lower region */
            q = Math.sqrt(-2*Math.log(p));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        } else if (p > HIGH) {
            /* Rational approximation for upper region */
            q  = Math.sqrt(-2*Math.log(1-p));
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        } else {
            /* Rational approximation for central region */
            q = p - 0.5;
            r = q*q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
        }
    }

    /** 
     * Vector-matrix multiplication
     *
     * @note b is treated as a row vector
     * 
     * @param a 
     * @param b 
     * 
     * @return 
     */
    Vector vm_mlt(Matrix a, Vector b) {
        if (a.m != b.n_dim)
            throw new IllegalArgumentException("Matrix and vector does not have the same dimension");

        Vector out = new Vector(b.n_dim);

        int i, j;

        for (j = 0; j < a.m;j++)
            if (b.V[j] != 0.0) {
                for (i = 0;i < a.n;i++)
                    out.V[i] += b.V[j]*a.A[j][i];
                }

        for (j = 0;j < a.n;j++) {
            double sum = 0.0;
            for (i = 0;i < a.m;i++)
                sum += b.V[i]*a.A[i][j];
            out.V[j] = sum;
        }

        return out;
    }

    /** 
     * Multivariate Normal Sampling
     * 
     * @param n_dim
     * @param mu mean vector
     * @param sigma deviation vector
     * @param rho correlation matrix (all values must be provided)
     * 
     * @return a sample
     */
    public Vector mvnsamp(int n_dim, Vector mu, Vector sigma, Matrix rho) {
        if (!rho_check(rho))
            throw new IllegalArgumentException("Illegal correlation matrix");

        Matrix cov = comp_cov(n_dim, sigma, rho);

        //System.out.println("cov");
        //cov.print(3,2);

        //cov.chol().CholeskyDecompositionInSitu(cov);
        Matrix chol = cov.chol().getL2();
        
        //System.out.println("chol");
        //chol.print(3,2);

        Vector z = new Vector(n_dim);

        for (int i = 0;i < n_dim;i++)
            z.V[i] = invnorm(Math.random());
        
        Vector x = vm_mlt(chol, z);

        for (int i = 0; i < n_dim;i++)
            x.V[i] += mu.V[i];
        
        return x;
    }

    /** 
     * sample bivariata correlated random normal variables
     * 
     * @param mu1 
     * @param sigma1 
     * @param mu2 
     * @param sigma2 
     * @param rho 
     * 
     * @return 
     */
    public Vector bvnsamp(double mu1, double sigma1,
                          double mu2, double sigma2, double rho) {
        double z1, z2;
        
        z1 = invnorm(Math.random());
        z2 = invnorm(Math.random());

        Vector v = new Vector(2);
        v.V[0] = mu1 + sigma1*z1;
        v.V[1] = mu2 + sigma2*(z1*rho+z2*Math.sqrt(1-rho*rho));
        
        return v;
    }
    
}
