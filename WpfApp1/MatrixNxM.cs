using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace WpfApp1
{
    public class MatrixNxM
    {
        double[,] data;
        public int rows;
        public int cols;

        public MatrixNxM(int rows, int cols)
        {
            this.rows = rows;
            this.cols = cols;
            this.data = new double[rows, cols];
        }

        public MatrixNxM(double[,] data)
        {
            this.rows = data.GetLength(0);
            this.cols = data.GetLength(1);
            this.data = new double[rows, cols];

            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    this.data[i, j] = data[i, j];
        }

        public void ZeroMatrix()
        {
            for( int i = 0 ; i < rows ; ++i )
                for( int j = 0 ; j < cols ; ++j )
                    data[i,j] = 0;
        }

        public void Identity()
        {
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                {
                    if (i == j)
                        data[i, j] = 1.0;
                    else
                        data[i, j] = 0.0;
                }
        }

        public MatrixNxM CopyMatrix()
        {
            MatrixNxM c = new MatrixNxM(rows, cols);

            for( int i = 0 ; i < rows ; ++i )
                for( int j = 0 ; j < cols ; ++j )
                    c[i,j] = data[i,j];

            return c;
        }

        public double this[int i, int j]
        {
            get{ return data[i,j]; }
            set{ data[i,j] = value; }
        }

        public double Determinant()
        {
            Debug.Assert( this.rows == this.cols );

            int N = this.rows;

            double[] a = new double[N * N];
            int[] pivot = new int [N];

            Func<int, int, int> IND = (i, j) => i + j * N;

            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    a[IND(i,j)] = this[i, j];
                }
            }

            int info = BurkardtLib.dge_fa(N, a, pivot);

            if (info != 0)
                throw new Exception(string.Format("The factorization failed on the {0} step", info));

            return BurkardtLib.dge_det(N, a, pivot);
        }

        public static MatrixNxM operator *(MatrixNxM a, MatrixNxM b) 
        {
            Debug.Assert(a.cols == b.rows);

            MatrixNxM c = new MatrixNxM(a.rows, b.cols);

	        c.ZeroMatrix();

            for (int i = 0; i < a.rows; ++i)
                for (int j = 0; j < b.cols; ++j)
				    for( int k = 0 ; k < a.cols ; ++k )
					    c[i,j] += a[i,k] * b[k,j];

		    return c;
	    }

        public static MatrixNxM operator *(MatrixNxM a, double b) 
        {
            MatrixNxM c = new MatrixNxM(a.rows, a.cols);

            for ( int i = 0 ; i < a.rows ; ++i)
			    for ( int j = 0 ; j < a.cols ; ++j )
				    c[i,j] = a[i,j] * b;

		    return c;
	    }

        public static MatrixNxM operator *(double b, MatrixNxM a) 
        {
            return a * b;
	    }

        public static MatrixNxM operator +(MatrixNxM a, MatrixNxM b) 
        {
            Debug.Assert(a.rows == b.rows && a.cols == b.cols);

            MatrixNxM c = new MatrixNxM(a.rows, a.cols);

            for (int i = 0; i < a.rows; ++i)
                for (int j = 0; j < a.cols; ++j)
                    c[i, j] = a[i, j] + b[i, j];

            return c;
        }

        public MatrixNxM Transpose()
        {
            return MatrixNxM.Transpose(this);
        }

        public static MatrixNxM Transpose(MatrixNxM a)
        {
            MatrixNxM c = new MatrixNxM(a.cols, a.rows);

            for (int i = 0; i < a.rows; ++i)
                for (int j = 0; j < a.cols; ++j)
                    c[j, i] = a[i, j];
            return c;
        }

        public void SVD(out MatrixNxM U, out MatrixNxM W, out MatrixNxM V)
        {
            double[,] u = new double[rows + 1, cols + 1];
            double[,] v = new double[cols + 1, cols + 1];
            double[] w = new double[cols + 1];

            for (int i = 1; i <= this.rows; ++i)
                for (int j = 1; j <= this.cols; ++j)
                    u[i, j] = this.data[i-1, j-1];

            NumericalRecipesLib.svdcmp(u, rows, cols, w, v);

            U = new MatrixNxM(rows, cols);
            V = new MatrixNxM(cols, cols);
            W = new MatrixNxM(cols, cols);

            for (int i = 1; i <= this.rows; ++i)
                for (int j = 1; j <= this.cols; ++j)
                    U[i-1, j-1] = u[i, j];

            for (int i = 1; i <= this.cols; ++i)
                for (int j = 1; j <= this.cols; ++j)
                    V[i - 1, j - 1] = v[i, j];

            W.ZeroMatrix();

            for (int i = 1; i <= this.cols; ++i)
                W[i - 1, i - 1] = w[i];
        }

        public MatrixNxM SolveSVD(MatrixNxM rhs)
        {
            Debug.Assert(rhs.rows == this.rows);
            Debug.Assert(this.rows == this.cols);

            double[,] u = new double[rows + 1, cols + 1];
            double[,] v = new double[cols + 1, cols + 1];
            double[] w = new double[cols + 1];
            double[] b = new double[rows + 1];
            double[] x = new double[cols + 1];

            const double illConditionedThreshold = 1e-6;

            for (int i = 1; i <= this.rows; ++i)
                for (int j = 1; j <= this.cols; ++j)
                    u[i, j] = this.data[i - 1, j - 1];

            NumericalRecipesLib.svdcmp(u, rows, cols, w, v);

            for (int i = 1; i <= cols; ++i)
                if (w[i] < illConditionedThreshold)
                    w[i] = 0.0;

            MatrixNxM res = new MatrixNxM(rhs.rows, rhs.cols);

            for (int j = 0; j < rhs.cols; ++j)
            {
                for (int i = 0; i < rhs.rows; ++i)
                {
                    b[i + 1] = rhs[i, j];
                }

                NumericalRecipesLib.svbksb(u, w, v, rows, cols, b, x);

                for (int i = 0; i < res.rows; ++i)
                {
                    res[i, j] = x[i + 1];
                }
            }

            return res;
        }

        public MatrixNxM Solve(MatrixNxM rhs)
        {
            Debug.Assert(rhs.rows == this.rows);
            Debug.Assert(this.rows == this.cols);

            int N = this.rows;
            int RHS_NUM = rhs.cols;

            double[] a = new double[N * (N + RHS_NUM)];
            Func<int, int, int> IND = (i, j) => i + j * N;

            for (int i = 0; i < this.rows; ++i)
            {
                for (int j = 0; j < this.cols; ++j)
                {
                    a[IND(i,j)] = this[i, j];
                }
                for (int j = 0; j < rhs.cols; ++j)
                {
                    a[IND(i, this.cols + j)] = rhs[i, j];
                }
            }

            int solution = BurkardtLib.dmat_solve(N, RHS_NUM, a);

            if (solution != 0)
                throw new Exception(string.Format("factorization failed on step {0}, and the solutions could not be computed.", solution));

            MatrixNxM res = new MatrixNxM(rhs.rows, rhs.cols);

            for (int i = 0; i < res.rows; ++i)
            {
                for (int j = 0; j < rhs.cols; ++j)
                {
                    res[i, j] = a[IND(i, this.cols + j)];
                }
            }

            return res;
        }

        public MatrixNxM SolveLU(MatrixNxM rhs)
        {
            Debug.Assert(rhs.rows == this.rows);
            Debug.Assert(this.rows == this.cols);

            int N = this.rows;

            double[] a = new double[N * N];
            double[] b = new double[N];
            int[] pivot = new int[N];

            Func<int, int, int> IND = (i, j) => i + j * N;

            for (int i = 0; i < this.rows; ++i)
            {
                for (int j = 0; j < this.cols; ++j)
                {
                    a[IND(i, j)] = this[i, j];
                }
            }

            int info = BurkardtLib.dge_fa(N, a, pivot);

            if (info != 0)
                throw new Exception(string.Format("The factorization failed on the {0} step", info));

            MatrixNxM res = new MatrixNxM(rhs.rows, rhs.cols);

            for (int j = 0; j < rhs.cols; ++j)
            {
                for (int i = 0; i < rhs.rows; ++i)
                {
                    b[i] = rhs[i, j];
                }

                BurkardtLib.dge_sl(N, a, pivot, b, 0);

                for (int i = 0; i < res.rows; ++i)
                {
                    res[i, j] = b[i];
                }
            }

            return res;
        }
    }

    public static class NumericalRecipesLib
    {
        //Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
        //v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
        //square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
        //No input quantities are destroyed, so the routine may be called sequentially with different b’s.
        public static void svbksb(double[,] u, double[] w, double[,] v, int m, int n, double[] b, double[] x)
        {
            Debug.Assert(u.GetLength(0) == m + 1);
            Debug.Assert(u.GetLength(1) == n + 1);
            Debug.Assert(v.GetLength(0) == n + 1);
            Debug.Assert(v.GetLength(1) == n + 1);
            Debug.Assert(w.Length == n + 1);
            Debug.Assert(b.Length == m + 1);
            Debug.Assert(x.Length == n + 1);

            int jj,j,i;
            double s;
            double[] tmp;

            tmp = new double[n + 1];
            for (j=1;j<=n;j++) {
                s=0.0;
                if (w[j] != 0.0) {
                    for (i = 1; i <= m; i++) s += u[i, j] * b[i];
                    s /= w[j];
                }
                tmp[j]=s;
            }
            for (j=1;j<=n;j++) {
                s=0.0;
                for (jj = 1; jj <= n; jj++) s += v[j, jj] * tmp[jj];
                x[j]=s;
            }
        }

        //Given a matrix A[1..m][1..n], this routine computes its singular value decomposition, A =
        //U·W·V^T. Thematrix U replaces A on output. The diagonal matrix of singular values W is output
        //as a vector w[1..n]. Thematrix V (not the transpose V^T ) is output as v[1..n][1..n].
        public static void svdcmp(double[,] a, int m, int n, double[] w, double[,] v)
        {
            Debug.Assert(a.GetLength(0) == m + 1);
            Debug.Assert(a.GetLength(1) == n + 1);
            Debug.Assert(v.GetLength(0) == n + 1);
            Debug.Assert(v.GetLength(1) == n + 1);
            Debug.Assert(w.Length == n + 1);

            //Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
            Func<double, double, double> PYTHAG = (valA, valB) =>
            {
                Func<double,double> SQR = X => X*X;

                double absa, absb;
                absa = Math.Abs(valA);
                absb = Math.Abs(valB);
                if (absa > absb) return absa*Math.Sqrt(1.0+SQR(absb/absa));
                else return (absb == 0.0 ? 0.0 : absb*Math.Sqrt(1.0+SQR(absa/absb)));
            };

            Func<double, double, double> SIGN = (valA, valB) => ((valB) >= 0.0 ? Math.Abs(valA) : -Math.Abs(valA));

            int flag,i,its,j,jj,k,l=0,nm=0;
            double c,f,h,s,x,y,z;
            double anorm=0.0, g=0.0, scale=0.0;
            double[] rv1;

            if (m < n) throw new Exception("SVDCMP: You must augment A with extra zero rows");
            rv1 = new double[n + 1];
            for (i=1;i<=n;i++)
            {
                l=i+1;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m) {
                    for (k = i; k <= m; k++) scale += Math.Abs(a[k, i]);
                    if (scale != 0.0) {
                        for (k=i;k<=m;k++) {
                            a[k,i] /= scale;
                            s += a[k,i]*a[k,i];
                        }
                        f=a[i,i];
                        g = -SIGN(Math.Sqrt(s),f);
                        h=f*g-s;
                        a[i,i]=f-g;
                        if (i != n) {
                            for (j=l;j<=n;j++) {
                                for (s=0.0,k=i;k<=m;k++) s += a[k,i]*a[k,j];
                                f=s/h;
                                for (k=i;k<=m;k++) a[k,j] += f*a[k,i];
                            }
                        }
                        for (k=i;k<=m;k++) a[k,i] *= scale;
                    }
                }
                w[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m && i != n) {
                    for (k = l; k <= n; k++) scale += Math.Abs(a[i,k]);
                    if (scale != 0.0) {
                        for (k=l;k<=n;k++) {
                            a[i,k] /= scale;
                            s += a[i,k]*a[i,k];
                        }
                        f=a[i,l];
                        g = -SIGN(Math.Sqrt(s), f);
                        h=f*g-s;
                        a[i,l]=f-g;
                        for (k=l;k<=n;k++) rv1[k]=a[i,k]/h;
                        if (i != m) {
                            for (j=l;j<=m;j++) {
                                for (s=0.0,k=l;k<=n;k++) s += a[j,k]*a[i,k];
                                for (k=l;k<=n;k++) a[j,k] += s*rv1[k];
                            }
                        }
                        for (k=l;k<=n;k++) a[i,k] *= scale;
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(w[i]) + Math.Abs(rv1[i])));
            }
            for (i=n;i>=1;i--) {
                if (i < n) {
                    if (g != 0.0) {
                        for (j=l;j<=n;j++)
                            v[j, i] = (a[i, j] / a[i, l]) / g;
                        for (j=l;j<=n;j++) {
                            for (s = 0.0, k = l; k <= n; k++) s += a[i, k] * v[k, j];
                            for (k = l; k <= n; k++) v[k, j] += s * v[k, i];
                        }
                    }
                    for (j = l; j <= n; j++) v[i, j] = v[j, i] = 0.0;
                }
                v[i, i] = 1.0;
                g=rv1[i];
                l=i;
            }
            for (i=n;i>=1;i--) {
                l=i+1;
                g=w[i];
                if (i < n)
                    for (j = l; j <= n; j++) a[i, j] = 0.0;
                if (g != 0.0) {
                    g=1.0/g;
                    if (i != n) {
                        for (j=l;j<=n;j++) {
                            for (s = 0.0, k = l; k <= m; k++) s += a[k, i] * a[k, j];
                            f = (s / a[i, i]) * g;
                            for (k = i; k <= m; k++) a[k, j] += f * a[k, i];
                        }
                    }
                    for (j = i; j <= m; j++) a[j, i] *= g;
                } else {
                    for (j = i; j <= m; j++) a[j, i] = 0.0;
                }
                ++a[i, i];
            }
            for (k=n;k>=1;k--) {
                for (its=1;its<=30;its++) {
                    flag=1;
                    for (l=k;l>=1;l--) {
                        nm=l-1;
                        if (Math.Abs(rv1[l]) + anorm == anorm)
                        {
                            flag=0;
                            break;
                        }
                        if (Math.Abs(w[nm]) + anorm == anorm) break;
                    }
                    if (flag != 0) {
                        c=0.0;
                        s=1.0;
                        for (i=l;i<=k;i++) {
                            f=s*rv1[i];
                            if (Math.Abs(f) + anorm != anorm)
                            {
                                g=w[i];
                                h=PYTHAG(f,g);
                                w[i]=h;
                                h=1.0/h;
                                c=g*h;
                                s=(-f*h);
                                for (j=1;j<=m;j++) {
                                    y = a[j, nm];
                                    z = a[j, i];
                                    a[j, nm] = y * c + z * s;
                                    a[j, i] = z * c - y * s;
                                }
                            }
                        }
                    }
                    z=w[k];
                    if (l == k) {
                        if (z < 0.0) {
                            w[k] = -z;
                            for (j = 1; j <= n; j++) v[j, k] = (-v[j, k]);
                        }
                        break;
                    }
                    if (its == 30) throw new Exception("No convergence in 30 SVDCMP iterations");
                    x=w[l];
                    nm=k-1;
                    y=w[nm];
                    g=rv1[nm];
                    h=rv1[k];
                    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                    g=PYTHAG(f,1.0);
                    f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                    c=s=1.0;
                    for (j=l;j<=nm;j++) {
                        i=j+1;
                        g=rv1[i];
                        y=w[i];
                        h=s*g;
                        g=c*g;
                        z=PYTHAG(f,h);
                        rv1[j]=z;
                        c=f/z;
                        s=h/z;
                        f=x*c+g*s;
                        g=g*c-x*s;
                        h=y*s;
                        y=y*c;
                        for (jj=1;jj<=n;jj++) {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }
                        z=PYTHAG(f,h);
                        w[j]=z;
                        if (z != 0.0) {
                            z=1.0/z;
                            c=f*z;
                            s=h*z;
                        }
                        f=(c*g)+(s*y);
                        x=(c*y)-(s*g);
                        for (jj=1;jj<=m;jj++) {
                            y = a[jj, j];
                            z = a[jj, i];
                            a[jj, j] = y * c + z * s;
                            a[jj, i] = z * c - y * s;
                        }
                    }
                    rv1[l]=0.0;
                    rv1[k]=f;
                    w[k]=x;
                }
            }
        }
    }


    public static class BurkardtLib
    {
        /// <summary>
        ///  Purpose:
        ///    DMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
        ///  Discussion:
        ///    The doubly dimensioned array A is treated as a one dimensional vector,
        ///    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*N]
        ///  Modified:
        ///    29 August 2003
        ///  Author:
        ///    John Burkardt
        /// </summary>
        /// <param name="n">the order of the matrix</param>
        /// <param name="rhs_num">the number of right hand sides. RHS_NUM must be at least 0.</param>
        /// <param name="a">Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
        /// to N the coefficient matrix, and in columns N+1 through N+RHS_NUM, the right hand sides.
        /// On output, the coefficient matrix area has been destroyed, while the right hand sides have
        /// been overwritten with the corresponding solutions.</param>
        /// <returns>0, the matrix was not singular, the solutions were computed;
        /// J, factorization failed on step J, and the solutions could not be computed.</returns>
        public static int dmat_solve(int n, int rhs_num, double[] a)
        {
            double apivot;
            double factor;
            int i;
            int ipivot;
            int j;
            int k;
            double temp;

            for (j = 0; j < n; j++)
            {
                //
                //  Choose a pivot row.
                //
                ipivot = j;
                apivot = a[j + j * n];

                for (i = j; i < n; i++)
                {
                    if (Math.Abs(apivot) < Math.Abs(a[i + j * n]))
                    {
                        apivot = a[i + j * n];
                        ipivot = i;
                    }
                }

                if (apivot == 0.0)
                {
                    return j;
                }
                //
                //  Interchange.
                //
                for (i = 0; i < n + rhs_num; i++)
                {
                    temp = a[ipivot + i * n];
                    a[ipivot + i * n] = a[j + i * n];
                    a[j + i * n] = temp;
                }
                //
                //  A(J,J) becomes 1.
                //
                a[j + j * n] = 1.0;
                for (k = j; k < n + rhs_num; k++)
                {
                    a[j + k * n] = a[j + k * n] / apivot;
                }
                //
                //  A(I,J) becomes 0.
                //
                for (i = 0; i < n; i++)
                {
                    if (i != j)
                    {
                        factor = a[i + j * n];
                        a[i + j * n] = 0.0;
                        for (k = j; k < n + rhs_num; k++)
                        {
                            a[i + k * n] = a[i + k * n] - factor * a[j + k * n];
                        }
                    }
                }
            }
            return 0;
        }

        /// <summary>
        ///  Purpose:
        ///    DGE_DET computes the determinant of a matrix factored by SGE_FA.
        ///  Discussion:
        ///    The doubly dimensioned array A is treated as a one dimensional vector,
        ///    stored by COLUMNS:  
        ///      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        ///    Entry A(I,J) is stored as A[I+J*N]
        ///  Modified:
        ///    04 September 2003
        ///  Author:
        ///    John Burkardt
        ///  Reference:
        ///    Dongarra, Bunch, Moler, Stewart,
        ///    LINPACK User's Guide,
        ///    SIAM, 1979
        /// </summary>
        /// <param name="n">the order of the matrix. N must be positive.</param>
        /// <param name="a">A[N*N], the LU factors computed by DGE_FA.</param>
        /// <param name="pivot">PIVOT[N], as computed by DGE_FA.</param>
        /// <returns>the determinant of the matrix.</returns>
        public static double dge_det(int n, double[] a, int[] pivot)
        {
            double det;
            int i;

            det = 1.0;

            for (i = 0; i < n; i++)
            {
                det = det * a[i + i * n];
                if (pivot[i] != i + 1)
                {
                    det = -det;
                }
            }

            return det;
        }

        /// <summary>
        ///  Purpose:
        ///    DGE_FA factors a general matrix.
        ///  Discussion:
        ///    DGE_FA is a simplified version of the LINPACK routine SGEFA.
        ///    The doubly dimensioned array A is treated as a one dimensional vector,
        ///    stored by COLUMNS:  
        ///      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        ///    Entry A(I,J) is stored as A[I+J*N]
        ///  Modified:
        ///    05 September 2003
        ///  Author:
        ///    John Burkardt
        ///  Reference:
        ///    Dongarra, Bunch, Moler, Stewart,
        ///    LINPACK User's Guide,
        ///    SIAM, 1979
        /// </summary>
        /// <param name="n">the order of the matrix. N must be positive.</param>
        /// <param name="a">Input/output, A[N*N], the matrix to be factored. On output, A contains an upper triangular 
        /// matrix and the multipliers which were used to obtain it.  The factorization can be written
        /// A = L * U, where L is a product of permutation and unit lower triangular matrices and U is 
        /// upper triangular.</param>
        /// <param name="pivot">Output, int PIVOT[N], a vector of pivot indices.</param>
        /// <returns>0, no singularity detected. nonzero, the factorization failed on the DGE_FA-th step.</returns>
        public static int dge_fa(int n, double[] a, int[] pivot)
        {
            int i;
            int ii;
            int info;
            int j;
            int k;
            int l;
            double t;

            info = 0;

            for (k = 1; k <= n - 1; k++)
            {
                //
                //  Find L, the index of the pivot row.
                //
                l = k;
                for (i = k + 1; i <= n; i++)
                {
                    if (Math.Abs(a[l - 1 + (k - 1) * n]) < Math.Abs(a[i - 1 + (k - 1) * n]))
                    {
                        l = i;
                    }
                }

                pivot[k - 1] = l;
                //
                //  If the pivot index is zero, the algorithm has failed.
                //
                if (a[l - 1 + (k - 1) * n] == 0.0)
                {
                    throw new Exception(string.Format("DGE_FA - Fatal error!. Zero pivot on step {0}", k));
                }
                //
                //  Interchange rows L and K if necessary.
                //
                if (l != k)
                {
                    t = a[l - 1 + (k - 1) * n];
                    a[l - 1 + (k - 1) * n] = a[k - 1 + (k - 1) * n];
                    a[k - 1 + (k - 1) * n] = t;
                }
                //
                //  Normalize the values that lie below the pivot entry A(K,K).
                //
                for (j = k + 1; j <= n; j++)
                {
                    a[j - 1 + (k - 1) * n] = -a[j - 1 + (k - 1) * n] / a[k - 1 + (k - 1) * n];
                }
                //
                //  Row elimination with column indexing.
                //
                for (j = k + 1; j <= n; j++)
                {
                    if (l != k)
                    {
                        t = a[l - 1 + (j - 1) * n];
                        a[l - 1 + (j - 1) * n] = a[k - 1 + (j - 1) * n];
                        a[k - 1 + (j - 1) * n] = t;
                    }

                    for (ii = k; ii < n; ii++)
                    {
                        a[ii + (j - 1) * n] = a[ii + (j - 1) * n] + a[ii + (k - 1) * n] * a[k - 1 + (j - 1) * n];
                    }
                }
            }

            pivot[n - 1] = n;

            if (a[n - 1 + (n - 1) * n] == 0.0)
            {
                throw new Exception(string.Format("DGE_FA - Fatal error!. Zero pivot on step {0}", n));
            }

            return info;
        }

        /// <summary>
        ///  Purpose:
        ///    DGE_SL solves a system factored by SGE_FA.
        ///  Discussion:
        ///    DGE_SL is a simplified version of the LINPACK routine SGESL.
        ///    The doubly dimensioned array A is treated as a one dimensional vector,
        ///    stored by COLUMNS:  
        ///      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        ///    Entry A(I,J) is stored as A[I+J*N]
        ///  Modified:
        ///    06 September 2003
        ///  Author:
        ///    John Burkardt
        /// </summary>
        /// <param name="n">the order of the matrix. N must be positive.</param>
        /// <param name="a">A[N*N], the LU factors from DGE_FA.</param>
        /// <param name="pivot">PIVOT[N], the pivot vector from DGE_FA.</param>
        /// <param name="b">Input/output, B[N]. On input, the right hand side vector.
        /// On output, the solution vector.</param>
        /// <param name="job">specifies the operation: 0, solve A * x = b. nonzero, solve A' * x = b.</param>
        public static void dge_sl(int n, double[] a, int[] pivot, double[] b, int job)
        {
            int i;
            int j;
            int k;
            int l;
            double t;
            //
            //  Solve A * x = b.
            //
            if (job == 0)
            {
                //
                //  Solve PL * Y = B.
                //
                for (k = 1; k <= n - 1; k++)
                {
                    l = pivot[k - 1];

                    if (l != k)
                    {
                        t = b[l - 1];
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }

                    for (i = k + 1; i <= n; i++)
                    {
                        b[i - 1] = b[i - 1] + a[i - 1 + (k - 1) * n] * b[k - 1];
                    }
                }
                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] = b[k - 1] / a[k - 1 + (k - 1) * n];
                    for (i = 1; i <= k - 1; i++)
                    {
                        b[i - 1] = b[i - 1] - a[i - 1 + (k - 1) * n] * b[k - 1];
                    }
                }
            }
            //
            //  Solve A' * X = B.
            //
            else
            {
                //
                //  Solve U' * Y = B.
                //
                for (k = 1; k <= n; k++)
                {
                    t = 0.0;
                    for (i = 1; i <= k - 1; i++)
                    {
                        t = t + b[i - 1] * a[i - 1 + (k - 1) * n];
                    }
                    b[k - 1] = (b[k - 1] - t) / a[k - 1 + (k - 1) * n];
                }
                //
                //  Solve ( PL )' * X = Y.
                //
                for (k = n - 1; 1 <= k; k--)
                {
                    t = 0.0;
                    for (i = k + 1; i <= n; i++)
                    {
                        t = t + b[i - 1] * a[i - 1 + (k - 1) * n];
                    }
                    b[k - 1] = b[k - 1] + t;

                    l = pivot[k - 1];

                    if (l != k)
                    {
                        t = b[l - 1];
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }
                }
            }
            return;
        }
    }
}
