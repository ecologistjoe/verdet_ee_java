package com.google.earthengine.examples.experimental;

import com.google.earthengine.api.base.AlgorithmBase;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

//import org.ejml.data.DMatrixRMaj;
//import org.ejml.dense.row.CommonOps_DDRM;

/**
 * Implements the VeRDET change detection:
 *   Patch-Based Forest Change Detection from Landsat Time Series
 *   M. Joseph Hughes, S. Douglas Kaylor, Daniel J. Hayes
 *   Forests 2017, 8(5), 166; https://doi.org/10.3390/f8050166
 * <p>
 * Original Matlab code:
 *   find_distubance.m, fit_piecewise_linear.m, tvr_1dm_any.m
 *
 * This implementation performs 1-dimensional, 1-st derivative TVR denoising 
 * followed by a simple merging and re-interpolation of a vector representing
 * equally-spaced (nominally yearly time series) data.
 * For default parameters, data values should be between 0 and 1.
 * This implementation does not perform any spatial denoising or segmentation.
 * Additionally, a simpler merging rule is used here than in fit_piecewise_linear.m
 * since the original rule was found to be non-significant in change detection.
 * <p>
 * @author Zhiqiang Yang, 04/1/2015
 * @author Joe Hughes, 2019-03-16
* <p>
 */


public final class Verdet {
  static class Args extends AlgorithmBase.ArgsBase {
    @Doc(help = "Regularization parameter for temporal segmentation. Larger "
            + "values produce more de-noising and fewer segments.")
    @Optional
    double alpha = 1/20.0;

    @Doc(help = "convergence tolerance")
    @Optional
    double tolerance = 0.0001;

    @Doc(help = "Maximum number of runs for convergence.")
    @Optional
    int maxRuns = 100;
  }
  private final Args args;
  private int size;

  public Verdet() {
    this(new Args());
  }

  public Verdet(Args args) {
    this.args = args;
  }

  /**
   * Compute the verdet scores.
   * @param a scores calculated in verdet
   */
  public double[] getResult(double[] f) {
    init(f.length);
    // Perform 1-dimensional 1st derivative of
    // Total Variation Regularization over the input
    // to denoise and generate approximate segmentations

    DenseMatrix64F slopes = tvr1D(f);
    int[] SegIdx = mergeSegments(slopes);
    DenseMatrix64F A = buildInterpolationMatrix(SegIdx);

    // Can loop this step with different fit-to-vertex indices
    double[] X = fitToVertices(A, f);

    return X;
  }

  public void init(int size) {
    if (this.size == size) {
      return;
    }
    this.size = size;
  }

  /**Accepts a matrix of slopes, representing the segmentation of
   * a time series, then merges neighboring segments that are
   * separated by a distance of less than alpha
  **/
  public int[] mergeSegments(DenseMatrix64F slopes) {

    final double thresh = args.alpha;
    
    // Identify near-same-sloped segments from denoised variable
    // This can't be folded into the next loop because we need to
    // know segment count before building the interpolation matrix
    int[] SegIdx = new int[size];
    double[] deltas = new double[size-1];
    int segCnt = 0;
    int lastSegIdx = 1;
    SegIdx[0] = 0;
    
    for (int i = 1; i < size-1; i++) {
      if( Math.abs(slopes.get(i+1) - slopes.get(lastSegIdx)) > thresh) {
        segCnt++;
        lastSegIdx = i+1;
      }
      SegIdx[i] = segCnt;
    }
    // The last node is always part of, and ends the current segment.
    SegIdx[size-1] = segCnt+1;
    return SegIdx;
  }

  // Make an interpolation matrix given a strictly-increasing array
  // where each entry is the segment number
  public DenseMatrix64F buildInterpolationMatrix(int [] SegIdx) {

    // Build a matrix to solve for the piecewise linear model
    // where segments as defined as mixtures of endpoint vertices
    DenseMatrix64F interpolater = new DenseMatrix64F(size, SegIdx[size-1]+1);
    int lastSegIdx = 0;
    interpolater.zero();
    for (int i = 1; i < size; i++) {
      if(SegIdx[i] != SegIdx[lastSegIdx]) {
        for (int j = lastSegIdx; j < i; j++) {
          double weight = (double) (j-lastSegIdx) / (i-lastSegIdx);
          interpolater.set(j, SegIdx[j], 1 - weight);
          interpolater.set(j, SegIdx[j]+1,  weight);
        }
        lastSegIdx = i;
      }
    }
    interpolater.set(size-1, SegIdx[size-1], 1.0);
    return interpolater;

  }

  // Given a piecewise linear model, fit those segments to a set of points
  // these points don't necessarily need to be the same that generated the
  // model. This, for example, allows segementation using one nbr, but
  // getting the data for the red band.
  private double[] fitToVertices(DenseMatrix64F interpolater, double[] f) {

    double [] Out = new double[size];
    
    DenseMatrix64F vertices = new DenseMatrix64F(interpolater.getNumCols(), 1);
    DenseMatrix64F ff = new DenseMatrix64F(size, 1);
    for(int i = 0; i< size; i++) {
      ff.set(i, f[i]);
    }
    CommonOps.solve(interpolater, ff, vertices);

    int j = 0;
    //Interpolate f between the solved vertices (Out = interpolater*vertices)
    for (int i=0; i<size; i++) {
      Out[i] = interpolater.get(i, j)   * vertices.get(j)
             + interpolater.get(i, j+1) * vertices.get(j+1);
      if(interpolater.get(i,j) == 0) {
       j++;
      }
    }

    return Out;
  }


  /**Accepts a time series of points (f) and returns a nearly-piecewise linear
   * approximation of that time series, where alpha controls the
   * complexity of the approximation.  Higher alpha means more denoising.
  **/
  private DenseMatrix64F tvr1D(double[] f) {

    /**
     * Original MATLAB code vectorized iterations using
     * a anti-differentiation matrix, A: A=tril(ones(size)).
     * from this A'A was calculated, and stored as AtA,
     * which is a (size x size) matrix with values =
     * [ size   size-1 size-2 ... 1
     *   size-1 size-1 size-2 ... 1
     *   size-2 size-2 size-2 ... 1
     *     ...
     *   1      1       1     ... 1 ]
     * ie, assuming 0-based indexing:
     *    AtA(i,j) = size-max(i,j)
     *
     * Since the java code is looped, we can replace
     * all instances of A and AtA with cumulative sums
     * and functions of size and iteration
    **/

    @SuppressWarnings("ConstantField")
    DenseMatrix64F L = new DenseMatrix64F(size, size);
    DenseMatrix64F Atf = new DenseMatrix64F(size, 1);
    DenseMatrix64F u = new DenseMatrix64F(size, 1);
    DenseMatrix64F u1 = new DenseMatrix64F(size, 1);

    // Establish an initial guess, and remove time=0 intercept
    for (int i = 0; i < size; i++) {
      u.set(i, f[i]-f[0]);
    }
    // Precalculate the target: A'*f - f[0]
    Atf.set(size-1, u.get(size-1));
    for (int i = size-2; i >= 0; i--) {
      Atf.set(i, Atf.get(i+1) + u.get(i));
    }
    
    //initialize some initial variables
    CommonOps.fill(u1, Double.POSITIVE_INFINITY);

    //initialize L to AtA
    for (int i=0; i < size; i++) {
      for (int j=0; j<size; j++) {
        if(i>j) {
          L.set(i,j, size-i);
        } else {
          L.set(i,j, size-j);
        }
      }
    }
    
    // iterate until convergence or maxRuns
    for (int i = 0; i < args.maxRuns; i++) {
    
      double prev = 0;
      double curr = 0;
      for (int j = 0; j < size-1; j++) {
        // Calculate inverse of deltas between appoximated slopes,
        // weighted by 1/alpha
        curr = args.alpha / (1e-6 + Math.abs(u.get(j+1) - u.get(j)));

        // Build (AtA)+L
        L.set(j,   j,   (size-j) + prev + curr);
        L.set(j,   j+1, (size-j-1) - curr);
        L.set(j+1, j,   (size-j-1) - curr);

        prev = curr;
      }
      L.set(size-1, size-1, 1+prev);

      // Update approximate solution (u) solving L*u = Atf for u
      CommonOps.solve(L, Atf, u);
      
      // Check for convergence after a few iterations
      if(i > 20) {
        CommonOps.subtract(u1, u, u1);
        if (CommonOps.elementMaxAbs(u1) <= args.tolerance) {
          break;
        }
        copy(u, u1);
      }
    }

    return u;
  }

  void copy(DenseMatrix64F a, DenseMatrix64F b) {
    for (int i = 0; i < a.numRows; i++) {
      for (int j = 0; j < a.numCols; j++) {
          b.set(i, j, a.get(i, j));
      }
    }
  }
  
  public static void main(String[] ar) {
    Verdet V = new Verdet();
    double[] test = {.82, 0.78, 0.77, 0.86, 0.94, 0.95, 0.70, 0.78, 0.61, 0.42, 0.28, 0.18, 0.10, 0.10, 0.12, 0.24, 0.39, 0.43, 0.50, 0.70};
    double[] result = V.getResult(test);
    System.out.println(Arrays.toString(result));      
  }
 
}
