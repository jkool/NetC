package utilities;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.jet.random.*;

/**
 * <p>
 * Title: Matrix Utilities
 * </p>
 * <p>
 * Description: A collection of methods for creating and manipulating matrices
 * of primitives.
 * </p>
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * <p>
 * Company: RSMAS, University of Miami
 * </p>
 * 
 * @author Johnathan Kool
 * @version 0.1
 */

public class MatrixUtilities {
	
	private static final float R_EARTH = 6378137f;

	/**
	 * Adds two matrices of equal dimension together.
	 * 
	 * @param a
	 *            double[][] - The first matrix
	 * @param b
	 *            double[][] - The second matrix
	 * @return double[][] - The result of adding the two matrices.
	 */

	public static double[][] add(double[][] a, double[][] b) {

		double[][] sum = new double[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				sum[i][j] = a[i][j] + b[i][j];

			}

		}

		return sum;

	}

	/**
	 * Adds two matrices of equal dimension together.
	 * 
	 * @param a
	 *            int[][] - The first matrix
	 * @param b
	 *            int[][] - The second matrix
	 * @return int[][] - The result of adding the two matrices.
	 */

	public static int[][] add(int[][] a, int[][] b) {

		int[][] sum = new int[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				sum[i][j] = a[i][j] + b[i][j];

			}

		}

		return sum;

	}

	/**
	 * Adds two vectors of equal dimension together.
	 * 
	 * @param a
	 *            int[][] - The first matrix
	 * @param b
	 *            int[][] - The second matrix
	 * @return int[][] - The result of adding the two matrices.
	 */

	public static long[] add(long[] a, int[] b) {

		long[] sum = new long[a.length];

		for (int i = 0; i < a.length; i++) {

			sum[i] = a[i] + b[i];
		}

		return sum;
	}
	
	public static float[] add(float[] a, float[] b) {

		float[] sum = new float[a.length];

		for (int i = 0; i < a.length; i++) {

			sum[i] = a[i] + b[i];
		}

		return sum;
	}
	
	public static double[] add(double[] a, double[] b) {

		double[] sum = new double[a.length];

		for (int i = 0; i < a.length; i++) {

			sum[i] = a[i] + b[i];
		}

		return sum;
	}

	/**
	 * Adds two vectors of equal dimension together.
	 * 
	 * @param a
	 *            int[][] - The first matrix
	 * @param b
	 *            int[][] - The second matrix
	 * @return int[][] - The result of adding the two matrices.
	 */

	public static long[] add(long[] a, long[] b) {

		long[] sum = new long[a.length];

		for (int i = 0; i < a.length; i++) {

			sum[i] = a[i] + b[i];
		}

		return sum;
	}

	/**
	 * Adds two matrices of equal dimension together.
	 * 
	 * @param a
	 *            int[][] - The first matrix
	 * @param b
	 *            int[][] - The second matrix
	 * @return int[][] - The result of adding the two matrices.
	 */

	public static long[][] add(long[][] a, long[][] b) {

		long[][] sum = new long[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				sum[i][j] = a[i][j] + b[i][j];

			}

		}

		return sum;

	}

	/**
	 * Combines two functions of standardizing and returning the cumulative
	 * function of the provided array.
	 * 
	 * @param da
	 *            double[]: The array for standardization and accumulation.
	 * @return double[]: The standardized, cumulative array.
	 */

	public static double[] cdf(double[] da) {

		return cumulative(norm1(da));

	}

	/**
	 * Generates a matrix with circulant form
	 * 
	 * @param vals -
	 *            a double array of values
	 * @param dim -
	 *            The dimension of the array
	 * @param shift -
	 *            a shift parameter for the array
	 * @return - The circulant matrix
	 */

	public static double[][] circulant(double[] vals, int dim) {

		double[][] out = new double[dim][dim];
		int vl = vals.length;
		int mod = (int) Math.floor((double) vl / 2);

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {

				if (i == j) {

					for (int k = 0; k < vl; k++) {

						if (mod - k > j) {

							out[i][dim + (j + k - mod)] = vals[k];
						}

						else if (j + k - mod >= dim) {
							out[i][-(dim - (j + k - mod))] = vals[k];
						}

						else {
							out[i][j + k - mod] = vals[k];
						}

					}
				}
			}
		}

		return out;

	}

	/**
	 * Generates a matrix with circulant form
	 * 
	 * @param vals -
	 *            a double array of values
	 * @param dim -
	 *            The dimension of the array
	 * @param shift -
	 *            a shift parameter for the array
	 * @return - The circulant matrix
	 */

	public static int[][] circulant(int[] vals, int dim) {

		int[][] out = new int[dim][dim];
		int vl = vals.length;
		int mod = (int) Math.floor((double) vl / 2);

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {

				if (i == j) {

					for (int k = 0; k < vl; k++) {

						if (mod - k > j) {

							out[i][dim + (j + k - mod)] = vals[k];
						}

						else if (j + k - mod >= dim) {
							out[i][-(dim - (j + k - mod))] = vals[k];
						}

						else {
							out[i][j + k - mod] = vals[k];
						}

					}
				}
			}
		}

		return out;

	}

	/**
	 * Compares two matrices for equivalency
	 * 
	 * @param a -
	 *            The first matrix
	 * @param b -
	 *            The second matrix
	 * @return - Indicates whether the two matrices are equivalent
	 */

	public static boolean compare(double[] a, double[] b) {

		if (a.length != b.length) {
			return false;
		}

		for (int i = 0; i < a.length; i++) {

			if (a[i] != b[i]) {
				return false;
			}
		}

		return true;

	}

	/**
	 * Compares two matrices for equivalency
	 * 
	 * @param a -
	 *            The first matrix
	 * @param b -
	 *            The second matrix
	 * @return - Indicates whether the two matrices are equivalent
	 */

	public static boolean compare(int[] a, int[] b) {

		if (a.length != b.length) {
			return false;
		}

		for (int i = 0; i < a.length; i++) {

			if (a[i] != b[i]) {
				return false;
			}
		}

		return true;

	}

	/**
	 * Concatenates (appends) two arrays horizontally (must have the same number
	 * of rows).
	 * 
	 * @param a -
	 *            The first array
	 * @param b -
	 *            The second array
	 * @return - The concatenated arrays
	 */

	public static double[][] concatenate(double[][] a, double[][] b) {
		int catlen = a[0].length + b[0].length;
		double[][] out = new double[a.length][catlen];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < catlen; j++) {

				if (j < a[0].length) {
					out[i][j] = a[i][j];
				} else {
					out[i][j] = b[i][j - a[0].length];

				}
			}
		}

		return out;

	}

	/**
	 * Concatenates (appends) two arrays horizontally (must have the same number
	 * of rows).
	 * 
	 * @param a -
	 *            The first array
	 * @param b -
	 *            The second array
	 * @return - The concatenated arrays
	 */

	public static long[][] concatenate(long[][] a, long[][] b) {
		int catlen = a[0].length + b[0].length;
		long[][] out = new long[a.length][catlen];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < catlen; j++) {

				if (j < a[0].length) {
					out[i][j] = a[i][j];
				} else {
					out[i][j] = b[i][j - a[0].length];

				}
			}
		}

		return out;

	}

	/**
	 * Returns a vector containing a constant value
	 * 
	 * @param length
	 *            int - the length of the desired vector
	 * @param value
	 *            double - an int value for the vector
	 * @return double[] - The constant value vector
	 */

	public static double[] constantVec(int length, double value) {

		double[] out = new double[length];
		for (int i = 0; i < out.length; i++) {

			out[i] = value;
		}

		return out;
	}

	/**
	 * Returns a vector containing a constant value
	 * 
	 * @param length
	 *            int - the length of the desired vector
	 * @param value
	 *            int - an int value for the vector
	 * @return int[] - The constant value vector
	 */

	public static int[] constantVec(int length, int value) {

		int[] out = new int[length];
		for (int i = 0; i < out.length; i++) {

			out[i] = value;
		}

		return out;

	}

	/**
	 * Returns a vector containing a constant value
	 * 
	 * @param length
	 *            int - the length of the desired vector
	 * @param value
	 *            long - an int value for the vector
	 * @return int[] - The constant value vector
	 */

	public static long[] constantVec(int length, long value) {

		long[] out = new long[length];
		for (int i = 0; i < out.length; i++) {

			out[i] = value;
		}

		return out;

	}

	/**
	 * Compute the cross product of two vectors
	 * 
	 * @param v1
	 *            The first vector
	 * @param v2
	 *            The second vector
	 * @param result
	 *            Where to store the cross product
	 **/
	public static double[] cross(double[] p1, double[] p2) {
		double[] result = new double[3];
		result[0] = p1[1] * p2[2] - p2[1] * p1[2];
		result[1] = p1[2] * p2[0] - p2[2] * p1[0];
		result[2] = p1[0] * p2[1] - p2[0] * p1[1];
		return result;
	}

	/**
	 * Compute the cross product of two vectors
	 * 
	 * @param v1
	 *            The first vector
	 * @param v2
	 *            The second vector
	 * @param result
	 *            Where to store the cross product
	 **/
	public static float[] cross(float[] p1, float[] p2) {
		float[] result = new float[3];
		result[0] = p1[1] * p2[2] - p2[1] * p1[2];
		result[1] = p1[2] * p2[0] - p2[2] * p1[0];
		result[2] = p1[0] * p2[1] - p2[0] * p1[1];
		return result;
	}
	
	/**
	 * Returns the cumulative form of the supplied array,
	 * 
	 * @param da
	 *            double[]: The array to be accumulated
	 * @return double[]: The cumulative form of the supplied array.
	 */

	public static double[] cumulative(double[] da) {

		double[] cm = new double[da.length];
		double sum = 0;

		for (int i = 0; i < da.length; i++) {

			sum += da[i];
			cm[i] = sum;

		}

		return cm;

	}

	/**
	 * Returns a matrix containing the cumulative probabilities of the
	 * transition matrix. The cumulative matrix is used in conjunction with a
	 * randomly generated number to allocate destinations.
	 * 
	 * @param matrix
	 *            double[][]: Converts a matrix of probabilities to cumulative
	 *            probabilities.
	 * @return double[][]: The cumulative probability matrix.
	 */

	public static double[][] cumulativeMatrix(double[][] matrix) {

		int matrixDim = matrix.length;
		double[][] cumulativeMatrix = new double[matrixDim][matrixDim];

		for (int i = 0; i < matrixDim; i++) {

			double val = 0;

			for (int j = 0; j < matrixDim; j++) {

				val += matrix[i][j];
				cumulativeMatrix[i][j] = val;

			}

		}

		return cumulativeMatrix;

	}

	/**
	 * Converts a vector into 2D matrix format
	 * 
	 * @param vec -
	 *            The input vector
	 * @return - The vector in 2D form.
	 */

	public static int[][] diag(int[] vec) {

		int len = vec.length;

		int[][] out = new int[len][len];

		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {

				if (i == j) {
					out[i][j] = vec[i];
				}
			}
		}

		return out;

	}

	/**
	 * Converts a vector into 2D matrix format
	 * 
	 * @param vec -
	 *            The input vector
	 * @return - The vector in 2D form.
	 */

	public static long[][] diag(long[] vec) {

		int len = vec.length;

		long[][] out = new long[len][len];

		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {

				if (i == j) {
					out[i][j] = vec[i];
				}
			}
		}

		return out;

	}

	public static float[] dilate(float[] arr, float multiplier){
		float[] out = new float[arr.length];
		for(int i = 0; i < arr.length; i++){
			out[i] = arr[i]*multiplier;
		}
		return out;
	}
	
	public static double[] dilate(double[] arr, double multiplier){
		double[] out = new double[arr.length];
		for(int i = 0; i < arr.length; i++){
			out[i] = arr[i]*multiplier;
		}
		return out;
	}
	
	/**
	 * Performs element-wise division of two arrays (vectors)
	 * 
	 * @param num -
	 *            the numerator vector
	 * @param denom -
	 *            the denominator vector
	 * @return - resulting vector
	 */

	public static double[] divide(long[] num, long[] denom) {

		int len = num.length;
		assert len == denom.length;
		if (len != denom.length) {
			throw new IllegalArgumentException(
					"Arrays must be of equal size for division");
		}

		double[] output = new double[len];

		for (int i = 0; i < len; i++) {

			output[i] = (double) num[i] / (double) denom[i];
		}

		return output;
	}

	/**
	 * Converts a matrix of doubles to integer values
	 * 
	 * @param a -
	 *            The matrix to be converted
	 * @return - The equivalent matrix as int values
	 */

	public static int[][] double2Int(double[][] a) {

		int[][] out = new int[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				out[i][j] = (int) a[i][j];

			}
		}

		return out;
	}

	/**
	 * Converts a matrix of doubles to long values.
	 * 
	 * @param a -
	 *            The matrix to be converted
	 * @return - The equivalent matrix as long values.
	 */

	public static long[][] double2Long(double[][] a) {

		long[][] out = new long[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				out[i][j] = (long) a[i][j];

			}

		}

		return out;
	}

	/**
	 * Generates a matrix of identical values where rows sum to 1
	 * 
	 * @param dim -
	 *            The dimensions of the matrix
	 * @return - The output matrix.
	 */

	public static double[][] evenDouble(int dim) {

		double[][] out = new double[dim][dim];

		double val = 1 / (double) dim;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				out[i][j] = val;
			}
		}
		return out;
	}

	/**
	 * Performs the actual work of generating an identity matrix (i.e. a matrix
	 * with values of 1 along the diagonal, top left to bottom right).
	 * 
	 * @param dim
	 *            int: The dimension of the matrix to be created
	 * @return double[][]: The identity matrix
	 */

	public static double[][] eyeAsDouble(int dim) {

		double[][] matrix = new double[dim][dim];

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				if (i == j) {
					matrix[i][j] = 1;
				}

			}

		}

		return matrix;

	}

	/**
	 * Generates an identity matrix (i.e. a matrix with values of 1 along the
	 * diagonal, top left to bottom right).
	 * 
	 * @param dim
	 *            int: The dimension of the matrix to be created
	 * @return int[][]: The identity matrix
	 */

	public static int[][] eyeAsInt(int dim) {

		int[][] matrix = new int[dim][dim];

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				if (i == j) {
					matrix[i][j] = 1;
				}

			}

		}

		return matrix;

	}

	/**
	 * Generates an identity matrix (i.e. a matrix with values of 1 along the
	 * diagonal, top left to bottom right).
	 * 
	 * @param dim
	 *            int: The dimension of the matrix to be created
	 * @return int[][]: The identity matrix
	 */

	public static long[][] eyeAsLong(int dim) {

		long[][] matrix = new long[dim][dim];

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				if (i == j) {
					matrix[i][j] = 1;
				}

			}

		}

		return matrix;
	}

	/**
	 * Guaranteed upper bound for a binary search of an array. If values are
	 * equivalent, default binarySearch can return any value randomly.
	 * 
	 * @param da -
	 *            The array of values to be searched
	 * @param target -
	 *            the target being searched for
	 * @return
	 */

	private static int gBinary_ub(double[] da, double target) {

		int idx = Arrays.binarySearch(da, target);

		if (idx == -1) {
			return -1;
		}

		// If the value is negative...

		if (idx < -1) {
			return Math.abs(idx) - 2;
		}

		while (Math.abs(idx) != da.length - 1 && da[idx] == da[idx + 1]) {
			idx++;
		}

		return idx;

	}

	/**
	 * Guaranteed upper bound for a binary search of an array. If values are
	 * equivalent, default binarySearch can return any value randomly.
	 * 
	 * @param ia -
	 *            The array of values to be searched
	 * @param target -
	 *            The target being searched for
	 * @return
	 */

	private static int gBinary_ub(int[] ia, int target) {

		int idx = Arrays.binarySearch(ia, target);

		// If the value is negative...

		if (idx == -1) {
			return -1;
		}

		if (idx < -1) {
			return Math.abs(idx) - 2;
		}

		while (Math.abs(idx) != ia.length - 1 && ia[idx] == ia[idx + 1]) {
			idx++;
		}

		return idx;

	}

	/**
	 * Guaranteed upper bound for a binary search of an array. If values are
	 * equivalent, default binarySearch can return any value randomly.
	 * 
	 * @param la -
	 *            The array of values to be searched
	 * @param target -
	 *            The target being searched for
	 * @return
	 */

	private static int gBinary_ub(long[] la, long target) {

		int idx = Arrays.binarySearch(la, target);

		// If the value is negative...

		if (idx == -1) {
			return -1;
		}

		if (idx < -1) {
			return Math.abs(idx) - 2;
		}

		while (Math.abs(idx) != la.length - 1 && la[idx] == la[idx + 1]) {
			idx++;
		}

		return idx;

	}

	/**
	 * Retrieves a column vector of doubles from an integer matrix.
	 * 
	 * @param mtx
	 *            double[][] - The integer matrix from which the value is to be
	 *            retrieved.
	 * @param col
	 *            int - The identifier of the column to be retrieved.
	 * @return int[] - The column vector.
	 */

	public static double[] getColumn(double[][] mtx, int col) {

		double[] out = new double[mtx[0].length];

		for (int i = 0; i < mtx[0].length; i++) {

			out[i] = mtx[i][col];

		}

		return out;
	}

	/**
	 * Retrieves a column vector of ints from an integer matrix.
	 * 
	 * @param mtx
	 *            int[][] - The integer matrix from which the value is to be
	 *            retrieved.
	 * @param col
	 *            int - The identifier of the column to be retrieved.
	 * @return int[] - The column vector.
	 */

	public static int[] getColumn(int[][] mtx, int col) {

		int[] out = new int[mtx.length];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = mtx[i][col];

		}

		return out;
	}

	/**
	 * Retrieves a column vector of longs from an integer matrix.
	 * 
	 * @param mtx
	 *            long[][] - The long matrix from which the value is to be
	 *            retrieved.
	 * @param col
	 *            int - The identifier of the column to be retrieved.
	 * @return long[] - The column vector.
	 */

	public static long[] getColumn(long[][] mtx, int col) {

		long[] out = new long[mtx[0].length];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = mtx[i][col];

		}

		return out;
	}

	/**
	 * Returns the sum of a row of the transition matrix. Used in standardizing
	 * the matrix.
	 * 
	 * @param matrix
	 *            double[][]: The matrix containing the row to be summed.
	 * @param row
	 *            int: The row being summed
	 * @return double: The sum of the row
	 */

	public static double getRowSum(double[][] matrix, int row) {

		double sum = 0.0d;
		for (int i = 0; i < matrix[row].length; i++) {

			sum += matrix[row][i];
		}

		return sum;

	}

	/**
	 * Returns the sum of a row of the transition matrix. Used in standardizing
	 * the matrix.
	 * 
	 * @param matrix
	 *            long[][]: The matrix containing the row to be summed.
	 * @param row
	 *            int: The row being summed
	 * @return double: The sum of the row
	 */

	public static double getRowSum(long[][] matrix, int row) {

		double sum = 0.0d;
		for (int i = 0; i < matrix[row].length; i++) {

			sum += matrix[row][i];
		}

		return sum;
	}

	/**
	 * Converts a list of values into frequencies based on threshold
	 * values,similar to constructing a histogram
	 * 
	 * @param da -
	 *            the array of values
	 * @param th -
	 *            threshold values
	 * @return - the number of elements falling within each bounded class.
	 */

	public static int[] hist(double[] da, double[] th) {

		Arrays.sort(da);
		Arrays.sort(th);
		int[] out = new int[th.length];
		int acc = 0;

		// This will return the cumulative. Subtract by the previous entry

		for (int i = 0; i < th.length; i++) {

			out[i] = gBinary_ub(da, th[i]) + 1 - acc;
			acc = acc + out[i];
			System.out.print("");
		}

		return out;

	}

	/**
	 * Converts a list of values into frequencies based on threshold
	 * values,similar to constructing a histogram
	 * 
	 * @param ia -
	 *            the array of values
	 * @param th -
	 *            threshold values
	 * @return - the number of elements falling within each bounded class.
	 */

	public static int[] hist(int[] ia, int[] th) {

		Arrays.sort(ia);
		Arrays.sort(th);
		int[] out = new int[th.length];
		int acc = 0;

		// This will return the cumulative. Subtract by the previous entry

		for (int i = 0; i < th.length; i++) {

			out[i] = Math.abs(gBinary_ub(ia, th[i])) + 1 - acc;
			acc = acc + out[i];
			System.out.print("");
		}

		return out;

	}

	/**
	 * Converts a list of values into frequencies based on threshold
	 * values,similar to constructing a histogram
	 * 
	 * @param la -
	 *            the array of values
	 * @param th -
	 *            threshold values
	 * @return - the number of elements falling within each bounded class.
	 */

	public static long[] hist(long[] la, long[] th) {

		Arrays.sort(la);
		Arrays.sort(th);
		long[] out = new long[th.length];
		long acc = 0;

		// This will return the cumulative. Subtract by the previous entry

		for (int i = 0; i < th.length; i++) {

			out[i] = Math.abs(gBinary_ub(la, th[i])) + 1 - acc;
			acc = acc + out[i];
			System.out.print("");
		}

		return out;

	}

	/**
	 * Converts an array of integers to doubles.
	 * 
	 * @param a -
	 *            The matrix to be converted
	 * @return - The corresponding matrix as doubles.
	 */

	public static double[][] intToDouble(int[][] a) {

		double[][] out = new double[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				out[i][j] = a[i][j];
			}
		}

		return out;

	}

	/**
	 * Converts a matrix of integers into longs.
	 * 
	 * @param ia -
	 *            The matrix of integer values
	 * @return - The corresponding matrix as long values.
	 */

	public static long[] intToLong(int[] ia) {

		long[] out = new long[ia.length];
		for (int i = 0; i < ia.length; i++) {

			out[i] = ia[i];
		}

		return out;

	}

	/**
	 * Load a vector from a file
	 * 
	 * @param f -
	 *            Thefile object containing the vector.
	 * @return - An integer array containing the values from the file
	 */

	public static int[] intVectorFromFile(File f) {

		String ln;
		int[] out = new int[0];

		try {
			FileReader fr = new FileReader(f);
			LineNumberReader lnr = new LineNumberReader(fr);

			ln = lnr.readLine();

			while (ln != null) {

				out = new int[lnr.getLineNumber()];
				ln = lnr.readLine();

			}

			fr = new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

			ln = br.readLine();
			int count = 0;

			while (ln != null) {

				out[count] = Integer.parseInt(ln);
				ln = br.readLine();
				count++;
			}

			br.close();
			fr.close();

		} catch (FileNotFoundException ex) {
		} catch (IOException ex) {
		}

		return out;

	}

	/**
	 * Creates a vector of linearly spaced values
	 * 
	 * @param start -
	 *            Starting element value
	 * @param end -
	 *            Ending element value
	 * @param intervals -
	 *            Number of elements
	 * @return - The vector of linearly spaced values
	 */
	
	public static BigDecimal[] linspace(BigDecimal start, BigDecimal end, BigDecimal intervals) {

		BigDecimal range = end.subtract(start);
		BigDecimal spacing = range.divide(intervals.subtract(BigDecimal.ONE),8,BigDecimal.ROUND_HALF_EVEN);
		BigDecimal[] out = new BigDecimal[intervals.intValue()];

		out[0] = start;

		for (int i = 1; i < intervals.intValue(); i++) {

			BigDecimal val = spacing.add(out[i-1]);
			out[i] = val;

		}

		return out;

	}
	
	public static double[] linspace(double start, double end, int intervals) {

		double [] out = new double[intervals];
		BigDecimal[] bda = linspace(new BigDecimal(start), new BigDecimal(end), new BigDecimal(intervals));
		
		for(int i = 0; i < intervals; i++){
			out[i] = bda[i].doubleValue();
		}

		return out;
		
	}
	
	public static float[] linspace(float start, float end, int intervals) {

		float [] out = new float[intervals];
		BigDecimal[] bda = linspace(new BigDecimal(start), new BigDecimal(end), new BigDecimal(intervals));
		
		for(int i = 0; i < intervals; i++){
			out[i] = bda[i].floatValue();
		}

		return out;
		
	}
	
	public static int[] linspace(int start, int end, int intervals) {

		int [] out = new int[intervals];
		BigDecimal[] bda = linspace(new BigDecimal(start), new BigDecimal(end), new BigDecimal(intervals));
		
		for(int i = 0; i < intervals; i++){
			out[i] = bda[i].intValue();
		}

		return out;
		
	}
	

	/**
	 * Loads a transition matrix into the class based on a stored text file.
	 * 
	 * @param file
	 *            File: The file containing the ASCII representation of the
	 *            transition matrix.
	 * 
	 * @return double[][]:The ASCII matrix converted to double[][] form.
	 */

	public static double[][] loadASCIIMatrix(File file) {

		double[][] matrix = null;

		try {

			FileReader fr = new FileReader(file);
			BufferedReader b = new BufferedReader(fr);

			String rl = b.readLine();

			if (rl == null) {

				throw new java.lang.NullPointerException(
						"Supplied matrix is empty");
			}

			StringTokenizer stk = new StringTokenizer(rl);
			int numTokens = stk.countTokens();
			matrix = new double[numTokens][numTokens];

			int i = 0;

			while (rl != null) {

				int j = 0;

				while (stk.hasMoreTokens()) {

					matrix[i][j] = Double.parseDouble(stk.nextToken());
					j++;

				}

				rl = b.readLine();
				if (rl != null) {
					stk = new StringTokenizer(rl);
				}
				i++;

			}
		} catch (NumberFormatException ex) {
			ex.printStackTrace();
			System.out.println();
			System.out.println("Incorrect Number Format.  Exiting.");
			System.exit(-1);

		} catch (FileNotFoundException ex) {
			ex.printStackTrace();
			System.out.println();
			System.out.println("File Not Found.  Exiting.");
			System.exit(-1);

		} catch (IOException ex) {
			ex.printStackTrace();
			System.out.println();
			System.out
					.println("Exception occurred when reading the file.  Exiting.");
			System.exit(-1);

		}

		return matrix;

	}

	/**
	 * Casts a matrix of long values into integers
	 * 
	 * @param la -
	 *            The matrix of long values
	 * @return - The matrix as corresponding int values
	 */

	public static int[] longToInt(long[] la) {

		int[] out = new int[la.length];
		for (int i = 0; i < la.length; i++) {

			out[i] = (int) la[i];
		}

		return out;
	}

	/**
	 * Converts an array of integers to doubles.
	 * 
	 * @param a -
	 *            The matrix to be converted
	 * @return - The corresponding matrix as doubles.
	 */

	public static int[][] longToInt(long[][] a) {

		int[][] out = new int[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				out[i][j] = (int) a[i][j];
			}
		}

		return out;

	}

	/**
	 * Converts the matrix into String form
	 * 
	 * @param om -
	 *            The matrix to be converted
	 * @return - String representation of the matrix
	 */

	public static String matrixToString(double[][] om) {

		StringBuffer sb = new StringBuffer();
		sb.append("\n\n");

		for (double[] ov : om) {
			for (double i : ov) {

				sb.append(i + "  ");
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	/**
	 * Converts the matrix into String form
	 * 
	 * @param om -
	 *            The matrix to be converted
	 * @return - String representation of the matrix
	 */

	public static String matrixToString(int[][] om) {

		StringBuffer sb = new StringBuffer();
		sb.append("\n\n");

		for (int[] ov : om) {
			for (int i : ov) {

				sb.append(i + "  ");
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	/**
	 * Converts the matrix into String form
	 * 
	 * @param om -
	 *            The matrix to be converted
	 * @return - String representation of the matrix
	 */

	public static String matrixToString(long[][] om) {

		StringBuffer sb = new StringBuffer();
		sb.append("\n\n");

		for (long[] ov : om) {
			for (long i : ov) {

				sb.append(i + "  ");
			}
			sb.append("\n");
		}

		return sb.toString();
	}

	/**
	 * Performs matrix multiplication
	 * 
	 * @param a -
	 *            The left-hand matrix
	 * @param b -
	 *            The right-hand matrix
	 * @return - The matrix multiplication result.
	 */

	public static int[][] multiply(int[][] a, int[][] b) {

		int dim1 = a[0].length;
		int dim2 = b.length;

		assert a.length == b[0].length;

		int[][] out = new int[dim1][dim2];

		for (int i = 0; i < dim1; i++) {
			for (int j = 0; j < dim2; j++) {
				for (int k = 0; k < a.length; k++) {

					out[i][j] += a[i][j] * b[i][j];

				}
			}
		}

		return out;

	}

	/**
	 * Compute the cross product of two vectors
	 * 
	 * @param v1
	 *            The first vector
	 * @param v2
	 *            The second vector
	 * @param result
	 *            Where to store the cross product
	 **/
	public static double[] ncross(double[] p1, double[] p2) {
		double[] result = new double[3];
		result[0] = -(p1[1] * p2[2] - p2[1] * p1[2]);
		result[1] = -(p1[2] * p2[0] - p2[2] * p1[0]);
		result[2] = -(p1[0] * p2[1] - p2[0] * p1[1]);
		return result;
	}

	/**
	 * Compute the cross product of two vectors
	 * 
	 * @param v1
	 *            The first vector
	 * @param v2
	 *            The second vector
	 * @param result
	 *            Where to store the cross product
	 **/
	public static float[] ncross(float[] p1, float[] p2) {
		float[] result = new float[3];
		result[0] = -(p1[1] * p2[2] - p2[1] * p1[2]);
		result[1] = -(p1[2] * p2[0] - p2[2] * p1[0]);
		result[2] = -(p1[0] * p2[1] - p2[0] * p1[1]);
		return result;
	}

	public static float[] negative(float [] f){
		return dilate(f,-1f);
	}
	
	public static double magnitude(double[] f){
		double sum = 0;
		for(int i = 0; i < f.length; i++){
			sum+=(f[i]*f[i]);
		}
		return Math.sqrt(sum);
	}
	
	public static float magnitude(float[] f){
		float sum = 0;
		for(int i = 0; i < f.length; i++){
			sum+=(f[i]*f[i]);
		}
		return (float) Math.sqrt(sum);
	}
	
	public static float[] normalize(float[] f){
		float n = magnitude(f);
		float[] out = new float[f.length];
		for(int i = 0; i < f.length; i++){
			out[i] = f[i]/n;
		}
		return out;
	}
	
	public static double[] normalize(double[] f){
		double n = magnitude(f);
		double[] out = new double[f.length];
		for(int i = 0; i < f.length; i++){
			out[i] = f[i]/n;
		}
		return out;
	}
	
	/**
	 * Utility algorithm providing a simple console printout of a matrix.
	 * 
	 * @param matrix
	 *            double[][]: The matrix to be printed
	 */

	public static void printMatrix(double[][] matrix) {

		StringBuffer sb = new StringBuffer("");

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				sb.append(matrix[i][j]);
				sb.append("\t");

			}

			sb.append("\n");

		}

		System.out.println(sb.toString());

	}

	/**
	 * Utility algorithm providing a simple console printout of a matrix.
	 * 
	 * @param matrix
	 *            int[][]: The matrix to be printed
	 */

	public static void printMatrix(int[][] matrix) {

		StringBuffer sb = new StringBuffer("");

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				sb.append(matrix[i][j]);
				sb.append("\t");

			}

			sb.append("\n");

		}

		System.out.println(sb.toString());

	}

	/**
	 * Utility algorithm providing a simple console printout of a matrix.
	 * 
	 * @param matrix
	 *            double[][]: The matrix to be printed
	 */

	public static void printMatrix(long[][] matrix) {

		StringBuffer sb = new StringBuffer("");

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				sb.append(matrix[i][j]);
				sb.append("\t");

			}

			sb.append("\n");

		}

		System.out.println(sb.toString());

	}
	
	/**
	 * Rounds an array of doubles to long values.
	 * 
	 * @param daa -
	 *            The array to be rounded
	 * @return - The rounded array.
	 */

	public static long[][] round(double[][] daa) {
		long[][] out = new long[daa.length][daa[0].length];
		for (int i = 0; i < daa.length; i++) {
			for (int j = 0; j < daa[0].length; j++) {
				out[i][j] = Math.round(daa[i][j]);
			}
		}
		return out;
	}
	
	public static void saveToFile(int[][] iaa,String filename){
		
		try {
			PrintWriter pw = new PrintWriter(new File(filename));
				for(int i = 0; i < iaa.length;i++){
					for(int j = 0; j < iaa[i].length;j++){
					pw.write(iaa[i][j] + "\t"); 
					}
					pw.write("\n");
				}
			pw.flush();
			pw.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}

	/**
	 * Multiplies the supplied matrix by a scalar value (double)
	 * 
	 * @param scalar
	 *            double - The scalar value by which the matrix will be
	 *            multiplied.
	 * @param mtx
	 *            double[][] - The matrix to be multiplied.
	 * @return double[][] - The resulting output.
	 */

	public double[][] scalarX(double scalar, double[][] mtx) {

		double[][] out = new double[mtx.length][];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = new double[mtx[i].length];

			for (int j = 0; j < mtx[i].length; j++) {

				out[i][j] = scalar * mtx[i][j];
			}

		}

		return out;

	}

	/**
	 * Multiplies the supplied matrix by a scalar value (double)
	 * 
	 * @param scalar
	 *            double - The scalar value by which the matrix will be
	 *            multiplied.
	 * @param mtx
	 *            int[][] - The matrix to be multiplied.
	 * @return double[][] - The resulting output.
	 */

	public double[][] scalarX(double scalar, int[][] mtx) {

		double[][] out = new double[mtx.length][];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = new double[mtx[i].length];

			for (int j = 0; j < mtx[i].length; j++) {

				out[i][j] = scalar * mtx[i][j];
			}

		}

		return out;

	}

	/**
	 * Multiplies the supplied matrix by a scalar value (double)
	 * 
	 * @param scalar
	 *            double - The scalar value by which the matrix will be
	 *            multiplied.
	 * @param mtx
	 *            long[][] - The matrix to be multiplied.
	 * @return double[][] - The resulting output.
	 */

	public double[][] scalarX(double scalar, long[][] mtx) {

		double[][] out = new double[mtx.length][];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = new double[mtx[i].length];

			for (int j = 0; j < mtx[i].length; j++) {

				out[i][j] = scalar * mtx[i][j];
			}

		}

		return out;

	}

	/**
	 * Multiplies the supplied matrix by a scalar value (int)
	 * 
	 * @param scalar
	 *            int - The scalar value by which the matrix will be multiplied.
	 * @param mtx
	 *            int[][] - The matrix to be multiplied.
	 * @return int[][] - The resulting output.
	 */

	public int[][] scalarX(int scalar, int[][] mtx) {

		int[][] out = new int[mtx.length][];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = new int[mtx[i].length];

			for (int j = 0; j < mtx[i].length; j++) {

				out[i][j] = scalar * mtx[i][j];
			}

		}

		return out;

	}

	/**
	 * Multiplies the supplied matrix by a scalar value (int)
	 * 
	 * @param scalar
	 *            int - The scalar value by which the matrix will be multiplied.
	 * @param mtx
	 *            long[][] - The matrix to be multiplied.
	 * @return long[][] - The resulting output.
	 */

	public long[][] scalarX(long scalar, long[][] mtx) {

		long[][] out = new long[mtx.length][];

		for (int i = 0; i < mtx.length; i++) {

			out[i] = new long[mtx[i].length];

			for (int j = 0; j < mtx[i].length; j++) {

				out[i][j] = scalar * mtx[i][j];
			}

		}

		return out;

	}

	/**
	 * Standardizes an array of doubles to values between 0 and 1;
	 * 
	 * @param da
	 *            double[]: The array to be standardized
	 * @return double[]: The standardized array.
	 */

	public static double[] norm1(double[] da) {

		double sum = 0;

		for (double d : da) {

			sum += d;

		}

		if (sum == 0) {
			return new double[da.length];
		}

		DoubleMatrix1D x = new DenseDoubleMatrix1D(da);
		x.assign(cern.jet.math.Functions.div(sum));

		return x.toArray();

	}

	/**
	 * Standardizes an array of integers to values between 0 and 1;
	 * 
	 * @param intArray
	 *            double[]: The array to be standardized
	 * @return double[]: The standardized array.
	 */

	public static double[] norm1(int[] intArray) {

		double[] da = new double[intArray.length];

		for (int i = 0; i < intArray.length; i++) {

			da[i] = intArray[i];

		}

		return norm1(da);

	}

	/**
	 * Divides each matrix row by the sum of its elements, bounding the values
	 * by 0 and 1 (inclusive).
	 * 
	 * @param ia -
	 *            The matrix to be standardized.
	 * @return - The standardized matrix.
	 */

	public static double[][] norm1(int[][] ia) {

		double[][] out = new double[ia.length][ia[0].length];

		for (int i = 0; i < ia.length; i++) {

			out[i] = norm1(ia[i]);

		}

		return out;

	}

	/**
	 * Standardizes an array of longs to values between 0 and 1;
	 * 
	 * @param longArray
	 *            double[]: The array to be standardized
	 * @return double[]: The standardized array.
	 */

	public static double[] norm1(long[] longArray) {

		double[] da = new double[longArray.length];

		for (int i = 0; i < longArray.length; i++) {

			da[i] = longArray[i];

		}

		return norm1(da);

	}

	/**
	 * Divides each matrix row by the sum of its elements, bounding the values
	 * by 0 and 1 (inclusive).
	 * 
	 * @param ia -
	 *            The matrix to be standardized.
	 * @return - The standardized matrix.
	 */

	public static double[][] norm1(long[][] ia) {

		double[][] out = new double[ia.length][ia[0].length];

		for (int i = 0; i < ia.length; i++) {

			out[i] = norm1(ia[i]);

		}

		return out;

	}

	/**
	 * Creates a standardized matrix of random values.
	 * 
	 * @param dim
	 *            int: The dimension of the matrix to be created
	 * @return double[][]: The standardized random matrix.
	 */

	public static double[][] standardizedRandomMatrix(int dim) {

		double[][] matrix = new double[dim][dim];

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				matrix[i][j] = Uniform.staticNextDouble();

			}

		}

		return standardizeMatrix(matrix);

	}

	/**
	 * Performs the actual work of creating a matrix with equal distribution
	 * probabilities, i.e. all values are equal to 1 / the matrix dimension.
	 * 
	 * @param dim
	 *            int: The dimensions of the matrix to be created
	 * @return double[][]: The standardized uniform matrix.
	 */

	public static double[][] standardizedUniformMatrix(int dim) {

		double[][] matrix = new double[dim][dim];
		double val = 1.0d / dim;

		for (int i = 0; i < matrix.length; i++) {

			for (int j = 0; j < matrix[i].length; j++) {

				matrix[i][j] = val;

			}

		}

		return matrix;

	}

	/**
	 * Converts the transition matrix into a form where each array value is a
	 * probability between 0 and 1 inclusive, and the row elements sum to 1,
	 * subject to precision limits of the Java language.
	 * 
	 * @param matrix
	 *            double[][]: The matrix to be standardized
	 * @return double[][]: The standardized matrix
	 */

	public static double[][] standardizeMatrix(double[][] matrix) {

		double[][] stdMatrix = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {

			double sum = getRowSum(matrix, i);

			for (int j = 0; j < matrix[i].length; j++) {

				// If there are no entries in the row, we will assume the agent
				// stays where it is (mortality is handled elsewhere).

				if (sum == 0) {

					stdMatrix[i][i] = 1d;

					continue;
				}

				stdMatrix[i][j] = matrix[i][j] / sum;

			}
		}

		return stdMatrix;
	}

	/**
	 * 
	 * Converts the transition matrix into a form where each array value is a
	 * probability between 0 and 1 inclusive, and the row elements sum to 1,
	 * subject to precision limits of the Java language.
	 * 
	 * @param matrix
	 *            long[][]: The matrix to be standardized
	 * @return double[][]: The standardized matrix
	 */

	public static double[][] standardizeMatrix(long[][] matrix) {

		double[][] stdMatrix = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {

			double sum = getRowSum(matrix, i);

			for (int j = 0; j < matrix[i].length; j++) {

				// If there are no entries in the row, we will assume the agent
				// stays where it is (mortality is handled elsewhere).

				if (sum == 0) {

					stdMatrix[i][i] = 1d;

					continue;
				}

				stdMatrix[i][j] = matrix[i][j] / sum;

			}
		}

		return stdMatrix;

	}

	/**
	 * Subtracts two vectors of equal dimension.
	 * 
	 * @param a
	 *            long[] - The first matrix
	 * @param b
	 *            int[] - The second matrix
	 * @return long[] - The result of adding the two matrices.
	 */

	public static long[] subtract(long[] a, int[] b) {

		long[] sum = new long[a.length];

		for (int i = 0; i < a.length; i++) {

			sum[i] = a[i] - b[i];
		}

		return sum;
	}

	/**
	 * Subtracts two vectors of equal dimension.
	 * 
	 * @param a
	 *            long[] - The first matrix
	 * @param b
	 *            long[] - The second matrix
	 * @return long[] - The result of adding the two matrices.
	 */

	public static long[][] subtract(long[][] a, long[][] b) {

		long[][] difference = new long[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				difference[i][j] = a[i][j] - b[i][j];
			}
		}
		return difference;
	}
	
	/**
	 * Subtracts two vectors of equal dimension.
	 * 
	 * @param a
	 *            long[] - The first matrix
	 * @param b
	 *            long[] - The second matrix
	 * @return long[] - The result of adding the two matrices.
	 */

	public static float[] subtract(float[] a, float[] b) {

		float[] difference = new float[a.length];

		for (int i = 0; i < a.length; i++) {
				difference[i] = a[i] - b[i];
		}
		return difference;
	}
	
	/**
	 * Subtracts two vectors of equal dimension.
	 * 
	 * @param a
	 *            long[] - The first matrix
	 * @param b
	 *            long[] - The second matrix
	 * @return long[] - The result of adding the two matrices.
	 */

	public static double[] subtract(double[] a, double[] b) {

		double[] difference = new double[a.length];

		for (int i = 0; i < a.length; i++) {
				difference[i] = a[i] - b[i];
		}
		return difference;
	}

	/**
	 * Returns the sum of the supplied vector as a double
	 * 
	 * @param doubleArray
	 *            long[] - A vector of longs
	 * @return double - The sum of the vector.
	 */

	public static double sumAsDouble(double[] doubleArray) {
		double sum = 0;

		for (double d : doubleArray) {
			sum += d;
		}

		return sum;

	}

	/**
	 * Returns the sum of the supplied vector as a double
	 * 
	 * @param longArray
	 *            long[] - A vector of longs
	 * @return double - The sum of the vector.
	 */

	public static double sumAsDouble(long[] longArray) {
		double sum = 0;

		for (double d : longArray) {
			sum += d;
		}

		return sum;

	}

	/**
	 * Returns the sum of the supplied vector as a float
	 * 
	 * @param intArray
	 *            int[] - A vector of longs
	 * @return float - The sum of the vector.
	 */

	public static float sumAsFloat(int[] intArray) {
		float sum = 0;

		for (float f : intArray) {
			sum += f;
		}

		return sum;

	}

	/**
	 * Returns the sum of the supplied vector as a float
	 * 
	 * @param intArray
	 *            long[] - A vector of longs
	 * @return float - The sum of the vector.
	 */

	public static float sumAsFloat(long[] intArray) {
		float sum = 0;

		for (float f : intArray) {
			sum += f;
		}

		return sum;

	}

	/**
	 * Returns the sum of the supplied vector as an int
	 * 
	 * @param intArray
	 *            long[] - A vector of ints
	 * @return int - The sum of the vector.
	 */

	public static int sumAsInt(int[] intArray) {
		int sum = 0;

		for (int i : intArray) {
			sum += i;
		}

		return sum;

	}

	/**
	 * Returns the sum of the supplied vector as a long
	 * 
	 * @param longArray
	 *            long[] - A vector of longs
	 * @return long - The sum of the vector.
	 */

	public static long sumAsLong(long[] longArray) {
		long sum = 0;

		for (long l : longArray) {
			sum += l;
		}

		return sum;

	}

	/**
	 * Returns the transpose of the matrix
	 * 
	 * @param a -
	 *            The input matrix
	 * @return - The transpose of the matrix
	 */

	public static double[][] transpose(double[][] a) {

		double[][] out = new double[a[0].length][a.length];

		for (int i = 0; i < out.length; i++) {
			for (int j = 0; j < out[0].length; j++) {

				out[i][j] = a[j][i];

			}
		}
		return out;
	}

	/**
	 * Returns the transpose of the matrix
	 * 
	 * @param a -
	 *            The input matrix
	 * @return - The transpose of the matrix
	 */

	public static int[][] transpose(int[][] a) {

		int[][] out = new int[a[0].length][a.length];

		for (int i = 0; i < out.length; i++) {
			for (int j = 0; j < out[0].length; j++) {

				out[i][j] = a[j][i];

			}
		}
		return out;
	}

	/**
	 * Returns the transpose of the matrix
	 * 
	 * @param a -
	 *            The input matrix
	 * @return - The transpose of the matrix
	 */

	public static long[][] transpose(long[][] a) {

		long[][] out = new long[a[0].length][a.length];

		for (int i = 0; i < out.length; i++) {
			for (int j = 0; j < out[0].length; j++) {

				out[i][j] = a[j][i];

			}
		}
		return out;
	}

	/**
	 * Creates a matrix having uniformly distributed random values bounded by
	 * the parameters provided.
	 * 
	 * @param dim
	 *            int: The dimensions of the output matrix
	 * @param from
	 *            double: Lower boundary of possible values
	 * @param to
	 *            double: Upper boundary of possible values
	 * @return double[][]: The output matrix
	 */

	public static double[][] uniformRandomMatrixFromTo(int dim, double from, double to) {

		double[][] matrix = new double[dim][dim];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {

				matrix[i][j] = Uniform.staticNextDoubleFromTo(from, to);

			}

		}

		return matrix;

	}

	/**
	 * Convenience method for writing a double[][] matrix to a file.
	 * 
	 * @param matrix
	 *            double[][]: The matrix to be written to the file
	 * @param destination
	 *            String: The path and name of the destination file/
	 */

	public static void writeToFile(double[][] matrix, String destination) {

		try {
			PrintWriter pw = new PrintWriter(destination);
			StringBuffer strbf = new StringBuffer();

			for (int i = 0; i < matrix.length; i++) {

				for (int j = 0; j < matrix[i].length; j++) {

					strbf.append(matrix[i][j] + ", ");

				}

				strbf.setLength(strbf.length() - 2);
				strbf.append("\n");

			}

			pw.write(strbf.toString());

		}

		catch (FileNotFoundException ex) {
			ex.printStackTrace();
			System.out.println("Matrix output destination: " + destination
					+ "was not found.");
		}
	}
	
	public static int circular_search(float[] a, float x){
	    //instead of using the division op. (which surprisingly fails on big numbers)
	    //we will use the unsigned right shift to get the average
		int low = 0;
		int high = a.length-1;
		int mid = (low + high) >>> 1;
	    if(a[mid] == x){
	        return mid;
	    }
	    //a variable to indicate which half is sorted
	    //1 for left, 2 for right
	    int sortedHalf = 0;
	    if(a[low] <= a[mid]){
	        //the left half is sorted
	        sortedHalf = 1;
	        if(x <= a[mid] && x >= a[low]){
	            //the element is in this half
	            return Arrays.binarySearch(a, low, mid, x);
	        }
	    }
	    if(a[mid] <= a[high]){
	        //the right half is sorted
	        sortedHalf = 2;
	        if(x >= a[mid] && x<= a[high] ){
	            return Arrays.binarySearch(a, mid, high, x);
	        }
	    }
	    // repeat the process on the unsorted half
	    if(sortedHalf == 1){
	        //left is sorted, repeat the process on the right one
	        return circular_search(a, x);
	    }else{
	        //right is sorted, repeat the process on the left
	        return circular_search(a, x);
	    }
	}
	
	public static double[] lonlat2ceqd(double[] coords){
		double [] out = new double[coords.length];
		out[1] = Math.toRadians(coords[1])*R_EARTH; //lat
		out[0] = Math.cos(Math.toRadians(coords[1]))*Math.toRadians(coords[0])*R_EARTH; //lon
		if(coords.length==3){out[2]=coords[2];}
		return out;
	}

	// Cylindrical Equidistant Meters to longitude and latitude
	
	public static double[] ceqd2lonlat(double[] coords){
		double [] out = new double[coords.length];
		out[1] = Math.toDegrees((coords[1]/R_EARTH)); //lat
		out[0] = Math.toDegrees((coords[0]/R_EARTH)*1/(Math.cos(coords[1]))); //lon
		if(coords.length==3){out[2]=coords[2];}
		return out;
	}
}
