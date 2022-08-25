import java.util.ArrayList; // import the ArrayList class

// Had to learn about (Absorbing) Markov Chains and do a lot of research in order to complete this problem.
// Watched this video to help me figure this out https://www.youtube.com/watch?v=qhnFHnLkrfA

public class Solution {
    // HELPER FUNCTIONS //
    static private int[] decimaltofraction(double x) { // takes a decimal and converts it into a fraction array representation, with the 0 index being the numerator and the 1 index being the denominator, using the algorithm in https://begriffs.com/pdf/dec2frac.pdf
        double tolerancelevel = 1E-10;
        double p1 = 1;
        double q1 = 0;
        double p2 = 0;
        double q2 = 1;
        double b = x;
        while (Math.abs(x - (p1 / p2)) > (x * tolerancelevel)) { // while the difference between the number and the fraction is bigger than the tolerance level times the fraction
            double a = Math.floor(b);
            double temp = p1;
            p1 = (a * p1) + q1; // computes the next iteration numerator
            q1 = temp;
            temp = p2;
            p2 = (a * p2) + q2; // computes the next iteration denominator
            q2 = temp;
            b = 1 / (b - a); // the fraction oscillates below and above x
        }
        return new int[] {
            (int) p1, (int) p2
        }; // return the fraction array
    }

    public static int gcd(int a, int b) // finds the gcd of a and b
    { // using the euclidean algorithm
        if (b == 0) {
            return a;
        }
        return gcd(b, a % b);
    }

    public static int gcdmatrix(int[][] m) { // finds the greatest common divisior of the integers in the matrix, using gcd(a, b)
        int value = 0;
        int length = m.length;
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {

                value = gcd(value, m[i][j]);
            }
        }
        return value;
    }

    public static int gcdarray(int[] array) { // finds the greatest common divisior of the integers in the array, using gcd(a, b)
        int value = 0;
        for (int i = 0; i < array.length; ++i) {

            value = gcd(value, array[i]);

        }
        return value;
    }

    // inverse of matrix algorithm created with help from https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
    public static double determinant(double[][] m) { // recursively produces the determinant of a matrix m, using the cofactor expansion
        int length = m.length;
        double finalvalue = 0;
        if (length == 1) {
            return m[0][0];
        } else {
            for (int j = 0; j < length; ++j) {
                // finds the sub matrices to compute on
                double[][] matrix = new double[length - 1][length - 1];
                int c = 0;
                for (int a = 1; a < length; ++a) {
                    int d = 0;
                    for (int b = 0; b < length; ++b) {
                        if (b != j) {
                            matrix[c][d] = m[a][b];
                            ++d;
                        }
                    }
                    ++c;
                }
                if ((j % 2) == 0) { // if the indexes of the rows and columns sum to even, then add
                    finalvalue += (m[0][j] * determinant(matrix));
                } else { // otherwise subtract
                    finalvalue -= (m[0][j] * determinant(matrix));
                }
            }
        }
        return finalvalue;
    }

    public static double[][] inversematrix(double[][] m) { // produces the inverse of matrix m, using minors, cofactors, and the adjugate (a lot of computation)
        // first we create the matrix of minors
        int length = m.length;
        double[][] matrix = new double[length][length];
        // singleton matrix base case
        if (length == 1) {
            matrix[0][0] = (1 / m[0][0]);
            return matrix;

        }
        for (int i = 0; i < length; ++i) { // iteratively produces the matrix of minors usning the determinants
            for (int j = 0; j < length; ++j) {
                // finds the sub matrices to compute on
                double[][] submatrix = new double[length - 1][length - 1];
                int c = 0;
                for (int a = 0; a < length; ++a) {
                    if (a != i) {
                        int d = 0;
                        for (int b = 0; b < length; ++b) {
                            if (b != j) {
                                submatrix[c][d] = m[a][b];
                                ++d;
                            }
                        }
                        ++c;
                    }
                }
                matrix[i][j] = determinant(submatrix);
            }
        }
        // next, we create the matrix of cofactors by applying negative symbols where necessary
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                if (((i + j) % 2) == 0) { // if sum of indexes are even, do nothing
                } else { // otherwise, negate the values
                    matrix[i][j] *= -1;
                }
            }
        }
        // next, we adjugate, or swap the positions of the matrix over the diagonalizer
        for (int i = 0; i < length; ++i) {
            for (int j = i; j < length; ++j) {

                if (j > i) { // if it's upper diagonal part, swap values, otherwise do nothing
                    double temp = matrix[j][i];
                    matrix[j][i] = matrix[i][j];
                    matrix[i][j] = temp;
                }
            }
        }
        // finally, multiply the adjugate matrix by 1 / determinant of the original matrix
        double detoriginal = determinant(m);
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                matrix[i][j] = (matrix[i][j] / detoriginal);
            }
        }
        return matrix;
    }

    public static double[][] standardizer(double[][] m) { // produces the standard form of the matrix m (helper function)
        int numstates = m.length;
        double[] denomarray = new double[numstates]; // denominator array
        for (int i = 0; i < numstates; ++i) { // iteratively initializes the denominator for each state and places it in the array
            double sum = 0;
            for (int j = 0; j < numstates; ++j) {
                sum += m[i][j];

            }
            denomarray[i] = sum;
        }
        double[][] matrix = new double[numstates][numstates];
        int n = 0;
        int x = 0;
        for (int i = 0; i < numstates; ++i) { // first finds absorbing states (states that are terminals OR states with singleton loops), and makes the terminals singleton "1" loops and places them in the beginning rows
            if (denomarray[i] == 0) { // if all 0s (terminating)
                matrix[n] = m[i];
                matrix[n][i] = 1;
                ++n;
                ++x;
                continue;
            }
            boolean singletonloop = true;
            for (int j = 0; j < numstates; ++j) { // loops that checks if state has only singleton loops
                if ((m[i][j] != 0) && (i != j)) {
                    singletonloop = false;
                    break;
                }
            }
            if (singletonloop) {
                matrix[n] = m[i];
                matrix[n][i] = 1;
                ++n;
                ++x;
                continue;
            }
        }
        for (int i = 0; i < numstates; ++i) { // finds the rest of the states and inserts them
            if (denomarray[i] == 0) { // if all 0s (terminating)
                continue;
            }
            boolean singletonloop = true;
            for (int j = 0; j < numstates; ++j) { // loops that checks if state has only singleton loops
                if ((m[i][j] != 0) && (i != j)) {
                    singletonloop = false;
                    break;
                }
            }
            if (singletonloop) {
                continue;
            }
            matrix[n] = m[i];
            ++n;
        }
        double[][] matrixrefined = new double[numstates][numstates];
        int a = 0;
        int b = 0;
        ArrayList < Integer > temparray = new ArrayList < Integer > (); // Create an ArrayList object
        // now, we clean the matrix up to produce a standard form matrix
        int diagonalizer = 0;
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < numstates; ++j) {
                if (matrix[i][j] == 1) {
                    for (int k = 0; k < numstates; ++k) { // swap the columns
                        matrixrefined[k][diagonalizer] = matrix[k][j];
                    }
                    temparray.add(j);

                    ++diagonalizer;
                }
            }
        }
        for (int j = 0; j < numstates; ++j) {
            boolean temp = false;
            for (int k = 0; k < temparray.size(); ++k) {
                if (j == temparray.get(k)) {
                    temp = true;
                }
            }
            if (!temp) {
                for (int s = 0; s < numstates; ++s) { // get the rest of the columns
                    matrixrefined[s][diagonalizer] = matrix[s][j];
                }
                ++diagonalizer;
            }
        }
        for (int i = 0; i < numstates; ++i) { // iteratively initializes the denominator again for each state and places it in the array
            double sum = 0;
            for (int j = 0; j < numstates; ++j) {
                sum += matrixrefined[i][j];
            }
            denomarray[i] = sum;
        }
        for (int i = 0; i < numstates; ++i) {
            for (int j = 0; j < numstates; ++j) {
                matrixrefined[i][j] = (matrixrefined[i][j] / denomarray[i]);
            }
        }
        return matrixrefined;
    }

    public static double[][] getRmatrix(double[][] m) { // produces the R matrix of the standard form matrix (bottom left square matrix partitioned from identity)
        int n = 0;
        int length = m.length;
        for (int i = 0; i < length; ++i) { // row cutoff algorithm
            boolean allzeroes = true;
            for (int j = 0; j < length; ++j) { // checks if identity matrix
                if ((m[i][j] != 0) && (i != j)) {
                    allzeroes = false;
                    break;
                }
            }
            if (allzeroes) {
                ++n;
            }
        }

        double[][] matrix = new double[length - n][n]; // initialize the matrix
        int a = 0; // to build the matrix
        for (int i = n; i < length; i++) {
            for (int j = 0; j < n; j++) {
                matrix[a][j] = m[i][j]; // get the values
            }
            ++a;
        }
        return matrix;
    }

    public static double[][] getQmatrix(double[][] m) { // produces the Q matrix of the standard form matrix (bottom right square matrix partitioned from identity and R matrix)
        // note that this code is very similar to the code in getRmatrix
        int n = 0;
        int length = m.length;
        for (int i = 0; i < length; ++i) { // row cutoff algorithm
            boolean allzeroes = true;
            for (int j = 0; j < length; ++j) { // checks if identity matrix
                if ((m[i][j] != 0) && (i != j)) {
                    allzeroes = false;
                    break;
                }
            }
            if (allzeroes) {
                ++n;
            }
        }
        double[][] matrix = new double[length - n][length - n]; // initialize the matrix
        int a = 0; // to build the rows of the matrix
        for (int i = n; i < length; i++) {
            int b = 0; // to build the columns of the matrix
            for (int j = n; j < length; j++) {
                matrix[a][b] = m[i][j]; // get the values
                ++b;
            }
            ++a;
        }
        return matrix;
    }
    public static double[][] getfundamentalmatrix(double[][] Q) { // produces the fundamental matrix of the standard form matrix, using its Q matrix. It is the inverse of (Identity matrix - Q matrix)
        int length = Q.length;
        double[][] matrix = new double[length][length]; // initialize the matrix
        for (int i = 0; i < length; ++i) { // finds the (I - Q) matrix
            for (int j = 0; j < length; ++j) {
                if (i == j) {
                    matrix[i][j] = (1 - Q[i][j]);
                } else {
                    matrix[i][j] = -Q[i][j];
                }
            }
        }
        // now we need to find the inverse of it...
        // in order to do so, we use the helper function (long code)
        matrix = inversematrix(matrix);
        return matrix;
    }

    public static double[][] getFRmatrix(double[][] f, double[][] R) { // produces the fundamental matrix of the standard form matrix, multiplied by the R matrix, FR.
        int m = f.length; // fundamental matrix will have m x n matrix
        int n = f.length; // it's a square matrix so
        int o = R.length; // R matrix will have o x p matrix NOT NECESSARILY A SQUARE MATRIX
        int p = R[0].length; // width
        double[][] matrix = new double[m][p]; // initialize the matrix (m x n matrix times n x p matrix equals m x p matrix)
        for (int i = 0; i < m; ++i) { // for each index in the new matrix,
            for (int j = 0; j < p; ++j) {
                double value = 0;
                // take the sum of the products of each ith index in the row in f with its respective ith index in the column R
                for (int k = 0; k < n; ++k) {
                    value += (f[i][k] * R[k][j]);
                }
                matrix[i][j] = value;
            }
        }
        return matrix;
    }

    public static void getlimitingmatrix(double[][] standardizedmatrix,
        double[][] FR) { // changes the standardized matrix doubleo the limiting matrix using the FR matrix, by turning all elements in the bottom right of the identity matrix 0s,
        // and replaceing the R part with FR.
        int length = standardizedmatrix.length;
        // note that this code is very similar to the code in getRmatrix
        int n = 0;
        for (int i = 0; i < length; ++i) { // row cutoff algorithm
            boolean allzeroes = true;
            for (int j = 0; j < length; ++j) { // checks if identity matrix
                if ((standardizedmatrix[i][j] != 0) && (i != j)) {
                    allzeroes = false;
                    break;
                }
            }
            if (allzeroes) {
                ++n;
            }
        }
        int a = 0;
        for (int i = n; i < length; i++) {
            int b = 0;
            for (int j = 0; j < length; j++) {
                if (j >= n) {
                    standardizedmatrix[i][j] = 0;
                } else {
                    standardizedmatrix[i][j] = FR[a][b];
                    ++b;
                }
            }
            ++a;
        }
    }

    public static int[] solution(int[][] m) {
        // base case singleton sets
        if (m.length == 1) {
            int solution[] = new int[m.length + 1];
            solution[0] = 1;
            solution[1] = 1;
            return solution;
        }
        double[] denomarray = new double[m.length]; // denominator array for ease of life
        for (int i = 0; i < m.length; ++i) { // iteratively initializes the denominator for each state and places it in the array
            double sum = 0;
            for (int j = 0; j < m.length; ++j) {
                sum += m[i][j];
            }
            denomarray[i] = sum;
        }
        // 0 denominator in first array base case
        if (denomarray[0] == 0) {
            int numterminals = 0;
            for (int i = 0; i < m.length; ++i) {
                if (denomarray[i] == 0) {
                    ++numterminals;
                }
            }
            int solution[] = new int[numterminals + 1];
            for (int i = 0; i < numterminals; ++i) {
                if (i == 0) {
                    solution[i] = 1;
                } else {
                    solution[i] = 0;
                }
            }
            solution[numterminals] = 1;
            return solution;
        }
        int rows = m.length;
        int columns = m[0].length;
        double[][] doublematrix = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                doublematrix[i][j] = m[i][j];
            }
        }
        double tempmatrix[][] = standardizer(doublematrix);
        double Rmatrix[][] = getRmatrix(tempmatrix);
        double Qmatrix[][] = getQmatrix(tempmatrix);
        double fmatrix[][] = getfundamentalmatrix(Qmatrix);
        double FR[][] = getFRmatrix(fmatrix, Rmatrix);
        getlimitingmatrix(tempmatrix, FR); // <- turns the standardized matrix into a limiting matrix
        int n = 0;
        for (int i = 0; i < tempmatrix.length; ++i) { // row cutoff algorithm
            boolean allzeroes = true;
            for (int j = 0; j < tempmatrix.length; ++j) { // checks if identity matrix
                if ((tempmatrix[i][j] != 0) && (i != j)) {
                    allzeroes = false;
                    break;
                }
                if ((i == j) && (tempmatrix[i][j] == 0)) {
                    allzeroes = false;
                    break;
                }
            }
            if (allzeroes) {
                ++n;
            }
        }
        int[] numarray = new int[n];
        int[] denominarray = new int[n];
        for (int j = 0; j < n; ++j) { // place numerator and denominator into its respective index arrays
            int[] temp = decimaltofraction(tempmatrix[n][j]);
            numarray[j] = temp[0];
            denominarray[j] = temp[1];
        }
        //loop through the array to find GCD
        //use GCD to find the LCM
        int lcmfinal = denominarray[0];
        int gcdfinal = denominarray[0];
        for (int i = 1; i < denominarray.length; i++) {
            gcdfinal = gcd(denominarray[i], lcmfinal);
            lcmfinal = (lcmfinal * denominarray[i]) / gcdfinal;
        }
        for (int i = 0; i < numarray.length; ++i) // make all integers have the same denominator
        {
            numarray[i] *= (lcmfinal / denominarray[i]);
        }
        // finally, put the lcm at the end of the array
        int solution[] = new int[numarray.length + 1];
        for (int i = 0; i < numarray.length; ++i) {
            solution[i] = numarray[i];
        }
        solution[numarray.length] = lcmfinal;
        return solution;
    }
}