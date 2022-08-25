import java.math.BigInteger; // to store large numbers
public class Solution {
    public static String solution(int[] xs) {
        int numnegative = 0; // counts the number of negative elements (< 0).
        int minnegnumabsindex = -1; // the index containing the minimum negative number in absolute value
        BigInteger currentproduct = new BigInteger("1"); // stores the product as it goes through the algorithm
        boolean morethanonenegvalue = false; // boolean for more than one negative value in array
        boolean atleastonepositivevalue = false; // boolean for at least one positive value in array
        BigInteger currentnegativeproduct = new BigInteger("1"); // stores the product of the negative values as it iterates
        boolean atleastonezero = false; // boolean for at least one zero in array
        if (xs.length == 0) {
            currentproduct = new BigInteger("0");
            return currentproduct.toString();
        }
        for (int n = 0; n < xs.length; ++n) {
            if (xs[n] < 0) {
                int currentnegnumabsindex = n;
                if (minnegnumabsindex == -1) {
                    minnegnumabsindex = n;
                } else {
                    if ((-1 * xs[currentnegnumabsindex]) < (-1 * xs[minnegnumabsindex])) {
                        minnegnumabsindex = currentnegnumabsindex;
                    }
                }
                ++numnegative;
                currentnegativeproduct = currentnegativeproduct.multiply(new BigInteger(String.valueOf(xs[n])));
            } else if (xs[n] > 0) {
                atleastonepositivevalue = true;
                currentproduct = currentproduct.multiply(new BigInteger(String.valueOf(xs[n])));
            } else {
                atleastonezero = true;
            }
        }
        if (numnegative > 1) {
            morethanonenegvalue = true;
        }
        if ((!atleastonepositivevalue) && ((!morethanonenegvalue) && atleastonezero)) { // that means all elements are 0 except for possibly 1 negative value, so 0 is greatest
            currentproduct = new BigInteger("0");
            return currentproduct.toString();
        }
        if ((numnegative % 2) == 1) { // odd number of negative elements
            if (xs.length == 1) { // only one negative element
                String solution = String.valueOf(xs[0]);
                return solution;
            } else {
                currentnegativeproduct = currentnegativeproduct.divide(new BigInteger(String.valueOf(xs[minnegnumabsindex])));
                currentproduct = currentproduct.multiply(currentnegativeproduct);
                return currentproduct.toString();
            }
        } else { // even number of negative elements
            currentproduct = currentproduct.multiply(currentnegativeproduct);
            return currentproduct.toString();
        }
    }
}