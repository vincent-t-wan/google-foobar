// this code works if we are not dealing with extremely large numbers

import java.math.BigInteger; // to store large numbers
import java.util.ArrayList; // import the ArrayList class
public class Solution {
    public static int solution(String x) {
        ArrayList < Integer > bestarray = new ArrayList < Integer > (); // Create an ArrayList object
        bestarray.add(1);
        bestarray.add(0);
        int n = Integer.valueOf(x);
        if (n >= 2) {
            for (int i = 2; i < n + 2; ++i) { // minus and divide part to iterate the bestarray
                if ((i % 2) == 1) {
                    bestarray.add(bestarray.get(i - 1) + 1);
                } else {
                    bestarray.add(Math.min(bestarray.get(i - 1),
                        bestarray.get(i / 2)) + 1);
                }
            }
            for (int i = 2; i < n + 1; ++i) { // do it again because we need to add the plus part
                if ((i % 2) == 1) { // if odd
                    bestarray.set(i,
                        Math.min(bestarray.get(i - 1),
                            bestarray.get(i + 1)) + 1);
                } else {
                    bestarray.set(i,
                        Math.min(bestarray.get(i - 1),
                            Math.min(bestarray.get(i / 2),
                                bestarray.get(i + 1))) +
                        1);
                }
            }
        }
        return bestarray.get(n);
    }
}

// this code always works (neat algorithm)

import java.math.BigInteger; // to store large numbers
import java.util.ArrayList; // import the ArrayList class
public class Solution {
    public static int solution(String x) {
        ArrayList < Integer > bestarray = new ArrayList < Integer > (); // Create an ArrayList object
        bestarray.add(1);
        bestarray.add(0);
        BigInteger n = new BigInteger(x);
        int counter = 0;
        while (n.compareTo(new BigInteger("1")) != 0) { // loop until n = 1
            // base cases
            if (n.compareTo(new BigInteger("0")) == 0) {
                n = n.add(new BigInteger("1"));
                ++counter;
            } else if (n.compareTo(new BigInteger("2")) == 0) {
                n = n.subtract(new BigInteger("1"));
                ++counter;
            } else if (n.compareTo(new BigInteger("3")) == 0) {
                n = n.subtract(new BigInteger("1"));
                ++counter;
            } else { // recursive step
                BigInteger nmod4 = n.mod(new BigInteger("4"));
                if (nmod4.compareTo(new BigInteger("1")) == 0) { // if remainder is 1, then it is better to subtract 1, since adding 1 will only allow for 1 division by 2
                    n = n.subtract(new BigInteger("1"));
                } else if (nmod4.compareTo(new BigInteger("3")) == 0) { // if remainder is 3, then it is better to add 1, since subtracting 1 will only allow for 1 division by 2
                    n = n.add(new BigInteger("1"));
                } else { // otherwise, we are at even numbers, and so we divide by 2, since dividing by 2 will always be the fastest solution
                    n = n.divide(new BigInteger("2"));
                }
                ++counter;
            }
        }
        return counter;
    }
}