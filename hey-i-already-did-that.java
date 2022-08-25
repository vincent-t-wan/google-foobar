// Importing Arrays class from the utility class
import java.util.Arrays;
import java.util.ArrayList; // import the ArrayList class
public class Solution {
    public static int solution(String n, int b) {
        int i = 0;
        int k = n.length();
        int p = Integer.valueOf(n); // converts it into int
        ArrayList < String > v1 = new ArrayList < String > (); // Create an ArrayList object
        while (true) {
            String x;
            String y;
            // Converting input string to character array
            char tempArray[] = n.toCharArray();
            // Sorting temp array (ascending) using
            Arrays.sort(tempArray);
            y = String.valueOf(tempArray); // gives y the ascending string
            // Sorting temp array (descending) using
            String reversedString = "";
            for (int o = tempArray.length - 1; o >= 0; o--) {
                reversedString = reversedString + tempArray[o];
            }
            x = reversedString; // gives x the descending string
            char newn[] = new char[k]; // make new array
            int carry = 0;
            for (int j = k - 1; j >= 0; j--) { // subtracts each value individually 
                int current = (Character.getNumericValue(x.charAt(j)) - Character.getNumericValue(y.charAt(j))) - carry;
                if (current < 0) {
                    newn[j] = (char)(b + current + 48);
                    carry = 1;
                } else {
                    newn[j] = (char)(current + 48);
                    carry = 0;
                }
            }
            n = String.valueOf(newn); // replace
            v1.add(String.valueOf(newn)); // add new string num to arraylist v1
            if (i > 0) { // if the arraylist finds a match, then return the index difference.
                for (int a = v1.size() - 2; a >= 0; --a) {
                    if ((v1.get(a)).equals(v1.get(v1.size() - 1))) {
                        return (v1.size() - 1) - a;
                    }
                }
            }
            ++i;
        }
    }
}