public class Solution {
    public static int solution(int[] l) {
        int luckytriplecount = 0; // counts the number of lucky triples
        for (int i = 0; i < l.length - 2; ++i) { // iteratively finds these triples; O(n^3) time
            for (int j = i + 1; j < l.length - 1; ++j) {
                for (int k = j + 1; k < l.length; ++k) {
                    if ((l[j] % l[i] == 0) && (l[k] % l[j] == 0)) {
                        ++luckytriplecount;
                    }
                }
            }
        }
        return luckytriplecount;
    }
}