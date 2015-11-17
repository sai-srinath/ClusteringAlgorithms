
public class Jaccard {
	
	public static double jaccard_coeff(int[][] gmatrix,int[][] cmatrix,int size) {
		
		double jaccard=0;
		int one=0;
		int one_zero=0;
		int zero_one=0;
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++) {
				if(cmatrix[i][j]==gmatrix[i][j]) {
					if(cmatrix[i][j]==1)
						one++;
					}
				else {
					if(cmatrix[i][j]==1 && gmatrix[i][j]==0)
						one_zero++;
					else if(cmatrix[i][j]==0 && gmatrix[i][j]==1)
						zero_one++;
				}
			}
		}
		jaccard=(one)/(double)(one + one_zero + zero_one);
		return jaccard;
	}
}
