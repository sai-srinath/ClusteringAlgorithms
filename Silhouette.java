import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

public class Silhouette {
	
	double[][] distance_matrix;
	HashMap<Integer, ArrayList<Integer>> cmap;
	ArrayList<Integer> Results;
	int count=0;
	public Silhouette(double[][] matrix,HashMap<Integer,ArrayList<Integer>> map,ArrayList<Integer> a, int number) {
		this.distance_matrix=matrix;
		this.cmap=map;
		this.Results=a;
		this.count=number;
		
	}
public double calculateSilhouette(){
		
		double silhouette=0;
		double B_sum=0;
		double A=0;
		double B=0;
		int clusterID=0;
		ArrayList<Integer> genesInclusterID=null;
		for(int i=0;i< Results.size();i++){
			clusterID=Results.get(i);
			
			A=0;
			genesInclusterID=cmap.get(clusterID);

			for(int j=0;j<genesInclusterID.size();j++){
				int geneID=genesInclusterID.get(j);
				A+=distance_matrix[i][geneID];
			}
			if(genesInclusterID.size() > 1)
				A=A/(genesInclusterID.size() - 1);
			
			B=0;
			Set<Integer> clusterIDKeys=cmap.keySet();
			Iterator<Integer> clusterIDKeys_irtr=clusterIDKeys.iterator();
			while(clusterIDKeys_irtr.hasNext()){
				int other_clusterID=clusterIDKeys_irtr.next();
				B_sum=0;
				if(other_clusterID == clusterID) continue;
				ArrayList<Integer> genesInOtherclusterID=cmap.get(other_clusterID);
				for(int j=0;j<genesInOtherclusterID.size();j++){
					int geneID=genesInOtherclusterID.get(j);
					B_sum+=distance_matrix[i][geneID];
				}
				B_sum=B_sum/genesInOtherclusterID.size();

				if(B==0)
					B=B_sum;
				else
					B=Math.min(B,B_sum);
				}
			
			double S=(B-A)/Math.max(A,B);
			silhouette=silhouette + S;
			}
		return silhouette/count;
	}

}
