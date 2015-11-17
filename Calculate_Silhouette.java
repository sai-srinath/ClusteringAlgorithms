import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;


public class Calculate_Silhouette {
	static List<ArrayList<Double>> expr_values=new ArrayList<ArrayList<Double>>();
	static ArrayList<Integer> ground_truth=new ArrayList<Integer>();
	static int count=0;
	private static BufferedReader reader;
	public static void main(String[] args) throws NumberFormatException, IOException {
		Scanner in = new Scanner(System.in);
		System.out.println("Enter name of the file.");
		
		String name=in.next();
		String path="/home/nischala/workspace/proj2/src/";
		String fname=path+name;
		BufferedReader br = new BufferedReader(new FileReader(fname));
		ArrayList<Double> gene=null;
		try {
		    
		    String line;

		    while ((line = br.readLine())!= null) {
		    	gene=new ArrayList<Double>();
		    	count++;
		    	String columns[]=line.split("\t");
		    	ground_truth.add(Integer.parseInt(columns[1]));
		    	
		    	for (int i = 2; i < columns.length; i++) {
					gene.add(Double.parseDouble(columns[i]));
					
					
				}
		    	
		    	expr_values.add(gene);
		    	
		    }
		    
		} finally {
		    br.close();
		}
		
		
		BufferedReader read = new BufferedReader(new FileReader("/home/nischala/workspace/CalculateSilhouette/src/hierarchichal_cho_target.txt"));
		ArrayList<Integer> Results=new ArrayList<Integer>();
		try {
		    
		    String line1;

		    while ((line1 = read.readLine())!= null) {
		    	
		    	
					Results.add(Integer.parseInt(line1));
					
					
				
		    	
		    }
		    
		} finally {
		    read.close();
		}
		
		HashMap<Integer, ArrayList<Integer>> cmap = new HashMap<Integer, ArrayList<Integer>>();
		for(int c = 0;c < Results.size(); c++){
			ArrayList<Integer> genesInCluster = null;
			int cluster = Results.get(c);
			
			
			
			if(cmap.get(cluster) == null){
				genesInCluster = new ArrayList<Integer>();
			}else{
				genesInCluster = cmap.get(cluster);
			}
			
			genesInCluster.add(c);
			
			cmap.put(cluster, genesInCluster);
		}
		double distance_matrix[][]=new double[count][count];
		for(int d = 0; d < count; d++){
			ArrayList<Double> a = expr_values.get(d);
			
			for(int j = 0; j < count; j++){
				ArrayList<Double> b = expr_values.get(j);
				distance_matrix[d][j] = dist(a,b);
				//System.out.println(distance_matrix[d][j]);
			}
		}
		double sil=silhouette(distance_matrix,cmap,Results,count);
		System.out.println("Silhouette coeff is "+sil);
	}
	private static double dist(ArrayList<Double> d, ArrayList<Double> c)
    {
	 	double distance=0.0;
    	for(int i=0;i<d.size();i++) {
    		double dist = (Math.pow((c.get(i) - d.get(i)), 2));
    		distance+=dist;
    	}
    		
    	return Math.sqrt(distance);
    }
	public static double silhouette(double[][] matrix,HashMap<Integer,ArrayList<Integer>> cmap,ArrayList<Integer> cresults, int number) {
		double silhouette=0;
		double B_sum=0;
		double A=0;
		double B=0;
		int clusterID=0;
		ArrayList<Integer> genesInclusterID=null;
		for(int i=0;i< cresults.size();i++){
			clusterID=cresults.get(i);
			
			A=0;
			genesInclusterID=cmap.get(clusterID);
			
			for(int j=0;j < genesInclusterID.size();j++){
				int geneID=genesInclusterID.get(j);
				A+=matrix[i][geneID];
			}
			if(genesInclusterID.size() > 1)
				A=A/(genesInclusterID.size() - 1);
			
			//System.out.println(A);
			
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
					B_sum+=matrix[i][geneID];
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
		return silhouette/number;
	}
	}


