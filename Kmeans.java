import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeMap;

public class Kmeans {
	static ArrayList<Integer> ground_truth=new ArrayList<Integer>();
	static List<ArrayList<Double>> expr_values=new ArrayList<ArrayList<Double>>();
	static List<ArrayList<Double>> centroids=new ArrayList<ArrayList<Double>>();
	static HashMap<ArrayList<Double>, ArrayList<ArrayList<Double>>> last = new HashMap<ArrayList<Double>, ArrayList<ArrayList<Double>>>(); 
	//static ArrayList<Double> previous=new ArrayList<Double>();
	static List<ArrayList<Double>> previous;
	static TreeMap<Integer, ArrayList<Double>> labels=new TreeMap<Integer,ArrayList<Double>>();
	static HashMap<ArrayList<Double>,Integer> labels1=new HashMap<ArrayList<Double>,Integer>();
	static ArrayList<Integer> Results=new ArrayList<Integer>();
	static int number=6;
	static int count=0;
	static int size=0;
	static int iterations=0;
	public static void main(String[] args) throws IOException {
		System.out.println("Enter the number of clusters.");
		Scanner in = new Scanner(System.in);
		number = in.nextInt();
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
		    	size=columns.length-2;
		    
		    	for (int i = 2; i < columns.length; i++) {
					gene.add(Double.parseDouble(columns[i]));
					
					
				}
		    
		    	expr_values.add(gene);
		    
		    }
		    
		} finally {
		    br.close();
		}
		
		Random r = new Random();
		int i1; 
		while(centroids.size()<number) {
		for (int j = 0; j < number; j++) {
			i1 = r.nextInt(count);
			if (!centroids.contains(expr_values.get(i1))) {
				centroids.add(expr_values.get(i1));
			}
		}
		}
		
		
		
		kmeans_cluster();
		
		
		
		while(!check(previous)) {
			//System.out.println("enter");
		previous=new ArrayList<ArrayList<Double>>();
		previous.addAll(centroids);
		update_centroids();
		//System.out.println("next run");
		kmeans_cluster();
		
		
		}
		ArrayList<Double> inverse;
		int i=1;
		for(int j:labels.keySet()){
			inverse=labels.get(j);
			if(!labels1.containsKey(inverse)){
				labels1.put(inverse, i);
				i++;
			}
		}
		
		
		Writer w;
		int z=1;
		w=new FileWriter("/home/nischala/workspace/proj2/src/result.txt");
		for(int k: labels.keySet()){
			Results.add(labels1.get(labels.get(k)));
			w.write(String.valueOf(labels1.get(labels.get(k))));
			w.write("\r\n");
			System.out.println(k+" : "+labels1.get(labels.get(k)));
			z++;
		}
		//System.out.println(count);
		w.close();
		int[][] gmatrix=new int[count][count];
		int[][] cmatrix=new int[count][count];
		for (int u=0;u<count;u++) {
            for (int j=u;j<count;j++) {
            	if (ground_truth.get(u)==ground_truth.get(j))
                    gmatrix[u][j]=gmatrix[j][u]=1;
                else
                	gmatrix[u][j]=gmatrix[j][u]=0;
                if (Results.get(u)==Results.get(j))
                    cmatrix[u][j]=cmatrix[j][u]=1;
                else
                    cmatrix[u][j]=cmatrix[j][u]=0;
                }
        }
		//System.out.println(gmatrix);
		double jaccard=Jaccard.jaccard_coeff(gmatrix,cmatrix,count);
		System.out.println("External Validation: Jaccard Coefficient is "+jaccard);
		//test
		//System.out.println(Results.get(0));
		HashMap<Integer, ArrayList<Integer>> cmap = new HashMap<Integer, ArrayList<Integer>>();
		for(int c = 0;c < Results.size(); c++){
			ArrayList<Integer> genesInCluster = null;
			int cluster = Results.get(c);
			//System.out.println(Results.get(0));
			
			if(cmap.get(cluster) == null){
				genesInCluster = new ArrayList<Integer>();
			}else{
				genesInCluster = cmap.get(cluster);
			}
			
			genesInCluster.add(c);
			
			cmap.put(cluster, genesInCluster);
		}
		/*for(int t:cmap.keySet())
			System.out.println("key: "+t+" value "+cmap.get(t));*/
		//distance matrix
		//System.out.println(expr_values);
		double distance_matrix[][]=new double[count][count];
		for(int d = 0; d < count; d++){
			ArrayList<Double> a = expr_values.get(d);
			
			for(int j = 0; j < count; j++){
				ArrayList<Double> b = expr_values.get(j);
				distance_matrix[d][j] = dist(a,b);
				//System.out.println(distance_matrix[d][j]);
			}
		}
		/*for(int d = 0; d < count; d++){
			for(int j = 0; j < count; j++){
				System.out.println(distance_matrix[d][j]);
			}*/
			
		//}
		Silhouette s=new Silhouette(distance_matrix,cmap,Results,count);
		double silhouette=s.calculateSilhouette();
		System.out.println("Silhouette coefficient is: "+silhouette);
		
		
	}
	public static boolean check(List<ArrayList<Double>> list) {
		ArrayList<ArrayList<Double>> list1=new ArrayList<ArrayList<Double>>();
		boolean y=false;
		//if(iterations>50)
			//y=true;
		//else
			//y=false;
		list1.addAll(centroids);
		
		if(!(list==null)) {
		if(list.equals(list1)) {
			
			y= true;
		}
		else {
			y= false;
			
		}
		}
		return y;
		
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
	public static void kmeans_cluster() {
		ArrayList<Double> list=new ArrayList<Double>();
		HashMap<Double, Integer> map=new HashMap<Double, Integer>();
		last.clear();
		for(int i=0;i<expr_values.size();i++) {
			list.clear();
            map.clear();
			ArrayList<Double> tempData = expr_values.get(i);
			for(int j = 0; j < number; j++)
            {
				
				
                
				double distance = dist(tempData, centroids.get(j));
                
                
               
                map.put(distance, j);
                list.add(distance);
            }//j
			Collections.sort(list);
			
			int temp;
			temp = map.get(list.get(0));
			labels.put(i, centroids.get(temp));
			//labels1.put(i, centroids.get(temp));
			
			
			//System.out.println("Tree "+labels);
			if (!last.containsKey(centroids.get(temp))) {
				ArrayList<ArrayList<Double>> entities = new ArrayList<ArrayList<Double>>();
				entities.add(expr_values.get(i));
				last.put(centroids.get(temp), entities);
			} else {
				ArrayList<ArrayList<Double>> entities = last.get(centroids.get(temp));
				entities.add(expr_values.get(i));
				last.put(centroids.get(temp), entities);
			}
		 
		
		}//i
		//System.out.println("1:"+centroids);
		
	}//kmeans_cluster
	public static void update_centroids() {
		centroids.clear();
		double[] sum=new double[size];
	    double[] mean=new double[size];
	    
	    double totalInCluster;
		for(ArrayList<Double> s: last.keySet()) {
			for(int i=0;i<sum.length;i++) {
        		sum[i]=0;
        		mean[i]=0;
        	}
        	//totalX = 0;
            totalInCluster = 0;
            
            for(ArrayList<Double> v:last.get(s)) {
            	
            	for(int j=0;j<v.size();j++) {
            		
            		sum[j]=sum[j]+v.get(j);
            		
            		
            	}
            	totalInCluster++;
            	
            	
             } 
            ArrayList<Double> m= new ArrayList<Double>();
            for(int i=0;i<sum.length;i++) {
            	mean[i]=sum[i]/totalInCluster;
            	m.add(mean[i]);
            	//System.out.println(m);
            	
            }
            
            	centroids.add(m);
            	
            
        }//s
		
		
	}//update-centroids
	
}//class
