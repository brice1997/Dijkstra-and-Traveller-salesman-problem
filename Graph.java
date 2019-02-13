/**Author of the student code: 
 * 
 * Alex Brice Horimbere
 * 
 * abh2167
 * December 13, 2018
 * 
 * */


import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
    
    // Your code here...
    // 
    
    Random m = new Random(); 
    
    for(int i = 0; i < n; i++){
      int x = m.nextInt(100); 
      int y = m.nextInt(100); 
      
      Vertex a1 = new Vertex(i, x, y); 
      addVertex(a1); 
    }
      
    LinkedList<Vertex> tmp = new LinkedList<>(getVertices()); 
    
      
    for(int i = 0; i < tmp.size(); i++){
        Vertex b = tmp.get(i); 
        
        for(int j = i + 1; j < tmp.size(); j++){
            Vertex c = tmp.get(j); 
            addUndirectedEdge(b.name, c.name, 100); 
        }
    }
    
    computeAllEuclideanDistances(); // compute distances
  }

  public List<Edge> nearestNeighborTsp() {
      
      //nearest neighbor algorithm
      
      LinkedList<Vertex> vertices = new LinkedList<>(getVertices()); 
      int size = vertices.size(); 
      int count = 0; 
      
      
      LinkedList<Edge> path = new LinkedList<>(); 
      Vertex origin = vertexNames.get(0); 
      Vertex start = vertexNames.get(0); 
      start.known = true; 
      
      while(count < size){
          
       
          List<Edge> tmp = new LinkedList<>(); 
          tmp = start.adjacentEdges; 
          
          Edge currentMin = null; 
              
          currentMin = getLowestDistance(tmp); 
          
          if(currentMin != null) {
              path.add(currentMin); 
          
              currentMin.target.known = true; 
              start = currentMin.target;
          }
          else{
              for (int i = 0; i < tmp.size(); i++){
                  if(tmp.get(i).target.equals(origin)){
                      path.add(tmp.get(i)); 
                  }
              }
          }
          count++; 
      }
      
      return path; 
  }
    
    public Edge getLowestDistance(List<Edge> x){
        
        Vertex origin = vertexNames.get(0); 
        Edge theEdge = null; 
        double min = Double.MAX_VALUE; 
        for(int i = 0; i < x.size(); i++){
            double edgeCost = x.get(i).distance;
            
            if(!x.get(i).target.known && !x.get(i).target.equals(origin)){
                if(edgeCost < min){
                    min = edgeCost; 
                    theEdge = x.get(i);
                }
            }
        }
        return theEdge; 
    }
    

    
  
  public List<Edge> bruteForceTsp() {
      
      //brute force algorithm
       List<Vertex> vertices = new ArrayList<>(getVertices()); 
       List<String> tours = new ArrayList<>(); 
       
       tours = permutate(vertices.size()); 
      
      //Let's find the starting route which is the best route by default
      
      String bestTour = tours.get(0); 
      double bestCost = 0.0; 
      List<Edge> bestPath = new LinkedList<>(); 
      int count = 0; 
      
      for(int i = 0; i < bestTour.length()-1; i++){
          
          int firstKey = Character.getNumericValue(bestTour.charAt(i)); 
          int endKey = Character.getNumericValue(bestTour.charAt(i+1)); 
          
          Vertex source = vertexNames.get(firstKey); 
          Vertex destination = vertexNames.get(endKey); 
          
          List<Edge> adj = new ArrayList<>(); 
          
          if(source != null){
              adj = source.adjacentEdges; 
          }
          
          for(int j = 0; j < adj.size(); j++){
              if(adj.get(j).target.equals(destination)){
                  bestCost = bestCost + adj.get(j).distance; 
                  bestPath.add(adj.get(j)); 
              } 
          } 
      }
      
      tours.remove(bestTour); 
      
      while(tours.size() != 0){
          
          String currentTour = tours.get(0); 
          double cost = 0; 
          List<Edge> path = new LinkedList<>(); 
          
          for(int i = 0; i < currentTour.length(); i++){
              
          
          if(i == currentTour.length() - 1){
              int firstKey = Character.getNumericValue(currentTour.charAt(i)); 
              int endKey = Character.getNumericValue(currentTour.charAt(0)); 
              
              Vertex source = vertexNames.get(firstKey); 
              Vertex destination = vertexNames.get(endKey); 
              
              List<Edge> adj = new ArrayList<>(source.adjacentEdges); 
              
              for(int j = 0; j < adj.size(); j++){
                  if(adj.get(j).target.equals(destination)){
                      cost = cost + adj.get(j).distance; 
                      path.add(adj.get(j)); 
                  }
              }
          }
          else{
              int firstKey = Character.getNumericValue(currentTour.charAt(i)); 
          int endKey = Character.getNumericValue(currentTour.charAt(i+1)); 
          
          Vertex source = vertexNames.get(firstKey); 
          Vertex destination = vertexNames.get(endKey); 
          
          List<Edge> adj = new ArrayList<>(); 
           
          if(source != null){
              adj = source.adjacentEdges; 
          }
          
          for(int j = 0; j < adj.size(); j++){
              if(adj.get(j).target.equals(destination)){
                  cost = cost + adj.get(j).distance; 
                  path.add(adj.get(j)); 
              } 
           }
          } 
      }
          if(cost < bestCost){
              bestCost = cost; 
              bestPath = path; 
          }
          tours.remove(currentTour); 
      }
      
      return bestPath; 
     
      
  }
    
    public static List<String> permutate (int x){
        //permutations methods to use during brute force
        
        String y = generateString(x-1); 
        
        List<String> permutations = new ArrayList<String>(); 
        
        permutate(" ", y, permutations);
        
        return permutations; 
        
    }
        private static void permutate (String first, String mot, List<String> thePermutations){
            
            LinkedList<String> permutations = new LinkedList<>(); 
            if(mot.isEmpty()){
                thePermutations.add(first); 
            }
            else{
                for(int i = 0; i < mot.length(); i++){
                    permutate(first + mot.charAt(i), mot.substring(0, i) + mot.substring(i+1, mot.length()), thePermutations); 
                }
            }

        }
    
    public static String generateString(int num){
        
        if(num == 0){
            return Integer.toString(num); 
        }
        
        return generateString(num - 1) + Integer.toString(num); 
    }

  // STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
