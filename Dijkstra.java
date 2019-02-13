/**Author of the student code: 
 * Alex Brice Horimbere
 * abh2167
 * 
 * This program implements Dijkstra Algorithm
 **/


import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.lang.Math; 
import java.util.Iterator;
import java.util.ArrayList; 
import java.util.Comparator; 
import java.util.Collections; 


public class Dijkstra {

  // Keep a fast index to nodes in the map
  private Map<String, Vertex> vertexNames;

  /**
   * Construct an empty Dijkstra with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Dijkstra() {
    vertexNames = new HashMap<String, Vertex>();
  }

  /**
   * Adds a vertex to the dijkstra. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the dijkstra
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the dijkstra
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the dijkstra
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
  public void addEdge(String nameU, String nameV, Double cost) {
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
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(String nameU, String nameV, double cost) {
    addEdge(nameU, nameV, cost);
    addEdge(nameV, nameU, cost);
  }

  // STUDENT CODE STARTS HERE

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
        // TODO
        double a = Math.pow(vx-ux, 2); 
        double b = Math.pow(vy-uy, 2); 
      
        return Math.sqrt(a+b); 
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
        // TODO
        
      Collection<Vertex> vertices = getVertices(); 
      Iterator<Vertex> l = vertices.iterator();
      
      while(l.hasNext()){
          Vertex a = l.next(); 
          List<Edge> isAdjacent = new LinkedList<>(); ;
          isAdjacent = a.adjacentEdges; 
          Iterator<Edge> p = isAdjacent.iterator(); 
          while(p.hasNext()){
          
          Edge tmp = p.next(); 
          Vertex start = tmp.source; 
          Vertex end = tmp.target; 
          
          double u1 = (double) start.x; 
          double u2 = (double) end.x; 
          double v1 = (double) start.y; 
          double v2 = (double) end.y; 
          
          tmp.distance = computeEuclideanDistance(u1, v1, u2, v2); 
          }
      }
  }

  /**
   * Dijkstra's Algorithm. 
   * 
   * @param s
   *          (String) starting city name
   */
  public void doDijkstra(String s) {
        // TODO  
         
        ArrayList<Vertex> vertices = new ArrayList<Vertex>(getVertices()); 
      
        for(int i = 0; i < vertices.size(); i++){
            
            //iterate through the vertices
            //set distance of first element to 0
            //Set distance of other elts to infinity
            
            Vertex current = vertices.get(i);
            
            if(current.equals(getVertex(s))){
                current.distance = 0.0; 
            }
            else{
                current.distance = Double.POSITIVE_INFINITY;
                current.known = false; 
                current.prev = null; 
            }
        }
      
      
      
      while(vertices.size() != 0){
          
          //while there are still vertices to visit
          //Find the vertex with minimum distance
          
          //Visit its adjacent Edges
          //Decide which next vertex to visit
          
          
          Vertex currentMin = getMinDist(vertices); 
          
          //Set known to true
          currentMin.known = true;
          
          List<Edge> adjacencyList = currentMin.adjacentEdges; 
          
          for(int i = 0; i < adjacencyList.size(); i++){
              
              Edge tmp = adjacencyList.get(i); 
              Vertex nextCity = tmp.target; 
              Vertex currentCity = tmp.source; 
              
              if(!nextCity.known){
                  double cost = tmp.distance; 
                  
                  if(currentCity.distance + cost < nextCity.distance){
                      nextCity.distance = currentCity.distance + cost; 
                      nextCity.prev = currentCity; 
                  }
              }
          }
          //Remove currentMin from unknown vertices
          vertices.remove(currentMin); 
           
      }
      
  }
    
    public Vertex getMinDist(List<Vertex> x){
        
        //to help determine the vertex with min Distances
        
        Vertex min = null; 
        double minDist = Double.MAX_VALUE; 
        for(int i = 0; i < x.size(); i++){
            double dist = x.get(i).distance;
            
            if(dist < minDist){
                minDist = dist; 
                min = x.get(i); 
            }
        }
        return min; 
    }
            

  /**
   * Returns a list of edges for a path from city s to city t. This will be the
   * shortest path from s to t as prescribed by Dijkstra's algorithm
   * 
   * @param s
   *          (String) starting city name
   * @param t
   *          (String) ending city name
   * @return (List<Edge>) list of edges from s to t
   */
  public List<Edge> getDijkstraPath(String s, String t) {
      
      doDijkstra(s); 
      
      LinkedList<Vertex> vertexPath = new LinkedList<>();
      
      
      LinkedList<Edge> path = new LinkedList<>(); 
      
      Vertex destination = getVertex(t); 
      
      while(destination.prev != null){
          
          //add all vertexes starting from destination and descending
          
          vertexPath.addFirst(destination); 
          destination = destination.prev; 
      }
      
      vertexPath.addFirst(getVertex(s)); //add the first vertex to the path
      
      for(int i = 0; i < vertexPath.size()-1; i++){
          
          //Find all the edges 
          
          Vertex source = vertexPath.get(i); 
          Vertex nextCity = vertexPath.get(i+1); 
          
          List<Edge> tmp = source.adjacentEdges; 
          
          for(int j = 0; j < tmp.size(); j++){
              Edge x = tmp.get(j); 
              if(x.target.equals(nextCity)){
                  path.add(x); 
              }
          }
      }
      return path; //return the path of edges. 
 
  }

  // STUDENT CODE ENDS HERE

  /**
   * Prints out the adjacency list of the dijkstra for debugging
   */
  public void printAdjacencyList() {
    for (String u : vertexNames.keySet()) {
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


  /** 
   * A main method that illustrates how the GUI uses Dijkstra.java to 
   * read a map and represent it as a graph. 
   * You can modify this method to test your code on the command line. 
   */
  public static void main(String[] argv) throws IOException {
    String vertexFile = "cityxy.txt"; 
    String edgeFile = "citypairs.txt";

    Dijkstra dijkstra = new Dijkstra();
    String line;

    // Read in the vertices
    BufferedReader vertexFileBr = new BufferedReader(new FileReader(vertexFile));
    while ((line = vertexFileBr.readLine()) != null) {
      String[] parts = line.split(",");
      if (parts.length != 3) {
        vertexFileBr.close();
        throw new IOException("Invalid line in vertex file " + line);
      }
      String cityname = parts[0];
      int x = Integer.valueOf(parts[1]);
      int y = Integer.valueOf(parts[2]);
      Vertex vertex = new Vertex(cityname, x, y);
      dijkstra.addVertex(vertex);
    }
    vertexFileBr.close();

    BufferedReader edgeFileBr = new BufferedReader(new FileReader(edgeFile));
    while ((line = edgeFileBr.readLine()) != null) {
      String[] parts = line.split(",");
      if (parts.length != 3) {
        edgeFileBr.close();
        throw new IOException("Invalid line in edge file " + line);
      }
      dijkstra.addUndirectedEdge(parts[0], parts[1], Double.parseDouble(parts[2]));
    }
    edgeFileBr.close();

    // Compute distances. 
    // This is what happens when you click on the "Compute All Euclidean Distances" button.
    dijkstra.computeAllEuclideanDistances();
    
    // print out an adjacency list representation of the graph
    dijkstra.printAdjacencyList();

    // This is what happens when you click on the "Draw Dijkstra's Path" button.

    // In the GUI, these are set through the drop-down menus.
    String startCity = "SanFrancisco";
    String endCity = "Boston";

    // Get weighted shortest path between start and end city. 
    List<Edge> path = dijkstra.getDijkstraPath(startCity, endCity);
    
    System.out.print("Shortest path between "+startCity+" and "+endCity+": ");
    System.out.println(path);
  }

}
