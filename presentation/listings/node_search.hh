/// Abstraction over a data structure to search nodes efficiently in a contact detector.
/// Possible implementations could include an exhaustive search, a spatial grid, an octree, ...
class NodeSearcher {
public:
  virtual ~NodeSearcher() = 0;

  /// Count nodes close enough to some position, optionally excluding a specific node.
  /// Pass a negative number to exclude no node.
  virtual UInt countCloseNodes(const Vector<Real> & position, Real maxDistance, UInt excludedNode) = 0;

  /// Find node index closest to the position, optionally excluding a specific node.
  /// The distance to the closest node is stored in best_distance.
  /// Pass a negative number to exclude no node.
  virtual UInt findClosestNode(const Vector<Real> & position, Real * best_distance, UInt excludedNode) = 0;
};