#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

extern int gridWidth, gridHeight;
extern int ra_lowerleft_x, ra_lowerleft_y;
extern int ra_width, ra_height;
extern float alpha, beta, gammaVal, delta, via_cost, overflowPenalty;

struct Point {
  int x, y;
  bool operator==(const Point& other) const {
    return x == other.x && y == other.y;
  }
  bool operator!=(const Point& other) const {
    return x != other.x || y != other.y;
  }
};

struct PointHash {
  std::size_t operator()(const Point& p) const {
    return std::hash<int>()(p.x) ^ std::hash<int>()(p.y);
  }
};

struct BumpPair {
  Point source;
  Point target;
  int bumpId;
};

int gridWidth = 0, gridHeight = 0;
int ra_lowerleft_x = 0, ra_lowerleft_y = 0;
int ra_width = 0, ra_height = 0;
float alpha = 0.0;
float beta = 0.0;
float gammaVal = 0.0;
float delta = 0.0;
float via_cost = 0.0;
float overflowPenalty = 0.0;

struct Node {
  Point coord;
  float g = 0.0, h = 0.0;
  Point parent = {-1, -1};
  int chip = -1;
  int bumpId = -1;
  int gcellId = -1;
  int leftCapacity = 0;
  int bottomCapacity = 0;
  int leftUsed = 0;
  int bottomUsed = 0;
  std::string type = "M1";
  float M1Cost = 0.0;
  float M2Cost = 0.0;

  float GradualCost(
      const std::unordered_map<Point, Node, PointHash>& tempMap) const {
    if (parent == Point{-1, -1}) return gammaVal * M1Cost;

    const Node& parentNode = tempMap.at(parent);

    int dx = std::abs(coord.x - parent.x);
    int dy = std::abs(coord.y - parent.y);

    float wireInc = alpha * (dx + dy);

    int overflow = 0;
    if (dx > 0) {
      if (coord.x > parent.x) {
        int used = tempMap.at(coord).leftUsed + 1;
        int cap = tempMap.at(coord).leftCapacity;
        int ov = used - cap;
        overflow = (ov > 0) ? ov : 0;
      } else {
        int used = tempMap.at(parent).leftUsed + 1;
        int cap = tempMap.at(parent).leftCapacity;
        int ov = used - cap;
        overflow = (ov > 0) ? ov : 0;
      }
    } else if (dy > 0) {
      if (coord.y > parent.y) {
        int used = tempMap.at(coord).bottomUsed + 1;
        int cap = tempMap.at(coord).bottomCapacity;
        int ov = used - cap;
        overflow = (ov > 0) ? ov : 0;
      } else {
        int used = tempMap.at(parent).bottomUsed + 1;
        int cap = tempMap.at(parent).bottomCapacity;
        int ov = used - cap;
        overflow = (ov > 0) ? ov : 0;
      }
    }

    float OV_initial = overflow * overflowPenalty;
    float OV_increment = beta * OV_initial;

    bool isVia = (type != parentNode.type);
    float stepCost = 0.0f;

    if (isVia) {
      float parentCellCost =
          gammaVal *
          ((parentNode.type == "M1") ? parentNode.M1Cost : parentNode.M2Cost);

      float viaNodeCost =
          gammaVal * ((parentNode.M1Cost + parentNode.M2Cost) * 0.5f) +
          (delta * via_cost);
      float currentNodeCost = gammaVal * ((type == "M1") ? M1Cost : M2Cost);
      stepCost = wireInc + OV_increment - parentCellCost + viaNodeCost +
                 currentNodeCost;
    } else {
      float GCell_increment = gammaVal * ((type == "M1") ? M1Cost : M2Cost);
      stepCost = wireInc + OV_increment + GCell_increment;
    }

    return stepCost;
  }

  float HeuristicCost(const Point& goal) const {
    int dx = std::abs(coord.x - goal.x);
    int dy = std::abs(coord.y - goal.y);
    return alpha * (dx + dy);
  }
};

std::unordered_map<Point, Node, PointHash> nodeMap;

struct CompareNode {
  bool operator()(const Node& a, const Node& b) {
    return a.g + a.h > b.g + b.h;
  }
};

void readGridInfo(const std::string& gmpFile, const std::string& gclFile,
                  const std::string& cstFile,
                  std::unordered_map<Point, Node, PointHash>& nodeMap,
                  std::vector<BumpPair>& bumpPairs) {
  std::ifstream gmpReader(gmpFile);
  if (!gmpReader.is_open()) {
    std::cerr << "Failed to open GMP file: " << gmpFile << "\n";
    return;
  }

  std::string line;
  int chip_offset_x = 0, chip_offset_y = 0;
  int currentChip = 0;
  std::unordered_map<int, BumpPair> tempBumpPairs;

  while (std::getline(gmpReader, line)) {
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    if (line.empty()) continue;

    std::istringstream iss(line);
    std::string tag;
    iss >> tag;

    if (tag == ".ra") {
      if (!std::getline(gmpReader, line) || line.empty()) continue;
      std::istringstream raStream(line);
      if (!(raStream >> ra_lowerleft_x >> ra_lowerleft_y >> ra_width >>
            ra_height))
        continue;
    } else if (tag == ".g") {
      if (!std::getline(gmpReader, line) || line.empty()) continue;
      std::istringstream gStream(line);
      if (!(gStream >> gridWidth >> gridHeight)) continue;
    } else if (tag == ".c") {
      if (!(iss >> chip_offset_x >> chip_offset_y)) {
        if (!std::getline(gmpReader, line) || line.empty()) continue;
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        std::istringstream nextLineStream(line);
        if (!(nextLineStream >> chip_offset_x >> chip_offset_y)) continue;
      }
      chip_offset_x += ra_lowerleft_x;
      chip_offset_y += ra_lowerleft_y;
      currentChip++;
    } else if (tag == ".b") {
      while (std::getline(gmpReader, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) break;

        std::istringstream bStream(line);
        int bumpId = 0, bump_x = 0, bump_y = 0;
        if (!(bStream >> bumpId >> bump_x >> bump_y)) continue;

        int real_bump_x = chip_offset_x + bump_x;
        int real_bump_y = chip_offset_y + bump_y;
        Point coord = {real_bump_x, real_bump_y};

        if (nodeMap.find(coord) == nodeMap.end()) {
          nodeMap[coord] = Node();
          nodeMap[coord].coord = coord;
        }

        nodeMap[coord].chip = currentChip;
        nodeMap[coord].bumpId = bumpId;

        if (tempBumpPairs.find(bumpId) == tempBumpPairs.end()) {
          tempBumpPairs[bumpId] = BumpPair{coord, {-1, -1}, bumpId};
        } else {
          tempBumpPairs[bumpId].target = coord;
          bumpPairs.push_back(tempBumpPairs[bumpId]);
          tempBumpPairs.erase(bumpId);
        }
      }
    }
  }
  gmpReader.close();

  std::sort(
      bumpPairs.begin(), bumpPairs.end(),
      [](const BumpPair& a, const BumpPair& b) { return a.bumpId < b.bumpId; });

  std::ifstream gclReader(gclFile);
  if (!gclReader.is_open()) {
    std::cerr << "Failed to open GCL file: " << gclFile << "\n";
    return;
  }

  const int xGrids = ra_width / gridWidth;
  const int yGrids = ra_height / gridHeight;

  while (std::getline(gclReader, line)) {
    if (line.find(".ec") != std::string::npos) {
      break;
    }
  }

  std::vector<std::pair<int, int>> capacities;
  while (std::getline(gclReader, line)) {
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    int leftCap, bottomCap;
    if (!(iss >> leftCap >> bottomCap)) continue;
    capacities.push_back({leftCap, bottomCap});
  }
  gclReader.close();

  int currentGCell = 0;
  for (int y = 0; y < yGrids; y++) {
    for (int x = 0; x < xGrids; x++) {
      int realX = ra_lowerleft_x + (x * gridWidth);
      int realY = ra_lowerleft_y + (y * gridHeight);
      Point coord = {realX, realY};

      currentGCell++;
      auto& cap = capacities[currentGCell - 1];

      if (nodeMap.find(coord) == nodeMap.end()) {
        nodeMap[coord] = Node();
        nodeMap[coord].coord = coord;
      }

      nodeMap[coord].gcellId = currentGCell;
      nodeMap[coord].leftCapacity = cap.first;
      nodeMap[coord].bottomCapacity = cap.second;
    }
  }

  std::ifstream cstReader(cstFile);
  if (!cstReader.is_open()) {
    std::cerr << "Failed to open CST file: " << cstFile << "\n";
    return;
  }

  std::string currentLayer = "";
  int currentY = ra_lowerleft_y;

  while (std::getline(cstReader, line)) {
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    if (line.empty()) continue;

    std::istringstream iss(line);
    std::string tag;
    iss >> tag;

    if (tag == ".alpha") {
      iss >> alpha;
    } else if (tag == ".beta") {
      iss >> beta;
    } else if (tag == ".gamma") {
      float value;
      iss >> value;
      gammaVal = value;
    } else if (tag == ".delta") {
      iss >> delta;
      if (delta < 1.0) {
        delta = 2.0 * delta;
      }
    } else if (tag == ".v") {
      std::getline(cstReader, line);
      std::istringstream viaStream(line);
      viaStream >> via_cost;
    } else if (tag == ".l") {
      if (currentLayer.empty()) {
        currentLayer = "M1";
      } else {
        currentLayer = "M2";
      }
      currentY = ra_lowerleft_y;
    } else {
      std::istringstream costStream(line);
      int currentX = ra_lowerleft_x;
      float cost;
      while (costStream >> cost) {
        Point coord = {currentX, currentY};
        if (nodeMap.find(coord) == nodeMap.end()) {
          nodeMap[coord] = Node();
          nodeMap[coord].coord = coord;
        }
        if (currentLayer == "M1") {
          nodeMap[coord].M1Cost = cost;
          overflowPenalty = std::max(overflowPenalty, cost);
        } else {
          nodeMap[coord].M2Cost = cost;
          overflowPenalty = std::max(overflowPenalty, cost);
        }
        currentX += gridWidth;
      }
      currentY += gridHeight;
    }
  }
  cstReader.close();
}

void outputPath(const std::deque<Node>& path, int bumpId,
                std::ofstream& outFile) {
  if (path.empty()) return;
  outFile << "n" << bumpId << "\n";

  size_t i = 0;
  Point currentPos = path[0].coord;
  bool lastVia = false;

  while (i < path.size() - 1) {
    const Node& current = path[i];
    size_t j = i;

    while (j < path.size() - 1 && path[j + 1].type == current.type) {
      j++;
    }

    if (current.coord == currentPos && current.coord == path[j].coord) {
      if (!lastVia) {
        outFile << "via\n";
        lastVia = true;
      }
    } else {
      if (current.type == "M1") {
        if (currentPos != path[j].coord) {
          outFile << "M1 " << currentPos.x << " " << currentPos.y << " "
                  << path[j].coord.x << " " << path[j].coord.y << "\n";
          lastVia = false;
        }
      } else {
        if (currentPos != path[j].coord) {
          outFile << "M2 " << currentPos.x << " " << currentPos.y << " "
                  << path[j].coord.x << " " << path[j].coord.y << "\n";
          lastVia = false;
        }
      }
      if (j < path.size() - 1) {
        outFile << "via\n";
        lastVia = true;
      }
    }
    currentPos = path[j].coord;
    i = j + 1;
  }

  if (i < path.size()) {
    if (currentPos == path.back().coord) {
      if (!lastVia) {
        outFile << "via\n";
      }
    } else {
      outFile << path.back().type << " " << currentPos.x << " " << currentPos.y
              << " " << path.back().coord.x << " " << path.back().coord.y
              << "\n";
    }
  }

  outFile << ".end\n";
}

std::deque<Node> AstarSearch(const Point& start, const Point& goal) {
  auto tempMap = nodeMap;
  std::priority_queue<Node, std::vector<Node>, CompareNode> openSet;
  std::unordered_map<Point, bool, PointHash> closedSet;

  Node startNode;
  startNode.coord = start;
  startNode.parent = {-1, -1};
  startNode.type = "M1";
  startNode.M1Cost = tempMap[start].M1Cost;
  startNode.M2Cost = tempMap[start].M2Cost;
  startNode.leftCapacity = tempMap[start].leftCapacity;
  startNode.bottomCapacity = tempMap[start].bottomCapacity;
  startNode.leftUsed = tempMap[start].leftUsed;
  startNode.bottomUsed = tempMap[start].bottomUsed;
  startNode.g = startNode.GradualCost(tempMap);
  startNode.h = startNode.HeuristicCost(goal);

  openSet.push(startNode);

  while (!openSet.empty()) {
    Node current = openSet.top();
    openSet.pop();

    if (closedSet[current.coord]) {
      continue;
    }

    if (current.parent.x != -1 && current.parent.y != -1) {
      if (current.coord.x != current.parent.x) {
        if (current.coord.x > current.parent.x) {
          tempMap[current.coord].leftUsed++;
        } else {
          tempMap[current.parent].leftUsed++;
        }
      } else if (current.coord.y != current.parent.y) {
        if (current.coord.y > current.parent.y) {
          tempMap[current.coord].bottomUsed++;
        } else {
          tempMap[current.parent].bottomUsed++;
        }
      }
    }

    closedSet[current.coord] = true;
    tempMap[current.coord] = current;

    if (current.coord == goal) {
      std::deque<Node> path;
      Node pathNode = current;

      if (current.type == "M1") {
        while (!(pathNode.coord == start)) {
          path.push_front(pathNode);
          pathNode = tempMap[pathNode.parent];
        }
        path.push_front(pathNode);
      } else {
        Node m1End = pathNode;
        m1End.type = "M1";

        while (!(pathNode.coord == start)) {
          path.push_front(pathNode);
          pathNode = tempMap[pathNode.parent];
        }
        path.push_front(pathNode);
        path.push_back(m1End);
      }

      for (size_t i = 0; i < path.size() - 1; i++) {
        const Node& current = path[i];
        const Node& next = path[i + 1];

        if (next.coord.x != current.coord.x) {
          if (next.coord.x > current.coord.x) {
            nodeMap[next.coord].leftUsed++;
          } else {
            nodeMap[current.coord].leftUsed++;
          }
        } else if (next.coord.y != current.coord.y) {
          if (next.coord.y > current.coord.y) {
            nodeMap[next.coord].bottomUsed++;
          } else {
            nodeMap[current.coord].bottomUsed++;
          }
        }

        nodeMap[current.coord].type = current.type;
        nodeMap[current.coord].g = current.g;
        nodeMap[current.coord].parent = current.parent;

        if (i == path.size() - 2) {
          nodeMap[next.coord].type = next.type;
          nodeMap[next.coord].g = next.g;
          nodeMap[next.coord].parent = next.parent;
        }
      }

      return path;
    }

    for (const auto& dir : std::vector<std::pair<int, int>>{{-gridWidth, 0},
                                                            {gridWidth, 0},
                                                            {0, -gridHeight},
                                                            {0, gridHeight}}) {
      Point newCoord = {current.coord.x + dir.first,
                        current.coord.y + dir.second};

      if (newCoord.x < ra_lowerleft_x ||
          newCoord.x >= ra_lowerleft_x + ra_width ||
          newCoord.y < ra_lowerleft_y ||
          newCoord.y >= ra_lowerleft_y + ra_height || closedSet[newCoord]) {
        continue;
      }

      Node neighbor;
      neighbor.coord = newCoord;
      neighbor.parent = current.coord;
      neighbor.type = (dir.first != 0) ? "M2" : "M1";

      neighbor.M1Cost = tempMap[newCoord].M1Cost;
      neighbor.M2Cost = tempMap[newCoord].M2Cost;
      neighbor.leftCapacity = tempMap[newCoord].leftCapacity;
      neighbor.bottomCapacity = tempMap[newCoord].bottomCapacity;
      neighbor.leftUsed = tempMap[newCoord].leftUsed;
      neighbor.bottomUsed = tempMap[newCoord].bottomUsed;

      float gradualCost = neighbor.GradualCost(tempMap);
      neighbor.g = current.g + gradualCost;
      neighbor.h = neighbor.HeuristicCost(goal);

      openSet.push(neighbor);
    }
  }

  std::cout << "No path found!\n";
  return std::deque<Node>();
}

int calculateManhattanDistance(const Point& p1, const Point& p2) {
    return std::abs(p1.x - p2.x) + std::abs(p1.y - p2.y);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <gmp_file> <gcl_file> <cst_file> <lg_file>\n";
        return 1;
    }

    std::string gmpFile = argv[1];
    std::string gclFile = argv[2];
    std::string cstFile = argv[3];
    std::string outputFile = argv[4];

    std::vector<BumpPair> bumpPairs;
    readGridInfo(gmpFile, gclFile, cstFile, nodeMap, bumpPairs);

    std::sort(bumpPairs.begin(), bumpPairs.end(),
        [](const BumpPair& a, const BumpPair& b) {
            int distA = calculateManhattanDistance(a.source, a.target);
            int distB = calculateManhattanDistance(b.source, b.target);
            return distA > distB; 
    });

    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Failed to create output file: " << outputFile << "\n";
        return 1;
    }

    for (const auto& pair : bumpPairs) {


        std::deque<Node> path = AstarSearch(pair.source, pair.target);

        if (!path.empty()) {
            outputPath(path, pair.bumpId, outFile);
        } else {
            std::cout << "Failed to route Bump " << pair.bumpId << "\n";
        }
    }
  

    outFile.close();
    return 0;
}