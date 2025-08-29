#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <utility>
#include <fstream>
#include <algorithm>   // for std::find_if
#include <stdexcept>   // for std::runtime_error
#include <cmath>       // for std::ceil

struct DieSize {
    int lowerLeftX; 
    int lowerLeftY; 
    int upperRightX; 
    int upperRightY; 
};

DieSize globalDieSize;

struct Cell {
    std::string name;
    int x, y;
    int width;
    int height;
    bool fixed;
    int new_x = -1;
    int new_y = -1;
    int placerow_idx = -1;

    // bool operator<(const Cell& other) const {
    //     return x < other.x || (x == other.x && y < other.y);
    // }
};

struct SubRow {
    int x_start;    
    int x_end;     
    int freespace; 

    std::set<Cell> cells;

    bool operator<(const SubRow& other) const {
        return x_start < other.x_start;
    }

    // void addCell(const Cell& cell) {
    //     if (cell.x >= x_start && cell.x + cell.width <= x_end) {
    //         cells.insert(cell);
    //         freespace -= cell.width;
    //     }
    // }

    // void removeCell(const Cell& cell) {
    //     auto it = cells.find(cell);
    //     if (it != cells.end()) {
    //         freespace += it->width; 
    //         cells.erase(it);  
    //     }
    // }
};

struct FixedCellInfo {
    std::string name;    
    int x, y;               
    int new_x, new_y;       
    int width, height;      
    int placerow_idx = -1; 
};

struct PlacementRow {
    int startX;      
    int startY;      
    int siteWidth;   
    int siteHeight;  
    int totalNumOfSites;
    int row_start;   
    int row_end;      
    std::set<SubRow> subRows;
};

void parseLGFile(
    const std::string& filename,
    double& alpha, double& beta,
    std::vector<PlacementRow>& rows,
    std::vector<Cell>& cells) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "Alpha") {
            iss >> alpha;
        } else if (token == "Beta") {
            iss >> beta;
        } else if (token == "PlacementRows") {
            PlacementRow row;
            iss >> row.startX >> row.startY >> row.siteWidth >> row.siteHeight >> row.totalNumOfSites;

            row.row_start = row.startX;
            row.row_end = row.startX + row.siteWidth * row.totalNumOfSites;

            rows.push_back(row); 
        } else if (token == "DieSize") {
            iss >> globalDieSize.lowerLeftX >> globalDieSize.lowerLeftY 
                >> globalDieSize.upperRightX >> globalDieSize.upperRightY;
        } else {
            Cell cell;
            cell.name = token;
            std::string type;
            iss >> cell.x >> cell.y >> cell.width >> cell.height >> type;

            cell.fixed = true;
            cells.push_back(cell);
        }
    }

    file.close();
}

void initializeSubRows(
    std::vector<PlacementRow>& rows,
    const std::vector<Cell>& cells,
    std::vector<FixedCellInfo>& fixedCellDetails) {

    for (PlacementRow& row : rows) {
        SubRow initialSubRow;
        initialSubRow.x_start = row.row_start;
        initialSubRow.x_end = row.row_end;
        initialSubRow.freespace = initialSubRow.x_end - initialSubRow.x_start;
        row.subRows.insert(initialSubRow);
    }
    std::sort(rows.begin(), rows.end(), [](const PlacementRow& a, const PlacementRow& b) {
    return a.startY < b.startY;
    });

    for (const Cell& cell : cells) {
        if (cell.fixed) {
            FixedCellInfo fixedCell;
            fixedCell.name = cell.name;
            fixedCell.x = cell.x;
            fixedCell.y = cell.y;
            fixedCell.width = cell.width;
            fixedCell.height = cell.height;

            for (PlacementRow& row : rows) {
                int rowStartY = row.startY;
                int rowEndY = row.startY + row.siteHeight;

                if (fixedCell.y < rowEndY && fixedCell.y + fixedCell.height > rowStartY) {
                    fixedCell.placerow_idx = (row.startY - rows[0].startY) / row.siteHeight;

                    std::set<SubRow> newSubRows;
                    for (const SubRow& subRow : row.subRows) {
                        if (fixedCell.x + fixedCell.width <= subRow.x_start || fixedCell.x >= subRow.x_end) {
                            newSubRows.insert(subRow);
                        } else if (fixedCell.x > subRow.x_start && fixedCell.x + fixedCell.width < subRow.x_end) {
                            SubRow leftSubRow = {subRow.x_start, fixedCell.x, fixedCell.x - subRow.x_start};
                            SubRow rightSubRow = {fixedCell.x + fixedCell.width, subRow.x_end, subRow.x_end - (fixedCell.x + fixedCell.width)};
                            newSubRows.insert(leftSubRow);
                            newSubRows.insert(rightSubRow);
                        } else if (fixedCell.x == subRow.x_start && fixedCell.x + fixedCell.width < subRow.x_end) {
                            SubRow adjustedSubRow = {fixedCell.x + fixedCell.width, subRow.x_end, subRow.x_end - (fixedCell.x + fixedCell.width)};
                            newSubRows.insert(adjustedSubRow);
                        } else if (fixedCell.x > subRow.x_start && fixedCell.x + fixedCell.width == subRow.x_end) {
                            SubRow adjustedSubRow = {subRow.x_start, fixedCell.x, fixedCell.x - subRow.x_start};
                            newSubRows.insert(adjustedSubRow);
                        }
                    }
                    row.subRows = newSubRows;
                }
            }

            fixedCellDetails.push_back(fixedCell);
        }
    }
}


void removeCell(std::vector<PlacementRow>& rows, const Cell& cellToRemove) {
    for (PlacementRow& row : rows) {
        if (cellToRemove.y < row.startY || cellToRemove.y >= row.startY + row.siteHeight) {
            continue;
        }

        bool leftMerged = false, rightMerged = false;
        auto leftIt = row.subRows.end(), rightIt = row.subRows.end();

        for (auto it = row.subRows.begin(); it != row.subRows.end(); ++it) {
            const SubRow& subRow = *it;

            if (subRow.x_end == cellToRemove.x) {
                leftMerged = true;
                leftIt = it;
            }
            if (subRow.x_start == cellToRemove.x + cellToRemove.width) {
                rightMerged = true;
                rightIt = it;
            }
        }

        if (leftMerged && rightMerged) {
            SubRow mergedSubRow = *leftIt;
            mergedSubRow.x_end = rightIt->x_end;
            mergedSubRow.freespace = mergedSubRow.x_end - mergedSubRow.x_start;

            row.subRows.erase(leftIt);
            row.subRows.erase(rightIt);

            row.subRows.insert(mergedSubRow);
            
        } else if (leftMerged) {
            SubRow updatedSubRow = *leftIt;
            updatedSubRow.x_end = cellToRemove.x + cellToRemove.width;
            updatedSubRow.freespace = updatedSubRow.x_end - updatedSubRow.x_start;

            row.subRows.erase(leftIt);

            row.subRows.insert(updatedSubRow);
            
        } else if (rightMerged) {
            SubRow updatedSubRow = *rightIt;
            updatedSubRow.x_start = cellToRemove.x;
            updatedSubRow.freespace = updatedSubRow.x_end - updatedSubRow.x_start;

            row.subRows.erase(rightIt);

            row.subRows.insert(updatedSubRow);
            
        } else {
            SubRow newSubRow;
            newSubRow.x_start = cellToRemove.x;
            newSubRow.x_end = cellToRemove.x + cellToRemove.width;
            newSubRow.freespace = newSubRow.x_end - newSubRow.x_start;

            row.subRows.insert(newSubRow);
       }

        break;
    }
}

bool insertCell(std::vector<PlacementRow>& rows, std::vector<Cell>& cells, Cell& cellToInsert, std::vector<Cell>& movedCells) {
    int closestRowIdx = -1;
    int y = cellToInsert.y;
    int height = cellToInsert.height;
    int width = cellToInsert.width;
    int siteHeight = rows[0].siteHeight; // 假設所有行的 siteHeight 相同

    // 找到最接近的行
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].startY <= y && y < rows[i].startY + siteHeight) {
            closestRowIdx = static_cast<int>(i);
            break;
        }
    }

    // 如果 y 小於所有 rows 的 startY，使用第一行
    if (closestRowIdx == -1 && y < rows[0].startY) {
        closestRowIdx = 0;
    }

    // 如果 y + height 超過 DieSize 上邊界，進行調整
    if (y + height > globalDieSize.upperRightY) {
        while (closestRowIdx >= 0 &&
               rows[closestRowIdx].startY + height > globalDieSize.upperRightY) {
            --closestRowIdx;
        }
        if (closestRowIdx < 0) {
            std::cout << "Cell " << cellToInsert.name
                      << " cannot be placed: exceeds DieSize upper boundary.\n";
            return false;
        }
    }

    // 計算需要的行數
    int rowsNeeded = (height + siteHeight - 1) / siteHeight;

    // 從 closestRowIdx 開始查找
    int startRowIdx = (closestRowIdx - 1) != 0 ? closestRowIdx - 1 : 0;
    do {
        for (const SubRow& subRow : rows[closestRowIdx].subRows) {
            int x_start = subRow.x_start;
            int x_end = subRow.x_end;

            // 每次嘗試向右移動 width / 4 的位置
            int moveStep = std::ceil(static_cast<double>(width) / 10.0);

            while (x_start + width <= x_end) {
                bool canInsert = true;

                // 檢查垂直方向是否可以插入
                for (int r = closestRowIdx; r < closestRowIdx + rowsNeeded; ++r) {
                    if (static_cast<std::size_t>(r) >= rows.size()) { // 行超出範圍
                        canInsert = false;
                        break;
                    }

                    bool rowHasSpace = false;
                    for (const SubRow& rowSubRow : rows[r].subRows) {
                        if (rowSubRow.x_start <= x_start &&
                            x_start + width <= rowSubRow.x_end &&
                            rowSubRow.freespace >= width) {
                            rowHasSpace = true;
                            break;
                        }
                    }
                    if (!rowHasSpace) {
                        canInsert = false;
                        break;
                    }
                }

                if (canInsert) {
                    // 切割 SubRow
                    for (int r = closestRowIdx; r < closestRowIdx + rowsNeeded; ++r) {
                        PlacementRow& row = rows[r];
                        std::set<SubRow> newSubRows;
                        for (const SubRow& subRow : row.subRows) {
                            if (x_start + width <= subRow.x_start || x_start >= subRow.x_end) {
                                newSubRows.insert(subRow);
                            } else if (x_start == subRow.x_start && x_start + width == subRow.x_end) {
                                continue;
                            } else if (x_start > subRow.x_start && x_start + width < subRow.x_end) {
                                SubRow leftSubRow = {subRow.x_start, x_start, x_start - subRow.x_start};
                                SubRow rightSubRow = {x_start + width, subRow.x_end, subRow.x_end - (x_start + width)};
                                newSubRows.insert(leftSubRow);
                                newSubRows.insert(rightSubRow);
                            } else if (x_start <= subRow.x_start && x_start + width < subRow.x_end) {
                                SubRow adjustedSubRow = {x_start + width, subRow.x_end, subRow.x_end - (x_start + width)};
                                newSubRows.insert(adjustedSubRow);
                            } else if (x_start > subRow.x_start && x_start + width >= subRow.x_end) {
                                SubRow adjustedSubRow = {subRow.x_start, x_start, x_start - subRow.x_start};
                                newSubRows.insert(adjustedSubRow);
                            }
                        }
                        row.subRows = newSubRows;
                    }

                    // 插入元件
                    cellToInsert.x = x_start;
                    cellToInsert.y = rows[closestRowIdx].startY;
                    cellToInsert.fixed = true;
                    cells.push_back(cellToInsert);

                    std::cout << "Cell " << cellToInsert.name << " successfully placed at ["
                              << cellToInsert.x << ", " << cellToInsert.y << "]\n";
                    return true;
                }

                // 向右移動
                x_start += moveStep;
            }
        }

        // 移動到下一行
        closestRowIdx = (closestRowIdx + 1) % rows.size();

    } while (closestRowIdx != startRowIdx);

    // 如果無法插入，返回 false
    std::cout << "Cell " << cellToInsert.name << " cannot be placed.\n";
    return false;
}

void processBankingInstruction(
    const std::string& instruction,
    std::vector<PlacementRow>& rows,
    std::vector<Cell>& cells,
    const std::string& outputFileName) {
    
    std::istringstream iss(instruction);
    std::string token;

    // Parse "Banking_Cell:" prefix
    iss >> token;
    if (token != "Banking_Cell:") {
        std::cerr << "Invalid instruction: " << instruction << "\n";
        return;
    }

    // Parse FFs to remove
    std::vector<std::string> ffsToRemove;
    while (iss >> token && token != "-->") {
        ffsToRemove.push_back(token);
    }

    // Parse merged FF name and its coordinates
    std::string mergedFFName;
    int x, y, w, h;
    iss >> mergedFFName >> x >> y >> w >> h;

    // Remove specified FFs
    for (const auto& ffName : ffsToRemove) {
        auto it = std::find_if(cells.begin(), cells.end(), [&ffName](const Cell& c) {
            return c.name == ffName;
        });

        if (it != cells.end()) {
            std::cout << "Removing FF: " << it->name << "\n";
            removeCell(rows, *it);
            cells.erase(it);       
        } else {
            std::cerr << "Warning: FF " << ffName << " not found in placement!\n";
        }
    }

    // Create and insert merged FF
    Cell mergedCell;
    mergedCell.name = mergedFFName;
    mergedCell.x = x;
    mergedCell.y = y;
    mergedCell.width = w;
    mergedCell.height = h;
    mergedCell.fixed = true; 

    std::vector<Cell> movedCells;
    std::ofstream outFile(outputFileName, std::ios::app);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open output file: " << outputFileName << "\n";
        return;
    }

    if (insertCell(rows, cells, mergedCell, movedCells)) { 

        outFile << mergedCell.x << " " << mergedCell.y << "\n";

        outFile << movedCells.size() << "\n";

        if (!movedCells.empty()) {
            for (const auto& moved : movedCells) {
                outFile << moved.name << " " << moved.x << " " << moved.y << "\n";
            }
        }
    } else {
        std::cerr << "Failed to place merged cell: " << mergedCell.name << "\n";
        outFile << mergedCell.x << " " << mergedCell.y << "\n";
        outFile << "0\n"; 
    }

    outFile.close();
}

int main(int argc, char* argv[]) {
    if (argc < 4) { 
        std::cerr << "Usage: " << argv[0] << " <lg_file> <banking_file> <output_file>\n";
        return 1;
    }

    double alpha = 0, beta = 0;
    std::vector<PlacementRow> rows;
    std::vector<Cell> cells;
    std::vector<FixedCellInfo> fixedCellDetails;

    try {
        parseLGFile(argv[1], alpha, beta, rows, cells);

        initializeSubRows(rows, cells, fixedCellDetails);
        std::ifstream bankingFile(argv[2]);
        if (!bankingFile.is_open()) {
            throw std::runtime_error("Unable to open banking file");
        }
        std::string instruction;
        while (std::getline(bankingFile, instruction)) {
            processBankingInstruction(instruction, rows, cells, argv[3]); 
        }

        bankingFile.close();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
