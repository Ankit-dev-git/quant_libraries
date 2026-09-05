#include "DataContainerCreation.h"

void savedatacontainer(Grid& dc, const std::string& filename)
{
	std::ofstream file(filename);
	for (int i = 0; i < dc.rows; ++i) {
		for (int j = 0; j < dc.cols; ++j) {
			file << dc(i, j);
			if (j < dc.cols - 1)
				file << ",";
		}
		file << "\n";
	}
}
